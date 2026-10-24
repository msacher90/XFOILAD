"""Pool of worker processes for the XFOIL evaluations.

The XFOIL library holds one global state, so parallel evaluations need
separate processes.  Each worker is a ``python -m xfoilad._worker``
subprocess with its own copy of the DLL; tasks and results go through a
``multiprocessing.connection`` channel (named pipe on Windows).  Using
plain subprocesses instead of ``multiprocessing`` means that user scripts
need no ``if __name__ == "__main__":`` guard, and a crash of the library
only kills one worker (the task is reported as failed and the worker is
restarted).
"""
import atexit
import os
import subprocess
import sys
import tempfile
import time
import uuid
from multiprocessing.connection import Client

_PKG_PARENT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _address():
    name = f"xfoilad-{os.getpid()}-{uuid.uuid4().hex[:12]}"
    if sys.platform.startswith("win"):
        return r"\\.\pipe" + "\\" + name
    return os.path.join(tempfile.gettempdir(), name + ".sock")


class _Worker:
    def __init__(self, startup_timeout=120.0):
        self.address = _address()
        self.authkey = os.urandom(16)
        env = dict(os.environ)
        env["PYTHONPATH"] = _PKG_PARENT + os.pathsep + env.get("PYTHONPATH", "")
        env["XFOILAD_WORKER_ADDRESS"] = self.address
        env["XFOILAD_WORKER_AUTHKEY"] = self.authkey.hex()
        self.proc = subprocess.Popen([sys.executable, "-m", "xfoilad._worker"], env=env,
                                     stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
        t0 = time.time()
        while True:                      # the worker creates the listener: connect when it is up
            try:
                self.conn = Client(self.address, authkey=self.authkey)
                break
            except OSError:
                if self.proc.poll() is not None:
                    raise RuntimeError(f"XFOILAD worker exited at startup (code {self.proc.returncode})")
                if time.time() - t0 > startup_timeout:
                    self.proc.kill()
                    raise RuntimeError("XFOILAD worker did not start")
                time.sleep(0.05)

    def close(self):
        try:
            self.conn.send(None)
            self.conn.close()
        except OSError:
            pass
        try:
            self.proc.wait(timeout=5)
        except subprocess.TimeoutExpired:
            self.proc.kill()


class WorkerPool:
    """Fixed-size pool; ``map`` runs a list of ``problem.run_point`` tasks."""

    def __init__(self, n_workers):
        self.n_max = max(1, int(n_workers))
        self.workers = []
        atexit.register(self.close)

    def _ensure(self, n):
        while len(self.workers) < min(n, self.n_max):
            self.workers.append(_Worker())

    def map(self, tasks):
        self._ensure(len(tasks))
        out = [None] * len(tasks)
        nw = len(self.workers)
        for start in range(0, len(tasks), nw):
            batch = list(range(start, min(start + nw, len(tasks))))
            sent = []
            for k, i in enumerate(batch):
                try:
                    self.workers[k].conn.send(tasks[i])
                    sent.append((k, i))
                except OSError:
                    out[i] = self._crashed(k, tasks[i])
            for k, i in sent:
                try:
                    res = self.workers[k].conn.recv()
                except (EOFError, OSError):
                    res = self._crashed(k, tasks[i])
                if "exception" in res:
                    raise RuntimeError(f"XFOILAD worker error: {res['exception']}")
                out[i] = res
        return out

    def _crashed(self, k, task):
        """The library crashed in worker k: restart it, report the task as failed."""
        try:
            self.workers[k].proc.kill()
        except OSError:
            pass
        self.workers[k] = _Worker()
        return dict(result=None, values={}, grads={}, continuation=False,
                    error=f"XFOIL crashed at operating point {task['point'].name}")

    def close(self):
        for w in self.workers:
            w.close()
        self.workers = []
