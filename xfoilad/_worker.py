"""Worker process of ``pool.WorkerPool``: runs ``problem.run_point`` tasks
with its own copy of the XFOIL library.  Started by the pool, not by users."""
import os
from multiprocessing.connection import Listener


def main():
    address = os.environ["XFOILAD_WORKER_ADDRESS"]
    authkey = bytes.fromhex(os.environ["XFOILAD_WORKER_AUTHKEY"])
    from .problem import run_point
    from .solver import XFoilAD
    solver = XFoilAD()
    with Listener(address, authkey=authkey) as listener:
        with listener.accept() as conn:
            while True:
                try:
                    task = conn.recv()
                except (EOFError, OSError):
                    break
                if task is None:
                    break
                try:
                    out = run_point(solver, task)
                except Exception as e:          # reported to the main process
                    out = dict(exception=f"{type(e).__name__}: {e}")
                conn.send(out)


if __name__ == "__main__":
    main()
