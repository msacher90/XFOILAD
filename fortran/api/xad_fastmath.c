/***********************************************************************
 *  XFOILAD - faster math functions on Windows (MinGW-w64 builds only)
 *
 *  The libm of MinGW-w64 (libmingwex) computes log, atan2, pow, exp, ...
 *  in x87 extended precision: log 23 ns, atan2 37 ns, pow 69 ns per call on
 *  an i7-1260P, against 4.4, 20 and 10 ns for the functions of the Windows
 *  runtime ucrtbase.dll (correctly rounded results in both cases, a few
 *  values differ in the last bit).  The library is linked with
 *  -Wl,--wrap=log,... (build.sh): every call goes to __wrap_<f>, which
 *  calls the ucrtbase function, or the MinGW one when ucrtbase.dll is not
 *  available.
 ***********************************************************************/
#include <windows.h>

typedef double (*fun1)(double);
typedef double (*fun2)(double, double);

double __real_log(double);
double __real_exp(double);
double __real_tanh(double);
double __real_atan2(double, double);
double __real_pow(double, double);

static fun1 p_log = __real_log, p_exp = __real_exp, p_tanh = __real_tanh;
static fun2 p_atan2 = __real_atan2, p_pow = __real_pow;
static int ready = 0;

static void load(void)
{
    HMODULE h = LoadLibraryA("ucrtbase.dll");
    if (h) {
        FARPROC f;
        if ((f = GetProcAddress(h, "log"))) p_log = (fun1)(void *)f;
        if ((f = GetProcAddress(h, "exp"))) p_exp = (fun1)(void *)f;
        if ((f = GetProcAddress(h, "tanh"))) p_tanh = (fun1)(void *)f;
        if ((f = GetProcAddress(h, "atan2"))) p_atan2 = (fun2)(void *)f;
        if ((f = GetProcAddress(h, "pow"))) p_pow = (fun2)(void *)f;
    }
    ready = 1;
}

double __wrap_log(double x) { if (!ready) load(); return p_log(x); }
double __wrap_exp(double x) { if (!ready) load(); return p_exp(x); }
double __wrap_tanh(double x) { if (!ready) load(); return p_tanh(x); }
double __wrap_atan2(double y, double x) { if (!ready) load(); return p_atan2(y, x); }
double __wrap_pow(double x, double y) { if (!ready) load(); return p_pow(x, y); }
