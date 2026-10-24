#!/bin/sh
# Builds xfoilad.dll (primal + adjoint XFOIL) into ../xfoilad/lib
# Requires gfortran/gcc (MSYS2 mingw64 on Windows, or any gfortran on Linux).
#   sh build.sh            -> build from generated/ (Tapenade-generated adjoint code)
# Optimization flags: XFOILAD_OPT (default -O3 -march=native: the library
# then runs on this processor type only; use XFOILAD_OPT=-O3 for a library
# that is copied to other machines).
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
[ -d /c/msys64/mingw64/bin ] && export PATH=/c/msys64/mingw64/bin:$PATH

OPT=${XFOILAD_OPT:-"-O3 -march=native"}
B=build
mkdir -p $B ../xfoilad/lib
FFLAGS="$OPT -fPIC -fdefault-real-8 -ffree-line-length-none -fno-range-check -w -J $B"

gcc $OPT -fPIC -c adfirstaid/adStack.c -o $B/adStack.o
gfortran $OPT -fPIC -fallow-argument-mismatch -w -c adfirstaid/adBuffer.f -o $B/adBuffer.o
OBJS="$B/adStack.o $B/adBuffer.o"

case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
        LIB=xfoilad.dll
        # math functions of ucrtbase.dll instead of the (slow) MinGW libm, see api/xad_fastmath.c
        gcc $OPT -c api/xad_fastmath.c -o $B/xad_fastmath.o
        OBJS="$OBJS $B/xad_fastmath.o"
        LDF="-shared -static -Wl,--wrap=log,--wrap=exp,--wrap=tanh,--wrap=atan2,--wrap=pow" ;;
    Darwin*) LIB=libxfoilad.dylib; LDF="-shared" ;;
    *)       LIB=libxfoilad.so; LDF="-shared -fPIC" ;;
esac

for f in handad/xad_hand_store generated/i_pindex_p generated/i_blpar_p generated/i_xfoil_b generated/i_xbl_b \
         generated/s_xbl_b generated/s_xfoil_b generated/m_userio_b generated/m_spline_b \
         generated/m_xsolve_b generated/m_xutils_b generated/m_xpanel_b generated/m_xblsys_b \
         generated/m_xbl_b generated/m_xgeom_b generated/m_xgdes_b generated/s_xoper_b \
         generated/m_xoper_AD_b generated/xad_core_b generated/xad_zero_adj api/xad_api; do
    o=$B/$(basename $f).o
    echo "  FC $f.f90"
    gfortran $FFLAGS -c $f.f90 -o $o
    OBJS="$OBJS $o"
done

gfortran $LDF $OPT -o ../xfoilad/lib/$LIB $OBJS
echo "built ../xfoilad/lib/$LIB"
