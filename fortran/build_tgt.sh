#!/bin/sh
# Builds xfoilad_tgt.dll (primal + tangent XFOIL) into ../xfoilad/lib
#   sh build_tgt.sh            -> build from generated_tgt/ (Tapenade-generated tangent code)
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"
[ -d /c/msys64/mingw64/bin ] && export PATH=/c/msys64/mingw64/bin:$PATH

OPT=${XFOILAD_OPT:-"-O3 -march=native"}     # see build.sh
B=build_tgt
mkdir -p $B ../xfoilad/lib
FFLAGS="$OPT -fPIC -fdefault-real-8 -ffree-line-length-none -fno-range-check -w -J $B"

OBJS=""
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*)
        LIB=xfoilad_tgt.dll
        gcc $OPT -c api/xad_fastmath.c -o $B/xad_fastmath.o      # ucrtbase math, see build.sh
        OBJS="$B/xad_fastmath.o"
        LDF="-shared -static -Wl,--wrap=log,--wrap=exp,--wrap=tanh,--wrap=atan2,--wrap=pow" ;;
    Darwin*) LIB=libxfoilad_tgt.dylib; LDF="-shared" ;;
    *)       LIB=libxfoilad_tgt.so; LDF="-shared -fPIC" ;;
esac
for f in generated_tgt/i_pindex_p generated_tgt/i_blpar_p generated_tgt/i_xfoil_d generated_tgt/i_xbl_d \
         generated_tgt/s_xbl_d generated_tgt/s_xfoil_d generated_tgt/m_userio_d generated_tgt/m_spline_d \
         generated_tgt/m_xsolve_d generated_tgt/m_xutils_d generated_tgt/m_xpanel_d generated_tgt/m_xblsys_d \
         generated_tgt/m_xbl_d generated_tgt/m_xgeom_d generated_tgt/m_xgdes_d generated_tgt/s_xoper_d \
         generated_tgt/m_xoper_AD_d generated_tgt/xad_core_d generated_tgt/xad_zero_tgt api/xad_api_tgt; do
    o=$B/$(basename $f).o
    echo "  FC $f.f90"
    gfortran $FFLAGS -c $f.f90 -o $o
    OBJS="$OBJS $o"
done

gfortran $LDF $OPT -o ../xfoilad/lib/$LIB $OBJS
echo "built ../xfoilad/lib/$LIB"
