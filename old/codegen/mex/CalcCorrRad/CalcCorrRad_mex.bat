@echo off
set MATLAB=C:\PROGRA~1\MATLAB\R2016a
set MATLAB_ARCH=win64
set MATLAB_BIN="C:\Program Files\MATLAB\R2016a\bin"
set ENTRYPOINT=mexFunction
set OUTDIR=.\
set LIB_NAME=CalcCorrRad_mex
set MEX_NAME=CalcCorrRad_mex
set MEX_EXT=.mexw64
call "C:\PROGRA~1\MATLAB\R2016a\sys\lcc64\lcc64\mex\lcc64opts.bat"
echo # Make settings for CalcCorrRad > CalcCorrRad_mex.mki
echo COMPILER=%COMPILER%>> CalcCorrRad_mex.mki
echo COMPFLAGS=%COMPFLAGS%>> CalcCorrRad_mex.mki
echo OPTIMFLAGS=%OPTIMFLAGS%>> CalcCorrRad_mex.mki
echo DEBUGFLAGS=%DEBUGFLAGS%>> CalcCorrRad_mex.mki
echo LINKER=%LINKER%>> CalcCorrRad_mex.mki
echo LINKFLAGS=%LINKFLAGS%>> CalcCorrRad_mex.mki
echo LINKOPTIMFLAGS=%LINKOPTIMFLAGS%>> CalcCorrRad_mex.mki
echo LINKDEBUGFLAGS=%LINKDEBUGFLAGS%>> CalcCorrRad_mex.mki
echo MATLAB_ARCH=%MATLAB_ARCH%>> CalcCorrRad_mex.mki
echo BORLAND=%BORLAND%>> CalcCorrRad_mex.mki
echo OMPFLAGS= >> CalcCorrRad_mex.mki
echo OMPLINKFLAGS= >> CalcCorrRad_mex.mki
echo EMC_COMPILER=lcc64>> CalcCorrRad_mex.mki
echo EMC_CONFIG=debug>> CalcCorrRad_mex.mki
"C:\Program Files\MATLAB\R2016a\bin\win64\gmake" -B -f CalcCorrRad_mex.mk
