#!/bin/bash
#Usage: make target CORE=[core] [options]
#Example targets:
#    ifort
#    gfortran
#    xlf
#    pgi
#    intel-xd2000 :: cptec/inpe ian xd2000
#Availabe Cores:
#    atmosphere
#    init_atmosphere
#    landice
#    ocean
#    seaice
#    sw
#    test
#Available Options:
#    DEBUG=true    - builds debug version. Default is optimized version.
#    USE_PAPI=true - builds version using PAPI for timers. Default is off.
#    TAU=true      - builds version using TAU hooks for profiling. Default is off.
#    AUTOCLEAN=true    - forces a clean of infrastructure prior to build new core.
#    GEN_F90=true  - Generates intermediate .f90 files through CPP, and builds with them.
#    TIMER_LIB=opt - Selects the timer library interface to be used for profiling the model. Options are:
#                    TIMER_LIB=native - Uses native built-in timers in MPAS
#                    TIMER_LIB=gptl - Uses gptl for the timer interface instead of the native interface
#                    TIMER_LIB=tau - Uses TAU for the timer interface instead of the native interface
#    OPENMP=true   - builds and links with OpenMP flags. Default is to not use OpenMP.
#    OPENACC=true  - builds and links with OpenACC flags. Default is to not use OpenACC.
#    USE_PIO2=true - links with the PIO 2 library. Default is to use the PIO 1.x library.
#    PRECISION=single - builds with default single-precision real kind. Default is to use double-precision.
#    SHAREDLIB=true - generate position-independent code suitable for use in a shared library. Default is false.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MONAN_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SCRIPTS_DIR="${MONAN_ROOT}/scripts"
SOURCE_DIR="${MONAN_ROOT}/sources/MONAN-Model_1.4.3-rc"
EXECS_DIR="${MONAN_ROOT}/execs"

echo "-=+*#@#*@#*+=--=+*#@#*@#*+=--=+*#@#*@#*+=--=+*#@#*@#*+=-"
echo ""
echo "Script Directory: ${SCRIPT_DIR}"
echo "MONAN Root Directory: ${MONAN_ROOT}"
echo "Scripts Directory: ${SCRIPTS_DIR}"
echo "Source Directory: ${SOURCE_DIR}"
echo "Executables Directory: ${EXECS_DIR}"
echo ""
echo "-=+*#@#*@#*+=--=+*#@#*@#*+=--=+*#@#*@#*+=--=+*#@#*@#*+=-"
sleep 4

cd "${SCRIPTS_DIR}"
. "${SCRIPTS_DIR}/setenv.bash"
cd "${SOURCE_DIR}"

rm -rf "${SOURCE_DIR}/default_inputs/" "${SOURCE_DIR}/src/core_atmosphere/physics/physics_wrf/files"
rm -f  "${SOURCE_DIR}/"stream_list.* "${SOURCE_DIR}/"streams.* "${SOURCE_DIR}/"namelist.*
rm -f  "${SOURCE_DIR}/"make*.output.atmosphere "${SOURCE_DIR}/"make*.output.init_atmosphere
rm -fr "${SOURCE_DIR}/src/core_atmosphere/inc" "${SOURCE_DIR}/src/core_init_atmosphere/inc"
DATE_TIME_NOW=$(date +"%Y%m%d%H%M%S")

export NETCDF=${NETCDFDIR}
export PNETCDF=${PNETCDFDIR}
# PIO is not necessary for version 8.* If PIO is empty, MPAS Will use SMIOL
#export PIO=${PIODIR}
export PIO=

MAKE_OUT_FILE="make_${DATE_TIME_NOW}_.output.atmosphere"

make clean CORE=atmosphere
make -j 8 intel-xd2000 CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
#make -j 8 intel-xd2000 CORE=atmosphere OPENMP=true USE_PIO2=true PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}

#make -j 8 intel-xd2000 CORE=atmosphere OPENMP=true USE_PIO2=false PRECISION=single OPTIMIZATION_LEVEL=O1 FFLAGS_OPT=-O1 CFLAGS_OPT=-O1 CXXFLAGS_OPT=-O1 2>&1 | tee ${MAKE_OUT_FILE}


#CR: TODO: put verify here if executable was created ok
mv "${SOURCE_DIR}/atmosphere_model" "${EXECS_DIR}"
mv "${SOURCE_DIR}/build_tables" "${EXECS_DIR}"
cp "${SOURCE_DIR}/VERSION.txt" "${EXECS_DIR}/MONAN-VERSION.txt"
cp "${SOURCE_DIR}/GF_ConvPar_nml" "${SCRIPTS_DIR}"
make clean CORE=atmosphere

MAKE_OUT_FILE="make_${DATE_TIME_NOW}_.output.init_atmosphere"

make clean CORE=init_atmosphere
make -j 8 intel2-xd2000 CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}
#make -j 8 intel2-xd2000 CORE=init_atmosphere OPENMP=true USE_PIO2=true PRECISION=single 2>&1 | tee ${MAKE_OUT_FILE}

#make -j 8 intel-xd2000 CORE=init_atmosphere OPENMP=true USE_PIO2=false PRECISION=single OPTIMIZATION_LEVEL=O1 FFLAGS_OPT=-O1 CFLAGS_OPT=-O1 CXXFLAGS_OPT=-O1 2>&1 | tee ${MAKE_OUT_FILE}


mv "${SOURCE_DIR}/init_atmosphere_model" "${EXECS_DIR}"
make clean CORE=init_atmosphere


if [ -s "${EXECS_DIR}/init_atmosphere_model" ] && [ -e "${EXECS_DIR}/atmosphere_model" ]; then
    echo ""
    echo -e "\033[1;32m==>\033[0m Files init_atmosphere_model and atmosphere_model generated Successfully in ${EXECS_DIR} !"
    echo
else
    echo -e "\033[1;31m==>\033[0m !!! An error occurred during build. Check output"
    exit -1
fi

