cmake_minimum_required(VERSION 3.16)

project(CosmiCalCpp LANGUAGES CXX)

# ------------------------------------------------------------------------------
# C++ standard & flags
# ------------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Equivalent to: -fopenmp -fPIC -g3 -funroll-loops
add_compile_options(-fopenmp -fPIC -g3 -funroll-loops)

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
set(PROJECT_ROOT ${CMAKE_SOURCE_DIR})

# External libs
set(AUX_DIR "${PROJECT_ROOT}/external_libs")
set(SOURCE_DIR "${PROJECT_ROOT}/src")
set(INCLUDE_DIR "${PROJECT_ROOT}/include")

# FFTW
#set(DIR_FFTW /usr/local/lib)
#set(DIR_FFTW_INC /usr/local/include)

# CFITSIO
set(DIR_CFITSIO_LIB "${AUX_DIR}/cfitsio")

# HEALPix
set(DIR_HEALPIX_INC
	"${AUX_DIR}/Healpix_3.31/src/cxx/generic_gcc/include"
)

set(DIR_HEALPIX_LIB
	"${AUX_DIR}/Healpix_3.31/src/cxx/generic_gcc/lib"
)

# ------------------------------------------------------------------------------
# Python (equivalent to python3-config)
# ------------------------------------------------------------------------------
find_package(Python3 COMPONENTS Development REQUIRED)

# ------------------------------------------------------------------------------
# OpenMP
# ------------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

# ------------------------------------------------------------------------------
# Sources
# ------------------------------------------------------------------------------
set(SRCS
    ${SOURCE_DIR}/NumericalMethods.cpp
    ${SOURCE_DIR}/FileOutput.cpp
    ${SOURCE_DIR}/ScreenOutput.cpp
    ${SOURCE_DIR}/Params.cpp
    ${SOURCE_DIR}/CosmologicalFunctions.cpp
    ${SOURCE_DIR}/CoordinateSystem.cpp
    ${SOURCE_DIR}/Miscelanious.cpp
    ${SOURCE_DIR}/FftwFunctions.cpp
    ${SOURCE_DIR}/Statistics.cpp
    ${SOURCE_DIR}/BiasFunctions.cpp
    ${SOURCE_DIR}/PowerSpectrumF.cpp
    ${SOURCE_DIR}/PowerSpectrumTH.cpp
    ${SOURCE_DIR}/CorrelationFunction.cpp
    ${SOURCE_DIR}/Astrophysics.cpp
    ${SOURCE_DIR}/HaloTools.cpp
    ${SOURCE_DIR}/GnuplotC.cpp
    ${SOURCE_DIR}/Cwclass.cpp
    ${SOURCE_DIR}/DensityProfiles.cpp
    ${SOURCE_DIR}/Hod.cpp
    ${SOURCE_DIR}/McmcFunctions.cpp
    ${SOURCE_DIR}/ScalingRelations.cpp
    ${SOURCE_DIR}/BiasMappingTechnique.cpp
    ${SOURCE_DIR}/CosmoLib.cpp
    ${SOURCE_DIR}/AngularPowerSpectrumF.cpp
    ${SOURCE_DIR}/Catalogue.cpp
)

# ------------------------------------------------------------------------------
# Create shared library
# ------------------------------------------------------------------------------
add_library(CosmiCalCpp SHARED ${SRCS})

# ------------------------------------------------------------------------------
# Include directories
# ------------------------------------------------------------------------------
target_include_directories(CosmiCalCpp PRIVATE
    "${INCLUDE_DIR}"
    "${DIR_HEALPIX_INC}"
    "${DIR_FFTW_INC}"
    "${Python3_INCLUDE_DIRS}"
)

# ------------------------------------------------------------------------------
# Link directories
# ------------------------------------------------------------------------------
target_link_directories(CosmiCalCpp PRIVATE
    ${DIR_FFTW}
    ${DIR_HEALPIX_LIB}
    ${DIR_CFITSIO_LIB}
)

# ------------------------------------------------------------------------------
# Libraries (ORDER MATTERS — same as Makefile)
# ------------------------------------------------------------------------------
target_link_libraries(CosmiCalCpp PRIVATE
    # Healpix + related
    healpix_cxx
    cxxsupport
    sharp
    c_utils
    fftpack
    cfitsio
    gomp

    # Python
    ${Python3_LIBRARIES}

    # FFTW
#    fftw3
#    fftw3f
#    fftw3f_omp
#    fftw3_omp

    # GSL + others (keep math LAST like your Makefile)
    gsl
    gslcblas
    openblas
    util
    boost_iostreams
    boost_system
    boost_filesystem
    m
)

# ------------------------------------------------------------------------------
# OpenMP (clean way)
# ------------------------------------------------------------------------------
target_link_libraries(CosmiCalCpp PRIVATE OpenMP::OpenMP_CXX)

# ------------------------------------------------------------------------------
# Output name (libCosmiCalCpp.so)
# ------------------------------------------------------------------------------
set_target_properties(CosmiCalCpp PROPERTIES
    OUTPUT_NAME "CosmiCalCpp"
)
