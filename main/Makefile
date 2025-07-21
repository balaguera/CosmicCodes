##########################################################################################################################################################
# Copyright (c) 2013-2024 A Balaguera
##################################################################################
# Makefile to compile COSMICATLASS
##################################################################################
CXX = g++
CXXFLAGS = --std=c++17 -fopenmp -fPIC -g3 -unroll
#-Wno-deprecated -Wno-write-strings -Wextra -msse -msse2 -msse4.1 -mfpmath=sse 
# Use CFLAGS to point out warnings
#CXXFLAGS += -Wall -Wextra -pedantic -Wno-deprecated -Wno-write-strings -Wall -Wextra -Wunused-variable -Wsign-compare -Wmisleading-indentation
##################################################################################
# Use /net/vaina if the library is to be used in diode, diosa or dardo
#AUX_DIR = ../Lib/
AUX_DIR = ../external_libs/
##################################################################################
# math and gsl libraries
LIBS_GSL = -lgsl -lgslcblas -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lm
#########################################################################
PYTHON = /usr/include/python3.8/
#########################################################################
# fftw libraries
dir_FFTW =  /usr/local/lib
dir_FFTW_INC = /usr/local/include/
LIBS_FFTW = -L$(dir_FFTW) -lfftw3 -lfftw3f -lfftw3f_omp -lfftw3_omp 
#########################################################################
# CFITSIO libraries
dir_CFITSIO_LIB=$(AUX_DIR)cfitsio/
#########################################################################
# healpix libraries
dir_HEALPIX_INC=$(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/include/
dir_HEALPIX_LIB=$(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/lib/
LIBS_HEALPIX = $(LIBS_GSL) -L$(dir_CFITSIO_LIB) -L$(dir_HEALPIX_LIB) -lhealpix_cxx -lcxxsupport -lsharp -lc_utils  -lfftpack -lcfitsio -lgomp
#########################################################################
#CCFITS libraries3
dir_CCFITS_LIB = $(AUX_DIR)CCfits-2.6/lib/
dir_CCFITS_INC = $(AUX_DIR)CCfits-2.6/include/CCfits/ 
LIBS_CCFITS =  -L$(dir_CFITSIO_LIB)  -L$(dir_HEALPIX_LIB) -L$(dir_CCFITS_LIB) -lCCfits
#########################################################################
#CosmiCalc Library
LCOSMO_DIR = ../lib
LCOSMO = -Wl,-rpath,$(LCOSMO_DIR) -L$(LCOSMO_DIR) -lCosmiCalCpp 
#########################################################################
HEADER_DIR = ../headers/
#########################################################################
CXXFLAGS += -I$(dir_HEALPIX_INC) -I$(dir_FFTW_INC) -I$(HEADER_DIR)
#########################################################################
# These dependencies are compiled by make using the implicit variables and CXXFLAGS
# Here I did not write a rule for compiling every dependency. Make does it when asking these dependencies
.PHONY: clean
#########################################################################

#########################################################################
# COSMOLIB
OBJS_CL = ../codes/cosmolib.o  ../source/CosmoLib.o
SRCS_CL = $(OBJS_CL:.o=.cpp)
TARGET_CLib = cosmolib.exe
#########################################################################
# BMT
OBJS = ../codes/bmt.o 
SRCS = $(OBJS:.o=.cpp)
TARGET_BMT = bmt.exe
#########################################################################
# power_spectrum
OBJS_PS = ../codes/power_spectrum.o 
SRCS_PS = $(OBJS_PS:.o=.cpp)
TARGET_PS = power.exe
#########################################################################
# power_spectrum
OBJS_SB = ../codes/secbias.o 
SRCS_SB = $(OBJS_SB:.o=.cpp)
TARGET_SB = secbias.exe
#########################################################################
# power_spectrum
OBJS_GT = ../codes/galtools.o ../source/Galaxy.o
SRCS_GT = $(OBJS_GT:.o=.cpp)
TARGET_GT = gtools.exe
#########################################################################
#########################################################################
# power_spectrum
OBJS_HT = ../codes/analysis_halo_catalogs.o 
SRCS_HT = $(OBJS_HT:.o=.cpp)
TARGET_HT = htools.exe

#########################################################################
# power_spectrum
OBJS_CW = ../codes/cosmic_web.o 
SRCS_CW = $(OBJS_CW:.o=.cpp)
TARGET_CW = cosmic_web.exe


#########################################################################
#########################################################################
# par
OBJS_PA = ../codes/paramsr.o 
SRCS_PA = $(OBJS_PA:.o=.cpp)
TARGET_PA = par
#########################################################################

# make will execute this rule by default
all: $(TARGET_BMT) 
##################################################################################
bmt: $(TARGET_BMT)
power: $(TARGET_PS)
secbias: $(TARGET_SB)
gtools: $(TARGET_GT)
htools: $(TARGET_HT)
cosmicweb: $(TARGET_CW)
cosmolib: $(TARGET_CLib)
params: $(TARGET_PA)
#########################################################################
$(TARGET_CLib): $(OBJS_CL) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling CosmoLib"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_CL) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_CLib) $(LIBS_GSL)
	@echo "Done!"
#########################################################################
$(TARGET_CW): $(OBJS_CW) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling Cosmicweb"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_CW) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_CW) $(LIBS_GSL)
	@echo "Done!"

#########################################################################
$(TARGET_BMT): $(OBJS) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling BMT"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_BMT) $(LIBS_GSL)
	@echo "Done!"
#########################################################################
$(TARGET_PS): $(OBJS_PS) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling Power"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_PS) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_PS) $(LIBS_GSL)
	@echo "Done!"
#########################################################################
$(TARGET_SB): $(OBJS_SB) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling Power"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_SB) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_SB) $(LIBS_GSL)
	@echo "Done!"
#########################################################################
$(TARGET_GT): $(OBJS_GT) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling Gtools"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_GT)  $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_GT) $(LIBS_GSL)
	@echo "Done!"
#########################################################################
$(TARGET_HT): $(OBJS_HT) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling Htools"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_HT)  $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_HT) $(LIBS_GSL)
	@echo "Done!"

#########################################################################
$(TARGET_PA): $(OBJS_PA) 
	@echo "*******************************************"
	@echo "*******************************************"
	@echo "Compiling params"
	@date  
	@echo "*******************************************"
# Este orden es de llamado a las librerias es importantisimo:
	@echo "*******************************************"         
	$(CXX) $(CXXFLAGS) $(SRCS_PA)  $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) -o $(TARGET_PA) $(LIBS_GSL)
	@echo "Done!"

clean:
	@echo "Cleaning:"
	rm -f  core.* *o *exe *~ 
#########################################################################
help:
	@echo "·········································································································"
	@echo "·········································································································"
	@echo "CosmiCalc compiling options"
	@echo "-make bmt: creates executable for the Bias Mapping Technique"
	@echo "	  execute with ./bmt.exe -b parameter_file to run bmt tasks."
	@echo "	  execute with ./bmt.exe -d to check preprocessor directives."
	@echo "	  execute with ./bmt.exe -i parameter_file to check input parameters directives."	
	@echo "-make power: creates executable for the power spectrum and related operations"
	@echo "	  execute with ./power.exe -p parameter_file to measure power (see documentation)"
	@echo "	  execute with ./power.exe -w parameter_file to measure window matrix"
	@echo "	  execute with ./power.exe -g parameter_file to obtain GRF."
	@echo "	  execute with ./power.exe -s parameter_file for low pass filter."
	@echo "	  execute with ./power.exe -m parameter_file for marked correlation function"			
	@echo "-make secbias: creates executable for the analysis of secondary bias"
	@echo "	  execute with ./secbias.exe -s parameter_file for low pass filter"		
	@echo "-make gtools: creates executable for the analysis of galaxy samples "
	@echo "	  execute with ./gtools.exe parameter_file (this is different from otopns 1-4). This is still under merging stages: merged with class::Catalog"		
	@echo "-make htools: creates executable for the analysis of halo catalogs "
	@echo "	  execute with ./htools.exe -c parameter_file to analyze catalog" 
	@echo "	  execute with ./htools.exe -b parameter_file to assign individual bias to input catalog"
	@echo "	  execute with ./htools.exe -s parameter_file to measure bias as a function of properteis"
	@echo "	  execute with ./htools.exe -h parameter_file to get galaxy catalogs using hod"
	@echo "	  execute with ./htools.exe -m parameter_file to build a mock catalog based on simulation (snapshot) using some dNdz" 	
	@echo "-make cosmolib: Cosmology functions and statistical tools (bias, mass function, power spectrum)"
	@echo "	  execute with ./cosmolib.exe -c parameter_file to get theoretical outputs for different cosmological functions at fixed redshift"
	@echo "	  execute with ./cosmolib.exe -h parameter_file to get theoretical outputs for HOD predictions (under constuction)"
	@echo "	  execute with ./cosmolib.exe -z parameter_file to get theoretical outputs for redshift dependent quantities" 		
	@echo "-make cosmicweb: Cosmic web analysis" 
	@echo "	  execute with ./cosmicweb.exe -f parameter_file to get all information on cosmnic web from DM and/or tracer field"
	@echo "·········································································································"
	@echo "·········································································································"

