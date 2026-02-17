##########################################################################################################################################################
# Copyright (c) 2013-2024 Andres Balaguera Antolinez
##################################################################################
CXX = g++
CXXFLAGS = --std=c++17 -fopenmp -fPIC -g3 -funroll-loops
DEBUGG   =
##################################################################################
AUX_DIR = ../external_libs/
# math and gsl libraries
PYTHON_INC := $(shell python3-config --embed --includes)
PYTHON_LIB := $(shell python3-config --embed --ldflags)
LIBS_GSL = -lgsl -lgslcblas -lopenblas -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lm $(PYTHON_LIB)
##################################################################################
# fftw libraries
DIR_FFTW = /usr/local/lib
DIR_FFTW_INC =  /usr/local/include/
LIBS_FFTW = -L$(DIR_FFTW) -lfftw3 -lfftw3f -lfftw3f_omp -lfftw3_omp 
##################################################################################
DIR_CFITSIO_LIB =  $(AUX_DIR)cfitsio/
##################################################################################
DIR_HEALPIX_INC :=  $(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/include
DIR_HEALPIX_LIB :=  $(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/lib/
LIBS_HEALPIX    := -L$(DIR_CFITSIO_LIB) -L$(DIR_HEALPIX_LIB) -lhealpix_cxx -lcxxsupport -lsharp -lc_utils  -lfftpack -lcfitsio -lgomp 
##################################################################################
# -Wall -Wextra to check all warnings
##################################################################################
#C++ Source files
SOURCE_DIR = ../source/
SRCS =  $(SOURCE_DIR)NumericalMethods.cpp $(SOURCE_DIR)FileOutput.cpp  $(SOURCE_DIR)ScreenOutput.cpp $(SOURCE_DIR)Params.cpp $(SOURCE_DIR)CosmologicalFunctions.cpp\
	$(SOURCE_DIR)CoordinateSystem.cpp $(SOURCE_DIR)Miscelanious.cpp $(SOURCE_DIR)FftwFunctions.cpp $(SOURCE_DIR)Statistics.cpp $(SOURCE_DIR)BiasFunctions.cpp\
	$(SOURCE_DIR)PowerSpectrumF.cpp $(SOURCE_DIR)PowerSpectrumTH.cpp  $(SOURCE_DIR)CorrelationFunction.cpp $(SOURCE_DIR)Astrophysics.cpp $(SOURCE_DIR)Catalog.cpp $(SOURCE_DIR)GnuplotC.cpp\
	$(SOURCE_DIR)Cwclass.cpp $(SOURCE_DIR)DensityProfiles.cpp $(SOURCE_DIR)Hod.cpp $(SOURCE_DIR)McmcFunctions.cpp $(SOURCE_DIR)ScalingRelations.cpp $(SOURCE_DIR)BiasMappingTechnique.cpp\
	$(SOURCE_DIR)CosmoLib.cpp 
#$(SOURCE_DIR)Galaxy.cpp
##################################################################################
#C++ Header files
HEADER_DIR = ../headers/
HDRS =  $(HEADER_DIR)CosmiCalcLIB.h  $(HEADER_DIR)Constants.h $(HEADER_DIR)constants.h $(HEADER_DIR)cosmo_parameters.h $(HEADER_DIR)def.h $(HEADER_DIR)Structures.h\
	$(HEADER_DIR)NumericalMethods.h  $(HEADER_DIR)FileOutput.h $(HEADER_DIR)ScreenOutput.h $(HEADER_DIR)ScalingRelations.h $(HEADER_DIR)Astrophysics.h\
	$(HEADER_DIR)Params.h $(HEADER_DIR)CosmologicalFunctions.h $(HEADER_DIR)FftwFunctions.h $(HEADER_DIR)CoordinateSystem.h\
	$(HEADER_DIR)Miscelanious.h  $(HEADER_DIR)lss_vector.h $(HEADER_DIR)FftwFunctions.h $(HEADER_DIR)BiasFunctions.h $(HEADER_DIR)Statistics.h\
	$(HEADER_DIR)PowerSpectrumTH.h $(HEADER_DIR)CorrelationFunction.h $(HEADER_DIR)PowerSpectrumF.h $(HEADER_DIR)Catalog.h $(HEADER_DIR)GnuplotC.h $(HEADER_DIR)Cwclass.h $(HEADER_DIR)CosmoLib.h\
	$(HEADER_DIR)DensityProfiles.h $(HEADER_DIR)Hod.h $(HEADER_DIR)McmcFunctions.h $(HEADER_DIR)BiasMappingTechnique.h	
	#$(HEADER_DIR)Galaxy.h
#$(HEADER_DIR)ClFunctions.h
##################################################################################
OBJS = $(SRCS:.cpp=.o)

INCLUDE :=  -I$(DIR_HEALPIX_INC) -I$(HEADER_DIR) -I$(DIR_FFTW_INC) $(PYTHON_INC)
##################################################################################
#Set library name
TARGET = libCosmiCalCpp.so
##################################################################################

all: $(TARGET)

.PHONY: clean

$(TARGET): $(OBJS)
# Este orden es de llamado a las librerias es importantisimo: GSL , con su -lm, debe quedar de ultimo
	$(CXX) $(CXXFLAGS) -shared  $(LIBS_HEALPIX) $(PYTHON_LIB) $(LIBS_FFTW) $(OBJS)  $(LIBS_GSL) -o $(TARGET)
	@echo  Successfully compiled the CosmiCalCpp.so library.

clean:
	@echo "Cleaning:"
	rm -f  core.* *o *exe *~ ../source/*o  ../source/*~  ../headers/*~
##################################################################################

#Esto evita tener que escribir 30 veces la compilación de todos los .cpp en SOURCE_DIR
$(SOURCE_DIR)%.o: $(SOURCE_DIR)%.cpp $(HDRS)
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@


#$(SOURCE_DIR)ClFunctions.o: $(SOURCE_DIR)ClFunctions.cpp $(HDRS)
#	@echo "Compiling ClFunctions.cpp"
#	$(CXX) $(CXXFLAGS) $(INCLUDE)  -c $< -o $@ $(LIBS_GSL)

#$(SOURCE_DIR)Galaxy.o: $(SOURCE_DIR)Galaxy.cpp #
#	@echo "Compiling Galaxy.cpp"
#	$(CXX) $(CXXFLAGS) $(INCLUDE)   $(LIBS_HEALPIX)  $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)


