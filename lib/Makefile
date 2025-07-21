##########################################################################################################################################################
# Copyright (c) 2013-2024 Andres Balaguera Antolinez
##################################################################################
CXX = gcc
CXXFLAGS = --std=c++17 -fopenmp -fPIC -g3 -unroll 
DEBUGG   =
##################################################################################
AUX_DIR = ../external_libs/
# math and gsl libraries
LIBS_GSL = -lgsl -lgslcblas -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lm
PYTHON = /usr/include/python3.8/
##################################################################################
# fftw libraries
DIR_FFTW = /usr/local/lib
DIR_FFTW_INC =  /usr/local/include/
LIBS_FFTW = -L$(DIR_FFTW) -lfftw3 -lfftw3f -lfftw3f_omp -lfftw3_omp 
##################################################################################
DIR_CFITSIO_LIB =  $(AUX_DIR)cfitsio/
##################################################################################
DIR_HEALPIX_INC =  $(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/include/
DIR_HEALPIX_LIB =  $(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/lib/
LIBS_HEALPIX    = -L$(DIR_CFITSIO_LIB) -L$(DIR_HEALPIX_LIB) -lhealpix_cxx -lcxxsupport -lsharp -lc_utils  -lfftpack -lcfitsio -lgomp 
CXXFLAGS  += -I$(DIR_HEALPIX_INC) 
##################################################################################
# -Wall -Wextra to check all warnings
##################################################################################
#C++ Source files
SOURCE_DIR = ../source/
SRCS =  $(SOURCE_DIR)NumericalMethods.cpp $(SOURCE_DIR)FileOutput.cpp  $(SOURCE_DIR)ScreenOutput.cpp $(SOURCE_DIR)Params.cpp $(SOURCE_DIR)CosmologicalFunctions.cpp\
	$(SOURCE_DIR)CoordinateSystem.cpp $(SOURCE_DIR)Miscelanious.cpp $(SOURCE_DIR)FftwFunctions.cpp $(SOURCE_DIR)Statistics.cpp $(SOURCE_DIR)BiasFunctions.cpp\
	$(SOURCE_DIR)PowerSpectrumF.cpp $(SOURCE_DIR)PowerSpectrumTH.cpp $(SOURCE_DIR)Astrophysics.cpp $(SOURCE_DIR)Catalog.cpp $(SOURCE_DIR)GnuplotC.cpp\
	$(SOURCE_DIR)Cwclass.cpp $(SOURCE_DIR)DensityProfiles.cpp $(SOURCE_DIR)Hod.cpp $(SOURCE_DIR)McmcFunctions.cpp $(SOURCE_DIR)ScalingRelations.cpp $(SOURCE_DIR)BiasMappingTechnique.cpp
#(SOURCE_DIR)ClFunctions.cpp $(SOURCE_DIR)Galaxy.cpp
##################################################################################
#C++ Header files
HEADER_DIR = ../headers/
HDRS =  $(HEADER_DIR)CosmiCalcLIB.h  $(HEADER_DIR)Constants.h $(HEADER_DIR)constants.h $(HEADER_DIR)cosmo_parameters.h $(HEADER_DIR)def.h $(HEADER_DIR)Structures.h\
	$(HEADER_DIR)NumericalMethods.h  $(HEADER_DIR)FileOutput.h $(HEADER_DIR)ScreenOutput.h $(HEADER_DIR)ScalingRelations.h $(HEADER_DIR)Astrophysics.h\
	$(HEADER_DIR)Params.h $(HEADER_DIR)CosmologicalFunctions.h $(HEADER_DIR)FftwFunctions.h $(HEADER_DIR)CoordinateSystem.h\
	$(HEADER_DIR)Miscelanious.h  $(HEADER_DIR)lss_vector.h $(HEADER_DIR)FftwFunctions.h $(HEADER_DIR)BiasFunctions.h $(HEADER_DIR)Statistics.h\
	$(HEADER_DIR)PowerSpectrumTH.h $(HEADER_DIR)PowerSpectrumF.h $(HEADER_DIR)Catalog.h $(HEADER_DIR)GnuplotC.h $(HEADER_DIR)Cwclass.h $(HEADER_DIR)CosmoLib.h\
	$(HEADER_DIR)DensityProfiles.h $(HEADER_DIR)Hod.h $(HEADER_DIR)McmcFunctions.h  $(HEADER_DIR)BiasMappingTechnique.h
#$(HEADER_DIR)ClFunctions.h
##################################################################################
OBJS = $(SRCS:.cpp=.o)

INCLUDE =  -I$(HEADER_DIR) -I$(DIR_FFTW_INC)
##################################################################################
#Set library name
TARGET = libCosmiCalCpp.so
##################################################################################

all: $(TARGET)

.PHONY: clean

$(TARGET): $(OBJS)
# Este orden es de llamado a las librerias es importantisimo: GSL , con su -lm, debe quedar de ultimo
	$(CXX) $(CXXFLAGS) -shared  $(LIBS_HEALPIX) $(LIBS_FFTW) $(OBJS)  $(LIBS_GSL) -o $(TARGET)
	@echo  Successfully compiled the CosmiCalCpp.so library.

clean:
	@echo "Cleaning:"
	rm -f  core.* *o *exe *~ ../source/*o  ../source/*~  ../headers/*~
##################################################################################
$(SOURCE_DIR)NumericalMethods.o: $(SOURCE_DIR)NumericalMethods.cpp $(HDRS) 
	@echo "Compiling NumericalMethods.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_HEALPIX) $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)FileOutput.o: $(SOURCE_DIR)FileOutput.cpp $(HDRS) 
	@echo "Compiling FileOutput.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)ScreenOutput.o: $(SOURCE_DIR)ScreenOutput.cpp $(HDRS)  
	@echo "Compiling ScreenOutput.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_GSL) $(LIBS_HEALPIX) $(LIBS_FFTW) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)Params.o: $(SOURCE_DIR)Params.cpp $(HDRS)
	@echo "Compiling Params.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS_GSL)


$(SOURCE_DIR)Miscelanious.o: $(SOURCE_DIR)Miscelanious.cpp $(HDRS)  
	@echo "Compiling Miscelanious.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_HEALPIX) $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)ScalingRelations.o: $(SOURCE_DIR)ScalingRelations.cpp $(HDRS)
	@echo "Compiling ScalingRelations.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_HEALPIX) $(LIBS_FFTW) -c $< -o $@  $(LIBS_GSL)

$(SOURCE_DIR)CosmologicalFunctions.o: $(SOURCE_DIR)CosmologicalFunctions.cpp $(HDRS) 
	@echo "Compiling CosmologicalFunctions.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@  $(LIBS_GSL)

$(SOURCE_DIR)CoordinateSystem.o: $(SOURCE_DIR)CoordinateSystem.cpp $(HDRS)
	@echo "Compiling CoordinateSystem.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  $(LIBS_GSL) $(LIBS_HEALPIX) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)Astrophysics.o: $(SOURCE_DIR)Astrophysics.cpp $(HDRS)
	@echo "Compiling Astrophysics.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS_GSL)


$(SOURCE_DIR)FftwFunctions.o: $(SOURCE_DIR)FftwFunctions.cpp $(HDRS)
	@echo "Compiling FFtwFunctions.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_FFTW) $(LIBS_HEALPIX) -c $< -o $@  $(LIBS_GSL)

$(SOURCE_DIR)Statistics.o: $(SOURCE_DIR)Statistics.cpp $(HDRS)
	@echo "Compiling Statistics.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  $(LIBS_HEALPIX)  $(LIBS_FFTW) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)BiasFunctions.o: $(SOURCE_DIR)BiasFunctions.cpp $(HDRS)
	@echo "Compiling BiasFunctions.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_FFTW) -c $< -o $@

$(SOURCE_DIR)PowerSpectrumTH.o: $(SOURCE_DIR)PowerSpectrumTH.cpp $(HDRS)
	@echo "Compiling PowerSpectrumTH.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  $(LIBS_HEALPIX)  $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)PowerSpectrumF.o: $(SOURCE_DIR)PowerSpectrumF.cpp $(HDRS)
	@echo "Compiling PowerSpectrumF.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  $(LIBS_HEALPIX)  $(LIBS_FFTW) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)Catalog.o: $(SOURCE_DIR)Catalog.cpp $(HDRS)
	@echo "Compiling Catalog.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_HEALPIX)  $(LIBS_FFTW) -c $< -o $@ $(LIBS_GSL) 

$(SOURCE_DIR)Cwclass.o: $(SOURCE_DIR)Cwclass.cpp $(HDRS)
	@echo "Compiling Cwclass.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)GnuplotC.o: $(SOURCE_DIR)GnuplotC.cpp $(HDRS)
	@echo "Compiling GnuplotC.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)DensityProfiles.o: $(SOURCE_DIR)DensityProfiles.cpp $(HDRS)
	@echo "Compiling DensityProfiles.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)Hod.o: $(SOURCE_DIR)Hod.cpp $(HDRS)
	@echo "Compiling Hod.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@ $(LIBS_GSL)

$(SOURCE_DIR)McmcFunctions.o: $(SOURCE_DIR)McmcFunctions.cpp $(HDRS)
	@echo "Compiling McmcFunctions.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE)  -c $< -o $@ $(LIBS_GSL)


$(SOURCE_DIR)BiasMappingTechnique.o: $(SOURCE_DIR)BiasMappingTechnique.cpp $(HDRS)
	@echo "Compiling BMT.cpp"
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS_HEALPIX)  $(LIBS_FFTW) -c $< -o $@ $(LIBS_GSL) 


#$(SOURCE_DIR)ClFunctions.o: $(SOURCE_DIR)ClFunctions.cpp $(HDRS)
#	@echo "Compiling ClFunctions.cpp"
#	$(CXX) $(CXXFLAGS) $(INCLUDE)  -c $< -o $@ $(LIBS_GSL)

#$(SOURCE_DIR)Galaxy.o: $(SOURCE_DIR)Galaxy.cpp 
#	@echo "Compiling Galaxy.cpp"
#	$(CXX) $(CXXFLAGS) $(INCLUDE)   $(LIBS_HEALPIX)  $(LIBS_FFTW)  -c $< -o $@ $(LIBS_GSL)


