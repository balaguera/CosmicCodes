##########################################################################################################################################################
# Copyright (c) 2013-2024 A Balaguera
GREEN=\033[0;32m
BLUE=\033[0;34m
YELLOW=\033[0;33m
RED=\033[0;31m

##########################################################################
# COSMICATLASS Makefile - clean, object-based compilation
# Copyright (c) 2013-2024 A Balaguera
##########################################################################

# Compiler
CXX = g++
CXXFLAGS = --std=c++17 -fopenmp -fPIC -g3 -funroll-loops

# Paths
AUX_DIR = ../external_libs/
SOURCE_DIR = ../source/
HEADER_DIR = ../headers/
LCOSMO_DIR = ../lib
dir_CFITSIO_LIB=$(AUX_DIR)cfitsio/
# Libraries
PYTHON_INC := $(shell python3-config --embed --includes)
PYTHON_LIB := $(shell python3-config --embed  --ldflags)
LIBS_GSL     = -lgsl -lgslcblas -lopenblas -lutil -lboost_iostreams -lboost_system -lboost_filesystem -lm $(PYTHON_LIB)
LIBS_FFTW    = -L/usr/local/lib -lfftw3 -lfftw3f -lfftw3f_omp -lfftw3_omp
LIBS_HEALPIX =  -L$(dir_CFITSIO_LIB) -L$(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/lib/ -lhealpix_cxx -lcxxsupport -lsharp -lc_utils -lfftpack -lcfitsio -lgomp
LCOSMO       = -Wl,-rpath,$(LCOSMO_DIR) -L$(LCOSMO_DIR) -lCosmiCalCpp

# Include paths
INCLUDE = -I$(AUX_DIR)Healpix_3.31/src/cxx/generic_gcc/include/ -I/usr/local/include/ -I$(HEADER_DIR) $(PYTHON_INC)

# Pattern rule to compile any .cpp -> .o
%.o: %.cpp $(HEADER_DIR)/*.h
	@echo "Compiling $<"
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

# Executables and object files
TARGET_CLib = cosmolib.exe
OBJS_CL     = ../codes/cosmolib.o 

TARGET_BMT  = bmt.exe
OBJS_BMT    = ../codes/bmt.o 

TARGET_PS   = power.exe
OBJS_PS     = ../codes/power_spectrum.o

TARGET_SB   = secbias.exe
OBJS_SB     = ../codes/secbias.o

TARGET_GT   = gtools.exe
OBJS_GT     = ../codes/galtools.o 

TARGET_HT   = htools.exe
OBJS_HT     = ../codes/analysis_halo_catalogs.o

TARGET_CW   = cosmicweb.exe
OBJS_CW     = ../codes/cosmic_web.o

TARGET_PA   = params.exe
OBJS_PA     = ../codes/paramsr.o

# Default target
all: $(TARGET_BMT)

##########################################################################
# Linking rules
$(TARGET_CLib): $(OBJS_CL)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_CL) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_BMT): $(OBJS_BMT)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_BMT) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_PS): $(OBJS_PS)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_PS) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_SB): $(OBJS_SB)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_SB) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_GT): $(OBJS_GT)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_GT) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_HT): $(OBJS_HT)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_HT) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_CW): $(OBJS_CW)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_CW) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

$(TARGET_PA): $(OBJS_PA)
	@echo "Linking $@"
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS_PA) $(LCOSMO) $(LIBS_HEALPIX) $(LIBS_FFTW) $(LIBS_GSL)

##########################################################################
# Aliases
cosmolib: $(TARGET_CLib)
bmt:     $(TARGET_BMT)
power:   $(TARGET_PS)
secbias: $(TARGET_SB)
gtools:  $(TARGET_GT)
htools:  $(TARGET_HT)
cosmicweb: $(TARGET_CW)
params:    $(TARGET_PA)

##########################################################################
# Clean
clean:
	@echo "Cleaning..."
	rm -f $(OBJS_CL) $(OBJS_BMT) $(OBJS_PS) $(OBJS_SB) $(OBJS_GT) \
	      $(OBJS_HT) $(OBJS_CW) $(OBJS_PA) *.exe core.* *~ 

##########################################################################
.PHONY: all clean cosmolib bmt power secbias gtools htools cosmicweb params

#########################################################################
help:
	@echo "$(GREEN)CosmicCodes compiling options"
	@echo ""
	@echo "$(BLUE)make params: $(GREEN)creates executable for the params code "
	@echo "	  $(GREEN)execute with $(BLUE)./params.exe -n name_of_parameter $(GREEN)to get information on a given parameter of the json parameter file"		
	@echo "	  $(GREEN)execute with $(BLUE)./params.exe -s parameter_file.json $(GREEN)to get information of the set of parameters present in the json parameter file"		
	@echo ""
	@echo "$(BLUE)make cosmolib: $(GREEN)creates executable for the Cosmolib code "
	@echo "	  $(GREEN)execute with $(BLUE)./cosmolib.exe --hmodel parameter_file.json $(GREEN)to get halo model outputs (with plots)"
	@echo ""
	@echo "$(BLUE)make bmt: $(GREEN)creates executable for the Bias Mapping Technique"
	@echo "	  $(GREEN)execute with $(BLUE)./bmt.exe -b parameter_file.json"
	@echo "	  $(GREEN)execute with $(BLUE)./bmt.exe -d to check preprocessor directives"
	@echo ""
	@echo "$(BLUE)make power: $(GREEN)creates executable for the power spectrum and related operations"
	@echo "	  $(GREEN)execute with $(BLUE)./power.exe -p parameter_file.json $(GREEN)to measure power (see documentation)"
	@echo "	  $(GREEN)execute with $(BLUE)./power.exe -w parameter_file.json $(GREEN)to measure window matrix"
	@echo "	  $(GREEN)execute with $(BLUE)./power.exe -g parameter_file.json $(GREEN)to obtain GRF"
	@echo "	  $(GREEN)execute with $(BLUE)./power.exe -s parameter_file.json $(GREEN)for low pass filter"		
	@echo ""
	@echo "$(BLUE)make secbias: $(GREEN)creates executable for the analysis of secondary bias"
	@echo "	  $(GREEN)execute with $(BLUE)./secbias.exe -s parameter_file.json" 
	@echo ""
	@echo "$(BLUE)make gtools: $(GREEN)creates executable for the analysis of galaxy samples "
	@echo "	  $(GREEN)execute with $(BLUE)./gtools.exe parameter_file.json"
	@echo "$(RED)Note: This is still under merging stages: merged with class::Catalog"		
	@echo ""
	@echo "$(BLUE)make htools: $(GREEN)creates executable for the analysis of halo catalogs "
	@echo "	  $(GREEN)execute with $(BLUE)./htools.exe parameter_file.json" 
	@echo "$(RED)Note: This is still under merging stages: merged with class::Catalog"		
	@echo ""
	@echo "$(BLUE)make cosmicweb: $(GREEN)creates executable for the cosmic web analysis ased on tracer catalogues"
	@echo "	  $(GREEN)execute with $(BLUE)./cosmicweb.exe parameter_file.json"		
	@echo ""
	@echo "$(BLUE)make params: $(GREEN)creates executable $(GREEN)for the params code "
	@echo "	  $(GREEN)execute with $(BLUE)./params.exe -n name_of_parameter $(GREEN)to get information on a given parameter of the json parameter file"		
	@echo "	  $(GREEN)execute with $(BLUE)./params.exe -s parameter_file.json $(GREEN)to get information of the set of parameters present in the json parameter file"		
	@echo ""

