## Description
Library in C++ to perform calculations in Cosmology (theoretical and observational). 
Theoretical cosmology: prediction for 1 and 2 point statistics (abundance, power spectrum in 2D and 3D).
Observational: measurements of abundance (e.g. luminosity functions, power spectrum, correlation function), marked statistics, halo properties.

## Documentation
Documentation is available in dOxygen and pdf. It can be requested via e-mail to abalant@gmail.com.


## Build
Follow the follow instructions:

$mkdir build

$cd build

$cmake ..

$cmake --build . -- VERBOSE=1


## Usage
In examples/ there are codes with examples if the diferent usages of the library. These are compiled with the Makefile
and executed with the following options:

-make bmt: creates executable for the computation of bias related propertioes of dark matter tracers
	execute with ./bmt.exe -bmt parameter_file.json to run the Bias Mapping Technique.
    execute with ./bmt.exe -ibias parameter_file.json to assign individual bias to input tracers.
	execute with ./bmt.exe -pbias parameter_file.json to produce plots of bias vs halo properties.
	execute with ./bmt.exe -sbias parameter_file.json to perform an analysis of secondary bias.



-make power: creates executable for the power spectrum and related operations
	  execute with ./power.exe -power parameter_file to measure power (see documentation)
	  execute with ./power.exe -window parameter_file to measure window matrix
	  execute with ./power.exe -grf parameter_file to obtain GRF.
	  execute with ./power.exe -lowpassfilter parameter_file for low pass filter.

-make angpower: creates executable for the measurement of angular power spectrum"
      execute with ./hgaps.exe -cl parameter_file.json 
	  to measure power (see documentation)"


-make gtools: creates executable for the analysis of galaxy samples 
	  execute with ./gtools.exe parameter_file (this is different from otopns 1-4). 

-make htools: creates executable for the analysis of halo catalogs 
	  execute with ./htools.exe -c parameter_file to analyze catalog
	  execute with ./htools.exe -b parameter_file to assign individual bias to input catalog
	  execute with ./htools.exe -s parameter_file to measure bias as a function of properteis
	  execute with ./htools.exe -h parameter_file to get galaxy catalogs using hod
	  execute with ./htools.exe -m parameter_file to build a mock catalog based on simulation (snapshot) using some dNdz

-make cosmolib: Cosmology functions and statistical tools (bias, mass function, power spectrum)
	  execute with ./cosmolib.exe -c parameter_file. to get theoretical outputs for different cosmological functions at fixed redshift
	  execute with ./cosmolib.exe -h parameter_file. to get theoretical outputs for HOD predictions (under constuction)
	  execute with ./cosmolib.exe -z parameter_file. to get theoretical outputs for redshift dependent quantities
