Library in C++ to perform calculations in Cosmology (theoretical and observational). 
Theoretical cosmology: prediction for 1 and 2 point statistics (abundance, power spectrum in 2D and 3D).
Observational: measurements of abundance (e.g. luminosity functions, power spectrum, correlation function), marked statistics, halo properties.

Documentation is available in dOxygen and pdf. It can be requested via e-mail to abalant@gmail.com.

The maike file in /main/ will react under the command make help as
CosmiCalc compiling options

-make bmt: creates executable for the Bias Mapping Technique
	  execute with ./bmt.exe -b parameter_file to run bmt tasks.
	  execute with ./bmt.exe -d to check preprocessor directives.
	  execute with ./bmt.exe -i parameter_file to check input parameters directives.

-make power: creates executable for the power spectrum and related operations
	  execute with ./power.exe -p parameter_file to measure power (see documentation)
	  execute with ./power.exe -w parameter_file to measure window matrix
	  execute with ./power.exe -g parameter_file to obtain GRF.
	  execute with ./power.exe -s parameter_file for low pass filter.
	  execute with ./power.exe -m parameter_file for marked correlation function

-make secbias: creates executable for the analysis of secondary bias
	  execute with ./secbias.exe -s parameter_file for low pass filter

-make gtools: creates executable for the analysis of galaxy samples 
	  execute with ./gtools.exe parameter_file (this is different from otopns 1-4). This is still under merging stages: merged with class::Catalog

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


Doxygen documentation available at html/index.html

