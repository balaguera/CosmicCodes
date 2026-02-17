A C++ library to perform calculations in Cosmology (theoretical and observational). 
Theoretical cosmology: prediction for 1 and 2 point statistics (abundance, power spectrum in 2D and 3D).
Observational: measurements of abundance (e.g. luminosity functions, power spectrum, correlation function), marked statistics, halo properties.

Documentation is available in Doxygen and pdf. 
See html/index.html


It can be requested via e-mail to abalant@gmail.com.

The maike file in /main/ will react under the command make help as
CosmiCalc compiling options

make params: creates executable for the params code
execute with ./params.exe -n name_of_parameter to get information on a given parameter of the json parameter file
execute with ./params.exe -s parameter_file.json to get information of the set of parameters present in the json parameter file

make cosmolib: creates executable for the Cosmolib code
execute with ./cosmolib.exe --hmodel parameter_file.json to get halo model outputs (with plots)

make bmt: creates executable for the Bias Mapping Technique
execute with ./bmt.exe -b parameter_file.json
execute with ./bmt.exe -d to check preprocessor directives

make power: creates executable for the power spectrum and related operations
execute with ./power.exe -power parameter_file.json to measure power (see documentation)
execute with ./power.exe -window parameter_file.json to measure window matrix
execute with ./power.exe -grf parameter_file.json to obtain GRF
execute with ./power.exe -lowpassfilter parameter_file.json for low pass filter

make secbias: creates executable for the analysis of secondary bias
execute with ./secbias.exe -s parameter_file.json

make gtools: creates executable for the analysis of galaxy samples
execute with ./gtools.exe parameter_file.json
Note: This is still under merging stages: merged with class::Catalog

make htools: creates executable for the analysis of halo catalogs
execute with ./htools.exe parameter_file.json
Note: This is still under merging stages: merged with class::Catalog

make cosmicweb: creates executable for the cosmic web analysis based on tracer catalogues
execute with ./cosmicweb.exe -ibias parameter_file.json to assign individual bias to tracers
execute with ./cosmicweb.exe -cwc parameter_file.json to perform cosmic-web analysis

make params: creates executable for the params code
execute with ./params.exe -n name_of_parameter to get information on a given parameter of the json parameter file
execute with ./params.exe -s parameter_file.json to get information of the set of parameters present in the json parameter file
