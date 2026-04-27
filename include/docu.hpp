// DOCUMENTATION FOR DOXYGEN COSMICATLAS
/**
* @mainpage
*
* @section intro_sec Introduction
 CosmiCCodes++: A C++ library for calculations in cosmology with general and detailed applications.<br>

*@subsection Applications
The cosmicatlass code is aimed to perform general as well as very specific tasks in the context of large scale structure of the Universe. Among the detailed taks we can identify:<br>
<b>Generation of halo mock catalogs based on the bias mapping method \b<br>
<b>Measurement of 3D power spectrum of dark matter tracers \b <br>
<b>Measurement of angular power spectrum and tomographic analysis \b <br>
<b>Generation of random catalogs for clustering analysis \b <br>
<b>Generation of lognormal catalogs for tomographic analysis\b<br>
<b>Parameter constraints with MonteCarlo Markov chains \b<br>
<b>Measurement of abundance of dark matter tracers \b<br>
<b>Analysis of secondary halo bias \b<br>
<tr>
<td>


* @section pap Related publications
<a href="https://ui.adsabs.harvard.edu/abs/2023arXiv231112991B/abstract">Secondary halo bias </a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2023A%26A...673A.130B/abstract">Mock comparison project for DESI</a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2021arXiv210707917S/abstrac">Hydro-BAM II</a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2021ApJ...921...66S/abstract">Hydro-BAM I</a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200511598K/abstract">Connection to PT</a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.493..586P/abstract">Comparison with other methods</a><br>
<a href="http://adsabs.harvard.edu/abs/2019arXiv190606109B">One simulation to have them all</a><br>
<a href="http://adsabs.harvard.edu/abs/2019MNRAS.483L..58B">The bias Assigment Method</a><br>
<a href="https://ui.adsabs.harvard.edu/abs/2012hcxa.confE..71B/abstract">Power Spectrum analysis</a><br>


*@section dw Download the code
The code will be soon made available. It is temporarly a private directory at  <a href="https://github.com/balaguera/CosmicCodes">GitHub</a> <br>
@note: Code will be subject to GNU-type licence.


*@section ex Executing the Code 

* @subsection rr Running CosmicAtlas
This file is compiled with several options. The main function is in the file 
\verbatim Main/cosmicatlass.cpp \endverbatim
The different options supported at the current version are summarized by the otion -h:
*@code
Options:         
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
			execute with ./secbias.exe -s parameter_file 
		-make gtools: creates executable for the analysis of galaxy samples 
			execute with ./gtools.exe parameter_file (this is different from otopns 1-4). This is still under merging stages: merged with class::Catalog
		-make htools: creates executable for the analysis of halo catalogs 
			execute with ./htools.exe -c parameter_file to analyze catalog
			execute with ./htools.exe -b parameter_file to assign individual bias to input catalog
			execute with ./htools.exe -s parameter_file to measure bias as a function of properteis
			execute with ./htools.exe -h parameter_file to get galaxy catalogs using hod
		-make cosmolib: Cosmology functions and statistical tools (bias, mass function, power spectrum)
			execute with ./cosmolib.exe -s parameter_file. to get theoretical outputs for different cosmological functions at fixed redshift
			execute with ./cosmolib.exe -z parameter_file. to get theoretical outputs for HOD predictions (under constuction)
*@endcode


*@subsection spar Parameter file
The input paramter file contains a list of parameters (numbers, strings, lists) needed to execute the code.

@note Not all the parameters in the parameter file are requested for the different options.

* @subsection rc Running Cosmolib
The cosmolib is devoted to calculate predictions of different cosmological observables such as the following:
Mass function <br>
Halo bias<br>
Linear matter power spectrum <br>
Non linear matter power spectrum<br>
Linear and non-liner galaxy power spectrum based on HOD<br>
Galaxy and dark matter correlation function<br>
Angular power spectrum <br>

*@subsubsection ssex Execution
code has its one parameter file. 
The user can also compile the cosmic_lib code
*@code
make clean;
make cosmiclib
./cosmiclib.exe parameter_cosmiclib.ini
*@endcode
</td>



@note: Library dependencies <br>
GSL <br>
Gnuplot<br>
Healpix<br>



 * @section _doc Documentation
A detailed documentation of the code can be found <a href="https://abalant.wixsite.com/abalan">Here</a> <br>


 * @section _author Author
 <a href="https://abalant.wixsite.com/abalan"> Andrés Balaguera-Antolínez</a> <br>
This code has been developed during working periods at different institutions:<br>
Instituto de Astrofísica de Canarias <br>
Univesita degli Studi Roma 3 <br>
Argelander Institut for Astronomy <br>
Max Planck Institute for Extraterrestrial Physics <br>
With the collaboration of many collabortors (co-authors) <br>
*/

 <td>
 <div class="column"></div>
 <img src="../Headers/flow.jpg" alt="BAM Flowchart" width="600" height="700">
</td>

 </tr>
 </table>



