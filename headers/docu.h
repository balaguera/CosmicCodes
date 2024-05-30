// DOCUMENTATION FOR DOXYGEN COSMICATLAS
/**
* @mainpage
*
* @section intro_sec Introduction
 CosmiCalc++: A C++ library for calculations in cosmology with general and detailed applications.<br>

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



*@section dw Download the code
The code will be soon made available 
@note: Code will be subject to GNU-type licence.


*@section ex Executing the Code 

* @subsection rr Running CosmicAtlas
This file is compiled with several options
*@code
make clean;
make bam
./cosmicatlass.exe -option parameter_bam.ini
*@endcode

The main function is in the file 
\verbatim Main/cosmicatlass.cpp \endverbatim
The different options supported at the current version are summarized by the otion -h:
*@code
Options:         
             -a for information on the author, no argument
	         -b parameter_file.ini: runs BAM
	         -c parameter_file.ini: analyzes tracer catalog
	         -d parameter_file.ini: shows preprocessor directives
	         -f parameter_file.ini: reads binarys form IC and ocmpute density field IC
	         -g parameter_file.ini: generates galaxy catalog from halo catalog using HOD
	         -h parameter_file.ini: help
	         -i parameter_file.ini: shows input parameters
	         -m parameter_file.ini: measures marked power spectrum
	         -p parameter_file.ini: measures power spectrum
	         -q parameter_file.ini: compute rings of bias
	         -s parameter_file.ini: applies a low-pass filter to input density field
	         -t parameter_file.ini: runs LPT (to be deprecated)
	         -u parameter_file.ini: analyzes input halo catalog and measurements of secondary bias
	         -w parameter_file.ini: computes window functions
	         -v parameter_file.ini: warnings 
	         -z parameter_file.ini: assignment of individual tracer bias
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
A detailed documentation of the code can be found <a href ="https://www.overleaf.com/read/kpyrhbyxdbmj#89834c">at this Overleaf link </a>. Check <a href = "../Headers/TalkBam.pdf">here</a> for a talk on BAM (Innsbrück 2020).





 * @section _author Author
 <a href="https://abalant.wixsite.com/abalan"> Andrés Balaguera-Antolínez</a> <br>
This code has been developed during working periods at <br>
Instituto de Astrofísica de Canarias <br>
Univesita degli Studi Roma 3 <br>
Argelander Institut for Astronomy <br>
Max Planck Institute for Extraterrestrial Physics <br>
With the collaboration of <br>
Ariel Sanchez<br>
Federico Marulli<br>
Daniele Tavagnacco<br>
Luca Tornatore<br>
Jeniffer Pollack<br>
Francisco Kitaura<br>
Carlos Hernandez-Monteagudo<br>
Marcos Pellejero<br>
Francesco Sinigaglia<br>




*/



 <td>
 <div class="column"></div>
 <img src="../Headers/flow.jpg" alt="BAM Flowchart" width="600" height="700">
</td>

 </tr>
 </table>



