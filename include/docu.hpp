
/**
 * @mainpage CosmiCCodes
 * @brief A \c C++ library to perform calculations in theoretical and observational cosmology.
*
 * @section contents Table of Contents
 *
 * - @ref intro_sec "Introduction"
 * - @ref dw "Download the Code and compilation"
 * - @ref exa "Applications and Examples"
 *   - @ref spar "Input Parameter File"
 *   - @ref gen "A detailed example"
 *   - @ref ssbias "Bias Analysis"
 *   - @ref sspower "Power Spectrum"
 *   - @ref ssapower "Angular Power Spectrum"
 *   - @ref sshtool "Halo Tools"
 *   - @ref ssgtool "Galaxy Tools"
 *   - @ref cwlib "Cosmic-Web Analysis"
 *   - @ref clib "Cosmological Functions"
 *
 * - @ref sfaqs "FAQS"
 * - @ref dep "Library Dependencies"
 * - @ref pap "Related Publications"
 * - @ref _author "Author"
 * - @ref _author_c "Contact"
 *
*
*
 * @section intro_sec Introduction
 *
 * The **CosmicCodes** is a \c C++ library which provides tools for the statistical analysis of the
 * large-scale structure of the Universe, using both N-body simulations or galaxy redshift surveys. Its main applications include:
 *
 * - Analysis of primary and secondary halo bias.
 * - Generation of halo mock catalogs using the Bias Assignment Method.
 * - Measurement of the three-dimensional power spectrum of dark-matter tracers.
 * - Measurement of angular power spectra and tomographic analysis.
 * - Generation of random catalogs for clustering analyses.
 * - Generation of lognormal catalogs for tomographic analyses.
 * - Cosmological parameter estimation using Markov Chain Monte Carlo methods.
 * - Measurement of abundance statistics, including mass and luminosity functions.
 * - Analysis of halo and galaxy catalogs.
 * - Cosmic-web analysis.
 *
 *
 * @section dw Download the code and compilation
 *
 * The source codecan bu downloaded from 
 * <a href="https://github.com/balaguera/CosmicCodes">GitHub</a>. 
 *
 * @note The code is distributed under a GNU-type license.
 *
 * The library is compiled using the CMake file available at <a href = "https://github.com/balaguera/CosmicCodes"> GitHub</a>. From the root directory of the project:
 *
 * \code{.sh}
 * mkdir build
 * cd build
 * cmake ..
 * cmake --build . -- VERBOSE=1
 * \endcode
 *
 *
 * @section exa Applications and Examples
 *
 * Example programs are available in the
 * \c examples/ directory. The corresponding executables can be compiled
 * using the Makefile located in \c main/.
 *
 * The available tools are described below.
 *
 *
 * @subsection spar Input Parameter File
 *
 * The input parameter file is written in \c JSON (JavaScript object notation) format, with a number of sections dedicated to 
 * either input catalogs of methods (calculations such as power spectrum, luminosity function, cosmological parameter constraints etc). 
 *
 * 
 * @note Go to \ref parameters_docu for detailed information on the parameters used in this library.
 * 
 * 
 * @subsection gen A detailed example
 * 
 * In the section \ref expl we show a more detail the usage of some classes of this library, as well as quick demonstrations and explanations.
 *   
 * 
 * 
 * @subsection ssbias Bias Analysis
 * This example aims at performing the analysis of cosmological bias using as input dark matter and dark matter tracer catalogs.
 * The <a href = "../../../examples/bmt.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make bmt
 * \endcode
 *
 * Available options:
 *
 * - **Assignment of individual bias to input tracers** (see e.g. @ref PowerSpectrumF::object_by_object_bias())
  *   \code{.sh}
 *   ./bmt.exe -ibias parameter_file.json
 *   \endcode
 *
 * - **Generate plots of bias as a function of halo properties**
  *   \code{.sh}
 *   ./bmt.exe -pbias parameter_file.json
 *   \endcode
 *
 * - **Perform a secondary-bias analysis**
  *   \code{.sh}
 *   ./bmt.exe -sbias parameter_file.json
 *   \endcode
 *
 * - **Bias Mapping Technique**
 *   \code{.sh}
 *   ./bmt.exe -bmt parameter_file.json
 *   \endcode
 *
 * @subsection sspower Power Spectrum
 *
 * The library offers the possibility to measure the 3D power spectrum from N-body simulations or redshift galaxy surveys. The <a href = "../../../examples/power_spectrum.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make power
 * \endcode
 *
 * Available options:
 *
 * - **Measure the power spectrum of an input catalog**
  *   \code{.sh}
 *   $ ./power.exe -power parameter_file
 *   \endcode
 *
 * - **Measure the window matrix using an input random catalog** (see e.g. @ref PowerSpectrumF::get_window_matrix_multipole())
  *   \code{.sh}
 *   $ ./power.exe -window parameter_file
 *   \endcode
 *
 * - **Generate Gaussian random fields from an input power spectrum** (see e.g. @ref PowerSpectrumF::get_GaussianRandomField())
  *   \code{.sh}
 *   $ ./power.exe -grf parameter_file
 *   \endcode
 *
 * - **Apply a low-pass filter to an input density field**
  *   \code{.sh}
 *   $ ./power.exe -lowpassfilter parameter_file
 *   \endcode
 * 
 * See @ref PowerSpectrumF::measure_power_spectrum_data() 
 *
 *
 * @subsection ssapower Angular Power Spectrum
 *
 * The library offers the possibility to perform tomographic analysis for galaxy redshift surveys using angular spectrum and redshift space distortions. The <a href = "../../../examples/power_spectrum.cpp">example code</a> can be compiled with
 * The <a href = "../../../examples/hgaps.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make angpower
 * \endcode
 * Available options:
 *
 * - **Measure the angular power spectrum from an input catalog and mask**
  *   \code{.sh}
 *   ./hgaps.exe -power parameter_file.json
 *   \endcode
 *
 * - **Measure cross-power spectra in different tomographic bins**
  *   \code{.sh}
 *   ./hgaps.exe -cross parameter_file.json
 *   \endcode
 *
 * @note The example to build log normal catalogs is under construction.
 * 
 *
 * @subsection sshtool Halo Tools
 * With \ref HaloTools we enclose different calculations and statistics that can be performed over a halo population read from an N-body simulation.
 * The <a href = "../../../examples/halo_analysis.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make htools
 * \endcode
 *
 * Available options:
 *
 * - **Analyze an input catalog**
 *
 *   \code{.sh}
 *   ./htools.exe -c parameter_file
 *   \endcode
 *
 * - **Assign individual bias to an input catalog**
 *
 *   \code{.sh}
 *   ./htools.exe -b parameter_file
 *   \endcode
 *
 * - **Measure halo bias as a function of selected properties**
 *
 *   \code{.sh}
 *   ./htools.exe -s parameter_file
 *   \endcode
 *
 * - **Build galaxy catalogs from an input halo catalog using an HOD model**
 *
 *   \code{.sh}
 *   ./htools.exe -h parameter_file
 *   \endcode
 *
 * - **Build a mock catalog from a simulation snapshot using an input dN/dz**
 *
 *   \code{.sh}
 *   ./htools.exe -m parameter_file
 *   \endcode
 *
 *
 * @subsection ssgtool Galaxy Tools
 * With \ref GalaxyTools we enclose different calculations and statistics that can be performed over a galaxy redshift survey.
 * The <a href = "../../../examples/galtools.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make gtools
 * \endcode
 *
 * The program is executed as:
 *
 * \code{.sh}
 * ./gtools.exe parameter_file
 * \endcode
 *
 * @note This module is currently under development.
 *
 *
 * @subsection cwlib Cosmic-Web Analysis
 * With this excample we aim at showing the approach to follow, using a halo population read from an N-body simulation, to perform cosmic-web analysis based on the todal field -or velocity shear- from the tracers or de dark matter particles in a N-body simulation.
 * The <a href = "../../../examples/cosmic_web.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make cosmicweb
 * \endcode
 *
 * To perform a cosmic-web analysis:
 *
 * \code{.sh}
 * ./cosmicweb.exe -cwc parameter_file.json
 * \endcode
 *
 *
 * @subsection clib Cosmological Functions
 *
 * The <a href = "../../../examples/cosmolib.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make cosmolib
 * \endcode
 *
 * Available options:
 *
 * - **Obtain theoretical predictions for cosmological functions at fixed redshift**
 *
 *   \code{.sh}
 *   ./cosmolib.exe -c parameter_file
 *   \endcode
 *
 * - **Obtain theoretical HOD predictions**
 *
 *   \code{.sh}
 *   ./cosmolib.exe -h parameter_file
 *   \endcode
 *
 * - **Obtain theoretical predictions for redshift-dependent quantities**
 *
 *   \code{.sh}
 *   ./cosmolib.exe -z parameter_file
 *   \endcode
 *
 * @note HOD predictions are currently under development.
 *

 * @section sfaqs FAQS
 *
 * - **Do I need to know C++ to use the library?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">Well...</span>
 *     <span class="tooltiptext">
 *      Sure! So if you are interested in Cosmolgy and c++, better be ready with some basic stuff about this programing language.       
 *     </span>
 *   </span>

 * - **Do I need to know cosmology to use the library?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">Well...</span>
 *     <span class="tooltiptext">
 *      Sure! But not to be an expert. One of the main goals of the library is that you get familiar with some basic concepts of observational cosmology whilye playing with it.        
 *     </span>
 *   </span>

 * - **Can I help?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">Well...</span>
 *     <span class="tooltiptext">
 *      Sure! If you have some module which can be implemented in a certain cosmological context, we can add it to the libray.        
 *     </span>
 *   </span>


 * - **How can I read an input catalogue?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">See how...</span>
 *     <span class="tooltiptext">
 *       The method @ref Catalogue::read_catalog_new()
 *       allocates the properties of the input ASCII catalogue in a number
 *       of vectors. As an example, the information of the mass of the i-th tracer is retrieved through methods of the form \c cat.mass_at(i). The section @ref TracerCatalogue must be enabled in
 *       the input parameter file. See more information in @ref Catalogue .
 *     </span>
 *   </span>
 *
 * - **How can I measure 3D power spectrum?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">See how...</span>
 *     <span class="tooltiptext">
 *       The class @ref PowerSpectrumF offers a number of possibilities. In the simplest case, youn can load an N-body 
 *       simulation and use \ref PowerSpectrumF::measure_power_spectrum_box; as a guide, follow
 *       @code
 *        string par_file=argv[2];
 *        Params params(par_file);
 *        PowerSpectrumF cPSF(params);
 *        cPSF.measure_power_spectrum_box();
 *       @endcode 
         See more information in @ref PowerSpectrumF and @ref FftwFunctions.
 *     </span>
 *   </span>


 * - **How can I measure the bias of dark matter tracers?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">See how...</span>
 *     <span class="tooltiptext">
 *       The class @ref PowerSpectrumF offers two aproaches. The first uses the standard estimationr through ratios of power spectra, \ref PowerSpectrumF::get_lss_bias(). 
 *       The second option is bsed on the assignment of individual tracer bias, using  \ref PowerSpectrumF::object_by_object_bias().
 *       See more information in @ref PowerSpectrumF.
 *     </span>
 *   </span>


 * - **How can I measure halo abundance?**
 *   <span class="tooltip">
 *     <span class="tooltip-trigger">See how...</span>
 *     <span class="tooltiptext">
 *       Once you have read the parameter file and the catalogue, create an instance of the class \ref HaloTools and use the method @ref HaloTools::get_property_function()
         See more information in @ref HaloTools::analyze_cat().
 *     </span>
 *   </span>


 * - **How can I generate a tracer field?**

 * - **How can I measure luminosity function?**




 * @section dep Library Dependencies
 *
 * The library currently depends on:
 *
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c GSL</span>
 *     <span class="tooltiptext">
        Used everywhere for mathematical calculations (integrations, interpolation, matrix algebra)
 *     </span>
 *   </span>
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c FFTW</span>
 *     <span class="tooltiptext">
 *     Unavoidable for Fourier analysis.
 *     </span>
 *   </span>
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c HealPix</span>
 *     <span class="tooltiptext">
 *     Used for angular analysis in harmonic space (angular power spectrum, Fourier-Bessel decomposition) 
 *     </span>
 *   </span>
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c Gnuplot</span>
 *     <span class="tooltiptext">
*      Used to show some results in codes going through large loops.
 *     </span>
 *   </span>
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c Python</span>
 *     <span class="tooltiptext">
 *      Used to show basic results at the end of some computation.       
 *     </span>
 *   </span>
 *   </span>
 * -   <span class="tooltip">
 *     <span class="tooltip-trigger">\c JSON</span>
 *     <span class="tooltiptext">
        Format to print out plotting information and input parameter file.
 *     </span>
 *   </span>


 *
 *
 *
 * @section pap Related Publications
 * 
 *  This library has been developed for -or used- in a number of research papers in the field of observational and theoretical cosmology. Among others, 
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2023arXiv231112991B/abstract">
 *   Secondary halo bias through cosmic time. I. Scaling relations and the connection with the cosmic web </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2026A%26A...711A.267A/abstract">
 *   The galaxy bias profile of cosmic voids: A comparison of void finders  </a>
 *
 * - <a href = "https://ui.adsabs.harvard.edu/abs/2026A%26A...711A.186R/abstract">Predicting galaxy bias using machine learning</a> 
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2025A%26A...700A.207G/abstract">Individual halo bias in models of f(R) gravity 
 *    </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2025A%26A...700A.207G/abstract">Power spectrum of the CODEX clusters  </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2023A%26A...673A.130B/abstract">
 *   DESI mock challenge. Halo and galaxy catalogues with the bias assignment method </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2021arXiv210707917S/abstract">
 *   Mapping the Three-dimensional Lyα Forest Large-scale Structure in Real and Redshift Space </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2021ApJ...921...66S/abstract">
 *   The Bias from Hydrodynamic Simulations: Mapping Baryon Physics onto Dark Matter Fields</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200511598K/abstract">
 *   The cosmic web connection to the dark matter halo distribution through gravity</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.493..586P/abstract">
 *   The bias of dark matter tracers: assessing the accuracy of mapping techniques </a>
 *
 * - <a href="http://adsabs.harvard.edu/abs/2019arXiv190606109B">
 *    One simulation to have them all: performance of the Bias Assignment Method against N-body simulations</a>
 *
 * - <a href="http://adsabs.harvard.edu/abs/2019MNRAS.483L..58B"> BAM: bias assignment method to generate mock catalogues </a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2025A%26A...700A.207G/abstract"> Extracting cosmological information from the angular power spectrum of the 2MASS Photometric Redshift catalogue </a>
 
 * - <a href="https://ui.adsabs.harvard.edu/abs/2012hcxa.confE..71B/abstract">
 *   The REFLEX II Galaxy Cluster sample: mock catalogues and clustering analysis</a>
 *
 * 
 *
 * @section _author Author
 *
 * <a href="https://abalant.wixsite.com/abalan">
 * Andrés Balaguera-Antolínez (PhD) </a>
 *
 * The development of the code has benefited from the contributions of
 * numerous collaborators and co-authors, during research positions at:
 *
 * - Instituto de Astrofísica de Canarias (2017–2024)
 * - Università degli Studi Roma Tre (2013–2016)
 * - Argelander-Institut für Astronomie (2010–2013)
 * - Max Planck Institute for Extraterrestrial Physics (2006–2010)
 *
 * 
 * @section _author_c Contact
 * 
 * Questions related to the library, please write to abalant@gmail.com
 *
 * @image html flow.jpg "BAM Flowchart"
 */






