
/**
 * @mainpage CosmicCodes Documentation
 *
 * @section contents Table of Contents
 *
 * - @ref intro_sec "Introduction"
 * - @ref dw "Download the Code"
 * - @ref ex "Compilation"
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
 * - @ref dep "Library Dependencies"
 * - @ref _doc "Additional Documentation"
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
 * @section dw Download the Code
 *
 * The source codecan bu downloaded from 
 * <a href="https://github.com/balaguera/CosmicCodes">GitHub</a>. 
 *
 * @note The code is distributed under a GNU-type license.
 *
 *
 * @section ex Compilation
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
 *
 *   \code{.sh}
 *   $ ./power.exe -power parameter_file
 *   \endcode
 *
 * - **Measure the window matrix using an input random catalog** (see e.g. @ref PowerSpectrumF::get_window_matrix_multipole())
 *
 *   \code{.sh}
 *   $ ./power.exe -window parameter_file
 *   \endcode
 *
 * - **Generate Gaussian random fields from an input power spectrum** (see e.g. @ref PowerSpectrumF::get_GaussianRandomField())
 *
 *   \code{.sh}
 *   $ ./power.exe -grf parameter_file
 *   \endcode
 *
 * - **Apply a low-pass filter to an input density field**
 *
 *   \code{.sh}
 *   $ ./power.exe -lowpassfilter parameter_file
 *   \endcode
 *
 *
 * @subsection ssapower Angular Power Spectrum
 *
 * The <a href = "../../../examples/hgaps.cpp">example code</a> can be compiled with
 *
 * \code{.sh}
 * make angpower
 * \endcode
 *
 * Available options:
 *
 * - **Measure the angular power spectrum from an input catalog and mask**
 *
 *   \code{.sh}
 *   ./hgaps.exe -cl parameter_file.json
 *   \endcode
 *
 * - **Measure cross-power spectra in different tomographic bins**
 *
 *   \code{.sh}
 *   ./hgaps.exe -r parameter_file.json
 *   \endcode
 *
 *
 * @subsection sshtool Halo Tools
 *
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
 *
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
 *
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
 *
 * @section dep Library Dependencies
 *
 * The library currently depends on:
 *
 * - \c GSL
 * - \c Gnuplot
 * - \c HEALPix
 * - \c Python
 *
 *
 * @section _doc Documentation
 *
 * Additional documentation is available
 * <a href="https://abalant.wixsite.com/abalan">here</a>.
 *
 *
 * @section pap Related Publications
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2023arXiv231112991B/abstract">
 *   Secondary halo bias</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2023A%26A...673A.130B/abstract">
 *   Mock comparison project for DESI</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2021arXiv210707917S/abstract">
 *   Hydro-BAM II</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2021ApJ...921...66S/abstract">
 *   Hydro-BAM I</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2020arXiv200511598K/abstract">
 *   Connection to perturbation theory</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2020MNRAS.493..586P/abstract">
 *   Comparison with other methods</a>
 *
 * - <a href="http://adsabs.harvard.edu/abs/2019arXiv190606109B">
 *   One simulation to have them all</a>
 *
 * - <a href="http://adsabs.harvard.edu/abs/2019MNRAS.483L..58B">
 *   The Bias Assignment Method</a>
 *
 * - <a href="https://ui.adsabs.harvard.edu/abs/2012hcxa.confE..71B/abstract">
 *   Power-spectrum analysis</a>
 *
 *
 * @section _author Author
 *
 * <a href="https://abalant.wixsite.com/abalan">
 * Andrés Balaguera-Antolínez</a>
 *
 * This code has been developed during research positions at:
 *
 * - Instituto de Astrofísica de Canarias (2017–2024)
 * - Università degli Studi Roma Tre (2013–2016)
 * - Argelander-Institut für Astronomie (2010–2013)
 * - Max Planck Institute for Extraterrestrial Physics (2006–2010)
 *
 * The development of the code has benefited from the contributions of
 * numerous collaborators and co-authors.
 * 
 * @section _author_c Contact
 * 
 * Questions related to the library, please write to abalant@gmail.com
 *
 * @image html flow.jpg "BAM Flowchart"
 */






