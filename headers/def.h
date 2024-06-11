// ****************************************************************************************
// ****************************************************************************************
/**
 * @file def.h
 * @brief Pre-processor directives for cosmicatlas
 * @author Andres Balaguera-Antolínez
 * @version 1.0 2018
 * @date   2018-2022
 * */
// ****************************************************************************************
// ****************************************************************************************
// This order disables the use of the assert() function
// #define NDEBUG
/**
*@brief Use OMP parallelization
*/
#define _USE_OMP_
/**
*@brief Use vectorization with OMP parallelization
*/
//#define _USE_SIMD_OMP_
//#define _USE_SIMD_OMP_mas_
#define _NTHREADS_ static_cast<int>(10)
// ****************************************************************************************
// ****************************************************************************************
#define NEGATIVE_INT static_cast<int>(-1)
#define POSITIVE_INT static_cast<int>(1)
// ****************************************************************************************
// ****************************************************************************************
/**
 * @brief Mock Mode
 * @details Define MOCK_MODE if calibration or mock production is to be run with cosmicatlas.
 *  If not defined, then the option BIAS follows, for which the bias statistics is computed
*/
#define MOCK_MODE // If MOCK_MODE is undefined, then define BIAS
// ****************************************************************************************
#ifndef MOCK_MODE
#define BIAS_MODE
#endif
// ****************************************************************************************
/**
 * @brief Mock Mode
 * @details Define MOCK_MODE if calibration or mock production is to be run with cosmicatlas.
 *  If not defined, then the option BIAS follows, for which the bias statistics is computed
*/
#ifdef MOCK_MODE
//#define _DO_BiasMT_CALIBRATION_
#ifndef _DO_BiasMT_CALIBRATION_
#define _GET_BiasMT_REALIZATIONS_
#endif
#endif
// ****************************************************************************************
#ifdef _DO_BiasMT_CALIBRATION_
#define _WRITE_AVERAGE_KERNEL_AND_BIAS_
#endif
// ****************************************************************************************
// ****************************************************************************************
/**
 * @brief Mock Mode
*/
#define _FULL_VERBOSE_

// ****************************************************************************************
#ifdef _FULL_VERBOSE_
//#define _VERBOSE_FREEMEM_
#define _VERBOSE_POWER_
//#define _VERBOSE_OVERDENS_
#define _VERBOSE_FILE_
#define _VERBOSE_CAT_
#endif
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
/**
 * @brief To show plots live
*/
#define _USE_GNUPLOT_
//#define FG_COLOR "white"
#define FG_COLOR "black"
// ****************************************************************************************
#ifdef _USE_GNUPLOT_
/**
 * @brief Define to show the 2d plot of halo bias P(Nh|delta)
*/
//#define _SHOW_BIAS_
/**
 * @brief Define to show pdf whwn gnerating DM density fiwlds with LPT
*/
//#define _USE_GNUPLOT_PDF_
/**
 * @brief Define to show power spectrum of IC
*/
//#define _USE_GNUPLOT_SREL_
#define _USE_GNUPLOT_POWER_
#define _USE_GNUPLOT_POWER_PLOT_
//#define _USE_GNUPLOT_ABUNDANCE_PLOT_
//#define _USE_GNUPLOT_ABUNDANCE_V_PLOT_
//#define _USE_GNUPLOT_ABUNDANCE_RS_PLOT_
//#define _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_
#endif // endif use_gnuplot
// ****************************************************************************************
#undef _USE_PYTHON_
// ****************************************************************************************
//#define _MG_
//#define _SLICS_
#define _UNITSIM_
//#define _TNG_
//#define _TNG_GAL_
//#define _FLAGSHIP_
// ****************************************************************************************
#ifdef _UNITSIM_
/**
 * @brief This factor corrects th growth computed in BiasMT agains the growth from FastPM. The difference in the two growthds is of about 2 percent.
*/
#define _correct_shape_theoretica_power_
#define factor_growth static_cast<real_prec>(1.018018854)
#else
#define factor_growth 1.0
#endif
// ****************************************************************************************
// When enabled, BiasMT uses the abacus cosmoloy and applies a half cell shift to the binning for
// the interpolation of tracers using ONLY NGP 
// and at TET for the DM evolved with LPT
//#define _ABACUS_
// ****************************************************************************************
#ifdef _ABACUS_
/**
 * @brief This section converts the fiducial input for abacus halo masses (given as number of particles) to Ms/h 
*/
#define PART_TO_MASSES_ABACUS static_cast<double>(2108093042.95605)
#endif
// ****************************************************************************************
// With this definition, two built-in mass-ranges are adopted in PowerSpectrumF.cpp
//#define _JPAS_
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// PREPROCESSOR DIRECTIVES FOR BiasMT-LPT (COSMICATLAS)************************************
// ****************************************************************************************
// Preproc for fast efficient runniing for andres. These only act in this file.
#define mode_b  // to run the -b option, it undefs _power_and mass bins, and _read_binary. Below we set that mode b deals with vmax as observable and mode p with mass, but this is not fixed, for we might want to use mode p with vmax as observable.
// ****************************************************************************************
#ifndef mode_b
#define mode_p  // does the contrary
#endif
// ****************************************************************************************
#define _VERBOSE_POWER_
// ****************************************************************************************
#ifdef mode_p
#define _VERBOSE_POWER_
#define _VERBOSE_CATALOG_
#endif
// ****************************************************************************************
#ifdef _SLICS_
//#define _USE_SLICS_PK_BINNING_
#endif
// ****************************************************************************************
#define _DOWNSAMPLING_
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
// ****************************************************************************************
// The following preproc directives are meant to specify the bias properties
// ****************************************************************************************
/**
 * @brief  CWEB uses the cosmic web classification CWEB only
*/
//#define _USE_CWEB_
// ****************************************************************************************
/**
 * @brief  MKNOTS uses ONLY the mass of collapsing regions (knots)
*/
#define _USE_MKNOTS_

// ****************************************************************************************
/**
 * @brief  ikweb uses the Iweb with the mass of collapsing regions (knots)
*/
//#define _USE_IKWEB_
#ifdef _USE_IKWEB_
#define MODEL_THETA "IKWEB"
#endif
// ****************************************************************************************
/**
 * @brief  TWEB uses the cosmic web classification CWEB plus the mass of collapsing regions (knots)
*/
#define _USE_TWEB_

#ifdef _USE_TWEB_
#define MODEL_THETA "TWEB"
//#define _USE_TIDAL_ANISOTROPY_
#endif


// ****************************************************************************************
#ifdef _USE_TIDAL_ANISOTROPY_
/**
 * @brief  Interpolates the tidal anisotropy at each halo position in order to uses this to assign halo primary property
*/
//#define _USE_ANISOTROPY_AT_TRACER_POSITION_
#endif
// ****************************************************************************************
/**
 * @brief  tiweb uses the CWC with the mass of collapsing regions (low res knots)
 * @details If below _USE_CURVATURE_ is dTNGefined, the eigenvalues are taken from the curvture
*/
//#define _USE_TIWEB_
#ifdef _USE_TIWEB_
#define MODEL_THETA "TIWEB"
#endif
// ****************************************************************************************
/**
 * @brief  IWEB uses the invariants of the tidal field I2 and I3. This can be used together with TEB,CWEB and PWEB
*/
//#define _USE_IWEB_
// ****************************************************************************************
#if defined (_USE_IWEB_) || defined (_USE_IKWEB_) 
//#define _USE_INVARIANT_TIDAL_FIELD_I_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
#define _USE_INVARIANT_TIDAL_FIELD_II_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
#define _USE_INVARIANT_TIDAL_FIELD_III_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
#endif
// ****************************************************************************************
#if defined (_USE_TIWEB_)
//#define _USE_INVARIANT_TIDAL_FIELD_I_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
//#define _USE_INVARIANT_TIDAL_FIELD_II_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
//#define _USE_INVARIANT_TIDAL_FIELD_III_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
#define _USE_TIDAL_ANISOTROPY_
#endif
// ****************************************************************************************
/**
 * @brief  PWEB Uses the invariats of the field ð_i ð_j delta
*/
//#define _USE_PWEB_
#ifdef _USE_PWEB_
//#define _USE_INVARIANT_PWEB_I_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
//#define _USE_INVARIANT_PWEB_II_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
//#define _USE_INVARIANT_PWEB_III_  //  This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II or s²delta
#endif
// ****************************************************************************************
// ****************************************************************************************
/**
 * @brief  CWEB_V Uses the invariats of the velocity shear and the mass of v-collapsing regions
*/
//#define _USE_CWEB_V
// ****************************************************************************************
/**
 * @brief  TWEB_V Uses the Cosmic web classification based on the eigenvañues of the velocity shear and the mass of collapsing regions defined by the shear
*/
//#define _USE_TWEB_V_
// ****************************************************************************************
/**
 * @brief  IWEB_V Uses the invariats of the velocity shear
*/
//#define _USE_IWEB_V
#ifdef _USE_IWEB_V
#define _USE_INVARIANT_SHEAR_VFIELD_I_  // def or undef here
#define _USE_INVARIANT_SHEAR_VFIELD_II_  // def or undef here
#define _USE_INVARIANT_SHEAR_VFIELD_III_  // def or undef here
//#define _USE_INVARIANT_SHEAR_VFIELD_IV_  // def or undef here
#endif
// ****************************************************************************************
// ****************************************************************************************
/**
 * @brief  AWEB Uses the ANISOTROPY PARAMETERS (ellipticity, anisotropy, prolatnes) of the tidal field
*/
//#define _USE_AWEB_
#ifdef _USE_AWEB_
#define _USE_TIDAL_ANISOTROPY_
//#define _USE_ELLIPTICITY_
//#define _USE_PROLATNESS_
#endif
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
// ======================================================================================================================================
#ifdef _UNITSIM_
//#define _USE_HALO_PROL_ELL_
//#define _USE_HALO_ELLIPTICITY_
#define _USE_CONCENTRATION_
#endif
// ****************************************************************************************
// ****************************************************************************************
// ************************COSMOLOGICAL PARAMETERS*****************************************
/**
 * @brief If _USE_COSMO_PARS_ is defined, the code uses the cosmological parametrs from  namespaces.
 * @details If undef, BiasMT uses cosmo params from input parameter file.
 * @details See the file cosmological_parameters.h for the namespaces with sets of cosmological parameters
*/
#define _USE_COSMO_PARS_
// ****************************************************************************************
/**
 * @brief Namespaces for cosmological parameters
*/
#if defined _UNITSIM_ || defined _TNG_ || defined _TNG_GAL_
#define _USE_UNITSIM_COSMOLOGY_
#endif
#ifdef _SLICS_
#define _USE_SLICS_COSMOLOGY_
#endif
#ifdef _MG_
#define _USE_MG_COSMOLOGY_
#endif
#ifdef _ABACUS_
#define _USE_ABACUS_COSMOLOGY_
#endif
#if !defined _USE_ABACUS_COSMOLOGY_ && !defined _USE_SLICS_COSMOLOGY_ && !defined _USE_MG_COSMOLOGY_ && !defined _FLAGSHIP_
#define _USE_PLANCK_COSMOLOGY_
#endif

#ifdef _USE_SLICS_COSMOLOGY_
#define COSMOPARS Cosmo_parameters_SLICS
#endif
#if defined (_USE_PLANCK_COSMOLOGY_) || defined (_USE_UNITSIM_COSMOLOGY_)
#define COSMOPARS Cosmo_parameters_PLANCK
#endif
#ifdef _USE_MINERVA_COMSOLOGY_
#define COSMOPARS Cosmo_parameters_Minerva
#endif
#ifdef _USE_MG_COSMOLOGY_
#define COSMOPARS Cosmo_parameters_MG
#endif
#ifdef _USE_ABACUS_COSMOLOGY_
#define COSMOPARS Cosmo_parameters_abacus
#endif
#ifdef _FLAGSHIP_
#define COSMOPARS Cosmo_parameters_Flagship
#endif
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

#define _WRITE_LOG_FILE_// not implemented yet

// ***************************************************PRECISION************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
//
// /**
// * @brief Define double precision for cosmicatlas
// * @details Applies all over the code except for gsl-type defined variables/containers
// */
//
// ********************************************************************************************************************************************
/**
*@brief Dynamical way of populating bins of density, MK, etc with the reference number of cells
*/
#define _DYNAMICAL_SAMPLING_  //optimal, faster

// ************************************************************************************************************************************
/**
 * @brief Define double precision for cosmicatlas
 * @details Applies all over the code except for gsl-type defined variables/containers
*/
//#define DOUBLE_PREC
//******************************************************
#ifndef DOUBLE_PREC
/**
 * @brief Define single precision for cosmicatlas
 * @details Applies all over the code except for gsl-type defined variables/containers
 * @details Defined when DOUBLE_PREC is undefined
*/
#define SINGLE_PREC
#endif
// ************************************************************************************************************************************
// Define Output type for binary arrays. If we work in double, we can ask to print bin files in double or float.
// If working with float, outputs will be by default in float
//#define OUTPUT_PREC_DOUBLE
#ifndef OUTPUT_PREC_DOUBLE
#define OUTPUT_PREC_FLOAT
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
/**
 * @brief This directive (if !defined) encloses any function that we wish to deprecate. Applied in the .cpp files
*/
#define _DEPRECATED_

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#define _USE_HEALPIX_

//#define _USE_HEALPIX_DNDZ_

/**
 * @brief Precision fopr Healpix
*/
#define healpix_real double
// ************************************************************************************************************************************
/**
 * @brief Precision for GSL
*/
#define gsl_real double  // Even with single precision, gsl will use double precision for all its calculations


// ********************************************
/**
 * @brief Define if the BiasMT_CATS are to be written
*/
//#define _WRITE_BiasMT_CATALOGS_
// ********************************************
/**
 * @brief Define if the BiasMT  cats are to be written in binary files. If undef, ascii files are written.
*/
//#define _WRITE_BINARY_BiasMT_FORMAT_
/**
 * @brief
*/
#define _WRITE_COORDINATES_
#define _APPLY_PERIODIC_BC_
// ********************************************
/**
 * @brief Define if the BiasMT  cats are to be written in binary files. If undef, ascii files are written.
*/
#define _WRITE_VELOCITIES_

// ********************************************
/**
 * @brief Define if the BiasMT  cats are to be read in binary files. If undef, ascii files are read.
*/
#ifdef mode_p
//#define _READ_BINARY_BiasMT_FORMAT_
#endif
//#define _READ_BINARY_BiasMT_FORMAT_
// ********************************************
/**
 * @brief // When defined, this allows for a raw metadata content, showing the name of the variables, its units.
*  @details Not yet working. Leave undefined.
*/
//#define _OUTPUT_WITH_HEADERS_
// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
#ifdef SINGLE_PREC
// ********************************************
/**
 * @brief Precision for FFT operations
*/
#define fftwf_real float
// ********************************************
/**
 * @brief Effective Precision for cosmicatlas
*/
#define real_prec float
// ********************************************
/**
 * @brief Effective Precision for complex containers used in FFTW
*/
#define complex_prec fftwf_complex
// ********************************************
// **************************************************************************************
// **************************************************************************************
#define _PREC_OUTPUT_ 6
#else
#define _PREC_OUTPUT_ 12
#define fftw_real double
#define real_prec double
#define complex_prec fftw_complex
#endif
// ********************************************
// **************************************************************************************
#ifdef _USE_SLICS_PK_BINNING_
#define DELTA_K static_cast<real_prec>(0.01)
#define MIN_K static_cast<real_prec>(0.0)
#endif
// ****************************************************************************************
#if defined _ABACUS_ || defined (_MG_)
#define BIN_SHIFT static_cast<real_prec>(0.5) 
#else
#define BIN_SHIFT static_cast<real_prec>(0.0)
#endif
// **************************************************************************************
// ********************************************
// Number of properties of a default halo catalog from BiasMT.
/**
 * @brief Number of properties to be allocated in a BiasMT_MOCK CAT
*/
#define N_PROP_CAT static_cast<int>(10)
#define MIN_N_PROP_CAT static_cast<int>(6)
// ********************************************
#define BIG_NUMBER static_cast<double>(1e7)
// ********************************************

#define ULONG unsigned long
// ********************************************
#define LONG long
// ********************************************
#define UULONG unsigned long long
// ********************************************
#define ASCII "ascii"
// *************************************************************
/**
 * @brief Precision type of the input par file of the DM field
*/
#define PrecType_X float
// ************************************************************************************************************************************
/**
 * @brief Precision type of the input par file of the Tracer field
*/
#define PrecType_Y float

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************CONSTANTS***************************************************************
// ************************************************************************************************************************************
// Define some constant factors
#define REAL 0
#define IMAG 1
#define ONE static_cast<int>(1)
#define dONE static_cast<double>(1.0)
#define KNOT_MAX_SIZE static_cast<ULONG>(100000000)
#define V_MAX_SIZE static_cast<ULONG>(100000000)
#define LARGE_NUMBER static_cast<double>(1e20)
#define MINUS_LARGE_NUMBER static_cast<double>(-1e30)
#define NOCELL static_cast<int>(-2)
#define NOVEL static_cast<int>(-999)
// ************************************************************************************************************************************
// Identifyers for cosmic-web environments
#define I_KNOT static_cast<int>(1)
#define I_FILAMENT static_cast<int>(2)
#define I_SHEET static_cast<int>(3)
#define I_VOID static_cast<int>(4)
// ************************************************************************************************************************************
// Identifyers for  mass assignment schemes
#define I_NGP static_cast<int>(0)
#define I_CIC static_cast<int>(1)
#define I_TSC static_cast<int>(2)
#define I_PCS static_cast<int>(3)
// ************************************************************************************************************************************
// Identifyers for system of coordinates
#define I_CART static_cast<int>(0) // x, y, z
#define I_EQR static_cast<int>(1)  // Ra, dec, r
#define I_EQZ static_cast<int>(2) // Ra, dec, redshift
#define I_EQRZ static_cast<int>(3) // R
// ************************************************************************************************************************************
// Identifyers for Lines of sight
#define I_LOSX static_cast<int>(1) // Along x
#define I_LOSY static_cast<int>(2) // Along y
#define I_LOSZ static_cast<int>(3) // Along z
// ************************************************************************************************************************************
#define NMIN_X_ONECELL static_cast<int>(0)
#define NMIN_Y_ONECELL static_cast<int>(0)
// ************************************************************************************************************************************
#define num_1   static_cast<double>(1.)
#define num_2   static_cast<double>(2.)
#define num_0_1 static_cast<double>(0.1)
#define num_0_5 static_cast<double>(0.5)
// ********************************************
#define VEL_BIAS_POWER static_cast<real_prec>(1.000)  // THis number is used in the power spectrum to multiply the velocitiys and compute P(k) in redshift space, used in FFTwfunctions
// ********************************************
// Number of mdoes used to get an estimate of the average bias on large
#define N_MODES static_cast<int>(30)
// Initial mode
#define N_MODE_INI static_cast<int>(1)
// Initial index (0..Nft) used to measure the residuals from INITIAL_MODE_RESIDUALS up to Nft/2
#define INITIAL_MODE_RESIDUALS static_cast<int>(1)
// ********************************************
// NUmbers used in the class measuring power spectrum
//#define CHUNK 27
#define ic_rank static_cast<int>(3) // For FFTW in 3 dimensions
#define NUMBER_WEIGHTS static_cast<int>(4)
#define ZERO static_cast<int>(0)
// ********************************************
#define TOLERANCE_FACTOR static_cast<int>(10)
// ************************************************************************************************************************************
// ************************************************************FUNCTIONS***************************************************************
// ************************************************************************************************************************************
/**
*@brief Absolute value
*@return  Absolute value of variable x
*/
#define myfabs(x) (*((int*)&x + 1) &= 0x7FFFFFFF)
// *************************************
/**
*@brief Bias model
*@return Bias model to pupulate DM density fields with overdensity x
*/
#define bias_test(x,alpha, rhoep, ep) (pow(1+x, alpha)*exp(-pow( (1+x)/rhoep, ep)))
// *************************************
/**
*@brief |r|
*@return  modulus of vector with coordinates x, y, z
*/
#define _get_modulo(x,y,z) static_cast<real_prec>(sqrt(x*x+y*y+z*z))
// *************************************
/**
*@brief 100000|r|
*@return  modulus squared of vector with coordinates x, y, z multiplied by a large number
*/
#define _get_modulo_squared(x,y,z)  static_cast<ULONG>(100000.0*(x*x+y*y+z*z))
// *************************************
/**
*@brief |r|
*@return  modulus squared of vector with coordinates x, y, z
*/
#define _get_modulo_squared_n(x,y,z) static_cast<real_prec>(x*x+y*y+z*z)
// ************************************************************************************************************************************
/**
*@brief Mmax (sigma_v) for dark matter halos of the slics. 
*@return  For a given x=sigma_v, halos cannot have masses above this limit
*/
#define limit_mass_halos(x) (1.5e11*pow(x,1.1)+1.5e11*pow(x,2.8))
// ************************************************************************************************************************************
/**
*@brief x sigma_v(Mv) for dark matter halos of the slics.
*@return  For a given x=M, halos cannot have sigma_v below this limit
*/
#define limit_vmax_halos(x) ((7.02e-6)*pow(x,0.405)+(1.1e-5)*pow(x,0.415))
// ************************************************************************************************************************************
// ********************************************************GENERAL STUFF AND OMP*******************************************************
// ************************************************************************************************************************************
/**
*@brief This allows the binning function to place in the extremes values outside the ranges
*/
#define _BIN_ACCUMULATE_
// ************************************************************************************************************************************
//#define _SHOW_ISSUES_
// ************************************************************************************************************************************
// *************************************************************
/**
* @brief Parallelism with OMP.
* @details Define to test parallel regions. Once these are properñly paralelized, change the def there to _USE_OMP_
* @details  Used now in %%new_new not working in the first loop of that function when thereshold_randoms is used
*/
#ifdef _USE_OMP_
#undef _USE_OMP_TEST_ 
#endif
// *************************************************************
/**
* @brief Parallelism with OMP.
* @details Define to test parallel regions. Once these are properñly paralelized, change the def there to _USE_OMP_
* @details  // Used now in %%new_new in assignmetn of 2nd-type prop, e.g. mvir, spin, rs. Working
*/
#ifdef _USE_OMP_
#define _USE_OMP_TEST2_
#endif
// *************************************************************
/**
* @brief Parallelism with OMP.
* @details Define to test parallel regions in the Catalog class. Once these are properñly paralelized, change the def there to _USE_OMP_
* @details So far the regions under this definition are giving problems, for these contains arrays with indices being updated from different threads. Leave undefineds
* @default Leave uncomment. Some loops cannot not be easily parallelized. This situation is also found in BiasMT:makecat()
*/
#ifdef _USE_OMP_
#undef _USE_OMP_TEST_CAT_
#endif
// *************************************************************
#ifdef _USE_OMP_
#define OMPPARRAN
#define OMPPARRANRSD
#define OMPPARGAR
#endif
// ************************************************************************************************************************************
//#define _ADD_NOISE_
#ifdef _ADD_NOISE_
#define MEAN_NEW static_cast<real_prec>(1000.0)
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************LPT ****************************************************************************

/**
* @brief Use Phase Space mapping withinn LPT to generte DM density fields.
* @details Leave commented. It provides much better results in terms of the filamentary struture of the LSS.
*/
#define _USE_TET_
// ********************************************
/**
 * @brief Use LPT
 * @details Define if we want to produce DM fields with LPT. If undef, the input DM file must be provided in the parameter file.
*/
//#undef _USE_LPT_
//#define _USE_LPT_

// ****************************************NAE****
/**
 * @brief To apply BiasMT to the displacements.
*/
//#define _DISPLACEMENTS_
//#undef _DISPLACEMENTS_
// ********************************************
/**
 * @brief Get density fields form LPTT
*/
#ifdef _DISPLACEMENTS_
#define DISP_COORD 3
#define _WRITE_DISPLACEMENTS_
#endif
// This is used to multiply the delta_ic if not it is extrapolated to z=0 from unknnown z
#ifdef _DISPLACEMENTS_
#define FACTOR_IC 1
#else
#define FACTOR_IC 1
#endif
// If this directive is defined, we are to define LPT and DO_BiasMT_CALIBRATION_
// ********************************************
// If the displacements are to be used, we do not use LPT as preparation, we use it inside the calibration.
#ifdef _DISPLACEMENTS_
#undef _USE_LPT_
#endif
// ********************************************
#ifdef _USE_LPT_
/**
 * @brief Define if only a simulation with LPT is to be run. The code stops before BiasMT starts doing the
*/
//#undef _ONLY_LPT_
#define _ONLY_LPT_
// ********************************************
/**
 * @brief Get the interpolated velocity field for LPT
*/
//#define _GET_VELOCITIES_LPT_
// ********************************************
/**
 * @brief Get density fields form LPTT
*/
#define _GET_DENSITY_FIELDS_
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************************MAIN MODE OF  BiasMT CODE********************************************************
// ************************************************************************************************************************************
/**
 * @brief IF defined, BA uses DM (local bias).
 * @details This should be defined by default. Undef if the local density is not to be used in the bias.
	This only affects functions as get_bias and get_mock_grid. Leave always defined.
*/
#define _USE_DM_IN_BiasMT_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
/**
*@brief This is defiend in the case when, in Catalog::read_input_bin, several .dat (binary) files are read as input in positions and velocities,
and the interpolation of density fields on a grid are perfoemd at while reading these files. This also affects the function, get_density_CIC and NGC in massFunctions
since we should not initialize the delta arrays there, as they are being filled in parallel with the reading.
*/
#undef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
//#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
/**
*@brief Compute the NGP interpolated density field from an input catalog
*/
#undef _GET_NGP_DENS_FIELD_
/**
*@brief Compute the CIC interpolated density field from an input catalog
*/
#define _GET_CIC_DENS_FIELD_
/**
*@brief Compute the TSC interpolated density field from an input catalog
*/
#undef _GET_TSC_DENS_FIELD_
#undef _GET_VEL_FIELD_
//#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *******************************************************REFERENCE TRACER CATALOG ****************************************************
// ************************************************************************************************************************************
#ifdef MOCK_MODE
#ifndef _USE_LPT_
/*
*@brief  Use the reference catalog. 
*@details This is available under _DO_BiasMT_CALIBRATION_
 	and allows to generate the tracer density field from the refenrece.
	The reference catalog is in the Catalog_file variable of the param file. If undef, the code shall read TR_DENS_FIELD.txt
*/
//#define _READ_REF_CATALOG_
#endif
#endif
// *********************************************************************
/**
*@brief If defined: The code selects objects from the input ref catalog using a minimum MASS, even if the observable is other quantity such as the vmax
*/
#define _SET_CAT_WITH_MASS_CUT_
// *********************************************************************
#ifndef _SET_CAT_WITH_MASS_CUT_
// Define this when power spectrum (-m opiton) is to me measured in cuts of VMAX
// In this case undef  _USE_MASS_TRACERS_
// Recall to undef _SET_GLOBAL_MASS_CUT_ below to use the -m option
// For -c option, undefine, and let mass define the cut
#undef _SET_CAT_WITH_VMAX_CUT_
#undef _SET_CAT_WITH_RS_CUT_
#undef _SET_CAT_WITH_SPIN_CUT_
#endif

#if defined (_SET_CAT_WITH_MASS_CUT_) || defined (_SET_CAT_WITH_VMAX_CUT_) || defined (_SET_CAT_WITH_RS_CUT_) || defined (_SET_CAT_WITH_SPIN_CUT_)
#define _SET_CAT_WITH_CUT_
#endif

#define MASS_SCALE static_cast<double>(1e12)

// This allows to read (if already created from the ref cat) the mass density field of the tracers
// Define ALSO when _READ_REF_CATALOG_ is defined in order to allow the reading of the masses of the tracer catalog
#define _USE_MASS_TRACERS_

// ***********************************************

// Define if the information on the satellite fraction is encoded in the ref catalog and is to be used
//#define _USE_SAT_FRACTION_
// ***********************************************
// Define if the information of the tracer velocities iºs present and is to be used
#define _USE_VELOCITIES_TRACERS_
// ***********************************************
// ***********************************************
// Defining something as OBSERVABLE, asks the code
// to compute all possible diagnosis as a function of cuts or bins
// in that particular quantity.
// The user must indicate the position of the VMAX in the -ini file in the slot for i_mass
// and specify mins and max in logMMIN and logMMX
// ***********************************************
/**
 * @brief Use the tracer mass as first observable 
*/
#ifdef mode_p
#define _USE_MASS_AS_PRIMARY_OBSERVABLE_
#endif
// ***********************************************
/**
 * @brief Use the tracer Vmax as first observable.
 * @detauils If this is undefined, the code uses Mvir as primary observavble
*/
#define _USE_VMAX_AS_PRIMARY_OBSERVABLE_
// ***********************************************
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#define _USE_VMAX_TRACERS_
#endif
// ***********************************************
#ifndef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#define _USE_MASS_AS_PRIMARY_OBSERVABLE_
#endif
// ***********************************************
/**
*@brief  When defined, BiasMT gets the vmax-mass scaling relation P(Mass|Vmax) or P(Vmax|Mass) and uses it to assign masses/vmax
 *once vmax/mass has been assigned.
*@details  In the end we obtain the right vmax-mass relation, the correct vmax function and
 3% residuals (over all mass range) in the mass function. Doing this allows us to increase the number of vmax bins
 as is explicitely shown below in the definition of N_VMAX_BINS
 Inside this definition we can add delta and use P(MKvir|Vmax, delta) or I/T web dependencies to use
 P(Mass|Vmax, {\theta})
*/
#define test_vmax_mass
// ***********************************************
#ifdef test_vmax_mass
/** 
*@brief When defined, add the DM information to obtain the mass from the conditional probability distribution P(Mass|Vmax, delta)
*/
#define _add_dm_density_
// ***********************************************
/** 
*@brief When getting P(M|V,delta) this factor increases NX (rad from parameter file) to EXTRA_BINS*NX the number of dm bins used.
*@details Use it >1 when using e.g IWeb with NX=100 for vmax assignment. INstnaces of TWEB use to have NX=500, in that case set this to 1 
*/
#define EXTRA_DM_BINS  static_cast<int>(1)
#endif //endif for test_vmax_mass

// ***********************************************
// ***********************************************
// ***********************************************
// ***********************************************
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
#define _ASSIGN_VMAX_POST_
#ifndef N_VMAX_BINS
#define N_VMAX_BINS static_cast<ULONG>(400)
#endif
#endif

#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#define _ASSIGN_MASS_POST_
#endif

// ***********************************************
// ***********************************************
#ifdef _USE_MASS_TRACERS_
// ---------------------------------
/**
*@brief This option use the masses of halos in cells to assign other properties. Leave undef.
*/
#undef _USE_MASS_FIELD_
// ---------------------------------
/**
*@brief This iption speeds calcuilations onder tests with collapse of randoms. Leave undef
*/
#undef _test_mass_assign_
// ---------------------------------
#endif // end of ifdef _USE_MASS_TRACERS_
// ***********************************************
/**
*@brief This option allows the get the dist. of min separations in reference
*/
//#define _GET_DIST_MIN_SEP_REF_
// ***********************************************
// ***********************************************
/**
*@brief This option allows the get the dist. of min separations in mock.
*/
//#define  _GET_DIST_MIN_SEP_MOCK_
// ***********************************************
// ***********************************************
/**
*@brief Use the Log10 of the trcer property under study. Default: undef
*/
#undef _USE_LOG_MASS_
// ***********************************************

// Maximum of the property (log, or sqrt, or whatever is selected above)
#define MAX_PROPERTY 16
#define _COUNTS_ "COUNTS"
#define _DENSITY_ "DENSITY"
#define _SAT_FRACTION_ "sat_fraction"
#define _TRACER_ "TRACER"

// ---------------------------------
// Please fix this
#define _MASS_ "_MASS_"
#define _MASS_LOG_   //ueful to measaure mass (or vmax) function
#define MBINS  /// 
// ---------------------------------
#define _VMAX_ "_VMAX_"
// ---------------------------------
#define _RS_ "_RS_"
//#define _USE_RS_AS_DERIVED_OBSERVABLE_
#define _RS_LOG_
#define RSBINS
#define _USE_RS_TRACERS_
// ---------------------------------
#define _RVIR_ "_RVIR_"
// ---------------------------------
#define _CONCENTRATION_ "_CONCENTRATION_"
#ifndef _USE_RS_AS_DERIVED_OBSERVABLE_
#define _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
#define _CONCENTRATION_LOG_
#define CONCENTRATIONBINS
#endif
// ---------------------------------
#define _SPIN_ "_SPIN_"
#define _USE_SPIN_AS_DERIVED_OBSERVABLE_
#define _SPIN_LOG_
#define SPINBINS
#define _USE_SPIN_TRACERS_
// ---------------------------------
#define _SPIN_BULLOCK_ "_SPIN_BULLOCK_"
// ---------------------------------
#define _VRMS_ "_VRSM_"
// ---------------------------------
#define _VIRIAL_ "_VIRIAL_"

// If we waant to reconstruct another property, we have to set in params files min and max of that propert. Define mins and max in Params.h,cpp
// and add ifs in BiasMT. in function get_X_function_complement.


// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ***********************************************************GENERATION OF MOCK CATALOG **********************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#ifdef _GET_BiasMT_REALIZATIONS_
/**
* @brief  This option enables the use of >1 sets of kernel+bias at generating independent halo number counts.
* @details If _UNITSIM_ or MG is defined, the code uses the set normal+pair. If SLICS (or any other GRFIC SIM) the code uses the refs from input parameter file
*/
#define _USE_TWO_REFS_MOCKS_
// *********************************
/**
* @brief  This option enabels the code to start the prcedure to generate a halo catalog
*/
#define _GET_BiasMT_CAT_
#endif
// ***********************************************
// ***********************************************
/**
* @brief If defined, the code prints the correlation betweeen the halo mass and the properties used to measure the conditional mass functio (or make the list of halo masses in a theta bin)
* @default Undefined
*/
//#defined _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
#undef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
// ***********************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************************PROPERTY ASSIGNMENT*********************************************************
// ************************************************************************************************************************************
#ifdef _GET_BiasMT_CAT_
// ********************************
#if defined _USE_MASS_TRACERS_ || defined _USE_VMAX_TRACERS_
/**
*@brief This option enables the assignment of halo properties.
*@details The ifdefs above only ensures that the algorithm has information from the tracer catalog to learn from//
*/
#define _ASSIGN_PROPERTY_
#endif
// ********************************
/**
* @brief Assign OORDS and properties to the calibration obtained after S1.
* @details In this case, only coords and vels can be assigned to the calibrated halo number counts if ASSIGN_PROPERTY is undef
* @details Undef in case of assignment to mock or assignmet to reference or new reference
*/
//#define _ASSIGN_TO_CALIBRATION_
#endif // end of ifdef _GET_BiasMT_CAT_

// ********************************
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Compute the bias object to object.
* @details Called from BiasMT methods. Implemented as a PowerSpectrumF method.
* @details This uses the slot used also by  INVARIANT_SHEAR_VFIELD_III_, _USE_S3_. Used when no hybrid is useful
*/
//#define _USE_BIAS_OBJECT_TO_OBJECT_ 
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Compute the bias object to object.
* @details Called from BiasMT methods. Implemented as a PowerSpectrumF method.
* @details This uses the slot used also by  INVARIANT_SHEAR_VFIELD_III_, _USE_S3_. Used when no hybrid is useful
*/
//#define _USE_BIAS_OBJECT_TO_OBJECT_ 
///////////////////////////////////////////////////////////////////////////////////////////////////
#define _USE_MACH_NUMBER_
#define _USE_LOCAL_OVERDENSITY_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// INFORMATION USED TO ASSIGN PROPERTIES:            NEWIGHBOURS
// ************************************************************************************************************************************
// ************************************************************************************************************************************
/**
* @brief  Use neighbours as a proxy for propèrty assignment.
* @details This option allows for using NUMBER OF NEIGHBOURS, MINIMIM DISTANCE TO NEIGHBOURS AND LOCAL CLUSTERING.
  These are three option are to be used disjointly.
*/
//#define _USE_NEIGHBOURS_
// ***********************************************
// ***********************************************
/**
* @brief This number controls how far from the current cell the code must search for neighbours
* @details This is also used outside the realms of the directive _USE_NEIGHBOURS_
*/
#define N_CELLS_BACK_FORTH static_cast<int>(1)
//***************

#ifdef _USE_NEIGHBOURS_
//***************
/**
* @brief  Use the number of neighbours as a proxy for propèrty assignment.
*/
//#define _USE_NUMBER_OF_NEIGHBOURS_   // This is in conflict with the invariants of the tidal field. Will move this slot to be shared with one of the shear-invariants
//***************
/**
* @brief  Number of close tracers within a range SCALE_MAX_N_NEIGHBOUR
*/
//#define _USE_MIN_DISTANCE_TO_NEIGHBOURS_
//***************
/**
* @brief // Information for local clusering. In the end it is the same as using the number of neighbours
*/
//#define _USE_LOCAL_CLUSTERING_
//***************

//***************
#define MIN_LOCAL_CLUSTERING static_cast<real_prec>(0)
//***************
#define MAX_LOCAL_CLUSTERING static_cast<real_prec>(1.2)
//***************

// This number is used for _USE_MIN_DISTANCE_TO_NEIGHBOURS_ and _USE_LOCAL_CLUSTERING_
// these numbers have to be tuned, e.g for UNITSIM (N_NEIGH_MAX,SCALE_MAX_N_NEI)=(62,8), (28,4)
// Also we can fix N_NEIGHBOURS_MAX and if a tracer happens to have more neighb, put them in the last bin
// When local clustering is used, SCALE_MAX-- =8Mpc
// N_NEIGHBOURS_MAX acts as the number of bins when the number of is used
// and takes the slot of I_C1

#define N_BINS_MIN_DIST_TO_NEI static_cast<int>(200)
//***************
#define N_NEIGHBOURS_MAX static_cast<int>(10)
//***************
#undef _USE_LOG_DIST_
//***************
#ifdef _USE_LOG_DIST_
#define MAX_OF_MIN_SEPARATION static_cast<real_prec>(log10(20.0))
#define MIN_OF_MIN_SEPARATION static_cast<real_prec>(-2)
#else
#define MAX_OF_MIN_SEPARATION static_cast<real_prec>(15.0)
#define MIN_OF_MIN_SEPARATION static_cast<real_prec>(0)
#endif //end of ifdef _USE_LOG_DIST_
// ***********************************************
#endif  //end of ifdef _USE_NEIGHBOURS_
// ***********************************************
// ***********************************************
/**
* @brief This number defines the maximum separations to count neighbours
* @details This is also used outside the realms of the directive _USE_NEIGHBOURS_ in -c options
*/
#define SCALE_MAX_N_NEIGHBOURS static_cast<real_prec>(8.0)
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************SEPARATION WITHIN CELLS ***********************************************************
// ************************************************************************************************************************************
/**
* @brief Use minimum separation in cells as a proxy for property assignments
* @details When defined for mocks results are not that good. High Vmax are spoiled.
*/
//#define _USE_MIN_SEPARATIONS_IN_CELLS_
// ***********************************************
/**
* @brief Number of bins in the dist of min (or mean) separations in cells
* @details This is also used outside the realms of the directive  _USE_MIN_SEPARATIONS_IN_CELLS_ in -c options
*/
#define N_BINS_MIN_SEP_IN_CELLS  static_cast<ULONG>(8)  // use 20 for tests, 100 for optimal mass assignment
// ***********************************************
#ifndef _USE_MIN_SEPARATIONS_IN_CELLS_
// ***********************************************
/**
* @brief Use mean separation in cells as a proxy for property assignments
* @details When defined for mocks results are not that good. High Vmax are spoiled.
*/
//#define _USE_MEAN_SEPARATIONS_IN_CELLS_
#define MIN_MEAN_SEP_IN_CELLS static_cast<real_prec>(0) // set 0 y linear, -2 if log This is used when mass is assigned to mocks.
#define MAX_MEAN_SEP_IN_CELLS static_cast<real_prec>(10.0)
#define N_BINS_MEAN_SEP_IN_CELLS  static_cast<int>(2)  // use 20 for tests, 100 for optimal mass assignment
#define DELTA_MEAN_SEP  static_cast<real_prec>(MAX_MEAN_SEP_IN_CELLS-MIN_MEAN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MEAN_SEP_IN_CELLS)
// ***********************************************
// ***********************************************
/**
* @brief Use stdv (in Mpc/h, linear scale) of the distributions of eparation in cells as a proxy for property assignments
* @details When defined for mocks results are not that good. High Vmax are spoiled.
*/
//#define _USE_STDV_SEPARATIONS_IN_CELLS_  // this is always linear saling (not log) by default
#define MIN_STDV_SEP_IN_CELLS static_cast<real_prec>(0) // set 0 y linear, -2 if log This is used when mass is assigned to mocks.
#define MAX_STDV_SEP_IN_CELLS static_cast<real_prec>(10.0)
#define N_BINS_STDV_SEP_IN_CELLS  static_cast<int>(10)  // use 20 for tests, 100 for optimal mass assignment
#define DELTA_STDV_SEP  static_cast<real_prec>(MAX_STDV_SEP_IN_CELLS-MIN_STDV_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_STDV_SEP_IN_CELLS)
// ***********************************************
#endif // endif for #ifndef _USE_MIN_SEPARATIONS_IN_CELLS_
// ***********************************************
/**
*@brief Define if eparatios in ells are binned in log
*/
//#define _USE_LOG_MSIC_
// ***********************************************
// This information goes into the slot for I_C_BIN3
#ifdef _USE_LOG_MSIC_
#define MIN_SEP_IN_CELLS static_cast<real_prec>(-0.3) // set 0 y linear, -2 if log This is used when mass is assigned to mocks. This depends on REDSHIFT!!!,it is GRAVITY IN ACTION
#define MAX_SEP_IN_CELLS static_cast<real_prec>(log10(7.0))
#else
#define MIN_SEP_IN_CELLS static_cast<real_prec>(0) // set 0 y linear, -2 if log This is used when mass is assigned to mocks.
#define MAX_SEP_IN_CELLS static_cast<real_prec>(10.0)
#endif //endif for #ifdef _USE_LOG_MSIC_
// ***********************************************
#define DELTA_MIN_SEP  static_cast<real_prec>(MAX_SEP_IN_CELLS-MIN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS)
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************************REDUCED MASS **********************************************************
// ************************************************************************************************************************************
/**
*@brief Enables ot use the reduced mass of pairs. Not useful.
*/
#undef _GET_DIST_REDUCED_MASS_
// ***********************************************
/**
*@brief Parameters for the reduced mass distribution
*@details Left outside the definition of  _GET_DIST_REDUCED_MASS_ since these are also used in -c mode
*/
#define N_BINS_REDUCED_MASS static_cast<ULONG>(200)  // use 20 for tests, 100 for optimal mass assignment
#define MAX_REDUCED_MASS static_cast<ULONG>(300)  // use 20 for tests, 100 for optimal mass assignment
#define MIN_REDUCED_MASS static_cast<ULONG>(100)  // use 20 for tests, 100 for optimal mass assignment
#define DELTA_REDUCED_MASS    (MAX_REDUCED_MASS-MIN_REDUCED_MASS)/static_cast<real_prec>(N_BINS_REDUCED_MASS)   //200 is roughly th maximum of mu in the sample
#define N_BINS_DIST_EX static_cast<ULONG>(10)   //NUmber of bis in separation to get dist of reduced masses
// ************************************************************************************************************************************
// ************************************************************************************************************************************
/**
* @brief Use cross correlation Halos-DM in Conf space
* @details THis computes the cross corr bethween halo number count sna DM density in FOurier space
* @details transforms back to conf space, make a uni-trnasform and use it to assign properties
* @warning
*/
//#define _USE_CROSS_CORRELATION_CONF_SPACE_
//#undef _USE_CROSS_CORRELATION_CONF_SPACE_

// ************************************************************************************************************************************
// *****************************************************TRACERS IN CELLS **************************************************************
// ************************************************************************************************************************************
/**
* @brief  Use the number of tracers in cells as a proxy for primary property assigment.
* @warning This is in potential conflict with the invariants of tnhe velocity field, I.E, this can be defined unless those invariants are also used.
*/
//#define _USE_TRACERS_IN_CELLS_
// *************************
#ifdef _USE_TRACERS_IN_CELLS_
/**
* @brief This is the maximumn number of tracers in a cell used for assignment.
*/
#define N_TRACERS_IN_CELLS_MAX  static_cast<int>(36)
// *************************
/**
* @brief Number of bins in wihch the "number of tracers in cells" is to be representd for BiasMT
*/
#define N_BINS_TRACERS_IN_CELLS  static_cast<int>(6)
// *************************
/**
* @brief Delta for this cquantity
*/
#define DELTA_TRACERS_IN_CELLS static_cast<double>(N_TRACERS_IN_CELLS_MAX)/static_cast<real_prec>(N_BINS_TRACERS_IN_CELLS)
#endif // end _USE_TRACERS_IN_CELLS_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *****************************************************TOTAL MASS IN CELLS ***********************************************************
// ************************************************************************************************************************************
/**
* @brief Define if the total property (mass or vmax) in cells is to be used property assignment.
* @default Undefined
* @warning Not yet used. Have to be checkd.
*/
#undef _USE_TOTAL_MASS_IN_CELL_ // this causes problems with the mock.
// *************************
#ifdef _USE_TOTAL_MASS_IN_CELL_
// *************************
/**
* @brief Numbers of bins in the
* @warning Not yet used. Have to be checkd.
*/
#define N_BINS_TOTAL_MASS_IN_CELL  static_cast<int>(20)
// *************************
/**
* @brief mAXIMUM VALUE USED FOR TOTAL MASS IN CELL
* @warning Not yet used. Have to be checkd.
*/
#define MAX_TOTAL_MASS_IN_CELL static_cast<real_prec>(1e15)
#endif  // endif for _USE_TOTAL_MASS_IN_CELL_

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************************TESTING ********************************************************************
// Define this in case we want to reassign properties to the reference catalog, in
// order to see whether the bias, from mass cuts, is still that from the reference
// or if the mass assignment induces the wrong clustering
// This option does not allow to use any probabilistic approach yet
//********************************************************
//********************************************************
#ifdef _ASSIGN_PROPERTY_
// *************************
/**
* @brief Define whether more than one reference are used to learn from in the procedure of property assignment
* @details If _UNITSIM_ or MG is defined, the code uses the set normal+pair. If SLICS (or any other GRFIC SIM) the code uses the refs from input parameter file
* @details If unabled,the code uses only one refernece poinetd by the input parameter file
* @details This option might be suitable for cases in which multiscale is used
*/
//#define _USE_TWO_REFS_MOCKS_ASSIGNMENT_
// *************************
/**
* @brief If properties are to be assigned, define whether to assig to mock or to reference
* @details Assigning to reference aims at checking the 1 dn 2 statistics of reference when propas has been reassigned
* @details The tracer coordinates and velocities for the mock are taken from the refernece catalog used to learn from
* @warning If >1 references are used to learn from, the core reads the coords and vels from the last read catalog
*/
#ifndef _ASSIGN_TO_CALIBRATION_
#define _ASSIGN_PROPERTIES_TO_REFERENCE_
#endif
// *************************
/**
* @brief define this if the kernesl obtained from the calibration is to be used in the assignment
* @default Undefined
*/
//#define _KONVOLVE_PASSIGN_
#endif // end if for _ASSIGN_PROPERTY_
// *********************
/**
@brief define this if the kernesl obtained from the calibration is to be used in number count generation
    @default Defined, it has to be defined always!!
*/
#define KONVOLVE_TO_MOCK
//********************************************************
#if !defined (_ASSIGN_PROPERTIES_TO_REFERENCE_) && !defined (_ASSIGN_TO_CALIBRATION_)
/**
* @brief Enable the assignment of properties to BiasMT mock
* @details This is enabled only if no properties are to be assigne to ref or to calibration
* @details BiasMT-MOCK s understood as an independent set of tracers with phase-space coords
*/
#define _ASSIGN_PROPERTIES_TO_MOCK_
#endif
//********************************************************
//********************************************************
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
/**
* @brief We can assign properties to another reference
* @details The other reference is a simulatin with right coords, provided in parameter file
* @details If undefied, properties are assigned to references or mock, provided whether _ASSIGN_PROPERTIES_TO_REFERENCE_ is defiend or not.
*/
#define _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
#endif
//********************************************************
/**
 * @brief Assign primary property using the actual values from the reference tracer catalog(s).
 * @details This option uses the method  BiasMT::assign_tracer_mass_new_new()
 * @details If undefined, the code uses the statistical assignement based on mass bins via BiasMT::assign_tracer_mass_new()
 * @remarks the function  BiasMT::assign_tracer_mass_new() is marginally sensitive to the bins of mass in which the mass function os measured
 * @remarks The function  BiasMT::assign_tracer_mass_new_new() is still less sensitive for the latter keeps the list of masses in a bin of DM properties.
 * @warning This is planned to be the standard assignment. When other properties are to be assigned, we shall use a probabilistic approach
 * @default Defined
*/
#ifdef _ASSIGN_PROPERTY_
#define _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
#endif
//********************************************************
#ifdef _ASSIGN_PROPERTY_
/**
 * @brief When enabled, the DM is interpolated at the position of the halos and then binned.
*/
//#define _USE_DM_DENSITY_AT_HALO_POSITION_
#endif
//********************************************************
//********************************************************
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
#ifdef _MULTISCALE_
#define _JOIN_AND_SORT_
#endif // end of MULTISCALE
#endif // end of _USE_TWO_REFS_MOCKS_ASSIGNMENT_
// ************************************************************************************************************************************
// ***************************************************COLLAPSE RANDOMS ****************************************************************
// IF ASSIGN TO REFERENCE is defined, we do not need to use randoms neigher to collapse them
#if defined _ASSIGN_PROPERTIES_TO_MOCK_ || defined (_ASSIGN_TO_CALIBRATION_)
//****************************************
//****************************************
/**
* @brief Collapse randoms towards DM particles
  @details This option is implemented *after assigning* halo properties. This shouled be standard for two main reasons:
  first, allows for some property-dependent collapse, and second, dies not change the number count of the mock number field.
  Although this is not important for the independnt DM fields (as long as we use the tolerance factors)
  when testing the assigment with the field produced by the calibration we might not be able to get the 0% error in the vmax function
  because of collapsing before collecting the vmax values in the bins of theta.
* @default Use this option by default.
*/
//#define _COLLAPSE_RANDOMS_
//#define _COLLAPSE_RANDOMS_VELS_

#ifndef _COLLAPSE_RANDOMS_
/**
* @brief Collapse randoms towards DM particles before property assignment, at the same time when positions of DM are assigned.
  @details Undefined. Conditional to _COLLAPSE_RANDOMS_. Does not use any tracer property
*/
#define _COLLAPSE_RANDOMS_AUX_
#endif
//********************************************************
//********************************************************
#ifdef _COLLAPSE_RANDOMS_
/**
* @brief Collapse randoms using the information of the particle properties (masses) to compute radii and mimic exclusion 
*/
#define _COLLAPSE_RANDOMS_USING_EXCLUSION_
// ------------------------------------------------------------------------------------
/**
* @brief Search for closest DM particle to random in neighbour cells
*/
#undef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_

#endif  // end  _COLLAPSE_RANDOMS_
// ------------------------------------------------------------------------------------
// ************************************************************************************************************************************
/**
* @brief Apply a random dispersion to the poisitons of high vmax (PROP LEVEL) to mimic exclusion.
* @details This is performed right before the codes writes the catalog to a file, so this operation might alter the fiducial bias and numer counts.
*/
//#define _APPLY_GLOBAL_EXCLUSION_

// ************************************************************************************************************************************

#endif  // end  of assign_properties_to_mock

// ***********************************************
// ***********************************************
/**
*@brief This option assigns vmax from the reference and the mass linked to the vmax chosen. Default: undefind
*/
//#define _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
// ***********************************************
// ***********************************************
// ***********************************************
// ***********************************************
// If we do not assign properties to reference, but instead we use all other possible cases
// (i.e, assign to mock in general, to ref statistically, to other refs in general)
// then we need to define a statistical assignment to enable the definitiona nd allocation of containers type ABUNDANCE
#if defined _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_ && defined _ASSIGN_PROPERTIES_TO_REFERENCE_
#undef _USE_STATISTICAL_ASSIGNMENT_
#else
#define _USE_STATISTICAL_ASSIGNMENT_
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
/**
 * @brief Assign properties (vmax, or mass) using multiscale.
 * @details Multiscaling works only when reading the actual primary property values from the reference(s)
 * @note The procedure consists in a loop to identify cells with objects above the X (X=mass, Vmax) threshold Mth, imposing already positions for the masses
 * (we can avoid that threshold here, and then we would be counting all cells). Note that counter_high_mass and N_props_multiscale are merely consequence of the X (vmax or mass) function
 * The number of tracers in each level is computed from the reference.
 * With the number of tracers in each interval (level), we can pass to assign the largest X values onto the low resolution mesh.
 * The fact that we assig the larger X-values is represented in the ordering imposed in the bins of theta, such that by increasing the index
 * every time we assign a mass we effectively assign top-bottom in each theta-bin.
 * Note that this applies for vmax and Mass, as these display a positive correlation. If a halo property displays an anticorrelation (but still strong)
 * with the underlying dm field and hence with the mass, the assignment would have to be done from bottom-to-top.
 * @note Regimes of multi-scale mass assignment are identified by a number of tracers in a given interval and the
 * correspnding grid size for that interval.
 * @note We need to play witht the threshold of tolerances given in the parameter file.
 * The higher the number of proerties/bins we uses to allocate the reference properties, the harder is for
 * the algoritm to assing properties, even if it assigns to the same reference catalog.
 * Hence, even in that ideal case, some tolerance threshold is advised, both from the parameter file as from preproc-directives
*/
#define _MULTISCALE_
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#ifdef _GET_BiasMT_CAT_
// ***********************************************
#ifdef _MULTISCALE_
#define  _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
//#define _NEW_APPROACH_ASS_
#endif
// ***********************************************
/**
* @brief Define to compute the property limits in multiscale from the Nft set provioded in parameter file
*/
#define _USE_FIXED_MULTISCALE_LEVELS_
// ***********************************************
// ***********************************************
#ifndef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
/**
*@brief Mulsti-scale property assignment.
*@details When defined, the code looks for the number Vmax from the cumulative
	such that vmax above this number is to be assigned to DM tracers. Vmax below
	this number will be assigned to random placed tracers.
	When this is defined, MULTISCALE is not defined.
*@default Undefined
*/
#undef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
// ***********************************************

#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
/**
* @brief When defined, the code assigns values of vmax<Vtrehsold to random particles. 
* @details This alone does not work. It works better though than the _BOTTOM_RANDOM_ option.
*/ 
//#define _BOTTOM_RANDOM_

#ifndef _BOTTOM_RANDOM_
/**
* @brief When defined, the code assigns valaues of vmax>Vtrehsold to random particles.
* @details This alone doe snot work. 
*/
//#define _TOP_RANDOM_
#endif
#endif // end _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#endif // end _USE_MULTISCALE_PROPERTY_ASSIGNMENT_
#endif // end multiscale
// ***********************************************
// ***********************************************
// Mass threshold: for mass es above thes values, masses are assigned in a grid-exclusion based approach. THis aims at accounting for exlcusion
// For masses below this values, masses are taken form the allowed lists, in a tracer-based approach. This neglects exclusion

// Mass threshold above which a low-resolution cell is used in order to assign masses
// If this number is larger than the maximum tracer mass, PLEASE undef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_

#define PROP_THRESHOLD static_cast<real_prec>(1e11)  //this is now mass threshold_leve1
// ***********************************************
/**
 * @note This tolerance factor is applied to the multilevel approach level 0 (particle as a funciton of theta).
 * The CPU time for the assignment of properties to generate mocks depends on these TOLERANCE factors
 * Low values yield high speed, though low accuracy in the clustering signal.
 * Also, the higher the number of properties/bins insed for the assignment, the more flexible we need to be with the tolerance
*/
#ifndef _MULTISCALE_
#define TOLERANCE_FACTOR_L0b static_cast<real_prec>(1.0)
#define TOLERANCE_FACTOR_L0 static_cast<real_prec>(0.97)
#else
#define TOLERANCE_FACTOR_L0 static_cast<real_prec>(0.93)
// This tolerance facor is applied to the multilevel approach level 0 (particle randomly selected)
#define TOLERANCE_FACTOR_L0b static_cast<int>(0.98)
#endif
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Define this if the glocbal mass funciton is to be used in the assignment of mass
// For a coprrect mass assignment, leave it undef
//#define _USE_GLOBAL_MASS_FUNCTION_
// This avoids the convolution of the DM field with the kernel, in order to see whether the convolution makes an efect on the precision of the assignment
// ************************************************************************************************************************************
// **************************************************BOUNDARIES FOR HISTOGRAMS *******************************************************
// This is important when assignning masses to tracers, in order to use all inforamtion.
//#define _MODIFY_LIMITS_
// ************************************************************************************************************************************
// **************************************************EXCLUSION ************************************************************************
//#define _CORRECT_FOR_EXCLUSION_   // Toooooooo slowwwwwwwww
// This is used in order to alloca the list of masses in pairs between the mini sep of the reference and MAXIMUM_DISTANCE_EXCLUSION
#define MAXIMUM_DISTANCE_EXCLUSION static_cast<real_prec>(7.0) // 1.5 Mpc/h as maximum separation explored??
#define MINIMUM_DISTANCE_EXCLUSION static_cast<real_prec>(0.1) // 1.5 Mpc/h as maximum separation explored??
#define EXCLUSION_SCALE static_cast<real_prec>(2.5) // Correct masses of pairs within this sparation
// ************************************************************************************************************************************
#define _VEL_UNITS_MPC_PER_h_    // define if the velocities in the final cat are to be given in Mpc/h
// ************************************************************************************************************************************
// **********************************************************POWER SPECTRUM **********************************************************
// ************************************************************************************************************************************
// Directives for the section of the code measuring Pk,(called with -m at compillation time)
// This definition has to work together with the one in the input parameter file.
// This only works on the method PowerSpectrumF::compute_power_spectrum(book, bool)
// to get power spectrum reading a catalog
//#define _CHECK_PARTICLES_IN_BOX_
//#define _USE_WEIGHTS_IN_POWER_              //This must be defined if weighs, as specified in the parameter file, are to be used in the determination of power spectrum
//#define _GET_BISPECTRUM_NUMBERS_   // This doe snot compute bispctrum, but allows to compute numbers used in the shot noise and normalziation of that statistics. Verify that this is de fined when that option is demanded from parfile
// DEFINE when the option -p is to be used
//If _POWER_ is defiened, the code distinguishes between the observed quantity and the defining (setting) quantity.
//The observed quantity is that with respect to which bins or cuts are to be done in order to measure power
// while the defining quantity is that global cut used to define the full sample.
// That is, we can have a halo catalog with mass above some mcut, and still want to make bins in Vmax
// If _POWER_ is undef, i.e, for BiasMT,  we shall assume that the catalog is DEFINED by the minimum halo MASS
#ifdef mode_p
#define _POWER_  
#endif
// ******************************************************************************************
// ******************************************************************************************
// These three disjoint options are to be used with the -p compilling flag.
#ifdef mode_b
#define _USE_ALL_PK_
#endif

// ******************************************************************************************
// ******************************************************************************************

#ifdef mode_p
#define _USE_MASS_BINS_PK_   //aplies for bins in the property specified by _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_
//#define _USE_MASS_CUTS_PK_
#endif
// ******************************************************************************************
// ******************************************************************************************
/**
* @brief Define this when measuring power from a grid with uts in number counts: P(k;N>=)1
*/
//#define _NCUTS_POWER_
#define N_MAX_OCCUPATION 20
// ********************************************
#ifdef _USE_MASS_BINS_PK_
// Define which property is to be binned when measuring pwoer spectrum. If no binning is desired,
#define _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_
//#define _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_0
//#define _USE_RS_AS_OBSERVABLE_POWER_
//#define _USE_SPIN_AS_OBSERVABLE_POWER_
#endif
// ********************************************
#define _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
//#define _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
//#define _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_   // This should be instead of using Mass or Vmax for the -u option. If this is defined, this wins over use_mass_as_pŕim_prop and in the analysis Mvir will be understood as nu.

// ********************************************

#ifdef _USE_MASS_CUTS_PK_
// Define which property is to be binned when measuring pwoer spectrum. If no binning is desired,
//#define _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_
#define _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_
//#define _USE_RS_AS_OBSERVABLE_POWER_
//#define _USE_SPIN_AS_OBSERVABLE_POWER_
#endif


// This is meant to measure power spectrum when passing the structure tracer.Halo
// This is due to the fact that when one reads the cat and measures P(k), one can pass the minimum cut
#ifdef _USE_MASS_CUTS_PK_
#define _SET_GLOBAL_MASS_CUT_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_
#define MINIMUM_PROP_CUT_original static_cast<real_prec>(2e13) // put the largest of the mcuts used in power, for comparison purposes
#define MINIMUM_PROP_CUT static_cast<real_prec>(pow(MINIMUM_PROP_CUT_original, exponent_mass_tracer))
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_
#define MINIMUM_PROP_CUT static_cast<real_prec>(500.0)
#endif
#endif

// Define if the 2D power spectrum in paralllel and perpendicular wavevectors is to be computed and written
//#define _WRITE2DPOWER_

// Define if the the power spectrum will be measured from redshift space in case we want to shift positions with known velocity field.
// If we deal with a survey in equatorial coordinates, RSD are alrady included, so leave undefined.

// *******************************************************************************************************************************************
/**
* @brief Define this when measuring power in redshift space
* @detail Define also when building mocks. BiasMT will measure real and redshift space
*/
#ifndef _ONLY_LPT_
#define _REDSHIFT_SPACE_
#endif
// *******************************************************************************************************************************************
/**
* @brief Identify the line of sight in order to move particles in redshift space
*/
#define LOS 3
// *******************************************************************************************************************************************
/**
* @brief 
*/
#ifdef _REDSHIFT_SPACE_
#define _WRITE_MULTIPOLES_
#endif
#define _WRITE_MULTIPOLES_
// *******************************************************************************************************************************************
/**
* @brief If defined, the code selects redshisft bin interval given in the input parameter file in order to compute the power spectrum.
*/
#ifdef _POWER_
//#define _USE_REDSHIFT_BINS_
#endif
// *******************************************************************************************************************************************

// Define this to use the vectorized version of grid assignment, written by L. Tornatore.
// That version has some memmory issues depending on the computer the code is run at.
// If undef, the code uses my own version which is  efficient enough
// IMPORTANT: when using the standar assignment, I demand here to specify whether PSC
// is to be used, in order to speed the code. I will set a warning in the parameter  file
//#define _USE_VECTORIZED_GRID_ASSIGNMENT_
#ifndef _USE_VECTORIZED_GRID_ASSIGNMENT_
#define MAX_MAS_DEG 5
#define MAX_MAS_DEG_TSC 3
#define MAX_MAS_DEG_CIC 3
#define CHUNK MAX_MAS_DEG*MAX_MAS_DEG*MAX_MAS_DEG   // this is 5*5*5
#endif


#define KMAX_RESIDUALS_high static_cast<real_prec>(0.5)
#define KMAX_RESIDUALS_low static_cast<real_prec>(0.3)


// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *************************************************BiasMT GENERAL ***********************************************************************
// ********************************************

// avoid writing or computing files already existing
// Not all files are under this def, COMPLETE!
//#define _VERIFY_FILES_

// ********************************************
// Change coordiantes x<->y in the input DM  density field
//#define _EXCHANGE_X_Y_DENSITY_
// ********************************************
// Change coordiantes x<->y in the input DM  density field
//#define _EXCHANGE_X_Y_DENSITY_
// ********************************************
// Change coordiantes x<->y in the input DM velocity fields
//#define _EXCHANGE_X_Y_VEL_

// ********************************************
// Compute and write to putput files PDF
//#define _WRITE_PDF_

// ********************************************
// Define if the TR density field is to be written (bin) during the iteration
//#define _WRITE_ALL_FIELDS_
#ifndef _WRITE_ALL_FIELDS_
#define _WRITE_TR_DENSITY_FIELD_
//#undef _WRITE_TR_DENSITY_FIELD_
#endif
// ********************************************
// Define if the DM density field is to be written (bin)
//#define _WRITE_DM_DENSITY_FIELD_

// ********************************************
// Define this if DM spectra duwing iterations is to be produced

#ifndef _ONLY_POST_PROC_
#define _GET_POWER_REFS_
#endif


#define _GET_POWER_FROM_CATS_

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************

// ********************************************
// Define if several realiations are used to compute bias statistics from. Narems in this case are hard-coded
//#undef _SEVERAL_REAL_
// ********************************************
#ifdef BIAS_MODE
// This isused in the context of chaning the reference DM when bias are below one
// situation in which BiasMT seems to be out of order.

#define NUM_IN_LOG static_cast<real_prec>(1.000001)
// Define this if the Joint B is to be computed from several realizations of the DM field (alpt so far) and the reference
//#define _SEVERAL_REAL_BIAS_
#endif

// ********************************************
// Compute smoothed version of the BiasMT kernel.
// Verified that smoothing is best for unitsim and abacus
#define _SMOOTHED_KERNEL_
// ********************************************
// Define if smoothed is to be done only in the last iteration
#ifdef _SMOOTHED_KERNEL_
/**
* @brief 
**/
//#define _SMOOTHED_KERNEL_LAST_ITERATION_
#endif
// ********************************************
// ********************************************
// This is the number ised in log(NuUM_IN_LOG+delta). If
// Undesirable situations with log10 are to be avoidded, I suggest this number to be set to 2.
// specially when mocks are to be built. If plots are to be done after using the -n option, set 1
// tHIS NUMBER MUST BE THE SAME IN CALIBRATION AND MOCK PRODUCTION.
#ifdef MOCK_MODE
#define NUM_IN_LOG static_cast<real_prec>(2.)
#endif
// --------------------------------------------
#ifdef MOCK_MODE
// Define to run the test in which a Halo density field is created from the input DM density
// field using a particular halo bias, given by bias_test
//#define _RUN_TEST_
#endif
// ********************************************************************************************************************************************
/**
* @brief 
**/
//#define _UNDER_BIASED_
// ********************************************************************************************************************************************
/**
	@brief  Do rank ordering of the DM field to the new_ref identified in parameter file
*/
// --------------------------------------------
/**
* @brief 
**/
//#define _RANK_ORDERING_AB_INITIO_
// --------------------------------------------
#ifdef _RANK_ORDERING_AB_INITIO_
/**
* @brief 
**/
#define _RO_WITH_DELTA_
#endif
// ********************************************************************************************************************************************
/**
* @brief 
**/
//#define _RANK_ORDERING_MOCK_GEN_ //  not working
/**
* @b/rief 
**/
  //#define _RO_WITH_DELTA_MOCK_GEN_  // not working
// ********************************************************************************************************************************************
#ifdef _ONLY_POST_PROC_
// This option was meant to avoid the convolution at the post_proc step.
// However the code collapse, even when assigning to the calibration, 
// even if the same dm field is used consistenly.
// So, it seems that: if you force the ref tracers to follow a convolved DF,, that is ok
// but if you force the mock to follow the original DM, wrong.
// What if you use the ref with the original dm and the mock_number count with the convolved one? Wrong.
//	#define _DO_NOT_CONVOLVE_
#endif
//IF this is defined, the calibration of the Y tracer is done with the delta_x and the delta_complement, where
// the complement is passed from the input parameter file
//#define _USE_X_COMPLEMENT_I_
// ********************************************************************************************************************************************
//#define file_one_cell this->params._Output_directory()+"nmax_cell.txt"
#ifdef MOCK_MODE
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#ifdef _GET_BiasMT_REALIZATIONS_
#ifndef file_one_cell
#define file_one_cell this->params._Input_Directory_BIAS_KERNEL()+"nmax_cell.txt"
#endif
#else
#define file_one_cell this->params._Output_directory()+"nmax_cell.txt"
#endif  // end of ifdef _GET_BiasMT_REALIZATIONS_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// If no mocks is to be created, the code assumes that we watn to perform the calibration of kernel.
// ********************************************
#ifdef _GET_BiasMT_REALIZATIONS_
/**
* @brief Post_processing for property assignment
* @details This is meant to be used when the halo number counts have been already computed and we are just to assign coordinates and  properties 
**/
#define _ONLY_POST_PROC_
#undef _EXTRAPOLATE_VOLUME_
#endif
#ifdef _DO_BiasMT_CALIBRATION_
//#define _CALIBRATION_WITHOUT_SN_
#endif
// THE OPION use_tracer_hr is allowed to be on only if we are doing the calibration
#ifndef _EXTRAPOLATE_VOLUME_
//This option is not working, The kernel, when checked with UNITsim, goes above 1 and1
// no convergece is acquired.
//#define _USE_TRACER_HR_ // to be deprecated
#endif
// ********************************************
#ifdef _DO_BiasMT_CALIBRATION_
// This is useful when calibration over different ref catalogs are to be done
// Note taht the paths to input files are hard coded in BiasMT.cpp
//#define _SEVERAL_REAL_CAL_
#endif
#endif   // end of get_BiasMT_realizations 791
// ********************************************
/**
* @brief This index is used at the def of ratio (Pref/Pnew)**KERNEL_INDEX
  @details By default this number should be 1.0
*/
#define KERNEL_INDEX static_cast<real_prec>(1.0)
// ********************************************
// ************************************************************************************************************************************
// ***********************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** COSMIC WEB ****************************************************************************
/**
* @brief Cosmic-Web classification
  @details Define this if the CWC is to be used as variable to characterize the bias.
 * This does not define though whether CWC is computed or not!!!
 * If undef, BiasMT could still do the classification in order to use the eigenvalues
 * of the Tidal field as bias variables. This happens if USE_INVARIANTS is defnied
 * RECALL TO SET UP THE VARIABLE v_CWT_used in the parameter file
 * in case this variable is defined
*/
// ********************************************
// ********************************************
#ifdef MOCK_MODE
/**
*@brief Perform the COmicWeb classification to be used in BiasMT.
*/
#if defined (_USE_TWEB_) || defined (_USE_CWEB_) || defined (_USE_TIWEB_)
/**
*@brief Perform the Cosmi-cWeb classification to be used in BiasMT.
*/
#define _USE_CWC_
#ifdef _USE_CWC_
//#define _USE_CURVATURE_  // This allows to use Tkmodels but based on the eigenvalue sof the curvature, istead of those form the tidal field.
#endif
// ********************************************
/**
*@brief Perform the COmicWeb classification to be used in BiasMT when assigning masses to halos.
*/
#define _USE_CWC_mass_
// ********************************************
#endif
#endif
// ********************************************
#define _USE_CWC_HALO_ANALYSIS_
#define _ASSEMBLY_BIAS_MASS_ONLY_
// ********************************************
// ***********************************************************************************************************************************
/**
* @brief Cosmic-Web classification inside iterations
  @details  If defined, the CWC is done in every iteration of the loop
 * If not, the CWC is done only in the first iteration, and that particular classification
 * is used in the iterative process. Preference: defined. This applies also to the Invariants
 * or any other quantity obtained from the gravitational potential
*/

#ifdef _USE_CWC_
#define _USE_CWC_INSIDE_LOOP_  //TO BE DEPREDCATED
#endif 
#define _USE_CWC_INSIDE_LOOP_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** COLLAPSING REGIONS ********************************************************************
// ********************************************
/**
* @brief  Use the mass of large-scale collapsing regions
* @details Note that in this case the CWC classification is done, even if _USE_CWC_ is not defined
*/
#ifdef MOCK_MODE

#if defined (_USE_TWEB_) || defined (_USE_MKNOTS_) || defined (_USE_IKWEB_) || defined (_USE_TIWEB_)
//#define _USE_MASS_KNOTS_   // This is used when asking inside the bias function for each of the components
//#define _hydro_
//#define _WRITE_MKNOTS_
#endif
#endif
// -----------------------------------------------------------
/**
* @brief  Max for the log of the mass of SuperKnots found by the FoF
* @details This number cchanges with the definition of lambda threshold.
* @details For lambda <0, more regions are classified as knots and hence MKMAX MUST BE INCRESAED e.g. to 1e7
*/
#ifdef _hydro_
#define _DM_NEW_UNITS_   //define this if youn have dm with cgs units. If a DM is done with counts, undef
#define MKMAX 1e-27
#define MKMIN 1e-30
#else
#ifdef _USE_CURVATURE_
#define MKMAX 5e1  // When using the curvature, the total mass of objects indetified as c-knots s lower
#else
#if defined _UNITSIM_ || defined _ABACUS_ || defined _TNG_ || defined _TNG_GAL_ || defined _FLAGSHIP_
#define MKMAX static_cast<real_prec>(2e5)  // I used 1e3 for UNITSIM and ABACUS (lambda=0) with ALPT. For the original DM, use 2e5
#elif defined _MG_
#define MKMAX static_cast<real_prec>(2e5)  // I used 1e3 for UNITSIM (lambda=0) and 2e7 for lambda=-0.3, 1e7 for -0.1
#endif
#endif
#if defined _UNITSIM_ || defined _ABACUS_ || defined _TNG_ || defined _TNG_GAL_ || defined _FLAGSHIP_
#define MKMIN static_cast<real_prec>(0.)
#elif defined _MG_
#define MKMIN static_cast<real_prec>(1.)
#endif
#endif
// -----------------------------------------------------------
/**
* @brief  MIN for the log of the mass of SuperKnots found by the FoF
*/
// -----------------------------------------------------------
/**
* @brief This is meant in case we have CIC DM fields and we want to transform it to NGP
*/
//#define _CONVERT_CIC_TO_NGP_
// ***********************************************************************************************************************************
// ***********************************************************************************************************************************
// ********************************************INVARIANTS OF TIDAL FIELD
// **********************************************************************************************************************************
/**
* @brief If this is defined, the code uses as invariant_II the value of the first eigenvalue,
*  as invariant_III -> second eigenvalue and invariant_IV -> third eigenvalue.
* This option is then used then tidal field eigenvalues are to be used.
*/
//#define _USE_EIGENVALUES_
// **********************************************************************************************************************************
#ifdef MOCK_MODE
/**
* @brief Define if invariants are t be written to binary files
*/
//#define _WRITE_INVARIANTS_
// **********************************************************************************************************************************
/**
* @brief This is only meant to write the invariant in a binary file.
  @details The invariant I is the same delta, so it is already there
*/

#ifdef _USE_IWEB_
//#define _USE_INVARIANT_TIDAL_FIELD_I_
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_I_
#define _USE_EXPONENT_INVARIANT_I_
#define _USE_SIGN_INVARIANT_I_
#define _MAP_TO_INTERVAL_INV_I_
#define NEWMAX_INV_I static_cast<real_prec>(1.0)
#define NEWMIN_INV_I static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_I_
#define EXPONENT_INVARIANT_I  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_I  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
/**
* @brief Invariant II of the tidal field
*/
#ifdef _USE_IWEB_
#define _USE_INVARIANT_TIDAL_FIELD_II_
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_II_)
#define _USE_EXPONENT_INVARIANT_II_
#define _USE_SIGN_INVARIANT_II_
#define _MAP_TO_INTERVAL_INV_II_
#define NEWMAX_INV_II static_cast<real_prec>(1.0)
#define NEWMIN_INV_II static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_II_
#define EXPONENT_INVARIANT_II  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_II  static_cast<real_prec>(1.0)
#endif
#endif
// OR
//#define _USE_DELTA2_   // ð².  This uses the memory space of _USE_INVARIANT_TIDAL_FIELD_I_
//#define _WRITE_DELTA2_
// ************************************************************************************************************************************
/**
* @brief Invariant III of the tidal field
*/
#ifdef _USE_IWEB_
#define _USE_INVARIANT_TIDAL_FIELD_III_
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_III_)
#define _USE_EXPONENT_INVARIANT_III_
#define _USE_SIGN_INVARIANT_III_
#define _MAP_TO_INTERVAL_INV_III_
#define NEWMAX_INV_III static_cast<real_prec>(1.0)
#define NEWMIN_INV_III static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_III_
#define EXPONENT_INVARIANT_III  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_III  static_cast<real_prec>(1.0)
#endif
#endif

// OR
//#define _USE_DELTA3_   //  ð³  This uses the memory space of _USE_INVARIANT_TIDAL_FIELD_II_
//#define _WRITE_DELTA3_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
/**
* @brief Invariant IV of the tidal field
*/

#ifdef _USE_IWEB_
//#define _USE_INVARIANT_TIDAL_FIELD_IV_
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
#define _USE_EXPONENT_INVARIANT_IV_
#define _USE_SIGN_INVARIANT_IV_
#define _MAP_TO_INTERVAL_INV_IV_
#define NEWMAX_INV_IV static_cast<real_prec>(1.0)
#define NEWMIN_INV_IV static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_IV_
#define EXPONENT_INVARIANT_IV  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_IV  static_cast<real_prec>(1.0)
#endif
#endif

// The following options can be also defined disjointly.
// ************************************************************************************************************************************
// OR

#ifdef _USE_TIDAL_ANISOTROPY_
#define _USE_EXPONENT_TIDAL_ANISOTROPY_
#define _USE_SIGN_USE_TIDAL_ANISOTROPY_
#define _MAP_TO_INTERVAL_TIDAL_ANISOTROPY_
#define NEWMAX_TIDAL_ANISOTROPY static_cast<real_prec>(1.0)
#define NEWMIN_TIDAL_ANISOTROPY static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_TIDAL_ANISOTROPY_
#define EXPONENT_TIDAL_ANISOTROPY  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_USE_TIDAL_ANISOTROPY  static_cast<real_prec>(1.0)
#endif
#endif
//  OR
// ************************************************************************************************************************************



#ifdef _USE_ELLIPTICITY_
#define _USE_EXPONENT_ELLIPTICITY_
#define _USE_SIGN_ELLIPTICITY_
#define _MAP_TO_INTERVAL_ELLIPTICITY_
#define NEWMAX_ELLIPTICITY static_cast<real_prec>(1.0)
#define NEWMIN_ELLIPTICITY static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_ELLIPTICITY_
#define EXPONENT_ELLIPTICITY  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_ELLIPTICITY  static_cast<real_prec>(1.0)
#endif
#endif

// ************************************************************************************************************************************


#ifdef _USE_PROLATNESS_
#define _USE_EXPONENT_PROLATNESS_
#define _USE_SIGN_PROLATNESS_
#define _MAP_TO_INTERVAL_PROLATNESS_
#define NEWMAX_PROLATNESS static_cast<real_prec>(1.0)
#define NEWMIN_PROLATNESS static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_PROLATNESS_
#define EXPONENT_PROLATNESS  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_PROLATNESS  static_cast<real_prec>(1.0)
#endif
#endif
//  OR
// ************************************************************************************************************************************
//#define _USE_S2_   // s²  This takes the memory space from _USE_TIDAL_ANISOTROPY_
#ifdef _USE_S2_
#define _WRITE_S2_
#define _USE_EXPONENT_S2_
#define _USE_SIGN_S2_
#define _MAP_TO_INTERVAL_S2_
#define NEWMAX_S2 static_cast<real_prec>(1.0)
#define NEWMIN_S2 static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_S2_
#define EXPONENT_S2  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_S2  static_cast<real_prec>(1.0)
#endif
#endif


#endif // end if mock mode


// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ******************************************** VELOCITIES ****************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************
/**
* @brief  Use information from the DM peculiar velocity field.
  @details  in order to characterize the bias In particular, use 2 eigenvalues of the shear of the V field from DM
*/
//#define _USE_VELOCITIES_
// -----------------------------------------------------------------------------

#ifdef _USE_VELOCITIES_
#ifndef _USE_LPT_
#define _READ_VELOCITIES_  // if we use vels, we can read them from files
#endif
#endif

// ************************************************************************************************************************************
// ******************************************** V--WEB ********************************************************************************
// ************************************************************************************************************************************

#if defined (_USE_CWEB_V_)  || defined (_USE_TWEB_V_)
#define _USE_CWC_V_  // Use the V_ classification
#endif

#ifdef _USE_TWEB_V_
#define _USE_VEL_KNOTS_V_  //Use the velocity dispersion of the cells classified as knots
#endif


#define VKMIN  0.0
#define VKMAX  1e10


// ************************************************************************************************************************************
//#ifndef _GET_BiasMT_REALIZATIONS_
//#define _TEST_THRESHOLDS_RESIDUALS_   // Define only when the thresholds of the T and V eigenvalues are to be explored sistematicaly in thw first (raw) interation
//#endif

// ************************************************************************************************************************************
// ******************************************** V-INVARIANTS **************************************************************************
// ************************************************************************************************************************************
#ifdef _USE_IVWEB_V_
// --------------------------------
#ifdef _USE_EXPONENT_INVARIANT_VS_I_
//#define _USE_EXPONENT_INVARIANT_VS_I_
//#define _USE_SIGN_INVARIANT_VS_I_
#define _MAP_TO_INTERVAL_INV_SHEAR_I_
#define NEWMAX_INV_SHEAR_I static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_I static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_I_
#define EXPONENT_INVARIANT_VS_I  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_I  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
#define _USE_EXPONENT_INVARIANT_VS_II_
//#define _USE_SIGN_INVARIANT_VS_II_
#define _MAP_TO_INTERVAL_INV_SHEAR_II_
#define NEWMAX_INV_SHEAR_II static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_II static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_II_
#define EXPONENT_INVARIANT_VS_II  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_II  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
#define _USE_EXPONENT_INVARIANT_VS_III_
//#define _USE_SIGN_INVARIANT_VS_III_
#define _MAP_TO_INTERVAL_INV_SHEAR_III_
#define NEWMAX_INV_SHEAR_III static_cast<real_prec>(1.0)
#define NEWMIN_INV_SHEAR_III static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_VS_III_
#define EXPONENT_INVARIANT_VS_III  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_VS_III  static_cast<real_prec>(1.0)
#endif
#endif

// --------------------------------
#endif

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// *********************************************P-WEB**********************************************************************************
// Use the three invariants of the tensor defined as the second derivative of delta.
// The first invariant is nabla²\delta. We speciofy it here despite of being uised below. Below er use it in the context of bias terms
// appearing in PT.

#ifdef _USE_PWEB_
// --------------------------------
#ifdef _USE_INVARIANT_PWEB_I_
#define _WRITE_INVARIANT_PWEB_I_
#define _USE_EXPONENT_INVARIANT_PWEB_I_
#define _USE_SIGN_INVARIANT_PWEB_I_
#define _MAP_TO_INTERVAL_INVARIANT_PWEB_I_
#define NEWMAX_INVARIANT_PWEB_I static_cast<real_prec>(1.0)
#define NEWMIN_INVARIANT_PWEB_I static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_PWEB_I_
#define EXPONENT_INVARIANT_PWEB_I  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_PWEB_I  static_cast<real_prec>(1.0)
#endif
#endif

// *****************************************************************

#ifdef _USE_INVARIANT_PWEB_II_
#define _WRITE_INVARIANT_PWEB_II_
#define _USE_EXPONENT_INVARIANT_PWEB_II_
#define _USE_SIGN_INVARIANT_PWEB_II_
#define _MAP_TO_INTERVAL_INVARIANT_PWEB_II_
#define NEWMAX_INVARIANT_PWEB_II static_cast<real_prec>(1.0)
#define NEWMIN_INVARIANT_PWEB_II static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_PWEB_II_
#define EXPONENT_INVARIANT_PWEB_II  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_PWEB_II  static_cast<real_prec>(1.0)
#endif
#endif
// *****************************************************************

#ifdef _USE_INVARIANT_PWEB_III_
#define _WRITE_INVARIANT_PWEB_III_
#define _USE_EXPONENT_INVARIANT_PWEB_III_
#define _USE_SIGN_INVARIANT_PWEB_III_
#define _MAP_TO_INTERVAL_INVARIANT_PWEB_III_
#define NEWMAX_INVARIANT_PWEB_III static_cast<real_prec>(1.0)
#define NEWMIN_INVARIANT_PWEB_III static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_INVARIANT_PWEB_III_
#define EXPONENT_INVARIANT_PWEB_III  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_INVARIANT_PWEB_III  static_cast<real_prec>(1.0)
#endif
#endif

// --------------------------------
#endif


// ******************************************** BIAS TERMS ****************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#ifndef _USE_PWEB_
// These trhee bias terms use the memmory space of the the PWEB
// ************************************************************************************************************************************
//#define _USE_NABLA2DELTA_  //  Nabla ²ð   This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_I
#ifdef _USE_NABLA2DELTA_
#define _WRITE_NABLA2DELTA_
#define _USE_EXPONENT_NABLA2DELTA_
//#define _USE_SIGN_NABLA2DELTA_
#define _MAP_TO_INTERVAL_NABLA2DELTA_
#define NEWMAX_NABLA2DELTA static_cast<real_prec>(1.0)
#define NEWMIN_NABLA2DELTA static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_NABLA2DELTA_
#define EXPONENT_NABLA2DELTA  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_NABLA2DELTA  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
//#define _USE_S2DELTA_  //  s²ð   This takes the memory space from _USE_INVARANT_SHEAR_VFIELD_II
#ifdef _USE_S2DELTA_
#define _WRITE_S2DELTA_
#define _USE_EXPONENT_S2DELTA_
#define _USE_SIGN_S2DELTA_
#define _MAP_TO_INTERVAL_S2DELTA_
#define NEWMAX_S2DELTA static_cast<real_prec>(1.0)
#define NEWMIN_S2DELTA static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_S2DELTA_
#define EXPONENT_S2DELTA  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_S2DELTA  static_cast<real_prec>(1.0)
#endif
#endif
// ************************************************************************************************************************************
//#define _USE_S3_        //   s³    This tkes the memory space from _USE_INVARANT_SHEAR_VFIELD_III
#ifdef _USE_S3_
#define _WRITE_S3_
#define _USE_EXPONENT_S3_
#define _USE_SIGN_S3_
#define _MAP_TO_INTERVAL_S3_
#define NEWMAX_S3 static_cast<real_prec>(1.0)
#define NEWMIN_S3 static_cast<real_prec>(-1.0)
#ifdef _USE_EXPONENT_S3_
#define EXPONENT_S3  static_cast<real_prec>(1./9.)
#else
#define EXPONENT_S3  static_cast<real_prec>(1.0)
#endif
#endif
// ***********************************************************************************************************************************
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ********************************************
// When BiasMT finishes, generate cats with particle posistions according to LPT
// ********************************************
// If catalog is to be produced, def if we want velocities or not
#ifdef _GET_BiasMT_CAT_
#ifndef _ASSIGN_PROPERTIES_TO_REFERENCE_
#define _GET_VELOCITIES_
#endif
#define CATWITHV
#endif

// ********************************************
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
#ifdef _USE_NEIGHBOURS_
#undef _USE_NEIGHBOURS_
#endif
#define N_C_BIN1 static_cast<ULONG>(4)
#define REDUCT_N_C_BIN1 static_cast<real_prec>(4.0)  // With this factor we reducefrom N-C_BIN1 to N-C_BIN1 / REDUCT_N_C_BIN1 the number of bins in the mass assignment procedure.
#define C1_MIN  static_cast<real_prec>(-1.)
#define C1_MAX  static_cast<real_prec>(1.)
#define DELTA_C1 (C1_MAX-C1_MIN)/(static_cast<double>(N_C_BIN1))
#else
#define N_C_BIN1 static_cast<ULONG>(1)
#define C1_MIN -200
#define C1_MAX  200
#define DELTA_C1 static_cast<real_prec>(1.)
#endif
// ********************************************
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_)
#define N_C_BIN2 static_cast<ULONG>(20)
#define REDUCT_N_C_BIN2 static_cast<real_prec>(4.0)
#define C2_MIN  static_cast<real_prec>(-1.)
#define C2_MAX  static_cast<real_prec>(1.)
#define DELTA_C2 (C2_MAX-C2_MIN)/(static_cast<double>(N_C_BIN2))
#else
#define N_C_BIN2 static_cast<ULONG>(1)
#define C2_MIN -200
#define C2_MAX  200
#define DELTA_C2 static_cast<real_prec>(1.)
#endif
// ********************************************
// This slot is dedicaqtd for the tital anisotropy
#if defined (_USE_TIDAL_ANISOTROPY_)  || defined (_USE_S2_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
#define N_C_BIN3 static_cast<ULONG>(100)
#define C3_MIN  static_cast<real_prec>(-1.0)
#define C3_MAX  static_cast<real_prec>(1.0)
#define DELTA_C3 (C3_MAX-C3_MIN)/(static_cast<double>(N_C_BIN3))
#else
#define N_C_BIN3 static_cast<ULONG>(1)
#define C3_MIN -200
#define C3_MAX  200
#define DELTA_C3 static_cast<real_prec>(1.)
#endif

// ********************************************
// This slot is dedicaqtd for the invariant odf the vshear field OR the Nabla² delta term
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_) || defined (_USE_INVARIANT_PWEB_I_)
#define N_CV_BIN1 static_cast<ULONG>(50)
#define CV1_MIN  static_cast<real_prec>(-40)
#define CV1_MAX  static_cast<real_prec>(40)
#define DELTA_CV1 (CV1_MAX-CV1_MIN)/(static_cast<double>(N_CV_BIN1))
#else
#define N_CV_BIN1 static_cast<ULONG>(1)
#define CV1_MIN  static_cast<real_prec>(-40)
#define CV1_MAX  static_cast<real_prec>(40)
#define DELTA_CV1 static_cast<real_prec>(1.)

#endif
// ********************************************
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_) || defined (_USE_INVARIANT_PWEB_II_) || defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
#define N_CV_BIN2 static_cast<ULONG>(2)
#define CV2_MIN  static_cast<real_prec>(-1.)
#define CV2_MAX  static_cast<real_prec>(1.)
#define DELTA_CV2 (CV2_MAX-CV2_MIN)/(static_cast<double>(N_CV_BIN2))
#else
#define N_CV_BIN2 static_cast<ULONG>(1)
#define CV2_MIN  static_cast<real_prec>(-40)
#define CV2_MAX  static_cast<real_prec>(40)
#define DELTA_CV2 static_cast<real_prec>(1.)
#endif

// ********************************************
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_) || defined (_USE_INVARIANT_PWEB_III_)
#define N_CV_BIN3 static_cast<ULONG>(20)
#define CV3_MIN  static_cast<real_prec>(-40)
#define CV3_MAX  static_cast<real_prec>(40)
#define DELTA_CV3 (CV3_MAX-CV3_MIN)/(static_cast<double>(N_CV_BIN3))
#else
#define N_CV_BIN3 static_cast<ULONG>(1)
#define CV3_MIN  static_cast<real_prec>(-40)
#define CV3_MAX  static_cast<real_prec>(40)
#define DELTA_CV3 static_cast<real_prec>(1.)
#endif





// Transform the vels from km/smto Mpc/h  / sec
// If input velñs are already in Mpc/h/sec, undef
#define _VEL_KMS_TO_MPCHS_


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//IMPORTANT DEFINITIONS FOR LPT
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//ONLY GALAXY NUMBER COUNTS CALCULATION
#define GONLY
//ONLY DM FIELD CALCULATION

#define DMONLY
//#undef DMONLY

#define TET

///////////////////////////////////////////////7
//STRUCTURE FORMATION MODEL
////////////////////////////////////////////////

#undef MUSCLE
#undef MUSCLE2LPT
#undef PTCURL
#undef SMOOTHCURL

#undef LPT3A
#undef LPT3B


///////////////////////////////////////////////
//RSDs
////////////////////////////////////////////////

#define _USE_FOGS_
#undef NORSD


///#define _MOVE_DM_TO_REDSHIFT_SPACE_

///////////////////////////////////////////////
//BASIC DEFINITIONS
////////////////////////////////////////////////

#define SAMPRAN

// ********************************************
#define SAVEMEM


// ********************************************
//#define PARTRAN // distribute halos/galaxies inside cells using the position of the available DM particles



// ********************************************
#ifdef _USE_LPT_
#define SAVECOMP // safe calculations for RSD at the expense of writing more fields
#undef CELLBOUND
#endif
// ********************************************
#undef GFINDIFF// gradient with finite differences
#define  GFFT    // gradient with FFTs
// ********************************************
//#define  _USE_GFINDIFF_EIGENV_    // gradient with FFTs FOR eIGENVALUES: esto da dips y bias
#define  _USE_GFFT_EIGENV_    // gradient with FFTs  FOR EIGENVALUES: Mejor que usar gradin
// ********************************************
//#define _USE_ZERO_PADDING_POT_
#define _EXTRA_NFT_FACTOR_   static_cast<int>(2)
// ********************************************
#define LOGRO
///////////////////////////////////////////////
///////////////////////////////////////////////
//BASIC DEFINITIONS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FOURIER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define FOURIER_DEF_1  // only tested with DEF_2
#define FOURIER_DEF_2  // only tested with DEF_2
#define FFTW_OPTION FFTW_ESTIMATE
#define FORWARD FFTW_FORWARD
#define BACKWARD FFTW_BACKWARD
#define fourier_space true
#define real_space false
#define to_Fspace true
#define to_Rspace false
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define _CONVERT_MPC_per_H_TO_KM_per_SECOND_
#define cgs_Mpc num_1
#define cgs_sec num_1
#define cgs_km static_cast<real_prec>(0.3240779290e-19 * cgs_Mpc) // 1 kilometer in 1 Mpc
#define cgs_clight static_cast<real_prec>(0.9715611892e-14 * cgs_Mpc/cgs_sec)
#define eps static_cast<real_prec>(1.e-14)
#define epsint static_cast<real_prec>(1.0e-6) /* numerical accuracy for integrations */
#define NEVAL 1000    /* numerical integration a2com */
#define _MESS_ So.message_screen("ABA", __LINE__);
#define _ACHTUNG_ So.message_warning("ACHTUNG in line ", __LINE__)
#define GAUSS_LEGANDERE_NODES_INTEGRATION static_cast<ULONG>(500)
#define MAXIMUM_MULTIPOLE_POWER static_cast<int>(1)
#define MAXIMUM_MULTIPOLE_WINDOW_MATRIX static_cast<int>(5)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define if Y is in NGP and want to be converted to CIC in order to smooth contours.
//By default
//#define  _KONV_
#ifdef BIAS_MODE
// Use this to convert NGP to CIC when using number counts in the tracers and want to display the bias
//#define  _NGP2CIC_Y_
//Uset this to convolve the DM with a kernel, which is read from the output directory
//#define  _KONV_
//define if a mask containig redshifts is to be used for bias analysis involving reconstrunction
//#define _USE_REDSHIFT_MASK_
//#define _USE_BINARY_MASK_  // if the file binary_mask is 1 or 0, define this. If not, do not use it
#ifdef _USE_REDSHIFT_MASK_
//#define REDSHIFT_MIN static_cast<real_prec>(0.35)
//#define REDSHIFT_MAX static_cast<real_prec>(2.34)
//#define N_REDSHIFT_BINS static_cast<int>(9)
#define REDSHIFT_MIN static_cast<real_prec>(0)
#define REDSHIFT_MAX static_cast<real_prec>(5)
#define N_REDSHIFT_BINS static_cast<int>(1)
#define DELTA_Z static_cast<real_prec>((REDSHIFT_MAX-REDSHIFT_MIN)/static_cast<real_prec>(N_REDSHIFT_BINS))
#endif
#ifndef _USE_REDSHIFT_MASK_
#define N_REDSHIFT_BINS 1
#endif
#else
#define N_REDSHIFT_BINS 1
#endif
#define __use_new_loops_bias_
//#define _use_random_kernel_  //derpecated
//#define _SHOW_EMPTY_CELLS_
#define DENSITY "density"
#define OVERDENSITY "overdensity"
#define _PERIODIC_BC_MAS_
// The values are used to get cosmological functions to interpolate upon when using the power spectrum code,
#define Z_MAX static_cast<real_prec>(1.5)
#define Z_MIN static_cast<real_prec>(0)
// if this is defined a test chaning the Invariantes for input foriles from the hydro project is done
//#define _HYDROTEST_
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Define to measure power spectrum using large random set split in a number of ascii files
*/
//#define _USE_SEVERAL_RANDOM_FILES_
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_TIWEB_
#if defined (_USE_CWEB_) || defined (_USE_TWEB_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_CWEB_
#if defined (_USE_TWEB_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_TWEB_
#if defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_IWEB_
#if defined (_USE_IKWEB_) || defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_IKWEB_
#if defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_PWEB_
#if defined (_USE_CWEB_V_) || defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_CWEB_V_
#if defined (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_TWEB_V_
#if defined (_USE_IWEB_V_) || defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif

#elif defined _USE_IWEB_V_
#if defined  (_USE_AWEB_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif


#elif defined _USE_AWEB_
#if defined (_USE_CWEB_) || defined (_USE_TWEB_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_PWEB_) || defined (_USE_CWEB_V_) || defined  (_USE_TWEB_V_) || defined  (_USE_IWEB_V_) 
#define WARNING_MODELS true
#else
#define WARNING_MODELS false
#endif


#endif


// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// DEFINE SOME COLORS
#define _USE_COLORS_

// ****************************************************************************************
/**
 * @brief Reset Color
*/
#define RESET   "\033[0m"
#ifdef _USE_COLORS_
// ********************************************
/**
 * @brief Color Black
*/
#define BLACK   "\033[30m"      /* Black */
#else
#define BLACK   RESET     /* Black */
#endif
// ********************************************



#ifdef _USE_COLORS_
/**
 * @brief Color Red
*/
#define RED     "\033[31m"      /* Red */
#else
#define RED   RESET
#endif
// ********************************************

#ifdef _USE_COLORS_
/**
 * @brief Color Green
*/
#define GREEN   "\033[32m"      /* Green */
#else
#define GREEN   RESET
#endif
// ********************************************
/**
 * @brief Color Yellow
*/
#ifdef _USE_COLORS_
#define YELLOW  "\033[33m"      /* Yellow */
#else
#define YELLOW   RESET
#endif
// ********************************************
/**
 * @brief Color Blue
*/
#ifdef _USE_COLORS_
#define BLUE    "\033[34m"      /* Blue */
#else
#define BLUE   RESET
#endif
// ********************************************
/**
 * @brief Color Magenta
*/
#ifdef _USE_COLORS_
#define MAGENTA "\033[35m"
#else
#define MAGENTA   RESET
#endif

// ********************************************
/**
 * @brief Color Cyan
*/
#ifdef _USE_COLORS_
#define CYAN    "\033[36m"      /* Cyan */
#else
#define CYAN   RESET
#endif
// ********************************************
/**
 * @brief Color White
*/
#ifdef _USE_COLORS_
#define WHITE   "\033[37m"      /* White */
#else
#define WHITE   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color Boldblack
*/
#ifdef _USE_COLORS_
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#else
#define BOLDBLACK   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldRed
*/
#ifdef _USE_COLORS_
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#else
#define BOLDRED   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldGreen
*/
#ifdef _USE_COLORS_
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#else
#define BOLDGREEN   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldYellow
*/
#ifdef _USE_COLORS_
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#else
#define BOLDYELLOW   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldBlue
*/
#ifdef _USE_COLORS_
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#else
#define BOLDBLUE   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldMagenta
*/
#ifdef _USE_COLORS_
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#else
#define BOLDMAGENTA   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldCyan
*/
#ifdef _USE_COLORS_
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#else
#define BOLDCYAN   RESET     /* Black */
#endif
// ********************************************
/**
 * @brief Color BoldWhite
*/
#ifdef _USE_COLORS_
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
#else
#define BOLDWHITE   RESET     /* Black */
#endif
// ********************************************
// ********************************************
#define COLOR_DEFINED BOLDGREEN
#define COLOR_UNDEFINED RED
// ************************************************************************************************************************************
// ************************************************************************************************************************************

#ifdef BIAS_MODE
//#define _USE_CWC_
//#undef _USE_CWC_
#endif
// ************************************************************************************************************************************
/**
 * @brief Enable if the input Pk is in form of De
*/
//#define _USE_IC_INPUT_POWER_DELTA_  
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
//Chuncks in Neighbouring Cells
/**
 * @brief This is used when low  ram memmory is at hand, dividing the assessment of e.g. mach number in slices of the box
*/
#define _USE_CHUNCKS_NC_
// ************************************************************************************************************************************
/**
 * @brief even numbers preferd as long as Nft is even too. If this is to be Nft, better disable _USE_CHUNCKS_NC_
*/
#define chunck_factor static_cast<int>(4)
// ************************************************************************************************************************************
#define SN_TOLERANCE_HBIAS static_cast<real_prec>(1.0)
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#define _ONLY_CWT_AND_PRIMARY_HBIAS_
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
//#define USE_GALAXY_TOOLS
// ************************************************************************************************************************************
// ************************************************************************************************************************************

#ifdef USE_GALAXY_TOOLS
#ifdef SINGLE_PREC
#undef SINGLE_PREC
#undef real_prec
#undef complex_prec
#define DOUBLE_PREC
#define real_prec double
#define complex_prec fftw_complex
#endif
#endif
// ************************************************************************************************************************************
// ************************************************************************************************************************************





