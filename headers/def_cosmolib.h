/**
 * @file def.h
 * @brief Pre-processor directives for cosmicatlass
 * @author Andres Balaguera-Antolínez
 * @version   1.0
 * @date      2008-2023
 * */


// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// PREPROCESSOR DIRECTIVES FOR COSMO-LIB **************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************
// ****************************************************************************************

// ************************COSMOLOGICAL PARAMETERS*****************************************
// See the file cosmological_parameters.h for the namespaces with sets of cosmological parameters

/**
 * @brief If defined, the c ode  uses  the cosmological parametrs from  namespaces.-
 * @brief If undef, BAM uses cosmo para from input parameter file.
*/

#define _USE_COSMO_PARS_
/**
 * @brief Namespaces for cosmological parameters
*/

//#define _USE_UNITSIM_COSMOLOGY_
//#define _USE_PLANCK_COSMOLOGY_
#define _USE_SLICS_COSMOLOGY_
//#define _USE_MINERVA_COSMOLOGY_

#ifdef _USE_SLICS_COSMOLOGY_
#define COSMOPARS Cosmo_parameters_SLICS
#elif defined (_USE_PLANCK_COMSOLOGY_) || defined (_USE_UNITSIM_COSMOLOGY_)
#define COSMOPARS Cosmo_parameters_PLANCK
#elif defined (_USE_MINERVA_COMSOLOGY_)
#define COSMOPARS Cosmo_parameters_Minerva
#endif

// ****************************************************************************************
// ****************************************************************************************
// DEFINE SOME COLORS
// ****************************************************************************************
/**
 * @brief Reset Color
*/
#define RESET   "\033[0m"
// ********************************************
/**
 * @brief Color Black
*/
#define BLACK   "\033[30m"      /* Black */
// ********************************************
/**
 * @brief Color Red
*/
// ********************************************
#define RED     "\033[31m"      /* Red */
/**
 * @brief Color Green
*/
// ********************************************
#define GREEN   "\033[32m"      /* Green */
/**
 * @brief Color Yellow
*/
// ********************************************
#define YELLOW  "\033[33m"      /* Yellow */
/**
 * @brief Color Blue
*/
// ********************************************
#define BLUE    "\033[34m"      /* Blue */
/**
 * @brief Color Magenta
*/
// ********************************************
#define MAGENTA "\033[35m"      /* Magenta */
/**
 * @brief Color Cyan
*/
// ********************************************
#define CYAN    "\033[36m"      /* Cyan */
/**
 * @brief Color White
*/
// ********************************************
#define WHITE   "\033[37m"      /* White */
/**
 * @brief Color Boldblack
*/
// ********************************************
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
/**
 * @brief Color BoldRed
*/
// ********************************************
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
/**
 * @brief Color BoldGreen
*/
// ********************************************
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
/**
 * @brief Color BoldYellow
*/
// ********************************************
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
/**
 * @brief Color BoldBlue
*/
// ********************************************
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
/**
 * @brief Color BoldMagenta
*/
// ********************************************
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
/**
 * @brief Color BoldCyan
*/
// ********************************************
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
/**
 * @brief Color BoldWhite
*/
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ***************************************************PRECISION************************************************************************
// ************************************************************************************************************************************

//
// /**
// * @brief Define double precision for cosmicatlas
// * @details Applies all over the code except for gsl-type defined variables/containers
// */
#define DOUBLE_PREC

#ifndef DOUBLE_PREC
/**
 * @brief Define single precision for cosmicatlas
 * @details Applies all over the code except for gsl-type defined variables/containers
 * @details Defined when DOUBLE_PREC is undefined
*/
#define SINGLE_PREC
#endif


#define gsl_real double  // Even with single precision, gsl will use double precision for all its calculations
// ********************************************
#ifdef SINGLE_PREC
/**
 * @brief Precision at output
*/
#define _PREC_OUTPUT_ 4
// ********************************************
/**
 * @brief Precision for FFT operations
*/
#define fftwf_real float
// ********************************************

#ifndef real_prec
/**
 * @brief Effective Precision for cosmicatlas
*/
#define real_prec fftwf_real
#endif

/**
 * @brief Effective Precision for complex containers used in FFTW
*/
#define complex_prec fftwf_complex
#endif
#ifdef DOUBLE_PREC
#define _PREC_OUTPUT_ 6
#define fftw_real double
#define real_prec fftw_real
#define complex_prec fftw_complex
#endif
#define BIG_NUMBER 1e7
#define ULONG unsigned long
#define LONG long
#define UULONG unsigned long long
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
#define KNOT_MAX_SIZE static_cast<ULONG>(100000000)
#define V_MAX_SIZE 100000000
#define LARGE_NUMBER 1e20
#define MINUS_LARGE_NUMBER -100000000
#define NOCELL -2
#define NOVEL -999

#define num_1   static_cast<double>(1.)
#define num_2   static_cast<double>(2.)
#define num_0_1 static_cast<double>(0.1)
#define num_0_5 static_cast<double>(0.5)

// ************************************************************************************************************************************
// ************************************************************************************************************************************
// ************************************************************************************************************************************
#define _USE_OMP_



