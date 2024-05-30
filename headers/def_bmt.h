#ifndef _DEF_BMT_
#define _DEF_BMT_
/**
 * @file def_bmt.h
 * @brief Pre-processor directives for Bias Mapping Technique
 * @author Andres Balaguera-Antolínez
 * @version 1.0 2018
 * @date   2018-2024
 * */
///////////////////////////////////////////////////////////////////////////////////////////////////
#include "def.h"
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Parallelization for the reconstruction of properties
*/
#define _USE_OMP_re
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Parallelism with OMP.
 * @details Define to test parallel regions in the BiasMT, specially at the multi-level assignment
*/
#define _USE_OMP_L0_
#define _USE_OMP_LX_
// used to parallelize assignments in do{} loops in BiasMT.cpp.
// The current implementation is not what I waNt, for it allows the assignment to a bit more number of objects than what I have defined through the tolerance and the thresholds
// Nevertheless, given that in the case in which the DM are not that used in the calibration,the tolerance are <1, and that excess induced by the parallelization just brings
// the figures towarsd the original values, so it compensates.
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief To compute halo ref power spectrum in bins of prpperties
*/
#define _COMPUTE_REF_POWER_
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Get power from the cat constructed without any info of properties
*/
//#define _CAT_POWER_ONLY_COORDS_
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
/**
* @brief Test for Vmax (or Mvir) assignment mixing stat and rad_from_ref approaches
* @details This option uses the method  BiasMT::assign_tracer_mass_new_new()
*/
//#define _USE_HYBRID_ASSIGNMENT_NEW_
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_HYBRID_ASSIGNMENT_NEW_
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
/**
 * @brief Number of bins in the primary proprty to pply psudo-learning approach to reduce scatter
*/
#define pNbins static_cast<ULONG>(200)
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
/**
 * @brief Target scatter
*/
//#define SIGMA_REC static_cast<real_prec>(0.05)
#define SIGMA_REC static_cast<real_prec>(0.03)
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
/**
 * @brief Number of bins away from the diagonal (exclusion zone) aimimng at not swaping there
*/
#define Exclusion_bins static_cast<ULONG>(2)
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
/**
 * @brief  Define to Test to see if what one learns from _USE_HYBRID_ASSIGNMENT_NEW_ is useful
 * @details This option uses the method  BiasMT::assign_tracer_mass_new_new()
 * @details and the oputput from the algorithm enabled by _USE_HYBRID_ASSIGNMENT_NEW_
*/
//#define _CHECK_HYBRID_V2_
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
#endif // endif _USE_HYBRID_ASSIGNMENT_NEW_
#endif // endif _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _GET_BiasMT_CAT_
#ifndef _ASSIGN_PROPERTIES_TO_REFERENCE_
#define _USE_NP_VEL_
#define _CORRECT_VEL_HALOMASS_  // ADHOC correction to vels in tow intervals of halo mass
#define _ADD_RANDOM_VEL_
#define _CORRECT_MEAN_VELOCITIES_
#define _VEL_CORR_EXP_ static_cast<real_prec>(0.2)
#define _VEL_CORR_EXP_KNOTS_ static_cast<real_prec>(0.21)
#define _VEL_CORR_EXP_FILAMENTS_ static_cast<real_prec>(0.2)
#define _VEL_CORR_EXP_SHEETS_ static_cast<real_prec>(0.22)
#define _VEL_CORR_EXP_VOIDS_ static_cast<real_prec>(0.1)
#define _CORRECT_VEL_DENSITY_
#endif
#endif
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\'
/**
* @brief Number of bins in bias. Mins and max are computed from the reference catalog.
*/
#define N_BINS_BIAS static_cast<ULONG>(15)
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_MACH_NUMBER_
/**
* @brief Compute the local overdenisty
*/
#define N_BINS_LO static_cast<int>(10)
#define N_BINS_MACH static_cast<int>(10)
#else
#define N_BINS_LO static_cast<int>(1)
#define N_BINS_MACH static_cast<int>(1)
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////
/**
* @brief Define to get power in bins of properties from the mock and the reference
*/
//#define _GET_POWER_BMT_

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
#endif
