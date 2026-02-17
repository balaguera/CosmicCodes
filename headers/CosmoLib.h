// ************************************************************************************************************************************************************************** //
// ************************************************************************************************************************************************************************** //
/** @file CosmoLib.cpp
 *
 * @brief This file contains headers of the class CosmoLib
 * @details The class Cosmolib generates estaimates of a number of cosmlogical observables
 * @author Andres Balaguera Antolinez
 * @date 2007-2024
 */
// ************************************************************************************************************************************************************************** //
#ifndef _COSMOLIB_
#define _COSMOLIB_
// ************************************************************************************************************************************************************************** //
// Pre-proc directives for COSMO-LIB
//#define TIME
//#define _GET_EFF_BIAS_REDSHIFT_
//#define _GET_ONLY_MASS_FUNCTION_
//#define _GET_POWER_SPECTRUM_
//#define _GET_NL_POWER_SPECTRUM_
//#define _GET_NL_PT_POWER_SPECTRUM_

#ifdef _GET_NL_POWER_SPECTRUM_
//#define _GET_HM_POWER_SPECTRUM_
//#define _GET_HM_CORRELATION_FUNCTION_
#endif
// ************************************************************************************************************************************************************************** //
// ************************************************************************************************************************************************************************** //
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include "NumericalMethods.h"
# include "Params.h"
# include "CosmologicalFunctions.h"
# include "FileOutput.h"
# include "ScreenOutput.h"
# include "PowerSpectrumTH.h"
# include "CorrelationFunction.h"
# include "Astrophysics.h"
# include "ScalingRelations.h"
# include "Statistics.h"
# include "BiasFunctions.h"
# include "DensityProfiles.h"
# include "Hod.h"
//# include "GalaxyOperations.h"
//# include "AngularPowerSpectrumTH.h"
# include "McmcFunctions.h"
# include "Galaxy.h"

using namespace std;
using namespace Constants;

// ************************************************************************************************************************************************************************** //
// ************************************************************************************************************************************************************************** //

class CosmoLib{
 private:
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type Cosmology
     *  @return 
     */
    Cosmology Cf;
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type FileManager
     */
    FileOutput File;
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type Statistics
     */
    Statistics Cs;
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type PowerSpectrum
     */
    PowerSpectrum Ps;
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type CorrelationFunction
     */
    CorrelationFunctionTH SCf;
    //////////////////////////////////////////////////////////
    /**
     *  @brief Object of type ScreenOutput
     */
    ScreenOutput So;

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
public:
    //////////////////////////////////////////////////////////
    /**
     *  @brief Default constructor of the CosmoLib class
     */
  CosmoLib(){} 
  //////////////////////////////////////////////////////////
  /**
   *  @brief Constructor passing an object of type Params 
   *  @return 
   */
  CosmoLib(Params _param) : 
    params(_param) ,
    s_cosmo_par(_param.s_cosmo_pars),  
    Cs(_param._M_min_effective(), _param._M_max_effective(), _param._npoints_mf()), 
    Ps(_param.s_cosmo_pars,_param._kmin_integration() ,_param._kmax_integration() , _param._npoints_dp_k(),10, _param._M_min_mf(),_param._M_max_mf(), _param._npoints_mf()),
    SCf(_param)
   {
    this->Cf.set_cosmo_pars(this->s_cosmo_par);
    this->Cf.check_cosmo_pars();

    // Compute tables to perform integrals using GL weights
    this->Cs.compute_int_table_wavenumber(this->params._kmin_ps(),this->params._kmax_ps());

    // Compute tables to perform integrals using GL weights
    this->Cs.compute_int_table_mass(this->params._M_min_effective(), this->params._M_max_effective());
  }
   //////////////////////////////////////////////////////////
   /**
    *  @brief Default destructor
   */
  ~CosmoLib(){}
  //////////////////////////////////////////////////////////
  /**
   *  @brief Object of type Parms.
   */
   Params params;
   //////////////////////////////////////////////////////////
   /**
    *  @brief Strcucure of type s_CosmologicalParameters
    */
   s_CosmologicalParameters s_cosmo_par;
   //////////////////////////////////////////////////////////
   /**
    *  @brief Container for masses
    */
  vector<gsl_real> v_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for mass dispersion
   */
  vector<gsl_real> v_sigma_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for halo mass function
   */
  vector<gsl_real> v_mass_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for halo mass 
   */
  vector<gsl_real> v_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for effective halo mass 
   */
  vector<gsl_real> v_effective_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for effective halo mena number density
   */
  vector<gsl_real> v_effective_halo_mean_number_density;

  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for power spectrum
   */
  vector<gsl_real> v_k_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @brief   Container for non linear power spectrum 
   */
  vector<gsl_real> v_nl_power_spectrum_halofit;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for non linear power spectrum derived from perturbation theory
   */
  vector<gsl_real> v_nl_power_spectrum_pt;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for non linear power spectrum derived from EH fitting formulae
   */

  vector<gsl_real> v_l_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Auxiliary container
   */
  vector<gsl_real>pk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Auxiliary container
   */
  vector<gsl_real>kk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for separations in CF analysis
   */
  vector<real_prec> v_r_cf;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for non linear correlation function derived from halo fit
   */
  vector<real_prec> v_nl_correlation_function_halofit;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for non linear correlation function derived from perturbation theory

   */
  vector<real_prec> v_nl_correlation_function_pt;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for linear correlation function 

   */
  vector<real_prec> v_l_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  

   */
  vector<real_prec> v_r_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for density profile in real space

   */
  vector<real_prec> v_density_profile_r;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  
   */
  vector<gsl_real> v_k_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for density profile in Fourier space
   */
  vector<gsl_real> v_density_profile_k;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for the satellite-satellite power spectrum
   */
  vector<gsl_real> v_galaxy_power_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for the satellite-central power spectrum
   */
  vector<gsl_real> v_galaxy_power_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for the 2-h term power spectrum
   */
  vector<gsl_real> v_galaxy_power_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for galaxy-matter bias
   */
  vector<gsl_real> v_galaxy_matter_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for the total galaxy power spectrum
   */
  vector<gsl_real> v_galaxy_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for the satellite-satellite correlation function
   */
  vector<real_prec> v_galaxy_correlation_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for the satellite-central correlation function 
   */
  vector<real_prec> v_galaxy_correlation_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Container for the 2h term correlation function
   */
  vector<real_prec> v_galaxy_correlation_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief  Container for the total galaxy correlation function
   */
  vector<real_prec> v_galaxy_correlation;
  //////////////////////////////////////////////////////////
  /**
   *  @brief   Shows in the shell the cosmological parameters
   */
  void show_cosmo_struct();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute the elements of galaxy power spectrum in the framwork of the Halo Model.
   *  @details Observables are computed at a fixed redshift (given in the parameter file) and according to the HOD model selected in the parameter file.
   *  @return Halo mass function. 
   *  \image html halo_mass_function.png 
   *  @return Halo mass bias. 
   *  \image html halo_mass_bias.png
   *  @return File with the different components of halo model galaxy power spectrum, 
   *  \image html galaxy_power_halo_model.png
   */
  void get_hmodel();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute cosmological observables
   *  @details Observables are computed at a fixed redshift (given in the parameter file)
   */
  void get_cosmological_information();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute dNdz for galaxies
   *  @warning TO b finished.
   *  @return
   */
  void get_dndz_gal();

};

#endif 





