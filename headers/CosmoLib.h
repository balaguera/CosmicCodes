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
#define _GET_ONLY_MASS_FUNCTION_
#define _GET_POWER_SPECTRUM_
#define _GET_NL_POWER_SPECTRUM_
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
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */

    Cosmology Cf;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    FileOutput File;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    Statistics Cs;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    PowerSpectrum Ps;
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */
    ScreenOutput So;

// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------------------------------------------
public:
    //////////////////////////////////////////////////////////
    /**
     *  @brief get the value of the the private member type_of_lbinning
     *  @return type_of_lbinning
     */

  CosmoLib(){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   *  @brief Constructor passed the impot parameter file to
   *  @return 
   */
  CosmoLib(Params _param) : params(_param) , s_cosmo_par(_param.s_cosmo_pars) {
    this->Cf.set_cosmo_pars(this->s_cosmo_par);
    this->Cf.check_cosmo_pars();
    Statistics Csa(this->params._M_min_mf(), this->params._M_max_mf(), this->params._npoints_mf());
    this->Cs=Csa;
    PowerSpectrum Psa(this->s_cosmo_par,this->params._kmin_integration() ,this->params._kmax_integration() ,this->params._npoints_dp_k(),10, this->params._M_min_mf(),params._M_max_mf(),params._npoints_mf());
    this->Ps=Psa;
  }
   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */
  ~CosmoLib(){}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return 
   */
   Params params;
   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */
   s_CosmologicalParameters s_cosmo_par;
   //////////////////////////////////////////////////////////
   /**
    *  @brief get the value of the the private member type_of_lbinning
    *  @return type_of_lbinning
    */
  vector<gsl_real> v_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_sigma_mass;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_mass_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_effective_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_effective_halo_mean_number_density;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_k_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_nl_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_nl_power_spectrum_pt;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real> v_l_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  vector<gsl_real>pk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real>kk_aux;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_r_cf;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_nl_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_l_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_r_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_density_profile_r;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_k_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_density_profile_k;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_matter_bias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<gsl_real> v_galaxy_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation_1h_ss;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation_1h_sc;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation_2h;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */
  vector<real_prec> v_galaxy_correlation;
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning
   */

  void show_cosmo_struct();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute cosmological obserrvables
   *  @return
   */
  void get_cosmolib();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute dNdz for galaxies
   *  @return
   */
  void get_dndz_gal();
  //////////////////////////////////////////////////////////
  /**
   *  @brief This does the same as the constructor above, but as an object of the type CosmoLib
   *  @return
   */
  void set_par_file(string);


};

#endif 





