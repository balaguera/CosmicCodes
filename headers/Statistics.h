////////////////////////////////////////////////////////////////////////////
/**
 * @class
 * @file Statistics.h
 * @brief Header file for the class Statistics::
 * @title Functions related to interpolations of mass function nand related halo statistics
 * @author Andres Balaguera-Antolínez
 * @version  1.0
 * @date     2017-2019
 */
////////////////////////////////////////////////////////////////////////////
#ifndef _STATISTICS_
#define _STATISTICS_
# include "NumericalMethods.h"
# include "Miscelanious.h"
# include "PowerSpectrumTH.h"
# include "Astrophysics.h"
# include "BiasFunctions.h"
# include "Hod.h"

////////////////////////////////////////////////////////////////////////////
#define GSL_INT_SIZE_k 200
#define GSL_INT_SIZE_M 200
////////////////////////////////////////////////////////////////////////////

class Statistics{
 private:
    //////////////////////////////////////////////////////////
    /**
     * @brief
    */
   Cosmology cosmo;
   //////////////////////////////////////////////////////////
   /**
   * @brief
  */
   s_CosmologicalParameters s_cosmo_pars;
    //////////////////////////////////////////////////////////
    /**
     * @brief
    */
  static gsl_real isigma(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real isigma_nu(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real isigma(gsl_real ,void *, vector<gsl_real>, vector<gsl_real>);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real iAs2sigma8(gsl_real,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real dsigma_dR(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_effective_halo_mass_bias(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_mass_function(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_mean_galaxy_number_density(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */  
  gsl_interp_accel *acc_mass_sigma;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */  
  gsl_interp_accel *acc_mass_function;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */  
   gsl_spline *spline_mass_sigma;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */  
   gsl_spline *spline_mass_function;
 
 
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
public:
  //////////////////////////////////////////////////////////
  /**
   * @brief Default constructor
  */
  Statistics(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  Statistics(real_prec m1, real_prec m2, int Np){
     compute_int_table_mass(m1,m2, Np);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  ~Statistics(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief Vairnce of mass fluctuations
   * @arg mass, log10(M),
   * @arg reds, Cosmological redshift,
  */
  real_prec sigma_masa(real_prec mass, real_prec reds);
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Vairnce of mass fluctuations
   * @arg mass, log10(M),
   * @arg reds, Cosmological redshift,
  */
  real_prec peak_height(real_prec mass, real_prec reds, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sigma_masa(real_prec ,real_prec,s_CosmologicalParameters *, vector<gsl_real>, vector<gsl_real>);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec As2sigma8(s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  /* real_prec sigma_l(real_prec, real_prec, void*); */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mass_function_D(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mass_function_D(real_prec, real_prec, vector<gsl_real>&, vector<gsl_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec bias(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec effective_halo_mass_bias(real_prec , real_prec , s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec effective_halo_mean_number_density(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  //  real_prec cluster_mass_function(real_prec, real_prec, s_CosmologicalParameters *);
  /* real_prec scale_dependent_bias(real_prec, real_prec, void*); */
  /* real_prec delta_fof(real_prec, real_prec, void*); */
  /* real_prec mass_function_light_cone(real_prec, real_prec, void*); */
  /* real_prec mass_function_prediction(real_prec, real_prec, void*); */
  /* real_prec lum_func_prediction(real_prec, real_prec, void*); */
  /* real_prec occupancy_variance_numerator(real_prec, real_prec, void*);   /\*See Smith et al 2011*\/ */
  /* real_prec lum_func_reflex(real_prec, real_prec, void*); */
  /* real_prec baryon_mass_function_reionization(real_prec, real_prec, real_prec, void *); */

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density(real_prec redshift,vector<gsl_real>&, vector<gsl_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density(real_prec redshift);
  real_prec mean_galaxy_number_density(real_prec redshift, real_prec min_mass);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  real_prec mean_halo_number_density(s_CosmologicalParameters *scp);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void non_linear_scales(real_prec &,real_prec& ,real_prec& , real_prec& );
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_mass(real_prec, real_prec, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_wavenumber(real_prec, real_prec, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_mass(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_int_table_wavenumber(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  gsl_integration_glfixed_table *wfd;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  gsl_integration_glfixed_table *wf;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_K;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_K;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
   bool use_external_power;
   //////////////////////////////////////////////////////////
   /**
    * @brief
   */
   vector<gsl_real> kvector_external;
   //////////////////////////////////////////////////////////
   /**
    * @brief
   */
   vector<gsl_real> power_external;
   //////////////////////////////////////////////////////////
   void set_cosmo_pars(s_CosmologicalParameters spars)
   {
       this->s_cosmo_pars=spars;
       this->cosmo.set_cosmo_pars(spars);
   }
  //////////////////////////////////////////////////////////
  void set_spline_vars_sigma()
  {
    acc_mass_sigma = gsl_interp_accel_alloc();
    spline_mass_sigma = gsl_spline_alloc (gsl_interp_linear, this->s_cosmo_pars.v_mass.size());
    gsl_spline_init (spline_mass_sigma, &this->s_cosmo_pars.v_mass[0], &this->s_cosmo_pars.v_sigma_mass[0], this->s_cosmo_pars.v_mass.size());
  }
//////////////////////////////////////////////////////////
  void set_spline_vars_mfunction()
  {
    acc_mass_function = gsl_interp_accel_alloc();
    spline_mass_function = gsl_spline_alloc (gsl_interp_linear, this->s_cosmo_pars.v_mass.size());
    gsl_spline_init (spline_mass_sigma, &this->s_cosmo_pars.v_mass[0], &this->s_cosmo_pars.v_mass_function[0], this->s_cosmo_pars.v_mass.size());
  }


};


#endif

