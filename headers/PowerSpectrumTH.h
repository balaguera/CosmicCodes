//////////////////////////////////////////////////////////
/**
 * @file PowerSpectrumTH.h
 * @brief This file contains methods of the class PowerSpectrum
 * @details CLASS TO COMPUTE LINEAR AND NON LINEAR MATTER POWER SPECTRUM BASED ON THE EISENSTEIN AND HU FITTING FORMULAE AND THE HALO-FIT FROM SMITH ET AL.
 * @author Andres Balaguera Antolinez
 * @date 2007-2023
 */
//////////////////////////////////////////////////////////

#ifndef __POWER_SPECTRUMTH__
#define __POWER_SPECTRUMTH__

# include <iostream>
# include <math.h>
# include <cmath>
# include <cctype>
# include <stdio.h>
# include <fstream>

# include "NumericalMethods.h"
# include "CosmologicalFunctions.h"
# include "BiasFunctions.h"
# include "Hod.h"
# include "DensityProfiles.h"
# include "ScreenOutput.h"

using namespace std;

class PowerSpectrum{
 private:
 //////////////////////////////////////////////////////////
  /**
  * @brief
 */
  ScreenOutput So;
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
  static gsl_real dsigma_dR(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_nw(gsl_real ,void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit2(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit4(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  static gsl_real fun_aux_halo_fit_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit2_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real fun_aux_halo_fit4_dw(gsl_real , void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real Power_Spectrum_i(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real iP1loop(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real iFkernel(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  static gsl_real i_Galaxy_power_spectrum_h1_ss(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_Galaxy_power_spectrum_h1_sc(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_Galaxy_matter_bias(gsl_real, void *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  static gsl_real i_mean_galaxy_number_density(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec kmin_int;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec kmax_int;

 public:

  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  PowerSpectrum(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  PowerSpectrum(s_CosmologicalParameters _s_cosmo):s_cosmo_pars(_s_cosmo){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */

  PowerSpectrum(s_CosmologicalParameters _s_cosmo,real_prec _kmin, real_prec _kmax, int nk, int nm):s_cosmo_pars(_s_cosmo){
    this->kmin_int=_kmin;
    this->kmax_int=_kmax;
    compute_int_table_k_mu(_kmin,_kmax, nk, nm);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  PowerSpectrum(real_prec M1, real_prec M2, int Np){
     compute_int_table_mass(M1,M2, Np);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  PowerSpectrum(s_CosmologicalParameters scos, real_prec _kmin, real_prec _kmax, int nk, int nmu,real_prec M1, real_prec M2,real_prec Nm){
    this->kmin_int=_kmin;
    this->kmax_int=_kmax;
    this->s_cosmo_pars=scos;
    compute_int_table_k_mu(_kmin,_kmax, nk, nmu);
    compute_int_table_mass(M1,M2, Nm);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  ~PowerSpectrum(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec omhh; /* Omega_matter*h^2 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec obhh;   /* Omega_baryon*h^2 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_equality; /* Scale of equality, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec R_equality; /* Photon-baryon ratio at equality epoch */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sound_horizon;  /* Sound horizon at drag epoch, in Mpc */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_silk;   /* Silk damping scale, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_c;  /* CDM suppression */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_c;   /* CDM log shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_b;  /* Baryon suppression */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_b;   /* Baryon envelope shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec beta_node;  /* Sound horizon shift */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec k_peak;   /* Fit to wavenumber of first peak, in Mpc^-1 */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sound_horizon_fit;  /* Fit to sound horizon, in Mpc */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha_gamma;  /* Gamma suppression in approximate TF */
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_z(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_interpolated(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Q_Model_Matter_Power_Spectrum(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Primordial_Matter_Power_Spectrum(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_NW(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_DW(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_G_DW(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Linear_Matter_Power_Spectrum_G_NW(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(real_prec, real_prec, real_prec, real_prec, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
   real_prec Non_Linear_Matter_Power_Spectrum_PT(real_prec k);
   //////////////////////////////////////////////////////////
   /**
    * @brief
   */
  real_prec P1loop(real_prec k);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec sigma_masa(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec bias(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mass_function_D(real_prec, real_prec, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void mass_function_M_Z(vector<real_prec>, s_CosmologicalParameters *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<vector<real_prec> > MASS_FUNCTION_M_Z;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<vector<real_prec> > MASS_BIAS_M_Z;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec mean_galaxy_number_density();
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Non_Linear_Matter_Power_Spectrum_Halo_Fit_DW(real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_power_spectrum_h1_ss(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_power_spectrum_h1_sc(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec Galaxy_matter_bias(real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief Setting variables for Fitting Formula Eisenstein and Hu
   * @details  Set all the scalars quantities for Eisenstein & Hu 1997 fitting formula
  /* @param omega0hh The density of CDM and baryons, in units of critical dens, multiplied by the square of the Hubble constant, in units of 100 km/s/Mpc
  /* @param f_baryon The fraction of baryons to CDM
  /* @param Tcmb The temperature of the CMB in Kelvin.  Tcmb<=0 forces use of the COBE value of  2.728 K.
  /* @warning Units are always Mpc, never h^-1 Mpc.
  */
  void TFset_parameters(real_prec omega0hh , real_prec f_baryon, real_prec Tcmb);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFfit_onek(real_prec , real_prec &, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void TFfit_hmpc(real_prec, real_prec, real_prec , real_prec,
                  int, real_prec *, real_prec *, real_prec *, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFsound_horizon_fit(real_prec, real_prec , real_prec );
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFk_peak(real_prec , real_prec , real_prec );
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFnowiggles(real_prec, real_prec, real_prec,
                     real_prec , real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec TFzerobaryon(real_prec, real_prec , real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief Halo fit based on the revision by Takahashi
   * @details https://arxiv.org/pdf/1208.2701
  */
  void halo_fit(real_prec,  real_prec*, real_prec*, real_prec*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void halo_fit_integrals(real_prec*, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void halo_fit_z(real_prec, real_prec, real_prec,real_prec, real_prec, real_prec, real_prec *,real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void normalization(real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief NOrmalziation of Spectrum
  */
  real_prec normalization();
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void nl_scales_halo_fit(real_prec *, real_prec*, vector<real_prec>&, vector<real_prec>&, bool);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void kstar_integral(real_prec*);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void hf_aux(real_prec, real_prec,real_prec, real_prec, real_prec *, real_prec *, real_prec *, real_prec *,real_prec *,real_prec *,real_prec *, real_prec *);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @details We will use this function to compute the HaloFit power spectrum or any other function
     that makes use of the linear matter power spectrum. Therefore, we compute Plin first
     and then we just interpolated
  */
  void compute_int_table_k_mu(real_prec, real_prec, int,  int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  * @details These vectors are allocated after the computation of the linear power spectrum,
   and are meant to speed up the calculation of the non-linear power spectrum using halo fit.
  */
  void compute_int_table_mass(real_prec, real_prec, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_Pk_linear;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_kk;

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
  vector<gsl_real> v_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_sigma_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_mass_function;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> v_halo_mass_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_Mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> XX_mu;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  vector<gsl_real> WW_mu;
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
  /**
   * @brief
  */
  bool use_external_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void set_cosmo_pars(s_CosmologicalParameters spars){this->s_cosmo_pars=spars;}
};

#endif  
  



