// ================================================================================================
// This file describes the class Parameters, used to set all the parameters adopted in the analysis
// ================================================================================================

#ifndef __PARAMETERS__
#define __PARAMETERS__

#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <string>
using namespace std;



class Parameters 
{
  
 private :
  

  double k_min_integration;
  double k_max_integration;

  double om_matter;
  double om_cdm;
  double om_radiation;
  double om_baryons;
  double om_vac;
  double om_k;
  double f_baryon;
  double Hubble;
  double hubble;
  double w_eos;
  double N_eff;
  double sigma8;
  double A_s;
  double n_s;
  double alpha_s;
  double Tcmb;
  double RR;
  double M_reference;
  bool use_wiggles;
  bool fixed_redshift;
  
  double redshift;
  double redshift_min;
  double redshift_max;

  string mass_function_fit;
  
  double M_min_effective;
  double M_max_effective;
  double A_gas;
  double B_gas;
  double mstar;
  
  double Delta_SO;
  int hod_model;
  double muno_hod;
  double alpha_hod;
  double mmin_hod;
  double scatter_hod;

  double sigma_ln;
  double sigma_red;
  double missing_flux;

  string density_profile;


  double M_min_mf;
  double M_max_mf;
  string scale_mf;
  int n_points_mf;
  string mass_function_output_file;
  string halo_mass_bias_fit;
  string halo_mass_bias_output_file;
  string effective_halo_mass_bias_output_file;
  string effective_halo_mean_number_density_output_file;
  bool compute_output_linear_power_spectrum;
  bool compute_output_non_linear_power_spectrum;
  string scale_ps;
  double k_min_ps;
  double k_max_ps;
  double kstar;
  double GAL_BIAS;
  double Amc;
  int n_points_ps;
  string linear_matter_ps_output_file;
  string non_linear_matter_ps_halo_fit_output_file;
  string non_linear_matter_ps_pt_output_file;
  double coef_concentration_amp;
  double coef_concentration;


  string galaxy_power_spectrum_halo_model_output_file;
  string galaxy_correlation_function_halo_model_output_file;

  string scale_cf;
  double r_min_cf;
  double r_max_cf;
  int n_points_cf;
  string linear_matter_cf_output_file;
  string non_linear_matter_cf_halo_fit_output_file;

  bool compute_output_linear_correlation_function;
  bool compute_output_non_linear_correlation_function;
  bool compute_density_profile;
  double r_min_dp;
  double r_max_dp;
  string scale_dp_r;
  int n_points_dp_r;
  string density_profile_r_output_file;
  double k_min_dp;
  double k_max_dp;
  string scale_dp_k;
  int n_points_dp_k;
  string density_profile_k_output_file ;




 public:

  // Default constructor
  Parameters () {}
  // Constructor: the input is the parameter file
  Parameters (string &);
  // Destructor
  ~Parameters () {}

  // function to get private variables
  double _k_max_integration () {return k_max_integration;}
  double _k_min_integration () {return k_min_integration;}
  double _om_matter () {return om_matter;}
  double _om_cdm () {return om_cdm;}
  double _om_radiation () {return om_radiation;}
  double _om_baryons () {return om_baryons;}
  double _om_vac () {return om_vac;}
  double _om_k () {return om_k;}
  double _f_baryon () {return f_baryon;}
  double _Hubble () {return Hubble;}
  double _hubble () {return hubble;}
  double _n_s () {return n_s;}
  double _w_eos () {return w_eos;}
  double _N_eff () {return N_eff;}
  double _sigma8 () {return sigma8;}
  double _A_s () {return A_s;}
  double _alpha_s () {return alpha_s;}
  double _Tcmb() {return Tcmb;}
  double _RR() {return RR;}
  double _Delta_SO () {return Delta_SO;}
  bool _use_wiggles() {return use_wiggles;}
  bool _fixed_redshift() {return fixed_redshift;}
  double _redshift() {return redshift;}
  double _redshift_min() {return redshift_min;}
  double _A_gas() {return A_gas;}
  double _B_gas() {return B_gas;}
  double _mstar() {return mstar;}
  double _kstar() {return kstar;}
  double _GAL_BIAS() {return GAL_BIAS;}
  double _Amc() {return Amc;}

  double _coef_concentration_amp() {return coef_concentration_amp;}
  double _coef_concentration() {return coef_concentration;}


  string _mass_function_fit() {return mass_function_fit;}

  double _sigma_ln() {return sigma_ln;}
  double _sigma_red() {return sigma_red;}
  double _missing_flux() {return missing_flux;}
  string _density_profile(){return density_profile;}
  double _muno_hod() {return muno_hod;}
  double _mmin_hod() {return mmin_hod;}
  double _scatter_hod() {return scatter_hod;}
  int _hod_model() {return hod_model;}
  double _alpha_hod() {return alpha_hod;}


  string _galaxy_power_spectrum_halo_model_output_file(){return galaxy_power_spectrum_halo_model_output_file;}
  string _galaxy_correlation_function_halo_model_output_file(){return galaxy_correlation_function_halo_model_output_file;}
  double _M_min_effective() {return M_min_effective;}
  double _M_max_effective() {return M_max_effective;}
  double _M_min_mf() {return M_min_mf;}
  double _M_max_mf() {return M_max_mf;}
  string _scale_mf() {return scale_mf ;}
  int _n_points_mf() {return n_points_mf ;}
  string _mass_function_output_file() {return mass_function_output_file ;}
  string _halo_mass_bias_fit() {return halo_mass_bias_fit ;}
  string _halo_mass_bias_output_file() {return halo_mass_bias_output_file ;}
  string _effective_halo_mass_bias_output_file() {return effective_halo_mass_bias_output_file ;}

  string _effective_halo_mean_number_density_output_file() {return effective_halo_mean_number_density_output_file ;}
  bool _compute_output_linear_power_spectrum() {return compute_output_linear_power_spectrum ;}
  bool _compute_output_non_linear_power_spectrum() {return compute_output_non_linear_power_spectrum  ;}
  string _scale_ps() {return scale_ps ;}
  double _k_min_ps() {return k_min_ps ;}
  double _k_max_ps() {return k_max_ps ;}
  int _n_points_ps() {return n_points_ps ;}
  string _linear_matter_ps_output_file() {return linear_matter_ps_output_file  ;}
  string _non_linear_matter_ps_halo_fit_output_file() {return non_linear_matter_ps_halo_fit_output_file ;}
  string _non_linear_matter_ps_pt_output_file() {return non_linear_matter_ps_pt_output_file ;}

  string _scale_cf() {return scale_cf ;}
  double _r_min_cf() {return r_min_cf ;}
  double _r_max_cf() {return r_max_cf ;}
  int _n_points_cf() {return n_points_cf;}
  string _linear_matter_cf_output_file() {return linear_matter_cf_output_file ;}
  string _non_linear_matter_cf_halo_fit_output_file() {return non_linear_matter_cf_halo_fit_output_file ;}

  bool _compute_density_profile() {return compute_density_profile ;}
  double _r_min_dp() {return r_min_dp ;}
  double _r_max_dp() {return r_max_dp ;}
  string _scale_dp_r() {return scale_dp_r ;}
  int _n_points_dp_r() {return n_points_dp_r ;}
  string _density_profile_r_output_file() {return density_profile_r_output_file ;}
  double _k_min_dp() {return k_min_dp ;}
  double _k_max_dp() {return k_max_dp ;}
  string _scale_dp_k() {return scale_dp_k ;}
  int _n_points_dp_k() {return n_points_dp_k ;}
  string _density_profile_k_output_file(){return density_profile_k_output_file ;}
  bool _compute_output_linear_correlation_function(){return compute_output_linear_correlation_function;}
  bool _compute_output_non_linear_correlation_function(){return compute_output_non_linear_correlation_function;}


};

#endif
