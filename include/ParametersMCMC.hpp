// ================================================================================================
// This file describes the class Parameters, used to set all the parameters adopted in the analysis
// ================================================================================================

#ifndef __PARAMETERS_MCMC__
#define __PARAMETERS_MCMC__

# include <fstream>
# include <iostream>
# include <vector>
# include <algorithm>
# include <math.h>
# include <sstream>
# include <iostream>
# include <iomanip>
# include <string>
using namespace std;



class ParametersMCMC 
{
  
 private :
  

  double k_min_integration;
  double k_max_integration;

  double om_matter;
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
  int n_points_ps;
  string linear_matter_ps_output_file;
  string non_linear_matter_ps_halo_fit_output_file;


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
  ParametersMCMC () {};
  // Constructor: the input is the parameter file
  ParametersMCMC(string &);
  // Destructor
  ~ParametersMCMC () {};

  int _n_points_dp_k() {return n_points_dp_k ;};
  string _density_profile_k_output_file(){return density_profile_k_output_file ;};
  bool _compute_output_linear_correlation_function(){return compute_output_linear_correlation_function;};
  bool _compute_output_non_linear_correlation_function(){return compute_output_non_linear_correlation_function;};


};

#endif
