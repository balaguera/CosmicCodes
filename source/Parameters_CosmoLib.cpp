// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains the function used to read parameters from input file                 **
// Andres Balaguera Antolinez                                                              **
// abalant@gmail.com                                                                       **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
#ifdef SINGLE_PREC
#undef SINGLE_PREC
#define DOUBLE_PREC
#endif

# include "../headers/Parameters_CosmoLib.h"

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ParametersCosmolib::read_pars()
{
  // Read the parameters from the parameter file
  ifstream fin_parameters (this->parameters_file.c_str());
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<parameters_file<<"!"<<endl; exit(1); }
  // a line in parameter file has a form: par_name = par_value
  string line_in_file;
  string par_name;
  string equality;
  string par_value;
  // go trough lines in parameters file
  while (getline(fin_parameters,line_in_file)) 
  {
    // ignore lines starting with hashtag
    if (line_in_file[0] != '#' && line_in_file.empty()==0) {
      // read parameter name 
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;
      // check that second word is "="
      line_string >> equality;
      if (equality != "=") {
	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      // read parameter value 
      line_string >> par_value;
      if (par_value.empty()) {
        cout << "Value of " << par_name << " not specified in " << parameters_file << endl;
        cout << "Assuming a default value for " << par_name << endl;
        continue;
      }
      if (par_name == "k_min_integration") k_min_integration = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "k_max_integration") k_max_integration = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_matter") this->om_matter = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_cdm") this->om_cdm = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_radiation") this->om_radiation = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_baryons") this->om_baryons = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_vac") this->om_vac = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "om_k") this->om_k = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "Hubble") this->Hubble = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "hubble") this->hubble = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "n_s") this->n_s = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "w_eos") this->w_eos = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "f_baryon")this->f_baryon = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "N_eff") this->N_eff = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "sigma8") this->sigma8 = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "alpha_s") this->alpha_s = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "Tcmb") this->Tcmb = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "RR") this->RR = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "Delta_SO") this->Delta_SO = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "A_s") this->A_s = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "use_wiggles") this->use_wiggles = par_value.c_str();
      else if (par_name == "kstar") this->kstar =static_cast<real_prec>(atof( par_value.c_str()));
      else if (par_name == "Amc") this->Amc =static_cast<real_prec>(atof( par_value.c_str()));
      else if (par_name == "GAL_BIAS") GAL_BIAS = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "fixed_redshift"){
        if(par_value=="true")fixed_redshift=true;
        else if(par_value=="false")fixed_redshift=false;
      }
      else if (par_name == "Get_SO_from_BN"){
        if(par_value=="true")Get_SO_from_BN=true;
        else if(par_value=="false")Get_SO_from_BN=false;
      }
      else if (par_name == "file_power") file_power =par_value;
      else if (par_name == "redshift") redshift = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "redshift_min") redshift_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "redshift_max") redshift_max = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "nbins_redshift") nbins_redshift = static_cast<int>(atof(par_value.c_str()));
      else if (par_name == "Delta_SO") Delta_SO = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "A_gas") A_gas = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "B_gas") B_gas = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mstar") mstar = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "sigma_red") sigma_red = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "sigma_ln") sigma_ln = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "missing_flux") missing_flux = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mass_function_fit") mass_function_fit = par_value.c_str();
      else if (par_name == "halo_mass_bias_fit") halo_mass_bias_fit = par_value.c_str();
      else if (par_name == "density_profile")density_profile = par_value.c_str();
      else if (par_name == "hod_model")hod_model = atoi(par_value.c_str());
      else if (par_name == "muno_hod")muno_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "alpha_hod")alpha_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "alpha_A")alpha_A = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "alpha_B")alpha_B = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "alpha_C")alpha_C = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mmin_hod")mmin_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "scatter_hod")scatter_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mt_hod")mt_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "ms_hod")ms_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "Mstep_hod")Mstep_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "s_faint_hod")s_faint_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "s_bright_hod")s_bright_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "width_hod")width_hod = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "area_survey")area_survey = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "coef_concentration_amp")coef_concentration_amp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "coef_concentration")coef_concentration = static_cast<real_prec>(atof(par_value.c_str()));
      //the new ones
      else if (par_name == "galaxy_power_spectrum_halo_model_output_file")galaxy_power_spectrum_halo_model_output_file=par_value.c_str();
      else if (par_name == "galaxy_correlation_function_halo_model_output_file")galaxy_correlation_function_halo_model_output_file=par_value.c_str();
      else if (par_name == "M_min_effective")M_min_effective = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "M_max_effective")M_max_effective = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "M_min_mf")M_min_mf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "M_max_mf")M_max_mf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "scale_mf")scale_mf = par_value.c_str();
      else if (par_name == "n_points_mf")n_points_mf = atoi(par_value.c_str());
      else if (par_name == "mass_function_output_file")mass_function_output_file = par_value.c_str();
      else if (par_name == "halo_mass_bias_output_file")halo_mass_bias_output_file = par_value.c_str();
      else if (par_name == "effective_halo_mass_bias_output_file")effective_halo_mass_bias_output_file = par_value.c_str();
      else if (par_name == "effective_halo_mean_number_density_output_file") effective_halo_mean_number_density_output_file=par_value.c_str();
      else if (par_name == "compute_output_linear_correlation_function")
	   {
	     if(par_value=="true")compute_output_linear_correlation_function=true;
	     else if(par_value=="false")compute_output_linear_correlation_function=false;
	   }
      else if (par_name == "compute_output_non_linear_correlation_function"){
        if(par_value=="true")compute_output_non_linear_correlation_function=true;
        else if(par_value=="false")compute_output_non_linear_correlation_function=false;
      }
      else if (par_name == "compute_output_linear_power_spectrum"){
        if(par_value=="true")compute_output_linear_power_spectrum=true;
        else if(par_value=="false")compute_output_linear_power_spectrum=false;
      }
      else if (par_name == "use_file_power"){
        if(par_value=="true")use_file_power=true;
        else if(par_value=="false")use_file_power=false;
      }
      else if (par_name == "compute_output_non_linear_power_spectrum"){
        if(par_value=="true")compute_output_non_linear_power_spectrum=true;
        else if(par_value=="false")compute_output_non_linear_power_spectrum=false;
      }
      else if (par_name == "scale_ps")scale_ps = par_value.c_str();
      else if (par_name == "k_min_ps")k_min_ps = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "k_max_ps")k_max_ps = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "n_points_ps")n_points_ps = atoi(par_value.c_str());
      else if (par_name == "linear_matter_ps_output_file")linear_matter_ps_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_halo_fit_output_file")non_linear_matter_ps_halo_fit_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_pt_output_file")non_linear_matter_ps_pt_output_file = par_value.c_str();
      else if (par_name == "scale_cf")scale_cf = par_value.c_str();
      else if (par_name == "Output_directory")this->Output_directory = par_value.c_str();
      else if (par_name == "r_min_cf")r_min_cf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "r_max_cf")r_max_cf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "n_points_cf")n_points_cf = atoi(par_value.c_str());
      else if (par_name == "linear_matter_cf_output_file")linear_matter_cf_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_cf_halo_fit_output_file")non_linear_matter_cf_halo_fit_output_file = par_value.c_str();
      else if (par_name == "compute_density_profile"){
        if(par_value=="true")compute_density_profile=true;
        else if(par_value=="false")compute_density_profile=false;
      }
      else if (par_name == "r_min_dp")r_min_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "r_max_dp")r_max_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "scale_dp_r")scale_dp_r = par_value.c_str();
      else if (par_name == "n_points_dp_r")n_points_dp_r = atoi(par_value.c_str());
      else if (par_name == "density_profile_r_output_file") density_profile_r_output_file = par_value.c_str();
      else if (par_name == "scale_dp_k")scale_dp_k = par_value.c_str();
      else if (par_name == "k_min_dp")k_min_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "k_max_dp")k_max_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "n_points_dp_k")n_points_dp_k = atoi(par_value.c_str());
      else if (par_name == "density_profile_k_output_file") density_profile_k_output_file = par_value.c_str();
      else { cerr << "Unknown parameter " << par_name << endl; cerr << "Doing nothing" << endl;}
    }
  }
  fin_parameters.clear(); fin_parameters.close(); 
#ifdef _USE_COSMO_PARS_
  this->om_matter = COSMOPARS::Om_matter;
  this->om_radiation = COSMOPARS::Om_radiation;
  this->om_baryons = COSMOPARS::Om_baryons;
  this->om_cdm = COSMOPARS::Om_cdm;
  this->om_vac = COSMOPARS::Om_vac;
  this->om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->n_s = COSMOPARS::n_s;
  this->w_eos = COSMOPARS::w_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = true;
  this->RR = COSMOPARS::RR;
  this->alpha_s=COSMOPARS::alpha_s;
  this->Delta_SO=COSMOPARS::Delta_SO;
  this->f_baryon=COSMOPARS::f_baryon;
#endif

  // append directory names in front of filenames
  string dat = ".txt";
  string zr = "_redshift_"+ std::to_string(static_cast<float>(redshift));
  string ll = "_";
  string output = "../Output/";
  mass_function_output_file = output+mass_function_output_file+ll+mass_function_fit+zr+dat;
  halo_mass_bias_output_file = output+halo_mass_bias_output_file+ll+halo_mass_bias_fit+zr+dat;
  effective_halo_mass_bias_output_file = output+effective_halo_mass_bias_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  effective_halo_mean_number_density_output_file=output+effective_halo_mean_number_density_output_file+ll+mass_function_fit+zr+dat;
  galaxy_power_spectrum_halo_model_output_file=output+galaxy_power_spectrum_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  galaxy_correlation_function_halo_model_output_file=output+galaxy_correlation_function_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  linear_matter_cf_output_file= output+linear_matter_cf_output_file+zr+dat;
  non_linear_matter_cf_halo_fit_output_file=output+non_linear_matter_cf_halo_fit_output_file+zr+dat;
  linear_matter_ps_output_file=output+linear_matter_ps_output_file+zr+dat;
  non_linear_matter_ps_halo_fit_output_file=output+non_linear_matter_ps_halo_fit_output_file+zr+dat;
  non_linear_matter_ps_pt_output_file=output+non_linear_matter_ps_pt_output_file+zr+dat;
  density_profile_r_output_file=output+density_profile_r_output_file+ll+density_profile+zr+dat;
  density_profile_k_output_file=output+density_profile_k_output_file+ll+density_profile+zr+dat;
}

