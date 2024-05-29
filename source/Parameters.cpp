// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains the function used to read parameters from input file                 **
// Andres Balaguera Antolinez                                                              **
// abalant@gmail.com                                                                       **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************


# include "../Include/Parameters.h"
# include "../Include/NumericalMethods.h"



Parameters::Parameters (string &parameters_file)
{

  

  // ***********************************************************************
  // HERE WE NEED TO INITIALIZE THE PARAMETERS


  // ***********************************************************************

  // *******************************************
  // read the parameters from the parameter file
  // *******************************************

  // open the parameters file
  ifstream fin_parameters (parameters_file.c_str());
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<parameters_file<<"!"<<endl; exit(1); }
  
  // a line in parameter file has a form: par_name = par_value
  string line_in_file;
  string par_name;
  string equality;
  string par_value;
  
  // go trough lines in parameters file
  while (getline(fin_parameters,line_in_file)) {
   
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
      
      // depending on parameter name assign parameter value to a corresponding variable

      // global and I/O parameters
      
      
      // cosmological parameters

      if (par_name == "k_min_integration") k_min_integration = (double)atof(par_value.c_str());
      else if (par_name == "k_max_integration") k_max_integration = (double)atof(par_value.c_str());
      else if (par_name == "om_matter") om_matter = (double)atof(par_value.c_str());
      else if (par_name == "om_cdm") om_cdm = (double)atof(par_value.c_str());
      else if (par_name == "om_radiation") om_radiation = (double)atof(par_value.c_str());
      else if (par_name == "om_baryons") om_baryons = (double)atof(par_value.c_str());
      else if (par_name == "om_vac") om_vac = (double)atof(par_value.c_str());
      else if (par_name == "om_k") om_k = (double)atof(par_value.c_str());
      else if (par_name == "Hubble") Hubble = (double)atof(par_value.c_str());
      else if (par_name == "hubble") hubble = (double)atof(par_value.c_str());
      else if (par_name == "n_s") n_s = (double)atof(par_value.c_str());
      else if (par_name == "w_eos") w_eos = (double)atof(par_value.c_str());
      else if (par_name == "f_baryon")f_baryon = (double)atof(par_value.c_str());
      else if (par_name == "N_eff") N_eff = (double)atof(par_value.c_str());
      else if (par_name == "sigma8") sigma8 = (double)atof(par_value.c_str());
      else if (par_name == "A_s") A_s = (double)atof(par_value.c_str());
      else if (par_name == "alpha_s") alpha_s = (double)atof(par_value.c_str());
      else if (par_name == "Tcmb") Tcmb = (double)atof(par_value.c_str());
      else if (par_name == "RR") RR = (double)atof(par_value.c_str());
      else if (par_name == "use_wiggles") use_wiggles = par_value.c_str();
      else if (par_name == "kstar") kstar =(double)atof( par_value.c_str());
      else if (par_name == "Amc") Amc =(double)atof( par_value.c_str());

      else if (par_name == "GAL_BIAS") GAL_BIAS = (double)atof(par_value.c_str());


      else if (par_name == "fixed_redshift"){
        if(par_value=="true")fixed_redshift=true;
        else if(par_value=="false")fixed_redshift=false;
      }


      else if (par_name == "redshift") redshift = (double)atof(par_value.c_str());
      else if (par_name == "redshift_min") redshift_min = (double)atof(par_value.c_str());
      else if (par_name == "redshift_max") redshift_max = (double)atof(par_value.c_str());
      else if (par_name == "Delta_SO") Delta_SO = (double)atof(par_value.c_str());
      else if (par_name == "A_gas") A_gas = (double)atof(par_value.c_str());
      else if (par_name == "B_gas") B_gas = (double)atof(par_value.c_str());
      else if (par_name == "mstar") mstar = (double)atof(par_value.c_str());
      else if (par_name == "sigma_red") sigma_red = (double)atof(par_value.c_str());
      else if (par_name == "sigma_ln") sigma_ln = (double)atof(par_value.c_str());
      else if (par_name == "missing_flux") missing_flux = (double)atof(par_value.c_str());
      else if (par_name == "mass_function_fit") mass_function_fit = par_value.c_str();
      else if (par_name == "halo_mass_bias_fit") halo_mass_bias_fit = par_value.c_str();
      else if (par_name == "density_profile")density_profile = par_value.c_str();
      else if (par_name == "hod_model")hod_model = atoi(par_value.c_str());
      else if (par_name == "muno_hod")muno_hod = (double)atof(par_value.c_str());
      else if (par_name == "alpha_hod")alpha_hod = (double)atof(par_value.c_str());
      else if (par_name == "mmin_hod")mmin_hod = (double)atof(par_value.c_str());
      else if (par_name == "scatter_hod")scatter_hod = (double)atof(par_value.c_str());

      else if (par_name == "coef_concentration_amp")coef_concentration_amp = (double)atof(par_value.c_str());
      else if (par_name == "coef_concentration")coef_concentration = (double)atof(par_value.c_str());
      //the new ones
      else if (par_name == "galaxy_power_spectrum_halo_model_output_file")galaxy_power_spectrum_halo_model_output_file=par_value.c_str();
      else if (par_name == "galaxy_correlation_function_halo_model_output_file")galaxy_correlation_function_halo_model_output_file=par_value.c_str();
      else if (par_name == "M_min_effective")M_min_effective = (double)atof(par_value.c_str());
      else if (par_name == "M_max_effective")M_max_effective = (double)atof(par_value.c_str());
      else if (par_name == "M_min_mf")M_min_mf = (double)atof(par_value.c_str());
      else if (par_name == "M_max_mf")M_max_mf = (double)atof(par_value.c_str());
      else if (par_name == "scale_mf")scale_mf = par_value.c_str();
      else if (par_name == "n_points_mf")n_points_mf = atoi(par_value.c_str());
      else if (par_name == "mass_function_output_file")mass_function_output_file = par_value.c_str();
      else if (par_name == "halo_mass_bias_output_file")halo_mass_bias_output_file = par_value.c_str();
      else if (par_name == "effective_halo_mass_bias_output_file")effective_halo_mass_bias_output_file = par_value.c_str();
      else if (par_name == "effective_halo_mean_number_density_output_file") effective_halo_mean_number_density_output_file=par_value.c_str();

      else if (par_name == "compute_output_linear_correlation_function"){
        if(par_value=="true")compute_output_linear_correlation_function=true;
        else if(par_value=="false")compute_output_linear_correlation_function=false;
      }

      else if (par_name == "compute_output_non_linear_correlation_function"){
        if(par_value=="true")compute_output_non_linear_correlation_function=true;
        else if(par_value=="false")compute_output_non_linear_correlation_function=false;
      }

      else if (par_name == "compute_output_non_linear_power_spectrum"){
        if(par_value=="true")compute_output_non_linear_power_spectrum=true;
        else if(par_value=="false")compute_output_non_linear_power_spectrum=false;
      }


      else if (par_name == "scale_ps")scale_ps = par_value.c_str();
      else if (par_name == "k_min_ps")k_min_ps = (double)atof(par_value.c_str());
      else if (par_name == "k_max_ps")k_max_ps = (double)atof(par_value.c_str());
      else if (par_name == "n_points_ps")n_points_ps = atoi(par_value.c_str());
      else if (par_name == "linear_matter_ps_output_file")linear_matter_ps_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_halo_fit_output_file")non_linear_matter_ps_halo_fit_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_pt_output_file")non_linear_matter_ps_pt_output_file = par_value.c_str();
      else if (par_name == "scale_cf")scale_cf = par_value.c_str();
      else if (par_name == "r_min_cf")r_min_cf = (double)atof(par_value.c_str());
      else if (par_name == "r_max_cf")r_max_cf = (double)atof(par_value.c_str());
      else if (par_name == "n_points_cf")n_points_cf = atoi(par_value.c_str());
      else if (par_name == "linear_matter_cf_output_file")linear_matter_cf_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_cf_halo_fit_output_file")non_linear_matter_cf_halo_fit_output_file = par_value.c_str();
      

      else if (par_name == "compute_density_profile"){
        if(par_value=="true")compute_density_profile=true;
        else if(par_value=="false")compute_density_profile=false;
      }
      else if (par_name == "r_min_dp")r_min_dp = (double)atof(par_value.c_str());
      else if (par_name == "r_max_dp")r_max_dp = (double)atof(par_value.c_str());
      else if (par_name == "scale_dp_r")scale_dp_r = par_value.c_str();
      else if (par_name == "n_points_dp_r")n_points_dp_r = atoi(par_value.c_str());
      else if (par_name == "density_profile_r_output_file") density_profile_r_output_file = par_value.c_str();
      
      else if (par_name == "scale_dp_k")scale_dp_k = par_value.c_str();
      else if (par_name == "k_min_dp")k_min_dp = (double)atof(par_value.c_str());
      else if (par_name == "k_max_dp")k_max_dp = (double)atof(par_value.c_str());
      else if (par_name == "n_points_dp_k")n_points_dp_k = atoi(par_value.c_str());
      else if (par_name == "density_profile_k_output_file") density_profile_k_output_file = par_value.c_str();
      
      
      else { cerr << "Unknown parameter " << par_name << endl; cerr << "Doing nothing" << endl;}

    }
  }
  fin_parameters.clear(); fin_parameters.close(); 

  // append directory names in front of filenames
  string dat = ".dat";
  string zr = "_redshift_"+ std::to_string((float)redshift);
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

