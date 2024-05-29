// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains the function used to read parameters from input file                 **
// Andres Balaguera Antolinez                                                              **
// abalant@gmail.com                                                                       **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
# include "../Headers/ParametersMCMC.h"


ParametersMCMC::ParametersMCMC(string &parameters_file)
{

  

  // // *******************************************
  // // read the parameters from the parameter file
  // // *******************************************

  // // open the parameters file
  // ifstream fin_parameters (parameters_file.c_str());
  // if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<parameters_file<<"!"<<endl; exit(1); }
  
  // // a line in parameter file has a form: par_name = par_value
  // string line_in_file;
  // string par_name;
  // string equality;
  // string par_value;
  
  // // go trough lines in parameters file
  // while (getline(fin_parameters,line_in_file)) {
   
  //   // ignore lines starting with hashtag
  //   if (line_in_file[0] != '#' && line_in_file.empty()==0) {

  //     // read parameter name 
  //     stringstream line_string (line_in_file);
  //     line_string << line_in_file;
  //     line_string >> par_name;

  //     // check that second word is "="
  //     line_string >> equality;

  //     if (equality != "=") {
  // 	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
  // 	cerr << "Using a default value for " << par_name << endl; exit(1);
  //     }

  //     // read parameter value 
  //     line_string >> par_value;
      
  //     if (par_value.empty()) {
  // 	cout << "Value of " << par_name << " not specified in " << parameters_file << endl;
  // 	cout << "Assuming a default value for " << par_name << endl;
  // 	continue;
  //     }
      
  //     // depending on parameter name assign parameter value to a corresponding variable

  //     // global and I/O parameters
      
      
  //     // cosmological parameters

  //   }
  // }
  // fin_parameters.clear(); fin_parameters.close(); 
  // cout<<n_points_cf<<endl;
  
  // // append directory names in front of filenames
  // string dat = ".dat";
  // string zr = "_redshift_";// stdto_string(redshift);
  // string ll = "_";
  // string output = "../Output/";
  
  // mass_function_output_file = output+mass_function_output_file+ll+mass_function_fit+zr+dat;

  // halo_mass_bias_output_file = output+halo_mass_bias_output_file+ll+halo_mass_bias_fit+zr+dat;

  // effective_halo_mass_bias_output_file = output+effective_halo_mass_bias_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;

  // effective_halo_mean_number_density_output_file=output+effective_halo_mean_number_density_output_file+ll+mass_function_fit+zr+dat;


  // galaxy_power_spectrum_halo_model_output_file=output+galaxy_power_spectrum_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  // galaxy_correlation_function_halo_model_output_file=output+galaxy_correlation_function_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;

  // linear_matter_cf_output_file= output+linear_matter_cf_output_file+zr+dat;
  // non_linear_matter_cf_halo_fit_output_file=output+non_linear_matter_cf_halo_fit_output_file+zr+dat;

  // linear_matter_ps_output_file=output+linear_matter_ps_output_file+zr+dat;
  // non_linear_matter_ps_halo_fit_output_file=output+non_linear_matter_ps_halo_fit_output_file+zr+dat;

  // density_profile_r_output_file=output+density_profile_r_output_file+ll+density_profile+zr+dat;
  // density_profile_k_output_file=output+density_profile_k_output_file+ll+density_profile+zr+dat;


}

