#include "../Headers/ParametersF.h"

ParametersF::ParametersF (string &parameters_file)
{

  // ***************************
  // initialize I/O variables 
  // ***************************

  input_type= "catalog";
  dir_input = "../Input/";
  dir_output = "../Output/";
  file_catalogue = "galaxy_coordinates.dat";
  file_random = "random_coordinates.dat";
  delta_grid_file = "delta_grid_file";

  measure_cross=false;
  measure_cross_from_1=1;
  measure_cross_from_2=1;

  sys_of_coord_g = 0;
  i_coord1_g = 0;  
  i_coord2_g = 1;  
  i_coord3_g = 2;
  angles_units_g = "D";    
  use_random_catalog = true;
  use_random_catalog_cl = false;
  sys_of_coord_r = 0; 
  i_coord1_r = 0;
  i_coord2_r = 1;
  i_coord3_r = 2;
  i_weight1_r = 10;
  i_weight2_r = 10;
  i_weight3_r = 10;
  i_weight4_r = 10;
  i_weight1_g = 10;
  i_weight2_g = 10;
  i_weight3_g = 10;
  i_weight4_g = 10;

  i_property1_g = 20;
  i_property2_g = 20;
  i_property3_g = 20;
  i_property4_g = 20;
  i_property5_g = 20;
  i_property6_g = 20;

  i_property1_r = 20;
  i_property2_r = 20;
  i_property3_r = 20;
  i_property4_r = 20;
  i_property5_r = 20;
  i_property6_r = 20;

  use_weight1_r= false;
  use_weight2_r= false;
  use_weight3_r= false;
  use_weight4_r= false;

  use_weight1_g = false;
  use_weight2_g = false;
  use_weight3_g = false;
  use_weight4_g = false;
  new_los = false;
  new_Lbox = false;

  n_catalogues = 1;

  // ******************************************************
  // initialize variables for covariance of 2point statistics
  // ******************************************************

  use_one_data_for_covariance = false;
  use_one_random_for_covariance = true;


  // ******************************************************
  // initialize variables for power spectrum 
  // ******************************************************

  Nft = 256;           
  Lbox = 1000.;         
  mass_assignment_scheme = "CIC"; 
  MAS_correction = true;
  type_of_binning = "linear";
  k_bin_step = 0.5;
  N_log_bins = 10;           
  ndel_data = 1;     
  ndel_window = 1;   
  N_mu_bins = 100;          
  FKP_weight = false; 
  Pest = 100000;        
  FKP_error_bars = false;
  FKP_error_bars_exact = false;
  SN_correction = true; 
  nbar_tabulated = false;
  constant_depth = false;   
  N_z_bins = 100;
  redshift_min_sample = 0.0;
  redshift_max_sample = 2.0;
  N_dndz_bins = 100;
  new_N_dndz_bins = 60;
  Healpix_resolution = 3;

  // **********************************
  // initialize parameter for Yamamoto Direct sum
  // **********************************
  kmax_y_ds=0.2;


  // **********************************
  // initialize parameter for Bispectrum
  // **********************************
  kmax_bk=0.2;
  kmin_bk=0.0;
  use_fundamental_mode_as_kmin_bk=true;

  // **********************************
  // initialize parameters for angular power spectrum
  // **********************************
  cL_max=200;
  cL_min=2;
  type_of_input_file= "cat";
  number_of_redshift_bins = 2;
  number_of_Lbins = 20;


  // **********************************
  // initialize cosmological parameters
  // **********************************

  om_matter = 0.237;
  om_radiation = 0.0;
  om_baryons = 0.041;
  om_vac = 0.763;
  om_k = 0.0;
  Hubble = 100.0;
  hubble = 0.7;
  spectral_index = 0.9;
  w_eos = -1.0;
  N_eff = 3.03;
  sigma8 = 0.773;
  Tcmb = 2.73;


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
      if (par_name == "Statistics") statistics = par_value;
      else if (par_name == "Input_directory") dir_input = par_value;
      else if (par_name == "Output_directory") dir_output = par_value;
      else if (par_name == "Catalogue_file") file_catalogue = par_value;
      else if (par_name == "Input_type") input_type = par_value;
      else if (par_name == "delta_grid_file") delta_grid_file = par_value;
      else if (par_name == "delta_grid_file2") delta_grid_file2 = par_value;
      else if (par_name == "delta_grid_file3") delta_grid_file3 = par_value;
      else if (par_name == "delta_grid_file4") delta_grid_file4 = par_value;
      else if (par_name == "ngal_delta") ngal_delta = atof(par_value.c_str());
      else if (par_name == "measure_cross") {
        if(par_value == "true")this->measure_cross = true;
        else if(par_value == "false")this->measure_cross = false;
      }

      else if (par_name == "measure_cross_from_1")measure_cross_from_1 = atoi(par_value.c_str());
      else if (par_name == "measure_cross_from_2")measure_cross_from_2 = atoi(par_value.c_str());


      else if (par_name == "sys_of_coord_g") sys_of_coord_g = atoi(par_value.c_str());
      else if (par_name == "i_coord1_g") i_coord1_g = atoi(par_value.c_str()); 
      else if (par_name == "i_coord2_g") i_coord2_g = atoi(par_value.c_str());
      else if (par_name == "i_coord3_g") i_coord3_g = atoi(par_value.c_str());
      else if (par_name == "i_mean_density_g") i_mean_density_g = atoi(par_value.c_str());
      else if (par_name == "i_weight1_g") i_weight1_g = atoi(par_value.c_str());
      else if (par_name == "i_weight2_g") i_weight2_g = atoi(par_value.c_str());
      else if (par_name == "i_weight3_g") i_weight3_g = atoi(par_value.c_str());
      else if (par_name == "i_weight4_g") i_weight4_g = atoi(par_value.c_str());

      else if (par_name == "use_weight1_g") {
        if(par_value == "true")this->use_weight1_g = true;
        else if(par_value == "false")this->use_weight1_g = false;
      }

      else if (par_name == "use_weight2_g") {
	if(par_value == "true")use_weight2_g = true;
	else if(par_value == "false")use_weight2_g = false;
      }
      else if (par_name == "use_weight3_g") {
	if(par_value == "true")use_weight3_g = true;
	else if(par_value == "false")use_weight3_g = false;
      }
      else if (par_name == "use_weight4_g") {
	if(par_value == "true")use_weight4_g = true;
	else if(par_value == "false")use_weight4_g = false;
      }


      
      else if (par_name == "i_property1_g") i_property1_r = atoi(par_value.c_str());
      else if (par_name == "i_property2_g") i_property2_r = atoi(par_value.c_str());
      else if (par_name == "i_property3_g") i_property3_r = atoi(par_value.c_str());
      else if (par_name == "i_property4_g") i_property4_r = atoi(par_value.c_str());
      else if (par_name == "i_property5_g") i_property5_r = atoi(par_value.c_str());
      else if (par_name == "i_property6_g") i_property6_r = atoi(par_value.c_str());

      else if (par_name == "new_Lbox"){
	if (par_value == "yes"){new_Lbox = 1;}
	else if (par_value == "no"){new_Lbox = 0;}
      }
      
      else if (par_name == "n_catalogues") n_catalogues = atoi(par_value.c_str());
      else if (par_name == "use_one_data_for_covariance") use_one_data_for_covariance = (par_value == "true" || par_value == "True") ? true : false;
      else if (par_name == "use_one_random_for_covariance") use_one_random_for_covariance = (par_value == "true" || par_value == "True") ? true : false;

      else if (par_name == "angles_units_g") {
	if (par_value=="D") angles_units_g = "D";
	else if (par_value=="R") angles_units_g = "R";
	else {cout<<"Unknown angle units"<<endl; exit(1) ;}
      }
      else if (par_name == "use_random_catalog"){
	if (par_value=="false") use_random_catalog = false;
	else if (par_value=="true") use_random_catalog = true;
      }

      else if (par_name == "use_random_catalog_cl"){
	if (par_value=="no") use_random_catalog_cl = false;
	else if (par_value=="yes") use_random_catalog_cl = true;
      }

      else if (par_name == "sys_of_coord_r") sys_of_coord_r = atoi(par_value.c_str());
      else if (par_name == "i_coord1_r") i_coord1_r = atoi(par_value.c_str());
      else if (par_name == "i_coord2_r") i_coord2_r = atoi(par_value.c_str());
      else if (par_name == "i_coord3_r") i_coord3_r = atoi(par_value.c_str());
      else if (par_name == "i_mean_density_r") i_mean_density_r = atoi(par_value.c_str());
      else if (par_name == "i_weight1_r") i_weight1_r = atoi(par_value.c_str());
      else if (par_name == "i_weight2_r") i_weight2_r = atoi(par_value.c_str());
      else if (par_name == "i_weight3_r") i_weight3_r = atoi(par_value.c_str());
      else if (par_name == "i_weight4_r") i_weight4_r = atoi(par_value.c_str());

      else if (par_name == "i_property1_r") i_property1_r = atoi(par_value.c_str());
      else if (par_name == "i_property2_r") i_property2_r = atoi(par_value.c_str());
      else if (par_name == "i_property3_r") i_property3_r = atoi(par_value.c_str());
      else if (par_name == "i_property4_r") i_property4_r = atoi(par_value.c_str());
      else if (par_name == "i_property5_r") i_property5_r = atoi(par_value.c_str());
      else if (par_name == "i_property6_r") i_property6_r = atoi(par_value.c_str());


     else if (par_name == "i_mask_pixel") i_mask_pixel= atoi(par_value.c_str());
     else if (par_name == "i_mask_alpha") i_mask_alpha= atoi(par_value.c_str());
     else if (par_name == "i_mask_delta") i_mask_delta= atoi(par_value.c_str());
     else if (par_name == "i_mask_flag") i_mask_flag= atoi(par_value.c_str());

      else if (par_name == "use_weight1_r") {
	if(par_value == "true")use_weight1_r = true;
	else if (par_value == "false")use_weight1_r = false;
      }
      
      else if (par_name == "use_weight2_r") {
	if(par_value == "true")use_weight2_r = true;
	else if (par_value == "false")use_weight2_r = false;
      }
      else if (par_name == "use_weight3_r") {
	if(par_value == "true")use_weight3_r = true;
	else if (par_value == "false")use_weight3_r = false;
      }
      else if (par_name == "use_weight4_r") {
	if(par_value == "true")use_weight4_r = true;
	else if (par_value == "false")use_weight4_r = false;
      }

      else if (par_name == "angles_units_r"){
	if (par_value == "D") angles_units_r = "D";
	else if (par_value == "R") angles_units_r = "R";
	else { cout<<"Unknown angle units"<<endl; exit(1); }
      }
      else if (par_name == "Random_file") file_random = par_value;
      else if (par_name == "Name_survey") Name_survey = par_value;


      // parameters for the two-point correlation function

      // parameters for the power spectrum 
      else if (par_name == "new_los") {
	if (par_value == "false") new_los = false;
	else if (par_value == "true") new_los = true;
      }
      else if (par_name == "Nft") Nft = atoi(par_value.c_str());
      else if (par_name == "Lbox") Lbox = (double)atof(par_value.c_str());
      else if (par_name == "mass_assignment_scheme")
	{
	  if (par_value=="NGP") mass_assignment_scheme = "NGP";
	  else if (par_value=="CIC") mass_assignment_scheme = "CIC";
	  else if (par_value=="TSC") mass_assignment_scheme = "TSC";
	  else{cout<<"Unknown mass assigment scheme"<<endl;exit(1);}
	}
      else if (par_name == "MAS_correction"){
	if (par_value=="true") MAS_correction = true;
	else if (par_value=="false") MAS_correction = false;
      }
      else if (par_name == "type_of_binning"){
	if (par_value=="linear") type_of_binning = "linear";
	else if (par_value=="log") type_of_binning = "log";
	else{cout<<"Unknown type of binning"<<endl; exit(1);}
      }
      else if (par_name == "k_bin_step") k_bin_step = atof(par_value.c_str());
      else if (par_name == "N_log_bins") N_log_bins = atoi(par_value.c_str());
      else if (par_name == "ndel_data") ndel_data = atoi(par_value.c_str());
      else if (par_name == "ndel_window") ndel_window = atoi(par_value.c_str());
      else if (par_name == "N_mu_bins") N_mu_bins = atoi(par_value.c_str());
      else if (par_name == "FKP_weight"){
	if (par_value=="false") FKP_weight = false;
	else if (par_value=="true") FKP_weight = true;
      }
      else if (par_name == "Pest") Pest = atoi(par_value.c_str());
      else if (par_name == "SN_correction"){
	if (par_value=="false") SN_correction = false;
	else if (par_value=="true") SN_correction = true;
      }
      else if (par_name == "FKP_error_bars"){
	if (par_value=="false") FKP_error_bars = false;
	else if (par_value=="true") FKP_error_bars = true;
      }
      else if (par_name == "FKP_error_bars_exact"){
	if (par_value=="false") FKP_error_bars_exact = false;
	else if (par_value=="true") FKP_error_bars_exact = true;
      }
      else if (par_name == "nbar_tabulated"){
	if (par_value=="false") nbar_tabulated = false;
	else if (par_value=="true") nbar_tabulated = true;
      }
      else if (par_name == "constant_depth"){
	if (par_value=="no") constant_depth = false;
	else if (par_value=="yes") constant_depth = true;
      }
      else if (par_name == "N_z_bins") N_z_bins = atoi(par_value.c_str());
      else if (par_name == "redshift_min_sample") redshift_min_sample = (double)atof(par_value.c_str());
      else if (par_name == "redshift_max_sample") redshift_max_sample = (double)atof(par_value.c_str());
      else if (par_name == "N_dndz_bins") N_dndz_bins = atoi(par_value.c_str());
      else if (par_name == "new_N_dndz_bins") new_N_dndz_bins = atoi(par_value.c_str());
      else if (par_name == "area_survey") area_survey = (double)atof(par_value.c_str());
      else if (par_name == "Healpix_resolution") Healpix_resolution = atoi(par_value.c_str());
      else if (par_name == "file_power") file_power = par_value;  
      else if (par_name == "file_power_log") file_power_log = par_value;

      else if (par_name == "file_window") file_window = par_value;
      else if (par_name == "file_dndz") file_dndz = par_value;
      else if (par_name == "file_power2d") file_power2d = par_value;
      else if (par_name == "file_power2d_mk") file_power2d_mk = par_value;
      else if (par_name == "file_bispectrum") file_bispectrum = par_value;
      
      //output files for Yamamoto estimator of moments
      else if (par_name == "kmax_y_ds")kmax_y_ds = (double)atof(par_value.c_str());


      //Parameters for the bispectrum
      else if (par_name == "kmax_bk")kmax_bk = (double)atof(par_value.c_str());
      else if (par_name == "kmin_bk")kmin_bk = (double)atof(par_value.c_str());
      else if (par_name == "use_fundamental_mode_as_kmin_bk"){
	if (par_value=="no") use_fundamental_mode_as_kmin_bk = false;
	else if (par_value=="yes") use_fundamental_mode_as_kmin_bk = true;
      }      
      
      
      //Angular power spectrum
      else if (par_name == "cL_max")cL_max = atoi(par_value.c_str());
      else if (par_name == "cL_min")cL_min = atoi(par_value.c_str());
      else if (par_name == "number_of_redshift_bins")number_of_redshift_bins = atoi(par_value.c_str());

      else if (par_name == "type_of_zbins")type_of_zbins = par_value.c_str();


      else if (par_name == "zbin1_min")zbin1_min = (double)atof(par_value.c_str());
      else if (par_name == "zbin1_max")zbin1_max = (double)atof(par_value.c_str());

      else if (par_name == "zbin2_min")zbin2_min = (double)atof(par_value.c_str());
      else if (par_name == "zbin2_max")zbin2_max = (double)atof(par_value.c_str());

      else if (par_name == "zbin3_min")zbin3_min = (double)atof(par_value.c_str());
      else if (par_name == "zbin3_max")zbin3_max = (double)atof(par_value.c_str());



      else if (par_name == "number_of_Lbins") number_of_Lbins = atoi(par_value.c_str());
      
      else if (par_name == "type_of_input_file") type_of_input_file = par_value;  
      else if (par_name == "input_file_mask") input_file_mask = par_value;
      else if (par_name == "type_of_lbinning") type_of_lbinning = par_value;

      else if (par_name == "file_power_cl") file_power_cl = par_value;
      else if (par_name == "file_cl_window")file_cl_window = par_value;
      //FB power spectrum

      else if (par_name == "Boundary_conditions")Boundary_conditions = par_value.c_str();
      else if (par_name == "compute_new_kln")compute_new_kln = par_value.c_str();
      else if (par_name == "file_kln")file_kln = par_value.c_str();
      else if (par_name == "file_modes")file_modes = par_value.c_str();
      else if (par_name == "l_max_fb")l_max_fb = atoi(par_value.c_str());
      else if (par_name == "l_min_fb")l_min_fb = atoi(par_value.c_str());
      else if (par_name == "redshift_max_fb")redshift_max_fb = (double)atof(par_value.c_str());
      else if (par_name == "redshift_min_fb")redshift_min_fb = (double)atof(par_value.c_str());
      else if (par_name == "k_max_fb")k_max_fb = (double)atof(par_value.c_str());
      else if (par_name == "file_power_fb")file_power_fb = par_value.c_str();


      // cosmological parameters
      else if (par_name == "om_matter") om_matter = (double)atof(par_value.c_str());
      else if (par_name == "om_radiation") om_radiation = (double)atof(par_value.c_str());
      else if (par_name == "om_baryons") om_baryons = (double)atof(par_value.c_str());
      else if (par_name == "om_vac") om_vac = (double)atof(par_value.c_str());
      else if (par_name == "om_k") om_k = (double)atof(par_value.c_str());
      else if (par_name == "Hubble") Hubble = (double)atof(par_value.c_str());
      else if (par_name == "hubble") hubble = (double)atof(par_value.c_str());
      else if (par_name == "spectral_index") spectral_index = (double)atof(par_value.c_str());
      else if (par_name == "w_eos") w_eos = (double)atof(par_value.c_str());
      else if (par_name == "N_eff") N_eff = (double)atof(par_value.c_str());
      else if (par_name == "sigma8") sigma8 = (double)atof(par_value.c_str());
      else if (par_name == "Tcmb") Tcmb = (double)atof(par_value.c_str());

      else
      {
	//cerr << "Unknown parameter " << par_name << endl; cerr << "Doing nothing" << endl;
      }


    }
  }
  fin_parameters.clear(); fin_parameters.close(); 



  // Append directory names in front of filenames.
  // The files related to the Fourier space (and harmonic) space measurements
  // are appended suth that the final output file has the structure
  // STATISTICS_SAMPLENAME_SOME FEATURES_+(name ginven in parameter file).dat

  file_catalogue = dir_input + file_catalogue;
  file_random = dir_input + file_random;


}


