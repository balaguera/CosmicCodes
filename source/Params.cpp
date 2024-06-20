////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @class<Params>
 *  @file Params.cpp
 *  @brief Methods of the class Params
 *  @details Reads and administrates input parameters
 *  @author Andrés Balaguera-Antolínez,
 *  @date 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../headers/Params.h"
#include "../headers/cosmo_parameters.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define NOT_USED -1
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This should be independent of reading the parameter file.


void Params::explain_pars(string par_name)
{
  s_message_pars message_pars;
  message_pars.set_par_name(par_name);
  bool v_string=false;
  if(par_name=="N_X"){
    message_pars.set_par_description("Number of bins for dark matter");
    message_pars.set_par_option("Integer > 0");
    message_pars.set_par_default(this->NX);
    v_string=false;// if parameter is a string, false. if a number, true
  }
  if(par_name=="N_Y")
  {
    message_pars.set_par_description("Number of bins for dark matter halos. ");
    message_pars.set_par_option("Integer > 0. If the interpolation is NGP, this value is overriden and transformed to maximum occupation number in cells");
    v_string=false;
  }
  if(par_name=="Redshift" || par_name=="redshift")
  {
    message_pars.set_par_description("Cosmological redshift");
    message_pars.set_par_option("A number > 0 identifiying");
    message_pars.set_par_default(this->redshift);
    v_string=false;
  }
  if(par_name=="IC_index")
  {
    message_pars.set_par_description("Index to characterize the IC of some procedures. ");
    message_pars.set_par_option("Integer");
    v_string=false;
  }
  if(par_name=="Nlambdath")
  {
    message_pars.set_par_description("Number of thersholds used to define cosmic-web types. Used only if test with cwc are available");
    message_pars.set_par_option("Integer > 0");
    v_string=false;
  }
  if(par_name=="lambdath")
  {
    message_pars.set_par_description("Threshold used to define cosmic-web types.");
    message_pars.set_par_option("Integer > 0");
    message_pars.set_par_default(this->lambdath);
    v_string=false;
  }
  if(par_name=="realization")
  {
    message_pars.set_par_description("Index to characterize realizations when producing mocks");
    message_pars.set_par_option("Integer > 0");
    message_pars.set_par_default(this->realization);
    v_string=false;
  }
  if(par_name=="delta_X_min")
  {
    message_pars.set_par_description("Minimum vaue of overdensity for histogram (in linear scale)");
    message_pars.set_par_option("Number > 0");
    message_pars.set_par_default(this->delta_X_min);
    v_string=false;
  }
  if(par_name=="delta_X_max")
  {
    message_pars.set_par_description("Maximum vaue of overdensity for histogram (in linear scale)");
    message_pars.set_par_option("Number > 0");
    message_pars.set_par_default(this->delta_X_max);
    v_string=false;
  }
  if(par_name=="ldelta_X_max")
  {
    message_pars.set_par_description("Maximum vaue of overdensity for histogram (in log(1+delta))");
    message_pars.set_par_option("Number > 0");
    message_pars.set_par_default(this->ldelta_X_max);
    v_string=false;
  }
  if(par_name=="ldelta_X_min")
  {
    message_pars.set_par_description("Minimum vaue of overdensity for histogram (in log(1+delta))");
    message_pars.set_par_option("Number > 0");
    message_pars.set_par_default(this->ldelta_X_min);
    v_string=false;
  }
  if(par_name=="Statistics")
  {
    message_pars.set_par_description("Statistics to measure from input catalog");
    message_pars.set_par_option("For power spectrum: Pk_fkp (FKP estimator), Pk_yb (Yamamoto-Bianchi), Pk_ys(Yamamoto-Scoccimarro)");
    v_string=true;
  }



  message_pars.show(v_string);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Params::init_pars()
{
  this->NX = 200;
  this->NY = 100;
  this->n_sknot_massbin = 200;
  this->n_vknot_massbin = 200;
  this->NMASSbins = 1;
  this->NMASSbins_power = 1;
  this->iMAS_X = 0;
  this->iMAS_X_REF_PDF = 0;
  this->iMAS_X_NEW = 0;
  this->iMAS_Y = 0;
  this->lambdath = 0.0;
  this->lambdath_v = 0.0;
  this->realization = 1;
  this->redshift = 0.0;
  this->Initial_Redshift_DELTA=0;
  this->smscale = 0.0;
  this->delta_Y_min = 0.;
  this->delta_Y_max = 10;
  this->delta_X_min = 0;
  this->delta_X_max = 10.;
  this->ldelta_Y_min = 0;
  this->ldelta_Y_max = 10;
  this->ldelta_X_min = 0.;
  this->ldelta_X_max = 10.;
  this->Input_Directory_Y = "../Input";
  this->Name_Catalog_Y = "UNKNOWN_NAME";
  this->Name_Catalog_Y_new_ref = "UNKNOWN_NAME";
  this->Name_Catalog_Y_HR = "UNKNOWN_NAME";
  this->Name_Catalog_Y_MWEIGHTED = "UNKNOWN_NAME";
  this->Input_Directory_X = "../Input";
  this->Input_Directory_X_REF = "../Input";
  this->Input_Directory_X_NEW = "../Input";
  this->Input_Directory_X_new_ref = "../Input";
  this->Input_Directory_BIAS_KERNEL = "../Input";
  this->Input_Directory_BIAS_KERNEL_TWO = "../Input";
  this->XNAME = "XNAME";
  this->YNAME = "YNAME";
  this->Name_Catalog_X = "UNKNOWN_NAME";
  this->Name_VelFieldx_X = "UNKNOWN_NAME";
  this->Name_VelFieldy_X = "UNKNOWN_NAME";
  this->Name_VelFieldz_X = "UNKNOWN_NAME";
  this->Name_Catalog_X_REF_PDF = "UNKNOWN_NAME";
  this->Name_Catalog_X_NEW = "UNKNOWN_NAME";
  this->Name_Catalog_X_new_ref = "UNKNOWN_NAME";
  this->extra_info = "";
  this->Name_Property_X = "X_PREPERTY";
  this->Name_Property_Y = "Y_PROPERTY";
  this->Output_directory = "../output";
  this->Quantity= "QUANTITY";
  this->Scale_X= "linear";
  this->Scale_Y= "linear";
  this->N_iterations_Kernel = 1;
  this->Iteration_Kernel_average = 1;
  this->Comp_conditional_PDF = false;
  this->Apply_Rankordering = false;
  this->Apply_Rankordering_ab_initio = false;
  this->vel_units_g = "kmps";
  this->IC_index=0;
  this->Redefine_limits= false;
  this->Input_dir_cat_new_ref = "../Input/";
  this->Input_dir_cat_TWO = "../Input/";
  this->dir_output = "../Output/";
  this->file_catalogue = "cat.dat";
  this->file_catalogue_new_ref = "cat.dat";
  this->input_type = "catalog";
  this->input_type_two = "catalog";
  this->delta_grid_file = "delta_file";
  this->delta_grid_file2 = "delta_file";
  this->delta_grid_file3 = "delta_file";
  this->delta_grid_file4 = "delta_file";
  this->ngal_delta = 1;
  this->measure_cross = false;
  this->measure_cross_from_1 = 1;
  this->measure_cross_from_2 = 2;
  // ************************************************************************
  this->Get_Mstellar_function = false;
  this->Get_Luminosity_function = false;
  this->Get_Color_function = false;
  this->LF_estimator = "Vmax_o";
  this->Nbins_color = 0;
  this->Nbins_Mstellar = 0;
  this->Mstellar_min = 0;
  this->Mstellar_max = 1.0;
  this->Color_min = 0;
  this->Color_max = 1.0;
  this->Get_Color_Mag_plane = false;
  this->Get_Random_Catalog = false;
  this->Number_of_random_files = 1;
  // ************************************************************************
  this->sys_of_coord_dm = 1;
  this->N_dm_initial=1;
  this->N_dm_realizations=1;
  this->N_iterations_dm =0;
  this->iteration_ini =0;
  this->i_coord1_dm = NOT_USED;
  this->i_coord2_dm = NOT_USED;
  this->i_coord3_dm = NOT_USED;
  this->i_v1_dm = NOT_USED;
  this->i_v2_dm = NOT_USED;
  this->i_v3_dm = NOT_USED;
  this->i_mass_dm = NOT_USED;
  this->vel_units_g = "kmps";
  this->Redefine_limits=false;
  this->Comp_joint_PDF=false;
  this->Scale_Y="linear";
    // ************************************************************************
  this->sys_of_coord_g = 1;
  this->redshift_space_coords_g =false;
  this->i_coord1_g = NOT_USED;
  this->i_coord2_g = NOT_USED;
  this->i_coord3_g = NOT_USED;
  this->i_v1_g = NOT_USED;
  this->i_v2_g = NOT_USED;
  this->i_v3_g = NOT_USED;
  this->i_mean_density_g = NOT_USED;
  this->i_weight1_g = NOT_USED;
  this->i_weight2_g = NOT_USED;
  this->i_weight3_g = NOT_USED;
  this->i_weight4_g = NOT_USED;
  this->i_mass_g = NOT_USED;
  this->i_rs_g = NOT_USED;
  this->i_spin_g = NOT_USED;
  this->i_vmax_g = NOT_USED;
  this->i_vrms_g = NOT_USED;
  this->i_virial_g = NOT_USED;
  this->i_b_to_a_g = NOT_USED;
  this->i_c_to_a_g = NOT_USED;
  this->i_color_g=NOT_USED;
  this->i_stellar_mass_g=NOT_USED;
  this->i_abs_mag_g=NOT_USED;
  this->i_app_mag_g=NOT_USED;
  this->type_of_object="TRACER";
  this->use_weight1_g = false;
  this->use_weight2_g = false;
  this->use_weight3_g = false;
  this->use_weight4_g = false;
  this->new_Lbox = false;
  this->n_catalogues = 1;
  this->angles_units_g = "D";
  this->MASS_units=1;
  this->NMASSbins_power=1;
  this->NMASSbins=1;
  this->LOGMASSmin=0;
 // ************************************************************************
  this->use_random_catalog = false;
  this->use_random_catalog_cl = false;
  this->use_file_nbar=false;
  this->sys_of_coord_r = 1;
  this->i_coord1_r = NOT_USED;
  this->i_coord2_r = NOT_USED;
  this->i_coord3_r = NOT_USED;
  this->i_mean_density_r = NOT_USED;
  this->i_weight1_r = NOT_USED;
  this->i_weight2_r = NOT_USED;
  this->i_weight3_r = NOT_USED;
  this->i_weight4_r = NOT_USED;
  this->i_color_r=NOT_USED;
  this->i_stellar_mass_r=NOT_USED;
  this->i_abs_mag_r=NOT_USED;
  this->i_app_mag_r=NOT_USED;
  // ************************************************************************
  this->i_mask_pixel= NOT_USED;
  this->i_mask_alpha= NOT_USED;
  this->i_mask_delta= NOT_USED;
  this->i_mask_flag= NOT_USED;
  // ************************************************************************
  this->use_weight1_r = false;
  this->use_weight2_r = false;
  this->use_weight3_r = false;
  this->use_weight4_r = false;
  this->angles_units_r = "D";
  this->file_random = "random";
  this->Name_survey = "survey";
  this->Get_marked_power_spectrum =false;
  this->Get_power_spectrum =false;
  this->Get_cross_power_spectrum =false;
  this->Get_tracer_number_counts=false;
  this->Get_tidal_anisotropy_at_halo=false;
  this->Get_tracer_mass_field=false;
  this->Get_tracer_vmax_field=false;
  this->Get_tracer_spin_field=false;
  this->Get_tracer_local_mach_number=false;
  this->Get_tracer_local_dach_number=false;
  this->Get_tracer_local_dm_density=false;
  this->Get_pearson_coefficient=false;
  this->Get_spearman_coefficient=false;
  this->Get_cell_local_mach_number=false;
  this->Scale_mach_number = static_cast<real_prec>(2.0);
  this->NPROPbins_bam=1;
  this->Number_of_bins_equal_number_tracers = 1;
  this->Number_of_bins_equal_number_tracers_main_property = 1;
  this->set_bins_equal_number_tracers = false;
  this->set_bins_equal_number_tracers_main_property = false;
  this->kmax_tracer_bias = 0;
  this->kmin_tracer_bias = 0;
  this->kmax_tracer_qbias = 0;
  this->kmin_tracer_qbias = 0;
  // ************************************************************************
  this->lambdath_v=0;
  this->Number_of_chunks_new_dm=1;
  this->Name_VelFieldx_X= "file";
  this->Name_VelFieldy_X= "file";
  this->Name_VelFieldz_X= "file";
  this->Apply_Rankordering=false;
    // ************************************************************************
  // parameters for the power spectrum
  this->ndel_window=1;
  this->ndel_data=1;
  this->N_mu_bins=10;
  this->Pest = 20000;
  this->new_los=false;
  this->Nft = 1;
  this->Nft_HR = 1;
  this->Nft_low = 1;
  this->Nft_JK = 1;
  this->Nft_random_collapse = 1;
  this->Distance_fraction=1.0;
  this->Lbox = 100;
  this->Lbox_low = 100;
  this->mass_assignment_scheme = "NGP";
  this->MAS_correction = false;
  this->type_of_binning = "linear";
  this->DeltaKmin = 0;
  this->N_log_bins = 10;
  this->ndel_data = 1;
  this->ndel_window = 1;
  this->N_mu_bins = 10;
  this->FKP_weight = false;
  this->Pest = 1000.0;
  this->SN_correction = false;
  this->FKP_error_bars = false;
  this->FKP_error_bars_exact = false;
  this->nbar_tabulated = false;
  this->constant_depth = false;
  this->Nbins_redshift = 10;
  this->redshift_min_sample = 0;
  this->redshift_max_sample = 1;
  this->N_dndz_bins = 2;
  this->new_N_dndz_bins = 2;
  this->area_survey = 1000;
  this->Healpix_resolution = 4;
  this->file_power = "power";
  this->file_power_log = "power_log";
  this->file_window = "window";
  this->file_dndz = "dndz";
  this->file_power2d = "power2d";
  this->file_power2d_mk = "power2d_mk";
  this->file_bispectrum = "bispectrum";
  //output files for Yamamoto estimator of moments
  this->kmax_y_ds = 0.2;
  //Parameters for the bispectrum
  this->kmax_bk = 0.2;
  this->kmin_bk = 0.1;
  this->Nft_random_collapse =32;
  this->use_fundamental_mode_as_kmin_bk = false;
  this->get_distribution_min_separations=false;
  this->M_exclusion=1;
  this->use_vel_kernel=false;
  this->slength=10;
  this->slengthv=10;
  this->velbias=0;
  this->velbias_dm=0;
  this->velbias_random=0;
  this->dilute_dm_sample = false;
  this->fraction_dilute=1.0;
  this->masskernel_vel=0;
  this->iteration_ini = 0;
  this->get_distribution_min_separations=false;
  this->ic_alias_corrected = false;
  this->ic_input_type = DENSITY;
  this->velbias_random=0.0;
  this->velbias_dm =0.0;
  this->use_vel_kernel=false;
  this->vkernel_exponent=0.;
  this->smscale=0.;
  this->masskernel=0;
  this->masskernel_vel=0; 
  this->N_lines_binary = 0;
  this->Number_of_GRF =1;
  this->Kmax_FA = 0.;
  this->Generate_FA=false;
  this->NVIRIALbins_power=0;
  // ************************************************************************
 // Set of parameter to be deprecated asap
  this->sfac=0; // these are not requested in the parameter file
  this->biasepsilon =0; // these are not requested in the parameter file
  this->biasepsilon2 =0; // these are not requested in the parameter file
  this->biasrhoexp =0; // these are not requested in the parameter file
  this->biasrhoexp2 =0; // these are not requested in the parameter file
  this->deltathH=0; // these are not requested in the parameter file
  this->biasL=0; // these are not requested in the parameter file
 // ************************************************************************
  this->xllc=0; // these are not requested in the parameter file
  this->yllc=0; // these are not requested in the parameter file
  this->zllc=0; // these are not requested in the parameter file
  this->xobs=0; // these are not requested in the parameter file
  this->yobs=0; // these are not requested in the parameter file
  this->zobs=0; // these are not requested in the parameter file
  this->Nmean=0; // these are not requested in the parameter file
  this->devpois=0; // these are not requested in the parameter file
  this->biassign2=0; // these are not requested in the parameter file
  this->biassign=0; // these are not requested in the parameter file
 // ************************************************************************
  this->xmin=0;
  this->ymin=0;
  this->zmin=0;
  this->Xoffset = 0;
  this->Yoffset = 0;
  this->Zoffset = 0;
  this->xmax = 0;
  this->ymax = 0;
  this->zmax = 0;
  // ************************************************************************
  // cosmological parameters
  this->Get_SO_from_BN = false;
  this->om_matter = COSMOPARS::Om_matter;
  this->om_radiation = COSMOPARS::Om_radiation;
  this->om_baryons = COSMOPARS::Om_baryons;
  this->om_cdm = this->om_matter-this->om_baryons;
  this->om_vac = COSMOPARS::Om_vac;
  this->om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->spectral_index = COSMOPARS::spectral_index;
  this->wde_eos = COSMOPARS::wde_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = true;
  this->RR = COSMOPARS::RR;
  // Parameter for auto and marked correlation function
  this->rmin_cf=1;
  this->rmax_cf=100.0;
  this->mark="mass";
  this->use_low_pass_filter=false;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Params::read_pars(string file)
{
  ifstream fin_parameters (file.c_str());
  if (fin_parameters.fail()) 
  { 
    cerr <<"Error at opening the parameter file "<<file<<std::endl;
    cerr <<"Cosmicatlass stops here."<<std::endl; 
    cerr <<"Please check."<<std::endl; 
    cout<<std::endl;
    exit(1); 
  }
  string line_in_file;
  string par_name;
  string equality;
  string par_value;
  while (getline(fin_parameters,line_in_file))
    {
      if (line_in_file[0] != '#' && line_in_file.empty()==0)      // ignore lines starting with hashtag
       	{
	  // read parameter name
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;
          line_string >> equality;
          if (equality != "=")       // check that second word is "="
    	      {
            std::cerr << "Error in parameters file at " << par_name << " = " << "???" << std::endl;
            std::cerr << "Symbol '=' not properly written"<< std::endl; 
            std::cerr << "Please check the parameter file "<<file<< std::endl; 
            std::cerr << "Cosmicatlass stops here."<<file<< std::endl; 
            exit(1);
	        }
        line_string >> par_value; 	  // read parameter value
	 	    if (true==par_value.empty())
          {
            cerr << "Warning from parameter file: " <<endl;
            cerr << "Value of " << par_name << " not specified in " <<file << endl;
            cerr << "Assuming a default value for " << par_name << endl;
            continue;
          }
        else if (par_name == "NX")
          {
            this->NX = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NX));
          }
        else if (par_name == "NY")
          {
            this->NY = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NY));
          }
        else if (par_name == "NY_MASS")
          {
            this->NY_MASS = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NY_MASS));
          }
        else if (par_name == "NY_SAT_FRAC")
          {
            this->NY_SAT_FRAC = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NY_SAT_FRAC));
          }
        else if (par_name == "N_SKNOT_MASSBIN")
          {
            this->n_sknot_massbin = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->n_sknot_massbin));
          }
        else if (par_name == "N_VKNOT_MASSBIN")
          {
                this->n_vknot_massbin = atoi(par_value.c_str());
                this->parameter_number.push_back(make_pair(par_name, this->n_vknot_massbin));
          }
        else if (par_name == "iMAS_X")
          {
          this->iMAS_X = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->iMAS_X));
          }
          else if (par_name == "iMAS_X_REF_PDF")
          {
            this->iMAS_X_REF_PDF = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->iMAS_X_REF_PDF));
          }
          else if (par_name == "iMAS_X_NEW")
          {
            this->iMAS_X_NEW = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->iMAS_X_NEW));
          }
        else if (par_name == "iMAS_Y")
          {
            this->iMAS_Y = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->iMAS_Y));
          }
        else if (par_name == "Unitsim_plabel")
          {
            this->Unitsim_plabel = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->Unitsim_plabel));
          }
        else if (par_name == "Number_of_chunks_new_dm")
          {
            this->Number_of_chunks_new_dm = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->Number_of_chunks_new_dm));
          }
        else if (par_name == "lambdath")
          {
            this->lambdath = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->lambdath));
          }
        else if (par_name == "lambdath_v")
          {
            this->lambdath_v = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->lambdath_v));
          }
        else if (par_name == "Realization")
          {
            this->realization = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->realization));
          }
        else if (par_name == "Redshift")
          {
            this->redshift = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->redshift));
          }
        else if (par_name == "smscale")
          {
            this->smscale = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->smscale));
          }
        else if (par_name == "delta_Y_min")
          {
            this->delta_Y_min = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->delta_Y_min));
          }
        else if (par_name == "delta_Y_max")
          {
            this->delta_Y_max = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->delta_Y_max));
          }
        else if (par_name == "delta_X_min")
          {
            this->delta_X_min = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->delta_X_min));
          }
        else if (par_name == "delta_X_max")
          {
            this->delta_X_max = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->delta_Y_max));
          }
        else if (par_name == "ldelta_Y_min")
          {
            this->ldelta_Y_min = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->ldelta_Y_min));
          }
        else if (par_name == "ldelta_Y_max")
          {
            this->ldelta_Y_max = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->ldelta_Y_max));
          }
        else if (par_name == "ldelta_X_min")
          {
            this->ldelta_X_min = static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->ldelta_X_min));	      
          }
        else if (par_name == "ldelta_X_max")
          {
            this->ldelta_X_max =  static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->ldelta_X_max));
          }
        else if (par_name == "Input_Directory_Y")
          {
              this->Input_Directory_Y = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_Y)); 
          }
        else if (par_name == "Input_Directory_Y_TWO")
          {
              this->Input_Directory_Y_TWO = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_Y_TWO)); 
          }
        else if (par_name == "Name_Catalog_Y")
          {
              this->Name_Catalog_Y = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_Y)); 
          }
        else if (par_name == "Name_Catalog_Y_new_ref")
          {
              this->Name_Catalog_Y_new_ref = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_Y_new_ref)); 
          }
        else if (par_name == "Name_Catalog_Y_HR")
          {
              this->Name_Catalog_Y_HR = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_Y_HR)); 
            
          }
        else if (par_name == "Name_Catalog_Y_MWEIGHTED")
          {
              this->Name_Catalog_Y_MWEIGHTED = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_Y_MWEIGHTED)); 
          }
        else if (par_name == "Input_Directory_X")
          {
              this->Input_Directory_X = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_X)); 
          }
        else if (par_name == "Input_Directory_BIAS_KERNEL")
          {
              this->Input_Directory_BIAS_KERNEL = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_BIAS_KERNEL)); 
          }
        else if (par_name == "Input_Directory_BIAS_KERNEL_TWO")
          {
              this->Input_Directory_BIAS_KERNEL_TWO = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_BIAS_KERNEL_TWO)); 
	    }
      else if (par_name == "Input_Directory_X_REF")
        {
            this->Input_Directory_X_REF = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_X_REF)); 
        }
      else if (par_name == "Input_Directory_X_new_ref")
        {
            this->Input_Directory_X_new_ref = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_X_new_ref)); 
        }
      else if (par_name == "Input_Directory_X_REF_TWO")
        {
            this->Input_Directory_X_REF_TWO = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_X_REF_TWO)); 
        }
      else if (par_name == "Input_Directory_X_NEW")
        {
            this->Input_Directory_X_NEW = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Input_Directory_X_NEW)); 
        }
      else if (par_name == "XNAME")
        {
            this->XNAME = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->XNAME)); 
        }
      else if (par_name == "YNAME")
        {
            this->YNAME = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->YNAME)); 
        }
      else if (par_name == "Name_Catalog_X")
        {
            this->Name_Catalog_X = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_X)); 
        }
      else if (par_name == "Name_VelFieldx_X")
        {
            this->Name_VelFieldx_X = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_VelFieldx_X)); 
        }
      else if (par_name == "Name_VelFieldy_X")
        {
            this->Name_VelFieldy_X = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_VelFieldy_X)); 
        }
      else if (par_name == "Name_VelFieldz_X")
        {
            this->Name_VelFieldz_X = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_VelFieldz_X)); 
        }
      else if (par_name == "Name_redshift_mask")
        {
            this->Name_redshift_mask = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_redshift_mask)); 
        }
      else if (par_name == "Name_binary_mask")
        {
            this->Name_binary_mask = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_binary_mask)); 
        }
      else if (par_name == "IC_index")
        {
          this->IC_index = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->IC_index)); 
        }
      else if (par_name == "Name_Catalog_X_REF_PDF")
        {
            this->Name_Catalog_X_REF_PDF = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_X_REF_PDF)); 
        }
      else if (par_name == "Name_Catalog_X_NEW")
        {
            this->Name_Catalog_X_NEW = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_X_NEW)); 
        }
      else if (par_name == "Name_Catalog_X_new_ref")
        {
            this->Name_Catalog_X_new_ref = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_X_new_ref)); 
        }
      else if (par_name == "Name_Catalog_X_NEW_TWO")
        {
            this->Name_Catalog_X_NEW_TWO = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Catalog_X_NEW_TWO)); 
        }
      else if (par_name == "Name_Property_X")
        {
            this->Name_Property_X = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Property_X)); 
        }
      else if (par_name == "Name_Property_Y")
        {
            this->Name_Property_Y = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Name_Property_Y)); 
        }
      else if (par_name == "Output_directory")
        {
            this->Output_directory = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Output_directory)); 
        }
      else if (par_name == "Quantity")
        {
            this->Quantity= par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Quantity)); 
        }
      else if (par_name == "Scale_X")
        {
            this->Scale_X= par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Scale_X)); 
        }
      else if (par_name == "Scale_Y")
        {
            this->Scale_Y= par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->Scale_Y)); 
        }
      else if (par_name == "N_iterations_Kernel")
        {
          this->N_iterations_Kernel = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->N_iterations_Kernel)); 
        }
      else if (par_name == "Iteration_Kernel_average")
        {
          this->Iteration_Kernel_average = static_cast<ULONG>(atoi(par_value.c_str()));
          this->parameter_number.push_back(make_pair(par_name, this->Iteration_Kernel_average)); 
        }
      else if (par_name == "iteration_ini")
        {
          this->iteration_ini = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->iteration_ini)); 
        }
      else if (par_name == "EXTRA_INFO")
        {
          this->extra_info = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->extra_info)); 
        }
      else if (par_name == "MASS_units")
        {
          this->MASS_units = atof(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->MASS_units)); 
        }
      else if (par_name == "LOGMASSmin")
        {
          this->LOGMASSmin = atof(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->LOGMASSmin)); 
        }
      else if (par_name == "VMAXmax")
        {
          this->VMAXmax = atof(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->VMAXmax)); 
        }
      else if (par_name == "VMAXmin")
        {
	      this->VMAXmin = atof(par_value.c_str());
        this->parameter_number.push_back(make_pair(par_name, this->VMAXmin));
      }
          else if (par_name == "VRMSmax")
            {
              this->VRMSmax = atof(par_value.c_str());
              this->parameter_number.push_back(make_pair(par_name, this->VRMSmax));
            }
          else if (par_name == "VRMSmin")
            {
              this->VRMSmin = atof(par_value.c_str());
              this->parameter_number.push_back(make_pair(par_name, this->VRMSmin));
            }
          else if (par_name == "RSmax")
          {
            this->RSmax = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->RSmax)); 
          }
          else if (par_name == "RSmin")
          {
            this->RSmin = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->RSmin)); 
          }
          else if (par_name == "CONCENTRATIONmax")
          {
            this->CONCENTRATIONmax = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->CONCENTRATIONmax)); 
          }
        else if (par_name == "CONCENTRATIONmin")
          {
            this->CONCENTRATIONmin = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->CONCENTRATIONmin)); 
          }
        else if (par_name == "SPINmax")
          {
            this->SPINmax = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->SPINmax)); 
          }
        else if (par_name == "SPINmin")
          {
            this->SPINmin = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->SPINmin)); 
          }
        else if (par_name == "Number_of_bins_equal_number_tracers")
          {
            this->Number_of_bins_equal_number_tracers  = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->Number_of_bins_equal_number_tracers)); 
          }
        else if (par_name == "Number_of_bins_equal_number_tracers_main_property")
          {
            this->Number_of_bins_equal_number_tracers_main_property  = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->Number_of_bins_equal_number_tracers_main_property)); 
          }
        else if (par_name == "set_bins_equal_number_tracers")
          {
            if (par_value=="true")this->set_bins_equal_number_tracers = true;
            else this->set_bins_equal_number_tracers = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->set_bins_equal_number_tracers)); 
          }
          else if (par_name == "set_bins_equal_number_tracers_main_property")
            {
              if (par_value=="true")this->set_bins_equal_number_tracers_main_property = true;
              else this->set_bins_equal_number_tracers_main_property = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->set_bins_equal_number_tracers_main_property));
            }
          else if (par_name == "dilute_dm_sample")
          {
            if (par_value=="true")this->dilute_dm_sample = true;
            else this->dilute_dm_sample = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->dilute_dm_sample)); 
          }
        else if (par_name == "fraction_dilute")
          {
            this->fraction_dilute = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->fraction_dilute)); 
          }
        else if (par_name == "Nbins_hist")
          {
            this->Nbins_hist = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->Nbins_hist)); 
          }
              else if (par_name == "LOGMASSmax")
          {
            this->LOGMASSmax = atof(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->LOGMASSmax)); 
          }
        else if (par_name == "NMASSbins")
          {
            this->NMASSbins = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NMASSbins)); 
          }
        else if (par_name == "NMASSbins_mf")
          {
            this->NMASSbins_mf = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NMASSbins_mf)); 
          }
              else if (par_name == "NPROPbins_bam")
          {
            this->NPROPbins_bam = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NPROPbins_bam)); 
          }
              else if (par_name == "file_bin_x_coord" )
          {
            this->file_bin_x_coord = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->file_bin_x_coord)); 
          }
        else if (par_name == "file_bin_y_coord" )
          {
            this->file_bin_y_coord = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->file_bin_y_coord)); 
          }
        else if (par_name == "file_bin_z_coord" )
          {
            this->file_bin_z_coord = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->file_bin_z_coord)); 
          }
        else if (par_name == "N_lines_binary")
          {
            this->N_lines_binary = static_cast<ULONG>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->N_lines_binary)); 
          }
        else if (par_name == "Input_dir_cat")
          {
              this->Input_dir_cat = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_dir_cat)); 
          }
        else if (par_name == "Input_dir_cat_new_ref")
          {
              this->Input_dir_cat_new_ref = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_dir_cat_new_ref)); 
          }
        else if (par_name == "Input_dir_cat_TWO")
          {
              this->Input_dir_cat_TWO = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->Input_dir_cat_TWO)); 
          }
        else if (par_name == "Type_of_object")
          {
              this->type_of_object = par_value.c_str();
            this->parameter_string.push_back(make_pair(par_name, this->type_of_object)); 
          }
        else if (par_name == "Number_of_GRF")
          {
            this->Number_of_GRF = static_cast<ULONG>(atoi(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->Number_of_GRF)); 
          }
        else if (par_name == "Kmax_FA")
          {
            this->Kmax_FA= static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->Kmax_FA)); 
          }
        else if (par_name == "Scale_mach_number")
          {
            this->Scale_mach_number= static_cast<real_prec>(atof(par_value.c_str()));
            this->parameter_number.push_back(make_pair(par_name, this->Scale_mach_number)); 
          }
        else if (par_name == "Comp_conditional_PDF")
          {
            if (par_value=="true")this->Comp_conditional_PDF = true;
            else this->Comp_conditional_PDF = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Comp_conditional_PDF)); 
          }
        else if (par_name == "Normalize_IC_to_initial_redshift")
          {
            if (par_value=="true")this->Normalize_IC_to_initial_redshift = true;
            else this->Normalize_IC_to_initial_redshift = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Normalize_IC_to_initial_redshift)); 
          }
        else if (par_name == "Apply_Rankordering")
          {
            if (par_value=="true")this->Apply_Rankordering = true;
            else this->Apply_Rankordering = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Apply_Rankordering)); 
          }
        else if (par_name == "Apply_Rankordering_ab_initio")
          {
            if (par_value=="true")this->Apply_Rankordering_ab_initio = true;
            else this->Apply_Rankordering_ab_initio = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Apply_Rankordering_ab_initio)); 
          }
        else if (par_name == "Redefine_limits")
          {
            if (par_value=="true")this->Redefine_limits= true;
            else this->Redefine_limits = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Redefine_limits)); 
          }
        else if (par_name == "Write_Scatter_Plot")
          {
            if (par_value=="true")this->Write_Scatter_Plot= true;
            else this->Write_Scatter_Plot = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Write_Scatter_Plot)); 
          }
        else if (par_name == "Write_PDF_number_counts")
          {
            if (par_value=="true")this->Write_PDF_number_counts= true;
            else this->Write_PDF_number_counts = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Write_PDF_number_counts)); 
          }
         else if (par_name == "use_real_and_redshift_space")
          {
            if (par_value=="true")this->use_real_and_redshift_space= true;
            else this->use_real_and_redshift_space = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->use_real_and_redshift_space));
          }
         else if (par_name == "Get_tracer_number_counts")
          {
            if (par_value=="true")this->Get_tracer_number_counts= true;
            else this->Get_tracer_number_counts = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_number_counts)); 
          }
          else if (par_name == "Get_marked_power_spectrum")
            {
              if (par_value=="true")this->Get_marked_power_spectrum= true;
              else this->Get_marked_power_spectrum = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_marked_power_spectrum));
            }
           else if (par_name == "Get_power_spectrum")
            {
              if (par_value=="true")this->Get_power_spectrum= true;
              else this->Get_power_spectrum = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_power_spectrum));
            }
           else if (par_name == "Get_cross_power_spectrum")
            {
              if (par_value=="true")this->Get_cross_power_spectrum= true;
              else this->Get_cross_power_spectrum = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_cross_power_spectrum));
            }
              else if (par_name == "Get_pearson_coefficient")
            {
              if (par_value=="true")this->Get_pearson_coefficient= true;
              else this->Get_pearson_coefficient = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_pearson_coefficient));
            }
              else if (par_name == "Get_spearman_coefficient")
            {
              if (par_value=="true")this->Get_spearman_coefficient= true;
              else this->Get_spearman_coefficient = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_spearman_coefficient));
            }
              else if (par_name == "Get_tidal_anisotropy_at_halo")
            {
              if (par_value=="true")this->Get_tidal_anisotropy_at_halo= true;
              else this->Get_tidal_anisotropy_at_halo = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tidal_anisotropy_at_halo));
            }
          else if (par_name == "Get_peak_height_at_halo")
            {
              if (par_value=="true")this->Get_peak_height_at_halo= true;
              else this->Get_peak_height_at_halo = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_peak_height_at_halo));
            }
          else if (par_name == "Get_tracer_mass_field")
          {
            if (par_value=="true")this->Get_tracer_mass_field= true;
            else this->Get_tracer_mass_field = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_mass_field)); 
          }
        else if (par_name == "Get_tracer_vmax_field")
          {
            if (par_value=="true")this->Get_tracer_vmax_field= true;
            else this->Get_tracer_vmax_field = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_vmax_field)); 
          }
        else if (par_name == "Get_tracer_spin_field")
          {
            if (par_value=="true")this->Get_tracer_spin_field= true;
            else this->Get_tracer_spin_field = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_spin_field)); 
          }
              else if (par_name == "Get_tracer_local_mach_number")
                {
                  if (par_value=="true")this->Get_tracer_local_mach_number= true;
                  else this->Get_tracer_local_mach_number = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_local_mach_number)); 
          }
              else if (par_name == "Get_tracer_local_dach_number")
                {
                  if (par_value=="true")this->Get_tracer_local_dach_number= true;
                  else this->Get_tracer_local_dach_number = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_local_dach_number));
            }
              else if (par_name == "Get_cell_local_mach_number")
                {
                  if (par_value=="true")this->Get_cell_local_mach_number= true;
                  else this->Get_cell_local_mach_number = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_cell_local_mach_number)); 
                }
              else if (par_name == "Get_tracer_bias")
            {
              if (par_value=="true")this->Get_tracer_bias= true;
              else this->Get_tracer_bias = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_bias));
            }
          else if (par_name == "Get_tracer_quadratic_bias")
            {
              if (par_value=="true")this->Get_tracer_quadratic_bias= true;
              else this->Get_tracer_quadratic_bias = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_quadratic_bias));
            }
          else if (par_name == "Get_tracer_relative_bias")
            {
              if (par_value=="true")this->Get_tracer_relative_bias= true;
              else this->Get_tracer_relative_bias = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_relative_bias));
            }
          else if (par_name == "Get_tracer_local_dm_density")
            {
              if (par_value=="true")this->Get_tracer_local_dm_density= true;
              else this->Get_tracer_local_dm_density = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_local_dm_density));
            }
          else if (par_name == "Get_PCA")
            {
              if (par_value=="true")this->Get_PCA= true;
              else this->Get_PCA = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_PCA));
            }
          else if (par_name == "kmax_tracer_bias"){
              this->kmax_tracer_bias= static_cast<real_prec>(atof(par_value.c_str()));
              this->parameter_number.push_back(make_pair(par_name, this->kmax_tracer_bias));
          }
          else if (par_name == "kmin_tracer_bias"){
              this->kmin_tracer_bias= static_cast<real_prec>(atof(par_value.c_str()));
              this->parameter_number.push_back(make_pair(par_name, this->kmin_tracer_bias));
          }
              else if (par_name == "kmax_tracer_qbias"){
              this->kmax_tracer_qbias= static_cast<real_prec>(atof(par_value.c_str()));
              this->parameter_number.push_back(make_pair(par_name, this->kmax_tracer_qbias));
          }
              else if (par_name == "kmin_tracer_qbias"){
              this->kmin_tracer_qbias= static_cast<real_prec>(atof(par_value.c_str()));
              this->parameter_number.push_back(make_pair(par_name, this->kmin_tracer_qbias));
          }
           else if (par_name == "Get_local_overdensity")
            {
              if (par_value=="true")this->Get_local_overdensity= true;
              else this->Get_local_overdensity = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_local_overdensity));
            }
          else if (par_name == "Get_prop_function_tracer")
	    {
	      if (par_value=="true")this->Get_prop_function_tracer= true;
	      else this->Get_prop_function_tracer = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Get_prop_function_tracer)); 
	    }
	  else if (par_name == "Get_prop_function_tracer_cwt")
	    {
	      if (par_value=="true")this->Get_prop_function_tracer_cwt= true;
	      else this->Get_prop_function_tracer_cwt = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Get_prop_function_tracer_cwt)); 
	    }
	  else if (par_name == "Convert_Density_to_Delta_X")
	    {
	      if (par_value=="true")this->Convert_Density_to_Delta_X= true;
	      else this->Convert_Density_to_Delta_X = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Convert_Density_to_Delta_X)); 
	    }
	  else if (par_name == "Convert_Density_to_Delta_Y")
	    {
	      if (par_value=="true")this->Convert_Density_to_Delta_Y= true;
	      else this->Convert_Density_to_Delta_Y = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Convert_Density_to_Delta_Y)); 
	    }
	  else if (par_name == "Comp_joint_PDF")
	    {
	      if (par_value=="true")this->Comp_joint_PDF = true;
	      else this->Comp_joint_PDF = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Comp_joint_PDF)); 
	    }
	  else if (par_name == "write_files_for_histograms")
	    {
	      if (par_value=="true")this->write_files_for_histograms = true;
	      else this->write_files_for_histograms = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->write_files_for_histograms)); 
	    }
          else if (par_name == "use_low_pass_filter")
            {
              if (par_value=="true")this->use_low_pass_filter = true;
              else this->use_low_pass_filter = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_low_pass_filter)); 
            }
	  else if (par_name == "Generate_FA")
	    {
	      if (par_value=="true")this->Generate_FA = true;
	      else this->Generate_FA = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->Generate_FA)); 
	    }
#ifdef _USE_CWC_
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 't' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality; 	  // check that second character is "="
	      if (equality != "=")
      		{
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
      		}
	      string par_value0;
	      int ending=0;
	      int n_m_lims=0;
	      do
		{
		  line_string >> par_value0;  //read value
		  n_m_lims++;
		  if(par_value0=="END")
		    ending=1;
		  this->cwt_used.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
	      this->cwt_used.pop_back(); //delete the last element, which is END
            }
#endif
#ifdef _USE_CWC_V_
          else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'v' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality; 	  // check that second character is "="
              if (equality != "=")
                {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
              string par_value0;
              int ending=0;
              int n_m_lims=0;
              do{
                line_string >> par_value0;  //read value
                n_m_lims++;
                if(par_value0=="END")
                  ending=1;
                this->cwv_used.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->cwv_used.pop_back(); //delete the last element, which is END
            }
#endif
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'x' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality; 	  // check that second character is "="
	      if (equality != "=") {
	       	cerr << "Error in parameters file "<<file<<" at " << par_name << " = " << "???" << endl;
		cerr << "Please check parameter file "<< std::endl; 
      		exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
    		line_string >> par_value0;  //read value
  	   	if(par_value0=="END")ending=1;
    	 	this->output_at_iteration.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->output_at_iteration.pop_back(); //delete the last element, which is END
            }
#ifdef _USE_TWO_REFS_MOCKS_
	  
#ifdef _SLICS_
	  
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'l' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
	       	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		line_string >> par_value0;  //read value
		if(par_value0=="END")ending=1;
		this->list_new_dm_fields.push_back(atoi(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_new_dm_fields.pop_back(); //delete the last element, which is END
	    }
	  else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'g' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
    	   	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		line_string >> par_value0;  //read value
		if(par_value0=="END")ending=1;
		this->list_bias_references.push_back(atoi(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_bias_references.pop_back(); //delete the last element, which is END
	    }
	  
	  // These arrays are inizialized here and updated in teh Bam.pp file 
	  this->Number_of_references=this->list_bias_references.size();
	  this->files_bias_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
            this->files_bias_references[i]="/net/deimos/scratch/balaguera/Numerics/BAM/Output_SLICS/CAL_R"+to_string(this->list_bias_references[i])+"GAL_TWEB/Bam_Bias.dat";
	  this->files_kernel_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
            this->files_kernel_references[i]="/net/deimos/scratch/balaguera/Numerics/BAM/Output_SLICS/CAL_R"+to_string(this->list_bias_references[i])+"GAL_TWEB/Bam_Kernel.dat";
	  this->files_tracer_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_tracer_references[i]=this->Input_dir_cat+"Halos_SLICS_LOS"+to_string(this->list_bias_references[i])+".txt";
	  this->files_tracer_field_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_tracer_field_references[i]=this->Input_Directory_Y+"SLICS_HALOS_LOS"+to_string(this->list_bias_references[i])+"_Nres192_MAS0.dat";
	  this->files_dm_references.resize(this->Number_of_references);
	  for(int i=0;i<this->Number_of_references;++i)
	    this->files_dm_references[i]=this->Input_Directory_X_NEW+"densDMALPTrS20.0TETCICz1.041G192V505.0S"+to_string(this->list_bias_references[i])+".dat";
	  this->Number_of_new_mocks=this->list_new_dm_fields.size();
	  this->files_new_dm_fields.resize(this->Number_of_new_mocks);
	  for(int i=0;i<this->Number_of_new_mocks;++i)
	    this->files_new_dm_fields[i]=this->Input_Directory_X_NEW+"densDMALPTrS20.0TETCICz1.041G192V505.0S"+to_string(this->list_new_dm_fields[i])+".dat";
	  
#elif defined _UNITSIM_  || defined (_ABACUS_)
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->Number_of_references=2;
#else
	  this->Number_of_references=1;
#endif
	  this->files_bias_references.resize(this->Number_of_references);
	  this->files_bias_references[0]=this->Input_Directory_BIAS_KERNEL+"Bam_Bias.dat";
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->files_bias_references[1]=this->Input_Directory_BIAS_KERNEL_TWO+"Bam_Bias.dat";
#endif
	  this->files_kernel_references.resize(this->Number_of_references);
	  this->files_kernel_references[0]=this->Input_Directory_BIAS_KERNEL+"Bam_Kernel.dat";
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->files_kernel_references[1]=this->Input_Directory_BIAS_KERNEL_TWO+"Bam_Kernel.dat";
#endif
	  this->files_tracer_references.resize(this->Number_of_references);
	  this->files_tracer_references[0]=this->Input_dir_cat+this->file_catalogue;
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->files_tracer_references[1]=this->Input_dir_cat_TWO+this->file_catalogue;
#endif
	  this->files_dm_references.resize(this->Number_of_references);
	  this->files_dm_references[0]=this->Input_Directory_X_REF+this->Name_Catalog_X;
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->files_dm_references[1]=this->Input_Directory_X_REF_TWO+this->Name_Catalog_X;// FOR Unitsim, the files will have the smae name regardless normal or Inv
#endif
	  this->files_tracer_field_references.resize(this->Number_of_references);
	  this->files_tracer_field_references[0]=this->Input_Directory_Y+this->Name_Catalog_Y;
#if defined _USE_TWO_REFS_MOCKS_ASSIGNMENT_  || defined _USE_TWO_REFS_MOCKS_
	  this->files_tracer_field_references[1]=this->Input_Directory_Y_TWO+this->Name_Catalog_Y;// FOR Unitsim, the files will have the smae name regardless normal or Inv
#endif
	  this->Number_of_new_mocks=1;
	  this->list_new_dm_fields.resize(this->Number_of_new_mocks);
	  this->files_new_dm_fields.resize(this->Number_of_new_mocks);
	  this->files_new_dm_fields[0]=this->Input_Directory_X_NEW+this->Name_Catalog_X_NEW; 
#endif
#endif
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'm' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->MASSbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->MASSbins_min.pop_back(); //delete the last element, which is END
              this->NMASSbins_power=this->MASSbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->MASSbins_min));
            }
	  else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'n' && line_in_file[1] == '_'))
	      {
		stringstream line_string (line_in_file);
		line_string << line_in_file;
		line_string >> par_name;   // read the first character, the name of the parameter
		line_string >> equality;    // check that second character is "="
		if (equality != "=") {
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
		}
		string par_value0;
		int ending=0;
		do{
		  line_string >> par_value0;  //read value
		  if(par_value0=="END")ending=1;
		  this->MASSbins_max.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
		this->MASSbins_max.pop_back(); //delete the last element, which is END
                this->parameter_vectors.push_back(make_pair(par_name, this->MASSbins_max));
                if(this->NMASSbins_power!=this->MASSbins_max.size())
                  {
		    cerr << "Error in parameters file. MASSbins_max has not the same dimension as MASSbins_min."<< endl;
		    cerr << "Check parameter file" <<endl; exit(1);
                  }
              }
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'a' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->VMAXbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->VMAXbins_min.pop_back(); //delete the last element, which is END
              this->NVMAXbins_power=this->VMAXbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->VMAXbins_min));

            }
	  else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'b' && line_in_file[1] == '_'))
	      {
            stringstream line_string (line_in_file);
            line_string << line_in_file;
            line_string >> par_name;   // read the first character, the name of the parameter
            line_string >> equality;    // check that second character is "="
            if (equality != "=") {
              cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
              cerr << "Using a default value for " << par_name << endl; exit(1);
            }
            string par_value0;
            int ending=0;
            do{
              line_string >> par_value0;  //read value
              if(par_value0=="END")ending=1;
              this->VMAXbins_max.push_back(atof(par_value0.c_str()));
            }while(ending!=1);
            this->VMAXbins_max.pop_back(); //delete the last element, which is END
            this->parameter_vectors.push_back(make_pair(par_name, this->VMAXbins_max));
            if(this->NVMAXbins_power!=this->VMAXbins_max.size())
             {
               cerr << "Error in parameters file. The container VMAXbins_max has not the same dimension as VMAXbins_min."<< endl;
               cerr << "Check parameter file" <<endl; exit(1);
             }
          }
	  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'w' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->VRMSbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->VRMSbins_min.pop_back(); //delete the last element, which is END
              this->NVRMSbins_power=this->VRMSbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->VRMSbins_min));
            }
	  else if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'z' && line_in_file[1] == '_'))
	      {
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;   // read the first character, the name of the parameter
          line_string >> equality;    // check that second character is "="
          if (equality != "=") {
            cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
            cerr << "Using a default value for " << par_name << endl; exit(1);
          }
          string par_value0;
          int ending=0;
          do{
            line_string >> par_value0;  //read value
            if(par_value0=="END")ending=1;
            this->VRMSbins_max.push_back(atof(par_value0.c_str()));
          }while(ending!=1);
          this->VRMSbins_max.pop_back(); //delete the last element, which is END
                      this->parameter_vectors.push_back(make_pair(par_name, this->VRMSbins_max));
                      if(this->NVRMSbins_power!=this->VRMSbins_max.size())
                        {
              cerr << "Error in parameters file. The container VRMSbins_max has not the same dimension as VRMSbins_min."<< endl;
              cerr << "Check parameter file" <<endl; exit(1);
                  }
              }
      // ************************************************************************************************************

      else if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'f' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->VIRIALbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->VIRIALbins_min.pop_back(); //delete the last element, which is END
              this->NVIRIALbins_power=this->VIRIALbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->VIRIALbins_min));
            }
	  else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'g' && line_in_file[1] == '_'&& line_in_file[2]=='V'))
	      {
		stringstream line_string (line_in_file);
		line_string << line_in_file;
		line_string >> par_name;   // read the first character, the name of the parameter
		line_string >> equality;    // check that second character is "="
		if (equality != "=") {
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
		}
		string par_value0;
		int ending=0;
		do{
		  line_string >> par_value0;  //read value
		  if(par_value0=="END")ending=1;
		  this->VIRIALbins_max.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
		this->VIRIALbins_max.pop_back(); //delete the last element, which is END
        this->parameter_vectors.push_back(make_pair(par_name, this->VIRIALbins_max));

        if(this->NVIRIALbins_power!=this->VIRIALbins_max.size())
		  {
		    cerr << "Error in parameters file. The container VIRIALbins_max has not the same dimension as VIRIALbins_min."<< endl;
		    cerr<<"Expected "<<this->NVIRIALbins_power<<". Read "<<this->VIRIALbins_max.size()<<"  "<<this->VIRIALbins_min.size()<<endl;
		    cerr << "Check parameter file" <<endl; exit(1);
		  }
         }
      // ************************************************************************************************************
      // ************************************************************************************************************

          else if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'c' && line_in_file[1] == '_'))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->BTOAbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->BTOAbins_min.pop_back(); //delete the last element, which is END
              this->NBTOAbins_power=this->BTOAbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->BTOAbins_min));
            }
	  else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'd' && line_in_file[1] == '_'))
	      {
		stringstream line_string (line_in_file);
		line_string << line_in_file;
		line_string >> par_name;   // read the first character, the name of the parameter
		line_string >> equality;    // check that second character is "="
		if (equality != "=") {
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
		}
		string par_value0;
		int ending=0;
		do{
		  line_string >> par_value0;  //read value
		  if(par_value0=="END")ending=1;
		  this->BTOAbins_max.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
		this->BTOAbins_max.pop_back(); //delete the last element, which is END
                this->parameter_vectors.push_back(make_pair(par_name, this->BTOAbins_max));
                if(this->NBTOAbins_power!=this->BTOAbins_max.size())
                  {
		    cerr << "Error in parameters file. The container BTOAbins_max has not the same dimension as BTOAbins_min."<< endl;
		    cerr << "Check parameter file" <<endl; exit(1);
                  }
              }
      // ************************************************************************************************************
       else if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'h' && line_in_file[1] == '_'))
         {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
		this->CTOAbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->CTOAbins_min.pop_back(); //delete the last element, which is END
              this->NCTOAbins_power=this->CTOAbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->CTOAbins_min));
            }
	  else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'j' && line_in_file[1] == '_'))
	      {
		stringstream line_string (line_in_file);
		line_string << line_in_file;
		line_string >> par_name;   // read the first character, the name of the parameter
		line_string >> equality;    // check that second character is "="
		if (equality != "=") {
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
		}
		string par_value0;
		int ending=0;
		do{
		  line_string >> par_value0;  //read value
		  if(par_value0=="END")ending=1;
		  this->CTOAbins_max.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
		this->CTOAbins_max.pop_back(); //delete the last element, which is END
                this->parameter_vectors.push_back(make_pair(par_name, this->CTOAbins_max));
                if(this->NCTOAbins_power!=this->CTOAbins_max.size())
                  {
                    cerr << "Error in parameters file. The container CTOAbins_max has not the same dimension as CTOAbins_min."<< endl;
                    cerr << "Check parameter file" <<endl; exit(1);
                  }
              }
      // ************************************************************************************************************
          else
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'o' && line_in_file[1] == '_'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
		cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		line_string >> par_value0;  //read value
		if(par_value0=="END")ending=1;
		this->MASScuts.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->MASScuts.pop_back(); //delete the last element, which is END
	      this->NMASScuts_power=this->MASScuts.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->MASScuts));
	    }
      // ************************************************************************************************************
      else
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'v' && line_in_file[1] == '_' && line_in_file[2] == 'R'))
        {
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;   // read the first character, the name of the parameter
          line_string >> equality;    // check that second character is "="
          if (equality != "=") {
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
          }
          string par_value0;
          int ending=0;
          do{
        line_string >> par_value0;  //read value
        if(par_value0=="END")ending=1;
        this->RANDOMfiles.push_back(par_value0.c_str());
          }while(ending!=1);
          this->RANDOMfiles.pop_back(); //delete the last element, which is END
          this->NRANDOMfiles=this->RANDOMfiles.size();
          this->parameter_vector_string.push_back(make_pair(par_name, this->RANDOMfiles));
        }
     // ************************************************************************************************************
          // READ THE PARAMETER FOR TRACER RS:
          else
	    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'k' && line_in_file[1] == '_'))
	      {
		stringstream line_string (line_in_file);
		line_string << line_in_file;
		line_string >> par_name;   // read the first character, the name of the parameter
		line_string >> equality;    // check that second character is "="
		if (equality != "=") {
		  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		  cerr << "Using a default value for " << par_name << endl; exit(1);
		}
		string par_value0;
		int ending=0;
		do{
		  line_string >> par_value0;  //read value
		  if(par_value0=="END")ending=1;
		  this->RSbins_min.push_back(atof(par_value0.c_str()));
		}while(ending!=1);
		this->RSbins_min.pop_back(); //delete the last element, which is END
		this->NRSbins_power=this->RSbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->RSbins_min));
              }
        else
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'r' && line_in_file[1] == '_'))
        {
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;   // read the first character, the name of the parameter
          line_string >> equality;    // check that second character is "="
          if (equality != "=") {
            cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
            cerr << "Using a default value for " << par_name << endl; exit(1);
          }
          string par_value0;
          int ending=0;
          do{
            line_string >> par_value0;  //read value
            if(par_value0=="END")ending=1;
            this->RSbins_max.push_back(atof(par_value0.c_str()));
          }while(ending!=1);
          this->RSbins_max.pop_back(); //delete the last element, which is END
                  this->parameter_vectors.push_back(make_pair(par_name, this->RSbins_max));
          if(this->NRSbins_power!=this->RSbins_max.size())
            {
              cerr << "Error in parameters file. The container RSbins_max has not the same dimension as RSbins_min."<< endl;
              cerr << "Check parameter file" <<endl; exit(1);
            }
        }
       // ************************************************************************************************************
          // READ THE PARAMETER FOR TRACER Rvir:
          else
        if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'k' && line_in_file[1] == 'v'&& line_in_file[2] == '_'))
          {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
          this->RVIRbins_min.push_back(atof(par_value0.c_str()));
        }while(ending!=1);
        this->RVIRbins_min.pop_back(); //delete the last element, which is END
        this->NRVIRbins_power=this->RVIRbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->RVIRbins_min));
              }

        else
      if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'r' && line_in_file[1] == 'v' && line_in_file[2] == '_'))
        {
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;   // read the first character, the name of the parameter
      line_string >> equality;    // check that second character is "="
      if (equality != "=") {
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      string par_value0;
      int ending=0;
      do{
        line_string >> par_value0;  //read value
        if(par_value0=="END")ending=1;
        this->RVIRbins_max.push_back(atof(par_value0.c_str()));
      }while(ending!=1);
      this->RVIRbins_max.pop_back(); //delete the last element, which is END
      this->parameter_vectors.push_back(make_pair(par_name, this->RVIRbins_max));
      }
     // ************************************************************************************************************
      // ************************************************************************************************************
     // READ THE PARAMETER FOR TRACER CONCENTRATION:
          else
        if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'k' && line_in_file[1] == 'c'&& line_in_file[2] == '_'))
          {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
          this->CONCENTRATIONbins_min.push_back(atof(par_value0.c_str()));
        }while(ending!=1);
        this->CONCENTRATIONbins_min.pop_back(); //delete the last element, which is END
        this->NCONCENTRATIONbins_power=this->CONCENTRATIONbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->CONCENTRATIONbins_min));
              }

        else
      if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'r' && line_in_file[1] == 'c' && line_in_file[2] == '_'))
        {
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;   // read the first character, the name of the parameter
      line_string >> equality;    // check that second character is "="
      if (equality != "=") {
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      string par_value0;
      int ending=0;
      do{
        line_string >> par_value0;  //read value
        if(par_value0=="END")ending=1;
        this->CONCENTRATIONbins_max.push_back(atof(par_value0.c_str()));
      }while(ending!=1);
      this->CONCENTRATIONbins_max.pop_back(); //delete the last element, which is END
      this->parameter_vectors.push_back(make_pair(par_name, this->CONCENTRATIONbins_max));
      }
      // ************************************************************************************************************
          // READ THE PARAMETER FOR TRACER SPIN:
	      else
		if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'q' && line_in_file[1] == '_'))
		  {
		    stringstream line_string (line_in_file);
		    line_string << line_in_file;
		    line_string >> par_name;   // read the first character, the name of the parameter
		    line_string >> equality;    // check that second character is "="
		    if (equality != "=") {
		      cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		      cerr << "Using a default value for " << par_name << endl; exit(1);
		    }
		    string par_value0;
		    int ending=0;
		    do{
		      line_string >> par_value0;  //read value
		      if(par_value0=="END")ending=1;
		      this->SPINbins_min.push_back(atof(par_value0.c_str()));
		    }while(ending!=1);
		    this->SPINbins_min.pop_back(); //delete the last element, which is END
		    this->NSPINbins_power=this->SPINbins_min.size();
                    this->parameter_vectors.push_back(make_pair(par_name, this->SPINbins_min));
                  }
		else
		  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'p' && line_in_file[1] == '_'))
		    {
		      stringstream line_string (line_in_file);
		      line_string << line_in_file;
		      line_string >> par_name;   // read the first character, the name of the parameter
		      line_string >> equality;    // check that second character is "="
		      if (equality != "=") {
			cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
			cerr << "Using a default value for " << par_name << endl; exit(1);
		      }
		      string par_value0;
		      int ending=0;
		      do{
			line_string >> par_value0;  //read value
			if(par_value0=="END")ending=1;
			this->SPINbins_max.push_back(atof(par_value0.c_str()));
		      }while(ending!=1);
		      this->SPINbins_max.pop_back(); //delete the last element, which is END
                      this->parameter_vectors.push_back(make_pair(par_name, this->SPINbins_max));

		      if(this->NSPINbins_power!=this->SPINbins_max.size())
			{
			  cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as SPINbins_min."<< endl;
			  cerr << "Check parameter file" <<endl; exit(1);
			}
           }
      // ************************************************************************************************************
      // ************************************************************************************************************

      // READ THE PARAMETER FOR TRACER SPIN-BULLOCK:
      else
    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'q' && line_in_file[1] == 'b' && line_in_file[2] == '_'))
      {
        stringstream line_string (line_in_file);
        line_string << line_in_file;
        line_string >> par_name;   // read the first character, the name of the parameter
        line_string >> equality;    // check that second character is "="
        if (equality != "=") {
          cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
          cerr << "Using a default value for " << par_name << endl; exit(1);
        }
        string par_value0;
        int ending=0;
        do{
          line_string >> par_value0;  //read value
          if(par_value0=="END")ending=1;
          this->SPINBULLOCKbins_min.push_back(atof(par_value0.c_str()));
        }while(ending!=1);
        this->SPINBULLOCKbins_min.pop_back(); //delete the last element, which is END
        this->NSPINBULLOCKbins_power=this->SPINBULLOCKbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->SPINBULLOCKbins_min));
              }
    else
      if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'p' && line_in_file[1] == 'b'&& line_in_file[2] == '_'))
        {
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;   // read the first character, the name of the parameter
          line_string >> equality;    // check that second character is "="
          if (equality != "=") {
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
          }
          string par_value0;
          int ending=0;
          do{
        line_string >> par_value0;  //read value
        if(par_value0=="END")ending=1;
        this->SPINBULLOCKbins_max.push_back(atof(par_value0.c_str()));
          }while(ending!=1);
          this->SPINBULLOCKbins_max.pop_back(); //delete the last element, which is END
                  this->parameter_vectors.push_back(make_pair(par_name, this->SPINBULLOCKbins_max));

          if(this->NSPINBULLOCKbins_power!=this->SPINBULLOCKbins_max.size())
           {
          cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as SPINbins_min."<< endl;
          cerr << "Check parameter file" <<endl; exit(1);
            }
       }
  // ************************************************************************************************************
          // READ THE PARAMETER FOR TRACER MACH NUMBER:
		  else
                    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'e' && line_in_file[1] == '_' && line_in_file[2] == 'M'))
		      {
			stringstream line_string (line_in_file);
			line_string << line_in_file;
			line_string >> par_name;   // read the first character, the name of the parameter
			line_string >> equality;    // check that second character is "="
			if (equality != "=") {
			  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
			  cerr << "Using a default value for " << par_name << endl; exit(1);
			}
			string par_value0;
			int ending=0;
			do{
			  line_string >> par_value0;  //read value
			  if(par_value0=="END")ending=1;
			  this->MACHbins_min.push_back(atof(par_value0.c_str()));
			}while(ending!=1);
			this->MACHbins_min.pop_back(); //delete the last element, which is END
			this->NMACHbins_power=this->MACHbins_min.size();
                        this->parameter_vectors.push_back(make_pair(par_name, this->MACHbins_min));
                    }
		    else
                      if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'i' && line_in_file[1] == '_' && line_in_file[2] == 'M'))
			{
			  stringstream line_string (line_in_file);
			  line_string << line_in_file;
			  line_string >> par_name;   // read the first character, the name of the parameter
			  line_string >> equality;    // check that second character is "="
			  if (equality != "=") {
			    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
			    cerr << "Using a default value for " << par_name << endl; exit(1);
			  }
			  string par_value0;
			  int ending=0;
			  do{
			    line_string >> par_value0;  //read value
			    if(par_value0=="END")ending=1;
			    this->MACHbins_max.push_back(atof(par_value0.c_str()));
			  }while(ending!=1);
			  this->MACHbins_max.pop_back(); //delete the last element, which is END
                          this->parameter_vectors.push_back(make_pair(par_name, this->MACHbins_max));
                          if(this->NMACHbins_power!=this->MACHbins_max.size())
			    {
			      cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as MACHbins_min."<< endl;
			      cerr << "Check parameter file" <<endl; exit(1);
			    }
			}
// READ THE PARAMETER FOR TRACER DACH NUMBER:
        else
                  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'e' && line_in_file[1] == '_' && line_in_file[2] == 'D'))
            {
          stringstream line_string (line_in_file);
          line_string << line_in_file;
          line_string >> par_name;   // read the first character, the name of the parameter
          line_string >> equality;    // check that second character is "="
          if (equality != "=") {
            cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
            cerr << "Using a default value for " << par_name << endl; exit(1);
          }
          string par_value0;
          int ending=0;
          do{
            line_string >> par_value0;  //read value
            if(par_value0=="END")ending=1;
            this->DACHbins_min.push_back(atof(par_value0.c_str()));
          }while(ending!=1);
          this->DACHbins_min.pop_back(); //delete the last element, which is END
          this->NDACHbins_power=this->DACHbins_min.size();
                      this->parameter_vectors.push_back(make_pair(par_name, this->DACHbins_min));
                  }
          else
                    if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'i' && line_in_file[1] == '_' && line_in_file[2] == 'D'))
          {
            stringstream line_string (line_in_file);
            line_string << line_in_file;
            line_string >> par_name;   // read the first character, the name of the parameter
            line_string >> equality;    // check that second character is "="
            if (equality != "=") {
              cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
              cerr << "Using a default value for " << par_name << endl; exit(1);
            }
            string par_value0;
            int ending=0;
            do{
              line_string >> par_value0;  //read value
              if(par_value0=="END")ending=1;
              this->DACHbins_max.push_back(atof(par_value0.c_str()));
            }while(ending!=1);
            this->DACHbins_max.pop_back(); //delete the last element, which is END
                        this->parameter_vectors.push_back(make_pair(par_name, this->DACHbins_max));
                        if(this->NDACHbins_power!=this->DACHbins_max.size())
              {
                cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as DACHbins_min."<< endl;
                cerr << "Check parameter file" <<endl; exit(1);
              }
          }
        // READ THE PARAMETER FOR TRACER BIAS:
        else
          if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'e' && line_in_file[1] == 'b' && line_in_file[2] == '_' && line_in_file[3] == 'B' ))
            {
              stringstream line_string (line_in_file);
              line_string << line_in_file;
              line_string >> par_name;   // read the first character, the name of the parameter
              line_string >> equality;    // check that second character is "="
              if (equality != "=") {
                cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                cerr << "Using a default value for " << par_name << endl; exit(1);
              }
              string par_value0;
              int ending=0;
              do{
                line_string >> par_value0;  //read value
                if(par_value0=="END")ending=1;
                this->BIASbins_min.push_back(atof(par_value0.c_str()));
              }while(ending!=1);
              this->BIASbins_min.pop_back(); //delete the last element, which is END
              this->NBIASbins_power=this->BIASbins_min.size();
              this->parameter_vectors.push_back(make_pair(par_name, this->BIASbins_min));
          }
          else
            if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'i' && line_in_file[1] == 'b' && line_in_file[2] == '_' && line_in_file[3] == 'B'))
              {
                stringstream line_string (line_in_file);
                line_string << line_in_file;
                line_string >> par_name;   // read the first character, the name of the parameter
                line_string >> equality;    // check that second character is "="
                if (equality != "=") {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
                string par_value0;
                int ending=0;
                do{
                  line_string >> par_value0;  //read value
                  if(par_value0=="END")ending=1;
                  this->BIASbins_max.push_back(atof(par_value0.c_str()));
                }while(ending!=1);
                this->BIASbins_max.pop_back(); //delete the last element, which is END
                this->parameter_vectors.push_back(make_pair(par_name, this->BIASbins_max));
                if(this->NBIASbins_power!=this->BIASbins_max.size())
                  {
                    cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as MACHbins_min."<< endl;
                    cerr << "Check parameter file" <<endl; exit(1);
                  }
              }
          // READ THE PARAMETER FOR peak height:
          else
            if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'p' && line_in_file[1] == 'm' && line_in_file[2] == '_' && line_in_file[3] == 'P' ))
              {
                stringstream line_string (line_in_file);
                line_string << line_in_file;
                line_string >> par_name;   // read the first character, the name of the parameter
                line_string >> equality;    // check that second character is "="
                if (equality != "=") {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
                string par_value0;
                int ending=0;
                do{
                  line_string >> par_value0;  //read value
                  if(par_value0=="END")ending=1;
                  this->PHbins_min.push_back(atof(par_value0.c_str()));
                }while(ending!=1);
                this->PHbins_min.pop_back(); //delete the last element, which is END
                this->NPHbins_power=this->PHbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->PHbins_min));
            }
            else
              if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'p' && line_in_file[1] == 'x' && line_in_file[2] == '_' && line_in_file[3] == 'P'))
                {
                  stringstream line_string (line_in_file);
                  line_string << line_in_file;
                  line_string >> par_name;   // read the first character, the name of the parameter
                  line_string >> equality;    // check that second character is "="
                  if (equality != "=") {
                    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                    cerr << "Using a default value for " << par_name << endl; exit(1);
                  }
                  string par_value0;
                  int ending=0;
                  do{
                    line_string >> par_value0;  //read value
                    if(par_value0=="END")ending=1;
                    this->PHbins_max.push_back(atof(par_value0.c_str()));
                  }while(ending!=1);
                  this->PHbins_max.pop_back(); //delete the last element, which is END
                  this->parameter_vectors.push_back(make_pair(par_name, this->PHbins_max));
                  if(this->NPHbins_power!=this->PHbins_max.size())
                    {
                      cerr << "Error in parameters file. The container SPINbins_max has not the same dimension as MACHbins_min."<< endl;
                      cerr << "Check parameter file" <<endl; exit(1);
                    }
                }
          // READ THE PARAMETER FOR local clustering:
          else
            if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'e' && line_in_file[1] == 'l' && line_in_file[2] == '_' && line_in_file[3] == 'L' ))
              {
                stringstream line_string (line_in_file);
                line_string << line_in_file;
                line_string >> par_name;   // read the first character, the name of the parameter
                line_string >> equality;    // check that second character is "="
                if (equality != "=") {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
                string par_value0;
                int ending=0;
                do{
                  line_string >> par_value0;  //read value
                  if(par_value0=="END")ending=1;
                  this->LCbins_min.push_back(atof(par_value0.c_str()));
                }while(ending!=1);
                this->LCbins_min.pop_back(); //delete the last element, which is END
                this->NLCbins_power=this->LCbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->LCbins_min));
            }
            else
              if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'i' && line_in_file[1] == 'l' && line_in_file[2] == '_' && line_in_file[3] == 'L'))
                {
                  stringstream line_string (line_in_file);
                  line_string << line_in_file;
                  line_string >> par_name;   // read the first character, the name of the parameter
                  line_string >> equality;    // check that second character is "="
                  if (equality != "=") {
                    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                    cerr << "Using a default value for " << par_name << endl; exit(1);
                  }
                  string par_value0;
                  int ending=0;
                  do{
                    line_string >> par_value0;  //read value
                    if(par_value0=="END")ending=1;
                    this->LCbins_max.push_back(atof(par_value0.c_str()));
                  }while(ending!=1);
                  this->LCbins_max.pop_back(); //delete the last element, which is END
                  this->parameter_vectors.push_back(make_pair(par_name, this->LCbins_max));
                  if(this->NLCbins_power!=this->LCbins_max.size())
                    {
                      cerr << "Error in parameters file. The container LCbins_max has not the same dimension as LCbins_min."<< endl;
                      cerr << "Check parameter file" <<endl; exit(1);
                    }
                }

          // READ THE PARAMETER FOR halo tidal anisotropy
          else
            if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'e' && line_in_file[1] == 't' && line_in_file[2] == '_' && line_in_file[3] == 'T' ))
              {
                stringstream line_string (line_in_file);
                line_string << line_in_file;
                line_string >> par_name;   // read the first character, the name of the parameter
                line_string >> equality;    // check that second character is "="
                if (equality != "=") {
                  cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                  cerr << "Using a default value for " << par_name << endl; exit(1);
                }
                string par_value0;
                int ending=0;
                do{
                  line_string >> par_value0;  //read value
                  if(par_value0=="END")ending=1;
                  this->TAbins_min.push_back(atof(par_value0.c_str()));
                }while(ending!=1);
                this->TAbins_min.pop_back(); //delete the last element, which is END
                this->NTAbins_power=this->TAbins_min.size();
                this->parameter_vectors.push_back(make_pair(par_name, this->TAbins_min));
            }
            else
              if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'i' && line_in_file[1] == 't' && line_in_file[2] == '_' && line_in_file[3] == 'T'))
                {
                  stringstream line_string (line_in_file);
                  line_string << line_in_file;
                  line_string >> par_name;   // read the first character, the name of the parameter
                  line_string >> equality;    // check that second character is "="
                  if (equality != "=") {
                    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
                    cerr << "Using a default value for " << par_name << endl; exit(1);
                  }
                  string par_value0;
                  int ending=0;
                  do{
                     line_string >> par_value0;  //read value
                    if(par_value0=="END")ending=1;
                    this->TAbins_max.push_back(atof(par_value0.c_str()));
                  }while(ending!=1);
                  this->TAbins_max.pop_back(); //delete the last element, which is END
                  this->parameter_vectors.push_back(make_pair(par_name, this->TAbins_max));
                  if(this->NTAbins_power!=this->TAbins_max.size())
                    {
                      cerr << "Error in parameters file. The container TAbins_max has not the same dimension as TAbins_min."<< endl;
                      cerr << "Check parameter file" <<endl; exit(1);
                    }
                }
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
      if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'u' && line_in_file[1] == '_' && line_in_file[2] == 'P'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
		cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
		cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
		line_string >> par_value0;  //read value
		if(par_value0=="END")ending=1;
		this->list_Props_Threshold_MultiLevels.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_Props_Threshold_MultiLevels.pop_back(); //delete the last element, which is END
          this->parameter_vectors.push_back(make_pair(par_name, this->list_Props_Threshold_MultiLevels));
       }
      else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'q' && line_in_file[1] == '_'&& line_in_file[2] == 'N'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
	      if (equality != "=") {
            cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
            cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
            line_string >> par_value0;  //read value
            if(par_value0=="END")ending=1;
            this->list_Nft_MultiLevels.push_back(atoi(par_value0.c_str()));
            }while(ending!=1);
           this->list_Nft_MultiLevels.pop_back(); //delete the last element, which is END
       //    this->parameter_vectors.push_back(make_pair(par_name, this->list_Nft_MultiLevels));
        }
      else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 's' && line_in_file[1] == '_' && line_in_file[2] == 'P'))
	    {
	      stringstream line_string (line_in_file);
	      line_string << line_in_file;
	      line_string >> par_name;   // read the first character, the name of the parameter
	      line_string >> equality;    // check that second character is "="
          if (equality != "=")
          {
            cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
            cerr << "Using a default value for " << par_name << endl; exit(1);
	      }
	      string par_value0;
	      int ending=0;
	      do{
            line_string >> par_value0;  //read value
            if(par_value0=="END")ending=1;
            this->list_Props_Tolerance_MultiLevels.push_back(atof(par_value0.c_str()));
	      }while(ending!=1);
	      this->list_Props_Tolerance_MultiLevels.pop_back(); //delete the last element, which is END
              this->parameter_vectors.push_back(make_pair(par_name, this->list_Props_Tolerance_MultiLevels));
            }
	  	  // These arrays are inizialized here and updated in teh Bam.pp file 
	  this->Number_of_MultiLevels = this->list_Nft_MultiLevels.size();
#else
	  this->Number_of_MultiLevels = 1;
#endif
	  this->list_Ntracers_MultiLevels.resize(this->Number_of_MultiLevels,0);
	  
	  // global and I/O parameters for Power spectrum
	  if (par_name == "Statistics")
	    {
          this->statistics = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->statistics)); 
	    }
	  else if (par_name == "Catalogue_file")
	    {
          this->file_catalogue = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_catalogue)); 
	    }
	  else if (par_name == "Catalogue_file_new_ref")
	    {
          this->file_catalogue_new_ref = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_catalogue_new_ref)); 
	    }
	  else if (par_name == "Input_type")
	    {
          this->input_type = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->input_type)); 
	    }
	  else if (par_name == "Input_type_two")
	    {
          this->input_type_two = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->input_type_two)); 
	    }
	  else if (par_name == "delta_grid_file")
	    {
          this->delta_grid_file = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->delta_grid_file)); 
	    }
	  else if (par_name == "delta_grid_file2")
	    {
          this->delta_grid_file2 = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->delta_grid_file2)); 
	    }
	  else if (par_name == "delta_grid_file3")
	    {
          this->delta_grid_file3 = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->delta_grid_file3)); 
	    }
	  else if (par_name == "delta_grid_file4")
	    {
          this->delta_grid_file4 = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->delta_grid_file4)); 
	    }
	  else if (par_name == "ngal_delta")
	    {
	      this->ngal_delta = atof(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->ngal_delta)); 
	    }
	  else if (par_name == "measure_cross")
	    {
	      if(par_value == "true")this->measure_cross = true;
	      else if(par_value == "false")this->measure_cross = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->measure_cross)); 
	    }
	  else if (par_name == "get_window_matrix")
	    {
	      if(par_value == "true")this->get_window_matrix = true;
	      else if(par_value == "false")this->get_window_matrix = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->get_window_matrix)); 
	    }
	  else if (par_name == "redshift_space_coords_g")
	    {
	      if(par_value == "true")this->redshift_space_coords_g  = true;
	      else if(par_value == "false")this->redshift_space_coords_g  = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->redshift_space_coords_g)); 
	    }
          else if (par_name == "xmin")
            {
              this->xmin = atof(par_value.c_str());
              this->parameter_number.push_back(make_pair(par_name, this->xmin));
            }
          else if (par_name == "ymin")
            {
              this->ymin = atof(par_value.c_str());
              this->parameter_number.push_back(make_pair(par_name, this->ymin));
            }
          else if (par_name == "zmin")
            {
              this->zmin = atof(par_value.c_str());
              this->parameter_number.push_back(make_pair(par_name, this->zmin));
            }
          else if (par_name == "Nbins_cf")
	    {
	      this->Nbins_cf = static_cast<ULONG>(atoi(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Nbins_cf)); 		      
	    }
	  else if (par_name == "rmin_cf")
	    {
	      this->rmin_cf = atof(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->rmin_cf)); 		      
	    }
	  else if (par_name == "rmax_cf")
	    {
	      this->rmax_cf = atof(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->rmax_cf)); 		      
	    }
	  else if (par_name == "mark")
	    {
          this->mark = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->mark)); 
	    }
	  else if (par_name == "rbin_type")
	    {
          this->rbin_type = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->rbin_type)); 
	    }
	  else if (par_name == "measure_cross_from_1")
	    {
	      this->measure_cross_from_1 = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->measure_cross_from_1)); 
	    }
	  else if (par_name == "measure_cross_from_2")
	    {
	      this->measure_cross_from_2 = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->measure_cross_from_2)); 
	    }
	  else if (par_name == "sys_of_coord_g")
	    {
	      this->sys_of_coord_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->sys_of_coord_g)); 
	    }
          else if (par_name == "i_coord1_g")
	    {
	      this->i_coord1_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord1_g)); 
	    }
          else if (par_name == "i_coord2_g")
	    {
	      this->i_coord2_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord2_g)); 
	    }
          else if (par_name == "i_coord3_g")
	    {
	      this->i_coord3_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord3_g)); 
	    }
          else if (par_name == "i_mean_density_g")
	    {
	      this->i_mean_density_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mean_density_g)); 
	    }
          else if (par_name == "i_weight1_g")
	    {
	      this->i_weight1_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight1_g)); 
	    }
          else if (par_name == "i_weight2_g")
	    {
	      this->i_weight2_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight2_g)); 
	    }
          else if (par_name == "i_weight3_g")
	    {
	      this->i_weight3_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight3_g)); 
	    }
          else if (par_name == "i_weight4_g")
	    {
	      this->i_weight4_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight4_g)); 
	    }
	  else if (par_name == "vel_units_g")
	    {
          this->vel_units_g = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->vel_units_g)); 
	    }
          else if (par_name == "i_v1_g")
	    {
	      this->i_v1_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v1_g)); 
	    }
          else if (par_name == "i_v2_g")
	    {
	      this->i_v2_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v2_g)); 	      
	    }
          else if (par_name == "i_v3_g")
	    {
	      this->i_v3_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v3_g)); 	      
	    }
          else if (par_name == "i_mass_g")
	    {
	      this->i_mass_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mass_g)); 	      
	    }
          else if (par_name == "i_vmax_g")
	    {
	      this->i_vmax_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_vmax_g)); 	      
	    }
	  else if (par_name == "i_vrms_g")
	    {
	      this->i_vrms_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_vrms_g)); 	      	      
	    }
	  else if (par_name == "i_sf_g")
	    {
	      this->i_sf_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_sf_g)); 	      	      
	    }
          else if (par_name == "i_rs_g")
	    {
	      this->i_rs_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_rs_g)); 	      	      
	    }
      else if (par_name == "i_rvir_g")
        {
           this->i_rvir_g = atoi(par_value.c_str());
           this->parameter_number.push_back(make_pair(par_name, this->i_rvir_g));
        }
          else if (par_name == "i_virial_g")
	    {
	      this->i_virial_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_virial_g)); 	      
	    }
          else if (par_name == "i_spin_g")
	    {
	      this->i_spin_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_spin_g)); 	      
	    }
      else if (par_name == "i_spin_bullock_g")
        {
          this->i_spin_bullock_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_spin_bullock_g));
        }
      else if (par_name == "i_b_to_a_g")
	    {
	      this->i_b_to_a_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_b_to_a_g)); 	      
	    }
	  else if (par_name == "i_c_to_a_g")
	    {
	      this->i_c_to_a_g = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_c_to_a_g)); 	      
	    }

      else if (par_name == "i_color_g")
        {
          this->i_color_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_color_g));
        }
      else if (par_name == "i_stellar_mass_g")
        {
          this->i_stellar_mass_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_stellar_mass_g));
        }
      else if (par_name == "i_app_mag_g")
        {
          this->i_app_mag_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_app_mag_g));
        }
      else if (par_name == "i_abs_mag_g")
        {
          this->i_abs_mag_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_abs_mag_g));
        }

      else if (par_name == "Get_Mstellar_function")
        {
          if(par_value == "true")this->Get_Mstellar_function = true;
          else if(par_value == "false")this->Get_Mstellar_function = false;
          this->parameter_boolean.push_back(make_pair(par_name, this->Get_Mstellar_function));
        }
      else if (par_name == "Get_Luminosity_function")
        {
          if(par_value == "true")this->Get_Luminosity_function= true;
          else if(par_value == "false")this->Get_Luminosity_function = false;
          this->parameter_boolean.push_back(make_pair(par_name, this->Get_Luminosity_function));
        }
      else if (par_name == "Get_Color_function")
        {
          if(par_value == "true")this->Get_Color_function= true;
          else if(par_value == "false")this->Get_Color_function = false;
          this->parameter_boolean.push_back(make_pair(par_name, this->Get_Color_function));
        }
      else if (par_name == "Get_Color_Mag_plane")
        {
          if(par_value == "true")this->Get_Color_Mag_plane= true;
          else if(par_value == "false")this->Get_Color_Mag_plane = false;
          this->parameter_boolean.push_back(make_pair(par_name, this->Get_Color_Mag_plane));
        }
      else if (par_name == "Get_Random_Catalog")
        {
          if(par_value == "true")this->Get_Random_Catalog= true;
          else if(par_value == "false")this->Get_Random_Catalog= false;
          this->parameter_boolean.push_back(make_pair(par_name, this->Get_Random_Catalog));
        }
      else if (par_name == "Nbins_color")
        {
          this->Nbins_color = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Nbins_color));
        }
      else if (par_name == "Nbins_Mstellar")
        {
          this->Nbins_Mstellar = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Nbins_Mstellar));
        }
      else if (par_name == "Mstellar_min")
        {
          this->Mstellar_min = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Mstellar_min));
        }
      else if (par_name == "Mstellar_max")
        {
          this->Mstellar_max = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Mstellar_max));
        }
      else if (par_name == "Color_min")
        {
          this->Color_min = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Color_min));
        }
      else if (par_name == "Color_max")
        {
          this->Color_max = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Color_max));
        }
      else if (par_name == "DEC_max")
        {
          this->DEC_max = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->DEC_max));
        }
      else if (par_name == "DEC_min")
        {
          this->DEC_min = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->DEC_min));
        }
      else if (par_name == "RA_max")
        {
          this->RA_max = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->RA_max));
        }
      else if (par_name == "RA_min")
        {
          this->RA_min = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->RA_min));
        }
      else if (par_name == "mKmin")
        {
          this->mKmin = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->mKmin));
        }
      else if (par_name == "mKmax")
        {
          this->mKmax = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->mKmax));
        }
      else if (par_name == "sys_of_coord_dm")
	    {
	      this->sys_of_coord_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->sys_of_coord_dm)); 	      
	    }
          else if (par_name == "i_coord1_dm")
	    {
	      this->i_coord1_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord1_dm)); 	      
	    }
          else if (par_name == "i_coord2_dm")
	    {
	      this->i_coord2_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord2_dm)); 	      
	    }
          else if (par_name == "i_coord3_dm")
	    {
	      this->i_coord3_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord3_dm)); 	      
	    }
          else if (par_name == "i_v1_dm")
	    {
	      this->i_v1_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v1_dm)); 	      
	    }
          else if (par_name == "i_v2_dm")
	    {
	      this->i_v2_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v2_dm)); 	      
	    }
          else if (par_name == "i_v3_dm")
	    {
	      this->i_v3_dm = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_v3_dm)); 	      
	    }
	  
	  else if (par_name == "weight_vel_with_mass")
	    {
	      if (par_value=="true")this->weight_vel_with_mass = true;
	      else this->weight_vel_with_mass = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->weight_vel_with_mass)); 	      
	    }
	  else if (par_name == "weight_with_mass")
	    {
	      if (par_value=="true")this->weight_with_mass = true;
	      else this->weight_with_mass = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->weight_vel_with_mass)); 	      
	    }
	  else if (par_name == "use_weight1_g")
	    {
	      if(par_value == "true")this->use_weight1_g = true;
	      else if(par_value == "false")this->use_weight1_g = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight1_g)); 	      
	    }
	  else if (par_name == "use_weight2_g")
	    {
	      if(par_value == "true")this->use_weight2_g = true;
	      else if(par_value == "false")this->use_weight2_g = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight2_g)); 	      
	      
	    }
	  else if (par_name == "use_weight3_g")
	    {
	      if(par_value == "true")this->use_weight3_g = true;
	      else if(par_value == "false")this->use_weight3_g = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight3_g)); 	      
	    }
	  else if (par_name == "use_weight4_g")
	    {
	      if(par_value == "true")this->use_weight4_g = true;
	      else if(par_value == "false")this->use_weight4_g = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight4_g)); 	      
	    }
	  else if (par_name == "new_Lbox")
	    {
	      if (par_value == "true"){this->new_Lbox = true;}
	      else if (par_value == "false"){this->new_Lbox = false;}
	      this->parameter_number.push_back(make_pair(par_name, this->new_Lbox)); 	      
	    }
	  else if (par_name == "n_catalogues")
	    {
	      this->n_catalogues = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->n_catalogues)); 	      
	    }
	  else if (par_name == "angles_units_g")
	    {
	      if (par_value=="D") this->angles_units_g = "D";
	      else if (par_value=="R") this->angles_units_g = "R";
	      else {cout<<"Unknown angle units"<<std::endl; exit(1) ;}
	      this->parameter_string.push_back(make_pair(par_name, this->angles_units_g)); 	      
	    }
	  else if (par_name == "use_random_catalog")
	    {
	      if (par_value=="false") this->use_random_catalog = false;
	      else if (par_value=="true") this->use_random_catalog = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_random_catalog)); 	      
	    }
	  else if (par_name == "use_random_catalog_cl")
	    {
	      if (par_value=="false") this->use_random_catalog_cl = false;
	      else if (par_value=="true") this->use_random_catalog_cl = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_random_catalog_cl)); 	      
	    }
	  else if (par_name == "sys_of_coord_r")
	    {
	      this->sys_of_coord_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->sys_of_coord_r)); 	      
	    }
          else if (par_name == "i_coord1_r")
	    {
	      this->i_coord1_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord1_r)); 	      
	    }
          else if (par_name == "i_coord2_r")
	    {
	      this->i_coord2_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord2_r)); 	      	      
	    }
          else if (par_name == "i_coord3_r")
	    {
	      this->i_coord3_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_coord3_r)); 	      
	    }
          else if (par_name == "i_mean_density_r")
	    {
	      this->i_mean_density_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mean_density_r)); 	      	      
	    }
          else if (par_name == "i_weight1_r")
	    {
	      this->i_weight1_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight1_r)); 	      	      
	    }
          else if (par_name == "i_weight2_r")
	    {
	      this->i_weight2_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight2_r)); 	      	      
	    }
          else if (par_name == "i_weight3_r")
	    {
	      this->i_weight3_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight3_r)); 	      	      
	    }
          else if (par_name == "i_weight4_r")
	    {
	      this->i_weight4_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_weight4_r)); 	      	      
	    }
      else if (par_name == "i_color_r")
        {
          this->i_color_r = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_color_g));
        }
      else if (par_name == "i_stellar_mass_g")
        {
          this->i_stellar_mass_r = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_stellar_mass_r));
        }
      else if (par_name == "i_app_mag_r")
        {
          this->i_app_mag_g = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_app_mag_r));
        }
      else if (par_name == "i_abs_mag_r")
        {
          this->i_abs_mag_r = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->i_abs_mag_r));
        }
          else if (par_name == "i_mass_r")
	    {
	      this->i_mass_r = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mass_r)); 	      	      
	    }
	  else if (par_name == "get_distribution_min_separations")
	    {
	      if (par_value=="true")this->get_distribution_min_separations = true;
	      else this->get_distribution_min_separations = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->get_distribution_min_separations)); 	      	      
	    }
	  else if (par_name == "i_mask_pixel")
	    {
	      this->i_mask_pixel= atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mask_pixel)); 	      	      
	    }
	  else if (par_name == "i_mask_alpha")
	    {
	      this->i_mask_alpha= atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mask_alpha)); 	      	      
	    }
	  else if (par_name == "i_mask_delta")
	    {
	      this->i_mask_delta= atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mask_delta)); 	      	      
	    }
	  else if (par_name == "i_mask_flag")
	    {
	      this->i_mask_flag= atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->i_mask_flag)); 	      	      
	    }
	  else if (par_name == "use_weight1_r")
	    {
	      if(par_value == "true")this->use_weight1_r = true;
	      else if (par_value == "false")this->use_weight1_r = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight1_r)); 	      	      
	    }
	  else if (par_name == "use_weight2_r")
	    {
	      if(par_value == "true")this->use_weight2_r = true;
	      else if (par_value == "false")use_weight2_r = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight2_r)); 	      	      
	    }
	  else if (par_name == "use_weight3_r")
	    {
	      if(par_value == "true")use_weight3_r = true;
	      else if (par_value == "false")use_weight3_r = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight3_r)); 	      	      
	    }
	  else if (par_name == "use_weight4_r")
	    {
	      if(par_value == "true")use_weight4_r = true;
	      else if (par_value == "false")use_weight4_r = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_weight4_r)); 	      	      
	    }
	  else if (par_name == "angles_units_r")
	    {
	      if (par_value == "D") angles_units_r = "D";
	      else if (par_value == "R") angles_units_r = "R";
	      else { cout<<"Unknown angle units"<<std::endl; exit(1); }
	      this->parameter_string.push_back(make_pair(par_name, this->angles_units_r)); 	      	      
	    }
	  else if (par_name == "Random_file")
	    {
          file_random = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_random)); 	      	      
	    }
	  else if (par_name == "Name_survey")
	    {
          Name_survey = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->Name_survey)); 	      	      	      
	    }
          else if (par_name == "new_los")
	    {
	      if (par_value == "false") new_los = false;
	      else if (par_value == "true") new_los = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->new_los)); 	      	      	      
	    }
	  else if (par_name == "Nft")
	    {
	      this->Nft = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Nft)); 	      	      	      
	    }
    else if (par_name == "Nft_JK")
      {
        this->Nft_JK = atoi(par_value.c_str());
        this->parameter_number.push_back(make_pair(par_name, this->Nft));                         
      }
	  else if (par_name == "Nft_low")
	    {
	      this->Nft_low = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Nft_low)); 	      	      	      
	    }
	  else if (par_name == "Nft_HR")
	    {
	      this->Nft_HR = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Nft_HR)); 	      	      	      
	    }
	  else if (par_name == "Nft_random_collapse")
	    {
	      this->Nft_random_collapse = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Nft_random_collapse));
	    }
	  else if (par_name == "Lbox")
	    {
	      this->Lbox = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Lbox));
	    }
	  else if (par_name == "Distance_fraction")
	    {
	      this->Distance_fraction = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Distance_fraction));
	    }
	  else if (par_name == "velbias_random")
	    {
	      this->velbias_random = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->velbias_random));
	    }
	  else if (par_name == "velbias_dm")
	    {
	      this->velbias_dm = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->velbias_dm));
	    }
	  else if (par_name == "Lbox_low")
	    {
	      this->Lbox_low = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Lbox_low));
	    }
	  else if (par_name == "mass_assignment_scheme")
	    {
	      if (par_value=="NGP") this->mass_assignment_scheme = "NGP";
	      else if (par_value=="CIC") this->mass_assignment_scheme = "CIC";
	      else if (par_value=="TSC") this->mass_assignment_scheme = "TSC";
	      else if (par_value=="PCS") this->mass_assignment_scheme = "PCS";
	      else{cout<<"Unknown mass assigment scheme"<<std::endl;exit(1);}
	      this->parameter_string.push_back(make_pair(par_name, this->mass_assignment_scheme));
	    }
	  else if (par_name == "MAS_correction")
	    {
	      if (par_value=="true") MAS_correction = true;
	      else if (par_value=="false") MAS_correction = false;
	      this->parameter_boolean.push_back(make_pair(par_name, this->MAS_correction));
	    }
	  else if (par_name == "type_of_binning")
	    {
	      if (par_value=="linear") type_of_binning = "linear";
	      else if (par_value=="log") type_of_binning = "log";
	      else{cout<<"Unknown type of binning"<<std::endl; exit(1);}
	      this->parameter_string.push_back(make_pair(par_name, this->type_of_binning));
	    }
	  else if (par_name == "DeltaKmin")
	    {
	      this->DeltaKmin = atof(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->DeltaKmin));
	    }
	  else if (par_name == "N_log_bins")
	    {
	      this->N_log_bins = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->N_log_bins));
	    }
	  else if (par_name == "ndel_data")
	    {
	      this->ndel_data = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->ndel_data));
	    }
	  else if (par_name == "ndel_window")
	    {
	      this->ndel_window = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->ndel_window));
	    }
	  else if (par_name == "N_mu_bins")
	    {
	      this->N_mu_bins = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->N_mu_bins));
	    }
	  else if (par_name == "FKP_weight")
	    {
	      if (par_value=="false") this->FKP_weight = false;
	      else if (par_value=="true") this->FKP_weight = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->FKP_weight));
	    }
	  else if (par_name == "Pest")
	    {
	      this->Pest = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Pest));
	    }
	  else if (par_name == "SN_correction")
	    {
	      if (par_value=="false") this->SN_correction = false;
	      else if (par_value=="true") this->SN_correction = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->SN_correction));
	    }
	  else if (par_name == "FKP_error_bars")
	    {
	      if (par_value=="false") this->FKP_error_bars = false;
	      else if (par_value=="true") this->FKP_error_bars = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->FKP_error_bars));
	    }
	  else if (par_name == "FKP_error_bars_exact")
	    {
	      if (par_value=="false")this->FKP_error_bars_exact = false;
	      else if (par_value=="true") this->FKP_error_bars_exact = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->FKP_error_bars_exact));
	    }
	  else if (par_name == "nbar_tabulated")
	    {
	      if (par_value=="false") nbar_tabulated = false;
	      else if (par_value=="true") nbar_tabulated = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->nbar_tabulated));
	    }
      else if (par_name == "use_file_nbar")
        {
          if (par_value=="false") use_file_nbar = false;
          else if (par_value=="true") use_file_nbar = true;
          this->parameter_boolean.push_back(make_pair(par_name, this->use_file_nbar));
        }
      else if (par_name == "constant_depth")
	    {
	      if (par_value=="false") constant_depth = false;
	      else if (par_value=="true") constant_depth = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->constant_depth));
	    }
      else if (par_name == "Nbins_redshift")
	    {
          Nbins_redshift = atoi(par_value.c_str());
          this->parameter_number.push_back(make_pair(par_name, this->Nbins_redshift));
	    }
	  else if (par_name == "redshift_min_sample")
	    {
	      redshift_min_sample = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->redshift_min_sample));
	    }
	  else if (par_name == "redshift_max_sample")
	    {
	      redshift_max_sample = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->redshift_max_sample));
	    }
	  else if (par_name == "N_dndz_bins")
	    {
	      N_dndz_bins = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->N_dndz_bins));
	    }
	  else if (par_name == "new_N_dndz_bins")
	    {
	      new_N_dndz_bins = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->new_N_dndz_bins));
	    }
	  else if (par_name == "area_survey")
	    {
	      area_survey = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->area_survey));
	    }
      else if (par_name == "Number_of_random_files")
        {
          this->Number_of_random_files = static_cast<real_prec>(atof(par_value.c_str()));
          this->parameter_number.push_back(make_pair(par_name, this->Number_of_random_files));
        }
      else if (par_name == "Healpix_resolution")
	    {
	      this->Healpix_resolution = atoi(par_value.c_str());
	      this->parameter_number.push_back(make_pair(par_name, this->Healpix_resolution));
	    }
      else if (par_name == "file_nbar")
        {
          this->file_nbar = par_value.c_str();
          this->parameter_string.push_back(make_pair(par_name, this->file_nbar));
        }
      else if (par_name == "file_power")
	    {
          file_power = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_power));
	    }
	  else if (par_name == "file_power_log")
	    {
          this->file_power_log = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_power_log));
	    }
	  else if (par_name == "file_window")
	    {
          this->file_window = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_window));
	    }
	  else if (par_name == "file_dndz")
	    {
          this->file_dndz = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_dndz));
	    }
	  else if (par_name == "file_power2d")
	    {
          file_power2d = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_power2d));
	    }
	  else if (par_name == "file_power2d_mk")
	    {
          this->file_power2d_mk = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_power2d_mk));
	    }
	  else if (par_name == "file_bispectrum")
	    {
          this->file_bispectrum = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->file_bispectrum));
	    }
	  else if (par_name == "kmax_y_ds")
	    {
	      kmax_y_ds = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->kmax_y_ds));
	    }
	  //Parameters for the bispectrum
	  else if (par_name == "kmax_bk")
	    {
	      kmax_bk = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->kmax_bk));
	    }
	  else if (par_name == "kmin_bk")
	    {
	      kmin_bk = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->kmin_bk));
	    }
	  else if (par_name == "use_fundamental_mode_as_kmin_bk")
	    {
	      if (par_value=="false") use_fundamental_mode_as_kmin_bk = false;
	      else if (par_value=="true") use_fundamental_mode_as_kmin_bk = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_fundamental_mode_as_kmin_bk));
	    }
#ifndef _USE_COSMO_PARS_
	  else if (par_name == "om_matter")
	    {
	      this->om_matter = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->om_matter));
	    }
	  else if (par_name == "om_radiation")
	    {
	      this->om_radiation = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->om_radiation));
	    }
	  else if (par_name == "om_baryons")
	    {
	      this->om_baryons = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->om_baryons));
	    }
	  else if (par_name == "om_vac")
	    {
	      this->om_vac = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->om_vac));

	    }
	  else if (par_name == "f_baryon")
	    {
	      this->f_baryon = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->f_baryon));

	    }
	  else if (par_name == "om_k")
	    {
	      this->om_k = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->om_k));
	    }
	  else if (par_name == "Hubble")
	    {
	      this->Hubble = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Hubble));
	    }
	  else if (par_name == "hubble")
	    {
	      this->hubble = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->hubble));
	    }
	  else if (par_name == "spectral_index")
	    {
	      this->spectral_index = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->spectral_index));
	    }
	  else if (par_name == "wde_eos")
	    {
	      this->wde_eos = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->wde_eos));

	    }
	  else if (par_name == "N_eff")
	    {
	      this->N_eff = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->N_eff));
	    }
	  else if (par_name == "sigma8")
	    {
	      this->sigma8 = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->sigma8));
	    }
	  else if (par_name == "Tcmb")
	    {
	      this->Tcmb = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Tcmb));
	    }
	  else if (par_name == "RR")
	    {
	      this->RR = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->RR));
	    }
      else if (par_name == "Delta_SO")
      {
        this->Delta_SO = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Delta_SO));
      }
#endif
 
          //here now read the aprams from LPT
	  else if (par_name == "inputmode")
	    {
	      this->inputmode = static_cast<int>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->inputmode));

	    }
	  else if (par_name == "seed")
	    {
	      this->seed = static_cast<int>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->seed));
	    }
	  else if (par_name == "seed_ref")
	    {
	      this->seed_ref = static_cast<int>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->seed_ref));
	    }
	  else if (par_name == "runsim")
	    {
	      if (par_value=="false") this->runsim = false;
	      else if (par_value=="true") this->runsim = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->runsim));
	    }
	  else if (par_name == "runv")
	    {
	      if (par_value=="false") this->runv = false;
	      else if (par_value=="true") this->runv = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->runv));
	    }
          else if (par_name == "ic_alias_corrected")
	    {
	      if (par_value=="false") this->ic_alias_corrected = false;
	      else if (par_value=="true") this->ic_alias_corrected = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->ic_alias_corrected));
	    }
	  else if (par_name == "Initial_Redshift_TH_power_file")
	    {
	      this->Initial_Redshift_TH_power_file = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Initial_Redshift_TH_power_file));
	    }
	  else if (par_name == "Initial_Redshift_DELTA")
	    {
	      this->Initial_Redshift_DELTA = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Initial_Redshift_DELTA));
	    }
	  else if (par_name == "Initial_Redshift_SIM")
	    {
	      this->Initial_Redshift_SIM = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Initial_Redshift_SIM));
	    }
	  else if (par_name == "ic_power_file")
	    {
          this->ic_power_file = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->ic_power_file));
	    }
	  else if (par_name == "ic_WN_file")
	    {
          this->ic_WN_file = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->ic_WN_file));
	    }
	  else if (par_name == "ic_file")
	    {
          this->ic_file = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->ic_file));

	    }
	  else if (par_name == "ic_input_type")
	    {
          this->ic_input_type = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->ic_input_type));
	      
	    }
	  else if (par_name == "ic_WN_dir")
	    {
          this->ic_WN_dir = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->ic_WN_dir));
	      
	    }
	  else if (par_name == "dir")
	    {
          this->dir = par_value.c_str();
	      this->parameter_string.push_back(make_pair(par_name, this->dir));
	    }
	  else if (par_name == "use_ic_file")
	    {
	      if (par_value=="false") this->use_ic_file = false;
	      else if (par_value=="true") this->use_ic_file = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_ic_file));
	    }
	  
	  else if (par_name == "readPS")
	    {
	      if (par_value=="false") this->readPS = false;
	      else if (par_value=="true") this->readPS = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->readPS));
	    }
	  else if (par_name == "use_vel_kernel")
	    {
	      if (par_value=="false") this->use_vel_kernel = false;
	      else if (par_value=="true") this->use_vel_kernel = true;
	      this->parameter_boolean.push_back(make_pair(par_name, this->use_vel_kernel));
	    }
	  else if (par_name == "vkernel_exponent")
	    {
	      this->vkernel_exponent = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_boolean.push_back(make_pair(par_name, this->vkernel_exponent));
	    }
	  else if (par_name == "Structure_Formation_Model")
	    {
	      this->Structure_Formation_Model = (int)atof(par_value.c_str());
	      this->parameter_boolean.push_back(make_pair(par_name, this->Structure_Formation_Model));
	    }
	  else if (par_name == "masskernel")
	    {
	      this->masskernel = (int)atof(par_value.c_str());
	      this->parameter_boolean.push_back(make_pair(par_name, this->masskernel));
	    }
	  else if (par_name == "masskernel_vel")
	    {
	      this->masskernel_vel = (int)atof(par_value.c_str());
	      this->parameter_boolean.push_back(make_pair(par_name, this->masskernel_vel));
	    }
	  else if (par_name == "slength")
	    {
	      this->slength = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_boolean.push_back(make_pair(par_name, this->slength));
	    }
	  else if (par_name == "slengthv")
	    {
	      this->slengthv = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_boolean.push_back(make_pair(par_name, this->slengthv));

	    }
	  else if (par_name == "velbias")
	    {
	      this->velbias = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_boolean.push_back(make_pair(par_name, this->velbias));
	    }
	  else if (par_name == "vslength")
	    {
	      this->vslength = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_boolean.push_back(make_pair(par_name, this->vslength));
	    }
      /// ------Parameter for cosmolib----------------------------------------//
      else if (par_name == "kmin_integration") 
        this->kmin_integration = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "kmax_integration") 
        this->kmax_integration = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "kstar") 
        this->kstar =static_cast<real_prec>(atof( par_value.c_str()));
      else if (par_name == "Amc") 
        this->Amc =static_cast<real_prec>(atof( par_value.c_str()));
      else if (par_name == "GAL_BIAS") 
        this->GAL_BIAS = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "fixed_redshift")
      {
        if(par_value=="true")this->fixed_redshift=true;
        else if(par_value=="false")this->fixed_redshift=false;
      }
            else if (par_name == "use_wiggles")
      {
        if(par_value=="true")this->use_wiggles=true;
        else if(par_value=="false")this->use_wiggles =false;
      }
      else if (par_name == "Get_SO_from_BN"){
        if(par_value=="true")this->Get_SO_from_BN=true;
        else if(par_value=="false")this->Get_SO_from_BN=false;
      }
      else if (par_name == "file_power_th") 
        this->file_power_th =par_value;
      else if (par_name == "redshift_min") 
        redshift_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "redshift_max") 
        this->redshift_max = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "nbins_redshift") this->nbins_redshift = static_cast<int>(atof(par_value.c_str()));
      else if (par_name == "A_gas") this->A_gas = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "B_gas") this->B_gas = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mstar") this->mstar = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "sigma_red") this->sigma_red = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "sigma_ln") this->sigma_ln = static_cast<real_prec>(atof(par_value.c_str()));
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
      else if (par_name == "npoints_mf")npoints_mf = atoi(par_value.c_str());
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
      else if (par_name == "kmin_ps")kmin_ps = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "kmax_ps")kmax_ps = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "npoints_ps")npoints_ps = atoi(par_value.c_str());
      else if (par_name == "linear_matter_ps_output_file")linear_matter_ps_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_halo_fit_output_file")non_linear_matter_ps_halo_fit_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_ps_pt_output_file")non_linear_matter_ps_pt_output_file = par_value.c_str();
      else if (par_name == "scale_cf")scale_cf = par_value.c_str();
      else if (par_name == "Output_directory")this->Output_directory = par_value.c_str();
      else if (par_name == "rmin_cf")rmin_cf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "rmax_cf")rmax_cf = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "npoints_cf")npoints_cf = atoi(par_value.c_str());
      else if (par_name == "linear_matter_cf_output_file")linear_matter_cf_output_file = par_value.c_str();
      else if (par_name == "non_linear_matter_cf_halo_fit_output_file")non_linear_matter_cf_halo_fit_output_file = par_value.c_str();
      else if (par_name == "compute_density_profile")
      {
        if(par_value=="true")compute_density_profile=true;
        else if(par_value=="false")compute_density_profile=false;
      }
      else if (par_name == "rmin_dp")rmin_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "rmax_dp")rmax_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "scale_dp_r")scale_dp_r = par_value.c_str();
      else if (par_name == "npoints_dp_r")npoints_dp_r = atoi(par_value.c_str());
      else if (par_name == "density_profile_r_output_file") density_profile_r_output_file = par_value.c_str();
      else if (par_name == "scale_dp_k")scale_dp_k = par_value.c_str();
      else if (par_name == "kmin_dp")kmin_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "kmax_dp")kmax_dp = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "npoints_dp_k")npoints_dp_k = atoi(par_value.c_str());
      else if (par_name == "density_profile_k_output_file") density_profile_k_output_file = par_value.c_str();
    	}
   }
//---------------------------------------end of while-----------------
// --------------------------------------------------------------------
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
  if(this->Number_of_MultiLevels > this->list_Props_Tolerance_MultiLevels.size() || this->Number_of_MultiLevels < this->list_Props_Tolerance_MultiLevels.size())
  {
    cerr << "Error in parameters for Multilevel. Not the same amount of entries" << endl;
    cerr << "Please check input parameter file." << endl;
    cerr<<RED<<this->Number_of_MultiLevels<<"  "<<this->list_Props_Tolerance_MultiLevels.size()<<endl;
    exit(1);
  }
#endif
#ifdef _USE_TWO_REFS_MOCKS_
  if(Number_of_references<Number_of_new_mocks)
  {
    cout<<RED<<"Warning: NUmber of references smaller than number of mocks to build simulstaneously"<<RESET<<std::endl;
  }
#endif
  if(this->mass_assignment_scheme=="NGP") {
    this->mass_assignment = 0;
  }
  else if(this->mass_assignment_scheme=="CIC") {
    this->mass_assignment = 1;
  }
  else if(this->mass_assignment_scheme=="TSC") {
    this->mass_assignment = 2;
  }
  else if(this->mass_assignment_scheme=="PCS") {
    this->mass_assignment = 3;
  }
// --------------------------------------------------------------------
// --------------------------------------------------------------------
#ifdef _USE_COSMO_PARS_
  this->om_matter = COSMOPARS::Om_matter;
  this->om_radiation = COSMOPARS::Om_radiation;
  this->om_baryons = COSMOPARS::Om_baryons;
  this->om_cdm = COSMOPARS::Om_cdm;
  this->om_vac = COSMOPARS::Om_vac;
  this->om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->spectral_index = COSMOPARS::spectral_index;
  this->wde_eos = COSMOPARS::wde_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = COSMOPARS::use_wiggles;
  this->RR = COSMOPARS::RR;
  this->alpha_s=COSMOPARS::alpha_s;
  this->Delta_SO=COSMOPARS::Delta_SO;
  this->f_baryon=COSMOPARS::f_baryon;
#endif
  this->s_cosmo_pars.cosmological_redshift=this->redshift;
  this->s_cosmo_pars.Hubble=this->Hubble;
  this->s_cosmo_pars.hubble=this->hubble;
  this->s_cosmo_pars.Om_matter =this->om_matter;
  this->s_cosmo_pars.Om_cdm =this->om_cdm;
  this->s_cosmo_pars.Om_baryons =this->om_baryons;
  this->s_cosmo_pars.Om_radiation =this->om_radiation;
  this->s_cosmo_pars.Om_vac =this->om_vac;
  this->s_cosmo_pars.Om_k =this->om_k;
  this->s_cosmo_pars.spectral_index =this->spectral_index;
  this->s_cosmo_pars.wde_eos =this->wde_eos;
  this->s_cosmo_pars.N_eff =this->N_eff;
  this->s_cosmo_pars.sigma8 =this->sigma8;
  this->s_cosmo_pars.f_baryon =this->om_baryons/this->om_matter;
  this->s_cosmo_pars.use_wiggles =this->use_wiggles;
  this->s_cosmo_pars.RR =this->RR;
  this->s_cosmo_pars.Tcmb =this->Tcmb;
  this->s_cosmo_pars.kmin_int = 0.00001;
  this->s_cosmo_pars.kmax_int = 100.0;
  this->s_cosmo_pars.use_external_power=false;
  this->s_cosmo_pars.GAL_BIAS=this->GAL_BIAS;
  this->s_cosmo_pars.Amc=this->Amc;
  this->s_cosmo_pars.kstar=this->kstar;
  this->s_cosmo_pars.Delta_SO=this->Delta_SO;
  this->s_cosmo_pars.use_wiggles=this->use_wiggles;
  this->s_cosmo_pars.kmin_int=this->kmin_integration;
  this->s_cosmo_pars.kmax_int=this->kmax_integration;
  this->s_cosmo_pars.mass_function_fit=this->mass_function_fit;
  this->s_cosmo_pars.halo_mass_bias_fit=this->halo_mass_bias_fit;
  this->s_cosmo_pars.density_profile=this->density_profile;
  this->s_cosmo_pars.coef_concentration_amp=this->coef_concentration_amp;
  this->s_cosmo_pars.coef_concentration=this->coef_concentration;
  this->s_cosmo_pars.mmin_hod=this->mmin_hod;
  this->s_cosmo_pars.scatter_hod=this->scatter_hod;
  this->s_cosmo_pars.hod_model=this->hod_model;
  this->s_cosmo_pars.alpha_hod=this->alpha_hod;
  this->s_cosmo_pars.alpha_A=alpha_A;
  this->s_cosmo_pars.muno_hod=this->muno_hod;
  this->s_cosmo_pars.ms_hod=this->muno_hod;
  this->s_cosmo_pars.s_bright_hod=this->s_bright_hod;
  this->s_cosmo_pars.s_faint_hod=this->s_faint_hod;
  this->s_cosmo_pars.width_hod=this->width_hod;
  this->s_cosmo_pars.Mstep_hod=this->Mstep_hod;
// --------------------------------------------------------------------
  this->d1=this->Lbox/static_cast<real_prec>(this->Nft);		/* grid spacing x-direction */
  this->d2=this->Lbox/static_cast<real_prec>(this->Nft);		/* grid spacing y-direction */
  this->d3=this->Lbox/static_cast<real_prec>(this->Nft);		/* grid spacing z-direction */
  this->Xoffset = this->Lbox/2; /* Xoffset for a box*/
  this->Yoffset = this->Lbox/2;
  this->Zoffset = this->Lbox/2;
  this->xmax=this->xmin+this->Lbox; /*xmax for a box */
  this->ymax=this->ymin+this->Lbox;
  this->zmax=this->zmin+this->Lbox;
  this->d1_HR=this->Lbox/static_cast<real_prec>(this->Nft_HR);		/* grid spacing x-direction */
  this->d2_HR=this->Lbox/static_cast<real_prec>(this->Nft_HR);		/* grid spacing y-direction */
  this->d3_HR=this->Lbox/static_cast<real_prec>(this->Nft_HR);		/* grid spacing z-direction */
  this->d1_low=this->Lbox/static_cast<real_prec>(this->Nft_low);		/* grid spacing x-direction */
  this->d2_low=this->Lbox/static_cast<real_prec>(this->Nft_low);		/* grid spacing y-direction */
  this->d3_low=this->Lbox/static_cast<real_prec>(this->Nft_low);		/* grid spacing z-direction */
  this->d1_JK=this->Lbox/static_cast<real_prec>(this->Nft_JK);		/* grid spacing x-direction */
  this->d2_JK=this->Lbox/static_cast<real_prec>(this->Nft_JK);		/* grid spacing y-direction */
  this->d3_JK=this->Lbox/static_cast<real_prec>(this->Nft_JK);		/* grid spacing z-direction */
#ifdef _USE_MULTISCALE_LEVEL_4_
  this->d1_low=this->Lbox/static_cast<real_prec>(this->Nft_low_l4);		/* grid spacing x-direction */
  this->d2_low=this->Lbox/static_cast<real_prec>(this->Nft_low_l4);		/* grid spacing y-direction */
  this->d3_low=this->Lbox/static_cast<real_prec>(this->Nft_low_l4);		/* grid spacing z-direction */
#endif
#ifndef _USE_CWC_
#ifdef _USE_MASS_KNOTS_
  this->cwt_used.push_back(1);
#elif !defined  _USE_MASS_KNOTS_
  this->cwt_used.push_back(0);
#endif
#endif
  this->n_cwt=this->cwt_used.size();
  
#ifdef _USE_CWC_
  if(this->n_cwt==1)
    {
      if(this->cwt_used[0]==0)
    	{
	  cout<<RED<<"Warning: _USE_CWC_ is defined but not specified in parameter file"<<RESET<<std::endl;
	  exit(1);
	}
    }
#endif
#ifndef _USE_CW_V
#ifdef _USE_VEL_KNOTS_V_
  this->cwv_used.push_back(1);
#elif !defined  _USE_VEL_KNOTS_V_
  this->cwv_used.push_back(0);
#endif
#endif
 this->n_cwv=this->cwv_used.size();
#ifdef _USE_CWC_V_
  if(this->n_cwv==1)
    {
      if(this->cwv_used[0]==0){
	cout<<RED<<"Warning: _USE_CWC_V_ is defined but not specified in parameter file"<<RESET<<std::endl;
	exit(0);
      }
    }
#endif
#ifndef _USE_V_DISP_VEL
  this->n_vknot_massbin=1;
#endif
#ifndef _USE_MASS_KNOTS_
  this->n_sknot_massbin=1;
#endif
  // These lines forget about the intervals put in the parameter file for the analysis of properties
  // and resize them in case we want bins with equal number of tracers
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
  if(true==this->set_bins_equal_number_tracers_main_property)
    {
      this->MASSbins_max.clear();this->MASSbins_max.shrink_to_fit();this->MASSbins_max.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
      this->MASSbins_min.clear();this->MASSbins_min.shrink_to_fit();this->MASSbins_min.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
      this->NMASSbins_power= this->Number_of_bins_equal_number_tracers_main_property;
   }
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
  if(true==this->set_bins_equal_number_tracers_main_property)
    {
      this->VMAXbins_min.clear();this->VMAXbins_min.shrink_to_fit();this->VMAXbins_min.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
      this->VMAXbins_max.clear();this->VMAXbins_max.shrink_to_fit();this->VMAXbins_max.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
      this->NVMAXbins_power= this->Number_of_bins_equal_number_tracers_main_property;
   }
#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_ // iN THIS CASE PEAK HEIGHT WILL TAKE THE ROLE OF MASS
  if(true==this->set_bins_equal_number_tracers_main_property)
    {
  this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
  this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(this->Number_of_bins_equal_number_tracers_main_property+1);
  this->NPHbins_power= this->Number_of_bins_equal_number_tracers_main_property;
}
#endif
  if(true==this->set_bins_equal_number_tracers)
    {
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
      this->VMAXbins_min.clear();this->VMAXbins_min.shrink_to_fit();this->VMAXbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->VMAXbins_max.clear();this->VMAXbins_max.shrink_to_fit();this->VMAXbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NVMAXbins_power= this->Number_of_bins_equal_number_tracers;
      this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NPHbins_power= this->Number_of_bins_equal_number_tracers;
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
      this->MASSbins_max.clear();this->MASSbins_max.shrink_to_fit();this->MASSbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->MASSbins_min.clear();this->MASSbins_min.shrink_to_fit();this->MASSbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NMASSbins_power= this->Number_of_bins_equal_number_tracers;
      this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NPHbins_power= this->Number_of_bins_equal_number_tracers;

#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_ // iN THIS CASE PEAK HEIGHT WILL TAKE THE ROLE OF MASS
      this->VMAXbins_min.clear();this->VMAXbins_min.shrink_to_fit();this->VMAXbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->VMAXbins_max.clear();this->VMAXbins_max.shrink_to_fit();this->VMAXbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NVMAXbins_power= this->Number_of_bins_equal_number_tracers;
      this->MASSbins_max.clear();this->MASSbins_max.shrink_to_fit();this->MASSbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->MASSbins_min.clear();this->MASSbins_min.shrink_to_fit();this->MASSbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      this->NMASSbins_power= this->Number_of_bins_equal_number_tracers;
#endif
     //-------------------------
      if(this->i_rs_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->RSbins_min.clear();this->RSbins_min.shrink_to_fit();this->RSbins_min.resize(1);
          this->RSbins_max.clear();this->RSbins_max.shrink_to_fit();this->RSbins_max.resize(1);
          this->RSbins_max[0]=1e10;
          this->RSbins_min[0]=0;
          this->NRSbins_power=1;
      }
      else
      {
        this->NRSbins_power=this->Number_of_bins_equal_number_tracers;
        this->RSbins_max.clear();this->RSbins_max.shrink_to_fit();this->RSbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->RSbins_min.clear();this->RSbins_min.shrink_to_fit();this->RSbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      }
          //-------------------------
      if(this->i_rvir_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->RVIRbins_min.clear();this->RVIRbins_min.shrink_to_fit();this->RVIRbins_min.resize(1);
          this->RVIRbins_max.clear();this->RVIRbins_max.shrink_to_fit();this->RVIRbins_max.resize(1);
          this->RVIRbins_max[0]=1e10;
          this->RVIRbins_min[0]=0;
          this->NRVIRbins_power=1;
      }
      else
      {
          this->RVIRbins_max.clear();this->RVIRbins_max.shrink_to_fit();this->RVIRbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->RVIRbins_min.clear();this->RVIRbins_min.shrink_to_fit();this->RVIRbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
          this->NRVIRbins_power=this->Number_of_bins_equal_number_tracers;
      }
      //-------------------------
      if(this->i_rvir_g<0 || this->i_rs_g<0)// If we do not use Rvir or Rs, we cannot define concentration
      {
          this->CONCENTRATIONbins_min.clear();this->CONCENTRATIONbins_min.shrink_to_fit();this->CONCENTRATIONbins_min.resize(1);
          this->CONCENTRATIONbins_max.clear();this->CONCENTRATIONbins_max.shrink_to_fit();this->CONCENTRATIONbins_max.resize(1);
          this->CONCENTRATIONbins_max[0]=1e10;
          this->CONCENTRATIONbins_min[0]=0;
          this->NCONCENTRATIONbins_power=1;
      }
      else
      {
        this->CONCENTRATIONbins_max.clear();this->CONCENTRATIONbins_max.shrink_to_fit();this->CONCENTRATIONbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->CONCENTRATIONbins_min.clear();this->CONCENTRATIONbins_min.shrink_to_fit();this->CONCENTRATIONbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NCONCENTRATIONbins_power=this->Number_of_bins_equal_number_tracers;
      }
      //-------------------------
      if(this->i_spin_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->SPINbins_min.clear();this->SPINbins_min.shrink_to_fit();this->SPINbins_min.resize(1);
          this->SPINbins_max.clear();this->SPINbins_max.shrink_to_fit();this->SPINbins_max.resize(1);
          this->SPINbins_max[0]=1e10;
          this->SPINbins_min[0]=0;
          this->NSPINbins_power=1;
      }
      else
      {
        this->NSPINbins_power=this->Number_of_bins_equal_number_tracers;
        this->SPINbins_max.clear();this->SPINbins_max.shrink_to_fit();this->SPINbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->SPINbins_min.clear();this->SPINbins_min.shrink_to_fit();this->SPINbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      }
      //-------------------------
      if(this->i_spin_bullock_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->SPINBULLOCKbins_min.clear();this->SPINBULLOCKbins_min.shrink_to_fit();this->SPINBULLOCKbins_min.resize(1);
          this->SPINBULLOCKbins_max.clear();this->SPINBULLOCKbins_max.shrink_to_fit();this->SPINBULLOCKbins_max.resize(1);
          this->SPINBULLOCKbins_max[0]=1e10;
          this->SPINBULLOCKbins_min[0]=0;
          this->NSPINBULLOCKbins_power=1;
      }
        else
      {
          this->NSPINBULLOCKbins_power=this->Number_of_bins_equal_number_tracers;
          this->SPINBULLOCKbins_max.clear();this->SPINBULLOCKbins_max.shrink_to_fit();this->SPINBULLOCKbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->SPINBULLOCKbins_min.clear();this->SPINBULLOCKbins_min.shrink_to_fit();this->SPINBULLOCKbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      }
          //-------------------------
      if(this->i_virial_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->VIRIALbins_min.clear();this->VIRIALbins_min.shrink_to_fit();this->VIRIALbins_min.resize(1);
          this->VIRIALbins_max.clear();this->VIRIALbins_max.shrink_to_fit();this->VIRIALbins_max.resize(1);
          this->VIRIALbins_max[0]=1e10;
          this->VIRIALbins_min[0]=0;
          this->NVIRIALbins_power=1;
      }
      else
      {
          this->VIRIALbins_max.clear();this->VIRIALbins_max.shrink_to_fit();this->VIRIALbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->VIRIALbins_min.clear();this->VIRIALbins_min.shrink_to_fit();this->VIRIALbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
          this->NVIRIALbins_power= this->Number_of_bins_equal_number_tracers;
        }
    //-------------------------
      this->VRMSbins_max.clear();this->VRMSbins_max.shrink_to_fit();this->VRMSbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
      this->VRMSbins_min.clear();this->VRMSbins_min.shrink_to_fit();this->VRMSbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
      if(this->i_vrms_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->VRMSbins_min.clear();this->VRMSbins_min.shrink_to_fit();this->VRMSbins_min.resize(1);
          this->VRMSbins_max.clear();this->VRMSbins_max.shrink_to_fit();this->VRMSbins_max.resize(1);
          this->VRMSbins_max[0]=1e10;
          this->VRMSbins_min[0]=0;
          this->NVRMSbins_power=1;
      }
    else
      {
          this->VRMSbins_max.clear();this->VRMSbins_max.shrink_to_fit();this->VRMSbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->VRMSbins_min.clear();this->VRMSbins_min.shrink_to_fit();this->VRMSbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
          this->NVRMSbins_power= this->Number_of_bins_equal_number_tracers;
        }
      //-------------------------
      if(this->i_b_to_a_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->BTOAbins_min.clear();this->BTOAbins_min.shrink_to_fit();this->BTOAbins_min.resize(1);
          this->BTOAbins_max.clear();this->BTOAbins_max.shrink_to_fit();this->BTOAbins_max.resize(1);
          this->BTOAbins_max[0]=1e10;
          this->BTOAbins_min[0]=0;
          this->NBTOAbins_power=1;
      }
    else
      {
          this->BTOAbins_max.clear();this->BTOAbins_max.shrink_to_fit();this->BTOAbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->BTOAbins_min.clear();this->BTOAbins_min.shrink_to_fit();this->BTOAbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
          this->NBTOAbins_power= this->Number_of_bins_equal_number_tracers;
      }

      //-------------------------
      if(this->i_c_to_a_g<0)// If we do not use that property, we need to set fictitious limits
      {
          this->CTOAbins_min.clear();this->CTOAbins_min.shrink_to_fit();this->CTOAbins_min.resize(1);
          this->CTOAbins_max.clear();this->CTOAbins_max.shrink_to_fit();this->CTOAbins_max.resize(1);
          this->CTOAbins_max[0]=1e10;
          this->CTOAbins_min[0]=0;
          this->NCTOAbins_power=1;
      }
      else
      {
          this->CTOAbins_max.clear();this->CTOAbins_max.shrink_to_fit();this->CTOAbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
          this->CTOAbins_min.clear();this->CTOAbins_min.shrink_to_fit();this->CTOAbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
          this->NCTOAbins_power= this->Number_of_bins_equal_number_tracers;
        }
      //-------------------------
      if(true==this->Get_tracer_local_mach_number)
      {
         this->MACHbins_min.clear();this->MACHbins_min.shrink_to_fit();this->MACHbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
         this->MACHbins_max.clear();this->MACHbins_max.shrink_to_fit();this->MACHbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
         this->NMACHbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else
      {
          this->MACHbins_min.clear();this->MACHbins_min.shrink_to_fit();this->MACHbins_min.resize(1);
          this->MACHbins_max.clear();this->MACHbins_max.shrink_to_fit();this->MACHbins_max.resize(1);
          this->MACHbins_max[0]=1e10;
          this->MACHbins_min[0]=-1e10;
          this->NMACHbins_power=1;
      }
      //-------------------------
      if(true==this->Get_tracer_local_dach_number)
      {
         this->DACHbins_min.clear();this->DACHbins_min.shrink_to_fit();this->DACHbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
         this->DACHbins_max.clear();this->DACHbins_max.shrink_to_fit();this->DACHbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
         this->NDACHbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->DACHbins_min.clear();this->DACHbins_min.shrink_to_fit();this->DACHbins_min.resize(1);
          this->DACHbins_max.clear();this->DACHbins_max.shrink_to_fit();this->DACHbins_max.resize(1);
          this->DACHbins_max[0]=1e10;
          this->DACHbins_min[0]=-1e10;
          this->NDACHbins_power=1;
      }
     //-------------------------
      if(true==this->Get_tracer_bias)
      {
        this->BIASbins_min.clear();this->BIASbins_min.shrink_to_fit();this->BIASbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->BIASbins_max.clear();this->BIASbins_max.shrink_to_fit();this->BIASbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NBIASbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->BIASbins_min.clear();this->BIASbins_min.shrink_to_fit();this->BIASbins_min.resize(1);
          this->BIASbins_max.clear();this->BIASbins_max.shrink_to_fit();this->BIASbins_max.resize(1);
          this->BIASbins_max[0]=1e10;
          this->BIASbins_min[0]=-1e10;
          this->NBIASbins_power=1;
      }
      //-------------------------
      if(true==this->Get_tracer_relative_bias)
      {
        this->RBIASbins_min.clear();this->RBIASbins_min.shrink_to_fit();this->RBIASbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->RBIASbins_max.clear();this->RBIASbins_max.shrink_to_fit();this->RBIASbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NRBIASbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->RBIASbins_min.clear();this->RBIASbins_min.shrink_to_fit();this->RBIASbins_min.resize(1);
          this->RBIASbins_max.clear();this->RBIASbins_max.shrink_to_fit();this->RBIASbins_max.resize(1);
          this->RBIASbins_max[0]=1e10;
          this->RBIASbins_min[0]=-1e10;
          this->NRBIASbins_power=1;
      }
      //-------------------------
      if(true==this->Get_tracer_quadratic_bias)
      {
        this->QBIASbins_min.clear();this->QBIASbins_min.shrink_to_fit();this->QBIASbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->QBIASbins_max.clear();this->QBIASbins_max.shrink_to_fit();this->QBIASbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NQBIASbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->QBIASbins_min.clear();this->QBIASbins_min.shrink_to_fit();this->QBIASbins_min.resize(1);
          this->QBIASbins_max.clear();this->QBIASbins_max.shrink_to_fit();this->QBIASbins_max.resize(1);
          this->QBIASbins_max[0]=1e10;
          this->QBIASbins_min[0]=-1e10;
          this->NQBIASbins_power=1;
      }
      //-------------------------
      if(true==this->Get_local_overdensity)
      {
        this->LCbins_min.clear();this->LCbins_min.shrink_to_fit();this->LCbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->LCbins_max.clear();this->LCbins_max.shrink_to_fit();this->LCbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NLCbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->LCbins_min.clear();this->LCbins_min.shrink_to_fit();this->LCbins_min.resize(1);
          this->LCbins_max.clear();this->LCbins_max.shrink_to_fit();this->LCbins_max.resize(1);
          this->LCbins_max[0]=1e10;
          this->LCbins_min[0]=-1e10;
          this->NLCbins_power=1;
      }
      //-------------------------
      if(true==this->Get_tidal_anisotropy_at_halo)
      {
        this->TAbins_min.clear();this->TAbins_min.shrink_to_fit();this->TAbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->TAbins_max.clear();this->TAbins_max.shrink_to_fit();this->TAbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NTAbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->TAbins_min.clear();this->TAbins_min.shrink_to_fit();this->TAbins_min.resize(1);
          this->TAbins_max.clear();this->TAbins_max.shrink_to_fit();this->TAbins_max.resize(1);
          this->TAbins_max[0]=1e10;
          this->TAbins_min[0]=-1e10;
          this->NTAbins_power=1;
      }
      //-------------------------
      if(true==this->Get_tracer_local_dm_density)
      {
        this->LOCALDMbins_min.clear();this->LOCALDMbins_min.shrink_to_fit();this->LOCALDMbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->LOCALDMbins_max.clear();this->LOCALDMbins_max.shrink_to_fit();this->LOCALDMbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NLOCALDMbins_power= this->Number_of_bins_equal_number_tracers;
      }
    else{
          this->LOCALDMbins_min.clear();this->LOCALDMbins_min.shrink_to_fit();this->LOCALDMbins_min.resize(1);
          this->LOCALDMbins_max.clear();this->LOCALDMbins_max.shrink_to_fit();this->LOCALDMbins_max.resize(1);
          this->LOCALDMbins_max[0]=1e10;
          this->LOCALDMbins_min[0]=-1e10;
          this->NLOCALDMbins_power=1;
      }
      //-------------------------
#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
      if(true==this->Get_peak_height_at_halo)
      {
        this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(this->Number_of_bins_equal_number_tracers+1);
        this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(this->Number_of_bins_equal_number_tracers+1);
        this->NPHbins_power= this->Number_of_bins_equal_number_tracers;
      }
      else{
          this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(1);
          this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(1);
          this->PHbins_max[0]=1e10;
          this->PHbins_min[0]=-1e10;
          this->NPHbins_power=1;
      }
#endif
  } // closes if for quartiles
   if(false==this->Get_tracer_local_mach_number)
    {
       this->MACHbins_min.clear();this->MACHbins_min.shrink_to_fit();this->MACHbins_min.resize(1);
       this->MACHbins_max.clear();this->MACHbins_max.shrink_to_fit();this->MACHbins_max.resize(1);
       this->MACHbins_max[0]=1e10;
       this->MACHbins_min[0]=-1e10;
       this->NMACHbins_power=1;
    }
   if(false==this->Get_tracer_local_dach_number)
    {
       this->DACHbins_min.clear();this->DACHbins_min.shrink_to_fit();this->DACHbins_min.resize(1);
       this->DACHbins_max.clear();this->DACHbins_max.shrink_to_fit();this->DACHbins_max.resize(1);
       this->DACHbins_max[0]=1e10;
       this->DACHbins_min[0]=-1e10;
       this->NDACHbins_power=1;
    }
   if(false==this->Get_local_overdensity)
    {
       this->LCbins_min.clear();this->LCbins_min.shrink_to_fit();this->LCbins_min.resize(1);
       this->LCbins_max.clear();this->LCbins_max.shrink_to_fit();this->LCbins_max.resize(1);
       this->LCbins_max[0]=1e10;
       this->LCbins_min[0]=-1e10;
       this->NLCbins_power=1;
    }
   if(false==this->Get_tracer_bias)
    {
       this->BIASbins_min.clear();this->BIASbins_min.shrink_to_fit();this->BIASbins_min.resize(1);
       this->BIASbins_max.clear();this->BIASbins_max.shrink_to_fit();this->BIASbins_max.resize(1);
       this->BIASbins_max[0]=1e10;
       this->BIASbins_min[0]=-1e10;
       this->NBIASbins_power=1;
       this->RBIASbins_min.clear();this->RBIASbins_min.shrink_to_fit();this->RBIASbins_min.resize(1);
       this->RBIASbins_max.clear();this->RBIASbins_max.shrink_to_fit();this->RBIASbins_max.resize(1);
       this->RBIASbins_max[0]=1e10;
       this->RBIASbins_min[0]=-1e10;
       this->NRBIASbins_power=1;
    }
   if(false==this->Get_tidal_anisotropy_at_halo)
    {
       this->TAbins_min.clear();this->TAbins_min.shrink_to_fit();this->TAbins_min.resize(1);
       this->TAbins_max.clear();this->TAbins_max.shrink_to_fit();this->TAbins_max.resize(1);
       this->TAbins_max[0]=1e10;
       this->TAbins_min[0]=-1e10;
       this->NTAbins_power=1;
    }
   if(false==this->Get_tracer_local_dm_density)
   {
     this->LOCALDMbins_min.clear();this->LOCALDMbins_min.shrink_to_fit();this->LOCALDMbins_min.resize(1);
     this->LOCALDMbins_max.clear();this->LOCALDMbins_max.shrink_to_fit();this->LOCALDMbins_max.resize(1);
     this->LOCALDMbins_max[0]=1e10;
     this->LOCALDMbins_min[0]=-1e10;
     this->NLOCALDMbins_power=1;
   }

#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_ // iN THIS CASE PEAK HEIGHT WILL TAKE THE ROLE OF MASS
   if(false==this->Get_peak_height_at_halo)
    {
       this->PHbins_min.clear();this->PHbins_min.shrink_to_fit();this->PHbins_min.resize(1);
       this->PHbins_max.clear();this->PHbins_max.shrink_to_fit();this->PHbins_max.resize(1);
       this->PHbins_max[0]=1e10;
       this->PHbins_min[0]=-1e10;
       this->NPHbins_power=1;
    }
#endif
  //--------------------------------------------------------------------------
  // append directory names in front of filenames
  string dat = ".txt";
  string ll = "_";
  string zr = "_redshift_"+ std::to_string(static_cast<float>(redshift));
  string output = "../output_cosmolib/";
  this->mass_function_output_file = output+mass_function_output_file+ll+mass_function_fit+zr+dat;
  this->halo_mass_bias_output_file = output+halo_mass_bias_output_file+ll+halo_mass_bias_fit+zr+dat;
  this->effective_halo_mass_bias_output_file = output+effective_halo_mass_bias_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  this->effective_halo_mean_number_density_output_file=output+effective_halo_mean_number_density_output_file+ll+mass_function_fit+zr+dat;
  this->galaxy_power_spectrum_halo_model_output_file=output+galaxy_power_spectrum_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  this->galaxy_correlation_function_halo_model_output_file=output+galaxy_correlation_function_halo_model_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat;
  this->linear_matter_cf_output_file= output+linear_matter_cf_output_file+zr+dat;
  this->non_linear_matter_cf_halo_fit_output_file=output+non_linear_matter_cf_halo_fit_output_file+zr+dat;
  this->linear_matter_ps_output_file=output+linear_matter_ps_output_file+zr+dat;
  this->non_linear_matter_ps_halo_fit_output_file=output+non_linear_matter_ps_halo_fit_output_file+zr+dat;
  this->non_linear_matter_ps_pt_output_file=output+non_linear_matter_ps_pt_output_file+zr+dat;
  this->density_profile_r_output_file=output+density_profile_r_output_file+ll+density_profile+zr+dat;
  this->density_profile_k_output_file=output+density_profile_k_output_file+ll+density_profile+zr+dat;
 //--------------------------------------------------------------------------
   this->warnings();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Params::derived_pars(){
  this->NGRID = static_cast<ULONG>(this->Nft*this->Nft*this->Nft);
  this->NGRID_low = static_cast<ULONG>(this->Nft_low*this->Nft_low*this->Nft_low);
  this->NGRID_HR = static_cast<ULONG>(this->Nft_HR*this->Nft_HR*this->Nft_HR);
  this->NGRID_JK = static_cast<ULONG>(this->Nft_JK*this->Nft_JK*this->Nft_JK);
  this->NGRID_h = static_cast<ULONG>(this->Nft*this->Nft*(this->Nft/2+1));  // h is for half
  this->delta_x = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->delta_y = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->delta_z = this->Lbox/ (static_cast<double>(this->Nft)); 
  this->delta_x_low = this->Lbox_low/(static_cast<double>(this->Nft_low)); 
  this->delta_y_low = this->Lbox_low/(static_cast<double>(this->Nft_low));
  this->delta_z_low = this->Lbox_low/(static_cast<double>(this->Nft_low)); 
  this->deltak_x = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_y = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_z = 2.*M_PI/static_cast<double>(this->Lbox);
  this->deltak_x_low = 2.*M_PI/static_cast<double>(this->Lbox_low);
  this->deltak_y_low = 2.*M_PI/static_cast<double>(this->Lbox_low);
  this->deltak_z_low = 2.*M_PI/static_cast<double>(this->Lbox_low);
  this->deltak_0 = sqrt(pow(this->deltak_x,2)+pow(this->deltak_y,2)+pow(this->deltak_z,2))/sqrt(3.0); 
  this->deltak_0_low = sqrt(pow(this->deltak_x_low,2)+pow(this->deltak_y_low,2)+pow(this->deltak_z_low,2))/sqrt(3.0); 
  this->kmin     = this->DeltaKmin*this->deltak_0;
  this->kmax     = sqrt(pow(0.5*this->Nft*this->deltak_x,2)+pow(0.5*this->Nft*this->deltak_y,2)+pow(0.5*this->Nft*this->deltak_z,2))/sqrt(3.0);         
  this->DeltaK_data   = this->ndel_data*this->deltak_0; 
  this->DeltaK_data_low   = this->ndel_data*this->deltak_0_low; 
  this->DeltaK_window = this->ndel_window*this->deltak_0; 
  this->Nnp_data      = (this->type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_data/2); 
  this->Nnp_window    = (this->type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_window/ 2); 
  this->Deltal        = log10(this->kmax/this->kmin)/static_cast<double>(this->N_log_bins); 
  this->Deltamu       = 2.0/(static_cast<real_prec>(this->N_mu_bins));
#ifdef _USE_SLICS_PK_BINNING_
  this->kmin = MIN_K;
  this->deltak_0 = DELTA_K; 
  this->DeltaK_data = this->deltak_0; 
  this->Nnp_data      = (this->type_of_binning=="log" ? this->N_log_bins : static_cast<int>(floor((M_PI*this->Nft/this->Lbox)/DELTA_K))) ; //* The new number of modes is the Nyq freq divided by DeltaK imposed
  this->file_power+="_offibinning";
#endif

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Params::show_params()
{
  std::cout<<BOLDYELLOW;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"CosmiCodes                                                                       *"<<std::endl;
  std::cout<<"Input values of parameters in parameter file                              *"<<std::endl;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"******************NUMERICAL PARAMTERS**************************************"<<RESET<<std::endl;
  for(int i=0;i<parameter_number.size();++i)
    std::cout<<YELLOW<<parameter_number[i].first<<" = "<<BLUE<<parameter_number[i].second<<RESET<<endl;
  std::cout<<YELLOW<<"******************STRING PARAMETERS****************************************"<<RESET<<std::endl;
  for(int i=0;i<parameter_string.size();++i)
    std::cout<<YELLOW<<parameter_string[i].first<<" = "<<BLUE<<parameter_string[i].second<<RESET<<endl;
  std::cout<<YELLOW<<"******************BOOLEAN PARAMETERS***************************************"<<RESET<<std::endl;
  for(int i=0;i<parameter_boolean.size();++i)
    std::cout<<YELLOW<<parameter_boolean[i].first<<" = "<<BLUE<<parameter_boolean[i].second<<RESET<<endl;
  std::cout<<YELLOW<<"******************VECTOR PARAMETERS****************************************"<<RESET<<std::endl;
  for(int i=0;i<parameter_vectors.size();++i)
  {
   std::cout<<YELLOW<<parameter_vectors[i].first<<" = ";
   for(int j=0;j<parameter_vectors[i].second.size();++j)cout<<BLUE<<parameter_vectors[i].second[j]<<" "<<RESET;
   std::cout<<endl;
  }
  std::cout<<YELLOW<<"******************VECTOR-STRING PARAMETERS***********************************"<<RESET<<std::endl;
  for(int i=0;i<parameter_vector_string.size();++i)
  {
   std::cout<<YELLOW<<parameter_vector_string[i].first<<" = ";
   for(int j=0;j<parameter_vector_string[i].second.size();++j)cout<<BLUE<<parameter_vector_string[i].second[j]<<" "<<RESET;
   std::cout<<endl;
  }
  std::cout<<YELLOW<<"******************SUMMARY PARAMETERS***************************************"<<RESET<<std::endl;
  std::cout<<YELLOW<<"Number of numerical parameters = "<<parameter_number.size()<<endl;
  std::cout<<YELLOW<<"Number of string parameters = "<<parameter_string.size()<<endl;
  std::cout<<YELLOW<<"Number of boolean parameters = "<<parameter_boolean.size()<<endl;
  std::cout<<YELLOW<<"Number of vector parameters = "<<parameter_vectors.size()<<endl;
  std::cout<<YELLOW<<"Total = "<<parameter_boolean.size()+parameter_string.size()+parameter_number.size()+parameter_vectors.size()<<endl;
  std::cout<<YELLOW<<"***************************************************************************"<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Params::warnings()
{
#ifdef USE_GALAXY_TOOLS
    cout<<RED<<endl;
    cout<<"USE_GALAXY_TOOLS is defined"<endl;
    cout<<"This is in conflict with the current set up "<endl;
    cout<<RESET<<endl;
    throw std::invalid_argument("Wrong set up. Fix and recompile");
#endif
    if(this->type_of_object!="TRACER" && this->type_of_object!="TRACER_REF" && this->type_of_object!="TRACER_MOCK" && this->type_of_object!="TRACER_MOCK_ONLY_COORDS" && this->type_of_object!="TRACER_REF_ONLY_COORDS")
     {
        cout<<RED<<endl;
        cout<<"Parameter Type_of_object = "<<this->type_of_object<<endl;
        cout<<"does not match any of the expected options"<<endl;
        cout<<"\tTRACER"<<endl;
        cout<<"\tTRACER_REF"<<endl;
        cout<<"\tTRACER_MOCK"<<endl;
        cout<<"\tTRACER_MOCK_ONLY_COORDS"<<endl;
        cout<<"\tTRACER_REF_ONLY_COORDS"<<endl;
        cout<<RESET<<endl;
        throw std::invalid_argument("Wrong input in parameter file");
    }
    if(this->use_file_nbar==true && this->nbar_tabulated ==true)
    {
        cout<<RED<<endl;
        cout<<"Parameter use_file_nbar = "<<this->use_file_nbar<<endl;
        cout<<"is in conflict with parameter nbar_tabulated = "<<this->nbar_tabulated<<endl;
        cout<<RESET<<endl;
        throw std::invalid_argument("Wrong input in parameter file");
    }
   if(this->Get_tracer_bias ==false)
    {
#ifdef  _USE_BIAS_OBJECT_TO_OBJECT_
      cout<<CYAN<<endl;
      cout<<"Parameter Get_tracer_bias  = "<<this->use_file_nbar<<endl;
      cout<<"is in conflict with parameter request from preproc directive _USE_BIAS_OBJECT_TO_OBJECT_ "<<endl;
      string ans;
      cout<<"Do you wish to continue? (y/n)"<<endl; cin>>ans;
      if(ans=="Y" || ans == "y" || ans=="yes")      
        cout<<RESET<<endl;
      else
        throw std::invalid_argument("Conflict between parameters not resolved.");
#endif
    }

  if(this->i_v1_g<0  || this->i_v2_g<0 || this->i_v3_g<0)
  {
#ifdef _REDSHIFT_SPACE_
      cout<<RED<<endl;
      cout<<"The opreproc directive _REDSHIFT_SPACE_ enabled but there is no info on velocities available"<<endl;
      string ans;
      cout<<"Do you wish to continue? (y/n)"<<endl; cin>>ans;
      if(ans=="Y" || ans == "y" || ans=="yes")      
        cout<<RESET<<endl;
      else
        throw std::invalid_argument("Conflict between parameters not resolved.");
#endif
  }


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
