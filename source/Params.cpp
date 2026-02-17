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
#define NN "NoName"
#define ON "on"
#define OFF "off"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Params::collect_params_info(string& pname, string& pname_c, string& section, string&description, string&options){
  s_properties_params mp;
  mp.par_name=pname;
  mp.par_name_in_code=pname_c;
  mp.options=options;
  mp.loading_sections=section;
  mp.description=description;
  this->parameter_properties.push_back(mp);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Params::explain_pars(string par_name)
{
  std::cout<<GREEN<<"CosmicCodes. Paramater Info. "<<par_name<<std::endl;
  for(int i=0; i<this->parameter_properties.size();++i)
  {
    if(parameter_properties[i].par_name==par_name || parameter_properties[i].par_name_in_code==par_name)
      parameter_properties[i].show();
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
void Params::init_pars() // Inizialization of parameters
{

this->redshift = 0.0;
  this->Initial_Redshift_DELTA=0;
  this->Initial_Redshift_SIM=0;
  this->IC_index =0;
  this->realization = 1;
//-----------------------------------------------------------------------------------------
  this->NX = 200;
  this->iMAS_X = 0;
  this->XNAME = "XNAME";
  this->Name_Property_X = "X_PROPERTY";
  this->iMAS_X_REF_PDF = 0;
  this->iMAS_X_NEW = 0;
  this->Type_of_file_X = "bin";
  this->Input_Directory_X = NN;
  this->Input_Directory_X_REF = NN;
  this->Input_Directory_X_REF_TWO = NN;
  this->Input_Directory_X_NEW = NN; 
  this->Input_Directory_X_new_ref = NN;
  this->Name_Catalog_X = NN;
  this->Name_Catalog_X_NEW = NN;
  this->Name_Catalog_X_REF_PDF = NN;
  this->Name_Catalog_X_new_ref = NN;
  this->Convert_Density_to_Delta_X = false;
  this->Name_VelFieldx_X = NN;
  this->Name_VelFieldy_X = NN;
  this->Name_VelFieldz_X  = NN;
  this->Scale_X= "linear";
  this->delta_X_min = 0;
  this->delta_X_max = 10.;
  this->ldelta_X_min = 0.;
  this->ldelta_X_max = 10.;
//-----------------------------------------------------------------------------------------
  this->NY = 100;
  this->iMAS_Y = 0;
  this->YNAME = "YNAME";
  this->Name_Property_Y = "Y_PROPERTY";
  this->Input_Directory_Y = NN;
  this->Input_Directory_Y_TWO = NN;
  this->Name_Catalog_Y = NN;
  this->Name_Catalog_Y_new_ref = NN;
  this->Name_Catalog_Y_HR = NN;
  this->Name_Catalog_Y_MWEIGHTED = NN;
  this->Convert_Density_to_Delta_Y = false;
  this->Scale_Y= "linear";
  this->delta_Y_min = 0.;
  this->delta_Y_max = 10;
  this->ldelta_Y_min = 0;
  this->ldelta_Y_max = 10;
  this->NY_MASS = 1;
  this->NY_SAT_FRAC = 1;
//-----------------------------------------------------------------------------------------
  this->Redefine_limits= false;
//-----------------------------------------------------------------------------------------
// CWC analysis
  this->Nlambdath=1;
  this->lambdath = 0.0;
  this->lambdath_v = 0.0;
  this->n_sknot_massbin = 1;
  this->n_vknot_massbin = 1;
  this->NMASSbins = 1;
  this->NMASSbins_power = 1;
//-----------------------------------------------------------------------------------------
// Bias analysis
  this->Write_PDF_number_counts =false;
  this->Comp_conditional_PDF = false;
  this->Comp_joint_PDF = false;
  this->write_files_for_histograms=false;
//-----------------------------------------------------------------------------------------
//BMT parameters
  this->N_iterations_Kernel = 1;
  this->Iteration_Kernel_average = 1;
  this->iteration_ini=0;
  this->N_dm_realizations=1;
  this->N_dm_initial=1;
  this->N_iterations_dm=1;
  this->Apply_Rankordering = false;
  this->Apply_Rankordering_ab_initio = false;
  this->Number_of_chunks_new_dm=1;
  this->Input_Directory_BIAS_KERNEL = NN;
  this->Input_Directory_BIAS_KERNEL_TWO = NN;
//-----------------------------------------------------------------------------------------
  this->Output_directory = "../output";
//-----------------------------------------------------------------------------------------
// Power Spectrum
  this->statistics = "Pk_fkp";
  this->Name_survey = "Zorbas";
  this->type_of_object = "TRACER";
  this->input_type = "catalog";
  this->input_type_two = NN;
  this->ngal_delta = 1;
  this->Input_dir_cat = NN;
  this->Input_dir_cat_TWO = NN;
  this->file_catalogue = NN;
  this->file_catalogue_new_ref = "cat.dat";
  this->file_random = NN;
  this->Number_of_random_files = 1;
  // ************************************************************************
  this->delta_grid_file = "delta_file";
  this->delta_grid_file2 = "delta_file";
  this->delta_grid_file3 = "delta_file";
  this->delta_grid_file4 = "delta_file";
  this->measure_cross = false;
  this->measure_cross_from_1 = 1;
  this->measure_cross_from_2 = 2;

  // ************************************************************************

  this->extra_info = "";
  this->vel_units_g = "kmps";
  // ************************************************************************

  this->sys_of_coord_dm = 1;
  this->i_coord1_dm = NOT_USED;
  this->i_coord2_dm = NOT_USED;
  this->i_coord3_dm = NOT_USED;
  this->i_v1_dm = NOT_USED;
  this->i_v2_dm = NOT_USED;
  this->i_v3_dm = NOT_USED;
  this->i_mass_dm = NOT_USED;

  // ************************************************************************
  this->sys_of_coord_g = 1;
  this->redshift_space_coords_g = false;
  this->angles_units_g = "D";
  this->i_coord1_g = NOT_USED;
  this->i_coord2_g = NOT_USED;
  this->i_coord3_g = NOT_USED;
  this->i_v1_g = NOT_USED;
  this->i_v2_g = NOT_USED;
  this->i_v3_g = NOT_USED;
  this->vel_units_g = "kmps";
  this->i_mass_g = NOT_USED;
  this->i_vmax_g = NOT_USED;
  this->i_vrms_g = NOT_USED;
  this->i_rs_g = NOT_USED;
  this->i_spin_g = NOT_USED;
  this->i_spin_bullock_g = NOT_USED;
  this->i_virial_g = NOT_USED;
  this->i_b_to_a_g = NOT_USED;
  this->i_c_to_a_g = NOT_USED;
  this->i_sf_g = NOT_USED;
  this->i_color_g = NOT_USED;
  this->i_stellar_mass_g=NOT_USED;
  this->i_mean_density_g = NOT_USED;
  this->i_abs_mag_g=NOT_USED;
  this->i_app_mag_g=NOT_USED;
  this->i_weight1_g = NOT_USED;
  this->use_weight1_g = false;
  this->i_weight2_g = NOT_USED;
  this->use_weight2_g = false;
  this->i_weight3_g = NOT_USED;
  this->use_weight4_g = false;
  this->i_weight4_g = NOT_USED;
  this->use_weight3_g = false;
  this->weight_with_mass = false;
  this->weight_vel_with_mass =false;
  // ************************************************************************
  this->use_random_catalog = false;
  this->use_random_catalog_cl = false;
  this->nbar_tabulated = false;
  this->use_file_nbar=false;
  this->file_nbar = NN;
  this->sys_of_coord_r = 1;
  this->angles_units_r = "D";
  this->i_coord1_r = NOT_USED;
  this->i_coord2_r = NOT_USED;
  this->i_coord3_r = NOT_USED;
  this->i_mass_r = NOT_USED;
  this->i_color_r=NOT_USED;
  this->i_stellar_mass_r=NOT_USED;
  this->i_abs_mag_r=NOT_USED;
  this->i_app_mag_r=NOT_USED;
  this->i_mean_density_r = NOT_USED;
  this->i_weight1_r = NOT_USED;
  this->use_weight1_r = false;
  this->i_weight2_r = NOT_USED;
  this->use_weight2_r = false;
  this->i_weight3_r = NOT_USED;
  this->use_weight3_r = false;
  this->i_weight4_r = NOT_USED;
  this->use_weight4_r = false;
  // ************************************************************************
  this->constant_depth = false;
  this->redshift_min_sample = 0;
  this->redshift_min_sample = 10;
  this->area_survey = 1000;
  this->i_mask_pixel= NOT_USED;
  this->i_mask_alpha= NOT_USED;
  this->i_mask_delta= NOT_USED;
  this->i_mask_flag= NOT_USED;
  // ************************************************************************
  this->new_los=false;
  this->Lbox = 100;
  this->Lbox_low = 100;
  this->new_Lbox = false;
  this->Nft = 256;
  this->Nft_HR = 256;
  this->Nft_low = 256;
  this->Nft_JK = 16;
  this->mass_assignment_scheme = "NGP";
  this->MAS_correction = false;
  this->type_of_binning = "linear";
  this->N_log_bins = 10;
  this->DeltaKmin = 0;
  this->ndel_data=1;
  this->ndel_window=1;
  this->FKP_weight = false;
  this->Pest = 20000;
  this->N_mu_bins=10;
  this->SN_correction = false;
  this->FKP_error_bars = false;
  this->FKP_error_bars_exact = false;
  // ************************************************************************
  this->Nbins_redshift = 10;
  this->N_dndz_bins = 2;
  this->new_N_dndz_bins = 2;
  this->file_dndz = "dndz";
  // ************************************************************************
  this->file_power = "power";
  this->file_power_log = "power_log";
  this->file_window = "window";
  this->file_power2d = "power2d";
  this->file_power2d_mk = "power2d_mk";
  // ************************************************************************
  this->kmax_y_ds = 0.2;
  this->kmax_bk = 0.2;
  this->kmin_bk = 0.1;
  this->use_fundamental_mode_as_kmin_bk = false;
  this->file_bispectrum = "bispectrum";
  this->Nft_random_collapse = 32;
  this->kmax_tracer_bias = 0.7;
  this->kmin_tracer_bias = 0.01;
  this->kmax_tracer_qbias = 0;
  this->kmin_tracer_qbias = 0;
  this->Healpix_resolution = 4;
 // ************************************************************************
 // ************************************************************************
 // Measurements of Abundance (Mass funcito, luminosity function)
  this->Get_Mstellar_function = false;
  this->Get_Luminosity_function = false;
  this->Get_Color_function = false;
  this->LF_estimator = "Vmax_o";
  this->Get_Color_Mag_plane = false;
  this->Get_Random_Catalog = false;
  this->Nbins_color = 1;
  this->Nbins_Mstellar = 1;
  this->Mstellar_min = 1;
  this->Mstellar_max = 1;
  this->Color_min = 0;
  this->Color_max = 1.0;
  // ************************************************************************
  // ************************************************************************
  this->MASS_units=1;
  this->NMASSbins=1;
  this->NMASSbins_mf=1;
  this->NMASSbins_power=1;
  this->LOGMASSmin=0;
  this->LOGMASSmax=16;
  this->VMAXmin = 0;
  this->VMAXmax = 1000;
  this->RSmin = 0;
  this->RSmax = 1000;
  this->SPINmin = 0;
  this->SPINmax = 1000;
  this->NPROPbins_bam=1;
  this->set_bins_equal_number_tracers = false;
  this->Number_of_bins_equal_number_tracers_main_property = 1;
  this->set_bins_equal_number_tracers_main_property = false;
  this->Number_of_bins_equal_number_tracers = 1;
  this->Nbins_hist = 10;

  this->Get_pearson_coefficient=false;
  this->Get_marked_power_spectrum =false;
  this->Get_power_spectrum =false;
  this->Get_cross_power_spectrum =false;
  this->Get_tracer_number_counts=false;
  this->Get_tidal_anisotropy_at_halo=false;
  this->Get_tracer_mass_field=false;
  this->Get_tracer_vmax_field=false;
  this->Get_tracer_spin_field=false;
  this->Get_tracer_spin_bullock_field=false;
  this->Get_tracer_local_mach_number=false;
  this->Get_tracer_local_dach_number=false;
  this->Get_tracer_local_dm_density=false;
  this->Get_spearman_coefficient=false;
  this->Get_cell_local_mach_number=false;
  this->Scale_mach_number = static_cast<real_prec>(2.0);
  this->NPROPbins_bam=1;
  this->lmax_bias=0;
  this->assign_bias_to_full_sample = false;
  this->Distance_fraction=1.0;
  // ************************************************************************
  // ************************************************************************
  this->get_distribution_min_separations=false;
  this->M_exclusion=1;
  this->use_vel_kernel=false;
  this->slength=10;
  this->slengthv=10;
  this->velbias=0;
  this->velbias_dm=0;
  this->velbias_random=0;
  this->dilute_dm_sample = false;
  this->fraction_dilute=0.1;
  this->masskernel_vel=0;
  this->iteration_ini = 0;
  this->get_distribution_min_separations=false;
  this->ic_alias_corrected = false;
  this->ic_input_type = DENSITY;
  this->velbias_random=0.0;
  this->velbias_dm =0.0;
  this->use_vel_kernel=false;
  this->vkernel_exponent=0.;
  this->masskernel=0;
  this->masskernel_vel=0; 
  this->N_lines_binary = 0;
  this->Number_of_GRF =1;
  this->Kmax_FA = 0.;
  this->Generate_FA=false;
  this->NVIRIALbins_power=0;

  this->Nft_random_collapse = 1;


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
  this->Xoffset = 0;
  this->Yoffset = 0;
  this->Zoffset = 0;
  this->xmax = 0;
  this->ymax = 0;
  this->zmax = 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Params::forbid_unknown_keys(const json& j,const std::set<std::string>& allowed,const std::string& path) 
{
  if (!j.is_object()) {
        throw std::runtime_error(path + " must be a JSON object");
    }
    for (const auto& [key, value] : j.items()) {
        if (allowed.find(key) == allowed.end()) {
            throw std::runtime_error("Unknown parameter in " + path + ": " + key);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// This method reads parameters from json and defines their default values. SO we do not need the init_pars method
void Params::read_pars_json(std::string file){

  
  std::ifstream rfile(file);
    if (!rfile) {
        throw std::runtime_error("Cannot open config file: " + file);
    }
  json cfg;
  rfile>>cfg;

 // -------------------------------------------------------------
  string pname="any";    
  string pname_c="any";    
  string description = "Any";
  string options = "any";
 // -------------------------------------------------------------

  if (!cfg.contains("Simulation"))
    throw std::runtime_error("Simulation section missing in config");

  std::string this_section = "Simulation";
  const auto& Simulation=cfg.at(this_section);
  std::string status=Simulation.value("status", "false");
  if(ON==status)
    {  
      
      this->enabled_sections.push_back(this_section);
      this->forbid_unknown_keys(Simulation, 
        {"status",
        "redshift", 
        "initial_redshift_sim",
        "initial_redshift_delta",
        "ic_index",
        "realization_index",
        "length_side_box"
        },
        this_section);

      // ____________________________________________________________________________
      string pname="redshift";    
      string pname_c="redshift";    
      string description = "Cosmological redshift";
      string options = "float > 0";
      this->redshift=Simulation.value(pname,0);
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,  this->redshift));
      // ____________________________________________________________________________
      pname="initial_redshift_sim";    
      pname_c="Initial_Redshift_SIM";    
      description="Initial redshift of a given N-body simulation";
      options="float > 0";
      this->Initial_Redshift_SIM=static_cast<real_prec>(Simulation.value(pname, 99.));
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,  this->Initial_Redshift_SIM));
      // ____________________________________________________________________________
      pname="initial_redshift_delta";    
      pname_c="Initial_Redshift_DELTA";    
      description="Redshift at which the initial linear power spectrum is normalized to genete IC";
      options="float > 0";
      this->Initial_Redshift_DELTA=Simulation.value(pname, 99.);
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,  this->Initial_Redshift_DELTA));
      // ____________________________________________________________________________
      pname="ic_index";    
      pname_c="IC_index";    
      description="Index to characterize the IC of some procedures";
      options="int > 0";
      this->IC_index=static_cast<int>(Simulation.value(pname, 0));
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,this->IC_index));
      // ____________________________________________________________________________
      pname="realization_index";    
      pname_c="realization"; 
      description="Index to characterize realizations when producing mocks";
      options="int > 0";
      this->realization=Simulation.value(pname, 0);
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,this->realization));
      // ____________________________________________________________________________
      pname="length_side_box";    
      pname_c="L_box"; 
      description="Length of the side of the box in unots of Mpc/h";
      options="float > 0";
      this->Lbox=Simulation.value(pname, 1000.0);
      this->collect_params_info(pname,pname_c,this_section, description,options);
      this->parameter_number.push_back(make_pair(pname,this->Lbox));
 
    
  }
 // -------------------------------------------------------------
 // -------------------------------------------------------------

   const auto& BiasAnalysis=cfg.at("BiasAnalysis");
   this_section = "BiasAnalysis";
   status=BiasAnalysis.value("status", "false");
   if(ON==status)
   {  
    this->enabled_sections.push_back("BiasAnalysis");
    this->forbid_unknown_keys(BiasAnalysis,
      {
        "status",
        "output_directory",
        "DarkMatter",
        "Tracers",
        "General",
        "BiasMappingTechnique"
      },
      "BiasAnalysis"); 
   // ____________________________________________________________________________
    pname="output_directory";    
    pname_c="Output_directory"; 
    description="Path to output directory for Bias Analysis";
    options="string";
    this->Output_directory = BiasAnalysis.value("output_directory", "null") ;
    this->collect_params_info(pname,pname_c,this_section, description,options);
    this->parameter_string.emplace_back(pname, this->Output_directory);
    // -------------------------------------------------------------
    const auto& DarkMatter = BiasAnalysis.at("DarkMatter");
    this_section="DarkMatter";
    this->forbid_unknown_keys(DarkMatter, 
        {
          "nbins",
          "mas" ,
          "mas_ref_pdf",
          "mas_new" ,
          "type_of_file",
          "name_property",
          "name_aux",
          "input_directory",
          "input_directory_ref",
          "input_directory_ref_two",
          "input_directory_new",
          "input_directory_new_ref",
          "catalogue_file_mesh", 
          "catalogue_file_mesh_new",
          "catalogue_file_mesh_new_ref",
          "catalogue_file_mesh_ref_pdf",
          "convert_density_to_delta",
          "velfieldx_file_mesh",
          "velfieldy_file_mesh",
          "velfieldz_file_mesh",
          "scale",
          "delta_min",
          "delta_max",
          "logdelta_min",
          "logdelta_max"
      },
      "DarkMatter");

// ____________________________________________________________________________
    pname="nbins";    
    pname_c="NX"; 
    description="Number of bins used to perform bias (distribution) analysis for the dark matter density --or overdenisty-- field.";
    options="ULONG";
    this->NX = static_cast<ULONG>(DarkMatter.value(pname, 10));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->NX);

    pname = "mas";
    pname_c = "iMAS_X";
    description = "Mass assignment scheme used to construct the density field.";
    options = "INT";
    this->iMAS_X = static_cast<int>(DarkMatter.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->iMAS_X);

    // mas_ref_pdf
    pname = "mas_ref_pdf";
    pname_c = "iMAS_X_REF_PDF";
    description = "Mass assignment scheme for the reference PDF density field.";
    options = "INT";
    this->iMAS_X_REF_PDF = static_cast<int>(DarkMatter.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->iMAS_X_REF_PDF);

    // mas_new
    pname = "mas_new";
    pname_c = "iMAS_X_NEW";
    description = "Mass assignment scheme for the new density field.";
    options = "INT";
    this->iMAS_X_NEW = static_cast<int>(DarkMatter.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->iMAS_X_NEW);

    // delta_min
    pname = "delta_min";
    pname_c = "delta_X_min";
    description = "Minimum overdensity value used in the analysis.";
    options = "real prec";
    this->delta_X_min = static_cast<real_prec>(DarkMatter.value(pname, 0.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->delta_X_min);

    // delta_max
    pname = "delta_max";
    pname_c = "delta_X_max";
    description = "Maximum overdensity value used in the analysis.";
    options = "real prec";
    this->delta_X_max = static_cast<real_prec>(DarkMatter.value("delta_max", 10.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->delta_X_max);

    // logdelta_min
    pname = "logdelta_min";
    pname_c = "ldelta_X_min";
    description = "Minimum log-overdensity value used in the analysis.";
    options = "real prec";
    this->ldelta_X_min = static_cast<real_prec>(DarkMatter.value("logdelta_min", 0.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ldelta_X_min);

    // logdelta_max
    pname = "logdelta_max";
    pname_c = "ldelta_X_max";
    description = "Maximum log-overdensity value used in the analysis.";
    options = "real prec";
    this->ldelta_X_max = static_cast<real_prec>(DarkMatter.value(pname, 10.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ldelta_X_max);

    // STRING parameters (no static_cast needed)
    pname = "type_of_file";
    pname_c = "Type_of_file_X";
    description = "Type of input file containing the density field.";
    options = "STRING";
    this->Type_of_file_X = DarkMatter.value(pname, "bin");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Type_of_file_X);

    pname = "name_property";
    pname_c = "Name_Property_X";
    description = "Name of the physical property stored in the density field.";
    options = "STRING";
    this->Name_Property_X = DarkMatter.value(pname, "DENSITY");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Property_X);

    pname = "name_aux";
    pname_c = "XNAME";
    description = "Auxiliary name used to identify the density field.";
    options = "STRING";
    this->XNAME = DarkMatter.value(pname, "DM");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->XNAME);

    // input directories
    pname = "input_directory";
    pname_c = "Input_Directory_X";
    description = "Input directory containing the density field files.";
    options = "STRING";
    this->Input_Directory_X = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X);

    pname = "input_directory_ref";
    pname_c = "Input_Directory_X_REF";
    description = "Input directory containing the reference density field.";
    options = "STRING";
    this->Input_Directory_X_REF = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X_REF);

    pname = "input_directory_ref_two";
    pname_c = "Input_Directory_X_REF_TWO";
    description = "Second input directory for reference density field.";
    options = "STRING";
    this->Input_Directory_X_REF_TWO = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X_REF_TWO);

    pname = "input_directory_new";
    pname_c = "Input_Directory_X_NEW";
    description = "Input directory containing the new density field.";
    options = "STRING";
    this->Input_Directory_X_NEW = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X_NEW);

    pname = "input_directory_new_ref";
    pname_c = "Input_Directory_X_new_ref";
    description = "Input directory for the reference of the new density field.";
    options = "STRING";
    this->Input_Directory_X_new_ref = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X_new_ref);

    // catalog files
    pname = "catalogue_file_mesh";
    pname_c = "Name_Catalog_X";
    description = "Catalog file associated with the density mesh.";
    options = "STRING";
    this->Name_Catalog_X = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Catalog_X);

    pname = "catalogue_file_mesh_new";
    pname_c = "Name_Catalog_X_NEW";
    description = "Catalog file associated with the new density mesh.";
    options = "STRING";
    this->Name_Catalog_X_NEW = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Catalog_X_NEW);

    pname = "catalogue_file_mesh_new_ref";
    pname_c = "Name_Catalog_X_new_ref";
    description = "Reference catalog file for the new density mesh.";
    options = "STRING";
    this->Name_Catalog_X_new_ref = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Catalog_X_new_ref);

    pname = "catalogue_file_mesh_ref_pdf";
    pname_c = "Name_Catalog_X_REF_PDF";
    description = "Catalog file used to compute the reference PDF.";
    options = "STRING";
    this->Name_Catalog_X_REF_PDF = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Catalog_X_REF_PDF);

    // velocity fields
    pname = "velfieldx_file_mesh";
    pname_c = "Name_VelFieldx_X";
    description = "Velocity field X-component mesh file.";
    options = "STRING";
    this->Name_VelFieldx_X = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_VelFieldx_X);

    pname = "velfieldy_file_mesh";
    pname_c = "Name_VelFieldy_X";
    description = "Velocity field Y-component mesh file.";
    options = "STRING";
    this->Name_VelFieldy_X = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_VelFieldy_X);

    pname = "velfieldz_file_mesh";
    pname_c = "Name_VelFieldz_X";
    description = "Velocity field Z-component mesh file.";
    options = "STRING";
    this->Name_VelFieldz_X = DarkMatter.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_VelFieldz_X);

    // scale
    pname = "scale";
    pname_c = "Scale_X";
    description = "Scaling used for the density field (e.g. log or linear).";
    options = "STRING";
    this->Scale_X = DarkMatter.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Scale_X);

    // convert_density_to_delta
    pname = "convert_density_to_delta";
    pname_c = "Convert_Density_to_Delta_X";
    description = "Convert density field into overdensity field.";
    options = "BOOL";
    this->Convert_Density_to_Delta_X = DarkMatter.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Convert_Density_to_Delta_X);

    // -------------------------------------------------------------
    // -------------------------------------------------------------
    const auto& Tracers = BiasAnalysis.at("Tracers");
    this_section="Tracers";
    this->forbid_unknown_keys(Tracers, 
        {
        "nbins",
        "mas", 
        "type_of_file",
        "name_property",
        "name_aux",
        "input_directory",
        "input_directory_two",
        "catalogue_file_mesh",
        "catalogue_file_hr_mesh",
        "convert_density_to_delta",
        "scale",
        "delta_min",
        "delta_max" ,
        "logdelta_min",
        "logdelta_max",
        "nbins_mass",
        "nbins_sat_frac"
        },
      "Tracers");

          // nbins
      pname = "nbins";
      pname_c = "NY";
      description = "Number of bins used to perform bias (distribution) analysis for the tracer density --or overdensity-- field.";
      options = "ULONG";
      this->NY = static_cast<ULONG>(Tracers.value("nbins", 10));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NY);

      // mas
      pname = "mas";
      pname_c = "iMAS_Y";
      description = "Mass assignment scheme used to construct the tracer density field.";
      options = "INT";
      this->iMAS_Y = static_cast<int>(Tracers.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->iMAS_Y);

      // delta_min
      pname = "delta_min";
      pname_c = "delta_Y_min";
      description = "Minimum overdensity value used in the tracer analysis.";
      options = "real prec";
      this->delta_Y_min = static_cast<real_prec>(Tracers.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->delta_Y_min);

      // delta_max
      pname = "delta_max";
      pname_c = "delta_Y_max";
      description = "Maximum overdensity value used in the tracer analysis.";
      options = "real prec";
      this->delta_Y_max = static_cast<real_prec>(Tracers.value(pname, 10));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->delta_Y_max);

      // logdelta_min
      pname = "logdelta_min";
      pname_c = "ldelta_Y_min";
      description = "Minimum log-overdensity value used in the tracer analysis.";
      options = "real prec";
      this->ldelta_Y_min = static_cast<real_prec>(Tracers.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->ldelta_Y_min);

      // logdelta_max
      pname = "logdelta_max";
      pname_c = "ldelta_Y_max";
      description = "Maximum log-overdensity value used in the tracer analysis.";
      options = "real prec";
      this->ldelta_Y_max = static_cast<real_prec>(Tracers.value(pname, 10));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->ldelta_Y_max);

      // nbins_mass
      pname = "nbins_mass";
      pname_c = "NY_MASS";
      description = "Number of mass bins used for tracer analysis.";
      options = "ULONG";
      this->NY_MASS = static_cast<ULONG>(Tracers.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NY_MASS);

      // nbins_sat_frac
      pname = "nbins_sat_frac";
      pname_c = "NY_SAT_FRAC";
      description = "Number of bins used for satellite fraction analysis.";
      options = "ULONG";
      this->NY_SAT_FRAC = static_cast<ULONG>(Tracers.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NY_SAT_FRAC);

      // type_of_file
      pname = "type_of_file";
      pname_c = "Type_of_file_Y";
      description = "Type of input file containing the tracer density field.";
      options = "STRING";
      this->Type_of_file_Y = Tracers.value(pname, "bin");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Type_of_file_Y);

      // name_property
      pname = "name_property";
      pname_c = "Name_Property_Y";
      description = "Name of the physical property stored in the tracer field.";
      options = "STRING";
      this->Name_Property_Y = Tracers.value(pname, "DENSITY");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Name_Property_Y);

      // name_aux
      pname = "name_aux";
      pname_c = "YNAME";
      description = "Auxiliary name used to identify the tracer field.";
      options = "STRING";
      this->YNAME = Tracers.value(pname, "TR");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->YNAME);

      // input_directory
      pname = "input_directory";
      pname_c = "Input_Directory_Y";
      description = "Input directory containing the tracer density field files.";
      options = "STRING";
      this->Input_Directory_Y = Tracers.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Input_Directory_Y);

      // input_directory_two
      pname = "input_directory_two";
      pname_c = "Input_Directory_Y_TWO";
      description = "Second input directory containing tracer data.";
      options = "STRING";
      this->Input_Directory_Y_TWO = Tracers.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Input_Directory_Y_TWO);

      // catalog_file_mesh
      pname = "catalogue_file_mesh";
      pname_c = "Name_Catalog_Y";
      description = "Catalog file associated with the tracer density mesh.";
      options = "STRING";
      this->Name_Catalog_Y = Tracers.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Name_Catalog_Y);

      // catalog_file_hr_mesh
      pname = "catalogue_file_hr_mesh";
      pname_c = "Name_Catalog_Y_HR";
      description = "High-resolution catalog file associated with the tracer mesh.";
      options = "STRING";
      this->Name_Catalog_Y_HR = Tracers.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Name_Catalog_Y_HR);

      // scale
      pname = "scale";
      pname_c = "Scale_Y";
      description = "Scaling used for the tracer density field (e.g. log or linear).";
      options = "STRING";
      this->Scale_Y = Tracers.value(pname, "log");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Scale_Y);

      pname = "convert_density_to_delta";
      pname_c = "Convert_Density_to_Delta_Y";
      description = "Convert tracer density field into overdensity field.";
      options = "BOOL";
      this->Convert_Density_to_Delta_Y = Tracers.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Convert_Density_to_Delta_Y);

    // -------------------------------------------------------------
    // -------------------------------------------------------------
    const auto& BiasGeneral = BiasAnalysis.at("General");
    this_section="BiasGeneral";
    this->forbid_unknown_keys(BiasGeneral, 
        {
        "ngrid_ft",
        "redefine_limits",
        "write_pdf_number_counts",
        "compute_joint_pdf",
        "compute_conditional_pdf",
        "write_files_for_histograms"
        },
      "BiasGeneral");
    // redefine_limits

    pname = "ngrid_ft";
    pname_c = "Nft";
    description = "Grid resolution used for Fourier transforms.";
    options = "INT";
    this->Nft = BiasGeneral.value(pname, 256);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft);


    pname = "redefine_limits";
    pname_c = "Redefine_limits";
    description = "Redefine the limits of the density or overdensity range used in the analysis.";
    options = "BOOL";
    this->Redefine_limits = BiasGeneral.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Redefine_limits);

    // write_pdf_number_counts
    pname = "write_pdf_number_counts";
    pname_c = "Write_PDF_number_counts";
    description = "Write PDF number counts to output files.";
    options = "BOOL";
    this->Write_PDF_number_counts = BiasGeneral.value("write_pdf_number_counts", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Write_PDF_number_counts);

    // compute_joint_pdf
    pname = "compute_joint_pdf";
    pname_c = "Comp_joint_PDF";
    description = "Compute the joint probability distribution function.";
    options = "BOOL";
    this->Comp_joint_PDF = BiasGeneral.value("compute_joint_pdf", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Comp_joint_PDF);

    // compute_conditional_pdf
    pname = "compute_conditional_pdf";
    pname_c = "Comp_conditional_PDF";
    description = "Compute the conditional probability distribution function.";
    options = "BOOL";
    this->Comp_conditional_PDF = BiasGeneral.value("compute_conditional_pdf", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Comp_conditional_PDF);

    // write_files_for_histograms
    pname = "write_files_for_histograms";
    pname_c = "write_files_for_histograms";
    description = "Write auxiliary files used to construct histograms.";
    options = "BOOL";
    this->write_files_for_histograms = BiasGeneral.value("write_files_for_histograms", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->write_files_for_histograms);

    // -------------------------------------------------------------
    // -------------------------------------------------------------
    this_section="BiasMappingTechnique";
    const auto& BiasMappingTechnique = BiasAnalysis.at("BiasMappingTechnique");

    this->forbid_unknown_keys(BiasMappingTechnique, 
        {
        "number_iterations_kernel" ,
        "initial_iteration_kernel_average",
        "index_initial_iteration",
        "number_dark_matter_realizations_mock_generation",
        "number_iterations_preprocessing_dm",
        "write_outputs_at_iterations",
        "apply_rank_ordering",
        "apply_rank_ordering_ab_initio",
        "number_chunks_new_dm",
        "input_directory_BMT_kernel",
        "input_directory_BMT_kernel_second_calibration",
        "levels_main_property_multiscale_assignment",
        "mesh_resolutions_multiscale_assignment",
        "tolerance_property_multiscale_assignment",
        "mesh_size_closest_neighbour_collapse_towards_randoms",
        "distance_fraction_to_nearest_random_collapse",
        "mass_scale_exclusion",
        "use_velkernel_velocity_assignment",
        "scale_length_velkernel",
        "velocity_bias_random",
        "velocity_bias_dm"
        },
      "BiasMappingTechnique");
// number_iterations_kernel
    pname = "number_iterations_kernel";
    pname_c = "N_iterations_Kernel";
    description = "Number of iterations used in the bias-mapping kernel.";
    options = "ULONG";
    this->N_iterations_Kernel =
        static_cast<ULONG>(BiasMappingTechnique.value("number_iterations_kernel", 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_iterations_Kernel);

    // initial_iteration_kernel_average
    pname = "initial_iteration_kernel_average";
    pname_c = "Iteration_Kernel_average";
    description = "Initial iteration used to compute the kernel average.";
    options = "ULONG";
    this->Iteration_Kernel_average =
        static_cast<ULONG>(BiasMappingTechnique.value("initial_iteration_kernel_average", 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Iteration_Kernel_average);

    // index_initial_iteration
    pname = "index_initial_iteration";
    pname_c = "iteration_ini";
    description = "Index of the initial iteration for the bias-mapping procedure.";
    options = "ULONG";
    this->iteration_ini =
        static_cast<ULONG>(BiasMappingTechnique.value("index_initial_iteration", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->iteration_ini);

    // number_dark_matter_realizations_mock_generation
    pname = "number_dark_matter_realizations_mock_generation";
    pname_c = "N_dm_realizations";
    description = "Number of dark matter realizations used for mock generation.";
    options = "ULONG";
    this->N_dm_realizations =
        static_cast<ULONG>(BiasMappingTechnique.value(
            "number_dark_matter_realizations_mock_generation", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_dm_realizations);

    // number_iterations_preprocessing_dm
    pname = "number_iterations_preprocessing_dm";
    pname_c = "N_iterations_dm";
    description = "Number of preprocessing iterations applied to the dark matter field.";
    options = "ULONG";
    this->N_iterations_dm =
        static_cast<ULONG>(BiasMappingTechnique.value(
            "number_iterations_preprocessing_dm", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_iterations_dm);

    // n_dm_initial_mock
    pname = "n_dm_initial_mock";
    pname_c = "N_dm_initial";
    description = "Initial number of dark matter realizations.";
    options = "ULONG";
    this->N_dm_initial =
        static_cast<ULONG>(BiasMappingTechnique.value("n_dm_initial_mock", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_dm_initial);

    // write_outputs_at_iterations
    pname = "write_outputs_at_iterations";
    pname_c = "output_at_iteration";
    description = "List of iterations at which outputs are written.";
    options = "VECTOR<ULONG>";
    this->output_at_iteration =
        BiasMappingTechnique.value("write_outputs_at_iterations",
                                  std::vector<ULONG>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_uvectors.emplace_back(pname, this->output_at_iteration);

    // levels_main_property_multiscale_assignment
    pname = "levels_main_property_multiscale_assignment";
    pname_c = "list_Props_Threshold_MultiLevels";
    description = "Threshold levels of the main property for multiscale assignment.";
    options = "VECTOR<REAL>";
    this->list_Props_Threshold_MultiLevels =
        BiasMappingTechnique.value("levels_main_property_multiscale_assignment",
            std::vector<real_prec>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_vectors.emplace_back(
        pname, this->list_Props_Threshold_MultiLevels);

    // mesh_resolutions_multiscale_assignment
    pname = "mesh_resolutions_multiscale_assignment";
    pname_c = "list_Nft_MultiLevels";
    description = "Mesh resolutions used in multiscale assignment.";
    options = "VECTOR<INT>";
    this->list_Nft_MultiLevels =
        BiasMappingTechnique.value("mesh_resolutions_multiscale_assignment", std::vector<int>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_ivectors.emplace_back(
        pname, this->list_Nft_MultiLevels);

    // tolerance_property_multiscale_assignment
    pname = "tolerance_property_multiscale_assignment";
    pname_c = "list_Props_Tolerance_MultiLevels";
    description = "Tolerance levels used in multiscale property assignment.";
    options = "VECTOR<REAL>";
    this->list_Props_Tolerance_MultiLevels =
        BiasMappingTechnique.value("tolerance_property_multiscale_assignment",
           std::vector<real_prec>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_vectors.emplace_back(
        pname, this->list_Props_Tolerance_MultiLevels);

        // apply_rank_ordering
    pname = "apply_rank_ordering";
    pname_c = "Apply_Rankordering";
    description = "Apply rank ordering during bias mapping.";
    options = "BOOL";
    this->Apply_Rankordering =
        BiasMappingTechnique.value("apply_rank_ordering", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Apply_Rankordering);

    // apply_rank_ordering_ab_initio
    pname = "apply_rank_ordering_ab_initio";
    pname_c = "Apply_Rankordering_ab_initio";
    description = "Apply rank ordering from the initial iteration.";
    options = "BOOL";
    this->Apply_Rankordering_ab_initio =
        BiasMappingTechnique.value("apply_rank_ordering_ab_initio", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(
        pname, this->Apply_Rankordering_ab_initio);

    // use_velkernel_velocity_assignment
    pname = "use_velkernel_velocity_assignment";
    pname_c = "use_vel_kernel";
    description = "Enable velocity kernel during velocity assignment.";
    options = "BOOL";
    this->use_vel_kernel =
        BiasMappingTechnique.value("use_velkernel_velocity_assignment", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_vel_kernel);

    // number_chunks_new_dm
    pname = "number_chunks_new_dm";
    pname_c = "Number_of_chunks_new_dm";
    description = "Number of chunks used to process the new dark matter field.";
    options = "ULONG";
    this->Number_of_chunks_new_dm =
        static_cast<ULONG>(BiasMappingTechnique.value("number_chunks_new_dm", 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Number_of_chunks_new_dm);

    // mesh_size_closest_neighbour_collapse_towards_randoms
    pname = "mesh_size_closest_neighbour_collapse_towards_randoms";
    pname_c = "Nft_random_collapse";
    description = "Mesh size used for collapse towards nearest randoms.";
    options = "ULONG";
    this->Nft_random_collapse =
        static_cast<ULONG>(BiasMappingTechnique.value(
            "mesh_size_closest_neighbour_collapse_towards_randoms", 16));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft_random_collapse);

    // distance_fraction_to_nearest_random_collapse
    pname = "distance_fraction_to_nearest_random_collapse";
    pname_c = "Distance_fraction";
    description = "Fractional distance to the nearest random used to collapse dark matter tracers towards.";
    options = "real prec";
    this->Distance_fraction =
        static_cast<real_prec>(BiasMappingTechnique.value(
            "distance_fraction_to_nearest_random_collapse", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Distance_fraction);

    // mass_scale_exclusion
    pname = "mass_scale_exclusion";
    pname_c = "M_exclusion";
    description = "Mass scale used to exclude objects during bias mapping.";
    options = "real prec";
    this->M_exclusion =
        static_cast<real_prec>(BiasMappingTechnique.value(
            "mass_scale_exclusion", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->M_exclusion);

    // scale_lenght_velkernel
    pname = "scale_length_velkernel";
    pname_c = "slengthv";
    description = "Scale length of the velocity kernel.";
    options = "real prec";
    this->slengthv =
        static_cast<real_prec>(BiasMappingTechnique.value(
            "scale_lenght_velkernel", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->slengthv);

    // velocity_bias_random
    pname = "velocity_bias_random";
    pname_c = "velbias_random";
    description = "Velocity bias applied to random particles.";
    options = "real prec";
    this->velbias_random =
        static_cast<real_prec>(BiasMappingTechnique.value(
            "velocity_bias_random", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->velbias_random);

    // velocity_bias_dm
    pname = "velocity_bias_dm";
    pname_c = "velbias_dm";
    description = "Velocity bias applied to dark matter particles.";
    options = "real prec";
    this->velbias_dm =
        static_cast<real_prec>(BiasMappingTechnique.value(
            "velocity_bias_dm", 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->velbias_dm);


    
  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------
 
  const auto& CWCAnalysis=cfg.at("CWCAnalysis");
  status = CWCAnalysis.value("status", "off"); // default = "off"
  this_section="CWCAnalysis";
   if(ON==status)
  {  
    this->enabled_sections.push_back(this_section);
    this->forbid_unknown_keys(CWCAnalysis, 
        {
          "status", 
          "output_directory",
          "box_size",
          "min_x_coord",
          "min_y_coord",
          "min_z_coord",
          "ngrid_ft",
          "number_thresholds_cosmicweb_types",
          "lambda_threshold",
          "lambda_threshold_vshear",
          "nbins_mass_sknots",
          "nbins_mass_vknot",
          "cosmic_web_types_tidal" , 
          "cosmic_web_types_velshear"
        },
      "CWCAnalysis");


    // output_directory
    pname = "output_directory";
    pname_c = "Output_directory";
    description = "Output directory where CWC analysis results are written.";
    options = "STRING";
    this->Output_directory = CWCAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Output_directory);

    // box_size
    pname = "box_size";
    pname_c = "Lbox";
    description = "Size of the simulation box.";
    options = "real prec";
    this->Lbox = CWCAnalysis.value(pname, 1.);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Lbox);

    // min_x_coord
    pname = "min_x_coord";
    pname_c = "xmin";
    description = "Minimum X coordinate of the analysis volume.";
    options = "real prec";
    this->xmin = CWCAnalysis.value(pname, 0.0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->xmin);

    // min_y_coord
    pname = "min_y_coord";
    pname_c = "ymin";
    description = "Minimum Y coordinate of the analysis volume.";
    options = "real prec";
    this->ymin = CWCAnalysis.value(pname, 0.0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ymin);

    // min_z_coord
    pname = "min_z_coord";
    pname_c = "zmin";
    description = "Minimum Z coordinate of the analysis volume.";
    options = "real prec";
    this->zmin = CWCAnalysis.value(pname, 0.0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->zmin);

    // ngrid_ft
    pname = "ngrid_ft";
    pname_c = "Nft";
    description = "Grid resolution used for Fourier transforms.";
    options = "INT";
    this->Nft = CWCAnalysis.value(pname, 256);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft);

    // number_thresholds_cosmicweb_types
    pname = "number_thresholds_cosmicweb_types";
    pname_c = "Nlambdath";
    description = "Number of lambda thresholds used to define cosmic web types.";
    options = "INT";
    this->Nlambdath =
        CWCAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nlambdath);

    // lambda_threshold
    pname = "lambda_threshold";
    pname_c = "lambdath";
    description = "Lambda threshold used for tidal cosmic web classification.";
    options = "real prec";
    this->lambdath = CWCAnalysis.value(pname, 0.0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->lambdath);

    // lambda_threshold_vshear
    pname = "lambda_threshold_vshear";
    pname_c = "lambdath_v";
    description = "Lambda threshold used for velocity shear cosmic web classification.";
    options = "real prec";
    this->lambdath_v =
        CWCAnalysis.value(pname, 0.);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->lambdath_v);

    // nbins_mass_sknots
    pname = "nbins_mass_sknots";
    pname_c = "n_sknot_massbin";
    description = "Number of mass bins for saddle–knot objects.";
    options = "INT";
    this->n_sknot_massbin =
        CWCAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->n_sknot_massbin);

    // nbins_mass_vknots
    pname = "nbins_mass_vknots";
    pname_c = "n_vknot_massbin";
    description = "Number of mass bins for velocity-knot objects.";
    options = "INT";
    this->n_vknot_massbin =
        CWCAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->n_vknot_massbin);

    // cosmic_web_types_tidal
    pname = "cosmic_web_types_tidal";
    pname_c = "cwt_used";
    description = "List of tidal cosmic web types included in the analysis.";
    options = "VECTOR<ULONG>";
    this->cwt_used =
        CWCAnalysis.value(pname,
                          std::vector<ULONG>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_uvectors.emplace_back(pname, this->cwt_used);

    // cosmic_web_types_vshear
    pname = "cosmic_web_types_vshear";
    pname_c = "cwv_used";
    description = "List of velocity shear cosmic web types included in the analysis.";
    options = "VECTOR<ULONG>";
    this->cwv_used =
        CWCAnalysis.value(pname,
                          std::vector<ULONG>{});
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_uvectors.emplace_back(pname, this->cwv_used);




  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  const auto& FourierAnalysis=cfg.at("FourierAnalysis");
  status = FourierAnalysis.value("status", "off");
  this_section = "FourierAnalysis";
  if(ON== status)
  {  
    this_section = "FourierAnalysis";
    this->enabled_sections.push_back(this_section);
    this->forbid_unknown_keys(FourierAnalysis, 
        {
          "status", 
          "output_directory",
          "statistics",
          "name_survey",
          "type_of_object",
          "input_type",
          "input_type_two",
          "ntracers_for_delta_input" ,
          "number_of_random_files",
          "delta_grid_file" ,
          "delta_grid_file2",
          "delta_grid_file3",
          "delta_grid_file4",
          "measure_cross_power",
          "measure_cross_power_using_a",
          "measure_cross_power_using_b",
          "get_new_line_of_sight",
          "clustering_space",
          "box_size",
          "box_size_lowres",
          "get_new_box_size",
          "min_x_coord",
          "min_y_coord",
          "min_z_coord",
          "ngrid_ft",
          "ngrid_ft_highres",
          "ngrid_ft_lowres",   
          "ngrid_jk",   
          "output_interpolated_field",
          "show_interpolated_field",
          "mass_assignment_scheme",
          "use_correction_mass_assignment_scheme",
          "type_of_kshell_binning",
          "number_log_kbins",
          "delta_kmin",
          "ndel_kbinning_data",
          "ndel_kbinning_window",
          "use_fkp_weight",
          "estimated_power",
          "number_mu_bins",
          "use_shot_noise_correction",
          "compute_fkp_variance",
          "compute_fkp_variance_exact",
          "output_file_power", 
          "show_power_spectrum",
          "output_file_power_log", 
          "output_file_window",
          "output_file_power2d",
          "output_file_power2d_mk",
          "kmax_yamamoto_direct_sum",
          "kmax_bispectrum",
          "kmin_bispectrum",
          "use_fundamental_mode_as_kmin_bk",
          "output_file_bispectrum"
        },
      this_section);

    pname = "output_directory";
    pname_c = "Output_directory";
    description = "Output directory where Fourier analysis results are written.";
    options = "STRING";
    this->Output_directory =
        FourierAnalysis.value("output_directory", "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Output_directory);

    // statistics
    pname = "statistics";
    pname_c = "statistics";
    description = "Type of Fourier-space statistic to compute (e.g. power spectrum).";
    options = "STRING. Denotes the estimator of poẁer used for two or three point stats: Pk_fkp (FKP), Pk_ys (Yamamoto with Scoccimarro implementation): Pk_yb (Yamamoto with Bianchi's iomplementation); Pk_ds (direct sum)";
    this->statistics =
        FourierAnalysis.value(pname, "Pk_fkp");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->statistics);

    // name_survey
    pname = "name_survey";
    pname_c = "Name_survey";
    description = "Name of the survey associated with the Fourier analysis. To append to output files";
    options = "STRING";
    this->Name_survey =
        FourierAnalysis.value(pname, "Zorbas");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_survey);

    // type_of_object
    pname = "type_of_object";
    pname_c = "type_of_object";
    description = "Type of object used in the Fourier analysis (e.g. TRACER, DM, TRACER_REF, TRACER_MOCK, TRACER_MOCK_ONLY_COORDS, TRACER_REF_ONLY_COORDS).";
    options = "STRING";
    this->type_of_object =
        FourierAnalysis.value(pname, "TRACER");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->type_of_object);


    if(this->type_of_object!="TRACER" && this->type_of_object!="TRACER_REF" && this->type_of_object!="TRACER_MOCK" && this->type_of_object!="TRACER_MOCK_ONLY_COORDS" && this->type_of_object!="TRACER_REF_ONLY_COORDS")
     {
        cout<<RED<<endl;
        cout<<"Parameter type_of_object = "<<this->type_of_object<<endl;
        cout<<"does not match any of the expected options"<<endl;
        cout<<"\tTRACER"<<endl;
        cout<<"\tTRACER_REF"<<endl;
        cout<<"\tTRACER_MOCK"<<endl;
        cout<<"\tTRACER_MOCK_ONLY_COORDS"<<endl;
        cout<<"\tTRACER_REF_ONLY_COORDS"<<endl;
        cout<<RESET<<endl;
        throw std::invalid_argument("Wrong input in parameter file");
    }


        // ===============================
    // INPUT TYPES
    // ===============================
    pname = "input_type";
    pname_c = "input_type";
    description = "Type of input used for the Fourier analysis (e.g. catalog or density field interpolated on a grid).";
    options = "STRING (catalog, density_grid, overdensity_grid)";
    this->input_type = FourierAnalysis.value(pname, "catalog");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->input_type);

    pname = "input_type_two";
    pname_c = "input_type_two";
    description = "Secondary input type used for cross-power spectrum calculations.";
    options = "STRING";
    this->input_type_two = FourierAnalysis.value(pname, "catalog");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->input_type_two);

    // ===============================
    // INPUT COUNTS AND FILES
    // ===============================
    pname = "ntracers_for_delta_input";
    pname_c = "ngal_delta";
    description = "Number of tracers used to construct the input density contrast field.";
    options = "ULONG";
    this->ngal_delta = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ngal_delta);

    pname = "number_of_random_files";
    pname_c = "Number_of_random_files";
    description = "Number of random catalogue files used in the analysis.";
    options = "ULONG";
    this->Number_of_random_files = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Number_of_random_files);

    pname = "delta_grid_file";
    pname_c = "delta_grid_file";
    description = "Filename containing the input density grid.";
    options = "STRING";
    this->delta_grid_file = FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->delta_grid_file);

    pname = "delta_grid_file2";
    pname_c = "delta_grid_file2";
    description = "Second density grid file used for cross-power calculations.";
    options = "STRING";
    this->delta_grid_file2 = FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->delta_grid_file2);

    pname = "delta_grid_file3";
    pname_c = "delta_grid_file3";
    description = "Third density grid file used in multi-field analyses.";
    options = "STRING";
    this->delta_grid_file3 = FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->delta_grid_file3);

    pname = "delta_grid_file4";
    pname_c = "delta_grid_file4";
    description = "Fourth density grid file used in multi-field analyses.";
    options = "STRING";
    this->delta_grid_file4 = FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->delta_grid_file4);

    // ===============================
    // CROSS POWER OPTIONS
    // ===============================
    pname = "measure_cross_power";
    pname_c = "measure_cross";
    description = "Enable computation of cross power spectra.";
    options = "BOOL";
    this->measure_cross = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->measure_cross);

    pname = "measure_cross_power_using_a";
    pname_c = "measure_cross_from_1";
    description = "Index of the first field used in cross-power computation.";
    options = "INT";
    this->measure_cross_from_1 =
        static_cast<int>(FourierAnalysis.value(pname, 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->measure_cross_from_1);

    pname = "measure_cross_power_using_b";
    pname_c = "measure_cross_from_2";
    description = "Index of the second field used in cross-power computation.";
    options = "INT";
    this->measure_cross_from_2 =
        static_cast<int>(FourierAnalysis.value(pname, 2));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->measure_cross_from_2);

    // ===============================
    // GEOMETRY AND GRID
    // ===============================
    pname = "get_new_line_of_sight";
    pname_c = "new_los";
    description = "Recompute the line of sight for each Fourier mode.";
    options = "BOOL";
    this->new_los = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->new_los);

    pname = "clustering_space";
    pname_c = "clustering_space";
    description = "Space where the statistics is to be computed. Applies only for simulations";
    options = "STRING, real_space or redshift_space";
    this->clustering_space = FourierAnalysis.value(pname, "real_space");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->clustering_space);

    pname = "box_size_lowres";
    pname_c = "Lbox_low";
    description = "Box size used for low-resolution Fourier grids, in Mpc/h";
    options = "real prec";
    this->Lbox_low =  static_cast<real_prec>(FourierAnalysis.value(pname, 1.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Lbox_low);


    pname = "box_size";
    pname_c = "Lbox";
    description = "Box size used for Fourier grids, in Mpc/h";
    options = "real prec";
    this->Lbox =  static_cast<real_prec>(FourierAnalysis.value(pname, 1.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Lbox);


    pname = "get_new_box_size";
    pname_c = "new_Lbox";
    description = "Enable recomputation of the simulation box size.";
    options = "BOOL";
    this->new_Lbox = static_cast<real_prec>(FourierAnalysis.value(pname, false));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->new_Lbox);

    pname = "min_x_coord";
    pname_c = "xmin";
    description = "Minimum x-coordinate in Mpc/h";
    options = "BOOL";
    this->xmin = static_cast<real_prec>(FourierAnalysis.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->xmin);

    pname = "min_y_coord";
    pname_c = "ymin";
    description = "Minimum x-coordinate in Mpc/h";
    options = "BOOL";
    this->ymin = static_cast<real_prec>(FourierAnalysis.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ymin);


    pname = "min_z_coord";
    pname_c = "zmin";
    description = "Minimum x-coordinate in Mpc/h";
    options = "BOOL";
    this->zmin = static_cast<real_prec>(FourierAnalysis.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->zmin);


    // ===============================
    // FFT GRIDS
    // ===============================
    pname = "ngrid_ft";
    pname_c = "Nft";
    description = "Number of grid cells per dimension for Fourier transforms.";
    options = "ULONG";
    this->Nft = FourierAnalysis.value(pname, 128);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft);

    pname = "ngrid_ft_highres";
    pname_c = "Nft_HR";
    description = "Number of grid cells per dimension for high-resolution FFTs.";
    options = "ULONG";
    this->Nft_HR = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft_HR);

    pname = "ngrid_ft_lowres";
    pname_c = "Nft_low";
    description = "Number of grid cells per dimension for low-resolution FFTs.";
    options = "ULONG";
    this->Nft_low = FourierAnalysis.value("ngrid_ft_lowres", 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft_low);

    pname = "ngrid_jk";
    pname_c = "Nft_JK";
    description = "Number of grid cells used for jackknife resampling.";
    options = "ULONG";
    this->Nft_JK = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Nft_JK);

    // ===============================
    // MASS ASSIGNMENT
    // ===============================
    pname = "mass_assignment_scheme";
    pname_c = "mass_assignment_scheme";
    description = "Mass assignment scheme used to deposit particles onto the grid.";
    options = "STRING";
    this->mass_assignment_scheme =
        FourierAnalysis.value(pname, "NGP");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->mass_assignment_scheme);

    pname = "output_interpolated_field";
    pname_c = "output_interpolated_field";
    description = "Output in binary filw the interpolated field";
    options = "BOOL";
    this->output_interpolated_field =
        FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->output_interpolated_field);

    this->set_output_file_interpolated_field(this->Output_directory+ "interpolated_feld");


    pname = "show_interpolated_field";
    pname_c = "show_interpolated_field";
    description = "Show a slice of the interpolated field of the tracers used on the power spectrum measurements";
    options = "BOOL";
    this->show_interpolated_field =
        FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_interpolated_field);

    pname = "show_power_spectrum";
    pname_c = "show_power_spectrum";
    description = "Show measurements of power spectrum";
    options = "BOOL";
    this->show_power_spectrum =
        FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_power_spectrum);



    pname = "use_correction_mass_assignment_scheme";
    pname_c = "MAS_correction";
    description = "Apply correction for the mass assignment window function.";
    options = "BOOL";
    this->MAS_correction =
        FourierAnalysis.value("use_correction_mass_assignment_scheme", true);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->MAS_correction);

    // ===============================
    // BINNING AND WEIGHTS
    // ===============================
    pname = "type_of_kshell_binning";
    pname_c = "type_of_binning";
    description = "Type of k-shell binning used in Fourier space.";
    options = "STRING";
    this->type_of_binning =
        FourierAnalysis.value("type_of_kshell_binning", "linear");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->type_of_binning);

    pname = "number_log_kbins";
    pname_c = "N_log_bins";
    description = "Number of logarithmic k-bins.";
    options = "ULONG";
    this->N_log_bins = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_log_bins);

    pname = "delta_kmin";
    pname_c = "N_log_bins";
    description = "This parameter f_{k} determines the minimum value of Fourier wavenumber in order to define spherical shells. Minimum wavenumber is f * Delta, where Delta is the fundamental mode. See documentation for further details";
    options = "ULONG";
    this->DeltaKmin = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->DeltaKmin);


    pname = "ndel_kbinning_data";
    pname_c = "ndel_data";
    description = "This parameter nd determines the width of the spherical shells in Fourier space. If the fundamental mode is Dk, the width of spherical shells is dk = nd * Dk such that the arrays defining the spherical shells are defined as ki = kmin + (i + 0.5)*dk.";
    options = "ULONG";
    this->ndel_data = static_cast<real_prec>(FourierAnalysis.value(pname, 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ndel_data);

    pname = "ndel_kbinning_window",
    pname_c = "ndel_window";
    description = "This parameter nd determines the width of the spherical shells in Fourier space. If the fundamental mode is Dk, the width of spherical shells is dk = nd * Dk such that the arrays defining the spherical shells are defined as ki = kmin + (i + 0.5)*dk.";
    options = "ULONG";
    this->ndel_window = static_cast<real_prec>(FourierAnalysis.value(pname, 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->ndel_window);


    pname = "use_fkp_weight";
    pname_c = "FKP_weight";
    description = "Enable FKP weighting scheme.";
    options = "BOOL";
    this->FKP_weight = FourierAnalysis.value("use_fkp_weight", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->FKP_weight);


    pname = "estimated_power";
    pname_c = "Pest";
    description = "Estimated power for FKP weights.";
    options = "float";
    this->Pest = static_cast<real_prec>(FourierAnalysis.value(pname, 1e5));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Pest);


    pname = "number_mu_bins";
    pname_c = "N_mu_bins";
    description = "Number of logarithmic k-bins.";
    options = "ULONG";
    this->N_mu_bins = FourierAnalysis.value(pname, 1);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_mu_bins);

    pname = "use_shot_noise_correction";
    pname_c = "SN_correction";
    description = "Correct for Poisson shot-noise. Set this to false when dark matter is to be analyzed, or when the cross correlation function is to be measured.";
    options = "BOOL";
    this->SN_correction = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->SN_correction);


    pname = "compute_fkp_variance";
    pname_c = "FKP_error_bars";
    description = "Compute FKP error bars from the FKP variance";
    options = "BOOL";
    this->FKP_error_bars = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->FKP_error_bars);


    pname = "compute_fkp_variance_exact";
    pname_c = "FKP_error_bars_exact";
    description = "Compute FKP error bars from the FKP variance";
    options = "BOOL";
    this->FKP_error_bars_exact = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->FKP_error_bars_exact);

    // ===============================
    // OUTPUT FILES
    // ===============================
    pname = "output_file_power";
    pname_c = "file_power";
    description = "Output file for the power spectrum.";
    options = "STRING";
    this->file_power =
        FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_power);

    pname = "output_file_power_log";
    pname_c = "file_power";
    description = "Output file for the logs of power spectrum.";
    options = "STRING";
    this->file_power_log =
        FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_power_log);


    pname = "output_file_window";
    pname_c = "file_window";
    description = "Output file for the power spectrum.";
    options = "STRING";
    this->file_window =
        FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_window);


    pname = "output_file_power2d";
    pname_c = "file_power2d";
    description = "Output file for the 2d power spectrum.";
    options = "STRING";
    this->file_power2d =
        FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_power2d);


    pname = "output_file_power2d_mk";
    pname_c = "file_power2d_mk";
    description = "Output file for two dimensional power spectrum in polar coordinates.";
    options = "STRING";
    this->file_power2d_mk =
        FourierAnalysis.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_power2d_mk);


    pname = "output_file_bispectrum";
    pname_c = "file_bispectrum";
    description = "Output file for the bispectrum.";
    options = "STRING";
    this->file_bispectrum =
        FourierAnalysis.value("output_file_bispectrum", "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_bispectrum);

    //Other

    pname = "kmax_yamamoto_direct_sum";
    pname_c = "kmax_y_ds";
    description = "Maximum wavenumber for the direct sum approach to Yamamoto-Blake";
    options = "float";
    this->kmax_y_ds = static_cast<real_prec>(FourierAnalysis.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->kmax_y_ds);

    pname = "kmax_bispectrum";
    pname_c = "kmax_bk";
    description = "Maximum k-value for constructing k-bins for the Bispectrum estimator";
    options = "float";
    this->kmax_bk = static_cast<real_prec>(FourierAnalysis.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->kmax_bk);

    pname = "kmin_bispectrum";
    pname_c = "kmin_bk";
    description = "Minimum k-value for constructing k-bins for the Bispectrum estimator";
    options = "float";
    this->kmin_bk = static_cast<real_prec>(FourierAnalysis.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->kmin_bk);

    pname = "use_fundamental_mode_as_kmin_bk";
    pname_c = "use_fundamental_mode_as_kmin_bk";
    description = "These parameters is used to define the shells in k-space for bispectrum  measurements.";
    options = "BOOL";
    this->use_fundamental_mode_as_kmin_bk = FourierAnalysis.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_fundamental_mode_as_kmin_bk);

  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------

  const auto& DarkMatterCatalogue=cfg.at("DarkMatterCatalogue");
  status = DarkMatterCatalogue.value("status", "off"); 
  this_section = "DarkMatterCatalogue";
  if(ON==status)
  {
    this->enabled_sections.push_back("DarkMatterCatalogue");
    this->forbid_unknown_keys(DarkMatterCatalogue, 
        {
          "status", 
          "input_directory",
          "catalogue_file_mesh",
          "file_bin_x_coord",
          "file_bin_y_coord",
          "file_bin_z_coord",
          "number_of_tracers_in_binary",
          "use_low_pass_filter",
          "system_of_coordinates",
          "i_coord1",
          "i_coord2",
          "i_coord3",
          "i_v1",
          "i_v2",
          "i_v3",
          "i_mass"
        },
      "DarkMatterCatalogue");

    // ===============================
    // BINARY COORDINATE FILES
    // ===============================

// input directories
    pname = "input_directory";
    pname_c = "Input_Directory_X";
    description = "Input directory containing the density field files.";
    options = "STRING";
    this->Input_Directory_X = DarkMatterCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_Directory_X);

// catalog files
    pname = "catalogue_file_mesh";
    pname_c = "Name_Catalog_X";
    description = "Catalogue file associated with the density mesh.";
    options = "STRING";
    this->Name_Catalog_X =DarkMatterCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Name_Catalog_X);

    pname = "file_bin_x_coord";
    pname_c = "file_bin_x_coord";
    description = "Binary file containing the x-coordinate of dark matter tracers.";
    options = "STRING";
    this->file_bin_x_coord =
        DarkMatterCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_bin_x_coord);

    pname = "file_bin_y_coord";
    pname_c = "file_bin_y_coord";
    description = "Binary file containing the y-coordinate of dark matter tracers.";
    options = "STRING";
    this->file_bin_y_coord =
        DarkMatterCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_bin_y_coord);

    pname = "file_bin_z_coord";
    pname_c = "file_bin_z_coord";
    description = "Binary file containing the z-coordinate of dark matter tracers.";
    options = "STRING";
    this->file_bin_z_coord =
        DarkMatterCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_bin_z_coord);

    // ===============================
    // BINARY INPUT OPTIONS
    // ===============================
    pname = "number_of_tracers_in_binary";
    pname_c = "N_lines_binary";
    description = "Total number of dark matter tracers stored in the binary catalogue.";
    options = "ULONG";
    this->N_lines_binary =
        static_cast<ULONG>(
            DarkMatterCatalogue.value("number_of_tracers_in_binary", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_lines_binary);

    pname = "use_low_pass_filter";
    pname_c = "use_low_pass_filter";
    description = "Apply a low-pass filter to the dark matter density field.";
    options = "BOOL";
    this->use_low_pass_filter =
        DarkMatterCatalogue.value("use_low_pass_filter", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_low_pass_filter);

    // ===============================
    // COORDINATE SYSTEM AND COLUMN INDICES
    // ===============================
    pname = "system_of_coordinates";
    pname_c = "sys_of_coord_dm";
    description = "Coordinate system identifier used for the dark matter catalogue.";
    options = "INT";
    this->sys_of_coord_dm =
        static_cast<int>(
            DarkMatterCatalogue.value("system_of_coordinates", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->sys_of_coord_dm);

    pname = "i_coord1";
    pname_c = "i_coord1_dm";
    description = "Column index corresponding to the first spatial coordinate.";
    options = "INT";
    this->i_coord1_dm =
        static_cast<int>(
            DarkMatterCatalogue.value("i_coord1", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord1_dm);

    pname = "i_coord2";
    pname_c = "i_coord2_dm";
    description = "Column index corresponding to the second spatial coordinate.";
    options = "INT";
    this->i_coord2_dm =
        static_cast<int>(
            DarkMatterCatalogue.value("i_coord2", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord2_dm);

    pname = "i_coord3";
    pname_c = "i_coord3_dm";
    description = "Column index corresponding to the third spatial coordinate.";
    options = "INT";
    this->i_coord3_dm =
        static_cast<int>(
            DarkMatterCatalogue.value("i_coord3", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord3_dm);

    pname = "i_mass";
    pname_c = "i_mass_dm";
    description = "Column index corresponding to the dark matter particle mass.";
    options = "INT";
    this->i_mass_dm =
        static_cast<int>(
            DarkMatterCatalogue.value("i_mass", 0)
        );
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mass_dm);

  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  const auto& TracerCatalogue=cfg.at("TracerCatalogue");
  status = TracerCatalogue.value("status", "off");
  this_section = "TracerCatalogue";

  if(ON==status)
  {
   this->enabled_sections.push_back("TracerCatalogue");
   this->forbid_unknown_keys(TracerCatalogue, 
        {
          "status", 
          "input_directory",
          "input_directory_catalogue_new_ref",
          "catalogue_file",
          "catalogue_file_new_ref", 
          "system_of_coordinates",
          "catalogue_file_mesh",
          "redshift_space_coordinates_included",
          "angles_units",
          "velocity_units",
          "i_coord1",
          "i_coord2",
          "i_coord3",
          "i_v1",
          "i_v2",
          "i_v3",
          "i_mass",
          "i_vmax",
          "i_vrms",
          "i_rs",
          "i_spin",
          "i_spin_bullock",
          "i_virial",
          "i_rvir",
          "i_b_to_a",
          "i_c_to_a",
          "i_star_formation",
          "i_color",
          "i_stellar_mass",
          "i_mean_density",
          "i_abs_mag",
          "i_app_mag",
          "i_weight1",
          "use_weight1",
          "i_weight2",
          "use_weight2",
          "i_weight3",
          "use_weight3",
          "i_weight4",
          "use_weight4",
          "weight_with_mass",
          "weight_vel_with_mass"
        },
      "TracerCatalogue");

    // ===============================
    // INPUT DIRECTORIES AND FILES
    // ===============================
    pname = "input_directory";
    pname_c = "Input_dir_cat";
    description = "Input directory containing the tracer catalogue.";
    options = "STRING";
    this->Input_dir_cat =
        TracerCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_dir_cat);

    pname = "input_directory_catalogue_new_ref";
    pname_c = "Input_dir_cat_TWO";
    description = "Input directory for the reference tracer catalogue.";
    options = "STRING";
    this->Input_dir_cat_TWO =
        TracerCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_dir_cat_TWO);

    pname = "catalogue_file";
    pname_c = "file_catalogue";
    description = "Filename of the tracer catalogue.";
    options = "STRING";
    this->file_catalogue =
        TracerCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_catalogue);

    pname = "catalogue_file_new_ref";
    pname_c = "file_catalogue_new_ref";
    description = "Filename of a reference tracer catalogue. Used for BMT";
    options = "STRING";
    this->file_catalogue_new_ref =
        TracerCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_catalogue_new_ref);

    // catalog_file_mesh
      pname = "catalogue_file_mesh";
      pname_c = "Name_Catalog_Y";
      description = "Catalog file associated with the tracer density interpolated on a mesh.";
      options = "STRING";
      this->Name_Catalog_Y = TracerCatalogue.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Name_Catalog_Y);

    // ===============================
    // COORDINATE SYSTEM AND UNITS
    // ===============================
    pname = "system_of_coordinates";
    pname_c = "sys_of_coord_g";
    description = "Coordinate system identifier for tracer catalogue.";
    options = "INT";
    this->sys_of_coord_g =
        TracerCatalogue.value("system_of_coordinates", 0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->sys_of_coord_g);

    pname = "redshift_space_coordinates_included";
    pname_c = "redshift_space_coords_g";
    description = "Flag indicating whether redshift-space coordinates are included.";
    options = "BOOL. Set to false when reading a simulation. If reading coordinates from a galaxy redshift survey, RS are alrady in: set it to true. If power from simulations is to be measured in redshift space (which implies that velocities are avaiblae, otherwise the code ignores the petition), enable to the pre-proc directive REDSHIFT_SPACE";
    this->redshift_space_coords_g =
        TracerCatalogue.value("redshift_space_coordinates_included", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->redshift_space_coords_g);

    pname = "angles_units";
    pname_c = "angles_units_g";
    description = "Units used for angular coordinates (e.g. degrees or radians).";
    options = "STRING";
    this->angles_units_g =
        TracerCatalogue.value("angles_units", "D");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->angles_units_g);

    pname = "velocity_units";
    pname_c = "vel_units_g";
    description = "Units used for velocity components.";
    options = "STRING";
    this->vel_units_g =
        TracerCatalogue.value("velocity_units", "kmps");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->vel_units_g);

    // ===============================
    // COLUMN INDICES: COORDINATES & VELOCITIES
    // ===============================
    pname = "i_coord1";
    pname_c = "i_coord1_g";
    description = "Column index for the first spatial coordinate.";
    options = "INT";
    this->i_coord1_g =
        static_cast<int>(TracerCatalogue.value("i_coord1", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord1_g);

    pname = "i_coord2";
    pname_c = "i_coord2_g";
    description = "Column index for the second spatial coordinate.";
    options = "INT";
    this->i_coord2_g =
        static_cast<int>(TracerCatalogue.value("i_coord2", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord2_g);

    pname = "i_coord3";
    pname_c = "i_coord3_g";
    description = "Column index for the third spatial coordinate.";
    options = "INT";
    this->i_coord3_g =
        static_cast<int>(TracerCatalogue.value("i_coord3", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord3_g);

    pname = "i_v1";
    pname_c = "i_v1_g";
    description = "Column index for the first velocity component.";
    options = "INT";
    this->i_v1_g =
        static_cast<int>(TracerCatalogue.value("i_v1", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_v1_g);

    pname = "i_v2";
    pname_c = "i_v2_g";
    description = "Column index for the second velocity component.";
    options = "INT";
    this->i_v2_g =
        static_cast<int>(TracerCatalogue.value("i_v2", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_v2_g);

    pname = "i_v3";
    pname_c = "i_v3_g";
    description = "Column index for the third velocity component.";
    options = "INT";
    this->i_v3_g =
        static_cast<int>(TracerCatalogue.value("i_v3", -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_v3_g);

    if(this->clustering_space=="redshift_space" && (this->i_v1_g<0 || this->i_v2_g<0 || this->i_v3_g<0)){
      this->clustering_space="real_space";
      cout<<RED<<"Warning. SImulation requested in redshift space, but no velocities are available"<<endl;
      exit(1);
    }


    // ===============================
    // PHYSICAL AND GALAXY PROPERTIES
    // ===============================
    #define ADD_INT_PARAM(key, member, desc) \
      pname = key; pname_c = #member; description = desc; options = "INT"; \
      this->member = static_cast<int>(TracerCatalogue.value(key, -1)); \
      this->collect_params_info(pname, pname_c, this_section, description, options); \
      this->parameter_number.emplace_back(pname, this->member);

    ADD_INT_PARAM("i_mass", i_mass_g, "Column index for tracer mass.")
    ADD_INT_PARAM("i_vmax", i_vmax_g, "Column index for maximum circular velocity.")
    ADD_INT_PARAM("i_vrms", i_vrms_g, "Column index for RMS velocity.")
    ADD_INT_PARAM("i_rs", i_rs_g, "Column index for scale radius.")
    ADD_INT_PARAM("i_spin", i_spin_g, "Column index for spin parameter.")
    ADD_INT_PARAM("i_spin_bullock", i_spin_bullock_g, "Column index for Bullock spin.")
    ADD_INT_PARAM("i_virial", i_virial_g, "Column index for virial mass.")
    ADD_INT_PARAM("i_rvir", i_rvir_g, "Column index for virial radius.")
    ADD_INT_PARAM("i_b_to_a", i_b_to_a_g, "Column index for b/a axis ratio.")
    ADD_INT_PARAM("i_c_to_a", i_c_to_a_g, "Column index for c/a axis ratio.")
    ADD_INT_PARAM("i_star_formation", i_sf_g, "Column index for star formation rate.")
    ADD_INT_PARAM("i_color", i_color_g, "Column index for color.")
    ADD_INT_PARAM("i_stellar_mass", i_stellar_mass_g, "Column index for stellar mass.")
    ADD_INT_PARAM("i_mean_density", i_mean_density_g, "Column index for mean density.")
    ADD_INT_PARAM("i_abs_mag", i_abs_mag_g, "Column index for absolute magnitude.")
    ADD_INT_PARAM("i_app_mag", i_app_mag_g, "Column index for apparent magnitude.")
    ADD_INT_PARAM("i_weight1", i_weight1_g, "Column index for first weight.")
    ADD_INT_PARAM("i_weight2", i_weight2_g, "Column index for second weight.")
    ADD_INT_PARAM("i_weight3", i_weight3_g, "Column index for third weight.")
    ADD_INT_PARAM("i_weight4", i_weight4_g, "Column index for fourth weight.")

    #undef ADD_INT_PARAM

    // ===============================
    // WEIGHT FLAGS
    // ===============================
    pname = "use_weight1";
    pname_c = "use_weight1_g";
    description = "Enable first tracer weight.";
    options = "BOOL";
    this->use_weight1_g =
        TracerCatalogue.value("use_weight1", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight1_g);

    pname = "use_weight2";
    pname_c = "use_weight2_g";
    description = "Enable second tracer weight.";
    options = "BOOL";
    this->use_weight2_g =
        TracerCatalogue.value("use_weight2", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight2_g);

    pname = "use_weight3";
    pname_c = "use_weight3_g";
    description = "Enable third tracer weight.";
    options = "BOOL";
    this->use_weight3_g =
        TracerCatalogue.value("use_weight3", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight3_g);

    pname = "use_weight4";
    pname_c = "use_weight4_g";
    description = "Enable fourth tracer weight.";
    options = "BOOL";
    this->use_weight4_g =
        TracerCatalogue.value("use_weight4", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight4_g);

    pname = "weight_with_mass";
    pname_c = "weight_with_mass";
    description = "Weight tracers by mass.";
    options = "BOOL";
    this->weight_with_mass =
        TracerCatalogue.value("weight_with_mass", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->weight_with_mass);

    pname = "weight_vel_with_mass";
    pname_c = "weight_vel_with_mass";
    description = "Weight tracer velocities by mass.";
    options = "BOOL";
    this->weight_vel_with_mass =
        TracerCatalogue.value("weight_vel_with_mass", false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->weight_vel_with_mass);

  }
  // -------------------------------------------------------------
  // -------------------------------------------------------------
  const auto& RandomCatalogue=cfg.at("RandomCatalogue");
  this_section="RandomCatalogue";
  status=RandomCatalogue.value("status", "off");
  if(ON==status)
  {
    this->enabled_sections.push_back("RandomCatalogue");
    this->forbid_unknown_keys(RandomCatalogue, 
        {
          "status", 
          "use_random_catalogue",
          "input_directory",
          "random_catalogue_file", 
          "nbar_tabulated",
          "use_external_file_nbar",
          "file_nbar",
          "sys_of_coord",
          "angles_units",
          "i_coord1",
          "i_coord2",
          "i_coord3",
          "i_mass",
          "i_color",
          "i_stellar_mass",
          "i_app_mag",
          "i_abs_mag",
          "i_mean_density",
          "i_weight1",
          "use_weight1",
          "i_weight2",
          "use_weight2",
          "i_weight3",
          "use_weight3",
          "i_weight4",
          "use_weight4",
          "compute_window_matrix"
        },
      "RandomCatalogue");


    // use_random_catalog
    pname = "use_random_catalogue";
    pname_c = "use_random_catalogue";
    description = "Use a random catalogue for the analysis.";
    options = "BOOL, true or false; if set to false, the code interprets that the input catalog is a simulation.";
    this->use_random_catalog = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_random_catalog);

    pname = "input_directory";
    pname_c = "Input_dir_cat_random";
    description = "Input directory containing the random catalogue.";
    options = "STRING";
    this->Input_dir_cat_random =
        RandomCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Input_dir_cat_random);

    
    // file_random
    pname = "random_catalogue_file";
    pname_c = "file_random";
    description = "Path to the random catalogue file.";
    options = "STRING";
    this->file_random = RandomCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_random);

    // nbar_tabulated
    pname = "nbar_tabulated";
    pname_c = "nbar_tabulated";
    description = "Whether the number density is tabulated in bins.";
    options = "BOOL";
    this->nbar_tabulated = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->nbar_tabulated);

    // use_file_nbar
    pname = "use_external_file_nbar";
    pname_c = "use_file_nbar";
    description = "Use an external file for the number density nbar.";
    options = "BOOL";
    this->use_file_nbar = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_file_nbar);

    // file_nbar
    pname = "file_nbar";
    pname_c = "file_nbar";
    description = "Path to the external nbar file.";
    options = "STRING";
    this->file_nbar = RandomCatalogue.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_nbar);

    // angles_units_r
    pname = "angle_units";
    pname_c = "angles_units_r";
    description = "Units for angles (D for degrees).";
    options = "STRING";
    this->angles_units_r = RandomCatalogue.value(pname, "D");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->angles_units_r);

    // i_coord1_r
    pname = "i_coord1";
    pname_c = "i_coord1_r";
    description = "Column index for X coordinates in the random catalogue.";
    options = "INT";
    this->i_coord1_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord1_r);

    // i_coord2_r
    pname = "i_coord2";
    pname_c = "i_coord2_r";
    description = "Column index for Y coordinates in the random catalogue.";
    options = "INT";
    this->i_coord2_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord2_r);

    // i_coord3_r
    pname = "i_coord3";
    pname_c = "i_coord3_r";
    description = "Column index for Z coordinates in the random catalogue.";
    options = "INT";
    this->i_coord3_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_coord3_r);

    // i_mass_r
    pname = "i_mass";
    pname_c = "i_mass_r";
    description = "Column index for mass in the random catalogue.";
    options = "INT";
    this->i_mass_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mass_r);

    // i_color_r
    pname = "i_color";
    pname_c = "i_color_r";
    description = "Column index for color in the random catalogue.";
    options = "INT";
    this->i_color_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_color_r);

    // i_stellar_mass_r
    pname = "i_stellar_mass";
    pname_c = "i_stellar_mass_r";
    description = "Column index for stellar mass in the random catalogue.";
    options = "INT";
    this->i_stellar_mass_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_stellar_mass_r);

    // i_mean_density_r
    pname = "i_mean_density";
    pname_c = "i_mean_density_r";
    description = "Column index for mean density in the random catalogue.";
    options = "INT";
    this->i_mean_density_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mean_density_r);

    // i_abs_mag_r
    pname = "i_abs_mag";
    pname_c = "i_abs_mag_r";
    description = "Column index for absolute magnitude in the random catalogue.";
    options = "INT";
    this->i_abs_mag_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_abs_mag_r);

    // i_app_mag_r
    pname = "i_app_mag";
    pname_c = "i_app_mag_r";
    description = "Column index for apparent magnitude in the random catalogue.";
    options = "INT";
    this->i_app_mag_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_app_mag_r);

    // i_weight1_r
    pname = "i_weight1";
    pname_c = "i_weight1_r";
    description = "Column index for weight1 in the random catalogue.";
    options = "INT";
    this->i_weight1_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_weight1_r);

    // use_weight1_r
    pname = "use_weight1";
    pname_c = "use_weight1_r";
    description = "Whether to use weight1 in the random catalogue.";
    options = "BOOL";
    this->use_weight1_r = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight1_r);

    // i_weight2_r
    pname = "i_weight2";
    pname_c = "i_weight2_r";
    description = "Column index for weight2 in the random catalogue.";
    options = "INT";
    this->i_weight2_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_weight2_r);

    // use_weight2_r
    pname = "use_weight2";
    pname_c = "use_weight2_r";
    description = "Whether to use weight2 in the random catalogue.";
    options = "BOOL";
    this->use_weight2_r = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight2_r);

    // i_weight3_r
    pname = "i_weight3";
    pname_c = "i_weight3_r";
    description = "Column index for weight3 in the random catalogue.";
    options = "INT";
    this->i_weight3_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_weight3_r);

    // use_weight3_r
    pname = "use_weight3";
    pname_c = "use_weight3_r";
    description = "Whether to use weight3 in the random catalogue.";
    options = "BOOL";
    this->use_weight3_r = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight3_r);

    // i_weight4_r
    pname = "i_weight4";
    pname_c = "i_weight4_r";
    description = "Column index for weight4 in the random catalogue.";
    options = "INT";
    this->i_weight4_r = static_cast<int>(RandomCatalogue.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_weight4_r);

    // use_weight4_r
    pname = "use_weight4";
    pname_c = "use_weight4_r";
    description = "Whether to use weight4 in the random catalogue.";
    options = "BOOL";
    this->use_weight4_r = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_weight4_r);

    // get_window_matrix
    pname = "compute_window_matrix";
    pname_c = "get_window_matrix";
    description = "Compute the window matrix for the random catalogue.";
    options = "BOOL";
    this->get_window_matrix = RandomCatalogue.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->get_window_matrix);

  }

  // -------------------------------------------------------------
  // -------------------------------------------------------------
  const auto& SurveyProperties=cfg.at("SurveyProperties");
  status=SurveyProperties.value("status", "off");
  if(ON==status)
  {
    this->enabled_sections.push_back("SurveyProperties");

     this->forbid_unknown_keys(SurveyProperties, 
        {
          "status", 
          "constant_sky_depth",
          "redshift_min_sample",
          "redshift_max_sample",
          "area_survey_sky",
          "MaskProperties",
          "RedshiftHistograms"
        },
      "SurveyProperties");
    
    // constant_depth
    pname = "constant_sky_depth";
    pname_c = "constant_depth";
    description = "Whether the survey has a constant sky depth.";
    options = "BOOL";
    this->constant_depth = SurveyProperties.value(pname, true);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->constant_depth);

    // redshift_min_sample
    pname = "redshift_min_sample";
    pname_c = "redshift_min_sample";
    description = "Minimum redshift in the sample.";
    options = "real prec";
    this->redshift_min_sample = static_cast<real_prec>(SurveyProperties.value(pname, 0.001));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->redshift_min_sample);

    // redshift_max_sample
    pname = "redshift_max_sample";
    pname_c = "redshift_max_sample";
    description = "Maximum redshift in the sample.";
    options = "real prec";
    this->redshift_max_sample = static_cast<real_prec>(SurveyProperties.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->redshift_max_sample);

    // area_survey
    pname = "area_surveyed_sky";
    pname_c = "area_survey";
    description = "Total area of the surveyed sky in square degrees.";
    options = "real prec";
    this->area_survey = static_cast<real_prec>(SurveyProperties.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->area_survey);

  // -------------------------------------------------------------

   const auto& MaskProperties=SurveyProperties.at("MaskProperties");

   this->forbid_unknown_keys(MaskProperties, 
        {
          "i_mask_pixel",
          "i_mask_alpha",
          "i_mask_delta",
          "i_mask_flag"},
    "MaskProperties");

    // i_mask_pixel
    pname = "i_mask_pixel";
    pname_c = "i_mask_pixel";
    description = "Index of the mask pixel.";
    options = "real prec";
    this->i_mask_pixel = static_cast<real_prec>(MaskProperties.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mask_pixel);

    // i_mask_alpha
    pname = "i_mask_alpha";
    pname_c = "i_mask_alpha";
    description = "Index of the mask alpha coordinate.";
    options = "real prec";
    this->i_mask_alpha = static_cast<real_prec>(MaskProperties.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mask_alpha);

    // i_mask_delta
    pname = "i_mask_delta";
    pname_c = "i_mask_delta";
    description = "Index of the mask delta coordinate.";
    options = "real prec";
    this->i_mask_delta = static_cast<real_prec>(MaskProperties.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mask_delta);

    // i_mask_flag
    pname = "i_mask_flag";
    pname_c = "i_mask_flag";
    description = "Index of the mask flag.";
    options = "real prec";
    this->i_mask_flag = static_cast<real_prec>(MaskProperties.value(pname, -1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->i_mask_flag);

  // -------------------------------------------------------------
    const auto& RedshiftHistograms=SurveyProperties.at("RedshiftHistograms");

    this->forbid_unknown_keys(RedshiftHistograms, 
          {
            "number_bins_dndz",
            "new_number_bins_dndz",
            "output_file_dndz" 
          },
      "RedshiftHistograms");

        // N_dndz_bins
    pname = "number_bins_dndz";
    pname_c = "N_dndz_bins";
    description = "Number of bins for the dN/dz histogram.";
    options = "real prec";
    this->N_dndz_bins = static_cast<real_prec>(RedshiftHistograms.value(pname, 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->N_dndz_bins);

    // new_N_dndz_bins
    pname = "new_number_bins_dndz";
    pname_c = "new_N_dndz_bins";
    description = "New number of bins for the dN/dz histogram.";
    options = "real prec";
    this->new_N_dndz_bins = static_cast<real_prec>(RedshiftHistograms.value(pname, 1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->new_N_dndz_bins);

    // file_dndz
    pname = "output_file_dndz";
    pname_c = "file_dndz";
    description = "Output file for the dN/dz histogram.";
    options = "STRING";
    this->file_dndz = RedshiftHistograms.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_dndz);

  }

  //---------------------------------------------------------------------------------------------
 //---------------------------------------------------------------------------------------------
  const auto& GalaxyAnalysis=cfg.at("GalaxyAnalysis");
  status = GalaxyAnalysis.value("status", "off");
  if(ON==status)
  {
     this->enabled_sections.push_back("GalaxyAnalysis");
     this->forbid_unknown_keys(GalaxyAnalysis, 
          {
            "status", 
            "output_directory",
            "measure_stellar_mass_function" ,
            "measure_luminosity_function" ,
            "measure_color_distribution",
            "luminosity_function_estimator",
            "compute_color_magnitude plane",
            "generate_random_catalogue",
            "number_bins_color",
            "min_color",
            "max_color",
            "number_bins_stellar_mass",
            "min_stellar_mass",
            "max_stellar_mass"
          },
      "GalaxyAnalysis");

      // Output_directory
      pname = "output_directory";
      pname_c = "Output_directory";
      description = "Directory for galaxy analysis outputs.";
      options = "STRING";
      this->Output_directory = GalaxyAnalysis.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Output_directory);

      // Get_Mstellar_function
      pname = "measure_stellar_mass_function";
      pname_c = "Get_Mstellar_function";
      description = "Compute the stellar mass function for galaxies.";
      options = "BOOL";
      this->Get_Mstellar_function = GalaxyAnalysis.value(pname, false);
      if (this->i_stellar_mass_g < 0 && this->Get_Mstellar_function) {
          std::cout << "Warning: stellar mass is not present in input galaxy catalog" << std::endl;
          this->Get_Mstellar_function = false;
      }
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_Mstellar_function);

      // Get_Luminosity_function
      pname = "measure_luminosity_function";
      pname_c = "Get_Luminosity_function";
      description = "Compute the galaxy luminosity function.";
      options = "BOOL";
      this->Get_Luminosity_function = GalaxyAnalysis.value(pname, false);
      if ((this->i_app_mag_g < 0 || this->i_abs_mag_g < 0) && this->Get_Luminosity_function) {
          std::cout << "Warning: luminosity or magnitude is not present in input galaxy catalog" << std::endl;
          this->Get_Luminosity_function = false;
      }
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_Luminosity_function);

      // Get_Color_function
      pname = "measure_color_distribution";
      pname_c = "Get_Color_function";
      description = "Compute the galaxy color distribution.";
      options = "BOOL";
      this->Get_Color_function = GalaxyAnalysis.value(pname, false);
      if (this->i_color_g < 0 && this->Get_Color_function) {
          std::cout << "Warning: color is not present in input galaxy catalog" << std::endl;
          this->Get_Color_function = false;
      }
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_Color_function);

      // LF_estimator
      pname = "luminosity_function_estimator";
      pname_c = "LF_estimator";
      description = "Estimator to use for the luminosity function.";
      options = "STRING";
      this->LF_estimator = GalaxyAnalysis.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->LF_estimator);

      // Get_Color_Mag_plane
      pname = "compute_color_magnitude_plane";
      pname_c = "Get_Color_Mag_plane";
      description = "Compute the color-magnitude plane.";
      options = "BOOL";
      this->Get_Color_Mag_plane = GalaxyAnalysis.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_Color_Mag_plane);

      // Get_Random_Catalog
      pname = "generate_random_catalogue";
      pname_c = "Get_Random_Catalog";
      description = "Generate a random galaxy catalog.";
      options = "BOOL";
      this->Get_Random_Catalog = GalaxyAnalysis.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_Random_Catalog);

      // Nbins_color
      pname = "number_bins_color";
      pname_c = "Nbins_color";
      description = "Number of bins for galaxy colors.";
      options = "ULONG";
      this->Nbins_color = static_cast<ULONG>(GalaxyAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Nbins_color);

      // Color_min
      pname = "min_color";
      pname_c = "Color_min";
      description = "Minimum color value.";
      options = "real prec";
      this->Color_min = static_cast<real_prec>(GalaxyAnalysis.value(pname, 0.1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Color_min);

      // Color_max
      pname = "max_color";
      pname_c = "Color_max";
      description = "Maximum color value.";
      options = "real prec";
      this->Color_max = static_cast<real_prec>(GalaxyAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Color_max);

      // Nbins_Mstellar
      pname = "number_bins_stellar_mass";
      pname_c = "Nbins_Mstellar";
      description = "Number of bins for stellar mass.";
      options = "ULONG";
      this->Nbins_Mstellar = static_cast<ULONG>(GalaxyAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Nbins_Mstellar);

      // Mstellar_min
      pname = "min_stellar_mass";
      pname_c = "Mstellar_min";
      description = "Minimum stellar mass.";
      options = "real prec";
      this->Mstellar_min = static_cast<real_prec>(GalaxyAnalysis.value(pname, 0.1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Mstellar_min);

      // Mstellar_max
      pname = "max_stellar_mass";
      pname_c = "Mstellar_max";
      description = "Maximum stellar mass.";
      options = "real prec";
      this->Mstellar_max = static_cast<real_prec>(GalaxyAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Mstellar_max);
  }  
 //---------------------------------------------------------------------------------------------
 //---------------------------------------------------------------------------------------------
  const auto& HaloAnalysis=cfg.at("HaloAnalysis");
  status=HaloAnalysis.value("status", "off");
  this_section = "HaloAnalysis";
  if(ON==status)
  {
    this->enabled_sections.push_back("HaloAnalysis");
    this->forbid_unknown_keys(HaloAnalysis, 
          {
            "status", 
            "output_directory",
            "mass_units",
            "number_bins_histograms",
            "number_mass_bins",
            "number_mass_bins_mass_function",
            "number_mass_bins_power_spectrum",
            "min_logmass",
            "max_logmass",
            "min_vmax",
            "max_vmax",
            "min_rs",
            "max_rs",
            "min_spin",
            "max_spin",
            "number_sec_property_bins_assignment_bmt",
            "set_bins_equal_number_tracers_main_property",
            "number_of_bins_equal_number_tracers_main_property", 
            "set_bins_equal_number_tracers",
            "number_of_bins_equal_number_tracers",
            "mass_cuts",
            "min_massbins", 
            "max_massbins",
            "min_vmaxbins",
            "max_vmaxbins",
            "min_rsbins",
            "max_rsbins" ,
            "min_spinbins",
            "max_spinbins",
            "min_vrmsbins", 
            "max_vrmsbins" ,
            "min_virialbins",
            "max_virialbins", 
            "min_btoabins", 
            "max_btoabins", 
            "min_ctoabins", 
            "max_ctoabins", 
            "min_machbins", 
            "max_machbins",
            "min_dachbins", 
            "max_dachbins", 
            "min_biasbins", 
            "max_biasbins", 
            "min_lcbins", 
            "max_lcbins", 
            "min_tabins", 
            "max_tabins", 
            "min_phbins",
            "max_phbins",
            "use_real_and_redshift_space",
            "compute_pearson_correlation_coefficient",
            "compute_spearman_coefficient",
            "perform_principal_component_analysis",
            "measure_marked_power_spectrum",
            "measure_power_spectrum" ,
            "measure_cross_power_spectrum",
            "get_tracer_number_counts",
            "get_tracer_mass_field",
            "get_tracer_vmax_field",
            "get_tracer_spin_field",
            "get_tracer_prop_function",
            "get_tracer_prop_function_cwt",
            "get_cell_local_mach_number",
            "get_tracer_local_mach_number",
            "get_tracer_local_dach_number",
            "get_tracer_local_dm_density",
            "scale_mach_number",
            "get_local_overdensity",
            "get_tidal_anisotropy_at_halo",
            "get_peak_height_at_halo",
            "get_distribution_minimum_separations",
            "IndividualTracerBias"
          },
      "HaloAnalysis");

      // Output_directory
      pname = "output_directory";
      pname_c = "Output_directory";
      description = "Directory for halo analysis outputs.";
      options = "STRING";
      this->Output_directory = HaloAnalysis.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Output_directory);

      // MASS_units
      pname = "mass_units";
      pname_c = "MASS_units";
      description = "Units for halo mass.";
      options = "real prec";
      this->MASS_units = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->MASS_units);

      // Nbins_hist
      pname = "number_bins_histograms";
      pname_c = "Nbins_hist";
      description = "Number of bins for histograms.";
      options = "ULONG";
      this->Nbins_hist = static_cast<ULONG>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Nbins_hist);

      // NMASSbins
      pname = "number_mass_bins";
      pname_c = "NMASSbins";
      description = "Number of mass bins for halos.";
      options = "ULONG";
      this->NMASSbins = static_cast<ULONG>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NMASSbins);

      // NMASSbins_mf
      pname = "number_mass_bins_mass_function";
      pname_c = "NMASSbins_mf";
      description = "Number of mass bins for the mass function.";
      options = "ULONG";
      this->NMASSbins_mf = static_cast<ULONG>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NMASSbins_mf);

      // NMASSbins_power
      pname = "number_mass_bins_power_spectrum";
      pname_c = "NMASSbins_power";
      description = "Number of mass bins for the power spectrum.";
      options = "ULONG";
      this->NMASSbins_power = static_cast<ULONG>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NMASSbins_power);

      // LOGMASSmin
      pname = "min_logmass";
      pname_c = "LOGMASSmin";
      description = "Minimum halo log-mass.";
      options = "real prec";
      this->LOGMASSmin = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->LOGMASSmin);

      // LOGMASSmax
      pname = "max_logmass";
      pname_c = "LOGMASSmax";
      description = "Maximum halo log-mass.";
      options = "real prec";
      this->LOGMASSmax = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->LOGMASSmax);

      // VMAXmin
      pname = "min_vmax";
      pname_c = "VMAXmin";
      description = "Minimum vmax.";
      options = "real prec";
      this->VMAXmin = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->VMAXmin);

      // VMAXmax
      pname = "max_vmax";
      pname_c = "VMAXmax";
      description = "Maximum vmax.";
      options = "real prec";
      this->VMAXmax = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->VMAXmax);

      // RSmin
      pname = "min_rs";
      pname_c = "RSmin";
      description = "Minimum scale radius rs.";
      options = "real prec";
      this->RSmin = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->RSmin);

      // RSmax
      pname = "max_rs";
      pname_c = "RSmax";
      description = "Maximum scale radius rs.";
      options = "real prec";
      this->RSmax = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->RSmax);

      // SPINmin
      pname = "min_spin";
      pname_c = "SPINmin";
      description = "Minimum halo spin.";
      options = "real prec";
      this->SPINmin = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->SPINmin);

      // SPINmax
      pname = "max_spin";
      pname_c = "SPINmax";
      description = "Maximum halo spin.";
      options = "real prec";
      this->SPINmax = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->SPINmax);

      // NPROPbins_bam
      pname = "number_sec_property_bins_assignment_bmt";
      pname_c = "NPROPbins_bam";
      description = "Number of bins for secondary property assignment (BMT).";
      options = "real prec";
      this->NPROPbins_bam = static_cast<real_prec>(HaloAnalysis.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->NPROPbins_bam);


      pname = "measure_power_spectrum";
      pname_c = "Get_power_spectrum";
      description = "Get power spectrum estimates for halo analysis";
      options = "BOOL";
      this->Get_power_spectrum = static_cast<real_prec>(HaloAnalysis.value(pname, false));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Get_power_spectrum);

      // Boolean flags
      #define ADD_BOOL_PARAMa(key, member, desc) \
          pname = key; pname_c = #member; description = desc; options = "BOOL"; \
          this->member = HaloAnalysis.value(key, false); \
          this->collect_params_info(pname, pname_c, this_section, description, options); \
          this->parameter_boolean.emplace_back(pname, this->member);

      ADD_BOOL_PARAMa("set_bins_equal_number_tracers_main_property", set_bins_equal_number_tracers_main_property, "Set bins to equal number of tracers for main property");
      ADD_BOOL_PARAMa("set_bins_equal_number_tracers", set_bins_equal_number_tracers, "Set bins to equal number of tracers");
      ADD_BOOL_PARAMa("use_real_and_redshift_space", use_real_and_redshift_space, "Use both real and redshift space data for analysis of halo clustering.");
      ADD_BOOL_PARAMa("compute_pearson_correlation_coefficient", Get_pearson_coefficient, "Compute Pearson correlation coefficient");
      ADD_BOOL_PARAMa("compute_spearman_correlation_coefficient", Get_spearman_coefficient, "Compute Spearman correlation coefficient in halo analysis");
      ADD_BOOL_PARAMa("perform_principal_component_analysis", Get_PCA, "Perform PCA analysis");
      ADD_BOOL_PARAMa("measure_marked_power_spectrum", Get_marked_power_spectrum, "Measure marked power spectrum within the halo analysis");
      ADD_BOOL_PARAMa("measure_cross_power_spectrum", Get_cross_power_spectrum, "Measure cross power spectrum");
      ADD_BOOL_PARAMa("get_tracer_number_counts", Get_tracer_number_counts, "Get tracer number counts");
      ADD_BOOL_PARAMa("get_tracer_mass_field", Get_tracer_mass_field, "Get tracer mass field");
      ADD_BOOL_PARAMa("get_tracer_vmax_field", Get_tracer_vmax_field, "Get tracer vmax field");
      ADD_BOOL_PARAMa("get_tracer_spin_field", Get_tracer_spin_field, "Get tracer spin field");
      ADD_BOOL_PARAMa("get_tracer_prop_function", Get_prop_function_tracer, "Get tracer property function");
      ADD_BOOL_PARAMa("get_tracer_prop_function_cwt", Get_prop_function_tracer_cwt, "Get tracer property function (CWT)");
      ADD_BOOL_PARAMa("get_cell_local_mach_number", Get_cell_local_mach_number, "Get cell-local Mach number");
      ADD_BOOL_PARAMa("get_tracer_local_mach_number", Get_tracer_local_mach_number, "Get tracer-local Mach number");
      ADD_BOOL_PARAMa("get_tracer_local_dach_number", Get_tracer_local_dach_number, "Get tracer-local DACH number");
      ADD_BOOL_PARAMa("get_tracer_local_dm_density", Get_tracer_local_dm_density, "Get tracer-local DM density");
      ADD_BOOL_PARAMa("get_local_overdensity", Get_local_overdensity, "Get local overdensity");
      ADD_BOOL_PARAMa("get_tidal_anisotropy_at_halo", Get_tidal_anisotropy_at_halo, "Get tidal anisotropy at halo");
      ADD_BOOL_PARAMa("get_peak_height_at_halo", Get_peak_height_at_halo, "Get peak height at halo");
      ADD_BOOL_PARAMa("get_distribution_minimum_separations", get_distribution_min_separations, "Get distribution of minimum separations");

      // Number flags
      #define ADD_NUMBER_PARAM(key, member, desc, typecast) \
          pname = key; pname_c = #member; description = desc; options = "real prec"; \
          this->member = typecast(HaloAnalysis.value(key, 0)); \
          this->collect_params_info(pname, pname_c, this_section, description, options); \
          this->parameter_number.emplace_back(pname, this->member);

      ADD_NUMBER_PARAM("scale_mach_number", Scale_mach_number, "Scaling factor for Mach number", static_cast<real_prec>);

      // Vector parameters
      #define ADD_VECTOR_PARAM(key, member, desc) \
          pname = key; pname_c = #member; description = desc; options = "VECTOR"; \
          this->member = HaloAnalysis[key].get<std::vector<real_prec>>(); \
          this->collect_params_info(pname, pname_c, this_section, description, options); \
          this->parameter_vectors.emplace_back(pname, this->member);

      ADD_VECTOR_PARAM("mass_cuts", MASScuts, "Mass cuts for analysis");
      ADD_VECTOR_PARAM("min_massbins", MASSbins_min, "Minimum mass bin values");
      ADD_VECTOR_PARAM("max_massbins", MASSbins_max, "Maximum mass bin values");
      ADD_VECTOR_PARAM("min_vmaxbins", VMAXbins_min, "Minimum vmax bin values");
      ADD_VECTOR_PARAM("max_vmaxbins", VMAXbins_max, "Maximum vmax bin values");
      ADD_VECTOR_PARAM("min_rsbins", RSbins_min, "Minimum rs bin values");
      ADD_VECTOR_PARAM("max_rsbins", RSbins_max, "Maximum rs bin values");
      ADD_VECTOR_PARAM("min_spinbins", SPINbins_min, "Minimum spin bin values");
      ADD_VECTOR_PARAM("max_spinbins", SPINbins_max, "Maximum spin bin values");
      ADD_VECTOR_PARAM("min_spinbins", VRMSbins_min, "Minimum VRMS bin values");
      ADD_VECTOR_PARAM("max_spinbins", VRMSbins_max, "Maximum VRMS bin values");
      ADD_VECTOR_PARAM("min_virialbins", VIRIALbins_min, "Minimum virial bin values");
      ADD_VECTOR_PARAM("max_virialbins", VIRIALbins_max, "Maximum virial bin values");
      ADD_VECTOR_PARAM("min_btoabins", BTOAbins_min, "Minimum b/a bin values");
      ADD_VECTOR_PARAM("max_btoabins", BTOAbins_max, "Maximum b/a bin values");
      ADD_VECTOR_PARAM("min_ctoabins", CTOAbins_min, "Minimum c/a bin values");
      ADD_VECTOR_PARAM("max_ctoabins", CTOAbins_max, "Maximum c/a bin values");
      ADD_VECTOR_PARAM("min_machbins", MACHbins_min, "Minimum Mach number bin values");
      ADD_VECTOR_PARAM("max_machbins", MACHbins_max, "Maximum Mach number bin values");
      ADD_VECTOR_PARAM("min_dachbins", DACHbins_min, "Minimum DACH number bin values");
      ADD_VECTOR_PARAM("max_dachbins", DACHbins_max, "Maximum DACH number bin values");
      ADD_VECTOR_PARAM("min_biasbins", BIASbins_min, "Minimum bias bin values");
      ADD_VECTOR_PARAM("max_biasbins", BIASbins_max, "Maximum bias bin values");
      ADD_VECTOR_PARAM("min_lcbins", LCbins_min, "Minimum LC bin values");
      ADD_VECTOR_PARAM("max_lcbins", LCbins_max, "Maximum LC bin values");
      ADD_VECTOR_PARAM("min_tabins", TAbins_min, "Minimum TA bin values");
      ADD_VECTOR_PARAM("max_tabins", TAbins_max, "Maximum TA bin values");
      ADD_VECTOR_PARAM("min_phbins", PHbins_min, "Minimum PH bin values");
      ADD_VECTOR_PARAM("max_phbins", PHbins_max, "Maximum PH bin values");

 #undef ADD_BOOL_PARAM
 #undef ADD_REAL_PARAM
 #undef ADD_VECTOR_PARAM

//------------------------------------------     
    const auto & IndividualTracerBias = HaloAnalysis.at("IndividualTracerBias");
    this->forbid_unknown_keys(IndividualTracerBias, 
          {
            "get_individual_tracer_bias",
            "assign_bias_to_full_sample",
            "fraction_to_dilute_sample",
            "get_tracer_relative_bias",
            "get_tracer_quadratic_bias",
            "get_tracer_bias_squared",
            "get_tracer_bias_multipoles",
            "max_kvalue_tracer_bias" ,
            "min_kvalue_tracer_bias" ,
            "min_kvalue_tracer_qbias",
            "max_kvalue_tracer_qbias",
            "max_multipole_tracer_lmbias"
          },
      "IndividualTracerBias");

    
    pname = "get_individual_tracer_bias";
    pname_c = "Get_tracer_bias";
    description = "Compute individual tracer bias.";
    options = "BOOL";
    this->Get_tracer_bias = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Get_tracer_bias);

    if(this->Get_tracer_bias == false)
    {
    #ifdef _USE_BIAS_OBJECT_TO_OBJECT_
        cout << CYAN << endl;
        cout << "Parameter Get_tracer_bias = " << this->Get_tracer_bias << endl;
        cout << "is in conflict with parameter request from preproc directive _USE_BIAS_OBJECT_TO_OBJECT_" << endl;
        string ans;
        cout << "Do you wish to continue? (y/n)" << endl; 
        cin >> ans;
        if(ans=="Y" || ans=="y" || ans=="yes")      
            cout << RESET << endl;
        else
            throw std::invalid_argument("Conflict between parameters not resolved.");
    #endif
    }

    // assign_bias_to_full_sample
    pname = "assign_bias_to_full_sample";
    pname_c = "assign_bias_to_full_sample";
    description = "Assign bias to full tracer sample. If set to true, a plot with the intepolated bias field is shown.";
    options = "BOOL";
    this->assign_bias_to_full_sample = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->assign_bias_to_full_sample);

    // fraction_dilute
    pname = "fraction_to_dilute_sample";
    pname_c = "fraction_dilute";
    description = "Fraction of sample to dilute. ";
    options = "real prec number in the range [0, 1]";
    this->fraction_dilute = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->fraction_dilute);

    // Get_tracer_relative_bias
    pname = "get_tracer_relative_bias";
    pname_c = "Get_tracer_relative_bias";
    description = "Compute relative bias for tracers.";
    options = "BOOL";
    this->Get_tracer_relative_bias = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Get_tracer_relative_bias);

    // Get_tracer_quadratic_bias
    pname = "get_tracer_quadratic_bias";
    pname_c = "Get_tracer_quadratic_bias";
    description = "Compute quadratic bias for tracers.";
    options = "BOOL";
    this->Get_tracer_quadratic_bias = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Get_tracer_quadratic_bias);

    // Get_tracer_bias squared
    pname = "get_tracer_bias_squared";
    pname_c = "Get_tracer_bias_squared";
    description = "Compute bias squared for tracers.";
    options = "BOOL";
    this->Get_tracer_bias_squared = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Get_tracer_bias_squared);


    // Get_tracer_bias multipoles
    pname = "get_tracer_bias_multipoles";
    pname_c = "Get_tracer_bias_multipoles";
    description = "Compute multipola bias for tracers.";
    options = "BOOL";
    this->Get_tracer_bias_multipoles = IndividualTracerBias.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Get_tracer_bias_multipoles);

    // kmax_tracer_bias
    pname = "max_kvalue_tracer_bias";
    pname_c = "kmax_tracer_bias";
    description = "Maximum k-value for tracer bias.";
    options = "real prec";
    this->kmax_tracer_bias = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmax_tracer_bias);

    // kmin_tracer_bias
    pname = "min_kvalue_tracer_bias";
    pname_c = "kmin_tracer_bias";
    description = "Minimum k-value for tracer bias.";
    options = "real prec";
    this->kmin_tracer_bias = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmin_tracer_bias);

    // kmax_tracer_qbias
    pname = "max_kvalue_tracer_qbias";
    pname_c = "kmax_tracer_qbias";
    description = "Maximum k-value for tracer quadratic bias.";
    options = "real prec";
    this->kmax_tracer_qbias = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.1));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmax_tracer_qbias);

    // kmin_tracer_qbias
    pname = "min_kvalue_tracer_qbias";
    pname_c = "kmin_tracer_qbias";
    description = "Minimum k-value for tracer quadratic bias.";
    options = "real prec";
    this->kmin_tracer_qbias = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmin_tracer_qbias);

    // lmax_bias
    pname = "max_multipole_tracer_lmbias";
    pname_c = "lmax_bias";
    description = "Maximum multipole for tracer bias.";
    options = "real prec";
    this->lmax_bias = static_cast<real_prec>(IndividualTracerBias.value(pname, 0.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->lmax_bias);

  }
//------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------
  const auto& CorrelationFunction=cfg.at("CorrelationFunction");
  status=CorrelationFunction.value("status", "off");
  if(ON==status)
  {
    this->enabled_sections.push_back("CorrelationFunction");
    this->forbid_unknown_keys(CorrelationFunction, 
          {
            "status", 
            "output_directory" ,
            "number_bins_separation",
            "min_separation",
            "max_separation",
            "binning_type",
            "property_as_mark"
                    },
      "CorrelationFunction");

            // Output_directory
      pname = "output_directory";
      pname_c = "Output_directory";
      description = "Directory for correlation function outputs.";
      options = "STRING";
      this->Output_directory = CorrelationFunction.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Output_directory);

      // Nbins_cf
      pname = "number_bins_separation";
      pname_c = "Nbins_cf";
      description = "Number of separation bins.";
      options = "ULONG";
      this->Nbins_cf = static_cast<ULONG>(CorrelationFunction.value(pname, 1));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Nbins_cf);

      // rmin_cf
      pname = "min_separation";
      pname_c = "rmin_cf";
      description = "Minimum separation for correlation function.";
      options = "real prec";
      this->rmin_cf = static_cast<real_prec>(CorrelationFunction.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->rmin_cf);

      // rmax_cf
      pname = "max_separation";
      pname_c = "rmax_cf";
      description = "Maximum separation for correlation function.";
      options = "real prec";
      this->rmax_cf = static_cast<real_prec>(CorrelationFunction.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->rmax_cf);

      // rbin_type
      pname = "binning_type";
      pname_c = "rbin_type";
      description = "Binning type for correlation function.";
      options = "STRING";
      this->rbin_type = CorrelationFunction.value(pname, "linear");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->rbin_type);

      // mark
      pname = "property_as_mark";
      pname_c = "mark";
      description = "Property to use as mark in correlation function.";
      options = "STRING";
      this->mark = CorrelationFunction.value(pname, "linear");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->mark);

  }
    //------------------------------------------------------------------------------------------------------

  const auto &AngularPowerSpectrum = cfg.at("AngularPowerSpectrum");
  status =AngularPowerSpectrum.value("status", "off");
  if(ON==status) // missing avoiding unkown keys
  {
      this->enabled_sections.push_back("AngularPowerSpectrum");

      pname = "output_directory";
      pname_c = "Output_directory";
      description = "Directory for angular power spectrum outputs.";
      options = "STRING";
      this->Output_directory = cfg["AngularPowerSpectrum"].value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Output_directory);

      // use_random_catalog_cl
      pname = "use_random_catalog";
      pname_c = "use_random_catalog_cl";
      description = "Use random catalog for angular power spectrum computation.";
      options = "BOOL";
      this->use_random_catalog_cl = AngularPowerSpectrum.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->use_random_catalog_cl);

  }
    //------------------------------------------------------------------------------------------------------
  const auto &FourierBessel = cfg.at("FourierBessel");
  status = FourierBessel.value("status", "off");
  if(ON==status)
  {
      this->enabled_sections.push_back("FourierBessel");
      pname = "output_directory";
      pname_c = "Output_directory";
      description = "Directory for Fourier-Bessel outputs.";
      options = "STRING";
      this->Output_directory = cfg["FourierBessel"].value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->Output_directory);

      // use_random_catalog
      pname = "use_random_catalog";
      pname_c = "use_random_catalog";
      description = "Use random catalog for Fourier-Bessel computation.";
      options = "BOOL";
      this->use_random_catalog = FourierBessel.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->use_random_catalog);
  }


  //------------------------------------------------------------------------------------------------------
  const auto &InitialConditions = cfg.at("InitialConditions");
  status = InitialConditions.value("status", "off");
  this_section="InitialConditions";
  if(ON==status)
  {
    this->enabled_sections.push_back("InitialConditions");
    this->forbid_unknown_keys(InitialConditions, 
          {
             "status", 
              "inputmode",
              "seed",
              "seed_ref",
              "diffcosmorz" ,
              "ic_power_file",
              "ic_white_noise_directory",
              "ic_white_noise_file",
              "ic_file",
              "use_ic_file",
              "ic_input_type",
              "normalize_ic_to_initial_redshift", 
              "initial_redshift_theoretical_power_file",
              "input_directory",
              "transf",
              "read_theoretical_power_spectrum",
              "masskernel" , 
              "masskernel_vel",
              "smoothing_length",
              "velbias",
              "vslength",
              "GaussianRandomField"          },
      "InitialConditions");


      // --------------------- InitialConditions ---------------------

      // input_mode
      pname = "input_mode";
      pname_c = "inputmode";
      description = "Mode of initial conditions input.";
      options = "INT";
      this->inputmode = static_cast<int>(InitialConditions.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->inputmode);

      // seed
      pname = "seed";
      pname_c = "seed";
      description = "Random seed for initial conditions.";
      options = "INT";
      this->seed = static_cast<int>(InitialConditions.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->seed);

      // seed_ref
      pname = "seed_ref";
      pname_c = "seed_ref";
      description = "Reference seed for paired simulations.";
      options = "INT";
      this->seed_ref = static_cast<int>(InitialConditions.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->seed_ref);

      // diffcosmorz
      pname = "diffcosmorz";
      pname_c = "diffcosmorz";
      description = "Use different cosmology in z-direction.";
      options = "BOOL";
      this->diffcosmorz = InitialConditions.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->diffcosmorz);

      // ic_power_file
      pname = "ic_power_file";
      pname_c = "ic_power_file";
      description = "Initial conditions power spectrum file.";
      options = "STRING";
      this->ic_power_file = InitialConditions.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->ic_power_file);

      // ic_white_noise_directory
      pname = "ic_white_noise_directory";
      pname_c = "ic_WN_dir";
      description = "Directory containing IC white noise files.";
      options = "STRING";
      this->ic_WN_dir = InitialConditions.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->ic_WN_dir);

      // ic_white_noise_file
      pname = "ic_white_noise_file";
      pname_c = "ic_WN_file";
      description = "IC white noise file name.";
      options = "STRING";
      this->ic_WN_file = InitialConditions.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->ic_WN_file);

      // ic_file
      pname = "ic_file";
      pname_c = "ic_file";
      description = "Initial conditions file to use directly.";
      options = "STRING";
      this->ic_file = InitialConditions.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->ic_file);

      // use_ic_file
      pname = "use_ic_file";
      pname_c = "use_ic_file";
      description = "Flag to use ic_file instead of generating ICs.";
      options = "BOOL";
      this->use_ic_file = InitialConditions.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->use_ic_file);

      // ic_input_type
      pname = "ic_input_type";
      pname_c = "ic_input_type";
      description = "Type of initial conditions input (delta, potential, etc.).";
      options = "STRING";
      this->ic_input_type = InitialConditions.value(pname, "delta");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->ic_input_type);

      // normalize_ic_to_initial_redshift
      pname = "normalize_ic_to_initial_redshift";
      pname_c = "Normalize_IC_to_initial_redshift";
      description = "Normalize ICs to initial redshift.";
      options = "BOOL";
      this->Normalize_IC_to_initial_redshift = InitialConditions.value(pname, true);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->Normalize_IC_to_initial_redshift);

      // initial_redshift_theoretical_power
      pname = "initial_redshift_thoeretical_power";
      pname_c = "Initial_Redshift_TH_power_file";
      description = "Initial redshift for theoretical power normalization.";
      options = "real prec";
      this->Initial_Redshift_TH_power_file = static_cast<real_prec>(
          InitialConditions.value(pname, 100.));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->Initial_Redshift_TH_power_file);

      // input_directory
      pname = "input_directory";
      pname_c = "dir";
      description = "Directory for initial condition files.";
      options = "STRING";
      this->dir = InitialConditions.value(pname, "null");
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_string.emplace_back(pname, this->dir);

      // read_theoretical_power_spectrum
      pname = "read_theoretical_power_spectrum";
      pname_c = "readPS";
      description = "Read theoretical power spectrum from file.";
      options = "BOOL";
      this->readPS = InitialConditions.value(pname, false);
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_boolean.emplace_back(pname, this->readPS);

      // masskernel
      pname = "masskernel";
      pname_c = "masskernel";
      description = "Mass assignment kernel for ICs.";
      options = "INT";
      this->masskernel = static_cast<int>(InitialConditions.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->masskernel);

      // masskernel_vel
      pname = "masskernel_vel";
      pname_c = "masskernel_vel";
      description = "Velocity assignment kernel for ICs.";
      options = "INT";
      this->masskernel_vel = static_cast<int>(InitialConditions.value(pname, 0));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->masskernel_vel);

      // smoothing_lenght
      pname = "smoothing_lenght";
      pname_c = "slength";
      description = "Smoothing length for ICs.";
      options = "real prec";
      this->slength = static_cast<real_prec>(InitialConditions.value(pname, 1.));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->slength);

      // velbias
      pname = "velbias";
      pname_c = "velbias";
      description = "Velocity bias for ICs.";
      options = "real prec";
      this->velbias = static_cast<real_prec>(InitialConditions.value(pname, 1.));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->velbias);

      // vslength
      pname = "vslength";
      pname_c = "vslength";
      description = "Velocity smoothing length.";
      options = "real prec";
      this->vslength = static_cast<real_prec>(InitialConditions.value(pname, 1.));
      this->collect_params_info(pname, pname_c, this_section, description, options);
      this->parameter_number.emplace_back(pname, this->vslength);

    //------------------------------------------------------------------------------------------------------
    const auto& GaussianRandomField=InitialConditions.at("GaussianRandomField");

    this->forbid_unknown_keys(GaussianRandomField, 
          {
            "max_kvalue_fa",
            "number_of_gaussian_random_field",
            "generate_fixed_amplitude_initial_conditions"
          },
      "GaussianRandomField");

    pname = "max_kvalue_fa";
    pname_c = "Kmax_FA";
    description = "Maximum wavenumber up to which fixed amplitude fluctuations are to be produced. For k>this value, fluctuations are drawn from gaussian distributions.";
    options = "real_prec";
    this->Kmax_FA = static_cast<real_prec>(GaussianRandomField.value(pname, 1.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Kmax_FA);


    pname = "number_of_gaussian_random_field";
    pname_c = "Number_of_GRF";
    description = "Number of GRF to generate";
    options = "int";
    this->Number_of_GRF = static_cast<real_prec>(GaussianRandomField.value(pname, 1.));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->Number_of_GRF);

    pname = "generate_fixed_amplitude_initial_conditions";
    pname_c = "Generate_FA";
    description = "Generate fixed amplitude initial conditions";
    options = "BOOLEAN";
    this->Generate_FA = static_cast<real_prec>(GaussianRandomField.value(pname, true));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->Generate_FA);

  }
    //------------------------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------------------------

  const auto &CosmologicalLibrary = cfg.at("CosmologicalLibrary");
  this_section="CosmologicalLibrary";
  status=CosmologicalLibrary.value("status", "off");
  if(ON==status)
  {
    this->enabled_sections.push_back("CosmologicalLibrary");
    this->forbid_unknown_keys(CosmologicalLibrary, 
       {
          "status", 
          "output_directory",
          "cosmological_redshift",
          "fixed_cosmological_redshift",
          "min_redshift",
          "max_redshift",
          "number_values_redshift",
          "MassFunctionTH",
          "HaloBiasTH",
          "DarkMatterHaloDensityProfile",
          "HOD",
          "PowerSpectrumTH", 
          "CorrelationFunctionTH",
          "AngularPowerSpectrumTH",
          "CosmologicalParameters"
        },
      "CosmologicalLibrary");

// --------------------- CosmologicalLibrary ---------------------

    // output_directory
    pname = "output_directory";
    pname_c = "Output_directory";
    description = "Directory for cosmological library outputs.";
    options = "STRING";
    this->Output_directory = CosmologicalLibrary.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->Output_directory);
    this->set_cosmo_output(this->Output_directory);
    // cosmological_redshift
    pname = "cosmological_redshift";
    pname_c = "redshift";
    description = "Redshift for cosmological calculations.";
    options = "real prec";
    this->redshift = static_cast<real_prec>(CosmologicalLibrary.value(pname, 0.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->redshift);

    // fixed_cosmological_redshift
    pname = "fixed_cosmological_redshift";
    pname_c = "fixed_redshift";
    description = "Flag to use a fixed cosmological redshift.";
    options = "BOOL";
    this->fixed_redshift = CosmologicalLibrary.value(pname, true);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->fixed_redshift);

    // min_redshift
    pname = "min_redshift";
    pname_c = "redshift_min";
    description = "Minimum redshift for calculations.";
    options = "real prec";
    this->redshift_min = static_cast<real_prec>(CosmologicalLibrary.value(pname, 0.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->redshift_min);

    // max_redshift
    pname = "max_redshift";
    pname_c = "redshift_max";
    description = "Maximum redshift for calculations.";
    options = "real prec";
    this->redshift_max = CosmologicalLibrary.value(pname, 1.0);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->redshift_max);

    // number_values_redshift
    pname = "number_values_redshift";
    pname_c = "nbins_redshift";
    description = "Number of redshift bins or values.";
    options = "INT";
    this->nbins_redshift = static_cast<ULONG>(CosmologicalLibrary.value(pname, 10));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->nbins_redshift);

    //-----------
    int index_rel=0;
    const auto &MassFunctionTH = CosmologicalLibrary.at("MassFunctionTH");
    this_section="MassFunctionTH";
    this->forbid_unknown_keys(MassFunctionTH, 
         {
            "mass_function_fit",
            "min_mass_mf"  , 
            "max_mass_mf"  , 
            "min_effective_mass" ,
            "max_effective_mass" ,
            "mass_scaling" ,
            "number_values_mass_function" , 
            "use_external_power_spectrum_file" ,
            "external_linear_power_spectrum_file",
            "show_plot"
        },
      "MassFunctionTH");

    // --------------------- MassFunctionTH ---------------------

    // mass_function_fit
    pname = "mass_function_fit";
    pname_c = "mass_function_fit";
    description = "Choice of theoretical mass function fit (e.g., Tinker).";
    options = "STRING";
    this->mass_function_fit = MassFunctionTH.value(pname, "Tinker");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->mass_function_fit);

    // min_mass_mf
    pname = "min_mass_mf";
    pname_c = "M_min_mf";
    description = "Minimum halo mass for mass function.";
    options = "real prec";
    this->M_min_mf = MassFunctionTH.value(pname, 1e7);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->M_min_mf);

    // max_mass_mf
    pname = "max_mass_mf";
    pname_c = "M_max_mf";
    description = "Maximum halo mass for mass function.";
    options = "real prec";
    this->M_max_mf = MassFunctionTH.value(pname, 1e16);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->M_max_mf);

    // min_effective_mass
    pname = "min_effective_mass";
    pname_c = "M_min_effective";
    description = "Minimum effective mass for mass function calculations.";
    options = "real prec";
    this->M_min_effective = MassFunctionTH.value(pname, 1e16);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->M_min_effective);

    // max_effective_mass
    pname = "max_effective_mass";
    pname_c = "M_max_effective";
    description = "Maximum effective mass for mass function calculations.";
    options = "real prec";
    this->M_max_effective = MassFunctionTH.value(pname, 1e16);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->M_max_effective);

    // mass_scaling
    pname = "mass_scaling";
    pname_c = "scale_mf";
    description = "Scaling type for mass function (e.g., linear, log).";
    options = "STRING";
    this->scale_mf = MassFunctionTH.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->scale_mf);

    // number_values_mass_function
    pname = "number_values_mass_function";
    pname_c = "npoints_mf";
    description = "Number of mass bins or values for the mass function.";
    options = "INT";
    this->npoints_mf = static_cast<ULONG>(MassFunctionTH.value(pname, 10));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->npoints_mf);

    // output_file_mass_function
    pname = "output_file_mass_function";
    pname_c = "mass_function_output_file";
    description = "Filename for output mass function.";
    options = "STRING";
    this->mass_function_output_file = MassFunctionTH.value(pname, "mass_function");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->mass_function_output_file);

    // use_external_power_spectrum_file
    pname = "use_external_power_spectrum_file";
    pname_c = "use_file_power";
    description = "Flag to use an external linear power spectrum file.";
    options = "BOOL";
    this->use_file_power = MassFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->use_file_power);

    // external_linear_power_spectrum_file
    pname = "external_linear_power_spectrum_file";
    pname_c = "file_power";
    description = "Path to external linear power spectrum file.";
    options = "STRING";
    this->file_power = MassFunctionTH.value(pname, "null");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->file_power);

    // external_linear_power_spectrum_file
    pname = "show_plot";
    pname_c = "show_mass_function_plot";
    description = "Show the mass function plot with python";
    options = "BOOL";
    this->show_mass_function_plot = MassFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_mass_function_plot);




    //-----------
    const auto &HaloBiasTH = CosmologicalLibrary.at("HaloBiasTH");
    this_section="HaloBiasTH";
    this->forbid_unknown_keys(HaloBiasTH, 
         {"halo_bias_fit",
          "show_plot"
         },
      "HaloBiasTH");

// external_linear_power_spectrum_file
    pname = "show_plot";
    pname_c = "show_halo_mass_bias_plot";
    description = "Show the mass function plot with python";
    options = "BOOL";
    this->show_halo_mass_bias_plot = HaloBiasTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_halo_mass_bias_plot);


    this->halo_mass_bias_fit = HaloBiasTH.value("halo_bias_fit", "Tinker");
    this->parameter_string.push_back(make_pair("halo_bias_fit", this->halo_mass_bias_fit));
    this->halo_mass_bias_output_file = HaloBiasTH.value("output_file_halo_bias", "halo_bias");//leave it with standard file name (but giving a change to set value in parf)
    this->effective_halo_mass_bias_output_file = HaloBiasTH.value("output_file_effective_bias", "effective_halo_bias");
    this->effective_halo_mean_number_density_output_file = HaloBiasTH.value("output_file_effective_halo_mean_number_density", "effective_halo_mean_number_density");
    //-----------

    const auto &DarkMatterHaloDensityProfile = CosmologicalLibrary.at("DarkMatterHaloDensityProfile");
    this_section="DarkMatterHaloDensityProfile";
    this->forbid_unknown_keys(DarkMatterHaloDensityProfile, 
         {
            "compute_density_profile" ,
            "density_profile_model",
            "min_radius_dp" ,
            "max_radius_dp",
            "radius_dp_scaling",
            "number_rvalues_density_profile", 
            "coef_concentration_amp",
            "coef_concentration" ,
            "show_plot_r",
            "min_kvalue_dp" ,
            "max_kvalue_dp"  ,
            "kval_dp_scaling" ,
            "number_kvalues_density_profile",
            "show_plot_k"
         },
      "DarkMatterHaloDensityProfile");

// --------------------- DarkMatterHaloDensityProfile ---------------------

    // compute_density_profile
    pname = "compute_density_profile";
    pname_c = "compute_density_profile";
    description = "Flag to compute the dark matter halo density profile.";
    options = "BOOL";
    this->compute_density_profile = DarkMatterHaloDensityProfile.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_density_profile);

    // density_profile_model
    pname = "density_profile_model";
    pname_c = "density_profile";
    description = "Choice of dark matter halo density profile model (e.g., nfw).";
    options = "STRING";
    this->density_profile = DarkMatterHaloDensityProfile.value(pname, "nfw");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->density_profile);

    // min_radius_dp
    pname = "min_radius_dp";
    pname_c = "rmin_dp";
    description = "Minimum radius for density profile computation.";
    options = "real prec";
    this->rmin_dp = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 0.01));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->rmin_dp);

    // max_radius_dp
    pname = "max_radius_dp";
    pname_c = "rmax_dp";
    description = "Maximum radius for density profile computation.";
    options = "real prec";
    this->rmax_dp = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 10.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->rmax_dp);

    // radius_dp_scaling
    pname = "radius_dp_scaling";
    pname_c = "scale_dp_r";
    description = "Scaling type for radius in density profile (linear or log).";
    options = "STRING";
    this->scale_dp_r = DarkMatterHaloDensityProfile.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->scale_dp_r);

    // number_rvalues_density_profile
    pname = "number_rvalues_density_profile";
    pname_c = "npoints_dp_r";
    description = "Number of radius values for the density profile.";
    options = "INT";
    this->npoints_dp_r = static_cast<int>(DarkMatterHaloDensityProfile.value(pname, 100));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->npoints_dp_r);

    // density_profile_r_output_file
    pname = "density_profile_r_output_file";
    pname_c = "density_profile_r_output_file";
    description = "Output file name for density profile as a function of radius.";
    options = "STRING";
    this->density_profile_r_output_file = DarkMatterHaloDensityProfile.value(pname, "density_profile_r");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->density_profile_r_output_file);

    // coef_concentration_amp
    pname = "coef_concentration_amp";
    pname_c = "coef_concentration_amp";
    description = "Amplitude coefficient for halo concentration.";
    options = "real prec";
    this->coef_concentration_amp = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->coef_concentration_amp);

    // coef_concentration
    pname = "coef_concentration";
    pname_c = "coef_concentration";
    description = "Halo concentration coefficient.";
    options = "real prec";
    this->coef_concentration = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->coef_concentration);

    // min_kvalue_dp
    pname = "min_kvalue_dp";
    pname_c = "kmin_dp";
    description = "Minimum k-value for density profile in Fourier space.";
    options = "real prec";
    this->kmin_dp = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 0.01));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmin_dp);

    // max_kvalue_dp
    pname = "max_kvalue_dp";
    pname_c = "kmax_dp";
    description = "Maximum k-value for density profile in Fourier space.";
    options = "real prec";
    this->kmax_dp = static_cast<real_prec>(DarkMatterHaloDensityProfile.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmax_dp);

    // kval_dp_scaling
    pname = "kval_dp_scaling";
    pname_c = "scale_dp_k";
    description = "Scaling type for k-values in density profile (linear or log).";
    options = "STRING";
    this->scale_dp_k = DarkMatterHaloDensityProfile.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->scale_dp_k);

    // number_kvalues_density_profile
    pname = "number_kvalues_density_profile";
    pname_c = "npoints_dp_k";
    description = "Number of k-values for density profile.";
    options = "INT";
    this->npoints_dp_k = static_cast<int>(DarkMatterHaloDensityProfile.value(pname, 100));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->npoints_dp_k);

    this->density_profile_k_output_file = DarkMatterHaloDensityProfile.value("density_profile_k_output_file", "density_profile_k");
    

    pname = "show_plot_r";
    pname_c = "show_density_profile_r";
    description = "Show the halo density profile with python";
    options = "BOOL";
    this->show_density_profile_r_plot= DarkMatterHaloDensityProfile.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_density_profile_r_plot);

    pname = "show_plot_k";
    pname_c = "show_density_profile_k";
    description = "Show the halo density in Fourier profile with python";
    options = "BOOL";
    this->show_density_profile_k_plot = DarkMatterHaloDensityProfile.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_density_profile_k_plot);

    //------------
    index_rel=0;
    const auto &HOD = CosmologicalLibrary.at("HOD");
    this_section="HOD";
    this->forbid_unknown_keys(HOD, 
         {
            "hod_model",
            "muno_hod",
            "alpha_hod",
            "mmin_hod",
            "scatter_hod"
         },
      "HOD");


    // --------------------- HOD Parameters ---------------------

    // hod_model
    pname = "hod_model";
    pname_c = "hod_model";
    description = "Halo Occupation Distribution (HOD) model. Choices are 0:, 1, 2, 3";
    options = "INT";
    this->hod_model = static_cast<int>(HOD.value(pname, 0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->hod_model);

    // mass_one_hod
    pname = "mass_one_hod";
    pname_c = "muno_hod";
    description = "Characteristic mass M1 for HOD.";
    options = "real prec";
    this->muno_hod = static_cast<real_prec>(HOD.value(pname, 1e12));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->muno_hod);

    // alpha_hod
    pname = "alpha_hod";
    pname_c = "alpha_hod";
    description = "Power-law slope alpha in the HOD model.";
    options = "real prec";
    this->alpha_hod = static_cast<real_prec>(HOD.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->alpha_hod);

    // mmin_hod
    pname = "mmin_hod";
    pname_c = "mmin_hod";
    description = "Minimum halo mass for hosting galaxies in HOD.";
    options = "real prec";
    this->mmin_hod = static_cast<real_prec>(HOD.value(pname, 1e12));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->mmin_hod);

    // scatter_hod
    pname = "scatter_hod";
    pname_c = "scatter_hod";
    description = "Log-normal scatter in the HOD relation.";
    options = "real prec";
    this->scatter_hod = static_cast<real_prec>(HOD.value(pname, 0.4));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->scatter_hod);


    //-----------
    const auto &PowerSpectrumTH = CosmologicalLibrary.at("PowerSpectrumTH");
    this_section="PowerSpectrumTH";

    this->forbid_unknown_keys(PowerSpectrumTH, 
         {
            "compute_linear_power_spectrum",
            "compute_non_linear_power_spectrum_halofit",
            "compute_non_linear_power_spectrum_pt",
            "compute_halo_model_galaxy_power_spectrum",
            "min_kvalue_power_spectrum",
            "max_kvalue_power_spectrum",
            "kval_scaling",
            "min_kvalue_for_integration",
            "max_kvalue_for_integration",
            "number_of_values_power_spectrum",
            "show_plot"        
         },
      "PowerSpectrumTH");

// --------------------- Power Spectrum Parameters ---------------------

    // compute_linear_power_spectrum
    pname = "compute_linear_power_spectrum";
    pname_c = "compute_linear_power_spectrum";
    description = "Compute the linear matter power spectrum.";
    options = "BOOL";
    this->compute_linear_power_spectrum =
        PowerSpectrumTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_linear_power_spectrum);

    // compute_nonlinear_power_spectrum
    pname = "compute_non_linear_power_spectrum_halofit";
    pname_c = "compute_non_linear_power_spectrum_halofit";
    description = "Compute the non-linear matter power spectrum using halo fit algorithm";
    options = "BOOL";
    this->compute_non_linear_power_spectrum_halofit =
        PowerSpectrumTH.value(pname, false); // note key typo preserved
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_non_linear_power_spectrum_halofit);

    // compute_nonlinear_power_spectrum using pt
    pname = "compute_non_linear_power_spectrum_pt";
    pname_c = "compute_non_linear_power_spectrum_pt";
    description = "Compute the non-linear matter power spectrum using eulerian perturbation theory. This is available when the halo fit power is computed.";
    options = "BOOL";
    this->compute_non_linear_power_spectrum_pt = PowerSpectrumTH.value(pname, false); // note key typo preserved
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_non_linear_power_spectrum_pt);

    // compute_nonlinear_power_spectrum using pt
    pname = "compute_halo_model_galaxy_power_spectrum";
    pname_c = "compute_halo_model_galaxy_power_spectrum";
    description = "Compute the non-linear galaxy power spectrum using Halo Model.";
    options = "BOOL";
    this->compute_halo_model_galaxy_power_spectrum = PowerSpectrumTH.value(pname, false); // note key typo preserved
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_halo_model_galaxy_power_spectrum);



    // min_kvalue_power_spectrum
    pname = "min_kvalue_power_spectrum";
    pname_c = "kmin_ps";
    description = "Minimum k-value for power spectrum calculation.";
    options = "real prec";
    this->kmin_ps =
        static_cast<real_prec>(PowerSpectrumTH.value(pname, 0.001));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmin_ps);

    // max_kvalue_power_spectrum
    pname = "max_kvalue_power_spectrum";
    pname_c = "kmax_ps";
    description = "Maximum k-value for power spectrum calculation.";
    options = "real prec";
    this->kmax_ps =
        static_cast<real_prec>(PowerSpectrumTH.value(pname, 1.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->kmax_ps);

    // kval_scaling
    pname = "kval_scaling";
    pname_c = "scale_ps";
    description = "Scaling of k-values (linear or log).";
    options = "STRING";
    this->scale_ps =
        PowerSpectrumTH.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->scale_ps);

    // number_of_values_power_spectrum
    pname = "number_of_values_power_spectrum";
    pname_c = "npoints_ps";
    description = "Number of k-values for the power spectrum.";
    options = "ULONG";
    this->npoints_ps =
        static_cast<ULONG>(PowerSpectrumTH.value(pname, 10));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->npoints_ps);

    // output_file_linear_power_spectrum
    pname = "output_file_linear_power_spectrum";
    pname_c = "linear_matter_ps_output_file";
    description = "Output file for linear matter power spectrum.";
    options = "STRING";
    this->linear_matter_ps_output_file =PowerSpectrumTH.value(pname, "linear_matter_power_spectrum_EH");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->linear_matter_ps_output_file);

    // output_file_non_linear_hm_power_spectrum
    pname = "output_file_non_linear_halo_fit_power_spectrum";
    pname_c = "non_linear_matter_ps_halo_fit_output_file";
    description = "Output file for nonlinear matter power spectrum from halo model.";
    options = "STRING";
    this->non_linear_matter_ps_halo_fit_output_file =
        PowerSpectrumTH.value(pname, "non_linear_matter_power_spectrum_hf");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->non_linear_matter_ps_halo_fit_output_file);

    // output_file_non_linear_pt_power_spectrum
    pname = "output_file_non_linear_pt_power_spectrum";
    pname_c = "non_linear_matter_ps_pt_output_file";
    description = "Output file for nonlinear matter power spectrum from perturbation theory.";
    options = "STRING";
    this->non_linear_matter_ps_pt_output_file =
        PowerSpectrumTH.value(pname, "non_linear_matter_power_spectrum_pt");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->non_linear_matter_ps_pt_output_file);

    // galaxy_power_spectrum_halo_model_output_file
    pname = "galaxy_ps_halo_model_output_file";
    pname_c = "galaxy_ps_halo_model_output_file";
    description = "Output file for galaxy power spectrum from halo model.";
    options = "STRING";
    this->galaxy_ps_halo_model_output_file =
        PowerSpectrumTH.value(pname, "galaxy_ps_halo_model");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->galaxy_ps_halo_model_output_file);

 // external_linear_power_spectrum_file
    pname = "show_plot";
    pname_c = "show_power_spectrum_plot";
    description = "Show the halo model galaxy power spectrum";
    options = "BOOL";
    this->show_power_spectrum_plot = PowerSpectrumTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->show_power_spectrum_plot);

//-----------
    const auto &CorrelationFunctionTH = CosmologicalLibrary.at("CorrelationFunctionTH");
    this_section="CorrelationFunctionTH";

    this->forbid_unknown_keys(CorrelationFunctionTH, 
         {
            "compute_linear_correlation_function",
            "compute_non_linear_correlation_function_halofit",
            "compute_non_linear_correlation_function_pt",
            "compute_halo_model_galaxy_correlation_function",
            "min_r_value_correlation_function",
            "max_r_value_correlation_function",
            "separation_scaling",
            "number_of_values_correlation_function"
         },
      "CorrelationFunctionTH");

   // --------------------- CorrelationFunctionTH ---------------------

    // compute_linear_correlation_function
    pname = "compute_linear_correlation_function";
    pname_c = "compute_linear_correlation_function";
    description = "Compute linear correlation function.";
    options = "BOOL";
    this->compute_linear_correlation_function = CorrelationFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_linear_correlation_function);

    // compute_non_linear_correlation_function
    pname = "compute_non_linear_correlation_function_halofit";
    pname_c = "compute_non_linear_correlation_function_halofit";
    description = "Compute non-linear correlation function using the halo fit power spectrum";
    options = "BOOL";
    this->compute_non_linear_correlation_function_halofit = CorrelationFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_non_linear_correlation_function_halofit);

    // compute_non_linear_correlation_function
    pname = "compute_non_linear_correlation_function_pt";
    pname_c = "compute_non_linear_correlation_function_pt";
    description = "Compute non-linear correlation function using the halo fit power spectrum";
    options = "BOOL";
    this->compute_non_linear_correlation_function_pt = CorrelationFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_non_linear_correlation_function_pt);

    // compute_non_linear_correlation_function
    pname = "compute_halo_model_galaxy_correlation_function";
    pname_c = "compute_halo_model_galaxy_correlation_function";
    description = "Compute galaxy correlation function using the halo model";
    options = "BOOL";
    this->compute_halo_model_galaxy_correlation_function = CorrelationFunctionTH.value(pname, false);
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_boolean.emplace_back(pname, this->compute_halo_model_galaxy_correlation_function);


    // min_r_value_correlation_function
    pname = "min_r_value_correlation_function";
    pname_c = "rmin_cf";
    description = "Minimum separation r for correlation function.";
    options = "real prec";
    this->rmin_cf = static_cast<real_prec>(CorrelationFunctionTH.value(pname, 0.01));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->rmin_cf);

    // max_r_value_correlation_function
    pname = "max_r_value_correlation_function";
    pname_c = "rmax_cf";
    description = "Maximum separation r for correlation function.";
    options = "real prec";
    this->rmax_cf = static_cast<real_prec>(CorrelationFunctionTH.value(pname, 100.0));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->rmax_cf);

    // separation_scaling
    pname = "separation_scaling";
    pname_c = "scale_cf";
    description = "Scaling type for separation bins.";
    options = "STRING";
    this->scale_cf = CorrelationFunctionTH.value(pname, "log");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->scale_cf);

    // number_of_values_correlation_function
    pname = "number_of_values_correlation_function";
    pname_c = "npoints_cf";
    description = "Number of points in correlation function.";
    options = "real prec";
    this->npoints_cf = static_cast<real_prec>(CorrelationFunctionTH.value(pname, 10));
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_number.emplace_back(pname, this->npoints_cf);



// output_file_linear_power_spectrum
    pname = "output_file_linear_correlation_function";
    pname_c = "linear_matter_cf_output_file";
    description = "Output file for linear matter correlation function";
    options = "STRING";
    this->linear_matter_cf_output_file =CorrelationFunctionTH.value(pname, "linear_matter_correlation_function_EH");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->linear_matter_ps_output_file);

    // output_file_non_linear_hm_power_spectrum
    pname = "output_file_non_linear_halo_fit_correlation function";
    pname_c = "non_linear_matter_cf_halo_fit_output_file";
    description = "Output file for nonlinear matter power spectrum from halo model.";
    options = "STRING";
    this->non_linear_matter_cf_halo_fit_output_file =
      CorrelationFunctionTH.value(pname, "non_linear_matter_correlation_function_hf");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->non_linear_matter_ps_halo_fit_output_file);


  // galaxy_power_spectrum_halo_model_output_file
    pname = "galaxy_cf_halo_model_output_file";
    pname_c = "galaxy_cf_halo_model_output_file";
    description = "Output file for galaxy correlation function from halo model.";
    options = "STRING";
    this->galaxy_cf_halo_model_output_file =
        CorrelationFunctionTH.value(pname, "galaxy_cf_halo_model");
    this->collect_params_info(pname, pname_c, this_section, description, options);
    this->parameter_string.emplace_back(pname, this->galaxy_cf_halo_model_output_file);

    //-----------
    const auto &CosmologicalParameters = CosmologicalLibrary.at("CosmologicalParameters");
    this_section="CosmologicalParameters";
    this->forbid_unknown_keys(CosmologicalParameters, 
         {
          "Hubble",
          "om_matter",
          "om_cdm", 
          "om_radiation",
          "om_baryons",
          "om_vac",
          "om_k",
          "spectral_index",
          "wde_eos",
          "hubble" ,
          "N_eff" ,
          "sigma8",
          "A_s",
          "alpha_s",
          "Tcmb",
          "RR",
          "use_wiggles",
          "kstar",
          "GAL_BIAS", 
          "Amc",
          "Get_SO_from_BN",
          "Delta_SO"
        },
      this_section);

// --------------------- CosmologicalLibrary ---------------------

    #define ADD_REAL_PARAM(lib, key, var, def, desc) \
    pname = key; \
    pname_c = var; \
    description = desc; \
    options = "real prec"; \
    this->var = static_cast<real_prec>(lib.value(key, def)); \
    this->collect_params_info(pname, pname_c, this_section, description, options); \
    this->parameter_number.emplace_back(pname, this->var);

    #define ADD_BOOL_PARAM(lib, key, var, def, desc) \
    pname = key; \
    pname_c = var; \
    description = desc; \
    options = "BOOL"; \
    this->var = lib.value(key, def); \
    this->collect_params_info(pname, pname_c, this_section, description, options); \
    this->parameter_boolean.emplace_back(pname, this->var);

    // Real parameters
    ADD_REAL_PARAM(CosmologicalParameters, "Hubble", Hubble, 0, "Hubble parameter.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_matter", om_matter, COSMOPARS::Om_matter, "Matter density parameter Omega_m.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_cdm", om_cdm, COSMOPARS::Om_cdm, "Cold dark matter density Omega_cdm.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_radiation", om_radiation, COSMOPARS::Om_radiation, "Radiation density Omega_r.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_baryons", om_baryons, COSMOPARS::Om_baryons, "Baryon density Omega_b.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_vac", om_vac, COSMOPARS::Om_vac, "Vacuum energy density Omega_Lambda.");
    ADD_REAL_PARAM(CosmologicalParameters, "om_k", om_k, COSMOPARS::Om_k, "Curvature density Omega_k.");
    ADD_REAL_PARAM(CosmologicalParameters, "spectral_index", spectral_index, COSMOPARS::spectral_index, "Primordial spectral index n_s.");
    ADD_REAL_PARAM(CosmologicalParameters, "wde_eos", wde_eos, COSMOPARS::wde_eos, "Dark energy equation of state w.");
    ADD_REAL_PARAM(CosmologicalParameters, "hubble", hubble, COSMOPARS::hubble, "Dimensionless Hubble parameter h.");
    ADD_REAL_PARAM(CosmologicalParameters, "N_eff", N_eff, COSMOPARS::N_eff, "Effective number of neutrino species.");
    ADD_REAL_PARAM(CosmologicalParameters, "sigma8", sigma8, COSMOPARS::sigma8, "Amplitude of matter fluctuations sigma_8.");
    ADD_REAL_PARAM(CosmologicalParameters, "A_s", A_s, 0, "Primordial scalar amplitude.");
    ADD_REAL_PARAM(CosmologicalParameters, "alpha_s", alpha_s, COSMOPARS::alpha_s, "Running of scalar spectral index.");
    ADD_REAL_PARAM(CosmologicalParameters, "Tcmb", Tcmb, COSMOPARS::Tcmb, "CMB temperature.");
    ADD_REAL_PARAM(CosmologicalParameters, "RR", RR, COSMOPARS::RR, "Reference scale RR.");
    ADD_REAL_PARAM(CosmologicalParameters, "kstar", kstar, 0, "Characteristic k value.");
    ADD_REAL_PARAM(CosmologicalParameters, "GAL_BIAS", GAL_BIAS, 1, "Galaxy bias factor.");
    ADD_REAL_PARAM(CosmologicalParameters, "Amc", Amc, 0, "Amplitude of power for mode coupling (perturbation theory)");
    ADD_REAL_PARAM(CosmologicalParameters, "Delta_SO", Delta_SO, COSMOPARS::Delta_SO, "Spherical overdensity Delta_SO.");

    // Boolean parameters
    ADD_BOOL_PARAM(CosmologicalParameters, "use_wiggles", use_wiggles, COSMOPARS::use_wiggles, "Include BAO wiggles.");
    ADD_BOOL_PARAM(CosmologicalParameters, "Get_SO_from_BN", Get_SO_from_BN, false, "Compute SO from BN.");
  #undef ADD_BOOL_PARAM
  #undef ADD_REAL_PARAM
 
}


}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define _USE_JSON_
void Params::read_pars(string file)
{

#ifndef _USE_JSON_

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
        else if (par_name == "lmax_bias")
          {
            this->lmax_bias  = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->lmax_bias)); 
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
          else if (par_name == "assign_bias_to_full_sample")
          {
            if (par_value=="true")this->assign_bias_to_full_sample = true;
            else this->assign_bias_to_full_sample = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->assign_bias_to_full_sample)); 
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
        else if (par_name == "number_mass_bins")
          {
            this->NMASSbins = atoi(par_value.c_str());
            this->parameter_number.push_back(make_pair(par_name, this->NMASSbins)); 
          }
        else if (par_name == "number_mass_bins_mass_function")
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
          else if (par_name == "Get_tracer_bias_field")
          {
            if (par_value=="true")this->Get_tracer_bias_field= true;
            else this->Get_tracer_bias_field = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_bias_field)); 
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
        else if (par_name == "Get_tracer_spin_bullock_field")
          {
            if (par_value=="true")this->Get_tracer_spin_bullock_field= true;
            else this->Get_tracer_spin_bullock_field = false;
            this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_spin_bullock_field)); 
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
           else if (par_name == "Get_tracer_mean_number_density")
            {
              if (par_value=="true")this->Get_tracer_mean_number_density= true;
              else this->Get_tracer_mean_number_density = false;
              this->parameter_boolean.push_back(make_pair(par_name, this->Get_tracer_mean_number_density));
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
      else if (par_name == "Delta_SO")
      {
        this->Delta_SO = static_cast<real_prec>(atof(par_value.c_str()));
	      this->parameter_number.push_back(make_pair(par_name, this->Delta_SO));
      }
      else if (par_name == "Get_SO_from_BN")
      {
        if(par_value== "true")this->Get_SO_from_BN=true;
        else if(par_value== "false")this->Get_SO_from_BN=false;
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
      else if (par_name == "compute_linear_correlation_function")
	     {
	       if(par_value=="true")compute_linear_correlation_function=true;
	       else if(par_value=="false")compute_linear_correlation_function=false;
	    }
      else if (par_name == "compute_non_linear_correlation_function"){
        if(par_value=="true")compute_non_linear_correlation_function=true;
        else if(par_value=="false")compute_non_linear_correlation_function=false;
      }
      else if (par_name == "compute_linear_power_spectrum"){
        if(par_value=="true")compute_linear_power_spectrum=true;
        else if(par_value=="false")compute_linear_power_spectrum=false;
      }
      else if (par_name == "use_file_power"){
        if(par_value=="true")use_file_power=true;
        else if(par_value=="false")use_file_power=false;
      }
      else if (par_name == "compute_non_linear_power_spectrum"){
        if(par_value=="true")compute_non_linear_power_spectrum=true;
        else if(par_value=="false")compute_non_linear_power_spectrum=false;
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
#endif

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
  this->f_baryon=COSMOPARS::f_baryon;
  this->kmin_integration=0.00001;
  this->kmax_integration=500.000;
  if(true ==this->Get_SO_from_BN)
    this->Delta_SO=COSMOPARS::Delta_SO;
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
  this->s_cosmo_pars.f_baryon = this->f_baryon;// this->om_baryons/this->om_matter;
  this->s_cosmo_pars.use_wiggles =this->use_wiggles;
  this->s_cosmo_pars.RR =this->RR;
  this->s_cosmo_pars.Tcmb =this->Tcmb;
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
  this->set_mass_function_output_file(this->cosmo_output+mass_function_output_file+ll+mass_function_fit+zr+dat);
  this->set_halo_mass_bias_output_file(this->cosmo_output+halo_mass_bias_output_file+ll+halo_mass_bias_fit+zr+dat);
  this->set_effective_halo_mass_bias_output_file(this->cosmo_output+effective_halo_mass_bias_output_file+ll+halo_mass_bias_fit+ll+mass_function_fit+zr+dat);
  this->set_effective_halo_mean_number_density_output_file(this->cosmo_output+effective_halo_mean_number_density_output_file+ll+mass_function_fit+zr+dat);
  this->set_linear_matter_cf_output_file(this->cosmo_output+linear_matter_cf_output_file+zr+dat);
  this->set_non_linear_matter_cf_halo_fit_output_file(this->cosmo_output+non_linear_matter_cf_halo_fit_output_file+zr+dat);
  this->set_linear_matter_ps_output_file(this->cosmo_output+linear_matter_ps_output_file+zr+dat);
  this->set_non_linear_matter_ps_halo_fit_output_file(this->cosmo_output+non_linear_matter_ps_halo_fit_output_file+zr+dat);
  this->set_non_linear_matter_ps_pt_output_file(this->cosmo_output+non_linear_matter_ps_pt_output_file+zr+dat);
  this->set_density_profile_r_output_file(this->cosmo_output+density_profile_r_output_file+ll+density_profile+zr+dat);
  this->set_galaxy_ps_halo_model_output_file(this->cosmo_output+galaxy_ps_halo_model_output_file+zr+dat);
  this->set_density_profile_k_output_file(this->cosmo_output+density_profile_k_output_file+ll+density_profile+zr+dat);
  this->set_galaxy_cf_halo_model_output_file(this->cosmo_output+galaxy_cf_halo_model_output_file+zr+dat);
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
  std::cout<<CYAN;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"CosmiCodes                                                                 *"<<std::endl;
  std::cout<<__PRETTY_FUNCTION__<<"                                                      *"<<std::endl;

 if(this->enabled_sections.size()==0){
    std::cout<<"No active sections in json file "<<this->_par_file()<<endl;
    std:: cout<<"Enable at least one section to run the code, genio."<<endl;
    std::cout<<"***************************************************************************"<<std::endl;
    std::cout<<"***************************************************************************"<<std::endl;
 }
 else{
  std::cout<<"Enabled Sections in JSON parameter file "<< this->_par_file()<<":"<<endl;
  for(int i=0;i<this->enabled_sections.size();++i)
    cout<<GREEN<<enabled_sections[i]<<RESET<<endl;
 
  std::cout<<"***************************************************************************"<<std::endl;
  std::cout<<"Input values of parameters in parameter file                              *"<<std::endl;
   std::cout<<CYAN<<"NUMERICAL PARAMTERS"<<RESET<<std::endl;
  for(int i=0;i<parameter_number.size();++i)
    std::cout<<GREEN<<parameter_number[i].first<<" = "<<BLUE<<parameter_number[i].second<<RESET<<endl;
  std::cout<<CYAN<<"STRING PARAMETERS"<<RESET<<std::endl;
  for(int i=0;i<parameter_string.size();++i)
    std::cout<<GREEN<<parameter_string[i].first<<" = "<<BLUE<<parameter_string[i].second<<RESET<<endl;
  std::cout<<CYAN<<"BOOLEAN PARAMETERS"<<RESET<<std::endl;
  for(int i=0;i<parameter_boolean.size();++i)
    std::cout<<GREEN<<parameter_boolean[i].first<<" = "<<BLUE<<parameter_boolean[i].second<<RESET<<endl;
  std::cout<<CYAN<<"VECTOR PARAMETERS"<<RESET<<std::endl;
  for(int i=0;i<parameter_vectors.size();++i)
  {
   std::cout<<GREEN<<parameter_vectors[i].first<<" = ";
   for(int j=0;j<parameter_vectors[i].second.size();++j)cout<<BLUE<<parameter_vectors[i].second[j]<<" "<<RESET;
   std::cout<<endl;
  }
  std::cout<<CYAN<<"VECTOR-STRING PARAMETERS"<<RESET<<std::endl;
  for(int i=0;i<parameter_vector_string.size();++i)
  {
   std::cout<<GREEN<<parameter_vector_string[i].first<<" = ";
   for(int j=0;j<parameter_vector_string[i].second.size();++j)cout<<BLUE<<parameter_vector_string[i].second[j]<<" "<<RESET;
   std::cout<<endl;
  }
  std::cout<<GREEN<<"******************SUMMARY PARAMETERS***************************************"<<RESET<<std::endl;
  std::cout<<GREEN<<"Number of numerical parameters = "<<parameter_number.size()<<endl;
  std::cout<<GREEN<<"Number of string parameters = "<<parameter_string.size()<<endl;
  std::cout<<GREEN<<"Number of boolean parameters = "<<parameter_boolean.size()<<endl;
  std::cout<<GREEN<<"Number of vector parameters = "<<parameter_vectors.size()<<endl;
  std::cout<<GREEN<<"Total = "<<parameter_boolean.size()+parameter_string.size()+parameter_number.size()+parameter_vectors.size()<<endl;
  std::cout<<GREEN<<"***************************************************************************"<<RESET<<std::endl;
 }
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


}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
