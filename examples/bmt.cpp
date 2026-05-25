////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for BMT.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
////////////////////////////////////////////////////////////////////////////
# include "../include/CosmiCalcLIB.hpp"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void ind_bias(Params &params)
{
    ScreenOutput So;
    FileOutput File;

    fs::path filetr = fs::path(params._Input_dir_cat()) / params._file_catalogue(); 

    Catalogue cat(params, "TRACER"); // Feed params into the Catlaog class and define an object of that type
    cat.read_catalog_new(filetr.string()); // Read the asciii file

    HaloTools htools(params, cat); 

    vector<real_prec> dm_field(params._NGRID(), 0); //COntainer for dm
    fs::path filedm = fs::path(params._Input_Directory_X()) / params._Name_Catalog_X(); 
    File.read_array(filedm.string(), dm_field); // Reading the binary file

    PowerSpectrumF power(params);
    string file_cat_new=params._Output_directory()+"diluted_cat_bias.txt";

    vector<real_prec> bias_field; // Define container to allocate the bias on a mesh
    vector<real_prec> tr_field_counts; // Cointainer for counts
    tr_field_counts.resize(params._NGRID(), 0); // Cointainer for counts
    htools.get_density_field_grid(_COUNTS_,tr_field_counts);  // Computing counts

    if(params._get_cwc_properties())
    {
      vector<real_prec> dm_alpha(params._NGRID(), 0); //Container for TA computed from the tidal field of dm
      vector<real_prec> tr_alpha(params._NGRID(), 0); //Container for TA computed from the tidal field of dm

      Cwclass cwclass_dm(params); //CWC type for dark matter
      cwclass_dm.get_tidal_anisotropy(dm_field, dm_alpha); //Compute tidal anisotropy from dm
      htools.assign_cwc_to_tracers(cwclass_dm.CWClass);
      htools.assign_tidal_anisotropy_to_tracers(dm_alpha, true); // 1 indic

      Cwclass cwclass_tr(params); //CWC type for tracers
      cwclass_tr.get_tidal_anisotropy(tr_field_counts, tr_alpha); //Compute tidal anisotropy from dm
      htools.assign_tidal_anisotropy_to_tracers(tr_alpha, false);//  0 for dm
      htools.assign_idgrid_to_tracers();
    }

    if(params._assign_bias_to_full_sample())
     {
       power.object_by_object_bias(cat, dm_field);

       if(params._Get_tracer_bias_squared())
        power.object_by_object_bias_squared(cat, dm_field, tr_field_counts); // to compute b² for each tracer

       if(params._Get_tracer_bias_multipoles())
         power.object_by_object_bias_lm(cat, dm_field, params._lmax_bias()); // to compute b² for each tracer

       bias_field.resize(params._NGRID(), 0);
       htools.get_density_field_grid(_BIAS_,tr_field_counts,bias_field); //Get halo bias averaged on a mesh.
       fs::path fileb = fs::path(params._Output_directory()) / "ls_bias"; //Defile output file
       cat.clear_mem();// Release memmory
       File.write_array(fileb.string() ,bias_field); //Write the bias on the mesh to poutput file (binary)
       fileb.replace_extension(".dat");
       bias_field.clear();bias_field.shrink_to_fit();// Release memmory

       string json_file_plotf ="plot_file_bias_field.json";
       std::ofstream jfilef(json_file_plotf);
       json ja;
       ja["output_file"] = fileb.string();
       ja["show_bias_field"] = true;
          ja["Lbox"] = params._Lbox();
          ja["Nft"] = params._Nft();
          ja["sample"] = params._Name_survey();
          ja["name"] = "Halo effective bias";
          ja["Initial_slice"] = static_cast<int>(floor(params._Nft()/2.));
          ja["Nslices"] = static_cast<int>(floor(params._Nft()/10));
          jfilef<<ja.dump(4);
          jfilef.close();
          system("python3 ../python/cosmolib_plots.py plot_file_bias_field.json &");
        }
      else
        {

          htools.select_random_subsample(params._fraction_dilute());      // Dilute the sample: this generates an object of type Catalogue, called catalogue_random_subsample
          htools.catalogue_random_subsample.set_params(params);           // Feed the classs catalog_random_subsample with params

          power.object_by_object_bias(htools.catalogue_random_subsample, dm_field);

          if (params._Get_tracer_relative_bias())
            So.message_warning("The relative bias still to be implemented in this example");
          
          if(params._Get_tracer_bias_squared())
              power.object_by_object_bias_squared(htools.catalogue_random_subsample, dm_field, tr_field_counts); // to compute b² for each tracer

          if(params._Get_tracer_bias_multipoles())
              power.object_by_object_bias_lm(htools.catalogue_random_subsample, dm_field, params._lmax_bias()); // to compute b² for each tracer
  
          if(params._print_catalogue())
              htools.catalogue_random_subsample.print_catalogue(file_cat_new); 
         }

      htools.catalogue_random_subsample.clear_mem();// Release memmory

      if(params._i_mass_g()<0)
        {
          So.message_warning("Input catalog does not contain mass information");
          exit(1);
        }
 
      string json_file_plot ="plot_file_individual_bias.json";
      std::ofstream jfile(json_file_plot);
      json j;
      j["output_new_catalog"] = file_cat_new;
      j["show_individual_halo_mass_bias"] = true;
      j["column_mass"] = 0;
      j["column_bias"] = 1;
      jfile<<j.dump(4);
      jfile.close();
      system("python3 ../Python/cosmolib_plots.py plot_file_individual_bias.json &");
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

  void bias_prop(Params &params)
  {
      FileOutput File;
      ScreenOutput So;

      int Nqrts = 4 + 1; // Number of qaurtiles plus the first
      So.welcome_message_c();
  
      Catalogue cat(params, "TRACER");
      cat.read_catalog_new(params._Input_dir_cat() + params._file_catalogue());
  
      HaloTools htools(params, cat); 
      vector<double> dm_field(params._NGRID(), 0);
      File.read_array_t<double>(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field);
  
      vector<real_prec> dm_field_f(params._NGRID(), 0);

      for (size_t i = 0; i < dm_field.size(); i++)
        dm_field_f[i] = static_cast<real_prec>(dm_field[i]);
  
      dm_field.clear();
      dm_field.shrink_to_fit();
  
      PowerSpectrumF power(params);
      power.object_by_object_bias(cat, dm_field_f);

      string primary_prop = "_MASS_";
      vector<s_info_in_bins> secondary_bias_info_ref(Nqrts); // container for seconday bias: 2  the number of quartiles to use

      int Nbins = N_BINS_BIAS;
      s_info_in_bins bias_info;
      bias_info.allocate_all(Nbins);
      bias_info.name_info = primary_prop;

      htools.get_mean_bias_relation(bias_info); // Get bias_prop relation from catalogue

      string file_o = params._Output_directory() + "halo_bias.txt";
      ofstream bal;
      So.message_screen("Writing bias in file", file_o);
      bal.open(file_o.c_str());
      for (auto i = 0; i < Nbins; ++i)
        if (bias_info.vq1[i] > 0)
          bal << bias_info.vbin[i] << "\t" << bias_info.vq1[i] << "\t" << bias_info.vq3[i] << endl;
      bal.close();

      string secondary_prop = "_CONCENTRATION_";
      vector<s_info_in_bins> secondary_bias_info(Nqrts); // container for seconday bias: 2  the number of quartiles to use
      for (int i = 0; i < secondary_bias_info.size(); ++i)
        secondary_bias_info[i].allocate_all(Nbins);

      for (int i = 0; i < secondary_bias_info.size(); ++i)
        secondary_bias_info[i].name_info = primary_prop;

      for (int i = 0; i < secondary_bias_info.size(); ++i)
        secondary_bias_info[i].name_info_sec = secondary_prop;

      htools.get_mean_secondary_bias_relation(secondary_bias_info); // Get bias_prop relation from reference

      file_o = params._Output_directory() + "halo_bias_" + secondary_prop + ".txt";

      for (int iq = 1; iq < Nqrts; iq++) // loop over quartiles
      {
        string bbfile = file_o + "_quartile" + to_string(iq);
        ofstream sbal;
        So.message_screen("Writing bias in file", bbfile);
        sbal.open(bbfile.c_str());
        for (int i = 0; i < Nbins; ++i)
          if (secondary_bias_info[iq].vq1[i] > 0)
            sbal << secondary_bias_info[0].vbin[i] << "\t" << secondary_bias_info[iq].vq1[i] << "\t" << secondary_bias_info[iq].vq3[i] << "\t" << secondary_bias_info[iq].vq1[i] << "\t" << secondary_bias_info[iq].vq3[i] << endl;
        sbal.close();
      }



  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  ScreenOutput So(start_all, "logfile_bmt.log");

  if (getenv("RUNNING_IN_XTERM") == nullptr) {

        std::string cmd = "RUNNING_IN_XTERM=1 xterm -hold -e ";

        // Add program name
        cmd += argv[0];

        // Add all original arguments
        for (int i = 1; i < argc; ++i) {
            cmd += " ";
            cmd += argv[i];
        }

        system(cmd.c_str());
        return 0;
    }

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  string par_file = argv[2];
  Params params(par_file);
 
  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"bmt", required_argument, 0, 'b'},
        {"ibas",  required_argument, 0, 'i'},
        {"pbias",  required_argument, 0, 's'},
        {"sbias",  required_argument, 0, 'x'},
        {0, 0, 0, 0}
    };
  

  if(!params.input_sections.DarkMatterCatalogue)
    {
     throw std::runtime_error("Section DarkMatterCatalogue is not enabled in parameter file"); 
     exit(0);
  }

  if(!params.input_sections.TracerCatalogue)
    {
     throw std::runtime_error("Section TracerCatalogue is not enabled in parameter file"); 
     exit(0);
  }

  if(!params.input_sections.BiasAnalysis)
    {
     throw std::runtime_error("Section BiasAnalysis is not enabled in parameter file"); 
     exit(0);
  }

  while ((c = getopt_long(argc, argv, "b:i:",long_options, &option_index)) != -1) 
   {
      switch (c)
       {

         case 'b': // Run BiasMT
          {
            So.set_params(params);
            BiasMT BiasMT(params,true);
            BiasMT.set_So(So);
            BiasMT.execute();
          }

        case 'i': // assign individual bias to input tracer
          {
              ind_bias(params);
              break;
          }

        case 's': // Compute bias as a funciton of intrinsic properties
          {
            bias_prop(params);
            break;
          }

        case 'x':  // SECONDARY BIAS ANALYSIS
          {
            So.welcome_message_c();
            params.set_mass_assignment_scheme("CIC");
            PowerSpectrumF power(params);
            // power.halo_bias_analysis("redshit_space");// this argument is overriden by parameter use_real_and_redshift_space, read from parameter file
            break;
          }

     default:
        break;

      }
  }

  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
