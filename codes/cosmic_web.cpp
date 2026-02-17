////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
#include "../headers/def.h"
#include "CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_hist(Params params, vector<real_prec>&dm_alpha, vector<real_prec>&tr_alpha, vector<real_prec>&alpha_bias, vector<ULONG>&cwc, bool fitype)
{

  ScreenOutput So;
  int nb=200;
  FileOutput Fo;
  string ftype = fitype == 0 ? "dm" : "tr";
  real_prec bmin=-2; real_prec bmax=2; real_prec delta=(bmax-bmin)/static_cast<real_prec>(nb);
  vector<real_prec>ibins(nb,0); 

#ifdef _USE_OMP_ 
#pragma omp parallel for 
#endif
  for(int i=0; i<nb; ++i)
    ibins[i]=bmin+(i+0.5)*delta;

   for (int icwt=0; icwt<=4;++icwt)
    {
      vector<real_prec>hist(nb,0); 
      vector<real_prec>hist_alphah(nb,0); 
      vector<real_prec>hist_alphadm(nb,0); 

       if(icwt==0)
          {
             for(ULONG i=0;i<dm_alpha.size();++i)
               {
                 int ibin=get_bin(log10(alpha_bias[i]),bmin,nb,delta,0);
                 if(log10(alpha_bias[i])<bmax && log10(alpha_bias[i])>=bmin)
                    hist[ibin]+=1./static_cast<real_prec>(dm_alpha.size());

                  int abin=get_bin(log10(tr_alpha[i]),bmin,nb,delta,0);
                  if(log10(tr_alpha[i])<bmax && log10(tr_alpha[i])>=bmin)
                      hist_alphah[abin]+=1./static_cast<real_prec>(dm_alpha.size());
                      
                  int dbin=get_bin(log10(dm_alpha[i]),bmin,nb,delta,0);
                  if(log10(dm_alpha[i])<bmax && log10(dm_alpha[i])>=bmin)
                      hist_alphadm[dbin]+=1./static_cast<real_prec>(dm_alpha.size());
                  }
             }
          else
          {
            for(ULONG i=0;i<dm_alpha.size();++i)
            {
              if(cwc[i]==icwt)
                {
                  int ibin=get_bin(log10(alpha_bias[i]),bmin,nb,delta,0);
                  if(log10(alpha_bias[i])<bmax && log10(alpha_bias[i])>=bmin)
                    hist[ibin]+=1./static_cast<real_prec>(dm_alpha.size());

                  int abin=get_bin(log10(tr_alpha[i]),bmin,nb,delta,0);
                  if(log10(tr_alpha[i])<bmax && log10(tr_alpha[i])>=bmin)
                    hist_alphah[abin]+=1./static_cast<real_prec>(dm_alpha.size());

                  int dbin=get_bin(log10(dm_alpha[i]),bmin,nb,delta,0);
                  if(log10(dm_alpha[i])<bmax && log10(dm_alpha[i])>=bmin)
                    hist_alphadm[dbin]+=1./static_cast<real_prec>(dm_alpha.size());
                   
                }
            }
        }
        So.message_screen("Writting ascii file for histogram of bias_alpha. CWC is done accorging to ", ftype);
        string fileo=params._Output_directory()+"alpha_hbias_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist);

        So.message_screen("Writting ascii file for histogram of tidal anisotroipy based on tracers. CWC is done accorging to ", ftype);
        fileo=params._Output_directory()+"alphah_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist_alphah);

        So.message_screen("Writting ascii file for histogram of tidal anisotroipy based on dm. CWC is done accorging to ", ftype);
        fileo=params._Output_directory()+"alphadm_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist_alphadm);
      }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void s(string s)
{
  std::cout<<BOLDGREEN;  
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\tCOSMICATLASS"<<endl;
  std::cout<<"\t\n\tCosmologicalCATalogs for LArge Scale Structure"<<endl;
  std::cout<<"\t\n\tHow to run: "<<s<<"\t [-option] [argument]"<<endl;
  std::cout<<"\t\n\tOptions: "<<endl;
  std::cout<<"\t         -a for information on the author (no argument). "<<endl;
  std::cout<<"\t         -b parameter_file.ini (with argument) To perform cosmic-web analysis and bias assignment"<<endl;
  std::cout<<"\t         -h parameter_file.ini (no argum,ent): Help"<<endl;
  std::cout<<"\t         -i parameter_file.ini (with argument): Shows input pars"<<endl;
  std::cout<<"\t\n\tArgument is a parameter file"<<endl;
  std::cout<<"\tConsult ../Headers/def.h for pre-procesor directives"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<RESET<<endl;
  std::cout<<RESET<<endl;                                       
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  auto start_time = chrono::high_resolution_clock::now();
  time_t start_all;
  time(&start_all);

  char temp;
  string par_file;

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec. and compile again.");
#endif
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile = "logfile.log";
  ScreenOutput So(start_all, logfile);
  if (argc == 1)
  {
    So.usage(argv[0]);
    exit(1);
  }

  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"ibias", required_argument, 0, 'i'},
        {"cwc", required_argument, 0, 'c'},
        {0, 0, 0, 0}
    };


  FileOutput File;     // File management
  So.welcome_message_c();
  par_file = argv[2];  // Get parameter file
  Params params(par_file);  // Read parameter file


  while ((temp = getopt(argc, argv, "i:c:")) != -1)
  {
    // -------------------------------------------------------------------------------------------------
    if ( temp  == 'i') // assign individual bias to input tracer
    {

      Catalog cat(params); // Feed params into the Catlaog class and define an object of that type
      cat.read_catalog(params._Input_dir_cat() + params._file_catalogue(), 0); // Read the asciii file

      vector<real_prec> dm_field(params._NGRID(), 0); //COntainer for dm
      File.read_array(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field); // Reading the binary file

      PowerSpectrumF power(params);
      string file_cat_new=params._Output_directory()+"diluted_cat_bias.txt";

      vector<real_prec> bias_field; // Define container to allocate the bias on a mesh
      vector<real_prec> tr_field_counts; // Cointainer for counts

      if(true==params._assign_bias_to_full_sample())
      {
        tr_field_counts.resize(params._NGRID(), 0); // Cointainer for counts
        cat.get_density_field_grid(_COUNTS_,tr_field_counts);  // Computing counts
        power.object_by_object_bias(cat.Halo, dm_field);
        if(true==params._Get_tracer_bias_squared())
          power.object_by_object_bias_squared(cat.Halo, dm_field, tr_field_counts); // to compute b² for each tracer
        if(true==params._Get_tracer_bias_multipoles())
          power.object_by_object_bias_lm(cat.Halo, dm_field, params._lmax_bias()); // to compute b² for each tracer


        bias_field.resize(params._NGRID(), 0);
        cat.get_density_field_grid(_BIAS_,tr_field_counts,bias_field); //Get halo bias averaged on a mesh.
        string fileb=params._Output_directory()+"ls_bias"; //Defile output file
        File.write_array(fileb,bias_field); //Write the bias on the mesh to poutput file (binary)
        bias_field.clear();bias_field.shrink_to_fit();// Release memmory
        string json_file_plotf ="plot_file_bias_field.json";
        std::ofstream jfilef(json_file_plotf);
        json ja;
        ja["output_file"] = fileb+".dat";
        ja["show_bias_field"] = true;
        ja["Lbox"] = params._Lbox();
        ja["Nft"] = params._Nft();
        ja["sample"] = params._Name_survey();
        ja["name"] = "Halo effective bias";
        ja["Initial_slice"] = static_cast<int>(floor(params._Nft()/2.));
        ja["Nslices"] = 70;
        jfilef<<ja.dump(4);
        jfilef.close();
        system("python3 ../python/cosmolib_plots.py plot_file_bias_field.json &");

      }
      else
      {
        cat.select_random_subsample(params._fraction_dilute());      // Dilute the sample:
        power.object_by_object_bias(cat.Halo_random_subsample, dm_field);

       if(true==params._Get_tracer_bias_squared())
          power.object_by_object_bias_squared(cat.Halo_random_subsample, dm_field, tr_field_counts); // to compute b² for each tracer

        if(true==params._Get_tracer_bias_multipoles())
          power.object_by_object_bias_lm(cat.Halo_random_subsample, dm_field, params._lmax_bias()); // to compute b² for each tracer

        So.message_warning("The function print catalog must be made automatic, not to be hard-coded.");
        print_catalog(cat.Halo_random_subsample, file_cat_new, false); // set true of blm are to be written, 
      }


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
      system("python3 ../python/cosmolib_plots.py plot_file_individual_bias.json &");


      // -------------------------------------------------------------------------------------------------
      // -------------------------------------------------------------------------------------------------

    }
    else if ( temp  == 'c') // Cosmic web analysis
    {


      //Read dm field
      vector<real_prec> dm_field(params._NGRID(), 0); //COntainer for dm
      File.read_array(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field); // Reading the binary file


      // Read catalog
      Catalog cat(params); // Feed params into the Catlaog class and define an object of that type
      cat.read_catalog(params._Input_dir_cat() + params._file_catalogue(), 0); // Read the asciii file

      //----------------------------------------------
      // Get density field from tracers
      //----------------------------------------------
 
      bool get_cwc=false;
      if(false==get_cwc){
        So.message_screen("============================================================");
        So.message_screen("NO CWC quantities are to be computed. Check line", __LINE__);
        So.message_screen("============================================================");
      }

        vector<real_prec> tr_field_counts; // Cointainer for counts
//      if(get_cwc)//commented for the test of computing b²
      {
        tr_field_counts.resize(params._NGRID(), 0); // Cointainer for counts
        cat.get_density_field_grid(_COUNTS_,tr_field_counts);  // Computing counts
      }  

 
 
      //-----------------------------
      Cwclass cwclass(params); // for halos
      vector<real_prec> tr_alpha; // Container for TA computed from the halo tidal field

      if(get_cwc)
      {
 
      // Define Cwc cobject to perform cosmic-web analysis based on the halo density field.  
      // Get tidal anisotropy (TA hereafter) from tracers
      tr_alpha.resize(params._NGRID(), 0); // Container for TA computed from the halo tidal field
      cwclass.get_tidal_anisotropy(tr_field_counts, tr_alpha); // Compute the TA from halo field --tr_field-- and allocate in tr_alpha

      So.message_screen("Writting binary file tidal anisotropy based in tracer tidal filed in the mesh");
      string fileo=params._Output_directory()+"alpha_tr";
      File.write_array(fileo,tr_alpha); //Write tidal anisotropy based on the tracers

      So.message_screen("Writting binary file for lambda_1 from tracer filed in the mesh");
      fileo=params._Output_directory()+"l1_tr";
      File.write_array(fileo,cwclass.lambda1); //Write first eigenvalue of halo tidal field

      So.message_screen("Writting binary file for lambda_2 from tracer filed in the mesh");
      fileo=params._Output_directory()+"l2_tr";
      File.write_array(fileo,cwclass.lambda2); //Write second eigenvalue of halo tidal field

      So.message_screen("Writting binary file for lambda_3 from tracer filed in the mesh");
      fileo=params._Output_directory()+"l3_tr";
      File.write_array(fileo,cwclass.lambda3); //Write third eigenvalue of halo tidal field
      
      So.message_screen("Writting binary file the cosmic-web classificcation based on tracers in the mesh");
      fileo=params._Output_directory()+"cwc_tr";
      File.write_array(fileo,cwclass.CWClass);  //Write the cwc of each cell, computed in from tracers in CWClass.
      }

      //-----------------------------
      //-----------------------------
      // Get tidal anisotropy from dm

      vector<real_prec> dm_alpha(params._NGRID(), 0); //Coitainer for TA computed from the tidal field of dm
      Cwclass cwclass_dm(params); //CWC type for dark matter
      if(get_cwc)
      {
        cwclass_dm.get_tidal_anisotropy(dm_field, dm_alpha); //Compute tidal anisotropy from dm
        string fileo=params._Output_directory()+"alpha_dm";
        File.write_array(fileo,dm_alpha); // Write TA from dm to binary
      }
      // --------------------------------------------------
      // Assign ID and cwc to halos in cat.Halo
      cat.assign_idgrid_to_tracers();
      if(get_cwc)
      {
      cat.assign_cwc_to_tracers(cwclass_dm.CWClass);
      cat.assign_tidal_anisotropy_to_tracers(tr_alpha, false);
      cat.assign_tidal_anisotropy_to_tracers(dm_alpha, true);
      }
      // --------------------------------------------------
      //Assign bias to tracers and get the bias field:
     // Define an object of type Power Spectrum, passing the parametes through the class variable param
      PowerSpectrumF power(params);
//      string file_cat_new=params._Output_directory()+"diluted_cat_blm.txt";
//      string file_cat_new=params._Output_directory()+"diluted_cat_blm_bsquared.txt";
      string file_cat_new=params._Output_directory()+"diluted_cat_bias.txt";
      if(true==params._assign_bias_to_full_sample())
      {
        power.object_by_object_bias(cat.Halo, dm_field);
        power.object_by_object_bias_lm(cat.Halo, dm_field, params._lmax_bias());
        cat.select_random_subsample_bl(params._fraction_dilute(), file_cat_new);
      }
      else{
      // Dilute the sample:
        cat.select_random_subsample(params._fraction_dilute());

        //Assign individual bias to objects in cat.Halo. That bias uses the dm_field:
        power.object_by_object_bias(cat.Halo_random_subsample, dm_field);

        //Assign individual bias^2 to objects in cat.Halo. That bias uses the dm_field:
        if(true==params._Get_tracer_bias_squared())
          power.object_by_object_bias_squared(  cat.Halo_random_subsample, dm_field, tr_field_counts); // to compute b² for each tracer

        //Assign individual harmonic-based bias to objects in cat.Halo. That bias uses the dm_field:
        if(true==params._Get_tracer_bias_multipoles())
            power.object_by_object_bias_lm(cat.Halo_random_subsample, dm_field, params._lmax_bias());  
        
            // Print catalog: do it below, after cwc is done
        print_catalog(cat.Halo_random_subsample, file_cat_new, true); // set true of blm are to be written, 
      }

      vector<real_prec> bias_field; // Define container to allocate the bias on a mesh
      if(true==params._assign_bias_to_full_sample())
      {
        bias_field.resize(params._NGRID(), 0);
        cat.get_density_field_grid(_COUNTS_,tr_field_counts,bias_field); //Get halo bias averaged on a mesh.
        string fileb=params._Output_directory()+"ls_bias"; //Defile output file
        File.write_array(fileb,bias_field); //Write the bias on the mesh to poutput file (binary)
        bias_field.clear();bias_field.shrink_to_fit();// Release memmory
      }
  // --------------------------------------------------

    if(get_cwc)
      {

      vector<real_prec> alpha_bias(params._NGRID(), 0); // Container for the bias--ratio-- of TA

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<dm_field.size();++i) // get the ratio between TA from tr and TA from dm
      alpha_bias[i]=tr_alpha[i]/dm_alpha[i];

      So.message_screen("Writting binary file for alpha_tr / alpha _dm in the mesh");
      string fileo=params._Output_directory()+"bias_alpha";
      File.write_array(fileo,alpha_bias); //Write the ratio of TA's

      So.message_screen("Writting binary file for eigenvalue lambda_1_dm in the mesh");
      fileo=params._Output_directory()+"l1_dm";
      File.write_array(fileo,cwclass_dm.lambda1); //Write first eigenvalue of DM tidal field

      So.message_screen("Writting binary file for eigenvalue lambda_2_dm in the mesh");
      fileo=params._Output_directory()+"l2_dm";
      File.write_array(fileo,cwclass_dm.lambda2);//Write second eigenvalue of DM tidal field

      So.message_screen("Writting binary file for eigenvalue lambda_3_dm in the mesh");
      fileo=params._Output_directory()+"l3_dm";
      File.write_array(fileo,cwclass_dm.lambda3);//Write third eigenvalue of DM tidal field

      So.message_screen("Writting binary file the cosmic-web classificcation based on dm in the mesh");
      fileo=params._Output_directory()+"cwc_dm";
      File.write_array(fileo,cwclass_dm.CWClass);  //Write the cwc of each cell, computed in from tracers in CWClass.


      // Histograms here are done classifying CW according to dm.
      get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass_dm.CWClass, 0);

      // Histograms here are done classifying CW according to tracers.
      get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass.CWClass, 1);

      //Relesas memmory (not requested as memmory is clean when this if ends)   
     dm_alpha.clear(); dm_alpha.shrink_to_fit();
     dm_field.clear(); dm_field.shrink_to_fit();
     tr_field_counts.clear(); tr_field_counts.shrink_to_fit();

      }
     // End of the '-b' option

    }
    else if ('?' == temp)
    {
      cout << endl;
      cout << "Argument not recognized." << endl;
      cout << "Please run ./cosmicatlass.exe -h." << endl;
      exit(1);
    }
  }

 auto end_time = chrono::high_resolution_clock::now();
 real_prec difft=chrono::duration<real_prec>(end_time-start_time).count();
 cout<<YELLOW<<"Time "<<difft<<" secs "<<RESET<<endl;

  exit(0);
  return 0;
}
// #################################################################################
// ##################################################################################
// #################################################################################
