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
        string fileo=params._Output_directory()+"alpha_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist);

        fileo=params._Output_directory()+"alphah_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist_alphah);

        fileo=params._Output_directory()+"alphadm_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+"_"+ftype+".txt";
        Fo.write_to_file(fileo,ibins,hist_alphadm);
      }



}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////





////////////////////////////////////////////////////////////////////////////
void usage(string s)
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
int main(int argc, char *argv[])
{
  //  system("clear");
  time_t start_all;
  time(&start_all);
  char temp;
  string par_file;

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec evecd da NaAaAAAArywhere. Plese desable and compile again.");
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
  while ((temp = getopt(argc, argv, "hadi:b:")) != -1)
  {
    if (temp == 'h')
      usage(argv[0]);

    else if (temp == 'a') // Show authors
      So.author();

    else if (temp == 'i') // displays input parameters
    {
      par_file = argv[2];
      Params params(par_file);
      params.show_params();
    }
    else if (temp == 'd') // displays input parameters
    {
      So.show_preproc();
    }
    else if ( temp  == 'b') // assign bias and do CW analysus using tracers and DM.
    {


      FileOutput File;     // File management
      So.welcome_message_c();
      par_file = argv[2];  // Get parameter file
      Params params(par_file);  // Read parameter file

      //Read dm field
      vector<real_prec> dm_field(params._NGRID(), 0); //COntainer for dm
      File.read_array(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field); // Readin the bunary file


      // Read catalog
      Catalog cat(params); // Feed params into the Catlaog class and define an object of that type
      cat.read_catalog(params._Input_dir_cat() + params._file_catalogue(), 0); // Read the asciii file

      //----------------------------------------------
      // Get density field from tracers
      //----------------------------------------------
 
      vector<real_prec> tr_field_counts(params._NGRID(), 0); // Cointainer for counts
      cat.get_density_field_grid(_COUNTS_,tr_field_counts);  // Computing counts

      //Assign bias to tracers and get the bias field:

      // Definr an object of type Power Spectrum, passing the parametes through the class variable param
      PowerSpectrumF power(params);

      //Assign individual bias to objhects in cat.Halo. That bias uses the dm_field.
      power.object_by_object_bias(cat.Halo, dm_field);

      vector<real_prec> bias_field(params._NGRID(), 0); // Define container to allocate the bias on a mesh
       
      cat.get_density_field_grid(_COUNTS_,tr_field_counts,bias_field); //Get halo bias averaged on a mesh.

      string fileb=params._Output_directory()+"ls_bias_p"+to_string(params._unitsim_plabel()); //Defile output file

      File.write_array(fileb,bias_field); //Write the bias on the mesh to poutput file (binary)

      bias_field.clear();bias_field.shrink_to_fit();// Release memmory
 

      vector<real_prec> tr_field(params._NGRID(), 0); //container for
      cat.get_density_field_grid(_iBIAS_,tr_field);  
      cat.get_density_field_grid(_COUNTS_,tr_field);
      string fileo=params._Output_directory()+"field_hbias_p"+to_string(params._unitsim_plabel());

      File.write_array(fileo,tr_field);


      //-----------------------------
      // Define Cwc cobject to perform cosmic-web analysis based on the halo density field.  
      Cwclass cwclass(params); // for halos

      // Get tidal anisotropy (TA hereafter) from tracers
      vector<real_prec> tr_alpha(params._NGRID(), 0); // Container for TA computed from the halo tidal field
      cwclass.get_tidal_anisotropy(tr_field, tr_alpha); // Compute the TA halo field --tr_field-- and allocate in tr_alpha

      fileo=params._Output_directory()+"alpha_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,tr_alpha); //Write tidal anisotropy

      fileo=params._Output_directory()+"l1_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda1); //Write first eigenvalue of halo tidal field


      fileo=params._Output_directory()+"l2_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda2); //Write second eigenvalue of halo tidal field

      fileo=params._Output_directory()+"l3_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda3); //Write third eigenvalue of halo tidal field

      fileo=params._Output_directory()+"cwc_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.CWClass);  //Write the cwc of each cell, computed in from tracers in CWClass.


      //-----------------------------
      //-----------------------------
      // Get tidal anisotropy from dm

      vector<real_prec> dm_alpha(params._NGRID(), 0); //Coitainer for TA computed from the tidal field of dm
      Cwclass cwclass_dm(params); //CWC type for dark matter

      cwclass_dm.get_tidal_anisotropy(dm_field, dm_alpha); //Compute tidal anisotropy from dm

      fileo=params._Output_directory()+"alpha_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,dm_alpha); // Write TA from dm to binary

      vector<real_prec> alpha_bias(params._NGRID(), 0); // Container for the bias--ratio-- of TA

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<dm_field.size();++i) // get the ratio between TA from tr and TA from dm
      alpha_bias[i]=tr_alpha[i]/dm_alpha[i];


      fileo=params._Output_directory()+"bias_alpha_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,alpha_bias); //Write the ratio of TA's

      fileo=params._Output_directory()+"l1_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda1); //Write first eigenvalue of DM tidal field

      fileo=params._Output_directory()+"l2_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda2);//Write second eigenvalue of DM tidal field

      fileo=params._Output_directory()+"l3_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda3);//Write third eigenvalue of DM tidal field

      fileo=params._Output_directory()+"cwc_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.CWClass);  //Write the cwc of each cell, computed in from tracers in CWClass.


      // Histograms here are done classifying CW according to dm.
      get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass_dm.CWClass, 0);

      // Histograms here are done classifying CW according to tracers.
      get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass.CWClass, 1);

     //Relesas memmory (not requested as memmory is clean when this if ends)   
     dm_alpha.clear(); dm_alpha.shrink_to_fit();
     dm_field.clear(); dm_field.shrink_to_fit();
     tr_field.clear(); tr_field.shrink_to_fit();
     
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
  exit(0);
  return 0;
}
// #################################################################################
// ##################################################################################
// #################################################################################
