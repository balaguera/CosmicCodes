////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
#include "../include/def.hpp"
#include "../include/CosmiCalcLIB.hpp"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_hist(Params params, vector<real_prec>&dm_alpha, vector<real_prec>&tr_alpha, vector<real_prec>&alpha_bias, vector<WebType>&cwc, bool fitype)
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
              if(static_cast<int>(cwc[i])==icwt)
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


void cw_analysis(Params &params)
{
  FileOutput File;
  ScreenOutput So;
  //Read dm field


  vector<real_prec> dm_field(params._NGRID(), 0); //Container for dm
  File.read_array(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field); // Reading the binary file

          // Read catalog
  Catalogue cat(params, "TRACER");
  cat.read_catalog_new(params._Input_dir_cat() + params._file_catalogue()); // Read the asciii file
  HaloTools htools(params, cat); // Feed params into the Catlaog class and define an object of that type

          // Get density field from tracers

  So.message_screen("============================================================");
  So.message_screen("NO CWC quantities are to be computed. Check line", __LINE__);
  So.message_screen("============================================================");

          vector<real_prec> tr_field_counts; // Cointainer for counts
          tr_field_counts.resize(params._NGRID(), 0); // Cointainer for counts
          htools.get_density_field_grid(_COUNTS_,tr_field_counts);  // Computing counts

          Cwclass cwclass(params); // for halos
          vector<real_prec> tr_alpha; // Container for TA computed from the halo tidal field

          // Define Cwc cobject to perform cosmic-web analysis based on the halo density field.  
          // Get tidal anisotropy (TA hereafter) from tracers
          tr_alpha.resize(params._NGRID(), 0); // Container for TA computed from the halo tidal field
          cwclass.get_tidal_anisotropy(tr_field_counts, tr_alpha); // Compute the TA from halo field --tr_field-- and allocate in tr_alpha
        
          const fs::path out_dir = params._Output_directory();

          std::vector<std::pair<std::string, std::vector<real_prec>&>> tracer_outputs = { //Using strucure bending
              {"alpha_tr", tr_alpha},
              {"l1_tr",    cwclass.lambda1},
              {"l2_tr",    cwclass.lambda2},
              {"l3_tr",    cwclass.lambda3}
          };

         for (const auto& [fname, data] : tracer_outputs)
          {
            So.message_screen("Writing binary file " + fname + " in the mesh");
            fs::path fileo = out_dir / fname;
            File.write_array(fileo.string(), data);
          }
          So.message_screen("Writing binary file for cosmic-web classification based on DM in the mesh");
          fs::path cwc_file = out_dir / "cwc_tr";
          File.write_array(cwc_file.string(), cwclass.CWClass);
 
          // Get tidal anisotropy from dm
          vector<real_prec> dm_alpha(params._NGRID(), 0); //Coitainer for TA computed from the tidal field of dm
          Cwclass cwclass_dm(params); //CWC type for dark matter
          cwclass_dm.get_tidal_anisotropy(dm_field, dm_alpha); //Compute tidal anisotropy from dm
          fs::path fileo=fs::path(params._Output_directory()) / "alpha_dm";
          File.write_array(fileo.string(),dm_alpha); // Write TA from dm to binary

          // Assign ID and cwc to halos in cat.Halo
          htools.assign_idgrid_to_tracers();
          htools.assign_cwc_to_tracers(cwclass_dm.CWClass);
          htools.assign_tidal_anisotropy_to_tracers(tr_alpha, false);
          htools.assign_tidal_anisotropy_to_tracers(dm_alpha, true);

          vector<real_prec> alpha_bias(params._NGRID(), 0); // Container for the bias--ratio-- of TA
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(size_t i=0;i<alpha_bias.size();++i) // get the ratio between TA from tr and TA from dm
              alpha_bias[i]=tr_alpha[i]/dm_alpha[i];

          So.message_screen("Writting binary file for alpha_tr / alpha _dm in the mesh");
          fileo= fs::path(params._Output_directory()) / "bias_alpha";
          File.write_array(fileo.string(),alpha_bias); //Write the ratio of TA's
          

              // Lambda eigenvalues
          std::vector<std::pair<std::string, std::vector<real_prec>&>> eigenfiles = {
               {"l1_dm", cwclass_dm.lambda1},
               {"l2_dm", cwclass_dm.lambda2},
               {"l3_dm", cwclass_dm.lambda3}
          };
              // Write each eigenvalue
          for (const auto& [fname, data] : eigenfiles)
            {
              So.message_screen("Writing binary file for eigenvalue " + fname + " in the mesh");
              fs:: path file_path = out_dir / fname;
              File.write_array(file_path.string(), data);
            }

              // Write cosmic web classification
          So.message_screen("Writing binary file for cosmic-web classification based on DM in the mesh");
          cwc_file = out_dir / "cwc_dm";
          File.write_array(cwc_file.string(), cwclass_dm.CWClass);

          // Histograms here are done classifying CW according to dm.
          get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass_dm.CWClass, 0);

              // Histograms here are done classifying CW according to tracers.
          get_hist(params, dm_alpha, tr_alpha, alpha_bias,cwclass.CWClass, 1);

}




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


  time_t start_all;
  time(&start_all);
  ScreenOutput So(start_all, "logfile_cw.log");

 if(argc==1)
  {
    So.usage(argv[0]);
    exit(1);
  }

  string  par_file = argv[2];
  So.welcome_message_c();
  Params params(par_file);
  if (getenv("RUNNING_IN_XTERM") == nullptr) {
    std::string cmd = "RUNNING_IN_XTERM=1 gnome-terminal -hold -e ";
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

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec. and compile again.");
#endif

  if(false==params.input_sections.CWCAnalysis)
    {
     throw std::runtime_error("Section BiasAnalysis is not enabled in parameter file"); 
     exit(0);
    }

  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"tcw", required_argument, 0, 't'},
        {"vcw", required_argument, 0, 'v'},
        {0, 0, 0, 0}
    };

 
  while ((c = getopt_long(argc, argv, "t:v:",long_options, &option_index)) != -1) 
   {
      switch (c) {

        case 't': // Cosmic web analysis
        {
          cw_analysis(params);
          break;
      }
      case 'v':
      {

        break;
      }
      default:
        break;
 }

}

  exit(0);
  return 0;
}
// #################################################################################
// ##################################################################################
// #################################################################################
