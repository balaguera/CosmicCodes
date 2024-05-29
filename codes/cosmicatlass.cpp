////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
# include "CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
//  system("clear");
  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_BiasMT;

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec everywhere. Plese desable and compile again.");
#endif
  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);
  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hadi:b:w:z:p:c:g:u:f:s:m:e:q:v")) != -1)
    {
      if(temp=='h')
        So.usage(argv[0]);

      else if (temp=='a')  // Show authors
        So.author();

      else if(temp=='i') // displays input parameters
        {
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          params.show_params();
        }
      else if(temp=='d') // displays input parameters
        {
          So.show_preproc();
        }
      else if(temp=='b') // Run BiasMT
        {
#if !defined _POWER__
//          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params Par(par_file_BiasMT);
          So.set_params(Par);
          BiasMT BiasMT(Par,true);
          BiasMT.set_So(So);
          BiasMT.BiasMT_execute();
#else
          So.message_warning("_POWER_ is defined in def.h.");
         exit(0);
#endif
          }

      else if(temp=='p')   // to mesure Power spectrum
        {
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          PowerSpectrumF cPSF(params);
          if(params._measure_cross()==true && params._input_type()!="catalog")
           cPSF.compute_cross_power_spectrum_grid(true, params._delta_grid_file(), params._delta_grid_file2(), true);
          else if(params._measure_cross()==true && params._input_type()=="catalog")
            cPSF.compute_power_spectrum(false,false);
          else if(params._measure_cross()==false && ( params._input_type()=="density_grid" || params._input_type()=="delta_grid") )
             cPSF.compute_power_spectrum_grid();
          else if(params._measure_cross()==false && params._input_type()=="catalog")
            cPSF.compute_power_spectrum(true,false); // This is the default function
//          cPSF.compute_power_spectrum(); // Thi is made for test in jpas mesuring power in bins of z, c, and ms

#ifdef _USE_REDSHIFT_BINS_
          if(params._sys_of_coord_g()==0)
          {
              So.message_warning("Asking for redshift bins while providing cartessian coordinates. Check input parameter file and def.h");
              exit(0);
          }
#endif
      }
      else if(temp=='m')   // to mesure marked power or corr function
        {
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          PowerSpectrumF cPSF(params);
          cPSF.compute_marked_correlation_function();
#ifdef _USE_REDSHIFT_BINS_
          if(params._sys_of_coord_g()==0)
          {
              So.message_warning("Asking for redshift bins while providing cartessian coordinates. Check input parameter file and def.h");
              exit(0);
          }
#endif


          So.message_time(start_all);
        }
      else if(temp=='c')   // to read input tracer catalg and analyze it. 
        {

          So.message_screen("cosmicatlass running under the option -c");
#ifdef _REDSHIFT_SPACE_
          So.message_warning("Option _REDSHIFT_SPACE_ in def.h is enabled. If velocities are not to be used, please desable it to save space. Code ends here.");
          So.message_screen("After editing the def.h file, please do make clean; make -j BiasMT;");
//          exit(0);
#endif
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog catalog(params); // this is redundant
          catalog.analyze_cat(true);
          So.message_time(start_all);
        }

      else if(temp=='s')   // to read an input density field and reduce resilution using a low-pass filter
        {
          So.message_screen("Cosmicatlass running with option -s (low pass filter)");
          So.message(start_all);
          par_file_BiasMT = argv[2];
          FileOutput File;
          Params params(par_file_BiasMT);
          params.set_MAS_correction(false);
          ULONG Nft_LR=params._Nft();
          ULONG Nft_HR=params._Nft_HR();
          ULONG Ngrid_HR=params._NGRID_HR();
          // ************************************************** //
          //* Read hr and get power *//
          string input_f=params._Input_Directory_X()+params._Name_Catalog_X();
          vector<real_prec>in_f(Ngrid_HR,0);
          File.read_array(input_f,in_f);
          params.set_Nft(Nft_HR);
          params.set_Name_survey("HR");
          PowerSpectrumF power_f(params);
          power_f.compute_power_spectrum_grid(in_f, true);
          // ************************************************** //
          params.set_Nft(Nft_LR);
          vector<real_prec>out_f(params._NGRID(),0);
          vector<real_prec>out_av(params._NGRID(),0);

          if("delta_grid"==params._input_type())
          {
            real_prec mean_HR=static_cast<real_prec>(Ngrid_HR)/pow(params._Lbox(),3);
#pragma omp parallel for
            for(ULONG i=0; i<in_f.size();++i)
                in_f[i]=mean_HR*(1.0+in_f[i]);
            }
          low_pass_filter(Nft_HR,Nft_LR,params._masskernel(),false,in_f,out_f, params._Lbox()); // low-pass filter
          average(Nft_HR,Nft_LR, params._Lbox(), in_f,out_av); // sown-sampling in configuration space.


          in_f.clear();in_f.shrink_to_fit();
        // Get ready for low*reslution in configuration space
          string output_fav=params._Input_Directory_X_NEW()+params._Name_Catalog_X_NEW()+"_averaged";
          File.write_array(output_fav,out_av);
          out_f.clear();out_f.shrink_to_fit();
          // ********************************************//
          //* Power of the LR overdensity*//
          params.set_Name_survey("LR");
          PowerSpectrumF power(params);
          get_overdens(out_f,out_f);
          power.compute_power_spectrum_grid(out_f, true);
          string output_f=params._Input_Directory_X_NEW()+params._Name_Catalog_X_NEW();
          File.write_array(output_f,out_f);
          out_f.clear();out_f.shrink_to_fit();
          // ********************************************//
      }
      else if(temp=='f')   // convert the xyz binary files from fastpm to a density field. Mianly used for IC
        {
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog tracer(params);
          tracer.read_catalog_bin();
          }
      else if(temp=='g')   // Read halo catalog and generate galaxy catalogs using an HOD
        {

          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog tracer(params);
          tracer.halos2galaxies_HOD();
        }
      else if('u'==temp)  // SECONDARY BIAS
        {
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          params.set_mass_assignment_scheme("CIC");
	  
          // Go to power spectrum and read the catalog therefrom.
          PowerSpectrumF power(params);
          power.compute_power_spectrum_halo_bias("redshit_space");// this argument is overriden by parameter use_real:_and_redshift_space from parameter file
//          power.compute_power_spectrum("redshift");
        }
      else if('z'==temp)  // assign bias to tracers
        {
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog cat(params);
          cat.read_catalog_bin_tng();
          vector<real_prec>dm_field(params._NGRID(),0);
          File.read_array(params._Input_Directory_X()+params._Name_Catalog_X(), dm_field);
          PowerSpectrumF power(params);
          power.object_by_object_bias(cat.Halo,  dm_field);

      }
      else if('q'==temp)  // to draw warnings
        {
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          PowerSpectrumF cPSF(params);
          cPSF.GetSuperClusters("under", "_RELATIVE_BIAS_");
//          cPSF.GetSuperClusters("_RELATIVE_BIAS_");
//          cPSF.GetSuperClusters("_TIDAL_ANISOTROPY_");
        }

      else if('v'==temp)  // to draw warnings
        {
              So.show_warnings();
        }
      else if('w'==temp)  // to draw warnings 
        {
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          PowerSpectrumF cPSF(params);
          cPSF.compute_window_function();
        }

      else if('e'==temp)  // to generate GRF
        {
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          PowerSpectrumF cPSF(params);
          cPSF.get_GaussianRandomField();
        }

      else if('?'==temp)
       {
         cout<<endl;
         cout<<"Argument not recognized."<<endl;
         cout<<"Please run ./cosmicatlass.exe -h."<<endl;
         exit(1);
       }

    }
  exit(0);
  return 0;

}
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################


