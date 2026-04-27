////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
# include "../include/CosmiCalcLIB.hpp" 
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


  ScreenOutput So(start_all,"logfile_pk.log");

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  string par_file=argv[2];
  Params params(par_file);

  if (getenv("RUNNING_IN_XTERM") == nullptr) {
      std::string cmd = "RUNNING_IN_XTERM=1 xterm -bg black -fg white -hold -e ";
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
  // Verify that section in json file is active: 
  if(false==params.input_sections.FourierAnalysis)
    {
      throw std::runtime_error("Section Fourier Analysis is not enabled in parameter file"); 
      exit(0);
   }

   int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"power", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"lpfilter",  required_argument, 0, 's'},
        {"grf",      no_argument,       0, 'g'},
        {"m",      required_argument,       0, 'm'},
        {0, 0, 0, 0}
    };

  while ((c = getopt_long(argc, argv, "p:w:s:g:m:",long_options, &option_index)) != -1)
    {

      switch (c) {
 
      case 'p':   // to mesure Power spectrum
        {

          So.message_screen("Power Spectrum analysis");
          string json_file_plot ="plot_file_power_spectrum_measurement.json";
          std::ofstream jfile(json_file_plot);
          json j;

          //Hay que hacer un enum con CATALOGUE, DENSITY_GRID, DELTA_GRID

          PowerSpectrumF cPSF(params);
          //For input deltas to measure cross power spectrum:
          if(params._measure_cross() && params._input_type()!="catalog")
           cPSF.measure_cross_power_spectrum_grid(true, params._delta_grid_file(), params._delta_grid_file2(), true);

           //For input catalog to measure cross power spectrum:
//          else if(params._measure_cross() && params._input_type()=="catalog")
 //           cPSF.measure_power_spectrum(false); HAY QUE HACE RESTA BIEN

            // For input fields to measure power (auto power)
          else if(!params._measure_cross() && ( params._input_type()=="density_grid" || params._input_type()=="delta_grid") )
             cPSF.measure_power_spectrum_grid();

             // For input catalogs to measure auto power
          else if(!params._measure_cross() && params._input_type()=="catalog" && params._use_random_catalog())
            cPSF.measure_power_spectrum_data(); // This is the default function

          else if(!params._measure_cross() && params._input_type()=="catalog" && !params._use_random_catalog())
            cPSF.measure_power_spectrum_box(); // This is the default function


          if(params._show_power_spectrum())
          {
            j["header"]="k vec, p0, p2, p4, sigma, nmodes";
            j["power_spectrum_measurement"]=true;
            j["kmin"] = cPSF._kvector_data(0);
            j["kmax"] = 2.0;
            j["sample"] = params._Name_survey();
            j["redshift"] = params._redshift();
            j["statistics"] = params._statistics();
            j["output_file"] = cPSF._file_power();
            j["clustering_space"] = params._clustering_space();
              
            cout<<"Generating json file "<<json_file_plot<<" for plotting"<<endl;  
            jfile.close();
            system("python3 ../python/cosmolib_plots.py plot_file_power_spectrum_measurement.json &");
          }

          break;
        }

      case 's':   // to read an input density field and reduce resilution using a low-pass filter
       {
          So.message_screen("CosmiCalc running with option -s (low pass filter)");
          So.message(start_all);
          FileOutput File;
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
          power_f.measure_power_spectrum_grid(in_f, true);
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
          average(Nft_HR,Nft_LR, params._Lbox(), in_f,out_av); // down-sampling in configuration space.

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
          power.measure_power_spectrum_grid(out_f, true);
          string output_f=params._Input_Directory_X_NEW()+params._Name_Catalog_X_NEW();
          File.write_array(output_f,out_f);
          out_f.clear();out_f.shrink_to_fit();
          // ********************************************//
          break;
        }
      case 'w':  // window
        {
          PowerSpectrumF cPSF(params);
          cPSF.measure_window_function();
          break;  
        }
      case 'g' :  // to generate GRF
        {
          if(false==params.input_sections.InitialConditions)
            {
              throw std::runtime_error("Section InitialConditions is not enabled in parameter file"); 
              exit(0);
            }
          So.message(start_all);
          PowerSpectrumF cPSF(params);
          cPSF.get_GaussianRandomField();
          string json_file_plot ="plot_file_power_spectrum_grf.json";
          std::ofstream jfile(json_file_plot);
          json j;
          j["power_spectrum_grf"]="k, P_grf, P_th_mesh, P_fa, P_fa-mixed";
          j["power_spectrum_grf"]=true;
          j["kmin"] = 0.005;
          j["kmax"] = 2.0;
          j["sample"] = "GRF";
          j["redshift"] = params._redshift();
          j["output_file"] = cPSF._file_power_grf();
          j["input_file"] = params._dir()+params._ic_power_file();
          j["clustering_space"] = "GRF";
          jfile<<j.dump(4);
          cout<<"Generating json file "<<json_file_plot<<" for plotting"<<endl;  
          jfile.close();
          system("python3 ../python/cosmolib_plots.py plot_file_power_spectrum_grf.json &");
    
          json_file_plot ="plot_file_grf.json";
          std::ofstream jfilen(json_file_plot); 
          json jn;
          jn["Lbox"]=params._Lbox();
          if (params._Generate_FA()==true)  
            jn["fa"]=true;
          if (params._Generate_FA()==false)  
            jn["grf"]=true;
          jn["Nft"] = params._Nft() ;
          jn["Nslices"] = 20;
          jn["Initial_slice"] = 20;
          jn["sample"] = params._Generate_FA()==true ? "Fixed amplitude IC ": "Gaussian random field";
          jn["kmax_fa"] = params._Kmax_FA();
          jn["redshift"] = params._redshift();
          jn["output_file_grf"] = cPSF._output_file_grf()+".dat";
          jn["output_file_fa"] = cPSF._output_file_fa()+".dat";
          jn["output_file_mixed"] = cPSF._output_file_mixed_grf_fa()+".dat";
          jfilen<<jn.dump(4);
          cout<<"Generating json file "<<json_file_plot<<" for plotting"<<endl;  
          jfilen.close();
          system("python3 ../python/cosmolib_plots.py plot_file_grf.json &");
          break;
        }
      case 'm':  // toCOmput Marked correlation function
        {
          So.message(start_all);
          string par_file = argv[2];
          Params params(par_file);
          PowerSpectrumF cPSF(params);
        //  cPSF.measure_marked_correlation_function();
          break;
        }
      default:
          break;

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


