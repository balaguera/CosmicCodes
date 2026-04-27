////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for Estimates of angular power.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
////////////////////////////////////////////////////////////////////////////
# include "../include/CosmiCalcLIB.hpp"
////////////////////////////////////////////////////////////////////////////
using namespace std;


int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  ScreenOutput So(start_all, "logfile_hgaps.log");

  if(argc==1)
  {
    So.usage(argv[0]);
    exit(1);
  }

  string  par_file = argv[2];
  Params params(par_file);
  AngularPowerSpectrum Cl(params);

/*
  if (getenv("RUNNING_IN_XTERM") == nullptr)
   {

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
        */

    // Verify that section in json file is active: 
  if(false==params.input_sections.AngularPowerSpectrum)
    {
      throw std::runtime_error("Section AngularPowerSpectrum is not enabled in parameter file"); 
      exit(0);
    }

  int option_index = 0;
  int c;
  static struct option long_options[] = {
        {"power", required_argument, 0, 'p'},
        {"cross",  required_argument, 0, 'r'},
        {0, 0, 0, 0}
    };

   while ((c = getopt_long(argc, argv, "p:r:",long_options, &option_index)) != -1) 
   {
      switch (c) {

        case 'p' :
         {
            So.message_screen("Angular Power Spectrum analysis");
            Cl.measure_angular_power();
            std::string json_file_plot ="plot_file_angular_power_spectrum_measurement.json";
            std::ofstream jfile(json_file_plot);
            json j;
            j["header"]="lmin, lmid, lmax, power, var";
            j["sky_fraction"]=Cl._sky_fraction();
            j["angular_power_spectrum_measurement"]=true;
            j["show_tomographic_maps"]=params._generate_fits_files();
            j["statistics"] = params._statistics();
            j["sample"] = params._Name_survey();
            j["lmin"] = params._Lmin();
            j["lmax"] = params._Lmax();
            j["number_of_zbins"] = params._n_z_bins_tomography();
            j["number_of_lbins"] = params._N_L_bins();
            for(int i=0;i<params._n_z_bins_tomography();++i)
             {
               string ifile="output_file_"+to_string(i);
               j[ifile] = Cl.output_files_hgaps(i);
             }          
            if(true==params._generate_fits_files())
             {
               j["mask_map"] = params._Output_directory()+"MASK.fits";
               for(int i=0;i<params._n_z_bins_tomography();++i)
                { 
                  string ifile_map="output_file_map"+to_string(i);
                  j[ifile_map] = Cl.output_files_hgaps_maps(i);
                }          
             }
            if(params._compute_mixing_matrix())
             {
               for(int i=0;i<params._N_L_bins();++i)
                {
                  j["show_mixing_matrix"] = true;
                  string ifile="output_file_mixing_matrix_lbin"+to_string(i);
                  j[ifile] = Cl.output_files_mixing_matrix(i);
                }
             }
            jfile<<j.dump(4);
            cout<<"Generating json file "<<json_file_plot<<" for plotting"<<endl;  
            jfile.close();
            system("python3 ../python/cosmolib_plots.py plot_file_angular_power_spectrum_measurement.json ");
            cout<<"============================================================================"<<endl;
            cout<<RED<<"Missing TASKS in HGAPS: check ARF, and z bins in plots"<<RESET<<endl;
            cout<<"============================================================================"<<endl;

            break;
           }
        case 'r':
        {
            Cl.get_cross_Cl();
            break;
        }

        default:
            break;  
            
        }
      }
  

  exit(0);
  return 0;
}

