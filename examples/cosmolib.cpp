////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/* Example of code obtaining different cosmologiccal observables
 * Andres Balaguera,
 * 2007-2023
*/
#include "../include/CosmiCalcLIB.hpp"
#include "../include/CosmoLib.hpp"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  ScreenOutput So(start_all,"logfile_cosmo.log");
  
  if(argc==1)
  {
    So.usage(argv[0]);
    exit(1);
  }

  string par_file = argv[2];
  Params param (par_file);
  CosmoLib Clib(param);

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
 
    // Verify that section in json file is active: 
  if(false==param.input_sections.CosmologicalLibrary)
     {
        throw std::runtime_error("Section CosmologicalLibrary is not enabled in parameter file"); 
        exit(0);
   }

  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"cosmology", required_argument, 0, 'c'},
        {"halomodel", required_argument, 0, 'm'},
        {"redshift",  required_argument, 0, 'z'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

   while ((c = getopt_long(argc, argv, "c:m:z:h",long_options, &option_index)) != -1)
    {
        switch (c) {

            case 'c':
            {
                So.message_screen("Cosmology");
                Clib.get_cosmological_information();
                system("python3 ../python/cosmolib_plots.py plot_file_cosmoinfo.json &");
                break;
            }
            case 'm':
            {
                So.message_screen("Halo Model");
                Clib.get_hmodel();
                system("python3 ../python/cosmolib_plots.py plot_file.json ");
                break;
            }
            case 'z':
            {
                std::cout << "Redshift: " << optarg << std::endl;
                Clib.get_dndz_gal();// partially developerd
                break;
            }
            case 'h':
            { 
                std::cout << "Help selected\n";
                break;
            }
            default:
                break;
        }
    }
    
   exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



