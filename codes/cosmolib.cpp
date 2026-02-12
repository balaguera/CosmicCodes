////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/* Example of code obtaining different cosmologiccal observables
 * Andres Balaguera,
 * 2007-2023
*/
#include "../headers/CosmiCalcLIB.h"
#include "../headers/CosmoLib.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
  time_t start_all;
  time(&start_all);
  string par_file = argv[2];
  char temp;

  Params param (par_file);
  string logfile="logfile_cosmo.log";
  ScreenOutput So(start_all, logfile);
  CosmoLib Clib(param);


int option_index = 0;
    int c;

    static struct option long_options[] = {
        {"cosmology", required_argument, 0, 'c'},
        {"halomodel", required_argument, 0, 'm'},
        {"redshift",  required_argument, 0, 'z'},
        {"help",      no_argument,       0, 'h'},
        {0, 0, 0, 0}
    };

  if(argc==1)
  {
    So.usage(argv[0]);
    exit(1);
  }
  else{
    while ((c = getopt_long(argc, argv, "c:z:h",
                            long_options, &option_index)) != -1) {

        switch (c) {

            case 'c':
                std::cout << "Cosmology: " << optarg << std::endl;
                Clib.get_cosmological_information();
                system("python3 ../python/cosmolib_plots.py plot_file_cosmoinfo.json");
                break;
            case 'm':
                std::cout << "Halomodel: " << optarg << std::endl;
                Clib.get_hmodel();
                system("python3 ../python/cosmolib_plots.py plot_file.json");
                break;

            case 'z':
                std::cout << "Redshift: " << optarg << std::endl;
                Clib.get_dndz_gal();// partially developerd
                break;

            case 'h':
                std::cout << "Help selected\n";
                break;

            default:
                break;
        }
    }
    }


  exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



