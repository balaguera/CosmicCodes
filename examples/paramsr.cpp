////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for BMT.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
# include "../include/CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
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

  Params params ("default_all_on.json");

  if (getenv("RUNNING_IN_XTERM") == nullptr) {
      std::string cmd = "RUNNING_IN_XTERM=1 xterm  -bg black -fg white  -hold -e ";
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
 

  int option_index = 0;
  int c;
  static struct option long_options[] = {
        {"par", required_argument, 0, 'p'},
        {"show", required_argument, 0, 's'},
        {0, 0, 0, 0}
  };


  while ((c = getopt_long(argc, argv, "p:s:",long_options, &option_index)) != -1)
    {

      switch (c) {

        case 'p': // displays input parameters
          {
            std::cout<<"=========================="<<endl;
            std::cout<<"Parameters of CosmicCodes"<<endl;
        news:
            std::cout<<"=========================="<<endl;
            std::cout<<"Enter parameter name"<<endl;
            std:string par_name;
            cin>>par_name;

            params.explain_pars(par_name);
            goto news;
            break;
          }
        case 's' : // displays input parameters
          {
            string par_file = argv[2];
            Params params(par_file);
            params.show_params();
            break;
          }
      }
    }

  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
