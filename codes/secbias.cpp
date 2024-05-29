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
  while((temp =  getopt(argc, argv, "hadi:s:")) != -1)
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
          So.show_preproc();
      else if('s'==temp)  // SECONDARY BIAS ANALYSIS
        {
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          params.set_mass_assignment_scheme("CIC");
          PowerSpectrumF power(params);
          power.halo_bias_analysis("redshit_space");// this argument is overriden by parameter use_real_and_redshift_space, read from parameter file
        }
      else if('v'==temp)  // to draw warnings
        {
          So.show_warnings();
       }
      else if('?'==temp)
       {
         cout<<endl;
         cout<<"Argument not recognized."<<endl;
         cout<<"Please run ./cosmicatlass.exe -h."<<endl;
         exit(1);
       }
    }

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


