////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for BMT.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
# include "CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
//  system("clear");
  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_BiasMT;

  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);
  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "n:s:")) != -1)
    {
      if(temp=='h')
        So.usage(argv[0]);

      else if(temp=='n') // displays input parameters
        {
          string par_name = argv[2];
          Params params("default.json");
          params.explain_pars(par_name);
        }
      else if(temp=='s') // displays input parameters
        {
          string par_file = argv[2];
          Params params(par_file);
          params.show_params();
        }
    }
  exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
