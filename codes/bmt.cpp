////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for BMT.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
# include "CosmiCalcLIB.h"
# include "BiasMappingTechnique.h"
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

  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);
  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  while((temp =  getopt(argc, argv, "hadi:b:")) != -1)
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
          par_file_BiasMT = argv[2];
          Params Par(par_file_BiasMT);
          cout<<"here "<<par_file_BiasMT<<endl;
          So.set_params(Par);
          BiasMT BiasMT(Par,true);
          BiasMT.set_So(So);
          BiasMT.execute();
       }
    }
  exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
