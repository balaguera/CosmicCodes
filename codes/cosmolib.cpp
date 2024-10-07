////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/* Example of code obtaining different cosmologiccal observables
 * Andres Balaguera,
 * 2007-2023
*/
#include "CosmiCalcLIB.h"
#include "CosmoLib.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
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
  if(argc==1)
  {
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "c:z:h:")) != -1)
    {
      if(temp=='c') // Run BiasMT
        {
          Clib.get_cosmolib();
       }
      else if(temp=='h') // Run BiasMT
        {
          Clib.get_dndz_gal();// partially developerd
       }
      else if(temp=='z') // Get cosmology functions as as funtion f redshift
        {
          Clib.get_dndz_gal();// partially developerd
       }
    }
  exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



