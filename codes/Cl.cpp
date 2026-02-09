// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************

/* Example of code obtaining different cosmologixcal observables
 * Andres Balaguera,
 * 2007-2019
*/
// **************************************************************************
// **************************************************************************
// **************************************************************************
// **************************************************************************
# include <iostream>
# include <fstream>
# include <vector>
# include <algorithm>
# include <math.h>
# include <sstream>
# include <iomanip>
# include <string>
# include "../Headers/ClFunctions.h"
# include "../Headers/AngularPowerSpectrum.h"

using namespace std;


int main(int argc, char *argv[]){
  string par_file = argv[1];

  string par_file1 = argv[2];
  Cl_FUNCTIONS cCl(par_file1);
  cCl.get_my_power();

  Cl_FUNCTIONS cl_cb;
  AngularPowerSpectrum cl_th;
  return 0;
}



