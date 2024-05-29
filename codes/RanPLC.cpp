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
# include "../Headers/Galaxy.h"

using namespace std;

// **********************************************************************************************************
//***********************************************************************************************************
struct s_bins_info{
    int index;
    real_prec min;
    real_prec max;
    int Nbins;
    real_prec delta;
    string type;
    void get_delta()
    {
     this->delta=(this->max-this->min)/static_cast<real_prec>(this->Nbins);
    };
    void show(){
      cout<<"Index = "<<this->index<<endl;
      cout<<"Min = "<<this->min<<endl;
      cout<<"Max = "<<this->max<<endl;
      cout<<"Nbins = "<<this->Nbins<<endl;
      cout<<"Delta = "<<this->delta<<endl;
      cout<<"Type of binning = "<<this->type<<endl;
    }
};


// **********************************************************************************************************
//***********************************************************************************************************
class GALAXY{








int main(int argc, char *argv[]){
  string par_file = argv[1];

  GALAXY Gal;

  Gal.read_pars(par_file);
  Gal.set_file_names();
  Gal.output_dir="/home/andres/data/Numerics/Pinocchio/plc/";

  int nn;
  auto iz_vls=0;

  if(Gal.redshift_type=="p")
      nn=100;
  if(Gal.redshift_type=="s")
      nn=1000;


 Gal.total_area=4.*M_PI/8.0;
 Gal.set_pars(iz_vls);

 Gal.get_cosmo(nn,true);


 Gal.set_vec();

 Gal.read_input_cats("g");
 Gal.gp.i_ra=Gal.i_ra;
 Gal.gp.i_dec=Gal.i_dec;
 Gal.gp.i_zs=Gal.i_zs;
 Gal.i_z=Gal.i_zs;

 Gal.phi_max=360.;//Gal.get_maxa(Gal.i_ra);
 Gal.phi_min=0.;//Gal.get_mina(Gal.i_ra);
 Gal.dec_max=90;//-Gal.get_maxa(Gal.i_dec);
 Gal.dec_min=45;//90-Gal.get_mina(Gal.i_dec);

 Gal.get_dndz();// these two are enough to tabulate in the nw catalog the nbar
 Gal.get_nbar();
 Gal.bspline(1,20);



 s_bins_info bx;
 bx.index=Gal.i_zs;
 bx.min=Gal.z_min;
 bx.max=Gal.z_max;
 bx.Nbins=Gal.N_bin_z_low_res;
 bx.type="linear";
 bx.get_delta();
 bx.show();

 s_bins_info by;
 by.index=Gal.i_mass;
 by.min=Gal.lmass_min;
 by.max=Gal.lmass_max;
 by.Nbins=Gal.N_bin_lmass;
 by.type="log";
 by.get_delta();
 by.show();

 Gal.get_P_X_Y(&bx,&by, "data");

 Gal.get_new_cat();



 Gal.prop.clear();
 Gal.prop.shrink_to_fit();

 Gal.get_random_cat(Gal.output_dir+"Random.txt", 5*Gal.NGAL);
 Gal.get_P_X_Y(&bx,&by, "random");
 Gal.prop.clear();
 Gal.prop.shrink_to_fit();


 return 0;
}



