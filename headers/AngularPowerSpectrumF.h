#ifndef _CL_FUNCTIONS_
#define _CL_FUNCTIONS_
//////////////////////////////////////////////////////////
/**
 *  @class Cl_FUNCTIONS.h
 *  @brief The class Cl_FUNCTIONS_
 *
 */
//////////////////////////////////////////////////////////
# include "Catalog.h"
# include "Galaxy.h"
using namespace std;
/////////////////////////////////////////////////////////
// Note that when we do averages over the m-modes,
// we do the following:
// (Sum_m Jlm) / (2l+1) from -l to +l is splitted
// in Jl0/(2l+1)+ sum_(m=1) Jlm/(l+0.5)
// Since we put all in the same loop, we make theinf (m==) to divide things by 2
//////////////////////////////////////////////////////////
class Cl_FUNCTIONS{

 protected:
 private:
  ParamsCl params;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  Catalog catalog;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  Catalog catalog_random;
  //////////////////////////////////////////////////////////
  /**
  * @brief Number of columns in Catalog cat
  */
  int n_columns_gal;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  int n_columns_ran;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of columns in the mask
   */
  int n_columns_mask;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type pointing, used in HealPix operations
  */
  pointing point;
  //////////////////////////////////////////////////////////
  /**
   * @brief RMS of galaxies in pixels
   */
  real_prec rms_nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of rings, depends on nside
   */
  int Nrings;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of redshift bins, set default zero
   */
  int nzbins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index identifying redhsift bin
   */
  int IZ;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the minimum redshifts in the redshift bins
   */
  vector<real_prec> z_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the maximum redshifts in the redshift bins
   */
  vector<real_prec> z_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the Wigner 3J symbols, used
   * to compute the mixing matrix Rll in it-s approximated
   * version,
   *  in the VLS
   */
  vector<vector<vector<real_prec> > > Wigner3J;
  //////////////////////////////////////////////////////////
   Gnuplot gp_power;
   //////////////////////////////////////////////////////////
   /**
    * @brief Object of type Cosmology
   */
   Cosmology Cosmo;
   //////////////////////////////////////////////////////////
   /**
    * @brief NUmber of used galaxies, i.e, not masked
   */
   real_prec mean_weight;
public:
  //////////////////////////////////////////////////////////
  /**
   * @brief Constructor
   */
  Cl_FUNCTIONS(){
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);
  }

  Cl_FUNCTIONS(ParamsCl _pars):params(_pars){
     std::cout.precision(6);
     std::cout.setf(ios::showpoint);
     std::cout.setf(ios::scientific);

#ifdef _USE_GNUPLOT_
     this->gp_power<<"set border linewidth 1.5\n";
     this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
     this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
     this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
     this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
#endif
     this->catalog.set_params(this->params);

     string file_Jlm_init_fixed=this->params._output_file_Jlm()+".txt";
     this->params.set_output_file_Jlm(file_Jlm_init_fixed);

  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Default destructor
   */
  ~Cl_FUNCTIONS(){}

  //////////////////////////////////////////////////////////
  /**
  *    @brief Inpput/Output
  */
  FILE_MANAGER<int> Fmi;
  //////////////////////////////////////////////////////////
  /**
   * @brief Fmd
  */
  FILE_MANAGER<real_prec> Fmd;
  //////////////////////////////////////////////////////////
  /**
   * @brief Fmd
  */
  FILE_MANAGER<float> Fmf;
  //////////////////////////////////////////////////////////
  /**
   * @brief Random generator
  */
  gsl_rng * r;
  //////////////////////////////////////////////////////////
  /**
  * @brief  Type of code running. Do not touch it.
  */
  string code;

  /**
   * @brief
  */
  int ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec zmin_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec zmax_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec zmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec zmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec MKmax_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec MKmin_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  real_prec area_pixel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Expected number of galaxies in one pixel
  */
  real_prec mean_number_galaxies_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  real_prec mean_number_galaxies; //expected number of galaxies. Ngal/Area
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  real_prec rms_ngal;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  real_prec shot_noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> theta_new_pix;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Mean_ngal_pix; // mean surface density of galaxies per pixel
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Mean_ngal; //mean surface density of galaxies in sample
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Shot_Noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<vector<real_prec> > Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector to the Alms in Magnitude or Redshift bins
   */
  vector<vector<vector<complex<real_prec> > > > Blm;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the Window function
  */
  vector<real_prec> Wl;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the Window function
  */
  vector<real_prec> VShot_Noise;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containing the l-modes
  */
  vector<int> lvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Clvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> lbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> lbin_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> lbin_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Clbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> Clbin_meas;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  vector<real_prec> eClbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of angular modes
  */
  vector<int> nmodes;
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power of the mask in lbins
  */
  vector<real_prec> Wlbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Mixing matrix Rll
  */
  vector<vector<real_prec> > R;
  //////////////////////////////////////////////////////////
  /**
   * @brief MIxing matrix is bins, l-lbin
  */
  vector<vector<real_prec> > Rll_bin;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used galaxies, i.e, not masked
  */
  ULONG ngal_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used random objects, i.e, not masked
  */
  ULONG nran_used;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of used random objects, i.e, not masked
  */
  ULONG nran;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of observed pixels from the mask
  */
  ULONG n_observed_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief NUmber of pixels in the mask
  */
  ULONG n_total_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  ULONG n_pixels;

  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void set_mean_ngal_pix();

  //////////////////////////////////////////////////////////
  /**
   * @brief   Index used to pindown the redshift label
  */
  void set_IZ(int newiz){this->IZ=newiz;}
  //////////////////////////////////////////////////////////
  /**
   * @brief   Define the type of L-bins
  */
  void set_Lbins();
  //////////////////////////////////////////////////////////
  /**
   * @brief  COMpute and set Healpix related quantities from the mask
  */
  void set_healpix_pars();

  //////////////////////////////////////////////////////////
  /**
   * @brief  Array for the ALm in redshift bins
  */
  void set_BLMZ();
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the mixing matrix
  */
  void get_mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief Write the parameters in screen
  */
  void write_pars_cl_estimator();
  //////////////////////////////////////////////////////////
  /**
   * @brief COmputes the Wigner symbols
  */
  void W3J();
   //////////////////////////////////////////////////////////
  /**
   * @brief  Get the ALms from the map, our way
  */
  void Map2Alm(Healpix_Map<real_prec>, Alm<xcomplex <real_prec> >&,  Alm<xcomplex <real_prec> >&, arr<arr<real_prec> >& );
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the Alm from the map
  */
  void Map2Alm(Healpix_Map<real_prec>, Alm<xcomplex <real_prec> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get the map from the Alms , using Healpix
  */
  void Alm2Map(Alm<xcomplex <real_prec> >&, Healpix_Map<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get Alm from the random map
  */
  void Map2Alm_ran(Healpix_Map<real_prec>,Healpix_Map<real_prec>, Alm<xcomplex <real_prec> >&,  Alm<xcomplex <real_prec> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get the function Jlm used n the D-like estimator
  */
  void Map2Jlm(Healpix_Map<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Map2Ilm_Jlm(Alm<xcomplex <real_prec> >&);
  //////////////////////////////////////////////////////////
  /**
   * @brief   .
  */
  void Alm2Cl(Alm<xcomplex <real_prec> >&, vector<int>&, vector<real_prec>&);
  /**
   * @brief   .
  */
  //////////////////////////////////////////////////////////
  /**
   * @brief   COnvert Alm to Cl when partial sky coverage is prensent
  */
  void Alm2Cl(string est, Alm<xcomplex <real_prec> >&, Alm<xcomplex <real_prec> >&, arr<arr<real_prec> >&, vector<int>&, vector<real_prec>&, vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  COnvert Alm to Cl for full sky
  */
  void Alm2Cl(int, int, Alm<xcomplex <real_prec> >&,  vector<real_prec> &);  // the one is
  void Alm2Cl_sn(int, Alm<xcomplex <real_prec> >&, Alm<xcomplex <real_prec> >&,  vector<real_prec> &);  // the one is this for sn
  //////////////////////////////////////////////////////////
  /**
   * @brief   Get Gaussian errors on the estimateos of Cl
  */
  void get_eCl(int, int, vector<real_prec>,vector<real_prec>,vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Get mixing matrix
  */
  void Mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @brief ANgular power in bins
  */
  void get_Cl_bins(vector<real_prec>, vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Gaussian errors in bins
  */
  void get_eCl_bins(vector<real_prec>, vector<real_prec>&);
   //////////////////////////////////////////////////////////
  /**
   * @brief Get the Mixing matix in bins
  */
  void Rll_bins();
  //////////////////////////////////////////////////////////
  /**
   * @brief Main member function doing all the job
  */
  void get_angular_power();

  //////////////////////////////////////////////////////////
  /**
   * @brief Main member function doing all the job
  */
  void get_Cl_bias();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_alm_gaussian(int, vector<real_prec>&, Alm<xcomplex<real_prec> >&);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_pixel_window();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec>pixel_window;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_min_clth;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec z_max_clth;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void read_mixing_matrix_lbins();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void read_mixing_matrix_l();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_lg_map(Healpix_Map<real_prec>, Healpix_Map<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  vector<real_prec> nbar_photo;
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void get_cross_Cl();
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */

  void set_mean_weight(real_prec new_ww){this->mean_weight=new_ww;}
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  real_prec _mean_weight(){return this->mean_weight;}
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
  *  @brief
  */
  void Get_shot_noise(Alm<complex<real_prec>>&Ilm);

};



#endif


