//////////////////////////////////////////////////////////
/** 
 *  @class AngularPowerSpectrum.h
 *  @brief Header file for the class AngularPowerSpectrum
 *  @file AngularPowerSpectrumF.h
 *  @author Andres Balaguera-Antolínez
 */
//////////////////////////////////////////////////////////


#ifndef _AngularPowerSpectrum_
#define _AngularPowerSpectrum_

# include <boost/tuple/tuple.hpp>
# include "Params.hpp"
# include "FileOutput.hpp"
# include "HaloTools.hpp"
# include "CoordinateSystem.hpp"
# include "FftwFunctions.hpp"
# include "Cwclass.hpp"
# include "ScreenOutput.hpp"
# include "McmcFunctions.hpp"
# include "GnuplotC.hpp"
using namespace std;

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////


class AngularPowerSpectrum{
  
 protected:

 private:
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type params
   */
   Params params;
//////////////////////////////////////////////////////////
  /**
  * @private
  * @brief Object of type FileOutput
  */
  FileOutput File;
//////////////////////////////////////////////////////////
  /** 
  * @private
  * @brief Object of type ScreenOutput 
  */
  ScreenOutput So;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of columns in the mask
   */
  int n_columns_mask;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type structure pointing, used in HealPix to determine the Healpix index from the angular coordinates and viceverza
  */
  pointing point;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of rings for sky healpixelization, depends on nside.
   */
  int Nrings;
//////////////////////////////////////////////////////////
   /**
   * @private
   * @brief Index identifying redshift bin within a tomographic analysis.
   */
  int index_tomographic_zbin;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector containing the minimum redshifts in the redshift bins defined for tomographic analysis.
   */
  vector<healpix_real> z_min;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector containing the maximum redshifts in the redshift bins defined for tomographic analysis.
   */
  vector<healpix_real> z_max;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Container for the Wigner 3J symbols, used to compute the mixing matrix Rll in it-s approximated version,
   */
  vector<vector<vector<healpix_real> > > Wigner3J;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Gnuplot
   */
  Gnuplot gp_power;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Cosmology
   */
   Cosmology Cosmo;
//////////////////////////////////////////////////////////
   /**
   * @private
   * @brief Mean number of used galaxies, weighted or not weighted.
   */
   healpix_real mean_weight;

   healpix_real sky_fraction;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////




   public:
  
  /**
   * @brief Constructor
   */
  AngularPowerSpectrum(){}

  AngularPowerSpectrum(Params &_pars):params(_pars){
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

     string file_Jlm_init_fixed=this->params._output_file_jlm()+".txt";
     this->params.set_output_file_jlm(file_Jlm_init_fixed);

  }
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Default destructor
   */
  ~AngularPowerSpectrum(){}

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  POisson shot noise correction applied to the angular power spectrum .
  */
  healpix_real shot_noise;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Container for the mean surface density of galaxies per pixel, for each redshift bin.
  */
  vector<healpix_real> Mean_ngal_pix; 
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Container for the mean number of galaxies, for each redshift bin. .
  */
  vector<healpix_real> Mean_ngal; 
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Container for the Poisson shot noise contribution for each redshift bin..
  */
  vector<healpix_real> Shot_Noise;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Container for the matrix elements Jlm
  */
  vector<vector<healpix_real> > Jlm;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Complex container for the harmonic decompoisiotn in different redshift bins.
   */
  vector<vector<vector<complex<healpix_real> > > > Blm;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Vector containing the Window function
  */
  vector<healpix_real> Wl;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Vector containing the l-modes
  */
  vector<healpix_real> lvec;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Vector containing the angular power spectrum
  */
  vector<healpix_real> Clvec;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   .Vector containing bin in multipoles
  */
  vector<healpix_real> lbin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   .Vector containing minimum valoe of multipole in each l-bin.
  */
  vector<healpix_real> lbin_min;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   .Vector containing maximim valoe of multipole in each l-bin.
  */
  vector<healpix_real> lbin_max;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Container for the power spectrum in bins of multipols
  */
  vector<healpix_real> Clbin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   .Container for the variance of the power spectrum in bins of multipols
  */
  vector<healpix_real> eClbin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the number of angular modes in l-bins
  */
  vector<int> nmodes;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Angular power of the mask in lbins
  */
  vector<healpix_real> Wlbin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Mixing matrix Rll
  */
  vector<vector<healpix_real> > R;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Mixing matrix in l-lbin
  */
  vector<vector<healpix_real> > Rll_bin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Number of pixels according to the NSIDE input parameter, or read from the survey mask..
  */
  ULONG n_pixels;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Method to resize vectors of this class.
  */
  void set_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   . Method to resize vectors of this class, associated to mean galaxy density in pixels and or redshift bins.
  */
  void set_vectors_mean_ngal_pix();

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Index used to pindown the redshift label
  */
  void set_index_tomographic_zbin(int newiz){this->index_tomographic_zbin=newiz;}
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Define the type of L-bins
  */
  void set_Lbins();
  //////////////////////////////////////////////////////////
  /***
   * @public
   * @brief  COMpute and set Healpix related quantities from the mask
  */
  void set_healpix_pars();

  //////////////////////////////////////////////////////////
  /***
   * @public
   * @brief  Array for the ALm in redshift bins
  */
  void set_BLMZ();
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Get the mixing matrix
  */
  void get_mixing_matrix(HaloTools &ht);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Write the parameters in screen
  */
  void write_pars_cl_estimator(HaloTools &ht);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief COmputes the Wigner symbols
  */
  void W3J();
   //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Get the ALms from the map, our way
  */
  void Map2Alm(HaloTools &ht, Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&,  Alm<xcomplex <healpix_real> >&, arr<arr<healpix_real> >& );
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Get the Alm from the map
  */
  void Map2Alm(HaloTools &ht, Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Get the map from the Alms , using Healpix
   * @note // Note that when we do averages over the m-modes,
   * we do the following:
   * (Sum_m Jlm) / (2l+1) from -l to +l is splitted
   * in Jl0/(2l+1)+ sum_(m=1) Jlm/(l+0.5)
   * Since we put all in the same loop, we make theinf (m==) to divide things by 2
  */
  void Alm2Map(Alm<xcomplex <healpix_real> >&, Healpix_Map<healpix_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @param ht: HaloTools 
   * @brief  Get Alm from the random map
  */
  void Map2Alm_ran(Healpix_Map<healpix_real>,Healpix_Map<healpix_real>, Alm<xcomplex <healpix_real> >&,  Alm<xcomplex <healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Method to compute the elements Jlm used n the D-like estimator
   * @param ht: HaloTools 
   * @param  Map: Input Healpix map
  */
  void Map2Jlm(HaloTools &ht, Healpix_Map<healpix_real>&Map);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Method to compute the elements Jlm and Ilm 
   * @param ht: HaloTools 
   * @param Alm: complex harmonic decomposiiton Alm
  */
  void Map2Ilm_Jlm(HaloTools &ht, Alm<xcomplex <healpix_real> >&Alm);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Method to convert Alm to Cl taking into account different redshift bins.
   * @param  a: Index of the first redshift bin
   * @param  b: Index of the second redshift bin
   * @param  Ilm: Complex harmonic decompositon Ilm of the mask
   * @param  Cl: Vector with estimates of angular power spectrum  
  */
  void Alm2Cl(HaloTools &ht, int a, int b, Alm<xcomplex <healpix_real> >&Ilm,  vector<healpix_real> &Cl);  
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  Method to transform Alm to Cl usef for the estimation of shot noise using randomized intances of the coordinates.
   * @param  a: Index of the first redshift bin
   * @param  Alm_r: Complex harmonic decompositon of the randomized field.
   * @param  Ilm: Complex harmonic decompositon Ilm of the mask
   * @param  Cl: Vector with estimates of angular power spectrum  
  */
  void Alm2Cl_sn(HaloTools &ht, int a, Alm<xcomplex <healpix_real> >&Alm_r, Alm<xcomplex <healpix_real> >&Ilm,  vector<healpix_real> &);  // the one is this for sn
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief   Get Gaussian errors on the estimates of Cl when measured in tomographic bins
   * @param  a: Index of the first redshift bin
   * @param  b: Index of the second redshift bin
   * @param  Cli: Raw measurements in first bin 
   * @param  Clj: Raw measurements in second bin 
   * @return  eCl Vector with estimates of uncertainty of angular power spectrum  
  */
  void get_eCl(HaloTools & ht, int a , int b, vector<healpix_real>&Cli,vector<healpix_real>&Clj,vector<healpix_real>&eCl);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief  This method computes the mixing matrix, usefl for likelihood analysis.
  */
  void Mixing_matrix();
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief CMethod to convert raw spectrum to bineed power spectrum
   * @param  Cl: Vector with estimates of angular power spectrum  
   * @return Clb: Vector with binned estimates of angular power spectrum  
  */
  void get_Cl_bins(vector<healpix_real>&Cl, vector<healpix_real>&Clb);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Gaussian errors in bins
  */
  void get_eCl_bins(vector<healpix_real>&, vector<healpix_real>&);
   //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Get the Mixing matix in bins
  */
  void Rll_bins();
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief This is the main method doing the job.
  */
  void measure_angular_power();
  //////////////////////////////////////////////////////////
  /**
   * @public
  *  @brief Generate gaussian random coefficients in harmonic space. Used for the production of log-normal catalogs.
  */
  void get_alm_gaussian(int, vector<healpix_real>&, Alm<xcomplex<healpix_real> >&);
  //////////////////////////////////////////////////////////
  /**
   * @public
  *  @brief Method to read external mixing matrix written in bins of multipoles, as generated by this same code.
  */
  void read_mixing_matrix_lbins();
  //////////////////////////////////////////////////////////
  /**
   * @public
  *  @brief Method to read external mixing matrix written in multipoles, as generated by this same code.
  */
  void read_mixing_matrix_l();
  //////////////////////////////////////////////////////////
  /**
  * @public
  * @brief Method to produce log normal maps.
  * @param map: Input map 
  * @return lmap: LOg normal map
  */
  void get_lg_map(Healpix_Map<healpix_real>&map, Healpix_Map<healpix_real>&lmap);
  //////////////////////////////////////////////////////////
  /**
  * @public
  * @brief Method to mesure the cross angular power spectrum between two samples with different masks
  */
  void get_cross_Cl();
  //////////////////////////////////////////////////////////
  /**
  * @public
  *  @brief
  */
  void set_mean_weight(healpix_real new_ww){this->mean_weight=new_ww;}
  //////////////////////////////////////////////////////////
  /**
  * @public
  *  @brief
  */
  healpix_real _mean_weight(){return this->mean_weight;}
  //////////////////////////////////////////////////////////
  /**
  * @public
  *  @brief Method to retrieve string parameters
  */
  string output_files_hgaps(int i){return this->params.output_files_hgaps[i];}
  //////////////////////////////////////////////////////////
  /**
  * @public
  * @brief Method to retrieve string parameters
  */
  string output_files_hgaps_maps(int i){return this->params.output_files_hgaps_maps[i];}
  //////////////////////////////////////////////////////////
  /**
  * @public
  * @brief Method to retrieve string parameters
  */
  string output_files_mixing_matrix(int i){return this->params.output_files_mixing_matrix[i];}
  //////////////////////////////////////////////////////////
  /**
  * @public
  * @brief Method to retrieve the i-th element of the container for the raw power spectrum.
  */
  real_prec _Cl(int i){return this->Clvec[i];}

  real_prec _sky_fraction(){return this->sky_fraction;}
};



#endif


