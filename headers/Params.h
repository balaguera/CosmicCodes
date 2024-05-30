////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @class<Params>
 *  @file Params.h
 *  @brief Headers for class Params
 *  @author Andrés Balaguera-Antolínez,
 *  @date 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef __PARAMS__
#define __PARAMS__
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <typeinfo>
using namespace std;
#include "NumericalMethods.h"  // def.h is included in the NumericalMethods.h file
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class Params
{

private :
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  vector<pair<string, real_prec>> parameter_number;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  vector<pair<string, string>> parameter_string;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  vector<pair<string, bool>> parameter_boolean;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  vector<pair<string,vector<real_prec>>> parameter_vectors;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  vector<pair<string,vector<string>>> parameter_vector_string;

  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string statistics;
  //////////////////////////////////////////////////////////
  /**
   *  @name 
   */
  real_prec Initial_Redshift_DELTA;
  //////////////////////////////////////////////////////////
  /**
   *  @name
   */
  real_prec Initial_Redshift_SIM;
  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string Input_dir_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the a new reference catalog aimed to be assigned properties is located
   **/
  string Input_dir_cat_new_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief input directory where the tracers, dm and random catalogues are stored
   **/
  string Input_dir_cat_TWO;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string ic_WN_dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string ic_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  unsigned long ngal_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string delta_grid_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  string delta_grid_file2;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  string delta_grid_file3;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string delta_grid_file4;
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
   **/
  int measure_cross_from_1;
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
   **/

  int measure_cross_from_2;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string Name_redshift_mask;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string Name_binary_mask;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  string type_of_object;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  int IC_index;
  //////////////////////////////////////////////////////////
  /**
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
   **/
  string file_random;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  bool use_ic_file;
  //////////////////////////////////////////////////////////
  /**
   *  @note coordinate system of the object catalogue
   *
   *  0 &rarr; positions are given in Cartesian coordinates (X, Y, Z)
   *
   *  1 &rarr; positions are given in equatorial coordinates (RA, Dec, r)
   *
   *  2 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *
   *  3 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *
   *  @note with option 2, the code uses the set of cosmological
   *  parameters to transform redshift z to comoving distance; with
   *  option 3, the redshift is directly used as radial coordinate
   */
  int sys_of_coord_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Column with the infomtion of the number of sub_structures of each tracer (if parent halo)
   **/
  int i_sf_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_mass_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_mass_r;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  int sys_of_coord_dm;
  //////////////////////////////////////////////////////////

  /// the column where the first object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
   **/
  int i_coord1_dm;

  //////////////////////////////////////////////////////////
  /// the column where the second object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
   **/
  int i_coord2_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  int Number_of_chunks_new_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  bool weight_vel_with_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minium apparent magnitude
   **/
  real_prec mKmin;

  //////////////////////////////////////////////////////////
  /**
   * @brief Minium apparent magnitude
   **/
  real_prec mKmax;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec MASS_units;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec LOGMASSmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec LOGMASSmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/

  real_prec VMAXmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec VMAXmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/

  real_prec VRMSmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec VRMSmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec RSmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec RSmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec CONCENTRATIONmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec CONCENTRATIONmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec SPINmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec SPINmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec VIRIALmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec VIRIALmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec BTOAmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec BTOAmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec CTOAmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec CTOAmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool use_real_and_redshift_space;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool Get_marked_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool Get_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool Get_cross_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool Get_pearson_coefficient;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool Get_spearman_coefficient;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  int Number_of_MultiLevels;
  //////////////////////////////////////////////////////////
  /// the column where the third object coordinate (according to the system of coordinates of the catalog) is written
  /**
   * @brief
   **/
  int i_coord3_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_v1_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_v2_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_v3_dm;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int i_weight1_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the second object weight is written

   **/
  int i_weight2_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the third object weight is written
   **/
  int i_weight3_g;

  //////////////////////////////////////////////////////////
  /**
   * @brief The column where the third object weight is written
   **/
  int i_weight4_g;

  //////////////////////////////////////////////////////////
  /**
   *  @brief units of angles in the object catalogue: D &rarr; degrees; R &rarr; radians
   */
  string angles_units_g;

  //////////////////////////////////////////////////////////
  /**
   *  @brief units of angles in the object catalogue: D &rarr; degrees; R &rarr; radians
   */
  string angles_units_r;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  real_prec M_exclusion;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Deltal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Deltamu;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int Number_of_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  int Number_of_new_mocks;
  //////////////////////////////////////////////////////////
  /**
   *  @brief coordinate system of the object catalogue
   *  0 &rarr; positions are given in Cartesian coordinates (X, Y, Z)
   *  1 &rarr; positions are given in equatorial coordinates (RA, Dec, r)
   *  2 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *  3 &rarr; positions are given in pseudo-equatorial coordinates (RA, Dec, z)
   *  @note with option 2, the code uses the set of cosmological
   *  parameters to transform redshift z to comoving distance; with
   *  option 3, the redshift is directly used as radial coordinate
   */
  int sys_of_coord_r;

  //////////////////////////////////////////////////////////
  /// the column where the first random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord2_r;

  //////////////////////////////////////////////////////////
  /// the column where the third random coordinate (according to the system of coordinates of the catalog) is written
  int i_coord3_r;

  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the color of tracer is located
   */
  int i_color_r ;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_stellar_mass_r;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_app_mag_r ;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_abs_mag_r ;
  //////////////////////////////////////////////////////////
  /// the column where the first random weight is written
  int i_weight1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random weight is written
  int i_weight2_r;
  //////////////////////////////////////////////////////////
  /// the column where the third random weight is written
  int i_weight3_r;

  //////////////////////////////////////////////////////////
  /// the column where the fourth random weight is written
  int i_weight4_r;

  //////////////////////////////////////////////////////////
  /// the column where the random mean number density is written
  int i_mean_density_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the first weight; false &rarr; do not use the first weight
  bool use_weight1_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the second weight; false &rarr; do not use the second weight
  bool use_weight2_r;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool use_weight3_r;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool use_weight4_r;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool get_distribution_min_separations;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec kmin_tracer_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec kmax_tracer_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec kmax_tracer_qbias;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec kmin_tracer_qbias;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool Get_Luminosity_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool Get_Color_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  int Nbins_color;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool Get_Mstellar_function;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec Mstellar_min;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec Mstellar_max;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec Color_min;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec Color_max;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  int Nbins_Mstellar;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  std::string LF_estimator;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool Get_Color_Mag_plane;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool Get_Random_Catalog;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /// Number of grid cells /per dimension for the Discrete Fourier Trasnform


  ULONG Nft_low;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  ULONG Nft_HR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  ULONG Nft_JK;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Lbox_low;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec vkernel_exponent;
  //////////////////////////////////////////////////////////
  /// true &rarr; compute a new Lbox; false &rarr; the code uses Lbox
  bool new_Lbox;

  //////////////////////////////////////////////////////////
  /// the power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false)
  bool use_random_catalog;


  //////////////////////////////////////////////////////////
  /// the angular power spectrum is computed with the help of a random catalogue (if use_random_catalog = true), or from a simulation with known mean (if use_random_catalog = false)
  bool use_random_catalog_cl;



  //////////////////////////////////////////////////////////
  ///Type of binning in Fourier space: linear/log
  string type_of_binning;

  //////////////////////////////////////////////////////////
  /// Number of log-spaced bins in Fourier space
  int N_log_bins;

  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bins in Foureir space
   * @details The size of k-bins is ndel times the fundamental mode
   */
  int ndel_data;


  //////////////////////////////////////////////////////////
  /// Ratio between the shell-width and the fundamental mode for window
  int ndel_window;

  //////////////////////////////////////////////////////////
  /// Number of mu-bins for P(k,mu)
  int N_mu_bins;

  //////////////////////////////////////////////////////////
  ///Use FKP weights (yes/no)
  bool FKP_weight;

  //////////////////////////////////////////////////////////
  /// Estimated power for FKP weights
  real_prec Pest;


  //////////////////////////////////////////////////////////
  /// Compute FKP error bars? (yes/no)
  bool FKP_error_bars;

  //////////////////////////////////////////////////////////
  /// Compute error bars following FKP exact formula(yes/no)
  /// If this is no, and the previous is yes
  /// the code uses the Veff approximation for the variance.
  bool FKP_error_bars_exact;

  //////////////////////////////////////////////////////////
  /// Is the mean number density tabulated?
  bool nbar_tabulated;
  //////////////////////////////////////////////////////////
  /// Is the mean number density tabulated?
  bool use_file_nbar;

  //////////////////////////////////////////////////////////
  /// Has the sample a constant depth?
  bool constant_depth;

  //////////////////////////////////////////////////////////
  ///Number of redshift bins to measure dNdz
  int Nbins_redshift;

  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample. Used when dNdz ahs to be measured
  real_prec redshift_min_sample;

  //////////////////////////////////////////////////////////
  /// Maximum redshift of the sample. Used when dNdz ahs to be measured
  real_prec redshift_max_sample;

  //////////////////////////////////////////////////////////
  ///  Number of dNdz bins to measure
  int N_dndz_bins;

  //////////////////////////////////////////////////////////
  /// Number of redshift bins withoin which the measuerd dNdz will be smoothed
  int new_N_dndz_bins;

  //////////////////////////////////////////////////////////
  /// Area of the survey.
  real_prec area_survey;

  //////////////////////////////////////////////////////////
  /// Resolution Healpix for pixelization. Used when no nbar is tabulated and dNdz is to be computed from a non-constant depth sample
  int Healpix_resolution;

  //////////////////////////////////////////////////////////
  /// output file for the redshift distribution
  string file_dndz;
  //////////////////////////////////////////////////////////
  /// output file for the redshift distribution
  string file_nbar;

  //////////////////////////////////////////////////////////
  /// Redefine a new line of sight
  bool new_los;


  //////////////////////////////////////////////////////////
  /// output log file for the FKP power spectrum
  string file_power_log;

  //////////////////////////////////////////////////////////
  /// output file for the power spectrum of thw window function
  string file_window;
  //////////////////////////////////////////////////////////
  /// output file for the 2d power spectrum in cartesian coordinates
  string file_power2d;

  //////////////////////////////////////////////////////////
  /// output file for the 2d power spectrum in polar coordinates
  string file_power2d_mk;

  //////////////////////////////////////////////////////////
  /// Maximum k value for the direct sum approach to Yamamoto-*Blake
  real_prec kmax_y_ds;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool use_vel_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int n_catalogues;

  //////////////////////////////////////////////////////////
  /**
   * @brief name of the output file used to store the bispectrum
   **/
  string file_bispectrum;

  //////////////////////////////////////////////////////////
  /// These parameters is used to define the shells in k-space
  /**
   * @brief
   **/
  bool use_fundamental_mode_as_kmin_bk;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmin_bk;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_x;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_y;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_z;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_x_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_y_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec delta_z_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_x;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_y;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_z;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_x_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_y_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_z_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_0;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec deltak_0_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmin;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmin_low;

  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec kmax;

  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec DeltaK_data;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec DeltaK_data_low;
  /////////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec DeltaK_window;


  /////////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  real_prec kmax_bk;
  /**
   *  @name angular power spectrum
   */
  ///@{



  string input_file_mask;
  //////////////////////////////////////////////////////////
  /// the column where the pixel of the mask is written
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  int i_mask_pixel;

  //////////////////////////////////////////////////////////
  /// the column where the RA of the pixel of the mask is written
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  int i_mask_alpha;

  //////////////////////////////////////////////////////////
  /// the column where the DEC of the pixel of the mask is written
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  int i_mask_delta;

  //////////////////////////////////////////////////////////
  /// the column where the FLAG of the mask is written
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  int i_mask_flag;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  string mark;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  real_prec rmin_cf;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  real_prec rmax_cf;
  //////////////////////////////////////////////////////////
  /**
   * @brief   NUmber of bins for correlation function
   **/
  ULONG Nbins_cf;
  //////////////////////////////////////////////////////////
  /**
   * @brief   Maximum k-value for constructing k-bins

   **/
  string rbin_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief 

   **/
  bool dilute_dm_sample;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec fraction_dilute;
  //////////////////////////////////////////////////////////
  /**
   *  @name cosmological parameters
   */
  ///@{
  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  /**
   * @brief 

   **/
  real_prec om_matter;
  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  /**
   * @brief 

   **/
  real_prec om_cdm;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>rad</SUB>: the radiation density at z=0
  /**
   * @brief 

   **/
  real_prec om_radiation;

  ///////////////////////////////////////////////////y///////
  /// &Omega;<SUB>b</SUB>: the baryon density at z=0
  /**
   * @brief 

   **/
  real_prec om_baryons;
  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>DE</SUB>: the dark energy density at z=0
  /**
   * @brief 

   **/
  real_prec om_vac;
  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>k</SUB>: the density of curvature energy
  /**
   * @brief 

   **/
  real_prec om_k;
  //////////////////////////////////////////////////////////
  /// H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc]
  /**
   * @brief 

   **/
  real_prec Hubble;
  //////////////////////////////////////////////////////////
  /// \e h: the Hubble parameter, H<SUB>0</SUB>/100
  /**
   * @brief 

   **/
  real_prec hubble;
  //////////////////////////////////////////////////////////
  /// n<SUB>spec</SUB>: the primordial spectral index
  /**
   * @brief 

   **/
  real_prec spectral_index;
  //////////////////////////////////////////////////////////
  /// w<SUB>0</SUB>: the parameter of the dark energy equation of state
  /**
   * @brief 

   **/
  real_prec wde_eos;
  //////////////////////////////////////////////////////////
  /// N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
  /**
   * @brief 

   **/
  real_prec N_eff;
  //////////////////////////////////////////////////////////
  /// &sigma;<SUB>8</SUB>: the power spectrum normalization
  /**
   * @brief 

   **/
  real_prec sigma8;
  //////////////////////////////////////////////////////////
  /// T<SUB>CMB</SUB>: the present day CMB temperature [K]
  /**
   * @brief 

   **/
  real_prec Tcmb;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   **/
  real_prec A_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec specral_index;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec alpha_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec RR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Delta_SO;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool use_wiggles;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec f_baryon;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool fixed_redshift;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_x_coord;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_y_coord;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string file_bin_z_coord;
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void warnings();
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  ///@}


  /**
   *  @name BAM parameters
   */
  ///@{




  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string par_file;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID_HR;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID_h;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG NGRID_JK;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG N_lines_binary;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.
   */
  int NY;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X to construct BIAS
   */
  ULONG NX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.
   */
  ULONG NY_MASS;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in Y to construct BIAS
   * If the Y- density field is NGP, this values is transformed
   * to the maximum number of TR in one cell.

   */
  int NY_SAT_FRAC;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int Nlambdath;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_Y_TWO;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_Y_new_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  string Name_Catalog_Y_HR;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/

  string Name_Catalog_Y_MWEIGHTED;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_REF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_new_ref;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_REF_TWO;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_BIAS_KERNEL;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_BIAS_KERNEL_TWO;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Input_Directory_X_NEW;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string XNAME;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_Catalog_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_Catalog_X_new_ref;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldx_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldy_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief INput density field for dark matter
   */
  string Name_VelFieldz_X;
  //////////////////////////////////////////////////////////

  /*
   * @brief
   */
  string Name_Catalog_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_X_NEW;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Catalog_X_NEW_TWO;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Property_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string new_Name_Property_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string YNAME;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string Name_Property_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string new_Name_Property_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  string extra_info;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_X_NEW;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iMAS_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of the TR overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_Y_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of the TR overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_Y_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of the DM overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_X_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief Minimum value of the DM overdensity
   * @details If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec delta_X_min;
  /////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of the log(1+overdensity) for TR
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_Y_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of the log(1+overdensity) for TR
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_Y_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   * @brief Maximum value of the log(1+overdensity) for DM
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_X_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief Minimum value of the log(1+overdensity) for DM
   * @detail If requiested from parameter filw, this value is overloaded
   * by computing it from the input fields
   */
  real_prec ldelta_X_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief Identification of the input density field (density, delta)
   */
  string Quantity;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins to measure property function for interpolation
   */
  ULONG NMASSbins;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Number of bins used in get_X_function and assign_oproperty_new_new
   */
  ULONG NPROPbins_bam;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of mass bins to mesure power spectrum
   */
  int NMASSbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NVMAXbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NVRMSbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NCONCENTRATIONbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NSPINbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NSPINBULLOCKbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NMACHbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NDACHbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NBIASbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NRBIASbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NQBIASbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in local clustering
   */
  int NLCbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in tidal anisotropy
   */
  int NTAbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in tidal anisotropy
   */
  int NLOCALDMbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in peak height at halo
   */
  int NPHbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NRSbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NRVIRbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NVIRIALbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NBTOAbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int NCTOAbins_power;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of mass cuts to measure power spectrum
   */
  int NMASScuts_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief Cosmological redshift
   * @brief Read from parameter file
   */
  real_prec redshift;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec smscale;
  //////////////////////////////////////////////////////////
  /*
   * @brief Identification for a realization
   * @brief Read from parameter file
   */
  int realization;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Comp_conditional_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Comp_joint_PDF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool write_files_for_histograms;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Redefine_limits;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Convert_Density_to_Delta_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Convert_Density_to_Delta_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the TWEB classification
   */
  real_prec lambdath;

  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the TWEB classification
   */
  bool get_window_matrix;
  //////////////////////////////////////////////////////////
  /**
   * @brief CPONtainer for the minimum value pf log Mass for mass bins to be used in power spectrum
   */
  vector<real_prec> MASSbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief CPONtainer for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> MASSbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VMAXbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VMAXbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VRMSbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VRMSbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> RSbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of log Rs
   */
  vector<real_prec> RSbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> CONCENTRATIONbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> CONCENTRATIONbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> RVIRbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> RVIRbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> SPINbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> SPINbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> SPINBULLOCKbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> SPINBULLOCKbins_max;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VIRIALbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> VIRIALbins_max;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> BTOAbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> BTOAbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> CTOAbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> CTOAbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> MACHbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> MACHbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> DACHbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> DACHbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of object-by-object bias
   */
  vector<real_prec> BIASbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum  value of object-by-object bias
   */
  vector<real_prec> BIASbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of object-by-object bias
   */
  vector<real_prec> RBIASbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum  value of object-by-object bias
   */
  vector<real_prec> RBIASbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of object-by-object bias
   */
  vector<real_prec> QBIASbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum  value of object-by-object bias
   */
  vector<real_prec> QBIASbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of local_overdensity LC
   */
  vector<real_prec> LCbins_min;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum value of local_overdensity LC
   */
  vector<real_prec> LCbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of tidal anisotropy
   */
  vector<real_prec> TAbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum value of tidal anisotropy
   */
  vector<real_prec> TAbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of tidal anisotropy
   */
  vector<real_prec> LOCALDMbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum value of localoverdensiuty
   */
  vector<real_prec> LOCALDMbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the minimum value of peak height at halo
   */
  vector<real_prec> PHbins_min;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum value of peak height at halo
   */
  vector<real_prec> PHbins_max;
  //////////////////////////////////////////////////////////
  /**
   * @brief CPONtainer for the minimum value pf log Mass for mass bins
   */
  vector<real_prec> MASScuts;


  //////////////////////////////////////////////////////////
  /**
   * @brief CPONtainer for the minimum value pf log Mass for mass bins
   */
  vector<string> RANDOMfiles;
  int NRANDOMfiles;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<int> list_new_dm_fields;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<int> list_Nft_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<ULONG> list_Ntracers_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<real_prec> list_Props_Threshold_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<real_prec> list_Props_Tolerance_MultiLevels;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_new_dm_fields;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container with the LOS of new DM fields used in the case in which seveal LOS are built simultaneously
   */
  vector<string> files_dm_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<int> list_bias_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_bias_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_kernel_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_tracer_references;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  vector<string> files_tracer_field_references;


  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec lambdath_v;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Write_Scatter_Plot;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Write_PDF_number_counts;
  //////////////////////////////////////////////////////////
  /*
   * @brief Specify whether the histograms are done in log or linear scale for Y
   * @brief Read from parameter file
   */
  string Scale_Y;
  //////////////////////////////////////////////////////////
  /*
   * @brief Specify whether the histograms are done in log or linear scale for X
   */
  string Scale_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in the Mk info to build BIAS.
   */
  int n_sknot_massbin;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of bins in the Veldisp info to build BIAS.
   */
  int n_vknot_massbin;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of BAM iterations
   * @brief Read from parameter file
   */
  ULONG N_iterations_Kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of BAM iterations
   * @brief Read from parameter file
   */
  ULONG Iteration_Kernel_average;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int iteration_ini;

  //////////////////////////////////////////////////////////
  /*
   * @brief Number of DM that will be created by the approximated gravity solver
   * @ and on which the Kernel will act to create a halo mock density field.
   */
  int N_dm_realizations;
  //////////////////////////////////////////////////////////
  /*
   * @brief Initial label of the DM realizations
   */
  int N_dm_initial;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of iterations for the pre-processing of DM
   * Parameter set in the parameter file.
   * Default 0
   */
  int N_iterations_dm;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Apply_Rankordering;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Apply_Rankordering_ab_initio;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  vector<ULONG> cwt_used;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  vector<ULONG> cwv_used;

  //////////////////////////////////////////////////////////
  /*
   * @brief  Container specifying iterations at which some outputs are produced
   */
  vector<ULONG> output_at_iteration;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of CWT used, computed as the size of the
   * Bam::cwt_used container
   */

  ULONG n_cwt;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number of CWT based on the shear of velocity field used
   */

  ULONG n_cwv;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // PATCHY PARS
  /**
   * @brief
   **/
  string dataFileName;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int inputmode;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int seed;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  int seed_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool runsim;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool runv;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool diffcosmorz;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string ic_power_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool ic_alias_corrected;


  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string ic_input_type;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string ic_WN_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief Inidcates approximated gravity solver to evolve IC
   **/
  int Structure_Formation_Model;
 
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  string dir;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  bool readPS;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec slength;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec slengthv;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec vslength;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec velbias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec velbias_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec velbias_random;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasepsilon;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasrhoexp;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasone;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasL;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biassign;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biassign2;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasepsilon2;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec biasrhoexp2;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec devpois;
  //////////////////////////////////////////////////////////
  /**
   * @brief0
   **/
  real_prec deltathH;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  real_prec Nmean;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cs2;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cs3;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cst;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cpsi;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec cdeltas2;
 
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec sfac;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec ep;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec xllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec yllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec zllc;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec xobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec yobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec zobs;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  bool Normalize_IC_to_initial_redshift;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec Initial_Redshift_TH_power_file;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1_HR;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2_HR;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3_HR;

  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1_JK;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2_JK;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3_JK;


  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d1_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d2_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  //////////////////////////////////////////////////////////
  /*
   * @brief Threshold value for the V-web classification
   */
  real_prec d3_low;
  //////////////////////////////////////////////////////////
  /*
   * @brief Number if cells per dimention in the mesh. 
   */
  ULONG Nft;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */
  bool SN_correction;

  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string Output_directory;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int NMASSbins_mf;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  string vel_units_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  int masskernel;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  

  int masskernel_vel;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string dir_output;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of the box in Mpc/h
   * @brief Read from parameter file
   */
  real_prec Lbox;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  real_prec DeltaKmin;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string file_power_th;
  //////////////////////////////////////////////////////////
  /*
   * @brief To identif if RSD are included in coordinates
   */  
  bool redshift_space_coords_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for x-coordinate of tracer
   */  
  int i_coord1_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for y-coordinate of tracer
   */  
  int i_coord2_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for z-coordinate of tracer
   */  
  int i_coord3_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for x-component of tracer
   */  
  int i_v1_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Column for y-vel compònent of tracer
   */  
  int i_v2_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Column for z-vel component of tracer
   */  
  int i_v3_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Column for mas of tracer
   */  
  int i_mass_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Column for maximum circular velocity of tracer
   */  
  int i_vmax_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column index for velocity dispersion of tracer
   */  
  int i_vrms_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Column for Shape radius Rs of tracer
   */  
  int i_rs_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for Virial radius of tracer
   */
  int i_rvir_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for the dimensionless spin (Peebles by default) of tracer
   */  
  int i_spin_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief Column for the dimensionless spin (Bullock by default) of tracer
   */
  int i_spin_bullock_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the halo mean densityis written
   **/
  int i_mean_density_g;
  //////////////////////////////////////////////////////////
  ///
  /*
   * @brief the column where the object T/|W| ratio  (0.5=virialization)
   */
  int i_virial_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where ration b_to_ of the semiaxi¡s ls located
   */
  int i_b_to_a_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where ration c_to_ of the semiaxis ls located
   */
  int i_c_to_a_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the color of tracer is located
   */
  int i_color_g ;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_stellar_mass_g ;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_app_mag_g ;
  //////////////////////////////////////////////////////////
  /*
   * @brief The column where the stellar mass of tracer is located
   */
  int i_abs_mag_g ;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec RA_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec RA_max;

  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec DEC_min;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec DEC_max;
  //////////////////////////////////////////////////////////
  /*
   * @brief use the fourth weight; false &rarr; do not use the fourth weight
   */
  bool use_weight1_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief use the fourth weight; false &rarr; do not use the fourth weight
   */
  bool use_weight2_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief use the fourth weight; false &rarr; do not use the fourth weight
   */
  bool use_weight3_g;

  //////////////////////////////////////////////////////////
  /*
   * @brief use the fourth weight; false &rarr; do not use the fourth weight
   */
  bool use_weight4_g;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool weight_with_mass;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nft_random_collapse;
  //////////////////////////////////////////////////////////
  /*
   * @brief Fraction to the closest dm tracer to collapse a random tracer
   */  
  real_prec Distance_fraction;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string mass_assignment_scheme;
  //////////////////////////////////////////////////////////
  /**
   * @brief Integer for mass assignment
   */
  int mass_assignment;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string file_catalogue;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string file_catalogue_new_ref;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool MAS_correction;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  string Name_survey;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  int Number_of_random_files;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nnp_data;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  ULONG Nnp_window;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  int Unitsim_plabel;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_tracer_number_counts;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_tracer_mass_field;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_tracer_vmax_field;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_tracer_spin_field;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_local_mach_number;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_local_dach_number;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_local_overdensity;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_bias;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_relative_bias;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_quadratic_bias;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tracer_local_dm_density;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_PCA;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_tidal_anisotropy_at_halo;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_peak_height_at_halo;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Get_cell_local_mach_number;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec Scale_mach_number;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_prop_function_tracer;
  //////////////////////////////////////////////////////////
  /*
   * @brief 
   */  
  bool Get_prop_function_tracer_cwt;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool measure_cross;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool set_bins_equal_number_tracers;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool set_bins_equal_number_tracers_main_property;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int Number_of_bins_equal_number_tracers;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  int Number_of_bins_equal_number_tracers_main_property;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  ULONG Number_of_GRF;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  real_prec Kmax_FA;
  //////////////////////////////////////////////////////////
  /*
   * @brief Minimum value of x-coordinate
   */
  real_prec xmin;
  //////////////////////////////////////////////////////////
  /*
   * @brief Minimum value of y-coordinate
   */
  real_prec ymin;
  //////////////////////////////////////////////////////////
  /*
   * @brief Minimum value of z-coordinate
   */
  real_prec zmin;
  //////////////////////////////////////////////////////////
  /*
   * @brief Max value of x-coordinate
   */
  real_prec xmax;
  //////////////////////////////////////////////////////////
  /*
   * @brief Max value of y-coordinate
   */
  real_prec ymax;
  //////////////////////////////////////////////////////////
  /*
   * @brief Max value of z-coordinate
   */
  real_prec zmax;
  //////////////////////////////////////////////////////////
  /*
   * @brief (xmax+xmin)/2, mid point between min and max in x-coordinate. Not asked in param file
   */
  real_prec Xoffset;
  //////////////////////////////////////////////////////////
  /*
   * @brief (ymax+ymin)/2, mid point between min and max in x-coordinate. Not asked in param file
   */
  real_prec Yoffset;
  //////////////////////////////////////////////////////////
  /*
   * @brief (zmax+zmin)/2, mid point between min and max in x-coordinate. Not asked in param file
   */
  real_prec Zoffset;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool Generate_FA;
  //////////////////////////////////////////////////////////
  /*
   * @brief
   */
  bool use_low_pass_filter;
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////
  //-------------------now parameter for cosmolib ---------//
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec redshift_min;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec redshift_max;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int nbins_redshift;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string mass_function_fit;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmin_integration;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmax_integration;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec M_min_effective;
   //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
 real_prec M_max_effective;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec A_gas;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec B_gas;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec mstar;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int hod_model;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec muno_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec alpha_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec alpha_A;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec alpha_B;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec alpha_C;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec mmin_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec scatter_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec mt_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec ms_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec s_faint_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec s_bright_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec width_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Mstep_hod;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec sigma_ln;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec sigma_red;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec missing_flux;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string density_profile;
   //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec M_min_mf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec M_max_mf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string scale_mf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int npoints_mf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string mass_function_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string halo_mass_bias_fit;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string halo_mass_bias_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string effective_halo_mass_bias_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string effective_halo_mean_number_density_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool compute_output_linear_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool compute_output_non_linear_power_spectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string scale_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kstar;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec GAL_BIAS;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Amc;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string linear_matter_ps_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string non_linear_matter_ps_halo_fit_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string non_linear_matter_ps_pt_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec coef_concentration_amp;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec coef_concentration;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string galaxy_power_spectrum_halo_model_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string galaxy_correlation_function_halo_model_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string scale_cf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int npoints_cf;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string linear_matter_cf_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string non_linear_matter_cf_halo_fit_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool compute_output_linear_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool compute_output_non_linear_correlation_function;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool compute_density_profile;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec rmin_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec rmax_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string scale_dp_r;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int npoints_dp_r;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string density_profile_r_output_file;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmin_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmax_dp;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string scale_dp_k;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int npoints_dp_k;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string density_profile_k_output_file ;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool use_file_power;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string file_power;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool Get_SO_from_BN;

  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmin_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec kmax_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int npoints_ps;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int Nbins_hist;
  ///////////////////////////////////////////////////////////


  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************
  // **********************************************************************************

public:
  /**
   *  @brief default constructor
   *  @return object of class Parameters
   */

  Params(){}

  //  Params(){this->init_pars();}
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /**
   *  @brief constructor
   *  @param parameters_file parameter file
   *  @return object of class Parameters
   */

  Params(string _par_file):par_file (_par_file)
  {
    this->init_pars();
    this->read_pars(par_file);
    this->derived_pars();
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /**
   *  @brief default destructor
   *  @return none
   */
  ~Params(){}


  // These variables smust be made private, 
  // andnif thy need to be modified, de define a method, e.g.
  // void set_Nft(int new_Nft){this->Nft=new_Nft;}

  //////////////////////////////////////////////////////////
  /**
   * @brief Structure containing cosmological parameters
   **/

  s_CosmologicalParameters s_cosmo_pars;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief Initialize all the parameteres defined as private variables of the class Params.
   **/
  void init_pars();

  //////////////////////////////////////////////////////////
  /**
   * @brief Read input parameter file and allocate var as private variables of Params class.
   **/
  void read_pars(string );

  //////////////////////////////////////////////////////////
  /**
   * @brief
   **/
  void derived_pars();
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _NX(){return this->NX;}
  void set_NX(ULONG newNX){
#ifdef _FULL_VERBOSE_
    cout<<BLUE<<"Changing parameter NX from "<< this->NX<<" to "<<newNX<<endl;
#endif
    this->NX=newNX;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NY(){return this->NY;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NY_MASS(){return this->NY_MASS;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NY_SAT_FRAC(){return this->NY_SAT_FRAC;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _Nlambdath(){return this->Nlambdath;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get/set the value of the private member
   *  @return
   */
  string _Output_directory(){return this->Output_directory;}
  void set_Output_directory(string new_Output_directory){this->Output_directory=new_Output_directory;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_Y(){return this->Input_Directory_Y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_Y_TWO(){return this->Input_Directory_Y_TWO;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_Y(){return this->Name_Catalog_Y;}

  string _Name_Catalog_Y_new_ref(){return this->Name_Catalog_Y_new_ref;}

  string _Name_Catalog_Y_HR(){return this->Name_Catalog_Y_HR;}


  string _Name_Catalog_Y_MWEIGHTED(){return this->Name_Catalog_Y_MWEIGHTED;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_X(){return this->Input_Directory_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_BIAS_KERNEL(){return this->Input_Directory_BIAS_KERNEL;}
  //////////////////////////////////////////////////////////

  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_BIAS_KERNEL_TWO(){return this->Input_Directory_BIAS_KERNEL_TWO;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _Number_of_references(){return this->Number_of_references;}
  void set_Number_of_references(int new_number_of_references){this->Number_of_references=new_number_of_references;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  
  int _Number_of_new_mocks(){return this->Number_of_new_mocks;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_X_REF(){return this->Input_Directory_X_REF;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_X_REF_TWO(){return this->Input_Directory_X_REF_TWO;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_X_NEW(){return this->Input_Directory_X_NEW;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Input_Directory_X_new_ref(){return this->Input_Directory_X_new_ref;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _XNAME(){return this->XNAME;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_X(){return this->Name_Catalog_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_X_new_ref(){return this->Name_Catalog_X_new_ref;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_VelFieldx_X(){return this->Name_VelFieldx_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_VelFieldy_X(){return this->Name_VelFieldy_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_VelFieldz_X(){return this->Name_VelFieldz_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_X_REF_PDF(){return this->Name_Catalog_X_REF_PDF;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_X_NEW(){return this->Name_Catalog_X_NEW;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Catalog_X_NEW_TWO(){return this->Name_Catalog_X_NEW_TWO;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Property_X(){return this->Name_Property_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _YNAME(){return this->YNAME;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Name_Property_Y(){return this->Name_Property_Y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _iMAS_X(){return this->iMAS_X;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _iMAS_X_NEW(){return this->iMAS_X_NEW;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _iMAS_X_REF_PDF(){return this->iMAS_X_REF_PDF;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _iMAS_Y(){return this->iMAS_Y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _delta_Y_max(){return this->delta_Y_max;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _delta_Y_min(){return this->delta_Y_min;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _delta_X_max(){return this->delta_X_max;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _delta_X_min(){return this->delta_X_min;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _ldelta_Y_max(){return this->ldelta_Y_max;}
  void set_ldelta_Y_max(real_prec new_xmin){this->ldelta_Y_max=new_xmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _ldelta_Y_min(){return this->ldelta_Y_min;}
  void set_ldelta_Y_min(real_prec new_xmin){this->ldelta_Y_min=new_xmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _ldelta_X_max(){return this->ldelta_X_max;}
  void set_ldelta_X_max(real_prec new_xm){this->ldelta_X_max=new_xm;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _ldelta_X_min(){return this->ldelta_X_min;}
  void set_ldelta_X_min(real_prec new_xmin){this->ldelta_X_min=new_xmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Quantity(){return this->Quantity;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _NMASSbins(){return this->NMASSbins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _NPROPbins_bam(){return this->NPROPbins_bam;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
 
  ULONG _NMASSbins_mf(){return this->NMASSbins_mf;}
  void set_NMASSbins_mf(int new_NMASSbins_mf){this->NMASSbins_mf=new_NMASSbins_mf;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NMASSbins_power(){return this->MASSbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NVMAXbins_power(){return this->VMAXbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NVRMSbins_power(){return this->VRMSbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NRSbins_power(){return this->RSbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NRVIRbins_power(){return this->RVIRbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NCONCENTRATIONbins_power(){return this->CONCENTRATIONbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NSPINbins_power(){return this->SPINbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NSPINBULLOCKbins_power(){return this->SPINBULLOCKbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NMACHbins_power(){return this->MACHbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NLOCALDMbins_power(){return this->LOCALDMbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NDACHbins_power(){return this->DACHbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NBIASbins_power(){return this->BIASbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NRBIASbins_power(){return this->RBIASbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NQBIASbins_power(){return this->QBIASbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NLCbins_power(){return this->LCbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NTAbins_power(){return this->TAbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _NPHbins_power(){return this->PHbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NVIRIALbins_power(){return this->VIRIALbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NBTOAbins_power(){return this->BTOAbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  int _NCTOAbins_power(){return this->CTOAbins_max.size();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */

  //  int _NMASSbins_power(){return this->NMASSbins_power;}
  int _NMASScuts_power(){return this->MASScuts.size();}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _redshift(){return this->redshift;}
  void set_redshift(real_prec new_redshift){this->redshift=new_redshift;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _smscale(){return this->smscale;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _realization(){return this->realization;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Comp_conditional_PDF(){return this->Comp_conditional_PDF;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Comp_joint_PDF(){return this->Comp_joint_PDF;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _write_files_for_histograms(){return this->write_files_for_histograms;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Redefine_limits(){return this->Redefine_limits;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Convert_Density_to_Delta_X(){return this->Convert_Density_to_Delta_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Convert_Density_to_Delta_Y(){return this->Convert_Density_to_Delta_Y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  real_prec _lambdath(){return this->lambdath;}
  void set_lambdath(real_prec new_lambdath){this->lambdath=new_lambdath;}

  real_prec _lambdath_v(){return this->lambdath_v;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Write_Scatter_Plot(){return this->Write_Scatter_Plot;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Write_PDF_number_counts(){return this->Write_PDF_number_counts;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Scale_X(){return this->Scale_X;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _Scale_Y(){return this->Scale_Y;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _N_dm_initial(){return this->N_dm_initial;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _n_sknot_massbin(){return this-> n_sknot_massbin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _n_vknot_massbin(){return this-> n_vknot_massbin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _N_iterations_Kernel(){return this->N_iterations_Kernel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _Iteration_Kernel_average(){return this->Iteration_Kernel_average;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  void set_N_iterations_Kernel(int new_N_iterations_Kernel){this->N_iterations_Kernel=new_N_iterations_Kernel;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Apply_Rankordering(){return this->Apply_Rankordering;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _Apply_Rankordering_ab_initio(){return this->Apply_Rankordering_ab_initio;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  vector<ULONG> _cwt_used(){return this->cwt_used;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  vector<ULONG> _cwv_used(){return this->cwv_used;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  vector<ULONG> _output_at_iteration(){return this->output_at_iteration;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n_cwt
   *  @return
   */
  ULONG _n_cwt(){return this->n_cwt;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n_cwv
   *  @return
   */
  ULONG _n_cwv(){return this->n_cwv;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n_cwv
   *  @return
   */
  int _NRANDOMfiles(){return this->NRANDOMfiles;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n_cwv
   *  @return
   */
  string _RANDOMfile(int i){return this->RANDOMfiles[i];}



  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //   Parameters of COSMOLIB    ///////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  // function to get private variables
  /**
   *  @brief get the value of the private member statistics
   *  @return statistics
   */
  string _statistics () {return statistics;}
  void set_statistics (string new_stats) {this->statistics=new_stats;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Input_dir_cat
   *  @return Input_dir_cat
   */
  string _Input_dir_cat () {return this->Input_dir_cat;}
  string _Input_dir_cat_new_ref () {return this->Input_dir_cat_new_ref;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Input_dir_cat
   *  @return Input_dir_cat
   */
  string _Input_dir_cat_TWO () {return this->Input_dir_cat_TWO;}
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
   **/
  void set_measure_cross(bool new_meas){this->measure_cross=new_meas;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _dir_output(){return this->dir_output;}
  void set_dir_output(string new_dir_output){this->dir_output=new_dir_output;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _input_type () {return input_type;}
  void set_input_type(string new_input_type){input_type=new_input_type;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */

  string _input_type_two () {return input_type_two;}
  void set_input_type_two(string new_input_type){input_type_two=new_input_type;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  ULONG _ngal_delta() {return ngal_delta;}
  void set_ngal_delta(ULONG new_ngal_delta) {this->ngal_delta=new_ngal_delta;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file () {return delta_grid_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file2 () {return delta_grid_file2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file3 () {return delta_grid_file3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  string _delta_grid_file4 () {return delta_grid_file4;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  bool _measure_cross () {return measure_cross;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _measure_cross_from_1 () {return measure_cross_from_1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  @return
   */
  int _measure_cross_from_2 () {return measure_cross_from_2;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief set the value of the private member Input_dir_cat
   *  @param _Input_dir_cat input directory where the object and random
   *  catalogues are stored
   *  @return none
   */
  void set_Input_dir_cat (string _Input_dir_cat) {Input_dir_cat = _Input_dir_cat;}
  void set_Input_dir_cat_new_ref (string _Input_dir_cat) {Input_dir_cat_new_ref = _Input_dir_cat;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_catalogue
   *  @return file_catalogue
   */
  string _file_catalogue () {return file_catalogue;}
  void set_file_catalogue (string new_file_catalogue) {file_catalogue=new_file_catalogue;}


  string _file_catalogue_new_ref () {return file_catalogue_new_ref;}
  void set_file_catalogue_new_ref (string new_file_catalogue) {file_catalogue_new_ref=new_file_catalogue;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_random
   *  @return file_random
   */
  string _file_random () {return file_random;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_Lbox
   *  @return new_Lbox
   */
  bool _new_Lbox () {return new_Lbox;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_g
   *  @return sys_of_coord_g
   */
  int _sys_of_coord_g () {return sys_of_coord_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_g
   *  @return sys_of_coord_dm
   */
  int _sys_of_coord_dm () {return sys_of_coord_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_coord1_g () {return this->i_coord1_g;}
  void set_i_coord1_g(int new_i_coord1_g){this->i_coord1_g=new_i_coord1_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_g () {return i_coord2_g;}
  void set_i_coord2_g(int new_i_coord2_g){this->i_coord2_g=new_i_coord2_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_g () {return this->i_coord3_g;}
  void set_i_coord3_g(int new_i_coord3_g){this->i_coord3_g=new_i_coord3_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_v1_g () {return i_v1_g;}
  void set_i_v1_g(int new_i_v1_g){this->i_v1_g=new_i_v1_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_v2_g () {return i_v2_g;}
  void set_i_v2_g(int new_i_v2_g){this->i_v2_g=new_i_v2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_v3_g () {return this->i_v3_g;}
  void set_i_v3_g(int new_i_v3_g){this->i_v3_g=new_i_v3_g;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_mass_g () {return this->i_mass_g;}
  void set_i_mass_g (int new_i_mass_g) {this->i_mass_g=new_i_mass_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_vmax_g () {return this->i_vmax_g;}
  void set_i_vmax_g (int new_i_vmax_g) {this->i_vmax_g=new_i_vmax_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_vrms_g () {return this->i_vrms_g;}
  void set_i_vrms_g (int new_i_vrms_g) {this->i_vrms_g=new_i_vrms_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_sf_g () {return this->i_sf_g;}
  void set_i_sf_g (int new_i_sf_g) {this->i_sf_g=new_i_sf_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_rs_g () {return this->i_rs_g;}
  void set_i_rs_g (int new_i_rs_g) {this->i_rs_g=new_i_rs_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_rvir_g () {return this->i_rvir_g;}
  void set_i_rvir_g (int new_i_rvir_g) {this->i_rvir_g=new_i_rvir_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_virial_g () {return this->i_virial_g;}
  void set_virial_g (int new_i_virial_g) {this->i_virial_g=new_i_virial_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_spin_g () {return this->i_spin_g;}
  void set_i_spin_g (int new_i_spin_g) {this->i_spin_g=new_i_spin_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_spin_bullock_g () {return this->i_spin_bullock_g;}
  void set_i_spin_bullock_g (int new_i_spin_g) {this->i_spin_bullock_g=new_i_spin_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_b_to_a_g () {return this->i_b_to_a_g;}
  void set_i_b_to_a_g (int new_i_g) {this->i_b_to_a_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _i_c_to_a_g () {return this->i_c_to_a_g;}
  void set_i_c_to_a_g (int new_i_g) {this->i_c_to_a_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_color_g () {return this->i_color_g;}
  void set_i_color_g (int new_i_g) {this->i_color_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_stellar_mass_g () {return this->i_stellar_mass_g;}
  void set_i_stellar_mass_g (int new_i_g) {this->i_stellar_mass_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_app_mag_g () {return this->i_app_mag_g;}
  void set_i_app_mag_g (int new_i_g) {this->i_app_mag_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_abs_mag_g () {return this->i_abs_mag_g;}
  void set_i_abs_mag_g (int new_i_g) {this->i_abs_mag_g=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_mass_dm () {return this->i_mass_dm;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_g
   *  @return i_weight1_g
   */
  int _i_weight1_g () {return this->i_weight1_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_g
   *  @return i_weight2_g
   */
  int _i_weight2_g () {return this->i_weight2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_g
   *  @return i_weight3_g
   */
  int _i_weight3_g () {return this->i_weight3_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight4_g
   *  @return i_weight4_g
   */
  int _i_weight4_g () {return i_weight4_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_g
   *  @return use_weight1_g
   */
  bool _use_weight1_g () {return this->use_weight1_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight2_g
   *  @return use_weight2_g
   */
  bool _use_weight2_g () {return this->use_weight2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight3_g
   *  @return use_weight3_g
   */
  bool _use_weight3_g () {return this->use_weight3_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight4_g
   *  @return use_weight4_g
   */
  bool _use_weight4_g () {return this->use_weight4_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _weight_vel_with_mass () {return this->weight_vel_with_mass;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_coord1_dm () {return i_coord1_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_dm () {return i_coord2_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_dm () {return this->i_coord3_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_v1_dm () {return i_v1_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_v2_dm () {return i_v2_dm;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_v3_dm() {return this->i_v3_dm;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Name_survey
   *  @return Name_survey
   */
  string _Name_survey () {return this->Name_survey;}
  void set_Name_survey(string new_Name_survey){this->Name_survey=new_Name_survey;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_mean_density_g
   *  @return i_mean_density_g
   */
  int _i_mean_density_g () {return i_mean_density_g;}
  void set_i_mean_density_g (int in) {this->i_mean_density_g=in;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member angles_units_g
   *  @return angles_units_g
   */
  string _angles_units_g () {return angles_units_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog () {return use_random_catalog;}
  void set_use_random_catalog (bool new_use_random_catalog) {this->use_random_catalog = new_use_random_catalog ;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog_cl () {return use_random_catalog_cl;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_los
   *  @return new_los
   */
  bool _new_los () {return new_los;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_r
   *  @return sys_of_coord_r
   */
  int _sys_of_coord_r () {return sys_of_coord_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_r
   *  @return i_coord1_r
   */
  int _i_coord1_r () {return i_coord1_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_r
   *  @return i_coord2_r
   */
  int _i_coord2_r () {return i_coord2_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_r
   *  @return i_coord3_r
   */
  int _i_coord3_r () {return i_coord3_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_r
   *  @return i_weight1_r
   */
  int _i_weight1_r () {return i_weight1_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_r
   *  @return i_weight2_r
   */
  int _i_weight2_r () {return i_weight2_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_r
   *  @return i_weight3_r
   */
  int _i_weight3_r () {return i_weight3_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight4_r
   *  @return i_weight4_r
   */
  int _i_weight4_r () {return i_weight4_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_mass_r
   *  @return i_mass_r
   */
  int _i_mass_r () {return i_mass_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_r
   *  @return use_weight1_r
   */
  bool _use_weight1_r () {return use_weight1_r;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_r
   *  @return use_weight1_r
   */

  bool _weight_with_mass () {return this->weight_with_mass;}
  void set_weight_with_mass(bool new_weight_with_mass){this->weight_with_mass=new_weight_with_mass;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight2_r
   *  @return use_weight2_r
   */
  bool _use_weight2_r () {return use_weight2_r;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight3_r
   *  @return use_weight3_r
   */
  bool _use_weight3_r () {return use_weight3_r;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight4_r
   *  @return use_weight4_r
   */
  bool _use_weight4_r () {return use_weight4_r;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_mean_density_r
   *  @return i_mean_density_r
   */
  int _i_mean_density_r () {return i_mean_density_r;}
  void set_i_mean_density_r (int in) {this->i_mean_density_r=in;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_color_r () {return this->i_color_r;}
  void set_i_color_r (int new_i_g) {this->i_color_r=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_stellar_mass_r () {return this->i_stellar_mass_r;}
  void set_i_stellar_mass_r (int new_i_g) {this->i_stellar_mass_r=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_app_mag_r () {return this->i_app_mag_r;}
  void set_i_app_mag_r (int new_i_g) {this->i_app_mag_r=new_i_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _i_abs_mag_r () {return this->i_abs_mag_r;}
  void set_i_abs_mag_r (int new_i_g) {this->i_abs_mag_r=new_i_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member angles_units_r
   *  @return angles_units_r
   */
  string _angles_units_r () {return angles_units_r;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n_catalogues
   *  @return n_catalogues
   */
  int _n_catalogues () {return n_catalogues;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get / set the value of the private member Nft
   *  @return Nft
   */
  ULONG _Nft () {return this->Nft;}
  //////////////////////////////////////////////////////////

  void set_Nft(ULONG new_Nft){
    this->Nft=new_Nft;
    this->NGRID=new_Nft*new_Nft*new_Nft;
    this->NGRID_h=new_Nft*new_Nft*(new_Nft/2+1);
    this->Nnp_data   = (this->type_of_binning=="log" ? this->N_log_bins : this->Nft/this->ndel_data/2); 
  }
 
  //////////////////////////////////////////////////////////
 
  ULONG d_NGRID(){return this->Nft*this->Nft*this->Nft;}
 
  //////////////////////////////////////////////////////////
 
  ULONG d_NGRID_JK(){return this->Nft_JK*this->Nft_JK*this->Nft_JK;}
 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_low () {return this->Nft_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_HR () {return this->Nft_HR;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_JK () {return this->Nft_JK;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nft_random_collapse () {return this->Nft_random_collapse;}
  void set_Nft_random_collapse (int new_Nft_random_collapse) {this->Nft_random_collapse=new_Nft_random_collapse;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Distance_fraction() {return this->Distance_fraction;}
  void set_Distance_fraction(real_prec new_Distance_fraction){this->Distance_fraction=new_Distance_fraction;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nf   */
  ULONG _NGRID(){return this->NGRID;}
  void set_NGRID(ULONG new_NGRID){this->NGRID=new_NGRID;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID_HR(){return this->NGRID_HR;}
  void set_NGRID_HR(ULONG new_NGRID){this->NGRID_HR=new_NGRID;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID_JK(){return this->NGRID_JK;}
  void set_NGRID_JK(ULONG new_NGRID){this->NGRID_JK=new_NGRID;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute size of grid
   *  @return Nft
   */
  ULONG _NGRID_h(){return this->NGRID_h;}
  void set_NGRID_h(ULONG new_NGRID_h){this->NGRID_h=new_NGRID_h;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  ULONG _NGRID_low(){return this->NGRID_low;}
  void set_NGRID_low(ULONG new_NGRID_l){this->NGRID_low=new_NGRID_l;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Lbox
   *  @return Lbox
   */
  real_prec _Lbox () {return this->Lbox;}
  void set_Lbox(real_prec new_Lbox){this->Lbox=new_Lbox;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Lbox_low () {return this->Lbox_low;}
  void set_Lbox_low(real_prec new_Lbox_low){this->Lbox_low=new_Lbox_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member mass_assignment_scheme
   *  @return mass_assignment_scheme
   */
  string _mass_assignment_scheme () {return this->mass_assignment_scheme;}
  void set_mass_assignment_scheme (string new_mass_assignment_scheme) {this->mass_assignment_scheme=new_mass_assignment_scheme;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member mass_assignment_scheme
   *  @return mass_assignment_scheme
   */
  int _mass_assignment () {return this->mass_assignment;}
  void set_mass_assignment_scheme (int new_mass_assignment_scheme) {this->mass_assignment=new_mass_assignment_scheme;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _type_of_binning () {return this->type_of_binning;}
  void set_type_of_binning(string new_type_of_binning){this->type_of_binning=new_type_of_binning;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */

  real_prec _vkernel_exponent(){return this->vkernel_exponent;}

  //////////////////////////////////////////////////////////
  /**
   * @brief Identifies dark matter "DM" or tracer "TR"
   */
  void set_type_of_object(string new_type_of_object){this->type_of_object=new_type_of_object;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _extra_info () {return this->extra_info;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N_log_bins
   *  @return N_log_bins
   */
  int _N_log_bins () {return this->N_log_bins;}
  void set_N_log_bins(int new_N_log_bins){this->N_log_bins=new_N_log_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member ndel_data
   *  @return ndel_data
   */
  int _ndel_data () {return this->ndel_data;}
  void set_ndel_data(int new_ndel_data){this->ndel_data=new_ndel_data;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member ndel_window
   *  @return ndel_window
   */
  int _ndel_window () {return this->ndel_window;}
  void set_ndel_window(int new_ndel_window){this->ndel_data=new_ndel_window;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N_mu_bins
   *  @return N_mu_bins
   */
  int _N_mu_bins () {return this->N_mu_bins;}
  void set_N_mu_bins(int new_N_mu_bins){this->N_mu_bins=new_N_mu_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member MAS_correction
   *  @return MAS_correction
   */
  bool _MAS_correction () {return this->MAS_correction;}
  void set_MAS_correction(bool new_MAS_correction){this->MAS_correction=new_MAS_correction;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member FKP_weight
   *  @return FKP_weight
   */
  bool _FKP_weight () {return FKP_weight;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get/set the value of the private member SN_correction
   *  @return SN_correction
   */
  bool _SN_correction () {return SN_correction;}
  void set_SN_correction (bool new_sn_correction) {this->SN_correction=new_sn_correction;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member FKP_error_bars
   *  @return FKP_error_bars
   */
  bool _FKP_error_bars () {return FKP_error_bars;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member FKP_error_bars_exact
   *  @return FKP_error_bars_exact
   */
  bool _FKP_error_bars_exact () {return FKP_error_bars_exact;}
  ///////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Pest
   *  @return Pest
   */
  real_prec _Pest () {return this->Pest;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member nbar_tabulated
   *  @return nbar_tabulated
   */
  bool _nbar_tabulated () {return this->nbar_tabulated;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member nbar_tabulated
   *  @return nbar_tabulated
   */
  bool _use_file_nbar () {return this->use_file_nbar;}
  ///////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member constant_depth
   *  @return constant_depth
   */
  bool _constant_depth () {return constant_depth;}
  //////////////////////////////////////////////////////////

  /**
   *  @brief get the value of the private member N_z_bins
   *  @return N_z_bins
   */
  int _Nbins_redshift () {return Nbins_redshift;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member redshift_min_sample
   *  @return redshift_min_sample
   */
  real_prec _redshift_min_sample () {return redshift_min_sample;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member redshift_max_sample
   *  @return redshift_max_sample
   */
  real_prec _redshift_max_sample () {return redshift_max_sample;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N_dndz_bins
   *  @return N_dndz_bins
   */
  int _N_dndz_bins () {return N_dndz_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_N_dndz_bins
   *  @return new_N_dndz_bins
   */
  int _new_N_dndz_bins () {return new_N_dndz_bins;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Healpix_resolution
   *  @return Healpix_resolution
   */
  int _Healpix_resolution () {return Healpix_resolution;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_dndz
   *  @return file_dndz
   */
  string _file_dndz () {return this->file_dndz;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_dndz
   *  @return file_dndz
   */
  string _file_nbar () {return this->file_nbar;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_power
   *  @return file_power
   */
  string _file_power () {return this->file_power;}
  void set_file_power(string new_file_power){this->file_power=new_file_power;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_power_log
   *  @return file_power_log
   */
  string _file_power_log () {return file_power_log;}
  //////////////////////////////////////////////////////////

  /**
   *  @brief get the value of the private member file_window
   *  @return file_window
   */
  string _file_window () {return file_window;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_power2d
   *  @return file_power2d
   */
  string _file_power2d () {return file_power2d;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_power2d_mk
   *  @return file_power2d_mk
   */
  string _file_power2d_mk () {return file_power2d_mk;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmin_bk
   */
  real_prec _kmin_bk () {return kmin_bk;}
  void set_kmin_bk(real_prec new_kmin_bk){kmin_bk=new_kmin_bk;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmax_bk
   */
  real_prec _kmax_bk () {return kmax_bk;}
  void set_kmax_bk(real_prec new_kmax_bk){kmax_bk=new_kmax_bk;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_fundamental_mode_as_kmin_bk
   *  @return Nshells_bk
   */
  bool _use_fundamental_mode_as_kmin_bk () {return use_fundamental_mode_as_kmin_bk;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_bispectrum
   *  @return file_bispectrum
   */
  string _file_bispectrum () {return file_bispectrum;}



  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member kmax_y_ds
   *  @return kmax_y_ds
   */
  real_prec _kmax_y_ds () {return kmax_y_ds;}
  void set_kmax_y_ds (real_prec new_kmax_y_ds ){kmax_y_ds=new_kmax_y_ds;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _DeltaKmin () {return this->DeltaKmin;}
  void set_DeltaKmin(real_prec new_DeltaKmin){this->DeltaKmin=new_DeltaKmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member input_file_mask
   *  @return input_file_mask
   */
  string _input_file_mask(){return input_file_mask;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member i_mask_pixel
   *  @return i_mask_pixel
   */
  int _i_mask_pixel(){return i_mask_pixel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member i_mask_alpha
   *  @return  i_mask_alpha
   */
  int _i_mask_alpha(){return i_mask_delta;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member i_mask_delta
   *  @return  i_mask_delta
   */
  int _i_mask_delta(){return i_mask_delta;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member i_mask_flag
   *  @return  i_mask_flag
   */
  int _i_mask_flag(){return i_mask_flag;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member
   *  @return
   */
  string _Name_redshift_mask(){return this->Name_redshift_mask;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the the private member
   *  @return
   */
  string _Name_binary_mask(){return this->Name_binary_mask;}



  //////////////////////////////////////////////////////////
  /**
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
   **/
  string input_type;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  string input_type_two;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  ////COSMOLOGICAL PARAMETERS///////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////



  /**
   *  @brief get the value of the private member &Omega;<SUB>M</SUB>
   *  @return &Omega;<SUB>M</SUB>
   */
  real_prec _om_matter () { return om_matter; }

  /**
   *  @brief get the value of the private member &Omega;<SUB>M</SUB>
   *  @return &Omega;<SUB>M</SUB>
   */
  real_prec _om_cdm () { return om_cdm; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &Omega;<SUB>rad</SUB>: the radiation density at z=0
   *  @return &Omega;<SUB>rad</SUB>
   */
  real_prec _om_radiation () { return om_radiation; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &Omega;<SUB>b</SUB>: the baryon density at z=0
   *  @return &Omega;<SUB>b</SUB>
   */
  real_prec _om_baryons () { return om_baryons; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &Omega;<SUB>DE</SUB>: the dark energy density at z=0
   *  @return &Omega;<SUB>DE</SUB>
   */
  real_prec _om_vac () { return om_vac; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &Omega;<SUB>k</SUB>: the density of curvature energy
   *  @return &Omega;<SUB>k</SUB>
   */
  real_prec _om_k () { return om_k; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc]
   *  @return H<SUB>0</SUB>
   */
  real_prec _Hubble () { return Hubble; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member \e h: the Hubble parameter, H<SUB>0</SUB>/100
   *  @return \e h
   */
  real_prec _hubble () { return hubble; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n<SUB>spec</SUB>: the primordial spectral index
   *  @return n<SUB>spec</SUB>
   */
  real_prec _spectral_index () { return spectral_index; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n<SUB>spec</SUB>: f_baryon
   *  @return n<SUB>spec</SUB>
   */
  real_prec _f_baryon () { return this->f_baryon;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member n<SUB>spec</SUB>: f_baryon
   *  @return n<SUB>spec</SUB>
   */
  real_prec _A_s () { return this->A_s;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member w<SUB>0</SUB>: the parameter of the dark energy equation of state
   *  @return w<SUB>0</SUB>
   */
  real_prec _wde_eos () { return wde_eos; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
   *  @return N<SUB>eff</SUB>
   */
  real_prec _N_eff () { return N_eff; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member &sigma;<SUB>8</SUB>: the power spectrum normalization
   *  @return &sigma;<SUB>8</SUB>
   */
  real_prec _sigma8 () { return sigma8; }

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member T<SUB>CMB</SUB>: the present day CMB temperature [K]
   *  @return T<SUB>CMB</SUB>
   */
  real_prec _Tcmb () { return Tcmb; }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _use_wiggles () { return use_wiggles; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _redshift_space_coords_g () { return this->redshift_space_coords_g; }
  void  set_redshift_space_coords_g (bool new_redshift_space_coords_g) { this->redshift_space_coords_g=new_redshift_space_coords_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string _vel_units_g () { return this->vel_units_g; }
  void set_vel_units_g (string new_vel_units_g) { this->vel_units_g=new_vel_units_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _RR () { return RR; }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Params useful for Patchy. OJO QUE EN PARAMS.cpp aun no estan leyendo muchos de estos
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _inputmode(){return this->inputmode;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _dataFileName(){return this->dataFileName;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_power_file(){return this->ic_power_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_WN_file(){return this->ic_WN_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_file(){return this->ic_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_input_type(){return this->ic_input_type;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _use_ic_file(){return this->use_ic_file;}
  void set_use_ic_file(bool newi){this->use_ic_file=newi;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _ic_alias_corrected(){return this->ic_alias_corrected;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _ic_WN_dir(){return this->ic_WN_dir;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _dir(){return this->dir;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Structure_Formation_Model(){return this->Structure_Formation_Model;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _Normalize_IC_to_initial_redshift(){return this->Normalize_IC_to_initial_redshift;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _runsim(){return runsim;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _runv(){return runv;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _readPS(){return readPS;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _use_vel_kernel(){return this->use_vel_kernel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Nmean(){return Nmean;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cs2(){return cs2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cdeltas2(){return cdeltas2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cpsi(){return cpsi;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cst(){return cst;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _cs3(){return cs3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _sfac(){return sfac;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _ep(){return ep;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _devpois(){return devpois;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _deltathH(){return deltathH;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _N1(){return this->Nft;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _L1(){return this->Lbox;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasrhoexp(){return this->biasrhoexp;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasepsilon(){return this->biasepsilon;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasrhoexp2(){return this->biasrhoexp2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasepsilon2(){return this->biasepsilon2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biassign(){return this->biassign;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biassign2(){return this->biassign2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasone(){return this->biasone;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias(){return this->velbias;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias_dm(){return this->velbias_dm;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _velbias_random(){return this->velbias_random;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _biasL(){return this->biasL;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _xllc(){return this->xllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _yllc(){return this->yllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _zllc(){return this->zllc;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _xobs(){return this->xobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _yobs(){return this->yobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _zobs(){return this->zobs;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _seed(){return this->seed;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _seed_ref(){return this->seed_ref;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _masskernel(){return this->masskernel;}
  void set_masskernel(int new_masskernel){this->masskernel=new_masskernel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _masskernel_vel(){return this->masskernel_vel;}
  void set_masskernel_vel(int new_masskernel_vel){this->masskernel_vel=new_masskernel_vel;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _slength(){return this->slength;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _slengthv(){return this->slengthv;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _vslength(){return this->vslength;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_x_coord(){return this->file_bin_x_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_y_coord(){return this->file_bin_y_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _file_bin_z_coord(){return this->file_bin_z_coord;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _dilute_dm_sample(){return this->dilute_dm_sample;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _fraction_dilute(){return this->fraction_dilute;}



  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG _N_lines_binary(){return  this->N_lines_binary;}
  void set_N_lines_binary(ULONG newL){this->N_lines_binary=newL;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  bool _get_distribution_min_separations(){return this->get_distribution_min_separations;} 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_chunks_new_dm(){return  this->Number_of_chunks_new_dm;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d1(){return this->d1;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d2(){return this->d2;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d3(){return this->d3;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d1_HR(){return this->d1_HR;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d2_HR(){return this->d2_HR;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d3_HR(){return this->d3_HR;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d1_JK(){return this->d1_JK;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d2_JK(){return this->d2_JK;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _d3_JK(){return this->d3_JK;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d1_low(){return this->d1_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d2_low(){return this->d2_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d3_low(){return this->d3_low;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _M_exclusion(){return this->M_exclusion;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _LOGMASSmin(){return this->LOGMASSmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _LOGMASSmax(){return this->LOGMASSmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXmin(){return this->VMAXmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXmax(){return this->VMAXmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VRMSmin(){return this->VRMSmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VRMSmax(){return this->VRMSmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _RSmin(){return this->RSmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _RSmax(){return this->RSmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CONCENTRATIONmin(){return this->CONCENTRATIONmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CONCENTRATIONmax(){return this->CONCENTRATIONmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINmin(){return this->SPINmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINmax(){return this->SPINmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VIRIALmin(){return this->VIRIALmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VIRIALmax(){return this->VIRIALmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _BTOAmin(){return this->BTOAmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _BTOAmax(){return this->BTOAmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CTOAmin(){return this->CTOAmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CTOAmax(){return this->CTOAmax;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _MASS_units(){return this->MASS_units;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Initial_Redshift_DELTA(){return this->Initial_Redshift_DELTA;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Initial_Redshift_SIM(){return this->Initial_Redshift_SIM;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _Initial_Redshift_TH_power_file(){return this->Initial_Redshift_TH_power_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  int _IC_index(){return this->IC_index;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  string _files_new_dm_fields(int i){
    return this->files_new_dm_fields[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _type_of_object(){
    return this->type_of_object;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_dm_references(int i){
    return this->files_dm_references[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */

  string _files_tracer_references(int i){
    return this->files_tracer_references[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_tracer_field_references(int i){
    return this->files_tracer_field_references[i];
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _list_new_dm_fields(int i){
    return this->list_new_dm_fields[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_bias_references(int i){
    return this->files_bias_references[i];
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _MASSbins_min(int i){
    return this->MASSbins_min[i];}
  void set_MASSbins_min(int i,real_prec value){
    this->check_dim(i,this->MASSbins_min.size(), __PRETTY_FUNCTION__);
    this->MASSbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _MASSbins_max(int i){
    return this->MASSbins_max[i];}
  void set_MASSbins_max(int i,real_prec value){
    this->check_dim(i,this->MASSbins_max.size(), __PRETTY_FUNCTION__);
    this->MASSbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXbins_min(int i){
    return this->VMAXbins_min[i];}

  void set_VMAXbins_min(int i,real_prec value){
    this->check_dim(i,this->VMAXbins_min.size(), __PRETTY_FUNCTION__);
    this->VMAXbins_min[i]=value;
  }



  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VMAXbins_max(int i){
    return this->VMAXbins_max[i];}

  void set_VMAXbins_max(int i,real_prec value){
    this->check_dim(i,this->VMAXbins_max.size(), __PRETTY_FUNCTION__);
    this->VMAXbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VRMSbins_min(int i){
    return this->VRMSbins_min[i];}
  void set_VRMSbins_min(int i,real_prec value){
    this->check_dim(i,this->VRMSbins_min.size(), __PRETTY_FUNCTION__);
    this->VRMSbins_min[i]=value;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VRMSbins_max(int i){
    return this->VRMSbins_max[i];}
  void set_VRMSbins_max(int i,real_prec value){
    this->check_dim(i,this->VRMSbins_max.size(), __PRETTY_FUNCTION__);
    this->VRMSbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RVIRbins_max(int i){
    return this->RVIRbins_max[i];}
  void set_RVIRbins_max(int i,real_prec value){
    this->check_dim(i,this->RVIRbins_max.size(), __PRETTY_FUNCTION__);
    this->RVIRbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RVIRbins_min(int i){
    return this->RVIRbins_min[i];}
  void set_RVIRbins_min(int i,real_prec value){
    this->check_dim(i,this->RVIRbins_min.size(), __PRETTY_FUNCTION__);
    this->RVIRbins_min[i]=value;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RSbins_min(int i){
    return this->RSbins_min[i];}
  void set_RSbins_min(int i,real_prec value){
    this->check_dim(i,this->RSbins_min.size(), __PRETTY_FUNCTION__);
    this->RSbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RSbins_max(int i){
    return this->RSbins_max[i];}
  void set_RSbins_max(int i,real_prec value){
    this->check_dim(i,this->RSbins_max.size(), __PRETTY_FUNCTION__);
    this->RSbins_max[i]=value;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _CONCENTRATIONbins_max(int i){
    return this->CONCENTRATIONbins_max[i];}
  void set_CONCENTRATIONbins_max(int i,real_prec value){
    if(i>=this->CONCENTRATIONbins_max.size())
      { cout<<"Error in "<<__PRETTY_FUNCTION__<<endl; cout<<"Requested space at i ="<<i<<endl;}
    this->CONCENTRATIONbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _CONCENTRATIONbins_min(int i){
    return this->CONCENTRATIONbins_min[i];}
  void set_CONCENTRATIONbins_min(int i,real_prec value){
    this->check_dim(i,this->CONCENTRATIONbins_min.size(), __PRETTY_FUNCTION__);
    this->CONCENTRATIONbins_min[i]=value;
  }


  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINbins_max(int i){
    return this->SPINbins_max[i];}
  void set_SPINbins_max(int i,real_prec value){
    this->check_dim(i,this->SPINbins_max.size(), __PRETTY_FUNCTION__);
    this->SPINbins_max[i]=value;
  }
 
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _SPINbins_min(int i){
    return this->SPINbins_min[i];}
  void set_SPINbins_min(int i,real_prec value){
    this->check_dim(i,this->SPINbins_min.size(), __PRETTY_FUNCTION__);
    this->SPINbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _SPINBULLOCKbins_max(int i){
    return this->SPINBULLOCKbins_max[i];}
  void set_SPINBULLOCKbins_max(int i,real_prec value){
    this->check_dim(i,this->SPINBULLOCKbins_max.size(), __PRETTY_FUNCTION__);
    this->SPINBULLOCKbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _SPINBULLOCKbins_min(int i){
    return this->SPINBULLOCKbins_min[i];}
  void set_SPINBULLOCKbins_min(int i,real_prec value){
    this->check_dim(i,this->SPINBULLOCKbins_min.size(), __PRETTY_FUNCTION__);
    this->SPINBULLOCKbins_min[i]=value;
  }


  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VIRIALbins_min(int i){
    return this->VIRIALbins_min[i];}
  void set_VIRIALbins_min(int i,real_prec value){
    this->check_dim(i,this->VIRIALbins_min.size(), __PRETTY_FUNCTION__);
    this->VIRIALbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _VIRIALbins_max(int i){
    return this->VIRIALbins_max[i];}
  void set_VIRIALbins_max(int i,real_prec value){
    this->check_dim(i,this->VIRIALbins_max.size(), __PRETTY_FUNCTION__);
    this->VIRIALbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _BTOAbins_min(int i){
    return this->BTOAbins_min[i];}
  void set_BTOAbins_min(int i,real_prec value){
    this->check_dim(i,this->BTOAbins_min.size(), __PRETTY_FUNCTION__);
    this->BTOAbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _BTOAbins_max(int i){
    return this->BTOAbins_max[i];}
  void set_BTOAbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->BTOAbins_max.size(), __PRETTY_FUNCTION__);
    this->BTOAbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CTOAbins_min(int i){
    return this->CTOAbins_min[i];}
  void set_CTOAbins_min(int i,real_prec value){
    this->check_dim(i,this->CTOAbins_min.size(), __PRETTY_FUNCTION__);
    this->CTOAbins_min[i]=value;
  }

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _CTOAbins_max(int i){
    return this->CTOAbins_max[i];}
  void set_CTOAbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->CTOAbins_max.size(), __PRETTY_FUNCTION__);
    this->CTOAbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _MACHbins_max(int i){
    return this->MACHbins_max[i];}
  void set_MACHbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->MACHbins_max.size(), __PRETTY_FUNCTION__);
    this->MACHbins_max[i]=value;
  }
  real_prec _MACHbins_max_last(){return this->MACHbins_max.back();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _MACHbins_min(int i){
    return this->MACHbins_min[i];}

  void set_MACHbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->MACHbins_min.size(), __PRETTY_FUNCTION__);
    this->MACHbins_min[i]=value;
  }
  real_prec _MACHbins_min_begin(){return this->MACHbins_min.front();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _DACHbins_max(int i){
    return this->DACHbins_max[i];}

  void set_DACHbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->DACHbins_max.size(), __PRETTY_FUNCTION__);
    this->DACHbins_max[i]=value;
  }
  real_prec _DACHbins_max_last(){return this->DACHbins_max.back();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _DACHbins_min(int i){
    return this->DACHbins_min[i];
  }
  void set_DACHbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->DACHbins_min.size(), __PRETTY_FUNCTION__);
    this->DACHbins_min[i]=value;
  }
  real_prec _DACHbins_min_begin(){return this->DACHbins_min.front();}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _BIASbins_max(int i)
  {
    return this->BIASbins_max[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RBIASbins_max(int i)
  {
    return this->RBIASbins_max[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _QBIASbins_max(int i)
  {
    return this->QBIASbins_max[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_BIASbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->BIASbins_max.size(), __PRETTY_FUNCTION__);
    this->BIASbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_RBIASbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->RBIASbins_max.size(), __PRETTY_FUNCTION__);
    this->RBIASbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_QBIASbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->QBIASbins_max.size(), __PRETTY_FUNCTION__);
    this->QBIASbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _BIASbins_min(int i)
  {
    return this->BIASbins_min[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RBIASbins_min(int i)
  {
    return this->RBIASbins_min[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _QBIASbins_min(int i)
  {
    return this->QBIASbins_min[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_BIASbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->BIASbins_min.size(), __PRETTY_FUNCTION__);
    this->BIASbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_RBIASbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->RBIASbins_min.size(), __PRETTY_FUNCTION__);
    this->RBIASbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_QBIASbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->QBIASbins_min.size(), __PRETTY_FUNCTION__);
    this->QBIASbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _MASScuts(int i){
    return this->MASScuts[i];
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_LCbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->LCbins_min.size(), __PRETTY_FUNCTION__);
    this->LCbins_min[i]=value;
  }
  real_prec _LCbins_min(int i)
  {
    return this->LCbins_min[i];
  }

  //////////////////////////////////////////////////////////
  real_prec _LCbins_max(int i)
  {
    return this->LCbins_max[i];
  }
  void set_LCbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->LCbins_max.size(), __PRETTY_FUNCTION__);
    this->LCbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  real_prec _TAbins_min(int i)
  {
    return this->TAbins_min[i];
  }
  void set_TAbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->TAbins_max.size(), __PRETTY_FUNCTION__);
    this->TAbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  real_prec _TAbins_max(int i)
  {
    return this->TAbins_max[i];
  }
  void set_TAbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->TAbins_min.size(), __PRETTY_FUNCTION__);
    this->TAbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  real_prec _LOCALDMbins_min(int i)
  {
    return this->LOCALDMbins_min[i];
  }
  void set_LOCALDMbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->LOCALDMbins_min.size(), __PRETTY_FUNCTION__);
    this->LOCALDMbins_min[i]=value;
  }
  real_prec _LOCALDMbins_min_begin(){return this->LOCALDMbins_min.front();}

  //////////////////////////////////////////////////////////
  real_prec _LOCALDMbins_max(int i)
  {
    return this->LOCALDMbins_max[i];
  }
  void set_LOCALDMbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->LOCALDMbins_max.size(), __PRETTY_FUNCTION__);
    this->LOCALDMbins_max[i]=value;
  }
  real_prec _LOCALDMbins_max_last(){return this->LOCALDMbins_max.back();}
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  real_prec _PHbins_min(int i)
  {
    return this->PHbins_min[i];
  }

  void set_PHbins_min(int i,real_prec value)
  {
    this->check_dim(i,this->PHbins_min.size(), __PRETTY_FUNCTION__);
    this->PHbins_min[i]=value;
  }
  //////////////////////////////////////////////////////////
  real_prec _PHbins_max(int i)
  {
    return this->PHbins_max[i];
  }
  void set_PHbins_max(int i,real_prec value)
  {
    this->check_dim(i,this->PHbins_max.size(), __PRETTY_FUNCTION__);
    this->PHbins_max[i]=value;
  }
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _files_kernel_references(int i){return this->files_kernel_references[i];}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_MultiLevels(){return this->Number_of_MultiLevels;}


  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_Nft_MultiLevels(int i){return this->list_Nft_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_Ntracers_MultiLevels(int i){return this->list_Ntracers_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  ULONG get_d_ml(int i){return this->Lbox/static_cast<real_prec>(this->list_Nft_MultiLevels[i]);}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec get_PropThreshold_MultiLevels(int i){return this->list_Props_Threshold_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec get_Props_Tolerance_MultiLevels(int i){return this->list_Props_Tolerance_MultiLevels[i];}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void set_Props_Tolerance_MultiLevels(int i, real_prec newTa){this->list_Props_Tolerance_MultiLevels[i]=newTa;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void set_Ntracers_MultiLevels(int i, ULONG new_Nt){this->list_Ntracers_MultiLevels[i]=new_Nt;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void set_PropThreshold_MultiLevels(int il, real_prec val){this->list_Props_Threshold_MultiLevels[il]=val;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  void show_params();
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _get_window_matrix(){return this->get_window_matrix;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_tracer_number_counts(){return this->Get_tracer_number_counts;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tidal_anisotropy_at_halo(){return this->Get_tidal_anisotropy_at_halo;}
  /**
   *  @brief
   */
  bool _Get_peak_height_at_halo(){return this->Get_peak_height_at_halo;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec _Scale_mach_number(){return this->Scale_mach_number;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_tracer_mass_field(){return this->Get_tracer_mass_field;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_tracer_vmax_field(){return this->Get_tracer_vmax_field;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_tracer_spin_field(){return this->Get_tracer_spin_field;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_local_mach_number(){return this->Get_tracer_local_mach_number;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_local_dach_number(){return this->Get_tracer_local_dach_number;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_cell_local_mach_number(){return this->Get_tracer_local_mach_number;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_bias(){return this->Get_tracer_bias;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_relative_bias(){return this->Get_tracer_relative_bias;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_quadratic_bias(){return this->Get_tracer_quadratic_bias;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_tracer_local_dm_density(){return this->Get_tracer_local_dm_density;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_PCA(){return this->Get_PCA;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  bool _Get_local_overdensity(){return this->Get_local_overdensity;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_prop_function_tracer(){return this->Get_prop_function_tracer;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool _Get_prop_function_tracer_cwt(){return this->Get_prop_function_tracer_cwt;}



  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // Derived paramerters
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_x(){return this->delta_x;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_nyquist_frequency(){return 0.5*this->Nft*(this->deltak_0);}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_y(){return this->delta_y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_z(){return this->delta_z;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_x_low(){return this->delta_x_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_y_low(){return this->delta_y_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_delta_z_low(){return this->delta_z_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_x(){return this->deltak_x;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_y(){return this->deltak_y;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_z(){return this->deltak_z;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_x_low(){return this->deltak_x_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_y_low(){return this->deltak_y_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_z_low(){return this->deltak_z_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_deltak_0(){return this->deltak_0;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_kmin(){return this->kmin;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_kmin_low(){return this->kmin_low;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_kmax(){return this->kmax;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_DeltaK_data(){return this->DeltaK_data;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_DeltaK_data_low(){return this->DeltaK_data_low;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  real_prec _d_DeltaK_window(){return this->DeltaK_window;}
  
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _d_Nnp_window(){return this->Nnp_window;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _d_Nnp_data(){return this->Nnp_data;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _d_Deltal(){return this->Deltal;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _d_Deltamu(){return this->Deltamu;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  real_prec _rmin_cf(){return this->rmin_cf;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  real_prec _rmax_cf(){return this->rmax_cf;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _Nbins_cf(){return this->Nbins_cf;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  string _mark(){return this->mark;}
  /////////////////////////////////////////////////////////
  /**
   * @brief  This directive acts on the generation of overdensity fields from particl distribution
   * @detail If enabled, reads lists of particles and gnerate high res fields. A ow-pass filter is then applied to get
   * the mesh with the fiducual mesh size without aliasing corrections.
   */
  bool _use_low_pass_filter(){return this->use_low_pass_filter;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  string _rbin_type(){return this->rbin_type;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  int _unitsim_plabel(){return this->Unitsim_plabel;}

  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  int _Nbins_hist(){return this->Nbins_hist;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  bool _set_bins_equal_number_tracers(){return this->set_bins_equal_number_tracers;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _set_bins_equal_number_tracers_main_property(){return this->set_bins_equal_number_tracers_main_property;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_bins_equal_number_tracers(){return this->Number_of_bins_equal_number_tracers;}
  void set_Number_of_bins_equal_number_tracers(int newN){this->Number_of_bins_equal_number_tracers=newN;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  int _Number_of_bins_equal_number_tracers_main_property(){return this->Number_of_bins_equal_number_tracers_main_property;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  ULONG _Number_of_GRF(){return this->Number_of_GRF;}
  void set_Number_of_GRF(ULONG value){this->Number_of_GRF=value;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  real_prec _Kmax_FA(){return this->Kmax_FA;}
  void set_Kmax_FA(ULONG value){this->Kmax_FA=value;}
  ///////////////////////////////////////////////////////// 
  /**
   *  @brief 
   *  @return
   */
  bool _Generate_FA(){return this->Generate_FA;}
  void set_Generate_FA(bool value){this->Generate_FA=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _kmax_tracer_bias(){return this->kmax_tracer_bias;}
  void set_kmax_tracer_bias(real_prec new_value){this->kmax_tracer_bias=new_value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _kmin_tracer_bias(){return this->kmin_tracer_bias;}
  void set_kmin_tracer_bias(real_prec new_value){this->kmin_tracer_bias=new_value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _kmin_tracer_qbias(){return this->kmin_tracer_qbias;}
  void set_kmin_tracer_qbias(real_prec new_value){this->kmin_tracer_qbias=new_value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _kmax_tracer_qbias(){return this->kmax_tracer_qbias;}
  void set_kmax_tracer_qbias(real_prec new_value){this->kmax_tracer_qbias=new_value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _xmin(){return this->xmin;}
  void set_xmin(real_prec value){this->xmin=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _ymin(){return this->ymin;}
  void set_ymin(real_prec value){this->ymin=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _zmin(){return this->zmin;}
  void set_zmin(real_prec value){this->zmin=value;}

  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _xmax(){return this->xmax;}
  void set_xmax(real_prec value){this->xmax=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _ymax(){return this->ymax;}
  void set_ymax(real_prec value){this->ymax=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _zmax(){return this->zmax;}
  void set_zmax(real_prec value){this->zmax=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Xoffset(){return this->Xoffset;}
  void set_Xoffset(real_prec value){this->Xoffset=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Yoffset(){return this->Yoffset;}
  void set_Yoffset(real_prec value){this->Yoffset=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Zoffset(){return this->Zoffset;}
  void set_Zoffset(real_prec value){this->Zoffset=value;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _mKmax(){return this->mKmax;}
  void set_mKmax(real_prec value){this->mKmax=value;}
  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _mKmin(){return this->mKmin;}
  void set_mKmin(real_prec value){this->mKmin=value;}


  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  string _par_file(){return this->par_file;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _use_real_and_redshift_space(){return this->use_real_and_redshift_space;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_marked_power_spectrum(){return this->Get_marked_power_spectrum;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_power_spectrum(){return this->Get_power_spectrum;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_cross_power_spectrum(){return this->Get_cross_power_spectrum;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_pearson_coefficient(){return this->Get_pearson_coefficient;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_spearman_coefficient(){return this->Get_spearman_coefficient;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_Mstellar_function(){return this->Get_Mstellar_function;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_Luminosity_function(){return this->Get_Luminosity_function;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_Color_function(){return this->Get_Color_function;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_Color_Mag_plane(){return this->Get_Color_Mag_plane;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  bool _Get_Random_Catalog(){return this->Get_Random_Catalog;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _Number_of_random_files(){return this->Number_of_random_files;}
  void set_Number_of_random_files(int Ni){this->Number_of_random_files=Ni;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _Nbins_color(){return this->Nbins_color;}
  void set_Nbins_color(int Nc){this->Nbins_color=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  int _Nbins_Mstellar() {return this->Nbins_Mstellar;}
  void set_Nbins_Mstellar(int Nc){this->Nbins_Mstellar=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Mstellar_min(){return this->Mstellar_min;}
  void set_Mstellar_min(real_prec Nc){this->Mstellar_min=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Mstellar_max(){return this->Mstellar_max;}
  void set_Mstellar_max(real_prec Nc){this->Mstellar_max=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Color_min(){return this->Color_min;}
  void set_Color_min(real_prec Nc){this->Color_min=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _Color_max(){return this->Color_max;}
  void set_Color_max(real_prec Nc){this->Color_max=Nc;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RA_max(){return this->RA_max;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _RA_min(){return this->RA_min;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _DEC_max(){return this->DEC_max;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  real_prec _DEC_min(){return this->DEC_min;}
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
  void check_dim(int i, int j,string fun)
  {
    if(i>j)
      {
        cout<<"In "<<fun<<":"<<endl;
        cout<<"Requested memmory slot at i ="<<i<<endl;
        cout<<"Size of array: "<<j<<endl;
        throw std::invalid_argument("Dimension of array smaller than label of slot requested");
      }
  }
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
    void _resizeSPINBULLOCKbins(int new_size )
    {
        this->SPINBULLOCKbins_max.resize(new_size,0);
        this->SPINBULLOCKbins_min.resize(new_size,0);
    }
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
    void _resizeRSbins(int new_size )
    {
        this->RSbins_max.resize(new_size,0);
        this->RSbins_min.resize(new_size,0);
    }
  /////////////////////////////////////////////////////////
  /**
   *  @brief
   *  @return
   */
    void _resizeCONCENTRATIONbins(int new_size )
    {
        this->CONCENTRATIONbins_max.resize(new_size,0);
        this->CONCENTRATIONbins_min.resize(new_size,0);
    }

  // -----------------------------------------------------------------------------------------------
  /**
   *  @name input/output
   */
  real_prec _kmax_integration () {return this->kmax_integration;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kmin_integration () {return this->kmin_integration;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _Get_SO_from_BN(){return this->Get_SO_from_BN;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _fixed_redshift() {return fixed_redshift;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _redshift_min() {return redshift_min;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _redshift_max() {return redshift_max;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _nbins_redshift() {return nbins_redshift;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _A_gas() {return A_gas;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _B_gas() {return B_gas;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _mstar() {return mstar;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kstar() {return kstar;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _GAL_BIAS() {return GAL_BIAS;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _Amc() {return Amc;}
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _alpha_s() {return this->alpha_s;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */

  real_prec _coef_concentration_amp() {return coef_concentration_amp;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _coef_concentration() {return coef_concentration;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _mass_function_fit() {return this->mass_function_fit;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _sigma_ln() {return this->sigma_ln;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _sigma_red() {return this->sigma_red;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _missing_flux() {return this->missing_flux;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _density_profile(){return this->density_profile;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _muno_hod() {return this->muno_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _ms_hod() {return this->ms_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _mmin_hod() {return this->mmin_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _scatter_hod() {return this->scatter_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _hod_model() {return this->hod_model;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _alpha_hod() {return this->alpha_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _alpha_A() {return this->alpha_A;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _alpha_B() {return this->alpha_B;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _alpha_C() {return this->alpha_C;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _mt_hod() {return this->mt_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _s_faint_hod() {return this->s_faint_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _s_bright_hod() {return this->s_bright_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _width_hod() {return this->width_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _Delta_SO() {return this->Delta_SO;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _Mstep_hod() {return this->Mstep_hod;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _galaxy_power_spectrum_halo_model_output_file(){return galaxy_power_spectrum_halo_model_output_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _galaxy_correlation_function_halo_model_output_file(){return galaxy_correlation_function_halo_model_output_file;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _M_min_effective() {return M_min_effective;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _M_max_effective() {return M_max_effective;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _M_min_mf() {return M_min_mf;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _M_max_mf() {return M_max_mf;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _scale_mf() {return scale_mf ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _npoints_mf() {return npoints_mf ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _mass_function_output_file() {return mass_function_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _halo_mass_bias_fit() {return halo_mass_bias_fit ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _halo_mass_bias_output_file() {return halo_mass_bias_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _effective_halo_mass_bias_output_file() {return effective_halo_mass_bias_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _effective_halo_mean_number_density_output_file() {return effective_halo_mean_number_density_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _compute_output_linear_power_spectrum() {return compute_output_linear_power_spectrum ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _compute_output_non_linear_power_spectrum() {return compute_output_non_linear_power_spectrum  ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _scale_ps() {return scale_ps ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kmin_ps() {return kmin_ps ;}
  void set_kmin_ps(real_prec nk) {this->kmin_ps=nk;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kmax_ps() {return kmax_ps ;}
  void set_kmax_ps(real_prec nk) {this->kmax_ps=nk;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _npoints_ps() {return npoints_ps ;}
  void set_npoints_ps(int nk) {this->npoints_ps=nk;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _linear_matter_ps_output_file() {return linear_matter_ps_output_file  ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _non_linear_matter_ps_halo_fit_output_file() {return non_linear_matter_ps_halo_fit_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _non_linear_matter_ps_pt_output_file() {return non_linear_matter_ps_pt_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _scale_cf() {return scale_cf ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _npoints_cf() {return npoints_cf;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _linear_matter_cf_output_file() {return linear_matter_cf_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _non_linear_matter_cf_halo_fit_output_file() {return non_linear_matter_cf_halo_fit_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _compute_density_profile() {return compute_density_profile ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _rmin_dp() {return rmin_dp ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _rmax_dp() {return rmax_dp ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _scale_dp_r() {return scale_dp_r ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _npoints_dp_r() {return npoints_dp_r ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _density_profile_r_output_file() {return density_profile_r_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kmin_dp() {return kmin_dp ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _kmax_dp() {return kmax_dp ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _scale_dp_k() {return scale_dp_k ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  int _npoints_dp_k() {return npoints_dp_k ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _density_profile_k_output_file(){return density_profile_k_output_file ;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _compute_output_linear_correlation_function(){return compute_output_linear_correlation_function;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _compute_output_non_linear_correlation_function(){return compute_output_non_linear_correlation_function;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  bool _use_file_power(){return this->use_file_power;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  string _file_power_th(){return this->file_power_th;}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec _area_survey(){return this->area_survey;}




};



#endif
