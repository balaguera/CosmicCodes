/**
 *  @file ParametersF.h
 *
 *  @brief The class Parameters 
 *
 *  @author Federico Marulli and Andres Balaguera Antolinez
 *  @author Optimization and parallelization by Luca Tornatore
 *  This file defines the interface of the class Parameters, used to
 *  set all the parameters adopted in the analysis
 */


#ifndef __PARAMETERSF__
#define __PARAMETERSF__

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

/**
 *  @enum binning
 *  @brief What kind of binning to use
 *
 *  The most general way of specifying a binning is to give the vector
 *  of the edges of the bins, this is what ARBITRARY option stands
 *  for. It is possible that we will never need this option, we put it
 *  in here just in case. The other two options are LINEAR and LOG.
 */

/**
 *  @class Parameters Parameters.h "Lib/Headers/Parameters.h"
 *
 *  @brief The class Parameters
 *
 *  This class is used to handle objects of type <EM> Parameters </EM>
 */
class ParametersF 
{
  
 private :

  //////////////////////////////////////////////////////////
  /**
   *  @brief statisticts to compute 
   *
   *  possible options are
   *
   *  - Pk_fkp &rarr; power spectrum using the FKP estimator 
   *
   *  - Bk_fkp &rarr; bispectrum using the FKP estimator
   *
   *  - Pk_ys &rarr; moment decomposition of the 3D power spectrum
   *    using the Yamamoto-Blake estimator with the Hexadecapole as
   *    estimated by Scoccimarro. 
   *
   *  - Pk_yb &rarr; moment decomposition of the 3D power spectrum using the
   *    Yamamoto-Blake estimator with the Hexadecapole as estimated by
   *    Bianchi et al.  
   *
   *  - Pk_y_ds &rarr; moment decomposition of the 3d power spectrum
   *    using the Yamamoto-Blake estimator with direct sum over
   *    galaxies. (SLOW!)
   *
   *  - cl &rarr; angular power spectrum C<SUB>ll</SUB>
   *
   *  
   */  
  string statistics;
  //////////////////////////////////////////////////////////
  
  // ----- I/O variables -----
   
  /**
   *  @name input/output
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// input directory where the object and random catalogues are stored
  string dir_input;

  //////////////////////////////////////////////////////////
  /// output directory
  string dir_output;

  
  //////////////////////////////////////////////////////////
  /// name of the object catalogue
  string file_catalogue;
  

  unsigned long ngal_delta;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  **/
  string delta_grid_file;
  //////////////////////////////////////////////////////////
  /**
   * @brief name of the survey
  **/

  string delta_grid_file2;
  //////////////////////////////////////////////////////////
  /**
   * @brief name of the survey
  **/

  string delta_grid_file3;
  //////////////////////////////////////////////////////////
  /**
   * @brief name of the survey
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
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
  **/
  string input_type;

  
  //////////////////////////////////////////////////////////
  /**
   * @brief Type if input file, options
   * are "cat" meaning catalog, and "grid_delta" meaning
   * that the input is alread the delta in the grid
   **/
  string file_random;

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

  /// the column where the first object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord1_g;

  //////////////////////////////////////////////////////////
  /// the column where the second object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object coordinate (according to the system of coordinates of the catalog) is written
  int i_coord3_g;

  //////////////////////////////////////////////////////////
  /// the column where the first object weight is written
  int i_weight1_g;

  //////////////////////////////////////////////////////////
  /// the column where the second object weight is written
  int i_weight2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object weight is written
  int i_weight3_g;

  //////////////////////////////////////////////////////////
  /// the column where the fourth object weight is written
  int i_weight4_g;

  //////////////////////////////////////////////////////////
  /// the column where the object mean number density is written
  int i_mean_density_g;

  //////////////////////////////////////////////////////////
  /// the column where the first object property is written
  int i_property1_g;

  //////////////////////////////////////////////////////////
  /// the column where the second object property is written
  int i_property2_g;

  //////////////////////////////////////////////////////////
  /// the column where the third object property is written
  int i_property3_g;

  //////////////////////////////////////////////////////////
  /// the column where the fourth object property is written
  int i_property4_g;

  //////////////////////////////////////////////////////////
  /// the column where the fifth object property is written
  int i_property5_g;

  //////////////////////////////////////////////////////////
  /// the column where the sixth object property is written
  int i_property6_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the first weight; false &rarr; do not use the first weight
  bool use_weight1_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the second weight; false &rarr; do not use the second weight
  bool use_weight2_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the third weight; false &rarr; do not use the third weight
  bool use_weight3_g;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the fourth weight; false &rarr; do not use the fourth weight
  bool use_weight4_g;
  
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
  /// the column where the first random property is written
  int i_property1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random property is written
  int i_property2_r;

  //////////////////////////////////////////////////////////
  /// the column where the third random property is written
  int i_property3_r;

  //////////////////////////////////////////////////////////
  /// the column where the fourth random property is written
  int i_property4_r;

  //////////////////////////////////////////////////////////
  /// the column where the fifth random property is written
  int i_property5_r;

  //////////////////////////////////////////////////////////
  /// the column where the sixth random property is written
  int i_property6_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the first weight; false &rarr; do not use the first weight
  bool use_weight1_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the second weight; false &rarr; do not use the second weight
  bool use_weight2_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the third weight; false &rarr; do not use the third weight
  bool use_weight3_r;

  //////////////////////////////////////////////////////////
  /// true &rarr; use the fourth weight; false &rarr; do not use the fourth weight
  bool use_weight4_r;


  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

  /**
   *  @name FKP power spectrum
   */
  ///@{
  //////////////////////////////////////////////////////////  
  /// Number of grid cells /per dimension for the Discrete Fourier Trasnform
  int Nft;

  //////////////////////////////////////////////////////////
  ///  Lenght of the Fourier box in configuration space (in Mpc/h)
  double Lbox;

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
  /// Type of linear binning in Fourier space: 0.5 or 1.
  /// for 0.5, the bins are k_i= (i+0.5)*Delta. In this case the center of the
  /// first bin is at the mid point between the zero mode and the fundamental mode
  /// if Delta=delta.
  /// for 1, the bins are k_i= (i+1)*Delta. In this case the center of the
  /// the first bin is at the fundamental mode, if Delta=delta.
  double k_bin_step; 


  //////////////////////////////////////////////////////////
  /// Number of log-spaced bins in Fourier space
  int N_log_bins; 


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
  double Pest; 

  //////////////////////////////////////////////////////////
  ///  Use Poisson shot-noise correction (yes/no)
  bool SN_correction; 
  
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
  /// Has the sample a constant depth?
  bool constant_depth; 
  
  //////////////////////////////////////////////////////////
  ///Number of redshift bins to measure dNdz
  int N_z_bins; 
  
  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample. Used when dNdz ahs to be measured 
  double redshift_min_sample; 

  //////////////////////////////////////////////////////////
  /// Maximum redshift of the sample. Used when dNdz ahs to be measured 
  double redshift_max_sample; 
  
  //////////////////////////////////////////////////////////
  ///  Number of dNdz bins to measure 
  int N_dndz_bins; 

  //////////////////////////////////////////////////////////
  /// Number of redshift bins withoin which the measuerd dNdz will be smoothed
  int new_N_dndz_bins; 

  //////////////////////////////////////////////////////////
  /// Area of the survey.
  double area_survey; 

  //////////////////////////////////////////////////////////
  /// Resolution Healpix for pixelization. Used when no nbar is tabulated and dNdz is to be computed from a non-constant depth sample
  int Healpix_resolution; 

  //////////////////////////////////////////////////////////
  /// output file for the redshift distribution
  string file_dndz;   

  //////////////////////////////////////////////////////////
  /// Redefine a new line of sight
  bool new_los;
  
  //////////////////////////////////////////////////////////
  /// output file for the FKP power spectrum
  string file_power;
  
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
  double kmax_y_ds;
  
  ///@}
  
  

  /**
   *  @name bispectrum
   */
  ///@{
  
  //////////////////////////////////////////////////////////
  /// name of the output file used to store the bispectrum
  string file_bispectrum;

  //////////////////////////////////////////////////////////
  /// These parameters is used to define the shells in k-space
  bool use_fundamental_mode_as_kmin_bk;

  
  /////////////////////////////////////////////////////////////
  /// Minimum k-value for constructing k-bins
  double kmin_bk;

  /////////////////////////////////////////////////////////////
  /// Maximum k-value for constructing k-bins
  double kmax_bk;


  ///@}
  


  /**
   *  @name angular power spectrum
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// Maximum multupole for measurements
  int cL_max;

  //////////////////////////////////////////////////////////
  /// Minumum multupole for measurements
  int cL_min;

  //////////////////////////////////////////////////////////
  /// Galaxies in form of Catalog (cat) or healpix map (map)
  string type_of_input_file;


  //////////////////////////////////////////////////////////
  /// Type of zbins for Cl
  string type_of_zbins;


  //////////////////////////////////////////////////////////
  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin1_min;

  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin1_max;

  //////////////////////////////////////////////////////////
  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin2_min;

  //////////////////////////////////////////////////////////
  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin2_max;

  //////////////////////////////////////////////////////////
  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin3_min;

  //////////////////////////////////////////////////////////
  /// Redshift bins, value of the monomum redshift of zbin 1
  double zbin3_max;


  //////////////////////////////////////////////////////////
  /// Input file for the mask
  string input_file_mask;

  //////////////////////////////////////////////////////////
  /// Output file for the power spectrum of the mask
  string file_cl_window;

  //////////////////////////////////////////////////////////
  /// Output file for the power spectrum
  string file_power_cl;

  //////////////////////////////////////////////////////////
  /// Number of redshift bins to measure the cross power 
  int number_of_redshift_bins;

  //////////////////////////////////////////////////////////
  /// Number of bins in L for the output
  int number_of_Lbins;

  //////////////////////////////////////////////////////////
  /// the column where the pixel of the mask is written
  int i_mask_pixel;

  //////////////////////////////////////////////////////////
  /// the column where the RA of the pixel of the mask is written
  int i_mask_alpha;

  //////////////////////////////////////////////////////////
  /// the column where the DEC of the pixel of the mask is written
  int i_mask_delta;

  //////////////////////////////////////////////////////////
  /// the column where the FLAG of the mask is written
  int i_mask_flag;

  //////////////////////////////////////////////////////////
  /// log or linear binning
  string type_of_lbinning; 
  ///@}





  /**
   *  @name Fourier-Bessel decomposition
   */
  ///@{
  //////////////////////////////////////////////////////////
  /// Selected buondary conditions
  string Boundary_conditions;

  //////////////////////////////////////////////////////////
  ///  Compute k-spectra or read from file
  bool compute_new_kln;

  //////////////////////////////////////////////////////////
  /// External file  with allowed kln according to the selected boudanry conditions
  string file_kln;
  
  //////////////////////////////////////////////////////////
  /// File with the number of k-modes as a function of the multipole l
  string file_modes;
  
  //////////////////////////////////////////////////////////
  /// Minimum l-mode for the analysis
  int l_min_fb;

  //////////////////////////////////////////////////////////
  /// Maximum l-mode for the analysis. 
  int l_max_fb;

  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample used in the FB analysis
  double redshift_min_fb;

  //////////////////////////////////////////////////////////
  /// Maximum redshift of the sample used in the FB analysis
  double redshift_max_fb;

  //////////////////////////////////////////////////////////
  /// Output file for FB power spectrum 
  string file_power_fb;

  //////////////////////////////////////////////////////////
  /// Maximum wavenumber for the FB analysis
  double k_max_fb;
  ///@}



  /**
   *  @name Covariance of two-point statistics
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// Instructs the covariance code to use always the same data catalog
  bool use_one_data_for_covariance;

  //////////////////////////////////////////////////////////
  /// Instructs the covariance code to use always the same random catalog
  bool use_one_random_for_covariance;

  //////////////////////////////////////////////////////////
  /// number of catalogues (for covariances)
  int n_catalogues;

  ///@}

  /**
   *  @name cosmological parameters
   */
  ///@{

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>M</SUB>: the density of baryons, cold dark matter and massive neutrinos (in units of the critical density) at z=0
  double om_matter;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>rad</SUB>: the radiation density at z=0
  double om_radiation;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>b</SUB>: the baryon density at z=0
  double om_baryons;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>DE</SUB>: the dark energy density at z=0
  double om_vac;

  //////////////////////////////////////////////////////////
  /// &Omega;<SUB>k</SUB>: the density of curvature energy
  double om_k;

  //////////////////////////////////////////////////////////
  /// H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
  double Hubble;

  //////////////////////////////////////////////////////////
  /// \e h: the Hubble parameter, H<SUB>0</SUB>/100
  double hubble;

  //////////////////////////////////////////////////////////
  /// n<SUB>spec</SUB>: the primordial spectral index
  double spectral_index;

  //////////////////////////////////////////////////////////
  /// w<SUB>0</SUB>: the parameter of the dark energy equation of state
  double w_eos;

  //////////////////////////////////////////////////////////
  /// N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
  double N_eff;

  //////////////////////////////////////////////////////////
  /// &sigma;<SUB>8</SUB>: the power spectrum normalization
  double sigma8;

  //////////////////////////////////////////////////////////
  /// T<SUB>CMB</SUB>: the present day CMB temperature [K]
  double Tcmb;

  ///@}

  
 public:

  /**
   *  @brief default constructor
   *  @return object of class Parameters
   */
  ParametersF () {}
  

  /**
   *  @brief constructor 
   *  @param parameters_file parameter file
   *  @return object of class Parameters
   */
  ParametersF (string &);
  
  
  //////////////////////////////////////////////////////////
  /// Select MAS: NGP=nearest grid point. CIC=cloud in cell. TSC= triangular shape cloud
  string mass_assignment_scheme; 
  
  //////////////////////////////////////////////////////////
  /// Ratio between the shell-width and the fundamental mode
  int ndel_data; 
  
  
  //////////////////////////////////////////////////////////
  /// Correct for the MAS: yes/no
  bool MAS_correction; 
  
  //////////////////////////////////////////////////////////
  /**
   * @brief name of the survey
  **/
  string Name_survey;

  
  //////////////////////////////////////////////////////////
  /**
   *  @brief default destructor
   *  @return none
   */
  ~ParametersF () {}


  //////////////////////////////////////////////////////////
  /**
   *  @brief set the names of the output files
   *  @return none
   */
  void build_file_names ();

  

  //////////////////////////////////////////////////////////
  
  // function to get private variables
  /**
   *  @brief get the value of the private member statistics
   *  @return statistics
   */
  string _statistics () {return statistics;}
  //////////////////////////////////////////////////////////  

  /**
   *  @brief get the value of the private member dir_input
   *  @return dir_input
   */
  string _dir_input () {return dir_input;}
  //////////////////////////////////////////////////////////
  /**
   * @brief ask if cross is desidered
  **/
  bool measure_cross;

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _dir_output () {return dir_output;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member dir_output
   *  @return dir_output
   */
  string _input_type () {return input_type;}

  unsigned long _ngal_delta() {return ngal_delta;}
  string _delta_grid_file () {return delta_grid_file;}
  string _delta_grid_file2 () {return delta_grid_file2;}
  string _delta_grid_file3 () {return delta_grid_file3;}
  string _delta_grid_file4 () {return delta_grid_file4;}

  bool _measure_cross () {return measure_cross;}

  int _measure_cross_from_1 () {return measure_cross_from_1;}
  int _measure_cross_from_2 () {return measure_cross_from_2;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief set the value of the private member dir_input
   *  @param _dir_input input directory where the object and random
   *  catalogues are stored
   *  @return none
   */
  void set_dir_input (string _dir_input) {dir_input = _dir_input;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief set the value of the private member dir_output
   *  @param _dir_output output directory
   *  @return none
   */
  void set_dir_output (string _dir_output) {dir_output = _dir_output;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member file_catalogue
   *  @return file_catalogue
   */
  string _file_catalogue () {return file_catalogue;}
    
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
   *  @brief get the value of the private member i_coord1_g
   *  @return i_coord1_g
   */
  int _i_coord1_g () {return i_coord1_g;}

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_coord2_g
   *  @return i_coord2_g
   */
  int _i_coord2_g () {return i_coord2_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_g
   *  @return i_coord3_g
   */
  int _i_coord3_g () {return i_coord3_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_g
   *  @return i_weight1_g
   */
  int _i_weight1_g () {return i_weight1_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_g
   *  @return i_weight2_g
   */
  int _i_weight2_g () {return i_weight2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_g
   *  @return i_weight3_g
   */
  int _i_weight3_g () {return i_weight3_g;}

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
  bool _use_weight1_g () {return use_weight1_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight2_g
   *  @return use_weight2_g
   */
  bool _use_weight2_g () {return use_weight2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight3_g
   *  @return use_weight3_g
   */
  bool _use_weight3_g () {return use_weight3_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight4_g
   *  @return use_weight4_g
   */
  bool _use_weight4_g () {return use_weight4_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member Name_survey
   *  @return Name_survey
   */
  string _Name_survey () {return Name_survey;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_mean_density_g
   *  @return i_mean_density_g
   */
  int _i_mean_density_g () {return i_mean_density_g;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member angles_units_g
   *  @return angles_units_g
   */
  string _angles_units_g () {return angles_units_g;};

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog () {return use_random_catalog;};


  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_random_catalog
   *  @return use_random_catalog
   */
  bool _use_random_catalog_cl () {return use_random_catalog_cl;};
  

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member new_los
   *  @return new_los
   */
  bool _new_los () {return new_los;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member sys_of_coord_r
   *  @return sys_of_coord_r
   */
  int _sys_of_coord_r () {return sys_of_coord_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord1_r
   *  @return i_coord1_r
   */
  int _i_coord1_r () {return i_coord1_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord2_r
   *  @return i_coord2_r
   */
  int _i_coord2_r () {return i_coord2_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_coord3_r
   *  @return i_coord3_r
   */
  int _i_coord3_r () {return i_coord3_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight1_r
   *  @return i_weight1_r
   */
  int _i_weight1_r () {return i_weight1_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight2_r
   *  @return i_weight2_r
   */
  int _i_weight2_r () {return i_weight2_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight3_r
   *  @return i_weight3_r
   */
  int _i_weight3_r () {return i_weight3_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_weight4_r
   *  @return i_weight4_r
   */
  int _i_weight4_r () {return i_weight4_r;};

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property1_r
   *  @return i_property1_r
   */
  int _i_property1_r () {return i_property1_r;};

  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_property2_r
   *  @return i_property2_r
   */
  int _i_property2_r () {return i_property2_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property3_r
   *  @return i_property3_r
   */
  int _i_property3_r () {return i_property3_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property4_r
   *  @return i_property4_r
   */
  int _i_property4_r () {return i_property4_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property5_r
   *  @return i_property5_r
   */
  int _i_property5_r () {return i_property5_r;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property6_r
   *  @return i_property6_r
   */
  int _i_property6_r () {return i_property6_r;};

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property1_g
   *  @return i_property1_g
   */
  int _i_property1_g () {return i_property1_g;};
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property2_g
   *  @return i_property2_g
   */
  int _i_property2_g () {return i_property2_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property3_g
   *  @return i_property3_g
   */
  int _i_property3_g () {return i_property3_g;}
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property4_g
   *  @return i_property4_g
   */
  int _i_property4_g () {return i_property4_g;}

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member i_property5_g
   *  @return i_property5_g
   */
  int _i_property5_g () {return i_property5_g;}
  //////////////////////////////////////////////////////////  
  /**
   *  @brief get the value of the private member i_property6_g
   *  @return i_property6_g
   */
  int _i_property6_g () {return i_property6_g;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member use_weight1_r
   *  @return use_weight1_r
   */
  bool _use_weight1_r () {return use_weight1_r;}
 
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
   *  @brief get the value of the private member use_one_data_for_covariance
   *  @return use_one_data_for_covariance
   */
  bool _use_one_data_for_covariance () {return use_one_data_for_covariance;};

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_one_random_for_covariance
   *  @return use_one_random_for_covariance
   */
  bool _use_one_random_for_covariance () {return use_one_random_for_covariance;};

  /**
   *  @brief get the value of the private member Nft
   *  @return Nft
   */
  int _Nft () {return Nft;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Lbox
   *  @return Lbox
   */
  double _Lbox () {return Lbox;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member mass_assignment_scheme
   *  @return mass_assignment_scheme
   */
  string _mass_assignment_scheme () {return mass_assignment_scheme;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member type_of_binning
   *  @return type_of_binning
   */
  string _type_of_binning () {return type_of_binning;}


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member k_bin_step
   *  @return k_bin_step
   */
  double _k_bin_step () {return k_bin_step;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_log_bins
   *  @return N_log_bins
   */
  int _N_log_bins () {return N_log_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member ndel_data 
   *  @return ndel_data 
   */
  int _ndel_data () {return ndel_data;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member ndel_window
   *  @return ndel_window
   */
  int _ndel_window () {return ndel_window;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_mu_bins
   *  @return N_mu_bins
   */
  int _N_mu_bins () {return N_mu_bins;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member MAS_correction 
   *  @return MAS_correction 
   */
  bool _MAS_correction () {return MAS_correction;}
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_weight
   *  @return FKP_weight
   */
  bool _FKP_weight () {return FKP_weight;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member SN_correction
   *  @return SN_correction
   */
  bool _SN_correction () {return SN_correction;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_error_bars
   *  @return FKP_error_bars
   */
  bool _FKP_error_bars () {return FKP_error_bars;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member FKP_error_bars_exact
   *  @return FKP_error_bars_exact
   */
  bool _FKP_error_bars_exact () {return FKP_error_bars_exact;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Pest
   *  @return Pest
   */
  double _Pest () {return Pest;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member nbar_tabulated
   *  @return nbar_tabulated
   */
  bool _nbar_tabulated () {return nbar_tabulated;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member constant_depth
   *  @return constant_depth 
   */
  bool _constant_depth () {return constant_depth;};
  ////////////////////////////////////////////////////////// 
  
  /**
   *  @brief get the value of the private member N_z_bins
   *  @return N_z_bins
   */
  int _N_z_bins () {return N_z_bins;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member redshift_min_sample
   *  @return redshift_min_sample
   */
  double _redshift_min_sample () {return redshift_min_sample;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member redshift_max_sample
   *  @return redshift_max_sample
   */
  double _redshift_max_sample () {return redshift_max_sample;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N_dndz_bins
   *  @return N_dndz_bins
   */
  int _N_dndz_bins () {return N_dndz_bins;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member new_N_dndz_bins
   *  @return new_N_dndz_bins 
   */
  int _new_N_dndz_bins () {return new_N_dndz_bins;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member area_survey
   *  @return area_survey 
   */
  double _area_survey () {return area_survey;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Healpix_resolution
   *  @return Healpix_resolution 
   */
  int _Healpix_resolution () {return Healpix_resolution;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_dndz
   *  @return file_dndz
   */
  string _file_dndz () {return file_dndz;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power
   *  @return file_power 
   */
  string _file_power () {return file_power;};
  
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power_log
   *  @return file_power_log 
   */
  string _file_power_log () {return file_power_log;};
  ////////////////////////////////////////////////////////// 

  /**
   *  @brief get the value of the private member file_window
   *  @return file_window 
   */
  string _file_window () {return file_window;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power2d
   *  @return file_power2d 
   */
  string _file_power2d () {return file_power2d;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power2d_mk
   *  @return file_power2d_mk
   */
  string _file_power2d_mk () {return file_power2d_mk;};

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmin_bk
   */
  double _kmin_bk () {return kmin_bk;};	
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmin_bk;
   *  @return kmax_bk
   */
  double _kmax_bk () {return kmax_bk;};	


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member use_fundamental_mode_as_kmin_bk
   *  @return Nshells_bk
   */
  bool _use_fundamental_mode_as_kmin_bk () {return use_fundamental_mode_as_kmin_bk;};	



  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_bispectrum
   *  @return file_bispectrum
   */
  string _file_bispectrum () {return file_bispectrum;};	



  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member kmax_y_ds
   *  @return kmax_y_ds
   */
  double _kmax_y_ds () {return kmax_y_ds;};




 
  ////////////////////////////////////////////////////////// 
  
  /**
   *  @brief get the value of the the private member cL_max
   *  @return cL_max
   */
  int _cL_max() {return cL_max;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member cL_min
   *  @return cL_min
   */
  int _cL_min() {return cL_min;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_input_file
   *  @return type_of_input_file
   */
  string _type_of_input_file(){return type_of_input_file;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member input_file_mask
   *  @return input_file_mask
   */
  string _input_file_mask(){return input_file_mask;};
  ////////////////////////////////////////////////////////// 

  /**
   *  @brief get the value of the the private member file_cl_window
   *  @return file_cl_window
   */
  string _file_cl_window(){return file_cl_window;};
  ////////////////////////////////////////////////////////// 

  /**
   *  @brief get the value of the the private memberfile_power_cl
   *  @return file_power_cl
   */
  string _file_power_cl(){return file_power_cl;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member number_of_redshift_bins
   *  @return number_of_redshift_bins
   */
  int _number_of_redshift_bins() {return number_of_redshift_bins;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  string _type_of_lbinning() {return type_of_lbinning;};

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  string _type_of_zbins() {return type_of_zbins;};

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin1_min() {return zbin1_min;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin1_max() {return zbin1_max;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin2_min() {return zbin2_min;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin2_max() {return zbin2_max;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin3_min() {return zbin3_min;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member type_of_lbinning
   *  @return type_of_lbinning 
   */
  double _zbin3_max() {return zbin3_max;};


  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member number_of_Lbins
   *  @return type_of_lbinning 
   */
  int _number_of_Lbins() {return number_of_Lbins;};
  
  /** ////////////////////////////////////////////////////////// 
   *  @brief get the value of the the private member i_mask_pixel
   *  @return i_mask_pixel
   */
  int _i_mask_pixel(){return i_mask_pixel;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_alpha
   *  @return  i_mask_alpha
   */
  int _i_mask_alpha(){return i_mask_delta;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_delta
   *  @return  i_mask_delta
   */
  int _i_mask_delta(){return i_mask_delta;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the the private member i_mask_flag
   *  @return  i_mask_flag
   */
  int _i_mask_flag(){return i_mask_flag;};

  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member redshift_min_fb
   *  @return redshift_min_fb
   */
  double _redshift_min_fb () {return redshift_min_fb;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member redshift_max_fb
   *  @return redshift_max_fb
   */
  double _redshift_max_fb () {return redshift_max_fb;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member k_max_fb
   *  @return k_max_fb
   */
  double _k_max_fb () {return k_max_fb;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member l_max_fb
   *  @return l_max_fb
   */
  int _l_max_fb () {return l_max_fb;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member l_min_fb
   *  @return l_min_fb
   */
  int _l_min_fb () {return l_min_fb;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member Boundary_conditions
   *  @return Boundary_conditions
   */
  string _Boundary_conditions () {return Boundary_conditions;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private membercompute_new_kln
   *  @return compute_new_kln
   */
  bool _compute_new_kln () {return compute_new_kln;};
  
  /** ////////////////////////////////////////////////////////// 
   *  @brief get the value of the private member file_kln
   *  @return file_kln
   */
  string _file_kln () {return file_kln;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_modes
   *  @return file_modes
   */
  string _file_modes () {return file_modes;};
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member file_power_fb
   *  @return file_power_fb
   */
  string _file_power_fb () {return file_power_fb;}; 
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>M</SUB>
   *  @return &Omega;<SUB>M</SUB>
   */
  double _om_matter () { return om_matter; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>rad</SUB>: the radiation density at z=0
   *  @return &Omega;<SUB>rad</SUB>
   */
  double _om_radiation () { return om_radiation; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>b</SUB>: the baryon density at z=0
   *  @return &Omega;<SUB>b</SUB>
   */
  double _om_baryons () { return om_baryons; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>DE</SUB>: the dark energy density at z=0
   *  @return &Omega;<SUB>DE</SUB>
   */
  double _om_vac () { return om_vac; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &Omega;<SUB>k</SUB>: the density of curvature energy
   *  @return &Omega;<SUB>k</SUB>
   */
  double _om_k () { return om_k; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member H<SUB>0</SUB>: the Hubble constant at z=0 [km/sec/Mpc] 
   *  @return H<SUB>0</SUB>
   */
  double _Hubble () { return Hubble; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member \e h: the Hubble parameter, H<SUB>0</SUB>/100
   *  @return \e h
   */
  double _hubble () { return hubble; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member n<SUB>spec</SUB>: the primordial spectral index
   *  @return n<SUB>spec</SUB>
   */
  double _spectral_index () { return spectral_index; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member w<SUB>0</SUB>: the parameter of the dark energy equation of state
   *  @return w<SUB>0</SUB>
   */
  double _w_eos () { return w_eos; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member N<SUB>eff</SUB>: the effective number (for QED + non-instantaneous decoupling)
   *  @return N<SUB>eff</SUB> 
   */
  double _N_eff () { return N_eff; };
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member &sigma;<SUB>8</SUB>: the power spectrum normalization
   *  @return &sigma;<SUB>8</SUB>
   */
  double _sigma8 () { return sigma8; };

  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member T<SUB>CMB</SUB>: the present day CMB temperature [K]
   *  @return T<SUB>CMB</SUB>
   */
  double _Tcmb () { return Tcmb; };



};

#endif


