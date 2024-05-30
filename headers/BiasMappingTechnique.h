//////////////////////////////////////////////////////////
/**
 * @class<BiasMT>
 * @brief Header file for the class BiasMT::
 * @file BiasMappingTechnique.h
 * @title Bias Mapping Technique for mock catalogs
 * @author Andres Balaguera-Antolínez.
 * @details Based on the BiasAssignment Method developed by A Balaguera and F.S.Kitaura
 * @version   1.0
 * @date      2020
 * @details: This is an example of a main function to call BiasMT. A file called cosmicatlas.cpp
 */
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#ifndef _BiasMT__
#define _BiasMT__
# include <boost/tuple/tuple.hpp>
# include<iostream>
# include<fstream>
# include<vector>
# include<algorithm>
# include<math.h>
# include<sstream>
# include<iomanip>
# include<string>
//////////////////////////////////////////////////////////
# include "NumericalMethods.h"
# include "ScreenOutput.h"
# include "Params.h"
# include "PowerSpectrumF.h"
# include "PowerSpectrumTH.h"
# include "McmcFunctions.h"
# include "Cwclass.h"
# include "Catalog.h"
#ifdef _USE_LPT_
# include "LPT.h"
#endif
#include "def_bmt.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<prop_collection>
 * @details Strucure to allocate set of masses, coordinates and their ID (in the mesh) for tracers within a given bin of the multidimendional Theta (DM) properties
 * @details Used in the assignment of halo properties.
 * @details Example: Private variable in class::Bam
 * @code
 vector<prop_collection> dm_properties_bins;
 *@endcode
 */
struct prop_collection{
  /**
   *@brief  Masses (or any other property) in that particular theta_bin
   */
  vector<real_prec> tracer_properties;
#ifdef _USE_HYBRID_ASSIGNMENT_
  /**
   * @brief Statistics (mean and variance) of the distribition of tracer properites on each theta bin
   * @details The size of this container depends on the number of moments we need, usually two.
   */
  vector<real_prec> stats_tracer_properties;

  /**
   * @brief Statistics (mean and variance) of the distribition of tracer properites on each theta bin
   * @details The size of this container depends on the number of moments we need, usually two.
   */
#endif
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
  /**
   * @brief  Masses (or any other property) in that particular theta_bin
   **/
  vector<real_prec> tracer_properties_secondary;
#endif
  /**
   * @brief  Flag to pinpoint properties used
   **/
  vector<bool> used_property;
#ifdef _USE_HYBRID_ASSIGNMENT_NEW_
  /**
   * @brief  Flag to pinpoint properties used
   **/
  vector<vector<ULONG>> GalID_in_Vmax_bin;
#ifdef _CHECK_HYBRID_V2_
  /**
   * @brief  Flag to pinpoint properties used
   **/
  vector<real_prec> Pjj_upward;
  vector<real_prec> Pjj_downward;
  vector<real_prec> Pjj;
  vector<real_prec> Pfj_upward;
  vector<real_prec> Pfj_downward;
#endif
#endif
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
  /**
   * @brief  Flag to pinpoint properties used
   **/
  vector<bool> used_property_secondary;
#endif
#ifdef _CORRECT_EXCLUSION_
  /**
   *@brief  x-coordinates of the objects in that theta-bin */
  vector<real_prec> x_coord_bin_properties;
  /**
   *@brief  y-coordinates of the objects in that theta-bin */
  vector<real_prec> y_coord_bin_properties;
  /**
   *@brief  z-coordinates of the objects in that theta-bin */
  vector<real_prec> z_coord_bin_properties;
  /**
   *@brief  Grid ID of the objects in that theta-bin */
  vector<ULONG> GridID_bin_properties;
#endif
  /**
   *@brief  Object ID of the objects in that theta-bin
   **/
  vector<ULONG> GalID_bin_properties;
#if defined (_USE_TWO_REFS_MOCKS_ASSIGNMENT_) || defined (_NEW_APPROACH_ASS_) || defined (_ASSIGN_MASS_LINKED_TO_V_FROM_REF_)
  /**
   *@brief Index to identify the reference from which this tracer comes */
  vector<int> index_reference;
#endif
};

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class BiasMT{
private:
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
#ifdef _USE_LPT_
  /**
   * @private
   * @brief Object of type LPT, to be used in BiasMT
   */
  LPT lpt;
#endif
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Cwclass: for Cosmic web analysis.
   */
  Cwclass cwclass;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Cwclass, for Cosmic web analysis.
   */
  Cwclass cwclass_ref;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Params
   */
  Params params;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Params
   */
  Params params_original;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Cosmology
   */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief OBject of type McmcmFunctions. Used in _BIAS_MODE_
   */
  McmcFunctions mcmc;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Catalog used to allocate info of BiasMT mocks
   */
  Catalog tracer;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Catalog used to allocate info of BiasMT mocks
   */
  Catalog tracer_aux;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Catalog: used to allocate info of reeference catalog
   */
  Catalog tracer_ref;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Structure of type s_CosmoInfo
   */
  s_CosmoInfo s_cosmo_info;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Boolean type variable to specify whether we are in BiasMT mode (mock or bias) or something else (e.g, running BiasMT with the -c option)
   */
  bool BiasMT_mode;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Boolean used to solve the behavior of under-biased tracers.
   */
  bool aux_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Boolean used to solve the behavior of under-biased tracers
   */
  bool used_once;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type ofstream
   */
  ofstream sal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type ofstream
   */
  ofstream output_res;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Mean gas density
   */
  real_prec Mean_rho_gas;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Mean density of DM
   */
  real_prec Mean_density_X;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Mean density of DM tracer
   */
  real_prec Mean_density_Y;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Mean mass-weighted density of dark mattr tracer
   */
  real_prec Mean_density_Y_MASS;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Mean f-weighted trcer density. f  is the fraction of satelites
   */
  real_prec Mean_density_Y_SAT_FRAC;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Auxiliary variable
   */
  real_prec mean_aux;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Poisson shot noise of the halo field generated at each iteration
   */
  real_prec shot_noise_new;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Poisson shot noise of the reference halo field
   */
  real_prec shot_noise_ref;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of properties used to characterize the halo bias
   */
  int N_bias_properties;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  real_prec minimum_multiscale_property;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief function of class BiasMT. Load the class BiasMT with parameters
   */
  void warnings();
  //////////////////////////////////////////////////////////
  /**
   * @brief Load cosmological information at the input redshift
   */
  void get_cosmo();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Computes the DM density field at each iteration
   * by convolviong the initial DM field with the Kernel
   */
  void get_BiasMT_DM();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Computes the DM density field at eacfh iteration
   * by convolviong the initial DM field with the Kernel
   */
  void get_new_DM_field();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Computes the DM density field at eacfh iteration
   * by convolviong the initial DM field with the Kernel
   */
  void get_minimum_multiscale_property();

  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @private
   * @brief Compute the bias as the probability of one cell having a given number of
   * halos conditional to the properties of the DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * string counts gets the bias B(N|{theta}_dm), used to assign number counts in cells,
   * while string mass computs
   * the bias B(M|{theta}_dm), used to assign halo mass to a cell
   * @returns This function loads the class member container BiasMT::CWT_X_Y_hist

   */
  void get_BIAS(string);
#endif
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE

  /**
   * @brief Compute the abundance with respect to a X property from the refernece catalog, conditional to the properties of teh DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * This is an alternative to the mass-assignmetn approach given by the the calculation of B(M|{theta}), where M 
   * is the total mass of tracers in cell. Here instead we compute n(Mh|{\theta}) to assing halo masses
   * (i., e, on an object by object basis). This approach  works better. DOES NOT USE THE INFORMATION FROM V-CLASS
   * @implements Dark matter prperties and Vmax
   * @returns This function loads the class member container BiasMT::CWT_X_Y_hist
   */
  void get_scaling_relations_primary_property();
  //////////////////////////////////////////////////////////
  /**
   * @brief Same goal as get_scaling_relations_primary_property, used when two or more references are used to learn the distribution of properties from
   */
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
  void get_scaling_relations_primary_property_two(int );
#endif
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Compute the abundance with respect to a X property from the reference catalog, conditional to the properties of the DM density field.
   * @details This is accomplished by constructing an histogram in the different variables.
   * @details The histogram has the form P(Y|X1,X2,X3,X4,X5)
   * <table>
  <caption id="multi_row">Memmory slots</caption>
  <tr><td rowspan="1">Mvir <td>Vmax <td> Delta_dm<td> CWT <td> Mach5<td> Delta5
  <tr><th>Porperty to assign Y <th>SLOT X1 <th>SLOT X2<th>SLOT X3<th>SLOT X4<th>SLOT X5
  <tr><td rowspan="1">Cvir <td>Mvir <td> Delta_dm<td> CWT <td> Mach5<td> Delta5
  <tr><td rowspan="1">Spin <td>Mvir <td> Cvir<td> CWT <td> Mach5<td> Delta5 
  </table>
   * @returns This function loads the class member container BiasMT::ABUNDANCE_normalized to be used in the asignment.
   */
  void get_scaling_relations_secondary_property(string);
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief 
   */
  void extrapolate_kernel(vector<real_prec>&, vector<real_prec>&);
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @private
   * @brief Generate the tracer mock density field using BiasMT
   * @details The TR number density field is generating by sampling the conditional 
   * probability distribution obtained by BiasMT onto the Target DM density field.
   * @param new_dm_field; bool true/false to ask whether a new target DM fiueld is used to creat mock
   * @details In general, this TDMF is another realization of the same IC used to calibrate the kernel
   * @return The output is a TR number density field with its power spectrum
   * used in the iterative procedure of BiasMT.
   */
  void get_mock_grid(string property);
#ifdef _USE_TWO_REFS_MOCKS_
  /**
   * @private
   * @brief
   */
  void get_mock_grid_two(string property);
#endif
#endif
#ifdef _GET_BiasMT_CAT_
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief This method is used if _USE_NEW_MASS_ASSIGNMENT_ is defined
   * @details This property assignment explicitely uses the values of the properties taken from the reference and assign them to the mocks provided
   * the underlying properties of the DM field.
   * For Initial assignment = true, the initial assignment is done over a main property (such as Vmax which seems to be that correlating the most with DM)
   * Initial_assignment = false, the mass is to be assigned based on the assignment of Vmax done in a previous call of this function.
   * @note This method needs to be in agreement with what the method BiasMT::get_scaling_relations_primary_property() has done
   */
  void assign_tracer_property(bool initial_assignment, string h_prop);
  //////////////////////////////////////////////////////////
#endif
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief This function is used if _USE_NEW_MASS_ASSIGNMENT_ is undefined
   */
  void assign_tracer_mass_new();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Generate random positions of tracers
   * @details The function reads the nhumber of tracers per cell and assign random coordinates
   * @return Writes to a file three columns with x, y, z
   */
  void sample_mock();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Read density fields using as inputs
   * @param file_dm: path to the DM density field, (file_dm)
   * @param files_tr: path to the DM  tracer (files_tr) 
   * @param files_dm_ref_pdf: path to reference dm density field. Binary files
   * After the call of this functions, the class members
   * delta_Y, delta_X and delta_X_ini contain the information of the TR 
   * and DM DENSITY field
   */

#ifdef _USE_VELOCITIES_
  void read_BiasMT_files(string file_dm, string files_tr, string files_tr_mw, string file_y_mass, string file_y_sf, string files_dm_ref_pdf, string file_Vx, string file_Vy, string file_Vz);
#else
  void read_BiasMT_files(string file_dm, string files_tr, string files_tr_mw, string file_Y_mass, string file_Y_sat_frac, string files_dm_ref_pdf);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Analyze the density fields read by read_BiasMT_files()
   * @details performing some basic statistics and power spectrum. If defiened as pre-processor directive (RUN_TEST), 
   * a TEST is performed here.
   * @details This function redefines the values of delta_X_min etc in case it is requeitsd in parameter file. 
   * Else, input values in parameter file are used.
   * @details The fields BiasMT::delta_X, BiasMT::delta_X_ini and BiasMT::delta_Y are, if requiested, converted to OVERDENSITIES.
   */
  void analyze_input();
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Compute, if desired, or define from input file,
   * the min and max values of variables X and Y
   */
  void get_min_max_X_Y(); // to be derpecated
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Compute, if desired, or define from input file,
   * the min and max values of variables X and Y
   */
  void get_new_min_max_properties();

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Structure containing the minimum values of the different DM properties used in the measurement of BIAS
   */
  s_minimums s_mins;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Structure containing the maximum values of the different DM properties used in the measurement of BIAS
   */
  s_maximums s_maxs;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Structure containing the bin sizes values of the different DM properties used in the measurement of BIAS
   */
  s_Deltas s_deltas;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Construction of a halo number counts field from a bias model
   * @details Given the DM overdensity, and follopwing the model pre-defined in bias() (def.h), this computs
   * the halo number counts in a mesh (defined by private variables)
   * @return aux: Halo number counts in a mesh
   */
  void dark_matter_to_halos_analytical(const vector<real_prec>& in, vector<real_prec>&aux);
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Size of bin in Y (often tracer) used to construct the BIAS
   */
  real_prec DELTAY;

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Size of bin in Y (often tracer) used to construct the BIAS
   */
  real_prec DELTAY_MW;

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Size of bin in Y (often DM) used to construct the BIAS
   */
  real_prec DELTAX;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Container for the minimum halo-masses in bins of Mh
   */
  vector<real_prec> MBmin;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  bool read_BiasMT_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  bool dm_already_done;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  vector<prop_collection> dm_properties_bins;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<prop_collection> dm_properties_for_randoms_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties for mocks
   */

  vector<prop_collection> dm_properties_bins_mock;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the maximum halo-masses in bins of Mh
   */
  vector<real_prec> MBmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the X (generally DM) density field
   * During running time it is converted to overdensity
   */
  vector<real_prec>delta_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */

  vector<real_prec>Displacement_inicial;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  vector<real_prec>delta_dm_aux; // this is used for mass assignmet
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  vector<real_prec>delta_dm_aux_mem; // this is used for mass assignmet

  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the x-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Velx_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the y-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Vely_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the z-component of the velocity field of X (generally DM) density field
   */
  vector<real_prec>Velz_X;
  //////////////////////////////////////////////////////////
  /*
   * @brief  Container for the divergence of the DM velocity field
   */
  vector<real_prec>Divergence_VelField;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_shear_velfield();
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the X (generally DM) density field of reference
   * During running time it is converted to overdensity
   */
  vector<real_prec>delta_X_REF_PDF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector of structures. Contains the values of the masses that belongs to a ginven bin in the DM properties
   */
  //  real_prec kmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the X (generally DM) density field. This is used for the calibration procedure, 
   * in which the DM density field is convolved with a kernel to match the TR P(k). In that procedure, the 
   * container BiasMT::delta_X is modified, so we always use the original BiasMT::delta_X_ini to do the convolution.
   */
  vector<real_prec>delta_X_ini;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) density field
   * @brief  During running time it is converted to overdensity, if requested in parameter file
   */
  
  vector<real_prec>delta_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) MASS weighted density field
   */
  vector<real_prec>delta_Y_HR;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the Y (generally TR) MASS weighted density field
   */

  vector<real_prec>delta_Y_MASS;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>delta_Y_SAT_FRAC;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the number counts or other prop used inside the get_mock function
   */
  vector<real_prec>delta_Y_new;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE TR power spectrum
   */
  vector<real_prec>Power_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>weight_mh;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE mass weighted power spectrum
   */
  vector<real_prec>Power_REF_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mock power spectrum
   */
  vector<real_prec>Power_NEW;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mock mass-weighted power spectrum
   */
  vector<real_prec>Power_NEW_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the REFERENCE DM power spectrum
   */
  vector<real_prec>Power_DM_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BiasMT DM power spectrum
   */
  vector<real_prec>Power_DM_NEW;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BiasMT DM power spectrum convolved with the kernel
   */
  vector<real_prec>Power_DM_KONV;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the wave vectors
   */
  vector<real_prec>kvec;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of grid cells 
   * @details computed inside BiasMT::BiasMT(params _par) as (BiasMT::Nft)*(BiasMT::Nft)*(BiasMT::Nft)
   */
  ULONG NGRID_low_mass_assignment;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of Grid cells N*N*(N/2+1). Computed inside BiasMT (constructor)
   */
  ULONG NT;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  ULONG NTT;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  ULONG NTT_low;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of x=, delta_x or log10(num_in_log+delta_x). Computed inside BiasMT
   */
  real_prec Xmin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of x=, delta_x or log10(num_in_log+delta_x). Computed inside BiasMT
   */
  real_prec Xmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief Minimum value of Y=, delta_Y or log10(num_in_log+delta_Y). Computed inside BiasMT
   */
  real_prec Ymin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum value of y=, delta_t or log10(num_in_log+delta_y). Computed inside BiasMT
   */
  real_prec Ymax;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BiasMT
   */
  ULONG new_nbins_x;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BiasMT
   */
  ULONG new_nbins_y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of bins in X for BIAS. Computed inside BiasMT
   */
  ULONG new_nbins_y_MW;

  //////////////////////////////////////////////////////////
  /**
   * @brief Aux, refers to index of Mass bin, to be deprecated
   */
  int im;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index specifying the CW type. 
   * @details Set to zero when mock catalog is created
   */
  int tstruct;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  BIAS_NCOUNTS;
  //////////////////////////////////////////////////////////
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  AVERAGE_BIAS_NCOUNTS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Joint distribution normalized within each bin of DM density
   */
  vector<real_prec>  AVERAGE_BIAS_NCOUNTS_normalized;
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS only as a function of local density
   */
  vector<ULONG>  BIAS_LOCAL_DM;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  BIAS_SAT_FRACTION;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  ABUNDANCE;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<ULONG>  NCELLSperDMBIN;

  //////////////////////////////////////////////////////////
  /**
   * @brief Joint distribution normalized within each bin of DM density
   */
  vector<real_prec>  BIAS_NCOUNTS_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<real_prec>  BIAS_SAT_FRACTION_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BIAS in the form o Number counts, mass, and satellite fraction 
   */
  vector<real_prec>  ABUNDANCE_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mean values of Y in a bin of DM density
   */
  vector<real_prec> mean_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of tracers Y obtained from the input density field
   * @details Computed in BiasMT::analyze_input()
   */
  real_prec N_objects_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Number of DM particles obtained from the input density field
   * @details Computed in BiasMT::analyze_input()
   */
  real_prec N_objects_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Number of DM (reference) particles obtained from the input density field
   */
  real_prec N_objects_X_REF;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the PDF of tracers
   */
  vector<ULONG> PDF_NC_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the PDF of DM particles
   */
  vector<ULONG> PDF_NC_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of objects in mock 
   */ 
  ULONG Nobjects;

  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of X particles in one cell
   */
  int nmax_x_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int mean_number_y_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int nmax_y_onecell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Maximum number of Y particles in one cell
   */
  int nmax_y_sat_onecell;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the bins in X
   */
  vector<real_prec> X_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container for the bins in X
   */
  vector<real_prec> Y_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief This vector will keep the maximum number of particles in one cell for the different cwt asked in par file
   */
  vector<int>nmax_y_onecell_cwt;
  //////////////////////////////////////////////////////////
  /**
   * @brief 2d vector containing the number of tracers in each multi-scale level for each reference simulation
   */
  vector<vector<ULONG>> tracers_multi;

  //////////////////////////////////////////////////////////

  /**
   * @brief Generate the PDF of Halos given the DM density field
   * @details According to the type of CW requested in parameter file, this function generates, for each CWT
   * the joint and the conditional probablity distribution. This function is not used when creating a mock catalog.
   * This function is active as long as the pre-processor directive define BIAS is found
   * @details Can only be used after having called read_BiasMT_files() and analyze_input()
   * @return The output if this function are files containing
   * @return The mean Y , var Y for bins in X, for a given bin of MK and a given CWT
   * @return All values of Y in a bin of X for a given Mk and CWT
   */

  // For different realizations
#ifdef _SEVERAL_REAL_
  void get_pdf(int);
#elif !_SEVERAL_REAL_
  // Same but for the single files from the input par
  void get_pdf();
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Container loaded by the function get_fof_info()
   * containing the bin in Mk in which a given cell is found, according to the results of the FoF.
   */
  vector<ULONG> SKNOT_M_info;
  //////////////////////////////////////////////////////////
  /**
   * @brief Index of the BiasMT Iteration
   * @details Used to identify output files
   */
  ULONG iteration;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */

  int step_mass_assignemt;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the BiasMT isotropic Kernel in Fourier space,
   * @details Computed as the ratio between the reference power and the mock power, with size BiasMT::Nft/BiasMT::ndel
   */
  vector<real_prec>Kernel;
  //////////////////////////////////////////////////////////
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
  /**
   * @brief Container for the BiasMT isotropic Kernel in Fourier space, averged over the last iterations (contrlled by the parameter N_Kernel_Average)
   * @details Computed as the ratio between the reference power and the mock power, with size BiasMT::Nft/BiasMT::ndel
   */
  vector<real_prec>Average_Kernel;
  vector<real_prec>average_kernel_updated;
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>Kernel_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for power spectrum used in the construction of the BiasMT kernel
   * @details with size BiasMT::Nft/BiasMT::ndel
   */
  vector<real_prec>power_ratio_unsmoothed;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the reference DM field, used when 
   * rank ordering is applied to the DMF used to calibrate the Kernel
   * based on the PDF of a high resolution DM field.
   * This should be depracated, since we should start from one realization of the low res
   * and apply the kernel to ahnother realziation of the low res
   */
  vector<real_prec> pdf_ref;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the original DM field
   * This is modified iteration after iteration
   */
  vector<real_prec> pdf_ite;
  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the pdf of the original DM field
   * Copied from the temporary container pdf_in filled in the iteration 0
   */
  vector<real_prec> pdf_ini;
  //////////////////////////////////////////////////////////
  /**
   * @brief This function builds the kernel from the ratio
   * @brief of the ref_tr power and the mock power spectrum
   * @returns Container BiasMT::Kernel
   */
  void GetKernel(bool, real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @brief Convolution of the initial DM field with the BiasMT Kernel
   * @params in: original DM density field
   * @params out: convolved with kernel
   * @returns convolved DM with Kernel, to be used in the determination of the bias inside the BiasMT iterative process
   */
  void Konvolve(vector<real_prec> &in, vector<real_prec>&out);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void v_konvolve(vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_Fourier_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void sort_properties();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void sort_properties(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief Take the properties in structure allocating the properties in bins od theta and randomize their order
   */
   void randomize_properties();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_power_spectrum(string type);

  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @brief Read from parameter file
   
   */
  string new_Name_Property_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  string new_Name_Property_Y;
  //////////////////////////////////////////////////////////
  /**
   * @brief Read from parameter file
   */
  string Type_of_File_X;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  string Type_of_File_Y;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  s_CosmologicalParameters s_cosmo_pars;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool bin_accumulate_borders;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  bool use_iteration_ini;
  //////////////////////////////////////////////////////////

#ifdef _USE_GNUPLOT_
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_swap;
  //////////////////////////////////////////////////////////
  /**
   * @brief Objects of Gnuplot type, used to plot property_vs biasa at running time
   */
  Gnuplot gp_swap_b;
  //////////////////////////////////////////////////////////
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_swap_c;
  //////////////////////////////////////////////////////////
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_swap_m;
  //////////////////////////////////////////////////////////
  /**
   * @brief Objects of Gnuplot type, yused to plot @running time.
   */
  Gnuplot gp_swap_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_ratio;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_pdf;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance_v;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_abundance_spin;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  void plot_power(vector<s_info_in_bins>,vector<s_info_in_bins>, string, string);
  void plot_power(vector<s_info_in_bins>, string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void plot_scaling_relation_assignment(string);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void plot_scaling_relation_assignment();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void plot_scaling_relation_assignment_bias();
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Meausure the bias-prop relation 
   * @details Measrements are done for the new mock as well as for the reference, and shown with an instance of Gnuplot.
   * @details This uses methods from Catalog.
   */
  void get_mean_scaling_relation_assignment_bias(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_mean_scaling_relation_assignment_secondary_bias(string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec lss_halo_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  int number_of_modes_upgrade_kernel;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec average_ratio_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void swap_mock_properties_in_theta_bins(int);


  // *****************************************************************************************************************************************
  // *****************************************************************************************************************************************
  // *****************************************************************************************************************************************
  // *****************************************************************************************************************************************
  

public:
  //////////////////////////////////////////////////////////
  /**
   *  @brief default Class Constructor
   *  @return object of class BiasMT
   */
  BiasMT(){
    if(true==WARNING_MODELS){
      So.message_warning("More than one model for properties of DM has been allowed. Please check def.h. Cosmicatlass stops here");
      exit(1);
    }

  }
  
  //////////////////////////////////////////////////////////
  /**
   *  @brief Use this constructor to pass an object of type Params
   *  @brief This constructor initializes the value of BiasMT parameters
   *  @return object of class BiasMT
   */
  

  BiasMT(Params _par):params (_par){}


  BiasMT(Params _par, bool _BiasMT_mode):params (_par), BiasMT_mode (_BiasMT_mode)
  {
    sal.precision(8);
    sal.setf(ios::showpoint); 
    sal.setf(ios::scientific); 
    this->s_cosmo_pars=this->params.s_cosmo_pars;
    if(true==WARNING_MODELS){
      So.message_warning("More than one model for properties of DM has been allowed. Please check def.h. Cosmicatlass stops here");
      exit(1);
    }
    this->params_original=_par;

    this->counter_assigned_secondary=0;
    
#ifdef _USE_MULTISCALE_LEVEL_4_
    this->NGRID_low_mass_assignment=(params._Nft_low_l4())*(params._Nft_low_l4())*(params._Nft_low_l4());
#endif
    this->NTT=this->params._NGRID_h();
    this->NTT_low=this->params._Nft_low()*this->params._Nft_low()*(this->params._Nft_low()/2+1);
    
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _MULTISCALE_
    this->get_minimum_multiscale_property();
#endif
#endif
    

#ifdef _USE_GNUPLOT_
    this->gp_kernel<<"set border linewidth 2.0\n";
    this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";

    this->gp_power<<"set border linewidth 2.0\n";
    this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_power<<"set size square \n";

    this->gp_pdf<<"set border linewidth 2.0\n";
    this->gp_pdf<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_pdf<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_pdf<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";

    this->gp_abundance<<"set border linewidth 2.0\n";
    this->gp_abundance<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_abundance<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance<<"set size square \n";

    this->gp_abundance_v<<"set border linewidth 2.0\n";
    this->gp_abundance_v<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_abundance_v<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set size square \n";


    this->gp_abundance_rs<<"set border linewidth 2.0\n";
    this->gp_abundance_rs<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_rs<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_rs<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_rs<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_abundance_rs<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance_rs<<"set size square \n";

    this->gp_abundance_spin<<"set border linewidth 2.0\n";
    this->gp_abundance_spin<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_spin<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_spin<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_abundance_spin<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_abundance_spin<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance_spin<<"set size square \n";

#endif



    this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
    So.message_screen("Loading parameters for BiasMT. THIS HAS TO BE UPDATED> USE THEM DIRECTLY FROM THE PARAMS CLASS");
#endif

    this->used_once = true;


    // --------------------------------------------------------------------------------------------
    if(true==this->BiasMT_mode)
      {
#ifdef _GET_BiasMT_REALIZATIONS_
	ifstream nxo;
	nxo.open(file_one_cell);
#ifdef _FULL_VERBOSE_
	if(nxo.is_open())
	  So.message_screen("Reading maximum number of reference tracers in cells from file ",file_one_cell);
	else
	  {
	    So.message_screen("File with max number of tracers in one cell does not exist",file_one_cell);
	    So.message_screen("BiasMT stops here");
	    exit(0);
	  }
#endif
	int iaux; int iaux2;
	nxo>>iaux>>iaux2;
	this->nmax_y_onecell=iaux;
	this->mean_number_y_onecell=iaux2;
	nxo.close();

#ifdef _FULL_VERBOSE_
	if(iaux==0)
	  So.message_warning("Maximum number of cells = 0. End");
	else
	  So.message_screen("maximum number of reference tracers in cells =",iaux);
	So.message_screen("mean number of reference tracers in cells =",iaux2);
#endif
#endif
      }

    if(true==this->BiasMT_mode)

      this->dm_already_done=false;
    // Feed the structure for cosmological parameters
    // This is also done in LPT, so verify that these lines are also in its init par member
   


#ifdef _BIN_ACCUMULATE_
    this->bin_accumulate_borders = true;
#else
    this->bin_accumulate_borders = false;
#endif
    So.DONE();
    // Determine the number of properties used to characterize the halo bias
    this->N_bias_properties=1;  // This accounts for the dark amtter density.
    vector<string> bias_properties;
    bias_properties.push_back("Local density");
#ifdef _USE_MASS_KNOTS_
    this->N_bias_properties++;
    bias_properties.push_back("Knot-Mass");
#endif
#ifdef _USE_CWC_
    this->N_bias_properties++;
    bias_properties.push_back("T-WEB");
#endif

#ifdef _USE_VEL_KNOTS_V_
    this->N_bias_properties++;
    bias_properties.push_back("V-dispersion in knots");
#endif
#ifdef _USE_CWEB_V_
    this->N_bias_properties++;
    bias_properties.push_back("V-WEB");
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
    bias_properties.push_back("I-2");
#else
    bias_properties.push_back("ð²");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
    bias_properties.push_back("I-3");
#else
    bias_properties.push_back("ð³");
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_)|| defined (_USE_PROLATNESS_)|| defined (_USE_S2_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
    bias_properties.push_back("I-4");
#elif defined (_USE_TIDAL_ANISOTROPY_)
    bias_properties.push_back("Tidal anisotropy");
#elif defined (_USE_ELLIPTICITY_)
    bias_properties.push_back("Ellipticity");
#elif defined (_PROLATNESS_)
    bias_properties.push_back("Prolatness");
#else
    bias_properties.push_back("S²");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined(_USE_NABLA2DELTA_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
    bias_properties.push_back("IV-1");
#else
    bias_properties.push_back("Nabla² ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
    bias_properties.push_back("IV-2");
#else
    bias_properties.push_back("s²ð");
#endif
#endif
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
    this->N_bias_properties++;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
    bias_properties.push_back("IV-3");
#else
    bias_properties.push_back("s²");
#endif
#endif
#ifdef _FULL_VERBOSE_
#ifndef _ONLY_LPT_
#ifdef _DO_BiasMT_CALIBRATION_
    So.message_screen("****************************************************************************");
    So.message_screen("BiasMT is using", this->N_bias_properties,"properties to characterise the bias:");
    for(int i=0;i<bias_properties.size();++i)std::cout<<YELLOW<<bias_properties[i]<<RESET<<endl;
    So.message_screen("****************************************************************************");
    std::cout<<endl;
#endif
#endif
#endif
    this->step_mass_assignemt=0; //initialize it
  }
  // *****************************************************************************************************************************************
  // *****************************************************************************************************************************************
  // *****************************************************************************************************************************************

  //////////////////////////////////////////////////////////
  /**
   *  @brief Default constructor
   *  @return object of class BiasMT
   */
  ~BiasMT(){}

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Runs BiasMTS. Called from main function
   */
  void execute();
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Name of the file with mock number density field
   */
  string fnameMOCK;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  string file_residuals;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  real_prec residuals_it;
  real_prec residuals_abs;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<real_prec>residuals_power;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<real_prec>residuals_power_unsigned;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<int>it_power;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  GnuplotC gp;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   * In this function we asign coordinates from DM + random particles based on the mock number counts
   * @details Masses are also assigned and the collapse of randoms towards dm (to correct small scale clusterin) is also performed after that
   * Randoms are decided to be collapsed after having assigned mass, for the information on the mass can be used to model the
   * fraction of distance to aproach teh rans to their closest DM particle.
   * Catalog is then written with positions, velocities and masses.
   */
#ifdef MOCK_MODE
#ifdef _USE_OMP_
  void BiasMT_makecat(gsl_rng ** gBaseRand);
#else
  void BiasMT_makecat(gsl_rng * gBaseRand);
#endif

#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @public
   * @brief Collapse random tracers towards their closesr dark matter particle
   * @details The algorithm identifies the closest DM particle to each ranomd tracer and reduces its distance to that dm by a fraction f. See documentation.
   */
  void collapse_randoms();
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @public
   * @brief
   */
  void assign_tracer_mass();
#endif
  //////////////////////////////////////////////////////////
#ifdef _USE_LPT_
#ifdef MOCK_MODE
  /**
   * @public
   * @brief
   */
  void move_particles_lpt();
#endif
#endif
  //////////////////////////////////////////////////////////
#ifdef MOCK_MODE
  /**
   * @public
   * @brief
   */
  void correct_for_exclusion(ULONG lenght_dm_bin);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<s_nearest_cells> ncells_info;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<gsl_real>mass_bins;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  vector<gsl_real>mfunc;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec prop_min;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  real_prec prop_max;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  real_prec min_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  real_prec min_mean_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  void set_So(ScreenOutput new_So){this->So=new_So;}
  /////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  ULONG counter_assigned_secondary;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  void get_IC_from_particles();
};
#endif

