////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef CATALOG__
#define CATALOG__
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<Catalog>
 * @brief Header file for the class Catalog::
 * @file Catalog.h
 * @title Class in charge of reading and analysing input catalogs
 * @details The methods in this class are desinged to generate different statistics from a catalog of dark matter tracers
 * @author Andres Balaguera-Antolínez
 * @version 2.0
 * @date    2008-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <math.h>
#include <fstream>
#include "bstream.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <string.h>
#include <cassert>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <netinet/in.h>
#include "NumericalMethods.h"
#include "Params.h"
#include "Miscelanious.h"
#include "CoordinateSystem.h"
#include "FileOutput.h"
#include "DensityProfiles.h"
#include "Statistics.h"
#include "PowerSpectrumTH.h"
#include "Cwclass.h"
#include "ScreenOutput.h"
//////////////////////////////////////////////////////////
using namespace std;

//////////////////////////////////////////////////////////
class Catalog{
private:
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief   Object of type Class::Cosmology
   * @details Calculations of redshift dependent cosmological functions
   */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Class::Cwclass
   * @details  Information of cosmic-web classification
   */
  Cwclass cwclass;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Class::stats
   * @details Calculations of statistical properties of dark matter tracers
   */
  Statistics stats;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of the class PowerSpectrum
   */
  PowerSpectrum ps;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Class::FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Class:ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Structure for cosmological parameters
   */
  s_CosmologicalParameters s_cosmo_pars;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Structure for derived cosmological information
   */
  s_CosmoInfo s_cosmo_info;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Container for comoving distance, used to interpolate
   */
  vector<gsl_real>rcosmo;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Container for redshift, used to interpolate
   */
  vector<gsl_real>zcosmo;

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Mean mass of the sample
   */
  real_prec mean_Mass;

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing x coordinates
   */
  vector<real_prec> xgal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing y-coordinates
   */
  vector<real_prec> ygal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing z-coordinates
   */
  vector<real_prec> zgal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing x-momemtum
   */
  vector<real_prec> Pxgal; //P meant for vector properties such as velocity or momentum
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing y-momemtum
   */
  vector<real_prec> Pygal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing z-momemtum
   */
  vector<real_prec> Pzgal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing something
   */
  vector<real_prec> VPgal;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing Mass
   */
  vector<real_prec> Mass;  // other scalar property such as the mass
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing the input parameters
   */
  vector<real_prec> property;  // other scalar

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing the input parameters
   */
  vector<real_prec>ABUNDANCE;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing the input parameters
   */
  vector<real_prec>ABUNDANCE_normalized;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  vector containing the input parameters
   */
  vector<real_prec> weight;
  //////////////////////////////////////////////////////////
  /**
   * @brief  File
   */
  string Output_directory;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to divide the sample
   */
  int NMASSbins;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  int NMASSbins_mf;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  int NMBINS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  Gnuplot gp_abundance;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of mass bins to measure mass function
   */
  Gnuplot gp_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief Identifies dark matter "DM" or tracer "TR"
   */
  string type_of_object;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  s_params_box_mas box;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  s_params_box_mas box_HR;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  s_params_box_mas box_low;
  s_params_box_mas box_n;  //new, auxiliary
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of objects load by the read file or read_file_bin objects
   */
  ULONG NOBJS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of objects chosen as a random subsample
   */
  ULONG NOBJS_subsample;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec normal_b;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec normal_p;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec sg2;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec sg1;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec s_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec w_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_property_value;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG n_gal;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_density;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_vmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_vmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_vmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_vmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_rs;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_rvir;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_rvir;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_rvir;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_rvir;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_concentration;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_concentration;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_concentration;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_concentration;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_spin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_spin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_spin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_spin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_virial;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_virial;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_virial;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_virial;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_vrms;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_vrms;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_vrms;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_vrms;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_b_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_b_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_b_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_b_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_c_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_c_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_c_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_c_to_a;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_mach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_mach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_mach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_mach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_dach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_dach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_dach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec var_dach;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_bias;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_local_overdensity;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_local_overdensity;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_tidal_anisotropy;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_tidal_anisotropy;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_local_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_local_dm;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_peak_height;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec max_peak_height;
  //////////////////////////////////////////////////////////

  /**
   * @brief
   */
  vector<ULONG>Number_of_tracers_in_mass_bins;
  //////////////////////////////////////////////////////////

  /**
   * @brief
   */
  vector<ULONG>Number_of_tracers_in_vmax_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<ULONG>Number_of_tracers_in_ph_bins;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_secondary_intervals(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_cosmo();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  const gsl_rng_type *Tn_cat;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */

  gsl_rng *rn_cat;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
public:
  //////////////////////////////////////////////////////////
  /**
   * @brief Constructor
   */
  Catalog():logdeltaRS(0),logdeltaVMAX(0),logdeltaM(0),logdeltaM_low(0),var_prop(0),mean_Mass(0),mean_number_density(0),min_halo_separation(0),NOBJS(0),Ntracers_ran(0),Ntracers_dm(0), aux_flag(true)
  {

    this->normal_b=0;
    this->normal_p=0;
    this->s_g=0;
    this->sg1=0;
    this->sg2=0;
    this->w_g=0;
    this->n_gal=0;
    this->mean_property_value=0;
    this->type_of_object=this->params._type_of_object();

    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Constructor (overloaded)
   * @details
   * @arg Object of type Params with the parameters read from the input parameter file
   */
  Catalog(Params _params):params (_params), logdeltaRS(0),logdeltaVMAX(0),logdeltaM(0),logdeltaM_low(0),var_prop(0),mean_Mass(0),mean_number_density(0),min_halo_separation(0),NOBJS(0),Ntracers_ran(0),Ntracers_dm(0),aux_flag(true)
  {
    // This constructor can be used when the type is declared
    // and used directly. Otherwise, if a class has an object of
    // thi type, we define the object calling the default
    // constructor and calling set_params_catalog to assign values to this class, specifçically, to a structure

    this->normal_b=0;
    this->normal_p=0;
    this->s_g=0;
    this->sg1=0;
    this->sg2=0;
    this->w_g=0;
    this->n_gal=0;
    this->mean_property_value=0;
    this->type_of_object=this->params._type_of_object();

    this->Output_directory=this->params._Output_directory();
    this->box.masskernel=this->params._masskernel();
    this->box.Lbox=this->params._Lbox();
    this->box.Nft=this->params._Nft();
    this->box.NGRID=this->params._NGRID();
    this->box.d1=this->params._d1();
    this->box.d2=this->params._d2();
    this->box.d3=this->params._d3();
    this->box.d1_HR=this->params._d1_HR();
    this->box.d2_HR=this->params._d2_HR();
    this->box.d3_HR=this->params._d3_HR();
    this->box.min1=this->params._xmin();
    this->box.min2=this->params._ymin();
    this->box.min3=this->params._zmin();

    this->box_n.masskernel=this->params._masskernel();
    this->box_n.Lbox=this->params._Lbox();
    this->box_n.Nft=this->params._Nft_low();
    this->box_n.NGRID=this->params._NGRID_low();
    this->box_n.d1=this->params._d1_low();
    this->box_n.d2=this->params._d2_low();
    this->box_n.d3=this->params._d3_low();
    this->box_n.min1=this->params._xmin();
    this->box_n.min2=this->params._ymin();
    this->box_n.min3=this->params._zmin();

    this->box_JK.masskernel=this->params._masskernel();
    this->box_JK.Lbox=this->params._Lbox();
    this->box_JK.Nft= this->params._Nft_JK();
    this->box_JK.NGRID=this->params._NGRID_JK();
    this->box_JK.d1=this->params._d1_JK();
    this->box_JK.d2=this->params._d2_JK();
    this->box_JK.d3=this->params._d3_JK();
    this->box_JK.min1=this->params._xmin();
    this->box_JK.min2=this->params._ymin();
    this->box_JK.min3=this->params._zmin();


    this->s_cosmo_pars=this->params.s_cosmo_pars;
    this->Cosmo.set_cosmo_pars(this->s_cosmo_pars);

#ifdef _USE_MULTISCALE_LEVEL_4_
    this->box_low.masskernel=this->params._masskernel();
    this->box_low.Lbox=this->params._Lbox();
    this->box_low.Nft=this->params._Nft_low_l4();
    this->box_low.d1=this->params._d1_low();
    this->box_low.d2=this->params._d2_low();
    this->box_low.d3=this->params._d3_low();
    this->box.min1=this->params._xmin();
    this->box.min2=this->params._ymin();
    this->box.min3=this->params._zmin();
#endif
    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;

    gsl_rng_env_setup();
    gsl_rng_default_seed=1015;
    this->Tn_cat =  gsl_rng_mt19937 ;//gsl_rng_default;
    rn_cat = gsl_rng_alloc (this->Tn_cat);


  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Default constructor
   */
  ~Catalog(){}

  //////////////////////////////////////////////////////////
  /**
   * @brief Get minimum of property named prop
   * @arg string prop
   * @return Minimum from tracer.Halo[].prop
   */
  real_prec get_min(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Get maximum of property named prop
   * @arg string prop
   * @return Maximum from tracer.Halo[].prop
   */
  real_prec get_max(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_intervals_equal_number(string prop,vector<real_prec>&min_aux,vector<real_prec>&max_aux);
  //////////////////////////////////////////////////////////
  /**
   * @brief Main function used to analyze catalog
   * @param read true  if read catalog, false: do no read: do tasks from input dens fields
   * @note Check the function while develpment: there are methods expecting orders from parameter file
   * @warning Ascii files are expected. Biinary files alre also expected with data_scheme as in BAM
   */
  void analyze_cat(bool read);
  //////////////////////////////////////////////////////////
  /**
   * @brief Writes catalog with N_PROP columns (see def.h) to a binary file
   */
  void write_catalog_bin(string output_file);
  //////////////////////////////////////////////////////////
  /**
   * @brief Computes the variance from a give halo property
   */
  real_prec get_variance(real_prec, string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Pair to return mean and variance
   */
  pair<real_prec, real_prec> get_variance(string prop, bool);
  //////////////////////////////////////////////////////////
  /**
   * @brief Writes catalog with N_PROP columns (see def.h) to an ascii file
   */
  void write_catalog(string output_file);
  //////////////////////////////////////////////////////////
  /**
   * @brief Writes catalog with N_PROP columns (see def.h) to an ascii file
   */
#if defined _USE_HYBRID_ASSIGNMENT_NEW_ || defined _USE_HYBRID_ASSIGNMENT_NEW_NEW_
  void select_random_subsample(real_prec fraction);
#else
  void select_random_subsample(string output_file, bool write_to_file, real_prec fraction);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Selects random sub-sample from catalog
   * @arg fraction:  fraction of the total number of tracers to be selected
   * @arg fname: output file name for reduced catalog
   */
  void select_random_subsample(real_prec fraction, string fname);
  //////////////////////////////////////////////////////////
  /**
   * @brief Selects random sub-sample from catalog
   * @arg fraction:  fraction of the total number of tracers to be selected
   */
  void select_random_subsample(real_prec fraction, int, vector<real_prec>&, string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief Selects ranbdom sub-sample from catalog
   * @arg fraction:  fraction of the total number of tracers to be selected
   */
  void select_random_subsample_v(real_prec fraction);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Spherical overdensity of the tracer
   */
  void Get_SO_tracer();
  //////////////////////////////////////////////////////////
  /**
   * @brief Reads particle catalog in ascii format.
   */
#ifdef _USE_MASS_CUTS_PK_
  void read_catalog(string input_file, real_prec mcut);
#elif defined (_USE_MASS_BINS_PK_)
  void read_catalog(string input_file, real_prec m_min, real_prec m_max);
#else
  void  read_catalog(string input_file, real_prec aux);
#endif
  void  read_catalog(string input_file, int, int, int);

  //////////////////////////////////////////////////////////
  /**
   * @brief Read input catalog in binary format, with x, y, z (and velocities if requested) in three separate files of lenght n
   */
#ifdef _USE_VELOCITIES_
  void read_catalog_bin(ULONG n, string fx, string fy, string fz,string fvx, string fvy, string fvz); // reads binary from fastpm
#else
  void read_catalog_bin(); // reads binary from fastpm
  void read_catalog_bin_tng(); // reads binary from fastpm
#endif
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the interpolated prop-density field and write to output_file
   * @details prop can be: _MASS_, _COUNTS_. TBD: must be generalized
   */
  void get_density_field_grid(string prop, string output_file);
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the interpolated prop-density and write to vector
   */
  void get_density_field_grid(string prop,vector<real_prec>&);
  /**
   * @brief Get the interpolated prop-density and write to vector
   */
  void get_density_field_grid(s_params_box_mas , string, vector<real_prec>&);
  /**
   * @brief Get the interpolated prop-density and write to vector
   */  
  void get_density_field_grid(std::string prop, vector<real_prec>&deltaTR_counts,vector<real_prec>&out);
  //////////////////////////////////////////////////////////
  /**
   * @brief Measure the abundance as a funvtion of a property
   * @args prop: string denoting the property
   */
  void get_property_function(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Measure the abundance as a function of a property in different CWT
   * @args prop: string denoting the property
   */
  void get_property_function_cosmic_web_types(string );
  //////////////////////////////////////////////////////////
  /**
   * @brief This method defines bins for properties.
   * @details Whether the property is read or not, this method genertes arrays for the properties
   * @details These are for the min_i, max_i and the center on the bin_i. The size of the bin is controled
   * with the parameter parms._NMASSbins(). A delta in log is defeind for all properties.
   * @warning In BAM, when assining properties, the size of the bins to allocate the container ABUNDANCE is controled by the quantity this->params._NPROPbins_bam()
   */
  void define_property_bins();
  //////////////////////////////////////////////////////////
  /**
   * @brief Thi method defines the bins in properties
   */
  void define_property_bins(ULONG Nbins, string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Measure the distributio of reduced mass in cells
   */
  void get_distribution_reduced_mass_in_cell();
  //////////////////////////////////////////////////////////
  /**
   * @brief Meausure the distributiuon of minimum separations
   * @details This is computed in a spheres of radii params._scale_mach_number() around each tracer
   */
  void get_distribution_min_separations(vector<s_nearest_cells> &nearest_cells_info);
  //////////////////////////////////////////////////////////
  /**
   * @brief Get distribution of separation within cells

   */
  void get_stats_separation_in_cell();

  //////////////////////////////////////////////////////////
  /**
   * @brief This method converts a snapshot into a mock catalog
   * @details Converts cartesian positions into RA, dec, and Z using cosmological and pecular redshifts.
   * @parameter  zsel se tto true if a dNdz is to be imposed, false otherwise 
   */
  void snap_to_mock(bool zsel);
  //////////////////////////////////////////////////////////
  /**
   * @brief Determine the ID of the neighbour tracers to each tracer
   */

  void get_neighbour_tracers(vector<s_nearest_cells> &nearest_cells_info );
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_masses_of_pairs_in_min_separation_bin_in_theta_bin(real_prec min_sep, vector<s_mass_members> & dm_properties_bins);
  //////////////////////////////////////////////////////////
  /**
   * @brief This function computs occupation numbers in cells
   * @details This function computs occupation numbers in cells in bins of a given property
   */
  void get_pdf_vmax(string);
  //////////////////////////////////////////////////////////
  /**
   * @brief This function generate galaxy catlaogs usingn halos and HOD recipe
   */
  void get_sep_distribution(real_prec min_value_prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief This function generate galaxy catlaogs usingn halos and HOD recipe
   */
  void halos2galaxies_HOD();
  //////////////////////////////////////////////////////////
  /**
   * @brief This function generate galaxy catlaogs usingn halos and HOD recipe
   */
  void get_HOD_web_types(string );
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the halo occupation distribution (Number of tracers as a function of halo mass)
   * @arg type: central, satellite
   */
  void get_HOD(string type);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the local MACH numbr in each cell
   * @details This uses velocities and grid ID
   * @arg rel: true if relative. If set to false, it computes the absolute mach number
   */
  void get_local_mach_number(bool rel);
  //////////////////////////////////////////////////////////
  /**
   * @brief This methods computs the local MACH numbr in each objectl
   * @details Uses veocity of tracrs within a radius R=Scale from each tracer
   * @details To that, the algorithm first identifies neighbouring cells (and tracers inside) to the cells where each tracer is
   * @details and then loops over them, avoiding to loop over the full set of tracers
   */
  void get_local_mach_number(real_prec);
  /**
   * @brief Same as this->get_local_mach_number
   * but using chiunks of slices to improve memmoery allocation
   */
  void get_local_mach_number_chuncks(real_prec);

  //////////////////////////////////////////////////////////
  /**
   * @brief  Main container to keep the input catalog
   * @notes This memebr will be the one used to assign external catalog and make statistics
   */
  vector<s_Halo> Halo;
  //////////////////////////////////////////////////////////
  /**
   * @brief Random catalog built based on the properties of input catalog
   */
  vector<s_Halo> Random;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Container to keep the input catalog
   */
  vector<s_Halo> Halo_random_subsample;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> tracer;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> Halo_ranked;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> Halo_randomized;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> tracer_aux;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> galaxy;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> galaxy_central;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_Halo> galaxy_satellite;
  //////////////////////////////////////////////////////////
  /**1
   * @brief
   */
  real_prec logdeltaM;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec logdeltaCONCENTRATION;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in log Mass
   */
  real_prec logdeltaM_low;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in log Vmax
   */
  real_prec logdeltaVMAX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in log Rs
   */
  real_prec logdeltaRS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in log Spin
   */
  real_prec logdeltaSPIN;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in log deltaDM
   */
  real_prec deltaLOCALDM;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in Mach
   */
  real_prec deltaMACH;
  //////////////////////////////////////////////////////////
  /**
   * @brief Size of bin in Neighbour statististics
   */
  real_prec deltaDACH;
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves x-coordinate of i-th galaxy in catalog
   */

  real_prec _x(ULONG i){return this->Halo[i].coord1;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves y-coordinate of i-th galaxyin catalog
   */

  real_prec _y(ULONG i){return this->Halo[i].coord2;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves z-coordinate of i-th galaxy in catalog
   */
  real_prec _z(ULONG i){return this->Halo[i].coord3;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  Retreaves vx-coordinate of i-th galaxy in catalog
   */
  real_prec _vx(ULONG i){return this->Halo[i].vel1;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves vy-coordinate of i-th galaxy in catalog
   */
  real_prec _vy(ULONG i){return this->Halo[i].vel2;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves vy-coordinate of i-th galaxy in catalog
   */
  real_prec _vz(ULONG i){return this->Halo[i].vel3;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves mass of i-th galaxy in catalog
   */
  real_prec _mass(ULONG i){return this->Halo[i].mass;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves mass of i-th galaxy in catalog
   */
  ULONG _Number_of_tracers_in_mass_bins(int i)
  {
    return this->Number_of_tracers_in_mass_bins[i];
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves mass of i-th galaxy in catalog
   */
  ULONG _Number_of_tracers_in_vmax_bins(int i)
  {
    return this->Number_of_tracers_in_vmax_bins[i];
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Retreaves mass of i-th galaxy in catalog
   */

  ULONG _Number_of_tracers_in_ph_bins(int i)
  {
    return this->Number_of_tracers_in_ph_bins[i];
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  sample variance of the tracer property
   */
  real_prec var_prop;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> MBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> VMAXBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> RSBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> SPINBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> DACHBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> MACHBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> LOCALDMBmin; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> DACHBmax; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> MACHBmax; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> LOCALDMBmax; // +1 to include the full sample
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> CONCENTRATIONBmin; // +1 to include the full sample
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> MBmax;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> VMAXBmax;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> RSBmax;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> SPINBmax;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> CONCENTRATIONBmax;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> MBin;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> VMAXBin;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> RSBin;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> SPINBin;

  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<gsl_real> CONCENTRATIONBin;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<real_prec> Dist_Min_Sepatations;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<int>Number_of_neighbours;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<real_prec>local_overdensity;

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<real_prec>mach_number;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<real_prec>dach_number;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  int NCOLS;
  //////////////////////////////////////////////////////////
  /**
   * @brief Numbner density of the catalog
   */
  real_prec mean_number_density;

  //////////////////////////////////////////////////////////
  /**
   * @brief Container for the mass function of the catalog
   */
  vector<real_prec> mass_function;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Container for the vmax function of the catalog
   */
  vector<real_prec> vmax_function;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Container for the rs_function of the catalog
   */
  vector<real_prec> rs_function;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Container for the rs_function of the catalog
   */
  vector<real_prec> cvir_function;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Container for the spin_function of the catalog
   */
  vector<real_prec> s_function;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containig the minimum tracer separation in cells
   */
  vector<real_prec>min_separation_in_cell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containig the mean tracer separation in cells
   */
  vector<real_prec>mean_separation_in_cell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containig the mean tracer separation in cells
   */
  vector<real_prec>stdv_separation_in_cell;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector containig the stdev of tracer separation in cells
   */
  vector<real_prec>sigma_separation_in_cell;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_mean_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec min_stdv_halo_separation;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<s_dist_in_dmbins> masses_in_cells_min_sep;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Class containing the input parameters
   */
  Params params;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold;

  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold_multi_scale_1;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold_multi_scale_2;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold_multi_scale_3;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold_multi_scale_4;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 0 of multi-scaling
   */
  ULONG N_props_0;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 1 of multi-scaling
   */
  ULONG N_props_1;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 2 of multi-scaling
   */
  ULONG N_props_2;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 3 of multi-scaling
   */
  ULONG N_props_3;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Number of properties in the Level 4 of multi-scaling
   */
  ULONG N_props_4;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG Ntracers_ran;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG Ntracers_dm;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec fraction_tracer_from_dm;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec fraction_tracer_from_random;
  ///////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec Prop_threshold_rand_dm;
  ///////////////////////////////////////////////////////
  /**
   * @brief  Auxiliary
   */
  bool aux_flag;
  //////////////////////////////////////////////////////////
  /**
   * @brief Set type of tracer
   */
  void set_type_of_object(string new_type_of_object){this->type_of_object=new_type_of_object;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Retrieve type of tracer
   */
  string _type_of_object(){return this->type_of_object;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Set minimum value of bias
   */
  //////////////////////////////////////////////////////////
  void set_min_bias(real_prec new_b){this->min_bias=new_b;}
  /**
   * @brief Retrieve minimum value of bias
   */
  real_prec _min_bias(){return this->min_bias;}
  //////////////////////////////////////////////////////////

  void set_max_bias(real_prec new_b){this->max_bias=new_b;}
  real_prec _max_bias(){return this->max_bias;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG _NOBJS(){return this->NOBJS;}
  int _NCOLS(){return this->NCOLS;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_NOBJS(ULONG new_NOBJS) {this->NOBJS=new_NOBJS;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */

  ULONG _NOBJS_subsample(){return this->NOBJS_subsample;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_NOBJS_subsample(ULONG new_NOBJS) {this->NOBJS_subsample=new_NOBJS;}

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_Number_of_tracers_in_mass_bins(int im, ULONG Ntr) {this->Number_of_tracers_in_mass_bins[im]=Ntr;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_Number_of_tracers_in_vmax_bins(int im, ULONG Ntr) {this->Number_of_tracers_in_vmax_bins[im]=Ntr;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_Number_of_tracers_in_ph_bins(int im, ULONG Ntr) {this->Number_of_tracers_in_ph_bins[im]=Ntr;}

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_intervals_equal_number_aux(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_tracer_tidal_anisotropy(vector<real_prec>&tidal);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @author ABA
   */
  void get_peak_height_at_tracer();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @author ABA
   */
  real_prec pearson_correlation(string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @author ABA
   */
  real_prec spearman_correlation(string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @author ABA
   */
  void Get_Ranked_Props(string);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @author ABA
   */
  void Get_Ranked_Props(vector<s_Halo>&,  string);
  //////////////////////////////////////////////////////////

  /**
   * @brief Principal component analysis for halo properties
   * @input Vector of strings containing the names of used properties
   * @input Vector of bool identifying used/unsed properties
   * @input bool to write out downsampled catalog
   * @details Principal component analysis for halo properties. Identify each halo as the linear combination of all its properties in the form
   * \f$ h=\sum_i \lambda_i e_i\f$ where \f$e_i\f$ is the multidimensional basis of halo properties and \f$\lambda _i\f$ are the standarized halo properties, e.g, mass ->to (Mass-mean_Mass)/standard_dev_Mass.
   * The argument holds the name of halos properties and whether they have been used.
   * @author ABA
   */
  void PCA(vector<string>&, vector<bool>&, string, bool);
  //////////////////////////////////////////////////////////
  /**
   * @brief Function to set parameters
   */
  void set_params(Params new_params){
    this->params=new_params;
    this->Output_directory=params._Output_directory();
    this->box.masskernel=this->params._masskernel();
    this->box.Lbox=this->params._Lbox();
    this->box.Nft=this->params._Nft();
    this->box.NGRID=this->params._NGRID();
    this->box.d1=this->params._d1();
    this->box.d2=this->params._d2();
    this->box.d3=this->params._d3();
    this->box.min1=this->params._xmin();
    this->box.min2=this->params._ymin();
    this->box.min3=this->params._zmin();
    this->type_of_object=this->params._type_of_object();

    this->box_JK.masskernel=this->params._masskernel();
    this->box_JK.Lbox=this->params._Lbox();
    this->box_JK.Nft= this->params._Nft_JK();
    this->box_JK.NGRID=this->params._NGRID_JK();
    this->box_JK.d1=this->params._d1_JK();
    this->box_JK.d2=this->params._d2_JK();
    this->box_JK.d3=this->params._d3_JK();
    this->box_JK.min1=this->params._xmin();
    this->box_JK.min2=this->params._ymin();
    this->box_JK.min3=this->params._zmin();

    this->Cosmo.set_cosmo_pars(this->params.s_cosmo_pars);


    this->box_n.masskernel=this->params._masskernel();
    this->box_n.Lbox=this->params._Lbox();
    this->box_n.Nft=this->params._Nft_low();
    this->box_n.NGRID=this->params._NGRID_low();
    this->box_n.d1=this->params._d1_low();
    this->box_n.d2=this->params._d2_low();
    this->box_n.d3=this->params._d3_low();
    this->box_n.min1=this->params._xmin();
    this->box_n.min2=this->params._ymin();
    this->box_n.min3=this->params._zmin();
  }

  //////////////////////////////////////////////////////////
  /**
   * @brief Get Random catalog based
   * @details Get random catalg using the distributions read from the tracer catalog. Randoms are kept in the Random  vector of structures.
   *
   */
  void get_random_catalog();

  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  void get_scaling_relation_primary_property(string prop);
  //////////////////////////////////////////////////////////

  /**
   * @brief Assignment of a halo property based on mass and a second property prop_env
   */
  void get_scaling_relation_primary_property(string prop, string prop_env);
  //////////////////////////////////////////////////////////

  /**
   * @brief Assignment of a halo property based on mass, a second property prop_env and a third property prop_extra
   */
  void get_scaling_relation_primary_property(string prop, string prop_env, string prop_extra);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property based on mass, a second property prop_env and a third property prop_extra
   */
  void get_scaling_relation_bias();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void get_superclusters(string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type s_params_box_mas
   */
  s_params_box_mas box_JK;  //new, auxiliary
  ///////////////////////////////////////////////////////
  /**
   * @brief Transformation of coordinates of input catalogues to cartessian coordinates
   * @details returning the same input vector in which the three first columns correspond to
   * to the X,Y,Z coordinates and in the i_nbar (see input parameter file)
   * the mean number density is written.
   * @param s_d structure containing information related to the catalogue
   * @result Catalogue with position in cartesian coordinates and mean number density
   * tabulated in the corresponding column as stated in the parameter file
   */
#ifdef _USE_SEVERAL_RANDOM_FILES_
  void ang_to_cart_coordinates(s_data_structure *, int);
#else
  void ang_to_cart_coordinates(s_data_structure *);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  void get_mean_number_density(real_prec alpha, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
#ifdef  _USE_SEVERAL_RANDOM_FILES_
  void get_interpolated_density_field(bool marked, string property, int ); // meant for redshift space
#else
  void get_interpolated_density_field(bool marked, string property); // meant for redshift space
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief THis method computes bias from individual asignment as a function of a property
   */
  void get_mean_bias_relation(s_info_in_bins &);
  //////////////////////////////////////////////////////////
  /**
   * @brief THis method computes secondary bias from individual asignment as a function of a primary property
   */
  void get_mean_secondary_bias_relation(vector<s_info_in_bins>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Define intervals in property once the Nft levels are defined
   */
  void get_intervals_multiscale(string prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */

  void get_interpolated_density_field_real_space(bool marked, string property);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  void get_interpolated_density_field_real_and_redshift_space(bool marked, string property);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  vector <real_prec> field_external;
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  vector <real_prec> field_external_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  vector <real_prec> field_external_marked;
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of a halo property prop based on mass
   */
  vector <real_prec> field_external_marked_s;
  //////////////////////////////////////////////////////////
  /**
   * @brief These are properties used in the estimation of power spectrum
   */
  real_prec _normal_b(){return this->normal_b;}
  //////////////////////////////////////////////////////////
  /**
   * @brief These are properties used in the estimation of power spectrum
   */
  real_prec _normal_p(){return this->normal_p;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _sg2(){return this->sg2;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _sg1(){return this->sg1;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _s_g(){return this->s_g;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _w_g(){return this->w_g;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _mean_property_value(){return this->mean_property_value;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG _n_gal(){return this->n_gal;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _mean_density(){return this->mean_density;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Bring an external catalog and make it class member.
   * @notes By default we assign the new catalog to the class member Halo
   */
  void set_tracer_catalog(vector<s_Halo> &new_cat){this->Halo=new_cat;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_max_mach(real_prec mm){this->max_mach=mm;}
  real_prec _max_mach(){return this->max_mach;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_min_mach(real_prec mm){this->min_mach=mm;}
  real_prec _min_mach(){return this->min_mach;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_max_local_overdensity(real_prec mm){this->max_local_overdensity=mm;}
  real_prec _max_local_overdensity(){return this->max_local_overdensity;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void set_min_local_overdensity(real_prec mm){this->min_local_overdensity=mm;}
  real_prec _min_local_overdensity(){return this->min_local_overdensity;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _min_vmax(){return this->min_vmax;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _max_vmax(){return this->max_vmax;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _min_spin_bullock(){return this->min_spin_bullock;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _max_spin_bullock(){return this->max_spin_bullock;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _max_tidal_anisotropy(){return this->max_tidal_anisotropy;}
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec _min_tidal_anisotropy(){return this->min_tidal_anisotropy;}

};


#endif
