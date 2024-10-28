////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 * @brief File containing definitions of structures
 * @file Type_structures_def.h
 * @author ABA
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef _STRUCTURES__
#define _STRUCTURES__
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>
# include "def.h"  // Do not define it after gnuplot or python
using namespace std;
////////////////////////////////////////////////////////////////////////////
/**
 * @struct<s_CosmologicalParameters>
 * @brief The s_CosmoInfo struct
 * @details Auxiliary structure containg cosmological parameters
 */
struct s_CosmologicalParameters
{
  real_prec cosmological_redshift;
  /**
   *@brief  Energy density of total matter in units of the the critical density */
  real_prec Om_matter;     
  /**
   *@brief  Energy density of cold_dark_matter in units of the the critical density */
  real_prec Om_cdm;      
  /**
   *@brief  Energy density of radiation in units of the the critical density */
  real_prec Om_radiation;  
  /**
   *@brief  Energy density of baryons in units of the the critical density */
  real_prec Om_baryons;    
  /**
   *@brief  Energy density of dark energy in units of the the critical density */
  real_prec Om_vac;        
  /**
   *@brief  Energy density of curvature in units of the the critical density */
  real_prec Om_k;          
  /**
   *@brief  Baryon fraction*/
  real_prec f_baryon;
  /**
   *@brief  Hubble parameter in units of h km / s / Mpc */
  real_prec Hubble;        
  /**
   *@brief  dimensionless Hubble parameter*/
  real_prec hubble;        
  /**
   *@brief  Equation of state of dark energy */
  real_prec wde_eos;         
  /**
   *@brief  Effective number of relativistic particles*/
  real_prec N_eff;         
  /**
   *@brief  RMS of matter fluctuations at R = 8 Mpc/h*/
  real_prec sigma8;        
  /**
   *@brief   Amplitude of primordial ower spectrum*/
  real_prec A_s;            
  /**
   *@brief  Primordial spectral index*/
  real_prec spectral_index;           
  /**
   *@brief  primordial spectral index */
  real_prec alpha_s;           
  /**
   *@brief  defined from patchy*/
  real_prec D1;   
  /**
   *@brief  defined from patchy*/
  real_prec D2;   
  /**
   *@brief  CMB temperature*/
  real_prec Tcmb;  
  /**
   *@brief  Something related to a magnitude, used for vmax stuff*/
  real_prec Mabs;   
  /**
   *@brief  Something related to a magnitude, used for vmax stuff*/
  real_prec mlim;   
  /**
   *@brief  Comoving scale to compute the rms of mass fluctuation. Usually 8mpc/h */
  real_prec RR;           
  /**
   *@brief  Spherical overdensity*/
  real_prec Delta_SO;     
  /**
   *@brief  Use baryonic wiggles in the linear P(k) */
  bool use_wiggles;        
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_a;  
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_b;  
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_c;  
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_d;  
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_e;  
  /**
   *@brief   constant factor in the K-correction */
  real_prec K_index_f;  
  /**
   *@brief   constant factor in the K-correction*/
  real_prec K_index_g;  
  /**
   *@brief   constant factor in the K-correction*/
  real_prec K_index_h;  
  /**
   *@brief   constant factor in the K-correction*/
  real_prec K_index_i;  
  /**
   *@brief   constant factor in the K-correction*/
  real_prec K_index_j;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_a;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_b;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_c;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_d;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_e;  
  /**
   *@brief   constant factor in the e-correction */
  real_prec e_index_zstar;  
  /**
   *@brief  Beta factor used to mimic RSD*/
  real_prec beta_rsd;   
  /**
   *@brief   constant factor in the k-correction */
  real_prec k_index;  
  /**
   *@brief   related to the density evolution */
  real_prec d_index; 
  /**
   *@brief  Galaxy bias, used in Cl */
  real_prec GAL_BIAS; 
  /**
   *@brief  Exponent defining bias */
  real_prec alpha_BIAS;
  // aca van variables que estan asociadas a integracion en k, o
  // que van servir para poner más parámetros y asi pasar una sola
  // estructura a las rutinas de integracion
  real_prec kmin_int;
  real_prec kmax_int;
  string mass_function_fit;
  string halo_mass_bias_fit;
  string density_profile;
/**
   *@brief  values of integrals used in the HaloFit for intergaration wrt to z of P(k,z)
*/
  real_prec h4; 
  real_prec h2;
  real_prec M_max_mf;
  real_prec M_min_mf;
  real_prec kmin_ps;
  real_prec kmax_ps;
/**
   *@brief  use it for m*/
  real_prec aux_var1; 
/**
   *@brief  use it for m*/
  real_prec aux_var2; 
/**
   *@brief  use it for z*/
  real_prec aux_var3; 
/**
   *@brief  use it for k*/
  real_prec aux_var4;
/**
   *@brief  Non linaer mass scale*/
  real_prec Mnl;
/**
   *@brief  Non linear wavenumber*/
  real_prec knl;
/**
   *@brief  Non linear scale*/
  real_prec rnl;
/**
   *@brief  Non linear wavenumber from Halo-fit*/
  real_prec knl_hf;
  /**
   *@brief  Non linear scale from Halo-fit*/
  real_prec rnl_hf;
  /**
   *@brief  Non linear wavenumber from Halo-fit*/
  real_prec kstar; 
  real_prec M_min_effective;
  real_prec M_max_effective;
  // Q-model related parameter
  real_prec A_PS;
  real_prec Q_PS;
  real_prec Amc; //Amplitude of the non linear correctionto P(k) in PT
  // Redshift dependent quantities
  real_prec critical_density;
  real_prec density_contrast_top_hat;
  real_prec mean_matter_density;
  real_prec growth_factor;
  real_prec pk_normalization;
  // Integration with respect to Mass
  int n_points_mass;
  int n_points_mass_integration;
  /**
   *@brief  Vector used to speed up calculations by means of interpolation. Must come in gsl_real precision if used in interpolations */
  vector<gsl_real> v_mass; 
  /**
   *@brief  Container for sigma_mass */
  vector<gsl_real> v_sigma_mass; 
  /**
   *@brief  Container for abundance (mass) */
  vector<gsl_real> v_mass_function; 
  /**
   *@brief  Container for bias as a functin of mass */
  vector<gsl_real> v_halo_mass_bias; 
  /**
   *@brief  Container for wavenumbers in P(k) */
  vector<gsl_real> v_k_ps;  
  /**
   *@brief  Container for linear P(k) */
  vector<gsl_real> v_lin_power_spectrum; 
  /**
   *@brief  Container for non-linear P(k) */
  vector<gsl_real> v_nl_power_spectrum; 
  /**
   *@brief  Container for density profile in Fourier space */
  vector<gsl_real> v_density_profile_k; 
  /**
   *@brief  Container for density profile in configuration space */
  vector<real_prec> v_density_profile_r; 
  /**
   *@brief  Container for galaxy P(k) 1-halo term satelite-satellite */
  vector<gsl_real> v_galaxy_power_spectrum_1h_ss; 
  /**
   *@brief  Container for galaxy P(k) 1-halo term satelite-central */
  vector<gsl_real> v_galaxy_power_spectrum_1h_sc; 
  /**
   *@brief  Container for galaxy P(k) 2-halo term */
  vector<gsl_real> v_galaxy_power_spectrum_2h;

  bool use_K_correction;
  bool use_e_correction;
  // Vector used to speed up cosmo-calculations by means of interpolation
  // Must come in gsl_real precision
  /**
   *@brief  Container for comoving distance r(z) */
  vector<gsl_real>rv;  
  /**0
   *@brief  Container for growth factor g(z) */
  vector<gsl_real>gv; 
  /**
   *@brief  Container for  redshift z */
  vector<gsl_real>zv; 
  /**
   *@brief  Container for transversal separation distance tr(z) */
  vector<gsl_real>trv; 
  /**
   *@brief  Selection for HOD model*/
  int hod_model;  
  /**
   *@brief  HOD parameter*/
  real_prec mmin_hod;  
  /**
   *@brief  HOD parameter*/
  real_prec ms_hod;
  /**
   *@brief  HOD parameter*/
  real_prec alpha_hod; 
  /**
   *@brief  HOD parameter*/
  real_prec alpha_A;
  /**
   *@brief  HOD parameter*/
  real_prec alpha_B;
  /**
   *@brief  HOD parameter*/
  real_prec alpha_C;
  /**
   *@brief  HOD parameter*/
  real_prec scatter_hod;
  /**
   *@brief  HOD parameter*/
  real_prec width_hod;
  /**
   *@brief  HOD parameter*/
  real_prec s_faint_hod;
  /**
   *@brief  HOD parameter*/
  real_prec s_bright_hod;
  /**
   *@brief  HOD parameter*/
  real_prec Mstep_hod;
  /**
   *@brief  HOD parameter*/
  real_prec mt_hod;
  /**
   *@brief  HOD parameter*/
  real_prec muno_hod;
  /**
   *@brief  HOD parameter*/
  real_prec coef_concentration;
  /**
   *@brief  HOD parameter*/
  real_prec coef_concentration_amp;
  /**
   *@brief  Container for comoving distance r(z) */
  vector<gsl_real>kvector_external;
  /**0
   *@brief  Container for growth factor g(z) */
  vector<gsl_real>power_external;
  bool use_external_power;
};
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @struct<matrices>
 * @brief The matrices struct
 * @details Structure with 2d cointainers with c onariance matrices, used in FB analysis
 */
struct matrices{
  /**
   *@brief  Measure
   */
  vector<real_prec>Cmed;
  /**
   *@brief  Measure fluctuation
   */
  vector<real_prec>Dmed;
  /**
   *@brief  Matrix R
   */
  vector<vector<real_prec> >R;
  /**
   *@brief  Matrix V
   */
  vector<vector<real_prec> >V;
  /**
   *@brief  Noise matrix
   */
  vector<vector<real_prec> >N;
  /**
   *@brief  Inverse of covariance matrix
   */
  vector<vector<real_prec> >iCov;
  /**
   *@brief  Covariance matrix
   */
  vector<vector<real_prec> >nCov;
  /**
   *@brief Something else
   */
  vector<vector<real_prec> >Step;
  real_prec det_matrix;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<params_clth>
 * @brief The params_clth struct
 * @details  Auxiliary strucure for the class:ClFunctions computation of theoretical angular power spectrum
 */
struct params_clth{
  /**
   *@brief   comoving disance
   */
  real_prec r;
  /**
   *@brief  Wave number*/
  real_prec k;
  /**
   *@brief  Angular multipole
   */
  int l;
  string wtype;
  /**
   *@brief  */
  real_prec zmax_bin;
  /**
   *@brief  */
  real_prec zmin_bin;
  /**
   *@brief  */
  real_prec rmax_bin;
  /**
   *@brief  */
  real_prec rmin_bin;
  /**
   *@brief  */
  real_prec k_min_integration;
  /**
   *@brief  */
  real_prec k_max_integration;
  /**
   *@brief  */
  real_prec sigma_p_errors; 
  /**
   *@brief  */
  string pdf_zerrors; 
  /**
   *@brief  */
  real_prec zaux; 
  /**
   *@brief  */
  vector<real_prec>pk; 
  /**
   *@brief  */
  vector<real_prec>kv; 
  /**
   *@brief  */
  vector<real_prec>Fkernel; 
  /**
   *@brief  */
  vector<real_prec>rv;/** *@brief  */
  /**
   *@brief  */
  vector<real_prec>zv; 
  /**
   *@brief  */
  vector<real_prec>gv; 
  /**
   *@brief  */
  vector<real_prec>bias_zv; 
  /**
   *@brief  */
  vector<real_prec>gfv; 
  /**
   *@brief  */
  vector<real_prec>Hv; 
  /**
   *@brief  */
  vector<real_prec>dn_photo; 
  /**
   *@brief  */
  vector<real_prec>dn_spect; 
  /**
   *@brief  */
  vector<real_prec>dz; 
  /**
   *@brief  */
  vector<real_prec>klnv; 
  /**
   *@brief  */
  vector<real_prec>sigma_mass; 
  /**
   *@brief  */
  vector<real_prec>M_nl; 
  /**
   *@brief  */
  vector<real_prec>HF_i4; 
  /**
   *@brief  */
  vector<real_prec>HF_i2; 
  /**
   *@brief  */
  vector< vector<gsl_real> > sBessel; 
  /**
   *@brief  */
  vector<gsl_real> sxBessel; 
  /**
   *@brief  */
  string redshift; 
  /**
   *@brief  */
  real_prec pk_normalization; 
  /**
   *@brief  */
  real_prec spectral_index;  
  /**
   *@brief  */
  bool use_non_linear_pk; 
  /**
   *@brief  */
  string type_of_nl_power; 
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_aux>
 @brief Template structure s_aux
 @details This structure is created in order to avoid missconfusing of information among threads when using OMP
*/
template<typename T>
struct s_aux{
  /**
   *@brief  */
  s_CosmologicalParameters *scp_a; /** *@brief  */
  /**
   *@brief  */
  params_clth *s_clth; 
  /**
   *@brief  */
  real_prec raux; 
  /**
   *@brief  */
  int laux; 
  /**
   *@brief  */
  real_prec kaux; 
  /**
   *@brief  */
  real_prec zaux; 
  /**
   *@brief  */
  vector<gsl_real>v_aux; 
  /**
   *@brief  */
  vector<gsl_real>XX; 
  /**
   *@brief  */
  vector<gsl_real>WW; 
  /**
   *@brief  */
  vector<gsl_real>XX_mu; 
  /**
   *@brief  */
  vector<gsl_real>WW_mu; 
  /**
   *@brief  */
  vector<gsl_real>XX_z; 
  /**
   *@brief  */
  vector<gsl_real>WW_z; 
  T Ps; 
  /**
   *@brief  */
  real_prec kmin_int;
  /**
   *@brief  */
  real_prec kmax_int;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<A1>
 * @brief The A1 struct
 * @details Structure used in the passage from Cl to Ps to spped up calculations
 */
struct A1{
  /**
   *@brief  */
  s_CosmologicalParameters *s_cp; 
  /**
   *@brief  */
  vector<gsl_real>MASS; 
  /**
   *@brief  */
  vector<gsl_real>MASS_FUNCTION; 
  /**
   *@brief  */
  vector<gsl_real>MASS_BIAS; 
  /**
   *@brief  */
  real_prec aux_k; 
  /**
   *@brief  */
  real_prec aux_z; 
  /**
   *@brief  */
  real_prec aux_m; 

};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_Halo>
 * @brief The s_Halo struct
 * @details The full catalog is read and allocated in a container of this structure-type. <br>
 * The FileOutput::read_file() method reads the catalog in a container props[] with lenght NBOBS*NCOLS,
 * The method Catalog::read_catalog() assigns the information of the column i_x (e.g., X-coords) as:
 *@code
 vector<s_Halo> halo(NOBJS);
 for(int i=0;i<NOBJECTS;++i)
     halo[i].coord1 = prop[i_x+i*NCOLS]; // x-coordinate of a catalog with NOBJECTS
 *@endcode
 */
struct s_Halo
{
  // -----------------------------------
  /**
   *@brief  Coordinate 1 (X,r or z) depending on params._sys_of_coords_g/r()  */
  real_prec coord1;
  // -----------------------------------
  /**
   *@brief   Coordinate 2 (Y,phi or phi) params._sys_of_coords_g/r() */
  real_prec coord2;
  // -----------------------------------
  /**
   *@brief   Coordinate 1 (Z or theta, ) params._sys_of_coords_g/r() */
  real_prec coord3;
  // -----------------------------------
  /**
   *@brief   Redshift (Z or theta, ) params._sys_of_coords_g/r() */
  real_prec cosmological_redshift;
  // -----------------------------------
  /**
   *@brief   Redshift (Z or theta, ) params._sys_of_coords_g/r() */
  real_prec redshift;
  // -----------------------------------
  /**
   *@brief  Velocity component in  coord1  */
  real_prec vel1;
  // -----------------------------------
  /**
   *@brief   Velocity component in  coord2 */
  real_prec vel2;
  // -----------------------------------
  /**
   *@brief   Velocity component in  coord3 */
  real_prec vel3;
  // -----------------------------------
  /**
   *@brief  Tracer Mass, idally in Ms/h  */
  real_prec mass;
  // -----------------------------------
  /**
   *@brief  Mass Bin in which the current tracer mass is located  */
  int mass_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Mass of tracer rassigned,  used for mocks  */
  real_prec mass_assigned;
  // -----------------------------------
  /**
   *@brief  Mass of refernece tracer, used for mocks  */
  real_prec mass_parent;
  // -----------------------------------
  /**
   *@brief  Maximum circular velocity   */
  real_prec vmax;
  /**
   *@brief  Bin of Maximum circular velocity  */
  int vmax_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer Vrms   */
  real_prec vrms;
  int vrms_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer Mach number   */
  real_prec mach_number;
  // -----------------------------------
  /**
   *@brief  Bin of tracer Mach number
  */
  int mach_number_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer Vmax   */
  real_prec vmax_assigned;
  // -----------------------------------
  /**
   *@brief  Tracer Vmax   */
  real_prec vmax_parent;
  // -----------------------------------
  /**
   *@brief  Tracer poistion in the set 0..Nobjects ordererd by increasing order in Vmax
   */
  ULONG vmax_index;
  // -----------------------------------
  /**
   *@brief  Effective bias computed object by object
   */
  real_prec  bias;
  // -----------------------------------
  /**
   *@brief  Quadratic bias computed object by object
   */
  real_prec  qbias;
  // -----------------------------------
  /**
   *@brief  Relative bias computed object by object
   */
  real_prec  relative_bias;
  // -----------------------------------
  /**
   *@brief  Effective bias computed object by object with kmax_2
   */
  real_prec  bias2;
  // -----------------------------------
  /**
   *@brief  Effective bias computed object by object with kmax3
   */
  real_prec  bias3;
  // -----------------------------------
  /**
   *@brief  bias computed from a box in redshift space
   */
  real_prec  bias_rs;
  // -----------------------------------
  /**
   *@brief  Parfameter linekd to the individual bias in redshift space
   */
  real_prec  rs_factor;
  // -----------------------------------
  /**
   *@brief  Number of neighbours in a sphere of radius R
   */
  real_prec  number_neigh;
  // -----------------------------------
  /**
   *@brief Mean separation to local neihbours in a sphere of radius R
   */
  real_prec  mean_sep_local_neigh;
  // -----------------------------------
  /**
   *@brief  Tracer Scale radius Rs   */
  real_prec rs;
  // -----------------------------------
  /**
   *@brief  Tracer Virial radius   */
  real_prec rvir;
  // -----------------------------------
  /**
   *@brief  Tracer Concentration   */
  real_prec concentration;
  // -----------------------------------
  /**
   *@brief Bin of Rs  */
  int rs_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer virial ratio T/|W|   */
  real_prec virial;
  // -----------------------------------
  /**
   *@brief  BIn of virial */
  int virial_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer Spin (Peebles)  */
  real_prec spin;
  // -----------------------------------
  /**
   *@brief  Tracer Spin defined by Bullock  */
  real_prec spin_bullock;
  // -----------------------------------
  /**
   *@brief Bin in spin
 */
  int spin_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer b_to_a ratio   */
  real_prec b_to_a;
  // -----------------------------------
  /**
   *@brief Bin in b_to_a
*/
  int btoa_bin_tracer;
  // -----------------------------------
  /**
   *@brief  Tracer b_to_a ratio   */
  real_prec c_to_a;
  // -----------------------------------
  /**
   *@brief  Bin in c_to_a ratio   */
  int ctoa_bin_tracer;

  /**
   *@brief  mean number density at current tracer mass  */
  real_prec mean_number_density;

  // -----------------------------------
  /**
   *@brief  Color (for galaxies)   */
  real_prec color;
  // -----------------------------------
  /**
   *@brief  Stellar mass (for galaxies)  */
  real_prec stellar_mass;
  // -----------------------------------
  /**
   *@brief  Apparent manitude   */
  real_prec app_mag;
  // -----------------------------------
  /**
   *@brief  Absolute manitude   */
  real_prec abs_mag;
  // -----------------------------------
  /**
   *@brief   ID (in the mesh) where this tracer is located given the nomial resolution () */
  ULONG GridID;
  // -----------------------------------
  /**
   *@brief   ID (in the mesh) where this tracer is located given the nomial resolution () */
  ULONG GridID_n;
  // -----------------------------------
  /**
   *@brief   bin number in the Theta properties binning in which the GridID of the galaxy has fallen */
  ULONG galThetaID;
  // -----------------------------------
  /**
   *@brief   ID (in the file) of each tracer (..Nobjects)  */
  ULONG galID;
  // -----------------------------------
  /**
   // *@brief   ID wghich identifiesthe multiscale levels in which the tracer acqured the property  */
  int multi_scale_level;
  // -----------------------------------
#ifdef _USE_MULTISCALE_
  /**
   *@brief  ID (in the mesh) where this tracer is located given the resolution with l1  */
  ULONG GridID_l1;
#endif
  // -----------------------------------
#ifdef _USE_MULTISCALE_
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l2  */
  ULONG GridID_l2;
#endif
  // -----------------------------------
#ifdef _USE_MULTISCALE_
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l3  */
  ULONG GridID_l3;
#endif
  // -----------------------------------
#ifdef _USE_MULTISCALE_
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l4  */
  ULONG GridID_l4;
#endif
  // -----------------------------------
  /**
   *@brief   Number of substructures (if the tracer is parent halo) */
  int number_sub_structures;
  // -----------------------------------
  /**
   *@brief  Identify whether this is a random object (-1) or a real tracer from DM particle sin LPT (1)  */
  int identity;
  // -----------------------------------
  /**
   *@brief  Identify whether this has been observed or avoided (used in different processes inside the code)  */
  bool observed;
  // -----------------------------------
  /**
   *@brief  Identify whether this has been observed or avoided (used in different processes inside the code)  */
  bool observed_secondary;
  // -----------------------------------
  /**
   *@brief   Weight (if any)
*/
  real_prec weight1;
  // -----------------------------------
  /**
   *@brief   Weight (if any) */
  real_prec weight2;
  // -----------------------------------
  /**
   *@brief   Weight (if any) */
  real_prec weight3;
  // -----------------------------------
  /**
   *@brief   Weight (if any) */
  real_prec weight4;
  // -----------------------------------
  /**
   *@brief   Mean density (uf any) evaluated at the position of the tracer */
  real_prec mean_density;
  // -----------------------------------
  /**
   *@brief   */
  ULONG PropID;
  // -----------------------------------
  /**
     @brief  COsmic-web clasification of tracer */
  int gal_cwt;
  // -----------------------------------
  /**
     @brief Interpolation of the DM in the cell where the halo is found */
  real_prec local_dm;
  // -----------------------------------
  /**
     @brief Eigenvalues of DM at the cellwherre tracers are found */
  real_prec lambda1;
  // -----------------------------------
  /**
     @brief Eigenvalues of DM at the cellwherre tracers are found */
  real_prec lambda2;
  // -----------------------------------
  /**
     @brief Eigenvalues of DM at the cellwherre tracers are found */
  real_prec lambda3;
  // -----------------------------------
  /**
     @brief  */
  real_prec reference;
  // -----------------------------------
  /**
     @brief  Measuremet of the number of tracers around a tracer in asphere of R set in parameter file.*/
  real_prec local_overdensity;
  // -----------------------------------
  /**
     @brief  Analogous to mach number, for distances to tracer: mean/var of the separations around the tracer in a sphere on a scale given in parameter file */
  real_prec dach_number;
  /**
     @brief  Analogous to mach number, for distances to tracer: mean/var of the separations around the tracer in a sphere on a scale given in parameter file */
  real_prec number_of_neighbours;
  // -----------------------------------
  /**
     @brief  Tidal anisotropy evaluated at the tracer position */
  real_prec tidal_anisotropy;
  // -----------------------------------
  /**
     @brief  Peak height delta/sigma computed for each tracer based on their virial mass*/
  real_prec peak_height;
  // -----------------------------------
    /**
  @brief  Spherical overdensity associated to the mass and radii of the tracer
    */
  real_prec spherical_overdensity;
  // -----------------------------------
  /**
    @brief  Mass of the closest neighbour
  */
  real_prec mass_closest_neighbour;
  /**
    @brief  Distance to the closest neighbour
  */
  real_prec distance_closest_neighbour;
  /**
    @brief  Mass of the most massive neighbour
  */
  real_prec most_massive_neighbour;
  /**
    @brief  Distnace ot the most massive neighbour
  */
  real_prec distance_to_most_massive_neighbour;
  /**
    @brief  Spin of the clisest neighbour
  */
  real_prec spin_closest_neighbour;
  /**
    @brief  Concentration of the closet neighbour
  */
  real_prec concentration_closest_neighbour;
  /**
    @brief  ID of the closest neighbour
  */
  ULONG id_closest_neighbour;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<experiments>
 * @brief The experiments struct
 * @details Containers used in the MCMC analysis of Cl. Used in class::HGAP
 */
struct experiments{
  vector<vector<real_prec> > acc_par1;
  vector<real_prec>  weight_par1;
  vector<vector<real_prec> > acc_par2;
  vector<real_prec>  weight_par2;
  vector<vector<real_prec> > acc_par3;
  vector<real_prec>  weight_par3;
  vector<vector<real_prec> > acc_par4;
  vector<real_prec>  weight_par4;
  vector<vector<real_prec> > acc_par5;
  vector<real_prec>  weight_par5;
  vector<vector<real_prec> > acc_par6;
  vector<real_prec>  weight_par6;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_data_structure>
 * @brief The s_data_structure struct
 * #details Auxiliary structure used in the measrement of power spectrum class::PowerSpectrum
 */
struct s_data_structure{
  /**
   *@brief  Mean density of catalog*/
  real_prec mean_density;
  /**
   *@brief  Ask whether the redshift distribution has to be computed from the input catalog */
  bool compute_dndz;
  /**
   *@brief  Container for redshift used if nbar hast to be computed here*/
  vector<gsl_real> zz_v;
  /**
   *@brief  Container for redshift distribution*/
  vector<gsl_real> dndz_v;
  /**
   *@brief  2D-Container for redshift distribution in different HEALPIX pixels*/
  vector<gsl_real> zz_c;
  vector<gsl_real> rr_c;
};
////////////////////////////////////////////////////////////////////////////
//Structure containing the properties of the catalogs
/**
 *@brief
 * @struct<s_data_structure>
 * @brief The s_data_structure struct
 * #details Auxiliary structure used in the measrement of power spectrum class::PowerSpectrum
 */
struct s_data_structure_direct_sum{
  /**
   *@brief  Mean density of catalog*/
  vector<s_Halo>properties;
  /**
   *@brief  Mean density of catalog*/
  int n_columns;
  /**
   *@brief  Mean density of catalog*/
  string catalog;
  /**
   *@brief  Mean density of catalog*/
  real_prec mean_density;
  /**
   *@brief  Ask whether the redshift distribution has to be computed from the input catalog */
  bool compute_dndz;
  /**
   *@brief  Container for redshift used if nbar hast to be computed here*/
  vector<gsl_real> zz_v;
  /**
   *@brief  Container for redshift distribution*/
  vector<gsl_real> dndz_v;
  /**
   *@brief  Mean density of catalog*/
  vector<gsl_real> zz_c;
  /**
   *@brief  Mean density of catalog*/
  vector<gsl_real> rr_c;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_parameters_estimator>
 * @brief The s_parameters_estimator structure
 @details Auxiliary structure used in the measurement of power spectrum, fed in class::PowerSpectrum, used in class::FftwFunctions
*/
struct s_parameters_estimator{
  /**
   *@brief  Use random catalog
   */
  bool use_random_catalog;
  /**
   *@brief  Number of selected real objects from which the P(k) will be measured
   */
  ULONG number_of_objects;
  /**
   *@brief  number of selected random objects from which the P(k) will be measured
   */
  ULONG number_of_randoms;
  /**
   *@brief  weighted number of selected real objects from which the P(k) will be measured
   */
  real_prec w_number_of_objects;
  /**
   *@brief weighted number of selected random objects from which the P(k) will be measured
   */
  real_prec w_number_of_randoms;  
  /**
   *@brief   w_number_of_objects divided by w_number_of_randoms
   */
  real_prec alpha;             
  /**
   *@brief  Sum of the squared of the weights in the real cataloge
   */
  real_prec S_g;
  /**
   *@brief  Sum of the squared of the weights in the random cataloge
   */
  real_prec S_r;               
  /**
   *@brief  Noramlization of power spectrum
   */
  real_prec normalization;
  /**
   *@brief  Normalization of window function
   */
  real_prec normalization_window;
  /**
   *@brief  Shot noise
   */
  real_prec shot_noise;         
  /**
   *@brief  Shot noise for window
   */
  real_prec shot_noise_window;  
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_parameters_box>
 * @brief The s_parameters_box struct
 */
struct s_parameters_box{
  //  string mas;       //Mass Assignment scheme
  //  string ave;        //Type of average
  vector<gsl_real>zz;
  vector<gsl_real>rc;
  long npixels;
  long nside;
  string file_dndz;
};
////////////////////////////////////////////////////////////////////////////
/**
 * @struct<s_params_box_mas>
 * @brief Auxiliary structure used to pass box arguments to different functions
 */
struct s_params_box_mas{
    /**
     * @brief mass kernel: 0 for NGP, 1 for CIC, 2 for TSC.
     */
  int masskernel;
  /**
   * @brief Size of the mesh in each direction
   */
  ULONG Nft;
  /**
   * @brief Size of a possible high resolution mesh in each direction
   */
  ULONG Nft_HR;
  /**
   * @brief Size of the mesh, \f$ N_{grid}=N_{ft}^{3}\f$
   */
  ULONG NGRID;
  /**
   * @brief Lenght of the side box, in Mpch /h
   */
  real_prec Lbox;
  /**
   * @brief Spatial resolution in the x-dir, \f$ L_{box}/N_{ft}\f$
   */
  real_prec d1;
  /**
   * @brief Spatial resolution in the y-dir, \f$ L_{box}/N_{ft}\f$
   */
  real_prec d2;
  /**
   * @brief Spatial resolution in the z-dir, \f$ L_{box}/N_{ft}\f$
   */
  real_prec d3;
  /**
   * @brief Spatial (high) resolution in the x-dir, \f$ L_{box}/N_{ft-HR}\f$
   */
  real_prec d1_HR;
  /**
   * @brief Spatial (high) resolution in the y-dir, \f$ L_{box}/N_{ft-HR}\f$
   */
  real_prec d2_HR;
  /**
   * @brief Spatial (high) resolution in the z-dir, \f$ L_{box}/N_{ft-HR}\f$
   */
  real_prec d3_HR;
  /**
   * @brief Minumum in x-dir
   */
  real_prec min1;
  /**
   * @brief Minumum in y-dir
   */
  real_prec min2;
  /**
   * @brief Minumum in z-dir
   */
  real_prec min3;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_minimums>
 * @brief The s_minimums struct
 * @details Auxiliary structure to allocate minimum values of properties used to characterize the bias in class::Bam
 */
struct s_minimums{
  real_prec prop0;  // Tracer property, counts
  real_prec prop0_mass;  // Tracer property, mass
  real_prec prop0_sf;  // Tracer property, satellite_fraction
  real_prec prop1;   // ð
  real_prec prop2;   // CWC
  real_prec prop3;   // Mk
  real_prec prop4;  // Inv1 or ð²
  real_prec prop5;  //InvII ot ð³
  real_prec prop6;  //Tidal Anis or s²
  real_prec prop7;  //InvI Shear or Nabla²ð
  real_prec prop8;  //InvII shear or s²ð
  real_prec prop9;  //InvIII Shear or s³
  real_prec prop10;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_maximums>
 * @brief The s_maximums struct
 * @details Auxiliary structure to allocate  maximum values of properties used to characterize the bias in class::Bam
 */
struct s_maximums{
  real_prec prop0;
  real_prec prop0_mass;  // Tracer property, mass
  real_prec prop0_sf;  // Tracer property, satellite_fraction
  real_prec prop1;
  real_prec prop2;
  real_prec prop3;
  real_prec prop4;
  real_prec prop5;
  real_prec prop6;
  real_prec prop7;
  real_prec prop8;
  real_prec prop9;
  real_prec prop10;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_Deltas>
 * @details Auxiliary structure to allocate  bin sizes in the different properties used to characterize the bias in class::Bam
 */
struct s_Deltas{
  /**
   *@brief   Bin size for property0 */
  real_prec prop0;
  /**
   *@brief   Bin size for property0 using tracer mass */
  real_prec prop0_mass;  // Tracer property, mass
  /**
   *@brief   Bin size for property0 satellite fraction */
  real_prec prop0_sf;  // Tracer property, satellite_fraction
  /**
   *@brief   Bin size for property1 */
  real_prec prop1;
  /**
   *@brief   Bin size for property2 */
  real_prec prop2;
  /**
   *@brief   Bin size for property3 */
  real_prec prop3;
  /**
   *@brief   Bin size for property4 */
  real_prec prop4;
  /**
   *@brief   Bin size for property5 */
  real_prec prop5;
  /**
   *@brief   Bin size for property6 */
  real_prec prop6;
  /**
   *@brief   Bin size for property7 */
  real_prec prop7;
  /**
   *@brief   Bin size for property8 */
  real_prec prop8;
  /**
   *@brief   Bin size for property9 */
  real_prec prop9;
  /**
   *@brief   Bin size for property10 */
  real_prec prop10;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_mass_members>
 * @details Strucure to allocate set of masses, coordinates and their ID (in the mesh) for tracers within a given bin of the multidimendional Theta (DM) properties
 * @details Used in the assignment of halo properties.
 * @details Example: Private variable in class::Bam
 * @code
 vector<s_mass_members> dm_properties_bins;
 *@endcode
 */
struct s_mass_members{
  /**
   *@brief  Masses (or any other property) in that particular theta_bin
   */
  vector<real_prec> tracer_properties;
  /**
   * @brief  Flag to pinpoint properties used
   **/
  vector<bool> used_property;
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
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_nearest_cells struct>.
 @details This is used in NumericalMethods::get_neighbour_cells() 
 @details Structure containing the closest cell neighbours to a given cell
 @code
 vector<s_nearest_cells>nearest_cells_to_cell(NGRID);
 @endcode
*/
struct s_nearest_cells{
  /**
   *@brief  Container with dimension Number_of_neighbours cointaining the ID of neighbour cells*/
  vector<ULONG> close_cell;
  /**
   *@brief  Tracks
   *info of boundary conditions in the x-direction */
  vector<int> bc_x;
  /**
   *@brief  Tracks info of boundary conditions in the y-direction */
  vector<int> bc_y;
  /**
   *@brief  Tracks info of boundary conditions in the z-direction */
  vector<int> bc_z;
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_cell_info struct>
 @details Structure for each cell of the mesh, containing the coordinates, mass, other properties, tracer-index and the position of the cell in the Theta-histgrams.
 Used in class::Catalog
 @code
 vector<s_cell_info> cell_info_tr(NGRID);
 for(ULONG i=0;i<halo[i].size() ;++i)
 {
     ULONG ID=halo[i].GridID;
     cell_info_tr[ID].posx_p.push_back(halo[i].coord1);  // and so on
 }
 @endcode
 @details The size of the containers of this structure is the number of tracers in each cell
 @code
 vector<int>Ncounts(NGRID,0);
 for(ULONG i=0;i<NGRID ;++i)
    Ncounts[i] = cell_info_tr[i].posx.size();
 @endcode
*/
struct s_cell_info{
  /**
   *@brief  Container to allocate the coordinates1 of objects in a cell*/
  vector<real_prec> posx_p;
  /**
   *@brief  Container to allocate the x-vels  of objects in a cell*/
  vector<real_prec> velx_p;
  /**
   *@brief  Container to allocate the coordinates2 of objects in a cell*/
  vector<real_prec> posy_p;
  /**
   *@brief  Container to allocate the y-vels of objects in a cell*/
  vector<real_prec> vely_p;
  /**
   *@brief  Container to allocate the coordinates3 of objects in a cell*/
  vector<real_prec> posz_p;
  /**
   *@brief  Container to allocate the z-vels of objects in a cell*/
  vector<real_prec> velz_p;
  /**
   *@brief  Container to allocate the mass of objects in a cell*/
  vector<real_prec> mass;
  /**
   *@brief  Container to allocate an extra-property of objects in a cell*/
  vector<real_prec> property;
  /**
   *@brief  Container to allocate the id of each object in a cell*/
  vector<ULONG> gal_index;
  /**
   *@brief  Position in the multidimensional DM property Theta bin in which the cell has been identified */
  ULONG Theta_bin;
  /**
   *@brief  ID of the current cell*/
  ULONG cell_index;
  /**
   *@brief  CWC current cell*/
  int cwt;
  /**
   *@brief  Local DM */
  real_prec local_dm;
};
////////////////////////////////////////////////////////////////////////////
/**
 * @struct<s_cell_info_reduced>
 * @brief Structure desingned to contain the information of the ID of the tracers living in a cell
*/
struct s_cell_info_reduced{
  /**
   *@brief  Container to allocate the ID of each tracer in a given cell
 */
  vector<ULONG> gal_index;
  /**
   *@brief  Clear container
 */
  void clear_mem(){gal_index.clear();gal_index.shrink_to_fit();}
};
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @struct<s_dist_in_dmbins struct>
 * @details Strucure to allocate the list of masses in pairs that -within one bin of the dm properties- are separated by a distnace between m_min and MAXIMUM_DISTANCE_EXCLUSION
*/
struct s_dist_in_dmbins{
  /**
   *@brief  Mass (or other prop) of one element of pair */
  vector<real_prec> M1;
  /**
   *@brief  Mass (or other prop) of second element of pair */
  vector<real_prec> M2;
  vector<ULONG> mass_to_swap;
};
////////////////////////////////////////////////////////////////////////////
/**
 * @struct<s_power_in_bins>
 * @brief Structure to allocate information of a number of measurements that come in differnt containers
 * @details Examples are : Power Spectrum (nmodes, k, l=0, l=2, l=4) or the primary and secondary bias.
 * @details For the power spetrum, one can define a container of thys type to allocate e.g. power in differnet mass bins
 * @details For primary bias, it is enough to define a single structure. FOr seconday bias
 * @it is advised ot define a vector of thys type-structure to allocate the bias in the different quartiles.
 * 
 * @author ABA
 */
struct s_info_in_bins{
/**
* @brief A name to identify the content, e,e, BIAS, POWER
*/
  string name_info;
/**
* @brief A name to identify the content, e,e, BIAS, POWER
*/
  string name_info_sec;
/**
* @brief Cointainer if integers: 
* @details Examples are: i) number of modes is k-shells, or counts in a bin in a histogram
*/
  vector<int> i_ncbin;
/**
* @brief Vector of real_prec type
* @details Mid values characterizing a binned property 
*/
  vector<real_prec> vbin;
/**
* @brief Vector of real_prec type 
* @details Excamples can be : i) Monopole of power, primary bias
*/
  vector<real_prec> vq1;
/**
* @brief Vector of real_prec type 
* @details Excamples can be : i) Quadrupole of power, primary bias
*/
  vector<real_prec> vq2;
/**
* @brief Vector of real_prec type 
* @details Excamples can be : i) Hexadecapole of power, primary bias
*/
  vector<real_prec> vq3;
/**
* @brief Function to allocate all members
*/
  void allocate_all(int Nbins)
  {
    this->i_ncbin.resize(Nbins,0);
    this->vbin.resize(Nbins,0);
    this->vq1.resize(Nbins,0);
    this->vq2.resize(Nbins,0);
    this->vq3.resize(Nbins,0);
  }
/**
* @brief Function to determine the size of a member, vbin by default
*/
  int s_size()
  {
    int ans=0;
    if(this->vbin.size()==0)
      std::cerr<<"Error allocating structure at "<<__LINE__<<std::endl;
    else
      ans= this->vbin.size();
    return ans;
  }
};
#endif
  
