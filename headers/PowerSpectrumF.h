//////////////////////////////////////////////////////////
/**
 * @class <PowerSpectrumF>
 * @brief Header file for the class PowerSpectrum
 * @file PowerSpectrumF.h
 * @author Andres Balaguera-Antolínez
 * @callgraph
 *
 * @details This file defines the interface of the class PowerSpectrum, used to obtain the measurements of Power spectrum (3D, Angular).
 */
//////////////////////////////////////////////////////////
#ifndef __POWERSPECTRUM__
#define __POWERSPECTRUM__
# include "Params.h"
# include "FileOutput.h"
# include "Catalog.h"
# include "CoordinateSystem.h"
# include "FftwFunctions.h"
# include "Cwclass.h"
# include "ScreenOutput.h"
# include "McmcFunctions.h"
# include "GnuplotC.h"
using namespace std;
class PowerSpectrumF{
 private :
    //////////////////////////////////////////////////////////
/**
  * @private
  * @brief Object of type GnuplotC
  * @note Used to make plots in Gnuplot
*/
  GnuplotC gp;
  //////////////////////////////////////////////////////////
  /**
   * @private
  * @brief Object of type GnuplotC
  * @note Used to make plots in Gnuplot
   */
  Gnuplot gp_pdf;
  //////////////////////////////////////////////////////////
  /**
   * @private
  * @brief Object of type GnuplotC
  * @note Used to make plots in Gnuplot
   */
  Gnuplot gp_power;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Catalog, meant to allcoate the tracer catalog
   */
  Catalog tracer_cat;
   //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Catalog, meant to allcoate the random catalog
   */
  Catalog random_cat;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Object of type params.
   */
  Params params;
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions
   */
  Cosmology cosmology; // nicer name
  
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions
   */
  McmcFunctions mcmc;
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type FftwFunctions
   */
  FftwFunctions fftw_functions;
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of galaxies in catalog
   * @details This is read from the class member tracer_cat
  */
  long N_galaxy;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of randoms in catalog
   * @details This is read from the class member random_cat
   */
  long N_random;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type cwclass
   */
  Cwclass  cwclass;
  //////////////////////////////////////////////////////////
  /**
   * @private
  *  @brief Vector containing the k identifying spherical shells for the estimates of gal sample
   */
  vector<real_prec> kvector_data;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector containing the k identifying spherical shells for the estimates of the window function
   */
  vector<real_prec> kvector_window;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Vector containing the k identifying spherical shells for the estimates of power spectrum* in 2d cart
   */
  vector<real_prec> kvector_data2d;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  vector<real_prec> muvector;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Complex field to allocate FT of dm
   */
   complex_prec *delta_dm;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Container for the power spectrum in real spece
   */
  vector<real_prec> pk;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Monopole of power spectrum in redshift space
   */
  vector<real_prec> pk0;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Quadrupole of power spectrum in redhoift space
   */
  vector<real_prec> pk2;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Hexadecapole of power in redhsift space
   */
  vector<real_prec> pk4;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Power Spectrum of window function
   */
  vector<real_prec> pk_w;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Number of Fourier modes in k-bins
   */
  vector<int> modes_g; 
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief  2D power spectrum kperp, kparallel
   */
  vector < vector<real_prec> > pkk;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief 2D power spectrum, |k|, mu
   */
  vector < vector<real_prec> > pmk;
    //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Variance of the FKP estimator
   */
  vector<real_prec> sigma_fkp;
    //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief
   */
  vector<real_prec> sigma_y_l2;
    //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief
   */
  vector<real_prec> sigma_y_l4;
 //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief
   */
  vector<real_prec> kvector_data_b;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief
   */
  vector<real_prec> bispectrum;
 //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  vector<real_prec> sn_bispectrum;
   //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  vector<int> modes_tri;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string rest_file;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  ULONG nside;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_random; 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_dndz;   
  //////////////////////////////////////////////////////////
    /**
   *  @brief 
   */
  string file_power;  
  //////////////////////////////////////////////////////////
    /**
   *  @brief
   */
  string file_power_cross;
 //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_power_real_space;  
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_power_redshift_space;  
  //////////////////////////////////////////////////////////
    /**
   *  @brief
   */
  string file_power_marked;
 //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  string file_power_marked_real_space;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  string file_power_marked_redshift_space;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_MCF;  
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_power_log;  
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_power2d; 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_power2d_mk;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_window; 
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string file_bispectrum;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  real_prec mean_density;
 //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  bool use_random_catalog_cl;
 //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  int Nbins_r;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  real_prec rmin;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  real_prec rmax;
  //////////////////////////////////////////////////////////
  /**
   *  @brief 
   */
  string r_bin_type;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
 public :

 /**
   *  @brief Default constructor
   *  @return object of PowerSpectrum
   */
  PowerSpectrumF ():N_random(0),mean_density(0),Nbins_r(0),rmin(0),rmax(0){
      time_t time_bam;
      time(&time_bam);
      this->So.initial_time=time_bam;
      this->So.welcome_message();
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /**
   *  @brief Overloaded constructor
   *  @param params 
   */
  PowerSpectrumF (Params &_params): params(_params),N_random(0),mean_density(0),Nbins_r(0),rmin(0),rmax(0)
  {
    this->fftw_functions.set_params(this->params);
    this->set_output_filenames();
    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;
    this->gp_pdf<<"set border linewidth 1.3\n";
    this->gp_pdf<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_pdf<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_pdf<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_pdf<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_pdf<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->tracer_cat.set_params(this->params);
    this->random_cat.set_params(this->params);
    this->cosmology.set_cosmo_pars(this->params.s_cosmo_pars);
    this->tracer_cat.set_params(this->params);
    this->random_cat.set_params(this->params);
  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  PowerSpectrumF (Params &_params, bool resize): params(_params),N_random(0),mean_density(0),Nbins_r(0),rmin(0),rmax(0)
  {
    this->fftw_functions.set_params(this->params);
    this->set_output_filenames();
    this->cosmology.set_cosmo_pars(this->params.s_cosmo_pars);
    if (resize)
    {
      this->kvector_data.clear();
      this->kvector_data.shrink_to_fit();
      for(int i=0;i<this->params._d_Nnp_data();i++)
        this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
    }
  }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /**
   *  @brief Default destructor
   *  @return none
   */
  ~PowerSpectrumF () {}
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined (_USE_MASS_CUTS_PK_)  || defined (_USE_ALL_PK_)
  /**
   * @brief Add both galaxy and random catalog at once
  */
void add_catalogues(real_prec mcuts);
void add_catalogues(int, int, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief Add both galaxy catalog, aimed at improving mmemory management
  */
  void add_catalogues_gal(real_prec mcuts);
  void add_catalogues_gal(int, int, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief Add both galaxy catalog, aimed at improving mmemory management
  */
  void add_catalogues_ran(real_prec mcuts);
  void add_catalogues_ran(int, int, int);
#elif defined (_USE_MASS_BINS_PK_)
  void add_catalogues(real_prec m_min, real_prec m_max);
  void add_catalogues_gal(real_prec m_min, real_prec m_max);
  void add_catalogues_ran(real_prec m_min, real_prec m_max);
#endif

#ifdef _USE_SEVERAL_RANDOM_FILES_
  void add_catalogues_ran(real_prec mcuts, int);
  void add_catalogues_ran(int, int, int, int);
#ifdef _USE_MASS_BINS_PK_
  void add_catalogues_ran(real_prec m_min, real_prec m_max, int);
#endif
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void add_random_catalogue();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  */
  string file_data;
//////////////////////////////////////////////////////////
  /**
   * @brief 
  */
  real_prec alpha;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  */
  real_prec normal_power;
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  real_prec var_prop;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  */
  vector<real_prec>window_matrix;
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure to allocate power (l=0,2,4), kvector and nmodes
  */
  vector<s_info_in_bins> power_in_bins;
   //////////////////////////////////////////////////////////
  /**
   * @brief Measurements of power spectrum
   * @details This function generates the estiametes of power spectrum ()
   *  (FKP, Yamamoto and their multipole decomposition) and the Bispectrum
   *  (using FKP). This method is aimed for realistic surveys in which a random file is expected to be provided
   *  along with (if defined in parameter file) a file with smoothed nbar tabulated
   * @param verbose: use dto write on screen
   * @param verbose: true or false to whether cuts in mass is to be done, to be deprecated.
   * @details This method reads input files and get power without touching the data, i-.e, does not make any binning in any property. It accepts catalogs as well as density fields.
   * @note This methood has been adapted to the case in which the random catalogs are too large to be kept in memmory or they come in chunks of randoms
   * with the SAME volume. The steps performed fr teh random catalog in this case are
   * 1 Read catalog.
   * 2 Get ready for nbar.
   * 3 Transform to cartessian coord (assigning nbar if needed and searching for box side lenght, offsets and mins).
   * 4 Interpolate on the mesh.
   * 5 Update params.
   * 6 Get window matrix if requested.
   * 7 Release memmory from the catalog.
   * If random file is too large, a loop over the chunks is done.
   * The a new box size is requested, the size of the box is computed from the randoms, and in particular, from the first chunk, such that the parameters in the PowerSpectrumF and FfftwFunctions are updatred with
   * @code
     this->set_params(this->random_cat.params);
     this->fftw_functions.set_params(this->random_cat.params);
   * @endcode
   * @warning No variance for Yammoto estimator have been coded.
   * @author ABA
*/
  void compute_power_spectrum (bool verbose, bool mcuts);
  //////////////////////////////////////////////////////////
  /**
   * @brief Measurements of power spectrum
   * @details This function generates (only) the estimates of power spectrum ()
   *  (FKP, Yamamoto and their multipole decomposition). This method is aimed for realistic surveys in which a random file is expected to be provided
   *  along with (if defined in parameter file) a file with smoothed nbar tabulated
   * @details This method reads input files and get power inside loops of bins of different properties specified in the parameter file
   * make any binning in any property.
   * @note This method is only applicabble for input catalogs (not for input inteprolated fields). See the description in this->compute_power_spectrum (bool, bool) to understand the steps
   * @warning This method contains explicit definitions that must be in agreement with the method Catalog::read_catalog().
   * Measurement of power spectrum of tracer catalog, with a random catalog and in bins of galaxy properties
   * The method is meant to be used whern bins in galaxy properties are to be done to mesaure power:  loops over bins, read data, select, compute.
   * This is desiged in this way as for large catalogs, these cannot be kept in the ram, otherwise one can do read bin compute, as it happens for the secondary bias analysi
   * Output file, keep roor in case intervals are used
   * Assumes that sys_of_coord is either I_EQZ or I_EQRZ
   * No window matrix omputation available
   * @author ABA
*/
  void compute_power_spectrum ();
  //////////////////////////////////////////////////////////
 /**
  * @brief Measrements of th window matrix
  * @details This function generates the elements of the window matrix of power spectrum, for differnet multipoles
  * and computed using the random catalog. Must be called inside the method compute_power_spectrum(bool, bool).
  * @note  This function is to be read from the cosmicatlas.cpp main function: verbose is to show detailed information, mcut is to force the analysis to have a minium cut in a given property
  depending on w. This method can be accesed from compilation with the -w option.

*/
  void compute_window_function ();
  //////////////////////////////////////////////////////////
  /**
   * @brief Measurements of power spectrum
   * @details This method has as input the tracer catalog in a box to measure the power once it has been read and allocated before in a Halo vector.
   * @details  This is meant for the case in which NO RANDOMS ARE USED AND IS NOT READY TO USE MASS CUTS
   * @details Assumed sys_of_coords = 0 and FKP
*/
  void compute_power_spectrum (bool verbose, vector<s_Halo>&tracer);
  //////////////////////////////////////////////////////////
  /**
   * @brief Measurements of power spectrum
   * @details This function is to be called from inside the code where catalogs are allocated in a s_halo type structure.
   * @param tracer: a s_Halo type container with the halo catalog.
   * @param space: "real space" or ·"redshift space"
   * @param property : "_MASS_", or "_VMAX_". No more properties of the catalog are allowed at this function
*/
  void compute_power_spectrum(vector<s_Halo>&tracer, string space, string property);
  //////////////////////////////////////////////////////////
  /**
   * @brief Analysis of secondary bias
   * @details This function computs lss bias, pearson, marked power and power of tracer in bins of differnet halo properties
   * The ordering of the properties is fixed, but the name of the ouput dir must emphasize which is the main property somehow.
   * This is the main method used for seconday bias studies.
   * @param space: real_space, redshift_space
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
*/
  void halo_bias_analysis(string space);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
  */
  void compute_marked_correlation_function ();
  //////////////////////////////////////////////////////////
  /**
   * @brief
  */
  void compute_power_spectrum_grid();
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Power_Spectrum and/or the Bispectrum
   * @details This function generates the estiametes of power spectrum ()
   *  (FKP, Yamamoto and their multipole decomposition) and the Bispectrum
   *  (using FKP).
   *  @note Used specifically in case we supply grids 1+delta1 and 1+delta2 and we want to compute the auto and cross power threof
   */
  void compute_power_spectrum_grid(const vector<real_prec>&, bool);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Power_Spectrum from a grid density
   * @details explicitely dedicatd to the case in which tracer bias is to be computed   */
  void compute_power_spectrum_grid(const vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the Marked Power_Spectrum
   */
  void compute_marked_power_spectrum_grid(const vector<real_prec>&, const vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the cross power spectrum between fields X and Y
    * @param dm bool set to true if X is dark matter so no shot-noise is computed for that component
 */
  void compute_cross_power_spectrum_grid(bool dm, vector<real_prec>&X,vector<real_prec>&Y, bool get_cross);
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute the cross power spectrum between fields X and Y
    * @param dm bool set to true if X is dark matter so no shot-noise is computed for that component
 */
  void compute_cross_power_spectrum_grid(bool dm, string file_X,string file_Y,bool get_cross);
  //////////////////////////////////////////////////////////
  /**
 * @brief Large-scale bias from power spectra
 *
 * @details
 * If the power spectra involved come from Gaussian initial conditions (IC),
 * we need to use the error bars from the Gaussian approximation on large scales.
 *
 * The bias is computed using:
 * \f[
 * \langle b \rangle = \frac{\sum \left( b_k / \sigma_{b_k}^2 \right)}{\sum \left( 1 / \sigma_{b_k}^2 \right)}
 * \f]
 * and the error is:
 * \f[
 * \sigma^2 = \frac{1}{\sum \left( 1 / \sigma_{b_k}^2 \right)}
 * \f]
 * where \f$ \sigma_{b_k} \f$ is the error on the bias,
 * and the bias is defined as:
 * \f[
 * b = \sqrt{\frac{P_h}{P_{\mathrm{dm}}}}
 * \f]
 *
 * The error on the bias is given by:
 * \f[
 * \sigma_{b_k}^2 = \frac{1}{4} b_k^2 \left( \frac{\sigma_{k,h}^2}{P_h^2} + \frac{\sigma_{k,\mathrm{dm}}^2}{P_{\mathrm{dm}}^2} \right)
 * \f]
 * (errors added in quadrature), where:
 * \f[
 * \sigma_k^2 = P^2 \cdot \frac{2}{N_k}
 * \f]
 * (neglecting shot noise), and \f$ N_k \f$ is the number of modes in each k-bin. This yields:
 * \f[
 * \langle b \rangle = \frac{\sum (N_k / b_k^2)}{\sum (N_k / b_k^2)}
 * \f]
 * and:
 * \f[
 * \sigma_{b_k}^2 = \frac{b_k^2}{N_k} = \frac{1}{\sum (N_k / b_k^2)}
 * \f]
 * If we use fixed-amplitude fields, the average bias is computed as a simple mean,
 * and its error follows from the standard deviation.
 *
 * @note Ojo
 */

  pair<real_prec, real_prec> get_lss_bias(vector<real_prec>&,vector<real_prec>&,vector<real_prec>&, vector<int>&, int init, real_prec kmax);
  //////////////////////////////////////////////////////////
  /**
   * @brief Large-scale bias from cross-power
   */
  pair<real_prec, real_prec> get_lss_cross_bias(vector<real_prec>&,vector<real_prec>& ,vector<real_prec>&, vector<int>&, int init, real_prec kmax);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @param write sigma true/false.
   */
  void write_power_spectrum (bool);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @param write sigma yes/no.
   */
  void write_power_and_modes();
  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @param write sigma yes/no.
   */
  void write_power_and_modes(string);
  //////////////////////////////////////////////////////////
  /** 
   * @brief Write output of power spectrum and/or bispectrum 
   * @details computed in the function compute_power_spectrum()
   */
  void write_power_spectrum_grid(string);
  //////////////////////////////////////////////////////////
 /**
   * @brief Set parameters
   */
   void set_parameters_power();
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the values of the private variables that regard filenames - with one catalogue
   */
  void set_output_filenames();
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::kvector_data
   *  @param i index of the vector kvector_data
   *  @return PowerSpectrum::kvector_data[i]
   */
  real_prec _kvector_data (int i) { return kvector_data[i]; }
 //////////////////////////////////////////////////////////
  /**
   *  @brief get the size of the private member
   *  PowerSpectrum::kvector_data
   *  @return PowerSpectrum::kvector_data.size()
   */
  ULONG _kvector_data_size () { return kvector_data.size(); }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk0
   *  @param i index of the vector pk0
   *  @return PowerSpectrum::pk0[i]
   */
  real_prec _pk0 (int i) { return pk0[i]; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::modes
   *  @param i index of the vector pk0
   *  @return PowerSpectrum::pk0[i]
   */
  real_prec _nmodes_k (int i) { return modes_g[i]; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk2
   *  @param i index of the vector pk2
   *  @return PowerSpectrum::pk2[i]
   */
  real_prec _pk2 (int i) { return pk2[i]; }
  ////////////////////////////////////////////////////////// 
  /**
   *  @brief get the value of the private member 
   *  PowerSpectrum::pk4
   *  @param i index of the vector pk4
   *  @return PowerSpectrum::pk4[i]
   */
  real_prec _pk4 (int i) { return pk4[i]; }
  //////////////////////////////////////////////////////////
  /**
   *  @brief get the value of the private member
   *  PowerSpectrum::pk4
   *  @param i index of the vector pk4
   *  @return PowerSpectrum::pk4[i]
   */
  real_prec _sigma_fkp (int i) { return sigma_fkp[i]; }
  //////////////////////////////////////////////////////////
  /**
   * @brief Computes the window matrix 
   *@details  This function provides the window matrix which in reality is a tensor of 4 indices:
   * two indices for multipoles l, l': these can be from 0 to 4, usually 0,2,4
   * two indices for wavenumber modes, k and k': k runs over the bins used in the measurements, 
   * while k' runs oer the modes implemented in the Gauss -legendre integration.
   * This computation involves a double rum over the random catalog.
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
  */
  void get_window_matrix_multipole();
  //////////////////////////////////////////////////////////
  /**
   * @brief Compute GRF
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void get_GaussianRandomField();
  //////////////////////////////////////////////////////////
  /**
   * @brief Method to assign effective bias to tracers
   * @param tracer_cat Catalog of dark matter tracers
   * @param dm_field Dark matter density field intepolated on a mesh
   * @details The individual bias is computed as 
   * \f[
   * b^{(i)}_{hm}=\frac{\sum_{j,k_{j}<k_{max}}N^{j}_{k}\langle {\rm e}^{-i\textbf{k} \cdot \textbf{r}_{i}} \delta_{\mathrm{dm}}^{*}(\textbf{k}) \rangle_{k_{j}}}{\sum_{j,k_{j}<k_{max}} N^{j}_{k}P_{\rm dm}(k_{j})},
   * \f]
   * where:  
   * - \f$ N_k^j \f$ is the number of Fourier modes in the \f$ j \f$-th spherical shell,
   * - \f$ \delta_{\mathrm{dm}}(\mathbf{k}) \f$ is the Fourier transform of the dark matter density field,
   * - \f$ P_{\mathrm{dm}}(k_j) \f$ is the dark matter power spectrum.
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void object_by_object_bias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field);
  //////////////////////////////////////////////////////////
  /**
  *  @brief Method to assign effective bias and relative bias to tracers
   * @param tracer_cat Catalog of dark matter tracers
   * @param dm_field Dark matter density field intepolated on a mesh
   * @param tracer_field tracer density field intepolated on a mesh
   * @details The power spectrum of the DM has been computed before the call of this method  and is still allocated in the container this->pk0
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void object_by_object_bias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field,vector<real_prec>& tracer_field);
  //////////////////////////////////////////////////////////
  /**
   * @brief Assignment of the object-by-object relative_bias.
   * @details Included in compute_power_spectrum meant to do that task in a realiztic sample. This method is meant to assign rbias for a realistic sample which imvolves the computation of the window fnuction * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void object_by_object_rbias();
  //////////////////////////////////////////////////////////
  /**
   * @brief Method to assign large scale bias as a function of harmonic index from Harmonic decomposition in Fourier space.
   * @details See details in this->object_by_object_bias
   */
  void object_by_object_qbias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field);
  //////////////////////////////////////////////////////////
  /**
   * @brief Method to assign effective bias to tracers
   * @param tracer_cat Catalog of dark matter tracers
   * @param dm_field Dark matter density field intepolated on a mesh
   * @param lmax_bias Maximum multipole
   * @details The individual bias is computed as 
   * \f[
   * b_{\ell}^{(i)}= \frac{1}{2\ell +1} \sum_{m=-\ell}^{m=\ell} \frac{\sum_{j,k_{j}<k_{max}}N^{j}_{k}\langle {\rm e}^{-i\textbf{k} \cdot \textbf{r}_{i}} \delta_{\mathrm{dm}}^{*}(\textbf{k})Y_{\ell m}(\hat{\textbf{k}}) \rangle_{k_{j}}}{\sum_{j,k_{j}<k_{max}} N^{j}_{k}P_{\rm dm}(k_{j})},
   * \f]
   * where:  
   * - \f$ N_k^j \f$ is the number of Fourier modes in the \f$ j \f$-th spherical shell,
   * - \f$ \delta_{\mathrm{dm}}(\mathbf{k}) \f$ is the Fourier transform of the dark matter density field,
   * - \f$ Y_{\ell m}(\hat{\mathbf{k}}) \f$ are spherical harmonics,
   * - \f$ P_{\mathrm{dm}}(k_j) \f$ is the dark matter power spectrum.
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void object_by_object_bias_lm(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field, int lmax);
  //////////////////////////////////////////////////////////
  /**
   *  @brief Compute cross correlatio.
   */
  void get_cross_correlation_config_space(vector<real_prec>&, vector<real_prec>& , vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief Method to do several taks lined to the secondary bias project
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void GetSuperClusters(string, string);
  //////////////////////////////////////////////////////////
  /**
   * @brief Method to do several taks lined to the secondary bias project
   * @warning This can be imlemented as an external function (e.g. in Miscelaneous.cpp) that implements methods of the current class.
   */
  void GetSuperClusters(string);
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set parameters
   */
  void set_params(Params _pars){this->params=_pars;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the number of galaxies
   */
  void set_N_galaxy(ULONG new_N_galaxy){this->N_galaxy=new_N_galaxy;}
  //////////////////////////////////////////////////////////
  /**
   *  @brief Set the number of randoms
   */
  void set_N_random(ULONG new_N_random){this->N_random=new_N_random;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Write output of power spectrum and/or bispectrum
   * @details computed in the function compute_power_spectrum()
   */
  void write_power_spectrum ();
  //////////////////////////////////////////////////////////
  /**
   * @brief Write log filw with numbers used and derived in the measuremetn of power
   */
  void write_fftw_parameters();


};

#endif
