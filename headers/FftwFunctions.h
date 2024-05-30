////////////////////////////////////////////////////////////////////////////
/**
 * @class <FftwFunctions>
 * @brief    Header file for the class FftwFunctions::
 * @file     FftwFunctions.h
 * @title    Manipulation of functions in Fourier space
 * @author   Andres Balaguera Antolinez
 * @author   Federico Marulli & Jennifer Pollack
 * @author   Optimization and parallelization by Luca Tornatore
 * @version  1.
 * @date     2013-2024
 * @details  See documentation
*/
////////////////////////////////////////////////////////////////////////////
#ifndef __FFTW_FUNCTIONS__
#define __FFTW_FUNCTIONS__
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include <vector>
# include <algorithm>
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <healpix_data_io.h>
# include <healpix_base.h>
# include <healpix_data_io.h>
# include <alm_powspec_tools.h>
# include <alm_healpix_tools.h>
# include <omp.h>
# include <complex>
# include <fftw3.h>
# include "fftw_array.h"
# include "NumericalMethods.h"
# include "Miscelanious.h"
# include "CosmologicalFunctions.h"
# include "FileOutput.h"
# include "CoordinateSystem.h"
# include "ScreenOutput.h"
# include "Params.h"
using namespace std;
//////////////////////////////////////////////////////////
class FftwFunctions{
 private:
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Obejct of type ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Obejct of type Cosmology
   */
  Cosmology CosmoF;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Minimum X-coordinate of the sample in Mpc/h
   */
  real_prec Xmin;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Maximim X-coordinate of the sample in Mpc/h
   */
  real_prec Xmax;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Minimum Y-coordinate of the sample in Mpc/h
  */
  real_prec Ymin;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Maximum Y-coordinate of the sample in Mpc/h
  */
  real_prec Ymax;
  //////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Minimum Z-coordinate of the sample in Mpc/h
   */
  real_prec Zmin;
  /////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief Maximum Z-coordinate of the sample in Mpc/h
   */
  real_prec Zmax;
 //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Vectors used in the estimator of the FKP variance of power spectrum
   */
  vector<real_prec> SN;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Vectors used in the estimator of the FKP variance of power spectrum
   */
  vector<real_prec> Q;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of random objects
   */
  ULONG n_ran;
 //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Weighted number of galaxies
   */
  real_prec w_g;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Weighted number of random objects
     */
  real_prec w_r;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Alpha parameter
   * @details The ratio between the weighted number of galaxies and the weighted number of randoms
   * @note This is computed as \f$ \alpha = \frac{\sum_{i=gal}w_{fkp}(r_i)}{\sum_{j=ran}w_{fkp}(r_j)}\f$
   */
  real_prec alpha;

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Parameter used to compute shot nnoise in power spectrum
   */
  real_prec s_g;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Parameter used to compute shot nnoise in power spectrum
   */
  real_prec s_r;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Parameter used to compute shot noise in bispectrum
   */
  real_prec sr1;
  //////////////////////////////////////////////////////////
  /**
   * @private
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sr2;
//////////////////////////////////////////////////////////
  /**
   * @private
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sg1;
//////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute shot noise in bispectrum
  */
  real_prec sg2;
  //////////////////////////////////////////////////////////
  /**
   * @brief Used in the normalization in bispectrum
   */
  real_prec normal_p;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Normalization of window function
   */
  real_prec normal_window;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Used in the normalization in bispectrum
   */
  real_prec normal_b;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Normalization in bispectrum
   */
  real_prec normal_bispectrum;
//////////////////////////////////////////////////////////
  /**
   * @private
   *  @brief  Normalization in power spectrum
   */
  real_prec normal;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Poisson Shot noise for window function
   */
  real_prec shot_noise_window;
//////////////////////////////////////////////////////////
  /**
   * @brief Poisson Shot noise for bispectrum
   */
  real_prec shot_noise_b1;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Poisson Shot noise for bispectrum
   */
  real_prec shot_noise_b2;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  real_prec correction_MAS_exp;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Number of grid cells in each direction for Bispectrum fast
   */
  int sgrid;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Total number of grid cells used in arrays defined for the Bispectrum fast
   */
  int new_sgrid;
  //////////////////////////////////////////////////////////
   /**
   * @private
     * @brief Offset in X direction
     */
  real_prec Xoffset;
  //////////////////////////////////////////////////////////
    /**
   * @private
     * @brief Offset in Y direction
     */
  real_prec Yoffset;

  //////////////////////////////////////////////////////////
    /**
   * @private
     * @brief Offset in Z direction
     */
  real_prec Zoffset;
//////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  inverse of the number of grid cells in each direction for the DFT
   */
  real_prec rNft;
//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xx;
//////////////////////////////////////////////////////////
  /**
   * @brief Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_yy;
//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_zz;
//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xy;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_xz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_g_yz;

//////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxx;

//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyy;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzz;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxy;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xxz;
//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyx;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yyz;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzx;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zzy;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xyy;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xzz;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yzz;
//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_xyz;
//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_yxz;
//////////////////////////////////////////////////////////
   /**
   * @brief  Vectors for the galaxy catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector<real_prec> data_g_zxy;
//////////////////////////////////////////////////////////
    /**
   * @brief  Vectors for the random catalogue used by the Yamamoto-Blake (fftw-based) estimator
   */
  vector <real_prec> data_r;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykx;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum. Contains the MAS correction for each mode
   */
  vector<real_prec> Array_corr;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arrayky;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykz;

  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> Arraykk;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<ULONG> VecArray;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   * Number of modes per shell for the Bispectrum
   */
  vector<int> Bmodes;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<int> kkminID;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<real_prec> kbins_bk;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   */
  vector<int> Ngrids_bk;
  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */
  real_prec rdeltax;
  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */

  real_prec rdeltay;
  //////////////////////////////////////////////////////////
  /**
   * @brief inverse of deltax/y/z, used in grid_assignment
   */
  real_prec  rdeltaz;

  //////////////////////////////////////////////////////////
  /**
   * @brief Number of galaxies
   */
  ULONG n_gal;
  //////////////////////////////////////////////////////////
  /**
   * @brief Product of NFT^2
   */
  int Nft2;
  //////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute the normalization in power spectrum
  */
  real_prec normal_power;
  //////////////////////////////////////////////////////////
  /**
      @brief Parameter used to compute the normalization in power spectrum
  */
  real_prec normal_power_two;

    //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */
  vector<ULONG> ArrayID;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */
  vector<vector<real_prec>> iFT_output_delta_full;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum
   */
  vector<vector<real_prec>> iFT_output_triangles_full;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vector used in Bispectrum.
   */
  vector<vector<real_prec>> iFT_shot_noise_p1_cyc_sum_full;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors used in the Bispectrum
   * Assigns to each vector in K-space the shell in which it's found
   */
  vector<ULONG> Kbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation
   */
  complex_prec *data_out_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation
   */
  complex_prec *data_out_g_rss;

  //////////////////////////////////////////////////////////
  /**
   * @brief Complex vector for the output of the FFTW for the fluctuation, used for the crossed power
   */
  complex_prec *data_out_gp;
  //////////////////////////////////////////////////////////
  /**
      @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xx;
  //////////////////////////////////////////////////////////
  /**
   *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yy;
  //////////////////////////////////////////////////////////
 /**
      @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zz;
  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xy;
  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xz;
  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yz;
  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xxx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyy;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzz;

  //////////////////////////////////////////////////////////
 /**
  *  @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xxy;

  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
   */
  complex_prec *data_out_g_xxz;

  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yyz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzx;

  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zzy;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xyy;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xzz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yzz;

  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_xyz;
  //////////////////////////////////////////////////////////
 /**
  *   @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_yxz;
  //////////////////////////////////////////////////////////
 /**
  *    @brief Complex vector for the output of the FFTW used by the Yamamoto-Blake estimator
  */
  complex_prec *data_out_g_zxy;
  //////////////////////////////////////////////////////////
  /**
   *   @brief Complex vector for the output of the FFTW
   */
  complex_prec *data_out_r;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Complex vector for the output of the FFTW used in the FKP estimation of the variance
   */
  complex_prec *data_out_SN;
  //////////////////////////////////////////////////////////
    /**
   *  @brief Complex vector for the output of the FFTW used in the FKP estimation of the variance
   */
  complex_prec *data_out_Q;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y0;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y2;
  //////////////////////////////////////////////////////////
  /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_g_out_y4;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y0;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y2;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector< complex<real_prec> > data_r_out_y4;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_g_out_y2;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_g_out_y4;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_r_out_y2;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> SN_r_out_y4;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y0;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y2;
  //////////////////////////////////////////////////////////
   /**
   * @brief Vector for moments in Yamamoto direct sum approach
   */
  vector<real_prec> data_g_y4;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  real_prec mean_property_value;
  //////////////////////////////////////////////////////////
  /**
   * @brief FKP Varianc
   * @details Estimation of the variance of the power spectrum using the FKP estimator with the exact expression, equation
   * 2.4.6  of FKP paper
   */
  void do_fkp_error_bars_exact(vector<real_prec> &, vector<real_prec> &, vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief FKP Variance
   * @details Computes the variance in the power spectrum as an approximation to the original FKP  estimation, introducing
   * the definition of the effective volume. This is computed using the random catalogue.
   * @return Effective volume as a function of k
   * @author ABA
   */
  void do_fkp_error_bars_veff(s_data_structure_direct_sum *, vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type Params
   */
  Params params;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Auxiliary function to compute Yamammoto related quantities
   * @author ABA
   */
  void aux_func_yamamoto(ULONG lp, real_prec kx, real_prec ky, real_prec kz, real_prec icorr, real_prec &F1,real_prec &F2);
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
 public:

/**
   * @public
   * @brief Default constructor
   * @author ABA
   * @note  The outputs of the FFTW has dimensions \f$ N1*N2*(N3/2+1)\f$, where only half of frequencies in the third component \f$(kz)\f$
   * are stored. The other half of the third component can be found using Hermitian Symmetry (see below).
   * The original output of the FFTW has the ordering
   *  @code
     index:    o_i      =  0    1     2     3    ...   N/2     N/2+1       N/2+2     ...   N-1
     freq:     k_i      =  0   k_f   2k_f  3k_f  ...   k_N   -k_n+k_f   -k_n+2k_f   ...   -k_f
     freq:    q_i/k_f   =  0    1     2     3    ...   N/2    -N/2+1      -N/2+2     ...   -1        in units of k_f                                                                                       *
   * @endcode
   *  where k_f is the fundamental mode, N is the number of grid cells per dimension in the DFT and \f$ k_n= (N/2)k_f\f$ is the Nyquist frequency. Hence, the slot for the frequency \f$-k_n\f$ is not written in the output. We need to explicitely use the negative and positive values
   * to account for all posible configurations properly, specialy when computing the Bispectrum. For the Power spectrum it is not necessary for we do not expect to exploit information at the scales of the Nyquist freq. We therfore expand the loops over the wavenumbers by introducing one more
   * slot, such that the new ordering (i) and coordinates (q_i/k_f) read as
 * @code
    index:    i       =   0   1     2      j  ...  N/2,  N/2+1    N/2+2     N/2+3    ...        j            N-2     N-1    N
    freq:     k_i     =   0   k_f  2k_f  jk_f ...  k_N,  -k_N   -k_n+k_f  -k_n+2k_f  ...   -k_n+(j-1)k_f   -3k_f    -2k_f  -k_f
    freq:     q_i/k_f =   0   1     2      j  ...  N/2   -N/2   -N/2+1     -N/2+2    ...     -N/2+(j-1)      -3      -2     -1    in units of k_f                                                          *
 * @endcode
 * Hence, the new coordinates of the modes are
 * @code
     q_i/k_f  = i      for i<=N/2,
 * @endcode
 * and
 * @code
     q_i/k_f  =  i-(N+1)    for i>=N/2+1
 * @endcode
 * Similarly, to map the new index to the original output, we have
 * @code
     o_i  = i       for i<=N/2,
 * @endcode
 * and
 * @code
     o_i = i-1    for i>=N/2+1,
 * @endcode
 * such that for \f$ i=N/2\f$, \f$o_i=i\f$ and for \f$i= N/2+1\f$, o_i is also \f$o_i=i\f$. For the third (\f$kz\f$) component of the DFT we only have \f$N3/2+1\f$ elements displayed originally as
 * @code
     index     o_k       =  0   1   2   3   i ...  N3/2
     freq     q_k/f_f    =  0   1   2   3   i ...  N3/2
 * @endcode
 * We reorder the components as we did above for the x-y components. However, in this case, the amplitude
 * at a coordinate -k_i is obtained by Hermitian symmetry \f$ \delta(kx,ky,-kz)=\delta(-kx,-ky,kz)\f$.
 * That is, when we use the complex conjugate to obtain the negative z-plane, we need to reverse the sign of the kx and ky components.
 * Then, the reordered index is
 * @code
       i3        =  0   1  2   j ...  N3/2   N3/2+1   N3/2+2    N3/2+3  ...     j      ... N3-2  N3-1   N3
    q_i3/k_f     =  0   1  2   j ...  N3/2    N3/2    N3/2-1    N3/2-2  ...  N3-(j-1) ...    3    2     1
 * @endcode
 * that is,
 * @code
   q_i3/k_f = i     for i<=N/2,
 * @endcode
 * and
 * @code
   q_i3/k_f = N3-(i-1)   for i>=N/2+1.
 * @endcode
 * The index need in the function index_3d() is \f$|q_i / k_f|\f$.
   */
  FftwFunctions(){}
//////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Constructor overloaded
   * @param Inizialization of private variables
   */
  FftwFunctions(Params _params): n_ran(0),w_g(0),w_r(0),alpha(0),s_g(0),s_r(0),sr1(0), sr2(0),sg1(0),sg2(0),normal_p(0),normal_window(0.),normal_b(0),normal(0),shot_noise(0),shot_noise2(0),shot_noise_b1(0),shot_noise_b2(0),correction_MAS_exp(0), sgrid(0),new_sgrid(0),rNft(128),DeltaK_Bis(0), Zmin(0),Zmax(0), Ymin(0),Ymax(0),Xmin(0),Xmax(0),Nshells_bk(1), rdeltax(0),rdeltay(0),rdeltaz(0),n_gal(0),Xoffset(0), Yoffset(0), Zoffset(0)
  {
    this->set_params(_params);
    this->CosmoF.set_cosmo_pars(this->params.s_cosmo_pars);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Destructor
   * @details Free fftw vector
   */
  ~FftwFunctions(){
    free_fftw_vectors();
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> data_g;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> data_g_rss;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue
   */
  vector <real_prec> data_g_mw;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Vectors for the galaxy catalogue used in case a cross power spectrum is to be measured
   */
  vector <real_prec> data_gp;
  //////////////////////////////////////////////////////////
  /**
   * @brief Width of the spherical shell (in h/Mpc) for the window function
   */
  real_prec DeltaK_window;
  //////////////////////////////////////////////////////////
    /**
     * @brief  Poisson Shot noise
     */
  real_prec shot_noise;
  //////////////////////////////////////////////////////////
   /**
     * @brief  Poisson Shot noise
     */
  real_prec shot_noise2;
  //////////////////////////////////////////////////////////
  /**
   * @brief Number of shells for the estimation of Bk as done by Jennifer
   */
  int Nshells_bk;
  //////////////////////////////////////////////////////////
  /**
   * @brief Width of the spherical shell (in h/Mpc) used int he estimates of Bispectrum
   */
  real_prec DeltaK_Bis;
  //////////////////////////////////////////////////////////
  /**
   * @brief Resize containers
   * @detail Resizes and initializes the input and output vectors (real or complex) for the FFTW
   * @note For the implementation of the Yamamoto-like estimator, the allocation of the complex vectors to the method where they are Fourier transformed, such that after the evaluation of the FFTW,
   * the memmory of the input vector is released.
   * @author ABA
   */
  void resize_fftw_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief Constant density grid assignment
   * @details Function used when no random catalogue is implemented.
   * A vector is filled with the mean number density of
   * the simulation. Such mean number density is computed
   * from the information if the size of the box
   * as given in the parameter file, together with the
   * number of objects.
   * @param vol Volumen of the sample if known
   * @result Private class member containing the galaxy fluctuation interpolated in the mesh of size NFT.
   */
  void raw_sampling(real_prec vol);
  //////////////////////////////////////////////////////////
  /**
   * @briefGet the tracer fluctuation
   * @details Build the galaxy fluctuation by subtracting data and random catalogue with factor alpha \f$ \delta(\vec{x})=\rho_{gal}(\vec{x})-\alpha \rho_{ran}(\vec{x})\f$.
   * @author ABA
   */
  void get_fluctuation();
  //////////////////////////////////////////////////////////
  /**
   * @brief Build FKP fluctuation
   * @details Build the galaxy fluctuation when ujsing a simulation box in real and reshift space.
   * @param rss Is a ghost parameter
   * @author ABA
   */
  void get_fluctuation(bool rss);
  //////////////////////////////////////////////////////////
  /**
   * @brief Computes the parameters associated to the
   * FKP estimator
   * @details e.g., normalization, shot noise, etc.
   * and assign them to public/or private variables of this class.
   * @ parameter ur Use random catalogue (true/false)
   * @author ABA
   */
  void get_parameters_estimator(bool verbose);
  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space.
   * @details This function returns the spherical average estimate
   * of the monopole, quadrupole, hexadecapole, and the
   * window function of based on the FKP estiamtor.
   * For the spherical averages, we only use a quarter of the full FOURIER box,
   * using Hermitian symmetry F(kx, ky, -kz)=F(-kx, -ky, kz)* to
   * recover the  information in the negative frequencies explicitely.
   * When counting modes and power, we weight by a factor 2
   * order to account for the negative z-quadrant. Although
   * the factor 2 cancels out when computing the average power in each
   * shell this allows comparisons with codes using the full box.
   * This subroutine is ideal also for the multipole decomposition
   * in the modes l=0, 2 and 4. For other moments, the full excusrion
   * through Fourier space has to be done.
   * BINNING: the floor function ensures that we are using intervals
   * of the form [). The spherical shells are such that the first bin
   * starts at the zero frequency (althought that mode is exlcuded
   * for the power spectrum, not for the window function),
   * Note therefore that in the case ndel=1, the fist bin
   * will contain ONLY one Fourier mode, is the zero frequency.
   * Therefore, it is convinient to start with ndel=2
   * This function is called by get_power_spectrum_fkp().
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   * @return p_g2 quadrupole power spectrum
   * @return p_g4 hexadecapole power spectrum
   * @return p_r power spectrum of the window function
   * @return p_2d 2d power spectrum in cartesian coordinates
   * @return p_2s 2d power spectrum in polar coordinates
   * @return nm Number of modes in spherical shells
   * @author ABA
   */
  void power_spectrum_fkp(vector<real_prec>&p_g0, vector<real_prec>&p_g2, vector<real_prec>&p_g4,vector<real_prec>&p_r, vector< vector<real_prec> >&p_2c, vector< vector<real_prec> >&p_2s, vector<int>&nm);
  //////////////////////////////////////////////////////////
  /**
   * @brief Power spectrum estimation
   * @return p_g0 Monopole power spectrum
   * @return nm Number of modes in k-bins
   * @author ABA
   */
  void power_spectrum(vector<real_prec>&p_g0,  vector<int>&nm);
  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space for Bispectrum
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   */
  void power_spectrum_fkp_for_bispectrum(vector<real_prec>&p_g0);
  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space.
   * @details This function returns the spherical average estimate
   * of the monopole, quadrupole, hexadecapole, based on the Yamamoto estimator.
   * For the spherical averages, we only use a quarter of the full FOURIER box,
   * using Hermitian symmetry F(kx, ky, -kz)=F(-kx, -ky, kz)* to
   *  recover the  information in the negative frequencies explicitely.
   *  When counting modes and power, we weight by a factor 2
   *  order to account for the negative z-quadrant. Although
   *  the factor 2 cancels out when computing the average power in each
   *  shell this allows comparisons with codes using the full box.
   *  This subroutine is ideal also for the multipole decomposition
   *  in the modes l=0, 2 and 4. For other moments, the full excurion
   *  through Fourier space has to be done.
   *  This function is called by get_power_spectrum_yamammoto().
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   * @return p_g2 quadrupole power spectrum
   * @return p_g4 hexadecapole power spectrum
   * @return mod Number of modes in spherical shells
   */
  void power_spectrum_yamamoto(vector<real_prec>&p0, vector<real_prec>&p2, vector<real_prec>&p4,vector<int>&mod);
  //////////////////////////////////////////////////////////
  /**
   * @brief Shell average in Fourier space.
   * @details This function returns the spherical average estimate
   * This is a modifiction to power_spectrum_yamamoto()
   * in which we do only one octant and visit the other three (z>0) as we did for fkp.
   * @param s_b structure of type s_parameters_box
   * @return p_g0 monopole power spectrum
   * @return p_g2 quadrupole power spectrum
   * @return p_g4 hexadecapole power spectrum
   * @return mod Number of modes in spherical shells
   */
  void power_spectrum_yamamoto_new(vector<real_prec>&p0, vector<real_prec>&p2, vector<real_prec>&p4,vector<int>&mod);
  //////////////////////////////////////////////////////////
  /**
  * @brief Shell average in Fourier space.
  * @details This function returns the spherical average estimate
  * This is a modifiction to power_spectrum_yamamoto()
  * in which we do only one octant and visit the other three (z>0) as we did for fkp.
  * @warning This method will be deprecated as has been shown that the method <power_spectrum_yamamoto_new> provides the same results with the same way of traveling the Fourier box as the methods dedicated for FKP
  */
  void _power_spectrum_yamamoto(vector<real_prec>&p0, vector<real_prec>&p2, vector<real_prec>&p4,vector<int>&mod);
  //////////////////////////////////////////////////////////
  /**
   * @brief  Computes the outputs for the P(k)
   * @detail Computes the outputs for the P(k) to generate the estimates of shot noise for Bispectrum
  */
  void get_power_spectrum_for_bispectrum(vector<real_prec>&p0);
  //////////////////////////////////////////////////////////
  /**
   *@brief  Compute the multipole decomposition using FKP estimator
   * @param s_b structure of type s_params_box
   * @return p0 Monopole power spectrum
   * @return p2 Quadrupole power spectrum
   * @return p4 Hexadecapole power spectrum
   * @return pr power spectrum of the window function
   * @return p2d 2D power spectrum in cartesian coordinates
   * @return p2s 2D power spectrum in polar coordinates
   * @return nm Number of modes in spherical shells
  */
 void get_power_spectrum_fkp(vector<real_prec>&p0,vector<real_prec>&p2,vector<real_prec>&p4, vector<real_prec>&pr, vector<vector<real_prec> >&p2d,  vector<vector<real_prec> >&p2s,vector<int>&nm);
 //////////////////////////////////////////////////////////
 /**
  *@brief  Compute the multipole decomposition using FKP estimator
  * @details This method is used when real and redhsoft space power spectrum are to be computed. Hence, it is valif for a simulation.
  * @param s_b structure of type s_params_box
  * @return p0 Monopole power spectrum
  * @return p2 Quadrupole power spectrum
  * @return p4 Hexadecapole power spectrum
  * @return pr power spectrum of the window function
  * @return p2d 2D power spectrum in cartesian coordinates
  * @return p2s 2D power spectrum in polar coordinates
  * @return nm Number of modes in spherical shells
 */
 void get_power_spectrum_fkp(vector<real_prec>&p, vector<real_prec>&p0,vector<real_prec>&p2,vector<real_prec>&p4, vector<real_prec>&pr, vector<vector<real_prec> >&p2d,  vector<vector<real_prec> >&p2s,vector<int>&nm);
 //////////////////////////////////////////////////////////
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
 void cross_power_spectrum_fkp(vector<real_prec>&p0,vector<int>&nm, vector<real_prec>&corr);
#else
 void cross_power_spectrum_fkp(vector<real_prec>&p0,vector<int>&nm, bool);
#endif
 //////////////////////////////////////////////////////////
 /**
  *@brief  Compute the multipole decomposition using Yamamoto estimator
  *@details with the FFTW based scheme.
  * @param s_b structure of type s_params_box
  * @result p0 monopole power spectrum
  * @result p2 quadrupole power spectrum
  * @result p4 hexadecapole power spectrum
  * @result p_r power spectrum of the window function
  * @result nm Number of modes in spherical shells
  */
 void get_power_spectrum_yamamoto(vector<real_prec>&p0,vector<real_prec>&p2,vector<real_prec>&p4,vector<int>&nm);
 //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of the variance for the power spectrum
   * @details based on the FKP estimator. Selects between the
   * exact and the approximate expression for the var(P)
   * @param s_b structure of type s_params_box
   * @param s_d structure of type s_data_structure
   * @param kv k-spherical shells
   * @param pk power spectrum in k-spherical shells
   * @param nm number of modes in k-spherical shells
   * @result sig FKP variance of the power spectrum
   */
  void get_fkp_error_bars(s_data_structure_direct_sum *s_d, vector<real_prec> &kv, vector<real_prec> &pk, vector<int>&nm, vector<real_prec> &sig);
  //////////////////////////////////////////////////////////
  /**
   * @brief Count of triangles in Fourier space
   * @details
   * @author Jennifer Pollack & Emiliao Sefusati & AB
   */
  void bispectrum_fkp(char, vector<real_prec> &, vector<real_prec> &, vector<real_prec> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of Bispectrum
   * @details This function generates the estimates of Bispectrum based on the counts of triangles the normalization and the shot-noise corrections.
   * @param s_b structure of type s_parameters_box
   * @warning This method needs to be revised/updated
   */
  void get_bispectrum_fkp(char, vector<real_prec> &, vector<real_prec> &, vector< int > &);
  //////////////////////////////////////////////////////////
  /**
   * @brief Estimates of Bispectrum
   * @details  based on trick by Scoccimarro and Jennifer
   * @param s_b structure of type s_parameters_box
   * @note The binning in the method that computes the power spectrum for the shot noise must have the same binning as that used in the bispectrum. Check that please.
   * @author Jennifer Pollack & ABA
   */
  void get_bispectrum_fkp_fast(vector<real_prec> &, vector<real_prec> &, vector< int > &, string file);
  //////////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void define_kshells();
  //////////////////////////////////////////////////////////
  /**
   * @brief Get the inverse Fourier transform in each Fourier shell
   * @details  Construc the arrays with the inverse Fourier transforms in shells.
   * For each shell we allocate an array with dimension N*N*N, where N is the maximum
   * size of the Nft used, i.e, that of the very last shell.
   * @param Pass the definition of binning
   */
  void get_ift_shells_bispectrum();
  //////////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void loop_shells_bispectrum(vector<real_prec> &pk, vector<real_prec> &bispect, vector< int > &mod, string file);
///////////////////////////////////////////////////////
  /**
   * @brief Define k-shells for the estimation of bispectrum
   * @param s_b structure of type s_parameters_box
   */
  void construct_shells(int ngrid, int kmnid, int kmxid, vector<real_prec> &iFT_output_delta, vector<real_prec> &iFT_output_triangles,vector<real_prec> &iFT_output_p1_cyc_sum);
//////////////////////////////////////////////////////////
  /**
   * @brief Mapping of vectos from one octant to the octants in the kz>0 sub-volume of Fourier space
   * @param s_b structure of type s_parameters_box
   */

  void cellsym(int id, int ngrid,complex_prec *data_ks, complex_prec *data_dk,complex_prec *data_pk_sn);
  //////////////////////////////////////////////////////////
  /**
   * @brief Evaluates the DFT required for the estimates of
   * multipole decomposition
   * @details of the power spectrum using the Yamamoto-Blake estimator
   * obtained from a direct-sum approach. Sampling of the elements used in the Yamamoto et al. estimator            *
   * for the moments of the 3d power spectrum. We are generating an output with the same structure as that
   * given by the FFTW algorithm, i.e, allocating a complex vector
   * with positive and negative frequencies except for the third component
   * which has only the positive quadrant. Furthermore, the first two
   * components only have the positive Nyquist frequency.
   * The negative z components will be then
   * recovered using the Hermitian symmetry, as is done for the power spectrum
   * computed using the FFTW.
   * Here we are explicitely computing the multipoles l=0,2,4
   * @param s_d structure of type s_data_structure
   */
  void get_power_moments_fourier_grid_ds_yam(s_data_structure_direct_sum *s_d);

  //////////////////////////////////////////////////////////
  /**
   * @brief Compute shell-averaged multipole decomposition of the
   * power spectrum obtained from the  Yamamoto-Blake estimator
   * implemening a a direct-sum approach.
   * @result p0 monopole
   * @result p2 quadrupole
   * @result p4 hexadecapole
   * @result nm number of modes in spherical shell
   */
  void power_yam_1d_ds(vector<real_prec>&p0 ,vector<real_prec>&p2 ,vector<real_prec>&p4, vector<int> &nm);
  //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on screen
   */
  void write_fftw_parameters();
 //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on .log file
   * @param p void pointer
   * @param log_f log file
   */
  void write_fftw_parameters(void *p, string log_f);
  //////////////////////////////////////////////////////////
  /**
   * @brief Release memmory
   */
  void free_fftw_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief Write parameters on .log file
   * @param p void pointer
   * @param log_f log file
   */
  vector<real_prec>mass_cuts;
  //////////////////////////////////////////////////////////
  /**
   * @brief set values of tracer overdenisty
   */
  void set_data_g(ULONG i, real_prec value){
    this->data_g[i]=value;
  }
//////////////////////////////////////////////////////////
  /**
   * @brief set values of tracer overdenisty
   */
  void set_data_r(ULONG i, real_prec value){
    this->data_r[i]=value;
  }
//////////////////////////////////////////////////////////
  /**
   * @brief set values of tracer overdenisty
   */
  void set_data_gp(ULONG i, real_prec value){
    this->data_gp[i]=value;
  }
 //////////////////////////////////////////////////////////
  /**
   * @brief Retrieve values of tracer overdenisty
   */
  real_prec _data_g(ULONG i){
    return this->data_g[i];
  }
 //////////////////////////////////////////////////////////
  /**
   * @brief set values of tracer overdenisty
   */
  void set_data_g_rss(ULONG i, real_prec value){
    this->data_g_rss[i]=value;
  }
 //////////////////////////////////////////////////////////
 /**
   * @brief Retrieve values of real part of Fourier transform of tracer overdensity
   */
  real_prec _data_out_g_real(ULONG i){
    return this->data_out_g[REAL][i];
  }
//////////////////////////////////////////////////////////
  /**
   * @brief Retrieve values of imaginary part of Fourier transform of tracer overdensity
   */
  real_prec _data_out_g_imag(ULONG i){
    return this->data_out_g[IMAG][i];
  }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::n_ran
   * @return FftwFunctions::n_ran
   */
  ULONG _n_ran(){return n_ran; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::w_g
   * @return FftwFunctions::w_g
   */
  real_prec _w_g(){return w_g; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::w_r
   * @return FftwFunctions::w_r
   */
  real_prec _w_r(){return w_r; }
 //////////////////////////////////////////////////////////
    /**
   * @brief get the value of private member Number of FftwFunctions::alpha
   * @return FftwFunctions::alpha
   * @details The ratio between the weighted number of galaxies and the weighted number of randoms
   * @note This is computed as \f$ \alpha = \frac{\sum_{i=gal}w_{fkp}(r_i)}{\sum_{j=ran}w_{fkp}(r_j)}\f$
   */
  real_prec _alpha(){return alpha; }
  //////////////////////////////////////////////////////////
     /**
    * @brief set the value of private member Number of FftwFunctions::alpha
    */
  void set_alpha(real_prec nalpha){this->alpha=nalpha; }
  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::normal_power
   * @return FftwFunctions::normal_power
   */
  real_prec _normal_power(){return normal_power; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_normal_power(real_prec new_normal_power){this->normal_power=new_normal_power; }
  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of FftwFunctions::normal_power
   * @return FftwFunctions::normal_power
 */
  real_prec _normal_power_two(){return normal_power_two; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_normal_power_two(real_prec new_normal_power_two){this->normal_power_two=new_normal_power_two; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_normal_p(real_prec new_normal_power){this->normal_p=new_normal_power; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_normal_b(real_prec new_nor){this->normal_b=new_nor;}
  //////////////////////////////////////////////////////////
  /**
   *    @brief  Return normalization of window function
   */
  real_prec _normal_window(){return normal_window; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_shot_noise(real_prec val){this->shot_noise=val;}
  real_prec _shot_noise(){return this->shot_noise; }
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_shot_noise2(real_prec val){this->shot_noise2=val;}
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_shot_noise_window(real_prec val){this->shot_noise_window=val;}
  real_prec _shot_noise_window(){return this->shot_noise_window;}
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  ULONG _n_gal(){return n_gal; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  void set_n_gal(ULONG new_ngal){this->n_gal=new_ngal; }
  //////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  void set_n_ran(ULONG new_ng){this->n_ran=new_ng; }
//////////////////////////////////////////////////////////
  /**
   * @brief get the value of private member Number of  FftwFunctions::n_gal
   * @return FftwFunctions::n_gal
  */
  void set_s_r(real_prec val){this->s_r=val;}
  real_prec _s_r(){return this->s_r;}
  //////////////////////////////////////////////////////////
    /**
     * @brief get the value of private member Number of  FftwFunctions::n_gal
     * @return FftwFunctions::n_gal
    */
  void set_s_g(real_prec val){this->s_g=val;}
  real_prec _s_g(){return this->s_g;}
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_sg1(real_prec val){this->sg1=val;}
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_sg2(real_prec val){this->sg2=val;}
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_sr1(real_prec val){this->sr1=val;}
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_sr2(real_prec val){this->sr2=val;}
   //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
  * @brief
  * @return
  */
  void set_w_r(real_prec val){this->w_r=val;}
//////////////////////////////////////////////////////////
/**
  * @brief get the value of private member Number of  FftwFunctions::n_gal
  * @return FftwFunctions::n_gal
  */
  void set_w_g(real_prec val){this->w_g=val;}
//////////////////////////////////////////////////////////
/**
 * @brief set new value for the variable mean_property_value
*/
 void set_mean_property_value(real_prec val){this->mean_property_value=val;}
//////////////////////////////////////////////////////////
/**
 * @brief set a new version of an intance of an object of type Params
*/
  void set_params(Params par);
//////////////////////////////////////////////////////////

};
#endif
