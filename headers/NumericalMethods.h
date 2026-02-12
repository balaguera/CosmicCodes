////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This file contains functions for numerical analysis
 * @file NumericalMethods.h
 * @author Andres Balaguera Antolinez
 * @version 1.0
 * @date 2017-2024
 */
////////////////////////////////////////////////////////////////////////////
#ifndef __NMETHODS__
#define __NMETHODS__
# include <boost/tuple/tuple.hpp>
# define GNUPLOT_ENABLE_PTY
//# include "../gnuplot-iostream/gnuplot-iostream.h"   //http://stahlke.org/dan/gnuplot-iostream/
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include <sstream>
# include <fftw3.h>
# include <getopt.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_histogram.h>
# include <gsl/gsl_histogram2d.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_expint.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_result.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_complex.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_statistics_float.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sort_vector.h>
# include <gsl/gsl_sort_vector_float.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_poly.h>
# include <vector>
# include <numeric>
# include <algorithm>
# include "Structures.h"
# include "GnuplotC.h"
#ifdef _USE_OMP_
# include <omp.h>
#endif
#ifdef _USE_PYTHON_
# include <Python.h>
//https://www.codeproject.com/Articles/820116/Embedding-Python-program-in-a-C-Cplusplus-code
#endif
#ifdef _USE_HEALPIX_
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <healpix_data_io.h>
# include <healpix_base.h>
# include <alm_powspec_tools.h>
# include <alm_healpix_tools.h>
# include <xcomplex.h>
# include <sort_utils.h>
# include <sharp_cxx.h>
#endif
# include <unistd.h>
# include <gsl/gsl_dht.h>
# include "fftw_array.h"
# include "lss_vector.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////
#ifndef FFTW_ARRAY
#define FFTW_ARRAY
//#include<iostream>
/* This template allows to avoid writing fftw_free */

template<typename T> class fftw_array
{
 public:
  T *data;
  fftw_array(long size)
    {
#ifdef SINGLE_PREC
      data = (T *) fftwf_malloc(size*sizeof(T));
#endif
#ifdef DOUBLE_PREC
      data = (T *) fftw_malloc(size*sizeof(T));
#endif
    }
  ~fftw_array()
    {
#ifdef SINGLE_PREC
      fftwf_free(data);
#endif
#ifdef DOUBLE_PREC
      fftw_free(data);
#endif
    }
  /* implicit conversion: works for functions, but not for templates. For the latter case write: <array_name>.data */
  operator T*()
  {
    return data;
  }
};
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 *@struct<s_eigen_vector>
 *@brief Structure to allocate eigen values and eigenvectors
 *@details This is mainly used in PCA and Tidal field analysis
 */
struct s_eigen_vector{
    /**
     *@brief Eigenvalue
     */
    real_prec eigen_val;
    /**
     *@brief Components of the eigenvector associted to the eigenvalue
     */
    vector<real_prec> eigen_vec;
    /**
     *@brief Number of eigenvalues
     */
    int Neigen;
    /**
     *@brief Get Number of eigenvalues
    */
    int Nkbins(){return this->eigen_vec.size();}
};
////////////////////////////////////////////////////////////////////////////
/**
 *@brief Randomize elements of a vector
 *@param data input vector
 *@returns elements of input vector placed in random positions
 */
void randomize_vector(vector<ULONG>&data);
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 *@brief
 */
void get_scalar(string FNAME,vector<real_prec>&OUT,ULONG N1,ULONG N2,ULONG N3);
#endif
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *, gsl_real, gsl_real);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(int N, gsl_real (*function)(gsl_real, void *) ,void *,gsl_real ,gsl_real);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>&, vector<gsl_real>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>&, vector<gsl_real>&, bool);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>&, vector<gsl_real>&, real_prec, bool);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *,vector<gsl_real>&, vector<gsl_real>&, gsl_real);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
void gsl_get_GL_weights(gsl_real,gsl_real, vector<real_prec>&, vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical integration using GSL
 * @returns 
 */
void gsl_get_GL_weights(gsl_real LowLimit,gsl_real UpLimit,  gsl_integration_glfixed_table *wf, vector<gsl_real> &XX, vector<gsl_real>&WW);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical interpolation
 * @returns 
 */
real_prec gsl_inter_pointers(real_prec *, real_prec *, int, real_prec );
#ifdef SINGLE_PREC
real_prec gsl_inter_pointers(double *, double *, int, double );
#endif

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical interpolation
 * @returns 
 */
real_prec gsl_inter(gsl_real *, gsl_real *, int, gsl_real);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical interpolation
 * @returns 
 */
real_prec gsl_inter_new(const vector<gsl_real> &, const vector<gsl_real> &, gsl_real);
#ifdef SINGLE_PREC
real_prec gsl_inter_new(const vector<float> &, const vector<float> &, gsl_real);
#endif

////////////////////////////////////////////////////////////////////////////
/**
 * @brief Numerical interpolation
 * @returns 
 */
real_prec gsl_inter_new2(vector<real_prec>&, vector<vector<real_prec>>&, int, real_prec);

////////////////////////////////////////////////////////////////////////////
#ifdef SINGLE_PREC
/**
 * @brief Numerical interpolation
 * @returns 
 */
real_prec gsl_inter_new2(vector<gsl_real>&, vector<vector<gsl_real>>&, int, real_prec);
#endif
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
real_prec MAS_CIC_public(real_prec x);

////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel(gsl_real (*function)(gsl_real, void *), void *p, gsl_real aux_var, gsl_real LowLimit,gsl_real UpLimit);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf(gsl_real (*function)(gsl_real, void *), void *,  gsl_real, gsl_real);
///////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
real_prec gsl_integration_sin_kernel_lowlimit_inf2(gsl_real (*function)(gsl_real, void *), void *, gsl_real, gsl_real, gsl_integration_workspace *, gsl_integration_workspace *);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void gsl_bspline(vector<gsl_real> &, vector<gsl_real>&, vector<gsl_real> &, vector<gsl_real> &);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void indexx(vector<int>&arr, vector<int>&indx);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void indexx_ulong(vector<ULONG>&arr, vector<ULONG>&indx);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void sort_index(int , int , int , int *, int *, int *);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void matrix_inversion(vector< vector<real_prec> > &, vector< vector<real_prec> > &);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
real_prec linInterpol(ULONG N, real_prec L, real_prec delta, real_prec xpos, real_prec ypos, real_prec zpos, const vector<real_prec>&tab);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void matrix_det(vector< vector<real_prec> >&, real_prec &);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void get_eigen(vector<vector<real_prec> >&icov_masses, vector<real_prec>&masses);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void get_eigen(vector<real_prec> &icov_masses, vector<s_eigen_vector>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void get_eigen(vector<real_prec> &icov_masses, vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void get_det_matrix(vector<vector<real_prec> >&matriz, real_prec &determinant);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void sort_vectors(vector<ULONG>& v1,vector<ULONG>&v2, vector<ULONG>&v3, vector<ULONG>& v4, vector<ULONG>&v5, vector<ULONG>& v6,vector<ULONG>
                  & v7, vector<real_prec>& v8);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Real part of spherical harmonics
 */
real_prec real_sh(int l, int m, real_prec theta, real_prec phi);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Imaginary part of spherical harmonics
 */
real_prec imag_sh(int l, int m, real_prec theta, real_prec phi);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Spherical Bessel funciton of order l
 */
real_prec bessel(real_prec , int);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief First derivative of spherical Bessel funciton of order l
 */
real_prec dbessel(real_prec , int );
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Second derivative of spherical Bessel funciton of order l
 */
real_prec ddbessel(real_prec , int );
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Counter of lines in an reading file
 */
size_t m_countLines(istream& inStream);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 */
size_t m_getBuffer(istream& inStream);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
string dto_string (real_prec Number);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the mean on an input vector
 */
real_prec get_mean(const vector<real_prec> &ini);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the minimum of an input vector
 */
real_prec get_min_nm(const vector<real_prec> &ini,bool zero);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the minimum of an input vector
 */
real_prec get_min_nm(const vector<real_prec> &ini);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the maximum of an input vector
 */
real_prec get_max_nm(const vector<real_prec> &ini);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the variance of ele elements of an input vector
 */
real_prec get_var(const vector<real_prec> &ini);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Compute the variance of elements of an input vector
 */
real_prec get_var(real_prec mean, const vector<real_prec> &ini);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Performs a logarithmic smoothing of the input vector
 */
void log_smooth(vector<real_prec>&xin, vector<real_prec>&vin);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Performs a linear smoothing of the input vector
 */
void lin_smooth(vector<real_prec>&xin, vector<real_prec>&vin, int n);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void gradfindif(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3,vector<real_prec>in,vector<real_prec>&out,int dim);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void convert_ngp_to_cic(vector<real_prec>&in, vector<real_prec>&out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_12d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int v, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns, int Nv);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_11d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_11d(ULONG i, ULONG j, ULONG k, ULONG l,ULONG m,ULONG n,ULONG o, ULONG p, ULONG q,ULONG r, ULONG s, ULONG Nj, ULONG Nk, ULONG Nl, ULONG Nm, ULONG Nn, ULONG No, ULONG Np, ULONG Nq, ULONG Nr, ULONG Ns);

////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_10d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_9d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq);

////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_8d(int i, int j, int k, int l, int m, int n, int o, int p, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np);

////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_7d(int i, int j, int k, int l, int m, int n, int o, int Nj, int Nk, int Nl, int Nm, int Nn, int No);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
ULONG index_6d(int i, int j, int k, int l, int m, int n, int Nj, int Nk, int Nl, int Nm, int Nn);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 * @returns ULONG index
 */
ULONG index_5d(int i, int j, int k, int l, int m, int Nj, int Nk, int Nl, int Nm);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 * @returns ULONG index
 */
ULONG index_4d(int i, int j, int k, int l, int Nj, int Nk, int Nl);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Computes the 3D index from the label coordinates i, j, k
 * @returns ULONG index
 */
ULONG index_3d(ULONG i, ULONG j, ULONG k, int Nj, int Nk);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Computes the 3D index from the label coordinates i, j, k
 * @returns ULONG index
 */
ULONG index_3d(ULONG i, ULONG j, ULONG k, ULONG Nj, ULONG Nk);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Computes the 2D index from the label coordinates i, j
 * @returns ULONG index
 */
ULONG index_2d(ULONG i, ULONG j, ULONG Nj);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Given a mesh of N (per-dimension), this function returns the coordinates (indices) of the cell
 * @params index = C-ordered (or Fortran) index of a cell (0,N*N*N-1)
 * @params N = Number of cells per dimention
 * @returns indices XG, YG, ZG of the cell if i is in c/row-major order or
 * @returns indices ZG, YG, XG of the cell if i is in Fortran/column-major order
 */
void index2coords(ULONG index, ULONG N, ULONG &XG, ULONG &YG, ULONG &ZG );
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Given cartasian coordiantes and box info, this returns the ID of the cell where the object is located
 */
ULONG grid_ID(s_params_box_mas *params, const real_prec &x, const real_prec &y, const real_prec &z);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Real to Complex FFT
 * @param Nft Size of the mes per dimenion
 * @param in array wth input field
 * @returns out Complex array containing the Fourier transform of the field in
 */
void do_fftw_r2c(ULONG Nft, vector<real_prec>&in, complex_prec *out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Complex to real FFTW
 * @param Nft Size of the mes per dimenion
 * @param in Complex array
 * @returns out real array with the inverse Fourier transform
 */
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Real to Complex FFT in 1d
 * @param Nft Size of the mes per dimenion
 * @param in array wth input field
 * @param dir +1 for real to complex, -1 for inverse Fourier transform
 * @returns out Complex array containing the Fourier transform of the field in
 */
void do_fftw_1d(ULONG Nft, complex_prec *in, complex_prec *out, int  dir);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Real to Complex FFT in 1d
 * @param Nft Size of the mes per dimenion
 * @param in array wth input field
 * @returns out Complex array containing the Fourier transform of the field in
 */
void do_fftw_1d_r2c(ULONG Nft, vector<real_prec>&in, complex_prec *out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Real to Complex FFT in 1d
 * @param Nft Size of the mes per dimenion
 * @param in Complex array
 * @returns out real array with the inverse Fourier transform
 */
void do_fftw_1d_c2r(ULONG Nft,  complex_prec *in, vector<real_prec>&out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  3D FFTW
 */
void do_fftw_3d(ULONG Nft, bool direction, complex_prec *in, complex_prec *out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function applies a low-pass filter (top-hat) to a high resolution field on a mesh Nft_HR**3
 * @returns container with density field with resolution Nft_LRr**3
 */
void low_pass_filter(ULONG Nft_HR, ULONG Nft_LR, int imas, bool correct, vector<real_prec>&HR_field, vector<real_prec>&LR_field, real_prec Lbox);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function applies a low-pass filter (top-hat) to a high resolution field on a mesh Nft_HR**3
 * @returns container with density field with resolution Nft_LRr**3
 */
void low_pass_filter(ULONG Nft_HR, ULONG Nft_LR, int imas, bool correct, vector<real_prec>&HR_field, complex_prec *LR_FOURIER, real_prec Lbox);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function applies a low-pass filter (top-hat) to a high resolution field on a mesh Nft_HR**3
 * @returns container with density field with resolution Nft_LRr**3
 */
void average(real_prec Lside, ULONG Nft_HR, ULONG Nft_LR, vector<real_prec>&HR_field, vector<real_prec>&LR_field);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void get_cumulative(const vector<real_prec>&, vector<real_prec> &, unsigned long &);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 * @returns
 */
void getDensity_NGP(ULONG Nft, ULONG Nft_HR, real_prec L1, real_prec d1, real_prec delta_hr, const vector<real_prec>&we, vector<real_prec>&delta);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 */
void exchange_xy(int, const vector<real_prec>&,vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 */
void exchange_xz(int, const vector<real_prec>&,vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void create_GARFIELD_FIXED_AMP(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed);
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void kernelcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, real_prec smol, int filtertype, string output_dir);
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void calc_twolptterm(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, vector<real_prec>&phiv, vector<real_prec> &m2v);
#endif
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void calc_LapPhiv(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, complex_prec *philv,vector<real_prec>&LapPhiv,int index1,int index2);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  This function finds the index of the neighbour cells for all cells in a mesh
 */
void calc_LapPhiv(ULONG N1,real_prec L1, complex_prec *delta,vector<real_prec>&LapPhiv,int index1,int index2);
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void calc_curlcomp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&phiv, vector<real_prec> phiv2, vector<real_prec> &m2v, int comp);
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void calc_mu2term(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec>&m2v);
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void calc_Det(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &in, vector<real_prec> &out);
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
/**
 * @brief 
 */
void convcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, vector<real_prec> &in, vector<real_prec> &out, int filtertype,real_prec smol,string file_kernel);
#endif
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void convolvek(ULONG N1, vector<real_prec>&in, vector<real_prec> &kernel, vector<real_prec> &out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void sort_2vectors(vector<vector<ULONG> >& v1,vector< vector<ULONG> >&v2);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief sorts vectors v1 and also returns v2 sorted according to v1. For closest neightbopugh applications, this
 * returns the element vc=v2[0] where v2 is sorted
 */
void sort_1dvectors(vector<ULONG>& v1,vector<ULONG> &v2);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief sorts vectors v1 and also returns v2 sorted according to v1. For closest neightbopugh applications.
 * @return sorted vector v1, v2 and the element vc=v2[0] where v2 is sorted. Returning this value saves one loopkin the funciton.
 */
void sort_1dvectors_v2(vector<ULONG>& v1, vector<ULONG> &v2, ULONG &v1c, ULONG &v2c); // 2nd version of sort_1dvectors
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void sort_1dvectors_v3(vector<ULONG>& v1, vector<ULONG> &v2, ULONG &v2c); // 2nd version of sort_1dvectors
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void sort_1dvectors_iv2(vector<int>& v1, vector<int> &v2, int &v1c, int &v2c); // 2nd version of sort_1dvectors
////////////////////////////////////////////////////////////////////////////
/**
 * @brief 
 */
void swap_amp_fourier(ULONG, vector<real_prec>& v1, vector<real_prec> &v2); // 2nd version of sort_1dvectors
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 */
void reduce_resolution(real_prec L, ULONG Nft_h, ULONG Nft_low, vector<real_prec>&field_hr, vector<real_prec>&fieldlr);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief
 */
void smooth_distribution(vector<real_prec>&x, vector<real_prec>&F,vector<real_prec>&Fs);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Get the Pearson correlation coefficient between the elements of input containers
 */
real_prec pearson_correlation(vector<real_prec>&, vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Get the Pearson correlation coefficient between the elements of input containers
 */
real_prec pearson_correlation(real_prec, vector<real_prec>&, real_prec, vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief  Get minimum and maximum of elements in container
 */
void min_max_vector(vector<real_prec>vec, size_t &imin, size_t &imax);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Swap input values
 */
template <class Type>void SWAP_L(Type &a, Type &b)
{
  Type dum=a;
  a=b;
  b=dum;
}
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Auxiliary function
 */
template <class Type>void indexx(vector<Type>&arr, vector<Type>&indx)
{
  const int M=7,NSTACK=50;
  ULONG n = indx.size();
  long ir;
  long i,j,k;
  ULONG indxt;
  long jstack=-1;
  long l=0;
  Type a;
  int istack[NSTACK];
  
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
	       a=arr[indxt];
	       for (i=j-1;i>=l;i--) 
          {
	          if (arr[indx[i]] <= a) break;
       	  indx[i+1]=indx[i];
          }
	      indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_L<Type>(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
        	SWAP_L<Type>(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	      SWAP_L<Type>(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
       	SWAP_L<Type>(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
    	 do i++; while (arr[indx[i]] < a);
	     do j--; while (arr[indx[j]] > a);
	     if (j < i) break;
	       SWAP_L<Type>(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
	       cerr << "NSTACK too small in indexx." << endl;
      	exit(1);
      }
      if (ir-i+1 >= j-l) {
	      istack[jstack]=ir;
	      istack[jstack-1]=i;
      	ir=j-1;
      } 
      else {
	     istack[jstack]=j-1;
       istack[jstack-1]=l;
    	l=i;
     }
    }
  }
}
////////////////////////////////////////////////////////////////////////////
/**
* @brief This template functions calculates the rank of the vector v1
* @details The rank is the position of the lowest value in the array v1
*/
template <typename Type> void sort_1d_vectors(vector<Type>&v1, ULONG &rank)
{
   ULONG n=v1.size();
   vector<Type>iwksp(n,0);
   indexx<Type>(v1,iwksp);
   rank=iwksp[0];
}
////////////////////////////////////////////////////////////////////////////
/**
* @brief This template functions calculates the rank of the vector v1
* @detAIL The rank is the position of the lowest value in the array v1
*/
template <typename Type> void sort_1d_vectors(vector<Type>&v1, vector<Type>&v1_s)
{
   ULONG n=v1.size();
   vector<Type>iwksp(n,0);
   indexx<Type>(v1,iwksp);
   for(int i=0;i<n;++i)
       v1_s[i]=v1[iwksp[i]];
}
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Get maximum of elements in container
 * @params in input 1D container
 */
template<typename Type> Type get_max(const vector<Type> &in)
{

  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;

  double lka=-static_cast<Type>(LARGE_NUMBER);
  double lkb;
  for(ULONG i=0;i<in.size();++i)
    {
      lkb=max(static_cast<double>(in[i]), lka);
      lka=lkb;
    }
  return static_cast<Type>(lkb);
}
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Get maximum of elements in container
 * @params in input 2D container
 */
template<typename Type> Type get_max(const vector<vector<Type>> &in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;

  Type lka=-static_cast<Type>(LARGE_NUMBER);
  Type lkb;
  for(ULONG i=0;i<in[0].size();++i)
    for(ULONG j=0;j<in[0].size();++j)
      {
      lkb=max(static_cast<Type>(in[i][j]), lka);
      lka=lkb;
    }
  return lkb;
}
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Get minimum of elements in container
 * @params in input container
 */
template<typename Type> Type get_min(const vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  double lka=static_cast<double>(LARGE_NUMBER);
  double lkb;
  for(ULONG i=0;i<in.size();++i)
    {
      lkb=min(static_cast<double>(in[i]),lka);
      lka=lkb;
    }
  return static_cast<Type>(lka);
}
////////////////////////////////////////////////////////////////////////////
  /**
   * @brief Get minimum of elements in container
   * @params in input container
   * @params exclude_zero If true, the minimum excludes the value in=0
   */
template<typename Type> Type get_min(const vector<Type>&in, bool exclude_zero)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  double lka=static_cast<double>(LARGE_NUMBER);
  double lkb;
  if(exclude_zero){
  for(ULONG i=0;i<in.size();++i)
    {
      if(in[i]!=0)
        {
          lkb=min(static_cast<double>(in[i]),lka);
          lka=lkb;
            }
    }
  }
  else
  for(ULONG i=0;i<in.size();++i)
    {
          lkb=min(static_cast<double>(in[i]),lka);
          lka=lkb;
    }
  return static_cast<Type>(lka);
}
////////////////////////////////////////////////////////////////////////////
  /**
   * @brief  Add the elements of input container
   */
template<typename Type> Type get_sum(const  vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;

  double lkb=0;
  //#pragma omp for reduction(+:lkb)
  for(ULONG i=0;i<in.size();++i)
    lkb+=static_cast<double>(in[i]);
  return static_cast<Type>(lkb);
}
////////////////////////////////////////////////////////////////////////////
template <typename Type> void check_if_positive_number(Type a){
    if(a<0)
        throw std::domain_error( "Received negative value" );
}
////////////////////////////////////////////////////////////////////////////
template <typename Type> void check_if_negative_number(Type a){
    if(a>0)
        throw std::domain_error( "Received positive value" );
}
#endif
