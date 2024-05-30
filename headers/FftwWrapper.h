#ifndef _NMETHODS_
#define _NMETHODS_

# include "def.h"
# include "bstream.h"
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

#ifdef _USE_OMP_
#include <omp.h>
#endif

# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_statistics_double.h>
# include <gsl/gsl_histogram.h>
# include <gsl/gsl_histogram2d.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_expint.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_sort_vector.h>
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
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <unistd.h>
# include <fftw3.h>
# include <gsl/gsl_dht.h>
#include "fftw_array.h"
using namespace std;


//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_r2c(int Nft, vector<real_prec>in, complex_prec *out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out);
//##################################################################################
//##################################################################################
//##################################################################################
/**
 * @brief This function builds the kernel from the ratio
 * @returns Container Bam::Kernel
 */
void do_fftw_3d(int Nft, bool direction, complex_prec *in, complex_prec *out);


