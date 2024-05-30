#ifndef __ANGULAR_POWER_SPECTRUM__
#define __ANGULAR_POWER_SPECTRUM__

// CLASS TO COMPUTE LINEAR AND NON LINEAR MATTER angular POWER SPECTRUM
// BASED ON THE EISENSTEIN AND HU FITTING FORMULAE
// AND THE HALO-FIT FROM SMITH ET AL.


# include <iostream>
# include <math.h>
# include <cmath>
# include <cctype>
# include <stdio.h>
# include <fstream>
# include <omp.h>
# include "CosmologicalFunctions.h"
# include "NumericalMethods.h"
# include "PowerSpectrumTH.h"
# include "FileOutput.h"


using namespace std;


struct cl_params{
  real_prec r;
  real_prec k;
  int l;
  real_prec zmean;
  real_prec width;
  string wtype;
  vector<gsl_real>F;
  vector<gsl_real>rv;
  vector<gsl_real>zv;
  vector<gsl_real>gv;
  vector<gsl_real>Hv;
  vector<gsl_real>dn;
  vector<gsl_real>dz;
  vector< vector<gsl_real > > sBessel;
  vector<gsl_real> sxBessel;
  s_CosmologicalParameters scp;
};





class AngularPowerSpectrum{
  
 private:
  
  int nz;
  
  static gsl_real idndz(gsl_real, void *);
  static gsl_real iKernel(gsl_real, void *);    
  static gsl_real iCl(gsl_real, void *);
  static gsl_real iCl_limber(gsl_real, void *);
  
  real_prec k_min, k_max, z_min, z_max;
  int nlbins;
  int L_min,L_max;
  int L_max_meas;
  PowerSpectrum PS;
  


 public:
  
  AngularPowerSpectrum(){}
  ~AngularPowerSpectrum(){}
  
  AngularPowerSpectrum(int nlb, int nn,real_prec _kmin, real_prec _kmax,int _lmin, int _lmax, int _lmaxmeas){
    nlbins=nlb;
    nz=nn;
    k_min=_kmin;
    k_max=_kmax;
    L_min=_lmin;
    L_max=_lmax;
    L_max_meas=_lmaxmeas;

  }


  FileOutput Fmi;
  FileOutput Fmd;

  gsl_real pk_normalization;
  gsl_real normal;
  vector<gsl_real>rv;
  vector<gsl_real>zv;
  vector<gsl_real>gv;
  vector<gsl_real>Hv;
  vector<gsl_real>dn;
  vector<gsl_real>dz;


  vector<gsl_real>XX_z;
  vector<gsl_real>WW_z;

  vector<gsl_real>XX_k;
  vector<gsl_real>WW_k;
  
  vector< vector<gsl_real > > Bessel;
  vector<gsl_real> xBessel;
  vector< vector<gsl_real> > R;
  
  gsl_integration_glfixed_table *wf;
  gsl_integration_glfixed_table *wfd;
  
  
  vector<real_prec>cl_model;
  vector<real_prec>cl_meas;
  vector<real_prec>sigma_cl;
  
  void get_cl_meas(string);
  void get_dndz(string);
  //  void Pk_normalization(s_CosmologicalParameters);
  void compute_int_table();
  void free_gsl_table();
  real_prec dndz(real_prec, void *);
  real_prec window(real_prec, void *);
  real_prec get_Kernel(void*);
  real_prec get_cl_exact(void*);
  gsl_real get_cl_ilimber(void*);
  void get_cl_limber(cl_params);
  void set_mixingM(string);
  void get_normal(void*);
  void get_cosmo(Cosmology , s_CosmologicalParameters);
  void get_bessel();
  
};

#endif  
