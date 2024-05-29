# include "../Headers/CorrelationFunction.h"

// *******************************************************************
// *******************************************************************

real_prec CorrelationFunctionTH::Linear_Matter_Correlation_Function(s_CosmologicalParameters *scp, real_prec r){        /*this k comes in h/Mpc */

  real_prec kmin=(scp->kmin_int);
  real_prec factor=0.5/pow(M_PI,2);
  return factor*gsl_integration_sin_kernel_lowlimit_inf2(i_Linear_Matter_Correlation_Function, (void *)scp, r, kmin, this->w,this->wc )/r;
}


// *******************************************************************
gsl_real CorrelationFunctionTH::i_Linear_Matter_Correlation_Function(gsl_real k, void *p){        /*this k comes in h/Mpc */
  // Integrand related to the power spectrum in the correlation function
  // xi= (1/2pi²) / r * int dk k²(sin(kr)/k)P(k) se integra en k directamente
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps(*scp);
  return k*Ps.Linear_Matter_Power_Spectrum(k);
}

// *******************************************************************
// *******************************************************************

real_prec CorrelationFunctionTH::Non_Linear_Matter_Correlation_Function_Halo_Fit(s_CosmologicalParameters *scp, real_prec r){        /*this k comes in h/Mpc */
  real_prec kmin=(scp->kmin_int);
  //  real_prec kmax=(scp->kmax_int);
  real_prec factor=0.5/pow(M_PI,2);
  return factor*gsl_integration_sin_kernel_lowlimit_inf(i_Non_Linear_Matter_Correlation_Function_Halo_Fit, (void *)scp, r, kmin)/r;
}
// *******************************************************************

gsl_real CorrelationFunctionTH::i_Non_Linear_Matter_Correlation_Function_Halo_Fit(gsl_real k, void *p){        /*this k comes in h/Mpc */
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps(*scp);
  vector<gsl_real> v_k_ps = scp->v_k_ps;
  vector<gsl_real> v_nl_power_spectrum = scp->v_nl_power_spectrum;
  return k*gsl_inter_new(v_k_ps, v_nl_power_spectrum, k);
}

// *******************************************************************
// *******************************************************************

real_prec CorrelationFunctionTH::Galaxy_Correlation_Function_1h_ss(s_CosmologicalParameters *scp, real_prec r){        /*this k comes in h/Mpc */
  real_prec kmin=(scp->kmin_int);
  real_prec kmax=(scp->kmax_int);
  real_prec factor=0.5/pow(M_PI,2);
  //return factor*gsl_integration_sin_kernel_lowlimit_inf(i_Galaxy_Correlation_Function_1h_ss, (void *)scp, r, kmin)/r;
  return factor*gsl_integration_sin_kernel(i_Galaxy_Correlation_Function_1h_ss, (void *)scp, r, kmin,kmax)/r;
}

// *******************************************************************
// *******************************************************************


gsl_real CorrelationFunctionTH::i_Galaxy_Correlation_Function_1h_ss(gsl_real lk, void *p){        /*this k comes in h/Mpc */
 /*this k comes in h/Mpc */
  real_prec k=lk;//pow(10,lk);
  real_prec jacobian=1.;//log(10)*k;
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  vector<gsl_real> v_k_ps = scp->v_k_ps;
  vector<gsl_real> v_galaxy_power_spectrum_1h_ss = scp->v_galaxy_power_spectrum_1h_ss;
  return jacobian*k*gsl_inter_new(v_k_ps,v_galaxy_power_spectrum_1h_ss, k);
}


// *******************************************************************
// *******************************************************************

real_prec CorrelationFunctionTH::Galaxy_Correlation_Function_1h_sc(s_CosmologicalParameters *scp, real_prec r){        /*this k comes in h/Mpc */
  real_prec kmin=(scp->kmin_int);
  real_prec kmax=(scp->kmax_int);
  real_prec factor=0.5/pow(M_PI,2);
  //return factor*gsl_integration_sin_kernel_lowlimit_inf(i_Galaxy_Correlation_Function_1h_sc, (void *)scp, r, kmin)/r;
  return factor*gsl_integration_sin_kernel(i_Galaxy_Correlation_Function_1h_sc, (void *)scp, r, kmin, kmax)/r;
}

gsl_real CorrelationFunctionTH::i_Galaxy_Correlation_Function_1h_sc(gsl_real lk, void *p){        /*this k comes in h/Mpc */
  real_prec k=lk;//pow(10,lk);
  real_prec jacobian=1.;//log(10)*k;
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  vector<gsl_real> v_k_ps = scp->v_k_ps;
  vector<gsl_real> v_galaxy_power_spectrum_1h_sc = scp->v_galaxy_power_spectrum_1h_sc;
  return jacobian*k*gsl_inter_new(v_k_ps,v_galaxy_power_spectrum_1h_sc, k);
}

// *******************************************************************
// *******************************************************************

real_prec CorrelationFunctionTH::Galaxy_Correlation_Function_2h(s_CosmologicalParameters *scp, real_prec r){        /*this k comes in h/Mpc */
  CorrelationFunctionTH Cf;
  real_prec kmin=(scp->kmin_int);
  real_prec kmax=(scp->kmax_int);
  real_prec factor=0.5/pow(M_PI,2);
  //return factor*gsl_integration_sin_kernel_lowlimit_inf(i_Galaxy_Correlation_Function_2h, (void *)scp, r, kmin)/r;
  return factor*gsl_integration_sin_kernel(i_Galaxy_Correlation_Function_2h, (void *)scp, r, kmin, kmax)/r;
}

gsl_real CorrelationFunctionTH::i_Galaxy_Correlation_Function_2h(gsl_real lk, void *p){        /*this k comes in h/Mpc */
 /*this k comes in h/Mpc */
  real_prec k=lk;//pow(10,lk);
  real_prec jacobian=1.;//log(10)*k;
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps(*scp);
  vector<gsl_real> v_k_ps = scp->v_k_ps;
  vector<gsl_real> v_galaxy_power_spectrum_2h = scp->v_galaxy_power_spectrum_2h;
  return jacobian*k*gsl_inter_new(v_k_ps,v_galaxy_power_spectrum_2h, k);

}

