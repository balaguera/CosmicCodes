//////////////////////////////////////////////////////////
/**
 * @file PowerSpectrumTH.cpp
 * @author Andres Balaguera-Antolínez (ABA)
 * @version 1.0
 * @author Andres Balaguera Antolinez
 * @date  2012-2023
 */
//////////////////////////////////////////////////////////
# include "../headers/PowerSpectrumTH.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::compute_int_table_k_mu(real_prec kmin_integration, real_prec kmax_integration, int nss_k, int nss_mu){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW.resize(nss_k);
  this->XX.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(kmin_integration)),static_cast<gsl_real>(log10(kmax_integration)),this->wfd,this->XX,this->WW);
  // INtegrate with respect to mu
  wf =  gsl_integration_glfixed_table_alloc (nss_mu);
  this->WW_mu.resize(nss_mu);
  this->XX_mu.resize(nss_mu);
  gsl_get_GL_weights(-1.0, 1.0,this->wf,this->XX_mu,this->WW_mu);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::compute_int_table_mass(real_prec M_min_integration, real_prec M_max_integration, int nss_k){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW_Mass.resize(nss_k);
  this->XX_Mass.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(1.1*log10(M_min_integration)),static_cast<gsl_real>(0.9*log10(M_max_integration)),this->wfd,this->XX_Mass,this->WW_Mass);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_interpolated(real_prec k){        /*this k comes in h/Mpc */
   return gsl_inter_new(this->s_cosmo_pars.v_k_ps, this->s_cosmo_pars.v_lin_power_spectrum, k);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_PT(real_prec k){
    real_prec g   = this->s_cosmo_pars.growth_factor;
    // The growth factor is factorized out and multiplied at the end. Hence we pass here 1.0 in the third argument
    real_prec P1  = exp(-0.5*pow(k/this->s_cosmo_pars.kstar,2))*this->Linear_Matter_Power_Spectrum_z(k,1.0);
    real_prec P2 = this->s_cosmo_pars.Amc*this->P1loop(k);
     return pow(g,2)*P1+ pow(g,4)*P2;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::P1loop(real_prec k){
    s_aux<PowerSpectrum> ssa;
    ssa.kaux=k;
    ssa.scp_a=&this->s_cosmo_pars;
    ssa.WW_mu=this->WW_mu;
    ssa.XX_mu=this->XX_mu;
    ssa.kmax_int=this->s_cosmo_pars.kmax_int;
    ssa.kmin_int=this->s_cosmo_pars.kmin_int;
    return gsl_integration(iP1loop, (void *)&ssa,this->XX, this->WW)/(2.*M_PI*M_PI);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::iP1loop(gsl_real lq, void *p){
   struct s_aux<PowerSpectrum> * saux = (struct s_aux<PowerSpectrum> *)p;
   struct s_CosmologicalParameters * scpa= saux->scp_a;
   PowerSpectrum Ps;
   Ps.set_cosmo_pars(*scpa);
   real_prec q=pow(10,lq);
   // Need to redefine here the aux structure to pass it at iFkernel
   s_aux<PowerSpectrum> ssaa;
   ssaa.kaux=saux->kaux; //k
   ssaa.raux=q; //q
   ssaa.scp_a=scpa;  //sca cosmological parameters
   ssaa.kmax_int=scpa->kmax_int;
   ssaa.kmin_int=scpa->kmin_int;
   //-------------------------
   // INtegral with respect to mu = k*q, the limits are set according to the kmin and kmax of the |k| intergration, see https://arxiv.org/pdf/astro-ph/0609547.pdf
  // THis sewction can be avoided if we are sure that the lñimits in mu are -1 and 1. As there is a kmin and max, these limites change
  real_prec kq_min=sqrt((saux->kaux*saux->kaux+q*q-scpa->kmax_int*scpa->kmax_int)/(2*q*saux->kaux));
  real_prec kq_max=sqrt((saux->kaux*saux->kaux+q*q-scpa->kmin_int*scpa->kmin_int)/(2*q*saux->kaux));
  real_prec mu_max=min(static_cast<float>(1.0), static_cast<float>(kq_max));
  real_prec mu_min=max(static_cast<float>(-1.0), static_cast<float>(kq_min));
   gsl_integration_glfixed_table *wf;
   wf =  gsl_integration_glfixed_table_alloc (30);
   gsl_get_GL_weights(mu_min, mu_max,wf,saux->XX_mu,saux->WW_mu);
   //-------------------------
   real_prec IF= gsl_integration(iFkernel, (void *)&ssaa, saux->XX_mu, saux->WW_mu);
   real_prec power=Ps.Linear_Matter_Power_Spectrum_z(q, 1.0);
   return static_cast<gsl_real>((log(10.0)*q)*q*q*power*IF);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::iFkernel(gsl_real mu, void *p){
    struct s_aux<PowerSpectrum> * saux = (struct s_aux<PowerSpectrum> *)p;
    struct s_CosmologicalParameters * scpa= saux->scp_a;
    PowerSpectrum Ps;
    Ps.set_cosmo_pars(*scpa);
    real_prec k = saux->kaux;
    real_prec q = saux->raux;
  //  real_prec kq= q*q+k*k>=2.0*k*q*mu ? q*q+k*k-2.0*k*q*mu : 0;
    real_prec kq=q*q+k*k-2.0*k*q*mu;
  //  real_prec kkk= q*q+k*k-2.0*k*q*mu;
    real_prec F2 = (5./7.)+0.5*((mu*k-q)/q)*(1.+q*q/kq)+(2./7.)*pow(mu*k-q,2)/kq;
    // If we are computing auto power spectrum, Pss has the kernel F**2. If it is cross power, F2 .
//    real_prec F2d = F2*F2;// only for auto power
    // The Plin is evaluated at z = 0  and extrapolated with the growth factor
    real_prec power=Ps.Linear_Matter_Power_Spectrum_z(sqrt(kq),1.0);
    return static_cast<gsl_real>(F2*power);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Linear_Matter_Power_Spectrum(real_prec k){        /*this k comes in units of h/Mpc */
  // Linear matter power spectrum 
    real_prec baryon_piece, cdm_piece;
    real_prec kk=k*this->s_cosmo_pars.hubble;// convert to units of 1/Mpc
    real_prec tf_thisk = TFfit_onek(kk, baryon_piece, cdm_piece);
    real_prec Tfunc = true==this->s_cosmo_pars.use_wiggles ? pow(fabs(tf_thisk),2.): pow(fabs(cdm_piece),2.) ;
    return this->s_cosmo_pars.pk_normalization*Tfunc*Primordial_Matter_Power_Spectrum(k)*pow(this->s_cosmo_pars.growth_factor,2);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_z(real_prec k, real_prec gr){        /*this k comes in h/Mpc */
  real_prec tf_thisk, baryon_piece, cdm_piece;
  tf_thisk = TFfit_onek(k*this->s_cosmo_pars.hubble, baryon_piece, cdm_piece);
  real_prec Tfunc = true==this->s_cosmo_pars.use_wiggles ? pow(fabs(tf_thisk),2.): pow(fabs(cdm_piece),2.) ;
  return this->s_cosmo_pars.pk_normalization*Tfunc*this->Primordial_Matter_Power_Spectrum(k)*gr*gr;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Q_Model_Matter_Power_Spectrum(real_prec k){        /*this k comes in h/Mpc */
  real_prec ans=0;
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec alpha_s=this->s_cosmo_pars.alpha_s;
  real_prec A_ps=this->s_cosmo_pars.A_PS;
  real_prec Q_ps=this->s_cosmo_pars.Q_PS;
  tf_thisk = TFfit_onek(k*this->s_cosmo_pars.hubble, baryon_piece, cdm_piece);
  if(true==this->s_cosmo_pars.use_wiggles)
    ans=Primordial_Matter_Power_Spectrum(k)*pow(fabs(tf_thisk),2.);
  else if(!this->s_cosmo_pars.use_wiggles)
    ans=Primordial_Matter_Power_Spectrum(k)*pow(fabs(cdm_piece),2.);
  return this->s_cosmo_pars.pk_normalization*ans*(1.+k*Q_ps)/(1.+A_ps*k)*pow(this->s_cosmo_pars.growth_factor,2);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Primordial_Matter_Power_Spectrum(real_prec k){        /*this k comes in h/Mpc */
  return pow(k,this->s_cosmo_pars.spectral_index);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Linear_Matter_Power_Spectrum_NW(real_prec k){        /*this k comes in h/Mpc */
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec alpha_s=this->s_cosmo_pars.alpha_s;
  tf_thisk =  TFfit_onek(k*this->s_cosmo_pars.hubble, baryon_piece, cdm_piece);
  real_prec prim= Primordial_Matter_Power_Spectrum(k)*pow(fabs(cdm_piece),2.);
  return this->s_cosmo_pars.pk_normalization*prim*pow(this->s_cosmo_pars.growth_factor,2);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec  PowerSpectrum ::Linear_Matter_Power_Spectrum_DW(real_prec k){        /*this k comes in h/Mpc */
  real_prec G=exp(-0.5*pow(k/(this->s_cosmo_pars.kstar),2));
  return Linear_Matter_Power_Spectrum(k)*G+Linear_Matter_Power_Spectrum_NW(k)*(1.-G);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec  PowerSpectrum ::Linear_Matter_Power_Spectrum_G_NW(real_prec k){        /*this k comes in h/Mpc */
  return Linear_Matter_Power_Spectrum_NW(k)*exp(-0.5*pow(k/(this->s_cosmo_pars.kstar),2));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum ::normalization(real_prec &nm){
  if(true==this->s_cosmo_pars.use_wiggles)
    nm=pow(this->s_cosmo_pars.sigma8,2)/gsl_integration(fun,(void *)&this->s_cosmo_pars, -4,5);
  else 
    nm=pow(this->s_cosmo_pars.sigma8,2)/gsl_integration(fun_nw,(void *)&this->s_cosmo_pars, -4,5);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum ::normalization(){
  real_prec nm=0;
  if(this->s_cosmo_pars.use_wiggles==true)
    nm=pow(this->s_cosmo_pars.sigma8,2)/gsl_integration(fun,(void *)&this->s_cosmo_pars, -7,7);
  else
    nm=pow(this->s_cosmo_pars.sigma8,2)/gsl_integration(fun_nw,(void *)&this->s_cosmo_pars, -7,7);
  return nm;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum ::fun(gsl_real lk, void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps;
  Ps.set_cosmo_pars(*scp);
  Ps.TFset_parameters();
  real_prec tf_thisk, baryon_piece, cdm_piece;
  real_prec k=pow(10,lk);
  tf_thisk =  Ps.TFfit_onek(k*scp->hubble, baryon_piece, cdm_piece);
  return static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Primordial_Matter_Power_Spectrum(k)*pow(tf_thisk,2)*pow(window(k,scp->RR),2));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum ::fun_nw(gsl_real lk, void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps;
  Ps.set_cosmo_pars(*scp);
  Ps.TFset_parameters();
  real_prec k=pow(10,lk);
  real_prec tfnw=Ps.TFnowiggles(k);
  return static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Primordial_Matter_Power_Spectrum(k)*pow(tfnw,2)*pow(window(k,scp->RR),2));//*pow(this->s_cosmo_pars.growth_factor,2);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Compute Kastar at z=0
void  PowerSpectrum::kstar_integral(real_prec &ksta){
  ksta=pow((real_prec)((1./(6.*M_PI*M_PI))*gsl_integration(Power_Spectrum_i,(void *)&this->s_cosmo_pars,log10(this->s_cosmo_pars.kmin_int),log10(this->s_cosmo_pars.kmax_int))),(real_prec)-0.5);
}
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::Power_Spectrum_i(gsl_real lk, void *p){        /*this k comes in h/Mpc */
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps;
  Ps.set_cosmo_pars(*scp);
  Ps.TFset_parameters();
  real_prec k=pow(10,lk);
  return static_cast<gsl_real>((log(10.0)*k)*Ps.Linear_Matter_Power_Spectrum_z(k,1.0));
}
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_Halo_Fit(real_prec k){
  real_prec ql, p, p_dw;
  halo_fit(k,ql,p,p_dw);
  return p;
}
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::halo_fit(real_prec k, real_prec &ql,real_prec &pp,real_prec &pp_dw){
  real_prec wde_eos=this->s_cosmo_pars.wde_eos;
  real_prec fac   = (1./(2.0*M_PI*M_PI))*pow(k,3);
  real_prec dh,dh_dw, dql, dql_dw;
  Cosmology cCf(this->s_cosmo_pars);
  real_prec Omz=cCf.omega_matter(this->s_cosmo_pars.cosmological_redshift);
  real_prec Omv=cCf.omega_dark_energy(this->s_cosmo_pars.cosmological_redshift);
  real_prec y     = k/(this->s_cosmo_pars.knl_hf);
  real_prec f1a    = pow(Omz,-0.0732);  /// Omega matter aca es al redhift, OJO
  real_prec f2a    = pow(Omz,-0.1423);
  real_prec f3a    = pow(Omz,+0.0725);
  real_prec f1b    = pow(Omz,-0.0307);  /// Omega matter aca es al redhift, OJO
  real_prec f2b    = pow(Omz,-0.0585);
  real_prec f3b    = pow(Omz,+0.0743);
  real_prec frac =Omv/(1.-Omz);
  real_prec f1=frac*f1b+(1-frac)*f1a;
  real_prec f2=frac*f2b+(1-frac)*f2a;
  real_prec f3=frac*f3b+(1-frac)*f3a;
  real_prec kmin  = log10(this->s_cosmo_pars.kmin_int);
  real_prec kmax  = log10(this->s_cosmo_pars.kmax_int);
  /*for the wiggled power spectrum*/
  real_prec rnl_hf=this->s_cosmo_pars.rnl_hf;
  real_prec a, b, c, alpha, beta, gama, nu, mu;
  this->s_cosmo_pars.aux_var3=rnl_hf;
  real_prec integration_aux4=gsl_integration(fun_aux_halo_fit4,(void *)&this->s_cosmo_pars, kmin,kmax);
  real_prec integration_aux2=gsl_integration(fun_aux_halo_fit2,(void *)&this->s_cosmo_pars, kmin,kmax);
  real_prec cc    = 4.0*pow(rnl_hf,2)*integration_aux4+4*pow(rnl_hf,4)*pow(integration_aux2,2);
  real_prec neff  =-3.0+2.0*pow(rnl_hf,2.)*integration_aux2;
  hf_aux(neff,cc,wde_eos, Omv,a,b,c,alpha,beta,gama,mu,nu);
  real_prec dl  = fac*this->Linear_Matter_Power_Spectrum(k);
  dql   = dl*(pow(1.+dl, beta)/(1.+alpha*dl))*exp(-y/4.-y*y/8.);
  real_prec dhp   = a*pow(y,3.0*f1)/(1.+b*pow(y,f2)+pow(c*f3*y,3-gama));
  dh    = dhp/(1.+(mu/y)+nu*pow(y,-2));
  ql=dql/fac;
  pp=(dh+dql)/fac;
  pp_dw=(dh_dw+dql_dw)/fac;
}
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(real_prec k, real_prec z, real_prec g, real_prec h4, real_prec h2, real_prec kln){
  real_prec ql, p;
  halo_fit_z(k, z,g,h4, h2, kln, ql,p);
  return p;
}
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::halo_fit_z(real_prec k, real_prec z, real_prec gf, real_prec h4, real_prec h2, real_prec knl, real_prec &ql,real_prec &pp){
  real_prec fac   = (1./(2.0*M_PI*M_PI))*pow(k,3);
  Cosmology cCf(this->s_cosmo_pars);
  real_prec Omz=cCf.omega_matter(z);
  real_prec Omv=cCf.omega_dark_energy(z);
  real_prec rnl_hf= 1./knl;
  real_prec y     = k/knl;
  real_prec f1a    = pow(Omz,-0.0732);  /// Omega matter aca es al redhift, OJO
  real_prec f2a    = pow(Omz,-0.1423);
  real_prec f3a    = pow(Omz,+0.0725);
  real_prec f1b    = pow(Omz,-0.0307);  /// Omega matter aca es al redhift, OJO
  real_prec f2b    = pow(Omz,-0.0585);
  real_prec f3b    = pow(Omz,+0.0743);
  real_prec frac =  Omv/(1.-Omz); // =1 appears in the original paper of Takahashi. The original paper by R SMith suggest interpolation. CLASS uses frac like this:
  real_prec f1=frac*f1b+(1.-frac)*f1a;
  real_prec f2=frac*f2b+(1.-frac)*f2a;
  real_prec f3=frac*f3b+(1.-frac)*f3a;
 // Las cantidades h2 y h4 ya tienen incorporado el gf**2 respectivo, que viene desde el cl_model.
  real_prec cc    = 4.0*pow(rnl_hf,2)*h4+4.0*pow(rnl_hf,4)*pow(h2,2);
  real_prec neff  =-3.0+2.*pow(rnl_hf,2.)*h2;
  real_prec a, b, c, alpha, beta, gama, mu, nu;
  hf_aux(neff,cc,this->s_cosmo_pars.wde_eos, Omv,a,b,c,alpha,beta,gama,mu,nu);
  // Two-halo term Delta Q
  real_prec dl    = fac*this->Linear_Matter_Power_Spectrum_z(k,gf);
  real_prec dql   = dl*(pow(1.+dl, beta)/(1.+alpha*dl))*exp(-y/4.0-y*y/8.0);
  // One halo term Delta H
  real_prec dhp   = a*pow(y,3.0*f1)/(1.+b*pow(y,f2)+pow(c*f3*y, 3.-gama));
  real_prec dh    = dhp/(1.+nu*pow(y,-2));  // mu = 0 for the fit of Takahashi
  ql=dh/fac;
  pp=(dh+dql)/fac;
}
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::halo_fit_integrals(real_prec &inte4, real_prec &inte2 ){
  real_prec kmin  = log10(this->s_cosmo_pars.kmin_int);
  real_prec kmax  = log10(this->s_cosmo_pars.kmax_int);
  int Ni = 100;
  real_prec integ4, integ2;
  if (true==this->s_cosmo_pars.use_wiggles){
   inte4=gsl_integration(Ni,fun_aux_halo_fit4,(void *)&this->s_cosmo_pars, kmin,kmax);
   inte2=gsl_integration(Ni,fun_aux_halo_fit2,(void *)&this->s_cosmo_pars, kmin,kmax);
  }
  else{
      /*for the de-wigled power spectrum*/
    inte4=gsl_integration(Ni,fun_aux_halo_fit4_dw,(void *)&this->s_cosmo_pars,kmin,kmax);
    inte2=gsl_integration(Ni,fun_aux_halo_fit2_dw,(void *)&this->s_cosmo_pars, kmin,kmax);
   }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::hf_aux(real_prec index, real_prec cc,real_prec weos, real_prec Omv,real_prec &a, real_prec &b, real_prec &c, real_prec &alpha,real_prec &beta,real_prec &gama,real_prec &mu, real_prec &nu){
  /*Funciones auxiliares del halo fit*/
  a   = pow(10,1.5222+2.8553*index+2.3706*index*index+0.9903*index*index*index+0.2250*pow(index,4)-0.6038*cc+0.1749*Omv*(1+weos));
  b   = pow(10,-0.5642+0.5864*index+0.5716*index*index-1.5474*cc+0.2279*Omv*(1+weos));
  c   = pow(10,0.3698+2.0404*index+0.8161*index*index+0.5869*cc);
  alpha= abs(6.0835+1.3373*index-0.1959*index*index-5.5274*cc);
  beta = 2.0379-0.7354*index+0.3157*index*index+1.2490*pow(index,3)+0.3980*pow(index,4)-0.1682*cc;
  gama = 0.1971-0.0843*index+0.8460*cc;
  mu   = 0.0;         //pow(10,-3.5442+0.1908*index);
  nu   = pow(10,5.2105+3.6902*index);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PowerSpectrum::nl_scales_halo_fit(real_prec &knl_hf,real_prec &rnl_hf, vector<real_prec>&rr, vector<real_prec>&sums, bool silence){
  if(silence)
    So.message_screen("Computing non linear scales for halo fit using Halo Fit by S03 and revised by Takahashi et al 2012");
  int nr=rr.size();
  real_prec rinic;
  rinic = this->s_cosmo_pars.Om_cdm<0.07? 1e-10 : 1e-3  ;
  real_prec rfinal=2e2;
  // We have an issue here. If we use pragma omp for, we get the numbers, but not that fast.
  // If we use pragma omp parallel for, the code executes jobs in parallel through different threads
  // but the information encoded in the structure is messed up, so no good answer. Issue Open.
  // Solution: I created another structure (struct s_aux, in Type_def*h) with two members, the cosmological_parameters structure
  // and a real_prec. Then I define an object of this type withon the loop so I get sure that every index
  // has its own strucuture defined, as if it were private.
  fill(sums.begin(), sums.end(), 0);
  vector<gsl_real> sums_aux (sums.size(),0);
  vector<gsl_real> rr_aux (sums.size(),0);
#ifdef SINGLE_PREC
  for(int i=0;i<sums.size();++i)
    sums_aux[i]=static_cast<gsl_real>(sums[i]);
  for(int i=0;i<sums.size();++i)
    rr_aux[i]=static_cast<gsl_real>(rr[i]);
#else
  sums_aux=sums;
#endif
  for(int i=0;i<rr.size();++i)
  {
    real_prec Rr= pow(10, log10(rinic)+i*log10(rfinal/rinic)/((real_prec)nr-1.));
    s_aux<PowerSpectrum>ssa;
    ssa.raux=Rr;
    ssa.scp_a=&this->s_cosmo_pars;
    rr_aux[nr-1-i]=Rr;
    // Compute sigma**2 (R,z=0)
    sums_aux[nr-1-i]= gsl_integration(400,fun_aux_halo_fit,(void *)&ssa,log10(this->s_cosmo_pars.kmin_int),log10(this->s_cosmo_pars.kmax_int));
  }
  rnl_hf=gsl_inter_new(sums_aux,rr_aux,num_1);
  knl_hf=1./rnl_hf;
  silence=true;
  if(silence){
    s_aux<PowerSpectrum> ssb;
    ssb.raux=rnl_hf;
    ssb.scp_a=&this->s_cosmo_pars;
    real_prec ss_check=gsl_integration(500,fun_aux_halo_fit,(void *)&ssb,log10(this->s_cosmo_pars.kmin_int),log10(this->s_cosmo_pars.kmax_int));
    So.message_screen("\tCheck: Sigma(r) at r = Rnl",ss_check);
    So.message_screen("\tAbsolute error",100.0*abs(1- ss_check)," %");
  }
#ifdef SINGLE_PREC
  for(int i=0;i<sums.size();++i)
       sums[i]=static_cast<real_prec>(sums_aux[i]);
  for(int i=0;i<sums.size();++i)
       rr[i]=static_cast<real_prec>(rr_aux[i]);
#else
  sums=sums_aux;
  rr=rr_aux;
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::fun_aux_halo_fit(gsl_real lk, void *p){
   struct s_aux<PowerSpectrum> * saa = (struct s_aux<PowerSpectrum> *)p;
   s_CosmologicalParameters * scp = saa->scp_a ;
   PowerSpectrum Ps;
   Ps.set_cosmo_pars(*scp);
   Ps.TFset_parameters();
   real_prec k=pow(10,lk);
   real_prec r = saa->raux;
   real_prec power=Ps.Linear_Matter_Power_Spectrum(k);
   return static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*power*exp(-pow(k*r,2)));
 }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::fun_aux_halo_fit_dw(gsl_real lk, void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum Ps;
  Ps.set_cosmo_pars(*scp);
  Ps.TFset_parameters();
  real_prec k=pow(10,lk);
  real_prec r = scp->aux_var3;
  return  static_cast<gsl_real>((log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(k)*exp(-pow(k*r,2)));
 }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::fun_aux_halo_fit2(gsl_real lk, void *p){
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    PowerSpectrum Ps;
    Ps.set_cosmo_pars(*scp);
    Ps.TFset_parameters();
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum(k)*exp(-pow(k*r,2))*pow(k,2);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::fun_aux_halo_fit2_dw(gsl_real lk, void *p){
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    PowerSpectrum Ps;
    Ps.set_cosmo_pars(*scp);
    Ps.TFset_parameters();
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return  (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(k)*exp(-pow(k*r,2))*pow(k,2);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
  gsl_real PowerSpectrum::fun_aux_halo_fit4(gsl_real lk, void *p){
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    PowerSpectrum Ps;
    Ps.set_cosmo_pars(*scp);
    Ps.TFset_parameters();
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    return (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum(k)*exp(-pow(k*r,2))*pow(k,2)*(1.-pow(r*k,2));
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
  gsl_real PowerSpectrum::fun_aux_halo_fit4_dw(gsl_real lk, void *p){
    struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
    PowerSpectrum Ps;
    Ps.set_cosmo_pars(*scp);
    Ps.TFset_parameters();
    real_prec k=pow(10,lk);
    real_prec r = scp->aux_var3;
    real_prec ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*Ps.Linear_Matter_Power_Spectrum_DW(k)*exp(-pow(k*r,2))*k*k*(1.-pow(r*k,2));
    return  static_cast<gsl_real>(ans);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Galaxy_power_spectrum_h1_ss(real_prec k, real_prec z){
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.MASS_BIAS=this->v_halo_mass_bias;
    sA1.s_cp=&this->s_cosmo_pars;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration(i_Galaxy_power_spectrum_h1_ss, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 gsl_real PowerSpectrum::i_Galaxy_power_spectrum_h1_ss(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod;
    DensityProfiles Dp;
    real_prec M=pow(10,m);
    real_prec z=sA1->aux_z;
    real_prec k=sA1->aux_k;
    real_prec jacobian=(log(10.0)*M);
    real_prec uden=Dp.density_k(k,M,z, scp);
    return static_cast<gsl_real>(jacobian*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m)*pow(Shod.SATELLITE(M)*uden,2));
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
  real_prec PowerSpectrum::Galaxy_power_spectrum_h1_sc(real_prec k, real_prec z){
    this->s_cosmo_pars.aux_var4=k;
    this->s_cosmo_pars.aux_var2=z;
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.s_cp=&this->s_cosmo_pars;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration(i_Galaxy_power_spectrum_h1_sc, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real PowerSpectrum::i_Galaxy_power_spectrum_h1_sc(gsl_real m, void *p)
  {
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod(*scp);
    DensityProfiles Dp;
//    MARKS Smark;
    real_prec M=pow(10,m);
    real_prec z=sA1->aux_z;
    real_prec k=sA1->aux_k;
    real_prec jacobian=(log(10.0)*M);
    real_prec uden=Dp.density_k(k,M,z,scp);
    return static_cast<gsl_real>(jacobian*2.0*Shod.CENTRAL(M)*Shod.SATELLITE(M)*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m)*uden);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::Galaxy_matter_bias(real_prec k, real_prec z){
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION = this->v_mass_function;
    sA1.MASS_BIAS = this->v_halo_mass_bias;
    sA1.s_cp=&this->s_cosmo_pars;
    sA1.aux_k=k;
    sA1.aux_z=z;
    return gsl_integration(i_Galaxy_matter_bias, (void *)&sA1,this->XX_Mass, this->WW_Mass);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 gsl_real PowerSpectrum::i_Galaxy_matter_bias(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod(*scp);
    DensityProfiles Dp;
    real_prec M=pow(10,m);
    real_prec k=sA1->aux_k;
    real_prec z=sA1->aux_z;
    real_prec jacobian=(log(10.0)*M);
    real_prec uk=Dp.density_k(k,M,z,scp);
    return static_cast<gsl_real>(jacobian*(Shod.CENTRAL(M)+Shod.SATELLITE(M)*uk)*gsl_inter_new(sA1->MASS, sA1->MASS_BIAS,m)*gsl_inter_new(sA1->MASS, sA1->MASS_FUNCTION,m));
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec PowerSpectrum::mean_galaxy_number_density(real_prec redshift){
    this->s_cosmo_pars.aux_var3=redshift;
    //Integrate from the value mmin_hod
    return gsl_integration(i_mean_galaxy_number_density,(void *)&this->s_cosmo_pars,this->XX_Mass,this->WW_Mass);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 real_prec PowerSpectrum::mean_galaxy_number_density(){
    A1 sA1;
    sA1.MASS=this->v_mass;
    sA1.MASS_FUNCTION=this->v_mass_function;
    sA1.s_cp=&this->s_cosmo_pars;
    //Integrate from the value mmin_hod
    return gsl_integration(i_mean_galaxy_number_density,(void *)&sA1,this->XX_Mass,this->WW_Mass);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 gsl_real PowerSpectrum::i_mean_galaxy_number_density(gsl_real m, void *p){
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *scp = sA1->s_cp;
    HOD Shod(*scp);
    real_prec M=pow(10,m);
    real_prec jacobian=log(10)*M;
    return static_cast<gsl_real>(jacobian*gsl_inter_new(sA1->MASS,sA1->MASS_FUNCTION,m)*(Shod.SATELLITE(M)+Shod.CENTRAL(M)));
 }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/* ------------------------ FITTING FORMULAE ROUTINES ------------------- */
  /* There are two routines here.  TFset_parameters() sets all the scalar
  parameters, while TFfit_onek() calculates the transfer function for a
  given wavenumber k.  TFfit_onek() may be called many times after a single
  call to TFset_parameters() */
  /* Global variables -- We've left many of the intermediate results as
  global variables in case you wish to access them, e.g. by declaring
  them as extern variables in your main program. */
  /* Note that all internal scales are in Mpc, without any Hubble constants! */
 void PowerSpectrum ::TFset_parameters()
  {
    Cosmology cf(this->s_cosmo_pars);
    real_prec omhh = this->s_cosmo_pars.Om_matter*pow(this->s_cosmo_pars.hubble,2);
    real_prec obhh = this->s_cosmo_pars.Om_baryons*pow(this->s_cosmo_pars.hubble,2);
    if (this->s_cosmo_pars.f_baryon<=0.0 || omhh <=0.0)
     {
       cout<<"f_baryon= "<<this->s_cosmo_pars.f_baryon<<" , omega0hh= "<<omhh<<endl;
       throw std::invalid_argument("Illegal input in cosmological parameters");
     }
    if (this->s_cosmo_pars.Tcmb<=0.0) 
      this->s_cosmo_pars.Tcmb=2.728;	/* COBE FIRAS */
    this->k_equality = cf.k_equality();
    this->sound_horizon = cf.sound_horizon(); 
    this->k_silk = cf.k_Silk();  
    real_prec alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    real_prec alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    this->alpha_c = pow(alpha_c_a1,-this->s_cosmo_pars.f_baryon)*pow(alpha_c_a2,-pow(this->s_cosmo_pars.f_baryon,3));
    real_prec beta_c_b1 = 0.944/(1.+pow(458.0*omhh,-0.708));
    real_prec beta_c_b2 = pow(0.395*omhh, -0.0266);
    this->beta_c = 1.0/(1+beta_c_b1*(pow(1-this->s_cosmo_pars.f_baryon, beta_c_b2)-1));
    real_prec y = cf.redshift_equality()/(1.+ cf.drag_redshift());
    real_prec alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    this->alpha_b = 2.07*this->k_equality*this->sound_horizon*pow(1+cf.R_drag(),-0.75)*alpha_b_G;
    this->beta_node = 8.41*pow(omhh, 0.435);  
    this->beta_b = 0.5+this->s_cosmo_pars.f_baryon+(3.-2.*this->s_cosmo_pars.f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);
    this->k_peak = 2.5*3.14159*(1+0.217*omhh)/this->sound_horizon;
    this->sound_horizon_fit = 44.5*log(9.83/omhh)/sqrt(1+10.0*pow(obhh,0.75));
    this->alpha_gamma = 1-0.328*log(431.0*omhh)*this->s_cosmo_pars.f_baryon + 0.38*log(22.3*omhh)*pow(this->s_cosmo_pars.f_baryon,2);
     return;
  }
////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 /* Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1.
 *tf_baryon, *tf_cdm -- Input value not used; replaced on output if
             the input was not NULL. */
  /* Output: Returns the value of the full transfer function fitting formula.
          This is the form given in Section 3 of Eisenstein & Hu (1997).
        *tf_baryon -- The baryonic contribution to the full fit.
        *tf_cdm -- The CDM contribution to the full fit. */
  /* Notes: Units are Mpc, not h^-1 Mpc. */
real_prec PowerSpectrum ::TFfit_onek(real_prec k, real_prec &tf_baryon, real_prec &tf_cdm)
  {
    if (k==0.0) 
    {
      if (&tf_baryon!=NULL) tf_baryon = 1.0;
      if (&tf_cdm!=NULL) tf_cdm = 1.0;
      return 1.0;
    }
    else
    {  
      real_prec q = k/13.41/this->k_equality;
      real_prec xx = k*this->sound_horizon;
      real_prec T_c_ln_beta = log(2.718282+1.8*this->beta_c*q);
      real_prec T_c_ln_nobeta = log(2.718282+1.8*q);
      real_prec T_c_C_alpha = 14.2/this->alpha_c + 386.0/(1+69.9*pow(q,1.08));
      real_prec T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));
      real_prec T_c_f = 1.0/(1.0+pow(xx/5.4,4));
      real_prec T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*(q*q))+(1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*(q*q));
      real_prec s_tilde = this->sound_horizon*pow(1+pow(this->beta_node/xx,3),-1./3.);
      real_prec xx_tilde = k*s_tilde;
      real_prec T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*(q*q));
      real_prec T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+pow(xx/5.2,2))+this->alpha_b*exp(-pow(k/this->k_silk,1.4)))/(1.+pow(this->beta_b/xx,3));
      real_prec T_full = this->s_cosmo_pars.f_baryon*T_b+ (1.0-this->s_cosmo_pars.f_baryon)*T_c;
            /* Now to store these transfer functions */
      if (&tf_baryon!=NULL) tf_baryon = T_b;
      if (&tf_cdm!=NULL) tf_cdm = T_c;
      return T_full;
    }
  }
////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
  real_prec PowerSpectrum ::TFnowiggles(real_prec k_hmpc)
  /* k_hmpc -- Wavenumber in units of (h Mpc^-1). */
  /* Output: The value of an approximate transfer function that captures the
  non-oscillatory part of a partial baryon transfer function.  In other words,
  the baryon oscillations are left out, but the suppression of power below
  the sound horizon is included. See equations (30) and (31).  */
  /* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      Cosmology cf(this->s_cosmo_pars);
      real_prec k = k_hmpc*this->s_cosmo_pars.hubble;	/* Convert to Mpc^-1 */
      real_prec omhh = this->s_cosmo_pars.Om_matter*this->s_cosmo_pars.hubble*this->s_cosmo_pars.hubble;
      if (this->s_cosmo_pars.Tcmb<=0.0) this->s_cosmo_pars.Tcmb=2.728;	/* COBE FIRAS */
      real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
      real_prec k_equality = 0.0746*omhh/pow(theta_cmb,2);
      real_prec q = k/13.41/k_equality;
      cf.set_f_baryon(1.0); // set fraction of baryons to unity in the s_cosmo_pars structure of Cosmology
      real_prec xx = k*cf.TFsound_horizon_fit();
      real_prec alpha_gamma = 1-0.328*log(431.0*omhh)*this->s_cosmo_pars.f_baryon + 0.38*log(22.3*omhh)*pow(this->s_cosmo_pars.f_baryon,2);
      real_prec gamma_eff = omhh*(alpha_gamma+(1-alpha_gamma)/(1+pow(0.43*xx,4)));
      real_prec q_eff = q*omhh/gamma_eff;
      real_prec T_nowiggles_L0 = log(2.0*2.718282+1.8*q_eff);
      real_prec T_nowiggles_C0 = 14.2 + 731.0/(1+62.5*q_eff);
      return T_nowiggles_L0/(T_nowiggles_L0+T_nowiggles_C0*pow(q_eff,2));
  }
  ////////////////////////////////////////////////////////////////////////////
  // ======================= Zero Baryon Formula ======================
  real_prec PowerSpectrum ::TFzerobaryon(real_prec k_hmpc)
  /* Input: omega0 -- CDM density, in units of critical density
        hubble -- Hubble constant, in units of 100 km/s/Mpc
        Tcmb -- Temperature of the CMB in Kelvin; Tcmb<=0 forces use of
              COBE FIRAS value of 2.728 K
        k_hmpc -- Wavenumber in units of (h Mpc^-1). */
  /* Output: The value of the transfer function for a zero-baryon universe. */
  /* Note: If you prefer to use wavenumbers in units of Mpc^-1, use hubble -> 1
  and omega0 -> omega0*hubble^2. */
  {
      real_prec k = k_hmpc*this->s_cosmo_pars.hubble;	/* Convert to Mpc^-1 */
      real_prec omhh = this->s_cosmo_pars.Om_matter*this->s_cosmo_pars.hubble*this->s_cosmo_pars.hubble;
      if (this->s_cosmo_pars.Tcmb<=0.0)this->s_cosmo_pars.Tcmb=2.728;	/* COBE FIRAS */
      real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
      real_prec k_equality = 0.0746*omhh/pow(theta_cmb,2);
      real_prec q = k/13.41/k_equality;
      real_prec T_0_L0 = log(2.0*2.718282+1.8*q);
      real_prec T_0_C0 = 14.2 + 731.0/(1+62.5*q);
      return T_0_L0/(T_0_L0+T_0_C0*q*q);
  }
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////
  // These functions are also in Statistics, but in order to use them in the
  // Cl analysis with omp, I had to copy them here
  real_prec PowerSpectrum::sigma_masa(real_prec m, real_prec z, s_CosmologicalParameters *scp){
    /* Sigma as a function of mass */
    this->s_cosmo_pars.aux_var1 = m;  //log10M
    this->s_cosmo_pars.aux_var2 = z;
    real_prec ans=gsl_integration(isigma,(void *)scp,this->XX,this->WW);
    return sqrt(ans);
  }
  ////////////////////////////////////////////////////////////////////////////
 gsl_real PowerSpectrum::isigma(gsl_real lk,void *p){  /* Integrand for sigma^2 */
    real_prec k=pow(10,lk);
    s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
    Cosmology cf(*s_cp);
    PowerSpectrum ps;
    ps.set_cosmo_pars(*s_cp);
    real_prec m= s_cp->aux_var1;
    real_prec z= s_cp->aux_var2;
    real_prec M=pow(10,m);
    real_prec ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum(k)*pow(window(k,cf.rr(M,z)),2);
    return ans;
  }
 ////////////////////////////////////////////////////////////////////////////
  real_prec PowerSpectrum::mass_function_D(real_prec m, real_prec z,s_CosmologicalParameters  *scp)
  {
    Cosmology Cf(*scp);
    MASS_BIAS_FUNCTIONS mb;
    real_prec deltac = Cf.critical_overdensity(z);
    real_prec Smasa=gsl_inter_new(this->v_mass, this->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
  //  real_prec Smasa=sigma_masa(M,z,scp);
    real_prec M=pow(10,m);
    real_prec nu=pow(deltac/Smasa,2);
    real_prec mean_matter_density=Cf.mean_matter_density(z);
    real_prec com_density=(mean_matter_density)*pow(1+z,(real_prec)-3.0);
    A1 sA1;
    sA1.aux_m = m;
    sA1.s_cp=scp;
    sA1.aux_z = z;
    real_prec dsdr=gsl_integration(dsigma_dR,(void *)&sA1, this->XX,this->WW);
    return (com_density/M)*(mb.mass_function(nu,z, (void *)scp))*pow(Smasa,-2)*(Cf.rr(M,z)/(3.*M))*dsdr;
  }
  ////////////////////////////////////////////////////////////////////////////
  void PowerSpectrum::mass_function_M_Z(vector<real_prec>XZ, s_CosmologicalParameters  *scp){
   for(int i=0; i<this->v_mass.size();i++){
      this->v_mass[i]=(log10(this->s_cosmo_pars.M_min_effective)+i*(log10(this->s_cosmo_pars.M_max_effective)-log10(this->s_cosmo_pars.M_min_effective))/this->v_mass.size());
   }
   for(int i=0; i<this->v_mass.size();i++){  //ESTO HACE EL SIGMA MASS DE HM SCALES INNECESARIO
      this->v_sigma_mass[i]=this->sigma_masa(this->v_mass[i],0.0,scp);
   }
   this->MASS_FUNCTION_M_Z.resize(XZ.size());
   for(int i=0;i<XZ.size();++i)this->MASS_FUNCTION_M_Z[i].resize(this->v_mass.size(),0);
   this->MASS_BIAS_M_Z.resize(XZ.size());
   for(int i=0;i<XZ.size();++i)this->MASS_BIAS_M_Z[i].resize(this->v_mass.size(),0);
   for(int i=0;i<XZ.size();++i)
    {
      real_prec gg=gsl_inter_new(this->s_cosmo_pars.zv,this->s_cosmo_pars.gv,XZ[i]);
      for(int j=0; j<this->v_mass.size();j++)
       {
         this->v_sigma_mass[j]*=gg;
         this->MASS_FUNCTION_M_Z[i][j]=mass_function_D(this->v_mass[j],XZ[i],scp);
         this->MASS_BIAS_M_Z[i][j]=bias(this->v_mass[j],XZ[i],scp);
       }
    }
  }
  ////////////////////////////////////////////////////////////////////////////
  gsl_real PowerSpectrum::dsigma_dR(gsl_real lk, void *p){  /* integrand to calculate the derivative of sigma^2 with respect to  R(M)*/
    struct A1 * sA1= (struct A1 *)p;
    s_CosmologicalParameters *s_cp = sA1->s_cp;
    PowerSpectrum ps;
    ps.set_cosmo_pars(*s_cp);
    Cosmology Cf(*s_cp);
    real_prec m= sA1->aux_m;
    real_prec M=pow(10,m);
    real_prec z= sA1->aux_z;
    real_prec k=pow(10,lk);
    real_prec r =Cf.rr(M,z);
    real_prec y=r*k;
    return   (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum_z(k, 1.0)*2.0*fabs(window(k,r))*k*fabs(((3.*pow(y,-2))*sin(y)-9.*pow(y,-4)*(sin(y)-y*cos(y))));
  }
  ////////////////////////////////////////////////////////////////////////////
  real_prec PowerSpectrum::bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
    Cosmology Cf(*scp);
    MASS_BIAS_FUNCTIONS mb;
    real_prec Smasa=gsl_inter_new(this->v_mass, this->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
    real_prec deltac = Cf.critical_overdensity(z);
    real_prec nu=pow(deltac/Smasa,2.);
    return mb.dm_h_bias(z,nu,(void *)scp);
  }

