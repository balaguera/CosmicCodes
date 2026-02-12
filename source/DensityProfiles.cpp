////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** 
 * @file DensityProfiles.cpp
 * @brief This file contains methods of density profiles
 * @author Andres Balaguera Antolinez 2017-2023
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../headers/DensityProfiles.h"
//////////////////////////////////////////////////////////
real_prec DensityProfiles::density_r(real_prec r){
  real_prec ans=0;
  if(this->s_cosmo_pars.density_profile=="nfw"){
    ans=this->rhos/((r/this->rs)*pow(1.+r/this->rs,2));
  }
  if(this->s_cosmo_pars.density_profile=="einasto"){
    ans=this->rhos*exp(-(2./this->alpha)*pow(r/this->rs,this->alpha)+(2./this->alpha));
  }  
  return ans;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::DensityProfile_NFW(real_prec r){
  return this->rhos/((r/this->rs)*pow(1.+r/this->rs,2));
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::DensityProfile_NFW_prob(real_prec r,real_prec rmin){
  real_prec normal=this->DensityProfile_NFW(rmin);
  return this->rhos/((r/this->rs)*pow(1.+r/this->rs,2))/normal;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::density_rc(real_prec r, real_prec M, real_prec z, void *cp){
  /*Profiles convolved with the c-M distribution
    Evaluates the integral with respect to the natural log of the concentration
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)cp;
  s_cp->aux_var1=r;
  s_cp->aux_var2=M;
  s_cp->aux_var3=z;

  // Ojo con este truco: como solo puedo pasar una estructura para las 
  // las funciones que se van a integrar, 
  // relleno esa estructura con los elementos que necesito
  // de una segunda estructura que hasta aca pude pasar como argumento
  return gsl_integration(idensity_rc, (void *)s_cp,-2,3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gsl_real DensityProfiles::idensity_rc(gsl_real c, void *p){
  /*Here we write explicitely the density profiles as a function of the concentration c
    in order to convolve these with with mass-concentration distribution function
    Integration is done in the log10(c)  of the concentration (the c-m dist is in log(c)!)
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology Cf(*s_cp);
  DensityProfiles Dp;
  PowerSpectrum Ps;
  c=pow(10,c);
  real_prec r = s_cp->aux_var1;
  real_prec M = s_cp->aux_var2;
  real_prec z = s_cp->aux_var3;
  real_prec dcontrast=Cf.density_contrast_top_hat(z);
  real_prec deltac=Cf.critical_overdensity(z);
  string prof=s_cp->density_profile;
  real_prec rvir,rs,alpha,f,nu;
  rvir=pow((real_prec)(3.*M/(4.*M_PI*(dcontrast)*Cf.mean_matter_density(z))),(real_prec)1./3.);
   if(prof=="nfw"){ 
    f=(1./3.)*(dcontrast)*(Cf.mean_matter_density(z))*pow(c,3)/(log(1+c)-(c/(1.+c)));
    f=f/((c*r/rvir)*pow(1.+c*r/rvir,2));
  }
  if(prof=="einasto"){ 
    //    nu=pow(deltac/gsl_inter_pointers(mass_array_p,sigma_mass_p,nn,m),2);

    vector<gsl_real>v_sigma_mass=s_cp->v_sigma_mass;
    vector<gsl_real>v_mass=s_cp->v_mass;
    nu=pow(deltac/gsl_inter_new(v_mass,v_sigma_mass, M),2);

    alpha=0.155+0.0095*nu;
    f=M*pow(alpha/2,(-3+alpha)/alpha)*exp(-2/alpha)*pow(rvir/c,-3)/(2.*M_PI*gsl_sf_gamma_inc_P(3./alpha,(2./alpha)*pow(c,alpha))*gsl_sf_gamma_inc_P((real_prec)3./alpha,(real_prec)(2./alpha)*c)*gsl_sf_gamma(3./alpha));
    f=f*exp(-(2./alpha)*pow(c*r/rvir,alpha)+(2./alpha));
  }
  
  f=(log(10)*c)*f*Dp.mass_concentration_dis(c,M,z, p);
  return static_cast<gsl_real>(f);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::density_k(real_prec k){
  /*Por defecto dejamos la transformada del NFW en su forma analitica, aunque la numerica se puede habilitar
    comentando el if(de_label==1) y sacando el  density_profile_fourier fuera del if(den_label==2), dejando solo un if(de_label==1)
    para llamar los parametros de NFW
  */
  string prof=this->s_cosmo_pars.density_profile;
  real_prec uk;
  if(prof=="nfw"){
    real_prec cia=gsl_sf_Ci(k*this->rs);
    real_prec sia=gsl_sf_Si(k*this->rs);
    real_prec cib=gsl_sf_Ci(k*this->rs*(this->concentration+1));
    real_prec sib=gsl_sf_Si(k*this->rs*(this->concentration+1));
    // analytical, behaves better
    uk = (4.*M_PI*rhos*pow(this->rs,3))*( sin(k*this->rs)*(sib-sia) -(sin(this->concentration*k*this->rs)/((1.+this->concentration)*k*this->rs))+cos(k*this->rs)*(cib-cia));
    // numerical estimate:
//      uk=gsl_integration_sin_kernel(idensity_fourier,(void *)s_cp, k, 1e-5, rvir);
  }
  if(prof=="einasto"){
    uk=gsl_integration_sin_kernel(idensity_fourier,(void *)&this->s_cosmo_pars, static_cast<gsl_real>(k), 1e-5, static_cast<gsl_real>(this->rvir));
    uk=4.0*M_PI*uk/k;
  }
  return uk;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gsl_real DensityProfiles::idensity_fourier(gsl_real r, void *ip){
  /*Intergand used to convert rho(r) -> u(k)*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)ip;
  DensityProfiles Sd(*s_cp);
  Cosmology Cf(*s_cp);
  real_prec M=s_cp->aux_var1;
  real_prec z=s_cp->aux_var2;
  return static_cast<gsl_real>((Sd.density_r(r))*r);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::density_kc(real_prec k, real_prec M, real_prec z){
  /*Aca definimos rvir independientemente, solo como funcion de la masa y no de la concentracion como sale de las void_%%_parameres*/
  Cosmology Cf(this->s_cosmo_pars);
  real_prec dcontrast=Cf.density_contrast_top_hat(z);
  real_prec rvir=pow((real_prec)(3.*M/(4.*M_PI*(dcontrast)*Cf.mean_matter_density(z))),(real_prec)1./3.);
  real_prec uk=gsl_integration_sin_kernel(idensity_fourier,(void *)&this->s_cosmo_pars, k, 1e-5, rvir);
  return 4.0*M_PI*uk/k / M;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
gsl_real DensityProfiles::idensity_fourier_c(gsl_real r,void *ip){
  //   /*INtergand used to convert rho_c(r) -> u_c(k), where _c means profiles convolved with the c-m distribution function*/
  DensityProfiles Sd;
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)ip;
  real_prec M=s_cp->aux_var1;
  real_prec z=s_cp->aux_var2;
  return (Sd.density_rc(r,M,z,ip))*r;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DensityProfiles::einasto_parameters(real_prec M, real_prec z, real_prec *alpha, real_prec *r_vir, real_prec *c,real_prec *rhos,real_prec *rs){
  real_prec deltac=this->s_cosmo_pars.critical_density;
  real_prec dcontrast=this->s_cosmo_pars.density_contrast_top_hat;
  real_prec mean_matter_density=this->s_cosmo_pars.mean_matter_density;
  vector<gsl_real>v_sigma_mass=this->s_cosmo_pars.v_sigma_mass;
  vector<gsl_real>v_mass=this->s_cosmo_pars.v_mass;
  real_prec nu=pow(deltac/gsl_inter_new(v_mass,v_sigma_mass, M),2);
  *r_vir=pow((real_prec)(3.*M/(4.*M_PI*(dcontrast)*mean_matter_density)),(real_prec)1./3.);
  *alpha=0.155+0.0095*nu;
  *c    = 5.70*pow(M/(2e12),(real_prec)-0.092)*pow(1+z,(real_prec)-0.55); /*c_200*/
  *rs   = *r_vir/(*c);
  *rhos = M*pow((*alpha)/2,(-3+(*alpha))/(*alpha))*exp(-2/(*alpha))*pow((*rs),-3)/(2.*M_PI*gsl_sf_gamma_inc_P((real_prec)3./(*alpha),(real_prec)(2./(*alpha))*pow((*c),(*alpha)))*gsl_sf_gamma_inc_P((real_prec)3./(*alpha),(real_prec)(2./(*alpha))*((*r_vir)/(*rs)))*gsl_sf_gamma(3./(*alpha)));
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DensityProfiles::einasto_parameters(real_prec M){
  real_prec deltac=this->s_cosmo_pars.critical_density;
  real_prec dcontrast=this->s_cosmo_pars.density_contrast_top_hat;
  real_prec mean_matter_density=this->s_cosmo_pars.mean_matter_density;
  vector<gsl_real>v_sigma_mass=this->s_cosmo_pars.v_sigma_mass;
  vector<gsl_real>v_mass=this->s_cosmo_pars.v_mass;
  real_prec nu=pow(deltac/gsl_inter_new(v_mass,v_sigma_mass, M),2);
  this->rvir=pow((real_prec)(3.*M/(4.*M_PI*(dcontrast)*mean_matter_density)),(real_prec)1./3.);
  this->alpha=0.155+0.0095*nu;
  this->concentration = 5.70*pow(M/(2e12),(real_prec)-0.092)*pow(1+this->s_cosmo_pars.cosmological_redshift ,(real_prec)-0.55); /*c_200*/
  this->rs   = this->rvir/(this->concentration);
  this->rhos = M*pow((alpha)/2,(-3+(alpha))/(alpha))*exp(-2/this->alpha)*pow(this->rs,-3)/(2.*M_PI*gsl_sf_gamma_inc_P(3./this->alpha,(2./this->alpha)*pow(this->concentration,this->alpha))*gsl_sf_gamma_inc_P(3./this->alpha,2./this->alpha)*(this->rvir/this->rs)*gsl_sf_gamma(3./alpha));
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DensityProfiles::nfw_parameters(real_prec M, real_prec z, real_prec *fc, real_prec *r_vir,real_prec *c,real_prec *rhos,real_prec *rs){
  Cosmology Cf(this->s_cosmo_pars);
  *c    = (this->s_cosmo_pars.coef_concentration_amp/(1.+z))*pow(M/this->s_cosmo_pars.Mnl,this->s_cosmo_pars.coef_concentration);
  *fc   = log(1+(*c))-((*c)/(1.+(*c)));
  real_prec density_contrast_top_hat=Cf.density_contrast_top_hat(z);
  real_prec mean_mat_den=Cf.mean_matter_density(z);
  *rhos = (1./3.)*(density_contrast_top_hat)*(mean_mat_den)*pow(*c,3)/(*fc);
  real_prec rho_aux=(1./3.)*(density_contrast_top_hat)*(mean_mat_den)*pow(*c,3);
  *rs   = pow(M/(4.*M_PI*rho_aux),(real_prec)1./3.);
  *r_vir= (*c)*(*rs);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DensityProfiles::nfw_parameters(real_prec M){
  this->Cosmo.set_cosmo_pars(this->s_cosmo_pars);
  this->concentration    = (this->s_cosmo_pars.coef_concentration_amp/(1.+this->s_cosmo_pars.cosmological_redshift))*pow(M/this->s_cosmo_pars.Mnl,this->s_cosmo_pars.coef_concentration);
  this->fc   = log(1+this->concentration)-this->concentration/(1.+this->concentration);
  real_prec density_contrast_top_hat=this->Cosmo.density_contrast_top_hat(this->s_cosmo_pars.cosmological_redshift);
  real_prec mean_mat_den=this->Cosmo.mean_matter_density(this->s_cosmo_pars.cosmological_redshift);
  this->rhos = (1./3.)*(density_contrast_top_hat)*(mean_mat_den)*pow(this->concentration,3)/(this->fc);
  real_prec rho_aux=(1./3.)*(density_contrast_top_hat)*(mean_mat_den)*pow(this->concentration,3);
  this->rs   = pow(M/(4.*M_PI*rho_aux),1./3.);
  this->rvir= concentration*this->rs;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec DensityProfiles::mass_concentration_dis(real_prec lc,real_prec M,real_prec z,  void *p){
  /*Log-normal disperssion in the assignement of concentration for a given Mass
    lc= conentacion
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec al,fc,rv,mean_c,rhos,rs;
  string prof=s_cp->density_profile;
  real_prec sigma_cln=0.8;
  if(prof=="nfw")
    this->nfw_parameters(M,z, &fc,&rv,&mean_c,&rhos,&rs);
  if(prof=="einasto")
    this->einasto_parameters(M,z,&al,&rv,&mean_c,&rhos,&rs);
  return  (1./sqrt(2.*M_PI*sigma_cln*sigma_cln))*exp(-0.5*pow(log(lc/mean_c),2)/(pow(sigma_cln,2)))/lc;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
