//////////////////////////////////////////////////////////
/**
 * @file CosmologicalFunctions.cpp
 * @brief This file contains methods of the class Cosmology
 * @details The class Cosmology generates values of cosmological functions as a function of redshft and cosmological parameters 
 * @author Andres Balaguera Antolinez 
 * @date 2007-2023
 */
//////////////////////////////////////////////////////////

#include "CosmologicalFunctions.h"

////////////////////////////////////////////////////////////////////////////
real_prec scale_factor(real_prec z){
    return 1./(1.+z);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::comoving_sound_horizon(real_prec redshift)
{
  real_prec a_end=static_cast<gsl_real>(1./(1.+redshift));
  return (Constants::speed_light/sqrt(3.))*static_cast<gsl_real>(gsl_integration(i_rs, (void *)&this->s_cosmo_pars, 0, a_end));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::rr(real_prec M, real_prec z){
  return (1+z)*pow((real_prec)3.*M/(4.*M_PI*this->mean_matter_density(z)),(real_prec)1./3.);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::rr_lag(real_prec M, real_prec z){
  /* Lagrangian radius in Mpc/h
    here z is the mass in the units of this code, z=x*/
  return pow(3.*M/(4.*M_PI*this->mean_matter_density(z)),1./3.);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::SO(real_prec z, real_prec M, real_prec r){
  return 3.*M/(4.*M_PI*pow(r, 3)*this->mean_matter_density(z));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::critical_overdensity(real_prec z){
  return  (3./20.)*pow(12.*M_PI,2./3.);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::density_contrast_top_hat(real_prec z){
    real_prec x=this->omega_matter(z)-1.;
    return  (18.0*M_PI*M_PI+82.0*x-39.0*pow(x,2))/this->omega_matter(z);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::Hubble_function(real_prec redshift){
  return this->s_cosmo_pars.Hubble*sqrt(this->s_cosmo_pars.Om_matter*pow(1.+redshift,3.)+this->s_cosmo_pars.Om_radiation*pow(1.+redshift,4.)+this->s_cosmo_pars.Om_k*pow(1+redshift,2.)+this->s_cosmo_pars.Om_vac*pow(1+redshift,3.*(1.+this->s_cosmo_pars.wde_eos)));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::comoving_distance(real_prec redshift){
  return gsl_integration(Hinv, (void *)&this->s_cosmo_pars, 0.0, redshift)*Constants::speed_light;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::transverse_comoving_distance(real_prec redshift)
{
  real_prec fac=this->s_cosmo_pars.Hubble*sqrt(fabs(this->s_cosmo_pars.Om_k));
  real_prec ans;

  if(this->s_cosmo_pars.Om_k<0){
    ans= Constants::speed_light*sin(fac*this->comoving_distance(redshift)/Constants::speed_light)/fac;
  }
  else if(this->s_cosmo_pars.Om_k>0){
    ans= Constants::speed_light*sinh(fac*this->comoving_distance(redshift)/Constants::speed_light)/fac;
  }
  else{
    ans= this->comoving_distance(redshift);
  }
  return ans;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::inter_transverse_comoving_distance(real_prec redshift)
{
  real_prec fac=this->s_cosmo_pars.Hubble*sqrt(fabs(this->s_cosmo_pars.Om_k));
  real_prec ans;
  real_prec cd=gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.rv,redshift);
  if(this->s_cosmo_pars.Om_k<0){
    ans= Constants::speed_light*sin(fac*cd/Constants::speed_light)/fac;
  }
  if(this->s_cosmo_pars.Om_k==0){
    ans=cd;
  }
  if(this->s_cosmo_pars.Om_k>0){
    ans= Constants::speed_light*sinh(fac*cd/Constants::speed_light)/fac;
  }
  return ans;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::derivative_transverse_comoving_distance(real_prec redshift)
{
  real_prec fac=this->s_cosmo_pars.Hubble*sqrt(fabs(this->s_cosmo_pars.Om_k));
  real_prec cd=gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.rv, redshift);
  real_prec ans;
  if(this->s_cosmo_pars.Om_k<0){
    ans= (Constants::speed_light/this->Hubble_function(redshift))*cos(fac*cd/Constants::speed_light)/fac;
  }
  if(this->s_cosmo_pars.Om_k==0){
    ans= Constants::speed_light/this->Hubble_function(redshift);
  }
  if(this->s_cosmo_pars.Om_k>0){
    ans= (Constants::speed_light/this->Hubble_function(redshift))*cosh(fac*cd/Constants::speed_light)/fac;
  }
  return ans;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::comoving_angular_diameter_distance(real_prec redshift)
//This is proper
{
  return transverse_comoving_distance(redshift);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::proper_angular_diameter_distance(real_prec redshift)
//This is proper
{
  return transverse_comoving_distance(redshift)/(1.+redshift);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::inter_proper_angular_diameter_distance(real_prec redshift)
//This is proper
{
  return gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.trv, redshift)/(1.+redshift);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::luminosity_distance(real_prec redshift)
{
  return transverse_comoving_distance(redshift)*(1.+redshift);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::inter_luminosity_distance(real_prec redshift)
{
  return gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.trv, redshift)*(1.+redshift);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::mean_matter_density(real_prec redshift)
/*M */
{
  return (3.*this->s_cosmo_pars.Hubble*this->s_cosmo_pars.Hubble/(8.*M_PI*Gravitational_constant))*omega_matter(redshift)*(Mpc_to_km/Solar_mass);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::age_universe(real_prec redshift)
/*Age of the universe in years / h*/
{
  return gsl_integration(Hinva,(void *)&this->s_cosmo_pars,static_cast<gsl_real>(redshift),static_cast<gsl_real>(1e4))*(Constants::Mpc_to_km/Constants::years_to_sec);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::growth_factor(real_prec redshift){
  real_prec a=1./(1.+redshift);
  return (5./2.)*(this->s_cosmo_pars.Om_matter)*this->Hubble_function(redshift)*gsl_integration(gint,(void *)&this->s_cosmo_pars,-10,log10(a))*pow(this->s_cosmo_pars.Hubble,2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::growth_index(real_prec redshift){
  real_prec om=this->omega_matter(redshift);
  return pow(om, 5./9.);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::growth_index2(real_prec redshift){
  real_prec om=this->omega_matter(redshift);
  return 2.*pow(om, 6./11.);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::halo_dynamical_time(real_prec redshift){
  return (0.1/this->Hubble_function(redshift))*(Mpc_to_km/years_to_sec);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::omega_matter(real_prec redshift){
    // Om(z) = Om(z=0) * (1+z)³ / E(z)²
  return this->s_cosmo_pars.Om_matter*pow(1+redshift,3)*pow(Hubble_function(redshift)/this->s_cosmo_pars.Hubble,-2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::omega_radiation(real_prec redshift){
  return this->s_cosmo_pars.Om_radiation*pow(1+redshift,4)*pow(Hubble_function(redshift)/this->s_cosmo_pars.Hubble,-2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::omega_curvature(real_prec redshift){
  return this->s_cosmo_pars.Om_k*pow(1+redshift,2)*pow(Hubble_function(redshift)/this->s_cosmo_pars.Hubble,-2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::omega_dark_energy(real_prec redshift){
  return this->s_cosmo_pars.Om_vac*pow(1+redshift,3.*(1+this->s_cosmo_pars.wde_eos))*pow(Hubble_function(redshift)/this->s_cosmo_pars.Hubble,-2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: Distance_Modulus(real_prec redshift){
  // ***************************
  // DISTANCE MODULUS
  // ***************************
  real_prec dm=(redshift==0 ? 0 : 25.+5.0*log10(luminosity_distance(redshift)));
  return dm;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: inter_Distance_Modulus(real_prec redshift){
  // ***************************
  // DISTANCE MODULUS
  // ***************************
  return 25.+5.0*log10(inter_luminosity_distance(redshift));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: K_correction(real_prec redshift){
  // ***************************
  // K-CORRECTION
    real_prec a0=0.01529;
    real_prec a1=3.2838;
    real_prec a2=-13.3196;
    real_prec a3=120.51;
    real_prec a4=-391.612;
    real_prec a5=395.23;
    return a0+a1*redshift+a2*pow(redshift,2)+a3*pow(redshift,3)+a4*pow(redshift,4)+a5*pow(redshift,5);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: K_correction(real_prec redshift, real_prec color1, real_prec color2){
  // ***************************
  // K-CORRECTION
  // ***************************
  //  return -6.0*log10(1+redshift);
  real_prec a0=this->s_cosmo_pars.K_index_a;
  real_prec a1=this->s_cosmo_pars.K_index_b;
  real_prec a2=this->s_cosmo_pars.K_index_c;
  real_prec a3=this->s_cosmo_pars.K_index_d;
  real_prec a4=this->s_cosmo_pars.K_index_e;
  real_prec a5=this->s_cosmo_pars.K_index_f;
  real_prec a6=this->s_cosmo_pars.K_index_g;
  real_prec a7=this->s_cosmo_pars.K_index_h;
  real_prec a8=this->s_cosmo_pars.K_index_i;
  real_prec a9=this->s_cosmo_pars.K_index_j;
//  real_prec lred=log10(1.0+redshift);
  real_prec lred=redshift;
  real_prec color_poly= a0+a1*color1 + a2*color1*color1 + a3*pow(color1,3)+a4*pow(color1,4);
  real_prec z_poly =    a5*lred + a6*lred*lred + a7*pow(lred,3) + a8*pow(lred,4)+a9*pow(lred,5);
  return color_poly*z_poly;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: dK_correction_dz(real_prec redshift){
    real_prec a1=3.2838;
    real_prec a2=-13.3196;
    real_prec a3=120.51;
    real_prec a4=-391.612;
    real_prec a5=395.23;
    return this->s_cosmo_pars.use_e_correction ? a1+2.*a2*redshift+3.*a3*pow(redshift,2)+4.*a4*pow(redshift,3)+5.*a5*pow(redshift,4): 0. ;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: e_correction(real_prec redshift){
  return this->s_cosmo_pars.use_e_correction ?  this->s_cosmo_pars.e_index_a*pow(redshift,2)+this->s_cosmo_pars.e_index_b*redshift+this->s_cosmo_pars.e_index_c  :   0.0 ;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology:: de_correction_dz(real_prec redshift){
  return this->s_cosmo_pars.use_e_correction ? 2.*this->s_cosmo_pars.e_index_a*redshift+this->s_cosmo_pars.e_index_b:  0. ;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::zmax_old(){
  real_prec z_ini=0.01;
  real_prec z=z_ini;
  vector<real_prec> zmm;
  int i=0;
  do{
    i++;
    if(z<0)std::cerr<<"Negative reds    hift present"<<endl;
    real_prec tcd = gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.trv, z);   //Trnasverse D<(z)
    real_prec F  = this->s_cosmo_pars.Mabs-this->s_cosmo_pars.mlim+25.0+5.0*log10(tcd*(1+z))+this->e_correction(z)+this->K_correction(z) ; //Function F(z)
    real_prec dF=  (5.0/log(10.0))*( 1./(1.+z) + derivative_transverse_comoving_distance(z)/tcd)+ this->dK_correction_dz(z)+ this->de_correction_dz(z);  // Derivative dF/dz
    z -= F/dF;
    zmm.push_back(z);

  }while(fabs((zmm[zmm.size()-1]-zmm[zmm.size()-2])/zmm[zmm.size()-2])>ERROR_LIMIT_NR);
  return z;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::zmax_old(real_prec z_ini){
  real_prec Mabs=this->s_cosmo_pars.Mabs;
  real_prec mlim=this->s_cosmo_pars.mlim;
  real_prec z=z_ini;
  vector<real_prec> zmm;
  int i=0;
  do{
    i++;
    real_prec tcd = gsl_inter_new(this->s_cosmo_pars.zv, this->s_cosmo_pars.trv, z);
    real_prec F  = Mabs-mlim+25.0+5.0*log10(tcd*(1+z))+K_correction(z)+e_correction(z);
    real_prec dF=  (5.0/log(10.0))*(1./(1.+z)+derivative_transverse_comoving_distance(z)/tcd)+ this->dK_correction_dz(z)+this->de_correction_dz(z);
    z -= F/dF;
    zmm.push_back(z);
  }while(fabs(zmm[zmm.size()-1]-zmm[zmm.size()-2])>ERROR_LIMIT_NR);
return z;
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::zmax(){
  int status;
  int iter = 0, max_iter = 50;
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  real_prec x0, x = 0.005;
  gsl_function_fdf FDF; 
  FDF.f = &Froot;
  FDF.df = &dFroot;
  FDF.fdf = &F_dF;
  FDF.params = (void*)&this->s_cosmo_pars;
  T = gsl_root_fdfsolver_secant;
  s = gsl_root_fdfsolver_alloc (T);
  gsl_root_fdfsolver_set (s, &FDF, x);
  do{
    iter++;
    status = gsl_root_fdfsolver_iterate (s);
    x0 = x;
    x = gsl_root_fdfsolver_root (s);
    status = gsl_root_test_delta (x, x0, 0, 1e-3);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fdfsolver_free (s);
  return x;
}
////////////////////////////////////////////////////////////////////////////
void Cosmology::free_gsl_table(){
  gsl_integration_glfixed_table_free(this->wf);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::gsl_cosmo_integration(gsl_real (*function)(gsl_real, void *) ,void *p,gsl_real LowLimit,gsl_real UpLimit){
  gsl_function F;
  F.params   = p;  
  F.function = function;
  return gsl_integration_glfixed(&F,LowLimit,UpLimit,this->wf);
}    
////////////////////////////////////////////////////////////////////////////
void Cosmology::check_cosmo_pars(){
/*
   this->s_cosmo_pars.Om_matter = this->s_cosmo_pars.Om_cdm+this->s_cosmo_pars.Om_baryons;
   this->s_cosmo_pars.f_baryon  = this->s_cosmo_pars.Om_baryons/this->s_cosmo_pars.Om_matter;
   this->s_cosmo_pars.Om_vac    = 1.- this->s_cosmo_pars.Om_matter-this->s_cosmo_pars.Om_radiation-this->s_cosmo_pars.Om_k;
   */
}
////////////////////////////////////////////////////////////////////////////
void Cosmology::Comoving_distance_tabulated(real_prec z_min, real_prec z_max, vector<gsl_real> & zz, vector<gsl_real> & rc)
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0; i<zz.size(); i++)
    {
       zz[i]=static_cast<gsl_real>(z_min+i*(z_max-z_min)/(zz.size()-1.0));
       rc[i]=static_cast<gsl_real>(this->comoving_distance(zz[i]));
  }  
}
////////////////////////////////////////////////////////////////////////////
// Write the static functions
gsl_real gint(gsl_real lscale_factor, void *p){   /*Used to obtain D(redshift)*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  real_prec scale_f=pow(10,lscale_factor);
  real_prec red_shift=-1.0+1./static_cast<double>(scale_f);
  return static_cast<gsl_real>((scale_f*log(10.))*pow(scale_f*cf.Hubble_function(red_shift),-3));
}
////////////////////////////////////////////////////////////////////////////
gsl_real i_rs(gsl_real a, void *p){  /*Used to obtain rs(redshift)*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  real_prec credshift = 1./a -1.0;
  real_prec Om_photons= (s_cp->Om_radiation)/(1.+(7./8.)*pow(4./11.,4./3.)*s_cp->N_eff);
  real_prec Om_baryons= (s_cp->Om_baryons);
  real_prec Rg= (3.*Om_baryons)/(4.*Om_photons);
  return static_cast<gsl_real>(1./(a*a*cf.Hubble_function(credshift)*sqrt(1.+Rg*a)));
}
////////////////////////////////////////////////////////////////////////////
gsl_real Hinv(gsl_real redshift, void *p){  /*Used to obtain r(redshift)*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  return static_cast<gsl_real>(1./cf.Hubble_function(redshift));
}
////////////////////////////////////////////////////////////////////////////
gsl_real Hinva(gsl_real redshift, void *p){ /*Used to obtan Age(redshift)*/
   struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
   Cosmology cf(*s_cp);
   return 1./((1.+redshift)*cf.Hubble_function(redshift));
 }
////////////////////////////////////////////////////////////////////////////
gsl_real Froot(gsl_real z, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  gsl_real tcd=gsl_inter_new(s_cp->zv, s_cp->trv,z);
  gsl_real ans=s_cp->Mabs-s_cp->mlim+25.0+5.0*log10(tcd*(1+z))+cf.K_correction(z)+cf.e_correction(z);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
gsl_real dFroot(gsl_real z, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  real_prec tcd=gsl_inter_new(s_cp->zv, s_cp->trv,z);
  real_prec ans=(5.0/log(10.0))*( 1./(1.+z) + cf.derivative_transverse_comoving_distance(z)/tcd)+cf.dK_correction_dz(z)+cf.de_correction_dz(z);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
void F_dF(gsl_real z, void *p, gsl_real *y, gsl_real *dy){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  real_prec tcd=static_cast<real_prec>(gsl_inter_new(s_cp->zv, s_cp->trv, z));
  *y=s_cp->Mabs-s_cp->mlim+25.0+5.0*log10(tcd*(1+z))+cf.K_correction(z)+cf.e_correction(z);
  *dy= (5.0/log(10.0))*( 1./(1.+z) + cf.derivative_transverse_comoving_distance(z)/tcd)+ cf.dK_correction_dz(z)+cf.de_correction_dz(z);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::drag_redshift(){
  real_prec omega_b  = this->s_cosmo_pars.Om_baryons*pow(this->s_cosmo_pars.hubble,2);
  real_prec omega_m  = this->s_cosmo_pars.Om_matter*pow(this->s_cosmo_pars.hubble,2);
  real_prec z_drag_a  = 0.0783*pow(omega_b,-0.238)/(1.+39.5*pow(omega_b,0.763));
  real_prec z_drag_b  = 0.560/(1+21.1*pow(omega_b,1.81));
  return  1048.0*(1.+0.00124*pow(omega_b,-0.738))*(1.+z_drag_a*pow(omega_m,z_drag_b));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::R_drag(){
  real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
  real_prec obhh  = this->s_cosmo_pars.Om_baryons*pow(this->s_cosmo_pars.hubble,2);
  return 31.5*obhh/pow(theta_cmb,4)*(1000.0/(1.+this->drag_redshift()));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::R_equality(){
  real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
  real_prec obhh  = (this->s_cosmo_pars.Om_baryons)*pow(this->s_cosmo_pars.hubble,2);
  return 31.5*obhh/pow(theta_cmb,4)*(1000./this->redshift_equality());
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::redshift_equality(){
  real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
  real_prec omhh  = (this->s_cosmo_pars.Om_matter)*pow(this->s_cosmo_pars.hubble,2);
  return  2.50e4*omhh/pow(theta_cmb,4);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::k_equality(){
  real_prec theta_cmb = this->s_cosmo_pars.Tcmb/2.7;
  real_prec omhh  = (this->s_cosmo_pars.Om_matter)*pow(this->s_cosmo_pars.hubble,2);
  return   0.0746*omhh/pow(theta_cmb,2);
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::k_Silk(){
  real_prec omhh  = (this->s_cosmo_pars.Om_matter)*pow(this->s_cosmo_pars.hubble,2);
  real_prec obhh  = (this->s_cosmo_pars.Om_baryons)*pow(this->s_cosmo_pars.hubble,2);
  return  1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::sound_horizon(){
  real_prec Rdrag=this->R_drag();
  real_prec Req=this->R_equality();
  return 2./3./this->k_equality()*sqrt(6./Req)*log((sqrt(1+Rdrag)+sqrt(Rdrag+Req))/(1+sqrt(Req)));
}
////////////////////////////////////////////////////////////////////////////
real_prec Cosmology::TFsound_horizon_fit(){
  real_prec omhh = this->s_cosmo_pars.Om_matter*pow(this->s_cosmo_pars.hubble,2) ;
  real_prec sound_horizon_fit_mpc = 44.5*log(9.83/omhh)/sqrt(1.+10.0*pow(omhh*this->s_cosmo_pars.f_baryon,0.75));
  return sound_horizon_fit_mpc*this->s_cosmo_pars.hubble;
}
////////////////////////////////////////////////////////////////////////////
  /* Output: The approximate location of the first baryonic peak, in h Mpc^-1 */
  real_prec Cosmology::TFk_peak()
  {
    real_prec omhh = this->s_cosmo_pars.Om_matter*pow(this->s_cosmo_pars.hubble,2) ;
    real_prec k_peak_mpc = 2.5*3.14159*(1+0.217*omhh)/(this->sound_horizon()*this->s_cosmo_pars.hubble);//    or use COsmology::TFsound_horizon_fit(omhh,f_baryon,1.0);
    return k_peak_mpc/this->s_cosmo_pars.hubble;
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
