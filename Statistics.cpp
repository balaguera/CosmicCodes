//////////////////////////////////////////////////////////
/**
 *  @brief Statistics
 *  @file Statistics.cpp
 *  @brief Methods of the class Statistics
 *  @author Andres Balaguera-Antolínez (ABA)
 *  @details This file contains functions used to transform coordiante system in the catalogues
 */
//////////////////////////////////////////////////////////
#ifndef __STATISTICS__
#define __STATISTICS__
////////////////////////////////////////////////////////////////////////////
#define GSL_INT_SIZE_k 200
#define GSL_INT_SIZE_M 600
////////////////////////////////////////////////////////////////////////////
# include "../headers/Statistics.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Statistics::compute_int_table_mass(real_prec M_min_integration, real_prec M_max_integration, int nss_k){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW_Mass.resize(nss_k);
  this->XX_Mass.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(M_min_integration)),static_cast<gsl_real>(0.99*log10(M_max_integration)),this->wfd,this->XX_Mass,this->WW_Mass);

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Statistics::compute_int_table_mass(real_prec M_min_integration, real_prec M_max_integration){
  wfd =  gsl_integration_glfixed_table_alloc (GSL_INT_SIZE_M);
  this->WW_Mass.resize(GSL_INT_SIZE_M,0);
  this->XX_Mass.resize(GSL_INT_SIZE_M,0);
  real_prec lMIN=static_cast<gsl_real>(log10(M_min_integration));
  real_prec lMAX=static_cast<gsl_real>(0.99*log10(M_max_integration));
  gsl_get_GL_weights(lMIN,lMAX,this->wfd,this->XX_Mass,this->WW_Mass);

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Statistics::compute_int_table_wavenumber(real_prec k_min_integration, real_prec k_max_integration, int nss_k){
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  this->WW_K.resize(nss_k);
  this->XX_K.resize(nss_k);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(k_min_integration)),static_cast<gsl_real>(log10(k_max_integration)),this->wfd,this->XX_K,this->WW_K);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Statistics::compute_int_table_wavenumber(real_prec k_min_integration, real_prec k_max_integration){
  wfd =  gsl_integration_glfixed_table_alloc (GSL_INT_SIZE_k);
  this->WW_K.resize(GSL_INT_SIZE_k,0);
  this->XX_K.resize(GSL_INT_SIZE_k,0);
  gsl_get_GL_weights(static_cast<gsl_real>(log10(k_min_integration)),static_cast<gsl_real>(log10(k_max_integration)),this->wfd,this->XX_K,this->WW_K);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec  Statistics::As2sigma8(s_CosmologicalParameters *scp){
Cosmology cf;
  real_prec lkmin=log10(scp->kmin_int);
  real_prec lkmax=log10(scp->kmax_int);
  real_prec ans=sqrt(gsl_integration(iAs2sigma8,(void *)scp,lkmin,lkmax));
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::iAs2sigma8(gsl_real lk,void *p){  /* Integrand for sigma^2 */
  real_prec k=pow(10,lk);
  s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum ps;
  ps.set_cosmo_pars(*scp);
  Cosmology cf(*scp);
  gsl_real ans=(log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*ps.Linear_Matter_Power_Spectrum(k)*pow(window(k,scp->RR),2);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::sigma_masa(real_prec m, real_prec z){
  /* Sigma as a function of mass */
  real_prec lkmin=log10(this->s_cosmo_pars.kmin_int);
  real_prec lkmax=log10(this->s_cosmo_pars.kmax_int);
  this->s_cosmo_pars.aux_var1 = m;  //log10M
  this->s_cosmo_pars.aux_var2 = z;
//  real_prec ans=gsl_integration(200,isigma_nu,(void *)&this->s_cosmo_pars,lkmin,lkmax);
  real_prec ans=gsl_integration(isigma_nu,(void *)&this->s_cosmo_pars,this->XX_K,this->WW_K,false);
  return sqrt(ans/(2.*M_PI*M_PI));
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::peak_height(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  /* Sigma as a function of mass */
  real_prec lkmin=log10(scp->kmin_int);
  real_prec lkmax=log10(scp->kmax_int);
  scp->aux_var1 = m;  //log10M
  scp->aux_var2 = z;
  real_prec ans=gsl_integration(400,isigma_nu,(void *)scp,lkmin,lkmax);
  ans=sqrt(ans/(2.*M_PI*M_PI));
  real_prec delta_c=this->cosmo.critical_overdensity(z);
  return delta_c/ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::sigma_masa(real_prec m, real_prec z, s_CosmologicalParameters *scp, vector<gsl_real>XX, vector<gsl_real>WW){
  /* Sigma as a function of mass */
  scp->aux_var1 = m;
  scp->aux_var2 = z;
  real_prec ans=gsl_integration(isigma, (void *)scp, XX, WW,false);
  return sqrt(ans);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::isigma(gsl_real lk,void *p){  /* Integrand for sigma^2 */
  real_prec k=pow(10,lk);
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum ps;
  ps.set_cosmo_pars(*s_cp);
  Cosmology cf(*s_cp);
  ps.use_external_power=s_cp->use_external_power;
  real_prec m= s_cp->aux_var1;
  real_prec z= s_cp->aux_var2;
  real_prec M=pow(10,m);
  real_prec power=0;
  real_prec rr=cf.rr(M,z);
 // if(true==s_cp->use_external_power)
 //     power= pow(s_cp->growth_factor,2)*gsl_inter_new(s_cp->kvector_external, s_cp->power_external, k);
 //  else
//      power=ps.Linear_Matter_Power_Spectrum(k);
  power=ps.Linear_Matter_Power_Spectrum(k);
  real_prec ans=(log(10.0)*k)*pow(k,2)*power*pow(window(k,rr),2);
  return static_cast<gsl_real>(ans);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::isigma_nu(gsl_real lk,void *p){  /* Integrand for sigma^2 */
  real_prec k=pow(10,lk);
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum ps;
  ps.set_cosmo_pars(*s_cp);
  Cosmology cf(*s_cp);
  real_prec m= s_cp->aux_var1;
  real_prec z= s_cp->aux_var2;
  real_prec M=pow(10,m);
  real_prec rr=cf.rr_lag(M,z);
  real_prec ans=(log(10.0)*k)*pow(k,2)*(ps.Linear_Matter_Power_Spectrum(k))*pow(window(k,rr),2);
  return static_cast<gsl_real>(ans);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// real_prec Statistics::sigma_l(real_prec l, real_prec z, void *p){
//     /* Sigma as a function of X-ray luminosity for the REFLEX II sample, 
//        using R=Vmax**1/3, with Vmax as a function of Lx
//     */
//     struct two_pars tpar=  {l,z};
//     return sqrt(gsl_integration(isigma_l,&tpar,k_min,k_max));
//   }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mass_function_D(real_prec m, real_prec z)
{
  MASS_BIAS_FUNCTIONS mb;
//  real_prec lkmin=log10(this->s_cosmo_pars.kmin_int);
//  real_prec lkmax=log10(this->s_cosmo_pars.kmax_int);
  real_prec deltac=this->s_cosmo_pars.critical_density;
//  real_prec Smasa=gsl_inter_new(this->s_cosmo_pars.v_mass, this->s_cosmo_pars.v_sigma_mass,static_cast<gsl_real>(m));    //
  real_prec M=pow(10,m);
  real_prec Smasa=sigma_masa(m,z);
  real_prec nu=pow(deltac/Smasa,2);
  real_prec com_density=(this->s_cosmo_pars.mean_matter_density)*pow(1+z,(real_prec)-3.0);
  this->s_cosmo_pars.aux_var1 = m;
  this->s_cosmo_pars.aux_var2 = z;
  real_prec dsdr=gsl_integration(dsigma_dR,(void *)&this->s_cosmo_pars,this->XX_K,this->WW_K, false);
  return (com_density/static_cast<double>(M))*(mb.mass_function(nu,z, (void *)&this->s_cosmo_pars))*pow(Smasa,-2)*(this->cosmo.rr(M,z)/(3.*static_cast<double>(M)))*dsdr;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mass_function_D(real_prec m, real_prec z,vector<gsl_real>&XX, vector<gsl_real>&WW)
{
  MASS_BIAS_FUNCTIONS mb;
  real_prec deltac=this->s_cosmo_pars.critical_density;
  real_prec Smasa=gsl_inter_new(this->s_cosmo_pars.v_mass, this->s_cosmo_pars.v_sigma_mass,static_cast<gsl_real>(m));    //     sigma_masa(M,z,scp);
  real_prec M=pow(10,m);
  real_prec nu=pow(deltac/Smasa,2);
  real_prec com_density=(this->s_cosmo_pars.mean_matter_density)*pow(1+z,(real_prec)-3.0);
  this->s_cosmo_pars.aux_var1 = m;
  this->s_cosmo_pars.aux_var2 = z;
  real_prec dsdr=gsl_integration(dsigma_dR,(void *)&this->s_cosmo_pars,XX,WW,false);
  return (com_density/M)*(mb.mass_function(nu,z, (void *)&this->s_cosmo_pars))*pow(Smasa,-2)*(this->cosmo.rr(M,z)/(3.*M))*dsdr;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::dsigma_dR(gsl_real lk, void *p){  /* integrand to calculate the derivative of sigma^2 with respect to  R(M)*/
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  PowerSpectrum ps;
  ps.set_cosmo_pars(*s_cp);
  Cosmology cf;
  cf.set_cosmo_pars(*s_cp);
  real_prec m= s_cp->aux_var1;
  real_prec M=pow(10,m);
  real_prec z= s_cp->aux_var2;
  real_prec k=pow(10,lk);
  real_prec r =cf.rr(M,z);
  real_prec y=r*k;
  real_prec power=1.0;
  if(true==s_cp->use_external_power)
      power= pow(s_cp->growth_factor,2)*gsl_inter_new(s_cp->kvector_external, s_cp->power_external, k);
   else
      power=ps.Linear_Matter_Power_Spectrum(k);
  return  (log(10.0)*k)*(1./(2.*pow(M_PI,2)))*pow(k,2)*power*2.0*fabs(window(k,r))*k*fabs(((3.*pow(y,-2))*sin(y)-9.*pow(y,-4)*(sin(y)-y*cos(y))));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  MASS_BIAS_FUNCTIONS mb;
  real_prec Smasa=gsl_inter_new(scp->v_mass, scp->v_sigma_mass,m);    //     sigma_masa(M,z,scp);
  real_prec nu=pow((scp->critical_density)/Smasa,2.);
  return mb.dm_h_bias(z,nu,(void *)scp);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::effective_halo_mass_bias(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  // This computes the halo_mass bias of haloes with masses greater than M
  // note that since we are interpolating over quantities already computed, 
  // we do not need now explicitely the redshift
   scp->aux_var2=z;
  return gsl_integration(i_effective_halo_mass_bias,(void *)scp, m, log10(scp->M_max_effective))/gsl_integration(i_mass_function,(void *)scp, m, log10(scp->M_max_effective));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::i_effective_halo_mass_bias(gsl_real m, void *p){
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec M=pow(10,m);
  real_prec jacobian= log(10.0)*M;
  vector<gsl_real> v_mass = s_cp->v_mass;
  vector<gsl_real> v_mass_function = s_cp->v_mass_function;
  vector<gsl_real> v_halo_mass_bias = s_cp->v_halo_mass_bias;
  return jacobian*gsl_inter_new(v_mass,v_mass_function,m)*gsl_inter_new(v_mass,v_halo_mass_bias,m);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::effective_halo_mean_number_density(real_prec m, real_prec z, s_CosmologicalParameters *scp){
  // This computes the mean number density of objects with masses greater than M
scp->aux_var2=z;     
   // note that since we are interpolating over quantities already computed, 
  // we do not need now explicitely the redshift
  return gsl_integration(i_mass_function,(void *)scp, m, log10(scp->M_max_effective));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real Statistics::i_mass_function(gsl_real m, void *p){
  // I integrate wrt log10(M); foir this reason, this function expects m=lg10(M)
  // insted of M. If I integrate with respect to M, I should pass M to this function.
  // In that case, jacobian = 1
  real_prec M=pow(10,m);
  real_prec jacobian= log(10.0)*M;
  s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  vector<gsl_real> v_mass = s_cp->v_mass;
  vector<gsl_real> v_mass_function = s_cp->v_mass_function;
  return jacobian*gsl_inter_new(v_mass,v_mass_function,m);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// real_prec Statistics::cluster_mass_function(real_prec x, real_prec xt, real_prec z, s_CosmologicalParameters *scp){
//   return gsl_integration(i_cluster_mass_function_feedback,p,m_min,m_max);
// }

// real_prec Statistics::i_cluster_mass_function_feedback(real_prec m, void *p){

//   real_prec xt=(params->p1);
//   real_prec z=(params->p2);
//   Astrophysics ap;
//   return (log(10.0)*pow(10,m)*(M_reference))*(ap.Mobs_Mnb_distribution(m,0, p))*gsl_inter_pointers(mass_array_p,mass_functionpoints,nn,m);
// }
// real_prec Statistics::scale_dependent_bias(real_prec x, real_prec z, void *p){
//   real_prec r=x;
//   /*r en log(scale)*/
//   /*Tinker parametrization of scale-dependent bias in terms of the non-linear correlatin function*/
//   real_prec xi= gsl_inter_pointers(xRp,XI_NLp,nn,r);
//   //  real_prec f= ( r<*xRp ? 0 : pow(1.+1.17*xi,1.49)/pow(1.+0.69*xi,2.09));
//   real_prec f=  pow(1.+1.17*xi,1.49)/pow(1.+0.69*xi,2.09);
//   return f;
// }
// // ******************************************************************************
// // ******************************************************************************

// real_prec Statistics::delta_fof(real_prec x, real_prec z, void *p){
//   /*Calculo de Delta para halos FOF siguiendo los resultados de MOre et al 2011, en donde 
//     se asume un perfil NFW y una concentracion dada por la masa*/
//   real_prec al,fc,rv,c,rhos,rs,ans;
//   /*NOT ACCURATE. WARNING*/
//   //nfw_parameters(x,&fc,&rv,&mean_c,&rhos,&rs);
//   c=8.0;
//   fc   = log(1+(c))-((c)/(1.+(c)));
//   ans=3.*(0.652960)*pow(0.2,-3)*fc*(1+c)*pow(c,2);
//   return ans;
// }
// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::mass_function_light_cone(){
//   real_prec ans;  PowerSpectrum ps(*scp);

//   real_prec rcmax,rcmin;
//   real_prec zmax=0.2;
//   real_prec zmin=0.001;
//   Cosmology cf;
//   real_prec mfz[nzpoints+1], zp[nzpoints+1];
//   for(int i=0;i<=nzpoints;i++)zp[i]=0;
//   for(int i=0;i<=nzpoints;i++)mfz[i]=0;  
//   mf_z_generator(x,1.1*zmax,zp,mfz);   /*THE FACTOR 1.1 PREVENTS NANS WHEN ONE NEEDS TO INTERPOLATE...*/
//   mfz_p=&mfz[0];
//   zp_p= &zp[0];
//   rcmax=cf.comoving_distance(zmax);
//   rcmin=cf.comoving_distance(zmin);
//   struct my_parameters par ={0,0,0};
//   return (3./(pow(rcmax,3)-pow(rcmin,3)))*gsl_integration(i_mass_function_m_z,&par,zmin,zmax);
// }
// // ******************************************************************************
// // ******************************************************************************

// real_prec Statistics::mass_functionpointsrediction(real_prec x, real_prec z, void *p){  /*Mass function as a function of x=log10(M/masa_ns)*/
//   /*ACA LAS ENTRADAS SON MASAS; EL Z LO TOMAMOS EL PUNTERO*/
//   real_prec xmax=x;
//   real_prec xmin=z;
//   real_prec dxmax=(M_reference)*pow(10,xmax);
//   real_prec dxmin=(M_reference)*pow(10,xmin);
//   struct my_parameters par ={0,0,0};
//   return gsl_integration(imass_functionpointsrediction,&par,xmin,xmax)/(dxmax-dxmin); // las masas estan en unidades de 10 a la 14  cuando mido n(m) de LBASICC
// }
// // ******************************************************************************
// // ******************************************************************************


// real_prec Statistics::lum_func_prediction(real_prec x, real_prec z, void *p){
//   real_prec l=x;
//   struct my_parameters par=  {l};
//   return gsl_integration(ilum_func_prediction,&par,m_min,m_max);
// }
// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::occupancy_variance_numerator(real_prec x, real_prec z, void *p){   /*See Smith et al 2011*/
//   real_prec l=x;
//   struct my_parameters par=  {l};
//   return gsl_integration(i_occupancy_variance,&par,m_min,m_max);
// }
// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::lum_func_reflex(real_prec x, real_prec z, void *p){
//   real_prec L=exp(x);                      /*L en unidades de 10^44 erg/s/h^2*/
//   return ncero*pow(L/lstar,alpha_lf+1)*Qexponential(qq,-L/lstar)/L; 
// // ******************************************************************************
// // ******************************************************************************
// real_prec Statistics::i_baryon_mass_function_reionization(real_prec x, real_prec z, void *p){
//   // struct my_parameters * params = (struct my_parameters *)p;
//   // real_prec xb=(params->p1);
//   // real_prec z =(params->p2);
//   Astrophysics ap;
//   real_prec ans=(log(10)*pow(10,x)*(M_reference))*ap.baryon_mass_virial_mass_distribution(x,z,p)*gsl_inter_pointers(mass_array_p,mass_functionpoints,nn,x);
//   return  ans;
// }
// real_prec Statistics::baryon_mass_function_reionization(real_prec x, real_prec l, real_prec z, void *p){
//   struct my_parameters par=  {x,z};
//   return gsl_integration(i_baryon_mass_function_reionization,&par,m_min,m_max);
// }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Statistics::non_linear_scales(real_prec *knl,real_prec *Mnl,real_prec *rnl, real_prec *sigman){
  real_prec z=this->s_cosmo_pars.cosmological_redshift;
  real_prec den=this->s_cosmo_pars.mean_matter_density*pow(1+z,-3);
  real_prec sig2,fac2,ms2,df2,dm2;
  /*Looking for mass scale where nu=1 and scales where sigma(M)=1 with Newton-Rhapson algorithm*/
  ms2=-4;
  //  if(ms1<m_min_interp){
  if(ms2<-4){
    cout<<"Potential error in function non_linear_scales:"<<endl;
    cout<<"minimum value of mass smaller than mn_it. CHECK!"<<endl;
  }
  int nr=10;
  real_prec init_mass=this->s_cosmo_pars.v_mass[0];
  for(int i=0;i<nr;i++)
  {
    this->s_cosmo_pars.aux_var1 = log10(pow(10,ms2)*M_reference);
    df2=-gsl_integration(dsigma_dR,(void *)&this->s_cosmo_pars,log10(this->s_cosmo_pars.kmin_int),log10(this->s_cosmo_pars.kmax_int));
    fac2=(4.*M_PI*den*pow(this->cosmo.rr(pow(10,ms2)*M_reference,z),2))/(log(10.0)*(M_reference)*pow(10,ms2));
    if(log10(pow(10,ms2)*M_reference)>init_mass){
      sig2=gsl_inter_new(this->s_cosmo_pars.v_mass, this->s_cosmo_pars.v_sigma_mass, log10(pow(10,ms2)*M_reference));
      dm2=fac2*(pow(sig2,2)-1.0)/df2;
      ms2-=dm2;
    }
  }
  *Mnl=pow(10,ms2)*(M_reference);
  *knl=pow(6.*(pow(10,ms2)*(M_reference))/(M_PI*(this->cosmo.mean_matter_density(z))),-1./3.);
  *sigman=sigma_masa(log10(pow(10,ms2)*M_reference),z);
  *rnl=1./(*knl);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*
real_prec Statistics::mean_galaxy_number_density(real_prec redshift, s_CosmologicalParameters *scp){
  scp->aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration(30,i_mean_galaxy_number_density,(void *)scp,log10(scp->M_min_effective),log10(scp->M_max_effective));
}
*/
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mean_galaxy_number_density(real_prec redshift, vector<gsl_real>&XX, vector<gsl_real>&WW){
    this->s_cosmo_pars.aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration(i_mean_galaxy_number_density,(void *)&this->s_cosmo_pars,XX,WW);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mean_galaxy_number_density(real_prec redshift){
  this->s_cosmo_pars.aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration(i_mean_galaxy_number_density,(void *)&this->s_cosmo_pars,this->XX_Mass,this->WW_Mass,false);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mean_galaxy_number_density(real_prec redshift, real_prec min_mass){
  this->s_cosmo_pars.aux_var3=redshift;
  //Integrate from the value mmin_hod
  return gsl_integration(i_mean_galaxy_number_density,(void *)&this->s_cosmo_pars,this->XX_Mass,this->WW_Mass,log10(min_mass),false);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
gsl_real  Statistics::i_mean_galaxy_number_density(gsl_real m, void *p){
  struct s_CosmologicalParameters * scp= (struct s_CosmologicalParameters *)p;  
  HOD Shod(*scp);
  real_prec M=pow(10,m);
  real_prec jacobian=log(10.)*M;
  return jacobian*gsl_inter_new(scp->v_mass,scp->v_mass_function,m)*(Shod.SATELLITE(M)+Shod.CENTRAL(M));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec Statistics::mean_halo_number_density(s_CosmologicalParameters *scp){
  // Integrates the mass function wrt the mass to get the mean number density of objects
  return gsl_integration(i_mass_function,(void *)scp,log10(scp->M_min_effective),log10(scp->M_max_effective));
}
#endif
