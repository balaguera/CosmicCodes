# include "../Headers/AngularPowerSpectrum.h"



real_prec top_hat(real_prec x){
    real_prec ans=0;
    if(fabs(x)>0.5) ans=0.0;
  else if(fabs(x)==0.5)ans= 0.5;
  else if(fabs(x)<0.5)ans=1.0;
  return ans;
}


// ***************************************************************
void AngularPowerSpectrum::get_bessel(){

  int nbess=1000;

  xBessel.resize(nbess+1);


  //for(int i=0;i<xBessel.size();i++)xBessel[i]=pow(10,log10(0.000001*this->k_min)+i*(log10(this->k_max*1000.0/(0.00001*this->k_min))/nbess));

  for(int i=0;i<xBessel.size();i++)xBessel[i]=pow(10,-5+i*(5.+5.)/nbess);

Bessel.resize(this->L_max_meas+1);
  for(int i=0;i<Bessel.size();i++)Bessel[i].resize(nbess);
  
  for(int xi=0;xi<Bessel[0].size();xi++)Bessel[0][xi]=gsl_sf_bessel_jl(0,xBessel[xi]);
  for(int xi=0;xi<Bessel[0].size();xi++)Bessel[1][xi]=gsl_sf_bessel_jl(1,xBessel[xi]);

  cout<<"Computing Bessel "<<endl;  
  for(int li=1;li<30;li++)for(int i=0;i<xBessel.size();i++)Bessel[li][i]=gsl_sf_bessel_jl(li,xBessel[i]);
  for(int li=40;li<this->L_max_meas;li++)for(int i=0;i<xBessel.size();i++)Bessel[li+1][i]=(2*li+1)*Bessel[li][i]/xBessel[i] - Bessel[li-1][i];

  cout<<"Done "<<endl;  
}



// ***************************************************************


void AngularPowerSpectrum::get_cosmo(Cosmology cf, s_CosmologicalParameters scp){


  rv.resize(this->nz);
  zv.resize(this->nz);
  gv.resize(this->nz);
  Hv.resize(this->nz);
 
  ofstream da;
  da.open("cosmo_pars.dat");
  real_prec zmn=0.000001;
  for(int i=0;i<this->nz;i++){
    zv[i]=zmn+(1.0-zmn)*i/(real_prec)this->nz;
    rv[i]=cf.comoving_distance(zv[i],(void *) &scp);
    gv[i]=cf.growth_factor(zv[i],(void *)&scp);
    Hv[i]=cf.Hubble_function(zv[i],(void *)&scp);
    da<<zv[i]<<"  "<<Hv[i]<<"  "<<rv[i]<<"  "<<gv[i]<<endl;
  }
  da.close();
  //  cf.free_gsl_table();
}
// ***************************************************************

void AngularPowerSpectrum::compute_int_table(){
  int nss_z=500;
  wf =  gsl_integration_glfixed_table_alloc (nss_z);
  WW_z.resize(nss_z);
  XX_z.resize(nss_z);
  gsl_get_GL_weights(z_min,z_max,wf,XX_z,WW_z);

  int nss_k=800;
  wfd =  gsl_integration_glfixed_table_alloc (nss_k);
  WW_k.resize(nss_k);
  XX_k.resize(nss_k);
  gsl_get_GL_weights(log10(k_min),log10(k_max),wfd,XX_k,WW_k);
}

// ***************************************************************

void AngularPowerSpectrum::free_gsl_table(){
  gsl_integration_glfixed_table_free(this->wf);
  gsl_integration_glfixed_table_free(this->wfd);
}

// ***************************************************************

real_prec AngularPowerSpectrum::window(real_prec z, void *p)
{
  struct cl_params * s_p= (struct cl_params *)p;
  real_prec x = (z-s_p->zmean)/(2*s_p->width);
  real_prec ans;
  if(s_p->wtype=="tophat")ans=top_hat(x);
  else if(s_p->wtype=="gaussian")ans= exp(-2.0*(x*x));
  return ans;
}

// ***************************************************************

real_prec AngularPowerSpectrum::dndz(real_prec z, void *p)
{
  struct cl_params * s_p= (struct cl_params *)p;
  return gsl_inter_new(s_p->dz,s_p->dn,z)*window(z,p);
}


// ***************************************************************
gsl_real AngularPowerSpectrum::idndz(gsl_real z, void *p)
{
  AngularPowerSpectrum cl;
  return cl.dndz(z,p);
}
// ***************************************************************

void AngularPowerSpectrum::get_normal(void *p)
{
  normal=gsl_integration2(idndz, p, this->XX_z,this->WW_z);
}

// ***************************************************************

gsl_real AngularPowerSpectrum::iKernel(gsl_real z, void *p)
{
  AngularPowerSpectrum cl;
  struct cl_params * s_p= (struct cl_params *)p;
  gsl_real Bessel;
  real_prec x=s_p->k*gsl_inter_new(s_p->zv,s_p->rv,z);

  if(s_p->l<=5000)  Bessel=gsl_sf_bessel_jl(s_p->l,x);
  else Bessel=gsl_inter_new2(s_p->sxBessel,s_p->sBessel,s_p->l,x);
  
  return gsl_inter_new(s_p->zv,s_p->gv,z)*cl.dndz(z,p)*Bessel;
}

// ***************************************************************
real_prec AngularPowerSpectrum::get_Kernel(void *p){
  return gsl_integration2(iKernel, p, this->XX_z,this->WW_z);
}

// ***************************************************************
// ***************************************************************

gsl_real AngularPowerSpectrum::iCl(gsl_real kl, void *p){

  struct cl_params * sp= (struct cl_params *)p;
  s_CosmologicalParameters scp=sp->scp;
  real_prec k=pow(10,kl);

  //real_prec k=kl;
  struct cl_params * s_p= (struct cl_params *)p;
  PowerSpectrum PS;
  real_prec power=1;//PS.Linear_Matter_Power_Spectrum(scp, k); 
  return (2./M_PI);//*(log(10.0)*k)*pow(k*gsl_inter_new(s_p->kv,s_p->F,k),2.0)*power;
}

// ***************************************************************
// ***************************************************************
real_prec AngularPowerSpectrum::get_cl_exact(void *p){
  return gsl_integration2(iCl, p, this->XX_k, this->WW_k)/pow(normal,2.0);
}

// ***************************************************************
// ***************************************************************

gsl_real AngularPowerSpectrum::iCl_limber(gsl_real z, void *p){
  AngularPowerSpectrum cl;
  PowerSpectrum PS;

  struct cl_params * sp= (struct cl_params *)p;

  s_CosmologicalParameters scp=sp->scp;
  real_prec gg=gsl_inter_new(sp->zv,sp->gv,z);
  real_prec r=gsl_inter_new(sp->zv,sp->rv,z);
  real_prec HH=gsl_inter_new(sp->zv,sp->Hv,z)/speed_light;
  real_prec ks=(sp->l+0.5)/r;
  real_prec power=PS.Q_Model_Matter_Power_Spectrum(&scp, ks); 
  return pow(gg*cl.dndz(z,p)/r,2)*power*HH;
}

// ***************************************************************

gsl_real AngularPowerSpectrum::get_cl_ilimber(void *p){
  // return gsl_integration(iCl_limber, p, Z_min,Z_max)/pow(normal,2);
  return gsl_integration2(iCl_limber, p, this->XX_z,this->WW_z);
}

// ***************************************************************************************
// ***************************************************************************************

void AngularPowerSpectrum::set_mixingM(string file){

  this->R.resize(this->nlbins+1);

  for(int i=0;i<R.size();i++)R[i].resize(this->L_max_meas+1);

  vector< vector<real_prec> > pR;
  Fmd.read_file(file,pR);
  
  int ik=-1; 
  for(int i=0;i<this->nlbins;i++){
    for(int j=0;j<this->L_max_meas+1;j++){
      ++ik;
      this->R[i][j]=pR[ik][2];
    }
  }
  pR.clear();


}


// ***************************************************************************************
// ***************************************************************************************

// void AngularPowerSpectrum::Pk_normalization(s_CosmologicalParameters scp){
//   PS.normalization((void *)&scp, this->pk_normalization);
// }



// ***************************************************************************************

void  AngularPowerSpectrum::get_cl_meas(string cl_meas_file){
  vector< vector<real_prec> > cl_m;

  Fmd.read_file(cl_meas_file,cl_m);

  cl_meas.resize(cl_m.size(),0);
  sigma_cl.resize(cl_m.size(),0);

  cl_model.resize(this->nlbins,0);

  for(int i=0;i<cl_m.size();i++)this->cl_meas[i]=cl_m[i][3];
  for(int i=0;i<cl_m.size();i++)this->sigma_cl[i]=cl_m[i][4];
  cl_m.clear();
}


void  AngularPowerSpectrum::get_dndz(string dndz_file){
  vector< vector<real_prec> > dnz;
  Fmd.read_file(dndz_file,dnz);
  dn.resize(dnz.size(),0);
  dz.resize(dnz.size(),0);
  for(int i=0;i<dnz.size();i++)this->dz[i]=dnz[i][0];
  for(int i=0;i<dnz.size();i++)this->dn[i]=dnz[i][1];
  dnz.clear();
  this->z_min=this->dz[0];
  this->z_max=this->dz[dz.size()-1];
}


// ***************************************************************************************
void AngularPowerSpectrum::get_cl_limber(cl_params sp){

  vector<real_prec> clr(this->L_max_meas+1,0);
  
  //  omp_set_num_threads(omp_get_max_threads());
  //#pragma omp parallel for
  for(int l=this->L_min;l<this->L_max_meas+1;l++)
    {
      sp.l=l;
      clr[l]=static_cast<real_prec>(get_cl_ilimber((void *)&sp)/pow(this->normal,2));
    }
  
  
  
  for(int i=0;i<this->nlbins;++i)this->cl_model[i]=0;
  
  //#pragma omp parallel for
  for(int i=0;i<this->nlbins;++i){
    real_prec aux=0;

    for(int lp=this->L_min;lp<this->L_max_meas+1;lp++)aux+=this->R[i][lp];

    for(int lp=this->L_min;lp<this->L_max_meas+1;lp++)
      this->cl_model[i]+=(this->R[i][lp]*clr[lp]/aux);
    
  }


  clr.clear();


}
