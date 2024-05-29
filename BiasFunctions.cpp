////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/** @file BiasFunctions.cpp
 *  @brief BiasFunctions 
 *  @author: Andrés Balaguera-Antolínez
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "BiasFunctions.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MASS_BIAS_FUNCTIONS::mass_function(real_prec nu, real_prec z, void *p){
  /*A esto llamo la cantidad nu*f(nu)*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  string mb=s_cp->mass_function_fit;
  gsl_real A_par[10]={0};
  gsl_real a_par[10]={0};
  gsl_real b_par[10]={0};
  gsl_real c_par[10]={0};
  gsl_real D[10]={0};
  D[0]=0;
  D[1]=log10(200);
  D[2]=log10(300);
  D[3]=log10(400);
  D[4]=log10(600);
  D[5]=log10(800);
  D[6]=log10(1200);
  D[7]=log10(1600);
  D[8]=log10(2400);
  D[9]=log10(3200);
  
  A_par[0]=0;
  A_par[1]=0.186;
  A_par[2]=0.200;
  A_par[3]=0.212;
  A_par[4]=0.218;    
  A_par[5]=0.248;
  A_par[6]=0.255;
  A_par[7]=0.260;
  A_par[8]=0.260;
  A_par[9]=0.260;
  
  a_par[0]=0;
  a_par[1]=1.47;
  a_par[2]=1.52;
  a_par[3]=1.56;
  a_par[4]=1.61;    
  a_par[5]=1.87;
  a_par[6]=2.13;
  a_par[7]=2.30;
  a_par[8]=2.53;
  a_par[9]=2.66;
  
  b_par[0]=0;
  b_par[1]=2.57;
  b_par[2]=2.25;
  b_par[3]=2.05;
  b_par[4]=1.87;
  b_par[5]=1.59;
  b_par[6]=1.51;
  b_par[7]=1.46;
  b_par[8]=1.44;
  b_par[9]=1.41;
  
  c_par[0]=0;
  c_par[1]=1.19;
  c_par[2]=1.27;
  c_par[3]=1.34;
  c_par[4]=1.45;
  c_par[5]=1.58;
  c_par[6]=1.80;
  c_par[7]=1.97;
  c_par[8]=2.24;
  c_par[9]=2.44;
  
  gsl_real deltac = cf.critical_overdensity(z);

  real_prec q,pp,AA;
  real_prec a,b,c;
  real_prec ans=0;
  if(mb=="Press_Schechter")ans=sqrt(nu/(2.*M_PI))*exp(-nu/2.0);
  if(mb=="Sheth_Tormen")
  {
    pp=0.3;
    AA=1./(1.+pow(2,-pp)*gsl_sf_gamma(0.5-pp)/sqrt(M_PI));
    q=0.75;
    ans=AA*sqrt(q*nu/(2.*M_PI))*(1.+pow(q*nu,-pp))*exp(-q*nu/2.0);
  }
  if(mb=="Jenkins")
      ans= 0.5*0.315*exp(-pow((real_prec)fabs(log(sqrt(nu)/(deltac))+0.61),(real_prec)3.8));
  if(mb=="Warren")
      ans= 0.5*0.7234*(pow((real_prec)(deltac)/sqrt(nu),(real_prec)-1.625)+0.2538)*exp(-1.1982*nu*pow(deltac,(real_prec)-2));
  if(mb=="Pillepich")
      ans= 0.5*0.6853*(pow((real_prec)(deltac)/sqrt(nu),(real_prec)-1.868)+0.3324)*exp(-1.2266*nu*pow(deltac,(real_prec)-2));
  if(mb=="MICE")
      ans= 0.5*0.5800*(pow((real_prec)(deltac)/sqrt(nu),(real_prec)-1.37)+0.3)*exp(-1.036*nu*pow(deltac,(real_prec)-2));
  if(mb=="Tinker")
   {
    AA=gsl_inter(D,A_par,9,log10((s_cp->Delta_SO)));
    AA*=pow(1+z,-0.14);
    a=gsl_inter(D,a_par,9,log10((s_cp->Delta_SO)));
    a*=pow(1+z,-0.06);
    b=gsl_inter(D,b_par,9,log10((s_cp->Delta_SO)));
    b*=pow(1+z,-pow(10, -pow(0.75/log10((s_cp->Delta_SO)/75.0),1.2)));
    c=gsl_inter(D,c_par,9,log10((s_cp->Delta_SO)));
    ans= 0.5*AA*(pow(deltac/(sqrt(nu)*b),-a)+1.0)*exp(-c*nu*pow(deltac,-2));
   }
  if("Watson"==mb)
   {
    b=1./1.406;
    ans= 0.282*(pow(deltac/(sqrt(nu)*b),-2.163)+1.0)*exp(-1.2010*nu*pow(deltac,-2));
   }
  return ans;
}




// **********************************************************************************
// **********************************************************************************
/*HALO-MASS BIAS*/
real_prec MASS_BIAS_FUNCTIONS::dm_h_bias(real_prec z, real_prec nu, void *s_p){

  real_prec y,a,b,c,A,B,C;

  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)s_p;
  Cosmology cf(*s_cp);
  string mb=s_cp->halo_mass_bias_fit;

  real_prec n_nu=sqrt(nu);

  real_prec deltac = cf.critical_overdensity(z);

  real_prec dcontrast=s_cp->Delta_SO;//    cf.density_contrast_top_hat(z,s_cp);
  real_prec ans=0;
  if(mb=="Peak_Background")ans= nu/(deltac);
  if(mb=="Mo_White")ans= 1.+(0.75*nu-1)/(deltac); 
  if(mb=="Sheth_Tormen")
  {
    a=0.707;
    b=0.35; //valores de Tinker 2005 0.35; valor original Sheth-Tormen 0.5
    c=0.80; //valores de Tinker 2005 0.80 ; valor original Sheth-Tormen 0.6
    ans= 1.+(1./((deltac)*sqrt(a)))*(sqrt(a)*(a*nu)+b*sqrt(a)*pow(a*nu,1.-c)-pow(a*nu,c)/(pow(a*nu,c)+b*(1.-c)*(1.-0.5*c)));
  }
  if(mb=="Tinker") //https://arxiv.org/pdf/1001.3162.pdf
  {
    y=log10(dcontrast); 
    A=1.0+0.24*y*exp(-pow(4./y,4));
    a=0.44*y-0.88;
    B=0.183;
    b=1.5;
    C=0.019 + 0.107*y + 0.19*exp(-pow(4./y,4));
    c=2.4;
    ans= 1.0-A*pow(n_nu,a)/(pow(n_nu,a)+pow(deltac,a)) +  B*pow(n_nu,b) + C*pow(n_nu,c);
  }
  
  if(mb=="Pillepich")
     ans = 0.647-0.540*pow(deltac/sqrt(nu),-1.0)+1.614*pow(deltac/sqrt(nu),-2.0);
  

  return ans;
}


