//////////////////////////////////////////////////////////
/**
 *  @brief Astrophysics
 *  @file Astrophysics.cpp
 *  @brief Methods of the class Astrophysics
 *  @author Andres Balaguera-Antolínez (ABA)
 */
//////////////////////////////////////////////////////////
# include "Astrophysics.h"
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
real_prec Astrophysics::mass2temp(real_prec x, real_prec z, void *p){
    Astrophysics ap;
    /////////////////////////////////////////////////////////////////
    /// Temperature in Kelvin as a function of circular velocity 
    /// (given in km/s). Mass in units of this code
    /////////////////////////////////////////////////////////////////

    return 0.5*mean_mol_weight*proton_mass*pow(vc(x,z,p),(real_prec)2.0)/Boltzmann_constant;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


real_prec Astrophysics::cooling_function(real_prec T, void*p){
  /*
    Cooling function as a function of the temperature in kelvins, taken from Robison and Silk 2000, after equation 5
    in units of kg m⁵ /s³ (these are J m³ /s)
    An isothermal gas density profile is assumed in order to give an analytial prescription
    x denotes temperature 
  */
  return (4e-36)*pow(T*Boltzmann_constant/kev,(real_prec)-0.5);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::filtering_mass(real_prec z, void *p){
  /*Function implemented to introduce Reionization effects on the Halo Mass function, following Gnedin
    and the Fitting formula of Kravtstov et al 2004
    Filtering mass as a function of redshift in the Code units*/
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec zr=7.0;
  real_prec zo=8.0;
  real_prec alpha_r=6.0;
  real_prec fa;
  real_prec a=1./(1.+z);
  real_prec ar=1./(1.+zr);
  real_prec ao=1./(1.+zo);
  real_prec m_mol_w=0.59;
  real_prec Jeans_Mass=2.5*sqrt(s_cp->Om_matter)*pow(m_mol_w,(real_prec)-1.5)*(1.e11)/(M_reference);
  if(z>=zo){
    fa=(3.*a/((2.+alpha_r)*(5.+2.0*alpha_r)))*pow(a/ao,alpha_r);
  }
  if(z<zo && z>=zr ){
    fa=(3./a)*(pow(1.+zo,-2)*(pow(2.+alpha_r,(real_prec)-1.0)-((2*pow(a*(1+zo),(real_prec)-0.5))/(5.+2.*alpha_r)))+0.1*a*a-0.1*pow(1.+zo,(real_prec)-2.0)*(5.0-4.*pow(a*(1.+zo),(real_prec)-0.5)));
  }
  if(z<zr){
    fa=(3./a)*(pow((real_prec)ao,(real_prec)2)*(pow(2.+alpha_r,(real_prec)-1.0)-((2*pow(a/ao,(real_prec)-0.5))/(5.+2.*alpha_r)))+0.1*pow(ar,(real_prec)2)*(5.0-4.*pow(a/ar,(real_prec)-0.5))-0.1*pow(ao,(real_prec)2)*(5.0-4.*pow(a/ao,(real_prec)-0.5))+(1./3.)*(a*ar)-(1./3.)*pow(ar,(real_prec)2)*(3.-2.*pow(a/ar,(real_prec)-0.5)));
  }
  return Jeans_Mass*fa;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


real_prec Astrophysics::baryon_fraction(real_prec x, real_prec z, void *p){
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  real_prec f_baryon=s_cp->Om_baryons/s_cp->Om_matter;
  real_prec Mass=pow(10,x)*(M_reference);
  real_prec ans= (f_baryon*(s_cp->Om_matter-s_cp->Om_baryons)/s_cp->Om_matter)/pow(1+0.26*filtering_mass(z,p)/Mass,(real_prec)3.0);
  return ans;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::cooling_radius(real_prec x, real_prec l, real_prec z, void *p){
  /*
    Cooling radius determined from the assumptions on White and Frenk 1999, using equations 16 and 17
    in units of Mpc/h
    Up to now I'm not using this expression explicitely anywhere
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(* s_cp);
  real_prec time_scale=cf.halo_dynamical_time(z)*years_to_sec;  /*Halo dynamical time in seconds/h. This quantity is suggested by Crotton, see their discussion after eq. 4*/
  real_prec fbary=baryon_fraction(x,z,p);
  real_prec temp=0.5*mean_mol_weight*proton_mass*pow(vc(x,z,p),(real_prec)2.0)/Boltzmann_constant;  /*en Kelvin*/
  real_prec ans=(1.e-8)*(49.*fbary*s_cp->Om_baryons)*cooling_function(temp, p)/(192.*M_PI*Constants::Gravitational_constant*pow(proton_mass,(real_prec)2));
  /*esta tiene unidades de km² / s : el factor 1e-8 viene del factor 1/(1000)⁵ que sale de pasar las uniades de metros a kilometros, ya que la function de cooling tiene unidades de metros*/
  //  ans=ans*((2./3.)/(Hubble*pow(1+(*red),3./2.)))/Mpc_to_km
  //  ans=sqrt(ans);
  ans=(ans/s_cp->hubble)*(time_scale); /*hago las unidades de tiempo en ans⁻¹ [sec/h]*/
  ans=sqrt(ans);               /*en unidades de km*/
  ans=ans*s_cp->hubble/Constants::Mpc_to_km;           /*en unidades de Mpc/h*/
  return ans;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::vc(real_prec x, real_prec z, void *p){
  /*
    Circular velocity in a dark matter halo observed at the same time of collapse
    as a function of the comoving lagrangian radius rr(x)
    in a matter dominated universe. Taken from  Robinson and Silk 2000
    En unidades de km/s
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  real_prec ans=(1./sqrt(2))*(s_cp->Hubble)*cf.rr(x, z)*sqrt(s_cp->Om_matter)*sqrt(1+z)*pow((cf.density_contrast_top_hat(z))/(cf.omega_matter(z)),(real_prec)1./6.);
  return ans; 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::cooling_mass_rate(real_prec x, real_prec z, void *p){
  /*Cooling mass rate, equation 20 from White and Frenk 1991    in units of (Solar_masses/h) /(years/h)  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  Cosmology cf(*s_cp);
  //return (3./4.)*f_baryon*omegabaryon*Hubble*pow(1+(*red),3./2.)*cooling_radius(x)*pow(vc(x),2)*years_to_sec/Gravitational_constant/Solar_mass;
  real_prec temp=mass2temp(x,z,p);
  /*Equation 8 of Robinson and Silk*/
  //  real_prec time_scale=cf.age_universe();    /*In years/h, Original suggestion from White and Frenk*/
  real_prec time_scale=cf.halo_dynamical_time(z); /*In years/h, Suggested by Crotton, see their discussion after eq. 4*/
  real_prec fbary=baryon_fraction(x,z,p);
  return 0.0045*s_cp->hubble*(s_cp->Om_matter*pow((real_prec)fbary,(real_prec)1.5)/sqrt(s_cp->Hubble*(Constants::years_to_sec/Constants::Mpc_to_km)*time_scale))*sqrt(cooling_function(temp,p)/1e-36)*pow(vc(x,z,p),2);
  /*el factor (years_to_sec/Mpc_to_km) convierte las unidades de Ho a km/s/MPc/h a Mpc/Year /Mpc/h*/
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::infall_mass_rate(real_prec x, real_prec z, void *p){
  //equation 7 of Robison and Silk, in Ms/h /(yrs/h)
  real_prec fbary=baryon_fraction(x,z,p);
  return 0.00024*Constants::sfr_efficiency*fbary*pow(vc(x,z,p),3);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::star_formation_rate(real_prec x, real_prec l, real_prec z, void *p){
  /*infall massa accretion rate in units Solar masses/h /years/h, equation 1 from White and Frenk 1990, 
    the term in the denominator is from their equation 23*/
  //  Minf=hubble*0.15*sfr_efficiency*f_baryon*pow(vc(x),3)/Gravitational_constant*years_to_sec/Solar_mass;
  //  return (vc(x)<min_circular_vel ? 0. : DMIN(infall_mass_rate(x),cooling_mass_rate(x))/(1.+0.02*pow(700./vc(x),2))); 
  return cooling_mass_rate(x,z,p); 
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::baryon_mass_virial_mass_distribution(real_prec xxb, real_prec xx, real_prec z, void *p){
  real_prec mf     =   pow(10,filtering_mass(z,p))*(M_reference);   /*Filtering mass*/
  real_prec mv     =   pow(10,xx)*(M_reference);                   /*DMH mass*/
  real_prec mg_mean=   mv*baryon_fraction(xxb,z,p);                /*Mean Mvir-Mg relation*/
  real_prec mg     =   pow(10,xxb)*(M_reference);                  /*Gas mass*/
  real_prec sigma_mass_reion=mf/(3.*mv);                       /*Intrinsic dispersion*/
  /*Log-normal distribution, Gnedin 2000:*/
  return (1./(sqrt(2.0*M_PI)*sigma_mass_reion))*exp(-0.5*pow(sigma_mass_reion,-2)*pow(log(mg/mg_mean)+0.5*pow(sigma_mass_reion,-2),2));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::total_mass_nb_mass_relation(real_prec x, real_prec z, void *cp, void *ap){
  /*recibe la masa de simulaciones (nb)y de acuerdo a la fraccion de bariones
    como function de la masa total fg(M)g enera la relacion Mt(Mnb) mediante 
    el metodo Newton-Rhapson.
    solucionando Mt=Mnb*((1-fc)/(1-fg))  
    Retorna la masa total en unidades de Ms/h
    Lo que llamamos la masa de las simulaciones es la masa total Mnb=Mdm+fc*Mnb  = Mdm/(1-fc) donde fc
    es la fraccion cosmical de bariones. Esta masa sobre-estima la masa  real de un cluster con un modelo de fgas realista.
  */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)cp;
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  real_prec f_baryon=s_cp->Om_baryons/s_cp->Om_matter;
  real_prec xnb=x;
  real_prec Mass_nb=pow(10,xnb)*(M_reference);
  real_prec masita,Fg,fg,dFg;
  real_prec alpha = (s_ap->A_gas)*(s_cp->f_baryon);
  real_prec beta  = (s_ap->B_gas);
  masita=Mass_nb;
  for(int j=1;j<=20;j++){ /*Newton-Rhapson loop*/
    fg    = alpha*pow(masita/s_ap->mstar,beta);
    Fg    = masita-Mass_nb*(1-s_cp->f_baryon)/(1-fg);
    dFg   = 1.-(1-s_cp->f_baryon)*alpha*beta*pow(masita/s_ap->mstar,beta-1.0)*pow(1.-fg,-2.0)*Mass_nb/s_ap->mstar;
    masita-= Fg/dFg;
  }
  masita=(fg>=(s_cp->f_baryon)? Mass_nb: masita);
  return masita;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//real_prec Astrophysics::Mobs_Mnb_distribution(real_prec x_nb, real_prec xt, real_prec z, void *p){
  // real_prec xx_nb=(M_reference)*pow(10,x_nb);
  // real_prec xxt=(M_reference)*pow(10,xt);
  
  // ************************************************************************************
  // Commented ver si se deja o no
  // ************************************************************************************
  
  
  // /*Do MC realizations to characterize the PDF */
  // real_prec pdf[3];
    // pdf[0]=0;
    // nb2m_mc(x_nb,pdf);
    // real_prec xxt_mean=pdf[1];  /*mean of ln(Mtot(Mnb)) */
    // real_prec ss      =pdf[2];
    // ss=sqrt(pow(scatter_Mobs_Mtot,2)+ss);

    // /*Use fit of sigma_mm from total_mass_nb*cpp and NR solver to avoid the MC realizations*/
    // // Astrophysics ap(x_nb,xt,z);
    // // real_prec xxt_mean=log(ap.total_mass_nb_mass_relation());
    // // real_prec ss=sqrt(pow(scatter_Mobs_Mtot,2)+pow(sigma_mm(x_nb),2));

    // // //Si queremos ignorar los efectos de fgas
    // // real_prec xxt_mean=log(xx_nb);
    // // real_prec ss=scatter_Mobs_Mtot;
    
    // real_prec ans=(1./(sqrt(2.0*M_PI)*ss))*exp(-0.5*pow(ss,-2)*pow(log(xxt)-xxt_mean-log(bias_mass),2.0))/xxt; /*log-normal*/
    // //cout<<ss<<"  "<<xxt_mean<<endl;
    // return ans;
  // return 1.0;
  // }
 

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::mass_lum_distribution(real_prec x, real_prec l, real_prec z, void *p, void *ap){
  /*PROBABILITY DISTRIBUTION FUNCTION FOR P(L|M)*/
  /*Log-normal distribution of mass-luminosities with intrinsec scatter. 
    Constant factors are not included for they calncel in the normalization
    x is the log10(mass/1e12)
    l is the natural log of the luminosity in units of  10^44 erg/sec h^-2 
    Nota importante: esta distribucion es log-normal, de modo que esta definida con p(ln L|M)d \ln L = (p(ln L|M)/L)dL */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  real_prec y         = pow(10,x)*(M_reference)/(1e14); //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  ScalingRelations ml;
  return      exp(-0.5*pow(l-ml.MOCKS_M2L(y,p),2)/(pow(s_ap->sigma_ln,2)))/exp(l);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::mass_lum_distribution_errors(real_prec x, real_prec l, real_prec z, void *cp, void *ap){
  /*this is the convolution of p(L|M) with another log-normal distribution with sigma for the flux errors
    Constant factors are not included for they cancel in the normalization
    x is the log10(mass/1e12)
    l is the natural log of the luminosity in units of  10^44 erg/sec h^-2 
    include missing flux correction in the factor missing flux */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  real_prec y =  pow(10,x)*(M_reference)/(1e14); //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  ScalingRelations ml;
  l =  log(exp(l)/s_ap->missing_flux);
  return exp(-0.5*pow(l-ml.MOCKS_M2L(y,cp),2)/pow(s_ap->sigma_red,2))/exp(l);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::G(real_prec x, real_prec l, real_prec z, void *p, void *ap){
  /*Analytic integration of the function  mass_lum_distribution_errors for objects with lum greater than L. 
    No coloco los factores multiplicativos pues se cancelan */
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  real_prec factores     =  1./sqrt(2.*s_ap->sigma_red*s_ap->sigma_red);
  real_prec y=pow(10,x)*(M_reference)/(1e14);             //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  ScalingRelations ml;
  return gsl_sf_erf(factores*(l-ml.MOCKS_M2L(y,p)));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Astrophysics::G_Lmin_Lmax(real_prec x, real_prec l, real_prec z, void *cp, void *ap){
    /*Analytic integration of the function  mass_lum_distribution_errors for objects with lum greater than L. 
      and smaller than a fixed value Lmax.
      No coloco los factores multiplicativos pues se cancelan en las funciones de normalizacion*/
  struct s_astrophysical_parameters * s_ap= (struct s_astrophysical_parameters *)ap;
  real_prec factores     =  1./sqrt(2.*s_ap->sigma_red*s_ap->sigma_red);
  real_prec y=pow(10,x)*(M_reference)/(1e14);             //ML RELATIONS RECIVE MASAS EN UNIDADES 1E14/h
  real_prec lmin=l;
  real_prec lmax=z;
  real_prec interval=(exp(lmax)-exp(lmin));
  ScalingRelations ml;
  return 0.5*(gsl_sf_erf(factores*(lmax-ml.MOCKS_M2L(y,cp)))-gsl_sf_erf(factores*(lmin-ml.MOCKS_M2L(y,cp))))/interval;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


