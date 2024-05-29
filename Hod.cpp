////////////////////////////////////////////////////////////////////////////
/**
 * @brief This file contains methods for the HOD class 
 * @file HOD.cpp
 * @author Andres Balaguera Antolinez
 * @version 1.0
 * @date 2017-2024
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "../headers/Constants.h"
# include "../headers/Hod.h"
////////////////////////////////////////////////////////////////////////////
real_prec HOD::CENTRAL(real_prec M){
  real_prec f;
  real_prec mmin=this->s_cosmo_pars.mmin_hod;
  real_prec sm=this->s_cosmo_pars.scatter_hod;
  switch(this->s_cosmo_pars.hod_model){
    case(1):
      f=(M<mmin? 0.0 : 1.0);
      break;
    case(2):
      f=exp(-mmin/M);
      break;
    case(3):
      f=0.5*(1.+gsl_sf_erf(log10(M/mmin)/sm));
      break;
    case(4):
      f=0.5*(1.+gsl_sf_erf(log10(M/Mmin_hod(M))/sigma_hod(M)));
    break;
  }
  return f;
}
////////////////////////////////////////////////////////////////////////////
real_prec HOD::SATELLITE(real_prec M){
  real_prec f=1;
  real_prec kappa=0.137;
  real_prec mmin=this->s_cosmo_pars.mmin_hod;
  real_prec muno=this->s_cosmo_pars.muno_hod;
  real_prec Mstep=this->s_cosmo_pars.Mstep_hod;
  real_prec al=this->s_cosmo_pars.alpha_hod;
  switch(this->s_cosmo_pars.hod_model){
    case(1):
      f=pow(M/muno,al)*this->CENTRAL(M);
      break;
    case(2):
      f=pow(M/muno,al)*this->CENTRAL(M);
      break;
    case(3):
      f=(M<mmin? 0.0 : pow((M-kappa*mmin)/muno,al))   ;
      break;
    case(4):
      f= M<Mstep? 0.0 : pow((M-Mstep)/muno,al)*this->CENTRAL(M);
  }
  return f;
}
////////////////////////////////////////////////////////////////////////////
real_prec HOD::Mmin_hod(real_prec Mh)
{
    real_prec Ms=Mhalo_to_Mstellar(Mh);
    return this->s_cosmo_pars.muno_hod;
}
////////////////////////////////////////////////////////////////////////////
real_prec HOD::sigma_hod(real_prec Mh)
{
    real_prec s_faint=this->s_cosmo_pars.s_faint_hod;
    real_prec s_bright=this->s_cosmo_pars.s_bright_hod;
    real_prec lMs = log10(this->Mhalo_to_Mstellar(Mh));
    real_prec Mstep = log10(this->s_cosmo_pars.Mstep_hod);
    real_prec width = this->s_cosmo_pars.width_hod;
    real_prec sigma=s_faint + (s_bright-s_faint)/(1.+exp((Mstep - lMs)*width));
    return sigma;
}
// **************************************|****************************
real_prec HOD::Mhalo_to_Mstellar(real_prec Mh)
{
    real_prec alpha=this->alpha_par();
    real_prec M_s= this->s_cosmo_pars.ms_hod;                          //
    real_prec M_t = this->s_cosmo_pars.mt_hod; // transition mass
    real_prec s_mass = M_s*pow(Mh/M_t,alpha)*exp(1.-M_t/Mh);
    return s_mass;
}
// ******************************************************************
real_prec HOD::alpha_par()
{

    real_prec ms= log10(this->s_cosmo_pars.ms_hod);
    real_prec al = log10(this->s_cosmo_pars.alpha_C+pow(s_cosmo_pars.alpha_A*ms,this->s_cosmo_pars.alpha_B));
    return al;
}
// ******************************************************************

