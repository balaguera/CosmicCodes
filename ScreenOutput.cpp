////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This file contains functions for numerical analysis
 * @file ScreenOutoput.h
 * @author Andres Balaguera Antolinez
 * @version 1.0
 * @date 2017-2024
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "../headers/ScreenOutput.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::show_warnings(){
#ifdef mode_p
  if(VEL_BIAS_POWER!=1.)
    this->message_warning("Parameter VEL_BIAS_POWER in def.h different of 1. Change if no vel bias in the power spectrum is to be used.");
#endif
#ifdef _REDSHIFT_SPACE_
    this->message_warning("Power to be measured in redshift space.");
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message(string mess){std::cout<<mess<<endl;}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<"\t"<<YELLOW<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tPOWER SPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP (+FFTW) ESTIMATOR            *"<<endl;
 std::cout<<"\tVERSION 1.3                                                                        *"<<endl;
 std::cout<<"\tFor more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tStarting time and date "<<ctime (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_c(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<"\t"<<YELLOW<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tANALYSIS OF DARK MATTER TRACERS - ASSEMBLY BIAS                                    *"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tStarting time and date "<<ctime (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_cl(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<"\t"<<YELLOW<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tANGULAR POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                *"<<endl;
 std::cout<<"\tUSING HEALPIX DECOMPOSITION AND PEEBLES ESTIMATOR                                  * "<<endl;
 std::cout<<"\tVERSION 1.0                                                                        *"<<endl ;
 std::cout<<"\tFor more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tStarting time and date"<<ctime (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_fb(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<"\t"<<YELLOW<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t3D POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                     *"<<endl;
 std::cout<<"\tUSING THE FOURIER BESSEL DECOMPOSITION                                             * "<<endl;
 std::cout<<"\tVERSION 1.0                                                                        *"<<endl ;
 std::cout<<"\tFor more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tStarting time "<<ctime (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_yama(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<"\t"<<YELLOW<<endl;
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<"\tMULTIPOLE DECOMPOSITION OF THE THREE DIMENSIONAL POWER SPECTRUM OF COSMOLOGICAL   *"<<endl;
  std::cout<<"\tMASS TRACERS USING YAMAMOTO ESTIMATOR (FFTW-based)                                *"<<endl;
  std::cout<<"\tVERSION 1.1                                                                       *"<<endl;
  std::cout<<"\tFor more documentation see README/README.pdf                                      *"<<endl;
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<"\tStarting time "<<ctime (&rawtime);
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_bispectrum(){
 std::cout<<"\t"<<YELLOW<<endl;
 time_t rawtime;
 time (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tBISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl;
 std::cout<<"\tVERSION 1.1                                                                        *"<<endl;
 std::cout<<"\tDocumentation and warnings in readme.ps                                            *"<<endl;
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<"\tStarting time "<<ctime (&rawtime);
 std::cout<<"\t************************************************************************************"<<endl;
 std::cout<<RESET;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_bispectrum_fast()
{
  time_t rawtime;
  time (&rawtime);
  std::cout<<"\t"<<YELLOW<<endl;
  std::cout<<"\t************************************************************************************"<<endl;
  std::cout<<"\t************************************************************************************"<<endl;
  std::cout<<"\tBISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl;
  std::cout<<"\tAND INVERSE FFTW TRICK                                                             *"<<endl;
  std::cout<<"\tVERSION 1.0                                                                        *"<<endl;
  std::cout<<"\tDocumentation and warnings in readme.ps                                            *"<<endl;
  std::cout<<"\t************************************************************************************"<<endl;
  std::cout<<"\tStarting time "<<ctime (&rawtime);
  std::cout<<"\t************************************************************************************"<<endl;
  std::cout<<RESET;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_interrupt()
{
  std::cout<<"\t"<<YELLOW<<endl;
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<"\tCode interrupted at"<<endl;
  time_t rawtime;
  time ( &rawtime );
  std::cout<<ctime (&rawtime);
  std::cout<<"\t***********************************************************************************"<<endl;
  std::cout<<"\t"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::enter(string text)
{
#ifdef _FULL_VERBOSE_
//  string time_now=ctime(&this->initial_time);
 // std::cout<<BOLDGREEN<<"\t"<<time_now<<"\t=============> Going to:"<<GREEN<<text<<RESET<<endl;
  std::cout<<BOLDGREEN<<"\t=============> Going to:"<<GREEN<<text<<RESET<<endl;
 // time(&this->initial_time);
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::leaving(string text)
{
#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<"\t==================================> leaving"<<GREEN<<text<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::ending_message()
{
  std::cout<<YELLOW<<endl;
  time_t rawtime;
  time ( &rawtime );
  std::cout<<"\tEnding date:"<<ctime (&rawtime)<<endl;
  std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::write_cosmo_parameters(void *p, void *pa)
{
  std::cout<<endl;  
  this->message_screen("\tCosmological Parameters");
#ifdef _USE_UNITSIM_COSMOLOGY_
  this->message_screen("\tUsing cosmological parameters from UNITsim (Planck 16)");
#endif
#ifdef _USE_SLICS_COSMOLOGY_
  this->message_screen("\tUsing cosmological parameters from SLICS (Planck 16)");
#endif
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  struct s_CosmoInfo * s_ci= (struct s_CosmoInfo *)pa;
  this->message_screen("Input redshift z =", s_cp->cosmological_redshift);
  this->message_screen("Omega Matter =", s_cp->Om_matter);
  this->message_screen("Omega CDM =", s_cp->Om_cdm);
  this->message_screen("Omega Vacuum =", s_cp->Om_vac);
  this->message_screen("Omega Baryons =", s_cp->Om_baryons);
  this->message_screen("Omega Curvature =", s_cp->Om_k);
  this->message_screen("Omega Radiation =", s_cp->Om_radiation);
  this->message_screen("DE eos =",s_cp->wde_eos);
  this->message_screen("Hubble parameter h =", s_cp->Hubble);
  this->message_screen("hubble parameter H =", s_cp->hubble);
  this->message_screen("Spectra al Index =", s_cp->spectral_index);
  this->message_screen("Sigma 8 =", s_cp->sigma8);
  this->message_screen("Primordial amplitude =", s_cp->A_s);
  this->message_screen("Running index =", s_cp->alpha_s);
  this->message_screen("Use wiggles =",s_cp->use_wiggles);
  this->message_screen("Mean CMB temperature =",s_cp->Tcmb);
  std::cout<<endl;
  this->message_screen("Derived quantities");
  cout<<endl;
  this->message_screen("Scale factor a =", 1./(s_cp->cosmological_redshift+1.));
  this->message_screen("Growth factor g(z)=D1(z)/D1(0) =",s_ci->growth_factor);
  this->message_screen("Growth factor D2(z) =",s_ci->D2);
  this->message_screen("Growth index f(z) =",s_ci->growth_index);
  this->message_screen("H(z) =",s_ci->Hubble_parameter, "(km/s)/(Mpc/h)");
  this->message_screen("Comoving distance r(z) =",s_ci->comoving_distance, "Mpc/h");
  this->message_screen("Comoving AD distance d(z) =", s_ci->comoving_angular_diameter_distance, "Mpc/h");
  this->message_screen("Comoving sound horizon rs(z) =",s_ci->comoving_sound_horizon, "Mpc/h");
  this->message_screen("Mean matter density =",s_ci->mean_matter_density, "(Ms/h)/(Mpc/h)^3");
  this->message_screen("Age of the Universe =",s_ci->age_universe, "Years/h");
  this->message_screen("Omega matter (z) =",s_ci->omega_matter);
#ifdef _USE_LPT_
  this->message_screen("Normalization of Linear Matter P(k,z=0) =",s_cp->pk_normalization);
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::write_cosmo_parameters(void *p)
{
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  std::cout<<BOLDCYAN<<"\tInput Cosmological parameters"<<RESET<<endl;
  this->message_screen("\tOmega matter =", s_cp->Om_matter);
  this->message_screen("\tOmega vac =", s_cp->Om_vac);
  this->message_screen("\tOmega baryons =", s_cp->Om_baryons);
  this->message_screen("\tOmega curv =", s_cp->Om_k);
  this->message_screen("\tHubble par =", s_cp->Hubble);
  this->message_screen("\tSpectral index =", s_cp->spectral_index);
  this->message_screen("\tSigma8 =", s_cp->sigma8);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::error_ncolumns(string fname){
  std::cout<<RED;
  std::cout<<fname<<"\t with less than 4 columns"<<endl;
  char yn;
  ScreenOutput So;
  std::cout<<"\tContinue (y/n)?"<<endl; cin>>yn;
  if(yn !='y' && yn !='n')
    {
      std::cout<<"\tPlease answer  yes (y) or not (n)"<<endl;
      std::cout<<"\tContinue ?"<<endl; cin>>yn;
    }
  else{
    if(yn !='y'){
    }
    So.message_interrupt();
  }
  std::cout<<RESET;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::comp_time(time_t start, unsigned long full, unsigned long step){
  // THIS IS TIME CONSUMING
  double fraction=100.0*((double)(step))/((double)full);
  time_t end;
  time(&end);
  double lapse=difftime(end,start);
  if(lapse<=60)std::cout<<"\t\r "<<fraction<<"\t % completed. Time elapsed: "<<lapse<<"\tsecs \r";std::cout.flush();
  if(lapse>60)std::cout<<"\t\r "<<fraction<<"\t % completed. Time elapsed: "<<lapse/60.<<"\tmins \r";std::cout.flush();
  if(lapse>3600)std::cout<<"\t\r "<<fraction<<"\t % completed. Time elapsed: "<<lapse/3600.<<"\thrs \r";std::cout.flush();
  if (fraction==25.0) std::cout <<"\t..25%  completed"<<endl;
  if (fraction==50.0) std::cout <<"\t..50%  completed"<<endl;
  if (fraction==75.0) std::cout <<"\t..75%  completed"<<endl;
  if (fraction==100.0)std::cout <<"\t..100% completed"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss)
{
  std::cout<<"\t"<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<"\t"<<GREEN<<ss<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss,string sa)
{
  std::cout<<"\t"<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<GREEN<<"\t"<<ss<<RESET<<" "<<sa<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss, int line)
{
  std::cout<<"\t"<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<GREEN<<"\t"<<ss<<"\t"<<line<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning_ini(int line, string fun,string file, string war)
{
  std::cout<<"\t"<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<RED<<"\tfrom line "<<line<<"\t in file "<<file<<"\t, function "<<fun<<RESET<<std::endl;
  std::cout<<GREEN<<war<<RESET<<std::endl;
  std::cout<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning_ini(int line, string fun,string file, string war, int ii)
{
  std::cout<<"\t"<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<RED<<"\tWarning from line "<<line<<"\tin function "<<fun<<"\tfrom file "<<file<<RESET<<std::endl;
  std::cout<<GREEN<<"\t"<<war<<"\t"<<ii<<RESET<<std::endl;
  std::cout<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error(string ss)
{
  std::cout<<"\t"<<CYAN<<"ERROR"<<RESET<<std::endl;
  std::cerr<<CYAN<<ss<<RESET<<std::endl;
  std::cerr<<RED<<"\t...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error_file(string ss)
{
  std::cerr<<CYAN<<"\tFile "<<ss<<"\t does not exist."<<RESET<<std::endl;
  std::cerr<<GREEN<<"\tPlease check parameter file and paths."<<endl;
  std::cerr<<GREEN<<"\tCosmicatlass ends here."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput:: message_error(string ss, double a)
{
  std::cerr<<CYAN<<ss<<"\t"<<a<<RESET<<std::endl;
  std::cerr<<RED<<"\t...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error(string ss, double a, string sa)
{
  std::cerr<<CYAN<<ss<<"\t"<<a<<"\t"<<sa<<RESET<<std::endl;
  std::cerr<<RED<<"\t...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss)
{
  std::clog<<"\t"<<YELLOW<<ss<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////

void ScreenOutput::message_screen_flush(string ss, real_prec s2, string sa, real_prec s3)
{
  std::cout<<"\t"<<"\r"<<YELLOW<<ss<<RESET<<"\t"<<s2<<"\t"<<sa<<"\t"<<s3<<RESET; std::cout.flush();
}
////////////////////////////////////////////////////////////////////////////

void ScreenOutput::message_screen_flush(string ss, real_prec s2, string sa, real_prec s3, string s4, real_prec s5)
{
  std::cout<<"\t"<<"\r"<<YELLOW<<ss<<RESET<<"\t"<<s2<<"\t"<<sa<<"\t"<<s3<<"\t"<<s4<<"\t"<<s5<<RESET; std::cout.flush();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen_flush(string ss, int s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<"\t"<<"\r"<<YELLOW<<ss<<RESET<<"\t"<<s2<<RESET;std::cout.flush();
#else
    std::cout<<RESET<<ss<<"\t"<<s2<<RESET<<std::endl;
#endif
}
void ScreenOutput::message_screen_flush(string ss, ULONG s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<"\t"<<"\r"<<YELLOW<<ss<<RESET<<"\t"<<s2<<RESET;std::cout.flush();
#else
    std::cout<<RESET<<ss<<"\t"<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen_flush(string ss, real_prec s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<"\t\r"<<YELLOW<<ss<<RESET<<"\t"<<s2<<RESET;
  std::cout.flush();
#else
    std::cout<<RESET<<ss<<"\t"<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double s2)
{
#ifdef _FULL_VERBOSE_
    std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<s2<<RESET<<std::endl;
#else
    std::cout<<"\t"<<RESET<<ss<<"\t"<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa,double s2)
{
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<i<<"\t"<<sa<<"\t  "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa, ULONG s2){
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<i<<"\t"<<sa<<"\t  "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa, ULONG s2,string s3){
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<i<<"\t"<<sa<<"\t"<<s2<<"\t"<<s3<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int s2)
{
  std::clog<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, ULONG s2)
{
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, string s2){
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double a, string s2, string f){
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<a<<"\t"<<s2<<"\t"<<f<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, string s2)
{
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<d2<<" "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, ULONG d2, string s2)
{
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<d2<<"\t"<<s2<<RESET<<std::endl;
}
void ScreenOutput::message_screen(string ss, int d2, string s2)
{
  std::cout<<"\t"<<YELLOW<<ss<<RESET<<"\t"<<CYAN<<d2<<"\t"<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, time_t time)
{
  double lapse=difftime(time,this->initial_time);
  if(lapse<60)
    std::clog<<RED<<ss<<"\t"<<d2<<"\t(T_init + "<<lapse<<" secs)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<RED<<ss<<"\t"<<d2<<"\t(T_init + "<<lapse/60<<" mins)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<RED<<ss<<"\t"<<d2<<"\t(T_init + "<<lapse/3600<<" hrs)"<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, time_t time, time_t time2)
{
  double lapse=difftime(time,this->initial_time);
  double lapse2=difftime(time,time2);

  if(lapse<60)
    std::clog<<BOLDCYAN<<ss<<"\t"<<d2<<"\t(T_init +"<<lapse<<" secs, time spent in last step "<<lapse2<<"\t s)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<BOLDCYAN<<ss<<"\t"<<d2<<"\t(T_init +"<<lapse/60<<" mins, time spent in last step "<<lapse2<<"\t s)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<BOLDCYAN<<ss<<"\t"<<d2<<"\t(T_init +"<<lapse/3600<<" hrs, time spent in last step "<<lapse2<<"\t s)"<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, string s2, string s3, double d3)
{
  std::clog<<GREEN<<ss<<CYAN<<"\t"<<s2<<"\t"<<BLUE<<s3<<"\t"<<CYAN<<d3<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time(time_t start_all)
{
  time_t end;
  time(&end);
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
#endif
  double lapse=difftime(end,start_all);
  string units_time="secs";
  if(lapse>60)
    {
      lapse/=60.0;
      units_time="mins";
    }
  else if(lapse>3060.0)
      units_time="hrs";
  std::cout<<__DATE__<<"\t"<<__TIME__<<endl;
  std::cout<<GREEN<<"\tComputing Time: "<<lapse<<"\t"<<units_time<<RESET<<endl;
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<"\t**************************************************************************"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time_mock(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  std::cout<<"\t"<<YELLOW<<"\tMock DF generated in : "<<lapse<<"\tsecs"<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time2(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  std::cout<<"\t"<<YELLOW<<"Computing Time: "<<CYAN<<lapse<<"\tsecs"<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::usage(string s)
{
  std::cout<<BOLDGREEN;  
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\tCOSMICATLASS"<<endl;
  std::cout<<"\t\n\tCosmologicalCATalogs for LArge Scale Structure"<<endl;
  std::cout<<"\t\n\tHow to run: "<<s<<"\t [-option] [argument]"<<endl;
  std::cout<<"\t\n\tOptions: "<<endl;
  std::cout<<"\t         -a for information on the author, no argument"<<endl;
  std::cout<<"\t         -b parameter_file.ini: Runs BiasMT"<<endl;
  std::cout<<"\t         -c parameter_file.ini: Analyzes tracer catalog"<<endl;
  std::cout<<"\t         -d parameter_file.ini: Shows preprocessor directives"<<endl;
  std::cout<<"\t         -f parameter_file.ini: Reads binarys form IC and cmpute density field IC"<<endl;
  std::cout<<"\t         -g parameter_file.ini: Generates galaxy catalog from halo catalog using HOD"<<endl;
  std::cout<<"\t         -h parameter_file.ini: Help"<<endl;
  std::cout<<"\t         -i parameter_file.ini: Shows input pars"<<endl;
  std::cout<<"\t         -m parameter_file.ini: Measures marked power spectrum"<<endl;
  std::cout<<"\t         -p parameter_file.ini: Measures Statistic"<<endl;
  std::cout<<"\t         -q parameter_file.ini: SuperClusters"<<endl;
  std::cout<<"\t         -s parameter_file.ini: Applies a low-pass filter to input density field"<<endl;
  std::cout<<"\t         -t parameter_file.ini: Runs LPT"<<endl;
  std::cout<<"\t         -u parameter_file.ini: Analyzes input halo catalog and measurements of secondary bias"<<endl;
  std::cout<<"\t         -w parameter_file.ini: Computes window functions. Under construction"<<endl;
  std::cout<<"\t         -v parameter_file.ini: Warnings "<<endl;
  std::cout<<"\t         -z parameter_file.ini: Assign,ment of individual tracer bias"<<endl;
  std::cout<<"\tConsult ../Headers/def.h for pre-procesor directives"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<RESET<<endl;
  std::cout<<RESET<<endl;                                       
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::author(){
  std::cout<<BOLDGREEN;  
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<BOLDCYAN;
  std::cout<<"\t*   ****  **  |*** *      *  *   ****    *   ***** *       *    |***  *"<<endl;
  std::cout<<"\t*  |     *  *  *   * *  * *  *  *       * *    *   *      * *    *    *"<<endl;
  std::cout<<"\t*  |     *  *   *| *  **  *  *  *      *****   *   *     *****    *|  *"<<endl;
  std::cout<<"\t*   ****  **  **** *      *  *   **** *     *  *   **** *     * ****  *"<<RESET<<endl;
  std::cout<<BOLDGREEN;  
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\tCOSMICATLASS                                                          *"<<endl;
  std::cout<<"\tAndres Balaguera-Antolinez (balaguera@iac.es)                         *"<<endl;
  std::cout<<"\tInstituto de Astrofisica de Canarias                                  *"<<endl;
  std::cout<<"\t2011-2023                                                             *"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<"\t***********************************************************************"<<endl;
  std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void  ScreenOutput::message(time_t start_all)
{
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\tCOSMICATLASS"<<endl;
  std::cout<<"\tCosmological CATalogs for LArge Scale Structure"<<endl;
  std::cout<<"\tIAC 2017-2023"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\tCompilation began on "<<__DATE__<<","<<__TIME__<<endl;
#ifdef _USE_OMP_
  std::cout<<"\tUsing OMP with "<<_NTHREADS_<<" threads"<<endl;
#endif
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<RESET<<endl;
#else
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  std::cout<<"\tBiasMT. Realization "<<this->params._realization()<<"\t p-index (Unitsim) "<<this->params._unitsim_plabel()<<endl;
#endif
  std::cout<<"\t*****************************************************************"<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  log<<"\t*****************************************************************"<<endl;
  log<<"\t*****************************************************************"<<endl;
  log<<"\tBiasMT. Launched at "<<__DATE__<<"\t at "<<__TIME__<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  log<<"\tBiasMT. Realization "<<this->params._realization()<<"\t p-index (Unitsim) "<<this->params._unitsim_plabel()<<endl;
#endif
  log<<"\t*****************************************************************"<<endl;
  log<<"\t*****************************************************************"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void  ScreenOutput::message_BiasMT(time_t start_all)
{
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\tCOSMICATLASS"<<endl;
  std::cout<<"\tCosmological CATalogs for LArge Scale Structure"<<endl;
  std::cout<<"\tIAC 2017-2023"<<endl;
  std::cout<<"\t****************************************************************************************"<<endl;
  std::cout<<"\tCompilation began on "<<__DATE__<<","<<__TIME__<<endl;
#ifdef _USE_OMP_
  std::cout<<"\tUsing OMP with "<<_NTHREADS_<<" threads"<<endl;
#endif
  std::cout<<"\t*****************************************************************"<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  std::cout<<"\tBiasMT. Realization "<<this->params._realization()<<"\t p-index (Unitsim) "<<this->params._unitsim_plabel()<<endl;
#endif
  std::cout<<RESET<<endl;
#else
  std::cout<<"\t*****************************************************************"<<endl;
  std::cout<<"\t*****************************************************************"<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  std::cout<<"\tBiasMT. Realization "<<this->params._realization()<<"\t p-index (Unitsim) "<<this->params._unitsim_plabel()<<endl;
#endif
  std::cout<<"\t*****************************************************************"<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  log<<"\t*****************************************************************"<<endl;
  log<<"\t*****************************************************************"<<endl;
  log<<"\tBiasMT. Launched at "<<__DATE__<<"\t at "<<__TIME__<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  log<<"\tBiasMT. Realization "<<this->params._realization()<<"\t p-index (Unitsim) "<<this->params._unitsim_plabel()<<endl;
#endif
  log<<"\t*****************************************************************"<<endl;
  log<<"\t*****************************************************************"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_output_file(string s, ULONG l)
{
#ifdef _FULL_VERBOSE_
  std::cout<<YELLOW<<"\tWritting output file "<<CYAN<<s<<" with "<<l<<" lines"<<RESET<<endl;
#endif
  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  if(log.is_open())
  {
     log<<"\tWritting output file "<<s<<" with "<<l<<" lines"<<endl;
     log.close();
  }
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_output_file(string s, int nlines, int ncols)
{
#ifdef _FULL_VERBOSE_
  std::cout<<"\t"<<YELLOW<<"Writting output file "<<CYAN<<s<<" with "<<nlines<<" lines and "<<ncols<<" columns"<<RESET<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(),ios_base::app);
  log<<"\tWritting output file "<<s<<" with "<<nlines<<" lines and "<<ncols<<" columns"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::DONE()
{
#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::show_preproc()
{
  message_screen("******************************");
  message_screen("******************************");
  message_screen("******************************");
  std::cout<<WHITE<<"\n\tPRE-PROCESSOR DIRECTIVES IN BiasMT"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<"\tDEFINED"<<RESET;
  std::cout<<COLOR_UNDEFINED<<"\t  UNDEFINED"<<RESET<<endl;
  message_screen("******************************");
  message_screen("PRECISION");
#ifdef SINGLE_PREC
  std::cout<<COLOR_DEFINED<<"\tSINGLE_PREC"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\tSINGLE_PREC"<<RESET<<RESET<<endl;
#endif
#ifdef DOUBLE_PREC
  std::cout<<COLOR_DEFINED<<"\tDOUBLE_PREC"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\tDOUBLE_PREC"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OMP");
#ifdef _USE_OMP_
  std::cout<<COLOR_DEFINED<<"\t_USE_OMP_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_OMP_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("GNUPLOT");
#ifdef _USE_GNUPLOT_
  std::cout<<COLOR_DEFINED<<"\t_USE_GNUPLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_GNUPLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_
  std::cout<<COLOR_DEFINED<<"\t_SHOW_BIAS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_SHOW_BIAS_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
  std::cout<<COLOR_DEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_V_PLOT_
  std::cout<<COLOR_DEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_V_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_V_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_RS_PLOT_
  std::cout<<COLOR_DEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_RS_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_RS_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_
  std::cout<<COLOR_DEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("VERBOSE");
#ifdef _USE_COLORS_
  std::cout<<COLOR_DEFINED<<"\t_USE_COLORS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_COLORS_"<<RESET<<RESET<<endl;
#endif
#ifdef _FULL_VERBOSE_
  std::cout<<COLOR_DEFINED<<"\t_FULL_VERBOSE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_FULL_VERBOSE_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("IO FORMATS");
#ifdef _WRITE_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _READ_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<"\t_READ_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_READ_BINARY_BiasMT_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _OUTPUT_WITH_HEADERS_
  std::cout<<COLOR_DEFINED<<"\t_OUTPUT_WITH_HEADERS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_OUTPUT_WITH_HEADERS_"<<RESET<<RESET<<endl;
#endif
#ifdef _WRITE_COORDINATES_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_COORDINATES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_COORDINATES_"<<RESET<<RESET<<endl;
#endif
#ifdef _WRITE_VELOCITIES_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_VELOCITIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_VELOCITIES_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("REFERENCE SIMULATION");
#ifdef _ABACUS_
  std::cout<<COLOR_DEFINED<<"\t_ABACUS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ABACUS_"<<RESET<<RESET<<endl;
#endif
#ifdef _MINERVA_
  std::cout<<COLOR_DEFINED<<"\t_MINERVA_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_MINERVA_"<<RESET<<RESET<<endl;
#endif
#ifdef _SLICS_
  std::cout<<COLOR_DEFINED<<"\t_SLICS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_SLICS_"<<RESET<<RESET<<endl;
#endif
#ifdef _UNITSIM_
  std::cout<<COLOR_DEFINED<<"\t_UNITSIM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_UNITSIM_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("COSMOLOGICAL PARAMETERS");
#ifdef _USE_COSMO_PARS_
  std::cout<<COLOR_DEFINED<<"\t_USE_COSMO_PARS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_COSMO_PARS_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_UNITSIM_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<"\t_USE_UNITSIM_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_UNITSIM_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_PLANCK_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<"\t_USE_PLANCK_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_PLANCK_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_SLICS_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<"\t_USE_SLICS_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_SLICS_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OPERATING MODE");
#ifdef BIAS_MODE
  std::cout<<COLOR_DEFINED<<"\tBIAS_MODE"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\tBIAS_MODE"<<RESET<<endl;
#endif
#ifdef mode_p
  std::cout<<COLOR_DEFINED<<"\tmode_p"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\tmode_p"<<RESET<<endl;
#endif
#ifdef mode_b
  std::cout<<COLOR_DEFINED<<"\tmode_b"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\tmode_b"<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("BiasMT");
/*
#ifdef MOCK_MODE
#ifdef _DO_BiasMT_CALIBRATION_
  std::cout<<COLOR_DEFINED<<"\t_DO_BiasMT_CALIBRATION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_DO_BiasMT_CALIBRATION_"<<RESET<<endl;
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
  std::cout<<COLOR_DEFINED<<"\t_GET_BiasMT_REALIZATIONS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_GET_BiasMT_REALIZATIONS_"<<RESET<<endl;
#endif
#ifdef _USE_VELOCITIES_
  std::cout<<COLOR_DEFINED<<"\t_USE_VELOCITIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_VELOCITIES_"<<RESET<<endl;
#endif
#endif
*/
  std::cout<<endl;
message_screen("******************************");
message_screen("LPT");
#ifdef _USE_LPT_
  std::cout<<COLOR_DEFINED<<"\t_USE_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_LPT_"<<RESET<<endl;
#endif
#ifdef _ONLY_LPT_
  std::cout<<COLOR_DEFINED<<"\t_ONLY_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_ONLY_LPT_"<<RESET<<endl;
#endif
#ifdef _POWER_BOOST_ALPT_
  std::cout<<COLOR_DEFINED<<"\t_POWER_BOSST_ALPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_POWER_BOOST_ALPT_"<<RESET<<endl;
#endif
#ifdef _GET_VELOCITIES_LPT_
  std::cout<<COLOR_DEFINED<<"\t_GET_VELOCITIES_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_GET_VELOCITIES_LPT_"<<RESET<<endl;
#endif
#ifdef _GET_DENSITY_FIELDS_
  std::cout<<COLOR_DEFINED<<"\t_GET_DENSITY_FIELDS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_GET_DENSITY_FIELDS_"<<RESET<<endl;
#endif

std::cout<<endl;
message_screen("******************************");
message_screen("DISPLACEMENTS");

#ifdef _DISPLACEMENTS_
  std::cout<<COLOR_DEFINED<<"\t_DISPLACEMENTS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_DISPLACEMENTS_"<<RESET<<endl;
#endif
std::cout<<endl;
  message_screen("******************************");
message_screen("CALIBRATION");

#ifdef _SMOOTHED_KERNEL_
  std::cout<<COLOR_DEFINED<<"\t_SMOOTHED_KERNEL_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_SMOOTHED_KERNEL_"<<RESET<<endl;
#endif
#ifdef _MODIFY_LIMITS_
  std::cout<<COLOR_DEFINED<<"\t_MODIFY_LIMITS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_MODIFY_LIMITS_"<<RESET<<endl;
#endif
#ifdef _RANK_ORDERING_AB_INITIO_
  std::cout<<COLOR_DEFINED<<"\t_RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#endif
#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _RHO_WITH_DELTA_
  std::cout<<COLOR_DEFINED<<"\t_RO_WITH_DELTA_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_RO_WITH_DELTA_"<<RESET<<endl;
#endif
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("DARK MATTER (THETA) MODELS AND PROPERTIES");
#ifdef _USE_DM_IN_BiasMT_
  std::cout<<COLOR_DEFINED<<"\t_USE_DM_IN_BiasMT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_DM_IN_BiasMT_"<<RESET<<endl;
#endif

#ifdef _USE_TWEB_
  std::cout<<COLOR_DEFINED<<"\t_USE_TWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_TWEB_"<<RESET<<endl;
#endif

#ifdef _USE_IWEB_
  std::cout<<COLOR_DEFINED<<"\t_USE_IWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_IWEB_"<<RESET<<endl;
#endif

#ifdef _USE_IKWEB_
  std::cout<<COLOR_DEFINED<<"\t_USE_IKWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_IKWEB_"<<RESET<<endl;
#endif


#ifdef _USE_CWC_
  std::cout<<COLOR_DEFINED<<"\t_USE_CWC_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_CWC_"<<RESET<<endl;
#endif
#ifdef _USE_CWC_INSIDE_LOOP_
  std::cout<<COLOR_DEFINED<<"\t_USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_KNOTS_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_KNOTS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_KNOTS_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  std::cout<<COLOR_DEFINED<<"\t_USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  std::cout<<COLOR_DEFINED<<"\t_USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_TIDAL_ANISOTROPY_
  std::cout<<COLOR_DEFINED<<"\t_USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  std::cout<<COLOR_DEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_I_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_I_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  std::cout<<COLOR_DEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_II_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  std::cout<<COLOR_DEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_III_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_INVARIANT_SHEAR_VFIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_S2_
  std::cout<<COLOR_DEFINED<<"\t_USE_S2_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_S2_"<<RESET<<endl;
#endif
#ifdef _USE_S3_
  std::cout<<COLOR_DEFINED<<"\t_USE_S3_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_S3_"<<RESET<<endl;
#endif
#ifdef _USE_S2DELTA_
  std::cout<<COLOR_DEFINED<<"\t_USE_S2DELTA_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_S2DELTA_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA2_
  std::cout<<COLOR_DEFINED<<"\t_USE_DELTA2_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_DELTA2_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA3_
  std::cout<<COLOR_DEFINED<<"\t_USE_DELTA3_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_DELTA3_"<<RESET<<endl;
#endif
#ifdef _USE_NABLA2DELTA_
  std::cout<<COLOR_DEFINED<<"\t_USE_NABLA2DELTA_"<<BOLDGREEN<<"\tdefined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_NABLA2DELTA_"<<RESET<<endl;
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("NUMBER COUNTS");
#ifdef  _USE_TWO_REFS_MOCKS_
  std::cout<<COLOR_DEFINED<<"\t _USE_TWO_REFS_MOCKS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t _USE_TWO_REFS_MOCKS_"<<RESET<<endl;
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("CATALOG");
#ifdef _GET_BiasMT_CAT_
  std::cout<<COLOR_DEFINED<<"\t_GET_BiasMT_CAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_GET_BiasMT_CAT_"<<RESET<<endl;
#endif
#ifdef _READ_REF_CATALOG_
  std::cout<<COLOR_DEFINED<<"\t_READ_REF_CATALOG_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_READ_REF_CATALOG_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_TRACERS_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_TRACERS_"<<RESET<<endl;
#else
 std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_TRACERS_"<<RESET<<endl;
#endif
#ifdef _USE_VELOCITIES_TRACERS_
  std::cout<<COLOR_DEFINED<<"\t_USE_VELOCITIES_TRACERS_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_VELOCITIES_TRACERS_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_DENSITY_
  std::cout<<COLOR_DEFINED<<"\t_CORRECT_VEL_DENSITY_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_CORRECT_VEL_DENSITY_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_HALOMASS_
  std::cout<<COLOR_DEFINED<<"\t_CORRECT_VEL_HALOMASS_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_CORRECT_VEL_HALOMASS_"<<RESET<<endl;
#endif
#ifdef _USE_VMAX_AS_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_VMAX_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_VMAX_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#endif


#ifdef _USE_MASS_AS_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#endif


#ifdef _SET_CAT_WITH_MASS_CUT_
  std::cout<<COLOR_DEFINED<<"\t_SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#endif
#ifdef _WRITE_BiasMT_CATALOGS_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#endif
#ifdef _WRITE_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#endif
 std::cout<<endl;
message_screen("******************************");
message_screen("ASSIGNMENT OF HALO PROPERTIES");
#ifdef _ONLY_POST_PROC_
  std::cout<<COLOR_DEFINED<<"\t_ONLY_POST_PROC_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ONLY_POST_PROC_"<<RESET<<endl;
#endif
#ifdef  _USE_TWO_REFS_MOCKS_ASSIGNMENT_
  std::cout<<COLOR_DEFINED<<"\t_USE_TWO_REFS_MOCKS_ASSIGNMENT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_TWO_REFS_MOCKS_ASSIGNMENT_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_TO_CALIBRATION_
  std::cout<<COLOR_DEFINED<<"\t_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_PROPERTY_
  std::cout<<COLOR_DEFINED<<"\t_ASSIGN_PROPERTY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ASSIGN_PROPERTY_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_MASS_POST_
  std::cout<<COLOR_DEFINED<<"\t_ASSIGN_MASS_POST_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ASSIGN_MASS_POST_"<<RESET<<endl;
#endif
#ifdef _ONLY_POST_PROC_
#ifdef _DO_NOT_CONVOLVE_
  std::cout<<COLOR_DEFINED<<"\t_DO_NOT_CONVOLVE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_DO_NOT_CONVOLVE_"<<RESET<<endl;
#endif
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
  std::cout<<COLOR_DEFINED<<"\t_ASSIGN_PROPERTIES_TO_REFERENCE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ASSIGN_PROPERTIES_TO_REFERENCE_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
  std::cout<<COLOR_DEFINED<<"\t_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_"<<RESET<<endl;
#endif
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  std::cout<<COLOR_DEFINED<<"\t_USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#endif

#ifdef _USE_STATISTICAL_ASSIGNMENT_
  std::cout<<COLOR_DEFINED<<"\t_USE_STATISTICAL_ASSIGNMENT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_STATISTICAL_ASSIGNMENT_"<<RESET<<endl;
#endif
#ifdef _MULTISCALE_
  std::cout<<COLOR_DEFINED<<"\t_MULTISCALE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_MULTISCALE_"<<RESET<<endl;
#endif
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
  std::cout<<COLOR_DEFINED<<"\t_USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_"<<RESET<<endl;
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<"\t_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#endif
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  std::cout<<COLOR_DEFINED<<"\t_USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  std::cout<<COLOR_DEFINED<<"\t_USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#endif
#ifdef _USE_LOCAL_CLUSTERING_
  std::cout<<COLOR_DEFINED<<"\t_USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#endif
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
  std::cout<<COLOR_DEFINED<<"\t_USE_CROSS_CORRELATION_CONF_SPACE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_CROSS_CORRELATION_CONF_SPACE_"<<RESET<<endl;
#endif

#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<"\t_USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif

#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<"\t_USE_MEAN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MEAN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif

#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<"\t_USE_STDV_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_STDV_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  std::cout<<COLOR_DEFINED<<"\t_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_ "<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_"<<RESET<<endl;
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _TOP_RANDOM_
  std::cout<<COLOR_DEFINED<<"\t_TOP_RANDOM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_TOP_RANDOM_"<<RESET<<endl;
#endif
#ifdef _BOTTOM_RANDOM_
  std::cout<<COLOR_DEFINED<<"\t_BOTTOM_RANDOM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_BOTTOM_RANDOM_"<<RESET<<endl;
#endif
#endif
#ifdef _COLLAPSE_RANDOMS_
  std::cout<<COLOR_DEFINED<<"\t_COLLAPSE_RANDOMS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_COLLAPSE_RANDOMS_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_AUX_
  std::cout<<COLOR_DEFINED<<"\t_COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_VELS_
  std::cout<<COLOR_DEFINED<<"\t_COLLAPSE_RANDOMS_VELS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_COLLAPSE_RANDOMS_VELS_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
  std::cout<<COLOR_DEFINED<<"\t_COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#endif
#ifdef _APPLY_GLOBAL_EXCLUSION_
  std::cout<<COLOR_DEFINED<<"\t_APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_DENSITY_
  std::cout<<COLOR_DEFINED<<"\t_CORRECT_VEL_DENSITY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_CORRECT_VEL_DENSITY_"<<RESET<<endl;
#endif
#ifdef _USE_NP_VEL_
  std::cout<<COLOR_DEFINED<<"\t_USE_NP_VEL_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_NP_VEL_"<<RESET<<endl;
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<"\t_USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_MACH_NUMBER_
  std::cout<<COLOR_DEFINED<<"\t_USE_MACH_NUMBER_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t _USE_MACH_NUMBER_"<<RESET<<endl;
#endif
#ifdef _USE_LOCAL_OVERDENSITY_
  std::cout<<COLOR_DEFINED<<"\t_USE_LOCAL_OVERDENSITY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_LOCAL_OVERDENSITY_"<<RESET<<endl;
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
  std::cout<<COLOR_DEFINED<<"\t_USE_BIAS_OBJECT_TO_OBJECT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_BIAS_OBJECT_TO_OBJECT_"<<RESET<<endl;
#endif

  std::cout<<endl;
  message_screen("******************************");
  message_screen("POWER SPECTRUM");
#ifdef _POWER_
  std::cout<<COLOR_DEFINED<<"\t_POWER_"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<"\t        This directive is devoted for -p execution (i.e, measure powr spectgrum)_"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<"\t        For other options (e.g., -b), please undefine it."<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_POWER_"<<RESET<<RESET<<endl;
#endif

#ifdef _USE_ALL_PK_
  std::cout<<COLOR_DEFINED<<"\t_USE_ALL_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_ALL_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_CUTS_PK_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_CUTS_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_CUTS_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_BINS_PK_
  std::cout<<COLOR_DEFINED<<"\t_USE_MASS_BINS_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_MASS_BINS_PK_"<<RESET<<endl;
#endif

#ifdef _REDSHIFT_SPACE_
  std::cout<<COLOR_DEFINED<<"\t_REDSHIFT_SPACE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_REDSHIFT_SPACE_"<<RESET<<endl;
#endif

#ifdef _WRITE_MULTIPOLES_
  std::cout<<COLOR_DEFINED<<"\t_WRITE_MULTIPOLES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_WRITE_MULTIPOLES_"<<RESET<<endl;
#endif
#ifdef _USE_SEVERAL_RANDOM_FILES_
  std::cout<<COLOR_DEFINED<<"\t_USE_SEVERAL_RANDOM_FILES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_USE_SEVERAL_RANDOM_FILES_"<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OTHER ACTIONS");
#ifdef _BIN_ACCUMULATE_
  std::cout<<COLOR_DEFINED<<"\t_BIN_ACCUMULATE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_BIN_ACCUMULATE_"<<RESET<<endl;
#endif



#ifdef _DYNAMICAL_SAMPLING_
  std::cout<<COLOR_DEFINED<<"\t_DYNAMICAL_SAMPLING_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"\t_DYNAMICAL_SAMPLING_"<<RESET<<endl;
#endif



  message_screen("******************************");
   message_screen("******************************");
   message_screen("******************************");
}
