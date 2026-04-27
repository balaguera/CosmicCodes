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
# include "../include/ScreenOutput.hpp"
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
 std::cout<<" "<<YELLOW<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP (+FFTW) ESTIMATOR            *"<<endl;
 std::cout<<" VERSION 1.3                                                                        *"<<endl;
 std::cout<<" For more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" Starting time and date "<<ctime (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_c(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<" "<<YELLOW<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ANALYSIS OF DARK MATTER TRACERS - ASSEMBLY BIAS                                    *"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" Starting time and date "<<ctime (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_cl(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<" "<<YELLOW<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ANGULAR POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                *"<<endl;
 std::cout<<" USING HEALPIX DECOMPOSITION AND PEEBLES ESTIMATOR                                  * "<<endl;
 std::cout<<" VERSION 1.0                                                                        *"<<endl ;
 std::cout<<" For more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" Starting time and date"<<ctime (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_fb(){
 time_t rawtime;
 time (&rawtime);
 std::cout<<" "<<YELLOW<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" 3D POWER SPECTRUM OF COSMOLOGICAL MASS TRACERS                                     *"<<endl;
 std::cout<<" USING THE FOURIER BESSEL DECOMPOSITION                                             * "<<endl;
 std::cout<<" VERSION 1.0                                                                        *"<<endl ;
 std::cout<<" For more documentation see README/README.pdf                                       *"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" Starting time "<<ctime (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_yama(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<" "<<YELLOW<<endl;
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<" MULTIPOLE DECOMPOSITION OF THE THREE DIMENSIONAL POWER SPECTRUM OF COSMOLOGICAL   *"<<endl;
  std::cout<<" MASS TRACERS USING YAMAMOTO ESTIMATOR (FFTW-based)                                *"<<endl;
  std::cout<<" VERSION 1.1                                                                       *"<<endl;
  std::cout<<" For more documentation see README/README.pdf                                      *"<<endl;
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<" Starting time "<<ctime (&rawtime);
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<RESET<<endl;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_bispectrum(){
 std::cout<<" "<<YELLOW<<endl;
 time_t rawtime;
 time (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" BISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl;
 std::cout<<" VERSION 1.1                                                                        *"<<endl;
 std::cout<<" Documentation and warnings in readme.ps                                            *"<<endl;
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<" Starting time "<<ctime (&rawtime);
 std::cout<<" ************************************************************************************"<<endl;
 std::cout<<RESET;

}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::welcome_message_bispectrum_fast()
{
  time_t rawtime;
  time (&rawtime);
  std::cout<<" "<<YELLOW<<endl;
  std::cout<<" ************************************************************************************"<<endl;
  std::cout<<" ************************************************************************************"<<endl;
  std::cout<<" BISPECTRUM OF COSMOLOGICAL MASS TRACERS USING FKP ESTIMATOR                        *"<<endl;
  std::cout<<" AND INVERSE FFTW TRICK                                                             *"<<endl;
  std::cout<<" VERSION 1.0                                                                        *"<<endl;
  std::cout<<" Documentation and warnings in readme.ps                                            *"<<endl;
  std::cout<<" ************************************************************************************"<<endl;
  std::cout<<" Starting time "<<ctime (&rawtime);
  std::cout<<" ************************************************************************************"<<endl;
  std::cout<<RESET;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_interrupt()
{
  std::cout<<" "<<YELLOW<<endl;
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<" Code interrupted at"<<endl;
  time_t rawtime;
  time ( &rawtime );
  std::cout<<ctime (&rawtime);
  std::cout<<" ***********************************************************************************"<<endl;
  std::cout<<" "<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::enter(string text)
{
#ifdef _FULL_VERBOSE_ENTER_
//  string time_now=ctime(&this->initial_time);
 // std::cout<<BOLDGREEN<<" "<<time_now<<" =============> Going to:"<<GREEN<<text<<RESET<<endl;
  std::cout<<BOLDGREEN<<" =============> Going to:"<<GREEN<<text<<RESET<<endl;
 // time(&this->initial_time);
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::leaving(string text)
{
#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<" ==================================> leaving"<<GREEN<<text<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::ending_message()
{
  std::cout<<YELLOW<<endl;
  time_t rawtime;
  time ( &rawtime );
  std::cout<<" Ending date:"<<ctime (&rawtime)<<endl;
  std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::write_cosmo_parameters(void *p, void *pa)
{
  std::cout<<endl;  
  this->message_screen(" Cosmological Parameters");
#ifdef _USE_UNITSIM_COSMOLOGY_
  this->message_screen(" Using cosmological parameters from UNITsim (Planck 16)");
#endif
#ifdef _USE_SLICS_COSMOLOGY_
  this->message_screen(" Using cosmological parameters from SLICS (Planck 16)");
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
  std::cout<<BOLDCYAN<<" Input Cosmological parameters"<<RESET<<endl;
  this->message_screen(" Omega matter =", s_cp->Om_matter);
  this->message_screen(" Omega vac =", s_cp->Om_vac);
  this->message_screen(" Omega baryons =", s_cp->Om_baryons);
  this->message_screen(" Omega curv =", s_cp->Om_k);
  this->message_screen(" Hubble par =", s_cp->Hubble);
  this->message_screen(" Spectral index =", s_cp->spectral_index);
  this->message_screen(" Sigma8 =", s_cp->sigma8);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::error_ncolumns(string fname){
  std::cout<<RED;
  std::cout<<fname<<"  with less than 4 columns"<<endl;
  char yn;
  ScreenOutput So;
  std::cout<<" Continue (y/n)?"<<endl; cin>>yn;
  if(yn !='y' && yn !='n')
    {
      std::cout<<" Please answer  yes (y) or not (n)"<<endl;
      std::cout<<" Continue ?"<<endl; cin>>yn;
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
  if(lapse<=60)std::cout<<" \r "<<fraction<<"  % completed. Time elapsed: "<<lapse<<" secs \r";std::cout.flush();
  if(lapse>60)std::cout<<" \r "<<fraction<<"  % completed. Time elapsed: "<<lapse/60.<<" mins \r";std::cout.flush();
  if(lapse>3600)std::cout<<" \r "<<fraction<<"  % completed. Time elapsed: "<<lapse/3600.<<" hrs \r";std::cout.flush();
  if (fraction==25.0) std::cout <<" ..25%  completed"<<endl;
  if (fraction==50.0) std::cout <<" ..50%  completed"<<endl;
  if (fraction==75.0) std::cout <<" ..75%  completed"<<endl;
  if (fraction==100.0)std::cout <<" ..100% completed"<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss)
{
  std::cout<<" "<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<" "<<GREEN<<ss<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss,string sa)
{
  std::cout<<" "<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<GREEN<<" "<<ss<<RESET<<" "<<sa<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning(string ss, int line)
{
  std::cout<<" "<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<GREEN<<" "<<ss<<" "<<line<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning_ini(int line, string fun,string file, string war)
{
  std::cout<<" "<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<RED<<" from line "<<line<<"  in file "<<file<<" , function "<<fun<<RESET<<std::endl;
  std::cout<<GREEN<<war<<RESET<<std::endl;
  std::cout<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_warning_ini(int line, string fun,string file, string war, int ii)
{
  std::cout<<" "<<GREEN<<"WARNING"<<RESET<<std::endl;
  std::cout<<RED<<" Warning from line "<<line<<" in function "<<fun<<" from file "<<file<<RESET<<std::endl;
  std::cout<<GREEN<<" "<<war<<" "<<ii<<RESET<<std::endl;
  std::cout<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error(string ss)
{
  std::cout<<" "<<CYAN<<"ERROR"<<RESET<<std::endl;
  std::cerr<<CYAN<<ss<<RESET<<std::endl;
  std::cerr<<RED<<" ...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error_file(string ss)
{
  std::cerr<<CYAN<<" File "<<ss<<"  does not exist."<<RESET<<std::endl;
  std::cerr<<GREEN<<" Please check parameter file and paths."<<endl;
  std::cerr<<GREEN<<" Cosmicatlass ends here."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput:: message_error(string ss, double a)
{
  std::cerr<<CYAN<<ss<<" "<<a<<RESET<<std::endl;
  std::cerr<<RED<<" ...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_error(string ss, double a, string sa)
{
  std::cerr<<CYAN<<ss<<" "<<a<<" "<<sa<<RESET<<std::endl;
  std::cerr<<RED<<" ...Exiting..."<<endl;
  exit (1);
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss)
{
  std::clog<<" "<<YELLOW<<ss<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////

void ScreenOutput::message_screen_flush(string ss, real_prec s2, string sa, real_prec s3)
{
  std::cout<<" "<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<" "<<sa<<" "<<s3<<RESET; std::cout.flush();
}
////////////////////////////////////////////////////////////////////////////

void ScreenOutput::message_screen_flush(string ss, real_prec s2, string sa, real_prec s3, string s4, real_prec s5)
{
  std::cout<<" "<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<" "<<sa<<" "<<s3<<" "<<s4<<" "<<s5<<RESET; std::cout.flush();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen_flush(string ss, int s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<" "<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<RESET;std::cout.flush();
#else
    std::cout<<RESET<<ss<<" "<<s2<<RESET<<std::endl;
#endif
}
void ScreenOutput::message_screen_flush(string ss, ULONG s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<" "<<"\r"<<YELLOW<<ss<<RESET<<" "<<s2<<RESET;std::cout.flush();
#else
    std::cout<<RESET<<ss<<" "<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen_flush(string ss, real_prec s2)
{
#ifdef _FULL_VERBOSE_
  std::cout<<" \r"<<YELLOW<<ss<<RESET<<" "<<s2<<RESET;
  std::cout.flush();
#else
    std::cout<<RESET<<ss<<" "<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double s2)
{
#ifdef _FULL_VERBOSE_
    std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
#else
    std::cout<<" "<<RESET<<ss<<" "<<s2<<RESET<<std::endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa,double s2)
{
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<i<<" "<<sa<<"   "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa, ULONG s2){
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<i<<" "<<sa<<"   "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int i, string sa, ULONG s2,string s3){
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<i<<" "<<sa<<" "<<s2<<" "<<s3<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, int s2)
{
  std::clog<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, ULONG s2)
{
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, string s2){
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double a, string s2, string f){
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<a<<" "<<s2<<" "<<f<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, string s2)
{
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<d2<<" "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, ULONG d2, string s2)
{
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<d2<<" "<<s2<<RESET<<std::endl;
}
void ScreenOutput::message_screen(string ss, int d2, string s2)
{
  std::cout<<" "<<YELLOW<<ss<<RESET<<" "<<CYAN<<d2<<" "<<s2<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, time_t time)
{
  double lapse=difftime(time,this->initial_time);
  if(lapse<60)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse<<" secs)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse/60<<" mins)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<RED<<ss<<" "<<d2<<" (T_init + "<<lapse/3600<<" hrs)"<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, double d2, time_t time, time_t time2)
{
  double lapse=difftime(time,this->initial_time);
  double lapse2=difftime(time,time2);

  if(lapse<60)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init +"<<lapse<<" secs, time spent in last step "<<lapse2<<"  s)"<<RESET<<std::endl;
  else if(lapse>60)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init +"<<lapse/60<<" mins, time spent in last step "<<lapse2<<"  s)"<<RESET<<std::endl;
  else if(lapse>3600)
    std::clog<<BOLDCYAN<<ss<<" "<<d2<<" (T_init +"<<lapse/3600<<" hrs, time spent in last step "<<lapse2<<"  s)"<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_screen(string ss, string s2, string s3, double d3)
{
  std::clog<<GREEN<<ss<<CYAN<<" "<<s2<<" "<<BLUE<<s3<<" "<<CYAN<<d3<<RESET<<std::endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time(time_t start_all)
{
  time_t end;
  time(&end);
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
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
  std::cout<<__DATE__<<" "<<__TIME__<<endl;
  std::cout<<GREEN<<" Computing Time: "<<lapse<<" "<<units_time<<RESET<<endl;
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
  std::cout<<CYAN<<" **************************************************************************"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time_mock(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  std::cout<<" "<<YELLOW<<" Mock DF generated in : "<<lapse<<" secs"<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_time2(time_t start_all)
{
  time_t end;
  time(&end);
  double lapse=difftime(end,start_all);
  std::cout<<" "<<YELLOW<<"Computing Time: "<<CYAN<<lapse<<" secs"<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::usage(string s)
{
  std::cout<<BOLDGREEN;  
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" COSMICATLASS"<<endl;
  std::cout<<" \n CosmologicalCATalogs for LArge Scale Structure"<<endl;
  std::cout<<" \n How to run: "<<s<<"  [-option] [argument]"<<endl;
  std::cout<<" \n Options: "<<endl;
  std::cout<<"          -a for information on the author, no argument"<<endl;
  std::cout<<"          -b parameter_file.ini: Runs BiasMT"<<endl;
  std::cout<<"          -c parameter_file.ini: Analyzes tracer catalog"<<endl;
  std::cout<<"          -d parameter_file.ini: Shows preprocessor directives"<<endl;
  std::cout<<"          -f parameter_file.ini: Reads binarys form IC and cmpute density field IC"<<endl;
  std::cout<<"          -g parameter_file.ini: Generates galaxy catalog from halo catalog using HOD"<<endl;
  std::cout<<"          -h parameter_file.ini: Help"<<endl;
  std::cout<<"          -i parameter_file.ini: Shows input pars"<<endl;
  std::cout<<"          -m parameter_file.ini: Measures marked power spectrum"<<endl;
  std::cout<<"          -p parameter_file.ini: Measures Statistic"<<endl;
  std::cout<<"          -q parameter_file.ini: SuperClusters"<<endl;
  std::cout<<"          -s parameter_file.ini: Applies a low-pass filter to input density field"<<endl;
  std::cout<<"          -t parameter_file.ini: Runs LPT"<<endl;
  std::cout<<"          -u parameter_file.ini: Analyzes input halo catalog and measurements of secondary bias"<<endl;
  std::cout<<"          -w parameter_file.ini: Computes window functions. Under construction"<<endl;
  std::cout<<"          -v parameter_file.ini: Warnings "<<endl;
  std::cout<<"          -z parameter_file.ini: Assign,ment of individual tracer bias"<<endl;
  std::cout<<" Consult ../Headers/def.h for pre-procesor directives"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<RESET<<endl;
  std::cout<<RESET<<endl;                                       
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::author(){
  std::cout<<BOLDGREEN;  
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<BOLDCYAN;
  std::cout<<" *   ****  **  |*** *      *  *   ****    *   ***** *       *    |***  *"<<endl;
  std::cout<<" *  |     *  *  *   * *  * *  *  *       * *    *   *      * *    *    *"<<endl;
  std::cout<<" *  |     *  *   *| *  **  *  *  *      *****   *   *     *****    *|  *"<<endl;
  std::cout<<" *   ****  **  **** *      *  *   **** *     *  *   **** *     * ****  *"<<RESET<<endl;
  std::cout<<BOLDGREEN;  
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" COSMICATLASS                                                          *"<<endl;
  std::cout<<" Andres Balaguera-Antolinez (abalant@gmail.com)                         *"<<endl;
  std::cout<<" 2011-2023                                                             *"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<" ***********************************************************************"<<endl;
  std::cout<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
void  ScreenOutput::message(time_t start_all)
{
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN;
  std::cout<<" *********************************************************"<<endl;
  std::cout<<" *********************************************************"<<endl;
  std::cout<<" *********************************************************"<<endl;
  std::cout<<" COSMICATLASS"<<endl;
  std::cout<<" Cosmological CATalogs for LArge Scale Structure"<<endl;
  std::cout<<" IAC 2017-2023"<<endl;
  std::cout<<" *********************************************************"<<endl;
  std::cout<<" Compilation began on "<<__DATE__<<","<<__TIME__<<endl;
#ifdef _USE_OMP_
  std::cout<<" Using OMP with "<<_NTHREADS_<<" threads"<<endl;
#endif
  std::cout<<" *********************************************************"<<endl;
  std::cout<<RESET<<endl;
#else
  std::cout<<" *********************************************************"<<endl;
  std::cout<<" *********************************************************"<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  log<<" *****************************************************************"<<endl;
  log<<" *****************************************************************"<<endl;
  log<<" BiasMT. Launched at "<<__DATE__<<"  at "<<__TIME__<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  log<<" BiasMT. Realization "<<this->params._realization()<<endl;
#endif
  log<<" *****************************************************************"<<endl;
  log<<" *****************************************************************"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void  ScreenOutput::message_BiasMT(time_t start_all)
{
#ifdef _FULL_VERBOSE_
  std::cout<<CYAN;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" COSMICATLASS"<<endl;
  std::cout<<" Cosmological CATalogs for LArge Scale Structure"<<endl;
  std::cout<<" IAC 2017-2023"<<endl;
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" Compilation began on "<<__DATE__<<","<<__TIME__<<endl;
#ifdef _USE_OMP_
  std::cout<<" Using OMP with "<<_NTHREADS_<<" threads"<<endl;
#endif
  std::cout<<" *****************************************************************"<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  std::cout<<" BiasMT. Realization "<<this->params._realization()<<"  p-index (Unitsim)"<<endl;
#endif
  std::cout<<RESET<<endl;
#else
  std::cout<<" *****************************************************************"<<endl;
  std::cout<<" *****************************************************************"<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  std::cout<<" BiasMT. Realization "<<this->params._realization()<<"  p-index (Unitsim) "<<endl;
#endif
  std::cout<<" *****************************************************************"<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  log<<" *****************************************************************"<<endl;
  log<<" *****************************************************************"<<endl;
  log<<" BiasMT. Launched at "<<__DATE__<<"  at "<<__TIME__<<endl;
#if (defined _GET_BiasMT_REALIZATIONS_ || defined _USE_LPT_ ) && defined _UNITSIM_
  log<<" BiasMT. Realization "<<this->params._realization()<<"  p-index (Unitsim) "<<endl;
#endif
  log<<" *****************************************************************"<<endl;
  log<<" *****************************************************************"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_output_file(string s, ULONG l)
{
#ifdef _FULL_VERBOSE_
  std::cout<<YELLOW<<" Writting output file "<<CYAN<<s<<" with "<<l<<" lines"<<RESET<<endl;
#endif
  ofstream log; log.open(this->logfile.c_str(), ios_base::app);
  if(log.is_open())
  {
     log<<" Writting output file "<<s<<" with "<<l<<" lines"<<endl;
     log.close();
  }
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::message_output_file(string s, int nlines, int ncols)
{
#ifdef _FULL_VERBOSE_
  std::cout<<" "<<YELLOW<<"Writting output file "<<CYAN<<s<<" with "<<nlines<<" lines and "<<ncols<<" columns"<<RESET<<endl;
#endif

  ofstream log; log.open(this->logfile.c_str(),ios_base::app);
  log<<" Writting output file "<<s<<" with "<<nlines<<" lines and "<<ncols<<" columns"<<endl;
  log.close();
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::DONE()
{
#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<"      ["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
void ScreenOutput::show_preproc()
{
  message_screen("******************************");
  message_screen("******************************");
  message_screen("******************************");
  std::cout<<WHITE<<"\n PRE-PROCESSOR DIRECTIVES IN BiasMT"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<" DEFINED"<<RESET;
  std::cout<<COLOR_UNDEFINED<<"   UNDEFINED"<<RESET<<endl;
  message_screen("******************************");
  message_screen("PRECISION");
#ifdef SINGLE_PREC
  std::cout<<COLOR_DEFINED<<" SINGLE_PREC"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" SINGLE_PREC"<<RESET<<RESET<<endl;
#endif
#ifdef DOUBLE_PREC
  std::cout<<COLOR_DEFINED<<" DOUBLE_PREC"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" DOUBLE_PREC"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OMP");
#ifdef _USE_OMP_
  std::cout<<COLOR_DEFINED<<" _USE_OMP_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_OMP_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("GNUPLOT");
#ifdef _USE_GNUPLOT_
  std::cout<<COLOR_DEFINED<<" _USE_GNUPLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_GNUPLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_
  std::cout<<COLOR_DEFINED<<" _SHOW_BIAS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _SHOW_BIAS_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
  std::cout<<COLOR_DEFINED<<" _USE_GNUPLOT_ABUNDANCE_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_GNUPLOT_ABUNDANCE_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_V_PLOT_
  std::cout<<COLOR_DEFINED<<" _USE_GNUPLOT_ABUNDANCE_V_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_GNUPLOT_ABUNDANCE_V_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_RS_PLOT_
  std::cout<<COLOR_DEFINED<<" _USE_GNUPLOT_ABUNDANCE_RS_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_GNUPLOT_ABUNDANCE_RS_PLOT_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_
  std::cout<<COLOR_DEFINED<<" _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("VERBOSE");
#ifdef _USE_COLORS_
  std::cout<<COLOR_DEFINED<<" _USE_COLORS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_COLORS_"<<RESET<<RESET<<endl;
#endif
#ifdef _FULL_VERBOSE_
  std::cout<<COLOR_DEFINED<<" _FULL_VERBOSE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _FULL_VERBOSE_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("IO FORMATS");
#ifdef _WRITE_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<" _WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _READ_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<" _READ_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _READ_BINARY_BiasMT_FORMAT_"<<RESET<<RESET<<endl;
#endif
#ifdef _OUTPUT_WITH_HEADERS_
  std::cout<<COLOR_DEFINED<<" _OUTPUT_WITH_HEADERS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _OUTPUT_WITH_HEADERS_"<<RESET<<RESET<<endl;
#endif
#ifdef _WRITE_COORDINATES_
  std::cout<<COLOR_DEFINED<<" _WRITE_COORDINATES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_COORDINATES_"<<RESET<<RESET<<endl;
#endif
#ifdef _WRITE_VELOCITIES_
  std::cout<<COLOR_DEFINED<<" _WRITE_VELOCITIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_VELOCITIES_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("REFERENCE SIMULATION");
#ifdef _ABACUS_
  std::cout<<COLOR_DEFINED<<" _ABACUS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ABACUS_"<<RESET<<RESET<<endl;
#endif
#ifdef _MINERVA_
  std::cout<<COLOR_DEFINED<<" _MINERVA_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _MINERVA_"<<RESET<<RESET<<endl;
#endif
#ifdef _SLICS_
  std::cout<<COLOR_DEFINED<<" _SLICS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _SLICS_"<<RESET<<RESET<<endl;
#endif
#ifdef _UNITSIM_
  std::cout<<COLOR_DEFINED<<" _UNITSIM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _UNITSIM_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("COSMOLOGICAL PARAMETERS");
#ifdef _USE_COSMO_PARS_
  std::cout<<COLOR_DEFINED<<" _USE_COSMO_PARS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_COSMO_PARS_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_UNITSIM_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<" _USE_UNITSIM_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_UNITSIM_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_PLANCK_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<" _USE_PLANCK_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_PLANCK_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
#ifdef _USE_SLICS_COSMOLOGY_
  std::cout<<COLOR_DEFINED<<" _USE_SLICS_COSMOLOGY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_SLICS_COSMOLOGY_"<<RESET<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OPERATING MODE");
#ifdef BIAS_MODE
  std::cout<<COLOR_DEFINED<<" BIAS_MODE"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" BIAS_MODE"<<RESET<<endl;
#endif
#ifdef mode_p
  std::cout<<COLOR_DEFINED<<" mode_p"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" mode_p"<<RESET<<endl;
#endif
#ifdef mode_b
  std::cout<<COLOR_DEFINED<<" mode_b"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" mode_b"<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("BiasMT");
/*
#ifdef MOCK_MODE
#ifdef _DO_BiasMT_CALIBRATION_
  std::cout<<COLOR_DEFINED<<" _DO_BiasMT_CALIBRATION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _DO_BiasMT_CALIBRATION_"<<RESET<<endl;
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
  std::cout<<COLOR_DEFINED<<" _GET_BiasMT_REALIZATIONS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _GET_BiasMT_REALIZATIONS_"<<RESET<<endl;
#endif
#ifdef _USE_VELOCITIES_
  std::cout<<COLOR_DEFINED<<" _USE_VELOCITIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_VELOCITIES_"<<RESET<<endl;
#endif
#endif
*/
  std::cout<<endl;
message_screen("******************************");
message_screen("LPT");
#ifdef _USE_LPT_
  std::cout<<COLOR_DEFINED<<" _USE_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_LPT_"<<RESET<<endl;
#endif
#ifdef _ONLY_LPT_
  std::cout<<COLOR_DEFINED<<" _ONLY_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_ONLY_LPT_"<<RESET<<endl;
#endif
#ifdef _POWER_BOOST_ALPT_
  std::cout<<COLOR_DEFINED<<" _POWER_BOSST_ALPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _POWER_BOOST_ALPT_"<<RESET<<endl;
#endif
#ifdef _GET_VELOCITIES_LPT_
  std::cout<<COLOR_DEFINED<<" _GET_VELOCITIES_LPT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _GET_VELOCITIES_LPT_"<<RESET<<endl;
#endif
#ifdef _GET_DENSITY_FIELDS_
  std::cout<<COLOR_DEFINED<<" _GET_DENSITY_FIELDS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _GET_DENSITY_FIELDS_"<<RESET<<endl;
#endif

std::cout<<endl;
message_screen("******************************");
message_screen("DISPLACEMENTS");

#ifdef _DISPLACEMENTS_
  std::cout<<COLOR_DEFINED<<" _DISPLACEMENTS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _DISPLACEMENTS_"<<RESET<<endl;
#endif
std::cout<<endl;
  message_screen("******************************");
message_screen("CALIBRATION");

#ifdef _SMOOTHED_KERNEL_
  std::cout<<COLOR_DEFINED<<" _SMOOTHED_KERNEL_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _SMOOTHED_KERNEL_"<<RESET<<endl;
#endif
#ifdef _MODIFY_LIMITS_
  std::cout<<COLOR_DEFINED<<" _MODIFY_LIMITS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _MODIFY_LIMITS_"<<RESET<<endl;
#endif
#ifdef _RANK_ORDERING_AB_INITIO_
  std::cout<<COLOR_DEFINED<<" _RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _RANK_ORDERING_AB_INITIO_"<<RESET<<endl;
#endif
#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _RHO_WITH_DELTA_
  std::cout<<COLOR_DEFINED<<" _RO_WITH_DELTA_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _RO_WITH_DELTA_"<<RESET<<endl;
#endif
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("DARK MATTER (THETA) MODELS AND PROPERTIES");
#ifdef _USE_DM_IN_BiasMT_
  std::cout<<COLOR_DEFINED<<" _USE_DM_IN_BiasMT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_DM_IN_BiasMT_"<<RESET<<endl;
#endif

#ifdef _USE_TWEB_
  std::cout<<COLOR_DEFINED<<" _USE_TWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_TWEB_"<<RESET<<endl;
#endif

#ifdef _USE_IWEB_
  std::cout<<COLOR_DEFINED<<" _USE_IWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_IWEB_"<<RESET<<endl;
#endif

#ifdef _USE_IKWEB_
  std::cout<<COLOR_DEFINED<<" _USE_IKWEB_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_IKWEB_"<<RESET<<endl;
#endif


#ifdef _USE_CWC_
  std::cout<<COLOR_DEFINED<<" _USE_CWC_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_CWC_"<<RESET<<endl;
#endif
#ifdef _USE_CWC_INSIDE_LOOP_
  std::cout<<COLOR_DEFINED<<" _USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_CWC_INSIDE_LOOP_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_KNOTS_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_KNOTS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MASS_KNOTS_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  std::cout<<COLOR_DEFINED<<" _USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_INVARIANT_TIDAL_FIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  std::cout<<COLOR_DEFINED<<" _USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_INVARIANT_TIDAL_FIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_TIDAL_ANISOTROPY_
  std::cout<<COLOR_DEFINED<<" _USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_TIDAL_ANISOTROPY_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  std::cout<<COLOR_DEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_I_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_I_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  std::cout<<COLOR_DEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_II_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_II_"<<RESET<<endl;
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  std::cout<<COLOR_DEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_III_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_INVARIANT_SHEAR_VFIELD_III_"<<RESET<<endl;
#endif
#ifdef _USE_S2_
  std::cout<<COLOR_DEFINED<<" _USE_S2_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_S2_"<<RESET<<endl;
#endif
#ifdef _USE_S3_
  std::cout<<COLOR_DEFINED<<" _USE_S3_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_S3_"<<RESET<<endl;
#endif
#ifdef _USE_S2DELTA_
  std::cout<<COLOR_DEFINED<<" _USE_S2DELTA_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_S2DELTA_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA2_
  std::cout<<COLOR_DEFINED<<" _USE_DELTA2_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_DELTA2_"<<RESET<<endl;
#endif
#ifdef _USE_DELTA3_
  std::cout<<COLOR_DEFINED<<" _USE_DELTA3_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_DELTA3_"<<RESET<<endl;
#endif
#ifdef _USE_NABLA2DELTA_
  std::cout<<COLOR_DEFINED<<" _USE_NABLA2DELTA_"<<BOLDGREEN<<" defined"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_NABLA2DELTA_"<<RESET<<endl;
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("NUMBER COUNTS");
#ifdef  _USE_TWO_REFS_MOCKS_
  std::cout<<COLOR_DEFINED<<"  _USE_TWO_REFS_MOCKS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"  _USE_TWO_REFS_MOCKS_"<<RESET<<endl;
#endif
std::cout<<endl;
message_screen("******************************");
message_screen("CATALOG");
#ifdef _GET_BiasMT_CAT_
  std::cout<<COLOR_DEFINED<<" _GET_BiasMT_CAT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _GET_BiasMT_CAT_"<<RESET<<endl;
#endif
#ifdef _READ_REF_CATALOG_
  std::cout<<COLOR_DEFINED<<" _READ_REF_CATALOG_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _READ_REF_CATALOG_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_TRACERS_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_TRACERS_"<<RESET<<endl;
#else
 std::cout<<COLOR_UNDEFINED<<" _USE_MASS_TRACERS_"<<RESET<<endl;
#endif
#ifdef _USE_VELOCITIES_TRACERS_
  std::cout<<COLOR_DEFINED<<" _USE_VELOCITIES_TRACERS_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_VELOCITIES_TRACERS_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_DENSITY_
  std::cout<<COLOR_DEFINED<<" _CORRECT_VEL_DENSITY_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _CORRECT_VEL_DENSITY_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_HALOMASS_
  std::cout<<COLOR_DEFINED<<" _CORRECT_VEL_HALOMASS_"<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _CORRECT_VEL_HALOMASS_"<<RESET<<endl;
#endif
#ifdef _USE_VMAX_AS_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_VMAX_AS_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_VMAX_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_VMAX_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#endif


#ifdef _USE_MASS_AS_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MASS_AS_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MASS_AS_PRIMARY_OBSERVABLE_"<<RESET<<endl;
#endif


#ifdef _SET_CAT_WITH_MASS_CUT_
  std::cout<<COLOR_DEFINED<<" _SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _SET_CAT_WITH_MASS_CUT_"<<RESET<<endl;
#endif
#ifdef _WRITE_BiasMT_CATALOGS_
  std::cout<<COLOR_DEFINED<<" _WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#endif
#ifdef _WRITE_BINARY_BiasMT_FORMAT_
  std::cout<<COLOR_DEFINED<<" _WRITE_BiasMT_CATALOGS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_BINARY_BiasMT_FORMAT_"<<RESET<<endl;
#endif
 std::cout<<endl;
message_screen("******************************");
message_screen("ASSIGNMENT OF HALO PROPERTIES");
#ifdef _ONLY_POST_PROC_
  std::cout<<COLOR_DEFINED<<" _ONLY_POST_PROC_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ONLY_POST_PROC_"<<RESET<<endl;
#endif
#ifdef  _USE_TWO_REFS_MOCKS_ASSIGNMENT_
  std::cout<<COLOR_DEFINED<<" _USE_TWO_REFS_MOCKS_ASSIGNMENT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_TWO_REFS_MOCKS_ASSIGNMENT_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_TO_CALIBRATION_
  std::cout<<COLOR_DEFINED<<" _ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ASSIGN_TO_CALIBRATION_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_PROPERTY_
  std::cout<<COLOR_DEFINED<<" _ASSIGN_PROPERTY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ASSIGN_PROPERTY_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_MASS_POST_
  std::cout<<COLOR_DEFINED<<" _ASSIGN_MASS_POST_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ASSIGN_MASS_POST_"<<RESET<<endl;
#endif
#ifdef _ONLY_POST_PROC_
#ifdef _DO_NOT_CONVOLVE_
  std::cout<<COLOR_DEFINED<<" _DO_NOT_CONVOLVE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _DO_NOT_CONVOLVE_"<<RESET<<endl;
#endif
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
  std::cout<<COLOR_DEFINED<<" _ASSIGN_PROPERTIES_TO_REFERENCE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ASSIGN_PROPERTIES_TO_REFERENCE_"<<RESET<<endl;
#endif
#ifdef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
  std::cout<<COLOR_DEFINED<<" _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_"<<RESET<<endl;
#endif
#ifdef _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_
  std::cout<<COLOR_DEFINED<<" _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_PROPERTY_ASSIGNMENT_READING_REF_PROPERTIES_"<<RESET<<endl;
#endif

#ifdef _USE_STATISTICAL_ASSIGNMENT_
  std::cout<<COLOR_DEFINED<<" _USE_STATISTICAL_ASSIGNMENT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_STATISTICAL_ASSIGNMENT_"<<RESET<<endl;
#endif
#ifdef _MULTISCALE_
  std::cout<<COLOR_DEFINED<<" _MULTISCALE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _MULTISCALE_"<<RESET<<endl;
#endif
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
  std::cout<<COLOR_DEFINED<<" _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_"<<RESET<<endl;
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<" _USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_TRACERS_IN_CELLS_"<<RESET<<endl;
#endif
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  std::cout<<COLOR_DEFINED<<" _USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_NUMBER_OF_NEIGHBOURS_"<<RESET<<endl;
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  std::cout<<COLOR_DEFINED<<" _USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MIN_DISTANCE_TO_NEIGHBOURS_"<<RESET<<endl;
#endif
#ifdef _USE_LOCAL_CLUSTERING_
  std::cout<<COLOR_DEFINED<<" _USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_LOCAL_CLUSTERING_"<<RESET<<endl;
#endif
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
  std::cout<<COLOR_DEFINED<<" _USE_CROSS_CORRELATION_CONF_SPACE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_CROSS_CORRELATION_CONF_SPACE_"<<RESET<<endl;
#endif

#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<" _USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MIN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif

#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<" _USE_MEAN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MEAN_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif

#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
  std::cout<<COLOR_DEFINED<<" _USE_STDV_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_STDV_SEPARATIONS_IN_CELLS_"<<RESET<<endl;
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  std::cout<<COLOR_DEFINED<<" _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_ "<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_"<<RESET<<endl;
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _TOP_RANDOM_
  std::cout<<COLOR_DEFINED<<" _TOP_RANDOM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _TOP_RANDOM_"<<RESET<<endl;
#endif
#ifdef _BOTTOM_RANDOM_
  std::cout<<COLOR_DEFINED<<" _BOTTOM_RANDOM_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _BOTTOM_RANDOM_"<<RESET<<endl;
#endif
#endif
#ifdef _COLLAPSE_RANDOMS_
  std::cout<<COLOR_DEFINED<<" _COLLAPSE_RANDOMS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _COLLAPSE_RANDOMS_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_AUX_
  std::cout<<COLOR_DEFINED<<" _COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _COLLAPSE_RANDOMS_AUX_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_VELS_
  std::cout<<COLOR_DEFINED<<" _COLLAPSE_RANDOMS_VELS_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _COLLAPSE_RANDOMS_VELS_"<<RESET<<endl;
#endif
#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
  std::cout<<COLOR_DEFINED<<" _COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _COLLAPSE_RANDOMS_USING_EXCLUSION_"<<RESET<<endl;
#endif
#ifdef _APPLY_GLOBAL_EXCLUSION_
  std::cout<<COLOR_DEFINED<<" _APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _APPLY_GLOBAL_EXCLUSION_"<<RESET<<endl;
#endif
#ifdef _CORRECT_VEL_DENSITY_
  std::cout<<COLOR_DEFINED<<" _CORRECT_VEL_DENSITY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _CORRECT_VEL_DENSITY_"<<RESET<<endl;
#endif
#ifdef _USE_NP_VEL_
  std::cout<<COLOR_DEFINED<<" _USE_NP_VEL_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_NP_VEL_"<<RESET<<endl;
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_RS_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  std::cout<<COLOR_DEFINED<<" _USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_SPIN_AS_DERIVED_OBSERVABLE_"<<RESET<<endl;
#endif
#ifdef _USE_MACH_NUMBER_
  std::cout<<COLOR_DEFINED<<" _USE_MACH_NUMBER_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<"  _USE_MACH_NUMBER_"<<RESET<<endl;
#endif
#ifdef _USE_LOCAL_OVERDENSITY_
  std::cout<<COLOR_DEFINED<<" _USE_LOCAL_OVERDENSITY_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_LOCAL_OVERDENSITY_"<<RESET<<endl;
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
  std::cout<<COLOR_DEFINED<<" _USE_BIAS_OBJECT_TO_OBJECT_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_BIAS_OBJECT_TO_OBJECT_"<<RESET<<endl;
#endif

  std::cout<<endl;
  message_screen("******************************");
  message_screen("POWER SPECTRUM");
#ifdef _POWER_
  std::cout<<COLOR_DEFINED<<" _POWER_"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<"         This directive is devoted for -p execution (i.e, measure powr spectgrum)_"<<RESET<<endl;
  std::cout<<COLOR_DEFINED<<"         For other options (e.g., -b), please undefine it."<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _POWER_"<<RESET<<RESET<<endl;
#endif

#ifdef _USE_ALL_PK_
  std::cout<<COLOR_DEFINED<<" _USE_ALL_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_ALL_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_CUTS_PK_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_CUTS_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MASS_CUTS_PK_"<<RESET<<endl;
#endif

#ifdef _USE_MASS_BINS_PK_
  std::cout<<COLOR_DEFINED<<" _USE_MASS_BINS_PK_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_MASS_BINS_PK_"<<RESET<<endl;
#endif

#ifdef _REDSHIFT_SPACE_
  std::cout<<COLOR_DEFINED<<" _REDSHIFT_SPACE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _REDSHIFT_SPACE_"<<RESET<<endl;
#endif

#ifdef _WRITE_MULTIPOLES_
  std::cout<<COLOR_DEFINED<<" _WRITE_MULTIPOLES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _WRITE_MULTIPOLES_"<<RESET<<endl;
#endif
#ifdef _USE_SEVERAL_RANDOM_FILES_
  std::cout<<COLOR_DEFINED<<" _USE_SEVERAL_RANDOM_FILES_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _USE_SEVERAL_RANDOM_FILES_"<<RESET<<endl;
#endif
  std::cout<<endl;
  message_screen("******************************");
  message_screen("OTHER ACTIONS");
#ifdef _BIN_ACCUMULATE_
  std::cout<<COLOR_DEFINED<<" _BIN_ACCUMULATE_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _BIN_ACCUMULATE_"<<RESET<<endl;
#endif



#ifdef _DYNAMICAL_SAMPLING_
  std::cout<<COLOR_DEFINED<<" _DYNAMICAL_SAMPLING_"<<RESET<<endl;
#else
  std::cout<<COLOR_UNDEFINED<<" _DYNAMICAL_SAMPLING_"<<RESET<<endl;
#endif



  message_screen("******************************");
   message_screen("******************************");
   message_screen("******************************");
}
