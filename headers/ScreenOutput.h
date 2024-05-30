#ifndef _SCREEN_OUTPUT_
#define _SCREEN_OUTPUT_
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "CosmologicalFunctions.h"
# include "Params.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
class ScreenOutput{
 private:
 Params params;
 public:
  ScreenOutput(){
    cout.precision(8);
    cout.setf(ios::showpoint);
  //  cout.setf(ios::scientific);
  }
 ScreenOutput(time_t _init_t, string _logfile): initial_time(_init_t), logfile(_logfile){

    cout.precision(8);
    cout.setf(ios::showpoint);
//    cout.setf(ios::scientific);

 }
 ~ScreenOutput(){}
  string logfile;
  time_t initial_time;
  void show_warnings();
  void message(string);
  void welcome_message();
  void welcome_message_cl();
  void welcome_message_c();
  void welcome_message_fb();
  void welcome_message_yama();
  void welcome_message_bispectrum();
  void welcome_message_bispectrum_fast();
  void message_interrupt();
  void enter(string);  
  void leaving(string);  
  void error_ncolumns(string);
  void ending_message();
  void write_parameters_estimator(void *);
  void write_cosmo_parameters(void *, void *);
  void write_cosmo_parameters(void *);
  void write_yam_parameters(void *, bool);
  void comp_time(time_t, unsigned long, unsigned long);
  void message_screen(string ss, double d2, time_t time);
  void message_screen(string ss, double d2, time_t time, time_t time2);
  void message_screen_flush(string ss, int s2);
  void message_screen_flush(string ss, ULONG s2);
  void message_screen_flush(string ss, real_prec s2);
  void message_screen_flush(string ss, real_prec s2, string sa, real_prec s3);
  void message_screen_flush(string ss, real_prec s2, string sa, real_prec s3, string, real_prec);
  void message_warning(string ss);
  void message_warning(string ss,string sa);
  void message_warning(string ss, int);
  void message_warning_ini(int, string,string, string);
  void message_warning_ini(int, string,string, string, int);
  void message_error(string ss);
  void message_error(string ss, double a);
  void message_error_file(string ss);
  void message_error(string ss, double a, string sa);
  void message_screen(string ss);
  void message_screen(string ss, int i, string sa,double s2);
  void message_screen(string ss, int i, string sa,ULONG s2);
  void message_screen(string ss, int i, string sa,ULONG s2,string s3);
  void message_screen(string ss, int s2);
  void message_screen(string ss, ULONG s2);
  void message_screen(string ss, double s2);
  void message_screen(string ss, string s2);
  void message_screen(string ss, double d2, string s2);
  void message_screen(string ss, ULONG d2, string s2);
  void message_screen(string ss, int d2, string s2);
  void message_screen(string ss, string s2, string s3, double d3);
  void message_time2(time_t start_all);
  void message_time_mock(time_t start_all);
  void message_time(time_t start_all);
  void message_screen(string ss, double a, string s2, string f);
  void usage(string );
  void author();
  void message(time_t);
  void message_BiasMT(time_t);
  void message_output_file(string s, ULONG l);
  void message_output_file(string s, int l, int n);
  void DONE();
  void show_preproc();
  void set_params(Params par){this->params=par;}
};

#endif
