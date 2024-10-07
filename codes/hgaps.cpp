// ***********************************************************************
// ESTIMATES OF ANGULAR POWER SPECTRUM USING HARMONIC ANALYSIS           *
// VERSION 1.0                                                           *
// v1.1 USES HEALPIX 3.30                                                *
// AND GSL 1.13                                                          *                                                     *
// Andres Balaguera Antolinez                                            *
// Roma3, IAC                                                            *
// 2014-2017                                                             *
// ***********************************************************************
// ***********************************************************************
# include "../headers/AngularPowerSpectrumF.h"

// ********************************-***********************************************************************************************************************
void message(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"ESTIMATES OF ANGULAR POWER SPECTRUM                   *"<<endl;
  std::cout<<"v1.0                                                  *"<<endl;
  std::cout<<"abalant@R3, IAC 2014-2017                             *"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"Starting time and date"<<std::endl;
  std::cout<<ctime (&rawtime)<<std::endl;
  std::cout<<"*******************************************************"<<endl;
}

// *******************************************************************************************************************************************************
void message_mock(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"MOCK CATALOGS FOR ANGULAR POWER SPECTRUM              *"<<endl;
  std::cout<<"v1.0                                                  *"<<endl;
  std::cout<<"abalant@IAC 2019                                      *"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"Starting time and date"<<std::endl;
  std::cout<<ctime (&rawtime)<<std::endl;
  std::cout<<"*******************************************************"<<endl;
}

// *******************************************************************************************************************************************************

void usage(string s){
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<RED<<"  Clgal a code to measure angular power spectrum of cosmological mass tracers"<<endl;
  cout<<RED<<"  Usage "<<s<<" [-option] [argument]"<<endl;
  cout<<RED<<"  Options:     -h for more information "<<endl;
  cout<<RED<<"               -a for information on the author"<<endl;
  cout<<RED<<"               -p followed by parameter_file.ini to measure angular power spectrum  "<<endl;
  cout<<RED<<"               -m followed by parameter_file.ini to create mock catalog based on log-normal  "<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;}

void author(){
  cout<<CYAN<<"                                                                                 "<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" Copyright. 2017-2023. Andres Balaguera Antolinez"<<endl;
  cout<<" Code developed for the power spectrum analysis of the 2MPZ sample"<<endl;
  cout<<" This code is public. If you want to use it for your research,     "<<endl;
  cout<<" contact the author at balaguera@iac.es"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<CYAN<<"                                                                                 "<<endl;
}


// *******************************************************************************************************************************************************

void message_time(time_t start_all){
    time_t end;
    time(&end);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
    double lapse=difftime(end,start_all);
    std::cout<<CYAN<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
    time(&start_all);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
}




// ######################################################################
// ######################################################################
// ######################################################################
void gal2equ(double b, double l, double *ra, double *dec){
    double delta_p=CFACTOR*27.128336;                         /*Declination of galactic north pole in radians*/
    double l_n   =CFACTOR*122.932;                           /*Galactic longitude of the celestial pole, in radians*/
    double delta=sin(b*CFACTOR)*sin(delta_p)+cos(b*CFACTOR)*cos(delta_p)*cos(CFACTOR*l-l_n);
    *dec=asin(delta);
    double alpha1=sin(CFACTOR*l-l_n);
    double alpha2=cos(CFACTOR*l-l_n)*sin(delta_p)-tan(CFACTOR*b)*cos(delta_p);
    *ra=atan2(alpha1, alpha2)+12.25;
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

int main(int argc, char *argv[]){
  time_t start_all;
  char temp;
  time(&start_all);

  string par_file;
  if(argc==1){
    usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hap:k:d:")) != -1){
    if(temp=='h')
      {
  	usage(argv[0]);
      }
    else if (temp=='a')
    	author();
    else if(temp=='p')
      {
      	message();
    	  par_file = argv[2];
          ParamsCl params(par_file);
          Cl_FUNCTIONS cCl(params);
          cCl.get_angular_power();
	  message_time(start_all);
      }
    else if(temp=='k')
      {
      	message_mock();
	par_file = argv[2];
	ParamsCl params(par_file);
	Cl_FUNCTIONS cCl(params);
	cCl.get_cross_Cl();
	message_time(start_all);
      	break;
      }
    else if(temp=='d')
      {
        message_mock();
        ParamsCl params(par_file);
        Cl_FUNCTIONS cCl(params);
        cCl.get_Cl_bias();
        message_time(start_all);
        break;
      }
  }
  return 0;
}

