/* Code to read IC's from the SLIC simulation. Andrés Balaguera-Antolínez IAC 2020 */


// This option is only available if FFTW is propery working

#define USE_SUPERSAMPLING
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_histogram.h>
# include <gsl/gsl_histogram2d.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_expint.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_statistics_float.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_sort_vector.h>
# include <gsl/gsl_sort_vector_float.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_poly.h>

#include <ctime>
#include <cmath>
#include <cctype>
#include <string>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <omp.h>

#ifdef USE_SUPERSAMPLING
#include <fftw3.h>
#endif

#define gsl_real double
#define USE_OMP  // Use OMP


using namespace std;
// Variables used in FFTW
#define REAL 0
#define IMAG 1
#define ic_rank 3
// Variables used in Supersampling and FFTW
#define ss_factor static_cast<int>(2)  // Used for ssuper-sampling
#define MAX_MAS_DEG static_cast<int>(3)
#define MAX_MAS_DEG_PSC static_cast<int>(4)


// Offsets provided by Joachim
#define x_offset static_cast<real_prec>(0)
#define y_offset static_cast<real_prec>(0)
#define z_offset static_cast<real_prec>(0)

#define i_xh  0  // column label for the x-coordinate of halos
#define i_yh  1  // column label for the y-coordinate of halos
#define i_zh  2  // column label for the z-coordinate of halos
#define i_mass  16  // column label for the z-coordinate of halos

#define Min_Mass (2e11)
#define mass_factor (5.2175e8)

#define RED     "\033[31m"      /* Red */
#define RESET   "\033[0m"

#define ULONG unsigned long
//#define DOUBLE_PREC
#ifdef DOUBLE_PREC
#define fftw_real double
#define real_prec fftw_real
#define complex_prec fftw_complex
#else
#define SINGLE_PREC
#define real_prec float
#define fftwf_real real_prec
#define complex_prec fftwf_complex
#endif
//#define _GET_INTERPOLATED_FIELD_FROM_HALO_CAT_

#define _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_

#define num_1 static_cast<int>(1)
#define column_x static_cast<int>(0)
#define column_y static_cast<int>(1)
#define column_z static_cast<int>(2)
#define num_0_5 static_cast<real_prec>(0.5)
#define Ncolumns static_cast<int>(6)  // Number of properties of IC (x,y,z,vx,vy,vz)
#define Number_headers static_cast<int>(1)  // Number of properties of IC (x,y,z,vx,vy,vz)
#define Nres_sim static_cast<int>(4)  //Number  of sub-cells (per dim) in which the simulation volume has ben divided
#define Lside_sim static_cast<real_prec>(3072.0)  //Lenght of simulation volume
#define Lside static_cast<real_prec>(505.0)  //Comoving lenght (in Mpc/h) of cosmological volume
#define redshift (1.04)   // Redshift

#ifdef _GET_INTERPOLATED_FIELD_FROM_HALO_CAT_
#undef _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#endif

// **********************************************************************************************************************
ULONG get_nobjects(vector<real_prec> &density_field)
{
  real_prec ans=0;
#pragma omp parallel for reduction(+:ans)
  for(ULONG i=0;i< density_field.size();++i)
    ans+=static_cast<real_prec>(density_field[i]);

  return ans;
}
// **********************************************************************************************************************
void get_overdens(vector<real_prec>&in, vector<real_prec>&out,real_prec &mean)
{

  cout<<"Converting to Overdensity"<<endl;
  mean=0;

#pragma omp parallel for reduction(+:mean)
  for(ULONG i=0;i< in.size();++i)
    mean+=static_cast<real_prec>(in[i]);

  mean/=static_cast<double>(in.size());
//  mean/=static_cast<double>(Lside*Lside*Lside);
#pragma omp parallel for
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/mean-1.0;
  
}


// **********************************************************************************************************************
size_t m_getBuffer(istream& inStream) {
    inStream.clear();
    inStream.seekg(0, inStream.beg);
    string line = "";
    while (false == inStream.eof() && ((true == line.empty()) || ('#' == line[0]))) {
      getline(inStream, line);
    }
    size_t numLines = 3*1024*1024/line.size();
    numLines += numLines&0x1;
    inStream.clear();
    inStream.seekg(0, inStream.beg);
    return numLines;
  }
// **********************************************************************************************************************
size_t m_countLines(istream& inStream){
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  size_t counter = -1;
  string line = "";
  while (false == inStream.eof()) {
    getline(inStream, line);
    counter++;
 }
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  return counter;
}
// **********************************************************************************************************************
ULONG read_file(string fname, vector< real_prec > &prop_ob, int NTHREADS)
{
  int items_per_line = 0;
  unsigned nElem = 0;
  cout<<"Reading input file "<<fname<<endl;
  ifstream inputf (fname.c_str(), ifstream::in);
  if (!inputf) {
   std::cerr<<"Error in opening the input file " << fname << "!" << endl; exit(1);
  }
  nElem=m_countLines(inputf);
  cout<<"Number of lines = "<<static_cast<int>(nElem)<<endl;
  string firstline;
  real_prec aa;
  if(!inputf.eof()) {
    getline (inputf, firstline);
    stringstream ss(firstline);
    while(ss>>aa)
      items_per_line++;
  }
  inputf.seekg (0, inputf.beg);
  cout<<"Number of columns = "<<items_per_line<<endl;
  prop_ob.resize(nElem * items_per_line);
  size_t bufSize = min(static_cast<int>(m_getBuffer(inputf)), static_cast<int>(nElem));
  vector<string> tmplines(bufSize, "");
  int line_idx = 0;
  int ii;
  size_t doneElem = 0;
  size_t load = 0; // buffer to load
  volatile size_t nfail = 0; // bad catalog lines
  string line = ""; //line to read
#pragma omp parallel num_threads(NTHREADS)
  {
    while ( (doneElem < nElem) && (false == inputf.eof()) ){
#pragma omp single
      {
        load = min(static_cast<int>(bufSize), static_cast<int>(nElem - doneElem));
        for (size_t i=0; i<load; ++i)
          getline(inputf, tmplines[i]);
        nfail = 0;
      }
      // -- process buffer - parallel -
#pragma omp for reduction (+:nfail)
      for (size_t i=0; i<load; ++i)
        {
          // skip empty lines (none of them should be empty)
          if (true == tmplines[i].empty())
            {
              nfail ++;
              continue;
            }
          int jj = 0;
          stringstream ss(tmplines[i]);
          real_prec bb;
          while (ss>>bb) prop_ob[(doneElem+i)*items_per_line + jj++] = static_cast<real_prec>(bb);
        }
#pragma omp single
      {
        size_t effective = load - nfail;
        doneElem += effective;
      }
    }
  }
  inputf.clear();
  inputf.close();
  return static_cast<unsigned long>(nElem);
}



// **********************************************************************************************************************
inline ULONG index_3d(int i, int j, int k, int Nj, int Nk)
{
  return static_cast<ULONG>(k)+static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
// **********************************************************************************************************************
void index2coords(int N, int index, int  &XG, int &YG, int &ZG )
{// Get the index in each direction : F-order cells (N,index, X,Y,Z). -Corder (N,index, Z,Y,X)
  XG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  YG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  ZG = index;
}

// **********************************************************************************************************************
// **********************************************************************************************************************
#ifdef USE_SUPERSAMPLING

void do_fftw_r2c(int Nft, vector<real_prec>in, complex_prec *out)
{
  int *n=(int *)fftw_malloc(ic_rank*sizeof(float));
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#ifdef SINGLE_PREC
  fftwf_plan plan_r2c=fftwf_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftwf_execute(plan_r2c);
  fftwf_destroy_plan(plan_r2c);
#endif
#ifdef DOUBLE_PREC
  fftw_plan plan_r2c=fftw_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftw_execute(plan_r2c);
  fftw_destroy_plan(plan_r2c);
#endif
  fftw_free(n);
}

//##################################################################################
//##################################################################################
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out)
{
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
  int *n=(int *)malloc(ic_rank*sizeof(float));
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef SINGLE_PREC
  fftwf_plan plan_c2r=fftwf_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftwf_execute(plan_c2r);
  fftwf_destroy_plan(plan_c2r);
#endif
#ifdef DOUBLE_PREC
  fftw_plan plan_c2r=fftw_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftw_execute(plan_c2r);
  fftw_destroy_plan(plan_c2r);
#endif

  free(n);

#ifdef USE_OMP
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;i++)  // normalize the c2r transform
    out[i]/=static_cast<double>(NGRID);
}

#endif

// **********************************************************************************************************************
// **********************************************************************************************************************
// **********************************************************************************************************************
inline real_prec MAS_TSC(real_prec x)
{
  x=fabs(x);
  real_prec ans;
  if(x<0.5)
      ans= (0.75-x*x);
  if(x>=0.5 && x<1.5)
  {
    real_prec r = 1.5 - x;
    r *= r;
    ans= (0.5*r);
  }
    else
      ans= 0;
    return ans;
}

// **********************************************************************************************************************
inline real_prec MAS_PCS(real_prec x)
{
   real_prec ans;
   real_prec y=fabs(x);
   if(y>=0 && y<1)
     ans=(1./6.)*(4.0- 6.0*x*x + 3.*pow(y,3));
   else if(y>=1. && y<2.)
     ans= (1./6.)*pow(2.0-y, 3);
   else
    ans=0.;
  return ans;
}
// **********************************************************************************************************************
void getDensity_PCS(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec deltax, real_prec deltay, real_prec deltaz, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{

#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif
  real_prec rdelta_x=1./deltax;
  real_prec rdelta_y=1./deltay;
  real_prec rdelta_z=1./deltaz;
  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];
  ULONG counter_a=0;
  for(ULONG ip=0;ip<xp.size();++ip)
    {
      real_prec xpos=xp[ip];
      real_prec ypos=yp[ip];
      real_prec zpos=zp[ip];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      int xc = static_cast<ULONG>(floor((xpos-min1)*rdelta_x)); // indices of the cell of the particle
      int yc = static_cast<ULONG>(floor((ypos-min2)*rdelta_y));
      int zc = static_cast<ULONG>(floor((zpos-min3)*rdelta_z));
      xc = static_cast<ULONG>(fmod(real_prec(xc),real_prec(N1)));
      yc = static_cast<ULONG>(fmod(real_prec(yc),real_prec(N2)));
      zc = static_cast<ULONG>(fmod(real_prec(zc),real_prec(N3)));
      real_prec xx  = deltax*(static_cast<real_prec>(xc)+0.5);
      real_prec yy  = deltay*(static_cast<real_prec>(yc)+0.5);
      real_prec zz  = deltaz*(static_cast<real_prec>(zc)+0.5);
      real_prec xxf = deltax*(static_cast<real_prec>(xc)+1.5);
      real_prec yyf = deltay*(static_cast<real_prec>(yc)+1.5);
      real_prec zzf = deltaz*(static_cast<real_prec>(zc)+1.5);
      real_prec xxff = deltax*(static_cast<real_prec>(xc)+2.5);
      real_prec yyff = deltay*(static_cast<real_prec>(yc)+2.5);
      real_prec zzff = deltaz*(static_cast<real_prec>(zc)+2.5);
      real_prec xxbb = deltax*(static_cast<real_prec>(xc)-1.5);
      real_prec yybb = deltay*(static_cast<real_prec>(yc)-1.5);
      real_prec zzbb = deltaz*(static_cast<real_prec>(zc)-1.5);
      real_prec xxb = deltax*(static_cast<real_prec>(xc)-0.5);
      real_prec yyb = deltay*(static_cast<real_prec>(yc)-0.5);
      real_prec zzb = deltaz*(static_cast<real_prec>(zc)-0.5);
      int Xb=(xc==0 ? N1: xc);
      int Yb=(yc==0 ? N2: yc);
      int Zb=(zc==0 ? N3: zc);
      int Xf=(xc==N1-1 ? -1: xc);
      int Yf=(yc==N2-1 ? -1: yc);
      int Zf=(zc==N3-1 ? -1: zc);
      int Xbb=(xc==0 ? N1: xc);
      int Ybb=(yc==0 ? N2: yc);
      int Zbb=(zc==0 ? N3: zc);
      if(xc!=0)
        Xbb=(xc==1 ? N1+1: xc);
      if(yc!=0)
       Ybb=(yc==1 ? N2+1: yc);
      if(zc!=0)
       Zbb=(zc==1 ? N3+1: zc);
      int Xff=(xc==N1-1 ? -1: xc);
      int Yff=(yc==N2-1 ? -1: yc);
      int Zff=(zc==N3-1 ? -1: zc);
      if(xc!=N1-1)
         Xff=(xc==N1-2 ? -2: xc);
      if(yc!=N2-1)
        Yff=(yc==N2-2 ? -2: yc);
      if(zc!=N3-1)
        Zff=(zc==N3-2 ? -2: zc);
      vector<int> i_idx = {Xbb-2, Xb-1, xc, Xf+1, Xff+2};
      vector<int> j_idx = {Ybb-2, Yb-1, yc, Yf+1, Xff+2};
      vector<int> k_idx = {Zbb-2, Zb-1, zc, Zf+1, Zff+2};
      vector<real_prec> MAS_xx=
        {
          MAS_PCS((xxbb - xpos)*rdelta_x),
          MAS_PCS((xxb  - xpos)*rdelta_x),
          MAS_PCS((xx   - xpos)*rdelta_x),
          MAS_PCS((xxf  - xpos)*rdelta_x),
          MAS_PCS((xxff - xpos)*rdelta_x)
      };
      vector<real_prec> MAS_yy =
        {
          MAS_PCS((yybb - ypos)*rdelta_y),
          MAS_PCS((yyb  - ypos)*rdelta_y),
          MAS_PCS((yy   - ypos)*rdelta_y),
          MAS_PCS((yyf  - ypos)*rdelta_y),
          MAS_PCS((yyff - ypos)*rdelta_y)
        };
      vector<real_prec> MAS_zz =
        {
          MAS_PCS((zzbb - zpos)*rdelta_z),
          MAS_PCS((zzb  - zpos)*rdelta_z),
          MAS_PCS((zz   - zpos)*rdelta_z),
          MAS_PCS((zzf  - zpos)*rdelta_z),
          MAS_PCS((zzff - zpos)*rdelta_z)
        };
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
      for(int ih=0;ih<MAX_MAS_DEG_PSC;++ih)
        for(int jh=0;jh<MAX_MAS_DEG_PSC;++jh)
          for(int kh=0;kh<MAX_MAS_DEG_PSC;++kh)
#ifdef USE_OMP
#pragma omp atomic update
#endif
            delta[index_3d(i_idx[ih],j_idx[jh], k_idx[kh], N3, N2) ] += MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
    counter_a++;
  }
  ULONG count=0;
#pragma omp parallel for reduction(+:count)
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// **********************************************************************************************************************
void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec deltax, real_prec deltay, real_prec deltaz, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{

#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif

  real_prec rdelta_x=1./deltax;
  real_prec rdelta_y=1./deltay;
  real_prec rdelta_z=1./deltaz;
  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];

  ULONG counter_a=0;
  for(ULONG ip=0;ip<xp.size();++ip)
    {
      real_prec xpos=xp[ip];
      real_prec ypos=yp[ip];
      real_prec zpos=zp[ip];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      int xc = static_cast<ULONG>(floor((xpos-min1)*rdelta_x)); // indices of the cell of the particle
      int yc = static_cast<ULONG>(floor((ypos-min2)*rdelta_y));
      int zc = static_cast<ULONG>(floor((zpos-min3)*rdelta_z));
      xc = static_cast<ULONG>(fmod(real_prec(xc),real_prec(N1)));
      yc = static_cast<ULONG>(fmod(real_prec(yc),real_prec(N2)));
      zc = static_cast<ULONG>(fmod(real_prec(zc),real_prec(N3)));
      real_prec xx  = deltax*(static_cast<real_prec>(xc)+0.5);
      real_prec yy  = deltay*(static_cast<real_prec>(yc)+0.5);
      real_prec zz  = deltaz*(static_cast<real_prec>(zc)+0.5);
      real_prec xxf = deltax*(static_cast<real_prec>(xc)+1.5);
      real_prec yyf = deltay*(static_cast<real_prec>(yc)+1.5);
      real_prec zzf = deltaz*(static_cast<real_prec>(zc)+1.5);
      real_prec xxb = deltax*(static_cast<real_prec>(xc)-0.5);
      real_prec yyb = deltay*(static_cast<real_prec>(yc)-0.5);
      real_prec zzb = deltaz*(static_cast<real_prec>(zc)-0.5);
      int Xb=(xc==0 ? N1: xc);
      int Yb=(yc==0 ? N2: yc);
      int Zb=(zc==0 ? N3: zc);
      int Xf=(xc==N1-1 ? -1: xc);
      int Yf=(yc==N2-1 ? -1: yc);
      int Zf=(zc==N3-1 ? -1: zc);
      vector<int>i_idx{Xb-1,xc,Xf+1};
      vector<int>j_idx{Yb-1,yc,Yf+1};
      vector<int>k_idx{Zb-1,zc,Zf+1};
      vector<real_prec>MAS_xx{MAS_TSC((xxb- xpos)*rdelta_x),MAS_TSC((xx - xpos)*rdelta_x),MAS_TSC((xxf- xpos)*rdelta_x)};
      vector<real_prec>MAS_yy{MAS_TSC((yyb- ypos)*rdelta_y),MAS_TSC((yy - ypos)*rdelta_y),MAS_TSC((yyf- ypos)*rdelta_y)};
      vector<real_prec>MAS_zz{MAS_TSC((zzb- zpos)*rdelta_z),MAS_TSC((zz - zpos)*rdelta_z),MAS_TSC((zzf- zpos)*rdelta_z)};
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
      for(int ih=0;ih<MAX_MAS_DEG;++ih)
        for(int jh=0;jh<MAX_MAS_DEG;++jh)
          for(int kh=0;kh<MAX_MAS_DEG;++kh)
#ifdef USE_OMP
#pragma omp atomic update
#endif
            delta[index_3d(i_idx[ih],j_idx[jh], k_idx[kh], N3, N2)] += MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
    counter_a++;
  }
  ULONG count=0;
#pragma omp parallel for reduction(+:count)
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// ***************************************************************************************************************************************************************

void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{
#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#ifdef USE_OMP
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif  


  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];

  ULONG counter_a=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:counter_a)
#endif
  for (ULONG ig=0; ig<xp.size(); ++ig)
    {   
      real_prec xpos=xp[ig]-num_0_5*d1;
      real_prec ypos=yp[ig]-num_0_5*d2;
      real_prec zpos=zp[ig]-num_0_5*d3;
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1));
      ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2));
      ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3));
      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
      ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
      ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
      ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
      real_prec xc = static_cast<real_prec>(i);
      real_prec yc = static_cast<real_prec>(j);
      real_prec zc = static_cast<real_prec>(k);
      real_prec dx = (xpos-min1)/d1-xc;
      real_prec dy = (ypos-min2)/d2-yc;
      real_prec dz = (zpos-min3)/d3-zc;
      real_prec tx = num_1-dx;
      real_prec ty = num_1-dy;
      real_prec tz = num_1-dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(j+N2*i)]    += tx*ty*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(j+N2*ii)]   += dx*ty*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(jj+N2*i)]   += tx*dy*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(j+N2*i)]   += tx*ty*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(jj+N2*ii)]  += dx*dy*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(j+N2*ii)]  += dx*ty*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(jj+N2*i)]  += tx*dy*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(jj+N2*ii)] += dx*dy*dz;
      counter_a++;
    }
  ULONG count=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// **********************************************************************************************************************
// **********************************************************************************************************************
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, vector<real_prec>&delta)
{
#ifndef _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
#endif

  ULONG count_previous=0;
#pragma omp parallel for reduction(+:count_previous)
  for(ULONG ip=0;ip<delta.size();++ip)
    count_previous+=delta[ip];
  ULONG counter_a=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:counter_a)
#endif
  for (ULONG ig=0; ig<xp.size(); ++ig)
    {
      real_prec xpos=xp[ig];
      real_prec ypos=yp[ig];
      real_prec zpos=zp[ig];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
      ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2));
      ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3));
      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
#pragma omp atomic update
      delta[k+N3*j+N3*N2*i]++;
      counter_a++;
    }

  ULONG count=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count)
#endif
  for(ULONG ip=0;ip<delta.size();++ip)
    count+=delta[ip];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid  = "<<count<<endl;
}
// ***************************************************************************************************************************************************************
void read_bin_file(string input_file,  int Nbytes_header, int Nbytes_data, int &Nparts, vector<float>&prop)
{
  cout<<"Reading input file "<<input_file<<endl;
  ifstream input;
  input.open(input_file.c_str(), ios::binary| ios::in);
  if(!input)
    {
      cout<<RED<<"File not found. Code exits here."<<RESET<<endl;
      exit(0);
    }
  input.read((char *)&Nparts,Nbytes_header);
  while(!input.eof())
    {
      float data_cat;
      input.read((char *)&data_cat,Nbytes_data);
      prop.push_back(data_cat);
    }
  input.close();
}

// ***************************************************************************************************************************************************************
void read_bin_array(string input_file, vector<real_prec>&out)
{
  cout<<"Reading input file "<<input_file<<endl;
  ifstream input(input_file.c_str(), ios::binary| ios::in);

  if(!input)
    {
      cout<<RED<<"File not found. Code exits here."<<RESET<<endl;
      exit(0);
    }

  vector<float>prop;
  while(!input.eof())
    {
      float data_cat;
      input.read((char *)&data_cat,sizeof(float));
      prop.push_back(data_cat);
    }
  for(ULONG i=0;i<out.size();++i)
      out[i]=static_cast<real_prec>(prop[i]);
  prop.clear();prop.shrink_to_fit();

  input.close();
}

// ***************************************************************************************************************************************************************
void write_to_binary(string FNAME, vector<real_prec>&out)
{
  cout<<"Writting in binary file "<<FNAME<<endl;
  ofstream outf(FNAME,ios::out | ios::binary);
  for(ULONG i=0;i<out.size();++i)
    {
      float new_out=static_cast<float>(out[i]);
      outf.write((char *)&new_out,sizeof(float));
    }
  outf.close();
}

// ***************************************************************************************************************************************************************
#ifdef USE_SUPERSAMPLING
void supersampling(int Nft_HR, int Nft_LR, int imas, vector<real_prec>&HR_field, vector<real_prec>&LR_field)
{
  cout<<"Using supersampling "<<endl;
  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  ULONG NGRID_LR=LR_field.size();  

  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }

  /*
  real_prec mean_d=0;
  get_overdens(HR_field, HR_field, mean_d);
  
  cout<<"Mean number density = "<<mean_d<<endl;
  */
  
    complex_prec *HR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_HR*sizeof(real_prec));
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][REAL]=0;
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][IMAG]=0;

   do_fftw_r2c(Nft_HR,HR_field,HR_FOURIER);



   
   vector<real_prec> correction(Nft_HR,0);
   //#pragma omp parallel for
   for(int i = 0 ; i < Nft_HR; ++i )
    {
      real_prec coords=static_cast<real_prec>( i<=Nft_HR/2? i: i-(Nft_HR));
      real_prec xx=coords*M_PI/static_cast<real_prec>(Nft_HR);
      correction[i]= i==0 ? 1.0 : pow(sin(xx)/static_cast<real_prec>(xx),imas+1) ;
    }
   

   vector<int> low_pass_filter(Nft_HR,0);
   //#pragma omp parallel for
   for(int i = 0 ; i < Nft_HR; ++i )
     if(i <= Nft_LR/2 || i >= 1+Nft_HR-Nft_LR/2)
       low_pass_filter[i]= 1.0;

  
   real_prec we=0;
#pragma omp parallel for collapse(3) reduction(+:we)
  for(int i=0;i<Nft_HR;++i)
    for(int j=0;j< Nft_HR;++j)
      for(int k=0;k<Nft_HR/2+1;++k)
	{
	  real_prec corr=correction[i]*correction[j]*correction[k];
	  real_prec filt=low_pass_filter[i]*low_pass_filter[j]*low_pass_filter[k];
	  ULONG index_hr=index_3d(i,j,k,Nft_HR,Nft_HR/2+1);

	  if(i==Nft_HR/2 && j==Nft_HR/2 && k==Nft_HR/2)
	    cout<<HR_FOURIER[index_hr][IMAG]<<endl;


	  HR_FOURIER[index_hr][REAL] *= filt/corr;
	  HR_FOURIER[index_hr][IMAG] *= filt/corr;
	  we+=filt;
        }

  vector<int> findex(Nft_LR,0);
#pragma omp parallel for
  for(ULONG i = 0 ; i < Nft_LR; ++i ){
    findex[i]= i < Nft_LR/2 +1 ? i :   i + Nft_HR-Nft_LR ;
  }
  
  complex_prec *LR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_LR*sizeof(real_prec));
#pragma omp parallel for
  for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][IMAG]=0;
  
  

  for(int i=0;i<Nft_LR;++i)
    {
      int i_hr=findex[i];
      for(int j=0;j< Nft_LR;++j)
	{
	  int j_hr=findex[j];
	  for(int k=0;k<Nft_LR/2+1;++k)
	    {
	      int k_hr=findex[k];
	      ULONG index_lr=index_3d(i,   j,   k,   Nft_LR, Nft_LR/2+1);
	      ULONG index_hr=index_3d(i_hr,j_hr,k_hr,Nft_HR ,Nft_HR/2+1);
	      LR_FOURIER[index_lr][REAL] = HR_FOURIER[index_hr][REAL];
	      LR_FOURIER[index_lr][IMAG] = HR_FOURIER[index_hr][IMAG];
	    }
	}
    }
  findex.clear();findex.shrink_to_fit();
  
  
  
  /*    
	for (ULONG i=0 ; i<Nft_LR/2+1;++i)
	for (ULONG j=0 ; j<Nft_LR/2+1;++j)
	for (ULONG k=0 ; k< Nft_LR/2 +1 ;++k)
	{  
	ULONG ii = Nft_LR - i;
	ULONG jj = Nft_LR - j;
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
	}
	

  for (ULONG i=1 ; i<Nft_LR;++i)
    for (ULONG j=1 ; j<Nft_LR/2;++j)
      {  
	ULONG ii = Nft_LR - i;
	ULONG jj = Nft_LR - j;
	ULONG k=Nft_LR/2;
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
	
      }
  

  for (ULONG i=1 ; i<Nft_LR;++i)
    for (ULONG j=1 ; j<Nft_LR/2;++j)
      {  
	ULONG k=0;
	ULONG ii = Nft_LR - i;
	ULONG jj = Nft_LR - j;
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
      }

      
  for (ULONG i=1 ; i<Nft_LR;++i)
    for (ULONG k=1 ; k<Nft_LR/2+1;++k)
      {  
	ULONG j=0;
	ULONG ii = Nft_LR - i;
	ULONG jj = j;
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
      }
   
  for (ULONG j=1 ; j<Nft_LR;++j)
    for (ULONG k=1 ; k<Nft_LR/2+1;++k)
      {  
	ULONG i=0;
	ULONG ii = i;
	ULONG jj = Nft_LR - j;
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
	LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
      }
  
  
 
  
  for (ULONG i=1 ; i<Nft_LR/2;++i)
    {  
      ULONG k=Nft_LR/2;
      ULONG j=Nft_LR/2;
      ULONG ii = Nft_LR - i;
      ULONG jj = Nft_LR - j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }
  

  for (ULONG i=1 ; i<Nft_LR/2;++i)
    {  
      ULONG j=Nft_LR/2;
      ULONG k=0;
      ULONG ii = Nft_LR - i;
      ULONG jj = j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }
  

  for (ULONG i=1 ; i<Nft_LR/2;++i)
    {  
      ULONG j=0;
      ULONG k=Nft_LR/2;
      ULONG ii = Nft_LR - i;
      ULONG jj = j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }
  


  for (ULONG j=1 ; j<Nft_LR/2;++j)
    {  
      ULONG i=0;
      ULONG k=Nft_LR/2;
      ULONG ii = i;
      ULONG jj = Nft_LR - j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }
  
 
  
  for (ULONG i=1 ; i<Nft_LR/2;++i)
    {  
      ULONG j=0;
      ULONG k=0;
      ULONG ii = Nft_LR - i;
      ULONG jj =j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }

  
  for (ULONG j=1 ; j<Nft_LR/2;++j)
    {  
      ULONG i=0;
      ULONG k=0;
      ULONG ii = i;
      ULONG jj = Nft_LR - j;
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][REAL] =  LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][REAL];
      LR_FOURIER[ index_3d(ii,jj,k, Nft_LR, Nft_LR/2+1)][IMAG] = -LR_FOURIER[index_3d(i,j,k, Nft_LR, Nft_LR/2+1)][IMAG];
    }
  */



  real_prec a, b;
  int ii, jj, kk;
  ULONG ind;
  
  ii=0; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=0; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=0; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=0; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=Nft_LR/2; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=Nft_LR/2; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=Nft_LR/2; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  ii=Nft_LR/2; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  
  


  do_fftw_c2r(Nft_LR ,LR_FOURIER,LR_field);

  
  /*

  real_prec factor=static_cast<real_prec>(NGRID_LR)/static_cast<real_prec>(2.*we) ;
#pragma omp parallel for
  for(ULONG i=0;i<NGRID_LR;i++)
    LR_field[i]=mean_d*(1+LR_field[i]/factor);
  */
  
  fftw_free(LR_FOURIER);
  fftw_free(HR_FOURIER);
}
#endif


real_prec gsl_inter_new(const vector<gsl_real> &xx, const vector<gsl_real> &yy, gsl_real xn){
  real_prec ans;
  int n=xx.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, &xx[0], &yy[0], n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans;
}

// ***************************************************************************
// ***************************************************************************

void subtract_cv(string file_IC, string file_Pkt, string file_H, ULONG Nft, int imas)
{
    ULONG Ngrid=Nft*Nft*Nft;
    ULONG NTT=Nft*Nft*(Nft/2+1);



    complex_prec *IC_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
    complex_prec *H_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)IC_FOURIER[i][REAL]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)IC_FOURIER[i][IMAG]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)H_FOURIER[i][REAL]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)H_FOURIER[i][IMAG]=0;



    // Theoretical redshift given at z=0.041
    real_prec ic_factor = 9017.46774669688 ; // multiply (divide) delta IC at z=120 (0) to take it to z=0 (120)
    real_prec factor =   1.04498241888233;  // multiply (divide) dektaIC to take it from at z=0.041 to z=0
    real_prec h_factor = 3224.62745724095; // takes the IC from z=120 to the halos redshift z=1.04
    real_prec norm=(Ngrid*Ngrid) / pow(Lside,3);  // Normalization of power spectrum. Multiply by this factor to un-normalize the theoretical power

    vector<real_prec>prop;
    int nlines_p=read_file(file_Pkt,prop,1);
    int ncols=(static_cast<ULONG>(prop.size()/nlines_p));
    vector<double>kv(nlines_p,0);
    vector<double>pv(nlines_p,0);
    for(int i=0;i<nlines_p;++i)
      {
        kv[i]=static_cast<double>(prop[ncols*i]);
        pv[i]=static_cast<double>(prop[1+ncols*i])*norm*factor/ic_factor;  //theoretical P(k) at z=120 unnormalized
      }
    prop.clear();prop.shrink_to_fit();



    real_prec deltak=2.*acos(-1.0)/Lside;
    vector<real_prec> coords(Nft,0);
    for(ULONG i=0;i<Nft ;++i)
      coords[i]=deltak*(i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));


    vector<int> nmodes(Nft/2,0);


    vector<real_prec> IC_field(Ngrid,0);
    vector<real_prec> IC_P(Nft/2,0);
    vector<real_prec> IC_PN(Nft/2,0);
    read_bin_array(file_IC,IC_field);
    real_prec icmean;
    get_overdens(IC_field,IC_field,icmean);
    cout<<"Mean IC number density = "<<icmean<<endl;

    vector<real_prec> H_field(Ngrid,0);
    vector<real_prec> H_P(Nft/2,0);
    vector<real_prec> H_PN(Nft/2,0);
    read_bin_array(file_H,H_field);
    real_prec hmean;
    get_overdens(H_field,H_field,hmean);
    cout<<"Mean h number density = "<<hmean<<endl;

    do_fftw_r2c(Nft,IC_field, IC_FOURIER);
    do_fftw_r2c(Nft, H_field,  H_FOURIER);

    real_prec eps=1;
    real_prec large_scale_bias = 1.422;

    for(int i=0;i<Nft;++i)
      for(int j=0;j< Nft;++j)
        for(int k=0;k<Nft/2+1;++k)
          {
            int lp=index_3d(i,j,k,Nft,Nft/2+1);
            if(lp>0)
            {
            real_prec kmod=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
            int ik=floor(kmod/deltak);

            real_prec h_re=H_FOURIER[lp][REAL];
            real_prec h_im=H_FOURIER[lp][IMAG];
            real_prec Amp_h = sqrt(pow(h_re,2)+pow(h_im,2));
            real_prec phi_h=atan2(H_FOURIER[lp][IMAG],H_FOURIER[lp][REAL]);

            real_prec ic_re=IC_FOURIER[lp][REAL];
            real_prec ic_im=IC_FOURIER[lp][IMAG];
            real_prec Amp_ic = sqrt(pow(ic_re,2)+pow(ic_im,2));
            real_prec phi_ic=atan2(IC_FOURIER[lp][IMAG],IC_FOURIER[lp][REAL]);

            real_prec Power_th = (kmod==0 ? 0 : gsl_inter_new(kv,pv,kmod));
            real_prec delta_th = sqrt(Power_th);


            if(ik<Nft/2)
              {
                    nmodes[ik]++;
                    H_P[ik] += pow(Amp_h,2);
                    IC_P[ik]+= pow(Amp_ic,2);
              }

            if(kmod<2.0){
              Amp_ic = delta_th;
              IC_FOURIER[lp][REAL] = Amp_ic*cos(phi_ic); // Amp_ic - alpha = delta_th
              IC_FOURIER[lp][IMAG] = Amp_ic*sin(phi_ic);

              Amp_h = large_scale_bias* sqrt(h_factor)* delta_th;
              H_FOURIER[lp][REAL] = Amp_h*cos(phi_h);
              H_FOURIER[lp][IMAG] = Amp_h*sin(phi_h);
            }

            if(ik<Nft/2)
              {
                H_PN[ik] += pow(Amp_h,2);
                IC_PN[ik]+= ow(Amp_ic,2);
              }
            }
          }




    // If the IC were to be built with the input power, the poower spectrum of the subtraction delta_iC - delta_th would be a white noise.
    // This is not the case where. The power of tha tdifferenc has structure coming from the displacements applied to the IC.

    ofstream sal;
    sal.open("power.txt");

    for(int i=0;i<Nft/2;++i)
      sal<<deltak*(i+0.5)<<"  "<<IC_P[i]/(nmodes[i]*norm)<<"  "<<H_P[i]/(nmodes[i]*norm)<<"   "<<IC_PN[i]/(nmodes[i]*norm)<<"  "<<H_PN[i]/(nmodes[i]*norm)<<endl;
    sal.close();

    IC_field.resize(Ngrid,0);
    do_fftw_c2r(Nft ,IC_FOURIER,IC_field);
    for(ULONG i=0; i<Ngrid;++i)
      IC_field[i]=icmean*(1+IC_field[i]);

    string output_den="SLICS_IC_LOS"+to_string(los)+"_Nres"+to_string(Nft)+"_MAS"+to_string(imas)+"ac_cvc.dat";
    write_to_binary(output_den,IC_field);

    H_field.resize(Ngrid,0);
    do_fftw_c2r(Nft ,H_FOURIER, H_field);
    for(ULONG i=0; i<Ngrid;++i)
      H_field[i]=hmean*(1.+H_field[i]);

    output_den="../HALOS/SLICS_HALOS_LOS"+to_string(los)+"_Nres"+to_string(Nft)+"_MAS0"+"cvc.dat";
    write_to_binary(output_den,H_field);

}










// **********************************************




void wn(string file_IC, string file_Pkt, ULONG Nft, int imas)
{
    ULONG Ngrid=Nft*Nft*Nft;
    ULONG NTT=Nft*Nft*(Nft/2+1);



    complex_prec *IC_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
    complex_prec *H_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)IC_FOURIER[i][REAL]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)IC_FOURIER[i][IMAG]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)H_FOURIER[i][REAL]=0;
  #pragma omp parallel for
    for(ULONG i=0;i<NTT;++i)H_FOURIER[i][IMAG]=0;



    // Theoretical redshift given at z=0.041
    real_prec ic_factor = 9017.46774669688 ; // multiply (divide) delta IC at z=120 (0) to take it to z=0 (120)
    real_prec factor =   1.04498241888233;  // multiply (divide) dektaIC to take it from at z=0.041 to z=0
    real_prec h_factor = 3224.62745724095; // takes the IC from z=120 to the halos redshift z=1.04
    real_prec norm=(Ngrid*Ngrid) / pow(Lside,3);  // Normalization of power spectrum. Multiply by this factor to un-normalize the theoretical power

    vector<real_prec>prop;
    int nlines_p=read_file(file_Pkt,prop,1);
    int ncols=(static_cast<ULONG>(prop.size()/nlines_p));
    vector<double>kv(nlines_p,0);
    vector<double>pv(nlines_p,0);
    for(int i=0;i<nlines_p;++i)
      {
        kv[i]=static_cast<double>(prop[ncols*i]);
        pv[i]=static_cast<double>(prop[1+ncols*i])*norm*factor/ic_factor;  //theoretical P(k) at z=120 unnormalized
      }
    prop.clear();prop.shrink_to_fit();



    real_prec deltak=2.*acos(-1.0)/Lside;
    vector<real_prec> coords(Nft,0);
    for(ULONG i=0;i<Nft ;++i)
      coords[i]=deltak*(i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));


    vector<int> nmodes(Nft/2,0);


    vector<real_prec> IC_field(Ngrid,0);
    vector<real_prec> WN_field(Ngrid,0);
    vector<real_prec> IC_P(Nft/2,0);
    vector<real_prec> WN_P(Nft/2,0);
    read_bin_array(file_IC,IC_field);
    real_prec icmean;
    get_overdens(IC_field,IC_field,icmean);
    cout<<"Mean IC number density = "<<icmean<<endl;


    do_fftw_r2c(Nft,IC_field, IC_FOURIER);

    real_prec eps=1;
    real_prec large_scale_bias = 1.422;

    for(int i=0;i<Nft;++i)
      for(int j=0;j< Nft;++j)
        for(int k=0;k<Nft/2+1;++k)
          {
            int lp=index_3d(i,j,k,Nft,Nft/2+1);
            if(lp>0)
            {
            real_prec kmod=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
            int ik=floor(kmod/deltak);


            real_prec ic_re=IC_FOURIER[lp][REAL];
            real_prec ic_im=IC_FOURIER[lp][IMAG];
            real_prec Amp_ic = sqrt(pow(ic_re,2)+pow(ic_im,2));
            real_prec phi_ic=atan2(IC_FOURIER[lp][IMAG],IC_FOURIER[lp][REAL]);

            real_prec Power_th = (kmod==0 ? 0 : gsl_inter_new(kv,pv,kmod));
            real_prec delta_th = sqrt(Power_th);


            if(ik<Nft/2)
              {
                    nmodes[ik]++;
                    IC_P[ik]+= pow(Amp_ic,2);
              }

            if(kmod<2.0){
              Amp_ic = delta_th;
              IC_FOURIER[lp][REAL] = Amp_ic*cos(phi_ic); // Amp_ic - alpha = delta_th
              IC_FOURIER[lp][IMAG] = Amp_ic*sin(phi_ic);

            }

            if(ik<Nft/2)
              {
                WN[ik]+= ow(Amp_ic,2);
              }
            }
          }




    // If the IC were to be built with the input power, the poower spectrum of the subtraction delta_iC - delta_th would be a white noise.
    // This is not the case where. The power of tha tdifferenc has structure coming from the displacements applied to the IC.

    ofstream sal;
    sal.open("power.txt");

    for(int i=0;i<Nft/2;++i)
      sal<<deltak*(i+0.5)<<"  "<<IC_P[i]/(nmodes[i]*norm)<<"  "<<H_P[i]/(nmodes[i]*norm)<<"   "<<IC_PN[i]/(nmodes[i]*norm)<<"  "<<H_PN[i]/(nmodes[i]*norm)<<endl;
    sal.close();

    IC_field.resize(Ngrid,0);
    do_fftw_c2r(Nft ,IC_FOURIER,IC_field);
    for(ULONG i=0; i<Ngrid;++i)
      IC_field[i]=icmean*(1+IC_field[i]);

    string output_den="WN_SLICS_IC_LOS"+to_string(los)+"_Nres"+to_string(Nft)+"_MAS"+to_string(imas)+"ac_cvc.dat";
    write_to_binary(output_den,WN_field);

}

// ***************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************
// ********************************************MAIN MAIN MAIN MAIN ***********************************************************************************************
// ***************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************
int main(int argc, char *argv[]){
  
  int los = 980;
  int Nft= 180;  // Number of cell/dimention high res
  int imas= 1;  // MAS
  string file_IC="SLICS_IC_LOS"+to_string(los)+"_Nres"+to_string(Nft)+"_MAS"+to_string(imas)+"ac.dat";
  string file_Pkt="CAMB_Z0.042_high_k_linear_WMAP9.txt";
  string file_H="../HALOS/SLICS_HALOS_LOS980_Nres180_MAS0.dat";

 // subtract_cv(file_IC, file_Pkt,filee_H, Nftm¡,imas);
  


 wn(file_IC, file_Pkt, Nft,imas);


}
















