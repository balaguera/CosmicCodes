/* Code to read IC's from the SLIC simulation. Andrés Balaguera-Antolínez IAC 2020 */


// This option is only available if FFTW is propery working

#define USE_DOWNSAMPLING
//#define USE_AVERAGE

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

#ifdef USE_DOWNSAMPLING
#include <fftw3.h>
#endif


#define USE_OMP  // Use OMP
#define _NTHREADS_ 16

using namespace std;
// Variables used in FFTW
#define REAL 0
#define IMAG 1
#define ic_rank 3
// Variables used in Supersampling and FFTW
#define ss_factor static_cast<int>(2)  // Used for ssuper-sampling
#define MAX_MAS_DEG static_cast<int>(3)
#define MAX_MAS_DEG_PSC static_cast<int>(4)


#define ULONG unsigned long
//#define DOUBLE_PREC
#ifdef DOUBLE_PREC
#define fftw_real double
#define real_prec fftw_real
#define complex_prec fftw_complex
#else
#define real_prec float
#define fftwf_real real_prec
#define complex_prec fftwf_complex
#endif

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
#define BLUE     "\033[34m"      /* Red */
#define RESET   "\033[0m"


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
#define Lside_sim static_cast<real_prec>(1000.0)  //Lenght of simulation volume
#define Lside static_cast<real_prec>(1000.0)  //Comoving lenght (in Mpc/h) of cosmological volume
#define redshift (1.04)   // Redshift

#ifdef _GET_INTERPOLATED_FIELD_FROM_HALO_CAT_
#undef _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#endif


template<typename T> class fftw_array
{
 public:
  T *data;
  
  fftw_array(long size)
    {
      
#ifdef SINGLE_PREC
      data = (T *) fftwf_malloc(size*sizeof(T));
#endif
#ifdef DOUBLE_PREC
      data = (T *) fftw_malloc(size*sizeof(T));
#endif
    }

  
  ~fftw_array()
    {
#ifdef SINGLE_PREC
      fftwf_free(data);
#endif
#ifdef DOUBLE_PREC
      fftw_free(data);
#endif
    }
  /* implicit conversion: works for functions, but not for templates. For the latter case write: <array_name>.data */
  operator T*()
  {
    return data;
  }
 
};



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
inline ULONG index_3d(ULONG i, ULONG j, ULONG k, ULONG Nj, ULONG Nk)
{
  return static_cast<ULONG>(k)+static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
// **********************************************************************************************************************
void index2coords(ULONG N, ULONG index, ULONG  &XG, ULONG &YG, ULONG &ZG )
{// Get the index in each direction : F-order cells (N,index, X,Y,Z). -Corder (N,index, Z,Y,X)
  XG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  YG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  ZG = index;
}

// **********************************************************************************************************************
// **********************************************************************************************************************
#ifdef USE_DOWNSAMPLING

void do_fftw_r2c(int Nft, vector<real_prec>&in, complex_prec *out)
{
  int *n=(int *)fftw_malloc(ic_rank*sizeof(float));
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif

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

#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif


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

  real_prec y=fabs(x);
  real_prec ans;
  if(y<0.5)
    ans= 0.75-pow(y,2);
  else if(y>=0.5 && y<1.5)
    ans = 0.5*pow(1.5 -y,2);
  else
    ans = 0.0;

  return ans;

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
  real_prec count_previous=0;

#pragma omp parallel for reduction (+:count_previous)
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
      int xc = static_cast<int>(floor((xpos-min1)*rdelta_x)); // indices of the cell of the particle
      int yc = static_cast<int>(floor((ypos-min2)*rdelta_y));
      int zc = static_cast<int>(floor((zpos-min3)*rdelta_z));
      // xc = static_cast<ULONG>(fmod(real_prec(xc),real_prec(N1)));
      // yc = static_cast<ULONG>(fmod(real_prec(yc),real_prec(N2)));
      // zc = static_cast<ULONG>(fmod(real_prec(zc),real_prec(N3)));
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
      
      for(int ih=0;ih<MAX_MAS_DEG;++ih)
        for(int jh=0;jh<MAX_MAS_DEG;++jh)
          for(int kh=0;kh<MAX_MAS_DEG;++kh)
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


#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

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
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec>&we, vector<real_prec>&delta, bool weight_with_field)
{


#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif


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

#pragma omp parallel for reduction(+:counter_a)
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
      // i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
      //      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
      //  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
      
      if(true==weight_with_field)
#pragma omp atomic update
        delta[k+N3*j+N3*N2*i]+=we[ig];
      else
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
  //  vector<float>prop;
  
  if(!input)
    
    {
      cout<<RED<<"File not found. Code exits here."<<RESET<<endl;
      exit(0);
    }
  
  
  while(!input.eof())
    {
      float data_cat;
      input.read((char *)&data_cat,sizeof(float));
      out.push_back(data_cat);
      //  cout<<data_cat<<endl;
    }
  
  /*  
  for(ULONG i=0;i<out.size();++i)
    out[i]=static_cast<real_prec>(prop[i]);
  prop.clear();prop.shrink_to_fit();
  */
  
  input.close();
}





// ***************************************************************************************************************************************************************
void read_bin_array_unitsim(string input_file, vector<real_prec>&out)
{
  cout<<"Reading input file "<<input_file<<" explicitely for UNITSIM density fields"<<endl;
  ifstream input(input_file.c_str(), ios::binary| ios::in);
  
  if(!input)
    {
      cout<<RED<<"File not found. Code exits here."<<RESET<<endl;
      exit(0);
    }
  
  cout.precision(12);
  input.seekg (0, ios::end);
  ULONG length = input.tellg();
  length/=sizeof(float);
  cout<<"Lenght "<<length<<endl;
  input.seekg (0, ios::beg);
  
  int Ndim;
  input.read((char *)&Ndim,sizeof(Ndim));
  cout<<"0 ="<<Ndim<<endl;
  input.read((char *)&Ndim,sizeof(Ndim));
  cout<<"1 ="<<Ndim<<endl;
  input.read((char *)&Ndim,sizeof(Ndim));
  cout<<"2 ="<<Ndim<<endl;
  input.read((char *)&Ndim,sizeof(Ndim));
  cout<<"3 ="<<Ndim<<endl;
  int Naux;
  input.read((char *)&Naux,sizeof(Naux));
  cout<<"4 ="<<Naux<<endl;

  ULONG Ng=static_cast<ULONG>(Ndim);
  ULONG Ngrid=(Ng*Ng)*Ng;
  cout<<"Ngrid= "<<Ndim<<"  "<<length-Ngrid<<endl;
  
  cout<<BLUE<<"Reading...."<<RESET<<endl;
  double counter=0;
  for(ULONG k=0; k<Ng; ++k)//Z is the slow varying component in the fortran (column) ordering
    for(ULONG j=0; j<Ng; ++j)
      for(ULONG ip=0; ip<=Ng+1; ++ip)// x is the fast varying cordinate in fortran
	{
	  if(ip>=1 && ip<=Ng )
	    {	
	      input.read((char *)&out[k+Ng*j+Ng*Ng*(ip-1)],sizeof(float));   // Here we change fortran to c like ordering
	      cout<<BLUE<<(100.0*(counter)/(static_cast<double>(Ngrid)))<<"  % \r";cout.flush();
	      counter++;
	    }
	  else
	    input.read((char *)&Naux,sizeof(Naux));  // this reads the Ng*Ng*4 which is an input of fortran unformatted at the beginning

	}
  input.close();
  
  double mean=0;
#pragma omp parallel for reduction(+:mean)
  for(ULONG i=0; i < out.size();++i)
    mean+=out[i];
  mean=static_cast<double>(mean/static_cast<double>(Ngrid));
  cout<<"MEAN = "<<mean<<endl;
  double var=0;
#pragma omp parallel for reduction(+:var)
  for(ULONG i=0; i < out.size();++i)
    var+=pow(out[i]-mean,2);
  var/=static_cast<double>(Ngrid);
  cout<<BLUE<<"Variance = "<<var<<endl;
  cout<<"Standard deviation = "<<sqrt(var)<<RESET<<endl;
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
#ifdef USE_DOWNSAMPLING


void downsampling(ULONG Nft_HR, ULONG Nft_LR, int imas, vector<real_prec>&HR_field, vector<real_prec>&LR_field)
{
  int NTHREADS=16;
  omp_set_num_threads(NTHREADS);
  cout<<"Applying low  pass filter "<<Nft_HR<<" to "<<Nft_LR<<endl;  

  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  ULONG NGRID_LR=LR_field.size();
  real_prec delta_box= Lside/static_cast<real_prec>(Nft_LR);
  real_prec delta_box_hr= Lside/static_cast<real_prec>(Nft_HR);
  real_prec delta_k= 2*M_PI/Lside;
  
  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }

  cout<<"delta_box = "<<delta_box<<"  delta_box_HR = "<<delta_box_hr<<" delta_k = "<<delta_k<<endl;


    complex_prec *HR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_HR*sizeof(real_prec));
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)
   {
       HR_FOURIER[i][REAL]=0;
       HR_FOURIER[i][IMAG]=0;
    }

   // DO FFT of the high-res density field

   do_fftw_r2c(Nft_HR,HR_field,HR_FOURIER);

    // Get the MAS correction for the H-Res density field
   vector<real_prec> correction(Nft_HR,1.0);
#pragma omp parallel for
   for(ULONG i = 1 ; i < Nft_HR; ++i )
   {
      int coords= i<Nft_HR/2+1? i: i-Nft_HR;
      real_prec xx=static_cast<real_prec>(coords)*M_PI/static_cast<real_prec>(Nft_HR);
      correction[i]= 1.0;//pow(sin(xx)/static_cast<real_prec>(xx), imas+1);  // applies directly to NGP
   }
   
   vector<real_prec> kmodes_lr(Nft_LR,0.0);
#pragma omp parallel for
   for(ULONG i = 0 ; i < Nft_LR; ++i )
     {
       ULONG coords= i < Nft_LR/2+1? i: i-Nft_LR;
       kmodes_lr[i]=static_cast<real_prec>(coords)*delta_k;
     }
   
   complex_prec *LR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_LR*sizeof(real_prec));
#pragma omp parallel for
   for(ULONG i=0;i<NTT_LR;++i)
     {
       LR_FOURIER[i][REAL]=0;
       LR_FOURIER[i][IMAG]=0;
    }
   // Get the FT of the Low-res density field
   real_prec alpha=static_cast<double>(Nft_HR)/static_cast<double>(Nft_LR);
   real_prec delta_shift= 0.5*(alpha-1.0)*delta_box_hr; //Yu
   
#pragma omp parallel for
   for(ULONG i=0;i<Nft_LR;++i)
     {
       ULONG i_hr = i < Nft_LR/2 +1 ? i : i+Nft_HR-Nft_LR;
       for(ULONG j=0;j< Nft_LR;++j)
	 {
	   ULONG j_hr = j < Nft_LR/2 +1 ? j : j+Nft_HR-Nft_LR;
	   for(ULONG k=0;k<Nft_LR/2+1;++k)
	     {
	       ULONG k_hr=k ; //< Nft_LR/2 +1 ? k : k+Nft_HR-Nft_LR;
	       ULONG index_lr=index_3d(i,   j,   k,   Nft_LR, Nft_LR/2+1);
	       ULONG index_hr=index_3d(i_hr,j_hr,k_hr,Nft_HR ,Nft_HR/2+1);
	       real_prec corr=correction[i_hr]*correction[j_hr]*correction[k_hr];
	       real_prec shift=delta_shift*(kmodes_lr[i] + kmodes_lr[j] + kmodes_lr[k]);
	       real_prec new_real=HR_FOURIER[index_hr][REAL]/corr;
	       real_prec new_imag=HR_FOURIER[index_hr][IMAG]/corr;
	       LR_FOURIER[index_lr][REAL] = (cos(shift)*new_real - sin(shift)*new_imag); //This follows the convention delta_2-> exp(iks)*delta_2
	       LR_FOURIER[index_lr][IMAG] = (cos(shift)*new_imag + sin(shift)*new_real);
	     }
	 }
     }
   
   fftw_free(HR_FOURIER);
   
   // Set Ny freq of the FT of the L-res density field to real:
   real_prec a, b;
   ULONG ii, jj, kk;
   ULONG ind;
   real_prec fac=1.0;
   
   ii=0; jj=0; kk=Nft_LR/2; // This requires no reflection in negative z
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=0; jj=Nft_LR/2; kk=0;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=0; jj=Nft_LR/2; kk=Nft_LR/2;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=Nft_LR/2; jj=0; kk=0;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=Nft_LR/2; jj=0; kk=Nft_LR/2;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=Nft_LR/2; jj=Nft_LR/2; kk=0;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   
   ii=Nft_LR/2; jj=Nft_LR/2; kk=Nft_LR/2;
   ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
   a  =  LR_FOURIER[ind][REAL];
   b  =  fac*LR_FOURIER[ind][IMAG];
   LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
   LR_FOURIER[ind][IMAG]=0.0;
   

   
   // Transform to L-res density field.
   do_fftw_c2r(Nft_LR ,LR_FOURIER,LR_field);
   fftw_free(LR_FOURIER);

}
#endif
// ***************************************************************************************************************************************************************

void average(ULONG Nft_HR, ULONG Nft_LR, vector<real_prec>&HR_field, vector<real_prec>&LR_field)
{
  cout<<"Averaging from high to low res"<<endl;
  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  
  ULONG NGRID_HR=static_cast<ULONG>(Nft_HR*Nft_HR*Nft_HR);
  ULONG NGRID_LR=static_cast<ULONG>(Nft_LR*Nft_LR*Nft_LR);
  
  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }
  
  real_prec delta=Lside/static_cast<real_prec>(Nft_LR);
  real_prec delta_hr=Lside/static_cast<real_prec>(Nft_HR);
  
  
  vector<real_prec>zcoord(NGRID_HR,0);
  vector<real_prec>ycoord(NGRID_HR,0);
  vector<real_prec>xcoord(NGRID_HR,0);
  
  for(ULONG i=0;i<NGRID_HR; ++i)
    {
      ULONG xc,yc,zc;
      index2coords(Nft_HR,i,xc,yc,zc);
      xcoord[i]=(static_cast<real_prec>(zc)+0.5)*delta_hr;
      ycoord[i]=(static_cast<real_prec>(yc)+0.5)*delta_hr;
      zcoord[i]=(static_cast<real_prec>(xc)+0.5)*delta_hr;
    }
  
  getDensity_NGP(Nft_LR,Nft_LR,Nft_LR,Lside,Lside,Lside,delta,delta,delta,0,0,0,xcoord,ycoord,zcoord,HR_field,LR_field,true);
  
  for(ULONG i=0;i<NGRID_LR; ++i)
    LR_field[i]/=static_cast<double>(pow(Nft_HR/Nft_LR,3));   //average                                                                                 
  
}

// ***************************************************************************************************************************************************************
// ********************************************MAIN MAIN MAIN MAIN ***********************************************************************************************
// ***************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************
int main(int argc, char *argv[]){
  
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  if(argc<5)
    {
      cout<<RED<<"Error: code expects 4  parameters. Only "<<argc-1<<" provided."<<RESET<<endl;
      exit(0);
    }
  cout<<RED<<"UNITSIM DM FIELD"<<RESET<<endl;
  int MAS=1;
  ULONG Nres_hr= atoi(argv[1]);  // Number of cell/dimention
  ULONG Nres_lr= atoi(argv[2]);  // Number of cell/dimention
  string file_in = argv[3];
  ULONG nIC  =atoi(argv[4]);
  
  if(Nres_lr<=0)
    {
      cout<<RED<<"Error: Mesh resolution (Nres) must be > 0."<<RESET<<endl;
      exit(0);
    }
  
  if(Nres_hr<=0)
    {
      cout<<RED<<"Error: Mesh resolution (Nres) must be > 0."<<RESET<<endl;
      exit(0);
    }
  if(MAS>3)
    {
      cout<<RED<<"Error: requested mass assignment scheme not valid."<<RESET<<endl;
      cout<<RED<<"Available: 0 (NGP), 1 (CIC), 2 (TSC), 3(PCS)"<<RESET<<endl;
      exit(0);
    }
  
  string sMAS;
  switch(MAS){case(0):sMAS="NGP";break; case(1):sMAS="CIC";break; case(2):sMAS="TSC";break; case(3):sMAS="PCS";break;}
  ULONG Ngrid=Nres_hr*Nres_hr*Nres_hr;
  cout<<BLUE<<"Allocating memory"<<endl;
  vector<real_prec>density_field(Ngrid,0);
  cout<<"DONE"<<RESET<<endl;
  string infile="/global/cscratch1/sd/abalant/UNITsim/NEW/fixedAmp_001/"+file_in;
  read_bin_array_unitsim(infile, density_field);
  cout<<"Size of array "<<density_field.size()<<"  Expected  "<<Ngrid<<endl;
  string i_den="/global/cscratch1/sd/abalant/UNITsim/NEW/fixedAmp_001/dm_field_p"+to_string(nIC)+"_Nres"+to_string(Nres_hr)+"_MAS"+to_string(MAS)+".dat";
  //  write_to_binary(i_den, density_field);
  ULONG Ngrid_lr=Nres_lr*Nres_lr*Nres_lr;
  vector<real_prec> LR_field(Ngrid_lr,0);

  cout<<density_field.size()<<"   "<<density_field[0]<<"  "<<density_field[10]<<endl;

  downsampling(Nres_hr,Nres_lr, 1+MAS, density_field, LR_field);
  string output_den="/global/cscratch1/sd/abalant/UNITsim/NEW/fixedAmp_001/dm_field_p"+to_string(nIC)+"_Nres"+to_string(Nres_lr)+"_MAS"+to_string(MAS)+"_downsampled.dat";
  cout<<LR_field.size()<<"   "<<LR_field[0]<<"  "<<LR_field[10]<<endl;
  write_to_binary(output_den,LR_field);
}

















