
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains some functions, based on GSL subroutines, used  in the measurements  **
// of the power spectrum                                                                   **
// Developer:                                                                              **
// Andres Balaguera Antolinez                                                              **
// abalant@gmail.com                                                                       **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
# include "../Headers/NumericalMethods.h"
# include "../Headers/FftwWrapper.h"
using namespace std;


//##################################################################################
//##################################################################################
//##################################################################################
void do_fftw_r2c(int Nft, vector<real_prec>in, fftw_complex *out)
{
  int *n=(int *)malloc(ic_rank*sizeof(float));
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
  fftw_plan plan_r2c=fftw_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftw_execute(plan_r2c);
  fftw_destroy_plan(plan_r2c);
  free(n);
}

//##################################################################################
//##################################################################################
void do_fftw_c2r(int Nft, fftw_complex *in, vector<real_prec>&out)
{
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
  int *n=(int *)malloc(ic_rank*sizeof(float));
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
  fftw_plan plan_c2r=fftw_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE); 
  fftw_execute(plan_c2r);
  fftw_destroy_plan(plan_c2r);
  free(n);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;i++)  // normalize the c2r transform
    out[i]/=static_cast<double>(NGRID);
  
}

//##################################################################################
void do_fftw_3d(int Nft, bool direction, complex_prec *in, complex_prec *out)
{


  // OJo que aca tengo una definicion de normali e la FT pre definida
  ULONG factor=Nft*Nft*Nft;
  
  int Ni1=static_cast<int>(Nft);
  int Ni2=static_cast<int>(Nft);
  int Ni3=static_cast<int>(Nft);
  
  fftw_plan fftp;
  if (direction == true)
    {
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,FORWARD,FFTW_OPTION);	
      fftw_execute(fftp);
    }
  else
    {
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,BACKWARD,FFTW_OPTION);	
      fftw_execute(fftp);
      
      for(ULONG i=0;i<factor;i++)
	{// normalize the c2r transform
	  out[i][REAL]/=static_cast<double>(factor);
	  out[i][IMAG]/=static_cast<double>(factor);
	}
      fftw_destroy_plan(fftp);
    }

  
}



