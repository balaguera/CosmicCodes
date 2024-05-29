////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/** @file Cwclass.cpp
 *
 *  @brief Cosmic web classification
 *  @author: Andrés Balaguera-Antolínez, Francisco-Shu Kitaura, IAC, 2017-2019
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "Cwclass.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int sign(real_prec x)
{
    int ans=1;
    if(x<0)
        ans=-1;
    return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_bias_terms(vector<real_prec>&delta)
{
#ifdef _FULL_VERBOSE_
    So.enter(__PRETTY_FUNCTION__);
#endif
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing bias terms");
#endif
  this->potential.resize(this->params._NGRID(),0);
  real_prec max_C1=0;
  real_prec min_C1=0;

#ifdef _USE_DELTA2_
  So.message_screen("Computing ð²");
  this->DELTA2.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID();++i)
    this->DELTA2[i]=delta[i]*delta[i];
  So.DONE();
#ifdef _WRITE_DELTA2_
  File.write_array(this->params._Output_directory()+"DELTA2", DELTA2);
#endif
  So.message_screen("DELTA2-term :");
  max_C1=get_max<real_prec>(DELTA2);
  So.message_screen("Maximum ð² =", max_C1);
  min_C1=get_min<real_prec>(DELTA2);
  So.message_screen("Minimum ð² =", min_C1);
#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif
#endif

#ifdef _USE_DELTA3_
  So.message_screen("Computing ð³");
  this->DELTA3.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID();++i)
    DELTA3[i]=delta[i]*delta[i]*delta[i];
  So.DONE();
#ifdef _WRITE_DELTA3_
  File.write_array(this->params._Output_directory()+"DELTA3", DELTA3);
#endif
  So.message_screen("DELTA3-term :");
  max_C1=get_max<real_prec>(DELTA3);
  So.message_screen("Maximum ð³  =", max_C1);
  min_C1=get_min<real_prec>(DELTA3);
  So.message_screen("Minimum ð³ =", min_C1);
#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif
#endif
  // Here we do not need the eigenvalues. Hence we do not resize these vectors
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing gravitational potential");
#endif
  PoissonSolver(this->params._Lbox(), this->params._Nft(),delta,this->potential);
  So.DONE();
#ifdef _USE_S2_
  So.message_screen("Computing s²");
  this->S2.resize(this->params._NGRID(),0);
  So.DONE();
#endif
#ifdef _USE_S3_
  So.message_screen("Computing s³");
  this->S3.resize(this->params._NGRID(),0);
  So.DONE();
#endif
#ifdef _USE_NABLA2DELTA_
  So.message_screen("Computing NABLA2DELTA term");
  this->N2D.resize(this->params._NGRID(),0);
#endif
#if defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_)
  EigenValuesTweb_bias(this->params._Nft(),this->params._Lbox(),delta, this->potential,this->S2, this->S3, this->N2D);
  So.DONE();
#endif
#ifdef _USE_S2_
#ifdef _WRITE_S2_
    File.write_array(this->params._Output_directory()+"S2_original", this->S2);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for (ULONG index=0;index<this->params._NGRID();++index)
   {
    real_prec x=this->S2[index];;
#ifdef _USE_EXPONENT_S2_
    int sign_x=1;
#ifdef _USE_SIGN_S2_
    sign_x=sign(x);
#endif
    this->S2[index]=sign_x*pow(abs(x), EXPONENT_S2);
#else
   this->S2[index]=x;
#endif
   }
#ifdef _MAP_TO_INTERVAL_S2_
real_prec xmin=get_min(this->S2);
real_prec xmax=get_max(this->S2);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for (ULONG index=0;index<this->params._NGRID();++index)
    this->S2[index]=(NEWMAX_S2-NEWMIN_S2)*(this->S2[index]-xmin)/(xmax-xmin)+NEWMIN_S2;
#endif
#ifdef _WRITE_S2_
if(0==this->step)
#ifdef _USE_EXPONENT_S2_
    File.write_array(this->params._Output_directory()+"S2", this->S2);
#endif
#endif
#endif
#ifdef _USE_S3_
#ifdef _WRITE_S3_
    this->File.write_array(this->params._Output_directory()+"S3_original", this->S3);
#endif

#pragma omp parallel for
    for (ULONG index=0;index<this->params._NGRID();++index)
   {
    real_prec x=this->S3[index];;
#ifdef _USE_EXPONENT_S3_
    int sign_x=1;
#ifdef _USE_SIGN_S3_
    sign_x=sign(x);
#endif
    this->S3[index]=sign_x*pow(abs(x), EXPONENT_S3);
#else
   this->S3[index]=x;
#endif
   }


#ifdef _MAP_TO_INTERVAL_S3_
    {
real_prec xmin=get_min(this->S3);
real_prec xmax=get_max(this->S3);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for (ULONG index=0;index<this->params._NGRID();++index)
    this->S3[index]=(NEWMAX_S3-NEWMIN_S3)*(this->S3[index]-xmin)/(xmax-xmin)+NEWMIN_S3;

#endif
    }
#ifdef _WRITE_S3_
if(0==this->step)
#ifdef _USE_EXPONENT_S3_
    File.write_array(this->params._Output_directory()+"S3", this->S3);
#endif
#endif
#endif
#ifdef _USE_S2DELTA_
  So.message_screen("Computing S²DELTA");
  this->S2DELTA.resize(this->params._NGRID(),0);
#pragma omp parallel for
for (ULONG index=0;index<this->params._NGRID();++index)
     this->S2DELTA[index]=this->S2[index]*delta[index];
#ifdef _WRITE_S2DELTA_
    File.write_array(this->params._Output_directory()+"S2DELTA_original", this->S2DELTA);
#endif
#pragma omp parallel for
    for (ULONG index=0;index<this->params._NGRID();++index)
     {
    real_prec x=this->S2DELTA[index];
#ifdef _USE_EXPONENT_S2DELTA_
    int sign_x=1;
#ifdef _USE_SIGN_S2DELTA_
    sign_x=sign(x);
#endif
    this->S2DELTA[index]=sign_x*pow(abs(x), EXPONENT_S2DELTA);
#else
   this->S2DELTA[index]=x;
#endif
   }
#ifdef _MAP_TO_INTERVAL_S2DELTA_
    {
    real_prec xmin=get_min(this->S2DELTA);
    real_prec xmax=get_max(this->S2DELTA);
#pragma omp parallel for
    for (ULONG index=0;index<this->params._NGRID();++index)
        this->S2DELTA[index]=(NEWMAX_S2DELTA-NEWMIN_S2DELTA)*(this->S2DELTA[index]-xmin)/(xmax-xmin)+NEWMIN_S2DELTA;
    }
#endif
#if defined (_WRITE_S2DELTA_ ) && defined (_USE_EXPONENT_S2DELTA_)
    if(0==this->step)
       File.write_array(this->params._Output_directory()+"S2DELTA", this->S2DELTA);
#endif
#endif
#ifdef _USE_NABLA2DELTA_
#ifdef _WRITE_NABLA2DELTA_
  File.write_array(this->params._Output_directory()+"NABLA2DELTA", N2D);
#endif
  /*
    So.message_screen("Nabla²ð:");
    max_C1=get_max<real_prec>(N2D);
    So.message_screen("Maximum NABLA²ð  =", max_C1);
    min_C1=get_min<real_prec>(N2D);
    So.message_screen("Minimum Nabla²ð =", min_C1);
    cout<<endl;
  */
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// This function gnerates the CWC, the Iweb  sor the Pweb eigenvalues
void Cwclass::get_CWC(vector<real_prec>&delta)
{
#ifdef _FULL_VERBOSE_
  So.enter(__PRETTY_FUNCTION__);
#endif
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  this->potential.resize(this->params._NGRID(),0);

#ifdef _USE_ZERO_PADDING_POT_
  So.message_screen("using zero-padding");
#endif
#if defined(_USE_PWEB_)  || defined(_USE_CURVATURE_)
  this->potential=delta;
#else

#ifndef _USE_GFFT_EIGENV_
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing gravitational potential");
#endif
// If the Tidal field is to be computed using FFT's, we pass delta directly and do not need to solve Poisson
  PoissonSolver(this->params._Lbox(), this->params._Nft(),delta,this->potential);
#endif
#endif
  So.DONE();
#ifdef _FULL_VERBOSE_
#if defined (_USE_PWEB_)  || defined (_USE_CURVATURE_)
  So.message_screen("Computing Eigenvalues of zeta field");
#else
  So.message_screen("Computing Eigenvalues of tidal field");
#endif
#endif
  this->lambda1.resize(this->params._NGRID(),0);
  this->lambda2.resize(this->params._NGRID(),0);
  this->lambda3.resize(this->params._NGRID(),0);
#ifdef _USE_GFFT_EIGENV_
  // If the Tidal field is to be computed using FFT's, we pass delta directly and do not need to solve Poisson
  EigenValuesTweb(this->params._Nft(),this->params._Lbox(),delta,this->lambda1,this->lambda2,this->lambda3);//or Pweb, in case this->potential is delta
#else
  EigenValuesTweb(this->params._Nft(),this->params._Lbox(),delta, this->potential,this->lambda1,this->lambda2,this->lambda3);//or Pweb, in case this->potential is delta
#endif



#ifdef _FULL_VERBOSE_
  So.DONE();
#endif

  /*
   real_prec max_C1;
   real_prec min_C1;
   max_C1=get_max<real_prec>(lambda1);
   So.message_screen("Maximum lambda1 =", max_C1);
   max_C1=get_min<real_prec>(lambda1);
   So.message_screen("Minimum lambda1 =", max_C1);
   max_C1=get_max<real_prec>(lambda2);
   So.message_screen("Maximum lambda2 =", max_C1);
   max_C1=get_min<real_prec>(lambda2);
   So.message_screen("Minimum lambda2 =", max_C1);
   max_C1=get_max<real_prec>(lambda3);
   So.message_screen("Maximum lambda3 =", max_C1);
   max_C1=get_min<real_prec>(lambda3);
   So.message_screen("Minimum lambda3 =", max_C1);
*/
// ---------------------------------------------------------------------------------
#if defined _USE_INVARIANT_TIDAL_FIELD_I_  || defined(_USE_INVARIANT_PWEB_I_)
#ifdef _USE_INVARIANT_TIDAL_FIELD_I_
  So.message_screen("Computing Invariant Tidal field I");
#else 
 So.message_screen("Computing Invariant Derivative field P1");
#endif
 this->Invariant_TF_I.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_TF_I[index]=invariant_field_I(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _WRITE_INVARIANT_TIDAL_FIELD_I_
    File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I_original", this->Invariant_TF_I);
#elif defined(_WRITE_INVARIANT_PWEB_I_)
    File.write_array(this->params._Output_directory()+"INVARIANT_PI_original", this->Invariant_TF_I);
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_I_) || defined (_USE_INVARIANT_PWEB_I_)  // the transformations so far only applies to the tidal field invariants
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
     {
      real_prec x= this->Invariant_TF_I[index];
#if defined (_USE_EXPONENT_INVARIANT_I_) || defined(_USE_EXPONENT__INVARIANT_PWEB_I_)
      int sign_x=1;
#if defined (_USE_SIGN_INVARIANT_I_) || defined(_USE_SIGN_INVARIANT_PWEB_I_)
      sign_x=sign(x);
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_I_      
      this->Invariant_TF_I[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_I);
#elif defined(_USE_INVARIANT_PWEB_I_ )
      this->Invariant_TF_I[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_PWEB_I);
#endif
#else
     this->Invariant_TF_I[index]=x;
#endif
     }

#if defined(_MAP_TO_INTERVAL_INV_I_) || defined(_MAP_TO_INTERVAL_INVARIANT_PWEB_I_)
  {
  real_prec xmin=get_min(this->Invariant_TF_I);
  real_prec xmax=get_max(this->Invariant_TF_I);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
#if defined(_MAP_TO_INTERVAL_INV_I_)
      this->Invariant_TF_I[index]=(NEWMAX_INV_I-NEWMIN_INV_I)*(this->Invariant_TF_I[index]-xmin)/(xmax-xmin)+NEWMIN_INV_I;
#elif defined (_MAP_TO_INTERVAL_INVARIANT_PWEB_I_)
    this->Invariant_TF_I[index]=(NEWMAX_INVARIANT_PWEB_I-NEWMIN_INVARIANT_PWEB_I)*(this->Invariant_TF_I[index]-xmin)/(xmax-xmin)+NEWMIN_INVARIANT_PWEB_I;
#endif    

}
#endif

#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EXPONENT_INVARIANT_I_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I", this->Invariant_TF_I);
#elif defined _USE_EXPONENT_INVARIANT_PWEB_I_
  File.write_array(this->params._Output_directory()+"INVARIANT_PI", this->Invariant_TF_I);
#endif

   if(this->step ==this->params._N_iterations_Kernel())
#if defined (_USE_EXPONENT_INVARIANT_I_) || defined(_USE_EXPONENT_PWEB_I_)
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_I_iteration"+to_string(this->step), this->Invariant_TF_I);
#endif

#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif
#endif
#endif
#endif

#if defined _USE_INVARIANT_TIDAL_FIELD_II_  || defined(_USE_INVARIANT_PWEB_II_)

#ifdef _FULL_VERBOSE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  So.message_screen("Computing Invariant Tidal field II");
#else
 So.message_screen("Computing Invariant Derivative field P2");
#endif
#endif
 
 this->Invariant_TF_II.resize(this->params._NGRID(),0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_TF_II[index]=invariant_field_II(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_I_original", this->Invariant_TF_II);
#else

#ifdef _WRITE_INVARIANT_TIDAL_FIELD_II_
    File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II_original", this->Invariant_TF_II);
#elif defined(_WRITE_INVARIANT_PWEB_II_)
    File.write_array(this->params._Output_directory()+"INVARIANT_PII_original", this->Invariant_TF_II);
#endif
#endif
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_ ) || defined (_USE_INVARIANT_PWEB_II_) || defined(_USE_EIGENVALUES_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
     {
      real_prec x=this->Invariant_TF_II[index];
#if defined (_USE_EXPONENT_INVARIANT_II_) || defined(_USE_EXPONENT_INVARIANT_PWEB_II_)
     int sign_x=1;
#if defined (_USE_SIGN_INVARIANT_II_) || defined(_USE_SIGN_INVARIANT_PWEB_II_)
      sign_x=sign(x);
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_II_      
      this->Invariant_TF_II[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_II);
#elif defined(_USE_INVARIANT_PWEB_II_ )
      this->Invariant_TF_II[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_PWEB_II);
#endif
#else
     this->Invariant_TF_II[index]=x;
#endif
     }

#if defined(_MAP_TO_INTERVAL_INV_II_) || defined(_MAP_TO_INTERVAL_INVARIANT_PWEB_II_)
 {
  real_prec xmin=get_min(this->Invariant_TF_II);
  real_prec xmax=get_max(this->Invariant_TF_II);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
#if defined(_MAP_TO_INTERVAL_INV_II_)
      this->Invariant_TF_II[index]=(NEWMAX_INV_II-NEWMIN_INV_II)*(this->Invariant_TF_II[index]-xmin)/(xmax-xmin)+NEWMIN_INV_II;
#elif defined (_MAP_TO_INTERVAL_INVARIANT_PWEB_II_)
    this->Invariant_TF_II[index]=(NEWMAX_INVARIANT_PWEB_II-NEWMIN_INVARIANT_PWEB_II)*(this->Invariant_TF_II[index]-xmin)/(xmax-xmin)+NEWMIN_INVARIANT_PWEB_II;
#endif

  }
#endif

#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#if defined (_USE_EXPONENT_INVARIANT_II_) || defined (_USE_EXPONENT_INVARIANT_PWEB_II_) || defined (_USE_EIGENVALUES_)
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_I", this->Invariant_TF_II);
#elif defined _USE_EXPONENT_INVARIANT_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II", this->Invariant_TF_II);
#elif defined _USE_EXPONENT_INVARIANT_PWEB_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_PII", this->Invariant_TF_II);
#endif

#endif

   if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_II_
     File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_II_iteration"+to_string(this->step), this->Invariant_TF_II);
#elif defined ( _USE_EXPONENT_INVARIANT_PWEB_II_)
   File.write_array(this->params._Output_directory()+"INVARIANT_PWEB_II_iteration"+to_string(this->step), this->Invariant_TF_II);
#endif

#ifdef _FULL_VERBOSE_
  cout<<endl;
#endif

#endif
  So.DONE();
#endif
#endif
#ifdef _HYDROTEST_
  this->File.read_array_t<PrecType_X>("/home/andres/data/Numerics/IACmocks/ANALYSIS/BAM/Output_Hydro/DensityGas_nohead.37.n128.dat",this->Invariant_TF_II);
  get_overdens(this->Invariant_TF_II, this->Invariant_TF_II);
  // if(this->step>0)
  //   this->Konvolve(Invariant_TF_II, Invariant_TF_II, "DELTA");
#pragma omp parallel for
  for(ULONG i = 0;i < this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->Invariant_TF_II[i] = this->Invariant_TF_II[i]<-1 ?  0 :  log10(NUM_IN_LOG+ static_cast<real_prec>(this->Invariant_TF_II[i]));
#endif

// ---------------------------------------------------------------------------------
#if defined _USE_INVARIANT_TIDAL_FIELD_III_  || defined(_USE_INVARIANT_PWEB_III_)

#ifdef _FULL_VERBOSE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  So.message_screen("Computing Invariant Tidal field III");
#else
 So.message_screen("Computing Invariant Derivative field P3");
#endif
#endif


  this->Invariant_TF_III.resize(this->params._NGRID(),0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_TF_III[index]=invariant_field_III(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_II_original", this->Invariant_TF_III);
#else

#ifdef _WRITE_INVARIANT_TIDAL_FIELD_III_
    File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III_original", this->Invariant_TF_III);
#elif defined(_WRITE_INVARIANT_PWEB_III_)
    File.write_array(this->params._Output_directory()+"INVARIANT_PIII_original", this->Invariant_TF_III);
#endif

  
#endif
#endif


#if defined (_USE_INVARIANT_TIDAL_FIELD_II_ )  || defined (_USE_INVARIANT_PWEB_III_) || defined(_USE_EIGENVALUES_) // TRANSOFMRATION ONLY AVAILABLE FOR TIDALFIELD EIGENVALUES

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     real_prec x=this->Invariant_TF_III[index];
#if defined (_USE_EXPONENT_INVARIANT_III_) || defined(_USE_EXPONENT_INVARIANT_PWEB_III_)
     int sign_x=1;
#if defined (_USE_SIGN_INVARIANT_III_) || defined(_USE_SIGN_INVARIANT_PWEB_III_)
     sign_x=sign(x);
#endif

#ifdef _USE_INVARIANT_TIDAL_FIELD_III_      
      this->Invariant_TF_III[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_III);
#elif defined(_USE_INVARIANT_PWEB_III_ )
      this->Invariant_TF_III[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_PWEB_III);
#endif
#else
     this->Invariant_TF_III[index]=x;
#endif
  }

#if defined(_MAP_TO_INTERVAL_INV_III_) || defined(_MAP_TO_INTERVAL_INVARIANT_PWEB_III_)
  {
    real_prec xmin=get_min(this->Invariant_TF_III);
    real_prec xmax=get_max(this->Invariant_TF_III);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG index=0;index<this->params._NGRID();++index)
#if defined(_MAP_TO_INTERVAL_INV_III_)
      this->Invariant_TF_III[index]=(NEWMAX_INV_III-NEWMIN_INV_III)*(this->Invariant_TF_III[index]-xmin)/(xmax-xmin)+NEWMIN_INV_III;
#elif defined (_MAP_TO_INTERVAL_INVARIANT_PWEB_III_)
    this->Invariant_TF_III[index]=(NEWMAX_INVARIANT_PWEB_III-NEWMIN_INVARIANT_PWEB_III)*(this->Invariant_TF_III[index]-xmin)/(xmax-xmin)+NEWMIN_INVARIANT_PWEB_III;
#endif
  }
#endif
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#if defined (_USE_EXPONENT_INVARIANT_III_) || defined (_USE_EXPONENT_INVARIANT_PWEB_III_) || defined (_USE_EIGENVALUES_)
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_II", this->Invariant_TF_III);
#elif defined _USE_EXPONENT_INVARIANT_III_
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III", this->Invariant_TF_III);
#elif defined _USE_EXPONENT_INVARIANT_PWEB_III_
      File.write_array(this->params._Output_directory()+"INVARIANT_PIII", this->Invariant_TF_III);
#endif

#endif
#endif

      if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_III_
	File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_III_iteration"+to_string(this->step), this->Invariant_TF_III);
#elif defined ( _USE_EXPONENT_INVARIANT_PWEB_III_)
      File.write_array(this->params._Output_directory()+"INVARIANT_PIII_iteration"+to_string(this->step), this->Invariant_TF_III);
#endif
    So.DONE();
#endif
#endif

#ifdef _HYDROTEST_
    this->File.read_array_t<PrecType_X>("/home/andres/data/Numerics/IACmocks/ANALYSIS/BAM/Output_Hydro/NumberDensityHI_nohead.37.n128.dat",this->Invariant_TF_III);
  get_overdens(this->Invariant_TF_III, this->Invariant_TF_III);
// if(this->step>0)
  //   this->Konvolve(Invariant_TF_III, Invariant_TF_III, "DELTA");
#pragma omp parallel for
  for(ULONG i = 0;i < this->params._NGRID() ;++i) //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->Invariant_TF_III[i] = this->Invariant_TF_III[i]<-1 ?  0 :  log10(NUM_IN_LOG+ static_cast<real_prec>(this->Invariant_TF_III[i]));
#endif
// ---------------------------------------------------------------------------------
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing Invariant Tidal field IV");
#endif
  this->Invariant_TF_IV.resize(this->params._NGRID(),0);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_TF_IV[index]=invariant_field_IV(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_III_original", this->Invariant_TF_IV);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_IV_original", this->Invariant_TF_IV);
#endif
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     real_prec x=Invariant_TF_IV[index];
#ifdef _USE_EXPONENT_INVARIANT_IV_
     int sign_x=1;
#ifdef _USE_SIGN_INVARIANT_IV_
     sign_x=sign(x);
#endif
     this->Invariant_TF_IV[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_IV);
#else
     this->Invariant_TF_IV[index]=x;
#endif
  }
#ifdef _MAP_TO_INTERVAL_INV_IV_
  {
  real_prec xmin=get_min(this->Invariant_TF_IV);
  real_prec xmax=get_max(this->Invariant_TF_IV);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_TF_IV[index]=(NEWMAX_INV_IV-NEWMIN_INV_IV)*(this->Invariant_TF_IV[index]-xmin)/(xmax-xmin)+NEWMIN_INV_IV;
}
#endif
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EXPONENT_INVARIANT_IV_
#ifdef _USE_EIGENVALUES_
      File.write_array(this->params._Output_directory()+"LAMBDA_III", this->Invariant_TF_IV);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_IV", this->Invariant_TF_IV);
#endif
  if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_INVARIANT_IV_
     File.write_array(this->params._Output_directory()+"INVARIANT_TIDAL_IV_iteration"+to_string(this->step), this->Invariant_TF_IV);
#endif
#endif
#endif
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_TIDAL_ANISOTROPY_
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing Tidal Anisotropy");
#endif
  this->Tidal_Anisotropy.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Tidal_Anisotropy[index]=tidal_anisotropy(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
      File.write_array(this->params._Output_directory()+"TIDAL_ANISOTROPY_original", this->Tidal_Anisotropy);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     real_prec x=this->Tidal_Anisotropy[index];
#ifdef _USE_EXPONENT_TIDAL_ANISOTROPY_
     int sign_x=1;
#ifdef _USE_SIGN_TIDAL_ANISOTROPY_
     sign_x=sign(x);
#endif
     this->Tidal_Anisotropy[index]=sign_x*pow(abs(x), EXPONENT_TIDAL_ANISOTROPY);
#else
     this->Tidal_Anisotropy[index]=x;
#endif
  }
#ifdef _MAP_TO_INTERVAL_TIDAL_ANISOTROPY_
  {
  real_prec xmin=get_min(this->Tidal_Anisotropy);
  real_prec xmax=get_max(this->Tidal_Anisotropy);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Tidal_Anisotropy[index]=(NEWMAX_TIDAL_ANISOTROPY-NEWMIN_TIDAL_ANISOTROPY)*(this->Tidal_Anisotropy[index]-xmin)/(xmax-xmin)+NEWMIN_TIDAL_ANISOTROPY;
}
#endif
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EXPONENT_TIDAL_ANISOTROPY_
      File.write_array(this->params._Output_directory()+"TIDAL_ANISOTROPY", this->Tidal_Anisotropy);
#endif
  if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_EXPONENT_USE_TIDAL_ANISOTROPY_
     File.write_array(this->params._Output_directory()+"TIDAL_ANISOTROPY_iteration"+to_string(this->step), this->Tidal_Anisotropy);
#endif
#endif
    So.DONE();
#endif
   // *********************************************************************************
  #ifdef _USE_ELLIPTICITY_
    So.message_screen("Computing Ellipticity");
    this->Ellipticity.resize(this->params._NGRID(),0);
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
    for (ULONG index=0;index<this->params._NGRID();++index)
        this->Ellipticity[index]=ellipticity(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
  #ifdef _WRITE_INVARIANTS_
    if(0==this->step)
        File.write_array(this->params._Output_directory()+"ELLIPTICITY_original", this->Ellipticity);
  #endif
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
    for (ULONG index=0;index<this->params._NGRID();++index)
     {
       real_prec x=this->Ellipticity[index];
  #ifdef _USE_EXPONENT_ELLIPTICITY_
       int sign_x=1;
  #ifdef _USE_SIGN_ELLIPTICITY_
       sign_x=sign(x);
  #endif
       this->Ellipticity[index]=sign_x*pow(abs(x), EXPONENT_ELLIPTICITY);
  #else
       this->Ellipticity[index]=x;
  #endif
    }
  #ifdef _MAP_TO_INTERVAL_ELLIPTICITY_
    {
    real_prec xmin=get_min(this->Ellipticity);
    real_prec xmax=get_max(this->Ellipticity);
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
    for (ULONG index=0;index<this->params._NGRID();++index)
        this->Ellipticity[index]=(NEWMAX_ELLIPTICITY-NEWMIN_ELLIPTICITY)*(this->Ellipticity[index]-xmin)/(xmax-xmin)+NEWMIN_ELLIPTICITY;
  }
  #endif
  #ifdef _WRITE_INVARIANTS_
    if(0==this->step)
  #ifdef _USE_EXPONENT_ELLIPTICITY_
        File.write_array(this->params._Output_directory()+"ELLIPTICITY", this->Ellipticity);
  #endif
    if(this->step ==this->params._N_iterations_Kernel())
  #ifdef _USE_ELLIPTICITY_
       File.write_array(this->params._Output_directory()+"ELLIPTICITY_iteration"+to_string(this->step), this->Ellipticity);
  #endif
  #endif
      So.DONE();
  #endif
    // ************************************************************************************
#ifdef _USE_PROLATNESS_
  So.message_screen("Computing Prolatness");
  this->Prolatness.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Prolatness[index]=prolatness(this->lambda1[index],this->lambda2[index],this->lambda3[index]);
#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
      File.write_array(this->params._Output_directory()+"PROLATNESS_original", this->Prolatness);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     real_prec x=this->Prolatness[index];
#ifdef _USE_EXPONENT_PROLATNESS_
     int sign_x=1;
#ifdef _USE_SIGN_PROLATNESS_
     sign_x=sign(x);
#endif
     this->Prolatness[index]=sign_x*pow(abs(x), EXPONENT_PROLATNESS);
#else
     this->Prolatness[index]=x;
#endif
  }
#ifdef _MAP_TO_INTERVAL_PROLATNESS_
  {
  real_prec xmin=get_min(this->Prolatness);
  real_prec xmax=get_max(this->Prolatness);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Prolatness[index]=(NEWMAX_PROLATNESS-NEWMIN_PROLATNESS)*(this->Prolatness[index]-xmin)/(xmax-xmin)+NEWMIN_PROLATNESS;
}
#endif

#ifdef _WRITE_INVARIANTS_
  if(0==this->step)
#ifdef _USE_EXPONENT_PROLATNESS_
      File.write_array(this->params._Output_directory()+"PROLATNESS", this->Prolatness);
#endif

  if(this->step ==this->params._N_iterations_Kernel())
#ifdef _USE_PROLATNESS_
     File.write_array(this->params._Output_directory()+"PROLATNESS_iteration"+to_string(this->step), this->Prolatness);
#endif
#endif
    So.DONE();
#endif
  this->CWClass.clear();
#if defined (_USE_CWC_) || defined (_USE_MASS_KNOTS_)
  this->CWClass.resize(this->params._NGRID(),0);
#endif

#ifdef _USE_CWC_
  ULONG nknots=0;
  ULONG nfilaments=0;
  ULONG nsheets=0;
  ULONG nvoids=0;
  ULONG nrest=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Extracting the Cosmic Web Classification:");
  So.message_screen("Lambda threshold = ", this->params._lambdath());
#endif

  /*

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nknots,nfilaments,nsheets,nvoids,nrest)
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
    {

          if (this->lambda1[index]>=this->params._lambdath() && this->lambda2[index]>=this->params._lambdath() && this->lambda3[index]>=this->params._lambdath())
	{
	  this->CWClass[index]=I_KNOT;
	  nknots++;
	}
      
      else if (this->lambda1[index]<this->params._lambdath() && lambda2[index]<this->params._lambdath() && lambda3[index]<this->params._lambdath())
	{
	  this->CWClass[index]=I_VOID;
	  nvoids++;
	}
      
      else if ((lambda1[index]<this->params._lambdath() && lambda2[index]<this->params._lambdath() && lambda3[index]>=this->params._lambdath()) || (lambda1[index]<this->params._lambdath() && lambda2[index]>=this->params._lambdath() && lambda3[index]<this->params._lambdath()) || (lambda1[index]>=this->params._lambdath() && lambda2[index]<this->params._lambdath() && lambda3[index]<this->params._lambdath()))
	{
	  this->CWClass[index]=I_SHEET;
	  nsheets++;
	}
      
      else if ((lambda1[index]<this->params._lambdath() && lambda2[index]>=this->params._lambdath() && lambda3[index]>=this->params._lambdath()) || (lambda1[index]>=this->params._lambdath() && lambda2[index]>this->params._lambdath() && lambda3[index]<this->params._lambdath()) || (lambda1[index]>this->params._lambdath() && lambda2[index]<this->params._lambdath() && lambda3[index]>=this->params._lambdath()))
	{
	  this->CWClass[index]=I_FILAMENT;
	  nfilaments++;
	}
      else
	nrest ++;
    }
 */


  // This is a more compact way of performing the classifictation
 vector<int>vff(4,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
    {
    vector<real_prec>cwc_v={this->lambda1[index], this->lambda2[index],this->lambda3[index]};
    int vtype=0;
    for(int it=0;it<cwc_v.size();++it)
        if(cwc_v[it]>this->params._lambdath())
            vtype++;
    this->CWClass[index]=4-vtype;
#pragma omp atomic
    vff[4-vtype-1]++;
  }
  So.DONE();
  this->volume_knots=static_cast<real_prec>(vff[0])/static_cast<real_prec>(this->params._NGRID()); // Volume in knots / Total volume
  this->volume_filaments=static_cast<real_prec>(vff[1])/static_cast<real_prec>(this->params._NGRID()); // Volume in knots / Total volume
  this->volume_sheets=static_cast<real_prec>(vff[2])/static_cast<real_prec>(this->params._NGRID()); // Volume in knots / Total volume
  this->volume_voids=static_cast<real_prec>(vff[3])/static_cast<real_prec>(this->params._NGRID()); // Volume in knots / Total volume
/*

  string ss = this->params._Output_directory()+"delta_cwc.txt";
  ofstream ssa; ssa.open(ss.c_str());
  for (ULONG index=0;index<this->params._NGRID();++index)
    ssa<<delta[index]<<"\t"<<CWClass[index]<<endl;
    ssa.close();
//    exit(0);

*/
//  File.write_array(this->params._Output_directory()+"CWC",this->CWClass);

/*

  File.write_array(this->params._Output_directory()+"CWC",this->CWClass);
  vector<real_prec>l1;
  vector<real_prec>l2;
  vector<real_prec>l3;
  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     if(I_VOID==this->CWClass[index])
      {
         l1.push_back(this->Invariant_TF_I[index]);
         l2.push_back(this->Invariant_TF_II[index]);
         l3.push_back(this->Invariant_TF_III[index]);
      }
  }
  File.write_array(this->params._Output_directory()+"eigenvalues1_voids",l1);
  File.write_array(this->params._Output_directory()+"eigenvalues2_voids",l2);
  File.write_array(this->params._Output_directory()+"eigenvalues3_voids",l3);
  l1.clear();l1.shrink_to_fit();
  l2.clear();l2.shrink_to_fit();
  l3.clear();l3.shrink_to_fit();

  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     if(I_SHEET==this->CWClass[index])
      {
         l1.push_back(this->Invariant_TF_I[index]);
         l2.push_back(this->Invariant_TF_II[index]);
         l3.push_back(this->Invariant_TF_III[index]);
      }
  }
  File.write_array(this->params._Output_directory()+"eigenvalues1_sheets",l1);
  File.write_array(this->params._Output_directory()+"eigenvalues2_sheets",l2);
  File.write_array(this->params._Output_directory()+"eigenvalues3_sheets",l3);
  l1.clear();l1.shrink_to_fit();
  l2.clear();l2.shrink_to_fit();
  l3.clear();l3.shrink_to_fit();

  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     if(I_FILAMENT==this->CWClass[index])
      {
         l1.push_back(this->Invariant_TF_I[index]);
         l2.push_back(this->Invariant_TF_II[index]);
         l3.push_back(this->Invariant_TF_III[index]);
      }
  }
  File.write_array(this->params._Output_directory()+"eigenvalues1_filaments",l1);
  File.write_array(this->params._Output_directory()+"eigenvalues2_filaments",l2);
  File.write_array(this->params._Output_directory()+"eigenvalues3_filaments",l3);
  l1.clear();l1.shrink_to_fit();
  l2.clear();l2.shrink_to_fit();
  l3.clear();l3.shrink_to_fit();

  for (ULONG index=0;index<this->params._NGRID();++index)
   {
     if(I_KNOT==this->CWClass[index])
      {
         l1.push_back(this->Invariant_TF_I[index]);
         l2.push_back(this->Invariant_TF_II[index]);
         l3.push_back(this->Invariant_TF_III[index]);
      }
  }
  File.write_array(this->params._Output_directory()+"eigenvalues1_knots",l1);
  File.write_array(this->params._Output_directory()+"eigenvalues2_knots",l2);
  File.write_array(this->params._Output_directory()+"eigenvalues3_knots",l3);
  l1.clear();l1.shrink_to_fit();
  l2.clear();l2.shrink_to_fit();
  l3.clear();l3.shrink_to_fit();

  */

#ifdef _FULL_VERBOSE_
  So.message_screen("Summary of T-web classification");
  So.message_screen("Knots (%) =", volume_knots*100.0);
  So.message_screen("Filaments (%) =", volume_filaments*100.0);
  So.message_screen("Sheets (%) =", volume_sheets*100.);
  So.message_screen("Voids (%) =", volume_voids*100.0);
#endif
  this->knots_fraction=volume_knots*100.0;
  this->filaments_fraction=volume_sheets*100.0;
  this->sheets_fraction=volume_filaments*100.0;
  this->voids_fraction=volume_voids*100.0;



#ifdef _VERBOSE_
#ifndef _TEST_THRESHOLDS_RESIDUALS_
  int index= (this->step <=this->params._N_iterations_Kernel())  ?  this->step  : this->step - (this->params._N_iterations_Kernel())+1;
  string label_aux = this->step <= this->params._N_iterations_Kernel() ? "_iteration": "_realization";
  string file_cwc=this->params._Output_directory()+"CWC"+label_aux+to_string(index)+".txt";
  ofstream fcwc; fcwc.open(file_cwc.c_str());
  So.message_screen("Writing CWC fractions in file", file_cwc);
  fcwc<<static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc.close();
  So.DONE();
#endif //* end for ifndef _TEST_THRESHOLDS_RESIDUALS_ *//


#ifdef _FULL_VERBOSE_
  if(nrest>0)
    So.message_screen("Unclassified (%) =",static_cast<real_prec>(nrest)*100./static_cast<real_prec>(this->params._NGRID()));
#endif
#endif //* end for _VERBOSE_

  // If we do not want to use the CWC but still use the infor from the knots, then
#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
#ifdef _FULL_VERBOSE_
  So.message_screen("Doing Cosmic Web Classification. Knots Only. Lambda threshold = ", this->params._lambdath());
#endif
#pragma omp parallel for
  for (ULONG index=0;index<this->params._NGRID();++index)
    if (lambda1[index]>this->params._lambdath() && lambda2[index]>this->params._lambdath() && lambda3[index]>this->params._lambdath())
      this->CWClass[index]=I_KNOT;
  So.DONE();
#endif
#endif  //* end for  #ifdef _USE_CWC_

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_CWC_V(vector<real_prec>&Vx,vector<real_prec>&Vy,vector<real_prec>&Vz)
{
#ifdef _FULL_VERBOSE_
    So.enter(__PRETTY_FUNCTION__);
    So.message_screen("V-WEB classification requested");
#endif

    // Transform the velocities to v / H(z), putting the units of H(z) with the cvel, such that the shear is dimensionless
/*
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0; i<this->params._NGRID();++i)
    {
       Vx[i]/=cvel;
       Vy[i]/=cvel;
       Vz[i]/=cvel;
    }

*/
  this->lambda1_vs.resize(this->params._NGRID(),0);
  this->lambda2_vs.resize(this->params._NGRID(),0);
  this->lambda3_vs.resize(this->params._NGRID(),0);
  this->Divergence_VelField.resize(this->params._NGRID(),0);
  So.message_screen("Computing Divergence and Eigenvalues of the shear of the velocity field");
  EigenValuesVweb(this->params._Nft(),params._Lbox(),Vx,Vy,Vz,this->Divergence_VelField, this->lambda1_vs,this->lambda2_vs,this->lambda3_vs);
  So.DONE();


#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
    So.message_screen("Invariant Shear Vfield I");
    this->Invariant_VS_I.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
     {
      real_prec x=invariant_field_I(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_I_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_I[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_I);
     }

#ifdef _FULL_VERBOSE_
  cout<<endl;
 #endif


#ifdef _MAP_TO_INTERVAL_INV_SHEAR_I_
  {
  real_prec xmin=get_min(this->Invariant_VS_I);
  real_prec xmax=get_max(this->Invariant_VS_I);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_VS_I[index]=(NEWMAX_INV_SHEAR_I-NEWMIN_INV_SHEAR_I)*(this->Invariant_VS_I[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_I;
}
#endif
#ifdef _USE_EXPONENT_INVARIANT_VS_I_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_I", this->Invariant_VS_I);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_I_original", this->Invariant_VS_I);
#endif
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
    So.message_screen("Invariant Shear Vfield II");
    this->Invariant_VS_II.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
     {
      real_prec x=invariant_field_II(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_II_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_II[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_II);
     }
#ifdef _FULL_VERBOSE_
  cout<<endl;
 #endif


#ifdef _MAP_TO_INTERVAL_INV_SHEAR_II_
  {
  real_prec xmin=get_min(this->Invariant_VS_II);
  real_prec xmax=get_max(this->Invariant_VS_II);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_VS_II[index]=(NEWMAX_INV_SHEAR_II-NEWMIN_INV_SHEAR_II)*(this->Invariant_VS_II[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_II;
}
#endif

#ifdef _USE_EXPONENT_INVARIANT_VS_II_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_II", this->Invariant_VS_II);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_II_original", this->Invariant_VS_II);
#endif

#endif


#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
    So.message_screen("Invariant Shear Vfield III");
    this->Invariant_VS_III.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
     {
      real_prec x=invariant_field_III(this->lambda1_vs[index],this->lambda2_vs[index],this->lambda3_vs[index]);
      int sign_x=1;
 #ifdef _USE_SIGN_INVARIANT_VS_III_
      sign_x=sign(x);
 #endif
      this->Invariant_VS_III[index]=sign_x*pow(abs(x), EXPONENT_INVARIANT_VS_III);
     }
#ifdef _FULL_VERBOSE_
  cout<<endl;
 #endif

#ifdef _MAP_TO_INTERVAL_INV_SHEAR_III_
  {
  real_prec xmin=get_min(this->Invariant_VS_III);
  real_prec xmax=get_max(this->Invariant_VS_III);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
      this->Invariant_VS_III[index]=(NEWMAX_INV_SHEAR_III-NEWMIN_INV_SHEAR_III)*(this->Invariant_VS_III[index]-xmin)/(xmax-xmin)+NEWMIN_INV_SHEAR_III;
}
#endif

  #ifdef _USE_EXPONENT_INVARIANT_VS_III_
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_III", this->Invariant_VS_III);
#else
      File.write_array(this->params._Output_directory()+"INVARIANT_SHEAR_III_original", this->Invariant_VS_III);
#endif
#endif
#ifdef _USE_CWC_V_
  this->CWClass_V.clear();
  this->CWClass_V.resize(this->params._NGRID(),0);
  ULONG nknots=0;
  ULONG nfilaments=0;
  ULONG nsheets=0;
  ULONG nvoids=0;
  ULONG nrest=0;
  So.message_screen("Extracting the Cosmic Web Classification:");
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nknots,nfilaments,nsheets,nvoids,nrest)
#endif
  for (ULONG index=0;index<this->params._NGRID();++index)
    {
      if (this->lambda1_vs[index]>this->params._lambdath()_v && this->lambda2_vs[index]>this->params._lambdath()_v && this->lambda3_vs[index]>this->params._lambdath()_v)
        {
          this->CWClass_V[index]=I_KNOT;
          nknots++;
        }

      else if (this->lambda1_vs[index]<this->params._lambdath()_v && lambda2_vs[index]<this->params._lambdath()_v && lambda3_vs[index]<this->params._lambdath()_v)
        {
          this->CWClass_V[index]=I_VOID;
          nvoids++;
        }

      else if ((lambda1_vs[index]<this->params._lambdath()_v && lambda2_vs[index]<this->params._lambdath()_v && lambda3_vs[index]>this->params._lambdath()_v) || (lambda1_vs[index]<this->params._lambdath()_v && lambda2_vs[index]>this->params._lambdath()_v && lambda3_vs[index]<this->params._lambdath()_v) || (lambda1_vs[index]>this->params._lambdath()_v && lambda2_vs[index]<this->params._lambdath()_v && lambda3_vs[index]<this->params._lambdath()_v))
        {
          this->CWClass_V[index]=I_SHEET;
          nsheets++;
        }

      else if ((lambda1_vs[index]<this->params._lambdath()_v && lambda2_vs[index]>this->params._lambdath()_v && lambda3_vs[index]>this->params._lambdath()_v) || (lambda1_vs[index]>this->params._lambdath()_v && lambda2_vs[index]>this->params._lambdath()_v && lambda3_vs[index]<this->params._lambdath()_v) || (lambda1_vs[index]>this->params._lambdath()_v && lambda2_vs[index]<this->params._lambdath()_v && lambda3_vs[index]>this->params._lambdath()_v))
        {
          this->CWClass_V[index]=I_FILAMENT;
          nfilaments++;
        }
      else
        nrest ++;
    }
  So.DONE();
#ifdef _VERBOSE_
  So.message_screen("Summary of V-web classification, Lambda threshold = ", this->params._lambdath()_v);
  So.message_screen("V_Knots (%) =",static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->params._NGRID()));
  So.message_screen("V_Filaments (%) =",static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->params._NGRID()));
  So.message_screen("V_Sheets (%) =",static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->params._NGRID()));
  So.message_screen("V_Voids (%) =",static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->params._NGRID()));
#endif


#ifndef _TEST_THRESHOLDS_RESIDUALS_
  int index= (this->step <=this->params._N_iterations_Kernel())  ?  this->step  : this->step - (this->params._N_iterations_Kernel())+1;
  string label_aux = this->step <= this->params._N_iterations_Kernel() ? "_iteration": "_realization";

  string file_cwc=this->params._Output_directory()+"CWC_V"+label_aux+to_string(index)+".txt";
  ofstream fcwc; fcwc.open(file_cwc.c_str());
  So.message_screen("Writing CWC-V fractions in file", file_cwc);
  fcwc<<static_cast<real_prec>(nknots)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nfilaments)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nsheets)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc<<static_cast<real_prec>(nvoids)*100./static_cast<real_prec>(this->params._NGRID())<<endl;
  fcwc.close();
  So.DONE();
#endif

  if(nrest>0)
    So.message_screen("Unclassified (%) =",static_cast<real_prec>(nrest)*100./static_cast<real_prec>(this->params._NGRID()));

  // If we do not want to use the CWC but still use the infor from the knots, then
#elif !defined _USE_CWC_
#ifdef _USE_MASS_KNOTS_
  So.message_screen("Doing Cosmic Web Classification. Knots Only. Lambda threshold = ", this->params._lambdath());
#pragma omp parallel for
  for (ULONG index=0;index<this->params._NGRID();++index)
    if (lambda1_vs[index]>this->params._lambdath() && lambda2_vs[index]>this->params._lambdath() && lambda3_vs[index]>this->params._lambdath())
      this->CWClass_V[index]=I_KNOT;
  So.DONE();
#endif
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_Mk_collapsing_regions(vector<real_prec>&in, real_prec Nmean)
{
  /*
    Routine taken from HADRON and adapted to BAM by A. Balaguera
  */
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

    // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,
#ifdef _FULL_VERBOSE_
  cout<<endl;
  So.message_screen("Getting FoF from Knots");
  So.message_screen("mean number density =", Nmean);
#endif
  vector<real_prec>rho(in.size(),0);

  if(Nmean>0)  // If Nob is zero, then we are passing a density field. If not, convert delta to density with the value of Nobs passed
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._NGRID();++i)
        rho[i]= Nmean*(num_1+in[i]);
    }
  else
    rho=in;
  

    ULONG aux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:aux)
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
      if(this->CWClass[i]==I_KNOT)
	aux++;

  // Initialize to -1 to then identify whether a cell it is used or not
  // This array will be finally used to allocate the bin in wehich a mass knot is found.
  this->SKNOT_M_info.clear();
  this->SKNOT_M_info.shrink_to_fit();
  this->SKNOT_M_info.resize(this->params._NGRID(), -1);
  
  // This vector contains the number of grids that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(KNOT_MAX_SIZE,0);  
  
  int k;
  real_prec knotmass;
  real_prec MKrange_max = -1e38;  // Maximum  and Min of Mass of DM in knots. Mass here is number of DM particles in the super knot
  real_prec MKrange_min = +1e38; // These two variables are updated when knot masses are computed
  
  for(ULONG i=0;i < this->params._NGRID(); ++i)
    {
      if(this->CWClass[i]==I_KNOT && rho[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
	{
	  if(this->SKNOT_M_info[i] == -1) // if this cell has not been already used in the fof, proceed
	    {
	      
	      knot[0]=i;            // The first cell classified as knot. knot is always inizialized like this for every cell
	      this->SKNOT_M_info[i]=1;    //Mark the first knot as visited
	      int n=1; // Starts from 1. It gets added 1 if one of the neighbour cells is a knot. Percolation.
	      
              for(int j=0;j<n;++j)  // j runs over the neighbour cells of the current one, up to n, when n is increased inside the "if" when a neighbour knot is found
		{
		  if( (k= knot[j]-1)>=0 && k < this->params._NGRID())
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
		      {
			knot[n] = k;   
			this->SKNOT_M_info[k] = 1; //mark this cell as already visited
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + 1) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] - this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE) 
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if((this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE) 
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] - this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		  
		  if( (k = knot[j] + this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if((this->CWClass[k]==I_KNOT && rho[k]>0) && this->SKNOT_M_info[k] == -1)
		      {
			knot[n] = k;
			this->SKNOT_M_info[k] = 1;
			n++;
			if(n >= KNOT_MAX_SIZE)
			  So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
		      }
		}
	      
	      knotmass = 0;
              for(int j = 0; j < n; j++) // Compute the mass of the SK by adding the mass of its n parts. 
		knotmass += rho[knot[j]];

	      for(int j = 0; j < n; j++) //Assign to all cells in this superknot the mass of the sp they are in.
		this->SKNOT_M_info[knot[j]] = knotmass;
	      
	      // Adjust the interval for the bins in KMass to the current values min and max of the variable knotmass
	      if(MKrange_min > knotmass) MKrange_min = knotmass;
	      if(MKrange_max < knotmass) MKrange_max = knotmass;
	    }
	}
      else
	this->SKNOT_M_info[i] = 0;
    }

 knot.clear();
 knot.shrink_to_fit();


#ifdef _WRITE_MKNOTS_
 if(this->step==0)
  File.write_array(this->params._Output_directory()+"MASS_KNOTS", this->SKNOT_M_info);
#endif

 real_prec number_log;
#ifdef _DM_NEW_UNITS_
 number_log=0.0;
#else
 number_log=1.0;
#endif

#ifdef _FULL_VERBOSE_
 So.message_screen("Min and Max from FoF: ");
 So.message_screen("log Minimum Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_min+number_log));
 So.message_screen("log Maximim Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_max+number_log));
#endif

    // *** TYhis region is commnted if we want the current limits to be used in each step of the iteration
#ifndef _MODIFY_LIMITS_
  MKrange_max = MKMAX;
  MKrange_min = MKMIN;
#ifdef _FULL_VERBOSE_
  So.message_screen("log Used minimum Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_min+number_log));
  So.message_screen("log Used maximim Knot-mass (in units of log_10 Mass of DM particle) =", log10(MKrange_max+number_log));
  So.message_screen("");
#endif
#endif


  real_prec deltaMK= static_cast<real_prec>((log10(MKrange_max+number_log)-log10(MKrange_min+number_log))/static_cast<real_prec>(this->params._n_sknot_massbin()));

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0; i < this->params._NGRID(); ++i)
    if(this->SKNOT_M_info[i] == -1) // If this cell was not used, exit
      this->SKNOT_M_info[i] = 0;
    else
     // Get the bin in Knot mass. At this point SKNOT_M_info[i] encodes the mass of the Super knot where this cell is included
     this->SKNOT_M_info[i] = get_bin(log10(SKNOT_M_info[i]+number_log),log10(MKrange_min+number_log), this->params._n_sknot_massbin(), deltaMK,true);
  So.DONE();
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_collapsing_regions(ULONG &Ntsk,vector<real_prec>&sk_index)
{
  
  // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,
  So.message_screen("Getting FoF from Knots");
  
  // Initialize to -1 to then identify whether a cell it is used or not
  // This array will be finally used to allocate the bin in wehich a mass knot is found.
  this->SKNOT_M_info.clear();
  this->SKNOT_M_info.shrink_to_fit();
  this->SKNOT_M_info.resize(this->params._NGRID(), -1);
  
  // This vector contains the number of grids-cells that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(KNOT_MAX_SIZE,0);  // id of cells in superknots
  
  int k;
  int Nsk=0;// Number of super-knots identified
  
  for(ULONG i=0;i < this->params._NGRID(); ++i)
    {
      if(this->CWClass[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
	{
	  if(this->SKNOT_M_info[i] == -1) // if has not been already used in the fof, proceed
	    {
	      knot[0]=i;            // The first cell classified as knot. Container "knot" is always inizialized like this for every cell, if not yet visited (according to the if above)
	      this->SKNOT_M_info[i]=1;    //Mark the the current cell as  visited
	      int n=1; // n denotes the number of cells in SuperKnot. Starts from 1 as the current cells (if not visited before) has been identified as a knot. It gets added 1 if one of the neighbour cells is a knot. 

	      for(int j=0;j<n;++j)  // j runs over the neighbour cells of the current one. The limit "n" is increased inside the "if" when a neighbour knot is found such that we keep on searching for more friends
	   	{
		  // The following lines are based on the decomposition I(i,j,k)=k+Nj+N²i
		  // such that
		  // I(i,j,k+-1) = I+-1
		  // I(i,j+-1,k) = I+-N
		  // I(i+-1,j,k) = I+-N²
		  
		  if( (k= knot[j]-1)>=0 && k < this->params._NGRID()) 		  // I(i,j,k-1) = I-1
		    {
		      if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
			{
			  knot[n] = k;   
			  this->SKNOT_M_info[k] = 1; //mark this cell as already visited
			  n++;
			}
		    }
		  if( (k = knot[j] + 1) >= 0 && k < this->params._NGRID())  // I(i,j,k+1) = I+1
		    {
		      if( this->CWClass[k] >0 && this->SKNOT_M_info[k] == -1)
			{
			  knot[n] = k;
			  this->SKNOT_M_info[k] = 1;
			  n++;
			}
		    }
		  if( (k = knot[j] - this->params._Nft()) >= 0 && k < this->params._NGRID())  // I(i,j-1,k) = I-N
		    {
		      if( this->CWClass[k]>0 && this->SKNOT_M_info[k] == -1)
			{
			  knot[n] = k;
			  this->SKNOT_M_info[k] = 1;
			  n++;
			}
		    }
		  if( (k = knot[j] + this->params._Nft()) >= 0 && k < this->params._NGRID()) // I(i,j+1,k) = I+N
		    {
		      if(this->CWClass[k]>0 && this->SKNOT_M_info[k] == -1)
			{
			  knot[n] = k;
			  this->SKNOT_M_info[k] = 1;
			  n++;
			}
		    }
		  if( (k = knot[j] - this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID()) // I(i-1,j,k) = I-N²
		    {
		      if( this->CWClass[k]>0 && this->SKNOT_M_info[k] == -1)
			{
			  knot[n] = k;
			  this->SKNOT_M_info[k] = 1;
			  n++;
		      }
		    }
		  if( (k = knot[j] + this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID()) // I(i+1,j,k) = I+N²
		    {
		      if( this->CWClass[k] >0 && this->SKNOT_M_info[k] == -1)
			{
			  knot[n] = k;
			  this->SKNOT_M_info[k] = 1;
			  n++;
		      }
		    }
		}
	      for(int j = 0; j < n; j++)
		sk_index[knot[j]]=Nsk;// allocate the ID of supoerknot in which the cell lives, e.g., sk_index[211]=5155 means that the ID mesh 211 lives in teh superknot 5155, and 211=knot[j=45] measn that the 45th cell in the Sknot 5155 has ID 211 in the mesh
	      Nsk++;
	    }
	}
    }
  
  knot.clear();
  knot.shrink_to_fit();
  Ntsk=Nsk;
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*
void Cwclass::get_collapsing_regions(ULONG &Ntsk,vector<real_prec>&sk_index)
{
  
  // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,
  So.message_screen("Getting FoF from Knots");
  
  // Initialize to -1 to then identify whether a cell it is used or not
  // This array will be finally used to allocate the bin in wehich a mass knot is found.
  this->SKNOT_M_info.clear();
  this->SKNOT_M_info.shrink_to_fit();
  this->SKNOT_M_info.resize(this->params._NGRID(), -1);
  
  // This vector contains the number of grids-cells that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(KNOT_MAX_SIZE,0);  // id of cells in superknots
  
  int k;
  int Nsk=0;// Number of super-knots identified
  
  for(ULONG ic=0;ic < this->params._Nft(); ++ic)
    for(ULONG jc=0;jc < this->params._Nft(); ++jc)
      for(ULONG kc=0;kc < this->params._Nft(); ++kc)
	{
	  
	  ULONG i=index_3d(ic,jc,kc, this->params._Nft(), this->params._Nft());
	  if(this->CWClass[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
	    {
	      if(this->SKNOT_M_info[i] == -1) // if has not been already used in the fof, proceed
		{
		  knot[0]=i;            // The first cell classified as knot. Container "knot" is always inizialized like this for every cell, if not yet visited (according to the if above)
		  this->SKNOT_M_info[i]=1;    //Mark the the current cell as  visited
		  int n=1; // n denotes the number of cells in SuperKnot. Starts from 1 as the current cells (if not visited before) has been identified as a knot. It gets added 1 if one of the neighbour cells is a knot. 
		  int ni=1;
		  int nj=1;
		  int nk=1;
		  
		  for(ULONG ip=0;ip<ni;++ip)
		    {
		      ULONG ji=index_3d(ic+ip,jc,kc, this->params._Nft(), this->params._Nft());	      
		      if( (k= knot[ji]-1)>=0 && k < this->params._NGRID())
			{
			  if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
			    {
			      knot[ni] = k;   
			      this->SKNOT_M_info[k] = 1; //mark this cell as already visited
			      ni++;
			    }
			}
		      if( (k= knot[ji]+1)>=0 && k < this->params._NGRID())
			{
			  if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
			    {
			      knot[ni] = k;   
			      this->SKNOT_M_info[k] = 1; //mark this cell as already visited
			      ni++;
			    }
			}
		      
		      for(ULONG jp=0;ip<nj;++jp)
			{
			  ULONG jj=index_3d(ic+ip,jc+jp,kc, this->params._Nft(), this->params._Nft());	      
			  if( (k= knot[jj]-1)>=0 && k < this->params._NGRID())
			    {
			      if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
				{
				  knot[nj] = k;   
				  this->SKNOT_M_info[k] = 1; //mark this cell as already visited
				  nj++;
				}
			    }
			  if( (k= knot[jj]+1)>=0 && k < this->params._NGRID())
			    {
			      if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
				{
				  knot[nj] = k;   
				  this->SKNOT_M_info[k] = 1; //mark this cell as already visited
				  nj++;
				}
			    }
			  for(ULONG kp=0;kp<nk;++kp)
			    {
			      ULONG jk=index_3d(ic+ip,jc+jp,kc+kp, this->params._Nft(), this->params._Nft());	      
			      if( (k= knot[jk]-1)>=0 && k < this->params._NGRID())
				{
				  if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
				    {
				      knot[nk] = k;   
				      this->SKNOT_M_info[k] = 1; //mark this cell as already visited
				      nk++;
				    }
				}
			      if( (k= knot[jk]+1)>=0 && k < this->params._NGRID())
				{
				  if( this->CWClass[k]>0 && this->SKNOT_M_info[k]== -1)   // If density is greater than zero and cell not yet used
				    {
				      knot[nk] = k;   
				      this->SKNOT_M_info[k] = 1; //mark this cell as already visited
				      nk++;
				    }
				}
			    }  
			}  
		    }
		  
		  for(int j = 0; j < n; j++)
		    sk_index[knot[j]]=Nsk;// allocate the ID of supoerknot in which the cell lives, e.g., sk_index[211]=5155 means that the ID mesh 211 lives in teh superknot 5155, and 211=knot[j=45] measn that the 45th cell in the Sknot 5155 has ID 211 in the mesh
		  Nsk++;
		}
	    }
	}
  knot.clear();
  knot.shrink_to_fit();
  Ntsk=Nsk;
  So.DONE();
}
*/
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_SigmaVel_collapsing_regions(const vector<real_prec>&in, vector<real_prec>&Vx, vector<real_prec>&Vy,vector<real_prec>&Vz, real_prec Nmean)
{
  /*
    Routine taken from HADRON
  */

    // If we are in _GET_REALIZATIONS_, the limits here must be fixed, since different realizations might have different min and max,
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  cout<<endl;
  So.message_screen("Getting FoF from V-Knots");
#endif

  real_prec Hubble_function=this->s_cosmo_info.Hubble_parameter;
  real_prec cvel=Hubble_function/(cgs_km/cgs_Mpc);
  // Transform the velocities to v / H(z), putting the units of H(z) with the cvel, such that the shear is dimensionless

/*
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for(int i=0; i<this->params._NGRID();++i)
  {
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
     Vx[i]/=cvel;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
     Vy[i]/=cvel;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
     Vz[i]/=cvel;
  }

  */


  vector<real_prec>rho(in.size(),0);


  //Now compute teh velocity dispersion
#ifdef _FULL_VERBOSE_
  So.message_screen("Computing Velocity dispersion");
#endif
  this->Dispersion_VelField.resize(this->params._NGRID(),0);
  real_prec mean_vel_x=0;
  real_prec mean_vel_y=0;
  real_prec mean_vel_z=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_vel_x,mean_vel_y,mean_vel_z)
#endif
  for(int i=0; i<this->params._NGRID();++i)
    {
      mean_vel_x+=Vx[i];
      mean_vel_y+=Vy[i];
      mean_vel_z+=Vz[i];
    }
  mean_vel_x/=static_cast<double>(this->params._NGRID());
  mean_vel_y/=static_cast<double>(this->params._NGRID());
  mean_vel_z/=static_cast<double>(this->params._NGRID());
  mean_vel_x=sqrt(pow(mean_vel_x,2)+pow(mean_vel_y,2)+pow(mean_vel_z,2));
  So.message_screen("Mean DM velocity = ", mean_vel_x, "Mpc/h");

  if(Nmean>0)  // If Nob is zero, then we are passing a density field. If not, convert delta to density with the value of Nobs passed
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._NGRID();++i)
        rho[i]= Nmean*(num_1+in[i]);
    }
  else
    rho=in;


    ULONG aux=0;
#pragma omp parallel for reduction(+:aux)
    for(ULONG i=0;i<this->params._NGRID();++i)
      if(this->CWClass_V[i]==I_KNOT)
        aux++;

  // Initialize to -1 to then identify whether it is used or not

  this->VDISP_KNOT_info.clear();
  this->VDISP_KNOT_info.resize(this->params._NGRID(), -1);

  // This vector contains the number of grids that belongs to a SKNOT. One expects
  // KNOT_MAX_SIZE superknots
  vector<ULONG>knot(V_MAX_SIZE,0);

  int k;
  real_prec v_disp;
  real_prec VKrange_max = -1e38;  // Maximum  and Min of Mass of DM in knots. Mass here is number of DM particles in the super knot
  real_prec VKrange_min = +1e38; // These two variables are updated when knot masses are computed

  for(ULONG i=0;i < this->params._NGRID(); ++i)
    {
      if(this->CWClass_V[i]==I_KNOT && rho[i]>0)  // If cell is classified as knot, proceed. This must be the classification of the REFERENCE density  field. The >0 indicates that there might be some cells classified as knots but empty
        {
          if(this->VDISP_KNOT_info[i] == -1) // if has not been already used in the fof, proceed
            {
              knot[0]=i;            // The first cell classified as knot. knot is always inizialized like this for every cell
              this->VDISP_KNOT_info[i]=1;    //Mark the first knot as visited
              int n=1; // Starts from 1. It gets added 1 if one of the neighbour cells is a knot. Percolation.

              for(int j=0;j<n;++j)  // j runs over the neighbour cells of the current one, up to n, when n is increased inside if a neirbour knot is found
                {
                  if( (k= knot[j]-1)>=0 && k < this->params._NGRID())
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k]== -1)   // If density is greater than zero and cell not yet used
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[n] = 1; //mark this cell as already visited
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + 1) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] - this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if((this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] - this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if( (this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }

                  if( (k = knot[j] + this->params._Nft()*this->params._Nft()) >= 0 && k < this->params._NGRID())
                    if((this->CWClass_V[k]==I_KNOT && rho[k]>0) && this->VDISP_KNOT_info[k] == -1)
                      {
                        knot[n] = k;
                        this->VDISP_KNOT_info[k] = 1;
                        n++;
                        if(n >= V_MAX_SIZE)
                          So.message_error("Error: size of knots is larger than KNOT_MAX_SIZE, ", n);
                      }
                }

              v_disp = 0;

              for(int j = 0; j < n; j++) // Compute the mass*sigma² of the SK by adding the mass of its n parts
                 {
                   int index=knot[j]; // identify the index of the cell
                   real_prec vel_knots=sqrt(pow(Vx[index],2)+pow(Vy[index],2)+pow(Vz[index],2));
                   v_disp += rho[index]* pow(vel_knots - mean_vel_x,2)/static_cast<double>(n); // m sigma²
                }

              for(int j = 0; j < n; j++) //Assign to all cells in this superknot the mass of the sp they are in.
                this->VDISP_KNOT_info[knot[j]] = v_disp;

              // Adjust the interval for the bins in KMass to the current values min and max of the variable knotmass
              if(VKrange_min > v_disp) VKrange_min = v_disp;
              if(VKrange_max < v_disp) VKrange_max = v_disp;
            }
        }
      else
        this->VDISP_KNOT_info[i] = 0;
    }
 So.DONE();


 real_prec number_log;
#ifdef _DM_NEW_UNITS_
 number_log=0.0;
#else
 number_log=1.0;
#endif

#ifdef _FULL_VERBOSE_
 So.message_screen("Min and Max of Kinetic from FoF regions: ");
  So.message_screen("Minimum log(1+Kinetic) =",  log10(VKrange_min+number_log));
  So.message_screen("Maximim log(1+Kinetic) =",  log10(VKrange_max+number_log));
#endif
    // *** TYhis region is commnted if we want the current limits to be used in each step of the iteration
#ifndef _MODIFY_LIMITS_
  VKrange_max = VKMAX;
  VKrange_min = VKMIN;
  So.message_screen("log Used minimum Kinetic =", log10(VKrange_min+number_log));
  So.message_screen("log Used maximum Kinetic =", log10(VKrange_max+number_log));
  So.message_screen("");
#endif


  real_prec deltaMK= static_cast<real_prec>((log10(VKrange_max+number_log)-log10(VKrange_min+number_log))/static_cast<real_prec>(this->params._n_vknot_massbin()));
  for(ULONG i = 0; i < this->params._NGRID(); i++)
    {
      if(this->CWClass_V[i]==I_KNOT && rho[i]>0)
        {
          if(this->VDISP_KNOT_info[i] == -1) // If this cell was not used, exit
            {
              this->VDISP_KNOT_info[i] = 0;
              continue;
            }
          // GetA the bin in Knot mass. At this point SKNOT_M_info[i] encodes the mass of the Super knot where this cell is included
          if(log10(VDISP_KNOT_info[i]+number_log)<=log10(VKrange_max+number_log) && log10(VDISP_KNOT_info[i]+number_log)>=log10(VKrange_min+number_log) )
            {
              int j = static_cast<int>(floor((log10(VDISP_KNOT_info[i]+number_log)-log10(VKrange_min+number_log))/deltaMK));
              if(j == this->params._n_sknot_massbin()) j--;
              // Assign the bin to the sknot_m_info to be used in the P(x,y,...)
              this->VDISP_KNOT_info[i] = j;
            }
        }
      else
        this->VDISP_KNOT_info[i] = 0;
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// This function helps to asses whether a cell should be counted (true) or not (false)
// according to the CW classification (member object) and the classifications
// requested in the input parameter file.
bool Cwclass::get_cell_classified(int sua, ULONG ig)
{
  bool res=false;
  if(this->cwt_used.size()==1 && this->cwt_used[0]==0)
    res=true;
  else
    {
      int classi=this->cwt_used[sua]; // This is the value passed fro the parameter file
      if(classi==0)
	res=true;
      
      else if(classi>0 && classi<5)  // This refers to the 4 typical classifications
	{
	  if(this->CWClass[ig]==classi)
	    res=true;
	}
      else         // This referes to combinations
	if(classi==234)
	  {
            if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET ||this->CWClass[ig]==I_VOID)
	      res=true;
	  }
	else if(classi==12)
	  {
            if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT)
	      res=true;
	  }
	else if(classi==23)
	  {
            if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
	      res=true;
	  }
	else if(classi==34)
	  {
            if(this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
	      res= true;
	  }
	else if(classi==123)
	  {
            if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
	      res= true;
	  }
    }
  return res;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/** 
 * @brief This function provides the index identifying the CWT classification
 * as proposed in the parameter file, for each cell, and according
 * to the CWClass array
 */

int Cwclass::get_Tclassification(ULONG ig)
{

  int ans=0;

  if(this->cwt_used.size()==1)
    ans=this->cwt_used[0];
  else
    {
      for(int ic=0;ic<this->cwt_used.size();++ic)// loop over the number of cwt required from the par file
	{
	  int cwc=this->cwt_used[ic];  // the cwt asked in par file. E.g, 1, 234
	  if(cwc<5)
	    {
	      if(this->CWClass[ig]==cwc)
		ans=ic;
	    }
	  else{
	    if(cwc==234)
	      {
                if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
		  ans=ic;
	      }
	    else if(cwc==12)
	      {
                if(this->CWClass[ig]==I_KNOT || this->CWClass[ig]==I_FILAMENT)
		  ans=ic;
	      }
	    else if(cwc==23)
	      {
                if(this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
		  ans=ic;
	      }
	    else if(cwc==34)
	      {
                if(this->CWClass[ig]==I_SHEET || this->CWClass[ig]==I_VOID)
		  ans=ic;
	      }
	    
	    else if(cwc==123)
	      {
                if(this->CWClass[ig]==I_KNOT ||  this->CWClass[ig]==I_FILAMENT || this->CWClass[ig]==I_SHEET)
		  ans=ic;
	      }
	    
	  }
	}
    }
  
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

/**
 * @brief This function provides the index identifying the CWT classification
 * as proposed in the parameter file, for each cell, and according
 * to the CWClass array
 */

int Cwclass::get_Vclassification(int ig)
{

  int ans=0;

  if(this->cwv_used.size()==1)
    ans=this->cwv_used[0];
  else
    {
      for(int ic=0;ic<this->cwv_used.size();++ic)// loop over the number of cwt required from the par file
        {
          int cwc=this->cwv_used[ic];  // the cwt asked in par file. E.g, 1, 234
          if(cwc<5)
            {
              if(this->CWClass_V[ig]==cwc)
                ans=ic;
            }
          else{
            if(cwc==234)
              {
                if(this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET || this->CWClass_V[ig]==I_VOID)
                  ans=ic;
              }
            else if(cwc==12)
              {
                if(this->CWClass_V[ig]==I_KNOT || this->CWClass_V[ig]==I_FILAMENT)
                  ans=ic;
              }
            else if(cwc==23)
              {
                if(this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET)
                  ans=ic;
              }
            else if(cwc==34)
              {
                if(this->CWClass_V[ig]==I_SHEET || this->CWClass_V[ig]==I_VOID)
                  ans=ic;
              }

            else if(cwc==123)
              {
                if(this->CWClass_V[ig]==I_KNOT ||  this->CWClass_V[ig]==I_FILAMENT || this->CWClass_V[ig]==I_SHEET)
                  ans=ic;
              }

          }
        }
    }

  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::Konvolve(vector<real_prec> &in, vector<real_prec>&out, string type)
{
  So.message_screen("Using Konvolve in Cwclass");
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef DOUBLE_PREC
  complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
  complex_prec *data_out= (complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID_h();++i){
    data_out[i][REAL]=0;
    data_out[i][IMAG]=0;
  }
  do_fftw_r2c(this->params._Nft(),in, data_out);
  double paux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:paux)
#endif
  for(ULONG i=0;i<this->params._NGRID();i++)
   paux+=static_cast<double>(in[i]);
/*
   if(true==isinf(paux)){
       So.message_warning("Not defined value found in container at function ",__PRETTY_FUNCTION__);
       So.message_warning("Code exits here");
       exit(0);
}
*/
#ifdef _EXTRAPOLATE_VOLUME_
  real_prec correction_factor = 1.0;
#else
  real_prec correction_factor = 1.00;
#endif
  double we=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG ind=0;ind< this->params._NGRID_h() ;++ind)
    {
      real_prec cor = this->Kernel[ind]*correction_factor;
      we+=cor;

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      data_out[ind][REAL]*=cor;

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      data_out[ind][IMAG]*=cor;
    }
  do_fftw_c2r(this->params._Nft(), data_out, out);
  // Here I have to correct for the normalization of the kernel
  // for the function  do_fftw_c2r returns the transform normalized by the this->params._NGRID(), so I divide by this->params._NGRID() and by multiply by 2 we
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID();i++)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
    out[i]/=(static_cast<double>(this->params._NGRID())/static_cast<double>(2.0*we));
  So.DONE();
  
#ifdef DOUBLE_PREC
  fftw_free(data_out);
#else
  fftwf_free(data_out);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_tidal_anisotropy(string input_file, string output_file){
    vector<real_prec>delta(this->params._NGRID(),0);
    this->File.read_array(input_file, delta);
    for(ULONG i=0;i<delta.size();++i)
    cout<<delta[i]<<endl;
    get_overdens(delta,delta);
    for(ULONG i=0;i<delta.size();++i)
        delta[i]=tidal_anisotropy(this->lambda1[i],this->lambda2[i],this->lambda3[i]);
    this->File.write_array(output_file, delta);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void Cwclass::get_tidal_anisotropy(vector<real_prec>&delta, string output_file){
    get_overdens(delta,delta);
    this->get_CWC(delta);
    for(ULONG i=0;i<delta.size();++i)
        delta[i]=tidal_anisotropy(this->lambda1[i],this->lambda2[i],this->lambda3[i]);
    this->File.write_array(output_file, delta);
}
