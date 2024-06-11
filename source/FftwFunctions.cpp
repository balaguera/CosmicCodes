////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @class<FftwFunctions>
 *  @file FftwFunctions.cpp
 *  @brief Methods of the class FftwFunctions
 *  @details Implementation of the methods of the class FftwFunctions,
 *  used to measure the 3D power spectrum and bispectrum (FKP and Yamamoto)
 *  @details Methods are based on the FFTW algrithm
 *  @author Andrés Balaguera-Antolínez,
 *  @author Jennifer Pollack wrote the original Fortran implementation of the Bispectrum. The original methods were first explained by Peter Schuecker in 2007.
 *  @date 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <string.h>
#ifdef _USE_OMP_
#include <omp.h>
#endif
#include<vector>
#include "FftwFunctions.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
#define ijk(i, j, k, nn1, nn2, nn3) (index_3d(i,j,k,nn2,nn3))
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::resize_fftw_vectors()
{
  this->data_g.resize(this->params._NGRID(),0.0);
  if(true==this->params._use_real_and_redshift_space())
      this->data_g_rss.resize(this->params._NGRID());
#ifdef _MASS_WEIGHT_POWER_
  this->data_g_mw.resize(this->params._NGRID());
#endif
  if(true==this->params._measure_cross() || true == this->params._use_real_and_redshift_space()) // this last one is becuase I needed a dirt/quick solution
    this->data_gp.resize(this->params._NGRID());
  if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_ysc")
    {
      this->data_g_xx.resize(this->params._NGRID(),0.0);
      this->data_g_yy.resize(this->params._NGRID(),0.0);
      this->data_g_zz.resize(this->params._NGRID(),0.0);
      this->data_g_xy.resize(this->params._NGRID(),0.0);
      this->data_g_xz.resize(this->params._NGRID(),0.0);
      this->data_g_yz.resize(this->params._NGRID(),0.0);
    }
  if(this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc" )
    {
      this->data_g_xx.resize(this->params._NGRID(),0.0);
      this->data_g_yy.resize(this->params._NGRID(),0.0);
      this->data_g_zz.resize(this->params._NGRID(),0.0);
      this->data_g_xy.resize(this->params._NGRID(),0.0);
      this->data_g_xz.resize(this->params._NGRID(),0.0);
      this->data_g_yz.resize(this->params._NGRID(),0.0);
      this->data_g_xxx.resize(this->params._NGRID(),0.0);
      this->data_g_yyy.resize(this->params._NGRID(),0.0);
      this->data_g_zzz.resize(this->params._NGRID(),0.0);
      this->data_g_xxy.resize(this->params._NGRID(),0.0);
      this->data_g_xxz.resize(this->params._NGRID(),0.0);
      this->data_g_yyx.resize(this->params._NGRID(),0.0);
      this->data_g_yyz.resize(this->params._NGRID(),0.0);
      this->data_g_zzx.resize(this->params._NGRID(),0.0);
      this->data_g_zzy.resize(this->params._NGRID(),0.0);
      this->data_g_xyy.resize(this->params._NGRID(),0.0);
      this->data_g_xzz.resize(this->params._NGRID(),0.0);
      this->data_g_yzz.resize(this->params._NGRID(),0.0);
      this->data_g_xyz.resize(this->params._NGRID(),0.0);
      this->data_g_yxz.resize(this->params._NGRID(),0.0);
      this->data_g_zxy.resize(this->params._NGRID(),0.0);
    }
    this->So.message_screen("Allocating memmory");
#ifdef DOUBLE_PREC
  this->data_out_g   =(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
  this->data_out_g   =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID_h();++i)data_out_g[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID_h();++i)data_out_g[i][IMAG]=0;
  if(true==this->params._measure_cross() || true == this->params._use_real_and_redshift_space())
  {
#ifdef DOUBLE_PREC
    data_out_gp   =(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
  data_out_gp   =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._NGRID_h();++i)data_out_gp[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID_h();++i)data_out_g[i][IMAG]=0;
}
  if(data_out_g==NULL)
  {
     this->So.message_screen("Allocation failed, check line", __LINE__);
     exit(0);
  }
  if(true==this->params._measure_cross() || true == this->params._use_real_and_redshift_space())
    if(data_out_gp==NULL)
    {
      this->So.message_screen("Allocation failed, check line", __LINE__);
      exit(0);
    }
  if(true==this->params._use_random_catalog())
    this->data_r.resize(this->params._NGRID(),0.0);
#ifdef DOUBLE_PREC
  data_out_r=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
  data_out_r=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
  if(this->params._statistics()=="Pk_y_ds")
  {
    data_g_out_y0.resize(this->params._NGRID_h(),0.0);
    data_g_out_y2.resize(this->params._NGRID_h(),0.0);
    data_g_out_y4.resize(this->params._NGRID_h(),0.0);
    data_r_out_y0.resize(this->params._NGRID_h(),0.0);
    data_r_out_y2.resize(this->params._NGRID_h(),0.0);
    data_r_out_y4.resize(this->params._NGRID_h(),0.0);
    SN_r_out_y2.resize(this->params._NGRID_h(),0.0);
    SN_r_out_y4.resize(this->params._NGRID_h(),0.0);
    SN_g_out_y2.resize(this->params._NGRID_h(),0.0);
    SN_g_out_y4.resize(this->params._NGRID_h(),0.0);
    data_g_y0.resize(this->params._NGRID_h(),0.0);
    data_g_y2.resize(this->params._NGRID_h(),0.0);
    data_g_y4.resize(this->params._NGRID_h(),0.0);
  }
  if(this->params._statistics()=="Pk_fkp" && this->params._FKP_error_bars()==true)
  {
    SN.resize(this->params._NGRID(),0);
    Q.resize(this->params._NGRID(),0);
  }
  if(this->params._statistics()=="Bk_fkp_fast")
    {
      this->Array_corr.resize(this->params._NGRID_h(),1.0);
      Arraykx.resize(this->new_sgrid,0);
      Arrayky.resize(this->new_sgrid,0);
      Arraykz.resize(this->new_sgrid,0);
      Arraykk.resize(this->new_sgrid,0);
      ArrayID.resize(this->new_sgrid,0);
      VecArray.resize(this->new_sgrid,0);
      Kbin.resize(this->new_sgrid,0);
      Bmodes.resize(this->Nshells_bk,0);
      kkminID.resize(this->Nshells_bk,0);
      kbins_bk.resize(this->Nshells_bk,0);
      Ngrids_bk.resize(this->Nshells_bk,0);
    }
  this->So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::raw_sampling(real_prec vol)
{
  So.enter(__PRETTY_FUNCTION__);
  this->normal_p=static_cast<real_prec>(this->n_gal)*static_cast<real_prec>(this->n_gal)/vol;
  this->normal_b=static_cast<real_prec>(this->n_gal)*pow(this->n_gal/vol,2);
  this->w_r=static_cast<real_prec>(this->n_gal);
  this->s_r=vol/static_cast<real_prec>(this->n_gal);
  this->sr1=this->n_gal*(this->n_gal/vol);
  this->sr2=static_cast<real_prec>(this->n_gal);
  this->n_ran=0;
  this->alpha=1.0;
  this->normal_power=normal_p;
  this->shot_noise=this->s_r;
  this->shot_noise_window=0;
  this->normal_window=1.0 ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_parameters_estimator(bool verbose)
{
  So.enter(__PRETTY_FUNCTION__);
  if(true==this->params._use_random_catalog())
    {
      this->alpha= static_cast<double>(this->w_g)/static_cast<double>(this->w_r);
      this->normal_power= this->alpha*this->normal_p;
      this->shot_noise= (alpha+1)*(s_r/normal_p);  // originally is alpha*(alpha+1)*(s_r/ normal_power);
      // this->shot_noise= s_g/ normal_power + pow( alpha,2)*( s_r/ normal_power);
      this->normal_window= normal_power/pow(alpha,2);
      this->shot_noise_window= (this->alpha*this->alpha)* s_r/ this->normal_power;
      this->normal_bispectrum= this->alpha* normal_b;
      this->shot_noise_b1= sr1/( this->alpha* this->normal_bispectrum);
      this->shot_noise_b2=(1.- this->alpha* this->alpha)* sr2/( this->alpha* this->normal_bispectrum);
      So.message_screen("** Ngal            =", this->n_gal);
      So.message_screen("** Nran            =", this->n_ran);
      So.message_screen("** Wg              =", this->w_g);
      So.message_screen("** Wr              =", this->w_r);
      So.message_screen("** Sg              =", this->s_g);
      So.message_screen("** Sr              =", this->s_r);
      So.message_screen("** alpha           =", this->alpha);
      So.message_screen("** Normalization   =", this->normal_power);
      So.message_screen("** Shot Noise      =", this->shot_noise);
    }
  else{
    this->alpha=1.0; //This alfa should be zero when using a box. We set it to one correspondingly with that case,
    // for all quantities are explicitely comnputed in that situation without assuming anything on this parameter
    this->normal= normal_power;
    this->shot_noise= s_r;
    this->normal_window=1.0;
    this->shot_noise_window=1.0; //care
    this->normal_bispectrum= normal_b;
    this->shot_noise_b1= this->sr1/(this->normal_b);
    this->shot_noise_b2= this->sr2/(normal_b);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_fluctuation(bool r_rss)
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
       this->data_g[i]-=static_cast<real_prec>(this->n_gal)/static_cast<real_prec>(this->params._NGRID());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
      this->data_g_rss[i]-=static_cast<real_prec>(this->n_gal)/static_cast<real_prec>(this->params._NGRID());
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_fluctuation()
{
  So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
  omp_set_num_threads(_NTHREADS_);
#endif

#ifdef _FULL_VERBOSE_POWER_
  So.message_screen("Computing galaxy fluctuation");
#endif
  // Define the coordinates of the cells within the grid
  // to be used in the Yamamoto estimator:
  vector<real_prec>cell_x;
  vector<real_prec>cell_y;
  vector<real_prec>cell_z;
  if(this->params._statistics()=="Pk_ybc"  || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_yb"  || this->params._statistics()=="Pk_ys")
  {
      cell_x.resize(this->params._NGRID(),0.0);
      cell_y.resize(this->params._NGRID(),0.0);
      cell_z.resize(this->params._NGRID(),0.0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
      for(int i=0;i<this->params._Nft();++i)
       for(int j=0;j<this->params._Nft();++j)
        for(int k=0;k<this->params._Nft();++k)
          {
            ULONG lp=index_3d(i ,j ,k, this->params._Nft(),this->params._Nft());
            cell_x[lp]=(i+0.5)*this->params._d_delta_x() + (this->params._Xoffset() - this->params._Lbox()*0.5);
            cell_y[lp]=(j+0.5)*this->params._d_delta_y() + (this->params._Yoffset() - this->params._Lbox()*0.5);
            cell_z[lp]=(k+0.5)*this->params._d_delta_z() + (this->params._Zoffset() - this->params._Lbox()*0.5);
          }
  }
  if(this->params._statistics()=="Pk_ds")
  {
    for(int i=0;i<this->params._NGRID_h();++i){
      this->data_g_out_y0[i] -= this->alpha*this->data_r_out_y0[i];
      this->data_g_out_y2[i] -= this->alpha*this->data_r_out_y2[i];
      this->data_g_out_y4[i] -= this->alpha*this->data_r_out_y4[i];
    }
  }
#ifdef _USE_OMP_
#pragma omp parallel
  {
#endif
      if(true==this->params._use_random_catalog())
      {
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
        for(ULONG i=0;i<this->params._NGRID();++i)
          this->data_g[i] -= this->alpha*this->data_r[i];  /// <----- FLUCTUATION
      }
    else
      {
#ifdef _MASS_WEIGHT_POWER_
#pragma omp for nowait
           for(int i=0;i<this->params._NGRID();++i)
               this->data_g[i]=this->data_g_mw[i]-data_g[i];
#else
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
          for(ULONG i=0;i<this->params._NGRID();++i)
                 this->data_g[i]-=static_cast<real_prec>(this->n_gal)/static_cast<real_prec>(this->params._NGRID());
#endif
      }
  // The following definitions are needed for estimators of multipoles different to FKP
  if(this->params._statistics()!="Pk_fkp")
  {
    if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_ysc")
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
      for(int i=0;i<this->params._NGRID();++i)
        {
          // Define the fluctuation F(r) interpolated
          real_prec delta_cc=this->data_g[i];
          real_prec absr2=cell_x[i]*cell_x[i]+cell_y[i]*cell_y[i]+cell_z[i]*cell_z[i];
          this->data_g_xx[i]=(cell_x[i]*cell_x[i]/absr2)*delta_cc;
          this->data_g_yy[i]=(cell_y[i]*cell_y[i]/absr2)*delta_cc;
          this->data_g_zz[i]=(cell_z[i]*cell_z[i]/absr2)*delta_cc;
          this->data_g_xy[i]=(cell_x[i]*cell_y[i]/absr2)*delta_cc;
          this->data_g_xz[i]=(cell_x[i]*cell_z[i]/absr2)*delta_cc;
          this->data_g_yz[i]=(cell_y[i]*cell_z[i]/absr2)*delta_cc;
         }
    else if(this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc")
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
      for(int i=0;i<this->params._NGRID();++i)
          {
              // Compute squared of the distance of the cell to the origin
              real_prec delta_cc=this->data_g[i];
              real_prec absr2=cell_x[i]*cell_x[i]+cell_y[i]*cell_y[i]+cell_z[i]*cell_z[i];
              real_prec iabsr4=1./(absr2*absr2);
              this->data_g_xx[i]=(cell_x[i]*cell_x[i]/absr2)*delta_cc;
              this->data_g_yy[i]=(cell_y[i]*cell_y[i]/absr2)*delta_cc;
              this->data_g_zz[i]=(cell_z[i]*cell_z[i]/absr2)*delta_cc;
              this->data_g_xy[i]=(cell_x[i]*cell_y[i]/absr2)*delta_cc;
              this->data_g_xz[i]=(cell_x[i]*cell_z[i]/absr2)*delta_cc;
              this->data_g_yz[i]=(cell_y[i]*cell_z[i]/absr2)*delta_cc;
              this->data_g_xxx[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_x[i])*delta_cc;
              this->data_g_yyy[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_y[i])*delta_cc;
              this->data_g_zzz[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_z[i])*delta_cc;
              this->data_g_xxy[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_y[i])*delta_cc;
              this->data_g_xxz[i]=iabsr4*(pow(cell_x[i],2)*cell_x[i]*cell_z[i])*delta_cc;
              this->data_g_yyx[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_x[i])*delta_cc;
              this->data_g_yyz[i]=iabsr4*(pow(cell_y[i],2)*cell_y[i]*cell_z[i])*delta_cc;
              this->data_g_zzx[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_x[i])*delta_cc;
              this->data_g_zzy[i]=iabsr4*(pow(cell_z[i],2)*cell_z[i]*cell_y[i])*delta_cc;
              this->data_g_xyy[i]=iabsr4*(pow(cell_x[i],2)*cell_y[i]*cell_y[i])*delta_cc;
              this->data_g_xzz[i]=iabsr4*(pow(cell_x[i],2)*cell_z[i]*cell_z[i])*delta_cc;
              this->data_g_yzz[i]=iabsr4*(pow(cell_y[i],2)*cell_z[i]*cell_z[i])*delta_cc;
              this->data_g_xyz[i]=iabsr4*(pow(cell_x[i],2)*cell_y[i]*cell_z[i])*delta_cc;
              this->data_g_yxz[i]=iabsr4*(pow(cell_y[i],2)*cell_x[i]*cell_z[i])*delta_cc;
              this->data_g_zxy[i]=iabsr4*(pow(cell_z[i],2)*cell_x[i]*cell_y[i])*delta_cc;
        }
    }
#ifdef _USE_OMP_
  }
#endif
  if(this->params._statistics()=="Pk_ybc"  || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_yb"  || this->params._statistics()=="Pk_ys")
  {
      cell_x.clear();cell_x.shrink_to_fit();
      cell_y.clear();cell_y.shrink_to_fit();
      cell_z.clear();cell_z.shrink_to_fit();
  }
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_power_spectrum_fkp(vector<real_prec> &power_g0, vector<real_prec> &power_g2,vector<real_prec> &power_g4,vector<real_prec> &power_r,  vector< vector<real_prec> >&power2d_cart,vector< vector<real_prec> >&power2d_spher,vector<int> &mod_g)
{
//#define ICONST
  So.enter(__PRETTY_FUNCTION__);
  do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
  if(true==this->params._measure_cross())
    if(true==this->params._use_random_catalog())
    do_fftw_r2c(this->params._Nft(),this->data_gp, this->data_out_gp);
  
  if(true==this->params._use_random_catalog())
  {
    do_fftw_r2c(this->params._Nft(),this->data_r, this->data_out_r);
#ifdef ICONST
    real_prec norm_ic=sqrt(pow(this->data_out_r[0][REAL],2)+pow(this->data_out_r[0][IMAG],2));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID_h();++i)
    {
      this->data_out_r[i][REAL]/=norm_ic;
      this->data_out_r[i][IMAG]/=norm_ic;
    }
#endif
  }

//#define EUTEST
#ifdef EUTEST
  ULONG NTT=this->params._Nft()*this->params._Nft()*(this->params._Nft()/2+1);
  this->So.message_screen("Reading densiy field for EUCLID test");
  vector<real_prec> field(2*NTT,0);
  this->File.read_array("/home/balaguera/data/Numerics/Euclid/Tests/data.dat",field);  
  vector<real_prec> field_r(2*NTT,0);
  this->File.read_array("/home/balaguera/data/Numerics/Euclid/Tests/random.dat",field_r);
  cout<<this->alpha<<endl;
#pragma omp parallel for
  for(ULONG i=0; i<NTT;++i)
  { 
    this->data_out_g[i][REAL]=static_cast<real_prec>(field[2*i]-this->alpha*field_r[2*i]);
    this->data_out_g[i][IMAG]=static_cast<real_prec>(field[2*i+1]-this->alpha*field_r[2*i+1]);
    this->data_out_r[i][REAL]=static_cast<real_prec>(field_r[2*i]);
    this->data_out_r[i][IMAG]=static_cast<real_prec>(field_r[2*i+1]);
  }
#endif
  this->power_spectrum_fkp(power_g0,power_g2,power_g4,power_r,power2d_cart,power2d_spher,mod_g);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_power_spectrum_fkp(vector<real_prec> &power_g, vector<real_prec> &power_g0, vector<real_prec> &power_g2,vector<real_prec> &power_g4,vector<real_prec> &power_r,  vector< vector<real_prec> >&power2d_cart,vector< vector<real_prec> >&power2d_spher,vector<int> &mod_g)
{
#ifdef _VERBOSE_POWER_
  So.enter(__PRETTY_FUNCTION__);
#endif
  if(true==this->params._use_random_catalog())
    do_fftw_r2c(this->params._Nft(),this->data_r, this->data_out_r);
  // Real space:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID_h();++i)
    {
      data_out_g[i][REAL]=0;
      data_out_g[i][IMAG]=0;
    }
  this->power_spectrum(power_g,mod_g); // the FFTW is done inside this method, OJO
  //Rss. We use the same this->data_out_g container
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._NGRID_h();++i)
    {
      data_out_g[i][REAL]=0;
      data_out_g[i][IMAG]=0;
    }
  do_fftw_r2c(this->params._Nft(),this->data_g_rss, this->data_out_g);
  this->power_spectrum_fkp(power_g0,power_g2,power_g4,power_r,power2d_cart,power2d_spher,mod_g);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_power_spectrum_yamamoto(vector<real_prec> &power_g0, vector<real_prec> &power_g2,vector<real_prec> &power_g4,vector<int> &mod_g)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  if(this->params._statistics()=="Pk_y_ds")
  {
    real_prec sn_aux =  this->params._SN_correction() ? 1.0: 0.0;
    real_prec r_normal_power = 1.0/this->normal_power;
    real_prec alpha2 = this->alpha * this->alpha;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<this->params._NGRID_h();++i){
      // Moment l=0: should be equal to the monopole from FKP
      // Given that this->shot_noise ya viene nomalizado, multiplico por normal_power
      // para cancelar y dividir todo el parentesis de nuevo por la normalizacion
      this->data_g_y0[i] = (norm(this->data_g_out_y0[i]) - sn_aux*this->shot_noise*this->normal_power) * r_normal_power;
      // Moment l=2
      this->data_g_y2[i]=(2.*2.0+1.)*(this->data_g_out_y2[i].real()*this->data_g_out_y0[i].real()+this->data_g_out_y2[i].imag()*this->data_g_out_y0[i].imag()
                                      - sn_aux*(this->SN_g_out_y2[i]+alpha2*SN_r_out_y2[i])) * r_normal_power;
      // Moment l=4
      this->data_g_y4[i]=(2.*4.0+1.)*(this->data_g_out_y4[i].real()*this->data_g_out_y0[i].real()+this->data_g_out_y4[i].imag()*this->data_g_out_y0[i].imag()
                                      - sn_aux*(this->SN_g_out_y4[i]+ alpha2*this->SN_r_out_y4[i])) * r_normal_power;
    }
    this->power_yam_1d_ds(power_g0, power_g2, power_g4,mod_g);
  }
  else
  {
    So.message_screen("Evaluating FFTW using ", omp_get_max_threads(), " threads");
    time_t start;
    time (&start);
    do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
    this->data_g.clear(); data_g.shrink_to_fit();  // Release memmory here
    if(this->params._statistics()=="Pk_ys"  || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc")
    {
        //----
#ifdef DOUBLE_PREC
     data_out_g_xx=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_xx=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_xx, this->data_out_g_xx);
      this->data_g_xx.clear(); data_g_xx.shrink_to_fit();  // Release memmory here
      //----
#ifdef DOUBLE_PREC
     data_out_g_yy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_yy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_yy, this->data_out_g_yy);
      this->data_g_yy.clear(); data_g_yy.shrink_to_fit();  // Release memmory here
      //----
#ifdef DOUBLE_PREC
     data_out_g_zz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_zz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_zz, this->data_out_g_zz);
      this->data_g_zz.clear(); data_g_zz.shrink_to_fit();  // Release memmory here
      //----
#ifdef DOUBLE_PREC
     data_out_g_xy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_xy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_xy, this->data_out_g_xy);
      this->data_g_xy.clear(); data_g_xy.shrink_to_fit();  // Release memmory here
      //----
#ifdef DOUBLE_PREC
     data_out_g_xz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_xz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_xz, this->data_out_g_xz);
      this->data_g_xz.clear(); data_g_xz.shrink_to_fit();  // Release memmory here
      //----
#ifdef DOUBLE_PREC
     data_out_g_yz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
     data_out_g_yz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
      do_fftw_r2c(this->params._Nft(),this->data_g_yz, this->data_out_g_yz);
      this->data_g_yz.clear(); data_g_yz.shrink_to_fit();  // Release memmory here
    }
    if(this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc")
      {
#ifdef DOUBLE_PREC
        data_out_g_xxx=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xxx=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_xxx, this->data_out_g_xxx);
        this->data_g_xxx.clear(); data_g_xxx.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_yyy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_yyy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_yyy, this->data_out_g_yyy);
        this->data_g_yyy.clear(); data_g_yyy.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_zzz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_zzz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_zzz, this->data_out_g_zzz);
        this->data_g_zzz.clear(); data_g_zzz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_xxy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xxy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
       do_fftw_r2c(this->params._Nft(),this->data_g_xxy, this->data_out_g_xxy);
        this->data_g_xxy.clear(); data_g_xxy.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_xxz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xxz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_xxz, this->data_out_g_xxz);
        this->data_g_xxz.clear(); data_g_xxz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_yyx=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_yyx=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_yyx, this->data_out_g_yyx);
        this->data_g_yyx.clear(); data_g_yyx.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_yyz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_yyz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_yyz, this->data_out_g_yyz);
        this->data_g_yyz.clear(); data_g_yyz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_zzx=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_zzx=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_zzx, this->data_out_g_zzx);
        this->data_g_zzx.clear(); data_g_zzx.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_zzy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_zzy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_zzy, this->data_out_g_zzy);
        this->data_g_zzy.clear(); data_g_zzy.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_xyy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xyy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_xyy, this->data_out_g_xyy);
        this->data_g_xyy.clear(); data_g_xyy.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_xzz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xzz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_xzz, this->data_out_g_xzz);
        this->data_g_xzz.clear(); data_g_xzz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_yzz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_yzz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_yzz, this->data_out_g_yzz);
        this->data_g_yzz.clear(); data_g_yzz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_xyz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_xyz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_xyz, this->data_out_g_xyz);
        this->data_g_xyz.clear(); data_g_xyz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_yxz=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_yxz=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_yxz, this->data_out_g_yxz);
        this->data_g_yxz.clear(); data_g_yxz.shrink_to_fit();  // Release memmory here
#ifdef DOUBLE_PREC
        data_out_g_zxy=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#else
        data_out_g_zxy=(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
        do_fftw_r2c(this->params._Nft(),this->data_g_zxy, this->data_out_g_zxy);
        this->data_g_zxy.clear(); data_g_zxy.shrink_to_fit();  // Release memmory here
      }
    if(this->params._use_random_catalog())
    {
      do_fftw_r2c(this->params._Nft(),this->data_r, this->data_out_r);
      this->data_r.clear(); data_r.shrink_to_fit();  // Release memmory here
    }
    power_spectrum_yamamoto_new(power_g0,power_g2,power_g4,mod_g);    // Compute shell-average //
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
void FftwFunctions::cross_power_spectrum_fkp(vector<real_prec> & powerk_g0,vector<int> & mod_g, vector<real_prec>&corr)
#else
void FftwFunctions::cross_power_spectrum_fkp(vector<real_prec> & powerk_g0,vector<int> & mod_g, bool corr_factor)
#endif
{
  this->So.enter(__PRETTY_FUNCTION__);
    this->So.message_warning("There is some shot noise subtraction missing in this function. Please check");
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
  do_fftw_r2c(this->params._Nft(),this->data_gp, this->data_out_gp);
  std::fill(mod_g.begin(),mod_g.end(),0);
  std::fill(powerk_g0.begin(),powerk_g0.end(),0);
  vector<real_prec>aux2(mod_g.size(),0);
  vector<real_prec>aux4(mod_g.size(),0);
  real_prec rDeltaK_data=0;
  real_prec rkmin=0;
  int my_s_box_ave=1;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
 complex_prec *CROSS=(complex_prec *)fftw_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
  my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear")
  {
    my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
  }
  else if(this->params._type_of_binning()=="log")
    {
      my_s_box_ave = s_box_log;
      rkmin = 1.0/this->params._d_kmin();
    }
  vector<real_prec>yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>zz_MAS_array(this->params._Nft()/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>sn_zz_MAS_array(this->params._Nft()/2+1,0);
  real_prec M_PI_rNft = M_PI * this->rNft;
  real_prec SN_aux=0;
  if(true==this->params._SN_correction())
    SN_aux=this->shot_noise;
  real_prec SN_aux2=0;
  if(true==this->params._SN_correction())
    SN_aux2=this->shot_noise2;
#pragma omp parallel num_threads(NTHREADS)
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;
#pragma omp sections
    {
#pragma omp section
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->params._Nft()/2;++j) {
          yy_MAS = j * M_PI_rNft;
          yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
          sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
        }
      }
#pragma omp section
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->params._Nft()/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
          sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
        }
      }
    }
#pragma omp barrier
#pragma omp for nowait
    for(int i=0;i<this->params._Nft()/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
        if(i > 0)
          {
            xx_MAS = i * M_PI_rNft;
            xx_MAS = sin(xx_MAS)/xx_MAS;
          }
        else xx_MAS = 1.0;
        sn_xx_MAS = sin(i * M_PI_rNft);
        i_deltak_x_2 = i * i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
        i_per_fact = i_deltak_x_2;
        if(my_s_box_ave == s_box_log)
          i_deltak_x_2 *= rkmin * rkmin;
        for(int j=0;j<this->params._Nft()/2;++j)
          {  // loop over octant ky>0. Half of what FFT gives
            yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];
            j_deltak_y_2 = j * j * this->params._d_deltak_y() * this->params._d_deltak_y();
            j_per_fact = i_per_fact + j_deltak_y_2;
            j_per_fact = sqrt(j_per_fact);
            if(my_s_box_ave == s_box_log)
              j_deltak_y_2 *= rkmin * rkmin;
            for(int k=0;k<=this->params._Nft()/2;++k)
              {  // loop over octant kz>0: this is all FFTW gives.
                zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];
                k_deltak_z_2 = k * k * this->params._d_deltak_z() * this->params._d_deltak_z();
                if(my_s_box_ave == s_box_log)
                  k_deltak_z_2 *= rkmin * rkmin;
                real_prec kv;
                int kmod_g, kmod_r;
                // Compute k-shell index
                if(my_s_box_ave == s_box_linear)
                  {
                    kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
                    kmod_g=(int)floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data));
                  }
                else
                  {
                    if(my_s_box_ave == s_box_log)
                      {
                        kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
                        if(kv!=0){
                          kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv-this->params._d_kmin()))/this->params._d_Deltal()))+1;
                          kmod_r=kmod_g;
                        }
                        else{kmod_g=0;kmod_r=0;}
                      }
                  }
                // Compute correction for MAS:
                real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
                // Compute index in c-order
                ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                // Compute thisgs related to 2d power spectrum and multipole
                // decomposition for the first octant
                // Define the angle bewtween k and los:
                // In the FKP we need to specify the LOS direction. set to z by default.
                if(lp!=0)
                  {
                    real_prec Pk= (this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                    real_prec Pk1= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;
                    real_prec Pk2= (this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
                    CROSS[lp][REAL]=(this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                    CROSS[lp][IMAG]=(-this->data_out_g[lp][REAL] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][REAL]) * icorr2;
#endif
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update
                        aux2[kmod_g]  += Pk1;
#pragma omp atomic update
                        aux4[kmod_g]  += Pk2;
#pragma omp atomic update
                        mod_g[kmod_g]++ ;
                      }
                  }
                //  *****************************************
                if(j>0  && k>0)
                  {
                    lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                    real_prec Pk= (this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                    real_prec Pk1= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;
                    real_prec Pk2= (this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
                        CROSS[lp][REAL]=(this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                        CROSS[lp][IMAG]=(-this->data_out_g[lp][REAL] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][REAL]) * icorr2;
#endif
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic  update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic  update
                        aux2[kmod_g]  += Pk1;
#pragma omp atomic  update
                        aux4[kmod_g]  += Pk2;

#pragma omp atomic update
                        mod_g[kmod_g]++ ;
                      }

                  }
               // ******************************************
                // Add negative frequencies in x:
                // ******************************************
                if(i>0  && (j>0 || k>0))
                  {
                    lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                    real_prec Pk= (this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG])* icorr2;
                    real_prec Pk1= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;
                    real_prec Pk2= (this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
                    CROSS[lp][REAL]=(this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                    CROSS[lp][IMAG]=(-this->data_out_g[lp][REAL] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][REAL]) * icorr2;
#endif
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update
                        aux2[kmod_g]  += Pk1;
#pragma omp atomic update
                        aux4[kmod_g]  += Pk2;
#pragma omp atomic update
                        mod_g[kmod_g]++;
                      }
                  }
                // ******************************************
                //  Add negative frequencies in x and y:
                // ******************************************
                if(i>0  && j>0  && k>0)
                  {
                    lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                    real_prec Pk= (this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                    real_prec Pk1= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;
                    real_prec Pk2= (this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
                    CROSS[lp][REAL]=(this->data_out_g[lp][REAL]*this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG]*this->data_out_gp[lp][IMAG]) * icorr2;
                    CROSS[lp][IMAG]=(-this->data_out_g[lp][REAL]*this->data_out_gp[lp][IMAG] + this->data_out_g[lp][IMAG]*this->data_out_gp[lp][REAL]) * icorr2;
#endif
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic update
                        powerk_g0[kmod_g] += Pk;
#pragma omp atomic update
                        aux2[kmod_g] += Pk1;
#pragma omp atomic update
                        aux4[kmod_g] += Pk2;
#pragma omp atomic update
                    mod_g[kmod_g] ++ ;
                      }
                  }
              }
          }
      } // END of parallel loop
  } // END of parallel regione

  // Subtract shot noise and normalize only he 1d power spectrum.
  // For the monopole:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<powerk_g0.size();++i)
    powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*sqrt(this->normal_power*this->normal_power_two))); // use this if corr is defined using shell averages of spectra
  if(corr_factor==true) // if correlation coefficient is requested,
  {

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<powerk_g0.size();++i)
    aux2[i]=(mod_g[i]== 0 ? 0 : aux2[i]/(static_cast<real_prec>(mod_g[i])*normal_power));
  // For the hexaecapole:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<powerk_g0.size();++i)
    aux4[i]=(mod_g[i]== 0 ? 0 : aux4[i]/(static_cast<real_prec>(mod_g[i])*normal_power_two));
  //      Get the cross power normalized:// use this if corr is defined using shell averages of spectra
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<powerk_g0.size();++i)
    powerk_g0[i]/=static_cast<real_prec>(sqrt(aux2[i]*aux4[i]));
    }
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
  real_prec deltak=this->params._d_DeltaK_data();
  vector<real_prec> coords(this->params._Nft(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._Nft() ;++i)
      coords[i]=deltak*(i<=this->params._Nft()/2? static_cast<real_prec>(i): -static_cast<real_prec>(this->params._Nft()-i));

     for(ULONG i=0; i< this->params._Nft();++i)
        for(ULONG j=0; j< this->params._Nft();++j)
          for(ULONG k=0; k< this->params._Nft()/2+1;++k)
           {
               ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               real_prec kv=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2));
               ULONG kmod=static_cast<ULONG>(floor(kv));
               real_prec acc2=(aux2[kmod]==0? 1. :aux2[kmod]);
               real_prec acc4=(aux4[kmod]==0? 1. :aux4[kmod]);
               CROSS[lp][IMAG]/=sqrt(acc2*acc4);
               CROSS[lp][REAL]/=sqrt(acc2*acc4);
          }
   do_fftw_c2r(this->params._Nft(), CROSS,corr);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for (ULONG index=0;index<corr.size();++index)
       corr[index]=pow(fabs(corr[index]),0.1);
   real_prec xmin=get_min(corr);
   real_prec xmax=get_max(corr);
    ULONG nb=30;
    real_prec del=2./nb;
   vector<real_prec>hiss(nb,0);;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for (ULONG index=0;index<corr.size();++index)
       corr[index]=(2.0)*(corr[index]-xmin)/(xmax-xmin)-1.0;
   for (ULONG index=0;index<corr.size();++index)
       hiss[get_bin(corr[index],-1.0,nb,del,false)]++;
   string ff="hist.txt";
   ofstream sa;sa.open(ff.c_str());
   for (ULONG index=0;index<nb;++index)
      sa<<-1.0+(index+0.5)*del<<"  "<<hiss[index]<<endl;
    sa.close();
#endif
  aux2.clear();aux2.shrink_to_fit();
   aux4.clear();aux4.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_spectrum_fkp(vector<real_prec> & powerk_g0,
                                       vector<real_prec> & powerk_g2,
                                       vector<real_prec> & powerk_g4,
                                       vector<real_prec> & powerk_r,
                                       vector< vector<real_prec> >&power2d_cart,
                                       vector< vector<real_prec> >&power2d_spher,
                                       vector<int> & mod_g)
{
    this->So.enter(__PRETTY_FUNCTION__);
    int NTHREADS = _NTHREADS_;

#ifdef _WRITE_2DPOWER_
  vector < vector<int> > mod_cart(power2d_cart.size(), vector<int>(power2d_cart[0].size(),0));
  vector < vector<int> > mod_spher(power2d_spher.size(), vector<int>(power2d_spher[0].size(),0));
#endif
  vector<int> mod_r(powerk_r.size(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<mod_g.size();++i)
    {
      mod_g[i]=0;
      powerk_g0[i]=0;
#ifdef _WRITE_MULTIPOLES_
      powerk_g2[i]=0;
      powerk_g4[i]=0;
#endif
    }
  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
  int my_s_box_ave;
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear")
  {
    my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
    rDeltaK_window = 1.0 / this->params._d_DeltaK_window();
  }
  else if(this->params._type_of_binning()=="log")
    {
      my_s_box_ave = s_box_log;
      rkmin = 1.0/this->params._d_kmin();
    }
  vector<real_prec>yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>zz_MAS_array(this->params._Nft()/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>sn_zz_MAS_array(this->params._Nft()/2+1,0);
  real_prec M_PI_rNft = M_PI * this->rNft;
  real_prec SN_aux=0;
  if(true==this->params._SN_correction())
    SN_aux=this->shot_noise;
  real_prec SN_aux2=0;
  if(true==this->params._SN_correction())
    SN_aux2=this->shot_noise2;
#ifdef _USE_OMP_
#pragma omp parallel num_threads(_NTHREADS_)
#endif
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;
#ifdef _USE_OMP_
#pragma omp sections
#endif
    {
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->params._Nft()/2;++j)
          {
            yy_MAS = j * M_PI_rNft;
            yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
            sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
          }
      }
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->params._Nft()/2;++k)
          {  // loop over octant kz>0: this is all FFTW gives.
            zz_MAS = k * M_PI_rNft;
            zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
            sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
          }
      }
    }
#ifdef _USE_OMP_
#pragma omp barrier
#endif
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
    for(int i=0;i<this->params._Nft()/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
        if(i > 0)
          {
            xx_MAS = i * M_PI_rNft;
            xx_MAS = sin(xx_MAS)/xx_MAS;
          }
        else
          xx_MAS = 1.0;
        sn_xx_MAS = sin(i * M_PI_rNft);
        i_deltak_x_2 = i * i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
        i_per_fact = i_deltak_x_2;
        if(my_s_box_ave == s_box_log)
          i_deltak_x_2 *= rkmin * rkmin;
        for(int j=0;j<this->params._Nft()/2;++j)
          {  // loop over octant ky>0. Half of what FFT gives
            yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];
            j_deltak_y_2 = j * j * this->params._d_deltak_y() * this->params._d_deltak_y();
            j_per_fact = i_per_fact + j_deltak_y_2;
            j_per_fact = sqrt(j_per_fact);
            if(my_s_box_ave == s_box_log)
              j_deltak_y_2 *= rkmin * rkmin;
            for(int k=0;k<=this->params._Nft()/2;++k)
              {  // loop over octant kz>0: this is all FFTW gives.
                zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];
                k_deltak_z_2 = k * k * this->params._d_deltak_z() * this->params._d_deltak_z();
                if(my_s_box_ave == s_box_log)
                  k_deltak_z_2 *= rkmin * rkmin;
                real_prec kv;
                int kmod_g, kmod_r;
                // Compute k-shell index
                if(my_s_box_ave == s_box_linear)
                  {
                    kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
                    kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                    kmod_r=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_window)));
                  }
                else
                  {
                    if(my_s_box_ave == s_box_log)
                      {
                        kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
                        if(kv!=0)
                          {
                            kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                            kmod_r=kmod_g;
                          }
                        else
                          {
                            kmod_g=0;kmod_r=0;
                          }
                      }
                  }
                // Compute correction for MAS:
                real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
                real_prec sn_corr=SN_aux*SN_correction_MAS(this->params._mass_assignment(), sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);
                // Compute index in c-order
                ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                // Compute thisgs related to 2d power spectrum and multipole
                // decomposition for the first octant
                // Define the angle bewtween k and los:
                // In the FKP we need to specify the LOS direction. set to z by default.
#if defined (_WRITE_2DPOWER_) || defined (_WRITE_MULTIPOLES_)
                real_prec mu=0;
                switch(LOS){
                case(I_LOSX):  mu = i*this->params._d_deltak_x()/kv; break;
                case(I_LOSY):  mu = j*this->params._d_deltak_y()/kv; break;
                case(I_LOSZ):  mu = k*this->params._d_deltak_z()/kv; break;
                }
#endif

#ifdef _WRITE_2DPOWER
                // Compute mu = kz/k.
                int i_par=0;
                switch(LOS)
                  {
                  case(I_LOSX):i_par = static_cast<int>(floor((float)(this->params._d_deltak_x()*i * rDeltaK_data))); break;
                  case(I_LOSY):i_par = static_cast<int>(floor((float)(this->params._d_deltak_y()*j * rDeltaK_data))); break;
                  case(I_LOSZ):i_par = static_cast<int>(floor((float)(this->params._d_deltak_z()*k * rDeltaK_data))); break;
                  }
                int i_per = static_cast<int>(floor((float)(j_per_fact * rDeltaK_data)));
#endif
#ifdef _WRITE_MULTIPOLES_
                // Compute index in mu-direction
                real_prec i_mu  = static_cast<int>(floor((float)((mu-(-1.0))/this->params._d_Deltamu())));
                i_mu  = (i_mu==this->params._N_mu_bins() ? 0: i_mu);
                // Compute Legendre functions neeeded for the Multipole decomposition:
                real_prec leg2=0.5*(3.*mu*mu-1);
                real_prec leg4=(35.*pow(mu,4)-30.*mu*mu+3.)/8.;
#endif
                //  *****************************************
                // Now, add modes in the first octant
                // Note that k>0 from here
                // Note also that Leg are invariant under the
                // shift os sigh kz-> - kz, and therefore we do not
                // need explicitely to account for the kz<0 region;
                // We can simply multiply by two both in the power as
                // in the number of counted modes. This factor two cancells
                // in the end, and therefore we don't employ it.
                //  *****************************************
                //Exclude the zero-frequency onLY when computing the P: for the W keep it:
                if(lp!=0)
                  {
                    real_prec Pk, Pk1, Pk2;
                    if(true==this->params._measure_cross())
                      {
                        Pk= (this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                        Pk1= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2;
                        Pk2= (this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                      }
                    else if(false==this->params._measure_cross())
                      {
                        Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                        Pk2= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2  ;//this is for l=2 and l=4, not corrected for shot-noise
                      }
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
                        if(true==this->params._measure_cross())
                          {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g0[kmod_g]  += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g2[kmod_g]  += Pk1;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g4[kmod_g]  += Pk2;
                        }
                        else
                          {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g2[kmod_g]  += leg2*Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g4[kmod_g]  += leg4*Pk2;
#endif
                          }
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_g[kmod_g]++ ;
                      }

#ifdef _WRITE_2DPOWER_
                    // 2d in Spherical coordinates:
                    if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        power2d_spher[i_mu][kmod_g]+= Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_spher[i_mu][kmod_g]++;
                      }

                    // 2d in Cartesian Coordinates:
                    if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        power2d_cart[i_par][i_per] += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_cart[i_par][i_per] ++ ;
                      }
#endif
                  }
                else
                  {  //Compute 1d spherical average for the window function only
                    if(this->params._use_random_catalog())
                      {
                        if(kmod_r<powerk_r.size()){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_r[kmod_r]  += (data_out_r[lp][IMAG] * data_out_r[lp][IMAG] + data_out_r[lp][REAL]*data_out_r[lp][REAL]) * icorr2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          mod_r[kmod_r]++ ;
                        }
                      }
                  }
                //  *****************************************
                //  Now, add modes in other octants:
                //  *****************************************
                //  Add negative frequencies in y:
                //  kv is the same, for ky->-ky.
                //  The same applies below for kx.
                //  Therefore, kmod_g is also the same.
                //  We only compute lp again to go to the right octant
                //  *****************************************
                if(j>0  && k>0)
                  {
                    lp=ijk(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);   // Go to the ky<0 octant

                    real_prec Pk, Pk1,Pk2;
                    if(true==this->params._measure_cross())
                      {
                        Pk= (this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                        Pk1= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;
                        Pk2= (this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL] + this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG]) * icorr2;
                        //Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
                      }
                    else if(false==this->params._measure_cross())
                      {
                        Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                        Pk2= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2  ;//this is for l=2 and l=4, not corrected for shot-noise
                      }
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
                        if(true==this->params._measure_cross())
                          {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g0[kmod_g]  += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g2[kmod_g]  += Pk1;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g4[kmod_g]  += Pk2;
                          }
                        else if(false==this->params._measure_cross())
                          {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g2[kmod_g]  += leg2*Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                            powerk_g4[kmod_g]  += leg4*Pk2;
#endif
                          }
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_g[kmod_g]++ ;
                      }
#ifdef _WRITE_2DPOWER_
                    // 2d in Spherical coordinates:
                    if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        power2d_spher[i_mu][kmod_g]+= Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_spher[i_mu][kmod_g]++;
                      }
                    // 2d in Cartesian Coordinates:
                    if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        power2d_cart[i_par][i_per] += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_cart[i_par][i_per] ++ ;
                      }
#endif
                    //Compute 1d spherical average for the window function only
                    if(this->params._use_random_catalog()){
                      if(kmod_r<powerk_r.size()){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_r[kmod_r]  += (data_out_r[lp][REAL]*data_out_r[lp][REAL] + data_out_r[lp][IMAG]*data_out_r[lp][IMAG]) * icorr2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_r[kmod_r]++ ;
                      }
                    }
                  }
                // ******************************************
                // Add negative frequencies in x:
                // ******************************************
                if(i>0  && (j>0 || k>0)){
                  lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);

                  real_prec Pk,Pk1,Pk2;
                  if(true==this->params._measure_cross())
                    {
                      Pk= (this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                      Pk1= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2;
                      Pk2= (this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                      //Pk=Pk/sqrt((Pk1-SN_aux)*(Pk2-SN_aux2));
                    }
                  else if(false==this->params._measure_cross()){
                      Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                      Pk2= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2  ;//this is for l=2 and l=4, not corrected for shot-noise
                  }
                  if(kmod_g<powerk_g0.size())
                    {
                      if(true==this->params._measure_cross())
                        {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g0[kmod_g]  += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g2[kmod_g]  += Pk1;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g4[kmod_g]  += Pk2;
                        }
                      else if(false==this->params._measure_cross()){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_

#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g2[kmod_g]  += leg2*Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g4[kmod_g]  += leg4*Pk2;
#endif
                      }
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_g[kmod_g]++;
                    }

#ifdef _WRITE_2DPOWER_
                  // 2d in Spherical coordinates:
                  if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      power2d_spher[i_mu][kmod_g]+= Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_spher[i_mu][kmod_g]++;
                    }
                  // 2d in Cartesian Coordinates:
                  if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      power2d_cart[i_par][i_per] += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_cart[i_par][i_per] ++ ;
                    }
#endif
                  //Compute 1d spherical average for the window function only
                  if(this->params._use_random_catalog())
                    {
                      if(kmod_r<powerk_r.size())
                        {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_r[kmod_r]  += (data_out_r[lp][IMAG]*data_out_r[lp][IMAG]+ data_out_r[lp][REAL]*data_out_r[lp][REAL]) * icorr2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          mod_r[kmod_r]++ ;
                        }
                    }
                }

                // ******************************************
                //  Add negative frequencies in x and y:
                // ******************************************
                if(i>0  && j>0  && k>0){
                  lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                  real_prec Pk,Pk1,Pk2;
                  if(true==this->params._measure_cross())
                   {
                     Pk= (this->data_out_g[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                     Pk1= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2;
                     Pk2= (this->data_out_gp[lp][IMAG] * this->data_out_gp[lp][IMAG] + this->data_out_gp[lp][REAL] * this->data_out_gp[lp][REAL]) * icorr2;
                   }
                  else if(false==this->params._measure_cross())
                  {
                    Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                    Pk2= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]) * icorr2  ;//this is for l=2 and l=4, not corrected for shot-noise
                  }
                  if(kmod_g<powerk_g0.size())
                    {
                      if(true==this->params._measure_cross())
                        {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g0[kmod_g] += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g2[kmod_g] += Pk1;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g4[kmod_g] += Pk2;

                        }
                      else if(false==this->params._measure_cross()){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g0[kmod_g] += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g2[kmod_g] += leg2*Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g4[kmod_g] += leg4*Pk2;
#endif

                      }
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_g[kmod_g] ++ ;
                    }

#ifdef _WRITE_2DPOWER_
                  // 2d in Spherical coordinates:
                  if(i_mu<power2d_spher.size() && kmod_g<power2d_spher[0].size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      power2d_spher[i_mu][kmod_g]+= Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_spher[i_mu][kmod_g]++;
                    }
                  // 2d in Cartesian Coordinates:
                  if(i_par<power2d_cart.size() && i_per<power2d_cart[0].size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      power2d_cart[i_par][i_per] += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_cart[i_par][i_per] ++ ;
                    }
#endif
                  //Compute 1d spherical average for the window function only

                  if(this->params._use_random_catalog())
                    {
                      if(kmod_r<powerk_r.size())
                        {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_r[kmod_r] += (data_out_r[lp][REAL]*data_out_r[lp][REAL]+data_out_r[lp][IMAG]*data_out_r[lp][IMAG]) * icorr2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          mod_r[kmod_r] ++ ;
                        }
                    }
                }
              }
          }
      } // END of parallel loop

  } // END of parallel regione
  // Subtract shot noise and normalize only he 1d power spectrum.
  // For the monopole:
  if(true==this->params._measure_cross()) // use the arrays g2 and g4 to get auto power needed for the correlation
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g0.size();++i)
        powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power)); // use this if corr is defined using shell averages of spectra

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g2.size();++i)
        powerk_g2[i]=(mod_g[i]== 0 ? 0 : powerk_g2[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

      // For the hexaecapole:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g4.size();++i)
        powerk_g4[i]=(mod_g[i]== 0 ? 0 : powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

      //      Get the cross power normalized:// use this if corr is defined using shell averages of spectra
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g0.size();++i)
        powerk_g0[i]/=static_cast<real_prec>(sqrt(powerk_g2[i]*powerk_g4[i]));
    }

  else
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g0.size();++i)
        powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g2.size();++i)
        powerk_g2[i]=(mod_g[i]== 0 ? 0 : 5.*powerk_g2[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

      // For the hexaecapole:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g4.size();++i)
        powerk_g4[i]=(mod_g[i]== 0 ? 0 : 9.*powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));
#endif
    }

  // For the quadrupole:
  // Also for the window function:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<powerk_r.size();++i)
    powerk_r[i]=(mod_r[i]== 0 ? 0 : powerk_r[i]/(static_cast<real_prec>(mod_r[i])*this->normal_window));


#ifdef _WRITE_2DPOWER_
  // And for the 2d estimates in cartesian coordinates
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(int i=0;i<power2d_cart.size();++i)
    for(int j=0;j<power2d_cart[0].size();++j)
      power2d_cart[i][j]=(mod_cart[i][j]== 0? 0. : power2d_cart[i][j]/((static_cast<real_prec>(mod_cart[i][j])*normal_power)-SN_aux));

#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(int i=0;i<power2d_spher.size();++i)
    for(int j=0;j<power2d_spher[0].size();++j)
      power2d_spher[i][j]=(mod_spher[i][j]==0? 0. : power2d_spher[i][j]/((static_cast<real_prec>((mod_spher[i][j]))*normal_power)-SN_aux));
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_spectrum(vector<real_prec> & powerk_g0,vector<int> & mod_g)
{
    int NTHREADS = _NTHREADS_;
    do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<mod_g.size();++i)
    {
      mod_g[i]=0;
      powerk_g0[i]=0;
    }
  real_prec rDeltaK_data;
  real_prec rkmin;
  int my_s_box_ave;
  my_s_box_ave = s_box_other;

  if(this->params._type_of_binning()=="linear") {
    my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
  }
  else if(this->params._type_of_binning()=="log")
    {
      my_s_box_ave = s_box_log;
      rkmin = 1.0/this->params._d_kmin();
    }
  vector<real_prec>yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>zz_MAS_array(this->params._Nft()/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>sn_zz_MAS_array(this->params._Nft()/2+1,0);

  real_prec M_PI_rNft = M_PI * this->rNft;

  real_prec SN_aux=0;
  if(true==this->params._SN_correction())
    SN_aux=this->shot_noise;

#pragma omp parallel num_threads(_NTHREADS_)
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;

#pragma omp sections
    {
#pragma omp section
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->params._Nft()/2;++j)
          {
            yy_MAS = j * M_PI_rNft;
            yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
            sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
          }
      }
#pragma omp section
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->params._Nft()/2;++k)
          {  // loop over octant kz>0: this is all FFTW gives.
            zz_MAS = k * M_PI_rNft;
            zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
            sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
          }
      }
    }
#pragma omp barrier
#pragma omp for nowait
    for(int i=0;i<this->params._Nft()/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
        if(i > 0)
          {
            xx_MAS = i * M_PI_rNft;
            xx_MAS = sin(xx_MAS)/xx_MAS;
          }
        else
          xx_MAS = 1.0;

        sn_xx_MAS = sin(i * M_PI_rNft);
        i_deltak_x_2 = i * i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
        i_per_fact = i_deltak_x_2;
        if(my_s_box_ave == s_box_log)
          i_deltak_x_2 *= rkmin * rkmin;

        for(int j=0;j<this->params._Nft()/2;++j)
          {  // loop over octant ky>0. Half of what FFT gives
            yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];

            j_deltak_y_2 = j * j * this->params._d_deltak_y() * this->params._d_deltak_y();
            j_per_fact = i_per_fact + j_deltak_y_2;
            j_per_fact = sqrt(j_per_fact);
            if(my_s_box_ave == s_box_log)
              j_deltak_y_2 *= rkmin * rkmin;

            for(int k=0;k<=this->params._Nft()/2;++k)
              {  // loop over octant kz>0: this is all FFTW gives.
                zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];

                k_deltak_z_2 = k * k * this->params._d_deltak_z() * this->params._d_deltak_z();
                if(my_s_box_ave == s_box_log)
                  k_deltak_z_2 *= rkmin * rkmin;

                real_prec kv;
                int kmod_g, kmod_r;

                // Compute k-shell index
                if(my_s_box_ave == s_box_linear)
                  {
                    kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
                    kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                  }
                else
                  {
                    if(my_s_box_ave == s_box_log)
                      {
                        kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
                        if(kv!=0)
                           kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                        else
                            kmod_g=0;
                      }
                  }


                // Compute correction for MAS:
                real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
                real_prec sn_corr=SN_aux*SN_correction_MAS(this->params._mass_assignment(), sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);

                // Compute index in c-order
                int lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);

                // Compute thisgs related to 2d power spectrum and multipole
                // decomposition for the first octant
                // Define the angle bewtween k and los:
                // In the FKP we need to specify the LOS direction. set to z by default.
                //Exclude the zero-frequency onLY when computing the P: for the W keep it:
                if(lp!=0)
                  {
                    real_prec Pk, Pk1, Pk2;
                    real_prec pk_raw=this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL];
                    Pk= (pk_raw - sn_corr*this->normal_power) * icorr2  ;

                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update
                        mod_g[kmod_g]++ ;
                      }
                  }
                if(j>0  && k>0)
                  {
                    lp=ijk(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);   // Go to the ky<0 octant
                    real_prec pk_raw=this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL];
                    real_prec Pk= (pk_raw - sn_corr*this->normal_power) * icorr2  ;
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#pragma omp atomic  update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update
                        mod_g[kmod_g]++ ;
                      }
                  }
                // ******************************************
                // Add negative frequencies in x:
                // ******************************************
                if(i>0  && (j>0 || k>0)){
                  lp=ijk(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                  real_prec pk_raw=this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL];
                  real_prec Pk= (pk_raw - sn_corr*this->normal_power) * icorr2  ;
                  if(kmod_g<powerk_g0.size())
                    {

#pragma omp atomic update
                        powerk_g0[kmod_g]  += Pk;
#pragma omp atomic update
                      mod_g[kmod_g]++;
                    }
                }

                // ******************************************
                //  Add negative frequencies in x and y:
                // ******************************************
                if(i>0  && j>0  && k>0){
                  lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                  real_prec pk_raw= this->data_out_g[lp][IMAG] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][REAL] ;
                  real_prec Pk= (pk_raw - sn_corr*this->normal_power) * icorr2  ;
                  if(kmod_g<powerk_g0.size())
                    {
#pragma omp atomic update
                        powerk_g0[kmod_g] += Pk;
#pragma omp atomic update
                      mod_g[kmod_g] ++ ;
                    }


                }
              }
          }
      } // END of parallel loop

  } // END of parallel regione
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<powerk_g0.size();++i)
        powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_spectrum_yamamoto_new(vector<real_prec> & powerk_g0,
                                       vector<real_prec> & powerk_g2,
                                       vector<real_prec> & powerk_g4,
                                       vector<int> & mod_g)
{
#ifdef _VERBOSE_POWER_
    this->So.enter(__PRETTY_FUNCTION__);
#endif

   int NTHREADS = _NTHREADS_;
   int nbink=mod_g.size();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<mod_g.size();++i)
    {
      mod_g[i]=0;
      powerk_g0[i]=0;
      powerk_g2[i]=0;
      powerk_g4[i]=0;
    }
  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
  int my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear") {
    my_s_box_ave = s_box_linear;
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
    rDeltaK_window = 1.0 / this->params._d_DeltaK_window();
  }
  else if(this->params._type_of_binning()=="log")
    {
      my_s_box_ave = s_box_log;
      rkmin = 1.0/this->params._d_kmin();
    }

  vector<real_prec>yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>zz_MAS_array(this->params._Nft()/2+1,0);
  vector<real_prec>sn_yy_MAS_array(this->params._Nft()/2,0);
  vector<real_prec>sn_zz_MAS_array(this->params._Nft()/2+1,0);

  real_prec M_PI_rNft = M_PI * this->rNft;

  real_prec SN_aux=0;
  if(true==this->params._SN_correction())
    SN_aux=this->shot_noise;

  real_prec SN_aux2=0;
  if(true==this->params._SN_correction())
    SN_aux2=this->shot_noise2;

#ifdef _USE_OMP_
#pragma omp parallel num_threads(_NTHREADS_)
#endif
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;
#ifdef _USE_OMP_
#pragma omp sections
#endif
    {
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        yy_MAS_array[0] = 1.0;
        sn_yy_MAS_array[0] = 0.0;
        for(int j=1;j<this->params._Nft()/2;++j)
          {
            yy_MAS = j * M_PI_rNft;
            yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
            sn_yy_MAS_array[j] = sin(j * M_PI_rNft);
          }
      }
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        zz_MAS_array[0] = 1.0;
        sn_zz_MAS_array[0] = 0.0;
        for(int k=1;k<=this->params._Nft()/2;++k)
          {  // loop over octant kz>0: this is all FFTW gives.
            zz_MAS = k * M_PI_rNft;
            zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
            sn_zz_MAS_array[k] = sin(k * M_PI_rNft);
          }
      }
    }
#ifdef _USE_OMP_
#pragma omp barrier
#endif
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
    for(int i=0;i<this->params._Nft()/2;++i)
      {  // loop over octant kx>0. Half of what FFT gives
        if(i > 0)
          {
            xx_MAS = i * M_PI_rNft;
            xx_MAS = sin(xx_MAS)/xx_MAS;
          }
        else
          xx_MAS = 1.0;
        sn_xx_MAS = sin(i * M_PI_rNft);
        i_deltak_x_2 = i * i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
        if(my_s_box_ave == s_box_log)
          i_deltak_x_2 *= rkmin * rkmin;
        for(int j=0;j<this->params._Nft()/2;++j)
          {  // loop over octant ky>0. Half of what FFT gives
            yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];
            j_deltak_y_2 = j * j * this->params._d_deltak_y() * this->params._d_deltak_y();
            if(my_s_box_ave == s_box_log)
              j_deltak_y_2 *= rkmin * rkmin;
            for(int k=0;k<=this->params._Nft()/2;++k)
              {  // loop over octant kz>0: this is all FFTW gives.
                zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];

                real_prec kz=k*this->params._d_deltak_z();

                k_deltak_z_2 = kz * kz ;
                if(my_s_box_ave == s_box_log)
                  k_deltak_z_2 *= rkmin * rkmin;

                real_prec kv=0;
                int kmod_g=0;

                // Compute correction for MAS:
                real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
                real_prec sn_corr=SN_aux*SN_correction_MAS(this->params._mass_assignment(), sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);

                // Compute index in c-order
                ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);

               //  *****************************************
                // Use informtion from the first octant:
                if(lp!=0)
                  {
                    real_prec kx=i*this->params._d_deltak_x();
                    real_prec ky=j*this->params._d_deltak_x();
                    kv=sqrt(kx*kx + ky*ky + kz*kz);
                    if(my_s_box_ave == s_box_linear)
                      kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                    else
                      if(my_s_box_ave == s_box_log)
                        if(kv!=0)
                          kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;

                    real_prec Pk2, Pk4;
                    real_prec Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                    aux_func_yamamoto(lp,kx,ky, kz,icorr2,Pk2,Pk4);
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g2[kmod_g]  += Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g4[kmod_g]  += Pk4;
#endif
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_g[kmod_g]++ ;
                      }
                  }
                //  *****************************************
                //  Now, add modes in other octants:
                //  Add negative frequencies in y:
                //  The same applies below for kx.
                //  Therefore, kmod_g is also the same.
                //  We only compute lp again to go to the right octant
                //  *****************************************
                if(j>0  && k>0)
                  {
                    real_prec kx=static_cast<real_prec>(i)*this->params._d_deltak_x();
                    real_prec ky=-static_cast<real_prec>(j)*this->params._d_deltak_x();
                    kv=sqrt(kx*kx+ ky*ky + kz*kz);
                    if(my_s_box_ave == s_box_linear)
                        kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                    else if(my_s_box_ave == s_box_log)
                      if(kv!=0)
                         kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                    lp=ijk(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);   // Go to the ky<0 octant
                    real_prec Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                    real_prec Pk2, Pk4;
                    aux_func_yamamoto(lp,kx,ky,kz,icorr2,Pk2,Pk4);
                    // 1d Spherical average
                    if(kmod_g<powerk_g0.size())
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g2[kmod_g]  += Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                          powerk_g4[kmod_g]  += Pk4;
#endif
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_g[kmod_g]++ ;
                      }
                  }
                // ******************************************
                // Add negative frequencies in x:
                // ******************************************
                if(i>0  && (j>0 || k>0))
                {
                    real_prec kx=-static_cast<real_prec>(i)*this->params._d_deltak_x();
                    real_prec ky=static_cast<real_prec>(j)*this->params._d_deltak_x();
                    kv=sqrt(kx*kx+ky*ky+kz*kz);
                    if(my_s_box_ave == s_box_linear)
                      {
                        kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                      }
                    else
                      {
                        if(my_s_box_ave == s_box_log)
                          if(kv!=0)
                            kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                      }

                  lp=ijk(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                  real_prec Pk2, Pk4;
                  real_prec Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                  aux_func_yamamoto(lp,kx, ky ,kz,icorr2,Pk2,Pk4);
                  // 1d Spherical average
                  if(kmod_g<powerk_g0.size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g2[kmod_g]  += Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g4[kmod_g]  += Pk4;
#endif
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_g[kmod_g]++ ;
                    }
                }

                // ******************************************
                //  Add negative frequencies in x and y:
                // ******************************************
                if(i>0  && j>0  && k>0)
                {
                    real_prec kx=-static_cast<real_prec>(i)*this->params._d_deltak_x();
                    real_prec ky=-static_cast<real_prec>(j)*this->params._d_deltak_y();
                    kv=sqrt(kx*kx+ ky*ky + kz*kz);
                    if(my_s_box_ave == s_box_linear)
                      {
                        kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                      }
                    else
                      {
                        if(my_s_box_ave == s_box_log)
                          if(kv!=0)
                            kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                      }

                  lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                  real_prec Pk= (this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG] + this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] - sn_corr*normal_power) * icorr2  ;
                  real_prec Pk2, Pk4;
                  aux_func_yamamoto(lp, kx,ky,kz,icorr2,Pk2,Pk4);
                  // 1d Spherical average
                  if(kmod_g<powerk_g0.size())
                    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g0[kmod_g]  += Pk;
#ifdef _WRITE_MULTIPOLES_
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g2[kmod_g]  += Pk2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g4[kmod_g]  += Pk4;
#endif
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                      mod_g[kmod_g]++ ;
                    }
                 }

             }
          }

    } // END of parallel loop
  } // END of parallel regione
  // Subtract shot noise and normalize only he 1d power spectrum.
  // For the monopole:
  for(int i=0;i<nbink;++i)
    powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));
  // For the quadrupole: no shot noise correction applies
  for(int i=0;i<nbink;++i)
    powerk_g2[i]=(mod_g[i]== 0 ? 0 : 5.*powerk_g2[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));
  if(this->params._statistics()=="Pk_ys" ||  this->params._statistics()=="Pk_ysc"){
    // For the hexadecapole: no shot noise correction applies. Scoccimarro recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : (35./2.)*powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*normal_power)- powerk_g2[i] - (7./2.)*(powerk_g0[i]+SN_aux));
  }
  else if(this->params._statistics()=="Pk_yb"||  this->params._statistics()=="Pk_ybc" ){
    //For the hexadecapole: no shot noise correction applies. Bianchi et al. recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : 9.0*powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*normal_power));
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_spectrum_yamamoto(vector<real_prec> & powerk_g0,
                                            vector<real_prec> & powerk_g2,
                                            vector<real_prec> & powerk_g4,
                                            vector<int> & mod_g)
{
  // Power spectrum and multipole decomposition
  // Using the Yamamoto estimator with the algorithm of Scocimarro
  // that permits the use of FFTW.
  // Loops over all Available modes in Fourier space
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#else
  int NTHREADS =1;
#endif
  int nbink=mod_g.size();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<nbink;++i)mod_g[i]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<nbink;++i)powerk_g0[i]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<nbink;++i)powerk_g2[i]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<nbink;++i)powerk_g4[i]=0;
  int my_s_box_ave;
  real_prec rkmin, rDeltaK_data, rDeltaK_window;
  my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear")
    my_s_box_ave = s_box_linear;
  else if(this->params._type_of_binning()=="log")
    my_s_box_ave = s_box_log;
  rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
  rDeltaK_window = 1.0 / this->params._d_DeltaK_window();
  rkmin = 1.0/this->params._d_kmin();
  real_prec M_PI_rNft = M_PI * this->rNft;
  vector<real_prec>yy_MAS_array(this->params._Nft()+1,0);
  vector<real_prec>zz_MAS_array(this->params._Nft()+1,0);
  vector<real_prec>_ky_array(this->params._Nft()+1,0);
  vector<real_prec>_ky2_array(this->params._Nft()+1,0);
  vector<real_prec>sn_yy_MAS_array(this->params._Nft()+1,0);
  vector<real_prec>sn_zz_MAS_array(this->params._Nft()+1,0);
  real_prec SN_aux;
  if(true==this->params._SN_correction())
    SN_aux=shot_noise;
  else
    SN_aux=0;

#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
  {
#endif
      real_prec _kx, _ky, kz;
    real_prec _kx2, _ky2, kz2;
    real_prec _kx_kmin2, _ky_kmin2;

    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec rxx_MAS, ryy_MAS, rzz_MAS;
    real_prec sn_xx_MAS, sn_yy_MAS, sn_zz_MAS;

    int q_i, q_j, q_k;

#ifdef _USE_OMP_
#pragma omp sections
#endif
    {
#ifdef _USE_OMP_
#pragma omp section
#endif
      for(int j = 0; j <= this->params._Nft(); j++)
        {
          q_j= (j<=this->params._Nft()/2? j: j-(this->params._Nft()+1));
          _ky_array[j] = (_ky = q_j*this->params._d_deltak_y());
          _ky2_array[j] = (_ky2 = _ky * _ky);
        }

#ifdef _USE_OMP_
#pragma omp section
#endif
{
      yy_MAS_array[0] = 1.0;
      sn_yy_MAS_array[0] = 0.0;
      for(int j = 1; j<= this->params._Nft(); j++)
        {
          q_j= (j<=this->params._Nft()/2? j: j-(this->params._Nft()+1));
          yy_MAS = q_j * M_PI_rNft;
          yy_MAS_array[j] = (yy_MAS = sin(yy_MAS)/yy_MAS);
          sn_yy_MAS_array[j] = sin(q_j * M_PI_rNft);
        }
}
#ifdef _USE_OMP_
#pragma omp section
#endif
    {
      zz_MAS_array[0] = 1.0;
      sn_zz_MAS_array[0] = 0.0;
      for(int k = 1; k <= this->params._Nft(); k++)
        {
          q_k= (k<=this->params._Nft()/2? k: k-(this->params._Nft()+1));
          zz_MAS = q_k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
          sn_zz_MAS_array[k] = sin(q_k * M_PI_rNft);
        }

       }

    }

#ifdef _USE_OMP_
#pragma omp barrier
#endif


#ifdef _USE_OMP_
#pragma omp for nowait
#endif
    for(int i=0;i<=this->params._Nft();++i)
      {
        //set coordinate in units of the fundamental mode
        q_i= (i<=this->params._Nft()/2? i: i-(this->params._Nft()+1));
        _kx = q_i*this->params._d_deltak_x();
        _kx2 = _kx * _kx;
        if(my_s_box_ave == s_box_log)
          _kx_kmin2 = _kx2 * rkmin * rkmin;
        if(i > 0) {
          xx_MAS = q_i * M_PI_rNft;
          xx_MAS = sin(xx_MAS)/xx_MAS;
        }
        else xx_MAS = 1.0;

        sn_xx_MAS = sin(q_i * M_PI_rNft);

            for(int j=0;j<=this->params._Nft();++j)
          {
            //set coordinate in units of the fundamental mode:
            _ky = _ky_array[j];
            _ky2 = _ky2_array[j];
            yy_MAS = yy_MAS_array[j];
            sn_yy_MAS=sn_yy_MAS_array[j];

            if(my_s_box_ave == s_box_log)
              _ky_kmin2 = _ky2 * rkmin * rkmin;

            for(int k=0;k<=this->params._Nft();++k)
              {
                //set coordinate in units of the fundamental mode:
                q_k= (k<=this->params._Nft()/2? k: k-(this->params._Nft()+1));

                zz_MAS = zz_MAS_array[k];
                sn_zz_MAS=sn_zz_MAS_array[k];
                // Note k>=Nft/2+1 is the region in which kz<0 (q_k = k-(Nft-1))
                // Then, in order to apply Hermitian symmetry, we use
                // delta(kx,ky,-kz)=delta(-kx,-ky,kz)*
                // when we cover this region we reverse the sign of kx and ky
                // Since kz is reversed (see the definiton of q_k) in this region, we also reverse it here
                // to make it positive.
                // Note that the terms assocaited to the multipole decomposition
                // are products of k-coordinates: these products do not change
                // when we apply Hermitian symetri, so there is no need
                // to remap coordinates.
                int kmod_g=0, kmod_r=0;
                // Find the indices assocaited to the original fftw array
                int o_i, o_j, o_k;
                real_prec factor,factor2;
                remap(this->params._Nft(), this->params._Nft(), this->params._Nft(), i, j, k, &o_i, &o_j, &o_k, &factor);
                factor2 = factor * factor;
                // Recall that the value of factor is +1 if kz>0 and
                // -1 if kz<0. We could multiply the imaginary part of the array
                // by this factor such that Hermitian symmetry follows.
                // However, note that factor comes squared,
                // then will be always one.
                // The coordinates however need to be reversed in the kz<0 region,
                // again, to apply HS. This because we need explicitely need the coordinates
                // in this estimator, not only the modulos of the wavevector.
                // This factor in the end appears on power of 4 though, i.e, 1.
                real_prec kx = _kx * factor;
                real_prec ky = _ky * factor;
                real_prec kx2 = _kx2 * factor2;
                real_prec ky2 = _ky2 * factor2;

                kz = factor*q_k*this->params._d_deltak_z();
                kz2 = kz * kz;

                real_prec kv, kv2;
                if(my_s_box_ave == s_box_linear){
                  kv2 = kx2 + ky2 + kz2;
                  kv = sqrt(kv2);
                  kmod_g=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_data)));
                  kmod_r=static_cast<int> (floor(static_cast<float>((kv-this->params._d_kmin()) * rDeltaK_window)));

                }
                else{
                  if(my_s_box_ave == s_box_log)
                    {
                      real_prec kx_kmin2 = _kx_kmin2 * factor2;
                      real_prec ky_kmin2 = _ky_kmin2 * factor2;

                      kv2 = kx_kmin2 + ky_kmin2 + kz2 * rkmin * rkmin;
                      if(kv2!=0){
                          kmod_g=(int)floor(static_cast<float>(log10(sqrt(kv))/this->params._d_Deltal()))+1;
                          kmod_r=kmod_g;
                      }
                      else{kmod_g=0;kmod_r=0;}
                    }
                }


                // Compute the c-ordered index associated to the original fftw array:
                ULONG lp=ijk(o_i,o_j,o_k,  this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);

                // Compute correction for MAS:
                real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
                real_prec icorr2 = 1.0/(corr*corr);
                real_prec F2=0;

                if(lp!=0)
                  { //Exclude zero-mode
                    // Subtract shot noise and normalize only the power spectrum.

                    real_prec sn_corr=SN_aux*SN_correction_MAS(this->params._mass_assignment(), sn_xx_MAS, sn_yy_MAS, sn_zz_MAS);
                    real_prec Pk=  (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL]+ this->data_out_g[lp][IMAG]*this->data_out_g[lp][IMAG] - sn_corr*normal_power) * icorr2;
                    real_prec fluc1=0;
                    aux_func_yamamoto(lp,kx,ky,kz,icorr2, fluc1,F2);

                    // 1d Spherical average:
                    if(kmod_g<nbink)
                      {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g0[kmod_g]  += Pk;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g2[kmod_g]  += fluc1;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        powerk_g4[kmod_g]  += F2;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                        mod_g[kmod_g]++;
                      }

                  }//closers lp!=0

            } // closes loop over k
          }//closes loop over j
      }//closes loop over i
  }


  // For the monopole:
  for(int i=0;i<nbink;++i)
    powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

  // For the quadrupole: no shot noise correction applies
  for(int i=0;i<nbink;++i)
    powerk_g2[i]=(mod_g[i]== 0 ? 0 : 5.*powerk_g2[i]/(static_cast<real_prec>(mod_g[i])*this->normal_power));

  if(this->params._statistics()=="Pk_ys" ||  this->params._statistics()=="Pk_ysc"){
    // For the hexadecapole: no shot noise correction applies. Scoccimarro recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : (35./2.)*powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*normal_power)- powerk_g2[i] - (7./2.)*(powerk_g0[i]+SN_aux));
  }
  else if(this->params._statistics()=="Pk_yb"||  this->params._statistics()=="Pk_ybc" ){
    //For the hexadecapole: no shot noise correction applies. Bianchi et al. recipie.
    for(int i=0;i<nbink;++i)
      powerk_g4[i]=(mod_g[i]== 0 ? 0 : 9.0*powerk_g4[i]/(static_cast<real_prec>(mod_g[i])*normal_power));
  }
}

// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************


void FftwFunctions::get_fkp_error_bars(s_data_structure_direct_sum *s_data, vector<real_prec> &kvector, vector<real_prec> &pk, vector<int> &modes, vector<real_prec> &sigma)
{
  //////////////////////////////////////////////////////////
  // Estimates of the variance for the power spectrum
  // Based on the FKP estimator. Selects between the
  // exact and the approximate expression for the var(P)
  //////////////////////////////////////////////////////////
  this->So.enter(__PRETTY_FUNCTION__);


  if(true==this->params._use_random_catalog())
    {
      if(true==this->params._FKP_error_bars_exact())
        {//Use the exact FKP formula
          for(int i=0;i<SN.size();++i)
            this->SN[i]=this->alpha*(1+this->alpha)*SN[i]/this->normal_power;
          for(int i=0;i<Q.size(); ++i)
            Q[i] =this->alpha*Q[i]/this->normal_power;
          do_fkp_error_bars_exact(this->Q, this->SN, pk, sigma);
          for(int i=0;i<sigma.size();++i)
            sigma[i]=sqrt(2.*sigma[i])/(2.*modes[i]);
          // FKP Formula. Multiply the number of modes by a facgor 2 to account also
          // for those kz<0 that were not explictely counted in the shell average
        }
      else
        { //If not the exact, use the approximation for large scales, involving the definition of the effective volume
          vector<real_prec> veff(sigma.size(),0);
          do_fkp_error_bars_veff(s_data, pk,veff);
          for(int i=0;i<veff.size();++i)
            {
              real_prec deltaVk=4.*M_PI*this->params._d_DeltaK_data()*pow(kvector[i],2)*(1.+(1./12.)*pow(this->params._d_DeltaK_data()/kvector[i],2))/(pow(2.*M_PI, 3));
              sigma[i]=sqrt(2./(deltaVk*this->alpha*veff[i]));
            }
        }
    }
  else
    {
      for(int i=0;i<sigma.size();++i)
        {
          real_prec deltaVk=4.*M_PI*this->params._d_DeltaK_data()*pow(kvector[i],2)*(1.+(1./12.)*pow(this->params._d_DeltaK_data()/kvector[i],2))/(pow(2.*M_PI, 3));
          sigma[i]=sqrt(2./(this->params._Lbox()*this->params._Lbox()*this->params._Lbox()*deltaVk))*(pk[i]+1./s_data->mean_density);
          // This has been compared to the variance from N-Body simulations
          // giving a good agreement, specially on large scales.
        }
    }
}

// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************
// ****************************************************************************************************************************




void FftwFunctions::do_fkp_error_bars_exact(vector<real_prec>&Q, vector<real_prec> &SN, vector<real_prec> &P, vector<real_prec> &sigma)
{
  //////////////////////////////////////////////////////////
  // Estimation of the variance of the power spectrum using
  // the FKP estimator with the exact expression, equation
  // 2.4.6  of FKP paper
  //////////////////////////////////////////////////////////
  this->So.enter(__PRETTY_FUNCTION__);

  real_prec temp1,temp2,perp1,perp2;
  real_prec kv, kvp;
  int ik, ikp, lp;

  /*Transform Q and SN*/

#ifdef DOUBLE_PREC
  this->data_out_SN=(complex_prec *)fftw_malloc(this->params._NGRID_h()*sizeof(real_prec));
  this->data_out_Q =(complex_prec *)fftw_malloc(this->params._NGRID_h()*sizeof(real_prec));
#else
  this->data_out_SN=(complex_prec *)fftwf_malloc(this->params._NGRID_h()*sizeof(real_prec));
  this->data_out_Q =(complex_prec *)fftwf_malloc(this->params._NGRID_h()*sizeof(real_prec));
#endif


  do_fftw_r2c(this->params._Nft(),this->SN, this->data_out_SN);
  do_fftw_r2c(this->params._Nft(),this->Q, this->data_out_Q);

  time_t start;
  time (&start);
  real_prec full=pow(static_cast<real_prec>(this->params._Nft())/2, 2)*pow(static_cast<real_prec>(this->params._Nft())/2, 2)*pow(static_cast<real_prec>(this->params._Nft())/2+1, 2)-1;
  long counter=0;

  /*sum over k' */
  for(int i=0;i<this->params._Nft()/2;++i)
    {
      for(int j=0;j<this->params._Nft()/2;++j)
        {
          for(int k=0;k<this->params._Nft()/2+1;++k)
            {
              if(this->params._type_of_binning()=="linear"){
                kv=sqrt(pow(i*this->params._d_deltak_x(),2)+pow(j*this->params._d_deltak_y(),2)+pow(k*this->params._d_deltak_z(),2));
                ik=(int)floor((float)(kv/this->params._d_DeltaK_data()));
              }
              else{
                if(this->params._type_of_binning()=="log"){
                  perp1=pow(this->params._d_deltak_x()*i/(this->params._d_kmin()),2)+pow(this->params._d_deltak_y()*j/(this->params._d_kmin()),2)+pow(this->params._d_deltak_z()*k/(this->params._d_kmin()),2);
                  if(perp1!=0){
                    ik= (int)floor( (float) (log10(sqrt(perp1))/this->params._d_Deltal()))+1;
                  }
                  else ik=0;
                }
              }

              /*sum over k'' */
              for(int l=0;l<this->params._Nft()/2;l++){
                for(int m=0;m<this->params._Nft()/2;m++){
                  for(int n=0;n<this->params._Nft()/2;n++){
                    counter++;
                    //	      So.comp_time(start, full, counter);

                    if(this->params._type_of_binning()=="linear"){
                      kvp=sqrt(pow(l*this->params._d_deltak_x(),2)+pow(m*this->params._d_deltak_y(),2)+pow(n*this->params._d_deltak_z(),2));
                      ikp=(int)floor((float)(kvp/this->params._d_DeltaK_data()));
                    }
                    else{
                      if(this->params._type_of_binning()=="log"){
                        kvp=pow(this->params._d_deltak_x()*l/(this->params._d_kmin()),2)+pow(this->params._d_deltak_y()*m/(this->params._d_kmin()),2)+pow(this->params._d_deltak_z()*n/(this->params._d_kmin()),2);
                        if(kvp!=0){
                          ikp= static_cast<int>(floor( (float) (log10(sqrt(perp1))/this->params._d_Deltal())))+1;
                        }
                        else ikp=0;
                      }
                    }

                    if(ik == ikp){
                      int I=fabs(l-i);
                      int J=fabs(m-j);
                      int K=fabs(n-k);
                      lp=ijk(I,J,K,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                      if(lp!=0){
                        temp1 = P[ik]*data_out_Q[lp][REAL]+this->data_out_SN[lp][REAL];
                        temp2 = P[ik]*data_out_Q[lp][IMAG]+data_out_SN[lp][IMAG];
                        sigma[ik] += (temp1*temp1+temp2*temp2);
                      }

                      if(J > 0 && K>0){
                        lp=ijk(I,this->params._Nft()-J,K,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                        temp1 = P[ik]*data_out_Q[lp][REAL]+data_out_SN[lp][REAL];
                        temp2 = P[ik]*data_out_Q[lp][IMAG]+data_out_SN[lp][IMAG];
                        sigma[ik] += (temp1*temp1+temp2*temp2);
                      }

                      if(I > 0 && (J>0 || K>0)){
                        lp=ijk(this->params._Nft()-I,J,K,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                        temp1 = P[ik]*data_out_Q[lp][REAL]+data_out_SN[lp][REAL];
                        temp2 = P[ik]*data_out_Q[lp][IMAG]+data_out_SN[lp][IMAG];
                        sigma[ik] += (temp1*temp1+temp2*temp2);
                      }

                      if(I > 0 && J > 0 && K>0){
                        lp=index_3d(this->params._Nft()-I,this->params._Nft()-J,K,this->params._Nft(),this->params._Nft()/2+1);
                        temp1 = P[ik]*data_out_Q[lp][REAL]+data_out_SN[lp][REAL];
                        temp2 = P[ik]*data_out_Q[lp][IMAG]+data_out_SN[lp][IMAG];
                        sigma[ik] += (temp1*temp1+temp2*temp2);
                      }
                    }
                  }
                }
              }
            }
        }
    }


#ifdef DOUBLE_PREC
  fftw_free(data_out_SN);
  fftw_free(data_out_Q);
#else
  fftwf_free(data_out_SN);
  fftwf_free(data_out_Q);
#endif


  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::do_fkp_error_bars_veff(s_data_structure_direct_sum *s_data, vector<real_prec> &P,vector<real_prec> &Veff)
{
  this->So.enter(__PRETTY_FUNCTION__);
  So.message_screen("Computing effetive volume");
  ULONG n_objects = s_data->properties.size() ;
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#pragma omp parallel
#endif
  {
    vector<real_prec>table_private(Veff.size(),0);
#ifdef _USE_OMP_
#pragma omp for collapse(2)
#endif
    for(int i = 0; i < n_objects; ++i)
      for(int k=0;k<Veff.size(); ++k)
        {
          real_prec nbar;
          nbar=s_data->properties[i].mean_density;
          table_private[k]+=nbar*1.0/((1+nbar*P[k])*(1+nbar*P[k]));  /*eff vol*/
       }
#ifdef _USE_OMP_
#pragma omp critical
#endif
    for(int k=0;k<Veff.size();++k)
      Veff[k]+=table_private[k];
 }
  return;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_power_moments_fourier_grid_ds_yam(s_data_structure_direct_sum *s_data){
  this->So.enter(__PRETTY_FUNCTION__);
  So.message_screen("Interpolating galaxy density field on a grid");
  real_prec fac      = M_PI/180.0;
  real_prec area     = this->params._area_survey()*pow(fac,2);        /*survey area in rad*/
  real_prec mean_den = s_data->mean_density;
  real_prec DELTAK   = this->params._d_DeltaK_data();
  // Identify columns in the corresponding catalog                               *
  int i_weight1      = (s_data->catalog=="data"? this->params._i_weight1_g() : this->params._i_weight1_r());
  int i_weight2      = (s_data->catalog=="data"? this->params._i_weight2_g() : this->params._i_weight2_r());
  int i_weight3      = (s_data->catalog=="data"? this->params._i_weight3_g() : this->params._i_weight3_r());
  int i_weight4      = (s_data->catalog=="data"? this->params._i_weight4_g() : this->params._i_weight4_r());
  bool use_weight1   = (s_data->catalog=="data"? this->params._use_weight1_g() : this->params._use_weight1_r());
  bool use_weight2   = (s_data->catalog=="data"? this->params._use_weight2_g() : this->params._use_weight2_r());
  bool use_weight3   = (s_data->catalog=="data"? this->params._use_weight3_g() : this->params._use_weight3_r());
  bool use_weight4   = (s_data->catalog=="data"? this->params._use_weight4_g() : this->params._use_weight4_r());
  // Arrays used in the interpolation of mean number density                     *
  vector<gsl_real> zzv =  (s_data->zz_v);
  vector<gsl_real> dndz = (s_data->dndz_v);
  int i_mean_density= (s_data->catalog=="data"? this->params._i_mean_density_g() : this->params._i_mean_density_r());
  string angles_units = (s_data->catalog=="data" ? this->params._angles_units_g() : this->params._angles_units_r());
  ULONG nlines=s_data->properties.size();
  int n_columns= s_data-> n_columns;
  long full=((real_prec)this->params._Nft()/2)*((real_prec)this->params._Nft()/2)*((real_prec)this->params._Nft()/2+1);
  long counter=0;
  real_prec total_weight=1.0;
  ULONG n_selected=0;
  real_prec normal_power=0;
  real_prec W_r=0;
  real_prec S_r_0=0;
  real_prec a11=0;
  real_prec a12=0;
  real_prec a13=0;
  real_prec a14=0;
  real_prec a15=0;
  real_prec a16=0;
  real_prec a21=0;
  real_prec a22=0;
  real_prec a23=0;
  real_prec a24=0;
  real_prec a25=0;
  real_prec a26=0;
  real_prec a31=0;
  real_prec a32=0;
  real_prec a33=0;
  real_prec a34=0;
  real_prec a35=0;
  real_prec a36=0;
  real_prec a41=0;
  real_prec a42=0;
  real_prec a43=0;
  real_prec a44=0;
  real_prec a45=0;
  real_prec a46=0;
  omp_set_num_threads(omp_get_max_threads());
  So.message_screen("Using", omp_get_max_threads(), "threads");
#pragma omp for
  // Loop over the grid
  for(int i=0;i<this->params._Nft()/2;i++){
    int i2=i+this->params._Nft()/2;
    int q_i2= (i==0? this->params._Nft()/2: i-this->params._Nft()/2);
    real_prec kx1=((real_prec)i)*params._d_deltak_x();
    real_prec kx2=q_i2*params._d_deltak_x();
    real_prec kx3=kx2;
    real_prec kx4=kx1;
    for(int j=0;j<this->params._Nft()/2;j++){
      int j2=j+this->params._Nft()/2;
      int q_j2=(j==0? this->params._Nft()/2 : j-this->params._Nft()/2);
      real_prec ky1=((real_prec)j)*params._d_deltak_y();
      real_prec ky2=ky1;
      real_prec ky3=q_j2*params._d_deltak_y();
      real_prec ky4=ky3;

      for(int k=0;k<=this->params._Nft()/2;k++){
        // counter++;
        // So.comp_time(start, full, counter);
        real_prec kz  = ((real_prec)k)*params._d_deltak_z();
        real_prec kk1=sqrt(kx1*kx1+ky1*ky1+kz*kz);   //Determine magnitude of k-vector
        real_prec kk2=sqrt(kx2*kx2+ky1*ky1+kz*kz);   //Determine magnitude of k-vector
        real_prec kk3=sqrt(kx2*kx2+ky3*ky3+kz*kz);   //Determine magnitude of k-vector
        real_prec kk4=sqrt(kx1*kx1+ky3*ky3+kz*kz);   //Determine magnitude of k-vector
        // RESETeando cada vex que inicia un loop sobre galaxias
        total_weight=1.0;
        n_selected=0;
        normal_power=0;
        W_r=0;
        S_r_0=0;
        a11=0;
        a12=0;
        a13=0;
        a14=0;
        a15=0;
        a16=0;
        a21=0;
        a22=0;
        a23=0;
        a24=0;
        a25=0;
        a26=0;
        a31=0;
        a32=0;
        a33=0;
        a34=0;
        a35=0;
        a36=0;
        a41=0;
        a42=0;
        a43=0;
        a44=0;
        a45=0;
        a46=0;
        //omp_set_num_threads(omp_get_max_threads());
        // std::cout<<RED<<"Warning: using one thread in Pk_y_ds"<<RESET<<std::endl;
        // omp_set_num_threads(1);
        real_prec *pprop;
        real_prec ow[4];
        for(ULONG Ig=0;Ig<nlines;++Ig){
          real_prec we;
          real_prec x=s_data->properties[i].coord1;
          real_prec y=s_data->properties[i].coord2;
          real_prec z=s_data->properties[i].coord3;
          real_prec nbar=s_data->properties[i].mean_density;
          real_prec rr= sqrt(x*x+y*y+z*z);
          // ***************************************************************
          // Compute weights for each galaxy
          // Compute the weights for each galaxy
          ow[0] = (use_weight1 && (i_weight1<n_columns))? s_data->properties[i].weight1 : 1.0;
          ow[1] = (use_weight2 && (i_weight2<n_columns))? s_data->properties[i].weight2 : 1.0;
          ow[2] = (use_weight3 && (i_weight3<n_columns))? s_data->properties[i].weight3 : 1.0;
          ow[3] = (use_weight4 && (i_weight4<n_columns))? s_data->properties[i].weight4 : 1.0;
          total_weight = ow[0] * ow[1] * ow[2] * ow[3];
          if(true==this->params._FKP_weight()){we=1.0/(1+this->params._Pest()*nbar);} else we=1.0;
          total_weight*=we;
          real_prec ptotal_weight2=total_weight*total_weight;
          if(kk1!=0){
            int lp1=ijk(i ,j ,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);//Find index as in FFTW
            int lp2=ijk(i2,j ,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);//Find index as in FFTW
            int lp3=ijk(i2,j2,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);//Find index as in FFTW
            int lp4=ijk(i ,j2,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);//Find index as in FFTW
            // *******************************************************************
            // Now compute r, k and mu and assign to the grid
            // first octant:
            if(kk1<=this->params._kmax_y_ds())
              {
                real_prec mu    = (x*kx1+y*ky1+z*kz)/(kk1*rr);   //cos of the angle between r and k
                real_prec cc    = total_weight*cos(kk1*rr*mu);
                real_prec ss    =-total_weight*sin(kk1*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
                a11+=cc;
                a12+=ss;
                a13+=cc*Leg_2;
                a14+=ss*Leg_2;
                a15+=cc*Leg_4;
                a16+=ss*Leg_4;
                if(s_data->catalog=="data"){
                  this->data_g_out_y0[lp1].real(a11);
                  this->data_g_out_y0[lp1].imag(a12);
                  this->data_g_out_y2[lp1].real(a13);
                  this->data_g_out_y2[lp1].imag(a14);
                  this->SN_g_out_y2[lp1]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_g_out_y4[lp1].real(a15);
                  this->data_g_out_y4[lp1].imag(a16);
                  this->SN_g_out_y4[lp1]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
                else{
                  this->data_r_out_y0[lp1].real(a11);
                  this->data_r_out_y0[lp1].imag(a12);
                  this->data_r_out_y2[lp1].real(a13);
                  this->data_r_out_y2[lp1].imag(a14);
                  this->SN_r_out_y2[lp1]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_r_out_y4[lp1].real(a15);
                  this->data_r_out_y4[lp1].imag(a16);
                  this->SN_r_out_y4[lp1]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
              }
            // // second octant:
            if(kk2<=this->params._kmax_y_ds())
              {
                real_prec mu    = (x*kx2+y*ky1+z*kz)/(kk2*rr);
                real_prec cc    = total_weight*cos(kk2*rr*mu);
                real_prec ss    =-total_weight*sin(kk2*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
                a21+=cc;
                a22+=ss;
                a23+=cc*Leg_2;
                a24+=ss*Leg_2;
                a25+=cc*Leg_4;
                a26+=ss*Leg_4;
                if(s_data->catalog=="data"){
                  this->data_g_out_y0[lp2].real(a21);
                  this->data_g_out_y0[lp2].imag(a22);
                  this->data_g_out_y2[lp2].real(a23);
                  this->data_g_out_y2[lp2].imag(a24);
                  this->SN_g_out_y2[lp2]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_g_out_y4[lp2].real(a25);
                  this->data_g_out_y4[lp2].imag(a26);
                  this->SN_g_out_y4[lp2]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
                else{
                  this->data_r_out_y0[lp2].real(a21);
                  this->data_r_out_y0[lp2].imag(a22);
                  this->data_r_out_y2[lp2].real(a23);
                  this->data_r_out_y2[lp2].imag(a24);
                  this->SN_r_out_y2[lp2]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_r_out_y4[lp2].real(a25);
                  this->data_r_out_y4[lp2].imag(a26);
                  this->SN_r_out_y4[lp2]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }

              }
            // Third octant
            if(kk3<=this->params._kmax_y_ds())
              {
                real_prec mu    = (x*kx2+y*ky3+z*kz)/(kk3*rr);
                real_prec cc    = total_weight*cos(kk3*rr*mu);
                real_prec ss    =-total_weight*sin(kk3*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
                a31+=cc;
                a32+=ss;
                a33+=cc*Leg_2;
                a34+=ss*Leg_2;
                a35+=cc*Leg_4;
                a36+=ss*Leg_4;
                if(s_data->catalog=="data"){
                  this->data_g_out_y0[lp3].real(a31);
                  this->data_g_out_y0[lp3].imag(a32);
                  this->data_g_out_y2[lp3].real(a33);
                  this->data_g_out_y2[lp3].imag(a34);
                  this->SN_g_out_y2[lp3]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_g_out_y4[lp3].real(a35);
                  this->data_g_out_y4[lp3].imag(a36);
                  this->SN_g_out_y4[lp3]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
                else{
                  this->data_r_out_y0[lp3].real(a31);
                  this->data_r_out_y0[lp3].imag(a32);
                  this->data_r_out_y2[lp3].real(a33);
                  this->data_r_out_y2[lp3].imag(a34);
                  this->SN_r_out_y2[lp3]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_r_out_y4[lp3].real(a35);
                  this->data_r_out_y4[lp3].imag(a36);
                  this->SN_r_out_y4[lp3]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
              }
            if(kk4<=this->params._kmax_y_ds())
              {
                // fourth octant:
                real_prec mu    = (x*kx1+y*ky3+z*kz)/(kk4*rr);   //cos of the angle between r and k
                real_prec cc    = total_weight*cos(kk4*rr*mu);
                real_prec ss    =-total_weight*sin(kk4*rr*mu);
                real_prec Leg_2 = 0.5*(3.*mu*mu-1.);
                real_prec Leg_4 = (1./8.)*(35*pow(mu,4)-30.*mu*mu+3);
                a41+=cc;
                a42+=ss;
                a43+=cc*Leg_2;
                a44+=ss*Leg_2;
                a45+=cc*Leg_4;
                a46+=ss*Leg_4;
                if(s_data->catalog=="data"){
                  this->data_g_out_y0[lp4].real(a41);
                  this->data_g_out_y0[lp4].imag(a42);
                  this->data_g_out_y2[lp4].real(a43);
                  this->data_g_out_y2[lp4].imag(a44);
                  this->SN_g_out_y2[lp4]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_g_out_y4[lp4].real(a45);
                  this->data_g_out_y4[lp4].imag(a46);
                  this->SN_g_out_y4[lp4]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
                else{
                  this->data_r_out_y0[lp4].real(a41);
                  this->data_r_out_y0[lp4].imag(a42);
                  this->data_r_out_y2[lp4].real(a43);
                  this->data_r_out_y2[lp4].imag(a44);
                  this->SN_r_out_y2[lp4]        += pow(total_weight,2)*Leg_2;     //Shot-noise for multipole
                  this->data_r_out_y4[lp4].real(a45);
                  this->data_r_out_y4[lp4].imag(a46);
                  this->SN_r_out_y4[lp4]        += pow(total_weight,2)*Leg_4;     //Shot-noise for multipoles
                }
              }
            // ***********************************************************
            // Parameters for the Pk
            n_selected++;
            W_r          +=total_weight;    //weighted number of selected objects
            S_r_0        +=ptotal_weight2;                  //Shot-noise for multipoles
            normal_power +=nbar*ptotal_weight2;            //normalization
            // **********************************************************************
          }
        }
      }
    }
  }
  if(s_data->catalog=="data"){
    this->n_gal=n_selected;
    this->w_g=W_r;
    this->s_g=S_r_0;
  }
  else{
    if(s_data->catalog=="random"){
      this->n_ran=n_selected;
      this->normal_p=normal_power;
      this->s_r=S_r_0;
      this->w_r=W_r;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_yam_1d_ds(vector<real_prec>&powerk0, vector<real_prec>&powerk2, vector<real_prec>&powerk4, vector<int> &mod){
  for(int i=0;i<mod.size();i++)mod[i]=0;
  for(int i=0;i<powerk0.size();i++)powerk0[i]=0;
  for(int i=0;i<powerk2.size();i++)powerk2[i]=0;
  for(int i=0;i<powerk4.size();i++)powerk4[i]=0;
  real_prec DELTAK=params._d_DeltaK_data();
  for(int i=0;i<=this->params._Nft();i++)
    {
      int q_i= (i<=this->params._Nft()/2? i: i-(this->params._Nft()+1)); //set coordinate in units of the fundamental mode
      for(int j=0;j<=this->params._Nft();j++)
        {
          int q_j= (j<=this->params._Nft()/2? j: j-(this->params._Nft()+1));  //set coordinate in units of the fundamental mode
          for(int k=0;k<=this->params._Nft();k++)
            {
              int q_k= (k<=this->params._Nft()/2? k: k-(this->params._Nft()+1)); //set coordinate in units of the fundamental mode
              real_prec kv; int kmod;
              if(this->params._type_of_binning()=="linear")
                {
                  kv=sqrt(pow(q_i*this->params._d_deltak_x(),2)+pow(q_j*this->params._d_deltak_y(),2)+pow(q_k*this->params._d_deltak_z(),2));
                  kmod=static_cast<int>(floor((float)((kv-this->params._d_kmin())/DELTAK)));
                }
              else{
                if(this->params._type_of_binning()=="log"){
                  kv=pow(this->params._d_deltak_x()*q_i/(this->params._d_kmin()),2)+pow(params._d_deltak_y()*q_j/(this->params._d_kmin()),2)+pow(params._d_deltak_z()*q_k/(this->params._d_kmin()),2);
                  if(kv!=0)
                    kmod=static_cast<int>(floor((float)( log10(sqrt(kv-this->params._d_kmin()))/this->params._d_Deltal())))+1;
                  else{kmod=0;}
                }
              }
              if(kv<=this->params._kmax_y_ds())
                {
                  int o_i, o_j, o_k;
                  real_prec factor;
                  remap(this->params._Nft(), this->params._Nft(), this->params._Nft(), i, j, k, &o_i, &o_j, &o_k, &factor);
                  ULONG lp=ijk(o_i,o_j,o_k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                  if(lp!=0)
                    {
                      if(kmod<powerk0.size())
                        {
                          powerk0[kmod]+=this->data_g_y0[lp];
                          powerk2[kmod]+=this->data_g_y2[lp];
                          powerk4[kmod]+=this->data_g_y4[lp];
                          mod[kmod]++ ;
                        }
                    }
                }
            }
        }
    }
  for(int i=0;i<powerk0.size();i++)powerk0[i]=(mod[i]==0 ? 0 : powerk0[i]/mod[i]);
  for(int i=0;i<powerk2.size();i++)powerk2[i]=(mod[i]==0 ? 0 : powerk2[i]/mod[i]);
  for(int i=0;i<powerk4.size();i++)powerk4[i]=(mod[i]==0 ? 0 : powerk4[i]/mod[i]);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::bispectrum_fkp(char d_w, vector<real_prec> &bispect, vector<real_prec> &sn_bispect, vector<real_prec> &p1, vector<int> &mod)
{
  //////////////////////////////////////////////////////////
  // Count of triangles in Fourier space
  // There is something wrong with the paralelization
  // in this function, Therefore it is commented.
  //////////////////////////////////////////////////////////


  real_prec Normal_delta= sqrt(normal_power);
  vector<int> modp(this->params._Nft(),0);
  int nnp= (d_w == 'd' ? this->params._d_Nnp_data() : this->params._d_Nnp_window());

  memset(&mod.front(), 0, sizeof(int) * this->params._Nft() * this->params._Nft() * this->params._Nft());
  memset(&bispect.front(), 0, sizeof(real_prec) * this->params._Nft() * this->params._Nft() * this->params._Nft());

  real_prec DELTAK=(d_w == 'd' ? this->params._d_DeltaK_data() : this->params._d_DeltaK_window());
  int ik, n_kmod, n_kmod2, n_kmod3;

  real_prec rDeltaK_data;
  real_prec rDeltaK_window;
  real_prec rkmin;
  int my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear")
    my_s_box_ave = s_box_linear;
  else if(this->params._type_of_binning()=="log")
    my_s_box_ave = s_box_log;

  if(my_s_box_ave == s_box_linear) {
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
    rDeltaK_window = 1.0 / this->params._d_DeltaK_window();
  }
  else if(my_s_box_ave == s_box_log)
    rkmin = 1.0/this->params._d_kmin();
real_prec full=pow((real_prec)this->params._Nft()+1, 2)*pow((real_prec)this->params._Nft()+1, 2)*pow((real_prec)this->params._Nft()+1, 2)-1;
  long counter=0;
  // #pragma omp parallel
  //   {
  real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
  real_prec xx_MAS, yy_MAS, zz_MAS;
  //  First loop k1
  //#pragma omp for nowait
  for(int i=0;i<=this->params._Nft();++i){
    int q_i= (i<=this->params._Nft()/2? i: i-(this->params._Nft()+1)); //set coordinate in units of the fundamental mode
    if(q_i > 0) {
      xx_MAS = q_i * M_PI / this->params._Nft();
      xx_MAS = sin(xx_MAS)/xx_MAS;
    } else xx_MAS = 1.0;
    i_deltak_x_2 = q_i * q_i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
    if(my_s_box_ave == s_box_log)
      i_deltak_x_2 *= rkmin * rkmin;
    for(int j=0;j<=this->params._Nft();++j){
      int q_j= (j<=this->params._Nft()/2? j: j-(this->params._Nft()+1)); //set coordinate in units of the fundamental mode
      if(q_j > 0) {
        yy_MAS = q_j * M_PI / this->params._Nft();
        yy_MAS = sin(yy_MAS)/yy_MAS;
      } else yy_MAS = 1.0;
      j_deltak_y_2 = q_j * q_j * this->params._d_deltak_y() * this->params._d_deltak_y() ;
      if(my_s_box_ave == s_box_log)
        j_deltak_y_2 *= rkmin * rkmin;
      for(int k=0;k<=this->params._Nft();++k){
        int q_k= (k<=this->params._Nft()/2? k: k-(this->params._Nft()+1)); //set coordinate in units of the fundamental mode
        if(q_k > 0) {
          zz_MAS = q_k * M_PI / this->params._Nft();
          zz_MAS = sin(zz_MAS)/zz_MAS;
        } else zz_MAS = 1.0;
        k_deltak_z_2 = q_k * q_k * this->params._d_deltak_z() * this->params._d_deltak_z() ;
        if(my_s_box_ave == s_box_log)
          k_deltak_z_2 *= rkmin * rkmin;

        real_prec kv, kvv, kvvv;
        int kmod, kmod2, kmod3;
        if(this->params._type_of_binning()=="linear"){
          kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
          kmod=(int)floor((float)((kv-this->params._d_kmin())/DELTAK));
        }
        else{
          if(this->params._type_of_binning()=="log"){
            kv=(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
            kmod=(kv==0? 0 : (int)floor((float)( log10(sqrt(kv))/this->params._d_Deltal()))+1);
          }
        }
        real_prec corr1 = (this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.00);
        int o_i, o_j, o_k;
        real_prec factor1;
        remap(this->params._Nft(), this->params._Nft(), this->params._Nft(), i, j, k, &o_i, &o_j, &o_k, &factor1);
        int lp=ijk(o_i,o_j,o_k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
        real_prec d1r   = this->data_out_g[lp][REAL] / corr1;
        real_prec d1i   = factor1*this->data_out_g[lp][IMAG] / corr1;
        real_prec p1s = pow(d1i,2)+pow(d1r,2);
        if(lp!=0){
          //	  #pragma omp atomic
          p1[kmod] += d1r*d1r + d1i*d1i;  //Compute power spectrum
          //#pragma omp atomic
          modp[kmod]++;
        }
        // *********************************************************************************
        // Second loop k2
        real_prec ii_deltak_x_2, jj_deltak_y_2, kk_deltak_z_2;
        real_prec xxx_MAS, yyy_MAS, zzz_MAS;
        for(int ii=0;ii<=this->params._Nft();++ii){
          int q_ii=(ii<=this->params._Nft()/2? ii: ii-(this->params._Nft()+1));
          if(q_ii > 0) {
            xxx_MAS = q_ii * M_PI / this->params._Nft();
            xxx_MAS = sin(xxx_MAS)/xxx_MAS;
          } else xxx_MAS = 1.0;
          ii_deltak_x_2 = q_ii * q_ii * this->params._d_deltak_x() * this->params._d_deltak_x() ;
          if(my_s_box_ave == s_box_log)
            ii_deltak_x_2 *= rkmin * rkmin;
          for(int jj=0;jj<=this->params._Nft();++jj){
            int q_jj=(jj<=this->params._Nft()/2? jj: jj-(this->params._Nft()+1));
            if(q_jj > 0) {
              yyy_MAS = q_jj * M_PI / this->params._Nft();
              yyy_MAS = sin(yyy_MAS)/yyy_MAS;
            } else yyy_MAS = 1.0;
            jj_deltak_y_2 = q_jj * q_jj * this->params._d_deltak_y() * this->params._d_deltak_y() ;
            if(my_s_box_ave == s_box_log)
              jj_deltak_y_2 *= rkmin * rkmin;
            for(int kk=0;kk<=this->params._Nft();++kk){
              int q_kk=(kk<=this->params._Nft()/2? kk: kk-(this->params._Nft()+1));
              if(q_kk > 0) {
                zzz_MAS = q_kk * M_PI / this->params._Nft();
                zzz_MAS = sin(zzz_MAS)/zzz_MAS;
              } else zzz_MAS = 1.0;
              kk_deltak_z_2 = q_kk * q_kk * this->params._d_deltak_z() * this->params._d_deltak_z() ;
              if(my_s_box_ave == s_box_log)
                kk_deltak_z_2 *= rkmin * rkmin;
              counter ++;
             if(my_s_box_ave==s_box_linear){
                kvv=sqrt(ii_deltak_x_2 + jj_deltak_y_2 + kk_deltak_z_2);
                kmod2=(int)floor((float)((kvv-this->params._d_kmin())/DELTAK));
              }
              else{
                if(my_s_box_ave==s_box_log){
                  kvv=(ii_deltak_x_2 + jj_deltak_y_2 + kk_deltak_z_2);
                  kmod2=(kvv==0? 0 : (int)floor((float)( log10(sqrt(kvv))/this->params._d_Deltal()))+1);
                }
              }
              int o_ii, o_jj, o_kk;
              real_prec factor2;
              remap(this->params._Nft(), this->params._Nft(), this->params._Nft(), ii, jj, kk, &o_ii, &o_jj, &o_kk, &factor2);
              int lp2=ijk(o_ii,o_jj,o_kk,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
              real_prec corr2 = (this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xxx_MAS, yyy_MAS, zzz_MAS): 1.000);
              real_prec d2r   = this->data_out_g[lp2][0]/corr2;
              real_prec d2i   = factor2*this->data_out_g[lp2][1]/corr2;
              real_prec p2s = pow(d2i,2)+pow(d2r,2);
              // set the components of the vector k3 in units of DELTAk from the components
              // of the vectors k1 and k2 to fullfil k3=-k1-k2
              int q_iii=-q_i-q_ii;
              int q_jjj=-q_j-q_jj;
              int q_kkk=-q_k-q_kk;
              // Now proceed if the triangle is in the grid:
              // the components of the third vector are required to be in the box,
              // that is, 0<= |q_i|<=n/2
              if(fabs(q_iii)<=this->params._Nft()/2 && fabs(q_jjj)<=this->params._Nft()/2 && fabs(q_kkk)<=this->params._Nft()/2){
                if(my_s_box_ave== s_box_linear){
                  kvvv=sqrt(q_iii*q_iii*params._d_deltak_x()*params._d_deltak_x() +
                            q_jjj*q_jjj*params._d_deltak_y()*params._d_deltak_y() +
                            q_kkk*q_kkk*params._d_deltak_z()*params._d_deltak_z());
                  //kmod3=(int)floor((float)(kvvv/DELTAK)); //OFFICIAL BINNING, RESTORE IT WHEN COMPARISON SUCCEDS
                  kmod3=(int)floor((float)((kvvv-this->params._d_kmin())/DELTAK)); //JENN AND CRIS BINNING,
                }
                else{
                  if(my_s_box_ave== s_box_log){
                    kvvv = (q_iii*q_iii*params._d_deltak_x()*params._d_deltak_x() +
                            q_jjj*q_jjj*params._d_deltak_y()*params._d_deltak_y() +
                            q_kkk*q_kkk*params._d_deltak_z()*params._d_deltak_z()) * rkmin*rkmin;
                    kmod3=(kvvv==0? 0 : (int)floor((float)( log10(sqrt(kvv))/this->params._d_Deltal()))+1);
                  }
                }
                // Find the position in the grid from the components of the k3 vector
                int iii = (q_iii>=0 ? q_iii: q_iii+this->params._Nft()+1);
                int jjj = (q_jjj>=0 ? q_jjj: q_jjj+this->params._Nft()+1);
                int kkk = (q_kkk>=0 ? q_kkk: q_kkk+this->params._Nft()+1);
                int o_iii, o_jjj, o_kkk;
                real_prec factor3;
                remap(this->params._Nft(), this->params._Nft(), this->params._Nft(), iii, jjj, kkk, &o_iii, &o_jjj, &o_kkk, &factor3);
                int lp3=ijk(o_iii,o_jjj,o_kkk,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
                if(lp!=0 &&  lp2!=0 && lp3!=0){			  // Exclude the zero-mode
                  real_prec corr3=(this->params._MAS_correction() ? correction_MAS(this->params._Nft(), this->params._mass_assignment(),q_iii,q_jjj,q_kkk): 1.000);
                  real_prec d3r = this->data_out_g[lp3][0]/corr3;
                  real_prec d3i = factor3*this->data_out_g[lp3][1]/corr3;
                  real_prec p3s = pow(d3i,2)+pow(d3r,2);
                  sort_index(kmod, kmod2, kmod3, &n_kmod, &n_kmod2, &n_kmod3);	//Sort
                  //#pragma omp atomic
                  ULONG ind=index_3d(n_kmod, n_kmod2, n_kmod3,this->params._Nft(), this->params._Nft());
                  mod[ind]++ ;
                  bispect[ind] +=(d1r*d2r*d3r- d1r*d2i*d3i- d1i*d2r*d3i- d1i*d2i*d3r);
                  sn_bispect[ind] += p1s+p2s+p3s;
                  // ************ Uncomment this to check the imaginary part ********************************//
                  //bispect[n_kmod][n_kmod2][n_kmod3] += (d1r*d2r*d3i + d1r*d2i*d3r + d1i*d2r*d3r - d1i*d2i*d3i);
                }
              }
            }
          }
        }
      }
    }
    //  } // END of outern parallel for
  } // END of parallel region
  //Compute average power in the spherical shells
  for(int i=0;i<this->params._Nft();++i)p1[i]= (modp[i]==0? 0. : p1[i]/modp[i]);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_power_spectrum_for_bispectrum(vector<real_prec> &power_g0)
{
  So.message_screen("Power spectrum for Bispectrum");
  do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
  power_spectrum_fkp_for_bispectrum(power_g0);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_bispectrum_fkp(char d_w, vector<real_prec> &bispect, vector<real_prec> &sn_bispect, vector< int > &mod)
{
  do_fftw_r2c(this->params._Nft(),this->data_g, this->data_out_g);
  vector<real_prec> power(this->params._Nft(),0);
  memset(&mod.front(), 0, sizeof(int) * this->params._Nft() * this->params._Nft() * this->params._Nft());
  memset(&bispect.front(), 0, sizeof(real_prec) * this->params._Nft() * this->params._Nft() * this->params._Nft());
  memset(&sn_bispect.front(), 0, sizeof(real_prec) * this->params._Nft() * this->params._Nft() * this->params._Nft());
  std::cout << CYAN << "Measuring bispectrum using brute-force approach ..." << RESET << endl;
  time_t start, stop;
  time (&start);
  bispectrum_fkp(d_w, bispect,sn_bispect, power, mod);
  time (&stop);
  std::cout << "done in " << difftime(stop,start) << " seconds" << endl;
  real_prec SN_aux;
  if(this->params._SN_correction())
    SN_aux=1.0;
  else
    SN_aux=0;
  //Estimates of the power spectrum using FKP, as done in the Pk main code.
  for(int i=0;i<power.size();++i)power[i]=power[i]/this->normal_power-shot_noise*SN_aux;
#pragma omp parallel
  {
    ULONG I;
#pragma omp for collapse(2) nowait
    for(int i=0;i<this->params._Nft();++i)
      for(int j=0;j<this->params._Nft();++j)
        for(int k=0;k<this->params._Nft();++k)
          if(i<=j && j<=k) {
            I = index_3d(i,j,k, this->params._Nft(), this->params._Nft());
            bispect[I]=(mod[I]== 0 ? 0 : bispect[I]/((real_prec)mod[I]*this->normal_b) - SN_aux*( (sn_bispect[I]/mod[I])*shot_noise_b1+shot_noise_b2));
          }
  }
#ifdef DOUBLE_PREC
  fftw_free(data_out_g);
#else
  fftwf_free(data_out_g);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::define_kshells(){
  // Getting ready for MAS correction
  real_prec xx_MAS, yy_MAS, zz_MAS; //for MAS correction
  real_prec *yy_MAS_array, *zz_MAS_array;
  real_prec *_ky_array, *_ky2_array;
  int q_i, q_j, q_k;
  yy_MAS_array = new real_prec[this->params._Nft()+1];
  zz_MAS_array = new real_prec[this->params._Nft()+1];
  real_prec M_PI_rNft = M_PI * this->rNft;
  yy_MAS_array[0] = 1.0;
  for(int j = 1; j<= this->params._Nft(); j++) {
    q_j= (j<=this->params._Nft()/2? j: j-(this->params._Nft()+1));
    yy_MAS = q_j * M_PI_rNft;
    yy_MAS_array[j] = (yy_MAS = sin(yy_MAS)/yy_MAS);
  }
  zz_MAS_array[0] = 1.0;
  for(int k = 1; k <= this->params._Nft(); k++) {
    q_k= (k<=this->params._Nft()/2? k: k-(this->params._Nft()+1));
    zz_MAS = q_k * M_PI_rNft;
    zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
  }
  // Given input parameter of kmin, kmax and the number of shells,
  // here the k-shells are defined
  int id=0;
  // This is a loop over the Fourier grid
  // that goes to the maximum value of k given as an
  // input, kmax_bk, in units of the fundamental mode.
  // From this value, the sgrid is computed in the set_pars method of this class.
  // Viewed as an quarter, the loop if done over the modes below the diagonal.
  // The other modes can be mapped using some Symmetries, designed by Jennifer.
  for(int i=0;i<this->sgrid;++i){
    if(i > 0) {
      xx_MAS = i * M_PI_rNft;
      xx_MAS = sin(xx_MAS)/xx_MAS;
    } else xx_MAS = 1.0;

    for(int j=i;j<this->sgrid;++j){
      yy_MAS = yy_MAS_array[j];
      for(int k=j;k<this->sgrid;++k){
        zz_MAS = zz_MAS_array[k];
        // Compute the squared of the magnitude, in units of the fundamental mode.
        real_prec kk= i*i+j*j+k*k;
        // Exclude the zero*mode. Be consistent with line 75 (set_pars())
        // in order to get the right dimension of the arrays
        if(kk>0){
         // Assign the MAS correction for this vector
          this->Array_corr[id]=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.000000);
         // Assign to the mode id its x-component in the following array
          this->Arraykx[id]=i;
         // Assign to the mode id its y-component in the following array:
          this->Arrayky[id]=j;
          // Assign to the mode id its z-component in the following array:
          this->Arraykz[id]=k;
         // Assign to the mode id its own original ID
          this->ArrayID[id]=id;
         // Assign to the mode id its magnitude squared in the following array:
          this->Arraykk[id]=kk;
         // Symmetries. These numbers are sacred.
          {
            if(i!=0 && i!=j && j!=k)this->VecArray[id]=48;
            else if(i!=0 && i!=j && j!=0 && j==k)this->VecArray[id]=25;
            else if(i!=0 && i==j && j!=k)this->VecArray[id]=24;
            else if(i==0 && i!=j && j!=k)this->VecArray[id]=23;
            else if(i==0 && j!=0 && j==k)this->VecArray[id]=12;
            else if(i!=0 && i==j && j==k)this->VecArray[id]=8;
            else if(i==0 && j==0 && k!=0)this->VecArray[id]=6;
            else if(i==0 && j==0 && k==0)this->VecArray[id]=0;
          }
         int kbin = (kk==0? 0 : (int)floor(((this->params._d_deltak_x()*sqrt(static_cast<real_prec>(kk))-this->params._d_kmin())/this->DeltaK_Bis)));
         //Assign to a mode identified with the id, the bin where it belongs to
          this->Kbin[id] = kbin;
         // Count the number of modes in the bin.
          if(kbin<this->Nshells_bk)
            this->Bmodes[kbin] ++;
          // Count the number of modes
          id++;
        }
      }
    }
  }
  // Define k bins
  for(int i=0;i<this->Nshells_bk;++i)
    this->kbins_bk[i]=this->params._d_kmin()+(i+0.5)*this->DeltaK_Bis;
  int nn_gg=0;
  for(int i=0;i<this->Nshells_bk;++i)
    {
      nn_gg=3*(int)((this->kbins_bk[i])/this->params._d_DeltaK_data());
      if(nn_gg%2!=0)nn_gg++;
      this->Ngrids_bk[i]=nn_gg;
    }
  // Sort array with the squared of the magnitude of the wavevector
  // and the other vectors, shuffled accordingly:
  sort_vectors(this->Arraykk, this->Arraykx,this->Arrayky,this->Arraykz,this->Kbin, this->VecArray, this->ArrayID, this->Array_corr);
  // lopp over sorted indices according to the k squarde value
  // Find, for each shell, the ID associated to the
  // vector with the smallest value of k**2 in that particular shell
  // This can be done easily because the vectors are already sorted.
  int idd=0;
  ArrayID[0]=0;
  //#pragma omp parallel for
  for(int i=1;i<this->new_sgrid;++i)
    {  //start from 1 to evaluate i-1 without segfault. That's why we lieave the id++ here instead of writing it below
      idd++;
      this->ArrayID[i]=i; //resort, ask Jennifer
      if(Kbin[i]!=Kbin[i-1]){
        int bin = Kbin[i];
        if(bin<this->Nshells_bk){
          this->kkminID[bin]=idd;
        }
      }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_ift_shells_bispectrum()
{
  So.message_screen("Construct Inverse Fourier transforms in shells");
  int Nshell_min = (floor)( ((this->params._kmin_bk()-this->params._d_kmin())/this->DeltaK_Bis));
  this->iFT_output_delta_full.resize(this->Nshells_bk);
  this->iFT_output_triangles_full.resize(this->Nshells_bk);
  this->iFT_shot_noise_p1_cyc_sum_full.resize(this->Nshells_bk);
  int Ngrid_max=Ngrids_bk[this->Nshells_bk-1];
  ULONG Nt_max=(Ngrid_max*Ngrid_max)*Ngrid_max;
  vector<real_prec> iFT_output_delta(Nt_max,0);
  vector<real_prec> iFT_output_triangles(Nt_max,0);
  vector<real_prec> iFT_output_p1_cyc_sum(Nt_max,0);
  // LOOP OVER SHELLS:
  for(int i=0;i<this->Nshells_bk;++i)
  {
    // Remap arrays and generate the iFT
    construct_shells(Ngrid_max, this->kkminID[i],  this->kkminID[i] +this->Bmodes[i]-1, iFT_output_delta, iFT_output_triangles, iFT_output_p1_cyc_sum);
    // For each shell allocate smae amount of memory
    // to the maximum, though loops will go as far
    // as the inverse FT needs.
    // For deltas
    this->iFT_output_delta_full[i].resize(Nt_max,0);
    // For number of triangles
    this->iFT_output_triangles_full[i].resize(Nt_max,0);
    // For the shot-noise
    this->iFT_shot_noise_p1_cyc_sum_full[i].resize(Nt_max,0);
    // Reassign arrays to matrices
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_output_delta_full[i][ik] = iFT_output_delta[ik]/Ngrid_max;
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_output_triangles_full[i][ik]= iFT_output_triangles[ik]/Ngrid_max;
    for(int ik=0;ik<Nt_max;++ik)
      this->iFT_shot_noise_p1_cyc_sum_full[i][ik]= iFT_output_p1_cyc_sum[ik]/Ngrid_max;
  }
  iFT_output_delta.shrink_to_fit();
  iFT_output_triangles.shrink_to_fit();
  iFT_output_p1_cyc_sum.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::loop_shells_bispectrum(vector<real_prec> &pk, vector<real_prec> &bispect, vector< int > &mod, string fname){
  // Next step: loop over shells
  std::cout<<CYAN<<"Loop over shells:"<<RESET<<endl;
  ofstream out;
  out.open(fname.c_str() , ios::out);
  out.precision(10);
  out.setf(ios::showpoint);
  out.setf(ios::scientific);
  cout.precision(12);
  cout.setf(ios::scientific);
  time_t start;
  time (&start);
  int Nshell_min = (floor)( ((this->params._kmin_bk()-this->params._d_kmin())/this->DeltaK_Bis));
  int Ngrid_max=Ngrids_bk[this->Nshells_bk-1];
  long Nt_max=(Ngrid_max*Ngrid_max)*Ngrid_max;
  real_prec SN_aux;
  if(this->params._SN_correction())SN_aux=1.0;
  else SN_aux=0;
  for(int i=Nshell_min;i<this->Nshells_bk;++i)
    {
      int nstep2 = floor(i/2);  //start at this position for the second loop
      int jstart=(nstep2==0 ? i: nstep2);
      for(int j=jstart;j<=i;++j)
        {
          // Prepare for loop over k:
          int nstep3=i-j-1;
          int kstart=(nstep3>j ? j: nstep3);
          if(nstep3<=0)
            kstart=Nshell_min;
          for(int k=kstart;k<=j;++k)
            {
              // Classify triangle bin as type [1] equilateral, [2]=isoceles, [3]Scalene
              int Triangle_type;
              if(i==j && j==k)
                Triangle_type=1;
              else if(i==j && j!=k)
                Triangle_type=2;
              else if(i!=j && j==k)
                Triangle_type=2;
              else if(i!=j && j!=k && i!=k)
                Triangle_type=3;
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              real_prec ntri=0;
              real_prec c_bispectrum=0;
              real_prec c_shot_noise=0;
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              switch(Triangle_type)
                {
                case(1): // Equilateral
                  for(int ik=0;ik<Nt_max;++ik)
                    ntri+=(this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];
                  if(ntri>0)
                    for(int ik=0;ik<Nt_max;++ik)
                      c_bispectrum+=(this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[i][ik])*this->iFT_output_delta_full[i][ik];
                  c_bispectrum/=(this->normal_b*(real_prec)ntri);

                  for(int ik=0;ik<Nt_max;++ik)
                    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];
                  c_shot_noise*=3.0/((real_prec)ntri);
                  break;
                case(2): // Isoceles
                  if(j==k)
                    {
                      for(int ik=0;ik<Nt_max;++ik)
                        ntri+=this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[j][ik]*this->iFT_output_triangles_full[j][ik];
                      if(ntri>0)
                        for(int ik=0;ik<Nt_max;++ik)
                          c_bispectrum+=this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[j][ik]*this->iFT_output_delta_full[j][ik];
                      for(int ik=0;ik<Nt_max;++ik)
                        c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];

                      for(int ik=0;ik<Nt_max;++ik)
                        c_shot_noise+=2.0*(this->iFT_shot_noise_p1_cyc_sum_full[j][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];

                    }
                  else if(k!=j)
                    { //i==j
                      for(int ik=0;ik<Nt_max;++ik)
                        ntri+=pow(this->iFT_output_triangles_full[i][ik],2)*this->iFT_output_triangles_full[k][ik];
                      if(ntri>0)
                        for(int ik=0;ik<Nt_max;++ik)
                          c_bispectrum+=pow(this->iFT_output_delta_full[i][ik],2)*this->iFT_output_delta_full[k][ik];

                      for(int ik=0;ik<Nt_max;++ik)
                        c_shot_noise+=2.0*(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];

                      for(int ik=0;ik<Nt_max;++ik)
                        c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[k][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[i][ik];

                    }
                  c_bispectrum/=(this->normal_b*static_cast<real_prec>(ntri));
                  c_shot_noise/=(real_prec)ntri;
                  ntri*=3; // Just to compare using the right number of modes in the full Fourier space!
                  break;
                  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                case(3): // Scaleno
                  for(int ik=0;ik<Nt_max;++ik)
                    ntri+=(this->iFT_output_triangles_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];
                  if(ntri>0)
                    for(int ik=0;ik<Nt_max;++ik)
                      c_bispectrum+=(this->iFT_output_delta_full[i][ik]*this->iFT_output_delta_full[j][ik])*this->iFT_output_delta_full[k][ik];
                  c_bispectrum/=(this->normal_b*(real_prec)ntri);


                  for(int ik=0;ik<Nt_max;++ik)
                    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[i][ik]*this->iFT_output_triangles_full[j][ik])*this->iFT_output_triangles_full[k][ik];

                  for(int ik=0;ik<Nt_max;++ik)
                    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[j][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[k][ik];

                  for(int ik=0;ik<Nt_max;++ik)
                    c_shot_noise+=(this->iFT_shot_noise_p1_cyc_sum_full[k][ik]*this->iFT_output_triangles_full[i][ik])*this->iFT_output_triangles_full[j][ik];


                  c_shot_noise/=(real_prec)ntri;

                  ntri*=6;// Just to compare usingthe right number of modes in the full Fourier space!

                  break;


                }
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

              //for(int ik=0;ik<Nt_n;++ik)cout<<ntri<<"\t"<<ik<<"\t"<<i<<"\t"<<this->iFT_output_triangles_full[i][ik]<<"\t"<<j<<" \t"<<this->iFT_output_triangles_full[j][ik]<<"\t"<<k<<"\t"<<this->iFT_output_triangles_full[k][ik]<<endl;
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              // Normalization of the Bispectrum: divide by a factor (Ngrid**3)*(Number of triangles)*(Normalization of the estimator)
              // Shot-noise correction
              // 	c_bispectrum-=(SN_aux*((pk[i]+pk[j]+pk[k])*this->shot_noise_b1+this->shot_noise_b2));
              c_bispectrum-= SN_aux*(c_shot_noise*this->shot_noise_b1+this->shot_noise_b2);
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              ULONG new_ind=index_3d(i,j,k,this->Nshells_bk,this->Nshells_bk);
              bispect[new_ind]=c_bispectrum;
              mod[new_ind]=ntri;
              //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              out<<i<<"\t"<<j<<"\t"<<k<<"\t"<<c_bispectrum<<"\t"<<c_shot_noise<<"\t"<<ntri<<endl;
              //cout<<i<<"\t"<<j<<"\t"<<k<<"\t"<<c_bispectrum<<"\t"<<c_shot_noise<<"\t"<<ntri<<endl;
            }
        }
    }
  out.close();
  std::cout<<CYAN<<"Output file "<<fname<<RESET<<std::endl;
  std::cout<<RED;
  time_t end;
  time (&end);
  real_prec diff = difftime(end,start);
  if (diff<60) std::cout <<"Lapse: "<<diff<<" secs"<<endl;
  else if (diff<3600) std::cout <<"Lapse: "<<diff/60<<" minutes"<<endl;
  else std::cout<<"Lapse: "<<diff/3600<<" hours"<<endl;
  std::cout<<RESET;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::construct_shells(int ngrid, int kmnid, int kmxid, vector<real_prec>&iFT_output_delta, vector<real_prec>&iFT_output_triangles, vector<real_prec>&iFT_output_p1_cyc_sum){
  long Ntt=ngrid*ngrid*(ngrid/2+1);
  // Inputs of the ift
  // Array used to get the count of triangles
#ifdef DOUBLE_PREC
  complex_prec *data_in_ks=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
#else
  complex_prec *data_in_ks=(complex_prec *)fftwf_malloc(2*Ntt*sizeof(real_prec));
#endif
  for(int i=0;i<Ntt;++i)data_in_ks[i][0]=0;
  for(int i=0;i<Ntt;++i)data_in_ks[i][1]=0;

  // Inputs of the ift
  // Array used to get the values of delta
#ifdef DOUBLE_PREC
  complex_prec *data_in_dk=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
#else
  complex_prec *data_in_dk=(complex_prec *)fftwf_malloc(2*Ntt*sizeof(real_prec));
#endif
  for(int i=0;i<Ntt;++i)data_in_dk[i][0]=0;
  for(int i=0;i<Ntt;++i)data_in_dk[i][1]=0;
  // Inputs of the ift
  // Array used to get shot noise power
#ifdef DOUBLE_PREC
  complex_prec *data_pk_shot_noise=(complex_prec *)fftw_malloc(2*Ntt*sizeof(real_prec));
#else
  complex_prec *data_pk_shot_noise=(complex_prec *)fftwf_malloc(2*Ntt*sizeof(real_prec));
#endif
  for(int i=0;i<Ntt;++i)data_pk_shot_noise[i][0]=0;
  for(int i=0;i<Ntt;++i)data_pk_shot_noise[i][1]=0;
 // Create 3d shell mesh to identify all vector components
  // This is a loop over the modes contained in this shell, from the
  // vector with the lowest value of k*k to the vector with the highest value.
  //  for(int i=kmnid;i<=kmxid;++i)cellsym(this->ArrayID[i],ngrid,data_in_ks, data_in_dk);
  for(int i=kmnid;i<=kmxid;++i)
    cellsym(i,ngrid,data_in_ks, data_in_dk, data_pk_shot_noise);
  // Do Inverse Fourier transforms:
  do_fftw_c2r(this->params._Nft(),data_in_ks, iFT_output_triangles);
  // // Do IFT of data_in_dk
  do_fftw_c2r(this->params._Nft(),data_in_dk, iFT_output_delta);
  // // Do IFT of data_pk_sn
  do_fftw_c2r(this->params._Nft(),data_pk_shot_noise, iFT_output_p1_cyc_sum);
#ifdef DOUBLE_PREC
  fftw_free(data_in_ks);
  fftw_free(data_in_dk);
  fftw_free(data_pk_shot_noise);
#else
  fftwf_free(data_in_ks);
  fftwf_free(data_in_dk);
  fftwf_free(data_pk_shot_noise);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::power_spectrum_fkp_for_bispectrum(vector<real_prec> & powerk_g0)
{
  // Shell average in Fourier space.
  // for the bispectrum
  int NTHREADS = 1;
#ifdef _USE_OMP_
  NTHREADS = _NTHREADS_;
#endif
  vector<int>mod_g(powerk_g0.size(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<mod_g.size();++i)mod_g[i]=0;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<mod_g.size();++i)powerk_g0[i]=0;
  real_prec rDeltaK_data;
  real_prec rkmin;
  int my_s_box_ave;
#define s_box_linear 1
#define s_box_log 2
#define s_box_other 0
  my_s_box_ave = s_box_other;
  if(this->params._type_of_binning()=="linear")  my_s_box_ave = s_box_linear;
  else if(this->params._type_of_binning()=="log")my_s_box_ave = s_box_log;

  if(my_s_box_ave == s_box_linear) {
    rDeltaK_data = 1.0 / this->params._d_DeltaK_data();
  }
  else if(my_s_box_ave == s_box_log)
    rkmin = 1.0/this->params._d_kmin();

  real_prec *yy_MAS_array, *zz_MAS_array;
  yy_MAS_array = new real_prec[this->params._Nft()/2];
  zz_MAS_array = new real_prec[this->params._Nft()/2+1];

  real_prec M_PI_rNft = M_PI * this->rNft;

  time_t start;
  time (&start);

#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
#endif
  {
    real_prec i_per_fact, j_per_fact;
    real_prec i_deltak_x_2, j_deltak_y_2, k_deltak_z_2;
    real_prec xx_MAS, yy_MAS, zz_MAS;
    real_prec rxx_MAS, ryy_MAS, rzz_MAS;

#ifdef _USE_OMP_
#pragma omp sections
#endif
    {
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        yy_MAS_array[0] = 1.0;
        for(int j=1;j<this->params._Nft()/2;++j) {
          yy_MAS = j * M_PI_rNft;
          yy_MAS_array[j] = sin(yy_MAS)/yy_MAS;
        }
      }
#ifdef _USE_OMP_
#pragma omp section
#endif
      {
        zz_MAS_array[0] = 1.0;
        for(int k=1;k<=this->params._Nft()/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = k * M_PI_rNft;
          zz_MAS_array[k] = (zz_MAS = sin(zz_MAS)/zz_MAS);
        }
      }
    }

#pragma omp barrier

#pragma omp for nowait  //reduction(+:powerk_g0, powerk_g2, powerk_g4, mod_g, power2d_spher, mod_spher, power2d_cart, mod_cart, powerk_r, mod_r)
    for(int i=0;i<this->params._Nft()/2;++i){  // loop over octant kx>0. Half of what FFT gives
      if(i > 0) {
        xx_MAS = i * M_PI_rNft;
        xx_MAS = sin(xx_MAS)/xx_MAS;
      } else xx_MAS = 1.0;

      i_deltak_x_2 = i * i * this->params._d_deltak_x() * this->params._d_deltak_x() ;
      i_per_fact = i_deltak_x_2;
      if(my_s_box_ave == s_box_log)
        i_deltak_x_2 *= rkmin * rkmin;

      for(int j=0;j<this->params._Nft()/2;++j){  // loop over octant ky>0. Half of what FFT gives
        yy_MAS = yy_MAS_array[j];

        j_deltak_y_2 = j * j * this->params._d_deltak_y() * this->params._d_deltak_y();
        j_per_fact = i_per_fact + j_deltak_y_2;
        j_per_fact = sqrt(j_per_fact);
        if(my_s_box_ave == s_box_log)
          j_deltak_y_2 *= rkmin * rkmin;

        for(int k=0;k<=this->params._Nft()/2;++k){  // loop over octant kz>0: this is all FFTW gives.
          zz_MAS = zz_MAS_array[k];

          k_deltak_z_2 = k * k * this->params._d_deltak_z() * this->params._d_deltak_z();
          if(my_s_box_ave == s_box_log)
            k_deltak_z_2 *= rkmin * rkmin;

          real_prec kv;
          int kmod_g;

          // Compute k-shell index
          if(my_s_box_ave == s_box_linear){
            kv=sqrt(i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2);
            kmod_g=(int)floor((float)((kv-this->params._d_kmin())*rDeltaK_data));
          }
          else{
            if(my_s_box_ave == s_box_log) {
              kv= i_deltak_x_2 + j_deltak_y_2 + k_deltak_z_2;
              if(kv!=0){
                kmod_g=(int)floor((float)( log10(sqrt(kv))/this->params._d_Deltal()))+1;
              }
              else{kmod_g=0;}
            }
          }
         // Compute correction for MAS:
          real_prec corr=(this->params._MAS_correction() ? correction_MAS(this->params._mass_assignment(), xx_MAS, yy_MAS, zz_MAS): 1.0);
          real_prec icorr2 = 1.0/(corr*corr);

          // Compute index in c-order
          int lp=ijk(i,j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);


          // Compute thisgs related to 2d power spectrum and multipole
          // decomposition for the first octant
          // Define the angle bewtween k and los:
          // In the FKP we need to specify the LOS direction. set to z by default.


          //Exclude the zero-frequency onLY when computing the P: for the W keep it:
          if(lp!=0){
            real_prec Pk= (this->data_out_g[lp][REAL] * this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG] * this->data_out_g[lp][IMAG]) * icorr2;

            // 1d Spherical average
            if(kmod_g<powerk_g0.size()){
#pragma omp atomic
              powerk_g0[kmod_g]  += Pk;
#pragma omp atomic
              mod_g[kmod_g]++ ;
            }

            //  *****************************************
            //  Now, add modes in other octants:
            //  Add negative frequencies in y:
            //  *****************************************
            if(j>0  && k>0){
              lp=ijk(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);   // Go to the ky<0 octant
              real_prec Pk= (this->data_out_g[lp][REAL]*this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG]*this->data_out_g[lp][IMAG]) * icorr2;

              // 1d Spherical average
              if(kmod_g<powerk_g0.size()){
#pragma omp atomic
                powerk_g0[kmod_g]  += Pk;
#pragma omp atomic
                mod_g[kmod_g]++ ;
              }
            }

            // ******************************************
            // Add negative frequencies in x:
            // ******************************************
            if(i>0  && (j>0 || k>0)){
              lp=ijk(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
              real_prec Pk= (this->data_out_g[lp][REAL]*this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG]*this->data_out_g[lp][IMAG]) * icorr2;
              if(kmod_g<powerk_g0.size()){
#pragma omp atomic
                powerk_g0[kmod_g]  += Pk;
#pragma omp atomic
                mod_g[kmod_g]++;
              }
            }

            // ******************************************
            //  Add negative frequencies in x and y:
            // ******************************************
            if(i>0  && j>0  && k>0){
              lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
              real_prec Pk= (this->data_out_g[lp][REAL]*this->data_out_g[lp][REAL] + this->data_out_g[lp][IMAG]*this->data_out_g[lp][IMAG]) * icorr2;
              if(kmod_g<powerk_g0.size()){
#pragma omp atomic
                powerk_g0[kmod_g] += Pk;
#pragma omp atomic
                mod_g[kmod_g] ++ ;
              }
            }
          }
        }
      }
    }
  }
 // Subtract shot noise and normalize only he 1d power spectrum.
  real_prec SN_aux=0;
  if(this->params._SN_correction())
    SN_aux=this->shot_noise;

  // For the monopole:
  for(int i=0;i<powerk_g0.size();++i)powerk_g0[i]=(mod_g[i]== 0 ? 0 : powerk_g0[i]/((real_prec)mod_g[i]*this->normal_power)-SN_aux);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::cellsym(int id, int ngrid, complex_prec *data_ks, complex_prec *data_dk,complex_prec *data_pk_shot_noise){
  int sx=this->Arraykx[id];
  int sy=this->Arrayky[id];
  int sz=this->Arraykz[id];
  int vt=this->VecArray[id];
  // MAS correction for each Fourier mode
  real_prec corr=1.0/this->Array_corr[id];
  int ngrid_two=ngrid/2+1;
  int n_ip, ip;
  switch(vt){
  case(0): //The zero mode
    data_ks[0][0]=1;
    data_ks[0][1]=0;
    data_dk[0][0]=this->data_out_g[0][0]*corr;
    data_dk[0][1]=this->data_out_g[0][1]*corr;
    data_pk_shot_noise[0][0]=pow(data_dk[0][0],2)+pow(data_dk[0][1],2);
    data_pk_shot_noise[0][1]=0;
    break;
  case(6)://Mapping vectors along the poles
    // Positive freqs;
    // 00z
    ip=ijk(0,0,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); // Original index. z-axis >0
    n_ip=ijk(0,0,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    ip=ijk(0,sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); //y-axis
    n_ip=ijk(0,sz,0,ngrid,ngrid,ngrid_two); //y-axis
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    ip=ijk(sz,0,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); //x-axis
    n_ip=ijk(sz,0,0,ngrid,ngrid,ngrid_two); //x-axis
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    // Negative freqs:
    // 0-z0
    ip  =ijk(0,this->params._Nft()-sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); //y-axis
    n_ip=ijk(0,ngrid-sz,0,ngrid,ngrid,ngrid_two); //y-axis
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    // -z00
    ip=ijk(this->params._Nft()-sz,0,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); //x-axis
    n_ip=ijk(ngrid-sz,0,0,ngrid,ngrid,ngrid_two); //x-axis
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    break;
  case(8): //x=y=z
    // 1: +++
    ip=ijk(sz,sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); // Original index. z-axis >0
    n_ip=ijk(sz,sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 1:   "<<sz<<"  "<<sz<<"   "<<sz<<endl;

    //2: +-+
    ip=ijk(sz,this->params._Nft()-sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); // Original index. z-axis >0
    n_ip=ijk(sz,ngrid-sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 2:   "<<sz<<"  "<<this->params._Nft()-sz<<"   "<<sz<<endl;

    //3: -++
    ip=ijk(this->params._Nft()-sz,sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); // Original index. z-axis >0
    n_ip=ijk(ngrid-sz,sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //cout<<"case 3:   "<<this->params._Nft()-sz<<"  "<<sz<<"   "<<sz<<endl;

    // 4: --+
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1); // Original index. z-axis >0
    n_ip=ijk(ngrid-sz,ngrid-sz,sz,ngrid,ngrid,ngrid_two); //New index. z-axis >0
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sz<<"   "<<sz<<endl;

    break;

  case(12):     // x=0, y=z

    // ++0
    ip=ijk(sz,sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 1:   "<<sz<<"  "<<sz<<"   "<<0<<endl;

    // +0+
    ip=ijk(sz,0,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,0,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 2:   "<<sz<<"  "<<0<<"   "<<sz<<endl;

    // 0++
    ip=ijk(0,sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 3:   "<<0<<"  "<<sz<<"   "<<sz<<endl;

    // 0-+
    ip=ijk(0,this->params._Nft()-sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,ngrid-sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<0<<"  "<<this->params._Nft()-sz<<"   "<<sz<<endl;


    // -0+
    ip=ijk(this->params._Nft()-sz,0,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,0,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 5:   "<<this->params._Nft()-sz<<"  "<<0<<"   "<<sz<<endl;



    // +-0
    ip=ijk(sz,this->params._Nft()-sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 6:   "<<sz<<"  "<<this->params._Nft()-sz<<"   "<<0<<endl;


    //-+0
    ip=ijk(this->params._Nft()-sz,sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 7:   "<<this->params._Nft()-sz<<"  "<<sz<<"   "<<0<<endl;


    //--0
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 8:   "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sz<<"   "<<0<<endl;


    break;

  case(23): //kx=0, kx!=ky, ky!=kz

    //1: 0yz
    ip=ijk(0,sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 1:   "<<0<<"  "<<sy<<"  "<<sz<<endl;



    //2: 0zy
    ip=ijk(0,sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 2:   "<<0<<"  "<<sz<<"  "<<sy<<endl;


    //3: y0z
    ip=ijk(sy,0,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,0,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 3:   "<<sy<<"  "<<0<<"  "<<sz<<endl;


    //4: yz0
    ip=ijk(sy,sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 4:   "<<sy<<"  "<<sz<<"  "<<0<<endl;


    //5: z0y
    ip=ijk(sz,0,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,0,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 5:   "<<sz<<"  "<<0<<"  "<<sy<<endl;

    //6: zy0
    ip=ijk(sz,sy,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sy,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 6:   "<<sz<<"  "<<sy<<"  "<<0<<endl;


    //7: -y-z0
    ip=ijk(this->params._Nft()-sy,this->params._Nft()-sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,ngrid-sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 7:   "<<this->params._Nft()-sy<<"  "<<this->params._Nft()-sz<<"  "<<0<<endl;


    //8: -z-y 0
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sy,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sy,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8:   "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sy<<"  "<<0<<endl;


    //9: 0-yz
    ip=ijk(0,this->params._Nft()-sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,ngrid-sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;

    //    cout<<"case 9:   "<<0<<"  "<<this->params._Nft()-sy<<"  "<<sz<<endl;


    //10: z-y0
    ip=ijk(sz,this->params._Nft()-sy,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sy,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10:   "<<sz<<"  "<<this->params._Nft()-sy<<"  "<<0<<endl;


    //11: -yz0
    ip=ijk(this->params._Nft()-sy,sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11:   "<<this->params._Nft()-sy<<"  "<<sz<<"  "<<0<<endl;


    //12: -y0z
    ip=ijk(this->params._Nft()-sy,0,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,0,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12:   "<<this->params._Nft()-sy<<"  "<<0<<"  "<<sz<<endl;

    //13: -z0y
    ip=ijk(this->params._Nft()-sz,0,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,0,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 13:   "<<this->params._Nft()-sz<<"  "<<0<<"  "<<sy<<endl;


    //14: y-z0
    ip=ijk(sy,this->params._Nft()-sz,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,ngrid-sz,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 14:   "<<sy<<"  "<<this->params._Nft()-sz<<"  "<<0<<endl;

    //15: 0-zy
    ip=ijk(0,this->params._Nft()-sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(0,ngrid-sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 15:   "<<0<<"  "<<this->params._Nft()-sz<<"  "<<sy<<endl;


    //16: -zy0
    ip=ijk(this->params._Nft()-sz,sy,0,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sy,0,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 16:   "<<this->params._Nft()-sz<<"  "<<sy<<"  "<<0<<endl;
    break;


  case(24): //x=y!=z not zero

    //1:xxz
    ip=ijk(sx,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1:   "<<sx<<"  "<<sx<<"  "<<sz<<endl;



    //2: xzx
    ip=ijk(sx,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2:   "<<sx<<"  "<<sz<<"  "<<sx<<endl;


    //3: zxx
    ip=ijk(sz,sx,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sx,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3:   "<<sz<<"  "<<sx<<"  "<<sx<<endl;


    //4: z-xx
    ip=ijk(sz,this->params._Nft()-sx,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sx,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4   "<<sz<<"  "<<this->params._Nft()-sx<<"  "<<sx<<endl;
    //5: x-xz
    ip=ijk(sx,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,ngrid-sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5   "<<sx<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;
    //6: -xzx
    ip=ijk(this->params._Nft()-sx,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6   "<<this->params._Nft()-sx<<"  "<<sz<<"  "<<sx<<endl;
    //7: -zxx
    ip=ijk(this->params._Nft()-sz,sx,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sx,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7   "<<this->params._Nft()-sz<<"  "<<sx<<"  "<<sx<<endl;
    //8: x-zx
    ip=ijk(sx,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8   "<<sx<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;
    //9: -xxz
    ip=ijk(this->params._Nft()-sx,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,sx, sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9   "<<this->params._Nft()-sx<<"  "<<sx<<"  "<<sz<<endl;

    //10: -z-xx
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sx,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sx,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10  "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sx<<"  "<<sx<<endl;
    //11: -x-xz
    ip=ijk(this->params._Nft()-sx,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,ngrid-sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->params._Nft()-sx<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;
    //12: -x-zx
    ip=ijk(this->params._Nft()-sx,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->params._Nft()-sx<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;
    break;
  case(25): //x!y y=z
    //1 : x z z
    ip=ijk(sx,sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1  "<<sx<<"  "<<sz<<"  "<<sz<<endl;
    //2 : z x z
    ip=ijk(sz,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2  "<<sz<<"  "<<sx<<"  "<<sz<<endl;
    //3 : z z x
    ip=ijk(sz,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3  "<<sz<<"  "<<sz<<"  "<<sx<<endl;
    //4 : z -x z
    ip=ijk(sz,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4  "<<sz<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;
    //5 : -z x z
    ip=ijk(this->params._Nft()-sz,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5  "<<this->params._Nft()-sz<<"  "<<sx<<"  "<<sz<<endl;



    //6 : x -z z
    ip=ijk(sx,this->params._Nft()-sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,ngrid-sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6  "<<sx<<"  "<<this->params._Nft()-sz<<"  "<<sz<<endl;



    //7 : -z z x
    ip=ijk(this->params._Nft()-sz,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7  "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;


    //8 : -x z z
    ip=ijk(this->params._Nft()-sx,sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8  "<<this->params._Nft()-sx<<"  "<<sz<<"  "<<sz<<endl;


    //9 : z -z x
    ip=ijk(sz,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9  "<<sz<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;

    //10 : -z -x z
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10  "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;


    //11 : -x -z z
    ip=ijk(this->params._Nft()-sx,this->params._Nft()-sz,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,ngrid-sz,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->params._Nft()-sx<<"  "<<this->params._Nft()-sz<<"  "<<sz<<endl;

    //12 : -z -zx
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->params._Nft()-sz<<"  "<<<this->params._Nft()-sz<<"  "<<sx<<endl;
    break;
  case(48): //x!=y!=z
    //1 :xyz
    ip=ijk(sx,sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 1  "<<sx<<"  "<<sy<<"  "<<sz<<endl;
    //2 :xzy
    ip=ijk(sx,sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 2  "<<sx<<"  "<<sz<<"  "<<sy<<endl;


    //3 :yxz
    ip=ijk(sy,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 3  "<<sy<<"  "<<sx<<"  "<<sz<<endl;


    //4 :yzx
    ip=ijk(sy,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 4  "<<sy<<"  "<<sz<<"  "<<sx<<endl;


    //5 :zxy
    ip=ijk(sz,sx,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sx,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 5  "<<sz<<"  "<<sx<<"  "<<sy<<endl;


    //6 :zyx
    ip=ijk(sz,sy,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,sy,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 6  "<<sz<<"  "<<sy<<"  "<<sx<<endl;
    //7 :x-yz
    ip=ijk(sx,this->params._Nft()-sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,ngrid-sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 7  "<<sx<<"  "<<this->params._Nft()-sy<<"  "<<sz<<endl;


    //8 :-xyz
    ip=ijk(this->params._Nft()-sx,sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 8  "<<this->params._Nft()-sx<<"  "<<sy<<"  "<<sz<<endl;


    //9 :-x-yz
    ip=ijk(this->params._Nft()-sx,this->params._Nft()-sy,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,ngrid-sy,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 9  "<<this->params._Nft()-sx<<"  "<<this->params._Nft()-sy<<"  "<<sz<<endl;


    //10 :x-zy
    ip=ijk(sx,this->params._Nft()-sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sx,ngrid-sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 10  "<<sx<<"  "<<this->params._Nft()-sz<<"  "<<sy<<endl;
    //11 :-xzy
    ip=ijk(this->params._Nft()-sx,sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 11  "<<this->params._Nft()-sx<<"  "<<sz<<"  "<<sy<<endl;
    //12 :-x-zy
    ip=ijk(this->params._Nft()-sx,this->params._Nft()-sz,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sx,ngrid-sz,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 12  "<<this->params._Nft()-sx<<"  "<<this->params._Nft()-sz<<"  "<<sy<<endl;
    //13 :y-xz
    ip=ijk(sy,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,ngrid-sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 13  "<<sy<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;
    //14 :-yxz
    ip=ijk(this->params._Nft()-sy,sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,sx,sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 14  "<<this->params._Nft()-sy<<"  "<<sx<<"  "<<sz<<endl;
    //15 :-y-xz
    ip=ijk(this->params._Nft()-sy,this->params._Nft()-sx,sz,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,ngrid-sx, sz,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 15  "<<this->params._Nft()-sy<<"  "<<this->params._Nft()-sx<<"  "<<sz<<endl;
    //16 :-yzx
    ip=ijk(this->params._Nft()-sy,sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 16  "<<this->params._Nft()-sy<<"  "<<sz<<"  "<<sx<<endl;
    //17 :y-zx
    ip=ijk(sy,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sy,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 17  "<<sy<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;
    //18 :-y-zx
    ip=ijk(this->params._Nft()-sy,this->params._Nft()-sz,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sy,ngrid-sz,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 18  "<<this->params._Nft()-sy<<"  "<<this->params._Nft()-sz<<"  "<<sx<<endl;
    //19 :z-xy
    ip=ijk(sz,this->params._Nft()-sx,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sx,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 19  "<<sz<<"  "<<this->params._Nft()-sx<<"  "<<sy<<endl;
    //20 :-zxy
    ip=ijk(this->params._Nft()-sz,sx,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sx,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 20  "<<this->params._Nft()-sz<<"  "<<sx<<"  "<<sy<<endl;
    //21 :-z-xy
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sx,sy,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sx,sy,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 21  "<<this->params._Nft()-sz<<"  "<<this->params._Nft()-sx<<"  "<<sy<<endl;
    //22 :-zyx
    ip=ijk(this->params._Nft()-sz,sy,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,sy,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 22  "<<this->params._Nft()-sz<<"  "<<sy<<"  "<<sx<<endl;
    //23 :z-yx
    ip=ijk(sz,this->params._Nft()-sy,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(sz,ngrid-sy,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    //    cout<<"case 23  "<<sz<<"  "<<this->params._Nft()-sy<<"  "<<sx<<endl;
    //24 :-z-yx
    ip=ijk(this->params._Nft()-sz,this->params._Nft()-sy,sx,this->params._Nft(),this->params._Nft(),this->params._Nft()/2+1);
    n_ip=ijk(ngrid-sz,ngrid-sy,sx,ngrid,ngrid,ngrid_two);
    data_ks[n_ip][0]=1;
    data_ks[n_ip][1]=0;
    data_dk[n_ip][0]=this->data_out_g[ip][0]*corr;
    data_dk[n_ip][1]=this->data_out_g[ip][1]*corr;
    data_pk_shot_noise[n_ip][0]=pow(data_dk[n_ip][0],2)+pow(data_dk[n_ip][1],2);
    data_pk_shot_noise[n_ip][1]=0;
    break;
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::get_bispectrum_fkp_fast(vector<real_prec> &pk, vector<real_prec> &bispect, vector<int> &mod, string file)
{
  // Define the shells:
  define_kshells();
 // Get the inverse FT in the defined shells
  get_ift_shells_bispectrum();
  // Loop pver k-shells and estimates of bispectrum
  loop_shells_bispectrum(pk, bispect, mod, file);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::write_fftw_parameters()
{
#ifdef _FULL_VERBOSE_
  So.message_screen("Input values and parameters of the FFT:");
  So.message_screen("Dimension of the grid =",this->params._Nft(),"³");
  So.message_screen("Lenght =",this->params._Lbox()," Mpc/h");
  if(this->params._type_of_binning()=="linear")
    {
      So.message_screen("Using linearly-spaced spherical shells");
      So.message_screen("MAS               =", this->params._mass_assignment_scheme());
      So.message_screen("eMAS              =", this->params._mass_assignment());
      So.message_screen("kmin              =",this->params._d_kmin()," h/Mpc");
      So.message_screen("Bin size for P(k) = ",this->params._d_DeltaK_data()," h/Mpc");
      So.message_screen("Bin size for W(k) = ",this->params._d_DeltaK_window()," h/Mpc");
      So.message_screen("Nyquist Frequency = ",this->params._d_nyquist_frequency()," h / Mpc");
    }
  else{
    So.message_screen("Using log-spaced spherical shells ");
    So.message_screen("Bin size in log =", this->params._d_Deltal());
  }
  std::cout<<RESET;
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::write_fftw_parameters(void *p, string fname)
{
  struct s_parameters_box * s_cp= (struct s_parameters_box *)p;
  ofstream out;
/*
  ifstream html(this->params._Output_directory()+"log_power.html");
  if (!html)
         cout << "Error opening"<<"  "<<this->params._Output_directory()+"log_power.html"<<endl;
 */
  out.open(fname.c_str() , ios::out);
  out.precision(12);
  out.setf(ios::showpoint);
  out.width(12);
  out<<"Selected options :"<<endl;
  out<<"Statistics: "<<this->params._statistics()<<endl;
  out<<"MAS: "<<(this->params._mass_assignment_scheme())<<endl;
  if(this->params._MAS_correction())
      out<<"MAS correction: enabled"<<endl;
  else
      out<<"MAS correction: disabled"<<endl;
  if(this->params._FKP_weight())
      out<<"Using FKP weights with Pest = "<<this->params._Pest()<<endl;
  else
      out<<"Using weights = 1"<<endl;
  if(this->params._FKP_error_bars())out<<"Computing FKP error bars"<<endl;
  else out<<"Estimate without error bars."<<endl;
  if(this->params._SN_correction())out<<"Shot-noise correction: enabled"<<endl;
  else out<<"Shot-noise correction: disabled.  "<<endl;
  out<<"******************************************************************************"<<endl;
  out<<"Sample information"<<endl;
  out<<"Number of objects = "<<this->n_gal<<endl;
  if(true==this->params._use_random_catalog())
    out<<"Number of random objects = "<<this->n_ran<<endl;
  out<<"******************************************************************************"<<endl;
  out<<"Fourier space information :"<<endl;
  out<<"Lbox = "<<this->params._Lbox()<<"  Mpc /h"<<endl;
  out<<"Nft = "<<this->params._Nft()<<"³"<<endl;
  if(this->params._statistics()=="Pk_y_ds")out<<"Number of shells to kmax = "<<this->sgrid<<std::endl;
  if(this->params._type_of_binning()=="linear")
    {
      out<<"Using linearly-spaced spherical shells "<<endl;
      out<<"Bin size for P(k) = "<<this->params._d_DeltaK_data()<<" h/Mpc"<<endl;
      out<<"Bin size for W(k) = "<<this->params._d_DeltaK_window()<<" h/Mpc"<<endl;
      out<<"Nyquist Frequency = "<<0.5*this->params._Nft()*(this->params._d_deltak_0())<<" h/Mpc"<<endl;
    }
  else
    {
      out<<"Using log-spaced spherical shells"<<endl;
      out<<"Bin size = "<<this->params._d_Deltal()<<endl;
    }
  out<<"****************************************************************************"<<endl;
  out<<"Parameters of the estimator: "<<endl;
  out<<"Normalization = "<<this->normal_power<<endl;
  out<<"Shot_noise (power) = "<<this->shot_noise<<"  (Mpc h^-1)^3"<<endl;
  if(this->params._use_random_catalog())
    out<<"Shot_noise (window) = "<<this->shot_noise_window<<endl;
  if(true==this->params._use_random_catalog())
    out<<"Mean number density (from weights) = "<<(this->n_ran-this->w_r)/(this->params._Pest()*this->w_r)<<" (Mpc h^-1)^(-3)"<<endl;
  out<<"Weighted number of objects = "<<this->w_g<<endl;
  if(true==this->params._use_random_catalog())
      out<<"Weighted number of random objects = "<<this->w_r<<endl;
  out<<"alpha = "<<this->alpha<<endl;
  out<<"Sum nw² galaxies = "<<this->s_g<<endl;
  if(this->params._use_random_catalog())out<<"Sum nw² random = "<<this->s_r<<endl;
  out<<"******************************************************************************"<<endl;
  time_t rawtime;
  time ( &rawtime );
  out<<"Date: "<<ctime (&rawtime)<<endl;
  out.close();
  So.message_screen("Log-file written in ",fname);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::aux_func_yamamoto(ULONG lp, real_prec kx, real_prec ky, real_prec kz, real_prec icorr2, real_prec &F1,real_prec &F2)
{
    real_prec kx2=kx*kx;
    real_prec ky2=ky*ky;
    real_prec kz2=kz*kz;
    real_prec kv2=kx2+ky2+kz2;
    real_prec F2_r=0.5*3.*((kx2*this->data_out_g_xx[lp][REAL])+ (ky2*this->data_out_g_yy[lp][REAL])
                          + (kz2*this->data_out_g_zz[lp][REAL])+(2.*kx*ky*this->data_out_g_xy[lp][REAL])
                          + (2.*kx*kz*this->data_out_g_xz[lp][REAL])+ (2.*ky*kz*this->data_out_g_yz[lp][REAL]))/kv2  - 0.5*this->data_out_g[lp][REAL];
    real_prec F2_i=0.5*3.*((kx2*this->data_out_g_xx[lp][IMAG]) + (ky2*this->data_out_g_yy[lp][IMAG])
                           + (kz2*this->data_out_g_zz[lp][IMAG]) + (2.*kx*ky*this->data_out_g_xy[lp][IMAG])
                           + (2.*kx*kz*this->data_out_g_xz[lp][IMAG]) + (2.*ky*kz*this->data_out_g_yz[lp][IMAG]))/kv2  - 0.5*this->data_out_g[lp][IMAG];

    F1 =  (F2_r*this->data_out_g[lp][REAL] + F2_i*this->data_out_g[lp][IMAG]) * icorr2;


    if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_ysc")
      F2= (F2_r * F2_r + F2_i*F2_i) * icorr2; //use this for Scoccimarro's version of Hexadecapole in the the Yamamoto-Blake estimator

    else if(this->params._statistics()=="Pk_yb"  ||  this->params._statistics()=="Pk_ybc")
      {
        real_prec F4_r = (35./8.)*((kx2*kx2)*this->data_out_g_xxx[lp][REAL]+(ky2*ky2)*this->data_out_g_yyy[lp][REAL]
                                   +(kz2*kz2)*this->data_out_g_zzz[lp][REAL]+4.*(kx2*kx)*ky*this->data_out_g_xxy[lp][REAL]
                                   +4.*(kx2*kx)*kz*this->data_out_g_xxz[lp][REAL]+4.*(ky2*ky)*kx*this->data_out_g_yyx[lp][REAL]
                                   +4.*(ky2*ky)*kz*this->data_out_g_yyz[lp][REAL]+4.*(kz2*kz)*kx*this->data_out_g_zzx[lp][REAL]
                                   +4.*(kz2*kz)*ky*this->data_out_g_zzy[lp][REAL]+6.*(kx2)*(ky2)*this->data_out_g_xyy[lp][REAL]
                                   +6.*(kx2)*(kz2)*this->data_out_g_xzz[lp][REAL]+6.*(ky2)*(kz2)*this->data_out_g_yzz[lp][REAL]
                                   +12.*kx*ky*kz*( kx* this->data_out_g_xyz[lp][REAL] + ky*this->data_out_g_yxz[lp][REAL]+kz*this->data_out_g_zxy[lp][REAL]))/(kv2*kv2)
          -(15./6.)*F2_r-(7./8.)*this->data_out_g[lp][REAL];
        real_prec F4_i = (35./8.)*((kx2*kx2)*this->data_out_g_xxx[lp][IMAG]+(ky2*ky2)*this->data_out_g_yyy[lp][IMAG]
                                   +(kz2*kz2)*this->data_out_g_zzz[lp][IMAG]+4.*(kx2*kx)*ky* this->data_out_g_xxy[lp][IMAG]
                                   +4.*(kx2*kx)*kz* this->data_out_g_xxz[lp][IMAG]+4.*(ky2*ky)*kx* this->data_out_g_yyx[lp][IMAG]
                                   +4.*(ky2*ky)*kz* this->data_out_g_yyz[lp][IMAG]+4.*(kz2*kz)*kx* this->data_out_g_zzx[lp][IMAG]
                                   +4.*(kz2*kz)*ky* this->data_out_g_zzy[lp][IMAG]+6.*(kx2)*(ky2)*this->data_out_g_xyy[lp][IMAG]
                                   +6.*(kx2)*(kz2)* this->data_out_g_xzz[lp][IMAG]+6.*(ky2)*(kz2)* this->data_out_g_yzz[lp][IMAG]
                                   +12.*kx*ky*kz*( kx* this->data_out_g_xyz[lp][IMAG] + ky*this->data_out_g_yxz[lp][IMAG]+kz*this->data_out_g_zxy[lp][IMAG]))/(kv2*kv2)
          -(15./6.)*F2_i-(7./8.)*this->data_out_g[lp][IMAG];

        F2=(F4_r*this->data_out_g[lp][REAL] + F4_i*this->data_out_g[lp][IMAG])*icorr2;
      }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::free_fftw_vectors()
{
#ifdef DOUBLE_PREC
  if(sizeof(this->data_out_g)>0)
    fftw_free(this->data_out_g);
  if(sizeof(this->data_out_r)>0)
    fftw_free(this->data_out_r);
  if(sizeof(this->data_out_g_rss)>0)
    fftw_free(this->data_out_g_rss);
   fftw_free(this->data_out_gp);
 //  fftw_free(this->data_out_g_gp);
#else
  if(sizeof(this->data_out_g)>0)
      fftwf_free(this->data_out_g);
  if(true==this->params._use_random_catalog())
    if(sizeof(this->data_out_r)>0)
      fftwf_free(this->data_out_r);
  if(this->params._measure_cross()==true)
      if(sizeof(this->data_out_gp)>0)
          fftwf_free(this->data_out_gp);
  if(this->params._use_real_and_redshift_space()==true)
      if(sizeof(this->data_out_g_rss)>0)
          fftwf_free(this->data_out_g_rss);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void FftwFunctions::set_params(Params _params){
    this->params=_params;
    this->rdeltax = static_cast<real_prec> (1.0/this->params._d_delta_x());
    this->rdeltay = this->rdeltax;
    this->rdeltaz = this->rdeltax;
    this->data_g.clear();data_g.shrink_to_fit();
    this->data_g_rss.clear();data_g_rss.shrink_to_fit();
    this->DeltaK_Bis    = this->params._d_DeltaK_data();// (this->params._kmax_bk()-this->params._kmin_bk())/((real_prec)this->Nshells_bk)
    this->Nshells_bk    = static_cast<int>(this->params._kmax_bk()/this->DeltaK_Bis);
    ULONG nff;
     if(this->params._statistics()!="Pk_y_ds")
       nff=this->params._Nft();
      // Get the number of grid-cells for the direct sum DFT as a function of kmax
     else if(this->params._statistics()=="Pk_y_ds")
       {
        nff=(int)(this->params._kmax_y_ds()*this->params._Lbox()/M_PI);
        if(nff%2!=0)nff++;
        this->params.set_Nft(nff);
        this->sgrid=nff;
      }
     // Get number of grid cells per dimension used for the bispectrum_fast given kmax
     int n_sgrid= static_cast<int>(this->params._kmax_bk()*this->params._Lbox()/M_PI);
     if(this->params._statistics()=="Bk_fkp_fast")
       {
        if(n_sgrid>this->params._Nft())
          {
           So.message_warning("Warning: kmax greater than Nyquist frequency");
           So.message_screen("setting kmax to the Nyquist");
           n_sgrid=this->params._Nft();
           this->params.set_kmax_bk(this->params._Nft()*M_PI/this->params._Lbox() );
         }
        if(n_sgrid%2!=0)n_sgrid++;
        this->sgrid=n_sgrid;
        if(this->params._kmin_bk()<2.*M_PI/this->params._Lbox() )
          {
            So.message_warning("Warning: kmin smaller than fundamental mode");
            So.message_screen("setting kmin as fundamental mode");
            this->params.set_kmin_bk(2.*M_PI/this->params._Lbox()) ;
          }
        }
      if(this->params._statistics()=="Bk_fkp_fast")
        {
          // Get effective number of grid cells
          // resulting from looping over half of the positive quadrant.
          // The result will be used to allocate memory for the arrays
          // used in the Bispectrum as performed by Jennifer.
          int new_sd=0;
          for(int i=0;i<this->sgrid;++i)
            for(int j=i;j<this->sgrid;++j)
              for(int k=j;k<this->sgrid;++k)
                if(i*i+j*j+k*k>0)
                  new_sd++;
          this->new_sgrid=new_sd;
        }
      // Initialize other private variables for the Pk
      this->Nft2=(int)(nff*nff);
      this->rNft =1.0 / this->params._Nft();
      this->params.set_NGRID(nff*nff*nff);
      this->params.set_NGRID_h(nff*nff*(nff/2+1));
      time_t time_bam;
      time(&time_bam);
      this->So.initial_time=time_bam;
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
