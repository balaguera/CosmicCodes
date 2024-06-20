////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/** @file BiasMT.cpp
 *  @brief Generation of mock catalogs of DM tracers
 *  based on the Bias Mapping Tecnique (ex-BAM) (Balaguera, Kitaura)
 *  @author: Andrés Balaguera-Antolínez
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#include "../headers/BiasMappingTechnique.h"
#include "../headers/def.h"
using namespace std;
#ifdef _USE_PYTHON_
void py_plot(vector<real_prec>& field)
{
  return;
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::warnings(){
#if !defined (_USE_DM_IN_BiasMT_)
      So.message_warning("Warning: _USE_DM_IN_BiasMT_ is an undefined pre-proc directive");
      So.message_warning("If this code is meant to run with DM (default option), please define it in def.h.");
      So.message_warning("Otherwise, type c to continue");
      string cont;
      cin>>cont;
      if (cont!="c" || cont!="C")
        So.message_warning("COSMICATLAS stops here.");
      exit(0);
#endif
#if !defined (_DO_BiasMT_CALIBRATION_) && defined (_MODIFY_LIMITS_)
    So.message_warning("Warning: _MODIFY_LIMITS_ is defined under directive _GET_BiasMT_REALIZATIONS_.");
    So.message_warning("This is a potential source of  bug/collapse of the code. Please check it.");
    exit(0);
#endif
#if defined (_USE_NABLA2DELTA_) && defined (_USE_INVARANT_SHEAR_VFIELD_I_)
    So.message_warning("Warning: _USE_NABLA2DELTA_ and _USE_INVARANT_SHEAR_VFIELD_I_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2DELTA_) && defined (_USE_INVARANT_SHEAR_VFIELD_II_)
    So.message_warning("Warning: _USE_S2DELTA_ and _USE_INVARANT_SHEAR_VFIELD_II_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S3_) && defined (_USE_INVARANT_SHEAR_VFIELD_III_)
    So.message_warning("Warning: _USE_S3_ and _USE_INVARANT_SHEAR_VFIELD_III_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2_) && defined (_USE_TIDAL_ANISOTROPY_)
    So.message_warning("Warning: _USE_S2_ and _USE_TIDAL_ANISOTROPY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_S2_) && defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
    So.message_warning("Warning: _USE_S2_ and _USE_INVARIANT_TIDAL_FIELD_IV_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_INVARIANT_TIDAL_FIELD_IV_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_TIDAL_ANISOTROPY_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_PROLATNESS_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_PROLATNESS_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif

#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_S2_) && defined (_USE_ELLIPTICITY_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_ELLIPTICITY_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_S2_) && defined (_USE_PROLATNESS_)
    So.message_warning("Warning: _USE_TIDAL_ANISOTROPY_ and _USE_PROLATNESS_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_DELTA3_) && defined (_USE_INVARIANT_TIDAL_FIELD_III_)
    So.message_warning("Warning: _USE_DELTA3_ and _USE_INVARIANT_TIDAL_FIELD_III_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
#if defined (_USE_DELTA2_) && defined (_USE_INVARIANT_TIDAL_FIELD_II_)
    So.message_warning("Warning: _USE_DELTA2_ and _USE_INVARIANT_TIDAL_FIELD_II_ are defined.");
    So.message_warning("These two properties use the smae memory slot in BAM. Please undefine one, according to your needs.");
    exit(0);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_cosmo()
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Getting cosmological derived-parameters");
  So.message_screen("Cosmological redshift z =",this->params._redshift());
#endif
  // Here we fill the structure s_cosmo_info with different cosmological quantities evaluated at the input redshift
  this->s_cosmo_info.scale_factor=1./(this->params._redshift()+1.);
  this->s_cosmo_info.critical_density=this->Cosmo.critical_overdensity(this->params._redshift());
  this->s_cosmo_info.density_contrast_top_hat=Cosmo.density_contrast_top_hat(this->params._redshift());
  this->s_cosmo_info.Hubble_parameter=this->Cosmo.Hubble_function(this->params._redshift());
  this->s_cosmo_info.comoving_distance=this->Cosmo.comoving_distance(this->params._redshift());
  this->s_cosmo_info.comoving_angular_diameter_distance=this->Cosmo.comoving_angular_diameter_distance(this->params._redshift());
  this->s_cosmo_info.mean_matter_density=this->Cosmo.mean_matter_density(this->params._redshift());
  this->s_cosmo_info.age_universe=Cosmo.age_universe(this->params._redshift());
  this->s_cosmo_info.comoving_sound_horizon=Cosmo.comoving_sound_horizon(this->params._redshift());
  this->s_cosmo_info.Delta_Vir = Cosmo.density_contrast_top_hat(this->params._redshift());
#if defined (_USE_LPT_ ) || defined (_DISPLACEMENTS_) // This is repeated here as it is in LPT
  this->s_cosmo_info.growth_factor=this->Cosmo.growth_factor(this->params._redshift(), (void *)&this->s_cosmo_pars)/static_cast<double>(this->Cosmo.growth_factor(0.0,(void *)&this->s_cosmo_pars));
  this->s_cosmo_info.growth_index=this->Cosmo.growth_index(this->params._redshift(), (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.growth_index2=this->Cosmo.growth_index2(this->params._redshift(), (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.halo_dynamical_time=Cosmo.halo_dynamical_time(this->params._redshift(), (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.omega_matter=Cosmo.omega_matter(this->params._redshift(), (void *)&this->s_cosmo_pars);
  this->s_cosmo_info.Distance_Modulus=Cosmo.Distance_Modulus(this->params._redshift(), (void *)&this->s_cosmo_pars);
  this->s_cosmo_pars.growth_factor=this->s_cosmo_info.growth_factor;
  real_prec fD2=static_cast<real_prec>(pow(this->s_cosmo_info.omega_matter,-1./143.));
  this->s_cosmo_info.D2=static_cast<real_prec>(-(3./7.)*pow(this->s_cosmo_info.growth_factor,2)*fD2);
  PowerSpectrum Pow;
  this->s_cosmo_pars.growth_factor=this->s_cosmo_info.growth_factor;
  this->s_cosmo_pars.pk_normalization=Pow.normalization((void *)&this->s_cosmo_pars);
#endif
#ifdef _USE_LPT_  // This is repeated here as it is in LPT
#ifdef _FULL_VERBOSE_
  this->So.write_cosmo_parameters((void *)&this->s_cosmo_pars, (void *)&this->s_cosmo_info);
#else
    this->So.message_screen("Cosmological redshift z =", this->s_cosmo_pars.cosmological_redshift);
#endif
#endif
  So.DONE();
 this->cwclass.s_cosmo_info=this->s_cosmo_info;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::set_Fourier_vectors()
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Inizializing vectors");
#endif
  ULONG NTT=(this->params._Nft())*(this->params._Nft())*(this->params._Nft()/2+1);
  this->Kernel.resize(NTT, 1.0); //initialize to unity
  this->Power_REF.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
  this->Power_NEW.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
  this->power_ratio_unsmoothed.resize(this->params._Nft()/2/this->params._ndel_data(), 10.0);
  this->kvec.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_MASS_KNOTS_
  this->SKNOT_M_info.resize( this->params._NGRID(), 0);
#endif
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_power_spectrum(string type)
{
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  So.message_screen("Measuring power spectrum of",type);
#endif
  this->params.set_mass_assignment_scheme("NGP");
  this->params.set_MAS_correction(false);
  this->params.set_dir_output(this->params._Output_directory());  //This line is important
  // in the case in which the output dir is changing. Otherwise, the O(k) will be written in the dir read initially by the class Params
  this->params.set_Name_survey(type);
  if(this->iteration <=this->params._N_iterations_Kernel())
    this->params.set_Name_survey(type+"_iteration"+to_string(this->iteration));
  if(type=="DM_REF" || type=="DM_REF_NEW")
    {
#ifdef _DISPLACEMENTS_
      this->params.set_input_type("delta_grid");
#else
      this->params.set_input_type("density_grid");
#endif
      this->params.set_SN_correction(false);
      this->params.set_MAS_correction(true);
      if(0==this->params._iMAS_X())
        {
          this->params.set_mass_assignment_scheme("NGP");
         this->params.set_MAS_correction(false);
        }
      if(1==this->params._iMAS_X()){
        this->params.set_mass_assignment_scheme("CIC");
      }
      if(11==this->params._iMAS_X())
           {
         this->params.set_mass_assignment_scheme("CIC");
         this->params.set_MAS_correction(false);
        }
      else if (2==this->params._iMAS_X())
        this->params.set_mass_assignment_scheme("TSC");
      else if (3==this->params._iMAS_X())
        this->params.set_mass_assignment_scheme("PSC");
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_X_ini, true);
      So.DONE();
      this->Power_DM_REF.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_DM_REF.size(); ++i)
        this->Power_DM_REF[i]=cPSF._pk0(i);
    }
  else if(type=="CROSS_TR_DM")
    {
      this->params.set_input_type("density_grid");
      this->params.set_SN_correction(true);
      this->params.set_mass_assignment_scheme("NGP");
      this->params.set_MAS_correction(false);
      this->params.set_measure_cross(true);
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_cross_power_spectrum_grid(true,this->delta_X, this->delta_Y, true);
      So.DONE();
      cPSF.write_power_and_modes();
    }
  else if(type=="DM_iteration" || type=="DM_real" || type=="DM_KONV"|| type=="DM_RO")
    {
      this->params.set_input_type("delta_grid");
      this->params.set_SN_correction(false);
      this->params.set_mass_assignment_scheme("CIC");
      this->params.set_MAS_correction(true);
      PowerSpectrumF cPSF(this->params);
      if(type=="DM_RO")
        {
              cPSF.compute_power_spectrum_grid(this->delta_X_ini, true);
             this->Power_DM_REF.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0;i<this->Power_DM_REF.size(); ++i)
            this->Power_DM_REF[i]=cPSF._pk0(i);
        }
      if(type=="DM_KONV")
        {
             cPSF.compute_power_spectrum_grid(this->delta_X,true);
             this->Power_DM_KONV.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0;i<this->Power_DM_KONV.size(); ++i)
            this->Power_DM_KONV[i]=cPSF._pk0(i);
        }
      else
        cPSF.compute_power_spectrum_grid(this->delta_X,true);
      So.DONE();
    }
  else if(type=="TR_MOCK" || type== "TR_MOCK_REALIZATION")//TR_MOCK refers to the itarations: TR_MOCK_REAL refers to the actuial mocks built from the learning process
    {
      if(this->params._Name_Property_Y()=="COUNTS")
        {
#ifdef _MG_
          this->params.set_SN_correction(false);
#else
          this->params.set_SN_correction(true);
#endif
          this->params.set_MAS_correction(false);
          this->params.set_mass_assignment_scheme("NGP");
        }
      else
        {
          this->params.set_SN_correction(false);
          this->params.set_mass_assignment_scheme("CIC");
          this->params.set_MAS_correction(true);
        }
      if(type== "TR_MOCK_REALIZATION")
        this->params.set_Name_survey("TR_MOCK");
#ifdef _CALIBRATION_WITHOUT_SN_
      if(_COUNTS_==this->params._Name_Property_Y() )
        if(type=="TR_MOCK")
#ifdef _MG_
            if(this->iteration<100)
#else
            if(this->iteration<2)
#endif
            this->params.set_SN_correction(false);
#endif
      this->params.set_input_type("density_grid");
#ifdef _DISPLACEMENTS_
      this->params.set_SN_correction(false);
      this->params.input_type="delta_grid";
#endif
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_Y_new,true);
      this->Power_NEW.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_NEW.size(); ++i)
        this->Power_NEW[i]=cPSF._pk0(i);  ///     +cPSF.shot_noise; Revisar esto: colocar arriba sn_cor =true, acá se los ponemos, pero lo escribimos con SN desde write_power_modes()
      So.DONE();
#ifndef _GET_BiasMT_REALIZATIONS_
      cPSF.write_power_and_modes(this->params._Output_directory()+"power_mock.txt");
#endif
      this->kvec.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->Power_NEW.size() ;++i)
        this->kvec[i]=cPSF._kvector_data(i);
    }
  else if(type=="TR_MOCK_CATb")
    {
      this->params.set_SN_correction(true);
      this->params.set_input_type("density_grid");
      this->params.set_mass_assignment_scheme("CIC");
      this->params.set_MAS_correction(true);
      PowerSpectrumF cPSF(this->params);
      cPSF.compute_power_spectrum_grid(this->delta_Y_new,true);
    }
  else if(type=="TR_MOCK_CAT")
    {
      this->params.set_SN_correction(false);
      if(this->params._Name_Property_Y()=="COUNTS")
        this->params.set_SN_correction(true);
      this->params.set_input_type("catalog");
      this->params.set_mass_assignment_scheme("TSC");
      this->params.set_MAS_correction(true);
#ifdef _USE_LPT_
      this->params.set_file_catalogue(this->lpt._fnameTRACERCAT());
#else
      this->params.set_file_catalogue("Auxfile");
#endif
      // Change this for the ordering in the param file might not be that of the one used by LPT to write the catalog
      this->params.set_i_coord1_g(0);
      this->params.set_i_coord2_g(1);
      this->params.set_i_coord3_g(2);
      /*      PowerSpectrumF cPSF(this->params);
              cPSF.tracer_cat=this->tracer;
              cPSF.compute_power_spectrum(false, false);
              So.DONE();
              //the argument false above forces to write explicitely here the write_power and modes:
              cPSF.write_power_and_modes();
      */
    }
  else if(type=="TR_REF")
    {
      this->params.set_SN_correction(false);
      if(_COUNTS_==this->params._Name_Property_Y())
        {
#ifdef _MG_
          this->params.set_SN_correction(false);
#else
          this->params.set_SN_correction(true);
#endif
          this->params.set_mass_assignment_scheme("NGP");
        this->params.set_MAS_correction(false);
        }
       else
         {
          this->params.set_mass_assignment_scheme("CIC");
           this->params.set_MAS_correction(false);
         }
      this->params.set_input_type("density_grid");
#ifdef _DISPLACEMENTS_
         this->params.set_SN_correction(false);
      this->params.input_type="delta_grid";
#endif
#ifdef _USE_TRACER_HR_
      this->params.set_Nft(this->params._Nft()_HR);
#endif
      PowerSpectrumF cPSF(this->params);
#ifdef _USE_TRACER_HR_
      cPSF.compute_power_spectrum_grid(this->delta_Y_HR,true);
#else
      cPSF.compute_power_spectrum_grid(this->delta_Y,true);
#endif
      // Return to the original value of Nft
      // and write ref power up the the NF of the Nft
#ifdef _USE_TRACER_HR_
      this->params.set_Nft(this->params._Nft());
#endif
      So.DONE();
      this->Power_REF.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._Nft()/2/this->params._ndel_data();++i)
        this->Power_REF[i]=cPSF._pk0(i);
      this->kvec.resize(this->params._Nft()/2/this->params._ndel_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._Nft()/2/this->params._ndel_data() ;++i)
        this->kvec[i]=cPSF._kvector_data(i);
    }
  // For the DM we do not take the kvectors. We take them from the calc of TR power
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define _use_more_chis_  //use this option to use the info from adjacen k-bins
void BiasMT::GetKernel(bool rejection, real_prec exponent)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _use_more_chis_
  real_prec weight_back=0.15;
  real_prec weight_forw=0.15;
  real_prec weight_central=0.7;
#endif
#ifdef _use_more_chis_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Applying Metropolis-Hasting algorithm to generate BAM-Kernel (three-modes based)");
#endif
#else
#ifdef _FULL_VERBOSE_
this->So.message_screen("Applying Metropolis-Hasting algorithm to generate BAM-Kernel (single-mode based)");
#endif
#endif
#ifdef _USE_GNUPLOT_
    // PLOT OF THE POWER SPECTRUM
     vector<pair<real_prec, real_prec> > xy_pts_ref;
     for(int i=1; i<kvec.size(); ++i)
       xy_pts_ref.push_back(std::make_pair(this->kvec[i], log10(this->Power_REF[i])));
     vector<pair<real_prec, real_prec> > xy_pts_new;
     for(int i=1; i<kvec.size(); ++i)
       xy_pts_new.push_back(std::make_pair(this->kvec[i], log10(this->Power_NEW[i])));
     vector<pair<real_prec, real_prec> > xy_pts_dm;
     for(int i=1; i<kvec.size(); ++i)
       xy_pts_dm.push_back(std::make_pair(this->kvec[i], log10(this->Power_DM_REF[i])));
     vector<pair<real_prec, real_prec> > xy_pts_dm_k;
     cout<<this->Power_DM_KONV.size()<<endl;
     if(this->iteration>1)
         for(int i=1; i<kvec.size(); ++i)
           xy_pts_dm_k.push_back(std::make_pair(this->kvec[i], log10(this->Power_DM_KONV[i])));
     this->gp_kernel<<"set size 0.3,0.5\n";
     this->gp_kernel<<"set origin 0.3,0.0\n";
     this->gp_kernel<<"set log x \n";
     this->gp_kernel<<"set grid \n";
     this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
     this->gp_kernel<<"set ylabel 'log P(k) [(Mpc / h)³]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
#ifdef _DISPLACEMENTS_
//     this->gp_kernel<<"plot"<<this->gp_kernel.file1d(xy_pts_ref) << "w l lw 3 lt 2 title 'Displacement',"<<this->gp_kernel.file1d(xy_pts_new) << " w l lw 3 lt 6 title 'New Displacement',"<<this->gp_kernel.file1d(xy_pts_dm)<< "w l lw 3 lt 3 title 'IC'"<<endl;
     this->gp_kernel<<"plot"<<this->gp_kernel.file1d(xy_pts_ref) << "w l lw 3 lt 2 title 'Ref. Displacement' ,"<<this->gp_kernel.file1d(xy_pts_new) << " w l lw 3 lt 6 title 'New Displacement'  "<<endl;
#else
     if(this->iteration>1)
         this->gp_kernel<<"plot"<<this->gp_kernel.file1d(xy_pts_ref) << "w l lw 4 lt 2 title 'Reference Halo' ,"<<this->gp_kernel.file1d(xy_pts_new) << " w l lw 4 lt 5 title 'Mock Halo' ,"<<this->gp_kernel.file1d(xy_pts_dm)<< "w l lw 4 lt 7 title 'DM', "<<this->gp_kernel.file1d(xy_pts_dm_k)<< "w l lw 4 lt 3 title 'DM Konv' "<<endl;
     else
         this->gp_kernel<<"plot"<<this->gp_kernel.file1d(xy_pts_ref) << "w l lw 4 lt 2 title 'Reference Halo' ,"<<this->gp_kernel.file1d(xy_pts_new) << " w l lw 4 lt 5 title 'Mock Halo' ,"<<this->gp_kernel.file1d(xy_pts_dm)<< "w l lw 4 lt 7 title 'DM'"<<endl;
#endif
     xy_pts_ref.clear();
     xy_pts_ref.shrink_to_fit();
     xy_pts_new.clear();
     xy_pts_new.shrink_to_fit();
     xy_pts_dm.clear();
     xy_pts_dm.shrink_to_fit();
     xy_pts_dm_k.clear();
     xy_pts_dm_k.shrink_to_fit();
#endif
  if(false==this->use_iteration_ini) // initialized as false in bamrunner if iteration_ini >0
    {
      ULONG nmodes_f=this->kvec.size();
      const gsl_rng_type *  T;
      gsl_rng * rng ;
      vector<real_prec> weight(nmodes_f,1.0);
      vector<real_prec> aux_power(nmodes_f,1.0);
      //#pragma omp parallel private(T, r)
      {
        gsl_rng_env_setup();
        ULONG seed=55457;//time(NULL);
        gsl_rng_default_seed=seed;
        T = gsl_rng_mt19937;//.gsl_rng_ranlux;
        rng = gsl_rng_alloc (T);
        // Select kernel according to the previous step. Update the previous kernel with the current step
        real_prec power_ratio=1.0;
        real_prec partial_ratio=1.00;
        real_prec deltak=2.*M_PI/this->params._Lbox();
        real_prec vol= pow(this->params._Lbox(),3);
        real_prec sfac=1./(2.*M_PI*deltak*vol); // this factor is 2 / (4pì delta_k V)
        int counter_mh=0;
        int counter_residuals=0;
        real_prec residuals=0;
        real_prec residuals_unsigned=0;
//#pragma omp parallel for
        for(ULONG i=0;i<nmodes_f;++i)
          {
            real_prec kmode=this->kvec[i];
#ifdef _DO_BiasMT_CALIBRATION_
            real_prec Power_ref=this->Power_REF[i];
            real_prec Power_new=this->Power_NEW[i];
#ifdef _use_more_chis_
            real_prec Power_ref_back= i > 0 ? this->Power_REF[i-1]: 0 ;
            real_prec Power_new_back= i > 0 ? this->Power_NEW[i-1]: 0 ;
            real_prec Power_ref_forw= i< nmodes_f ? this->Power_REF[i+1]:0;
            real_prec Power_new_forw= i< nmodes_f ? this->Power_NEW[i+1]:0;
#endif
#else
            real_prec Power_ref=this->Power_REF_MW[i];
            real_prec Power_new=this->Power_NEW_MW[i];
#endif
            // I make this distiction, for there are approx methods with negative dm power, which lead to "nan" under the square root
            // in the case of the kernel for the DM.
            if(exponent >=0 && exponent<1.0)
              power_ratio=(Power_new == 0 || Power_ref == 0) ? 1.0 : pow(fabs(static_cast<double>(Power_ref)/static_cast<double>(Power_new)), exponent);
            else
              {
                power_ratio = Power_new == 0. ? 1.0 : static_cast<double>(Power_ref)/static_cast<double>(Power_new);   //idelly for exponent = 1.0

                if(i< N_MODES)
                  partial_ratio+=power_ratio;
              }
            if(false==rejection)
              {
#ifndef _use_random_kernel_
                weight[i]=power_ratio;
#endif
                this->power_ratio_unsmoothed[i]=power_ratio;
              }
            else if(true==rejection)
              {
                real_prec kmode_squared=kmode*kmode;  // k²
                real_prec deltaVK_squared=kmode_squared+(1./12.)*deltak*deltak; // This is (1/Delta_K)²
                real_prec inv_deltaVK_squared=sfac/deltaVK_squared;
                // ***********************************
                // These applies if we use instead the likelihhod and take ratios of it.
                // The approach here resembles more a MCMC with  L=exp(-chi**2) for each kbin.
                // The variance at each k-bin is assumed to be Gaussian, and we neglect here shot-noise
//                real_prec sigma_squared =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_new*Power_new); // this is sigma²
#ifndef _DISPLACEMENTS_
                real_prec sigma_squared =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_ref+this->shot_noise_ref)*(Power_ref+this->shot_noise_ref); // this is sigma²
#else
                real_prec sigma_squared =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*Power_ref*Power_ref; // this is sigma²
#endif
#ifdef _use_more_chis_
                real_prec sigma_back_squared = Power_ref_back == 0? 1.0 : inv_deltaVK_squared*(Power_ref_back+this->shot_noise_ref)*(Power_ref_back+this->shot_noise_ref);
                real_prec sigma_forw_squared = Power_ref_forw == 0? 1.0 : inv_deltaVK_squared*(Power_ref_forw+this->shot_noise_ref)*(Power_ref_forw+this->shot_noise_ref);
                real_prec new_H =0.5*(pow(Power_ref- Power_new, 2)/sigma_squared)*weight_central   ;
                new_H+= 0.5*(pow(Power_ref_back- Power_new_back,2)/sigma_back_squared)*weight_back;
                new_H+= 0.5*(pow(Power_ref_forw- Power_new_forw,2)/sigma_forw_squared)*weight_forw;
#else
                // proposed likelihood in the present step
                real_prec sigma_squared_new =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_new+this->shot_noise_ref)*(Power_new+this->shot_noise_ref); // this is sigma²
                real_prec new_H = 0.5*(Power_ref- Power_new)*(Power_ref- Power_new)/sigma_squared_new;
#endif
                // power of the previous step obtained from the ratio ref/new saved in the previous iteration
                real_prec Power_old= static_cast<double>(Power_ref)/static_cast<double>(this->power_ratio_unsmoothed[i]);
                real_prec sigma_squared_old =  Power_ref == 0? 1.0 :  inv_deltaVK_squared*(Power_old+this->shot_noise_ref)*(Power_old+this->shot_noise_ref); // this is sigma²
#ifdef _use_more_chis_
                real_prec Power_old_back= i>0 ? static_cast<double>(Power_ref_back)/static_cast<double>(this->power_ratio_unsmoothed[i-1]): 0 ;
                real_prec Power_old_forw= i<nmodes_f? static_cast<double>(Power_ref_forw)/static_cast<double>(this->power_ratio_unsmoothed[i+1]):0 ;
                real_prec old_H = 0.5*(pow(Power_ref-Power_old,2)/sigma_squared)*weight_central;   // likelihood of the previous step
                old_H+= 0.5*(pow(Power_ref_back-Power_old_back,2)/sigma_back_squared)*weight_back;
                old_H+= 0.5*(pow(Power_ref_forw-Power_old_forw,2)/sigma_forw_squared)*weight_forw;
#else
                real_prec old_H = 0.5*pow(Power_ref-Power_old,2)/sigma_squared_old;   // likelihood of the previous step
#endif
                // likelihood of the previous step
  //                double ratio_diff= exp(-static_cast<double>(new_H/10.0))/exp(-static_cast<double>(old_H/10.0)); // ratios between lilekihoods
                double ratio_diff= static_cast<double>(old_H)/static_cast<double>(new_H);  //ratios between chi²
                // ***********************************  Selection Criteria*******************************************************
                double xran= gsl_rng_uniform (rng);
                if(xran  < min(1.0, ratio_diff))
                  {
                    counter_mh++;
                    weight[i]=power_ratio;
                    aux_power[i]=power_ratio;    //recorded just to print out
                  }
                else
                  {
#ifdef _use_random_kernel_
                    weight[i]=gsl_rng_uniform(rng);  //in the
//#else
//                    weight[i]=1.0; // this is not needed: if not accepted, the weight remains as inizialized above, i.e, =1.
#endif
                    aux_power[i]=this->power_ratio_unsmoothed[i]; //recorded just to print out
                  }
                // *************************************************************************************************************
                this->power_ratio_unsmoothed[i]=power_ratio;
                if(i>INITIAL_MODE_RESIDUALS)
                 {
                   residuals+=fabs(static_cast<real_prec>(power_ratio)-1.0);
                   residuals_unsigned+=static_cast<real_prec>(power_ratio)-1.0;
                   counter_residuals++;
                 }
            }
        }
        So.DONE();
        this->residuals_it=100.0*residuals/static_cast<real_prec>(counter_residuals);
        this->residuals_abs=100.0*residuals_unsigned/static_cast<real_prec>(counter_residuals);
#ifdef _FULL_VERBOSE_
        So.message_screen("Number of modes upgraded for Kernel =", counter_mh);
        So.message_screen("Residuals at this iteration (%) =",100.0*residuals/static_cast<real_prec>(counter_residuals));
        So.message_screen("Average of ratio P(k)_ref / P(k)_new =" , partial_ratio/static_cast<real_prec>(N_MODES));
        So.message_screen("Computed from fundamental mode up to k =", this->kvec[N_MODES]);
        std::cout<<endl;
#else
        std::cout<<BLUE<<"Iteration "<<this->iteration<<RESET<<endl;
        std::cout<<BLUE<<"Residuals = "<<100.0*residuals/static_cast<real_prec>(counter_residuals)<<" %  \r"<<RESET<<endl;
        std::cout<<BLUE<<"Residuals(unsigned) = "<<100.0*residuals_unsigned/static_cast<real_prec>(counter_residuals)<<" %  \r"<<RESET<<endl;
#endif
        this->number_of_modes_upgrade_kernel=counter_mh;
        this->average_ratio_power=partial_ratio/static_cast<real_prec>(N_MODES);
        residuals_power.push_back(100.0*residuals/static_cast<real_prec>(counter_residuals));
        residuals_power_unsigned.push_back(100.0*residuals_unsigned/static_cast<real_prec>(counter_residuals));
        it_power.push_back(this->iteration);
    }
#ifdef _DO_BiasMT_CALIBRATION_
  string kernel_file_or=this->params._Output_directory()+"RatioT_iteration"+to_string(this->iteration)+".txt";
  this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
  this->File.write_to_file(this->params._Output_directory()+"RatioT.txt", this->kvec,aux_power,power_ratio_unsmoothed);
#ifdef _USE_GNUPLOT_
  // PLOT OF THE RESIDUALS
  std::vector<std::pair<double, double> > xy_pts_r;
  std::vector<std::pair<double, double> > xy_pts_ru;
  for(int i=1; i<it_power.size(); ++i)
    xy_pts_r.push_back(std::make_pair(this->it_power[i], this->residuals_power[i]));
  for(int i=1; i<it_power.size(); ++i)
    xy_pts_ru.push_back(std::make_pair(this->it_power[i], this->residuals_power_unsigned[i]));
  this->gp_kernel<<"set size 0.3,0.5\n";
  this->gp_kernel<<"set origin 0.,0.5\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"unset log\n";
  this->gp_kernel<<"set yrange[-1.5:3.5]\n";
  this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
//  this->gp_kernel<<"set title 'Reference: "<<this->params._seed()<<", Iteration "<<this->iteration <<"\n";
  this->gp_kernel<<"set xlabel 'Iteration' font 'Times-Roman,15'\n";
  this->gp_kernel<<"set ylabel 'Residuals %' font 'Times-Roman,15'\n";
  this->gp_kernel<<"plot "<< gp_kernel.file1d(xy_pts_r)<< " w l lw 4 lt 9 title 'Absolute',"<<gp_kernel.file1d(xy_pts_ru)<< "w l lw 4 lt 3 title 'Relative'"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
  xy_pts_ru.clear();
  xy_pts_ru.shrink_to_fit();
  // PLOT OF THE RATIO
  for(int i=1; i<kvec.size(); ++i)
    xy_pts_r.push_back(std::make_pair(this->kvec[i], aux_power[i]));
  this->gp_kernel<<"set size 0.30,0.5\n";
  this->gp_kernel<<"set origin 0.0,0.0\n";
  this->gp_kernel<<"set yrange[*:*]\n";
  this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"unset log y\n";
  this->gp_kernel<<"set log x\n";
  this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ylabel 'Ratio T(k)' font 'Times-Roman,15'\n";
  this->gp_kernel<<"plot " << gp_kernel.file1d(xy_pts_r) << "w l lw 4 lt 6 title 'Ratio to Reference'"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
#endif  // endif for #ifdef _USE_GNUPLOT_
  aux_power.clear(); //to release memory before going out of scope
  aux_power.shrink_to_fit();
#else
  string kernel_file_or=this->params._Output_directory()+"RatioT_mass_assignment_iteration"+to_string(this->iteration)+".txt";
  this->File.write_to_file(kernel_file_or, this->kvec,aux_power,power_ratio_unsmoothed);
  aux_power.clear(); //to release memory before going out of scope
  aux_power.shrink_to_fit();
#ifdef _FULL_VERBOSE_
  So.message_screen("Updating BAM-Kernel");
#endif
#endif // endif for #ifdef _DO_BiasMT_CALIBRATION_
  vector<real_prec> coords(this->params._Nft(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<coords.size() ;++i)
    coords[i]= (i<=this->params._Nft()/2? static_cast<real_prec>(i): -static_cast<real_prec>(this->params._Nft()-i));
  vector<int>nmodes(this->kvec.size(), 0);
  vector<real_prec>kernel_updated(this->kvec.size(), 0);// to print out the shell-averaged kernel
#ifndef _use_random_kernel_
  if(true==rejection)
    {
#endif
#pragma omp parallel for collapse(3)
      for(ULONG i=0;i< this->params._Nft(); ++i)
        for(ULONG j=0;j< this->params._Nft(); ++j)
          for(ULONG k=0;k< this->params._Nft()/2+1; ++k)
            {
              ULONG ind=index_3d(i,j,k, this->params._Nft(), this->params._Nft()/2+1);
              real_prec kv=_get_modulo(coords[i]*this->params._d_deltak_x(),coords[j]*this->params._d_deltak_y(),coords[k]*this->params._d_deltak_z());
              int kmod=static_cast<int>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
              if(kmod<this->kvec.size())
               {
#ifdef _use_random_kernel_
                 this->Kernel[ind]=weight[i]; // Random kernel
#else
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                 this->Kernel[ind]*=weight[kmod]; //weight is 1 if no improvement, so the kernel is not updated;
#endif
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                  kernel_updated[kmod]+=this->Kernel[ind]; //weight is 1 of no improvement; the new kernel if it gets closer
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                  nmodes[kmod]++;
                }
            }
#ifndef _use_random_kernel_
    }
    else
      {
#pragma omp parallel for collapse(3)
        for(ULONG i=0;i< this->params._Nft(); ++i)
          for(ULONG j=0;j< this->params._Nft(); ++j)
            for(ULONG k=0;k< this->params._Nft()/2+1; ++k)
              {
                ULONG ind=index_3d(i,j,k, this->params._Nft(), this->params._Nft()/2+1);
                real_prec kv=_get_modulo(coords[i]*this->params._d_deltak_x(),coords[j]*this->params._d_deltak_y(),coords[k]*this->params._d_deltak_z());
                int kmod=static_cast<int>(floor( (kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
                if(kmod<this->kvec.size())
                  {
                    this->Kernel[ind] = weight[kmod];
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                    nmodes[kmod]++;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                    kernel_updated[kmod]+=this->Kernel[ind]; //weight is 1 of no improvement; the new kernel if it gets closer
                  }
              }
        }
#endif
// Get the shell-averaged version of the Kernel in FOurier space:
  for(ULONG i=0;i< kernel_updated.size(); ++i)
    kernel_updated[i]/=static_cast<real_prec>(nmodes[i]);
  nmodes.clear(); nmodes.shrink_to_fit();
//   Fix the value of the kernel at the fundamental mode by assigning to that mode an average of the first three values of kernel
 // kernel_updated[0]=1;
// if(this->iteration == this->params._N_iterations_Kernel())
  kernel_updated[0]=(kernel_updated[1]+kernel_updated[2]+kernel_updated[3])/3.0;  // USE THIS WHEN NON FIXED AMPLITUD IS USED
#ifdef _DO_BiasMT_CALIBRATION_
  string file_kernel=this->params._Output_directory()+"Kernel_Fourier_iteration"+to_string(this->iteration)+".txt";
  this->File.write_to_file(file_kernel, this->kvec, kernel_updated);
  this->File.write_to_file(this->params._Output_directory()+"Kernel_Fourier.txt", this->kvec, kernel_updated);
#ifdef _USE_GNUPLOT_
  // PLOT OF THE KERNEL
  for(int i=0; i<kvec.size(); ++i)
    xy_pts_r.push_back(std::make_pair(this->kvec[i], kernel_updated[i]));
  xy_pts_ru.clear();
  xy_pts_ru.shrink_to_fit();
  if(this->iteration>2)
      for(int i=0; i<kvec.size(); ++i)
        xy_pts_ru.push_back(std::make_pair(this->kvec[i], 0.5*sqrt(this->Power_NEW[i]/Power_DM_KONV[i])));
  std::vector<std::pair<double, double> > xy_pts_ra;
  for(int i=0; i<kvec.size(); ++i)
    xy_pts_ra.push_back(std::make_pair(this->kvec[i], sqrt(this->Power_REF[i]/Power_DM_REF[i])));
  std::vector<std::pair<double, double> > xy_pts_rb;
  if(this->iteration>2)
      for(int i=0; i<kvec.size(); ++i)
        xy_pts_rb.push_back(std::make_pair(this->kvec[i], 0.5*sqrt(this->Power_DM_REF[i]/Power_DM_KONV[i])));
#endif
  // **********************************************************************************************************************************************************
  // Update averaged version of Kernel
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
#ifdef _USE_GNUPLOT_
  std::vector<std::pair<double, double> > xy_pts_rz;
#endif
  this->average_kernel_updated.resize(this->kvec.size(), 0);// to print out the shell-averaged kernel
  if(this->iteration>=this->params._Iteration_Kernel_average())
    {
      for(int i=0; i<kvec.size(); ++i)
        this->average_kernel_updated[i]+=kernel_updated[i]/static_cast<real_prec>(this->params._N_iterations_Kernel()-this->params._Iteration_Kernel_average());
    }
  if(this->iteration == this->params._N_iterations_Kernel())
   {
      this->average_kernel_updated[0]=(this->average_kernel_updated[1]+this->average_kernel_updated[2]+this->average_kernel_updated[3]+this->average_kernel_updated[4])/4.0;  // USE THIS WHEN NON FIXED AMPLITUD IS USED
      this->File.write_to_file(this->params._Output_directory()+"Averaged_Kernel_Fourier.txt", this->kvec, average_kernel_updated);
    }
#ifdef _USE_GNUPLOT_
  if(this->iteration>=this->params._Iteration_Kernel_average())
    {
      for(int i=1; i<kvec.size(); ++i)
        xy_pts_rz.push_back(std::make_pair(this->kvec[i], this->average_kernel_updated[i]));
    }
#endif // endif for gnuplot
#endif // end for _WRITE_AVERAGE_KERNEL_AND_BIAS_
  // **********************************************************************************************************************************************************
#ifdef _USE_GNUPLOT_
  this->gp_kernel<<"set size 0.3,0.5\n";
  this->gp_kernel<<"set origin 0.3,0.5\n";
  this->gp_kernel<<"set border linecolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"unset log y\n";
  this->gp_kernel<<"set log x\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"'\n";
  this->gp_kernel<<"set ylabel 'Kernel K(k)' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"'\n";
  if(this->iteration>=this->params._Iteration_Kernel_average())
      this->gp_kernel<<"plot" << gp_kernel.file1d(xy_pts_r) << "w l lw 4 lt 0 title 'BAM Kernel',"<<gp_kernel.file1d(xy_pts_rz) << "w l lw 4 lt 3 title 'Average BAM Kernel',"<<gp_kernel.file1d(xy_pts_ru) << "w l lw 4 lt 2 title '(P^{ref}_{H}/P^{K}_{DM})^{1/2} /10' ,"<<gp_kernel.file1d(xy_pts_ra) << "w l lw 4 lt 3 title '(P^{ref}_{H}/P^{ref}_{DM})^{1/2}',"<<gp_kernel.file1d(xy_pts_rb) << "w l lw 4 lt 4 title '(P^{ref}_{DM}/P^{K}_{DM})^{1/2}/10'"<<endl;
  else
      this->gp_kernel<<"plot" << gp_kernel.file1d(xy_pts_r) << "w l lw 4 lt 0 title 'BAM Kernel',"<<gp_kernel.file1d(xy_pts_ru) << "w l lw 4 lt 2 title '(P^{ref}_{H}/P^{K}_{DM})^{1/2} /10' ,"<<gp_kernel.file1d(xy_pts_ra) << "w l lw 4 lt 3 title '(P^{ref}_{H}/P^{ref}_{DM})^{1/2}',"<<gp_kernel.file1d(xy_pts_rb) << "w l lw 4 lt 4 title '(P^{ref}_{DM}/P^{K}_{DM})^{1/2}/10'"<<endl;
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
  xy_pts_ru.clear();
  xy_pts_ru.shrink_to_fit();
  xy_pts_ra.clear();
  xy_pts_ra.shrink_to_fit();
  xy_pts_rz.clear();
  xy_pts_rz.shrink_to_fit();
  // **********************************************************************************************************************************************************
  // This section writes the information on the gnuplot screen
  real_prec counter_gp=-100.0;
  real_prec delta_space=4;
  int counter_ggp=0;
  this->gp_kernel<<"set size 0.25, 0.5 \n";
  this->gp_kernel<<"set origin 0.6,0.5 \n";
  this->gp_kernel<<"set border linecolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"unset grid\n";
  this->gp_kernel<<"unset tics\n";
  this->gp_kernel<<"unset log y\n";
  this->gp_kernel<<"unset log x\n";
  this->gp_kernel<<"set xrange[0:1] \n";
  this->gp_kernel<<"set yrange[-200:-100] \n";  //con este range evito que estos labels salgan en otras subplots
  this->gp_kernel<<"set xlabel '' \n";
  this->gp_kernel<<"set ylabel '' \n";
  counter_gp-=delta_space;
  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'BAM: CALIBRATION PROCESS' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'Times-Bold,12'\n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Reference = "<<this->params._realization()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Lbox = "<<this->params._Lbox()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Nft = "<<this->params._Nft()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Nyquist frequency = "<<this->params._d_nyquist_frequency()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Large-scale bias = "<<this->lss_halo_bias<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Maximum number of tracers in cells = "<<this->nmax_y_onecell<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'MODEL = "<<MODEL_THETA<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'NX = "<<this->params._NX()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Nmk = "<<this->params._n_sknot_massbin()<<"' at 0.1, "<<counter_gp<<"  textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Number of Cwt = "<<this->params._n_cwt()<<"' at 0.1, "<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Iteration "<<this->iteration<<"' at  0.1, "<<counter_gp<<" textcolor rgb 'orange' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  if(this->residuals_it>1)
      this->gp_kernel<<"set label "<< counter_ggp<<" 'Absolute Residuals = "<<floor(1000*this->residuals_it)/1000<<"' at 0.1, "<<counter_gp<<" textcolor rgb 'red' font 'arial,10' \n";
  else
      this->gp_kernel<<"set label "<< counter_ggp<<" 'Absolute Residuals = "<<floor(1000*this->residuals_it)/1000<<"' at 0.1, "<<counter_gp<<" textcolor rgb 'green' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Relative Residuals = "<<floor(1000*this->residuals_abs)/1000<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Knots = "<<this->cwclass._knots_fraction()<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;
  counter_ggp++;  this->gp_kernel<<"set label "<< counter_ggp<<" 'Filaments = "<<this->cwclass._filaments_fraction()<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Sheets = "<<this->cwclass._sheets_fraction()<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<" 'Voids = "<<this->cwclass._voids_fraction()<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<"  'Number of modes updated = "<<this->number_of_modes_upgrade_kernel<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<"  'Average rtr (power) = "<<this->average_ratio_power<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<"  'Number of tracers in reference = "<<this->N_objects_Y<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space;  counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<"  'Number of tracers in new mock = "<<this->Nobjects<<"' at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  counter_gp-=delta_space; counter_ggp++;
  this->gp_kernel<<"set label "<< counter_ggp<<"  'Date: "<<__DATE__<<"  "<<__TIME__<<" '  at 0.1,"<<counter_gp<<" textcolor rgb '"<<FG_COLOR<<"' font 'arial,10' \n";
  this->gp_kernel<<"plot 2.0 notitle "<<endl;
// **********************************************************************************************************************************************************
#endif // end for _USE_GNUPLOT_
#else  // else for ifdef _DO_BiasMT_CALIBRATION_
  string file_kernel=this->params._Output_directory()+"Kernel_Fourier_mass_assignment_iteration"+to_string(this->iteration)+".txt";
  this->File.write_to_file(file_kernel, this->kvec, kernel_updated, weight);
#endif
#ifdef _SMOOTHED_KERNEL_
#ifdef _SMOOTHED_KERNEL_LAST_ITERATION_
  if(this->iteration==this->params._N_iterations_Kernel())
    {
#endif
      lin_smooth(this->kvec, kernel_updated,8);
      lin_smooth(this->kvec, kernel_updated,4);
      lin_smooth(this->kvec, kernel_updated,2);
        //This line sets the kernek in the first mode to an average of the first three modes
      string kernel_file=this->params._Output_directory()+"Kernel_Fourier_smoothed_iteration"+to_string(this->iteration)+".txt";
#ifdef _FULL_VERBOSE_
      So.message_screen("Writting smoothed version of kernel in file", kernel_file);
#endif
      this->File.write_to_file(kernel_file, this->kvec, kernel_updated, weight);
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
      string a_kernel_file=this->params._Output_directory()+"Average_Kernel_Fourier_smoothed_iteration"+to_string(this->iteration)+".txt";
      this->File.write_to_file(a_kernel_file, this->kvec, this->average_kernel_updated, weight);
#endif
#ifdef _SMOOTHED_KERNEL_LAST_ITERATION_
    }
#endif
// ****************************************************************************************************************************
// Assign smoothed shell averaged Kernel to 3D kernel:
// ****************************************************************************************************************************
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
    this->Average_Kernel.resize(this->Kernel.size(),0);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
  for(ULONG i=0;i< this->params._Nft(); ++i)
    for(ULONG j=0;j< this->params._Nft(); ++j)
      for(ULONG k=0;k< this->params._Nft()/2+1; ++k)
        {
          ULONG ind=index_3d(i,j,k, this->params._Nft(), this->params._Nft()/2+1);
          real_prec kv=_get_modulo(coords[i]*this->params._d_deltak_x(),coords[j]*this->params._d_deltak_y(),coords[k]*this->params._d_deltak_z());
          int kmod=static_cast<int>(floor( (kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
          if(kmod<this->params._Nft()/2/this->params._ndel_data())
          {
            this->Kernel[ind]=kernel_updated[kmod];
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
            this->Average_Kernel[ind]=average_kernel_updated[kmod];
#endif
        }
    }
#endif // end of smooth
// ****************************************************************************************************************************
// All other vectors defined here and not released
// are destroyed her when going out of scope
// ****************************************************************************************************************************
#ifndef _TEST_THRESHOLDS_RESIDUALS_
      if(this->iteration == this->params._N_iterations_Kernel())
      {
        this->File.write_array(this->params._Output_directory()+"Bam_Kernel", this->Kernel);
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
        this->File.write_array(this->params._Output_directory()+"Average_BiasMT_Kernel", this->Kernel);
#endif
        }
#endif
  if(this->iteration == this->params._N_iterations_Kernel())
    {
#ifdef DOUBLE_PREC
     complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
     complex_prec *data_out= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif
// ****************************************************************************************************************************
// The original BAM code had only the real part assigned, with IMAG set to zero.
// ****************************************************************************************************************************
#pragma omp parallel for
    for(ULONG ind=0;ind< this->NTT ;++ind)
      {
        data_out[ind][REAL]=this->Kernel[ind];
        data_out[ind][IMAG]=0;
      }
    vector<real_prec>aux( this->params._NGRID(),0);
    do_fftw_c2r(this->params._Nft(), data_out,aux);
    this->File.write_array(this->params._Output_directory()+"Bam_Kernel_config_space", aux);
    aux.clear(); aux.shrink_to_fit();

#ifdef DOUBLE_PREC
    fftw_free(data_out);
#else
    fftwf_free(data_out);
#endif
    }
  this->use_iteration_ini=false;
  }
/*
else
    {
      auto out_it = std::find(std::begin(this->params._output_at_iteration()), std::end(this->params._output_at_iteration()), this->iteration);
      if((this->iteration==this->params._N_iterations_Kernel()) || (out_it != std::end(this->params._output_at_iteration())))
      {
#ifdef _FULL_VERBOSE_
        So.message_screen("Reading Kernel from iteration", this->params._iteration_ini());
#endif
        this->Kernel.clear();
        this->Kernel.resize(this->NTT, 0.0);
        this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat", this->Kernel);
        this->use_iteration_ini=true; // we set it true for we will need it for the bias
      }
   }*/
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::Konvolve(vector<real_prec> &in, vector<real_prec>&out)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _FULL_VERBOSE_
  if(this->iteration==0)
    So.message_screen("Generating new DM density field by convolution of input DM with input Kernel");
  else
    So.message_screen("Generating new DM density field by convolution of DM with updated mass-weighterd kernel ");
#endif
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Generating new DM density field by convolution of input DM with updated Kernel");
#endif
#endif
#ifdef DOUBLE_PREC
  complex_prec *data_out= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
  complex_prec *data_out= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->NTT;++i)
   {
     data_out[i][REAL]=0;
     data_out[i][IMAG]=0;
   }
  do_fftw_r2c(this->params._Nft(),in, data_out);
double paux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:paux)
#endif
  for(ULONG i=0;i< this->params._NGRID();i++)
     paux+=static_cast<double>(in[i]);
  if(true==isinf(paux))
   {
      So.message_warning("Not defined value found in container at function ",__PRETTY_FUNCTION__);
      So.message_warning("Line" ,__LINE__);
      So.message_warning("Code exits here");
      exit(0);
   }
#ifdef _EXTRAPOLATE_VOLUME_
  real_prec correction_factor = 1.0;
#else
  real_prec correction_factor = 1.0;
#endif
  double we=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG ind=0; ind< this->NTT ;++ind)
    {
      real_prec cor = this->Kernel[ind]*correction_factor;
      real_prec reald=data_out[ind][REAL]*cor;
      real_prec imagd=data_out[ind][IMAG]*cor;
      data_out[ind][REAL]=reald;
      data_out[ind][IMAG]=imagd;
      we+=cor;
    }
  // Do INverse Fourier transform:
  do_fftw_c2r(this->params._Nft(), data_out, out);
  // Correct for the normalization of the kernel
  // for the function  do_fftw_c2r returns the transform normalized by the NGRID, so I divide by NGRID and by multiply by 2 we
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< this->params._NGRID();i++)
    {
      real_prec out_s=out[i];
      out[i]=out_s/(static_cast<double>( this->params._NGRID())/static_cast<double>(2.0*we));
    }
  So.DONE();
#ifdef DOUBLE_PREC
  fftw_free(data_out);
#else
  fftwf_free(data_out);
#endif
/*
  if(this->params._iteration_ini()>0 && this->iteration!=this->params._iteration_ini())
    {
      auto out_it = std::find(std::begin(this->params._output_at_iteration()), std::end(this->params._output_at_iteration()), this->iteration);
      if (out_it != std::end(this->params._output_at_iteration()))
        {
          // Write the kernel interpolated to 3D in Fourier space.
#ifdef DOUBLE_PREC
          complex_prec *kern= (complex_prec *)fftw_malloc(2*this->NTT*sizeof(real_prec));
#else
          complex_prec *kern= (complex_prec *)fftwf_malloc(2*this->NTT*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->NTT;++i)
            {
              kern[i][REAL]=this->Kernel[i];
              kern[i][IMAG]=0;
            }

          vector<real_prec>aux( this->params._NGRID(),0);
          do_fftw_c2r(this->params._Nft(),kern,aux);
          // para no escribir lo que ha leído
          string file_kernel=this->params._Output_directory()+"3DKernel_iteration"+to_string(this->iteration);
          this->File.write_array(file_kernel, aux);
          aux.clear();aux.shrink_to_fit();
#ifdef DOUBLE_PREC
          fftw_free(kern);
#else
          fftwf_free(kern);
#endif
          }
   }

*/
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_new_min_max_properties()
{
  this->So.enter(__PRETTY_FUNCTION__);
  // These are the limits in the tracer info, which should remain constant
  // Only computed once in the first step (see below, in case we do not have -COUNTS:)
  // Limits on the reference
  this->s_mins.prop0=this->Ymin;
  this->s_maxs.prop0=this->Ymax;
  this->s_deltas.prop0=this->DELTAY;
  // This is only done once, as it is meant for the reference. THIS AHS TO BE REVISED, SINCE, WHEN WE USE CONTINIOUS VARIABLES, MION AND MAX ARE IN ANY CASE REFERRED TO BINS IN TH HISTOGRAMS
  // IN REALITY, this->Name_Property_Y MUST BE ALWAYS COUNTS, EVEN IF IT IS A CONTINIOUS QUANTITY, WE WILL COBNVERT IT TO A HISTOGRAM
  if(this->iteration==0)
    {
#ifdef _USE_MASS_FIELD_
      this->s_mins.prop0_mass=get_min(this->delta_Y_MASS);
      this->s_maxs.prop0_mass=get_max(this->delta_Y_MASS);
      this->s_deltas.prop0_mass=(this->s_maxs.prop0_mass-this->s_mins.prop0_mass)/static_cast<real_prec>(this->params._NY_MASS());
      So.message_screen("Maximim mass Tracer =", this->s_maxs.prop0_mass);
      So.message_screen("Minimum mass  Tracer =", this->s_mins.prop0_mass);
#endif
    }//close else of if _COUNTS_!=this->Name_Property_Y
#ifdef _MODIFY_LIMITS_
#ifdef _USE_DM_IN_BiasMT_
  // Gert new extremes for delta_X
  this->s_mins.prop1=get_min(this->delta_X);
  this->s_maxs.prop1=get_max(this->delta_X);
  this->s_deltas.prop1=(this->s_maxs.prop1-this->s_mins.prop1)/static_cast<real_prec>(this->new_nbins_x);
#ifdef _FULL_VERBOSE_
  if(this->params._Scale_X()=="linear"){
   So.message_screen("Minimum ð DM =", this->s_mins.prop1);
   So.message_screen("Maximum ð DM =", this->s_maxs.prop1);
  }
else{
  So.message_screen("Minimum log(2+ð) DM = ", this->s_mins.prop1);
  So.message_screen("Maximum log(2+ð) DM = ", this->s_maxs.prop1);
  }
  //  So.message_screen("N_bins  ð DM = ", this->new_nbins_x);
#endif
#endif // end use dm in bam
#ifdef _DISPLACEMENTS_
  this->s_mins.prop0=get_min(this->delta_Y);
  this->s_maxs.prop0=get_max(this->delta_Y);
  this->s_deltas.prop0=(this->s_maxs.prop0-this->s_mins.prop0)/static_cast<real_prec>(this->new_nbins_y);
  So.message_screen("Maximim displacement = ", this->s_maxs.prop0);
  So.message_screen("Minimum displacement =", this->s_mins.prop0);
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
  this->s_mins.prop4=get_min(this->cwclass.Invariant_TF_II);
  this->s_maxs.prop4=get_max(this->cwclass.Invariant_TF_II);
  this->s_deltas.prop4=(this->s_maxs.prop4-this->s_mins.prop4)/static_cast<real_prec>(N_C_BIN1);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximim InvTF2 = ", this->s_maxs.prop4);
  So.message_screen("Minimum InvTF2 = ", this->s_mins.prop4);
  So.message_screen("Delta InvTF2 = ", this->s_deltas.prop4);
#endif
#elif defined _USE_DELTA2_
  this->s_mins.prop4=get_min(this->cwclass.DELTA2);
  this->s_maxs.prop4=get_max(this->cwclass.DELTA2);
  this->s_deltas.prop4=(this->s_maxs.prop4-this->s_mins.prop4)/static_cast<real_prec>(N_C_BIN1);
  So.message_screen("Maximim ð² = ", this->s_maxs.prop4);
  So.message_screen("Minimum ð² = ", this->s_mins.prop4);
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
  this->s_mins.prop5=get_min(this->cwclass.Invariant_TF_III);
  this->s_maxs.prop5=get_max(this->cwclass.Invariant_TF_III);
  this->s_deltas.prop5=(this->s_maxs.prop5-this->s_mins.prop5)/static_cast<real_prec>(N_C_BIN2);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximim InvTF3 = ", this->s_maxs.prop5);
  So.message_screen("Minimum InvTF3 = ", this->s_mins.prop5);
  So.message_screen("Delta InvTF3 = ", this->s_deltas.prop5);
#endif
#elif defined _USE_DELTA3_
  this->s_mins.prop5=get_min(this->cwclass.DELTA3);
  this->s_maxs.prop5=get_max(this->cwclass.DELTA3);
  this->s_deltas.prop5=(this->s_maxs.prop5-this->s_mins.prop5)/static_cast<real_prec>(N_C_BIN2);
  So.message_screen("Maximim ð³ = ", this->s_maxs.prop5);
  So.message_screen("Minimum ð³ = ", this->s_mins.prop5);
#endif
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
  this->s_mins.prop6=get_min(this->cwclass.Invariant_TF_IV);
  this->s_maxs.prop6=get_max(this->cwclass.Invariant_TF_IV);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim InvTFIII = ", this->s_maxs.prop6);
  So.message_screen("Minimum InvTFII = ", this->s_mins.prop6);
  So.message_screen("Delta InvTFII = ", this->s_deltas.prop6);
#elif defined _USE_TIDAL_ANISOTROPY_
  this->s_mins.prop6=get_min(this->cwclass.Tidal_Anisotropy);
  this->s_maxs.prop6=get_max(this->cwclass.Tidal_Anisotropy);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim Tidal Anisotropy = ", this->s_maxs.prop6);
  So.message_screen("Minimum Tidal Anisotropy = ", this->s_mins.prop6);
#elif defined _USE_S2_
  this->s_mins.prop6=get_min(this->cwclass.S2);
  this->s_maxs.prop6=get_max(this->cwclass.S2);
  this->s_deltas.prop6=(this->s_maxs.prop6-this->s_mins.prop6)/static_cast<real_prec>(N_C_BIN3);
  So.message_screen("Maximim s² = ", this->s_maxs.prop6);
  So.message_screen("Minimum s² = ", this->s_mins.prop6);
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
  this->s_mins.prop7=get_min(this->cwclass.Invariant_VS_I);
  this->s_maxs.prop7=get_max(this->cwclass.Invariant_VS_I);
  this->s_deltas.prop7=(this->s_maxs.prop7-this->s_mins.prop7)/static_cast<real_prec>(N_CV_BIN1);
  So.message_screen("Maximim InvVS1 = ", this->s_maxs.prop7);
  So.message_screen("Minimum InvVS1 = ", this->s_mins.prop7);
#elif defined _USE_NABLA2DELTA_
  this->s_mins.prop7=get_min(this->cwclass.N2D);
  this->s_maxs.prop7=get_max(this->cwclass.N2D);
  this->s_deltas.prop7=(this->s_maxs.prop7-this->s_mins.prop7)/static_cast<real_prec>(N_CV_BIN1);
  So.message_screen("Maximim Nabla²ð = ", this->s_maxs.prop7);
  So.message_screen("Minimum Nabla²ð = ", this->s_mins.prop7);
#elif defined _USE_INVARIANT_PWEB_I_
  this->s_mins.prop7=get_min(this->cwclass.Invariant_TF_I);
  this->s_maxs.prop7=get_max(this->cwclass.Invariant_TF_I);
  this->s_deltas.prop7=(this->s_maxs.prop7-this->s_mins.prop7)/static_cast<real_prec>(N_CV_BIN1);
  So.message_screen("Maximim PwebI1 = ", this->s_maxs.prop7);
  So.message_screen("Minimum PwebI1= ", this->s_mins.prop7);
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
  this->s_mins.prop8=get_min(this->cwclass.Invariant_VS_II);
  this->s_maxs.prop8=get_max(this->cwclass.Invariant_VS_II);
  this->s_deltas.prop8=(this->s_maxs.prop8-this->s_mins.prop8)/static_cast<real_prec>(N_CV_BIN2);
  So.message_screen("Maximim InvVSII = ", this->s_maxs.prop8);
  So.message_screen("Minimum InvVSII = ", this->s_mins.prop8);
#elif defined _USE_S2DELTA
  this->s_mins.prop8=get_min(this->cwclass.S2DELTA);
  this->s_maxs.prop8=get_max(this->cwclass.S2DELTA);
  this->s_deltas.prop8=(this->s_maxs.prop8-this->s_mins.prop8)/static_cast<real_prec>(N_CV_BIN2);
  So.message_screen("Maximim s²ð = ", this->s_maxs.prop8);
  So.message_screen("Minimum s²ð = ", this->s_mins.prop8);
#elif defined _USE_INVARIANT_PWEB_II_
  this->s_mins.prop8=get_min(this->cwclass.Invariant_TF_II);
  this->s_maxs.prop8=get_max(this->cwclass.Invariant_TF_II);
  this->s_deltas.prop78(this->s_maxs.prop8-this->s_mins.prop8)/static_cast<real_prec>(N_CV_BIN2);
  So.message_screen("Maximim PwebI2 = ", this->s_maxs.prop8);
  So.message_screen("Minimum PwebI2= ", this->s_mins.prop8);
#endif
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
  this->s_mins.prop9=get_min(this->cwclass.Invariant_VS_III);  //not yet assigned to a container
  this->s_maxs.prop9=get_max(this->cwclass.Invariant_VS_III);
  this->s_deltas.prop9=(this->s_maxs.prop9-this->s_mins.prop9)/static_cast<real_prec>(N_CV_BIN3);
  So.message_screen("Maximim InvVSIII = ", this->s_maxs.prop9);
  So.message_screen("Minimum InvVSIII = ", this->s_mins.prop9);
#elif defined _USE_S3_
  this->s_mins.prop9=get_min(this->cwclass.S3);
  this->s_maxs.prop9=get_max(this->cwclass.S3);
  this->s_deltas.prop9=(this->s_maxs.prop9-this->s_mins.prop9)/static_cast<real_prec>(N_CV_BIN3);
  So.message_screen("Maximim s³ = ", this->s_maxs.prop9);
  So.message_screen("Minimum s³ = ", this->s_mins.prop9);
#elif defined _USE_INVARIANT_PWEB_III_
  this->s_mins.prop9=get_min(this->cwclass.Invariant_TF_III);
  this->s_maxs.prop9=get_max(this->cwclass.Invariant_TF_III);
  this->s_deltas.prop9=(this->s_maxs.prop8-this->s_mins.prop9)/static_cast<real_prec>(N_CV_BIN3);
  So.message_screen("Maximim PwebI3 = ", this->s_maxs.prop9);
  So.message_screen("Minimum PwebI3= ", this->s_mins.prop9);
#endif
#else    // if _MODIFY_LIMITS_is not defined,
  real_prec xmin_temp=get_min(this->delta_X);
  real_prec xmax_temp=get_max(this->delta_X);
#ifdef _FULL_VERBOSE_
  if(this->Xmax< xmax_temp)
  {
    So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Xmax in delta below the nominal value");
    std::cout<<"this->Xmax="<<this->Xmax<<"  Current="<<xmax_temp<<endl;
  }
  if(this->Xmin> xmin_temp)
    So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Xmin in delta_X above the nominal value");
#endif
  this->s_mins.prop1=this->Xmin;
  this->s_mins.prop4=C1_MIN;
  this->s_mins.prop5=C2_MIN;
  this->s_mins.prop6=C3_MIN;
  this->s_mins.prop7=CV1_MIN;
  this->s_mins.prop8=CV2_MIN;
  this->s_mins.prop9=CV3_MIN;
  this->s_maxs.prop1=this->Xmax;
  this->s_maxs.prop4=C1_MAX;
  this->s_maxs.prop5=C2_MAX;
  this->s_maxs.prop6=C3_MAX;
  this->s_maxs.prop7=CV1_MAX;
  this->s_maxs.prop8=CV2_MAX;
  this->s_maxs.prop9=CV3_MAX;
  this->s_deltas.prop1=this->DELTAX;
  this->s_deltas.prop4=DELTA_C1;
  this->s_deltas.prop5=DELTA_C2;
  this->s_deltas.prop6=DELTA_C3;
  this->s_deltas.prop7=DELTA_CV1;
  this->s_deltas.prop8=DELTA_CV2;
  this->s_deltas.prop9=DELTA_CV3;
#endif   // end of ifdef _MODIFY_LIMITS_
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void  BiasMT::get_min_max_X_Y()
{
  // This function redefines the values of delta_X_min etc in case. If this is not called, the input values are used
  // The filds delta X and delta Y ahave been already, if requiested, converted to overdensities
   this->So.enter(__PRETTY_FUNCTION__);
   if(true==this->params._Redefine_limits())
    {
      real_prec num_in_log_x = true==this->params._Convert_Density_to_Delta_X() ? NUM_IN_LOG: 0.;
      real_prec num_in_log_y = true==this->params._Convert_Density_to_Delta_Y() ? NUM_IN_LOG: 0.;
      vector<real_prec>AUX(this->delta_X.size());
      if(this->params._Scale_X()=="linear")
        {
          this->Xmin=get_min<real_prec>(this->delta_X);
          this->Xmax=get_max<real_prec>(this->delta_X);
        }
      else
        if(this->params._Scale_X()=="log")
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_X.size();++i)AUX[i]=log10(num_in_log_x+this->delta_X[i]);
          this->Xmin=get_min<real_prec>(AUX);
          this->Xmax=get_max<real_prec>(AUX);
        }
      if(this->params._Scale_Y()=="linear")
        {
          this->Ymin=get_min<real_prec>(this->delta_Y);
          this->Ymax=get_max<real_prec>(this->delta_Y);
        }
      else
        if(this->params._Scale_Y()=="log")
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_Y.size();++i)
            AUX[i]=log10(num_in_log_y+this->delta_Y[i]);
          this->Ymin=get_min<real_prec>(AUX);
          this->Ymax=get_max<real_prec>(AUX);
        }
      AUX.clear();  // to release memory before going out of scope
      AUX.shrink_to_fit();
    }
  else
    {
      if(this->params._Scale_X()=="log")
        {
          this->Xmin=this->params._ldelta_X_min();
          this->Xmax=this->params._ldelta_X_max();
        }
      else
        {
          this->Xmin=this->params._delta_X_min();
          this->Xmax=this->params._delta_X_max();
        }
      if(this->params._Scale_Y()=="log")
        {
          this->Ymin=this->params._ldelta_Y_min();
          this->Ymax=this->params._ldelta_Y_max();
        }
      else
        {
          this->Ymin=this->params._delta_Y_min();
          this->Ymax=this->params._delta_Y_max();
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_VELOCITIES_
void BiasMT::read_BiasMT_files(string file_X, string file_Y, string file_Y_HR, string file_Y_mass, string file_Y_sat_frac, string file_X_ref, string file_Vx, string file_Vy, string file_Vz)
#else
  void BiasMT::read_BiasMT_files(string file_X, string file_Y, string file_Y_HR, string file_Y_mass, string file_Y_sat_frac, string file_X_ref)
#endif
{
   this->So.enter(__PRETTY_FUNCTION__);


   int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);

#ifdef _FULL_VERBOSE_
  std::cout<<endl;
  this->So.message_screen("Reading input *reference* files");
  std::cout<<endl;
#endif
  ULONG NGRID_NEW;
  ULONG NGRID_NEW_HR;
#ifdef _EXTRAPOLATE_VOLUME_
  NGRID_NEW=static_cast<ULONG>(this->params._Nft_low()*this->params._Nft_low()*this->params._Nft_low());
#else
#ifdef _USE_TRACER_HR_
  NGRID_NEW= this->params._NGRID();
  NGRID_NEW_HR=static_cast<ULONG>(this->params._Nft_HR()*this->params._Nft_HR()*this->params._Nft_HR());
#else
  NGRID_NEW= this->params._NGRID();
  NGRID_NEW_HR= this->params._NGRID();
#endif
#endif
#ifndef _EXTRAPOLATE_VOLUME_
  this->delta_X.resize(NGRID_NEW,0);
  this->delta_X_ini.resize(NGRID_NEW,0);
  this->File.read_array(file_X,this->delta_X);
  if(this->params._Quantity()=="DELTA")
    for(ULONG i=0;i<this->params._NGRID();++i)
    this->delta_X[i]*=pow(256,3)/(this->params._Lbox(),3)*(1+this->delta_X[i]);
#ifdef _DISPLACEMENTS_
  for(ULONG i=0;i <  this->params._NGRID(); ++i)
    this->delta_X[i]*=FACTOR_IC;
#endif
#ifdef _CONVERT_CIC_TO_NGP_
  convert_cic_to_ngp(delta_X,delta_X);
#endif
#endif
#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.resize(NGRID_NEW_HR,0);
  this->File.read_array_t<PrecType_Y>(file_Y_HR, this->delta_Y_HR);
#endif
  //  this->delta_Y_new.resize( this->params._NGRID(),0); //Isn't this resized in get_mock?
  this->delta_Y.resize( this->params._NGRID(),0);
  this->File.read_array_t<PrecType_Y>(file_Y, this->delta_Y);
#ifdef _SHOW_EMPTY_CELLS_
  ULONG empty_cells=0;
  for(ULONG i=0;i <  this->params._NGRID(); ++i)
    if(this->delta_Y[i]==0)
      empty_cells++;
  So.message_screen("Number of empty cells from number density field =", empty_cells);
#endif

#ifdef _USE_MASS_TRACERS_
#ifdef _USE_MASS_FIELD_
  this->delta_Y_MASS.resize( this->params._NGRID(),0);
  this->File.read_array_t<PrecType_Y>(file_Y_mass, this->delta_Y_MASS);
#endif
#ifdef _USE_LOG_MASS_
  this->Mean_density_Y_MASS=1;
#else
  /*
    this->Mean_density_Y_MASS=get_mean(this->delta_Y_MASS);
    get_overdens(this->delta_Y_MASS,this->Mean_density_Y_MASS, this->delta_Y_MASS);
  */
#ifdef _USE_MASS_FIELD_
#ifdef _SHOW_EMPTY_CELLS_
  empty_cells=0;
  for(ULONG i=0;i <  this->params._NGRID(); ++i)
    {
      if(this->delta_Y_MASS[i]<1.0)
        empty_cells++;
      this->delta_Y_MASS[i]=log10(NUM_IN_LOG+ this->delta_Y_MASS[i]/MASS_SCALE);
    }
  So.message_screen("Number of empty cells from mass density field =", empty_cells);
#endif
  So.message_screen("Using log(2+total central mass in cells). Check line ",__LINE__);
#endif
#endif
#endif
#ifdef _NGP2CIC_Y_
  convert_ngp_to_cic(this->delta_Y, this->delta_Y);
#endif
  // this->So.message_screen("Minimum Number of Y", get_min(delta_Y));
  // this->So.message_screen("Maximum Number of Y", get_max(delta_Y));
  // VELOCITIES: TAKEN FROM LPT OR READ
#ifdef _USE_VELOCITIES_
  this->Velx_X.resize( this->params._NGRID(),0);
  this->Vely_X.resize( this->params._NGRID(),0);
  this->Velz_X.resize( this->params._NGRID(),0);
  this->File.read_array(file_Vx, this->Velx_X);
  this->File.read_array(file_Vy, this->Vely_X);
  this->File.read_array(file_Vz, this->Velz_X);
#endif

#ifdef _EXCHANGE_X_Y_DENSITY_
  So.message("Exchanging axis X-Y in density field");
  exchange_xy(this->params._Nft(),this->delta_X,this->delta_X);
#endif

#ifdef _USE_VELOCITIES_
#ifdef _EXCHANGE_X_Y_VEL_
  So.message_screen("Exchanging axis X-Y in velocity field");
  exchange_xy(this->params._Nft(),this->Velx_X,this->Velx_X);
  exchange_xy(this->params._Nft(),this->Vely_X,this->Vely_X);
  exchange_xy(this->params._Nft(),this->Velz_X,this->Velz_X);
  So.DONE();
#endif
#endif
#ifdef _KONV_
  vector<real_prec>kernel(this->NTT,0);
  this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat",kernel);
  convolvek(this->params._Nft(),this->delta_X, kernel,this->delta_X);
#endif
  // Initialize the DM density field with the input (REF) density field
  this->delta_X_ini=this->delta_X;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::analyze_input()
{
  So.enter(__PRETTY_FUNCTION__);
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#ifdef _FULL_VERBOSE_
#ifndef _DISPLACEMENTS_
    So.message_screen("Getting some statistics from input files: ");
    std::cout<<endl;
#endif
#endif
// Analyze input file for X variables:
  real_prec nmean_X=0.;
#ifndef _EXTRAPOLATE_VOLUME_
  nmean_X=get_nobjects(this->delta_X);
  this->N_objects_X=nmean_X;
#ifdef _FULL_VERBOSE_
#ifndef _DISPLACEMENTS_
  if(this->params._Name_Property_X()==_COUNTS_)
    So.message_screen("Total number of X objects =", nmean_X);
#endif
#endif
  nmean_X=static_cast<real_prec>(nmean_X)/static_cast<real_prec>( this->params._NGRID());
#ifdef _FULL_VERBOSE_
  if(_COUNTS_==this->params._Name_Property_X())
    So.message_screen("Mean number of X objects in cells =", nmean_X);
  else
    So.message_screen("Mean property X in cells =", nmean_X);
#endif

#ifdef _FULL_VERBOSE_
  if(_COUNTS_==this->params._Name_Property_X())
    So.message_screen("Mean number X density =", nmean_X*static_cast<real_prec>( this->params._NGRID())/pow(this->params._Lbox(),3), "(Mpc / h )⁻³");
  else
    So.message_screen("Mean property X density =", nmean_X*static_cast<real_prec>( this->params._NGRID())/pow(this->params._Lbox(),3));
#endif
  if(this->params._iMAS_X()==0) // if DM commes in NGP, we can make a pdf
    {
      int NPART=static_cast<int>(get_max<real_prec>(this->delta_X)); // Number of particles in cells
      this->PDF_NC_X.resize(NPART, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i< this->delta_X.size();++i)
        if(static_cast<int>(this->delta_X[i])<NPART)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          this->PDF_NC_X[static_cast<int>(this->delta_X[i])]++;

      if(true==this->params._Write_PDF_number_counts())
        {
          string fileX=this->params._Output_directory()+"PDF_NC"+"_X_"+this->params._XNAME()+"_"+this->params._Name_Property_X()+"_MASX"+to_string(this->params._iMAS_X())+"_Nft"+to_string(this->params._Nft())+"_SmoothingScale"+to_string(this->params._smscale())+"_z"+to_string(this->params._redshift())+"_LambdaTH"+to_string(this->params._lambdath())+"_CW"+to_string(this->tstruct)+"_CWclass.txt";
          this->File.write_to_file_i(fileX,this->PDF_NC_X);
        }
      this->nmax_x_onecell=NPART;
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum number of DM particles in one cell =", this->nmax_x_onecell);
      So.message_screen("Estimated Poisson signal-to-noise ratio from DM=", sqrt(this->nmax_x_onecell));
#endif
    }
#ifdef MOCK_MODE
#ifdef _GET_POWER_REFS_
  this->get_power_spectrum("DM_REF");  //gets power from delta_X
#endif
#endif
#endif  // END EXTRAPOLATE
#ifdef _USE_X_COMPLEMENT_I_
#endif
// *********************************************************************************************
// Analyze input file for Y variables
  real_prec nmean_Y=0;
  nmean_Y=get_nobjects(this->delta_Y);
  if(_COUNTS_==this->params._Name_Property_Y())
    {
      this->N_objects_Y=nmean_Y;
#ifdef _EXTRAPOLATE_VOLUME_
      So.message_screen("Note that these are properties of the reference simulation (i.e, smaller volume)");
#endif
#ifdef _FULL_VERBOSE_
      std::cout<<endl;
      So.message_screen("Total number of Y objects =", nmean_Y);
#endif
    }
  real_prec LLBOX;
  ULONG NNGRID;
#ifdef _EXTRAPOLATE_VOLUME_
  LLBOX=this->params._Lbox_low();
  NNGRID= static_cast<ULONG>(this->params._Nft_low()*this->params._Nft_low()*this->params._Nft_low());
#else
  LLBOX=this->params._Lbox();
  NNGRID=  this->params._NGRID();
#endif
  nmean_Y=static_cast<real_prec>(nmean_Y)/static_cast<real_prec>(NNGRID);
#ifdef _FULL_VERBOSE_
  if(_COUNTS_==this->params._Name_Property_Y())
    So.message_screen("Mean number of Y objects in cells =", nmean_Y);
  else
    So.message_screen("Mean property Y in cells =", nmean_Y);
  this->Mean_density_Y=nmean_Y;
  if(_COUNTS_==this->params._Name_Property_Y())
    So.message_screen("Mean Y number density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3), "(Mpc / h )⁻³");
  else
    {
      So.message_screen("Mean Y property-density =", nmean_Y*static_cast<real_prec>(NNGRID)/pow(LLBOX,3));
    }
#endif
  this->new_Name_Property_Y=this->params._Name_Property_Y();
  if(_COUNTS_==this->params._Name_Property_Y() && ZERO==this->params._iMAS_Y())
    {
      this->nmax_y_onecell=static_cast<int>(get_max<real_prec>(this->delta_Y)); // Maximum number of particles in cells
      // Here we get the pdf of the counts in order to measure the mean occupation number
      // thee value is to be allocated in the variable this->mean_number_x_onecell;
#ifdef _UNITSIM_
#ifdef _USE_TWO_REFS_MOCKS_
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Checking occupation number in paired reference");
#endif
      vector<real_prec>auxv(this->params._NGRID(),0);
      this->File.read_array(this->params._Input_Directory_Y_TWO()+this->params._Name_Catalog_Y(),auxv);
      int aux_nmax=static_cast<int>(get_max<real_prec>(auxv)); // Maximum number of particles in cells
      int mmax_cell=max(aux_nmax,this->nmax_y_onecell);
      this->So.message_screen("Maximum number of Y particles in one cell (paired) =",  mmax_cell);
      this->nmax_y_onecell=mmax_cell;
      auxv.clear(); auxv.shrink_to_fit();
#endif
#endif
      this->mean_number_y_onecell=static_cast<int>(get_mean_occupation_number<real_prec>(this->nmax_y_onecell,  this->delta_Y));
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Maximum number of Y particles in one cell =", this->nmax_y_onecell);
      this->So.message_screen("Average number of Y particles in one cell =", this->mean_number_y_onecell );
      this->So.message_screen("Estimated Maximum Poisson signal-to-noise ratio =", sqrt(this->nmax_y_onecell));
#endif
     // If we are using two refs to learn the bias from, as we do for the unitsim, here we choose the max occupation number as the max of those two.
#ifdef _USE_MASS_FIELD_
        this->So.message_screen("Maximum tracer mass in one cell =", get_max<real_prec>(this->delta_Y_MASS));
#endif
      ofstream nxo;
      nxo.open(file_one_cell);
      So.message_screen("Writting max number of objects in cell in file", file_one_cell);
      nxo<<this->nmax_y_onecell<<"\t"<<this->mean_number_y_onecell<<endl;
      nxo.close();
       if(true==this->params._Write_PDF_number_counts())
          {

            this->PDF_NC_Y.resize(this->nmax_y_onecell, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(ULONG i=0;i<NNGRID;++i)
              {
                int nob=static_cast<int>(this->delta_Y[i]);

#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                this->PDF_NC_Y[nob]++;
              }
            string fileY=this->params._Output_directory()+"PDF_NC_Y_REF_"+this->new_Name_Property_Y+"_MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
            this->File.write_to_file_i(fileY,this->PDF_NC_Y);
            this->PDF_NC_Y.clear();this->PDF_NC_Y.shrink_to_fit();
           }
    }
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
#endif
#ifdef MOCK_MODE
#ifndef _EXTRAPOLATE_VOLUME_
#ifdef _GET_POWER_REFS_
  this->get_power_spectrum("TR_REF");
#endif
#ifdef _USE_TRACER_HR_
  this->delta_Y_HR.clear();
  this->delta_Y_HR.shrink_to_fit();
#endif
#endif
#endif
#ifdef _GET_POWER_REFS_
//  this->get_power_spectrum("CROSS_TR_DM"); // NOT WORKING, MEMCHINK
#endif
#ifdef MOCK_MODE
#ifndef _GET_BiasMT_REALIZATIONS_
  real_prec lss_bias=0;
  int ncount=0;
  if(this->Power_DM_REF.size()==0 || this->Power_REF.size()==0)
    {
      So.message_warning("Container Power_DM_REF is not resized. Check preprocessor directive _GET_POWER_REFS_. BAM ends here");
      exit(1);
    }
#ifndef _DISPLACEMENTS_ // The following lines are not requested when BAM is applied to displacement
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:lss_bias, ncount)
#endif
  for(int i=0;i<N_MODES;++i)
    if(this->Power_DM_REF[i]!=0)
      {
        lss_bias += this->Power_REF[i]/this->Power_DM_REF[i];
               ncount++;
      }
  lss_bias/=static_cast<real_prec>(ncount);
  this->lss_halo_bias=sqrt(lss_bias);
#ifdef _FULL_VERBOSE_
  So.message_screen("Large-Scale Bias: P(k)_tracer / P(k)_dm  = ", sqrt(lss_bias));
  std::cout<<endl;
#endif
#endif
#ifdef _UNDER_BIASED_
  if(lss_bias<1)
    {
      this->used_once=false;
      So.message_screen("Treating under-biased tracers");
      std::cout<<endl;
    }
#endif
#endif
#endif
#ifdef _RUN_TEST_
  dark_matter_to_halos_analytical(this->delta_X, this->delta_Y);
  exit(0);
#endif
  // Assign the density field delta_X_ini to *density*  delta_X
  // Below, if delta_X is to be converted to overdensity, we reassign delta_X_ini to overdensity
  // BAM starts with the outcome of this loop. This line is important.
  this->delta_X_ini=this->delta_X;
#ifdef _USE_BINARY_MASK_
  vector<real_prec> binary_mask( this->params._NGRID(),0);
  this->File.read_array(this->params._Name_binary_mask(), binary_mask);
#endif
  // Convert density to delta for both the reference and the approximated field
  if(true==this->params._Convert_Density_to_Delta_X())
    {
      this->new_Name_Property_X="DELTA";
#ifdef _ADD_NOISE_
      // With this we add Poisson noise to the DM density field.
      // The noise is added by increasing the mean number density from that of the original filed
      // to a new valie MEAN_NEW, from wihch we construct a new density field
      // MEAN_NEW*1+delta) and use this value as mean to draw Poisson distributed values.
      // Note that, when measuring power spectrum, if we subtract the Poisson shot-noise from this,
      // we obtain the same input field.
      // Hence if we do that operation (subtract Poiss sn) we better use another distribution to add this
      // variance.
      gsl_rng_env_setup();
      gsl_rng_default_seed=75;
      const gsl_rng_type *  T= gsl_rng_ranlux;
      gsl_rng * r = gsl_rng_alloc (T);
      get_overdens(this->delta_X,  nmean_X, this->delta_X);  // get overdensity
      So.message_screen("Converting DENSITY -> Poisson(DENSITY):");
      gsl_real proba=0.8;
      real_prec expo=0.98;
        real_prec mean_new=0;
      for(ULONG i=0;i<delta_X.size();++i)
            mean_new+=pow(1.+delta_X[i],expo);
      mean_new/=static_cast<real_prec>(delta_X.size());

      for(ULONG i=0;i<delta_X.size();++i)
//         delta_X[i]=gsl_ran_poisson(r,MEAN_NEW*(1.+delta_X[i]));   //recover density
         delta_X[i]=gsl_ran_negative_binomial(r,proba, (MEAN_NEW/mean_new)*pow(1.+delta_X[i],expo));   //recover density
    //    delta_X[i]=(nmean_X/mean_new)*pow(1.+delta_X[i],1.8);   //recover density
      So.DONE();
       {
        this->params.set_input_type("density_grid");
        this->params.set_Name_survey("DM_Poisson");
        this->params.set_SN_correction(true);
        this->params.set_mass_assignment_scheme("CIC");
        this->params.set_MAS_correction(true);
        PowerSpectrumF dPSFa(this->params);
        dPSFa.compute_power_spectrum_grid(this->delta_X,true); //ojo, el argumento debe ser DENSIDAD
        dPSFa.write_power_and_modes();
      }
#endif
#ifdef _FULL_VERBOSE_
      So.message_screen("Converting DENSITY into DELTA for X:");
#endif
#ifndef _USE_BINARY_MASK_
#ifdef _ADD_NOISE_
      get_overdens(this->delta_X,  this->delta_X);
#else
      get_overdens(this->delta_X, nmean_X, this->delta_X);
#ifdef _GIVE_POWER_
      So.message_screen("Adding power on small scales to overdensity");
      give_power(this->params._Lbox(), this->params._Nft() , this->delta_X,this->delta_X);
      So.DONE();
      this->params.set_input_type("delta_grid");
      this->params.set_Name_survey("DM_extra");
      this->params.set_SN_correction(false);
      this->params.set_mass_assignment_scheme("CIC");
      this->params.set_MAS_correction(true);
      PowerSpectrumF dPSFa(this->params);
      dPSFa.compute_power_spectrum_grid(this->delta_X,true,false); //ojo, el argumento debe ser DENSIDAD
      dPSFa.write_power_and_modes();
#endif
#endif
#else
      get_overdens(this->delta_X, binary_mask, this->delta_X,true);
#endif
      // Reassign the delta_X_ini to *overdensities*  delta_X
      this->delta_X_ini=this->delta_X;
#ifdef _USE_BINARY_MASK_
      So.message_screen("Weighting Delta X with selection");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->delta_X.size();++i)
        delta_X[i]*=binary_mask[i];
      this->So.DONE();
#endif
    }
  else
    {
      this->new_Name_Property_X=this->params._Name_Property_X();
#ifdef _USE_BINARY_MASK_
      if(this->params._Name_Property_X()=="DELTA")
        {
          So.message_screen("Weighting Delta X with selection");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->delta_X.size();++i)
            delta_X[i]*=binary_mask[i];
          this->So.DONE();
        }
#endif
    }
// Convert the TRACER to delta
    if(true==this->params._Convert_Density_to_Delta_Y())
      {
#ifdef _FULL_VERBOSE_
        So.message_screen("Converting DENSITY into DELTA for Y:");
#endif
        this->new_Name_Property_Y="DELTA";

#ifdef _USE_BINARY_MASK_
        get_overdens(this->delta_Y, binary_mask,this->delta_Y, false);
        binary_mask.clear();
        binary_mask.shrink_to_fit();
#else
        get_overdens(this->delta_Y, this->delta_Y);
#endif
      }
  else // IF NET, THE NEW NAME REMAINS THE SAME
    this->new_Name_Property_Y=this->params._Name_Property_Y();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function converts the DM DELTA fields to log 2+delta and compute the CWC if requested
// The BAM kernel is initialized
void BiasMT::get_BiasMT_DM()
{
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
  real_prec num_in_log_x= NUM_IN_LOG;
  this->get_min_max_X_Y();// THIS CAN BE REPLACED by finding min and max at converting time
  vector<real_prec>xbaux(this->params._NX(), 0);
  vector<real_prec>pdf_in(this->params._NX(), 0); //contains the pdf of DM in each iteration, before updating it with teh convolution
  this->pdf_ref.resize(this->params._NX(), 0);
  this->pdf_ini.resize(this->params._NX(), 0);
#if defined _WRITE_PDF_ || defined (_USE_GNUPLOT_)
  for(int i=0;i<this->params._NX(); ++i)
    xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(this->params._NX()));
  if(this->iteration==0)
    {
#ifdef _FULL_VERBOSE_
          So.message_screen("Computing PDF from log(1+delta)  dark matter");
#endif
      calc_pdf("log",  this->params._NGRID(),this->params._NX(), this->Xmax, this->Xmin, this->delta_X, this->pdf_ref);
      string filex=this->params._Output_directory()+"pdf_X_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
      this->File.write_to_file(filex, xbaux,pdf_ref);
    }
#endif
  // *********************************************************************************************
  // At this point, we had already said delta_X_ini =delta_x WITHOUT transforming to logs
  // (we convert delta_X to log10(num_in_log+delta_X) in the next function get_BIAS)
#ifdef _UNDER_BIASED_
  if(false==this->used_once)
    {
      So.message_screen("Using TR as DM for b<1 tests. Trick to treat under-biased populations");
      So.message_screen("Test in line", __LINE__, ", function ",__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i = 0;i <  this->params._NGRID() ;++i)
        this->delta_X_ini[i] = (this->delta_X[i]<-1 ?  0 :  log10(num_in_log_x+ static_cast<real_prec>(this->delta_Y[i])));
      this->used_once=true;
    }
#endif
  // *********************************************************************************************
  // Do the CWC class for the first iteration, or in the last , when te target field is loaded
  if(this->iteration==0)
    {
#ifdef _RANK_ORDERING_AB_INITIO_
#ifdef _FULL_VERBOSE_
          So.message_screen("Ready to perform rank ordering of DM to target distribution");
          So.message_screen("From field at", this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X());
          So.message_screen("to pdf from field at", this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X_REF_PDF());
#endif
          this->delta_X_REF_PDF.resize( this->params._NGRID(),0);
          string file_X_ref_pdf=this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X_REF_PDF();
          this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_
          get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF);
#endif
          pdf_in=this->pdf_ref;  // This was already computed above, so we do not replicate it here
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of log(1+delta) DM target reference");
#endif
          calc_pdf("log",  this->params._NGRID(),this->params._NX(), this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
          So.DONE();
#ifdef _WRITE_PDF_
          string filex=this->params._Output_directory()+"pdf_X_REF_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X_REF_PDF())+rest_file;
          this->File.write_to_file(filex, xbaux,this->pdf_ref);
#endif
          this->delta_X_REF_PDF.clear();
          this->delta_X_REF_PDF.shrink_to_fit();
#ifdef _FULL_VERBOSE_
          So.message_screen("Executing rank ordering from DM to DM-target");
#endif
          rankorder(ZERO, xbaux, this->params._NX(),  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);
          So.DONE();
          this->get_power_spectrum("DM_RO"); // from here, Power_DM_REF is loaded
          // Get pdf of the rank-ordered and write it to file
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of delta DM rank-ordered");
#endif
          calc_pdf("log",  this->params._NGRID(),this->params._NX(), this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
          So.DONE();
#ifdef _WRITE_PDF_
          filex=this->params._Output_directory()+"pdf_rank_ordered_iteration"+to_string(this->iteration)+"_X_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
          this->File.write_to_file(filex, xbaux,pdf_in);
#endif
          real_prec lss_bias=0;
          int ncount=0;
          for(int i=0;i<N_MODES;++i)
           if(this->Power_DM_REF[i]>0)
             {
               lss_bias += this->Power_REF[i]/this->Power_DM_REF[i];
               ncount++;
             }
            lss_bias/=static_cast<real_prec>(ncount);
#ifdef _FULL_VERBOSE_
          So.message_screen("New large-Scale Bias: P(k)_tracer / P(k)_dm_RO  = ", sqrt(lss_bias));
          std::cout<<endl;
#endif
#endif // END RANK ORDERING AB INITIO
#ifdef _USE_CWC_
      this->cwclass.get_CWC(this->delta_X_ini);
#ifdef _USE_MASS_KNOTS_
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#else
#if defined (_USE_MASS_KNOTS_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_AWEB_) || defined (_USE_PWEB_) || defined (_USE_TIWEB_)
      this->cwclass.get_CWC(this->delta_X_ini);
#if defined (_USE_MASS_KNOTS_)
      this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#endif
#endif   //endif _USE_CWC_
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_DELTA2_)  || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_)) && (!defined (_USE_CWC_)) && (!defined (_USE_PWEB_))
      this->cwclass.get_bias_terms(this->delta_X_ini);
#endif
    } // end if step==0
  else    // steps >= 1
    {
// *********************************************************************************************
// Compute the kernel
// *********************************************************************************************
      this->GetKernel(true, KERNEL_INDEX);
// #ifdef _HYDROTEST_
//       // ********hydro test, *******************************
//       this->cwclass.Kernel.clear();
//       this->cwclass.Kernel.shrink_to_fit();
//       this->cwclass.Kernel.resize(this->NTT, 0.0);
//       this->cwclass.Kernel=this->Kernel;
//       // ********************************************
// #endif
// *********************************************************************************************
// Convolve the Dm with the kernel and output delta_X
// with this we always convolve the original overdensity field
// *********************************************************************************************
//       if (this->iteration==this->N_iterations_Kernel)
//            this->File.write_array(this->params._Output_directory()+"DM_DELTA_NOconvolved_iteration"+to_string(this->iteration), this->delta_X_ini);
      this->Konvolve(this->delta_X_ini, this->delta_X);
      // Get the PDF of the convolved field
#if defined _WRITE_PDF_  || defined (_USE_GNUPLOT_)
      this->pdf_ite.resize(this->params._NX(), 0);
#ifdef _FULL_VERBOSE_
      So.message_screen("Measuring pdf of log 1+ delta DM convolved");
#endif
      calc_pdf("log", this->params._NGRID(),this->params._NX(), this->Xmax, this->Xmin, this->delta_X, this->pdf_ite);
#ifdef _USE_GNUPLOT_
 vector<pair<real_prec, real_prec> > xy_pts_r;
 vector<pair<real_prec, real_prec> > xy_pts_ra;
 if(this->iteration>2)
  {
    for(int i=0; i<xbaux.size(); ++i)
      xy_pts_r.push_back(std::make_pair(xbaux[i], log10(pdf_ite[i])));
    for(int i=0; i<xbaux.size(); ++i)
      xy_pts_ra.push_back(std::make_pair(xbaux[i], log10(pdf_ref[i])));
    this->gp_kernel<<"set size 0.13,0.2\n";
    this->gp_kernel<<"set origin 0.34,0.06\n";
    this->gp_kernel<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set xrange[0:2] \n";
    this->gp_kernel<<"set yrange[-7:0] \n";  //con este range evito que estos labels salgan en otras subplots
    this->gp_kernel<<"set grid\n";
    this->gp_kernel<<"unset log y\n";
    this->gp_kernel<<"unset log x\n";
    this->gp_kernel<<"set xlabel 'log (1 + {/Symbol d})' enhanced font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_kernel<<"set ylabel 'log P({/Symbol d})' font 'Times-Roman,15'\n";
    this->gp_kernel<<"plot " << gp_kernel.file1d(xy_pts_ra) << "w l lw 3 lt 7 notitle, "<< gp_kernel.file1d(xy_pts_r) << "w l lw 4 lt 31 notitle "<<endl;
  }
  xy_pts_r.clear();
  xy_pts_r.shrink_to_fit();
  xy_pts_ra.clear();
  xy_pts_ra.shrink_to_fit();
#endif
  int index= (this->iteration <=this->params._N_iterations_Kernel())  ?  this->iteration  : this->iteration - (this->params._N_iterations_Kernel())+1;
  string label_aux = this->iteration <= this->params._N_iterations_Kernel() ? "_iteration": "_realization";
  string afilex=this->params._Output_directory()+"pdf_convolved"+label_aux+to_string(index)+"_X_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
  this->File.write_to_file(afilex, xbaux,this->pdf_ite);
#endif
#ifdef _USE_CWC_INSIDE_LOOP_
#ifdef _USE_CWC_
      this->cwclass.get_CWC(this->delta_X);
#ifdef _USE_MASS_KNOTS_
      this->cwclass.get_Mk_collapsing_regions(this->delta_X,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_) || defined (_USE_IWEB_)|| defined (_USE_AWEB_) || defined (_USE_PWEB_) || defined (_USE_IKWEB_) || defined (_USE_TIWEB_)
      this->cwclass.get_CWC(this->delta_X);
#if defined (_USE_MASS_KNOTS_)
      this->cwclass.get_Mk_collapsing_regions(this->delta_X,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#endif
#endif
#endif
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_DELTA2_) || defined (_USE_S2DELTA_) || defined (_USE_DELTA3_) || defined (_USE_S3_) ) && (!defined (_USE_CWC_))
      this->cwclass.get_bias_terms(this->delta_X);
#endif
#ifdef _WRITE_DM_DENSITY_FIELD_
      auto out_it = std::find(std::begin(this->params._output_at_iteration()), std::end(this->params._output_at_iteration()), this->iteration);
      if (this->iteration==this->params._N_iterations_Kernel() || out_it != std::end(this->params._output_at_iteration()))
        this->File.write_array(this->params._Output_directory()+"DM_DELTA_convolved_iteration"+to_string(this->iteration), this->delta_X);
      this->File.write_array(this->params._Output_directory()+"DM_DELTA_convolved", this->delta_X);
#endif
#ifdef _GET_POWER_REFS_
    this->get_power_spectrum("DM_KONV");
#endif
    }//end if step>0
  // At the end of this function, dela_X denotes an overdensity
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MOCK_MODE
// Same goal as get_scaling_relations_primary_property, used when two or more references are used to learn the distribution of properties from
// This function will be called from the makecat funciton as many times as references catags are to be used.
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
void BiasMT::get_scaling_relations_primary_property_two(int ifile)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
 int NTHREADS=_NTHREADS_;
 omp_set_num_threads(NTHREADS);
#endif
  //HEWRE WE HAVE TO USE tracer_ref to read the reference catalog
#ifdef _FULL_VERBOSE_
 if(ifile==0){
    So.message_screen("*************************************************************************");
    So.message_screen("*************************************************************************");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    So.message_screen("***Measuring conditional Vmax Function from reference tracer catalog*****");
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    So.message_screen("***Measuring conditional Mass Function from reference tracer catalog*****");
#endif
    So.message_screen("****as a function of DM properties of the reference DM density field*****");
    So.message_screen("*************************************************************************");
  std::cout<<endl;
  }
#endif
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_a=N_C_BIN1;
#ifdef _USE_TOTAL_MASS_IN_CELL_
  N_a = N_BINS_TOTAL_MASS_IN_CELL;
#endif
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_b=N_C_BIN2;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
  N_b=N_BINS_MIN_DIST_TO_NEI;
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_ || defined _USE_MEAN_SEPARATIONS_IN_CELLS_
  N_b = N_BINS_MIN_SEP_IN_CELLS;
#endif
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_c = N_C_BIN3;
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_v= N_CV_BIN1;
#ifdef _USE_TRACERS_IN_CELLS_
  N_v = N_BINS_TRACERS_IN_CELLS;
#endif
 // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG LENGHT_AB=N_v*N_a* N_b* N_c* N_CV_BIN2*N_CV_BIN3;
  LENGHT_AB*=this->params._n_sknot_massbin() * this->params._n_cwt() * this->params._n_vknot_massbin() * this->params._n_cwv()*this->params._NX();
  // ********************************************************************************
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  ULONG N_x = this->params._NPROPbins_bam();
  ULONG LENGHT_AB_ONE= LENGHT_AB*N_x;
#endif
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   real_prec lm_min=this->params._LOGMASSmin();
#endif
// *********************************************************************************************
// Get min and max of the differnet properties invovled
// *********************************************************************************************
  if(ifile==0)
    this->get_new_min_max_properties();
  this->tracer_ref.set_type_of_object("TRACER_REF");
  string newcat=this->params._files_tracer_references(ifile);
  // *********************************************************************************************
#ifdef _FULL_VERBOSE_
  So.message_screen("*************************************");
  So.message_screen("********Reading Reference TRACER*****");
  So.message_screen("*************************************");
#endif
  //read catalog passing as argument the file and the mininum mass requested
  // *********************************************************************************************
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_ref.read_catalog(newcat,pow(10,this->params._LOGMASSmin())*this->params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,this->params._VMAXmin());
#else
  this->tracer_ref.read_catalog(newcat,0);
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif
#endif
  if(ifile==this->params._Number_of_references()-1)// Measure the abundance of the last refernce,
      this->tracer_ref.get_property_function(this->params._Output_directory()+"tracer_ref_abundance.txt");
#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  So.message_screen("Identifying Neighbouring Cells (this is done once)");
  this->ncells_info.resize( this->params._NGRID());
  get_neighbour_cells(this->params._Nft(), N_CELLS_BACK_FORTH, this->ncells_info);
  So.DONE();
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
#endif
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
 if(ifile==0)
 {
#ifdef _FULL_VERBOSE_
  So.message_warning("Properties are to be assigned to last read reference catalog");
#endif
// This aims at assigning the coordinates of ref to mocks, in case we want to reassign prop to the refernece halos
#ifdef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_   // If we want to assign to a new reference sim,
  this->So.message_warning("!!! Properties are to be assigned to paired reference simulation");
  this->tracer_aux.set_params(this->params);
  this->tracer_aux.set_type_of_object("TRACER_REF");
  string newcat_aux=this->params._Input_dir_cat_new_ref()+this->params._file_catalogue_new_ref();
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer_aux.read_catalog(newcat_aux,pow(10,this->params._LOGMASSmin())*this->params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_aux.read_catalog(newcat_aux,params._VMAXmin());
#else
  this->tracer_aux.read_catalog(newcat_aux,0);
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_aux.read_catalog(newcat_aux,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_aux.read_catalog(newcat_aux,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif
#endif
     this->tracer.Halo.resize(this->tracer_aux._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer._NOBJS(); ++i)
       {
         this->tracer.Halo[i].coord1=this->tracer_aux.Halo[i].coord1;
         this->tracer.Halo[i].coord2=this->tracer_aux.Halo[i].coord2;
         this->tracer.Halo[i].coord3=this->tracer_aux.Halo[i].coord3;
         this->tracer.Halo[i].vel1=this->tracer_aux.Halo[i].vel1;
         this->tracer.Halo[i].vel2=this->tracer_aux.Halo[i].vel2;
         this->tracer.Halo[i].vel3=this->tracer_aux.Halo[i].vel3;
         this->tracer.Halo[i].GridID=this->tracer_aux.Halo[i].GridID;
       }
     this->tracer.set_NOBJS(this->tracer_aux._NOBJS());
#else  // else for  _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
      // Else assign to the *last reference read*.
     this->tracer.Halo.resize(this->tracer_ref._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_ref._NOBJS(); ++i)
       {
         this->tracer.Halo[i].coord1=this->tracer_ref.Halo[i].coord1;
         this->tracer.Halo[i].coord2=this->tracer_ref.Halo[i].coord2;
         this->tracer.Halo[i].coord3=this->tracer_ref.Halo[i].coord3;
         this->tracer.Halo[i].GridID=this->tracer_ref.Halo[i].GridID;
//         this->tracer.Halo[i].vel1=this->tracer_ref.Halo[i].vel1;
//         this->tracer.Halo[i].vel2=this->tracer_ref.Halo[i].vel2;
//         this->tracer.Halo[i].vel3=this->tracer_ref.Halo[i].vel3;
       }
     this->tracer.set_NOBJS(this->tracer_ref._NOBJS());
#endif
 }
#endif
  // *********************************************************************************************
  for(int i=0;i<this->params._Number_of_MultiLevels();++i)
    this->tracers_multi[i][ifile]=this->tracer_ref.params.get_Ntracers_MultiLevels(i);
 // *********************************************************************************************
// We need to load all the info of the tracers together in a single class member when using multilevel. WHen we do not se multilevel, all the info of all refs go into the dm_prop container
  // *****************************Deal with reference catalogs and fields ***************
  // Here we open the reference DM catalog: we have to
  // i ) get delta
  // ii) convolve with kernel from BAM
  // iII) get mins and max
  // if inside iterative mass assignment, convolve with kernel computed from mass power spectra
  // iv) do Iweb or Tweb classification
  // v) convert to log(num_in_log+delta)
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
#endif
  // ********************************************************************************
  this->cwclass_ref.s_cosmo_pars=this->s_cosmo_pars;
  // Ideally Read the reference only in the first iter  ation, all containers with *aux* are not meant to be kept in memmory, they are like dummy
#ifdef _FULL_VERBOSE_
  So.message_screen("*************************************");
  So.message_screen("********Reading Reference DM*********");
  So.message_screen("*************************************");
#endif
  this->delta_dm_aux_mem.clear();
  this->delta_dm_aux_mem.shrink_to_fit();
  this->delta_dm_aux_mem.resize( this->params._NGRID(),0);  //Keep this untouched, will be used along the iterations
  // ********************************************************************************
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_ // we need below the dm, so we save it before trnasforming to delta and then to log
 vector<real_prec> REF_DM( this->params._NGRID(),0);
 REF_DM=delta_dm_aux_mem;
#endif
  File.read_array(this->params._files_dm_references(ifile),this->delta_dm_aux_mem);
  this->mean_aux=get_mean(this->delta_dm_aux_mem);
  get_overdens(this->delta_dm_aux_mem,this->mean_aux, this->delta_dm_aux_mem);
  ULONG NXn=600;// This is as NX but higher to make pdf
#ifdef _RANK_ORDERING_MOCK_GEN_
#ifdef _FULL_VERBOSE_
       So.message_screen("Executing rank ordering from DM to DM-target");
#endif
   this->delta_X_REF_PDF.resize( this->params._NGRID(),0);
   string file_X_ref_pdf=this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X_REF_PDF(); // TBDep
   this->File.read_array_t<PrecType_X>(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_MOCK_GEN_
   get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF); //a better bispect is found if the rank ordering is done to the CIC of number counts, not the delta thereof
#endif
   vector<real_prec>xbaux(NXn, 0);
   vector<real_prec>pdf_in(NXn, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<xbaux.size(); ++i)
     xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(NXn));
   this->pdf_ref.resize(this->params._NX(), 0);
   calc_pdf("log",  this->params._NGRID(),NXn, this->Xmax, this->Xmin, delta_dm_aux_mem, pdf_in);
   calc_pdf("log",  this->params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
   rankorder(ZERO, xbaux, NXn,  this->Xmax, this->Xmin, delta_dm_aux_mem, pdf_in, this->pdf_ref);
   this->delta_X_REF_PDF.clear();
   this->delta_X_REF_PDF.shrink_to_fit();
   xbaux.clear();xbaux.shrink_to_fit();
   pdf_in.clear();pdf_in.shrink_to_fit();
#endif   // end for _RANK_ORDERING_MOCK_GEN_
// **********************************************************************************
#ifdef _KONVOLVE_PASSIGN_
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
  this->Kernel.clear();this->Kernel.resize(this->NTT,0);
  int n_refs=this->params._Number_of_references();
  for(int j=0; j<n_refs; ++j)  // Loop over nrefs-1 refereces. The first reference will be counted here, hence the function get_new_dm field is not used when this function is called
  {
    vector<real_prec>kernel_ghost(this->NTT,0);
    this->File.read_array(this->params._files_kernel_references(j), kernel_ghost); // FIX NAME OF GHOST DM, I have used one of the calibration
    for(ULONG i=0; i<this->NTT; ++i)
       this->Kernel[i]+=kernel_ghost[i]/static_cast<real_prec>(n_refs);
  }
  So.DONE();
#else
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat", this->Kernel);
#endif
  this->Konvolve(this->delta_dm_aux_mem,this->delta_dm_aux_mem);  // COnvolve with BAM Kernel
//  File.write_array(this->params._Output_directory()+"DM_field",dm_ref);
  this->Kernel.clear();
  this->Kernel.shrink_to_fit();
#endif   // NO NEED OF ELSE FOR THIS IF, FOR THE NAME NAME OF THE DELTAS BEFORE AND AFTER KONV HAS BEEN SET TO THE SAME STRING
#if defined _USE_CWC_ || defined _USE_CWC_mass_
      this->cwclass_ref.get_CWC(this->delta_dm_aux_mem);   //get the CW info
#ifdef _USE_MASS_KNOTS_
      this->cwclass_ref.get_Mk_collapsing_regions(this->delta_dm_aux_mem,this->mean_aux);  //get the MK info
#endif //use_mass_knots
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_)  ||  defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)|| defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
   this->cwclass_ref.get_CWC(this->delta_dm_aux_mem);   //get the CW info
#endif
#endif
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
  So.message_screen("Transforming delta_ref -> log10(2+delta_ref). Line ", __LINE__);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i <  this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->delta_dm_aux_mem[i] = (this->delta_dm_aux_mem[i]<-1 ? 0  :  log10(NUM_IN_LOG+ this->delta_dm_aux_mem[i]));
  So.DONE();
  // ********************************************************************************
#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD( this->params._NGRID(),0);
  this->File.read_array(this->params._files_tracer_field_references(ifile) ,REF_DEN_FIELD);
  int nmax=get_max<real_prec>(REF_DEN_FIELD);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum number of tracer in cell", nmax);
#endif
#endif
  // ********************************************************************************
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
    Params params_new=this->params;
    vector<real_prec> CROSS_CORR( this->params._NGRID(),0);
    PowerSpectrumF cpower(params_new);
#ifndef  _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD( this->params._NGRID(),0);
  this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),REF_DEN_FIELD);
#endif
  cpower.get_cross_correlation_config_space(REF_DM,REF_DEN_FIELD,CROSS_CORR);
   File.write_array("test", CROSS_CORR);
  REF_DM.clear(); REF_DM.shrink_to_fit();
#endif
#if !defined (_USE_TRACERS_IN_CELLS_) && defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
  REF_DEN_FIELD.clear();REF_DEN_FIELD.shrink_to_fit();
#endif
  // ********************************************************************************
#ifdef _USE_TOTAL_MASS_IN_CELL_
  vector<real_prec> REF_MASS_FIELD( this->params._NGRID(),0);
  this->tracer_ref.get_density_field_grid(_MASS_, REF_MASS_FIELD);
  real_prec mmax=get_max<real_prec>(REF_MASS_FIELD);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum mass of tracer in cell", mmax);
#endif
#endif
  // ********************************************************************************
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  this->tracer_ref.get_neighbour_tracers(this->ncells_info);
#endif
  // ********************************************************************************
#ifdef _GET_DIST_MIN_SEP_REF_
    this->tracer_ref.get_distribution_min_separations(this->ncells_info);
#endif
// ********************************************************************************
#ifdef _GET_DIST_REDUCED_MASS_
  this->tracer_ref.get_distribution_reduced_mass_in_cell();
#endif
  // ********************************************************************************
#if defined _USE_MIN_SEPARATIONS_IN_CELLS_
  this->tracer_ref.get_stats_separation_in_cell();
#if defined _ASSIGN_PROPERTIES_TO_REFERENCE_ || defined (_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_)
  this->min_halo_separation=this->tracer_ref.min_halo_separation;
#else
  this->min_halo_separation=MIN_SEP_IN_CELLS;
#endif
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-this->min_halo_separation)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
  this->tracer_ref.get_mean_separation_in_cell();
#if defined _ASSIGN_PROPERTIES_TO_REFERENCE_ || defined (_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_)
  this->min_mean_halo_separation=this->tracer_ref.min_halo_separation;
#else
  this->min_mean_halo_separation=MIN_MEAN_SEP_IN_CELLS;
#endif
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-MIN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#endif
  // ********************************************************************************
  if(ifile==0)
  {
#ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
    So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of Vmax (reference) in theta-bins");
#endif
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
    So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses (reference) in theta-bins");
#endif
#endif
    this->dm_properties_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
    this->dm_properties_bins.shrink_to_fit();
    this->dm_properties_bins.resize(LENGHT_AB);
    So.DONE();
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  this->dm_properties_for_randoms_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
  this->dm_properties_for_randoms_bins.shrink_to_fit();
  this->dm_properties_for_randoms_bins.resize(LENGHT_AB);
  So.DONE();
#endif
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB_ONE* (sizeof(ULONG))/(1e9), "Gb for probability distribution");
#endif
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  So.DONE();
#endif
}// closes if ifile==0
  // ********************************************************************************
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  real_prec mean_hmass=0;
  real_prec mean_ldm=0;
  real_prec m_dm=0;
  real_prec mean_ntr=0;
  real_prec m_ntr=0;
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  So.message_screen("**Measuring n(M|theta) from the reference  using ", this->tracer_ref._NOBJS()," objects");
#else
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
    So.message_screen("**Identifying reference masses in bins of {theta}_ref using ", this->tracer_ref._NOBJS()," objects");
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  So.message_screen("**Identifying reference values of Vmax in bins of {theta}_ref with ", this->tracer_ref._NOBJS()," objects from the reference catalog");
#endif
#endif
#endif
#ifdef  _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
    ULONG counter_vmax=0;
    ULONG counter_vmax_r=0;
#endif
   vector<s_cell_info> cell_info_tr( this->params._NGRID());
 // This can be improved by doing two loops. In the first loop, Abundance is computed in a prallel way,
  // while the structore allocating the thera_bin ofin which the id where the galaxy ig licves is allocated.
  // This can then be used in the second loop over grid, filling the dm_properties
  // *********************************************************************************************
// Do not parallelize. There are push_backs inside this loop.
  for(ULONG ig = 0; ig< this->tracer_ref._NOBJS() ; ++ig)     // Open loop over the number of reference tracer objects
    {
      ULONG I_Y=0;
#ifdef _USE_STATISTICAL_ASSIGNMENT_
     // When assigning via reading from the reference values, those values are allocated in the dm_prop structure and then
     // assigned: hence we do not need to ask in which bin of vmax the value of vmax of each tracer is located.
     // Therefore, the variable I_Y is not necessary and is set to zero.
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
      real_prec halo_prop=log10(this->tracer_ref.Halo[ig].mass/this->params._MASS_units());  // Using the logarithm of Mass
      I_Y= get_bin(halo_prop,lm_min,this->params._NMASSbins_mf(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     real_prec halo_prop=log10(this->tracer_ref.Halo[ig].vmax);   //Using the logarithn of Vmax
      I_Y =get_bin(halo_prop, log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#endif
#endif
    // Get ID of this reference dark matter tracers
      ULONG ID=this->tracer_ref.Halo[ig].GridID;
    // Get bin of DM ovedensity
#ifdef _USE_DM_DENSITY_AT_HALO_POSITION_
      real_prec xdm=linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->tracer_ref.Halo[ig].coord1,this->tracer_ref.Halo[ig].coord2,this->tracer_ref.Halo[ig].coord3,this->delta_dm_aux_mem);
      ULONG I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
      this->tracer_ref.Halo[ig].local_dm=xdm;
#else
      real_prec xdm  = static_cast<real_prec>(this->delta_dm_aux_mem[ID]);
      ULONG I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
      this->tracer_ref.Halo[ig].local_dm=xdm;
#endif
      ULONG I_CWT=0;
#ifdef _USE_CWC_
      I_CWT=this->cwclass_ref.get_Tclassification(ID);
      this->tracer_ref.Halo[ig].gal_cwt=this->cwclass_ref.CWClass[ID];
#endif
      ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
      I_MK= (this->cwclass_ref.cwt_used[I_CWT]== I_KNOT ? this->cwclass_ref.SKNOT_M_info[ID]: 0);
#endif
      ULONG I_CWV=0;
#ifdef _USE_CWC_V_
      I_CWV=cwclass_ref.get_Vclassification(ID);
#endif
      ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
      I_VK= (this->cwclass_ref.cwv_used[I_CWV]== I_KNOT ? this->cwclass_ref.VDISP_KNOT_info[ID]: 0);
#endif
      ULONG I_C1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
      I_CV1=get_bin( log10(REF_MASS_FIELD[ID]), this->params._LOGMASSmin(),N_BINS_TOTAL_MASS_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
      real_prec C1 = this->cwclass_ref.Invariant_TF_II[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
      real_prec C1 = this->cwclass_ref.DELTA2[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
      ULONG I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
      I_C2=this->tracer_ref.Number_of_neighbours[ig];
      if(I_C2==N_NEIGHBOURS_MAX)
        I_C2=N_NEIGHBOURS_MAX-1;
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
      real_prec C2 = this->cwclass_ref.Invariant_TF_III[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_
      real_prec min_sep = this->tracer_ref.min_separation_in_cell[ID];
      I_C2 = get_bin(min_sep,  MIN_SEP_IN_CELLS, N_b, delta_min_sep , this->bin_accumulate_borders);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
      real_prec min_sep = this->tracer_ref.min_mean_separation_in_cell[ID];
      I_C2 = get_bin(min_sep,  MIN_SEP_IN_CELLS, N_b, delta_mean_sep , this->bin_accumulate_borders);
#endif
      ULONG I_C3=0;
#if defined _USE_INVARIANT_TIDAL_FIELD_IV_
      real_prec C3 = this->cwclass_ref.Invariant_TF_IV[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_
       real_prec C3 = this->cwclass_ref.Tidal_Anisotropy[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
      real_prec C3 = this->cwclass_ref.S2[ID];             // s²
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
      ULONG I_CV1=0;
#ifdef _USE_TRACERS_IN_CELLS_
      I_CV1=get_bin(static_cast<int>(REF_DEN_FIELD[ID]),0,N_BINS_TRACERS_IN_CELLS,DELTA_TRACERS_IN_CELLS,true);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_I_)
      real_prec CV1 = invariant_field_I(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]); // Not yet assigned to container
      I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
      real_prec CV1 = this->cwclass_ref.N2D[ID];      // Nabla² ð
      I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif
      ULONG I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
      real_prec CV2 = invariant_field_II(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
      real_prec CV2 = this->cwclass_ref.S2DELTA[ID];         // s²ð
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif
      ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
      real_prec CV3 = invariant_field_III(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
      real_prec CV3 = this->cwclass_ref.S3[ID];                                   // s³
      I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
      mean_ldm+=xdm;
      mean_hmass+=halo_prop;
      m_dm+=halo_prop*xdm;
#ifdef _USE_TRACERS_IN_CELLS_
      mean_ntr+=REF_DEN_FIELD[ID];
      m_ntr+=halo_prop*REF_DEN_FIELD[ID];
#endif
#endif
#ifndef _BIN_ACCUMULATE_
      if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
        if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
          if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined(_USE_TIDAL_ANISOTROPY_)
            if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
              if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
                if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
                  if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                      {
#endif
                      ULONG index_ant=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, static_cast<ULONG>(this->params._n_cwt()), static_cast<ULONG>(this->params._n_sknot_massbin()), static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()),N_a, N_b ,N_c, N_v, N_CV_BIN2,N_CV_BIN3);
                      // The container ABUNDANCE will be used used to assign property to mock tracers below the threshold Xthreshold
#ifdef _USE_STATISTICAL_ASSIGNMENT_
                      ULONG index_prop_ab= index_2d(I_Y,index_ant,LENGHT_AB);  // In case that we assign the the paired-*referennce, we need this.If we assign to the same reference, we do not need to get abundance
                      this->ABUNDANCE[index_prop_ab]++;
#endif
#ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                      this->dm_properties_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].mass);
                      this->dm_properties_bins[index_ant].used_property.push_back(false);
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                     // If the ordered index vmax_index (ordered from bottom-to-top in vmax) for each reference tracer is higher than Nran, means that the vmax for this tracer wll be assigned to DMparticles at assignment campaign
                      // while the lowest this->tracer.Ntracers_ran vmax available will be assigned to tracers generated with random coordinates.
#ifdef _BOTTOM_RANDOM_
                      if(this->tracer_ref.Halo[ig].vmax_index >= this->tracer.Ntracers_ran)// This "if" guarranties that the number of vmax contained here are equal to the number of dm- assigned tracers.
#elif defined (_TOP_RANDOM_)                                                                                     // In the Bottom-random case, the dm particles get the highst Vmax values
                      if(this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_dm)  // In the Top_Random case , dm particles get the lowest vmax values
#endif
                       {
                           counter_vmax++;
#endif
                           this->dm_properties_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                           this->dm_properties_bins[index_ant].used_property.push_back(false);
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
                           this->dm_properties_bins[index_ant].tracer_properties_secondary.push_back(this->tracer_ref.Halo[ig].mass);
                           this->dm_properties_bins[index_ant].used_property_secondary.push_back(false);
#endif

#if defined (_ASSIGN_MASS_LINKED_TO_V_FROM_REF_)  || defined _NEW_APPROACH_ASS_ || defined _SORT_AND_JOIN_
                           this->dm_properties_bins[index_ant].index_reference.push_back(ifile+1);
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
                        }
                       else  // if this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_ran, we fill withthe value sof Vmax the vectors of structures for the randoms
                        {
                          this->dm_properties_for_randoms_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                          this->dm_properties_for_randoms_bins[index_ant].used_property.push_back(false);
                           counter_vmax_r++;
                        }
#endif // end  #if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
#endif // end ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#endif //  end defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#endif // end _USE_STATISTICAL_ASSIGNMENT_
#ifndef _BIN_ACCUMULATE_
                      }
#endif
    }
  this->So.DONE();
    //  So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"tracer_ref.Halo.clear() is disabled in order to compute correlation between masses");
     // this->tracer_ref.Halo.clear();
 // this->tracer_ref.Halo.shrink_to_fit();
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  mean_ldm/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  mean_hmass/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  m_dm/=static_cast<real_prec>(this->tracer_ref._NOBJS());
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
  So.message_screen("Correlation log(M) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Correlation log(Vmax) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#endif
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  mean_ntr/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  m_ntr/=static_cast<real_prec>(this->tracer_ref._NOBJS());
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-mass - N_tracers in cells  =", sqrt(fabs(m_ntr-mean_hmass*mean_ntr)));
#endif
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-Vmax - N_tracers in cells =", sqrt(fabs(m_ntr-mean_hmass*mean_ldm)));
  std::cout<<endl;
#endif
#endif
#endif
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  REF_DEN_FIELD.clear();
  REF_DEN_FIELD.shrink_to_fit();
#endif
#ifndef _USE_STATISTICAL_ASSIGNMENT_
  real_prec aux_n=-100;
  real_prec aux_m;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_bins[ih].tracer_properties.size()), aux_n);
       aux_n=aux_m;
      }
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax > Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax < Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif
#endif
  aux_n=-100;
  aux_m=0;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_for_randoms_bins[ih].tracer_properties.size()), aux_n);
       aux_n=aux_m;
      }
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax <= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax >= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif
#endif
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum number of tracers in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
  std::cout<<endl;
#endif
#endif
  ULONG aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in theta-containers:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_a+=this->dm_properties_bins[ih].used_property.size();// using masses_bin_properties leads to the same result.
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  So.message_screen("Number of tracers_dm in theta-containers =", aux_a);
  So.message_screen("check =", counter_vmax);
  So.message_screen("Expected =", this->tracer.Ntracers_dm);
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Number of tracers in theta-containers in =", aux_a);
  So.message_screen("Number of references used =", this->params._Number_of_references());
#endif
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  //Here we sum the tracers with Vmax below the threshold identified below which random tracers will take their vmax
  ULONG aux_b=0;
#pragma omp parallel for reduction(+:aux_b)
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_b+=this->dm_properties_for_randoms_bins[ih].used_property.size();
  if(ifile<this->params._Number_of_references()-1)
    So.message_screen("Partial Number of tracers_random in theta-containers =", aux_b);
  else
    So.message_screen("Total Number of tracers_random in theta-containers =", aux_b);
  So.message_screen("check =", counter_vmax_r);
  So.message_screen("Expected =", this->tracer.Ntracers_ran);
  aux_a+=aux_b;
#endif
  if(aux_a<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in vec<dm_props> dm_properties is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cout<<aux_a<<"  "<<this->tracer_ref._NOBJS()<<endl;
      exit(0);
    }
  So.DONE();
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  if(ifile==this->params._Number_of_references()-1)
   {
   aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in abundance:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE[ih];
  if(aux_a<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cout<<aux_a<<"  "<<this->tracer_ref._NOBJS()<<endl;
      exit(0);
    }
  So.DONE();
  }
#endif
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
if(ifile==this->params._Number_of_references()-1)
  {
#ifdef _FULL_VERBOSE_
    So.message_screen("Normalizing to get joint probability distribution");
#endif
   this->ABUNDANCE_normalized.clear();
   this->ABUNDANCE_normalized.shrink_to_fit();
   this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG idm = 0; idm < LENGHT_AB; ++idm)
     {
      long aux_a=-10000;
      long aux_b;

      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
          ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[indexa];
          aux_b=max(aux_a, AUX);
          aux_a=aux_b;
        }
      for(int tr_j = 0; tr_j < N_x; ++tr_j)
        {
          ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
         }
     }
   }
/*
  aux_a=0;
  So.message_screen("Checking sum in abundance_normalized:");
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE_normalized[ih];
   std::cout<<aux_a<<endl;
*/
#endif
#ifdef _MULTISCALE_
//  If we use MULTISCALE, we need to sort th properties allocated in the dm_properties structure
#ifdef _JOIN_AND_SORT_
    if(ifile==this->params._Number_of_references()-1)
        this->sort_properties();// Sort all tracers in dm_prop only when reading the last reference, i.e, join and sort
#endif
#endif
// ========================================================================================
#ifdef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _FULL_VERBOSE_
  So.DONE();
#endif
  if(ifile==this->params._Number_of_references()-1)
    {
#ifdef _VERBOSE_FREEMEM_
      So.message_screen("Freeing memory");
#endif
      this->ABUNDANCE.clear();
      this->ABUNDANCE.shrink_to_fit();
      So.DONE();
    }
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  //* Final check to verify that all tracers have been counted in
  ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->ABUNDANCE));
  if(ncells_used<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cout<<"Number of cells used = "<<ncells_used<<"  Number expected = "<<this->tracer_ref._NOBJS()<<endl;
    }
   So.message_screen("Assigning memmory space for conditional probability function ABUNDANCE_normalized");
   this->ABUNDANCE_normalized.shrink_to_fit();
   this->ABUNDANCE_normalized.clear();
   this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE, 0);

    ULONG LENGHT_AC  = N_CV_BIN1*N_CV_BIN2*N_CV_BIN3*N_a* N_b* N_c * this->params._n_sknot_massbin() * this->params._n_cwt()*this->params._n_vknot_massbin()*this->params._n_cwv()*this->params._NX();
    this->NCELLSperDMBIN.clear();
    this->NCELLSperDMBIN.shrink_to_fit();
    this->NCELLSperDMBIN.resize(LENGHT_AC, 0);
    So.DONE();
    So.message_screen("**Normalizing number counts and marginalizing with respect to Y-bins:");
#ifdef _USE_OMP_
#pragma omp parallel for collapse(11)
#endif
    for(int i=0;i < this->params._NX();++i)// loop sobre bines de dark matter
      for(int sua = 0; sua < this->params._n_cwt(); ++sua)
        for(int k = 0; k < this->params._n_sknot_massbin() ;  ++k)
          for(int vua = 0; vua < this->params._n_cwv(); ++vua)
            for(int vk = 0; vk < this->params._n_vknot_massbin() ;  ++vk)
              for(int l1 = 0; l1 < N_a; ++l1)
                for(int l2 = 0; l2 < N_b ; ++l2)
                  for(int l3 = 0; l3 < N_c; ++l3)
                    for(int lv1 = 0; lv1 < N_CV_BIN1; ++lv1)
                      for(int lv2 = 0; lv2 < N_CV_BIN2; ++lv2)
                        for(int lv3 = 0; lv3 < N_CV_BIN3; ++lv3)
                           {
                                ULONG index_l=index_11d(i,sua,k,vua,vk,l1,l2,l3,lv1,lv2,lv3,this->params._n_cwt(),this->params._n_sknot_massbin(),this->params._n_cwv(),this->params._n_vknot_massbin(),N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                long aux_a=-10000;
                                long aux_b;
                                for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
                                  {
                                    ULONG indexa=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1, lv2, lv3,this->params._NX(),this->params._n_cwt(), this->params._n_sknot_massbin(),this->params._n_cwv(), this->params._n_vknot_massbin(), N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                    long AUX=this->ABUNDANCE[indexa];
                                    aux_b=max(aux_a, AUX);
                                    aux_a=aux_b;
                                    ULONG index_h=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1,lv2,lv3,this->params._NX(),this->params._n_cwt(),this->params._n_sknot_massbin(),this->params._n_cwv(),this->params._n_vknot_massbin(),N_a, N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                    this->NCELLSperDMBIN[index_l]+=this->ABUNDANCE[index_h];
                                  }
                                for(int tr_j = 0; tr_j < this->params._NMASSbins_mf(); ++tr_j)
                                  {
                                    ULONG indexa=index_12d(tr_j,i,sua,k,vua,vk,l1,l2,l3,lv1, lv2, lv3,this->params._NX(),this->params._n_cwt(), this->params._n_sknot_massbin(),this->params._n_cwv(), this->params._n_vknot_massbin(), N_a,N_b,N_c, N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                                    long AUX=this->ABUNDANCE[indexa];
                                    this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
                                  }
                           }
    this->So.DONE();
#ifdef _FULL_VERBOSE_
    So.message_screen("Check on number of bins used...");
#endif
    real_prec aux_h=0;
#pragma omp parallel for reduction(+:aux_h)
    for(int ih=0;ih<this->ABUNDANCE_normalized.size();++ih)
      aux_h+=this->ABUNDANCE_normalized[ih];
    if(aux_h<=0)
      {
        So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts ill-defined. CosmicAtlas stops here.");
        exit(0);
      }
    else
      this->So.DONE();
#ifdef _VERBOSE_FREEMEM_
    So.message_screen("Freeing memmory");
#endif
//    this->ABUNDANCE.clear();
//    this->ABUNDANCE.shrink_to_fit();
    So.DONE();
    // We do not need here to write these files as they will be inmediatly used in the assign_tracer_mass
    //  this->File.write_array(this->params._Output_directory()+"Bam_Abundance", this->ABUNDANCE);
    // this->File.write_array(this->params._Output_directory()+"Bam_Abundance_Normalized", this->ABUNDANCE_normalized);
// Update a bif tracaer_ref class member to allcoate all properties of merged files
#endif  // end for #ifdef _USE_STATISTICAL_ASSIGNMENT_
}
#endif
//end class member function get_scaling_relations_primary_property_two()
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function gets the scaling relation P(V|theta)
// and must be in accordance with what BiasMT::get_scaling_relations_primary_property()
void BiasMT::get_scaling_relations_primary_property()
{
  this->So.enter(__PRETTY_FUNCTION__);
  //HEWRE WE HAVE TO USE tracer_ref to read the reference catalog
#ifdef _FULL_VERBOSE_
  So.message_screen("*************************************************************************");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  So.message_screen("***Measuring conditional Vmax Function from reference tracer catalog*****");
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
  So.message_screen("***Measuring conditional Mass Function from reference tracer catalog*****");
#endif
  So.message_screen("****as a function of DM properties of the reference DM density field*****");
  So.message_screen("*************************************************************************");
  cout<<endl;
#endif
  // ******************************************************************************
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  ULONG N_a=N_C_BIN1;
#ifdef _USE_TOTAL_MASS_IN_CELL_
  N_a = N_BINS_TOTAL_MASS_IN_CELL;
#elif defined _USE_STDV_SEPARATIONS_IN_CELLS_
   N_a= N_BINS_STDV_SEP_IN_CELLS ;
#endif
   // ******************************************************************************
  ULONG N_b=1;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
  N_b=N_BINS_MIN_DIST_TO_NEI;
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_
  N_b = N_BINS_MIN_SEP_IN_CELLS;
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
  N_b = N_BINS_MEAN_SEP_IN_CELLS;
#endif
  // ******************************************************************************
  ULONG N_c=1; // this applies to tidal anisotropy, inv tidal field iv, s2
#ifdef _USE_TIDAL_ANISOTROPY_
  N_c = N_C_BIN3;
#elif defined _USE_MACH_NUMBER_
  N_c = N_BINS_MACH;
#endif
  // ******************************************************************************
  ULONG N_v=1;
#ifdef _USE_TRACERS_IN_CELLS_
  N_v = N_BINS_TRACERS_IN_CELLS;
#endif
  // ******************************************************************************
  ULONG N_x=1;
#ifdef _USE_LOCAL_OVERDENSITY_
  N_x = N_BINS_LO;
#endif
  // ******************************************************************************
  ULONG N_bias=1;
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
   N_bias=N_BINS_BIAS;
#endif
  // ******************************************************************************
  ULONG LENGHT_AB=N_v*N_a* N_b* N_c*N_x*N_bias;
  LENGHT_AB*=this->params._n_sknot_massbin() * this->params._n_cwt() * this->params._n_vknot_massbin() * this->params._n_cwv()*this->params._NX();
  // ******************************************************************************
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  ULONG N_x = this->params._NPROPbins_bam();
  ULONG LENGHT_AB_ONE= LENGHT_AB*N_x;
#endif
  // ********************************************************************************
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   real_prec lm_min=this->params._LOGMASSmin();
#endif
   // ******************************************************************************
// Get min and max of the differnet properties invovled
  this->get_new_min_max_properties();
  this->tracer_ref.set_params(this->params);
  this->tracer_ref.set_type_of_object("TRACER_REF");
  string newcat=this->params._Input_dir_cat()+this->params._file_catalogue();
  //read catalog passing as argument the file and the mininum mass requested
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
      // Read catalog: if multilevel is enabled, this object updates params, so we must bring that params to bmt again below, at update_params
     this->tracer_ref.read_catalog(newcat,pow(10,this->params._LOGMASSmin())*this->params._MASS_units());   // READ INPUT CATALOG //
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,params._VMAXmin());
#else
  this->tracer_ref.read_catalog(newcat,0);
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_ref.read_catalog(newcat,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_ref.read_catalog(newcat,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif
#endif
 // *********************************************************************************************
update_params:
#ifdef _USE_FIXED_MULTISCALE_LEVELS_
  this->params=this->tracer_ref.params; // Update params
#endif
  // *********************************************************************************************
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_   // THis applies in case we want to assign properties to the same catalog where we learnt from (reference)
#ifdef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_   // If we want to assign to the paired version of the reference using the paired (e.g. using UNITSIM), then read the paired: we assume that the two files have the smae name at different dirs
  this->So.message_screen("!!! Properties are to be assigned to a new referencs simulation");
  this->tracer_aux.set_params(this->params);
  this->tracer_aux.set_type_of_object("TRACER_REF");
  string newcat_aux=this->params._Input_dir_cat_new_ref()+this->params._file_catalogue_new_ref();
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer_aux.read_catalog(newcat_aux,pow(10,this->params._LOGMASSmin())*this->params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_aux.read_catalog(newcat_aux,params._VMAXmin());
#else
  this->tracer_aux.read_catalog(newcat_aux,0);
#endif
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
  this->tracer_aux.read_catalog(newcat_aux,pow(10,params._LOGMASSmin())*params._MASS_units(),static_cast<real_prec>(BIG_NUMBER));
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
  this->tracer_aux.read_catalog(newcat_aux,params._VMAXmin(),static_cast<real_prec>(BIG_NUMBER));
#endif
#endif
    this->tracer.set_params(this->params);
    this->tracer.Halo.resize(this->tracer_aux._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_aux._NOBJS(); ++i)
       {
         this->tracer.Halo[i].coord1=this->tracer_aux.Halo[i].coord1;
         this->tracer.Halo[i].coord2=this->tracer_aux.Halo[i].coord2;
         this->tracer.Halo[i].coord3=this->tracer_aux.Halo[i].coord3;
         this->tracer.Halo[i].vel1=this->tracer_aux.Halo[i].vel1;
         this->tracer.Halo[i].vel2=this->tracer_aux.Halo[i].vel2;
         this->tracer.Halo[i].vel3=this->tracer_aux.Halo[i].vel3;
         this->tracer.Halo[i].GridID=this->tracer_aux.Halo[i].GridID;
         this->tracer.Halo[i].GridID_n=this->tracer_aux.Halo[i].GridID_n;
       }
    this->tracer.set_NOBJS(this->tracer_aux._NOBJS());
#else  // end of assign_to_new_ref. iF assign to the refernece directly
  this->So.message_warning("!!!Properties are to be assigned to reference simulation");
  this->tracer.Halo.resize(this->tracer_ref._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_ref._NOBJS(); ++i) // Copy the information of the reference catalog to the "new" catalog, except for the properties to be assigned
       {
         this->tracer.Halo[i].coord1=this->tracer_ref.Halo[i].coord1;
         this->tracer.Halo[i].coord2=this->tracer_ref.Halo[i].coord2;
         this->tracer.Halo[i].coord3=this->tracer_ref.Halo[i].coord3;
         this->tracer.Halo[i].vel1=this->tracer_ref.Halo[i].vel1;
         this->tracer.Halo[i].vel2=this->tracer_ref.Halo[i].vel2;
         this->tracer.Halo[i].vel3=this->tracer_ref.Halo[i].vel3;
         this->tracer.Halo[i].GridID=this->tracer_ref.Halo[i].GridID;
         this->tracer.Halo[i].GridID_n=this->tracer_ref.Halo[i].GridID_n;
         this->tracer.Halo[i].mass=0;
         this->tracer.Halo[i].vmax=0;
       }
    this->tracer.set_NOBJS(this->tracer_ref._NOBJS());
    this->tracer.set_type_of_object("TRACER_MOCK");
#endif
#endif
  this->tracer_aux.Halo.clear();  this->tracer_aux.Halo.shrink_to_fit();
  this->tracer_ref.get_property_function(this->params._Output_directory()+"tracer_ref_abundance.txt");
#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  So.message_screen("Identifying Neighbouring Cells (this is done once)");
  this->ncells_info.resize( this->params._NGRID());
  get_neighbour_cells(this->params._Nft(), N_CELLS_BACK_FORTH, this->ncells_info);
  So.DONE();
#endif
  // *****************************Deal with reference catalogs and fields ***************
  // Here we open the reference DM catalog: we have to
  // i ) get delta
  // ii) convolve with kernel from BAM
  // iII) get mins and max
  // if inside iterative mass assignment, convolve with kernel computed from mass power spectra
  // iv) do Iweb or Tweb classification
  // v) convert to log(num_in_log+delta)
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
#endif
  this->cwclass_ref.s_cosmo_pars=this->s_cosmo_pars;
  // Ideally Read the reference only in the first iter  ation, all containers with *aux* are not meant to be kept in memmory, they are like dummy
#ifdef _FULL_VERBOSE_
  So.message_screen("********Reading Reference DM*********");
#endif
  this->delta_dm_aux_mem.clear();
  this->delta_dm_aux_mem.shrink_to_fit();
  this->delta_dm_aux_mem.resize( this->params._NGRID(),0);  //Keep this untouched, will be used along the iterations
  File.read_array(this->params._Input_Directory_X()+this->params._Name_Catalog_X(),this->delta_dm_aux_mem);
     // ********************************************************************************
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_ // we need below the dm, so we save it before trnasforming to delta and then to log
   vector<real_prec> REF_DM( this->params._NGRID(),0);
   REF_DM=delta_dm_aux_mem;
#endif
    // ********************************************************************************
bias_pr:
   if(this->params._Get_tracer_bias())
    { // WE leave the param question separated from the prerpc, as we would want to use bias to plot but not to assign props
#ifndef _READ_BIAS_
     PowerSpectrumF power_new(this->params, true);
     power_new.object_by_object_bias(this->tracer_ref.Halo,this->delta_dm_aux_mem);
#else
    string out_bias=params._Output_directory()+"individual_bias.txt";
    ifstream bout; bout.open(out_bias.c_str());
    this->So.message_screen("Reading bias from  file ", out_bias);
    for(ULONG i=0;i<this->tracer_ref.Halo.size();++i)
      bout>>this->tracer_ref.Halo[i].bias>>this->tracer_ref.Halo[i].mach_number>>this->tracer_ref.Halo[i].local_overdensity>>this->tracer_ref.Halo[i].number_of_neighbours;
    bout.close();
#endif
#ifndef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_ref._NOBJS(); ++i)//Copy the bias from the "ref" to the "new" catalog. Since the bias is computed from the coords, both can have the same unless coords change
         this->tracer.Halo[i].bias=this->tracer_ref.Halo[i].bias;
#endif
     this->tracer_ref.set_min_bias(this->tracer_ref.get_min("_BIAS_"));
     this->tracer_ref.set_max_bias(this->tracer_ref.get_max("_BIAS_"));
     So.message_screen("Minimum bias",this->tracer_ref._min_bias());
     So.message_screen("Maximum bias",this->tracer_ref._max_bias());
   }
     // ********************************************************************************
#ifdef _USE_MACH_NUMBER_
#ifndef _READ_BIAS_
#ifdef _USE_CHUNCKS_NC_
    this->tracer_ref.get_local_mach_number_chuncks(this->params._Scale_mach_number());//here get mach and delta5
#else
    this->tracer_ref.get_local_mach_number(this->params._Scale_mach_number());//here get mach and delta5
#endif
#endif
#ifndef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_ // if this is not defined, measn that we might be assigning to the same reference so we just replicate it
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_ref._NOBJS(); ++i)//Copy the mach_number from the "ref" to the "new" catalog. Since mach is computed from the phase-space coords, both can have the same unless coords change
         this->tracer.Halo[i].mach_number=this->tracer_ref.Halo[i].mach_number ;
#endif
     this->tracer_ref.set_min_mach(this->tracer_ref.get_min("_MACH_"));
     this->tracer_ref.set_max_mach(this->tracer_ref.get_max("_MACH_"));
     So.message_screen("Minimum mach",this->tracer_ref._min_mach());
     So.message_screen("Maximum mach",this->tracer_ref._max_mach());
#endif
     // ********************************************************************************
#ifdef _USE_LOCAL_OVERDENSITY_
#ifndef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i<this->tracer_ref._NOBJS(); ++i)//Copy the local_delta from the "ref" to the "new" catalog. Since local_delta is computed from the phase-space coords, both can have the same unless coords change
         this->tracer.Halo[i].local_overdensity=this->tracer_ref.Halo[i].local_overdensity;
#endif
     this->tracer_ref.set_min_local_overdensity(this->tracer_ref.get_min("_LOCAL_OVERDENSITY_"));
     this->tracer_ref.set_max_local_overdensity(this->tracer_ref.get_max("_LOCAL_OVERDENSITY_"));
     So.message_screen("Minimum local-overd",this->tracer_ref._min_local_overdensity());
     So.message_screen("Maximum local-overd",this->tracer_ref._max_local_overdensity());
#ifdef _WRITE_BIAS_
    string out_bias=params._Output_directory()+"individual_bias.txt";
    ofstream bout; bout.open(out_bias.c_str());
    this->So.message_screen("Writing bias in  file ", out_bias);
    for(ULONG i=0;i<this->tracer_ref.Halo.size();++i)
      bout<<this->tracer_ref.Halo[i].bias<<"\t"<<this->tracer_ref.Halo[i].mach_number<<"\t"<<this->tracer_ref.Halo[i].local_overdensity<<"\t"<<this->tracer_ref.Halo[i].number_of_neighbours<<endl;
    bout.close();
#endif
#endif
     // ********************************************************************************
#if defined _USE_TIDAL_ANISOTROPY_ || defined _USE_TIDAL_ANISOTROPY_SEC_PROP_
       vector<real_prec>tidal(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._NGRID() ;++i)
        tidal[i]=tidal_anisotropy(this->cwclass.lambda1[i], this->cwclass.lambda2[i], this->cwclass.lambda3[i]);
      this->tracer_ref.get_tracer_tidal_anisotropy(tidal);
      So.message_screen("\tMin TA :",this->tracer_ref._min_tidal_anisotropy());
      So.message_screen("\tMax TA :",this->tracer_ref._max_tidal_anisotropy());
      tidal.clear();tidal.shrink_to_fit();
#endif
     // ********************************************************************************
     this->mean_aux=get_mean(this->delta_dm_aux_mem);
     get_overdens(this->delta_dm_aux_mem,this->mean_aux, this->delta_dm_aux_mem);
// **********************************************************************************
#if defined _USE_CWC_ || defined _USE_CWC_mass_
      this->cwclass_ref.get_CWC(this->delta_dm_aux_mem);   //get the CW info
#ifdef _USE_MASS_KNOTS_
      this->cwclass_ref.get_Mk_collapsing_regions(this->delta_dm_aux_mem,this->mean_aux);  //get the MK info
#endif //use_mass_knots
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_)  ||  defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_ANISOTROPY_AT_TRACER_POSITION_ ) || defined (_USE_ELLIPTICITY_) || defined (_USE_PROLATNESS_) || defined (_USE_PWEB_)
   this->cwclass_ref.get_CWC(this->delta_dm_aux_mem);   //get the CW info
#endif
#endif
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
  So.message_screen("Transforming delta_ref -> log10(2+delta_ref). Line ", __LINE__);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i <  this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
    this->delta_dm_aux_mem[i] = (this->delta_dm_aux_mem[i]<-1 ? 0  :  log10(NUM_IN_LOG+ this->delta_dm_aux_mem[i]));
  So.DONE();
#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD( this->params._NGRID(),0);
  this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),REF_DEN_FIELD);
//  this->tracer_ref.get_density_field_grid(_COUNTS_,REF_DEN_FIELD);
  int nmax=get_max<real_prec>(REF_DEN_FIELD);
  So.message_screen("Maximum number of tracer in cell", nmax);
#endif
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
    Params params_new=this->params;
    vector<real_prec> CROSS_CORR( this->params._NGRID(),0);
    PowerSpectrumF cpower(params_new);
#ifndef  _USE_TRACERS_IN_CELLS_
  vector<real_prec> REF_DEN_FIELD( this->params._NGRID(),0);
  this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),REF_DEN_FIELD);
#endif
  cpower.get_cross_correlation_config_space(REF_DM,REF_DEN_FIELD,CROSS_CORR);
   File.write_array("test", CROSS_CORR);
  REF_DM.clear(); REF_DM.shrink_to_fit();
#endif
#if !defined (_USE_TRACERS_IN_CELLS_) && defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
  REF_DEN_FIELD.clear();REF_DEN_FIELD.shrink_to_fit();
#endif
#ifdef _USE_TOTAL_MASS_IN_CELL_
  vector<real_prec> REF_MASS_FIELD( this->params._NGRID(),0);
  this->tracer_ref.get_density_field_grid(_VMAX_, REF_MASS_FIELD);
  real_prec mmax=get_max<real_prec>(REF_MASS_FIELD);
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum PROP of tracer in cell", mmax);
#endif
#endif
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  this->tracer_ref.get_neighbour_tracers(this->ncells_info);
#endif
#ifdef _GET_DIST_MIN_SEP_REF_
    this->tracer_ref.get_distribution_min_separations(this->ncells_info);
#endif
#ifdef _GET_DIST_REDUCED_MASS_
  this->tracer_ref.get_distribution_reduced_mass_in_cell();
#endif
#if defined _USE_MIN_SEPARATIONS_IN_CELLS_
  this->tracer_ref.get_stats_separation_in_cell();
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
  this->min_halo_separation=this->tracer_ref.min_halo_separation;
#else
  this->min_halo_separation=MIN_SEP_IN_CELLS;
#endif
  real_prec delta_min_sep = (MAX_SEP_IN_CELLS-MIN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
  this->tracer_ref.get_stats_separation_in_cell(); // tghis method also computes the men and stdev ojo
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
  this->min_mean_halo_separation= MIN_MEAN_SEP_IN_CELLS;// this->tracer_ref.min_mean_halo_separation;
#else
  this->min_mean_halo_separation=MIN_SEP_IN_CELLS;
#endif
  real_prec delta_mean_sep = DELTA_MEAN_SEP;//(MAX_MEAN_SEP_IN_CELLS-MIN_MEAN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MEAN_SEP_IN_CELLS);
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
  N_a= N_BINS_STDV_SEP_IN_CELLS ;
#if !defined _USE_MEAN_SEPARATIONS_IN_CELLS_
    this->tracer_ref.get_stats_separation_in_cell(); // this also works to get the mean and stedev
#endif
  real_prec delta_stdv_sep = (MAX_STDV_SEP_IN_CELLS-MIN_STDV_SEP_IN_CELLS)/static_cast<real_prec>(N_a);
#endif
  // ********************************************************************************
#ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB*(sizeof(real_prec)+sizeof(bool))/(1e9), "Gb for list of Vmax (reference) in theta-bins");
#endif
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB* (sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for list of masses (reference) in theta-bins");
#endif
#endif
  this->dm_properties_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
  this->dm_properties_bins.shrink_to_fit();
  this->dm_properties_bins.resize(LENGHT_AB);
  So.DONE();
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  this->dm_properties_for_randoms_bins.clear();            // Container of structures meant to keep the properties of the tracers found in a bin of Theta
  this->dm_properties_for_randoms_bins.shrink_to_fit();
  this->dm_properties_for_randoms_bins.resize(LENGHT_AB);
  So.DONE();
#endif
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _FULL_VERBOSE_
  So.message_screen("Allocating", LENGHT_AB_ONE* (sizeof(ULONG))/(1e9), "Gb for probability distribution");
#endif
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  So.DONE();
#endif
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  real_prec mean_hmass=0;
  real_prec mean_ldm=0;
  real_prec m_dm=0;
  real_prec mean_ntr=0;
  real_prec m_ntr=0;
#endif
#ifdef _USE_STATISTICAL_ASSIGNMENT_
  So.message_screen("**Measuring n(Vmax|{Theta}) from the reference using ", this->tracer_ref._NOBJS()," objects");
#else
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
    So.message_screen("**Identifying reference masses in bins of {theta}_ref using ", this->tracer_ref._NOBJS()," objects");
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  So.message_screen("**Identifying reference values of Vmax in bins of {Theta}_ref with ", this->tracer_ref._NOBJS()," objects from the reference catalog");
#endif
#endif
#endif
#ifdef  _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
    ULONG counter_vmax=0;
    ULONG counter_vmax_r=0;
#endif
#if defined (_MULTISCALE_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_) 
    vector<s_cell_info> cell_info_tr( this->params._NGRID());
#endif
 // This can be improved by doing two loops. In the first loop, Abundance is computed in a prallel way,
  // while the structore allocating the thera_bin ofin which the id where the galaxy ig licves is allocated.
  // This can then be used in the second loop over grid, filling the dm_properties
// Open loop over the number of tracer objects
// Do not parallelize when reading from the reference. There are push_backs inside this loop.
   for(ULONG ig = 0; ig< this->tracer_ref._NOBJS() ; ++ig)
     {
       ULONG I_Y=0;
#ifdef _USE_STATISTICAL_ASSIGNMENT_
     // When assigning via reading from the reference values, those values are allocated in the dm_prop structure and then
     // assigned: hence we do not need to ask in which bin of vmax the value of vmax of each tracer is located.
     // Therefore, the variable I_Y is not necessary and is set to zero.
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
      real_prec halo_prop=log10(this->tracer_ref.Halo[ig].mass/this->params._MASS_units());  // Using the logarithm of Mass
      I_Y= get_bin(halo_prop,lm_min,this->params._NMASSbins_mf(),this->tracer_ref.logdeltaM,this->bin_accumulate_borders);//ojo con el numero de bies cá, arreglarlo
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     real_prec halo_prop=log10(this->tracer_ref.Halo[ig].vmax);   //Using the logarithn of Vmax
      I_Y =get_bin(halo_prop, log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#endif
#endif
      ULONG ID=this->tracer_ref.Halo[ig].GridID; // Get ID of this reference dark matter tracers
#ifdef _USE_DM_DENSITY_AT_HALO_POSITION_
      // Get bin of DM ovedensity
      real_prec xdm=linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->tracer_ref.Halo[ig].coord1,this->tracer_ref.Halo[ig].coord2,this->tracer_ref.Halo[ig].coord3,this->delta_dm_aux_mem);
      ULONG I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
      this->tracer_ref.Halo[ig].local_dm=xdm;
#else
      real_prec xdm  = static_cast<real_prec>(this->delta_dm_aux_mem[ID]);
      ULONG I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
      this->tracer_ref.Halo[ig].local_dm=xdm;
#endif
      ULONG I_CWT=0;
#if defined _USE_CWC_ || defined (_USE_MASS_KNOTS_)
      I_CWT= static_cast<ULONG>(this->cwclass_ref.get_Tclassification(ID));
      this->tracer_ref.Halo[ig].gal_cwt=this->cwclass_ref.CWClass[ID];
#endif
      ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
      I_MK= static_cast<ULONG>(this->cwclass_ref.cwt_used[I_CWT]== I_KNOT ? this->cwclass_ref.SKNOT_M_info[ID]: 0);
#endif
      ULONG I_CWV=0;
#ifdef _USE_CWC_V_
      I_CWV=cwclass_ref.get_Vclassification(ID);
#endif
      ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
      I_VK= (this->cwclass_ref.cwv_used[I_CWV]== I_KNOT ? this->cwclass_ref.VDISP_KNOT_info[ID]: 0);
#endif
      ULONG I_C1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
      I_CV1=get_bin( log10(REF_MASS_FIELD[ID]), this->params._LOGMASSmin(),N_BINS_TOTAL_MASS_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
      real_prec C1 = this->cwclass_ref.Invariant_TF_II[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
      real_prec C1 = this->cwclass_ref.DELTA2[ID];
      I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
#if defined _USE_STDV_SEPARATIONS_IN_CELLS_ && !defined  (_USE_TOTAL_MASS_IN_CELL_)
       real_prec sep_std = this->tracer_ref.stdv_separation_in_cell[ID];
       I_C1 = get_bin(sep_std, MIN_STDV_SEP_IN_CELLS, N_a, delta_stdv_sep , this->bin_accumulate_borders);
#endif
      ULONG I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
      I_C2=this->tracer_ref.Number_of_neighbours[ig];
      if(I_C2==N_NEIGHBOURS_MAX)
        I_C2=N_NEIGHBOURS_MAX-1;
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
      real_prec C2 = this->cwclass_ref.Invariant_TF_III[ID];
      I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_
      real_prec min_sep = this->tracer_ref.min_separation_in_cell[ID];
      I_C2 = get_bin(min_sep,  MIN_SEP_IN_CELLS, N_b, delta_min_sep , this->bin_accumulate_borders);
#elif  defined _USE_MEAN_SEPARATIONS_IN_CELLS_
      real_prec min_sep = this->tracer_ref.mean_separation_in_cell[ID];
      I_C2 = get_bin(min_sep,  MIN_MEAN_SEP_IN_CELLS, N_b, delta_mean_sep , this->bin_accumulate_borders);
#endif
      ULONG I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
      real_prec C3 = this->cwclass_ref.Invariant_TF_IV[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_
#ifdef _USE_ANISOTROPY_AT_TRACER_POSITION_
      real_prec C3  =linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->tracer_ref.Halo[ig].coord1,this->tracer_ref.Halo[ig].coord2,this->tracer_ref.Halo[ig].coord3,this->cwclass_ref.Tidal_Anisotropy);
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#else
      I_C3= get_bin(this->cwclass_ref.Tidal_Anisotropy[ID], this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
#elif defined _USE_S2_
      I_C3= get_bin(this->cwclass_ref.S2[ID], this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_MACH_NUMBER_
      I_C3= get_bin(this->tracer_ref.Halo[ig].mach_number,this->tracer_ref._min_mach(), N_BINS_MACH, (this->tracer_ref._max_mach()-this->tracer_ref._min_mach())/static_cast<real_prec>(N_BINS_MACH), this->bin_accumulate_borders);
#endif
       ULONG I_CV1=0;
#ifdef _USE_TRACERS_IN_CELLS_
      I_CV1=get_bin(static_cast<int>(REF_DEN_FIELD[ID]),0,N_BINS_TRACERS_IN_CELLS,DELTA_TRACERS_IN_CELLS,true);
//      I_CV1=static_cast<int>(REF_DEN_FIELD[ID]);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_I_)
      real_prec CV1 = invariant_field_I(cwclass_ref.lambda1_vs[ID],cwc_ref.lambda2_vs[ID],cwc_ref.lambda3_vs[ID]); // Not yet assigned to container
      I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
      I_CV1= get_bin(this->cwclass_ref.N2D[ID], this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif
      ULONG I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
      real_prec CV2 = invariant_field_II(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
      real_prec CV2 = this->cwclass_ref.S2DELTA[ID];         // s²ð
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined _USE_CROSS_CORRELATION_CONF_SPACE_
      real_prec CV2 = CROSS_CORR[ID];         //cross correlation between dm and halo number counts
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined _USE_LOCAL_OVERDENSITY_
      I_CV2= get_bin(this->tracer_ref.Halo[ig].local_overdensity, this->tracer_ref._min_local_overdensity(), N_BINS_LO,  (this->tracer_ref._max_local_overdensity()-this->tracer_ref._min_local_overdensity())/static_cast<real_prec>(N_BINS_LO) ,this->bin_accumulate_borders);
#endif
      ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
      real_prec CV3 = invariant_field_III(cwclass_ref.lambda1_vs[ID],cwclass_ref.lambda2_vs[ID],cwclass_ref.lambda3_vs[ID]);// Not yet assigned to container
      I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
      I_CV3= get_bin(this->cwclass_ref.S3[ID], this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#elif defined _USE_BIAS_OBJECT_TO_OBJECT_
      I_CV3 = get_bin(this->tracer_ref.Halo[ig].bias, this->tracer_ref._min_bias(), N_BINS_BIAS,(this->tracer_ref._max_bias()-this->tracer_ref._min_bias())/static_cast<real_prec>(N_BINS_BIAS),this->bin_accumulate_borders);
#endif
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
      mean_ldm+=xdm;
      mean_hmass+=halo_prop;
      m_dm+=halo_prop*xdm;
#ifdef _USE_TRACERS_IN_CELLS_
      mean_ntr+=REF_DEN_FIELD[ID];
      m_ntr+=halo_prop*REF_DEN_FIELD[ID];
#endif
#endif
#ifndef _BIN_ACCUMULATE_
      if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
        if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
          if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_IV_) || defined(_USE_TIDAL_ANISOTROPY_)
            if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
              if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
                if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
                  if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                     {
#endif
                      ULONG index_ant=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, static_cast<ULONG>(this->params._n_cwt()), static_cast<ULONG>(this->params._n_sknot_massbin()), static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()),N_a, N_b ,N_c, N_v, N_x,N_bias);
                      // The container ABUNDANCE will be used in statistical mode
#ifdef  _USE_STATISTICAL_ASSIGNMENT_
                      this->ABUNDANCE[index_2d(I_Y,index_ant,LENGHT_AB)]++;
#endif
                      // This will be used for mocks and assign_to_ref for X>=Xthres level 1 only
#ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                      this->dm_properties_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].mass);
                      this->dm_properties_bins[index_ant].used_property.push_back(false);
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                     // If the ordered index vmax_index (ordered from bottom-to-top in vmax) for each reference tracer is higher than Nran, means that the vmax for this tracer wll be assigned to DMparticles at assignment campaign
                      // while the lowest this->tracer.Ntracers_ran vmax available will be assigned to tracers generated with random coordinates.
#ifdef _BOTTOM_RANDOM_
                      if(this->tracer_ref.Halo[ig].vmax_index >= this->tracer.Ntracers_ran)// This "if" guarranties that the number of vmax contained here are equal to the number of dm- assigned tracers.
#elif defined (_TOP_RANDOM_)                                                                                     // In the Bottom-random case, the dm particles get the highst Vmax values
                      if(this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_dm)  // In the Top_Random case , dm particles get the lowest vmax values
#endif
                       {
                           counter_vmax++;
#endif
                           this->dm_properties_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                           this->dm_properties_bins[index_ant].used_property.push_back(false);
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
                        }
                       else  // if this->tracer_ref.Halo[ig].vmax_index < this->tracer.Ntracers_ran, we fill withthe value sof Vmax the vectors of structures for the randoms
                        {
                          this->dm_properties_for_randoms_bins[index_ant].tracer_properties.push_back(this->tracer_ref.Halo[ig].vmax);
                          this->dm_properties_for_randoms_bins[index_ant].used_property.push_back(false);
                           counter_vmax_r++;
                        }
#endif // end  #if defined _TOP_RANDOM_ || defined (_BOTTOM_RANDOM_)
#endif // end ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#endif //  end defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#endif // end _USE_STATISTICAL_ASSIGNMENT_
#ifndef _BIN_ACCUMULATE_
                      }
#endif
    }
  this->So.DONE();
#ifdef _USE_GNUPLOT_SREL_
//      this->plot_scaling_relation_assignment("MASS_DENS");
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  REF_DEN_FIELD.clear();
  REF_DEN_FIELD.shrink_to_fit();
#endif
#ifdef _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
  mean_ldm/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  mean_hmass/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  m_dm/=static_cast<real_prec>(this->tracer_ref._NOBJS());
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
  So.message_screen("Correlation log(M) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Correlation log(Vmax) - DM density =", sqrt(fabs(m_dm-mean_hmass*mean_ldm)));
#endif
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  mean_ntr/=static_cast<real_prec>(this->tracer_ref._NOBJS());
  m_ntr/=static_cast<real_prec>(this->tracer_ref._NOBJS());
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-mass - N_tracers in cells  =", sqrt(fabs(m_ntr-mean_hmass*mean_ntr)));
#endif
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("Halo-Vmax - N_tracers in cells =", sqrt(fabs(m_ntr-mean_hmass*mean_ldm)));
  std::cout<<endl;
#endif
#endif
#endif
#endif // end for _GET_HALO_MASS_PROPERTIES_CORRELATIONS_
#ifndef _USE_STATISTICAL_ASSIGNMENT_
  real_prec aux_n=-100;
  real_prec aux_m=LARGE_NUMBER;
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_bins[ih].tracer_properties.size()), aux_n);
       aux_n=aux_m;
  }
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax > Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax < Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif // closes _BOTTOM_RANDOM_
#endif// closes   _FULL_VERBOSE_

  aux_n=-100;
  aux_m=0;
  for(int ih=0;ih<LENGHT_AB ;++ih)
     {
       aux_m=max(static_cast<real_prec>(this->dm_properties_for_randoms_bins[ih].tracer_properties.size()), aux_n);
       aux_n=aux_m;
      }
#ifdef _FULL_VERBOSE_
#ifdef _BOTTOM_RANDOM_
  So.message_screen("Maximum number of tracers with Vmax <= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#else
  So.message_screen("Maximum number of tracers with Vmax >= Vmax_threshold in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
#endif // closes _BOTTOM_RANDOM_
#endif// closes   _FULL_VERBOSE_
#else // else for #ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
  So.message_screen("Maximum number of tracers in a bin of DM properties {Theta} = ", static_cast<real_prec>(aux_m));
  std::cout<<endl;
#endif
#endif // closes  _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  ULONG aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in theta-containers:");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:aux_a)
#endif
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_a+=this->dm_properties_bins[ih].used_property.size();// using masses_bin_properties leads to the same result.
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  So.message_screen("Number of tracers_dm in theta-containers =", aux_a);
  So.message_screen("check =", counter_vmax);
  So.message_screen("Expected =", this->tracer.Ntracers_dm);
#else
#ifdef _FULL_VERBOSE_
  So.message_screen("Number of tracers in theta-containers in =", aux_a);
#endif
#endif // closes _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  //Here we sum the tracers with Vmax below the threshold identified below which random tracers will take their vmax
  ULONG aux_b=0;

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:aux_b)
#endif
  for(ULONG ih=0;ih<LENGHT_AB ;++ih)
    aux_b+=this->dm_properties_for_randoms_bins[ih].used_property.size();
  So.message_screen("Number of tracers_random in theta-containers =", aux_b);
  So.message_screen("check =", counter_vmax_r);
  So.message_screen("Expected =", this->tracer.Ntracers_ran);
  aux_a+=aux_b;
#endif // closes _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  try{
      check_if_positive_number<long>(this->tracer_ref._NOBJS()-aux_a);
     }
  catch( const std::domain_error& e )
  {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in vec<dm_props> dm_properties is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cerr << e.what() << endl;
  }
#else // else for _USE_STATISTICAL_ASSIGNMENT_
  ULONG baux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in abundance:");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:baux_a)
#endif
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     baux_a+=ABUNDANCE[ih];
  if(baux_a<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cout<<baux_a<<"  "<<this->tracer_ref._NOBJS()<<endl;
      exit(0);
    }
  So.DONE();
  //* Final check to verify that all tracers have been counted in
  ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->ABUNDANCE));
  if(ncells_used<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
      std::cout<<"Number of objects allocated = "<<ncells_used<<"  Number expected = "<<this->tracer_ref._NOBJS()<<endl;
    }
    else
      So.message_screen("Test on ABUNDANCE at line ", __LINE__, " passed. ");
#ifdef _FULL_VERBOSE_
  So.message_screen("Normalizing to get joint probability distribution");
#endif
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
  long aux_a,aux_b;
#ifdef _USE_OMP_
#pragma omp parallel for private(aux_a,aux_b)
#endif
  for(ULONG idm = 0; idm < LENGHT_AB; ++idm)
  {
    aux_a=-10000;
    for(ULONG tr_j = 0; tr_j < N_x; ++tr_j)
      {
        aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_j,idm,LENGHT_AB)]));
        aux_a=aux_b;
     }
    for(ULONG tr_j = 0; tr_j < N_x; ++tr_j)
     {
       ULONG indexa=index_2d(tr_j,idm,LENGHT_AB);
       long AUX=this->ABUNDANCE[indexa];
       this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
//           if(this->ABUNDANCE[indexa]>0) cout<<aux_b<<"  "<<AUX<<"  "<<tr_j<<" "<<idm<<endl;
        }
   }
  this->So.DONE();
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking conditional prob distribution:");
#endif
    real_prec aux_h=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:aux_h)
#endif
    for(ULONG ih=0;ih<this->ABUNDANCE_normalized.size();++ih)
      aux_h+=this->ABUNDANCE_normalized[ih];

    if(aux_h<=0)
      {
        So.message_screen("Negative/zero value not tollerated : ", aux_h);
        So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts ill-defined. CosmicAtlas stops here.");
        exit(0);
      }
    else
        So.message_screen("Test on ABUNDANCE_normalized at line ", __LINE__, " passed. ");
     this->So.DONE();
#ifdef _VERBOSE_FREEMEM_
    So.message_screen("Freeing memmory");
#endif
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    So.DONE();

#endif  // end for #ifndef _USE_STATISTICAL_ASSIGNMENT_
}
//end class member function get_scaling_relations_primary_property()
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_scaling_relations_secondary_property(string h_property)
//* This function is explicitely meant to be applid to assign masses/vmax after having assigned vmax/masses  *//
{
  So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  if(h_property==_MASS_)// Mvir is secondary if Vmax is primary observavble
   {
    So.message_screen("***Measuring conditional MASS Function from reference tracer catalog*****");
    So.message_screen("****as a function of Vmax and properties of the reference DM density field*****");
   }
  else if(h_property==_VMAX_) // Vmax is secondary if Mvir is primary observavble
   {
      So.message_screen("***Measuring conditional Vmax Function from reference tracer catalog*****");
      So.message_screen("****as a function of Mvir and properties of the reference DM density field*****");
    }
  else if(h_property==_CONCENTRATION_)
   {
      So.message_screen("***Measuring conditional Cvir Function from reference tracer catalog*****");
      So.message_screen("****as a function of M and Vmax *****");
    }
  else if(h_property==_RS_)
   {
      So.message_screen("***Measuring conditional Rs Function from reference tracer catalog*****");
      So.message_screen("****as a function of M and Vmax *****");
    }
  else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
   {
     So.message_screen("***Measuring conditional Spin Function from reference tracer catalog*****");
     So.message_screen("****as a function of M and Vmax *****");
   }
#endif
  ULONG extra_bins=1;
  if(h_property==_MASS_ || h_property==_VMAX_)
    extra_bins=EXTRA_DM_BINS;
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  So.message_warning("See documentation");
  // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
  // The histogram has the form P(Y|X1,X2,X3,X4,X5)
  ULONG Ntot=1; // Ntot is the total number of bins to use in the P(Y|Xi)
  //-----------------------
  // For X1 (see table in documentation):
  Ntot*=this->params._NPROPbins_bam();
  //-----------------------
  // For X2: 
  if(h_property==_MASS_ || h_property==_VMAX_ || h_property==_CONCENTRATION_ || h_property==_RS_)
    {    // for mass, vmax, cvir, we  use dm in IX_2
      if(EXTRA_DM_BINS>1)
      {
        real_prec oldNX=params._NX();
        this->params.set_NX(oldNX*extra_bins);
        this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->params._NX()); // redefine DELTA
        this->s_deltas.prop1=this->DELTAX;
      }
      Ntot*=this->params._NX();
   }
  else if( h_property==_SPIN_ || h_property==_SPIN_BULLOCK_){   // For spin, we use cvir 
    Ntot*=this->params._NPROPbins_bam();
  }
  //-----------------------
  //For X3:
  ULONG n_cwt=1;
#ifdef _USE_CWC_mass_
   n_cwt= static_cast<ULONG>(this->params._n_cwt());
#endif
  Ntot*=n_cwt;
  //-----------------------
  // For X4:
#ifdef _USE_MACH_NUMBER_
  Ntot*=N_BINS_MACH;
#endif
  //-----------------------
  // For X5:
#ifdef _USE_LOCAL_OVERDENSITY_
  Ntot*=N_BINS_LO;
#endif
  //-----------------------
  // For X6:
#ifdef _USE_TIDAL_ANISOTROPY_SEC_PROP_
    Ntot*=N_BINS_TA;
#endif
  //-----------------------
  ULONG Ntot_partial=Ntot;
  //-----------------------
  Ntot*=this->params._NPROPbins_bam(); // For the variable Y, done separately because elow we need explicitely to split (for noimalization)
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(Ntot,0);
  this->get_new_min_max_properties();
#ifdef _USE_TOTAL_MASS_IN_CELL_
  vector<real_prec> REF_MASS_FIELD( this->params._NGRID(),0);
  this->tracer_ref.get_density_field_grid(_MASS_, REF_MASS_FIELD);
#endif

warning:  
  So.message_warning("Using Vmax instead of M to assign Cvir and Rs");

#ifdef _USE_OMP_
#pragma omp parallel
  {
    vector<real_prec>abundance_parallel(this->ABUNDANCE.size(),0);
#pragma omp for nowait
#endif
    for(ULONG ig = 0; ig< this->tracer_ref._NOBJS() ; ++ig)
     {
       ULONG ID=this->tracer_ref.Halo[ig].GridID;       // Get ID of this tracer
      // -----------------------------------------------------------------------------------------
      // P(Y|X1,X2,X3,X4,X5) h_property refers to Y
       ULONG I_Y=0;//This index takes the place of the variable to assign (e.g. Mvir, Rs), so here it is read form ref
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
       if(h_property==_MASS_)// if we want to assign mass, learn from it
        I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),this->params._NPROPbins_bam(), (this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NPROPbins_bam()),this->bin_accumulate_borders);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
       if(h_property==_VMAX_)
        I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].vmax),log10(this->params._VMAXmin()),this->params._NPROPbins_bam(),this->tracer_ref.logdeltaVMAX,this->bin_accumulate_borders);
#endif
#ifdef  _USE_RS_AS_DERIVED_OBSERVABLE_  
        else if(h_property==_RS_ )
         I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].rs),log10(this->params._RSmin()),this->params._NPROPbins_bam(),this->tracer_ref.logdeltaRS,this->bin_accumulate_borders);
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
        else if( h_property == _CONCENTRATION_)
         I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].concentration),log10(this->params._CONCENTRATIONmin()),this->params._NPROPbins_bam(),log10(this->params._CONCENTRATIONmax()/this->params._CONCENTRATIONmin())/static_cast<double>(this->params._NPROPbins_bam()),this->bin_accumulate_borders);
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
        else if(h_property==_SPIN_)
          I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].spin), log10(this->params._SPINmin()),this->params._NPROPbins_bam(),this->tracer_ref.logdeltaSPIN,this->bin_accumulate_borders);
        else if( h_property==_SPIN_BULLOCK_)
          I_Y  = get_bin(log10(this->tracer_ref.Halo[ig].spin_bullock), log10(this->params._SPINmin()),this->params._NPROPbins_bam(),this->tracer_ref.logdeltaSPIN,this->bin_accumulate_borders);
#endif
      ULONG I_X1=0;
#if defined  _USE_RS_AS_DERIVED_OBSERVABLE_ || defined _USE_SPIN_AS_DERIVED_OBSERVABLE_  || defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
      // IF we neeed to assign RS, CVIR, or Spin, put Mass in the X1 slot
      if(h_property==_RS_ || h_property == _CONCENTRATION_ || h_property==_SPIN_  ||  h_property==_SPIN_BULLOCK_ ) // Using DELTA only allowed so far when assigning mass. FOr Rs and Spin we use so far P(Rs|M,Vmax) or P(S|M,Vmax)
        I_X1 =  get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),this->params._NPROPbins_bam(), (this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NPROPbins_bam())  ,this->bin_accumulate_borders);
        // Use this if we want to use Vmax (even if primary) to assign seconday, passing over Mvir assigned
//        I_X1 =  get_bin(log10(this->tracer_ref.Halo[ig].vmax),log10(this->params._VMAXmin()),this->params._NPROPbins_bam(), (log10(this->params._VMAXmax())-log10(this->params._VMAXmin()))/static_cast<double>(this->params._NPROPbins_bam())  ,this->bin_accumulate_borders);
#endif
      else{ // Unless we want to assign M(or Vmax) in which case we put Vmax (or M)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_  // Here this means that Vmax is already available and we can learn from it
        if(h_property==_MASS_)
         I_X1=get_bin(log10(this->tracer_ref.Halo[ig].vmax), log10(this->params._VMAXmin()),this->params._NPROPbins_bam(),log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(this->params._NPROPbins_bam()), true);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ // If mvir is primary  we can now assign vmax from P(Vmax|M...)
       I_X1=get_bin(log10(this->tracer_ref.Halo[ig].mass), this->params._LOGMASSmin(),this->params._NPROPbins_bam(),(this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NPROPbins_bam()), true);
#endif
      }
   // Get bin of DM ovedensity
   // We use this bin as I_X ->DM for P(M|DM,Vmax) , or I_X -> Mass for P(Rs (or spin)|M,Vmax). Vmax will take the index I_CV2. If add_dm_density undef, I_X=0 for P(M|Vmax)
      ULONG I_X2=0; 
      ULONG NbinsX2=0;
      if(h_property==_MASS_ || h_property==_VMAX_  || h_property==_CONCENTRATION_ || h_property==_RS_) // Using DELTA only allowed so far when assigning mass (or vmax) and concentration. FOr Rs we use so far P(Rs|M,Vmax)
        {
          real_prec xdm  = this->delta_dm_aux_mem[ID]; // this and its I-or T-web  was computed in member function  get_scaling_relations_primary_property()
          I_X2  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
          NbinsX2=this->params._NX();
        }
      else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_) // If we want to assign spin, we learn from concentration (already assigned)
        {
#if defined  _USE_RS_AS_DERIVED_OBSERVABLE_  // Recal that for assignment, Rs and Cvir cannot coexist
          I_X2  = get_bin(log10(this->tracer_ref.Halo[ig].rs),log10(this->params._RSmin()),this->params._NPROPbins_bam(),this->tracer_ref.logdeltaRS,this->bin_accumulate_borders);
#elif defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
          I_X2  = get_bin(log10(this->tracer_ref.Halo[ig].concentration),log10(this->params._CONCENTRATIONmin()),this->params._NPROPbins_bam(),log10(this->params._CONCENTRATIONmax()/this->params._CONCENTRATIONmin())/static_cast<double>(this->params._NPROPbins_bam()),this->bin_accumulate_borders);
#endif
          NbinsX2=this->params._NPROPbins_bam();
       }
    // If assignment is for secondary (Mas has been already assigned), Let IX1 takes the place of halo mass
      ULONG I_X3=0;
#ifdef _USE_CWC_mass_
       I_X3=this->cwclass_ref.get_Tclassification(ID);
#endif
      ULONG I_X4=0;
#ifdef _USE_MACH_NUMBER_
      I_X4= get_bin(this->tracer_ref.Halo[ig].mach_number,this->tracer_ref._min_mach(), N_BINS_MACH, (this->tracer_ref._max_mach()-this->tracer_ref._min_mach())/static_cast<real_prec>(N_BINS_MACH), this->bin_accumulate_borders);
#endif
      ULONG I_X5=0;
#ifdef _USE_LOCAL_OVERDENSITY_
      I_X5= get_bin(this->tracer_ref.Halo[ig].local_overdensity, this->tracer_ref._min_local_overdensity(), N_BINS_LO,  (this->tracer_ref._max_local_overdensity()-this->tracer_ref._min_local_overdensity())/static_cast<real_prec>(N_BINS_LO) ,this->bin_accumulate_borders);
#endif
      ULONG index_dm=0;
      ULONG I_X6=0;     
#ifdef _USE_TIDAL_ANISOTROPY_SEC_PROP_ 
      I_X6= get_bin(this->tracer_ref.Halo[ig].tidal_anisotropy, this->tracer_ref._min_tidal_anisotropy(), N_BINS_TA,  (this->tracer_ref._max_tidal_anisotropy()-this->tracer_ref._min_tidal_anisotropy())/static_cast<real_prec>(N_BINS_TA) ,this->bin_accumulate_borders);
      index_dm=index_6d(I_X1,I_X2,I_X3,I_X4,I_X5,I_X6, NbinsX2,n_cwt,N_BINS_MACH,N_BINS_LO,N_BINS_TA);
#else
      index_dm=index_5d(I_X1,I_X2,I_X3,I_X4,I_X5,NbinsX2,n_cwt,N_BINS_MACH,N_BINS_LO);
#endif
#ifndef _BIN_ACCUMULATE_
       if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
          {
#endif
            ULONG index_prop_b=index_2d(I_Y,index_dm,Ntot_partial); // I_Y is the bin of mass
#ifdef _USE_OMP_
            abundance_parallel[index_prop_b]++;
#endif
#ifndef _BIN_ACCUMULATE_
       }
#endif
    }
#ifdef _USE_OMP_
  {
#pragma omp critical
    for(ULONG i=0;i <this->ABUNDANCE.size();++i)
       this->ABUNDANCE[i]+=abundance_parallel[i];
   }
}// close parallel region
#endif
  this->So.DONE();
  ULONG aux_a=0;
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking number of tracers in abundance:");
#endif
#pragma omp parallel for reduction(+:aux_a)
  for(ULONG ih=0;ih<this->ABUNDANCE.size() ;++ih)
     aux_a+=ABUNDANCE[ih];
  if(aux_a<this->tracer_ref._NOBJS())
    {
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of tracers counted in ABUNDANCE(Y,X) is smaller than the input. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCode stops here");
      std::cout<<aux_a<<"  "<<this->tracer_ref._NOBJS()<<endl;
      exit(0);
    }
    else
  So.message_screen("Test on Abundance in line", __LINE__, ": passed");
// Free memmory only in the last property (spin) has been assigned
#if defined (_USE_RS_AS_DERIVED_OBSERVABLE_) || defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_ && !defined _USE_SPIN_AS_DERIVED_OBSERVABLE_
      this->delta_dm_aux_mem.clear();delta_dm_aux_mem.shrink_to_fit();
#elif defined  _USE_SPIN_AS_DERIVED_OBSERVABLE_ 
   if(h_property==_SPIN_ || h_property == _SPIN_BULLOCK_)
      this->delta_dm_aux_mem.clear();delta_dm_aux_mem.shrink_to_fit();
#endif
/*  // These lines are commented since there is no strickt rule to free memmory here for the container this->tracer_ref.Halo, as it willbe used below ehen comparing/writting bias to output files for comparison
{
#ifndef _USE_BIAS_OBJECT_TO_OBJECT_  // if bias is not to be used at all, we can clar memory here
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  if(h_property == _SPIN_ || h_property == _SPIN_BULLOCK_) // The tracer info canot be cleared until we have assigned Rs or the last variable
   {
    So.message_screen("**Freeing memmory from tracer_ref and auxiliary containers, line", __LINE__);
    this->tracer_ref.Halo.clear();
    this->tracer_ref.Halo.shrink_to_fit();
    this->So.DONE();
  }
#endif
#endif
}
*/
 this->ABUNDANCE_normalized.clear();
 this->ABUNDANCE_normalized.shrink_to_fit();
 this->ABUNDANCE_normalized.resize(Ntot,0);
#ifdef _FULL_VERBOSE_
  So.message_screen("Normalizing joint probability distribution");
#endif
#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
 {
#pragma omp parallel for
#endif
    for(ULONG idm = 0; idm < Ntot_partial; ++idm)
     {
     long aux_a=-10000;
     long aux_b=0;
     for(ULONG tr_j = 0; tr_j < this->params._NPROPbins_bam()  ; ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,Ntot_partial);
         long AUX=this->ABUNDANCE[indexa];
         aux_b=max(aux_a, AUX);
         aux_a=aux_b;
       }
     for(ULONG tr_j = 0; tr_j < this->params._NPROPbins_bam()  ; ++tr_j)
       {
         ULONG indexa=index_2d(tr_j,idm,Ntot_partial);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
   }
#ifdef _USE_OMP_
}
#endif
#ifdef _FULL_VERBOSE_
  So.message_screen("Checking conditional prob distribution:");
#endif
    real_prec aux_h=0;
    real_prec aux_j=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:aux_h,aux_j)
#endif
    for(ULONG ih=0;ih<this->ABUNDANCE_normalized.size();++ih)
      {
        aux_h+=this->ABUNDANCE_normalized[ih];
        if(this->ABUNDANCE_normalized[ih]<0)
            aux_j++;
    }
    if(aux_h<=0 || aux_j>0)
      {
        So.message_screen("Negative/zero value not tollerated : ", aux_h);
        So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts ill-defined. CosmicCode stops here.");
        exit(0);
      }
    else
        So.message_screen("Test on ABUNDANCE_normalized at line ", __LINE__, " passed. ");
    this->So.DONE();
 So.DONE();
 this->ABUNDANCE.clear();
 this->ABUNDANCE.shrink_to_fit();
}
#endif  // end of ifdef MOCK_MODE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MOCK_MODE
 void BiasMT::get_BIAS(string property)
 {
    So.enter(__PRETTY_FUNCTION__);
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
   real_prec num_in_log_x = true==this->params._Convert_Density_to_Delta_X() ? NUM_IN_LOG: 0.;
   real_prec num_in_log_y = true==this->params._Convert_Density_to_Delta_Y() ? NUM_IN_LOG: 0.;
#ifdef _DO_BiasMT_CALIBRATION_
#ifdef _FULL_VERBOSE_
   std::cout<<endl;
    So.message_screen("Measuring Bias: ");
    std::cout<<endl;
#endif
   if(this->params._Scale_Y()=="log")
     {
       if(this->iteration==0) // Do this only in the first step, sicne in the iteration process does not change the tracer density field.
        {
#ifdef _FULL_VERBOSE_
       So.message_screen("Transforming delta_Y->log10(num_in_log+delta_Y_ref). Line ", __LINE__);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i = 0; i <  this->params._NGRID(); ++i)
          this->delta_Y[i] = this->delta_Y[i]<-1 ? 0: log10(num_in_log_y+static_cast<real_prec>(this->delta_Y[i]));
        }
      }
   if(this->params._Scale_X()=="log")
     {
#ifdef _FULL_VERBOSE_
       So.message_screen("Transforming delta_X->log10(2+delta_ref). Line ", __LINE__);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i = 0;i <  this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA). tHIS IS DONE SINCE THE KONVOLUITION ALWAYS GOVES DELTAS, WE HAVE TO TRANSFORMTO LOG(1+DELTA)
         this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(num_in_log_x + static_cast<real_prec>(this->delta_X[i]));
      So.DONE();
     }
#endif  //end if DO_BiasMT_CALIBRATION
   this->new_nbins_x = this->params._NX();
 // ******************************************************************************
   // For the mocks we convert always X to delta, so no question mark here. xmin, xmax are not updated
   //  this->new_nbins_x = this->params._NX();
   if(_COUNTS_==property)
     {
       this->Ymin=NMIN_Y_ONECELL;
       this->Ymax=this->nmax_y_onecell;
       this->new_nbins_y = this->Ymax+1; // One bin per particle, plus 0
     }
   else
     {
       this->new_nbins_y = this->params._NY();
       if(this->params._Scale_Y()=="log")
        {
           this->Ymin=this->params._ldelta_Y_min();
           this->Ymax=this->params._ldelta_Y_max();
        }else
        {
           this->Ymin=this->params._delta_Y_min();
           this->Ymax=this->params._delta_Y_max();
         }
     }
   this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->params._NX());
   this->DELTAY=(this->Ymax-this->Ymin)/static_cast<real_prec>(this->new_nbins_y);
   // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
   ULONG NX_NEW=1;
#ifdef _USE_DM_IN_BiasMT_
   NX_NEW=this->params._NX();
#endif
   ULONG LENGHT_BIAS_NCOUNT=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->params._n_sknot_massbin() * this->params._n_cwt() * this->params._n_vknot_massbin() * this->params._n_cwv()* NX_NEW * this->new_nbins_y;
   ULONG LENGHT_BIAS_NCOUNT_aux=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3* N_C_BIN1 * N_C_BIN2* N_C_BIN3* this->params._n_sknot_massbin() * this->params._n_cwt() * this->params._n_vknot_massbin() * this->params._n_cwv()* NX_NEW;
#ifdef _USE_GNUPLOT_
   ULONG LENGHT_BIAS_LOCAL_DENSITY=NX_NEW * this->new_nbins_y;
#endif
   int count_arrays=1;
#ifdef _USE_MASS_TRACERS_
   count_arrays++;
#endif
#if !defined (_ASSIGN_PROPERTIES_TO_REFERENCE_)
#ifdef _DO_BiasMT_CALIBRATION_
#ifdef _FULL_VERBOSE_
   So.message_screen("Allocating", LENGHT_BIAS_NCOUNT*count_arrays*(sizeof(real_prec)+sizeof(ULONG))/(1e9), "Gb for bias");
#endif
#endif
// This section will be done in parallel only if we have more than two threads
  this->BIAS_NCOUNTS.clear();
  this->BIAS_NCOUNTS.shrink_to_fit();
  this->BIAS_NCOUNTS.resize(LENGHT_BIAS_NCOUNT, 0);
  this->BIAS_NCOUNTS_normalized.clear();
  this->BIAS_NCOUNTS_normalized.shrink_to_fit();
  this->BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);
#endif
#ifdef _SHOW_BIAS_
  this->BIAS_LOCAL_DM.clear();
  this->BIAS_LOCAL_DM.shrink_to_fit();
  this->BIAS_LOCAL_DM.resize(LENGHT_BIAS_NCOUNT, 0);
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
#ifndef _ONLY_POST_PROC_ // when post_proc is defiened, we do not need the bias for the numbr counts is already generated
#ifdef _FULL_VERBOSE_
   So.message_screen("Reading BIAS from calibration");
#endif
#ifndef _USE_TWO_REFS_MOCKS_
#ifndef __ASSIGN_PROPERTIES_TO_REFERENCE_
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
    ULONG ncounts_aux = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
#ifdef _FULL_VERBOSE_
    So.message_screen("Number of cells in BIAS = ", ncounts_aux);
#endif
    if(ncounts_aux!= this->params._NGRID())
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Missing cells in get_mock_grid(). Perhaps density range is not wide enough to contain all cells");
   else
    So.DONE();
   real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
   if(Kd<=0)
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Conditional probablity ill-defined.");

#endif  // end of ifndef _use_two:refs
#endif  // end of ifndef _MASS_ASSIGNMENT_TO_REFERENCE_
#endif  // endif for  ifndef _ONLY_POST_PROC_
#else  // else for ifdef _GET_BiasMT_REALIZATIONS_. The following lines are then allowed for CALIBRATION
   if(true==this->use_iteration_ini) //NOT WORKING YET
     {
       this->File.read_array(this->params._Output_directory()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
       this->File.read_array(this->params._Output_directory()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
       this->use_iteration_ini=false; //set it to false such that we will not read this again,
     }
   else
     {
       this->get_new_min_max_properties();
       // The outputs of this method (private variables) are also used in get_mock
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG ig = 0; ig<  this->params._NGRID() ; ++ig)
         {
           // Get number counts or the property tracer. If displacement, the definition of Delta is on line 12277 of Bam.cpp (bamrunner)
           real_prec halo_prop = this->delta_Y[ig];
           int I_Y=0;
           if(_COUNTS_== property)
             I_Y=static_cast<int>(halo_prop);
           else
             I_Y=get_bin(halo_prop,this->s_mins.prop0, this->new_nbins_y, this->s_deltas.prop0,this->bin_accumulate_borders);
#ifdef _USE_MASS_FIELD_   //get bin of tracer mass
           real_prec halo_mass_prop=this->delta_Y_MASS[ig];
           int I_Y_MASS=get_bin(halo_mass_prop,this->s_mins.prop0_mass, this->params._NY_MASS(), this->s_deltas.prop0_mass,this->bin_accumulate_borders);
#endif
           ULONG I_X=0;
#ifdef _USE_DM_IN_BiasMT_
           real_prec xdm    = this->delta_X[ig];
           I_X  = ((this->params._iMAS_X() == 0  && false==this->params._Convert_Density_to_Delta_X()) ? static_cast<int>(this->delta_X[ig]) : get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders));
#endif
           ULONG I_CWT=0;
#ifdef _USE_CWC_
           I_CWT=this->cwclass.get_Tclassification(ig);
#endif
           ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
           I_MK= (this->cwclass.cwt_used[this->cwclass.get_Tclassification(ig)]== I_KNOT ? this->cwclass.SKNOT_M_info[ig]: 0);
#endif
           ULONG I_CWV=0;
#ifdef _USE_CWC_V_
           I_CWV=this->cwclass.get_Vclassification(ig);
#endif
           ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
           I_VK= (this->cwclass.cwv_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[ig]: 0);
#endif
           ULONG I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
           real_prec C1 = this->cwclass.Invariant_TF_II[ig];
           I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
           real_prec C1 = this->cwclass.DELTA2[ig];
           I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
           ULONG I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
           real_prec C2 = this->cwclass.Invariant_TF_III[ig];
           I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, this->s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
           real_prec C2 = this->cwclass.DELTA3[ig];
           I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif
           ULONG I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
           real_prec C3 = this->cwclass.Invariant_TF_IV[ig];
           I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, this->s_deltas.prop6,this->bin_accumulate_borders);
#elif defined (_USE_TIDAL_ANISOTROPY_)
#ifdef _USE_ANISOTROPY_AT_TRACER_POSITION_
      real_prec C3  =linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->tracer_ref.Halo[ig].coord1,this->tracer_ref.Halo[ig].coord2,this->tracer_ref.Halo[ig].coord3,this->cwclass.Tidal_Anisotropy);
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#else
      real_prec C3 = this->cwclass.Tidal_Anisotropy[ID];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
#elif defined _USE_S2_
           real_prec C3 = this->cwclass.S2[ig];             // s²
           I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
           ULONG I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
           real_prec CV1 = this->cwclass.Invariant_VS_I[ig];
           I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, this->s_deltas.prop7, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
           real_prec CV1 = this->cwclass.N2D[ig];      // Nabla² ð
           I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
 #elif defined (_USE_INVARIANT_PWEB_I_)
         real_prec CV1 = this->cwclass.Invariant_TF_I[ig]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#endif
           ULONG I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
           real_prec CV2 = this->cwclass.Invariant_VS_II[ig];
           I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2,this->s_deltas.prop8, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
           real_prec CV2 = this->cwclass.S2DELTA[ig];         // s²ð
           I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_II_)
           real_prec CV2 = this->cwclass.Invariant_TF_II[ig]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
           I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif
           ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
           real_prec CV3 = this->cwclass.Invariant_VS_III[ig];
           I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3,this->s_deltas.prop9, this->bin_accumulate_borders);
#elif defined _USE_S3_
           real_prec CV3 = this->cwclass.S3[ig];                                   // s³
           I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_III_)
           real_prec CV3 = this->cwclass.Invariant_TF_III[ig]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
           I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif
#ifndef _BIN_ACCUMULATE_
           if(halo_prop <=this->s_maxs.prop0 && halo_prop >=this->s_mins.prop0)
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_DM_IN_BiasMT_
             if(xdm <=this->s_maxs.prop1 && xdm >=this->s_mins.prop1)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
               if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
                 if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_)  || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
                   if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
                     if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif

#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
                       if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
                         if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
#endif
                           {
                             ULONG index_dm=index_11d(I_X, I_CWT, I_MK, I_CWV,I_VK, I_C1, I_C2, I_C3, I_CV1, I_CV2, I_CV3, static_cast<ULONG>(this->params._n_cwt()), static_cast<ULONG>(this->params._n_sknot_massbin()), static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()),N_C_BIN1,N_C_BIN2,N_C_BIN3, N_CV_BIN1, N_CV_BIN2,N_CV_BIN3);
                             ULONG Index=index_2d(I_Y, index_dm, LENGHT_BIAS_NCOUNT_aux);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                             this->BIAS_NCOUNTS[Index]++; // Add the number of cells that fall in this bin of Theta Index
#ifdef _SHOW_BIAS_
                             this->BIAS_LOCAL_DM[index_2d(I_Y, I_X, this->params._NX())]++; // Add the number of cells that fall in this bin of Theta Index
#endif
                        }
         }
       this->So.DONE();
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
    if(this->iteration==this->params._Iteration_Kernel_average())
    {
      this->AVERAGE_BIAS_NCOUNTS.clear();
      this->AVERAGE_BIAS_NCOUNTS.shrink_to_fit();
      this->AVERAGE_BIAS_NCOUNTS.resize(LENGHT_BIAS_NCOUNT, 0);
    }
    if(this->iteration>=this->params._Iteration_Kernel_average())
    {
        for(ULONG i=0;i<this->AVERAGE_BIAS_NCOUNTS.size();++i)
          this->AVERAGE_BIAS_NCOUNTS[i]+=this->BIAS_NCOUNTS[i]/static_cast<real_prec>(this->params._N_iterations_Kernel()-this->params._Iteration_Kernel_average());
    }
#endif
// **********************************************************************************************************
#ifdef _SHOW_BIAS_    // Gnuplot section to plot the bias as a funciton of the loca density
/*
    vector<vector<real_prec>> image;
   for(int j=0; j<this->params._NX(); j++) {
        vector<real_prec> row;
        for(int i=0; i<new_nbins_y; i++)
            row.push_back(this->BIAS_LOCAL_DM[index_2d(i,j,this->params._NX())]);
        image.push_back(row);
    }
*/
   string data_aux="data_aux.txt";
   ofstream datasal; datasal.open(data_aux.c_str());
   for(int v=0; v<new_nbins_y; v++)
        for(int u=0; u<this->params._NX(); u++)
            if(u<this->params._NX()-1)
                datasal<<log10(this->BIAS_LOCAL_DM[index_2d(v, u, this->params._NX())])<<"\t";
            else
                datasal<<log10(this->BIAS_LOCAL_DM[index_2d(v, u, this->params._NX())])<<endl;
  datasal.close();
  this->gp_kernel<<"set size 0.34,0.54\n";
  this->gp_kernel<<"set origin 0.6,0.0\n";
  this->gp_kernel<<"set border linecolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_kernel<<"set xrange[0:2] \n";
  this->gp_kernel<<"set yrange[0:20] \n";  //con este range evito que estos labels salgan en otras subplots
  this->gp_kernel<<"unset log\n";
  this->gp_kernel<<"set grid\n";
  this->gp_kernel<<"set xlabel 'log(1+{/Symbol d})' font 'Times-Roman,15'\n";
  this->gp_kernel<<"set ylabel 'N_{halos}' font 'Times-Roman,15'\n";
  this->gp_kernel<<"xmax="<<this->s_maxs.prop1<<"\n";
  this->gp_kernel<<"xmin="<<this->s_mins.prop1<<"\n";
  this->gp_kernel<<"delta="<<this->s_deltas.prop1<<"\n";
  this->gp_kernel<<"unset colorbox\n";
   this->gp_kernel<<"set pm3d map\n";
  this->gp_kernel<<"pm3d interpolate 0,0 \n";
  //this->gp_kernel<<"plot" << this->gp_kernel.binFmt2d(image, "array") << " u (xmin+($1+0.5)*delta):2:3 with image"<<endl;
  // this->gp_kernel.send2d(image);
   this->gp_kernel<<"splot '"<<data_aux<<"' u (xmin+($1+0.5)*delta):2:3 matrix with image  title 'Halo Bias'"<<endl;
  this->BIAS_LOCAL_DM.clear(); this->BIAS_LOCAL_DM.shrink_to_fit();
#endif
// **********************************************************************************************************
// **********************************************************************************************************
       //* Final check done only in case we use counts-in-cells
       // This is only done when we gemnerate mocks with the same volume of the reference
#ifdef _FULL_VERBOSE_
       So.message_screen("Checking the number of grid-cells accounted for in BIAS_NCOUNTS:");
#endif
#ifndef _EXTRAPOLATE_VOLUME_
       ULONG ncells_used = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
       if(ncells_used< this->params._NGRID()){
         So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"The number of cells counted in BIAS(Y,X) is smaller than NGRID. Perhaps the range for some of the DMF properties *must* be increased to include all cells. CosmiCatlas stops here");
       }
       else
         this->So.DONE();
#endif
#ifdef _FULL_VERBOSE_
       So.message_screen("Normalizing using",NTHREADS,"threads");
#endif
       //     int nbins_y_temp = property==_COUNTS_ ?  this->new_nbins_y: this->new_nbins_y_MW;
       int nbins_y_temp = this->new_nbins_y;
#if !defined (_ASSIGN_PROPERTIES_TO_REFERENCE_) || defined (_DO_BiasMT_CALIBRATION_)
       this->BIAS_NCOUNTS_normalized.shrink_to_fit();
       this->BIAS_NCOUNTS_normalized.clear();
       this->BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);
#endif
long aux_a_full=-10000;
long aux_b_full;
ULONG aux_b_empty=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(ULONG i=0;i < LENGHT_BIAS_NCOUNT_aux;++i)// loop sobre bines de dark matter
    {
       long aux_a=-10000;
       long aux_b;
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->BIAS_NCOUNTS[inde];
           if(AUX==0)
            aux_b_empty++;
           aux_b=max(aux_a, AUX);
           aux_a=aux_b;
           aux_b_full=max(aux_a_full, AUX);
           aux_a_full=aux_b_full;
          }
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->BIAS_NCOUNTS[inde];
           this->BIAS_NCOUNTS_normalized[inde] = (aux_a<=0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
         }
          }
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
if(this->iteration==this->params._Iteration_Kernel_average())
{
  this->AVERAGE_BIAS_NCOUNTS_normalized.clear();
  this->AVERAGE_BIAS_NCOUNTS_normalized.shrink_to_fit();
  this->AVERAGE_BIAS_NCOUNTS_normalized.resize(LENGHT_BIAS_NCOUNT, 0);
  aux_a_full=-10000;
  aux_b_full;
  aux_b_empty=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(ULONG i=0;i < LENGHT_BIAS_NCOUNT_aux;++i)// loop sobre bines de dark matter
    {
       long aux_a=-10000;
       long aux_b;
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->AVERAGE_BIAS_NCOUNTS[inde];
           if(AUX==0)
            aux_b_empty++;
           aux_b=max(aux_a, AUX);
           aux_a=aux_b;
           aux_b_full=max(aux_a_full, AUX);
           aux_a_full=aux_b_full;
          }
       for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
         {
           ULONG inde=index_2d(tr_j, i, LENGHT_BIAS_NCOUNT_aux);
           long AUX=this->AVERAGE_BIAS_NCOUNTS[inde];
           this->AVERAGE_BIAS_NCOUNTS_normalized[inde] = (aux_a<=0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
         }
    }
}
#endif   // end for _WRITE_AVERAGE_KERNEL_AND_BIAS_
double effecive_nbins=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:effecive_nbins)
#endif
   for(ULONG i=0;i < this->BIAS_NCOUNTS.size();++i)// loop sobre bines de dark matter
    effecive_nbins+=this->BIAS_NCOUNTS[i];
  effecive_nbins/=static_cast<double>(aux_b_full);
#ifdef _FULL_VERBOSE_
   So.message_screen("Effective number of bins =",effecive_nbins);
   So.message_screen("Number of empty bins =",aux_b_empty);
   So.message_screen("Number of used bins =",this->BIAS_NCOUNTS.size()-aux_b_empty);
#endif
   this->So.DONE();
   real_prec aux_h=static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
   if(aux_h<=0)
     {
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Normalized Joint probability for number counts is ill-defined. CosmicAtlas stops here");
       exit(0);
     }
   aux_h=0;
#ifdef _UNDER_BIASED_
  if(0==this->iteration)
  {
      this->File.write_array(this->params._Output_directory()+"Bam_Raw_Bias", this->BIAS_NCOUNTS);
      this->File.write_array(this->params._Output_directory()+"Bam_Raw_Bias_Normalized", this->BIAS_NCOUNTS_normalized);
  }
#endif
#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->params._output_at_iteration()), std::end(this->params._output_at_iteration()), this->iteration);
//   if( (this->iteration==this->N_iterations_Kernel) || (out_it != std::end(this->params._output_at_iteration())))
   if( this->iteration==this->params._N_iterations_Kernel())  // Print arrays in the last step if the iteration procedure
     {

       this->File.write_array(this->params._Output_directory()+"Bam_Bias", this->BIAS_NCOUNTS);
       this->File.write_array(this->params._Output_directory()+"Bam_Bias_Normalized", this->BIAS_NCOUNTS_normalized);

#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
       this->File.write_array(this->params._Output_directory()+"Averaged_BiasMT_Bias", this->AVERAGE_BIAS_NCOUNTS);
       this->File.write_array(this->params._Output_directory()+"Averaged_BiasMT_Bias_Normalized", this->AVERAGE_BIAS_NCOUNTS_normalized);
#endif
       vector<real_prec>BIAS_AUX(LENGHT_BIAS_NCOUNT_aux,0);
#pragma omp parallel for collapse(2)
   for(ULONG i=0;i<LENGHT_BIAS_NCOUNT_aux;++i)
     for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
      {
       ULONG inde_l=index_2d(tr_j,i,LENGHT_BIAS_NCOUNT_aux);
#pragma omp atomic update
       BIAS_AUX[i]+=this->BIAS_NCOUNTS[inde_l];
     }
     this->File.write_array(this->params._Output_directory()+"Bam_Bias_AUX", BIAS_AUX);// print aux array from the last iteration
#ifdef _WRITE_AVERAGE_KERNEL_AND_BIAS_
   BIAS_AUX.clear();
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
   for(ULONG i=0;i<LENGHT_BIAS_NCOUNT_aux;++i)
     for(int tr_j = 0; tr_j < nbins_y_temp; ++tr_j)
      {
       ULONG inde_l=index_2d(tr_j,i,LENGHT_BIAS_NCOUNT_aux);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
       BIAS_AUX[i]+=this->AVERAGE_BIAS_NCOUNTS[inde_l];
     }
     this->File.write_array(this->params._Output_directory()+"Averaged_BiasMT_Bias_AUX", BIAS_AUX); // print aux array averaged
#endif

     BIAS_AUX.clear(); BIAS_AUX.shrink_to_fit();
    }  // closes  if( (this->iteration==this->N_iterations_Kernel) || (out_it != std::end(this->params._output_at_iteration())))
#endif

     }  // closes if(this->params._iteration_ini()>0 && true==this->use_iteration_ini) //NOT WORKING YET, mut be deprecated
#endif // endif for _GET_BiasMT_REALIZARTIONS_
 }
//end class member function
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // This function is used in mode _GET_BiasMT_REALIZATIONS_ and gets the approximated DMF
 // Action: get new alpt dm filed. Get bam kernel and convolve
 // This fuctionn doies not need to know if we use cwc insde loops, for those loope are meant for the calibration procedure.
 void BiasMT::get_new_DM_field()
 {
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_TEST_
   int NTHREADS=_NTHREADS_;   // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#else
   int NTHREADS=1;   // omp_get_max_threads();
#endif
#ifndef _ONLY_POST_PROC_
#ifdef _FULL_VERBOSE_
   std::cout<<endl;
   this->So.message_screen("*************************************************************************");
   this->So.message_screen("****Getting new DM density field in line", __LINE__, "of file",__FILE__);
   this->So.message_screen("*************************************************************************");
   std::cout<<endl;
#endif
#endif
   real_prec num_in_log_x = (true==this->params._Convert_Density_to_Delta_X() ? NUM_IN_LOG: 0.);
   string file_X;
   string file_Y;
   string file_Vx;
   string file_Vy;
   string file_Vz;
   this->delta_X.clear();
   this->delta_X.shrink_to_fit();
   this->delta_X.resize( this->params._NGRID(),0);
   this->delta_X_ini.clear();
   this->delta_X_ini.shrink_to_fit();
   this->delta_X_ini.resize( this->params._NGRID(),0);
   int n_realization=this->iteration - this->params._N_iterations_Kernel();
   real_prec nmean=0;
    ULONG NXn=600;// This is as NX but higher to make pdf
   vector<real_prec>xbaux(NXn, 0);
   vector<real_prec>pdf_in(NXn, 0);
#ifndef _USE_LPT_  //If we cannot run LPT, we need to read DM files from somewhere, so we read them from paths hard coded here
#if defined (_ASSIGN_PROPERTIES_TO_REFERENCE_) && !defined (_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_)
   this->File.read_array(this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X(), this->delta_X_ini);
#elif defined (_ASSIGN_PROPERTIES_TO_REFERENCE_) && defined (_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_)
   this->File.read_array(this->params._Input_Directory_X_new_ref()+this->params._Name_Catalog_X_new_ref(), this->delta_X_ini);
#elif defined  _ASSIGN_PROPERTIES_TO_MOCK_
   this->File.read_array(this->params._Input_Directory_X_NEW()+this->params._Name_Catalog_X_NEW(), this->delta_X_ini);
#endif
#ifdef _USE_VELOCITIES_
   file_Vx="../";
   file_Vy="../";
   file_Vz="../";
   this->Velx_X.resize( this->params._NGRID(),0);
   this->Vely_X.resize( this->params._NGRID(),0);
   this->Velz_X.resize( this->params._NGRID(),0);
   this->File.read_array(file_Vx, this->Velx_X);
   this->File.read_array(file_Vy, this->Vely_X);
   this->File.read_array(file_Vz, this->Velz_X);
#endif
   nmean=get_nobjects(this->delta_X_ini);
#ifdef _FULL_VERBOSE_
   So.message_screen("Total number of X objects (new) =", nmean);
#endif
   nmean/=static_cast<real_prec>( this->params._NGRID());
#ifdef _FULL_VERBOSE_
   So.message_screen("Average number of X objects (new) =", nmean);
   std::cout<<endl;
#endif
#ifdef _ASSIGN_PROPERTIES_TO_MOCK_  // just to save tmie while doing the test of mass assingment to tracers
#ifdef _GET_POWER_REFS_
   this->get_power_spectrum("DM_REF_NEW");
#endif
#endif
   get_overdens(this->delta_X_ini, nmean, this->delta_X_ini);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<NXn; ++i)
      xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(NXn));
#ifdef _FULL_VERBOSE_
     So.message_screen("Measuring pdf of log(NumInLogG+delta) of input NEW DM field");
#endif
     {
       pdf_in.resize(NXn, 0);
       string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
       calc_pdf("log",  this->params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
       string filex=this->params._Output_directory()+"pdf_NEW_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
#ifdef _WRITE_PDF_
       this->File.write_to_file(filex, xbaux,pdf_in);
#endif
       So.DONE();
     }
#ifdef _RANK_ORDERING_MOCK_GEN_
#ifdef _FULL_VERBOSE_
    So.message_screen("***Rank ordering to target field****");
    So.message_screen("Looking for TARGET reference field");
#endif

    this->delta_X_REF_PDF.resize( this->params._NGRID(),0);
    string file_X_ref_pdf=this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X_REF_PDF();
    this->File.read_array(file_X_ref_pdf, this->delta_X_REF_PDF); // Read target DM field for RO
#ifdef _RO_WITH_DELTA_MOCK_GEN_
    get_overdens(this->delta_X_REF_PDF, this->delta_X_REF_PDF);
#endif
#ifdef _FULL_VERBOSE_
    So.message_screen("Measuring pdf of log(1+delta) of TARGET reference field");
#endif
    this->pdf_ref.resize(NXn, 0);
    calc_pdf("log",  this->params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
    So.DONE();
#ifdef _WRITE_PDF_
     {
       string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
       string filex=this->params._Output_directory()+"pdf_TARGET_REF_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
       this->File.write_to_file(filex, xbaux,this->pdf_ref);
       So.DONE();
     }
#endif
#ifdef _FULL_VERBOSE_
       So.message_screen("Executing rank ordering from NEW DM to DM-target");
#endif
      //rankorder(ZERO, xbaux, NXn,  this->Xmax, this->Xmin, this->delta_X_ini, pdf_in, this->pdf_ref);
      swap_amp_fourier(this->params._Nft(),this->delta_X_REF_PDF,this->delta_X_ini);
      this->delta_X_REF_PDF.clear();
      this->delta_X_REF_PDF.shrink_to_fit();
      So.DONE();
#ifdef _WRITE_PDF_
       {
#ifdef _FULL_VERBOSE_
        So.message_screen("Measuring pdf of log(1+delta) of NEW DM filed Rank-ordered to TARGET");
#endif
        this->pdf_ref.resize(NXn, 0);
        string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
        calc_pdf("log",  this->params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_ini, this->pdf_ref);
        string filex=this->params._Output_directory()+"pdf_NEW_RO_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
        this->File.write_to_file(filex, xbaux,pdf_ref);
        So.DONE();
      }
#endif
#endif //end if apply rank ordering mock_gen
       // we do not need else for the output of rank ordering is again delta_X_ini
#else   //if we have to use LPT
   {
     const gsl_rng_type *rng_t;
#ifdef _USE_OMP_
     gsl_rng **gBaseRand;
#else
     gsl_rng *gBaseRand;
#endif
     gsl_rng_env_setup();
     rng_t = gsl_rng_mt19937;// gsl_rng_default;
#ifdef _USE_OMP_
     int nt=omp_get_max_threads();
     gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(nt)
#endif
     for(int i=0;i<nt;i++)
       {
         gBaseRand[i] = gsl_rng_alloc(rng_t);
         gsl_rng_set(gBaseRand[i],this->params._seed());
       }
     time_t start_LPT;
     time(&start_LPT);
     if(false==dm_already_done) //COmpute DM if it has not been done before through a call of this function
       {
#ifdef _USE_OMP_
         this->lpt.get_dm_field(gBaseRand);
#else
         this->lpt.get_dm_field(gBaseRand);
#endif
         this->dm_already_done=true;
         this->So.message_screen("LPT has created DMDF using ALPT");
 #ifdef _FULL_VERBOSE_
         std::cout<<endl;
 #endif
       }
     // LPT has written in files, so now we read them
     file_X=this->lpt._fnameDM()+".dat";
#ifdef _USE_TRACER_HR_
     file_Y=this->params._Input_Directory_Y()+this->params._Name_Catalog_Y_HR();
#else
     file_Y=this->params._Input_Directory_Y()+this->params._Name_Catalog_Y();
#endif
     this->delta_X_ini.resize( this->params._NGRID(),0);
#ifdef _USE_VELOCITIES_
     //here we have to read the vels of the newly created density field
     file_Vx=this->lpt.fnameVX+".dat";
     file_Vy=this->lpt.fnameVY+".dat";
     file_Vz=this->lpt.fnameVZ+".dat";
     this->Velx_X.resize( this->params._NGRID(),0);
     this->Vely_X.resize( this->params._NGRID(),0);
     this->Velz_X.resize( this->params._NGRID(),0);
     this->File.read_array(file_Vx, this->Velx_X);
     this->File.read_array(file_Vy, this->Vely_X);
     this->File.read_array(file_Vz, this->Velz_X);
#endif // end ifdef _USE_VELOCITIES_
     // read the just created density field
     this->File.read_array(file_X, this->delta_X_ini);
#ifdef _GET_POWER_REFS_
     this->get_power_spectrum("DM_REF_NEW");
#endif
     nmean=get_nobjects(this->delta_X_ini);
     So.message_screen("Total number of X objects (new) =", nmean);
     nmean/=static_cast<real_prec>( this->params._NGRID());
     So.message_screen("Average number of X objects (new) =", nmean);
     // Tansform input DF to Overdensity field
     get_overdens(this->delta_X_ini, this->delta_X_ini);
   } // end else uf use LPT
#endif
#ifndef _RANK_ORDERING_MOCK_GEN_
{
// If this RO was defined, the pdf of the delta_ini_new was already computed, so we avoid it here
   this->pdf_ini.resize(NXn, 0);
   string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<NXn; ++i)
     xbaux[i]=this->Xmin+static_cast<real_prec>(i+0.5)*(this->Xmax-this->Xmin)/(static_cast<real_prec>(NXn));
#ifdef _FULL_VERBOSE_
   So.message_screen("Measuring pdf of log(1+delta) DM: ");
#endif
    pdf_in.resize(NXn,0);
    calc_pdf("log", this->params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_ini, pdf_in);
   this->pdf_ini=pdf_in;
#ifdef _WRITE_PDF_
   string filex=this->params._Output_directory()+"pdf_NEW_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X())+rest_file;
   this->File.write_to_file(filex, xbaux,pdf_in);
#endif
   So.DONE();
  }
#endif
   // Now that the new field is loadad, convolve with kernel and then get the T or I WEB
#ifdef _USE_MASS_KNOTS_
   this->cwclass.SKNOT_M_info.resize( this->params._NGRID(), 0);
#endif
   // The order is important: if we do CWC inside each loop,
   // then we convolve the target DF with the updated kernel and then do the CWC
   // Otherwise one does first the CWC once and then convolves in each iteration
   this->delta_X.resize(this->delta_X_ini.size(),0);
#ifdef _KONVOLVE_TO_MOCK_
#ifdef _FULL_VERBOSE_
    So.message_screen("Reading Kernel ");
#endif
    this->Kernel.clear();
    this->Kernel.shrink_to_fit();
    this->Kernel.resize(this->NTT, 0.0);
#ifdef _USE_TWO_REFS_MOCKS_
    int n_refs=this->params._Number_of_references();
   for(int j=0; j<n_refs; ++j)  // Loop over nrefs-1 refereces. The first reference will be counted here, hence the function get_new_dm field is not used when this function is called
   {
     vector<real_prec>kernel_ghost(this->NTT,0);
     this->File.read_array(this->params._files_kernel_references(j), kernel_ghost); // FIX NAME OF GHOST DM, I have used one of the calibration
     for(ULONG i=0; i<this->NTT; ++i)
        this->Kernel[i]+=kernel_ghost[i]/static_cast<real_prec>(n_refs);
   }
  So.DONE();
#else
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Kernel.dat", this->Kernel);
#endif
    this->Konvolve(this->delta_X_ini, this->delta_X);
#else
    this->delta_X=this->delta_X_ini;
#endif
#ifdef _WRITE_DM_DENSITY_FIELD_
        this->File.write_array(this->params._Output_directory()+"DM_DELTA_convolved_realization", this->delta_X);
#endif
#ifdef _USE_CWC_
   if(this->params._n_cwt() > 1)
     {
#endif
// Open the kernel
#if defined(_USE_CWC_) || defined (_USE_CWC_mass_)
       this->cwclass.get_CWC(this->delta_X);   //get the CW info
#ifdef _USE_MASS_KNOTS_
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif //use_mass_knots
#elif !defined _USE_CWC_ || !defined _USE_CWC_mass_
#if defined (_USE_MASS_KNOTS_) || defined (_USE_IWEB_)  || defined (_USE_PWEB_) || defined (_USE_IKWEB_) || defined (_USE_TIWEB_)
       this->cwclass.get_CWC(this->delta_X);   //get the CW info
#if defined (_USE_MASS_KNOTS_)
       this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif // use_mass_knots
#endif // use_mass_knots || other XWEB model
#endif    // !use_cwc
          // If the terms in the bias expansion are to be used, then :
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_) ) && (!defined (_USE_CWC_))
       this->cwclass.get_bias_terms(this->delta_X);
#endif
#ifdef _USE_CWC_
   }
#endif // end ifdef use_cwc
#ifndef _USE_TWO_REFS_MOCKS_
#ifdef _GET_POWER_REFS_
     this->get_power_spectrum("DM_KONV"); // This has to be done *BEFORE* transfoming to LOg
#endif
#endif
     {// sec instance
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring pdf of log(1+delta) Konvolved DM New field ");
#endif
     pdf_in.resize(this->params._NX(),0);
     calc_pdf("log",  this->params._NGRID(),this->params._NX(), this->Xmax, this->Xmin, this->delta_X, pdf_in);
     So.DONE();
#ifdef _WRITE_PDF_
     string rest_file="_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
#ifdef _RANK_ORDERING_MOCK_GEN_
     string filex=this->params._Output_directory()+"pdf_NEW_RO_Konv_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X_REF_PDF())+rest_file;
#else
     string filex=this->params._Output_directory()+"pdf_NEW_Konv_"+this->params._XNAME()+"_MASX"+to_string(this->params._iMAS_X_REF_PDF())+rest_file;
#endif
     this->File.write_to_file(filex, xbaux,pdf_in);
#endif
      pdf_in.clear(); pdf_in.shrink_to_fit();
    }// sec instance
   if(this->params._Scale_X()=="log")
     {
#ifdef _FULL_VERBOSE_
       std::cout<<endl;
       // If the mass_iterative_nes is nit to be used, here we directly convert to log(numinlog+delta)
       So.message_screen("Transforming new delta->log10(2+delta)");
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i = 0;i <  this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
         this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(num_in_log_x + static_cast<real_prec>(this->delta_X[i]));
       So.DONE();
     }
   // Once the delta field is obtained, get the limits. Note that these limits mut be those of the last iteration
   // of the iterative procedure,
   // but it is not guaranteed that this is the case. So, if in _GET_BiasMT_REALIZATIONS_mde, Leave the limits fixed for the time being.
   this->get_new_min_max_properties();
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MOCK_MODE
#ifdef _GET_BiasMT_REALIZATIONS_
#define _achtung_inside_  // This check inside the mock loop whether the cond probability is >0, such that we do not fall in infinite while loops
//#define _new_app_  //define it,
//#define _assign_again_  // with this, thesame pdf is recycled to assign Nh to missing cells once the available refs are over.
#endif
 //#define _use_filled_cells_  Use this when cañlibrating from something else such as sat fraction or mass
 void BiasMT::get_mock_grid(string property)
 {
   So.enter(__PRETTY_FUNCTION__);
   bool silent=true;
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_
   So.message_screen("Function get_mock_grid is using",omp_get_max_threads(),"threads");
#endif
   std::cout<<endl;
#endif
    ULONG nx=1;
#ifdef _USE_DM_IN_BiasMT_
   // Define seed vector for each thread in the paralellized run
   nx=this->params._NX();  // Number of bins id DM
#endif
#ifdef _use_filled_cells_
   vector<bool>filled_cells( this->params._NGRID(), true);
   vector<real_prec> ncounts( this->params._NGRID(),0); // this is only used in case:_prop=2, whcih was halo mass in cells, deprecated
#endif
   string pname,fname;
   int case_prop=0; int ny=0;
   ULONG size_bias_array=0;
   if(_COUNTS_==property)
     {
       case_prop=1;
       size_bias_array=this->BIAS_NCOUNTS.size();
       ny = this->new_nbins_y;

       if(this->iteration <=this->params._N_iterations_Kernel())
         pname ="_iteration"+to_string(this->iteration);
       else
         pname= "_realization"+to_string(this->params._realization());
       fname=this->params._Output_directory()+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft());//+"_z"+to_string(this->params._redshift());
     }

   if(_DENSITY_==property)
     {
       case_prop=2;
       size_bias_array=this->BIAS_NCOUNTS.size();
       ny = this->new_nbins_y;
       fname=this->params._Output_directory()+"MOCK_TR"+pname+"_"+"MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift());
     }
   this->fnameMOCK=fname;
   // Initialize counter for those cells for which no available positions were found
   ULONG counter_orphan=0;
   ULONG counter_out_of_theta=0;
   ULONG counter_orphan_b=0;
   ULONG counter_dyn=0;
   ULONG Ntracer_dyn=0;
   ULONG Ntracer_glob=0;
   time_t start_mock;
   time(&start_mock);
   // Allocate memmory for the mock density field
   this->delta_Y_new.resize( this->params._NGRID(),0);
#ifdef _DYNAMICAL_SAMPLING_
   // Vector to allocate the Joint distribution updated after assigning Nhalos to a cell.
   vector<ULONG> X_Y_hist_dynamical(size_bias_array);
   vector<real_prec>  X_Y_hist_dynamical_normalized(size_bias_array);
#endif
   // Initialize these vectors with the original distribution in each density bin.
#if defined _USE_MASS_FIELD_
   switch(case_prop)
     {
     case(1):
#endif
#ifdef _DYNAMICAL_SAMPLING_
     X_Y_hist_dynamical=this->BIAS_NCOUNTS;
     X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
#ifdef _new_counter_
   vector<ULONG> new_X_Y_hist_dynamical(size_bias_array);
   vector<real_prec>  new_X_Y_hist_dynamical_normalized(size_bias_array);
   new_X_Y_hist_dynamical=this->BIAS_NCOUNTS;
   new_X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
#endif
#if defined _USE_SAT_FRAC_ || defined _USE_MASS_FIELD_
       break;
#endif
#ifdef _USE_MASS_FIELD_
     case(3):
       X_Y_hist_dynamical=this->BIAS_SAT_FRACTION;
       X_Y_hist_dynamical_normalized=this->BIAS_SAT_FRACTION_normalized;
       break;
#endif
#if defined _USE_MASS_FIELD_
   }
#endif
   // Vector allocating the Joint distribution normalized within each Den bin, after having assigned a value of Nhalos to a cell.
   // Initialize these vectors with the original Joint and the normalized distribution in each density bin.
#endif  // enf if def _DYNAMICAL_SAMPLING_
   ULONG lenght_bias=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_C_BIN1 * N_C_BIN2* N_C_BIN3 * this->params._n_sknot_massbin() * this->params._n_cwv() * this->params._n_vknot_massbin() * this->params._n_cwt() * nx;
   // Vector containing the total number if cells in a given density bin and CWT and KNOT mass
   this->NCELLSperDMBIN.clear();
   this->NCELLSperDMBIN.shrink_to_fit();
   this->NCELLSperDMBIN.resize(lenght_bias, 0);
   // These two loops replace the 11 loops below. tHIS MODIFICATIONB MUST BE DONE CCONSISTENTLYT EVERYWHERE WHERE BIAS IS LOADED
#ifdef _GET_BiasMT_REALIZATIONS_
   this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_AUX.dat", this->NCELLSperDMBIN);
#else
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
   for(ULONG i=0; i< lenght_bias ;++i)
       for(ULONG j=0; j< ny ; ++j)
         {
           ULONG index_h=index_2d(j,i,lenght_bias);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
           this->NCELLSperDMBIN[i]+= this->BIAS_NCOUNTS[index_h];
          }
#endif // enf ifdef get bam realizations
   {
#ifndef _EXTRAPOLATE_VOLUME_
     if(property==_COUNTS_) // We have only allowed the number counts bias to have all cells.
       {
#ifdef _FULL_VERBOSE_
         So.message_screen("Double-checking the number of grid-cells accounted for in BIAS_NCOUNTS:");
#endif
         ULONG KK = static_cast<ULONG>(get_nobjects(this->BIAS_NCOUNTS));
         if(KK!= this->params._NGRID())
           {
             So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Missing cells in get_mock_grid(). Perhaps density range is not wide enough to contain all cells", KK);
             exit(0);
            }
         else
             So.DONE();
#endif
         real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
         if(Kd<=0){
           So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Conditional probablity ill-defined. BAM stops here.", Kd);
           exit(0);
           }
         }
   }
   // Vector to allocate the number of cells in a given density bin during the mapping */
   vector<ULONG>Ncells_density_bin_new(size_bias_array , 0);
   // Vector containing the number of cells in a given density bin, updated everytime a cell has been assigend a value of Nhalos
   vector<ULONG> NCELLSperDMBIN_now(lenght_bias, 0);
  ULONG nt_partial=0;
  ULONG nt_partial_b=0;
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _FULL_VERBOSE_
   So.message_screen("**Generating new Mock (tracer) number counts.");
   std::cout<<endl;
#endif
#else
#ifdef _FULL_VERBOSE_
#ifdef _DISPLACEMENTS_
   So.message_screen("**Generating New Delta-Dispalacement  @ iteration ", this->iteration);
#else
   So.message_screen("**Generating Mock (tracer) number counts @ iteration ", this->iteration);
#endif
   std::cout<<endl;
#endif
#endif
#ifdef _new_app_
    string file_Pdf=this->params._Output_directory()+"PDF_NC_Y_REF_COUNTS_MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+".txt";
    So.message_screen("Reading PDF from reference halo catalog in file", file_Pdf);
    ifstream pdf_file; pdf_file.open(file_Pdf);
    vector<ULONG> pdf_ref(ny+1,0);
    vector<ULONG> pdf_temp(ny+1,0);
    int iaux;
    for(int i =0; i < ny +1;++i)
      pdf_file>>iaux>>pdf_ref[i];
    pdf_file.close();
#ifdef _FULL_VERBOSE_
    So.DONE();
#endif
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _USE_OMP_
    NTHREADS=1;
    omp_set_num_threads(NTHREADS);
#endif
#endif
#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_
        So.message_screen("Using",omp_get_max_threads(),"threads");
#endif
#endif
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRand;
    int jthread=0;
    gsl_rng_env_setup();
#ifdef _USE_OMP_
/*
        //initialize gsl random number generator
  const gsl_rng_type *rng_t;
  gsl_rng **gBaseRand;
  gsl_rng_env_setup();
  rng_t = gsl_rng_mt19937;// gsl_rng_default;
  int nt=NTHREADS;
  int jthread;
  gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));
#pragma omp parallel for num_threads(nt)
  for(int i=0;i<NTHREADS;i++)
    {
      gBaseRand[i] = gsl_rng_alloc(rng_t);
      gsl_rng_set(gBaseRand[i],seed);
    }
*/
   vector<ULONG>vseeds(NTHREADS,0);
   for(int i=0;i<vseeds.size();++i)
     vseeds[i]=35+static_cast<ULONG>(i)*565;
#endif
#if defined (_new_app_) || defined (_assign_again_)
   struct s_scalar_cell_info{
     ULONG theta_index;
#ifdef _new_app_
     ULONG delta_index;
#endif
     bool assigned;
     double auxvar;
   };
    vector<s_scalar_cell_info> cell_info( this->params._NGRID());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< this->params._NGRID();++i)
    {
      cell_info[i].assigned=false;
      cell_info[i].theta_index=0;
      cell_info[i].auxvar=0;
   }
#endif
#ifdef _FULL_VERBOSE_
      So.message_screen("Assigning Nh to cells");
#endif
#ifdef _USE_OMP_
#pragma omp parallel private (jthread, gBaseRand, rng_t)
   {
#endif
#ifdef _USE_OMP_
     jthread=omp_get_thread_num();
     gsl_rng_default_seed=vseeds[jthread];
#else
     gsl_rng_default_seed=35;
#endif
     rng_t = gsl_rng_mt19937;//_default;
     gBaseRand = gsl_rng_alloc (rng_t);
#ifdef _use_more_references_
     int count_refs=0;
#endif
#ifdef _assign_again_
#ifdef _use_more_references_
     ULONG NASS=0;
     for(int Irefs=0;Irefs<this->params._Number_of_references();++Irefs)
     {
#else
      do{
#endif
#endif
#ifdef _use_more_references_
   if(Irefs>0)
   {
   if(NASS>= this->params._NGRID())
     this->So.message_screen("All cells have been assigned now.");
#if defined(_FULL_VERBOSE_)
   std::cout<<endl;
   std::cout<<endl;
   So.message_screen("Using a new reference");
   So.message_screen("Number of reference to be used", this->params._Number_of_references());
   So.message_screen("Number of reference used", Irefs+1);
   std::cout<<endl;
#endif
   // This is the place to call for a second reference. Make a function (call two) which reads the bias
   // Read kernel and bias. Assign this->Kernel and this->BIAS_NCOUNTS and this->BIAS_NCOUNTS_normalized

   this->delta_X_ini.clear(); this->delta_X_ini.shrink_to_fit();this->delta_X_ini.resize( this->params._NGRID(),0);
   this->File.read_array(this->params._Input_Directory_X_NEW()+this->params._Name_Catalog_X_NEW(), this->delta_X_ini);
   real_prec nmean=get_nobjects(this->delta_X_ini);
   nmean/=static_cast<real_prec>( this->params._NGRID());
   get_overdens(this->delta_X_ini,nmean, this->delta_X_ini);
   this->delta_X.clear(); this->delta_X.shrink_to_fit();this->delta_X.resize( this->params._NGRID(),0);
   if(1==Irefs)
     this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_TWO()+"Bam_Kernel.dat", this->Kernel);// We could assume that the kernel here is very close to the first used as refenence and avoid to get again delta_X.
   else if(2==Irefs)
     this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_THREE()+"Bam_Kernel.dat", this->Kernel);
    this->Konvolve(this->delta_X_ini, this->delta_X); //Convolve with the new kernel
#ifdef _USE_CWC_
    this->cwclass.get_CWC(this->delta_X);   //get the CW info
#ifdef _USE_MASS_KNOTS_
    this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif //use_mass_knots
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_) ||  defined (_USE_IWEB_) || defined (_USE_AWEB_) || defined (_USE_PWEB_) || defined (_USE_IKWEB_) || defined (_USE_TIWEB_)
    this->cwclass.get_CWC(this->delta_X);   //get the CW info
#if defined (_USE_MASS_KNOTS_)
    this->cwclass.get_Mk_collapsing_regions(this->delta_X,nmean);  //get the MK info
#endif // use_mass_knots
#endif // use_mass_knots || other models
#endif    // !use_cwc
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_) ) && (!defined (_USE_CWC_))
    this->cwclass.get_bias_terms(this->delta_X);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i = 0;i <  this->params._NGRID() ;++i)  //TRANSFORM DELTA TO LOG10(NUM_IN_LOG + DELTA)
     this->delta_X[i] = this->delta_X[i]<-1 ?  0 : log10(NUM_IN_LOG + static_cast<real_prec>(this->delta_X[i]));
 if(1==Irefs){
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_TWO()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_TWO()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_TWO()+"Bam_Bias_AUX.dat", this->NCELLSperDMBIN);
  }
  else if(2==Irefs)
  {
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_THREE()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_THREE()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL_THREE()+"Bam_Bias_AUX.dat", this->NCELLSperDMBIN);
   }
  }
   Ncells_density_bin_new.clear();
   Ncells_density_bin_new.resize(size_bias_array,0);
   NCELLSperDMBIN_now.clear();
   NCELLSperDMBIN_now.resize(lenght_bias,0);
   X_Y_hist_dynamical=this->BIAS_NCOUNTS;
   X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
#endif // end of _use_more_references
#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan, counter_out_of_theta,Ntracer_dyn,counter_dyn)
#endif
     //Start loop over the cells
     for(ULONG i=0;i< this->params._NGRID();++i)        // First block, meant to construct halo number counts
       {
#ifdef _assign_again_
        if(false==cell_info[i].assigned)
        {
#endif
         ULONG I_X=0;
#ifdef _USE_DM_IN_BiasMT_
         real_prec dm = this->delta_X[i];
         // Get the bin in the delta dark matter (or log 1+delta) in each cell
         I_X  = static_cast<ULONG>(((this->params._iMAS_X() == 0  && false==this->params._Convert_Density_to_Delta_X()) ? static_cast<int>(this->delta_X[i]) : get_bin(dm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders)));
#endif
         // **********CWT
         ULONG I_CWT=0;
#ifdef _USE_CWC_
         I_CWT=static_cast<ULONG>(this->cwclass.get_Tclassification(i));
#endif
         // **********MK
         // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
         ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
         I_MK= static_cast<ULONG>((this->cwclass.cwt_used[this->cwclass.get_Tclassification(i)]== I_KNOT ? this->cwclass.SKNOT_M_info[i]: 0));
#endif
         // **********CW-V
         ULONG I_CWV=0;
#ifdef _USE_CWC_V_
         I_CWV=this->cwclass.get_Vclassification(i);
#endif
         // **********MV
         ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
         I_VK= (this->cwclass.cwv_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[i]: 0);
#endif
        // **********C1
        // Get the corresponding bin in the two invariants of the shear of the tidal field
         ULONG I_C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
         real_prec C1 = this->cwclass.Invariant_TF_II[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
         real_prec C1 = this->cwclass.DELTA2[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
        // **********C2
         ULONG I_C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
         real_prec C2 = this->cwclass.Invariant_TF_III[i];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
         real_prec C2 = this->cwclass.DELTA3[i];
         I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif
       ULONG I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
         real_prec C3 = this->cwclass.Invariant_TF_IV[i];
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined (_USE_TIDAL_ANISOTROPY_)
      real_prec C3 = this->cwclass.Tidal_Anisotropy[i];
      I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
         real_prec C3 = this->cwclass.S2[i];             // s²
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
         ULONG I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
         real_prec CV1 = this->cwclass.Invariant_VS_I[i];
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
         real_prec CV1 = this->cwclass.N2D[i];      // Nabla² ð
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_I_)
         real_prec CV1 = this->cwclass.Invariant_TF_I[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#endif
         ULONG I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
         real_prec CV2 = this->cwclass.Invariant_VS_II[i];
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2,this->s_deltas.prop8, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
         real_prec CV2 = this->cwclass.S2DELTA[i];         // s²ð
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_II_)
         real_prec CV2 = this->cwclass.Invariant_TF_II[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif
         ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
         real_prec CV3 = this->cwclass.Invariant_VS_III[i];
         I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3,this->s_deltas.prop9, this->bin_accumulate_borders);
#elif defined _USE_S3_
         real_prec CV3 = this->cwclass.S3[i];                                   // s³
         I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_III_)
         real_prec CV3 = this->cwclass.Invariant_TF_III[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
         I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif
#ifndef _BIN_ACCUMULATE_
         if(dm >=this->s_mins.prop1 && dm <=this->s_maxs.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
           if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
             if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
               if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                 if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                   if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                     if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
#endif
#ifdef _use_filled_cells_
                       if(true==filled_cells[i])
#endif
                       {
                           // Get the bin in DM
                           ULONG index_el=index_11d(I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,static_cast<ULONG>(this->params._n_cwt()),static_cast<ULONG>(this->params._n_sknot_massbin()),static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()),N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
#if defined (_new_app_) || defined (_assign_again_)
                            cell_info[i].theta_index=index_el;
#endif
#ifdef _new_app_
                            cell_info[i].delta_index=I_X;
#endif
#ifdef _achtung_inside_
                           // check if prob is zero for all possible valoes of Nh
                           // if aux_h=0, given that the quantity CWT_hist is always > 0, means that all values of prob in this bin are zero
                           // Hence no cells with this DM were found to have halos (for all values of ny).
                           // This is likely to happen when we apply the Bias and the Kernel to a different DM field in order to create a mock.
                                             // During the calibration procedure, this must not happen (by construction)
                           double aux_h=0;
                           for(int ih=0;ih<ny;++ih)
                               aux_h +=  static_cast<double>(this->BIAS_NCOUNTS[index_2d(ih,index_el, lenght_bias)]);//this applyes for cases 1 and 2
#ifdef _assign_again_
                            cell_info[i].auxvar=aux_h;
#endif

                           if(aux_h>0) //This applies in the case of generating mocks. This is likely the reason why parallelization in the mock production is not working, for it works for all cells......
                             {
#endif // end _achtung_inside_
                               // Number of cells in the density bin and CWT.
                               // This is always greater than zero, by construction. A fixed value inside the loop
                               ULONG N_available_positions_in_denbin_cwt_Mk = this->NCELLSperDMBIN[index_el];
                               // Number of used cells for the current bin of DM and CWT
                               ULONG N_used_positions_in_denbin_cwt_Mk = NCELLSperDMBIN_now[index_el];
                               // If the density bin to which this cell belongs to is already filled
                               // then go below to get a value of Nhalos selected according to the original distribution in this density bin
                               // setting a flag=false
                               bool flag = true;
                               if (N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
                                 flag=false;
                               // If the density bin still has available positions to be assigned, then proceed
                               if(true==flag)
                                 {
                                   bool cell_accepted=false;
                                   int halo_prop;
                                  // This density bin as available positions. Proceed until the cell is assigned a value of Nhalo
                                   while (false==cell_accepted)
                                     {
                                       // Throw one value of Nhalos in the range [0, nmax_y_onecell], where  nmax_y_onecell represents the
                                       // maximum number of tracers in one cell, read from the cell number counts.

                                       // Accept this number according to the normalized distribution of Nhalos in the corresponding density bin.
                                       // Everytime one cell is assigned this value of Nhalo, the normalized distribution is updated (below)
                                       // and computed extracting that position already assigned. This will ensure that more probability
                                       // is given to the remaining available positions.
                                       real_prec prob_ncounts=-10.0;
                                       real_prec ran = 10.0;
                                       ULONG index=static_cast<ULONG>(LARGE_NUMBER);
                                       // Probability for number counts
                                       while (prob_ncounts<ran)
                                         {
#if defined _USE_MASS_FIELD_
                                           if(1==case_prop || 2==case_prop)
#endif
                                              halo_prop= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins in the Y distirbution (id _DENSITY_)

#ifdef _USE_MASS_FIELD_
                                           else if(4==case_prop)// this is not expected
                                             halo_prop= static_cast<int>(floor(static_cast<real_prec>(ny)*gsl_rng_uniform(r)));                         // draw bins n the mass distribution
#endif
                                           index  = index_2d(halo_prop,index_el,lenght_bias );
                                           prob_ncounts = static_cast<real_prec>(X_Y_hist_dynamical_normalized[index]);
                                           ran   = gsl_rng_uniform(gBaseRand);
                                       }// closes while(prob_ncounts<ran)
                                       ULONG N_available_positions=this->BIAS_NCOUNTS[index];
                                       ULONG N_used_positions = Ncells_density_bin_new[index];
                                       // Now proceed to assign the value of Nhalos to the current cell i, only in case
                                       // that we have available positions in the current Den-Nh bin.
                                       // If there are no available positions, then cell_accepted will be false and we cannot
                                       // get out pf the while loop. A new Nhalo will be proposed.
                                       //For the other properties such as Mass in a cell or the fraction of satellits, we
                                       // we let the available number of cells be the driver criteria
                                       if(N_available_positions > N_used_positions)
                                         {
                                           this->delta_Y_new[i]=static_cast<real_prec>(halo_prop);
                                          Ntracer_dyn+=this->delta_Y_new[i];
                                          counter_dyn++;
#ifdef _assign_again_
                                          NASS++;
#endif
                                           // Claim the current cell as already assigned a value of Nhalos. This breaks the while and continues to the next cell.
                                           cell_accepted=true;
#if defined (_new_app_) || defined (_assign_again_)
                                           cell_info[i].assigned=true;
#endif
#ifdef _new_app_
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           pdf_temp[halo_prop]++;
#endif
                                           // Add one (i.e, the accepted) to count the number of positions used in the current den-N bin
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           Ncells_density_bin_new[index]++;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           NCELLSperDMBIN_now[index_el]++;
                                           // In the current Den-N bin, subtract one (i.e, the accepted) in order
                                           // to upgrade the distribution and give more weight to the remaining available positions
                                           // Attention, do not ask whether this is <0, for it is an unsigned long
                                           if(X_Y_hist_dynamical[index] >=1)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                             X_Y_hist_dynamical[index]--;
                                           else
                                             X_Y_hist_dynamical[index]=0;

                                           // Normalize the dynamical histogram,  in each Den bin, to the maximum.
                                           // The resulting distribution, in this den bin, will be used for the next cell
                                           // in order to accept or reject the new value of Nhalos,
                                           // such that more probability is given to the remaining available positions.

                                           long aux_a=-10000;
                                           long lkk;
                                           for(int j=0; j< ny ;++j)// loop over the number of tracers in a mesh cell
                                           {
                                             long AUX=X_Y_hist_dynamical[index_2d(j,index_el, lenght_bias)];
                                             lkk=max(aux_a, AUX);
                                             aux_a=lkk;
                                           }
                                           for(int j=0;j< ny ;++j)
                                           {
                                             ULONG iin=index_2d(j,index_el, lenght_bias);
                                             ULONG AUX=X_Y_hist_dynamical[iin];
                                             X_Y_hist_dynamical_normalized[iin]= (lkk == 0  ? 0.0 :  static_cast<double>(AUX)/static_cast<double>(lkk));
                                           }
                                         }   // end   if(N_available_positions > N_used_positions) .
                                     }// end  while (false==cell_accepted)
                                 } // end of if(true==flag)
#if !defined (_new_app_)  || !defined (_assign_again_)
                               else
                                 // If the density bin has been already filled, we call the full joint distribution in that bin to assign randomly the value Nhalo according to the other properties
                                 {
                                   real_prec prob=-1.0;
                                   real_prec ran=0.0;
                                   int Nhalos_orfan=0;
                                   while(prob<ran)
                                    {
                                      Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins of a density-like quantity
                                      prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index_el, lenght_bias)];
                                      ran = gsl_rng_uniform(gBaseRand);
                                     }
                                   this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
                                   Ntracer_glob+=this->delta_Y_new[i];
                                   counter_orphan++;
                                }// end else
#endif
#ifdef _achtung_inside_
                             }// end of if(aux_h>0), available only under  _achtung_inside_

                           else // if in the theta_bin the prob is **ALWAYS zero**, assign Poisson distributed number (or zero) with mean that of the reference;
                             {
                              // This cells can be assigned using the marginalized bias, TO BE DONE
                              // if(this->iteration==this->N_iterations_Kernel)
                                 this->delta_Y_new[i]=0;
                            //    this->delta_Y_new[i]=gsl_ran_poisson(gBaseRand,this->mean_number_y_onecell); //This gives a Poisson realization with mean=mean_part_in_cell. More realistic perhaps
#ifndef _new_app_
                                counter_out_of_theta++;
#endif
#ifdef _new_app_
                                cell_info[i].assigned=true;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                pdf_temp[0]++;
#endif //ende new_app
                               }
#endif // end achtung inside
            }// end if(dm is in the defined range)
#ifdef _assign_again_
        } // cloes if cell not assined
#endif
     }// end loop over cells
#ifdef _assign_again_
#if defined(_FULL_VERBOSE_)
    So.message_screen_flush("Number cells assigned = ",static_cast<int>(NASS));
#endif
#ifdef _use_more_references_
    if(NASS== this->params._NGRID())
      break;
  } // closes loop over references
#else
  }while(NASS<0.9* this->params._NGRID());// 0.9 is some threshold put by hand. THis while closes a do that I have omitted,
#endif
#ifdef _FULL_VERBOSE_
  std::cout<<endl;
#endif
#endif
 this->Nobjects=static_cast<ULONG>(get_nobjects(this->delta_Y_new));
 nt_partial=Ntracer_dyn; // keep rackof the number of traers assigned with dynamical
#if !defined(_assign_again_)
 nt_partial_b=Ntracer_glob; // keep track of the number of tracers assigned with global prob dist
#endif
#if defined _new_app_ || defined (_assign_again_)
 // Given that the biased mocks are assigned more objects (power is below) we can avoid that by reducing the nmax in each cell
 // usefd
 // get the number of tracers so far assigned
// get the ratio btween the number of objects assigned and the total number of tracers form the reference

#ifdef _new_app_
  this->delta_Y.resize( this->params._NGRID(),0);
  this->File.read_array_t<PrecType_Y>(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(), this->delta_Y);
  ULONG Nobjects_ref=static_cast<int>(get_nobjects(this->delta_Y));
/*
  ny= static_cast<int>(floor(ny*static_cast<double>(this->Nobjects)/static_cast<double>(Nobjects_ref)));
  So.message_screen("Reducing maximum occuipation number to", ny);
  So.message_screen("First loop done");
*/
 this->delta_Y.clear(); this->delta_Y.shrink_to_fit();
#endif // end _new_app_

#ifdef _assign_again_   // Here we assign randomly if twe do not manage to assign all with the do-whole loop
  So.message_screen("Start second loop: random");
  So.message_screen("NUmber of cells to assign: ",  this->params._NGRID()-NASS);
  counter_orphan_b=0;
#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan_b)
#endif
 for(ULONG i=0;i< this->params._NGRID();++i)      //Start loop over the cells to assign Nh following th eglobal pbias
   {
    if(cell_info[i].auxvar>0)
      if(false==cell_info[i].assigned)
       {
         ULONG index=cell_info[i].theta_index;
         int Nhalos_orfan;
         real_prec prob=-1.0;
         real_prec ran=0.0;
         while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
          {
            Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
            prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index, lenght_bias)];
            ran = gsl_rng_uniform(gBaseRand);
          }
          this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
          counter_orphan_b++;
          cell_info[i].assigned=true;
          Ntracer_glob+=this->delta_Y_new[i];
       }
  }
   So.message_screen("Second loop done");
#endif // end assign_again
     // One option here is to collect the cells that are not yet assigned and randomize them
     // Then assign number of tracers in cells until the number of tracers reaches a value close to the reference * Some Poisson dispersion.

#endif  // closes #if defined _new_app_ || defined (_assign_again_) in order to allow to close omp region.
#ifdef _USE_OMP_
     gsl_rng_free (gBaseRand);
     }// end parallelized region
#endif
     this->So.DONE();
#if defined _new_app_ || defined (_assign_again_)
     So.message_screen("Partial number of objects in new mock= ", nt_partial);
     So.message_screen("Start third loop"); // this is the same second part of the main loop, where the main prob dist is used.
#ifdef _new_app_
#ifndef _DYN_MISSING_
     const gsl_rng_type * rng_p;
     gsl_rng * gBaseRand_p;
     rng_p = gsl_rng_mt19937;//_default;
     gsl_rng_env_setup();
     gsl_rng_default_seed=125*this->params._realization();
     gBaseRand_p = gsl_rng_alloc (rng_p);
     ULONG Ntracer_new=static_cast<int>(get_nobjects(this->delta_Y_new));
//     ULONG Ntracer_ref_poiss = static_cast<ULONG>(gsl_ran_poisson(gBaseRand,Nobjects_ref));
     ULONG Ntracer_ref_poiss = static_cast<ULONG>(Nobjects_ref*(1+gsl_ran_gaussian(gBaseRand,0.002)));

     So.message_screen("Setting threshold at ",Ntracer_ref_poiss );
     So.message_screen("Number of reference tracers",Nobjects_ref);
     //do{
//       gsl_ran_shuffle(gBaseRand,&i_nac[0], i_nac.size(),sizeof(ULONG));
       //Start loop over the cells that are not yeat assigned
    vector<ULONG>i_nac;
   for(ULONG i=0;i< this->params._NGRID();++i)
      if(false==cell_info[i].assigned)
        i_nac.push_back(i);
    counter_dyn= this->params._NGRID()-i_nac.size();
    So.message_screen("Getting bias as a function of delta");
    ULONG AUX_p=0;
#pragma omp parallel for reduction(+:AUX_p)
    for(ULONG i=0; i<X_Y_hist_dynamical.size();++i)
        AUX_p+=X_Y_hist_dynamical[i];
    So.message_screen("Number of available bins =",AUX_p);
    So.message_screen("Number of cells without Nh=",i_nac.size());
    vector<real_prec> new_bias(this->params._NX()*ny,0);
    So.message_screen("Marginalizing.");
    for(int iy = 0; iy < ny; ++iy)
    {
      for(int i=0;i < this->params._NX();++i)// loop sobre bines de dark matter
       {
         ULONG index=index_2d(iy,i,this->params._NX());
         for(int sua = 0; sua < this->params._n_cwt(); ++sua)
           for(int k = 0; k < this->params._n_sknot_massbin() ;  ++k)
            for(int vua = 0; vua < this->params._n_cwv(); ++vua)
              for(int vk = 0; vk < this->params._n_vknot_massbin(); ++vk)
               for(int l1 = 0; l1 < N_C_BIN1; ++l1)
                 for(int l2 = 0; l2 < N_C_BIN2 ; ++l2)
                   for(int l3 = 0; l3 < N_C_BIN3; ++l3)
                     for(int lv1 = 0; lv1 < N_CV_BIN1; ++lv1)
                       for(int lv2 = 0; lv2 < N_CV_BIN2; ++lv2)
                         for(int lv3 = 0; lv3 < N_CV_BIN3; ++lv3)
                          {
                            ULONG index_el=index_12d(iy,i,sua,k,vua,vk,l1,l2,l3,lv1,lv2,lv3,this->params._NX(),this->params._n_cwt(),this->params._n_sknot_massbin(),this->params._n_cwv(), this->params._n_vknot_massbin(),N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                            new_bias[index]+=X_Y_hist_dynamical[index_el];
                          }
        }
    }
    So.DONE();
  // ahora normalizamos en bines de dm
    vector<real_prec> new_bias_normalized(this->params._NX()*ny,0);
    vector<real_prec> new_bias_normalized_dyn(this->params._NX()*ny,0);
    vector<real_prec> new_bias_dyn(this->params._NX()*ny,0);
    new_bias_dyn=new_bias;
    new_bias_normalized_dyn=new_bias_normalized;
    So.message_screen("Normalizing");
    for(int i=0;i < this->params._NX();++i)// loop sobre bines de dark matter
     {
       long auxp=-100000;
       long lkk;
       for(int iy = 0; iy < ny; ++iy)
        {
          ULONG index=index_2d(iy,i,this->params._NX());
          long aux=new_bias[index];
          lkk=max(aux,auxp);
          auxp=lkk;
        }
        for(int iy = 0; iy < ny; ++iy)
         {
          ULONG index=index_2d(iy,i,this->params._NX());
          ULONG newb=new_bias[index];
          new_bias_normalized[index]= (lkk==0? 0: static_cast<real_prec>(newb)/static_cast<real_prec>(lkk));
         }
     }
    So.DONE();
#endif
    // Loop over those cells not assigned. Runs untlil all cells have been assigned.
    counter_orphan=0;
#ifdef _DYN_MISSING_
     vector<real_prec>new_NCELLSperDMBIN_now(this->params._NX(),0);
     vector<real_prec>new_NCELLSperDMBIN(this->params._NX(),0);
     vector<real_prec>new_Ncells_density_bin(this->params._NX(),0);
     new_NCELLSperDMBIN=new_bias;
     So.message_screen("Dynamicall assigning number counts to missing cells. NOT WORKING");
     for(ULONG ic=0;ic<i_nac.size(); ++ic)
       {
         ULONG i=i_nac[ic];
         ULONG index_el= cell_info[i].delta_index;
         double aux_h=0;
         for(int ih=0;ih<ny;++ih)
           aux_h +=  static_cast<double>(new_bias[index_2d(ih,index_el, ny*this->params._NX())]);
         if(aux_h>0)
          {
            ULONG N_available_positions_in_denbin_cwt_Mk = new_NCELLSperDMBIN[index_el];
            ULONG N_used_positions_in_denbin_cwt_Mk = new_NCELLSperDMBIN_now[index_el];
            bool flag = true;
            if (N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
              flag=false;
            if(true==flag)
             {
              bool cell_accepted=false;
              int halo_prop;
              while (false==cell_accepted)
               {
                real_prec prob_ncounts=-10.0;
                real_prec ran = 10.0;
                ULONG index=static_cast<ULONG>(LARGE_NUMBER);
                while (prob_ncounts<ran)
                 {
                  halo_prop= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins in the Y distirbution (id _DENSITY_)
                  index  = index_2d(halo_prop,index_el,ny*NX );
                  prob_ncounts = static_cast<real_prec>(new_bias_normalized_dyn[index]);
                  ran   = gsl_rng_uniform(gBaseRand);
                }
               ULONG N_available_positions=new_bias[index];
               ULONG N_used_positions = new_Ncells_density_bin[index];
               if(N_available_positions > N_used_positions)
                {
                  this->delta_Y_new[i]=static_cast<real_prec>(halo_prop);
                  cell_accepted=true;
                  Ntracer_glob +=this->delta_Y_new[i];
                  new_Ncells_density_bin[index]++;
                  new_NCELLSperDMBIN_now[index_el]++;
                  counter_orphan++;
                  if(new_bias_dyn[index] >=1)
                    new_bias_dyn[index]--;
                  else
                    new_bias_dyn[index]=0;
                 long aux_a=-10000;
                 long lkk;
                 for(int j=0; j< ny ;++j)// loop over the number of tracers in a mesh cell
                  {
                   long AUX=new_bias_dyn[index_2d(j,index_el, ny*NX)];
                   lkk=max(aux_a, AUX);
                   aux_a=lkk;
                  }
                 for(int j=0;j< ny ;++j)
                  {
                    ULONG iin=index_2d(j,index_el, ny*NX);
                    ULONG AUX=new_bias_dyn[iin];
                    new_bias_normalized_dyn[iin]= (lkk == 0  ? 0.0 :  static_cast<double>(AUX)/static_cast<double>(lkk));
                  }
                }   // end   if(N_available_positions > N_used_positions) .
            }// end  while (false==cell_accepted)
          } // end of if(true==flag)
          else
           {
            real_prec prob=-1.0;
            real_prec ran=0.0;
            int Nhalos_orfan=0;
            while(prob<ran)
             {
              Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //d   cout<<endl;
                  //raw number counts in cells or bins of a density-like quantity
              prob = new_bias_normalized[index_2d(Nhalos_orfan,index_el, NX*ny)];
              ran = gsl_rng_uniform(gBaseRand);
             }
             this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
             Ntracer_glob+=this->delta_Y_new[i];
             counter_orphan++;
           }// end else
         }// end of if(aux_h>0), available only under  _achtung_inside_
         else // if in the theta_bin the prob is **ALWAYS zero**, assign Poisson distributed number (or zero) with mean that of the reference;
          {
           this->delta_Y_new[i]=0;  // This is radical. Might not allow for cosmic variance
           counter_out_of_theta++;
          }
     }
     nt_partial_b=Ntracer_glob;
#else
    for(ULONG i=0;i<i_nac.size();++i)
      {
        ULONG i_cell = i_nac[i]; // original ID of the not assigned cell
    //    ULONG index= cell_info[i_cell].theta_index;
        ULONG index= cell_info[i_cell].delta_index;// thiis is to use the bias as a function of delta
        int Nhalos_orfan=0;
        real_prec prob=-1.0;
        real_prec ran=0.0;
        while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
         {
           Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
//            ULONG Index_t=index_2d(Nhalos_orfan, index, lenght_bias);
           ULONG Index_t=index_2d(Nhalos_orfan, index, NX); // thiis is to use the bias as a function of delta
//            prob = this->BIAS_NCOUNTS_normalized[Index_t];
           prob = new_bias_normalized[Index_t];// thiis is to use the bias as a function of delta
           ran = gsl_rng_uniform(gBaseRand);
         }
         counter_orphan++;
         Ntracer_glob+=Nhalos_orfan;
         this->delta_Y_new[i_cell]= static_cast<real_prec>(Nhalos_orfan);
         cell_info[i_cell].assigned=true;
      }
      new_bias.clear(); new_bias.shrink_to_fit();
      new_bias_normalized.clear(); new_bias_normalized.shrink_to_fit();
      nt_partial_b=Ntracer_glob;
#endif // end dyn_assignment_
/*
// Assign Nh randomly to cells wihtout any distribution. This  brings problems with clustering, although gioves good number of tracers
// This  brings problems with clustering, although gioves good number of tracers
    i_nac.clear();i_nac.shrink_to_fit();
    for(ULONG i=0;i< this->params._NGRID();++i)
         i_nac.push_back(i);
     counter_out_of_theta= this->params._NGRID()-i_nac.size();
    if(i_nac.size()>0) // if we have cells to distribute. When dealing with the DM use for the calibration, this third part is not used for all cells are assigned.
      {
      do{
          gsl_ran_shuffle(gBaseRand,&i_nac[0], i_nac.size(),sizeof(ULONG));
         //Start loop over the cells that are not yeat assigned
          for(ULONG i=0;i<i_nac.size();++i) // loop ove rrandomized not assigned cells
            {
            ULONG i_cell = i_nac[i]; // original ID of the not assigned cell
            ULONG index= cell_info[i_cell].theta_index;
            int Nhalos_orfan=0;
            Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
            counter_orphan++;
            Ntracer_new+=Nhalos_orfan;
            this->delta_Y_new[i_cell]+= static_cast<real_prec>(Nhalos_orfan);
            if(Ntracer_new>Ntracer_ref_poiss)
                break;
         }
       Ntracer_new=static_cast<int>(get_nobjects(this->delta_Y_new));
      }while(Ntracer_new<Ntracer_ref_poiss);
    }
*/
    i_nac.clear();i_nac.shrink_to_fit();
/*
#ifdef _USE_OMP_
#pragma omp for reduction(+:counter_orphan)
#endif
     //Start loop over the cells
     for(ULONG i=0;i< this->params._NGRID();++i)
       {
        if(false==cell_info[i].assigned)
          {
            ULONG index=cell_info[i].theta_index;
            int Nhalos_orfan;
            real_prec prob=-1.0;
            real_prec ran=0.0;
            while(prob<ran) // This while aims at chosing the right Nh according to the full propability distribution
              {
                 Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells
                 prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index, lenght_bias)];
                 ran = gsl_rng_uniform(gBaseRand);
              }
              this->delta_Y_new[i]= static_cast<real_prec>(Nhalos_orfan);
              counter_orphan++;
           }
      }
*/
     So.message_screen("Third loop done");
    cell_info.clear();cell_info.shrink_to_fit();
#endif // end_new_app_
#endif // end_new_app_ ior new_assingn
    X_Y_hist_dynamical.clear();X_Y_hist_dynamical.shrink_to_fit();
    X_Y_hist_dynamical_normalized.clear();X_Y_hist_dynamical_normalized.shrink_to_fit();
#ifdef _new_counter_
    new_X_Y_hist_dynamical.clear();new_X_Y_hist_dynamical.shrink_to_fit();
    new_X_Y_hist_dynamical_normalized.clear();new_X_Y_hist_dynamical_normalized.shrink_to_fit();
#endif
    this->NCELLSperDMBIN.clear(); this->NCELLSperDMBIN.shrink_to_fit();
    NCELLSperDMBIN_now.clear();    NCELLSperDMBIN_now.shrink_to_fit();
    Ncells_density_bin_new.clear(); Ncells_density_bin_new.shrink_to_fit();
//     testa.close();
#ifdef _FULL_VERBOSE_
   So.message_time_mock(start_mock);
#endif
   // ******************************************** end assigning ncounts or other prop*******************************
#ifdef _FULL_VERBOSE_
   So.message_screen("Fraction of cells assigned with DYNAMICAL prob. dist =", 100.0*static_cast<double>(counter_dyn)/static_cast<double>( this->params._NGRID()), "%");
#ifdef _assign_again_
 if(counter_orphan_b>0)
   So.message_screen("Fraction of cells assigned with GLOBAL prob. dist =", 100.0*static_cast<double>(counter_orphan_b)/static_cast<double>( this->params._NGRID()), "%");
#else
 if(counter_orphan>0)
   So.message_screen("Fraction of cells assigned with GLOBAL prob. dist =", 100.0*static_cast<double>(counter_orphan)/static_cast<double>( this->params._NGRID()), "%");
#endif
 if(counter_out_of_theta>0)
     So.message_screen("Fraction of cells without representation in the ref bias=", 100.0*static_cast<real_prec>(counter_out_of_theta)/static_cast<real_prec>( this->params._NGRID()), "%");
#endif
   string file_info=this->params._Output_directory()+"info_cells_realization"+to_string(this->params._realization())+"_reference"+to_string(this->params._seed())+".txt";
   ofstream fi; fi.open(file_info.c_str());
   fi<<"# Fraction of cells (%): dynamical, global(orphan)"<<endl;
#ifdef _assign_again_
   fi<<100.0*static_cast<double>(counter_dyn)/static_cast<double>( this->params._NGRID())<<"  "<<100.0*static_cast<double>(counter_orphan_b)/static_cast<double>( this->params._NGRID())<<endl;
#else
   fi<<100.0*static_cast<double>(counter_dyn)/static_cast<double>( this->params._NGRID())<<"  "<<100.0*static_cast<double>(counter_orphan)/static_cast<double>( this->params._NGRID())<<endl;
#endif
   fi<<"# Number of tracers in each case"<<endl;
   fi<<nt_partial<<"   "<<nt_partial_b<<endl;
   fi<<"# Total"<<endl;
   fi<<nt_partial+nt_partial_b<<endl;
   fi.close();
   if(property==_COUNTS_)
     {
       this->Nobjects=static_cast<int>(get_nobjects(this->delta_Y_new));
#ifdef _FULL_VERBOSE_
       So.message_screen("Number of objects in new mock= ", this->Nobjects);
       So.message_screen("Mean number density = ", static_cast<real_prec>(this->Nobjects)/pow(this->params._Lbox(), 3), "(Mpc / h )⁻³");
#endif
#ifdef _DO_BiasMT_CALIBRATION_
#ifndef _EXTRAPOLATE_VOLUME_
#ifdef _FULL_VERBOSE_
       if(this->Nobjects < this->N_objects_Y)
         So.message_screen("->Less objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
       else if(this->Nobjects  > this->N_objects_Y)
         So.message_screen("->More objects assigned to mock. Fraction = ", fabs(100.0-100.0*(static_cast<real_prec>(this->Nobjects)/this->N_objects_Y)), "%");
#endif
#endif
#endif
       if(false==silent)
         {
           real_prec anmean=static_cast<real_prec>(this->Nobjects)/static_cast<real_prec>( this->params._NGRID());
#ifdef _FULL_VERBOSE_
           So.message_screen("<N> objects in mock =", anmean);
#endif
           int nyy=30;
           real_prec dmax=300.0;
           vector<real_prec>DELTA_HIST(nyy,0);
           real_prec delta_delta= (dmax+1.0)/static_cast<real_prec>(nyy);  //deltamax - deltamin
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0;i<  this->params._NGRID();++i)
             {
               real_prec delta_mock= (this->delta_Y_new[i]/static_cast<real_prec>(anmean) - 1.0);
               if(delta_mock>=-1.0 && delta_mock<=dmax)
                 {
                   int ID=static_cast<int>(floor(delta_mock +1.0)/delta_delta);
                   if(ID==nyy)ID--;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   DELTA_HIST[ID]++;
                 }
             }
           real_prec lkk=get_max<real_prec>(DELTA_HIST);
           for(int j=0;j< nyy ;++j)
             DELTA_HIST[j]= lkk == 0  ? 0.0 :  static_cast<real_prec>(DELTA_HIST[j])/static_cast<real_prec>(lkk);

           ofstream aja;
           aja.open("delta_mock_dist.txt");
           for (int i=0;i<DELTA_HIST.size();++i)aja<<i<<"\t"<<DELTA_HIST[i]<<endl;
           aja.close();
           DELTA_HIST.clear();
           DELTA_HIST.shrink_to_fit();
           vector<real_prec>AUX( this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0; i<  this->params._NGRID() ;++i)
             AUX[i]=this->delta_Y_new[i]/static_cast<real_prec>(anmean) - 1.0;
           lkk=get_max<real_prec>(AUX);
           So.message_screen("Maximum delta in mock =", lkk);
           lkk=get_min<real_prec>(AUX);
           So.message_screen("Minimum delta in mock =", lkk);
           AUX.clear();
           AUX.shrink_to_fit();
        }//closes if(false==silent)
       if(true==this->params._Write_PDF_number_counts()) 	  // Write the PDF of the outcome */
         {
               string fileY=this->params._Output_directory()+"PDF_NC"+"_Y_MOCK"+this->new_Name_Property_Y+"_MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+"_"+this->params._file_power()+".txt";
                So.message_screen("Writting PDF_NC for MOCKS");
               int NPART=ny; // Number of particles in cells
               this->PDF_NC_Y.clear();
               this->PDF_NC_Y.resize(NPART, 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
               for(ULONG i=0;i<  this->params._NGRID();++i)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 this->PDF_NC_Y[static_cast<int>(this->delta_Y_new[i])]++;  // Get pdf of tracer number counts
               this->File.write_to_file_i(fileY,this->PDF_NC_Y);// Write pdf to file
               vector<int> cells_with_one;
               for(int i=0;i<this->PDF_NC_Y.size();++i)
                 if(this->PDF_NC_Y[i]>0)
                   cells_with_one.push_back(i);
               So.message_screen("Maximum number of mock-objects in one cell =", get_max<int>(cells_with_one) );

               cells_with_one.clear();
               this->PDF_NC_X.clear();
               cells_with_one.shrink_to_fit();
               this->PDF_NC_X.shrink_to_fit();
                }
       } //closes if(true==this->Write_PDF_number_counts)
   else if (_DENSITY_==property)
     {
       // if we have Y as density, then we have to sample the bins in Y to rweturn values of a Y field
       gsl_rng_env_setup();
       gsl_rng_default_seed=236655;//time(NULL);
       const gsl_rng_type *  Tn= gsl_rng_ranlux;
       gsl_rng * rn = gsl_rng_alloc (Tn);
#ifdef _DISPLACEMENTS_
    So.message_screen("Building Displacement");
#endif
      if(true==this->params._Convert_Density_to_Delta_Y())// Here I am also assuming that if this is true, we also took the log10(num_in_log+delta)
       {
        for(ULONG i=0;i<  this->params._NGRID();++i) //tbc
              {
                 real_prec xr=gsl_rng_uniform (rn);
           real_prec aux_den=this->s_mins.prop0+(static_cast<int>(delta_Y_new[i])+xr)*this->s_deltas.prop0;   // Get value of log10(2+delta_y)
                 this->delta_Y_new[i]=this->Mean_density_Y*(1.0+pow(10,aux_den)-NUM_IN_LOG);                                        // Get nbar*(1+delta_y)
#ifdef _DISPLACEMENTS_
          this->delta_Y_new[i]+=this->Displacement_inicial[i];  // Construct psi new = delta Psi + psi_ref
#endif
          }
      }
      else{
        for(ULONG i=0;i<  this->params._NGRID();++i) //tbc
          {
           real_prec xr=gsl_rng_uniform (rn);
            int bin_prop=static_cast<int>(delta_Y_new[i]);
           this->delta_Y_new[i]=this->s_mins.prop0+(bin_prop+xr)*this->s_deltas.prop0;   // Get value of log10(2+delta_y)
#ifdef _DISPLACEMENTS_
           this->delta_Y_new[i]+=this->Displacement_inicial[i];  // Construct psi new = delta Psi + psi_ref
#endif
          }//closes for
      }// closes else
    So.DONE();
  }// closes  else if (_DENSITY_==property)
#ifdef  _EXTRAPOLATE_VOLUME_
#ifdef _WRITE_TR_DENSITY_FIELD_
   this->File.write_array(fname, this->delta_Y_new);
#endif
#endif
   // this is commented for only was applicable to DENSITY case. For the case in which we calibrate ncounts
   // and other properties, those extra fields are written directly as they come in the last iteration
   if(case_prop==4)
     {
       gsl_rng_env_setup();
       gsl_rng_default_seed=236655;//time(NULL);
       const gsl_rng_type *  Tn= gsl_rng_ranlux;
       gsl_rng * rn = gsl_rng_alloc (Tn);
#ifndef _USE_LOG_MASS_
       real_prec mass_min=pow(10,params._LOGMASSmin())*params._MASS_units();
        So.message_warning_ini(__LINE__,__PRETTY_FUNCTION__,__FILE__,"Calibrating from sat fraction. Please check definition of containers vector<int>ncounts and vector<bool>filled_cells");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<  this->params._NGRID();++i)
         {
           real_prec xr=gsl_rng_uniform (rn);
           real_prec aux=this->s_mins.prop0_mass+(static_cast<int>(delta_Y_new[i])+xr)*this->s_deltas.prop0_mass;   // Get value of log10(2+delta_y)
           real_prec mass_cell = (pow(10,aux)-NUM_IN_LOG)*MASS_SCALE;
           real_prec filled=1;
           int ncounts_m=1;
#ifdef _use_filled_cells_
           filled=static_cast<real_prec>(filled_cells[i]);
           ncounts_m=ncounts[i];
#endif
           this->delta_Y_new[i]=mass_cell*filled;
           int Neff = static_cast<int>(floor(this->delta_Y_new[i]/mass_min));
           // This weights help to upweight the mass at a cell such that, give the number counts it has, it can at least provide mass for all particles
           real_prec weight = (Neff >= ncounts_m ? 1.0 : static_cast<real_prec>(ncounts_m)/static_cast<real_prec>(Neff));
           this->delta_Y_new[i]*=weight;
         }//closes for
       gsl_rng_free (rn);
       So.message_screen("Minimum of new mass field = ", get_min(this->delta_Y_new));
       So.message_screen("Maximim of new mass field = ", get_max(this->delta_Y_new));
#endif
     } //closes if(case_prop==4)
#ifdef _use_filled_cells_
   ncounts.clear();
   ncounts.shrink_to_fit();
#endif
   // Write the density fields of properties of mock tracers:.
#ifndef _GET_BiasMT_REALIZATIONS_
#ifndef _TEST_THRESHOLDS_RESIDUALS_
   auto out_it = std::find(std::begin(this->params._output_at_iteration()), std::end(this->params._output_at_iteration()), this->iteration);
//   if((this->iteration==this->N_iterations_Kernel) || (out_it == std::end(this->params._output_at_iteration())))
    if(this->iteration==this->params._N_iterations_Kernel())
      this->File.write_array(fname, this->delta_Y_new);
#endif
#else
     this->File.write_array(fname, this->delta_Y_new);
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _USE_LPT_
   if(case_prop==1 || case_prop==2)
     this->lpt.set_fname_MOCK_NCOUNTS(fname);
   if(case_prop==4)
     this->lpt.set_fname_MOCK_MASS(fname);
//   if(case_prop==3)
//     this->lpt.set_fname_MOCK_NCOUNTS_SAT(fname);
#else
   So.message_warning("CHeck names in line",__LINE__);
#endif
#endif
   // Get power spectrum only number counts
   if(case_prop==1 || case_prop==2)
   {
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _USE_GNUPLOT_
    this->get_power_spectrum("TR_MOCK_REALIZATION");
    this->delta_Y.resize( this->params._NGRID(),0);
    this->File.read_array_t<PrecType_Y>(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(), this->delta_Y);
#ifdef _EXTRAPOLATE_KERNEL_
    this->params.set_Nft(this->params._Nft_low());
    this->params.set_Lbox_low(this->params._Lbox_low());
#endif
    this->get_power_spectrum("TR_REF");
    this->get_power_spectrum("TR_MOCK");
    this->delta_Y.clear();this->delta_Y.shrink_to_fit();
    this->gp_power<<"set border linewidth 2.\n";
   this->gp_power<<"set size square \n";
   this->gp_power<<"set log x \n";
   this->gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
   this->gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";

     vector<pair<real_prec, real_prec> > xy_pts_ref;
     vector<pair<real_prec, real_prec> > xy_pts_new;
     for(int i=0; i<kvec.size(); ++i)
       xy_pts_ref.push_back(std::make_pair(this->kvec[i], log10(this->Power_REF[i])));
     for(int i=0; i<kvec.size(); ++i)
       xy_pts_new.push_back(std::make_pair(this->kvec[i], log10(this->Power_NEW[i])));
     this->gp_power<<"plot "<<this->gp_power.file1d(xy_pts_ref) << "w l lw 3.3 lt 2 title 'Reference',"<<this->gp_power.file1d(xy_pts_new)<< " w l lw 3.3 lt 6 title 'Mock'"<<endl;
     this->gp_power<<"set terminal pngcairo\n";
     xy_pts_ref.clear();
     xy_pts_ref.shrink_to_fit();
     xy_pts_new.clear();
     xy_pts_new.shrink_to_fit();
#endif
#endif
#ifdef _EXTRAPOLATE_KERNEL_
    this->params.set_Nft(this->params._Nft());
    this->params.set_Lbox(this->params._Lbox());
#endif
   }// closes if
   this->shot_noise_new = pow(this->params._Lbox(),3)/static_cast<real_prec>(this->Nobjects);
   this->shot_noise_ref=this->shot_noise_new;  // This should be the case during the calibration
   this->delta_Y_new.clear();
   this->delta_Y_new.shrink_to_fit();
 }
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_TWO_REFS_MOCKS_
#define _achtung_inside_  // This check inside the mock loop whether the cond probability is >0, such that we do not fall in infinite while loops
//#define testz  // if this is defined, this does the same as get_grid_mock, i.-e, uses on reference
void BiasMT::get_mock_grid_two(string property)
 {
   So.enter(__PRETTY_FUNCTION__);
   bool silent=true;
#ifdef _FULL_VERBOSE_
   So.message_screen("**Generating new Mock (tracer) number counts based on two references");
   std::cout<<endl;
#ifdef _USE_OMP_
 int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
   So.message_screen("Function get_mock_grid is using",omp_get_max_threads(),"threads");
#endif
   std::cout<<endl;
#endif
    ULONG nx=1;
#ifdef _USE_DM_IN_BiasMT_
   // Define seed vector for each thread in the paralellized run
   nx=this->params._NX();  // Number of bins id DM
#endif
  ULONG lenght_bias=N_CV_BIN1 * N_CV_BIN2* N_CV_BIN3 * N_C_BIN1 * N_C_BIN2* N_C_BIN3 * this->params._n_sknot_massbin() * this->params._n_cwv() * this->params._n_vknot_massbin() * this->params._n_cwt() * nx;
  int n_refs=this->params._Number_of_references();
  int n_new_dm=this->params._Number_of_new_mocks();
  struct tracers_dfield{
   vector<real_prec>delta_field;// container for the number counts of the mocks to be simultaneously generated
   vector<ULONG>id_cells;//container for the id of the cells in each box
   string fname;
   bool file_exist;
   ULONG Ntracer_dyn;
   ULONG Nobjects;
   ULONG Ntracer_glob;
   ULONG counter_dyn;
   ULONG counter_orphan;
   ULONG mesh_id_counter;
   Cwclass cwclass_field;
   void get_nobjects_tracers(){this->Nobjects=get_nobjects(this->delta_field);}
  };
   string pname,fname;
   int case_prop=0;
   if(_COUNTS_==property)
            case_prop=1;
   if(_DENSITY_==property)
     case_prop=2;
   ULONG ny = this->new_nbins_y;
   ULONG size_bias_array= ny*lenght_bias;//    this->BIAS_NCOUNTS.size();
   this->BIAS_NCOUNTS.clear();
   this->BIAS_NCOUNTS.shrink_to_fit();
   this->BIAS_NCOUNTS.resize(size_bias_array,0);
   // ADDING THE BIAS FROM THE TWO references
  // THIS HAS TO BE OPTIMIEZ FROM THE PARAMETER FILE
#ifdef _FULL_VERBOSE_
  So.message_screen("Reading and merging reference bias");
#endif
 vector<string>bias_files;
  if(n_refs>1)
   {
     vector<real_prec> aux_vec;
     for(int il=0;il<n_refs;++il)
      {
       aux_vec.resize(size_bias_array,0);
  //     cout<<size_bias_array<<endl;
      this->File.read_array(this->params._files_bias_references(il), aux_vec);
/*
       ULONG ncheck=0;
#pragma omp parallel for reduction(+:ncheck)
       for(ULONG i=0;i<aux_vec.size();++i)
         ncheck+=aux_vec[i];
     if(ncheck!= this->params._NGRID() )
      {
         So.message_warning("Number of cells in BIAS not equal to total numbr of cells.");
         cout<<ncheck<<"   "<<this->params._NGRID()<<"   "<< aux_vec.size()<<endl;
       }
 */
     for(ULONG i=0;i<aux_vec.size();++i)
      this->BIAS_NCOUNTS[i]+=aux_vec[i];
    }
    aux_vec.clear();aux_vec.shrink_to_fit();
   }
  else
    this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias.dat", this->BIAS_NCOUNTS);
 bias_files.clear();bias_files.shrink_to_fit();
 So.DONE();
//Normalize the bias
#ifdef _FULL_VERBOSE_
   So.message_screen("Normalizing total bias");
#endif
   this->BIAS_NCOUNTS_normalized.clear();
   this->BIAS_NCOUNTS_normalized.shrink_to_fit();
   this->BIAS_NCOUNTS_normalized.resize(size_bias_array,0);
  long aux_a_full=-10000;
  long aux_b_full;
  ULONG aux_b_empty=0;
  if(n_refs>1)
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i < lenght_bias;++i)// loop sobre bines de dark matter
     {
        long aux_a=-10000;
        long aux_b;
        for(int tr_j = 0; tr_j < ny; ++tr_j)
          {
            ULONG inde=index_2d(tr_j, i, lenght_bias);
            long AUX=this->BIAS_NCOUNTS[inde];
            if(AUX==0)
              aux_b_empty++;
            aux_b=max(aux_a, AUX);
            aux_a=aux_b;
            aux_b_full=max(aux_a_full, AUX);
            aux_a_full=aux_b_full;
           }
         for(int tr_j = 0; tr_j < ny; ++tr_j)
          {
            ULONG inde=index_2d(tr_j, i, lenght_bias);
            long AUX=this->BIAS_NCOUNTS[inde];
            this->BIAS_NCOUNTS_normalized[inde] = (aux_a<=0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
          }
     }
   }
    else
       this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_Normalized.dat", this->BIAS_NCOUNTS_normalized);
  So.DONE();
  if(n_refs==1)
    {
      So.DONE();
      real_prec KK = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS));
      real_prec Kd = static_cast<real_prec>(get_nobjects(this->BIAS_NCOUNTS_normalized));
      if(Kd<=0){
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Conditional probablity ill-defined. BAM stops here.", Kd);
       exit(0);
      }
  }
  vector<ULONG> X_Y_hist_dynamical(size_bias_array);
  vector<real_prec>  X_Y_hist_dynamical_normalized(size_bias_array);
  X_Y_hist_dynamical=this->BIAS_NCOUNTS;
  X_Y_hist_dynamical_normalized=this->BIAS_NCOUNTS_normalized;
  this->NCELLSperDMBIN.clear();
  this->NCELLSperDMBIN.shrink_to_fit();
  this->NCELLSperDMBIN.resize(lenght_bias, 0);
#ifdef _FULL_VERBOSE_
   So.message_screen("Marginalizing");
#endif
   if(n_refs>1)
   {
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
   for(ULONG i=0; i< lenght_bias ;++i)
       for(ULONG j=0; j< ny ; ++j)
         {
           ULONG index_h=index_2d(j,i,lenght_bias);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
           this->NCELLSperDMBIN[i]+= this->BIAS_NCOUNTS[index_h];
          }
  }
  else
   this->File.read_array(this->params._Input_Directory_BIAS_KERNEL()+"Bam_Bias_AUX.dat",this->NCELLSperDMBIN);
  So.DONE();
   // Vector to allocate the number of cells in a given density bin during the mapping */
   vector<ULONG>Ncells_density_bin_new(size_bias_array , 0);
   // Vector containing the number of cells in a given density bin, updated everytime a cell has been assigend a value of Nhalos
   vector<ULONG> NCELLSperDMBIN_now(lenght_bias, 0);
// In this part we enable the possiboility to extapolate the kernel to the current L from  LOWER VALUE OF l.
#ifdef _FULL_VERBOSE_
   So.message_screen("Building Kernel");
#endif
   cout<<n_refs<<endl;
  this->Kernel.clear();this->Kernel.resize(this->NTT,0);
  for(int j=0; j<n_refs; ++j)  // Loop over nrefs-1 refereces. The first reference will be counted here, hence the function get_new_dm field is not used when this function is called
  {
#ifndef _EXTRAPOLATE_KERNEL_
    vector<real_prec>kernel_ghost(this->NTT,0);
    this->File.read_array(this->params._files_kernel_references(j), kernel_ghost); // FIX NAME OF GHOST DM, I have used one of the calibration
#else
    vector<real_prec>low_kernel_ghost(this->NTT_low,0);
    vector<real_prec>kernel_ghost(this->NTT,0);
    this->File.read_array(this->params._files_kernel_references(j), low_kernel_ghost); // FIX NAME OF GHOST DM, I have used one of the calibration
    this->extrapolate_kernel(low_kernel_ghost,kernel_ghost);
    low_kernel_ghost.clear();low_kernel_ghost.shrink_to_fit();
#endif

    for(ULONG i=0; i<this->NTT; ++i)
       this->Kernel[i]+=kernel_ghost[i]/static_cast<real_prec>(n_refs);
  }
  So.DONE();
#ifdef _FULL_VERBOSE_
   So.message_screen("Reading and computing theta properties for new DM fields");
#endif

   vector<int>list_aux(n_new_dm,0);
   vector<string>list_file_aux(n_new_dm);
  for(int j=0; j<n_new_dm; ++j)  // Loop over new dm fields
     list_aux[j]=this->params._list_new_dm_fields(j);
  for(int j=0; j<n_new_dm; ++j)  // Loop over new dm fields
     list_file_aux[j]=this->params._files_new_dm_fields(j);
   string file_aux_power=this->params._file_power();
  for(int ih=0;ih<this->params._Number_of_chunks_new_dm(); ++ih)
   {
      vector<tracers_dfield>dm_new(n_new_dm); // container to allocate the dm density fields and id of references
      vector<tracers_dfield>tracers(n_new_dm);  // container to allocate the fields and id of each new tracer genreated
#ifdef _FULL_VERBOSE_
   So.message_screen("Resizing structure for deltas");
#endif
  for(ULONG i=0;i<n_new_dm;++i)
    {
      tracers[i].delta_field.resize( this->params._NGRID(),0);
      dm_new[i].delta_field.resize( this->params._NGRID(),0);
      dm_new[i].id_cells.resize( this->params._NGRID(),0);
      tracers[i].Ntracer_dyn=0;
      tracers[i].Ntracer_glob=0;
      tracers[i].counter_dyn=0;
      tracers[i].counter_orphan=0;
      dm_new[i].mesh_id_counter=0;
      tracers[i].file_exist=false;
      for(ULONG j=0;j< this->params._NGRID();++j)
        dm_new[i].id_cells[j]=j;
    }
  So.DONE();
 if(ih>0)
     for(int j=0; j<n_new_dm; ++j)  // Loop over new dm fields
       list_aux[j]=list_aux[j]+n_new_dm;
    if(ih>0)
     for(int j=0; j<n_new_dm; ++j)  // Loop over new dm fields
#ifdef _SLICS_
       list_file_aux[j]=this->params._Input_Directory_X_NEW()+"densDMALPTrS20.0TETCICz1.041G192V505.0S"+to_string(list_aux[j])+".dat";
#elif defined _UNITSIM_ || defined _ABACUS_
       list_file_aux[j]=this->params._Input_Directory_X_NEW()+this->params._Name_Catalog_X_NEW();
#endif
    for(int j=0; j<n_new_dm; ++j)  // Loop over new dm fields
     {
      tracers[j].fname=this->params._Output_directory()+"MOCK_TR_Realization"+to_string(list_aux[j])+"_"+"MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft());//+"_z"+to_string(this->params._redshift());
      vector<real_prec>delta_ghost( this->params._NGRID(),0);
      string file_dm=list_file_aux[j];
      this->File.inStreamp.open(file_dm.data(), file_is_natural);
      tracers[j].file_exist = this->File.inStreamp.is_open();
      this->File.inStreamp.close();
      dm_new[j].cwclass_field.set_params(this->params);
     if(false==tracers[j].file_exist)
       {
        So.message_warning("File not found. Stop");
        exit(1);
       }
     else
       {
        this->File.read_array(list_file_aux[j], delta_ghost); // FIX NAME OF GHOST DM, I have used one of the calibration
        real_prec nmean=get_nobjects(delta_ghost);
        nmean/=static_cast<real_prec>( this->params._NGRID());
        get_overdens(delta_ghost,nmean, delta_ghost);
//andres this was just a test meant to account for discrepancies in growhts, leave commented
//        for(ULONG i=0;i<delta_ghost.size();++i)delta_ghost[i]*=1.033;
//
        this->Konvolve(delta_ghost, delta_ghost); //Convolve with the new kernel
//        for(ULONG i=0;i<this->Kernel.size();++i)cout<<this->Kernel[i]<<endl;
#if defined _USE_CWC_
        dm_new[j].cwclass_field.get_CWC(delta_ghost);   //get the CW info
#ifdef _USE_MASS_KNOTS_
        dm_new[j].cwclass_field.get_Mk_collapsing_regions(delta_ghost,nmean);  //get the MK info
#endif //use_mass_knots
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_) || defined (_USE_IWEB_) || defined (_USE_AWEB_) || defined (_USE_PWEB_) || defined (_USE_IKWEB_) || defined (_USE_TIWEB_)
        dm_new[j].cwclass_field.get_CWC(delta_ghost);   //get the CW info
#if defined (_USE_MASS_KNOTS_)
        dm_new[j].cwclass_field.get_Mk_collapsing_regions(delta_ghost,nmean);  //get the MK info
#endif // use_mass_knots
#endif // use_mass_knots || other models

#endif    // !use_cwc
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_) || defined (_USE_S3_) || defined (_USE_S2DELTA_) ) && (!defined (_USE_CWC_))
        dm_new[j].cwclass_field.get_bias_terms(delta_ghost);
#endif
        real_prec alpha_c=1.0;
#ifdef _SLICS_
        //        So.message_warning("For the slics, we are using a factor 0.95 to correct the DM delta: due to a overwritting of dm files.");
        alpha_c=1.0;
#endif
        for(ULONG i = 0; i <  this->params._NGRID() ;++i )  //TRANSFORM DELTATO LOG10(NUM_IN_LOG + DELTA)
// Andres, ojo, el 0.95 es a mano. LO he puesto cuando me di cuenta de que las dm de las calibraciones de la slics tenian un
          // espectro 0.954 mayor que las dm que usaba para generar mocks, que probablemente los generé despues de la calibracion con tora cosmologia.
          dm_new[j].delta_field[i] = delta_ghost[i] <-1? 0: log10(NUM_IN_LOG + alpha_c*static_cast<real_prec>(delta_ghost[i]));
       }
     delta_ghost.clear();
     delta_ghost.shrink_to_fit();
    } //closes loop over n_refs
// WE need to do this becuase the get_new_min_mass function expects this->delta_X
   this->delta_X.resize( this->params._NGRID(),0);
   this->delta_X=dm_new[0].delta_field;
   this->get_new_min_max_properties();
   this->delta_X.clear();delta_X.shrink_to_fit();
   So.DONE();
   time_t start_mock;
   time(&start_mock);
#ifdef _FULL_VERBOSE_
  So.message_screen("Assigning Nh to cells");
#endif
   const gsl_rng_type * rng_t;
   gsl_rng * gBaseRand;
   gsl_rng_env_setup();
   gsl_rng_default_seed=35;
   rng_t = gsl_rng_mt19937;//_default;
   gBaseRand = gsl_rng_alloc (rng_t);
  ULONG cell_counter=0;
   for(ULONG ii=0;ii<(n_new_dm)* this->params._NGRID();++ii)   // loop over (nrefs-1)*NGRID: we randomly choose mesh1 or mesh2, and in each mesh, we pick-up a random cell. NOte that we need to do the CWC of the second field as well and plage the inner with ifs()
    {
      int mesh_id=0; // Identify one of the n_refs
      ULONG i=ii;
      if(n_new_dm>1)
        {
         if(ii<n_new_dm)
            mesh_id=ii;
          else
           mesh_id= ii % (n_new_dm) ; // Identify one of the n_refs
//              mesh_id=gsl_rng_uniform_int(gBaseRand,n_new_dm);
          cell_counter=dm_new[mesh_id].mesh_id_counter;
          i=dm_new[mesh_id].id_cells[cell_counter];
      }
      if(true==tracers[mesh_id].file_exist)
      {
      ULONG I_X=0;
      real_prec dm = dm_new[mesh_id].delta_field[i];
      I_X  = get_bin(dm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
         // **********CWT
      ULONG I_CWT=0;
#ifdef _USE_CWC_
       I_CWT=static_cast<ULONG>(dm_new[mesh_id].cwclass_field.get_Tclassification(i));
#endif
        // **********MK
         // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
        ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
        I_MK= (dm_new[mesh_id].cwclass_field.cwt_used[dm_new[mesh_id].cwclass_field.get_Tclassification(i)]== I_KNOT ? dm_new[mesh_id].cwclass_field.SKNOT_M_info[i]: 0);
//        cout<<I_MK<<" "<<I_CWT<<"   "<<dm<<endl;
#endif
         // **********CW-V
        ULONG I_CWV=0;
#ifdef _USE_CWC_V_
         I_CWV=dm_new[mesh_id].cwclass_field.get_Vclassification(i);
#endif
         // **********MV
         ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
         I_VK= (dm_new[mesh_id].cwclass_field.cwv_used[I_CWV]== I_KNOT ? dm_new[mesh_id].cwclass_field.VDISP_KNOT_info[i]: 0);
#endif
         // **********C1
         // Get the corresponding bin in the two invariants of the shear of the tidal field
         ULONG I_C1=0;
         real_prec C1=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_II_
         C1 = dm_new[mesh_id].cwclass_field.Invariant_TF_II[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
         C1 = dm_new[mesh_id].cwclass_field.DELTA2[i];
         I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
         // **********C2
        ULONG I_C2=0;
        real_prec C2=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_III_
        C2 = dm_new[mesh_id].cwclass_field.Invariant_TF_III[i];
        I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_DELTA3_
        C2 = dm_new[mesh_id].cwclass_field.DELTA3[i];
        I_C2 = get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#endif
       ULONG I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
         real_prec C3 = dm_new[mesh_id].cwclass_field.Invariant_TF_IV[i];
        I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined (_USE_TIDAL_ANISOTROPY_)
         real_prec C3 = dm_new[mesh_id].cwclass_field.Tidal_Anisotropy[i];
       I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_S2_
         real_prec C3 = dm_new[mesh_id].cwclass_field.S2[i];             // s²
         I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
         ULONG I_CV1=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_I_
         real_prec CV1 = dm_new[mesh_id].cwclass_field.Invariant_VS_I[i];
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
         real_prec CV1 = dm_new[mesh_id].cwclass_field.N2D[i];      // Nabla² ð
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_I_)
         real_prec CV1 = dm_new[mesh_id].cwclass_field.Invariant_TF_I[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
         I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1,this->s_deltas.prop7, this->bin_accumulate_borders);
#endif
         ULONG I_CV2=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_II_
        real_prec CV2 = dm_new[mesh_id].Invariant_VS_II[i];
        I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2,this->s_deltas.prop8, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
         real_prec CV2 = dm_references[mesh_id].S2DELTA[i];         // s²ð
        I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_II_)
         real_prec CV2 = dm_new[mesh_id].cwclass_field.Invariant_TF_II[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
        I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#endif
         ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
        real_prec CV3 = dm_dm_references[j].cwclass_field[mesh_id].Invariant_VS_III[i];
        I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3,this->s_deltas.prop9, this->bin_accumulate_borders);
#elif defined _USE_S3_
       real_prec CV3 =dm_dm_references[j].cwclass_field[mesh_id].S3[i];                                   // s³
       I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#elif defined (_USE_INVARIANT_PWEB_II_)
        real_prec CV3 = dm_new[mesh_id].cwclass_field.Invariant_TF_III[i]; // When using PWEB, the iarrays contaning th nvariants of the tidal field are used to allocate the invariats of the Pfield
       I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#endif
#ifndef _BIN_ACCUMULATE_
         if(dm >=this->s_mins.prop1 && dm <=this->s_maxs.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
           if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
             if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_) || defined (_USE_INVARIANT_TIDAL_FIELD_IV_)
               if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                 if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                   if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                     if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
#endif
                     {
                           // Get the bin in DM
                           ULONG index_el=index_11d(I_X,I_CWT,I_MK,I_CWV,I_VK,I_C1,I_C2,I_C3,I_CV1,I_CV2,I_CV3,static_cast<ULONG>(this->params._n_cwt()),static_cast<ULONG>(this->params._n_sknot_massbin()),static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()),N_C_BIN1,N_C_BIN2,N_C_BIN3,N_CV_BIN1,N_CV_BIN2,N_CV_BIN3);
                           double aux_h=0;
                           for(int ih=0;ih<ny;++ih)
                            aux_h += static_cast<double>(this->BIAS_NCOUNTS[index_2d(ih,index_el,lenght_bias)]);
                           if(aux_h>0) //This applies in the case of generating mocks. This is likely the reason why parallelization in the mock production is not working, for it works for all cells......
                             {
                               ULONG N_available_positions_in_denbin_cwt_Mk = this->NCELLSperDMBIN[index_el];
                               ULONG N_used_positions_in_denbin_cwt_Mk = NCELLSperDMBIN_now[index_el];
                               bool flag = true;
                               if (N_available_positions_in_denbin_cwt_Mk == N_used_positions_in_denbin_cwt_Mk)
                               flag=false;
                               // If the density bin still has available positions to be assigned, then proceed
                               if(true==flag)
                                {
                                   bool cell_accepted=false;
                                  int halo_prop;
                                   // This density bin as available positions. Proceed until the cell is assigned a value of Nhalo
                                   while (false==cell_accepted)
                                     {
                                       // Throw one value of Nhalos in the range [0, nmax_y_onecell], where  nmax_y_onecell represents the
                                       // maximum number of tracers in one cell, read from the cell number counts.
                                       // Accept this number according to the normalized distribution of Nhalos in the corresponding density bin.
                                       // Everytime one cell is assigned this value of Nhalo, the normalized distribution is updated (below)
                                       // and computed extracting that position already assigned. This will ensure that more probability
                                       // is given to the remaining available positions.
                                       real_prec prob_ncounts=-10.0;
                                       real_prec ran = 10.0;
                                       ULONG index=static_cast<ULONG>(LARGE_NUMBER);
                                       // Probability for number counts
                                       while (prob_ncounts<ran)
                                         {
                                           halo_prop= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins in the Y distirbution (id _DENSITY_)
                                           index  = index_2d(halo_prop,index_el,lenght_bias );
                                           prob_ncounts = static_cast<real_prec>(X_Y_hist_dynamical_normalized[index]);
                                           ran   = gsl_rng_uniform(gBaseRand);
                                        }// closes while(prob_ncounts<ran)
                                       ULONG N_available_positions=this->BIAS_NCOUNTS[index];
                                       ULONG N_used_positions = Ncells_density_bin_new[index];
                                       if(N_available_positions > N_used_positions)
                                         {
                                             tracers[mesh_id].delta_field[i]=static_cast<real_prec>(halo_prop);
                                             tracers[mesh_id].Ntracer_dyn+=tracers[mesh_id].delta_field[i];
                                             tracers[mesh_id].counter_dyn++;
                                         // Claim the current cell as already assigned a value of Nhalos. This breaks the while and continues to the next cell.
                                           cell_accepted=true;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           Ncells_density_bin_new[index]++;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                           NCELLSperDMBIN_now[index_el]++;
                                           // -----------------------------------------------------------------
                                           // In the current Den-N bin, subtract one (i.e, the accepted) in order
                                           // to upgrade the distribution and give more weight to the remaining available positions
                                           // Attention, do not ask whether this is <0, for it is an unsigned long
                                           if(X_Y_hist_dynamical[index] >=1)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                                             X_Y_hist_dynamical[index]--;
                                           else
                                             X_Y_hist_dynamical[index]=0;
                                           // Normalize the dynamical histogram,  in each Den bin, to the maximum.
                                           // The resulting distribution, in this den bin, will be used for the next cell
                                           // in order to accept or reject the new value of Nhalos,
                                           // such that more probability is given to the remaining available positions.
                                           long aux_a=-10000;
                                           long lkk;
                                           for(int j=0; j< ny ;++j)// loop over the number of tracers in a mesh cell
                                           {
                                             long AUX=X_Y_hist_dynamical[index_2d(j,index_el, lenght_bias)];
                                             lkk=max(aux_a, AUX);
                                             aux_a=lkk;
                                           }
                                           for(int j=0;j< ny ;++j)
                                           {
                                             ULONG iin=index_2d(j,index_el, lenght_bias);
                                             ULONG AUX=X_Y_hist_dynamical[iin];
                                             X_Y_hist_dynamical_normalized[iin]= (lkk == 0  ? 0.0 :  static_cast<double>(AUX)/static_cast<double>(lkk));
                                           }
                                         }   // end   if(N_available_positions > N_used_positions) .
                                     }// end  while (false==cell_accepted)
                              } // end of if(true==flag)
                            else    // If the density bin has been already filled, we call the full joint distribution in that bin to assign randomly the value Nhalo according to the other properties
                              {
                                real_prec prob=-1.0;
                                real_prec ran=0.0;
                                int Nhalos_orfan=0;
                                while(prob<ran)
                                     {
                                        Nhalos_orfan= gsl_rng_uniform_int(gBaseRand,ny); //draw number counts in cells or bins of a density-like quantity
                                        prob = this->BIAS_NCOUNTS_normalized[index_2d(Nhalos_orfan,index_el, lenght_bias)];
                                        ran = gsl_rng_uniform(gBaseRand);
                                     }
                                    tracers[mesh_id].delta_field[i]= static_cast<real_prec>(Nhalos_orfan);
                                     tracers[mesh_id].Ntracer_glob+=tracers[mesh_id].delta_field[i];
                                    tracers[mesh_id].counter_orphan++;
                            }// end else
                }// end of if(aux_h>0), available only under  _achtung_inside_
              else
               tracers[mesh_id].delta_field[i]= 0;
      }
    dm_new[mesh_id].mesh_id_counter++;
  }
  }// end loop over cells
#ifdef _FULL_VERBOSE_
 std::cout<<endl;
#endif
  gsl_rng_free (gBaseRand);
  this->So.DONE();
#ifdef _FULL_VERBOSE_
   So.message_time_mock(start_mock);
#endif
  for(int il=0;il<n_new_dm;++il)
    {
      if(tracers[il].file_exist==true)
      {
      tracers[il].get_nobjects_tracers();
 #ifdef _FULL_VERBOSE_
      So.message_screen("Fraction of cells assigned with DYNAMICAL prob. dist =", 100.0*static_cast<double>(tracers[il].counter_dyn)/static_cast<double>( this->params._NGRID()), "%");
      So.message_screen("Fraction of cells assigned with GLOBAL  prob. dist =", 100.0*static_cast<double>(tracers[il].counter_orphan)/static_cast<double>( this->params._NGRID()), "%");
      So.message_screen("Partial number of objects in new mock= ", tracers[il].Nobjects);
      So.message_screen("Mean number density = ", static_cast<real_prec>(tracers[il].Nobjects)/pow(this->params._Lbox(), 3), "(Mpc / h )⁻³");
#endif
      string fool=this->params._Output_directory()+"Ntracers_realization"+to_string(this->params._seed())+".txt";
      ofstream ncc;
      ncc.open(fool.c_str());
      ncc<<tracers[il].Nobjects<<endl;
      ncc.close();
    }
  }
   delta_Y_new.resize( this->params._NGRID(),0);
  for(int il=0;il<n_new_dm;++il)
    {
      if(true==tracers[il].file_exist)
      {
       for(ULONG i=0;i< this->params._NGRID();++i)this->delta_Y_new[i]=tracers[il].delta_field[i];
       this->File.write_array(tracers[il].fname, this->delta_Y_new);
#ifdef _USE_LPT_
       if(case_prop==1 || case_prop==2)
           this->lpt.set_fname_MOCK_NCOUNTS(fname);
#endif
      }
    }
#ifdef _USE_GNUPLOT_
    this->delta_Y.resize( this->params._NGRID(),0);
    this->File.read_array_t<PrecType_Y>(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(), this->delta_Y);
    this->get_power_spectrum("TR_REF");
    this->delta_Y.clear();this->delta_Y.shrink_to_fit();
    this->gp_power<<"set border linewidth 2.\n";
    this->gp_power<<"set size square \n";
    this->gp_power<<"set log x \n";
    this->gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
    this->gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";
    vector<pair<real_prec, real_prec> > xy_pts_ref;
    for(int i=0; i<kvec.size(); ++i)
       xy_pts_ref.push_back(std::make_pair(this->kvec[i], log10(this->Power_REF[i])));
#endif
  int ilc=0;
  for(int il=0;il<n_new_dm;++il)
    {
     if(true==tracers[il].file_exist)
      {
        for(ULONG i=0;i< this->params._NGRID();++i)
        this->delta_Y_new[i]=tracers[il].delta_field[i];
//        this->params.set_file_power(file_aux_power+to_string(list_aux[il])); // this was ment for cases in which various dm are generated at once, but created conflict with names.
        this->params.set_file_power(file_aux_power);
        this->get_power_spectrum("TR_MOCK_REALIZATION");
#ifdef _USE_GNUPLOT_POWER_PLOT_    // PLOT OF THE POWER SPECTRUM
        vector<pair<real_prec, real_prec> > xy_pts_new;
       for(int i=0; i<kvec.size(); ++i)
         xy_pts_new.push_back(std::make_pair(this->kvec[i], log10(this->Power_NEW[i])));
       if(ilc==0){
           this->gp_power<<"plot [][1.5:5]"<<this->gp_power.file1d(xy_pts_ref) << "w l lw 3.3 lt 2 title 'Reference',"<<this->gp_power.file1d(xy_pts_new) << " w l lw 3.3 lt 6 title 'Mock "<<list_aux[il]<<"'"<<endl;
        }
       else{
         this->gp_power<<"replot "<<this->gp_power.file1d(xy_pts_new) << " w l lw 3.3 lt 6 title 'Mock "<<list_aux[il]<<"'"<<endl;
       }
       xy_pts_new.clear();
       xy_pts_new.shrink_to_fit();
#endif
       ilc++;
     }
   }// closes loop over references, opendd after the assignmeng
#ifdef _USE_GNUPLOT_POWER_PLOT_    // PLOT OF THE POWER SPECTRUM
  xy_pts_ref.clear();
  xy_pts_ref.shrink_to_fit();
#endif
  this->delta_Y_new.clear();
  this->delta_Y_new.shrink_to_fit();
   }  // closes loop over chunks
// fre memory after closing the loop over the set of 5-new dm files to sample uopon
  X_Y_hist_dynamical.clear();
  X_Y_hist_dynamical.shrink_to_fit();
  X_Y_hist_dynamical_normalized.clear();
  X_Y_hist_dynamical_normalized.shrink_to_fit();
  this->NCELLSperDMBIN.clear(); this->NCELLSperDMBIN.shrink_to_fit();
  NCELLSperDMBIN_now.clear();    NCELLSperDMBIN_now.shrink_to_fit();
  Ncells_density_bin_new.clear(); Ncells_density_bin_new.shrink_to_fit();
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _GET_BiasMT_CAT_
void BiasMT::assign_tracer_property(bool initial_assignment, string h_property)
{
    So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
  if(initial_assignment==true)
  {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     So.message_screen("*Assigning Vmax using the values read from the reference*");
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
      So.message_screen("*Assigning Halo Mass using the values read from the reference*");
#endif
    }
  So.message_screen("Number of tracers in reference = ", this->tracer_ref._NOBJS());
  So.message_screen("Number of tracers in mock      = ", this->tracer._NOBJS());
#endif
  ULONG cumulative_counter=0; // this counts all
  // *********************************************OMP stuff***********************************************************************
  int NTHREADS=_NTHREADS_;
  gsl_rng_env_setup();
  const gsl_rng_type *Tn;
  gsl_rng_default_seed=1015;
  Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (Tn);
  // ********************************************************************************************************************
  ULONG LENGHTdm=this->dm_properties_bins.size();
  //Vector of structures aimed at keeping the values of the tracer properties (vmax, or halos) which fall in a bin of {Theta}
  //This information will be modified below, so we make a copy ni order not to measure that again above in get_scaling_relations_primary_property()
  // ********************************************************************************************************************
   ULONG counter_masses_0=0;
#ifndef _MULTISCALE_
  vector<ULONG>number_in_theta_ref(LENGHTdm,0);
  // container to track the number of available masses in each theta-bin
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<LENGHTdm; ++i)
    number_in_theta_ref[i]=this->dm_properties_bins[i].tracer_properties.size();
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
  ULONG LENGHTdm_for_randoms=this->dm_properties_for_randoms_bins.size(); //In principe this number is equal to this->dm_properties_bins.size();
  vector<ULONG>number_in_theta_ref_for_randoms(LENGHTdm_for_randoms,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<LENGHTdm_for_randoms; ++i)
    number_in_theta_ref_for_randoms[i]=this->dm_properties_for_randoms_bins[i].tracer_properties.size();
#endif
#endif  // end of ifndef _MULTISCALE_
  // *******************************Mass, Vmax function stufff************************************************************************************
  int nbins_mf=0;
  real_prec lm_min, lm_max;
  gsl_interp_accel *acc;
  gsl_spline *spline ;
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  if(true==initial_assignment)
    nbins_mf=this->tracer_ref.vmax_function.size();
  else 
   {
     if(h_property==_MASS_)
      nbins_mf=this->tracer_ref.mass_function.size();
     else if (h_property==_RS_ )
      nbins_mf=this->tracer_ref.rs_function.size();
     else if ( h_property == _CONCENTRATION_)
      nbins_mf=this->tracer_ref.cvir_function.size();
     else if (h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
      nbins_mf=this->tracer_ref.s_function.size();
}
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
  if(true==initial_assignment)
      nbins_mf=this->tracer_ref.mass_function.size();
  else
  {
   if(h_property==_VMAX_)
     nbins_mf=this->tracer_ref.vmax_function.size();
   else if (h_property==_RS_ )
     nbins_mf=this->tracer_ref.rs_function.size();
   else if ( h_property == _CONCENTRATION_)
     nbins_mf=this->tracer_ref.concentration_function.size();
   else if (h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
     nbins_mf=this->tracer_ref.s_function.size();
  }
#endif // endif for #ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
   if(nbins_mf==0)
     {
       So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"Reference mass function has not be allocated and perhaps not measured. Code exits here.");
       exit(0);
     }
#ifdef _FULL_VERBOSE_
   So.message_screen("Preparing arrays for global X-function from reference: ");
#endif
   this->tracer.define_property_bins();
   this->mfunc.resize(nbins_mf,0);
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     if(true==initial_assignment)
      {
        this->tracer_ref.params.set_i_vmax_g(POSITIVE_INT); //any positive number, allows to define bins in vmax
        this->tracer_ref.define_property_bins();
        this->prop_min=this->tracer_ref.VMAXBin[0];
        this->prop_max=this->tracer_ref.VMAXBin[nbins_mf-1];
     }
     else if(false==initial_assignment)
      {
        if(h_property==_MASS_)
         {
           this->tracer_ref.params.set_i_mass_g(POSITIVE_INT); //any positive number, allows to define bins in mass
           this->tracer_ref.define_property_bins();
           this->prop_min=this->tracer_ref.MBin[0];
           this->prop_max=this->tracer_ref.MBin[nbins_mf-1];
         }
     }
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
   if(true==initial_assignment)
    {
      this->tracer_ref.params.set_i_mass_g(POSITIVE_INT); //any positive number, allows to define bins in mass
      this->tracer_ref.params.set_i_vmax_g(NEGATIVE_INT); //any positive number, allows to define bins in mass
      this->tracer_ref.define_property_bins();
      this->prop_min=this->tracer_ref.MBin[0];
      this->prop_max=this->tracer_ref.MBin[nbins_mf-1];
   }
   else if(false==initial_assignment)
   {
      if(h_property==_VMAX_)
       {
         this->tracer_ref.params.set_i_vmax_g(POSITIVE_INT); //any positive number, allows to define bins in vmax
         this->tracer_ref.params.set_i_mass_g(NEGATIVE_INT); //any positive number, allows to define bins in mass
         this->tracer_ref.define_property_bins();
         this->prop_min=pow(10,this->tracer_ref.VMAXBin[0]);
         this->prop_max=pow(10,this->tracer_ref.VMAXBin[nbins_mf-1]);
       }
   }
#endif
   if(false==initial_assignment)
   {
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
     if(h_property==_RS_ )
       {
          this->tracer_ref.params.set_i_rs_g(POSITIVE_INT); //any positive number, allows to define bins in rs
          this->tracer_ref.define_property_bins();
          this->prop_min=this->tracer_ref.RSBin[0];
          this->prop_max=this->tracer_ref.RSBin[nbins_mf-1];
        }
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
     if( h_property == _CONCENTRATION_)
       {
          this->tracer_ref.define_property_bins();
          this->prop_min=this->tracer_ref.CONCENTRATIONBin[0];
          this->prop_max=this->tracer_ref.CONCENTRATIONBin[nbins_mf-1];
       }
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
      else if (h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
        {
          if (h_property==_SPIN_)
          { 
            this->tracer_ref.params.set_i_spin_g(POSITIVE_INT);//any positive number, allows to define bins in spin
            this->tracer_ref.params.set_i_spin_bullock_g (NEGATIVE_INT);//any positive number, allows to define bins in spin
          }
          if (h_property==_SPIN_BULLOCK_)
          {
            this->tracer_ref.params.set_i_spin_g(NEGATIVE_INT);//any positive number, allows to define bins in spin
            this->tracer_ref.params.set_i_spin_bullock_g(POSITIVE_INT);//any positive number, allows to define bins in spin
          }
          this->prop_min=this->tracer_ref.SPINBin[0];
          this->tracer_ref.define_property_bins();
          this->prop_max=this->tracer_ref.SPINBin[nbins_mf-1];
//          this->prop_min= this->params._SPINmin();//  this->tracer_ref.SPINBin[0];
//          this->prop_max= this->params._SPINmax(); // this->tracer_ref.SPINBin[nbins_mf-1];

          cout<<RED<<this->prop_min<<" "<<this->prop_max<<endl;
        }
#endif
      }
      // I cannot use here the this->tracer MBmin and max for it might happen that the reference mass function has been measured with a different number of bins
      // different to the current this>tracer_ref_NMBINS
      vector<real_prec>delta_prop_aux(nbins_mf,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0; i< nbins_mf; ++i)
        {
          real_prec lmmin=log10(this->prop_min)+i*log10(this->prop_max/this->prop_min)/static_cast<real_prec>(nbins_mf);
          real_prec lmmax=log10(this->prop_min)+(i+1)*log10(this->prop_max/this->prop_min)/static_cast<real_prec>(nbins_mf);
          delta_prop_aux[i]=pow(10,lmmax)-pow(10,lmmin);
        }
#ifdef _USE_GNUPLOT_ABUNDANCE_V_PLOT_
   vector<pair<real_prec, real_prec> > xy_pts_ref_v;
#endif
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
   vector<pair<real_prec, real_prec> > xy_pts_ref_m;
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      if(true==initial_assignment)
      {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0; i< nbins_mf; ++i)
            this->mfunc[i]=this->tracer_ref.vmax_function[i]*delta_prop_aux[i];
        }
        else{
         if(h_property==_MASS_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0; i< nbins_mf; ++i)
           this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_prop_aux[i];
        }
#ifdef _USE_GNUPLOT_ABUNDANCE_V_PLOT_
   xy_pts_ref_v.clear();xy_pts_ref_v.shrink_to_fit();
   for(int i=0; i<this->tracer_ref.vmax_function.size(); ++i)
     xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.VMAXBin[i]+(i+0.5)*(this->tracer_ref.VMAXBin[i]-this->tracer_ref.VMAXBin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.vmax_function[i])));
#endif

#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
   if(true==initial_assignment)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i< nbins_mf; ++i)
       this->mfunc[i]=this->tracer_ref.mass_function[i]*delta_prop_aux[i];
     }
   else if(false==initial_assignment)
     if(h_property==_VMAX_)
       {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0; i< nbins_mf; ++i)
          this->mfunc[i]=this->tracer_ref.vmax_function[i]*delta_prop_aux[i];
     }
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
   xy_pts_ref_m.clear();xy_pts_ref_m.shrink_to_fit();
   for(ULONG i=0; i<this->tracer_ref.mass_function.size(); ++i)
       xy_pts_ref_m.push_back(std::make_pair(this->tracer_ref.MBmin[i]+(i+0.5)*(this->tracer_ref.MBmax[i]-this->tracer_ref.MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.mass_function[i])));
#endif
#endif  // end ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_ line 8503
   if(false==initial_assignment)// now do rs and spin
    {
       if(h_property==_RS_ )
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(int i=0; i< nbins_mf; ++i)
           {
            this->mfunc[i]=this->tracer_ref.rs_function[i]*delta_prop_aux[i];
           }
         }
       if(h_property == _CONCENTRATION_)
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(int i=0; i< nbins_mf; ++i)
           {
            this->mfunc[i]=this->tracer_ref.cvir_function[i]*delta_prop_aux[i];
           }
         }



         else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_ )
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0; i< nbins_mf; ++i)
          this->mfunc[i]=this->tracer_ref.s_function[i]*delta_prop_aux[i];
        }
      }
      delta_prop_aux.clear(); delta_prop_aux.shrink_to_fit();
      real_prec max_mf=static_cast<real_prec>(get_max(mfunc));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i< nbins_mf; ++i)
        {
          real_prec mf=this->mfunc[i];
          this->mfunc[i]=mf/max_mf;
        }
      acc = gsl_interp_accel_alloc ();
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    if(true==initial_assignment) // for true, we initially assign the primary property
      {
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.vmax_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.VMAXBin[0]), &(this->mfunc[0]), this->tracer_ref.vmax_function.size());
        lm_min=log10(this->params._VMAXmin());
        lm_max=log10(this->params._VMAXmax());
      }
      else
        if(h_property==_MASS_) // if false, we assign the secondary
         {
           spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
           gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
           lm_min=this->params._LOGMASSmin();
           lm_max=this->params._LOGMASSmax();
        }
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
      if(true==initial_assignment) // for true, we initially assign the primary property
        {
            spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.mass_function.size());
            gsl_spline_init (spline, &(this->tracer_ref.MBin[0]), &(this->mfunc[0]), this->tracer_ref.mass_function.size());
            lm_min=this->params._LOGMASSmin();
            lm_max=this->params._LOGMASSmax();
        }
      else if (false==initial_assignment)
        if(h_property==_VMAX_)// if false, we assign the secondary
         {
           spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.vmax_function.size());
           gsl_spline_init (spline, &(this->tracer_ref.VMAXBin[0]), &(this->mfunc[0]), this->tracer_ref.vmax_function.size());
           lm_min=log10(this->params._VMAXmin());
           lm_max=log10(this->params._VMAXmax());
        }
#endif //endif _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    if(false==initial_assignment)
    {
      if(h_property==_RS_ )
       {
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.rs_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.RSBin[0]), &(this->mfunc[0]), this->tracer_ref.rs_function.size());
        lm_min=log10(this->params._RSmin());
        lm_max=log10(this->params._RSmax());
       }
      if( h_property == _CONCENTRATION_)
       {
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.cvir_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.CONCENTRATIONBin[0]), &(this->mfunc[0]), this->tracer_ref.cvir_function.size());
        lm_min=log10(this->params._CONCENTRATIONmin());
        lm_max=log10(this->params._CONCENTRATIONmax());
       }
      else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_ )
       {
        spline    = gsl_spline_alloc (gsl_interp_linear, this->tracer_ref.s_function.size());
        gsl_spline_init (spline, &(this->tracer_ref.SPINBin[0]), &(this->mfunc[0]), this->tracer_ref.s_function.size());
        lm_min=log10(this->params._SPINmin());
        lm_max=log10(this->params._SPINmax());
       }
    }
  // ************************************************RENAME SOME BIN PROPPERTIES*****************************************
  ULONG N_a=N_C_BIN1; //THis will be used by std_sep_in cell as well
  if(true==initial_assignment) // This applies for the assignment of Vmax with the different techniques
    {
#ifdef _USE_TOTAL_MASS_IN_CELL_  //check this, things need to be fixed inside, related to paths
      vector<real_prec> MOCK_MASS_FIELD;
      N_a = N_BINS_TOTAL_MASS_IN_CELL;
      // Can happen that the used mock properties fall in a theta bin in which the reference has nothing, hence  dm_properties_bins[index_bins].tracer_properties.size()=0
      // This has been visible when using a total mass density field sampled from the calibration (as it must be) instead of using the reference
      // This is likely to happen if we do not use the sampled total mass field that correspond to the final bumber density field that one uses to assign positions
      // Hence, if we are assigning as test_mode to the reference, use the mass density field from the reference
      MOCK_MASS_FIELD.resize( this->params._NGRID(),0);
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
      So.message_warning_ini(__LINE__, __PRETTY_FUNCTION__, __FILE__,"I am passing the mass density field from the reference, in order toto do tests. Ideally this should be the one sampled in the last step of the calibration");
      this->tracer_ref.get_density_field_grid(_MASS_, MOCK_MASS_FIELD); // ojo que estoy cargando la ref, no el mock:esto debería ser el mock con la masa total dada por la calibración
#else
      File.read_array(this->params._Output_directory()+"MOCK_TR_MASS_iteration200_MASY0_Nft256_z1.124.dat",MOCK_MASS_FIELD);
#endif

      real_prec mmax=get_max<real_prec>(MOCK_MASS_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum mass of tracer in mock-cell", mmax);
#endif
      mmax=get_min<real_prec>(MOCK_MASS_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Minimum mass of tracer in mock-cell", mmax);
#endif
#endif //end totalmass in cell
    }
  // ************************************************RENAME SOME BIN PROPPERTIES********************************************************************
  ULONG N_b=N_C_BIN2;
  real_prec delta_min_sep =0;
  if(true==initial_assignment)
    {
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
      this->tracer.get_neighbour_tracers(this->ncells_info);
      N_b = N_NEIGHBOURS_MAX;
#elif defined (_USE_MIN_DISTANCE_TO_NEIGHBOURS_) || defined (_USE_LOCAL_CLUSTERING_)
      N_b=N_BINS_MIN_DIST_TO_NEI;
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_
     N_b= N_BINS_MIN_SEP_IN_CELLS;
     this->tracer.get_stats_separation_in_cell();
     delta_min_sep = (MAX_SEP_IN_CELLS-MIN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MIN_SEP_IN_CELLS);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
     N_b= N_BINS_MEAN_SEP_IN_CELLS;
     this->tracer.get_stats_separation_in_cell(); // this also works to get the mean and stedev
     delta_min_sep = DELTA_MEAN_SEP; //(MAX_MEAN_SEP_IN_CELLS-MIN_MEAN_SEP_IN_CELLS)/static_cast<real_prec>(N_BINS_MEAN_SEP_IN_CELLS);
#endif
    }
  // ******************************************************************************
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
  N_a=  N_BINS_STDV_SEP_IN_CELLS ;
#if !defined _USE_MEAN_SEPARATIONS_IN_CELLS_
    this->tracer.get_stats_separation_in_cell(); // this also works to get the mean and stedev
#endif
  real_prec delta_stdv_sep = DELTA_STDV_SEP;// (MAX_STDV_SEP_IN_CELLS-MIN_STDV_SEP_IN_CELLS)/static_cast<real_prec>(N_a);
#endif
  // ******************************************************************************
  ULONG N_c = 1; // this applies to tidal anisotropy, inv tidal field iv, s2
#ifdef _USE_TIDAL_ANISOTROPY_
  N_c = N_C_BIN3;
#elif defined _USE_MACH_NUMBER_
  N_c = N_BINS_MACH;
#endif
  // ******************************************************************************
  ULONG N_v= N_CV_BIN1;
#ifdef _USE_TRACERS_IN_CELLS_
  vector<real_prec> MOCK_DEN_FIELD;
#endif
  if(true==initial_assignment)  // Ww shall not use his info for post-mass assignment
    {
#ifdef _USE_TRACERS_IN_CELLS_
      N_v = N_BINS_TRACERS_IN_CELLS;
      MOCK_DEN_FIELD.resize( this->params._NGRID(),0);

#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
      this->File.read_array(this->fnameMOCK,MOCK_DEN_FIELD); // fname mock has been defined in bamrunner before calling makecat
#else
      this->File.read_array(this->fnameMOCK+".dat",MOCK_DEN_FIELD); // fname mock has been defined in bamrunner before calling makecat
#endif

      int nmax=get_max<real_prec>(MOCK_DEN_FIELD);
#ifdef _FULL_VERBOSE_
      So.message_screen("Maximum number of tracer in cell", nmax);
#endif
#endif // end USE_TRACERS_IN_CELLS_
    }
  vector<real_prec> CROSS_CORR;
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
  if(true==initial_assignment)
  {  // only applies for the asignment of first primary property
    Params params_new=this->params;
    vector<real_prec> NEW_DM( this->params._NGRID(),0);
    this->File.read_array(this->params._Input_Directory_X_REF()+this->params._Name_Catalog_X(), NEW_DM);
    CROSS_CORR.resize(this->params._NGRID(),0);
    PowerSpectrumF cpower(params_new);
#ifndef  _USE_TRACERS_IN_CELLS_
    vector<real_prec> MOCK_DEN_FIELD( this->params._NGRID(),0);
    this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),MOCK_DEN_FIELD);
#endif
    cpower.get_cross_correlation_config_space(NEW_DM,MOCK_DEN_FIELD,CROSS_CORR);
    NEW_DM.clear(); NEW_DM.shrink_to_fit();
  }
#endif
   ULONG N_x=1;
    if(true==initial_assignment)
#ifdef _USE_STATISTICAL_ASSIGNMENT_
          N_x = this->params._NPROPbins_bam();
#endif
#if defined (_USE_VMAX_AS_PRIMARY_OBSERVABLE_) || defined (_USE_MASS_AS_PRIMARY_OBSERVABLE_)
#ifdef _ASSIGN_MASS_POST_
  N_x=1;
  if(false==initial_assignment) // if initial_assignem t is false, we proceed to assign mass/vmax once the vmax/mass is already assigned. This will be false upon a second call of this function.
    N_x=this->params._NPROPbins_bam();
  else if(true==initial_assignment)
    {
      N_x=1;
#ifdef _USE_LOCAL_OVERDENSITY_
      N_x=N_BINS_LO;
#endif
    }
#elif defined _ASSIGN_VMAX_POST_
  if(false==initial_assignment) // if initial_assignem t is false, we proceed to assign mass/vmax once the vmax/mass is already assigned. This will be false upon a second call of this function.
    N_x = this->params._NPROPbins_bam();
#endif
#endif

  ULONG counter_fmf=0;
#if defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_ || defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
  if(true==initial_assignment) // This applies for the assignment of vmax or mass, reading from the reference catalogs, whether we are using or not multi-scale approach
    {
#endif
#ifdef _FULL_VERBOSE_
      So.message_screen("ASSIGNMENT, in function ",__PRETTY_FUNCTION__);
#endif
#if defined (_MULTISCALE_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_) 
// This container of structures can in principle hold the same info as the vector tracer[].galThetaID.
// but since it is useful when loops are performed over the grid (where we tracer[].galThetaID is not accesible)
// we use this:
    vector<s_cell_info> cell_info_tr( this->params._NGRID());
#endif
    // ******************************************************************************
   ULONG N_bias=1;
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
   N_bias=N_BINS_BIAS;
#endif
// ************************************************************************************
//  The following loop aims at collecting the values of Vmax in theta bins.
  counter_fmf=0;
  ULONG count_dm=0;
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
  ULONG count_dm_full=0;
  ULONG count_ran_full=0;
#endif
  // Note:
  // If the number of neighbours is to be used, we do a loop over the particles.
  // Also, when multiscale *is not* used, we do a loop over the particles.
  // Else (e.g. using multiscale or using neighbour info), we do a loop over the grid cells, since the some of
  // the used quantities are computed at the cells.
  // In both cases, parallelization can be done, taking care of defining the random objects inside parallel region
  // Start parallel region. Problems with this paralellization
#ifdef _USE_OMP_TEST_
  int jthread=0;
  const gsl_rng_type *Trn;
  gsl_rng *randa;
  vector<ULONG>vseeds(NTHREADS,0);
  for(int i=0;i<vseeds.size();++i)
   vseeds[i]=35+static_cast<ULONG>(i+14)*56045;
#ifdef _MULTISCALE_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Collecting info for multiscaling");
#endif
#elif defined _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Assinging with threhold in Vmax");
#endif
#endif
#pragma omp parallel private (jthread, Trn, randa)// this is causing problem
  {
    jthread=omp_get_thread_num();
    gsl_rng_default_seed=vseeds[jthread];
    Trn = gsl_rng_mt19937;//_default;
    randa= gsl_rng_alloc (Trn);
#endif
     // ********************************************************************************
#ifdef _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
    //If we are assigning properties to the same reference catalog, we use the same bias previously computed in label bias_pr, so we do not need to compute here agan
    // Or, if we are to assign coords to a new reference, we proceed to compute properties from that new target reference
#ifndef _READ_BIAS_
      if(true==this->params._Get_tracer_bias())
      {
        vector<real_prec> DM_DEN_FIELD( this->params._NGRID(),0);// Read the reference DM field
        this->File.read_array(this->params._Input_Directory_X_new_ref()+this->params._Name_Catalog_X_new_ref(),DM_DEN_FIELD);
        So.message_screen("\tMin DM field:",get_min_nm(DM_DEN_FIELD));
        So.message_screen("\tMax DM field:",get_max_nm(DM_DEN_FIELD));
        So.message_screen("\tMean DM field:",get_mean(DM_DEN_FIELD));
        PowerSpectrumF Psb(this->params, true);
        Psb.object_by_object_bias(this->tracer.Halo,DM_DEN_FIELD);
        DM_DEN_FIELD.clear();DM_DEN_FIELD.shrink_to_fit();
      }
#ifdef _USE_MACH_NUMBER_
      this->tracer.get_local_mach_number_chuncks(this->params._Scale_mach_number());
#endif
#ifdef _WRITE_BIAS_
      string out_bias=params._Output_directory()+"individual_bias_invPhase.txt";
      ofstream bout; bout.open(out_bias.c_str());
      this->So.message_screen("Writing bias in  file ", out_bias);
      for(ULONG i=0;i<this->tracer.Halo.size();++i)
        bout<<this->tracer.Halo[i].bias<<"\t"<<this->tracer.Halo[i].mach_number<<"\t"<<this->tracer.Halo[i].local_overdensity<<endl;
      bout.close();
#endif
#endif
#ifdef _READ_BIAS_
      string out_bias=params._Output_directory()+"individual_bias_invPhase.txt";
      ifstream bout; bout.open(out_bias.c_str());
      this->So.message_screen("Reading bias from  file ", out_bias);
      for(ULONG i=0;i<this->tracer.Halo.size();++i)
        bout>>this->tracer.Halo[i].bias>>this->tracer.Halo[i].mach_number>>this->tracer.Halo[i].local_overdensity;
      bout.close();
#endif 
#ifdef _USE_TIDAL_ANISOTROPY_SEC_PROP_
      vector<real_prec>tidal(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->params._NGRID() ;++i)
        tidal[i]=tidal_anisotropy(this->cwclass.lambda1[i], this->cwclass.lambda2[i], this->cwclass.lambda3[i]);
      this->tracer.get_tracer_tidal_anisotropy(tidal);
      So.message_screen("\tMin TA :",this->tracer._min_tidal_anisotropy());
      So.message_screen("\tMax TA :",this->tracer._max_tidal_anisotropy());
      tidal.clear();tidal.shrink_to_fit();
#endif
#endif // endf of _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
     // ********************************************************************************
#ifdef _USE_HYBRID_ASSIGNMENT_NEW_
    this->dm_properties_bins_mock.resize(this->dm_properties_bins.size());
#endif
     // ********************************************************************************
#ifdef _MULTISCALE_
#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:count_dm)
#endif
#if defined _USE_NUMBER_OF_NEIGHBOURS_ || defined _ASSIGN_MASS_POST_
    for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
#else
    for(ULONG id=0; id <  this->params._NGRID();++id)
#endif
#else   //else for ifdef _MULTISCALE_
#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:counter_fmf, count_dm)
#endif
    for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig) // Loop enabled if no multiscale is used. WE assign directly here
#endif    // endif for else to def _MULTISCALE_
        {
       // Get the cell ID where the mock tracer is located
         ULONG id=this->tracer.Halo[ig].GridID;
         ULONG I_X=0;
       // Get the bin in the delta dark matter (or log 1+delta) in each cell
#ifdef _USE_DM_DENSITY_AT_HALO_POSITION_
        real_prec xdm=linInterpol(this->params._Nft(),this->params._Lbox(),this->params._d1(),this->tracer.Halo[ig].coord1,this->tracer.Halo[ig].coord2,this->tracer.Halo[ig].coord3,this->delta_X);
        I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
#else
        real_prec xdm = static_cast<real_prec>(this->delta_X[id]);
        I_X  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
#endif
#if defined (_MULTISCALE_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
         real_prec xdm_b = static_cast<real_prec>(this->delta_X[id]);
         cell_info_tr[id].local_dm=pow(10, xdm_b)-NUM_IN_LOG; // let us write delta_dm
#endif
         ULONG I_CWT=0;
#if defined _USE_CWC_ || defined (_USE_MASS_KNOTS_)
         I_CWT=this->cwclass.get_Tclassification(id);
#if defined (_MULTISCALE_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
         cell_info_tr[id].cwt=this->cwclass.cwt_used[I_CWT];
#endif
#endif
// Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
         ULONG I_MK=0;
#ifdef _USE_MASS_KNOTS_
         I_MK= (this->cwclass.cwt_used[I_CWT]== I_KNOT ? this->cwclass.SKNOT_M_info[id]: 0);
#endif
         ULONG I_CWV=0;
#ifdef _USE_CWC_V_
         I_CWV=this->cwclass.get_Vclassification(id);
#endif
       // Get the corresponding bin on Knot-mass. The current cell is in a cluster with a given mass falling in that particular bin
        ULONG I_VK=0;
#ifdef _USE_VEL_KNOTS_V_
         I_VK= (this->cwclass.cwt_used[I_CWV]== I_KNOT ? this->cwclass.VDISP_KNOT_info[id]: 0);
#endif
       // Get the corresponding bin in the two invariants of the shear of the tidal field
        ULONG I_C1=0;
#ifdef _USE_TOTAL_MASS_IN_CELL_
        if(true==initial_assignment)
          I_C1=get_bin(log10(MOCK_MASS_FIELD[id]),this->params._LOGMASSmin(),N_BINS_TOTAL_MAS  S_IN_CELL,(log10(MAX_TOTAL_MASS_IN_CELL)-this->params._LOGMASSmin())/static_cast<double>(N_BINS_TOTAL_MASS_IN_CELL), true);
#elif defined _USE_INVARIANT_TIDAL_FIELD_II_
        real_prec C1 = this->cwclass.Invariant_TF_II[id];
        I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#elif defined _USE_DELTA2_
       real_prec C1 = this->cwclass.DELTA2[id];
       I_C1= get_bin(C1, this->s_mins.prop4, N_C_BIN1, s_deltas.prop4,this->bin_accumulate_borders);
#endif
#if defined _USE_STDV_SEPARATIONS_IN_CELLS_ && !defined  (_USE_TOTAL_MASS_IN_CELL_)
        if(true==initial_assignment)
         {
           real_prec sep_std = this->tracer.stdv_separation_in_cell[id];
           I_C1 = get_bin(sep_std, MIN_STDV_SEP_IN_CELLS, N_a, delta_stdv_sep , this->bin_accumulate_borders);
        }
#endif
       ULONG I_C2=0;
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
       if(true==initial_assignment)
         {
           I_C2=this->tracer.Number_of_neighbours[ig];
           if(I_C2>=N_NEIGHBOURS_MAX)
              I_C2=N_NEIGHBOURS_MAX-1;
          }
#elif defined _USE_INVARIANT_TIDAL_FIELD_III_
          real_prec C2 = this->cwclass.Invariant_TF_III[id];
          I_C2= get_bin(C2, this->s_mins.prop5, N_C_BIN2, s_deltas.prop5,this->bin_accumulate_borders);
#elif defined _USE_MIN_SEPARATIONS_IN_CELLS_
          if(true==initial_assignment)
           {
             real_prec min_sep = this->tracer.min_separation_in_cell[id];
             I_C2 = get_bin(min_sep, MIN_SEP_IN_CELLS, N_b, delta_min_sep , this->bin_accumulate_borders);
           }
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
          if(true==initial_assignment)
           {
             real_prec min_sep = this->tracer.mean_separation_in_cell[id];
             I_C2 = get_bin(min_sep, MIN_MEAN_SEP_IN_CELLS, N_b, delta_min_sep , this->bin_accumulate_borders);
           }
#endif
          ULONG I_C3=0;
#ifdef _USE_INVARIANT_TIDAL_FIELD_IV_
          real_prec C3 = this->cwclass.Invariant_TF_IV[id];
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_TIDAL_ANISOTROPY_

#ifdef  _USE_ANISOTROPY_AT_TRACER_POSITION_
          real_prec C3 = linInterpol(this->params._Nft(),this->params._Lbox(),this->params._d1(),this->tracer.Halo[ig].coord1,this->tracer.Halo[ig].coord2,this->tracer.Halo[ig].coord3,this->cwclass.Tidal_Anisotropy);
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#else
          real_prec C3 = this->cwclass.Tidal_Anisotropy[id];
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#endif
#elif defined _USE_S2_
          real_prec C3 = this->cwclass.S2[id];             // s²
          I_C3= get_bin(C3, this->s_mins.prop6, N_C_BIN3, s_deltas.prop6,this->bin_accumulate_borders);
#elif defined _USE_MACH_NUMBER_
          I_C3= get_bin(this->tracer.Halo[ig].mach_number, this->tracer_ref._min_mach(), N_BINS_MACH, (this->tracer_ref._max_mach()-this->tracer_ref._min_mach())/static_cast<real_prec>(N_BINS_MACH),this->bin_accumulate_borders);
#endif
          ULONG I_CV1=0;
#ifdef _USE_TRACERS_IN_CELLS_
          if(true==initial_assignment)
              I_CV1=get_bin(static_cast<int>(MOCK_DEN_FIELD[id]),0,N_BINS_TRACERS_IN_CELLS,DELTA_TRACERS_IN_CELLS,true);
#elif defined (_USE_INVARIANT_SHEAR_VFIELD_I_)
          real_prec CV1 = invariant_field_I(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]); // Not yet assigned to container
          I_CV1= get_bin(CV1, CV1_MIN, N_CV_BIN1,DELTA_CV1, this->bin_accumulate_borders);
#elif defined _USE_NABLA2DELTA_
          real_prec CV1 = this->cwclass.N2D[id];      // Nabla² ð
          I_CV1= get_bin(CV1, this->s_mins.prop7, N_CV_BIN1, s_deltas.prop7,this->bin_accumulate_borders);
#endif
         ULONG I_CV2=0;
//#ifdef _ASSIGN_MASS_POST_
         if(false==initial_assignment)
         {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
             I_CV2=get_bin(log10(this->tracer.Halo[ig].vmax), log10(this->params._VMAXmin()),N_x,log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(N_x), true);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
             I_CV2=get_bin(log10(this->tracer.Halo[ig].mass), this->params._LOGMASSmin(),N_x,(this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(N_x), true);
#endif
         }
         else{
#if defined _USE_INVARIANT_SHEAR_VFIELD_II_
         real_prec CV2 = invariant_field_II(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]);// Not yet assigned to container
         I_CV2= get_bin(CV2, CV2_MIN, N_CV_BIN2,DELTA_CV2, this->bin_accumulate_borders);
#elif defined _USE_S2DELTA_
         real_prec CV2 = this->cwclass.S2DELTA[id];         // s²ð
         I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
#elif defined _USE_CROSS_CORRELATION_CONF_SPACE_
//      real_prec CV2 = CROSS_CORR[id];         // cross correlation
      real_prec CV2  =linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->tracer.Halo[ig].coord1,this->tracer.Halo[ig].coord2,this->tracer.Halo[ig].coord3,CROSS_CORR);
      I_CV2= get_bin(CV2, this->s_mins.prop8, N_CV_BIN2, s_deltas.prop8,this->bin_accumulate_borders);
      //#endif // endif for _ASSIGN_MASS_POST_
#elif defined _USE_LOCAL_OVERDENSITY_
      I_CV2= get_bin(this->tracer.Halo[ig].local_overdensity, this->tracer_ref._min_local_overdensity(), N_BINS_LO,  (this->tracer_ref._max_local_overdensity()-this->tracer_ref._min_local_overdensity())/static_cast<real_prec>(N_BINS_BIAS) ,this->bin_accumulate_borders);
#endif
         }
        ULONG I_CV3=0;
#ifdef _USE_INVARIANT_SHEAR_VFIELD_III_
        real_prec CV3 = invariant_field_III(this->cwclass.lambda1_vs[id],this->cwclass.lambda2_vs[id],this->cwclass.lambda3_vs[id]);// Not yet assigned to container
        I_CV3= get_bin(CV3, CV3_MIN, N_CV_BIN3,DELTA_CV3, this->bin_accumulate_borders);
#elif defined _USE_S3_
        real_prec CV3 = this->cwclass.S3[id];                                   // s³
        I_CV3= get_bin(CV3, this->s_mins.prop9, N_CV_BIN3, s_deltas.prop9,this->bin_accumulate_borders);
#elif defined _USE_BIAS_OBJECT_TO_OBJECT_ // we use the limits computed from the reference
        I_CV3 = get_bin(this->tracer.Halo[ig].bias, this->tracer_ref._min_bias(),N_BINS_BIAS,(this->tracer_ref._max_bias()-this->tracer_ref._min_bias())/static_cast<real_prec>(N_BINS_BIAS),this->bin_accumulate_borders);
#endif
#ifndef _BIN_ACCUMULATE_
        if(xdm >=this->s_mins.prop1 && xdm <=this->s_maxs.prop1)
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_DELTA2_)
          if(C1<= this->s_maxs.prop4 && C1>= this->s_mins.prop4)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_TIDAL_FIELD_III_) || defined (_USE_DELTA3_)
            if(C2<= this->s_maxs.prop5 && C2>= this->s_mins.prop5)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_TIDAL_ANISOTROPY_) || defined (_USE_S2_) || (_USE_INVARIANT_TIDAL_FIELD_I_)
              if(C3<= this->s_maxs.prop6 && C3>= this->s_mins.prop6)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_I_) || defined (_USE_NABLA2DELTA_)
                if(CV1<= this->s_maxs.prop7 && CV1>= this->s_mins.prop7)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_II_) || defined (_USE_S2DELTA_)
                  if(CV2<= this->s_maxs.prop8 && CV2>= this->s_mins.prop8)
#endif
#endif
#ifndef _BIN_ACCUMULATE_
#if defined (_USE_INVARIANT_SHEAR_VFIELD_III_) || defined (_USE_S3_)
                   if(CV3<= this->s_maxs.prop9 && CV3>= this->s_mins.prop9)
#endif
                     {
#endif
                       ULONG index_bins = index_11d(I_X, I_CWT, I_MK,I_CWV,I_VK,I_C1,I_C2, I_C3, I_CV1,I_CV2, I_CV3, static_cast<ULONG>(this->params._n_cwt()), static_cast<ULONG>(this->params._n_sknot_massbin()),static_cast<ULONG>(this->params._n_cwv()), static_cast<ULONG>(this->params._n_vknot_massbin()), N_a, N_b, N_c,N_v,N_x,N_bias);
#if defined (_MULTISCALE_) || defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_)
                       cell_info_tr[id].Theta_bin=index_bins;// New_Theta-bin in which the cell ID has been identified. This is used below in multiscale DEDICATED SECTION
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                       if(this->tracer.Halo[ig].identity >0 )  // Select tracer associated to a DM particle
                                           {
                           count_dm_full++;
#endif
#ifndef _MULTISCALE_
#ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifdef  _USE_GLOBAL_MASS_FUNCTION_
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
                                             int N_props_left=0;
#endif //endif for _MASS_ASSIGNMENT_TO_REFERENCE_
#else // else for  _USE_GLOBAL_MASS_FUNCTION_
                                   ULONG N_tracers_in_theta_bin = this->dm_properties_bins[index_bins].tracer_properties.size();
                                   ULONG N_tracers_left =  number_in_theta_ref[index_bins]; // this container is updated below
#endif // endif for _USE_GLOBAL_MASS_FUNCTION_
      // If we have available properties in the theta-bin to assign, proceed
                                    if(N_tracers_left>0) //when we use the ref as a mock, this condition will be always satisfied by construction
                                     {
                                       bool flag=false;
                                       while(false == flag) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                                        {
#ifdef _USE_OMP_TEST_
                                           int i_halo_index= gsl_rng_uniform_int(randa,N_tracers_in_theta_bin);
#else
                                      ULONG i_halo_index= gsl_rng_uniform_int(rn,N_tracers_in_theta_bin);
#endif
                                      bool used_tracer = this->dm_properties_bins[index_bins].used_property[i_halo_index]; //true or false if the mass was already chosen or not
                                       if(false == used_tracer)// if the property of the tracer (e.g., vmax, mass) has not been used before, then assign that property to the current particle i
                                        {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                                           this->tracer.Halo[ig].vmax= this->dm_properties_bins[index_bins].tracer_properties[i_halo_index]; //prop_to_assign;
#ifdef  _USE_HYBRID_ASSIGNMENT_NEW_              //fill the structure with the new assigment in order to make the swap procedure below
                                           this->dm_properties_bins_mock[index_bins].tracer_properties.push_back(this->tracer.Halo[ig].vmax);
                                           // keep track of the id  of the tracers
                                           this->dm_properties_bins_mock[index_bins].GalID_bin_properties.push_back(ig);
#endif
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
                                          this->tracer.Halo[ig].mass= this->dm_properties_bins[index_bins].tracer_properties_secondary[i_halo_index]; //prop_to_assign;
                                          this->tracer.Halo[ig].observed_secondary=true; // Mark this tarcer as already assigned a property. This will be used below.
#endif // end
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                                          this->tracer.Halo[ig].mass= this->dm_properties_bins[index_bins].tracer_properties[i_halo_index];  //prop_to_assign;
#ifdef  _USE_HYBRID_ASSIGNMENT_NEW_              //fill the structure with the new assigment in order to make the swap procedure below
                                           this->dm_properties_bins_mock[index_bins].tracer_properties.push_back(this->tracer.Halo[ig].mass);
                                           // keep track of the id  of the tracers
                                           this->dm_properties_bins_mock[index_bins].GalID_bin_properties.push_back(ig);
#endif // end  _USE_HYBRID_ASSIGNMENT_NEW_
#endif // end _USE_MASS_AS_PRIMARY_OBSERVABLE_

                                           // this->tracer.Halo[ig].GridID=index_bins; // Mark this tarcer as already assigned a property. This will be used below.
                                          this->tracer.Halo[ig].observed=true; // Mark this tarcer as already assigned a property. This will be used below.
#ifndef _USE_STATISTICAL_ASSIGNMENT_
                                          this->dm_properties_bins[index_bins].used_property[i_halo_index] = true; //mark this tracer as already assigned a property. This is used withhin this loop.
#endif
                                           flag = true;
                                           counter_fmf++;
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                                      count_dm++;
#endif
                                         if(number_in_theta_ref[index_bins]>1)
#ifdef _USE_OMP_TEST_
#pragma omp atomic update
#endif
                                            number_in_theta_ref[index_bins]--;
                                        else
                                           number_in_theta_ref[index_bins]=0;
                                       }
                                     }
                                  }
// IF noo multiscale is to be used, and when available options (props) are left, we use the statistical approach P(M|theta)
// This is also needed only in case we are NOT assigning to the reference or not reading from it to assign to mock in general
                            else
                                    {
#endif // enf of _USE_STATISTICAL_ASSIGNMENT_
#ifdef _USE_STATISTICAL_ASSIGNMENT_
                              real_prec aux_h=0;//
                              for(ULONG iy=0;iy< this->params._NPROPbins_bam();++iy) // we only need to go up to the threshold mass
                                  aux_h+=this->ABUNDANCE_normalized[index_2d(iy,index_bins,LENGHTdm)];
                              if(aux_h>0)
                               {
                                 bool flag=false;
                                 while(false==flag)
                                  {
                                   real_prec prob=-10.0;
                                   real_prec ran=10.0;
                                   int i_halo_index=0;
                                   while(prob<ran)
                                    {
#ifdef _USE_OMP_TEST_
                                      i_halo_index= gsl_rng_uniform_int(randa,this->params._NMASSbins()); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#else
                                      i_halo_index= gsl_rng_uniform_int(rn,this->params._NPROPbins_bam()); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#endif
                                      prob = this->ABUNDANCE_normalized[index_2d(i_halo_index,index_bins,LENGHTdm)];
#ifdef _USE_OMP_TEST_
                                      ran = gsl_rng_uniform(randa);
#else
                                      ran = gsl_rng_uniform(rn);
#endif
                                   }
#ifdef _USE_OMP_TEST_
                                   real_prec xr = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                                   real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                                   //OJO ACA QUE NO ESTA CONTEPADO EL CASO EN QEU USEAMOS ASSIGN_MASSES_FROM_REF CUANDO VMAX ES MAIN PROPERTY
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                                   real_prec fraction_p= xr*log10(this->tracer_ref.VMAXBmax[i_halo_index]/this->tracer_ref.VMAXBmin[i_halo_index]);
                                   real_prec lp_halo = log10(this->tracer_ref.VMAXBmin[i_halo_index])+fraction_p ;
                                   this->tracer.Halo[ig].vmax = pow(10,lp_halo);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                                   real_prec fraction_p= xr*log10(this->tracer_ref.MBmax[i_halo_index]/this->tracer_ref.MBmin[i_halo_index]);
                                   real_prec lp_halo = log10(this->tracer_ref.MBmin[i_halo_index])+fraction_p ;
                                   this->tracer.Halo[ig].mass = pow(10,lp_halo);
#endif
                                   this->tracer.Halo[ig].observed=true; // Mark this tarcer as already assigned a property. This will be used below.
                                   counter_fmf++;
                                   flag=true;
                                } //closes while
                             }//closes aux_h>0
#endif // endif for #ifndef _USE_STATISTICAL_ASSIGNMENT_
#ifndef _USE_STATISTICAL_ASSIGNMENT_
                                    }// closes else ion line 9531
#endif
#endif  // end of ifndef MULTISCALE
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_  //this is only done for vmax so far, since we shall  use that property as main property
                     }
#endif
#ifndef _BIN_ACCUMULATE_
                  }  // end of ifs for theta-bins
#endif
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
#if defined(_FULL_VERBOSE_) && !defined (_USE_OMP_TEST_)
                   So.message_screen_flush("Number of dm-like tracers with assigned property = ",static_cast<int>(counter_fmf));
#endif
#else
#ifndef _MULTISCALE_
#if defined _FULL_VERBOSE_ && !defined (_USE_OMP_TEST_)
//                  So.message_screen_flush("Number of tracers with assigned property = ",static_cast<int>(counter_fmf));
#endif
#endif
#endif
              }// end loop over tracers or grid
#ifdef _FULL_VERBOSE_
            std::cout<<endl;
#endif
#ifdef _USE_OMP_TEST_
        }  // end of parallel region
#endif
#if !defined (_SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_) && !defined (_MULTISCALE_)
#ifdef _FULL_VERBOSE_
    this->So.DONE();
    So.message_screen("Number of tracers with assigned property = ",static_cast<int>(counter_fmf));
    std::cout<<endl;
#endif
#endif
    // ***********************
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
    vector<real_prec>left_overs_properties;
    vector<bool>left_overs_properties_used;

    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
      for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
        if(false==this->dm_properties_bins[i].used_property[j])
          {
            left_overs_properties.push_back(this->dm_properties_bins[i].tracer_properties[j]);
            left_overs_properties_used.push_back(false);
          }
#ifdef _FULL_VERBOSE_
       So.message_screen("Assigning randomly to remainig set -dm:");
       std::cout<<endl;
#endif
    for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig)
      if( this->tracer.Halo[ig].identity>0  && false==this->tracer.Halo[ig].observed)// proceed if this object is random
        {
          bool flag=false;
          while(false==flag){
            int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
            if(false==left_overs_properties_used[jj])
              {
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                this->tracer.Halo[ig].observed=true;
                left_overs_properties_used[jj]=true;
                count_dm++;
                flag=true;
              }
           }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property (remaining) = ",static_cast<int>(count_dm));
#endif
         }
#ifdef _FULL_VERBOSE_
    std::cout<<endl;
#endif
    left_overs_properties.clear();
    left_overs_properties.shrink_to_fit();
    left_overs_properties_used.clear();
    left_overs_properties_used.shrink_to_fit();
    // This part is meant for the random distributed tracers

#ifdef _FULL_VERBOSE_
    So.message_screen("Assigning property to random tracers. Expected ", this->tracer.Ntracers_ran);
#endif
    ULONG count_ran=0;
    // This loop will not work perfectly: The container this->dm_properties_for_randoms_bins[theta_bin].tracer_properties[]
    // has been filled with the reference tracers below a threshold Vmax, but this doest not guarranty that
    // if a object here is a random living in a theta bin with Nrans inside, we can authomatically find "Nrans"
    // in a theta bin of the this->dm_properties_for_randoms_bins[theta_bin].tracer_properties[]
    // There are some tracers marked as "not-observed" and "random" living in a theta_bin where the reference has none of these "randoms".
    // Thes objeects are theferefre not assigned a value of vmax. But the container this->dm_properties_for_randoms_bins has more,
    // so we will just assign vmax to the remainnhg orphan objects randomly from the left-overs.
    for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig)
      {
        if( this->tracer.Halo[ig].identity<0)// proceed if this object is random
          {
            count_ran_full++;
            ULONG id=this->tracer.Halo[ig].GridID;
            ULONG index_bins = cell_info_tr[id].Theta_bin;
            ULONG N_tracers_in_theta_bin_r = this->dm_properties_for_randoms_bins[index_bins].tracer_properties.size();
            ULONG N_tracers_left_r =  number_in_theta_ref_for_randoms[index_bins]; // this container is updated below
           // If we have available properties in the theta-bin
        //    std::cout<<ig<<" "<<this->tracer.Halo[ig].identity<<"   "<<this->tracer.Halo[ig].observed<<"   "<<N_tracers_left_r<<endl;
            if(N_tracers_left_r>0) //when we use the ref as a mock, this condition will be always satisfied by construction
              {
                bool flag=false;
                while(flag == false) // este while obliga a elegir una masa -no elegida antes- de entre las disponibles
                 {
                   int i_halo_index_r= gsl_rng_uniform_int(rn,N_tracers_in_theta_bin_r);
                   bool used_tracer = this->dm_properties_for_randoms_bins[index_bins].used_property[i_halo_index_r]; //true or false if the mass was already chosen or not
                   if(false == used_tracer)// if the property of the tracer (e.g., vmax, mass) has not been used before, then assign that property to the current particle i
                    {
                      this->tracer.Halo[ig].vmax= this->dm_properties_for_randoms_bins[index_bins].tracer_properties[i_halo_index_r]; //prop_to_assign;
                      this->tracer.Halo[ig].observed = true; // Mark this tarcer as already assigned a property. This will be used below.
                      this->dm_properties_for_randoms_bins[index_bins].used_property[i_halo_index_r] = true; //mark this tracer as already assigned a property. This is used within this loop.
                      flag = true;
                      counter_fmf++;
                      count_ran++;
                      number_in_theta_ref_for_randoms[index_bins]--;
                    }
                 }
           }
        }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property = ",static_cast<int>(count_ran));
#endif
      }
#ifdef _FULL_VERBOSE_
    std::cout<<endl;
#endif
        // Left overs:
    for(ULONG i = 0; i < this->dm_properties_for_randoms_bins.size(); ++i)
      for(ULONG j = 0; j < this->dm_properties_for_randoms_bins[i].used_property.size(); ++j)
        if(false==this->dm_properties_for_randoms_bins[i].used_property[j])
          {
            left_overs_properties.push_back(this->dm_properties_for_randoms_bins[i].tracer_properties[j]);
            left_overs_properties_used.push_back(false);
          }
    this->dm_properties_for_randoms_bins.clear();
    this->dm_properties_for_randoms_bins.shrink_to_fit();
#ifdef _FULL_VERBOSE_
       So.message_screen("Assigning randomly remainig set (random):");
       std::cout<<endl;
#endif
    for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig)
      if( this->tracer.Halo[ig].identity<0  && false==this->tracer.Halo[ig].observed)// proceed if this object is random
        {
          bool flag=false;
          while(false==flag){
            int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
            if(false==left_overs_properties_used[jj])
              {
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                this->tracer.Halo[ig].observed=true;
                left_overs_properties_used[jj]=true;
                count_ran++;
                flag=true;
              }
           }
#ifdef _FULL_VERBOSE_
       So.message_screen_flush("Number of ran-like tracers with assigned property (remaining) = ",static_cast<int>(count_ran));
#endif
         }
    left_overs_properties.clear();
    left_overs_properties.shrink_to_fit();
    left_overs_properties_used.clear();
    left_overs_properties_used.shrink_to_fit();
#ifdef _FULL_VERBOSE_
   std::cout<<endl;
    So.message_screen("Number of tracers_dm with assigned property = ",static_cast<int>(count_dm));
    So.message_screen("Number of tracers_dm with assigned property_full = ",static_cast<int>(count_dm_full));//to check, ok
    So.message_screen("Number of tracers_ran with assigned property = ",static_cast<int>(count_ran));
    So.message_screen("Number of tracers_ran with assigned property_full = ",static_cast<int>(count_ran_full));//to check, ok
#endif
#endif
#ifndef _MULTISCALE_
      number_in_theta_ref.clear();
      number_in_theta_ref.shrink_to_fit();
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
      number_in_theta_ref_for_randoms.clear();
      number_in_theta_ref_for_randoms.shrink_to_fit();
#endif
#endif
#ifdef _USE_TOTAL_MASS_IN_CELL_
      MOCK_MASS_FIELD.clear();
      MOCK_MASS_FIELD.shrink_to_fit();
#endif
// See documentation for multiscale:
// See documentation for multiscale:
// See documentation for multiscale:
#ifdef _MULTISCALE_
   int N_props_particle_level=0; // Particle level (=0)
   vector<int>number_in_cells_aux( this->params._NGRID(),0); // this is needed regardless the level defined
   ULONG ep_ncounts=this->tracer_ref._NOBJS()-this->tracer._NOBJS();
   if(this->tracer_ref._NOBJS()-this->tracer._NOBJS()<0)
       ep_ncounts*=-1.0;
   // This line is important. The number of tracers in each level is computed by the Catalog class and stored in the params of that clases, which hs been unidirectionally copied from Bam to Class
   // Hence we need to update the params object of the Bam class
//#ifdef _NEW_APPROACH_ASS_
   // Compute the number of tracers in each multiscale level as an average from all references
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
   for(int i=0;i<this->params._Number_of_MultiLevels();++i)
    {
       ULONG ntr=0;
       vector<int>auxa(this->params._Number_of_references(),0);
       for(int j=0;j<this->params._Number_of_references();++j)
//        ntr+=this->tracers_multi[i][j];
            auxa[j]=this->tracers_multi[i][j];
//           ntr/=static_cast<real_prec>(this->params._Number_of_references());
      ntr=get_min<int>(auxa);
      this->params.set_Ntracers_MultiLevels(i,static_cast<ULONG>(floor(ntr)));
    }
#else   // if not, read the numbers from the last load reference.
   for(int i=0;i<this->params._Number_of_MultiLevels();++i)
          this->params.set_Ntracers_MultiLevels(i,this->tracer_ref.params.get_Ntracers_MultiLevels(i));
#endif
    ULONG N_mlevels=0;
     for(int i=0;i<this->params._Number_of_MultiLevels();++i)
        N_mlevels+=this->params.get_Ntracers_MultiLevels(i);
     N_props_particle_level=this->tracer_ref._NOBJS()-ep_ncounts-N_mlevels;
#ifdef _FULL_VERBOSE_
   So.message_screen("Total to be assigned from grid-approach = ", N_mlevels);
   std::cout<<endl;
#endif
     if(this->tracer_ref._NOBJS()<N_mlevels)
      {
       So.message_warning("Error in line", __LINE__);
       throw std::invalid_argument( "ULONG defined variable was assigned a negative value" );
       So.message_warning("This might be caused by pre-proc declarations defining 'mass_cuts' or 'mass_bins'");
       So.message_screen("N_props=", N_props_particle_level);
       So.message_screen("Assigned =",N_mlevels);
       So.message_warning("BMT stops here to avoid error propagation");
       exit(0);
      }
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of reference masses to be assigned at particle level (level 0) = ",N_props_particle_level);
   std::cout<<endl;
   So.message_screen("Identifying galaxy-id to get X-values in the different levels");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG ig=0;ig < this->tracer._NOBJS(); ++ig)
       this->tracer.Halo[ig].observed=false;  // this will be used below, so let keep it inzialized like that
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG ig=0;ig < this->tracer._NOBJS(); ++ig)
      this->tracer.Halo[ig].multi_scale_level=-1;  // this will be used below, so let keep it inzialized like that
   // Since this has a push_bvack inside, I leave it without parallelization
   for(ULONG ig=0;ig < this->tracer._NOBJS(); ++ig)
     {
       ULONG id=this->tracer.Halo[ig].GridID;
       cell_info_tr[id].gal_index.push_back(ig);  // allocate the galaxy ID "ig" in the corresponding grid cell "id"
    }
    So.DONE();
    // NOTA: esta sección estaba puesta después del loop sobre los halos. Ha sido puesta acá para evitar un seg fault que sale con el uso del array e_cells
    // bajo la orden _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
    this->sort_properties();// SOrt all collected properties from the different references top-bottom on each theta bin
#ifdef _NEW_APPROACH_ASS_
/*
#ifdef _FULL_VERBOSE_
   So.message_screen("Re-ordering");
#endif
   // Her we shall order agan the values of siga in eac theta bin such that we do not take athe highst values because of using several refs.
   for(ULONG it=0; it < LENGHTdm; ++it)  // We now sort all properties in the different theta bins from top-to-bottom
     {
       int ll=this->dm_properties_bins[it].tracer_properties.size();// number of collected property-values in this theta bin
       if(ll>0)
        {
          vector<int>nprop_mul(this->params._Number_of_references(),0);
          for(int ip=0; ip<this->params._Number_of_references(); ++ip)
            for(int im=0; im< ll; ++im)// loop over the total number of tracers in this theta bin
              if(this->dm_properties_bins[it].index_reference[im]==ip+1)
                nprop_mul[ip]++; // get the number of tracers from echs ref in the current theta-bin
          vector<vector<real_prec>>prop_at_mul;
          prop_at_mul.resize(this->params._Number_of_references());
          ULONG icounter=0;
          for(int ip=0; ip<this->params._Number_of_references(); ++ip)
           for(int imm=0; imm< nprop_mul[ip]; ++imm)// loop over the total number of tracers in each realization at this theta bin
             {
               prop_at_mul[ip].push_back(this->dm_properties_bins[it].tracer_properties[icounter]);  //prop[index_reference][sigma]
               icounter++;
             }
          int imax=get_min<int>(nprop_mul);
          int ij=0;
          int im=0;
          int ip=0;
          while(im<imax)
           {
            if(nprop_mul[ip]>0)
            {
              this->dm_properties_bins[it].tracer_properties[im]=prop_at_mul[ip][ij];
              im++;
              ip++; // jump to the next reference
              if(ip==this->params._Number_of_references()-1)
                ip=0;// reset to start again from the first reference
              ij++; // junp to the next property
             }
          }
       }
   }
   So.DONE();
*/
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG ig=0;ig < this->tracer._NOBJS(); ++ig)// Loop over the mock tracers to which properties are to be assigned
     {
       ULONG id=this->tracer.Halo[ig].GridID;
       this->tracer.Halo[ig].galThetaID=cell_info_tr[id].Theta_bin;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
       number_in_cells_aux[id]++;    // this is the number (of mocks) counts in the ref grid and will be updated in every level
     }
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_aa
   // At this point, the container dm_porperies has -for each bin of theta, the values of vmax top-bottom sorted. In this part we randomize those values in for each theta bin *but only*  within each lavel.
   // This can hepl becuase we are now collecting properties form several references, and, since we are sorting and taking by order the vmax , for each theta bin we take alway the  "heavier"
   // ones and not use the lighter ones in that bin. BUT, somehow this spoils the n(Vmax) function. It is because of the shape of the abundance; if we randomize, the order will most likely select the low mass objects which are more abundant.
#ifdef _FULL_VERBOSE_
    this->So.message_screen("Randomizing within each level");
#endif
   for(ULONG it=0; it < LENGHTdm; ++it)
     {
       ULONG ll=this->dm_properties_bins[it].tracer_properties.size();
       if(ll>0)
        {
         ULONG llp=0;
         //First randomize for the first level
         vector<ULONG>aux_l;
          vector<real_prec>aux_p;
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
          vector<real_prec>aux_p_secondary;
#endif
          for(int ip=0; ip< ll; ++ip )
            if(this->dm_properties_bins[it].tracer_properties[ip]>=this->params.get_PropThreshold_MultiLevels(0))
              {
                 aux_l.push_back(ip);
                 aux_p.push_back(this->dm_properties_bins[it].tracer_properties[ip]);
            }
          if(aux_l.size()>0)
            {
              gsl_ran_shuffle(rn,&aux_l[0], aux_l.size() ,sizeof(ULONG));  // randomize within the level 0
              for(ULONG ipnew=0; ipnew< aux_l.size(); ++ipnew )  // reassign to dm_prop at this theta bin up tp the size of aux_level4
                {
                  real_prec prop_main=aux_p[aux_l[ipnew]];
                  if(prop_main>=this->params.get_PropThreshold_MultiLevels(0))
                  {
                    this->dm_properties_bins[it].tracer_properties[ipnew]=prop_main;
                      }
                }
             }
          llp=aux_l.size();
          aux_l.clear();aux_l.shrink_to_fit();
          aux_p.clear();aux_p.shrink_to_fit();
          aux_p_secondary.clear();aux_p_secondary.shrink_to_fit();
          // now for the other levels:
          if(this->params._Number_of_MultiLevels()>1)
            {
              for(int iml=1;iml<this->params._Number_of_MultiLevels();iml++)// Now randomize for the other levels if there are
              {
                vector<ULONG>aux_l;
                vector<real_prec>aux_p;
                vector<real_prec>aux_p_secondary;
                for(int ip=0; ip< ll; ++ip )
                  if(this->dm_properties_bins[it].tracer_properties[ip]>=this->params.get_PropThreshold_MultiLevels(iml) && this->dm_properties_bins[it].tracer_properties[ip] < this->params.get_PropThreshold_MultiLevels(iml-1))
                    {
                      aux_l.push_back(ip);
                      aux_p.push_back(this->dm_properties_bins[it].tracer_properties[ip]);
                    }
                if(aux_l.size()>0)
                 {
                  gsl_ran_shuffle(rn,&aux_l[0], aux_l.size() ,sizeof(ULONG));  // randomize within the level 4
                  for(int im=llp; im< llp+aux_l.size(); ++im )  // reassign to dm_prop at this theta bin up tp the size of aux_level4
                    {
                      real_prec prop_main=aux_p[aux_l[im-llp]];
                      if(prop_main>=this->params.get_PropThreshold_MultiLevels(iml) && prop_main < this->params.get_PropThreshold_MultiLevels(iml-1))
                      {
                         this->dm_properties_bins[it].tracer_properties[im]=prop_main;
                      }
                     }
                  }
                llp+=aux_l.size();
               }
             }
          }
       }
    So.DONE();
#endif  //endif _USE _MULTISCALE_PROPERTY_ASSIGNMENT_NEW_aa
//---------------------------------------------------------------------------------------------------------------------------------------
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
 ULONG cumulative_counter=0;
 vector<ULONG>count_aux_cell(LENGHTdm,0);  // This must be outside of the do-while loop
 cout<<endl;
 for(int iml=0; iml<this->params._Number_of_MultiLevels(); ++iml)   // Loop over the multiscale levels >1, split to avoid if's inside loops
  {
    real_prec min_prop_level=0;
    real_prec max_prop_level=0;
    if(iml==0)
     {
       min_prop_level=this->params.get_PropThreshold_MultiLevels(0);
       max_prop_level=LARGE_NUMBER;
     }
    else
     {
       min_prop_level= this->params.get_PropThreshold_MultiLevels(iml);
       max_prop_level= this->params.get_PropThreshold_MultiLevels(iml-1);
     }
    ULONG Ntracers_ml=this->params.get_Ntracers_MultiLevels(iml);
    ULONG NFT_ml=this->params.get_Nft_MultiLevels(iml);
    ULONG NGRID_MS=NFT_ml*NFT_ml*NFT_ml;// Grid size for the current level
    ULONG counter_av=0;
    ULONG counter_masses_a_tolerance=static_cast<ULONG>(floor(Ntracers_ml*static_cast<real_prec>(this->params.get_Props_Tolerance_MultiLevels(iml))));
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_av)
#endif
    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
        for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
            if((this->dm_properties_bins[i].tracer_properties[j]>=min_prop_level && this->dm_properties_bins[i].tracer_properties[j]<max_prop_level))
              counter_av++;
#ifdef _FULL_VERBOSE_
    So.message_screen("*********************************************************************");
    So.message_screen("Going through LEVEL",iml+1);
    So.message_screen("Minimum value property",min_prop_level);
    So.message_screen("Maximum value property",max_prop_level);
    So.message_screen("Number of requested properties (computed from one reference)=", Ntracers_ml);
    So.message_screen("Requested with tolerance factor=", counter_masses_a_tolerance);
    So.message_screen("Available number of properties (from dm_container) in this level=",counter_av);
    So.message_screen("Suggested grid size for this population=", static_cast<int>(floor(pow(Ntracers_ml,1./3.))));
    So.message_screen("Using Nft=",NFT_ml);
#endif
    vector<s_cell_info_reduced>cell_info_cell(NGRID_MS);
    get_low_res_id_from_high_res_id(this->params._Nft(), NFT_ml,cell_info_cell);
    ULONG counter_multiscale=0;
    do
     {
       ULONG N_cells=NGRID_MS;// cells_id_still_to_assign.size();
       ULONG counter_internal=0;
       ULONG id_ini;
#ifdef _USE_OMP_LX_
#pragma omp parallel for default(shared) private(id_ini) reduction(+:counter_internal)
#endif
       for(id_ini=0; id_ini< N_cells; ++id_ini) // loop over cells with resolution according to current level
        {
          ULONG N_cells_in_cells= cell_info_cell[id_ini].gal_index.size();  // Get the Number of high-res cells (i.e, the original ones) in one low-res cell
          int index_cell_in_low_res_cell= gsl_rng_uniform_int(rn,N_cells_in_cells); //Select randomy one ID of the high-re cells witin the low res cell:
          ULONG id=cell_info_cell[id_ini].gal_index[index_cell_in_low_res_cell];    //this chosen id will be used below to assign mass
          if(number_in_cells_aux[id]>0)// If there are tracers in the spatial cell id with no property assigned (this is reduced by one after every assignment), go
            {
             ULONG index_theta_bins = cell_info_tr[id].Theta_bin;      //Get the theta-bin in which this high-res cell has been classified. This has been calculated in the step PRop<Prop above in this function
             ULONG i_mass_halo_label=count_aux_cell[index_theta_bins]; //Get the -updated in the previous label- the stage at which, in the sorted list of prop, the assignemt is located
             ULONG N_mocks_in_cell=cell_info_tr[id].gal_index.size();  //Number of mock tracers in the cell id
             ULONG jk= gsl_rng_uniform_int(rn,N_mocks_in_cell);        //Choose a random integer in  [0,N_mocks_in_cell).
             ULONG ig=cell_info_tr[id].gal_index[jk];                  //Get the Galaxy ID of the mock tracer chosen from the cell id in the position jk of tracers without mass
             if(false==this->tracer.Halo[ig].observed)                 //If this mock tracer has not yet a property, proceed
               {
                if(this->dm_properties_bins[index_theta_bins].tracer_properties.size()>0)// this if must be added a theta bin filled from the mock might be empty in the ref
                 {
                  real_prec assigned_property = this->dm_properties_bins[index_theta_bins].tracer_properties[i_mass_halo_label];
                  bool used_prop = this->dm_properties_bins[index_theta_bins].used_property[i_mass_halo_label]; //true or false if the mass was already chosen or not
                  if(false==used_prop && (assigned_property>=min_prop_level &&  assigned_property< max_prop_level) )
                   {
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
                     if(true==initial_assignment)  // should we take this out of the loop??
#endif
                       this->tracer.Halo[ig].mass=assigned_property ;
#ifdef _ASSIGN_MASS_POST_
                     else
                       this->tracer.Halo[ig].vmax=assigned_property ;
#endif
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_  // THis applies if mass are to be obtained in as econd assignment campaign
                     if(true==initial_assignment)
                     {
#endif
                       this->tracer.Halo[ig].vmax=assigned_property ;
#ifdef _ASSIGN_MASS_POST_
                      }
                      else
                         this->tracer.Halo[ig].mass=assigned_property ;
#endif
#endif
//                     if(counter_internal< counter_masses_a_tolerance-counter_multiscale)
                       //We could exceed the desired value during the loop when tolerance is <1. But in that case we will end-up asiging the right properties, as if tolerance=1
                       //On the other hand, we will necer assign more than what there is available. Hence I omit any if here.
                       counter_internal++;
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
                      int auxx=0;
//                      if(iml==this->params._Number_of_MultiLevels()-1)
//                      auxx=gsl_rng_uniform_int(rn,3);
                      count_aux_cell[index_theta_bins]+=1+auxx;  // This keeps track of the X-properties used. This is IMPORTANT for it increases 1 as long as a prop is assigned, such that the following assignment will poin to the next less massive objectand will be used below for campaign at proceeding levels. However, if more than one ref is used, thismight bring problems for will select the jo
#else
#ifdef _USE_OMP_LX_
#pragma omp atomic
#endif
                      count_aux_cell[index_theta_bins]++;  // This keeps track of the X-properties used. This is IMPORTANT for it increases 1 as long as a prop is assigned, such that the following assignment will poin to the next less massive objectand will be used below for campaign at proceeding levels
#endif
#ifdef _USE_OMP_LX_
#pragma omp atomic
#endif
                      number_in_cells_aux[id]--; // Subtract one unity from the count of mock tracers in the id cell
                      this->tracer.Halo[ig].gal_cwt=cell_info_tr[id].cwt;
                      this->tracer.Halo[ig].local_dm=cell_info_tr[id].local_dm;
                      this->tracer.Halo[ig].observed=true;
                      this->tracer.Halo[ig].multi_scale_level=iml+1;
                      this->dm_properties_bins[index_theta_bins].used_property[i_mass_halo_label]=true;  //mark this mass as already assigned
                 }// closes if(false==used_prop
             }// closes if(this->dm_properties_bins[index_theta_bins].tracer_properties.size()>0
         } // closes if(false==this->tracer.Halo[ig].observed)
       }//closes if(number_in_cells_aux[id]>0)
   //    if(counter_multiscale == counter_masses_a_tolerance)
   //      break;
      }// end loop over cells

     counter_multiscale+=counter_internal;
#ifdef _FULL_VERBOSE_
     So.message_screen_flush("\tNumber of properties assigned =",static_cast<int>(counter_multiscale));
#endif
    }while(counter_multiscale < counter_masses_a_tolerance);   // closes do: do untill all tracers requested for this lavel have assigned propety
    cumulative_counter+=counter_multiscale;
    ULONG counter_masses_bb=0; // Number of tracers without assigned property below the current minimum of this level
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_masses_bb)
#endif
    for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
        for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
            if(false==this->dm_properties_bins[i].used_property[j] && (this->dm_properties_bins[i].tracer_properties[j]<=min_prop_level))//count what is below the minimum of the current level
              counter_masses_bb++;
#ifdef _FULL_VERBOSE_
    So.message_screen(". Partial number of properties assigned =",cumulative_counter);
    So.message_screen("Pending for assignment =",this->tracer.Halo.size()-cumulative_counter);
    So.message_screen("Number of properties available from the reference(s) container =",counter_masses_bb);
#endif
  } // closes Loop over the multiscale levels
#endif // end for_USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
// \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\INdividual assignment  for X<X threshold *****************************
#ifdef _FULL_VERBOSE_
  So.message_screen("*********************************************************************");
  So.message_screen("Individual assignment");
#endif
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
   So.message_screen("Available = ", N_props_particle_level);
#endif
   ULONG counter_from_ref=0; // Number of tracers read from the reference
   vector<ULONG>id_to_assign;
   ULONG counter_masses_b=0; // Number of tracers without assigned property
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_masses_b)
#endif
   for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
     if(this->tracer.Halo[ig].observed==false)
       counter_masses_b++;
#ifdef _FULL_VERBOSE_
   So.message_screen("Requested(obtained from non-observed in tracer.Halo) =", counter_masses_b);
   So.message_screen("Requested(obtained from N_props_particle_level) =",N_props_particle_level);
   So.message_screen("Requested(check Nobs - cumulative) =", this->tracer._NOBJS()-cumulative_counter);
#endif
   ULONG counter_masses_b_tolerance=static_cast<ULONG>(floor(counter_masses_b*static_cast<real_prec>(TOLERANCE_FACTOR_L0)));
#ifdef _FULL_VERBOSE_
   So.message_screen("Requested with tolerance factor = ", counter_masses_b_tolerance);
#endif
  this->get_minimum_multiscale_property();
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   int bin_threshold=get_bin(log10(this->minimum_multiscale_property), this->params._LOGMASSmin(),this->params._NMASSbins(),(this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NMASSbins()), false);
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
   int bin_threshold=get_bin(log10(this->minimum_multiscale_property), log10(this->params._VMAXmin()),this->params._NMASSbins(),log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(this->params._NMASSbins()), false);
#endif

   // count if there are still available masses from the reference:
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_from_ref)
#endif
   for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
       for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
           if(false==this->dm_properties_bins[i].used_property[j] && this->dm_properties_bins[i].tracer_properties[j]< this->minimum_multiscale_property)
             counter_from_ref++;

   // Since at some level there can be some tolerance, we can have, at this level, more tracers without assignment than tracers in the reference. Hence we ask
   // If all tolerance values are set to unity and we are assigning to the rference, the available number of references at this level is equal to the number of tracers without assignment.
   ULONG ntracers_this_level=min(counter_masses_b_tolerance, counter_from_ref);
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of tracers available from the reference(s) container at individual level =", counter_from_ref);
   So.message_screen("Number of properties to assign ",ntracers_this_level);
#endif
    this->randomize_properties(); // Randomize properties in the contained dm_properties: these have been previously sorted for the assignment
    ULONG counter_masses_0=0; // NUmber of tracers with property assigned: This will be updated inside the loop
#ifdef _MASS_ASSIGNMENT_TO_REFERENCE_
   if(N_props_particle_level>0)
#else
   if(counter_masses_b>0) // If we need to assign (because there are still tracer with observedc=false), proceed
#endif
     {
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
       do
        { // DO this part until all available tracers in the reference set are used. When sampling from the DM of the calibration, this do-loop is not needed as we can use assignment with 1pt abundance if needed
#endif
// ++++++++++++++++++++++++++++++++++++++++++
// This is an ideal space for paralelization: we can parallelize each loop over the tracers, collect the number of assigned props
// add them
// ++++++++++++++++++++++++++++++++++++++++++
          ULONG counter_internal=0;
          ULONG ig;
#ifdef _USE_OMP_L0_
#pragma omp parallel for default(shared) private(ig) reduction(+:counter_internal)
#endif
          for(ig=0; ig <this->tracer.Halo.size(); ++ig)// loop over the MOCK tracers
           {
             if(false==this->tracer.Halo[ig].observed)// if this tracer has not been assigned a property, proceed
                {
                  ULONG id=this->tracer.Halo[ig].GridID; // Retrieve ID of spatial cell where the tracer is located
                  ULONG index_bins = cell_info_tr[id].Theta_bin; // Retrieve the bin in theta where that cells has been identified
                  if(counter_from_ref>counter_masses_0) // counter_masses_0 is updated inside this loop. This "if" statment ensures that we try to use all refernece properties
                   {
                    ULONG N_props_in_bin = this->dm_properties_bins[index_bins].tracer_properties.size();// Get the number of tracers in that theta bin
                    if(N_props_in_bin>0)
                      {
                        int i_mass_halo_label= gsl_rng_uniform_int(rn,N_props_in_bin);// select randomly one mock halo in the corresponding theta bin
                        real_prec assigned_property=this->dm_properties_bins[index_bins].tracer_properties[i_mass_halo_label] ;
                        bool used_prop = this->dm_properties_bins[index_bins].used_property[i_mass_halo_label]; //true or false if the mass was already chosen or not
#ifndef _ASSIGN_TO_CALIBRATION_
                        if(false==used_prop && assigned_property <this->minimum_multiscale_property )   //If used_prop=false (the property has not been used before), then assign that mass to the current particle i
                         {
#endif
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                           this->tracer.Halo[ig].mass=assigned_property ;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                           this->tracer.Halo[ig].vmax=assigned_property ;
#endif
                           this->tracer.Halo[ig].observed=true;
                           this->tracer.Halo[ig].multi_scale_level=0;
                           this->dm_properties_bins[index_bins].used_property[i_mass_halo_label] = true; //mark this mass as already assigned
                           this->tracer.Halo[ig].gal_cwt=cell_info_tr[id].cwt;
                           this->tracer.Halo[ig].local_dm=cell_info_tr[id].local_dm;
                           counter_internal++;
                         }
                        else
                         {
                           this->tracer.Halo[ig].observed=false;
                           this->tracer.Halo[ig].multi_scale_level=0;
                        }
                     } //closes  if(N_props_in_bin>0)
                  } // closes if(counter_from_ref>counter_masses_0)
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                   So.message_screen_flush("Number of masses assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#endif
#endif
               }// closes if (false==observed)
            }// end loop over tracers
           counter_masses_0+=counter_internal;
#ifdef _FULL_VERBOSE_
           So.message_screen_flush("\tNumber of properties assigned  = ", counter_masses_0);
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
          }while(counter_masses_0 < ntracers_this_level); // this will be never exact, for the inner loop cannot comunicate outside the current nu ber of accepted `
#endif
        } // end if counter_masses_b > 0
        So.DONE();
#ifdef _FULL_VERBOSE_
#ifdef _USE_OMP_BARR_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     So.message_screen("Number of Vmax assigned   = ", counter_masses_0);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
     So.message_screen("Number of Masses assigned   = ", counter_masses_0);
#endif
#endif
#endif
     cumulative_counter+=counter_masses_0;
#ifdef _FULL_VERBOSE_
   std::cout<<endl;
   So.message_screen("Partial number of properties assigned = ",cumulative_counter);
   So.message_screen("Number of properties to assign = ",this->tracer._NOBJS()-cumulative_counter);
   std::cout<<endl;
#endif
   ULONG nc_test=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nc_test)
#endif
   for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
      if(false==this->tracer.Halo[ig].observed)
         nc_test++;
#ifdef _FULL_VERBOSE_
   So.message_screen("Number of 'un-observed' (not assigned) tracers = ", nc_test);
   So.message_screen("Number of properties to assign from unobserved tracers= ",this->tracer._NOBJS()-nc_test);
#endif
    // IF THERE RE STILL PARTICLES TO BE ASSIGNED PROPERTIES:
   if(this->tracer._NOBJS()-cumulative_counter>0)
    {
     if(true==initial_assignment)//Only valid for initial assignment
      {
        counter_from_ref=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:counter_from_ref)
#endif
        for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
         for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
          if(false==this->dm_properties_bins[i].used_property[j])
           counter_from_ref++;
//---------------------------------------------------------------------------------------------
// When sampling the same DM field used for the calibration, at this point there should
// be no available tracers in the reference from which properties can be read.
// If we are assigning to the DM used for calibration, things should end up here.
// However, at this point we still have tracers without propertes and lefts-overs from
// the reference container
// This is true also if we are assigning to the same ref if the tolerance for L0 is <1
// (which is better, as a full assignment is very slow). We can then assign the missing
// ones as randomly from the lefts overs of the refereece. In oder to do that, we collect
// in a single arraay the lefst overs. Also, here we do not set a limit as a function of
// the minumum multiscale property. This is of no impact for the reconstruction of the n(M)
// function though
//---------------------------------------------------------------------------------------------
#ifdef _FULL_VERBOSE_
        So.message_screen("CHECK: available number of properties from the reference container = ", counter_from_ref);
        So.message_screen("Assigning remaining properties randomly selecting from the previous number:");
#endif
        vector<real_prec>left_overs_properties;
        vector<bool>left_overs_properties_used;
        for(ULONG i = 0; i < this->dm_properties_bins.size(); ++i)
           for(ULONG j = 0; j < this->dm_properties_bins[i].used_property.size(); ++j)
              if(false==this->dm_properties_bins[i].used_property[j])
                {
                  left_overs_properties.push_back(this->dm_properties_bins[i].tracer_properties[j]  );
                  left_overs_properties_used.push_back(false);
                 }

#ifdef _FULL_VERBOSE_
        So.message_screen("CHECK: size of container for these objects = ", left_overs_properties.size());
        So.message_screen("Number of tracers to assign with tolerance L0b = ", static_cast<ULONG>(nc_test*TOLERANCE_FACTOR_L0b));
#endif

        if(nc_test*TOLERANCE_FACTOR_L0b >left_overs_properties.size())// In this case, we have more requested than available
         {
#ifdef _FULL_VERBOSE_
            So.message_screen("There are not enough properties to assign");
            So.message_screen("Changing the number to be assigned to the number of available refs");
#endif
           nc_test=left_overs_properties.size();
         }
        ULONG counter_masses_left=0;
        do{
         for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
           {
            if(false==this->tracer.Halo[ig].observed)
              {
               int jj = gsl_rng_uniform_int(rn,left_overs_properties.size());
               if(false==left_overs_properties_used[jj])
                 {
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                   this->tracer.Halo[ig].mass=left_overs_properties[jj] ;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                   this->tracer.Halo[ig].vmax= left_overs_properties[jj] ;
#endif
                   this->tracer.Halo[ig].observed=true;
                   this->tracer.Halo[ig].multi_scale_level=8;
                   left_overs_properties_used[jj]=true;
                   counter_masses_left++;
                   counter_masses_0++;
                 }
              }
           }
#ifdef _FULL_VERBOSE_
          So.message_screen_flush("\tAssigned from left-overs= ", counter_masses_left);
#endif
          if(counter_masses_left==static_cast<ULONG>(TOLERANCE_FACTOR_L0b*nc_test))
               break;
        }while(counter_masses_left < static_cast<ULONG>(TOLERANCE_FACTOR_L0b*nc_test));
#ifdef _FULL_VERBOSE_
      So.DONE();
      So.message_screen("Total Assigned at Level 0 = ", counter_masses_0);
#endif
      cumulative_counter+=counter_masses_left;
      left_overs_properties_used.shrink_to_fit();
      So.DONE();
     }// close if(true==initial_assignment)
  }//  close if(this->tracer._NOBJS()-cumulative_counter<0)
  N_props_particle_level=counter_masses_0;  // update this number. Asked below
   // END Level 0
#endif // end for _MULTISCALE_
#ifndef _USE_STATISTICAL_ASSIGNMENT_
   this->dm_properties_bins.clear();
   this->dm_properties_bins.shrink_to_fit();
   So.DONE();
#endif
#if defined (_USE_VMAX_AS_PRIMARY_OBSERVABLE_) || defined (_USE_MASS_AS_PRIMARY_OBSERVABLE_)
    }  // close if(true==initial_assignment):
  else  // else for if(true==initial_assignment):
#endif
   if (false==initial_assignment)  //to assign MASSES (or Vmax) based on the information of  Vmax (or Mass) already assigned
    {
#ifdef _FULL_VERBOSE_
      if(h_property==_MASS_)
        So.message_screen("Assigning Mvir.  Expected number of properties to assign = ", this->tracer._NOBJS());
      if(h_property==_VMAX_)
        So.message_screen("Assigning Vmax.  Expected number of properties to assign = ", this->tracer._NOBJS());
      else if(h_property==_RS_)
        So.message_screen("Assigning RS.  Expected number of properties to assign = ", this->tracer._NOBJS());
      else if(h_property==_CONCENTRATION_)
       So.message_screen("Assigning Cvir.  Expected number of properties to assign = ", this->tracer._NOBJS());
      else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_ )
        So.message_screen("Assigning SPIN.  Expected number of properties to assign = ", this->tracer._NOBJS());
#endif
      ULONG N_props_particle_level=this->tracer._NOBJS();
      ULONG Ntot=1;
      // For X1 (see table in documentation):
      Ntot*=this->params._NPROPbins_bam();
      // For X2: 
      if(h_property==_MASS_ || h_property==_VMAX_ || h_property==_CONCENTRATION_ || h_property==_RS_)  // ask because one can use MASS as primeary property
        Ntot*=this->params._NX(); // This has been already updaterd in get_scaling_relations_secondary_propert
      else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
        Ntot*=this->params._NPROPbins_bam();
      // For X3:
      ULONG n_cwt=1;
#ifdef _USE_CWC_mass_
      n_cwt= static_cast<ULONG>(this->params._n_cwt());
#endif
      Ntot*=n_cwt;
      // For X4:
#ifdef _USE_MACH_NUMBER_
      Ntot*=N_BINS_MACH;
#endif
      // For X5:
#ifdef _USE_LOCAL_OVERDENSITY_
      Ntot*=N_BINS_LO;
#endif
#ifdef _USE_TIDAL_ANISOTROPY_SEC_PROP_
  Ntot*=N_BINS_TA;
#endif
//*****************************************
      ULONG Ntot_partial=Ntot;// This must be after the last property advocated
//*****************************************

      Ntot*=this->params._NPROPbins_bam(); // Assign the space for the variable Y
// The following loops are meant to inizialize the tracer containers
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifndef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_ //If this is on, mases have been laready assigned
      if(h_property==_MASS_)
      {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
        {
         this->tracer.Halo[ig].observed=false;
         this->tracer.Halo[ig].mass=0;
        }
    }
#endif
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    if(h_property==_VMAX_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].vmax=0;
        }
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
  if(h_property==_RS_ )
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].rs=0;
        }
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
  if(h_property == _CONCENTRATION_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].concentration=0;
        }
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_ )
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
        {
          this->tracer.Halo[ig].observed=false;
          this->tracer.Halo[ig].spin=0;
          this->tracer.Halo[ig].spin_bullock=0;
        }
#endif
 //  The following section aims at using P(M|Vmax, theta) to assign M (or P(V|M.theta) to assogn V)
  vector<real_prec>aux_cond(Ntot_partial,0);
  ULONG  NBins=this->params._NPROPbins_bam();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Ntot_partial;++i)
    for(ULONG j=0;j<NBins;++j)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
         aux_cond[i]+=static_cast<double>(this->ABUNDANCE_normalized[index_2d(j,i,Ntot_partial)]);
    counter_masses_0=0;
    ULONG Ng_k=0;
    ULONG Ng_s=0;
    ULONG Ng_f=0;
    ULONG Ng_v=0;
    ULONG counter_masses_or=0;
    ULONG counter_masses_global=0;
    ULONG counter_masses_test=0;
#ifdef _USE_OMP_TEST2_
  vector<ULONG>vseeds(NTHREADS,0);
     for(int  i=0;i<vseeds.size();++i)
         vseeds[i]=565+static_cast<ULONG>(i)*565;
     int jthread=0;
     const gsl_rng_type *Trn;
     gsl_rng *randa;
#pragma omp parallel private(jthread, randa, Trn)
    {
      jthread=omp_get_thread_num();
      gsl_rng_default_seed=vseeds[jthread];
      Trn = gsl_rng_default;
      randa = gsl_rng_alloc (Trn);
#pragma omp parallel for reduction(+:counter_masses_0, counter_masses_or, counter_masses_global,Ng_k,Ng_f,Ng_s,Ng_v, counter_masses_test)
#endif
      for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
       {
         ULONG id=this->tracer.Halo[ig].GridID;
         ULONG I_X1=0;
#if defined  _USE_RS_AS_DERIVED_OBSERVABLE_ || defined _USE_SPIN_AS_DERIVED_OBSERVABLE_  || defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
      // IF we neeed to assign RS, CVIR, Spin, put Mass in the X1 slot
         if(h_property==_RS_ || h_property == _CONCENTRATION_ || h_property==_SPIN_  ||  h_property==_SPIN_BULLOCK_ ) // Using DELTA only allowed so far when assigning mass. FOr Rs and Spin we use so far P(Rs|M,Vmax) or P(S|M,Vmax)
            I_X1 =  get_bin(log10(this->tracer.Halo[ig].mass),this->params._LOGMASSmin(),this->params._NPROPbins_bam(), (this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NPROPbins_bam())  ,this->bin_accumulate_borders);
        // Use this if we want to use Vmax (even if primary) to assign seconday, passing over Mvir assigned
//        I_X1 =  get_bin(log10(this->tracer.Halo[ig].vmax),log10(this->params._VMAXmin()),this->params._NPROPbins_bam(), (log10(this->params._VMAXmax())-log10(this->params._VMAXmin()))/static_cast<double>(this->params._NPROPbins_bam())  ,this->bin_accumulate_borders);
#endif
        else
         { // Unless we want to assign M(or Vmax) in which case we put Vmax (or M)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_  // Here this means that Vmax is already available and we can learn from it
          if(h_property==_MASS_)
            I_X1=get_bin(log10(this->tracer.Halo[ig].vmax), log10(this->params._VMAXmin()),this->params._NPROPbins_bam(),log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<double>(this->params._NPROPbins_bam()), true);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ // If mvir is primary  we can now assign vmax from P(Vmax|M...)
            I_X1=get_bin(log10(this->tracer.Halo[ig].mass), this->params._LOGMASSmin(),this->params._NPROPbins_bam(),(this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<double>(this->params._NPROPbins_bam()), true);
#endif
         }
        ULONG I_X2=0;
        ULONG NbinsX2=1;
        if(h_property==_MASS_ || h_property==_VMAX_ || h_property==_CONCENTRATION_ || h_property==_RS_)
           {
             real_prec xdm = static_cast<real_prec>(this->delta_X[id]); // this has been read and kept from get_new_dm(): there the CWC has been performed allcoated in cwclass
             I_X2  = get_bin(xdm,this->s_mins.prop1,this->params._NX(),this->s_deltas.prop1,this->bin_accumulate_borders);
             NbinsX2=this->params._NX();
          }
#if defined  _USE_RS_AS_DERIVED_OBSERVABLE_ || defined _USE_SPIN_AS_DERIVED_OBSERVABLE_ ||  defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
         else if( h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
           {
#if defined  _USE_RS_AS_DERIVED_OBSERVABLE_  // Recal that for assignment, Rs and Cvir cannot coexist
            I_X2  = get_bin(log10(this->tracer.Halo[ig].rs),log10(this->params._RSmin()),this->params._NPROPbins_bam(),this->tracer.logdeltaRS,this->bin_accumulate_borders);
#elif defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
            I_X2  = get_bin(log10(this->tracer.Halo[ig].concentration),log10(this->params._CONCENTRATIONmin()),this->params._NPROPbins_bam(),log10(this->params._CONCENTRATIONmax()/this->params._CONCENTRATIONmin())/static_cast<double>(this->params._NPROPbins_bam()),this->bin_accumulate_borders);
#endif
            NbinsX2=this->params._NPROPbins_bam();
          }
#endif
         ULONG I_X3=0;
#ifdef _USE_CWC_mass_
          I_X3=static_cast<ULONG>(this->cwclass.get_Tclassification(id));
          if(true==initial_assignment)
            {
              if(I_X3 == I_KNOT)
                  Ng_k++;
              if(I_X3 == I_FILAMENT)
                  Ng_f++;
              if(I_X3 == I_SHEET)
                  Ng_s++;
              if(I_X3 == I_VOID)
                  Ng_v++;
            }
#endif
          ULONG I_X4=0;
#ifdef _USE_MACH_NUMBER_
          I_X4= get_bin(this->tracer.Halo[ig].mach_number,this->tracer_ref._min_mach(), N_BINS_MACH, (this->tracer_ref._max_mach()-this->tracer._min_mach())/static_cast<real_prec>(N_BINS_MACH), this->bin_accumulate_borders);
#endif
          ULONG I_X5=0;
#ifdef _USE_LOCAL_OVERDENSITY_
          I_X5= get_bin(this->tracer.Halo[ig].local_overdensity, this->tracer_ref._min_local_overdensity(), N_BINS_LO,  (this->tracer_ref._max_local_overdensity()-this->tracer._min_local_overdensity())/static_cast<real_prec>(N_BINS_LO) ,this->bin_accumulate_borders);
#endif
         ULONG index_dm=0; 
         ULONG I_X6=0;
#ifdef _USE_TIDAL_ANISOTROPY_SEC_PROP_ // if defined, onlyu applies to cvir, spin,. not to mass
         I_X6= get_bin(this->tracer.Halo[ig].tidal_anisotropy, this->tracer_ref._min_tidal_anisotropy(), N_BINS_TA,  (this->tracer_ref._max_tidal_anisotropy()-this->tracer_ref._min_tidal_anisotropy())/static_cast<real_prec>(N_BINS_TA) ,this->bin_accumulate_borders);
         index_dm=index_6d(I_X1,I_X2,I_X3,I_X4,I_X5,I_X6, NbinsX2,n_cwt,N_BINS_MACH,N_BINS_LO,N_BINS_TA);
#else
          index_dm=index_5d(I_X1,I_X2,I_X3,I_X4,I_X5,NbinsX2,n_cwt,N_BINS_MACH,N_BINS_LO);
#endif
        real_prec aux_h=aux_cond[index_dm];
         if(aux_h>0)//if all bins in abundance_normalized are zero, means no availability of props.
              {
                counter_masses_test++;
                real_prec prob=-10.0;
                real_prec ran=10.0;
                int i_halo_index=0;
                while(prob<ran)
                 {
#ifdef _USE_OMP_TEST2_
                      i_halo_index= gsl_rng_uniform_int(randa,NBins); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#else
                      i_halo_index= gsl_rng_uniform_int(rn,NBins); // we have chosen to use NMASSBIns also for RS, Vmax etc when these are to be assigned
#endif
                      ULONG index_or=index_2d(i_halo_index,index_dm,Ntot_partial);
                      prob = this->ABUNDANCE_normalized[index_or];
#ifdef _USE_OMP_TEST2_
                      ran = gsl_rng_uniform(randa);
#else
                      ran = gsl_rng_uniform(rn);
#endif
                  }
#ifdef _USE_OMP_TEST2_
                  real_prec xr = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                  real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                  if(_VMAX_ == h_property)
                    {
                      real_prec fraction_mass= xr*log10(this->tracer_ref.VMAXBmax[i_halo_index]/this->tracer_ref.VMAXBmin[i_halo_index]);
                      real_prec lmass_halo = log10(this->tracer_ref.VMAXBmin[i_halo_index])+fraction_mass ;
                      this->tracer.Halo[ig].vmax = pow(10,lmass_halo);
                  }
#elif defined (_USE_VMAX_AS_PRIMARY_OBSERVABLE_)
                  if(_MASS_ == h_property)
                    {
                      real_prec fraction_mass= xr*log10(this->tracer_ref.MBmax[i_halo_index]/this->tracer_ref.MBmin[i_halo_index]);
                      real_prec lprop_halo = log10(this->tracer_ref.MBmin[i_halo_index])+fraction_mass ;
                      this->tracer.Halo[ig].mass = pow(10,lprop_halo);
                    }
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
                   else if(_RS_== h_property)
                    {
                      real_prec fraction_mass= xr*log10(this->tracer_ref.RSBmax[i_halo_index]/this->tracer_ref.RSBmin[i_halo_index]);
                      real_prec lprop_halo = log10(this->tracer_ref.RSBmin[i_halo_index])+fraction_mass ;
                      this->tracer.Halo[ig].rs = pow(10,lprop_halo);
                   }
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
                   else if(_CONCENTRATION_== h_property)
                    {
                      real_prec fraction_mass= xr*log10(this->tracer_ref.CONCENTRATIONBmax[i_halo_index]/this->tracer_ref.CONCENTRATIONBmin[i_halo_index]);
                      real_prec lprop_halo = log10(this->tracer_ref.CONCENTRATIONBmin[i_halo_index])+fraction_mass ;
                      this->tracer.Halo[ig].concentration = pow(10,lprop_halo);
                    }
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
                    else if(_SPIN_== h_property)
                      {
                        real_prec fraction_mass= xr*log10(this->tracer_ref.SPINBmax[i_halo_index]/this->tracer_ref.SPINBmin[i_halo_index]);
                        real_prec lprop_halo = log10(this->tracer_ref.SPINBmin[i_halo_index])+fraction_mass ;
                        this->tracer.Halo[ig].spin = pow(10,lprop_halo);
                      }
                    else if(h_property==_SPIN_BULLOCK_)
                      {
                        real_prec fraction_mass= xr*log10(this->tracer_ref.SPINBmax[i_halo_index]/this->tracer_ref.SPINBmin[i_halo_index]);
                        real_prec lprop_halo = log10(this->tracer_ref.SPINBmin[i_halo_index])+fraction_mass ;
                        this->tracer.Halo[ig].spin_bullock = pow(10,lprop_halo);
                      }
#endif
                     this->tracer.Halo[ig].observed=true;
                     counter_masses_0++;
                     counter_masses_or++;
                }//closes if (aux_h>0)
#if !defined  _ASSIGN_PROPERTIES_TO_REFERENCE_ || defined _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_ || defined _ASSIGN_PROPERTIES_TO_MOCK_
                  else
                  {
                      real_prec prob=-10.0;
                      real_prec ran=10.0;
                      real_prec mass_tracer=0;
                    // here we assign mass according to the measured abundance
                      while(prob<ran)
                       {
#ifdef _USE_OMP_TEST2_
                        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                        real_prec xr = static_cast<real_prec>(gsl_rng_xuniform(rn));
#endif
                        real_prec fraction_mass=0;
                        fraction_mass = xr*log10(this->prop_max/this->prop_min);
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_ // if we arwe assigning mass here, the minium is the nominal
                        mass_tracer=pow(10,log10(this->prop_min)+fraction_mass);
#elif  defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ // if we are assigning vmax here, the minium depends on the mass
                        mass_tracer=pow(10,log10(limit_vmax_halos(this->tracer.Halo[ig].mass))+fraction_mass);
#endif
                        real_prec aux_mass= (mass_tracer <= this->prop_min? this->prop_min: mass_tracer);
                        aux_mass= (mass_tracer >= this->prop_max ? this->prop_max : mass_tracer);
                        if(aux_mass<=this->prop_min || aux_mass >= this->prop_max )
                          prob=0.0;
                        else
                          {
#ifdef _USE_OMP_TEST2_
                            ran  = static_cast<real_prec>(gsl_rng_uniform(randa));
#else
                            ran  = static_cast<real_prec>(gsl_rng_uniform(rn));
#endif
                            prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_mass, acc));
                          }//closes else
                        mass_tracer=aux_mass;
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_ // if this applies, we now assign masses
                        if(h_property==_MASS_)
                            this->tracer.Halo[ig].mass=mass_tracer;
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ // if this applies, we now assign vmax
                        if(h_property==_VMAX_)
                            this->tracer.Halo[ig].vmax=mass_tracer;
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
                        else if(h_property==_RS_)
                            this->tracer.Halo[ig].rs=mass_tracer;
#endif
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
                        else if(h_property==_CONCENTRATION_)
                          this->tracer.Halo[ig].concentration=mass_tracer;
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
                        else if(h_property==_SPIN_)
                          this->tracer.Halo[ig].spin=mass_tracer;
                        else if(h_property==_SPIN_BULLOCK_)
                          this->tracer.Halo[ig].spin_bullock=mass_tracer;
#endif
                      } // closes while(prop<ran)
                    this->tracer.Halo[ig].observed=true;
                    counter_masses_global++;
                    counter_masses_0++;
                   }// closs else
#endif   //end of ifndef  _ASSIGN_PROPERTIES_TO_REFERENCE_
#ifndef test_vmax_mass
            }
#endif
#ifdef _FULL_VERBOSE_
#ifndef _USE_OMP_TEST2_
          //andres commented when no parallel is on
          So.message_screen_flush("Number of properties assigned   = ",static_cast<int>(counter_masses_0));
#endif
#endif
     } // end loop over tracers
#ifdef _USE_OMP_TEST2_
     gsl_rng_free(randa);
    } // closes parallel region
#endif
        So.message_screen("Number of properties assigned (test) = ",static_cast<ULONG>(counter_masses_test));


#ifdef _FULL_VERBOSE_
    std::cout<<endl;
#ifdef _USE_OMP_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
          So.message_screen("Number of vmax's assigned   = ",static_cast<int>(counter_masses_0));
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        So.message_screen("Number of properties assigned = ",static_cast<ULONG>(counter_masses_0));
#endif
#endif
    if(h_property==_MASS_)
      So.message_screen("Number of Mvir values assigned (using joint prob-dist) = ",static_cast<ULONG>(counter_masses_or));
    else if(h_property==_RS_)
       So.message_screen("Number of Rs values assigned (using joint prob-dist) = ",static_cast<ULONG>(counter_masses_or));
    else if(h_property==_CONCENTRATION_)
       So.message_screen("Number of Cvir values assigned (using joint prob-dist) = ",static_cast<ULONG>(counter_masses_or));
    else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
      So.message_screen("Number of Spin values assigned (using joint prob-dist) = ",static_cast<ULONG>(counter_masses_or));
    So.message_screen("Number of properties assigned with global abundance = ",static_cast<ULONG>(counter_masses_global));
#ifdef _USE_CWC_mass_
    if(true==initial_assignment)
    {
      So.message_screen("Fraction of tracers in knots:",  Ng_k*100./static_cast<real_prec>(this->tracer._NOBJS()), "%");
      So.message_screen("Fraction of tracers in filaments:", Ng_f*100./static_cast<real_prec>(this->tracer._NOBJS()), "%");
      So.message_screen("Fraction of tracers in sheets:",  Ng_s*100./static_cast<real_prec>(this->tracer._NOBJS()), "%");
      So.message_screen("Fraction of tracers in voids:", Ng_v*100./static_cast<real_prec>(this->tracer._NOBJS()), "%");
    }
#endif
    So.DONE();
#endif
    N_props_particle_level=counter_masses_0;
  } // closes  else if (false==initial_assignment)
#ifdef _FULL_VERBOSE_
  So.message_screen("Freeing memmory in line",__LINE__);
#endif
  this->ABUNDANCE_normalized.clear(); this->ABUNDANCE_normalized.shrink_to_fit();
#ifndef test_vmax_mass
  this->dm_properties_bins.clear();
  this->dm_properties_bins.shrink_to_fit();
#endif
#ifdef _USE_TRACERS_IN_CELLS_
  MOCK_DEN_FIELD.clear();
  MOCK_DEN_FIELD.shrink_to_fit();
#endif
  So.DONE();
  ULONG assigned_properties_mf=0;
#ifndef _MULTISCALE_   // If no multiscale, assigment was already done for some.
  assigned_properties_mf=counter_fmf;
#else
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:assigned_properties_mf)// ask, after all the multilevel or other assignmetn, how many are there left without property
#endif
   for(ULONG ig=0; ig < this->tracer._NOBJS();++ig)
      if(false==this->tracer.Halo[ig].observed)
         assigned_properties_mf++;
#endif
  ULONG Ncount_ab=0;
  if(true==initial_assignment)  //this applies for vmax assignment when the global vmax function  is needed
    {
      if(assigned_properties_mf != this->tracer._NOBJS())
        {
#ifdef _FULL_VERBOSE_
          So.message_screen("Number of tracers with *no property* assigned so far = ", this->tracer._NOBJS()-assigned_properties_mf);
#endif
          ULONG nc_test=0;
#pragma parallel for reduction(+:nc_test)
          for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig)
            if(false==this->tracer.Halo[ig].observed)
              nc_test++;
#ifdef _FULL_VERBOSE_
          So.message_screen("Check: Number of tracers with *no property* assigned so far, from 'non-observed tracer' = ", nc_test);
          std::cout<<endl;
          So.message_screen("Assigning properties using global *property* function");
#endif
          NTHREADS=_NTHREADS_;
          omp_set_num_threads(NTHREADS);
#ifdef _USE_OMP_TEST_
          int jthread=0;
          vector<ULONG>vseeds(NTHREADS,0);
          for(int i=0;i<vseeds.size();++i)
             vseeds[i]=35+static_cast<ULONG>(i+14)*56045;
#endif
          const gsl_rng_type *Trn;
          gsl_rng *rna;
#ifdef _USE_OMP_TEST_
#pragma omp parallel private (jthread, Trn, rna)// this is causing problem
        {
              jthread=omp_get_thread_num();
              gsl_rng_default_seed=vseeds[jthread];
#else
             gsl_rng_default_seed=1555;
#endif
            Trn = gsl_rng_mt19937;//_default;
            rna= gsl_rng_alloc (Trn);
#ifdef _USE_OMP_TEST_
#pragma omp for reduction(+:Ncount_ab) // this parallelization was causing a problem
#endif
            for(ULONG ig=0; ig < this->tracer._NOBJS(); ++ig)
              {
                if(false==this->tracer.Halo[ig].observed || this->tracer.Halo[ig].vmax<=0)
                  {
                     real_prec prob=-10.0;
                     real_prec ran=10.0;
                     real_prec property_tracer=0;
                     while(prob<ran)
                        {
                           real_prec xr = static_cast<real_prec>(gsl_rng_uniform(rna));
                           real_prec fraction_prop = xr*log10(this->prop_max/this->prop_min);
                           property_tracer=pow(10,log10(this->prop_min)+fraction_prop);
                           real_prec aux_prop= (property_tracer < this->prop_min? this->prop_min : property_tracer);
                           aux_prop= (property_tracer > this->prop_max? this->prop_max: property_tracer);
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                           if(aux_prop<this->tracer_ref.VMAXBin[0] || aux_prop>this->tracer_ref.VMAXBin[nbins_mf-1])
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                           if(aux_prop<=this->tracer_ref.MBmin[0] || aux_prop>=this->tracer_ref.MBmin[nbins_mf-1])
#endif
                             prob=0.0;
                           else
                             {
                                ran  = static_cast<real_prec>(gsl_rng_uniform(rna));
                                prob = static_cast<real_prec>(gsl_spline_eval (spline, aux_prop, acc));
                             }
                            property_tracer=aux_prop;
                          }
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
                         this->tracer.Halo[ig].mass=property_tracer;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                            this->tracer.Halo[ig].vmax=property_tracer;
#endif
                          this->tracer.Halo[ig].observed=true;
                          this->tracer.Halo[ig].multi_scale_level=-1;
                          Ncount_ab++;
                      }//closes if
               }// closes loop over tracers
          gsl_rng_free(rna);
#ifdef _USE_OMP_TEST_
             }  // closes parallel region
#endif
          So.DONE();
#ifdef _FULL_VERBOSE_
      So.message_screen("Number of props assigned using global abundance function =",Ncount_ab);
      std::cout<<endl;
#endif
            }
    }
  gsl_rng_free (rn);
  cumulative_counter+=Ncount_ab;
  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  string fname_mass_function_Y_ref = this->params._Output_directory()+"tracer_ref_abundance.txt";
  string fname_mass_function_Y = this->params._Output_directory()+"tracer_mock_abundance_R"+to_string(this->params._realization())+".txt";
  this->tracer.aux_flag=false;
  if(true==initial_assignment)
    {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      this->tracer.params.set_i_vmax_g(POSITIVE_INT); //this allows the tracer to get the vmax function
      this->tracer.params.set_i_mass_g(NEGATIVE_INT); //no mass information here
      this->tracer.params.set_i_spin_bullock_g(NEGATIVE_INT); //no mass information here
      this->tracer.params.set_i_rs_g(NEGATIVE_INT); //no mass information here
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
      this->tracer.params.set_i_vmax_g(-5); //no mass information here
      this->tracer.params.set_i_mass_g(5); //this allows the tracer to get the mass function
#endif
      this->tracer.get_property_function(fname_mass_function_Y);
      this->tracer_ref.get_property_function(fname_mass_function_Y_ref);
#ifdef _USE_GNUPLOT_ABUNDANCE_V_PLOT_
      this->gp_abundance_v<<"set log x \n";
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      this->gp_abundance_v<<"set xlabel 'Vrms [km/s]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      this->gp_abundance_v<<"set ylabel 'log n(Vrms)' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
#elif defined  _USE_MASS_AS_PRIMARY_OBSERVABLE_
      this->gp_abundance_v<<"set xlabel 'Mass [Ms / h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      this->gp_abundance_v<<"set ylabel 'log n(M) h /Ms (h / Mpc)³]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
#endif
      vector<pair<real_prec, real_prec> > xy_pts_v;
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      for(ULONG i=0; i<this->tracer.vmax_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.VMAXBmin[i]+(i+0.5)*(this->tracer.VMAXBmax[i]-this->tracer.VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.vmax_function[i])));
      this->gp_abundance_v<<"plot [120:1000]"<<this->gp_abundance_v.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance_v.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock' "<<endl;
#elif defined  _USE_MASS_AS_PRIMARY_OBSERVABLE_
      for(int i=0; i<this->tracer.mass_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.MBmin[i]+(i+0.5)*(this->tracer.MBmax[i]-this->tracer.MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.mass_function[i])));
      this->gp_abundance_v<<"plot"<<this->gp_abundance_v.file1d(xy_pts_ref_m) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance_v.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock'"<<endl;
      xy_pts_ref_m.clear(); xy_pts_ref_m.shrink_to_fit();
#endif
      xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
      xy_pts_v.clear(); xy_pts_v.shrink_to_fit();
#endif
      real_prec residuals=0;
      ULONG count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
      for(ULONG i=0;i<this->params._NMASSbins_mf() ;++i)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        if(this->tracer.vmax_function[i]>0)
#elif defined  _USE_MASS_AS_PRIMARY_OBSERVABLE_
          if(this->tracer.mass_function[i]>0)
#endif
        {
          count_b++;
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
          residuals+= fabs(this->tracer_ref.vmax_function[i]/this->tracer.vmax_function[i] -1.0 );
#elif defined  _USE_MASS_AS_PRIMARY_OBSERVABLE_
          residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i] -1.0 );
#endif
        }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        So.message_screen("Residuals from vmax-function = ", residuals, "%");
#elif defined  _USE_MASS_AS_PRIMARY_OBSERVABLE_
        So.message_screen("Residuals from mass-function = ", residuals, "%");
#endif
        std::cout<<endl;
#endif
  }//closes  if(true==initial_assignment)
  else if(false==initial_assignment)
    {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     if(h_property==_MASS_)
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
  if(h_property==_VMAX_)
#endif
      {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        this->tracer.params.set_i_vmax_g(NEGATIVE_INT); //this avoids the measurement of the vmax function again
        this->tracer.params.set_i_mass_g(POSITIVE_INT); //this allow the tracer to get the mass function
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
        this->tracer.params.set_i_vmax_g(POSITIVE_INT); //this allowsthe measurement of the vmax function
        this->tracer.params.set_i_mass_g(NEGATIVE_INT); //this avoids the tracer to get the mass function
#endif
        this->tracer.params.set_i_rs_g(NEGATIVE_INT); //no rs information here
        this->tracer.params.set_i_spin_bullock_g(NEGATIVE_INT); //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
        this->tracer_ref.get_property_function(fname_mass_function_Y_ref);
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
        this->gp_abundance<<"set log x \n";
        this->gp_abundance<<"set border linewidth 2.2\n";
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
        this->gp_abundance<<"set xlabel 'Vrms [km / h]' font 'Times-Roman,15'\n";
        this->gp_abundance<<"set ylabel 'log n(Vrms) ' font 'Times-Roman,15'\n";
#elif defined  _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      this->gp_abundance<<"set xlabel 'Mass [Ms / h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      this->gp_abundance<<"set ylabel 'log n(M) ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
#endif
        vector<pair<real_prec, real_prec> > xy_pts_v;
        xy_pts_ref_v.clear();xy_pts_ref_v.shrink_to_fit();
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
        for(int i=0; i<this->tracer.vmax_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.VMAXBmin[i]+(i+0.5)*(this->tracer.VMAXBmax[i]-this->tracer.VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.vmax_function[i])));
        for(int i=0; i<this->tracer.vmax_function.size(); ++i)
          xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.VMAXBmin[i]+(i+0.5)*(this->tracer_ref.VMAXBmax[i]-this->tracer_ref.VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.vmax_function[i])));
        this->gp_abundance<<"plot"<<this->gp_abundance.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock'"<<endl;
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        for(int i=0; i<this->tracer.mass_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.MBmin[i]+(i+0.5)*(this->tracer.MBmax[i]-this->tracer.MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.mass_function[i])));
        for(int i=0; i<this->tracer.mass_function.size(); ++i)
          xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.MBmin[i]+(i+0.5)*(this->tracer_ref.MBmax[i]-this->tracer_ref.MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.mass_function[i])));
        this->gp_abundance<<"plot"<<this->gp_abundance.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock'"<<endl;
#endif
        xy_pts_v.clear(); xy_pts_v.shrink_to_fit();
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
#endif
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(int i=0;i<this->params._NMASSbins_mf() ;++i)
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
          if(this->tracer.vmax_function[i]>0)
#elif defined  _USE_VMAX_AS_PRIMARY_OBSERVABLE_
            if(this->tracer.mass_function[i]>0)
#endif
                {
                  count_b++;
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
              residuals+= fabs(this->tracer_ref.vmax_function[i]/this->tracer.vmax_function[i]-1.0);
#elif defined  _USE_VMAX_AS_PRIMARY_OBSERVABLE_
              residuals+= fabs(this->tracer_ref.mass_function[i]/this->tracer.mass_function[i]-1.0);
#endif
            }
          residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
          So.message_screen("Residuals from vmax-function = ", residuals, "%");
#elif defined  _USE_VMAX_AS_PRIMARY_OBSERVABLE_
          So.message_screen("Residuals from mass-function = ", residuals, "%");
#endif
#endif
      }
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
    else if(h_property==_RS_)
      {
        this->tracer.params.set_i_vmax_g(NEGATIVE_INT); //this allows the avoid the measurement of the vmax function again
        this->tracer.params.set_i_mass_g(NEGATIVE_INT); //this allow aviud the tracer to get the mass function
        this->tracer.params.set_i_rs_g(POSITIVE_INT); //this allow the tracer to get the mass function
        this->tracer.params.set_i_spin_g(NEGATIVE_INT); //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
#ifdef _USE_GNUPLOT_ABUNDANCE_RS_PLOT_
        this->gp_abundance_rs<<"set log x \n";
        this->gp_abundance_rs<<"set border linewidth 2.2\n";
        this->gp_abundance_rs<<"set xlabel 'Rs ' font 'Times-Roman,15'\n";
        this->gp_abundance_rs<<"set ylabel 'log n(Rs) ' font 'Times-Roman,15'\n";
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
        vector<pair<real_prec, real_prec> > xy_pts_v;
        for(ULONG i=0; i<this->tracer.rs_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.RSBmin[i]+(i+0.5)*(this->tracer.RSBmax[i]-this->tracer.RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.rs_function[i])));
        for(ULONG i=0; i<this->tracer_ref.rs_function.size(); ++i)
          xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.RSBmin[i]+(i+0.5)*(this->tracer_ref.RSBmax[i]-this->tracer_ref.RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.rs_function[i])));
        this->gp_abundance_rs<<"plot"<<this->gp_abundance_rs.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance_rs.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock' textcolor rgb '"<<FG_COLOR<<"'"<<endl;
        xy_pts_v.clear(); xy_pts_v.shrink_to_fit();
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
#endif
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(ULONG i=0;i<this->params._NMASSbins_mf() ;++i)
          if(this->tracer.rs_function[i]>0)
            {
              count_b++;
              residuals+= fabs(this->tracer_ref.rs_function[i]/this->tracer.rs_function[i]-1.0);
            }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from Rs-function = ", residuals, "%");
        std::cout<<endl;
#endif 
      }
#endif //closes _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
    else if(h_property==_CONCENTRATION_)
      {
        this->tracer.params.set_i_vmax_g(NEGATIVE_INT); //this allows the avoid the measurement of the vmax function again
        this->tracer.params.set_i_mass_g(NEGATIVE_INT); //this allow aviud the tracer to get the mass function
        this->tracer.params.set_i_rs_g(POSITIVE_INT); //this allow the tracer to get the mass function
        this->tracer.params.set_i_spin_g(NEGATIVE_INT); //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
#ifdef _USE_GNUPLOT_ABUNDANCE_CONCENTRATION_PLOT_
        this->gp_abundance_rs<<"set log x \n";
        this->gp_abundance_rs<<"set border linewidth 2.2\n";
        this->gp_abundance_rs<<"set xlabel 'Cvir ' font 'Times-Roman,15'\n";
        this->gp_abundance_rs<<"set ylabel 'log n(Cvir) ' font 'Times-Roman,15'\n";
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
        vector<pair<real_prec, real_prec> > xy_pts_v;
        for(ULONG i=0; i<this->tracer.rs_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.CONCENTRATIONBmin[i]+(i+0.5)*(this->tracer.CONCENTRATIONBmax[i]-this->tracer.CONCENTRATIONBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.rs_function[i])));
        for(ULONG i=0; i<this->tracer_ref.rs_function.size(); ++i)
          xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.CONCENTRATIONBmin[i]+(i+0.5)*(this->tracer_ref.CONCENTRATIONBmax[i]-this->tracer_ref.RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.rs_function[i])));
        this->gp_abundance_rs<<"plot"<<this->gp_abundance_rs.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance_rs.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock' textcolor rgb '"<<FG_COLOR<<"'"<<endl;
        xy_pts_v.clear(); xy_pts_v.shrink_to_fit();
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
#endif
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(ULONG i=0;i<this->params._NMASSbins_mf() ;++i)
          if(this->tracer.cvir_function[i]>0)
            {
              count_b++;
              residuals+= fabs(this->tracer_ref.cvir_function[i]/this->tracer.cvir_function[i]-1.0);
            }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from Cvir-function = ", residuals, "%");
        std::cout<<endl;
#endif 
      }
#endif //closes _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
      else if(h_property==_SPIN_ || h_property==_SPIN_BULLOCK_)
       {
        this->tracer.params.set_i_vmax_g(NEGATIVE_INT); //this allows the avoid the measurement of the vmax function again
        this->tracer.params.set_i_mass_g(NEGATIVE_INT); //this allow aviud the tracer to get the mass function
        this->tracer.params.set_i_rs_g(NEGATIVE_INT); //this allow the tracer to get the mass function
        if(h_property==_SPIN_)
          this->tracer.params.set_i_spin_g(POSITIVE_INT); //this allow the tracer to get the mass function
        if(h_property==_SPIN_BULLOCK_)
          this->tracer.params.set_i_spin_bullock_g(POSITIVE_INT); //this allow the tracer to get the mass function
        this->tracer.get_property_function(fname_mass_function_Y);
#ifdef _USE_GNUPLOT_ABUNDANCE_SPIN_PLOT_
        this->gp_abundance_spin<<"set log x \n";
        this->gp_abundance_spin<<"set border linewidth 2.2\n";
        this->gp_abundance_spin<<"set xlabel '{/Symbol l} ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
        this->gp_abundance_spin<<"set ylabel 'log n({/Symbol l}) ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
        vector<pair<real_prec, real_prec> > xy_pts_v;
        for(ULONG i=0; i<this->tracer.s_function.size(); ++i)
          xy_pts_v.push_back(std::make_pair(this->tracer.SPINBmin[i]+(i+0.5)*(this->tracer.SPINBmax[i]-this->tracer.SPINBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.s_function[i])));
        for(ULONG i=0; i<this->tracer_ref.s_function.size(); ++i)
          xy_pts_ref_v.push_back(std::make_pair(this->tracer_ref.SPINBmin[i]+(i+0.5)*(this->tracer_ref.SPINBmax[i]-this->tracer_ref.SPINBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer_ref.s_function[i])));
        this->gp_abundance_spin<<"plot"<<this->gp_abundance_spin.file1d(xy_pts_ref_v) << "w l lw 3 lt 2 title 'Reference', "<<this->gp_abundance_spin.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Mock'"<<endl;
        xy_pts_v.clear(); xy_pts_v.shrink_to_fit();
        xy_pts_ref_v.clear(); xy_pts_ref_v.shrink_to_fit();
#endif
        real_prec residuals=0;
        int count_b=0;
#pragma omp parallel for reduction(+:count_b,residuals)
        for(int i=0;i<this->params._NMASSbins_mf() ;++i)
          if(this->tracer.s_function[i]>0)
            {
              count_b++;
              residuals+= fabs(this->tracer_ref.s_function[i]/this->tracer.s_function[i]-1.0);
            }
        residuals/=static_cast<real_prec>(count_b)/100.0;
#ifdef _FULL_VERBOSE_
        So.message_screen("Residuals from Spin-function = ", residuals, "%");
        std::cout<<endl;
#endif
      }
#endif  // closes _USE_SPIN_AS_DERIVED_OBSERVABLE_
//  string rfile=this->params._Output_directory()+"ref_ncat.txt";
//  this->tracer_ref.select_random_subsample(0.2, rfile);
  string bfile=this->params._Output_directory()+"new_cat_with_dm_";
#ifdef  _USE_CWC_
     bfile+="cwc_";
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
     bfile+="indiv_bias_";
#endif
#ifdef _USE_MACH_NUMBER_
     bfile+="mach5_";
#endif
#ifdef  _USE_LOCAL_OVERDENSITY_
     bfile+="delta5_";
#endif
#ifdef  _USE_TIDAL_ANISOTROPY_SEC_PROP_
    bfile+="ta_";
#endif
#ifdef _MULTISCALE_
     bfile+="MS_";
#endif
    bfile+=".txt";
   if(h_property==_SPIN_BULLOCK_) // es la última
      this->tracer.select_random_subsample(0.2, bfile);

    } // closes else if(false==initial_assignment)
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::sample_mock()
   {
     So.enter(__PRETTY_FUNCTION__);
     this->So.message_screen("Sampling halo density field with random positions within cells");
     const gsl_rng_type *  T;
     gsl_rng * r ;
     gsl_rng_env_setup();
     gsl_rng_default_seed=1152;
     T = gsl_rng_ranlux;
     r = gsl_rng_alloc (T);
     string pname;
     if(this->iteration <=this->params._N_iterations_Kernel())
       pname ="_iteration"+to_string(this->iteration);
     else
       pname= "_realization"+to_string(this->params._realization());
     string fname=this->params._Output_directory()+"MOCK_TR"+pname+"_CAT_z"+to_string(this->params._redshift())+".txt";
     vector<real_prec>x_cart;
     vector<real_prec>y_cart;
     vector<real_prec>z_cart;
     real_prec delta=this->params._Lbox()/static_cast<real_prec>(this->params._Nft());
     for(ULONG i=0;i<this->delta_Y_new.size();++i)
       {
         int Ngal=this->delta_Y_new[i];
         int counter=0;
         ULONG x_min, y_min, z_min;
         index2coords(i,this->params._Nft(), x_min, y_min, z_min);
         x_min*=delta;
         y_min*=delta;
         z_min*=delta;
         while(counter<=Ngal)
           {
             x_cart.push_back(x_min+delta*gsl_rng_uniform (r));
             y_cart.push_back(y_min+delta*gsl_rng_uniform (r));
             z_cart.push_back(z_min+delta*gsl_rng_uniform (r));
             counter++;
           }
       }
     So.message_screen("Sampled with",x_cart.size(),"objects");
     this->File.write_to_file(fname, x_cart,y_cart,z_cart);
     So.DONE();
   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _SEVERAL_REAL_
void BiasMT::get_pdf(int real)
#elif !defined _SEVERAL_REAL_
void BiasMT::get_pdf()
#endif
 {
    So.enter(__PRETTY_FUNCTION__);
    this->So.message_screen("BAM in BIAS mode");
    this->So.message_screen("Computing statistics from input density fields");
    // The int sua indicates 0, 1, 2.,  ... n_cwt usados. SI tomamos knots y el resto, sua =0, 1
    string type_X=this->params._Scale_X();
    string type_Y=this->params._Scale_Y();
    // real_prec num_in_log_x = true==this->params._Convert_Density_to_Delta_X() ? NUM_IN_LOG: 0.;
    // real_prec num_in_log_y = true==this->params._Convert_Density_to_Delta_Y() ? NUM_IN_LOG: 0.;
    real_prec num_in_log_x=0;
    if(type_X=="log")
      num_in_log_x = NUM_IN_LOG;
    real_prec num_in_log_y=0;
    if(type_Y=="log")
      num_in_log_y = NUM_IN_LOG;
#ifdef _SEVERAL_REAL_
     string prop_X="_X"+this->params._XNAME()+"_ScaleX"+type_X+"_MASX"+to_string(this->params._iMAS_X())+"_Realization"+to_string(real);
     string prop_Y="Y"+this->params._YNAME()+"_ScaleY"+type_Y+"_MASY"+to_string(this->params._iMAS_Y());
#else
     string prop_X="_X"+this->params._XNAME()+"_ScaleX"+type_X+"_MASX"+to_string(this->params._iMAS_X());
     string prop_Y="Y"+this->params._YNAME()+"_ScaleY"+type_Y+"_MASY"+to_string(this->params._iMAS_Y())+"_"+this->params._extra_info();
#endif
#ifdef _USE_CWC_
     this->cwclass.get_CWC(this->delta_X_ini); // CWClass done with the full delta
#ifdef _USE_MASS_KNOTS_
     this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#elif !defined _USE_CWC_
#if defined (_USE_MASS_KNOTS_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_AWEB_) || defined (_USE_PWEB_) || defined (_USE_TIWEB_)
     this->cwclass.get_CWC(this->delta_X_ini); // CWClass done with the full delta
#ifdef _USE_MASS_KNOTS_
     this->cwclass.get_Mk_collapsing_regions(this->delta_X_ini,  static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#endif
#endif
#if (defined (_USE_NABLA2DELTA_) || defined (_USE_S2_)  || defined (_USE_S3_) || defined (_USE_S2DELTA_)) && (!defined (_USE_CWC_)) && !(defined (_USE_PWEB_))
     this->cwclass.get_bias_terms(this->delta_X_ini);
#endif
     vector<real_prec> delta_Y_ini( this->params._NGRID(),0);
     delta_Y_ini=this->delta_Y;
     vector<real_prec> delta_X_ini( this->params._NGRID(),0);
     delta_X_ini=this->delta_X;
     int sua=0;
#ifdef _USE_CWC_
     for(sua=0;sua<this->params._n_cwt();++sua) // Loop over the different CWC requested in the parameter file
       {
#endif
         this->tstruct=0;
         string file=this->params._Output_directory()+"2dbias"+prop_X+"_"+prop_Y+"Nft"+to_string(this->params._Nft());
#ifdef _USE_CWC_
         this->tstruct=this->cwclass.cwt_used[sua];
         file=this->params._Output_directory()+"2dbias"+prop_X+"_"+prop_Y+"Nft"+to_string(this->params._Nft())+"_CW"+to_string(this->tstruct);
#endif
//         string file=this->params._Output_directory()+"2dbias"+prop_X+"_"+prop_Y+"_Nft"+to_string(this->params._Nft())+"_z"+to_string(this->params._redshift())+"_LambdaTH"+to_string(lambdath)+"_CW"+to_string(this->tstruct);
         // So.message_screen("Computing variances for cosmic web type ",this->tstruct);
         // real_prec var_XX=gsl_stats_variance(&this->delta_X[0],1,  this->params._NGRID());
         // real_prec var_YY=gsl_stats_variance(&this->delta_Y[0],1,  this->params._NGRID());
         // real_prec corr_XY=gsl_stats_correlation(&this->delta_X[0],1, &this->delta_Y[0], 1,  this->params._NGRID());
         // So.message_screen("Variance DM =",var_XX);
         // So.message_screen("Variance TR =",var_YY);
         // So.message_screen("Correlation X-Y =",corr_XY);
         this->get_min_max_X_Y();// THIS CAN BE REPLACED by finding min and max at converting time
         if(this->tstruct!=0)
           {
             // Here we need to redefine the overdensities according to the mean for each CWT
#ifdef _FULL_VERBOSE_
             So.message_screen("Transforming to overdensities for CWT",this->tstruct);
#endif
             ULONG new_nobjects=0;
             real_prec beta=1.0;
             if(true==this->params._Convert_Density_to_Delta_Y())
               {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:new_nobjects)
#endif
                 for(ULONG i = 0; i <  this->params._NGRID(); ++i)
#ifdef _USE_CWC_
                   if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                     new_nobjects+=(this->N_objects_Y/static_cast<real_prec>( this->params._NGRID()))*(1.0+delta_Y_ini[i]);
                 beta=static_cast<real_prec>(new_nobjects)/static_cast<real_prec>(this->N_objects_Y);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                 for(ULONG i = 0; i <  this->params._NGRID(); ++i)
                   this->delta_Y[i]=(1./beta)*delta_Y_ini[i]-(beta-1)/beta;
                 So.message_screen("New nobs",  new_nobjects);
                 So.message_screen("Original nobs",  this->N_objects_Y);

                 So.message_screen("Factor beta for TR", beta);
               }
             if(true==this->params._Convert_Density_to_Delta_X())
               {
                 new_nobjects=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:new_nobjects)
#endif
                 for(ULONG i = 0; i <  this->params._NGRID(); ++i)
#ifdef _USE_CWC_
                   if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                     new_nobjects+=(this->N_objects_X/static_cast<real_prec>( this->params._NGRID()))*(1.0+delta_X_ini[i]);
                 beta=static_cast<real_prec>(new_nobjects)/static_cast<real_prec>(this->N_objects_X);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                 for(ULONG i = 0; i <  this->params._NGRID(); ++i)
                   this->delta_X[i]=(1./beta)*this->delta_X[i]-(beta-1)/beta;
                 So.message_screen("Factor beta for DM", beta);
               }
           }
         this->So.message_screen("Minimum of delta X", get_min(delta_X));
         this->So.message_screen("Maximum of delta X", get_max(delta_X));
         this->So.message_screen("Minimum of delta Y", get_min(delta_Y));
         this->So.message_screen("Maximum of delta Y", get_max(delta_Y));
         if(type_Y=="log")
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(ULONG i = 0; i <  this->params._NGRID(); ++i)
#ifdef _USE_CWC_
               if(true==this->cwclass.get_cell_classified(sua, i))
#endif
                 this->delta_Y[i]=  log10(num_in_log_y+static_cast<real_prec>(delta_Y_ini[i]));
#ifdef _USE_CWC_
               else
                 delta_Y[i]=NOCELL;
#endif
           }
#ifdef _USE_CWC_
         else
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i <  this->params._NGRID() ;++i)
               if(false!=this->cwclass.get_cell_classified(sua, i))
                 delta_Y[i]=NOCELL;
           }
#endif
         if(type_X=="log")
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i <  this->params._NGRID() ;++i)
               {
#ifdef _USE_CWC_
                 if(true==this->cwclass.get_cell_classified(sua, i))
                   {
#endif
                     delta_X[i]= log10(num_in_log_x+static_cast<real_prec>(delta_X_ini[i]));
#ifdef _USE_CWC_
                   }
                 else
                   delta_X[i]=NOCELL;
#endif
               }
           }
#ifdef _USE_CWC_
         else
           {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
             for(ULONG i = 0;i <  this->params._NGRID() ;++i)
               if(false!=this->cwclass.get_cell_classified(sua, i))
                 delta_X[i]=NOCELL;
           }
#endif
         if(type_Y=="log")
           {
             this->So.message_screen("Minimum of log10 1 + delta Y", get_min(delta_Y));
             this->So.message_screen("Maximum of log10 1 + delta Y", get_max(delta_Y));
           }
         if(type_X=="log")
           {
             this->So.message_screen("Minimum of log10 1 + delta X", get_min(delta_X));
             this->So.message_screen("Maximum of log10 1 + delta X", get_max(delta_X));
           }
         // ******************************************************************************
#ifdef _FULL_VERBOSE_
         So.message_screen("Getting X and Y bins used for the bias information");
#endif
         // ******************************************************************************
         // Aca redefinimos los bines en X y Y.
         // Hacemos una distinición importante. Cuando usamos NGP,
         // los bines estarán identificados con el número de particulas
         // Para otros tipos de interpolaciones, hacemos bines propiamente.
         // nmax_x_onecell is the maximum number of objects in one cell.
         if(this->params._iMAS_X()==0)
           {
             if(false==params._Convert_Density_to_Delta_X())
               if(type_X=="linear")
                 {
                   this->new_nbins_x = this->nmax_x_onecell+1;
                   this->Xmin=0;
                   this->Xmax=static_cast<int>(nmax_x_onecell);
                 }
               else{
                 this->new_nbins_x = this->params._NX();
                 this->Xmin=this->params._ldelta_X_min();
                 this->Xmax=this->params._ldelta_X_max();
               }
             else
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x =this->params._NX();
                     this->Xmin=this->params._delta_X_min();
                     this->Xmax=this->params._delta_X_max();
                   }
                 else{
                   this->new_nbins_x = this->params._NX();
                   this->Xmin=this->params._ldelta_X_min();
                   this->Xmax=this->params._ldelta_X_max();
                 }
               }
           }
         else
           {
             if(true==params._Convert_Density_to_Delta_X())
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x = this->params._NX();
                     this->Xmin=this->params._delta_X_min();
                     this->Xmax=this->params._delta_X_max();
                   }
                 else if(type_X=="log")
                   {
                     this->new_nbins_x = this->params._NX();
                     this->Xmin=this->params._ldelta_X_min();
                     this->Xmax=this->params._ldelta_X_max();
                   }
               }
             else
               {
                 if(type_X=="linear")
                   {
                     this->new_nbins_x = this->params._NX();
                     this->Xmin=this->params._delta_X_min();
                     this->Xmax=this->params._delta_X_max();
                   }
                 else if(type_X=="log")
                   {
                     this->new_nbins_x = this->params._NX();
                     this->Xmin=this->params._ldelta_X_min();
                     this->Xmax=this->params._ldelta_X_max();
                   }
               }
           }
         // ******************************************************************************
         if(this->params._iMAS_Y()==0)
           {
             if(false==params._Convert_Density_to_Delta_Y())
               if(type_Y=="linear")
                 {
                   this->new_nbins_y = this->nmax_y_onecell+1;
                   this->Ymin=0;
                   this->Ymax=static_cast<int>(nmax_y_onecell);
                 }
               else{
                 this->new_nbins_y = this->params._NY();
                 this->Ymin=this->params._ldelta_Y_min();
                 this->Ymax=this->params._ldelta_Y_max();
               }
             else
               {
                 if(type_Y=="linear")
                   {
                     this->new_nbins_y =this->params._NY();
                     this->Ymin=this->params._delta_Y_min();
                     this->Ymax=this->params._delta_Y_max();
                   }
                 else{
                   this->new_nbins_y = this->params._NY();
                   this->Ymin=this->params._ldelta_Y_min();
                   this->Ymax=this->params._ldelta_Y_max();
                 }
               }
           }
         else{  // if CIC or anything higher
           if(true==params._Convert_Density_to_Delta_Y())
             {
               if(type_Y=="linear")
                 {
                   this->Ymin=this->params._delta_Y_min();
                   this->Ymax=this->params._delta_Y_max();
                 }
               else if(type_Y=="log")
                 {
                   this->Ymin=this->params._ldelta_Y_min();
                   this->Ymax=this->params._ldelta_Y_max();
                 }
             }
           else
             {
               if(type_Y=="linear")
                 this->Ymax=this->params._delta_Y_max();
               else if(type_Y=="log")
                 this->Ymax=this->params._ldelta_Y_max();
               if(type_Y=="linear")
                 this->Ymin=this->params._delta_Y_min();
               else if(type_Y=="log")
                 this->Ymin=this->params._ldelta_Y_min();
           }
           this->new_nbins_y = this->params._NY();
         }
         So.DONE();
/*
         this->Ymax=get_max(this->delta_Y);
         this->Ymin=get_min(this->delta_Y);
         this->Xmax=get_max(this->delta_X);
         this->Xmin=get_min(this->delta_X);
*/
         // ******************************************************************************
         this->DELTAX=(this->Xmax-this->Xmin)/static_cast<real_prec>(this->params._NX());
         this->DELTAY=(this->Ymax-this->Ymin)/static_cast<real_prec>(this->new_nbins_y);
#ifdef _nomas_
         this->So.message_screen("Minimum  Y", this->Ymin);
         this->So.message_screen("Maximum  Y", this->Ymax);
         this->So.message_screen("Minimum  X", this->Xmin);
         this->So.message_screen("Maximum  X", this->Xmax);
#endif
         this->X_bins.resize(new_nbins_x,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(ULONG i = 0; i < new_nbins_x; ++i)
           X_bins[i]= (this->params._iMAS_X() == 0 && false==this->params._Convert_Density_to_Delta_X() && type_X=="linear") ? static_cast<real_prec>(i) :  this->Xmin+(i+0.5)*(this->Xmax-this->Xmin)/static_cast<real_prec>(new_nbins_x);
         ofstream sx; sx.open("xbins.txt");
         for(ULONG i = 0; i < new_nbins_x; ++i)
           sx<<i<<"  "<<X_bins[i]<<endl;
         sx.close();
         this->Y_bins.resize(new_nbins_y,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
//         for(ULONG i = 0; i < new_nbins_y; ++i)
 //          Y_bins[i]= (iMAS_Y == 0  && false==this->params._Convert_Density_to_Delta_Y()  && type_Y=="linear") ? static_cast<real_prec>(i) :  this->Ymin+(i+0.5)*(this->Ymax-this->Ymin)/static_cast<real_prec>(new_nbins_y);

         for(ULONG i = 0; i < new_nbins_y; ++i)
           Y_bins[i]=   this->Ymin+(i+0.5)*(this->Ymax-this->Ymin)/static_cast<real_prec>(new_nbins_y);
         ofstream sy; sy.open("ybins.txt");
         for(ULONG i = 0; i < new_nbins_y; ++i)sy<<i<<"  "<<Y_bins[i]<<endl;
         sy.close();
         // ************************RESIZE VECTORS FOR HISTOGRAMS*************************
         // Resize the vector to allocate the number of cells in the Den-Nh bins
         // This vectors will contain P(X,Y) only used for analyizing , not to create mocks
         BIAS_NCOUNTS.resize(this->params._n_sknot_massbin()*this->params._n_cwt()*new_nbins_x*new_nbins_y*N_REDSHIFT_BINS, 0);
         BIAS_NCOUNTS_normalized.resize(this->params._n_sknot_massbin()*this->params._n_cwt()*new_nbins_x*new_nbins_y*N_REDSHIFT_BINS, 0);
#ifdef _FULL_VERBOSE_
         So.message_screen("Computing BIAS(X,Y, CWT, ...) ");
#endif
#ifdef _USE_REDSHIFT_MASK_
         vector<real_prec> redshift_mask( this->params._NGRID(),0);
         this->File.read_array(this->params._Name_redshift_mask(),redshift_mask);
         std::cout<<"Max redshift in mask = "<<get_max(redshift_mask)<<endl;
         std::cout<<"Min redshift in mask = "<<get_min(redshift_mask)<<endl;
         vector<int> redshift_mask_bin( this->params._NGRID(),0);
         vector<int> ncells_zbin(N_REDSHIFT_BINS,0);
#endif
#ifdef _USE_BINARY_MASK_
         vector<real_prec> binary_mask( this->params._NGRID(),0);
         this->File.read_array(this->params._Name_binary_mask(), binary_mask);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(ULONG ig = 0; ig<  this->params._NGRID() ; ++ig)
           {
             int redshift_bin=0;
             real_prec redshift_mask_cell=0;
#ifdef _USE_REDSHIFT_MASK_
#ifdef _USE_BINARY_MASK_
             if(static_cast<real_prec>(binary_mask[ig])>0)
               {
#endif
                 redshift_mask_cell=redshift_mask[ig];

                 if(redshift_mask_cell<REDSHIFT_MAX && redshift_mask_cell>=REDSHIFT_MIN)
                   {
                     redshift_bin = static_cast<int>(floor((redshift_mask_cell-REDSHIFT_MIN)/DELTA_Z));
                     redshift_mask_bin[ig]=redshift_bin;


#pragma omp atomic update
                     ncells_zbin[redshift_bin]++;

#endif

                     int ict=0;
#ifdef _USE_MASS_KNOTS_
                     ict=this->cwclass.SKNOT_M_info[ig];
#endif
#ifdef _USE_CWC_
                     if(true==this->cwclass.get_cell_classified(sua, ig))
                       {
#endif
                         if(delta_X[ig]!= NOCELL && delta_Y[ig]!= NOCELL)   // This means, use cells being classified.
                           if((delta_X[ig]>=this->Xmin && delta_X[ig]<=this->Xmax) && (delta_Y[ig]>=this->Ymin && delta_Y[ig]<=this->Ymax))
                             {
                               int IX= get_bin(delta_X[ig], this->Xmin, this->params._NX(),this->DELTAX, false);
                               int IY= get_bin(delta_Y[ig], this->Ymin, this->params._NY(),this->DELTAY, false);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                              this->BIAS_NCOUNTS[index_5d(IX, IY, sua, ict, redshift_bin, new_nbins_y, this->params._n_cwt(), this->params._n_sknot_massbin(), N_REDSHIFT_BINS)]++;
                           }
#ifdef _USE_CWC_
                       }
#endif

#ifdef _USE_REDSHIFT_MASK_
                   }
#ifdef _USE_BINARY_MASK_
               }//close if binary mask
#endif

#endif
           }
#ifdef _USE_REDSHIFT_MASK_
         ULONG NCELLS=0;
         for(int i=0;i<ncells_zbin.size();++i)
           NCELLS+=ncells_zbin[i];

         std::cout<<"Used "<<NCELLS<<" out of "<< this->params._NGRID()<<"  ("<<100.*abs(( this->params._NGRID()-NCELLS)/static_cast<real_prec>( this->params._NGRID()))<<" %)"<<endl;

         redshift_mask.clear();
         redshift_mask.shrink_to_fit();
#endif
         // DO THIS ONLY FOR KNOTS, SINCE THIS IS A LOOP OVER THE MKNOTS
         //Now I normalize the 4d array such that I can also do contours for the knots in different bins of MK
         int iz=0;
#ifdef _USE_REDSHIFT_MASK_
         for(iz=0;iz<N_REDSHIFT_BINS;++iz)
           {
#endif
             int mk=0;
#ifdef _USE_MASS_KNOTS_
             for(mk=0;mk<this->params._n_sknot_massbin();++mk)
               {
#endif
                 vector<LONG>aux_v(new_nbins_x*new_nbins_y, 0);
                 for(int lx=0;lx<new_nbins_x;++lx)
                   for(int ly=0;ly<new_nbins_y;++ly)
                     aux_v[index_2d(lx,ly,new_nbins_y)]=this->BIAS_NCOUNTS[index_5d(lx,ly,sua,mk,iz,new_nbins_y, this->params._n_cwt(), this->params._n_sknot_massbin(), N_REDSHIFT_BINS)];
                 ULONG lkk=get_max<LONG>(aux_v);
                 vector<real_prec>aux_n(new_nbins_x*new_nbins_y, 0);
                 for(int i=0;i< aux_v.size() ;++i)
                   aux_n[i]= lkk==0 ? 0 : static_cast<real_prec>(aux_v[i])/static_cast<real_prec>(lkk);
                 aux_v.clear(); aux_v.shrink_to_fit();

                 So.message_screen("Writting the joint distribution P(X,Y)");
                 vector<vector<real_prec> > Vaux;
                 Vaux.resize(new_nbins_x);
                 for(int i=0;i<new_nbins_x;++i)
                   for(int j=0;j<new_nbins_y;++j)
                     Vaux[i].push_back(aux_n[index_2d(i,j,new_nbins_y)]);
                 aux_n.clear(); aux_n.shrink_to_fit();
                 string filek=file;
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                 filek+="_MKbin"+to_string(mk);
#endif

#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                 filek+="_zbin"+to_string(iz);
#endif

#if  defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                 filek+="_zbin"+to_string(iz)+"_MKbin"+to_string(mk);
#endif
                 this->File.write_to_file(filek+".txt",X_bins,Y_bins,Vaux);
                 mcmc.get_contour_levels(filek+"_contour_levels",Vaux);
                 Vaux.clear(); Vaux.shrink_to_fit();
#ifdef _USE_MASS_KNOTS_
               }
#endif
#ifdef _USE_REDSHIFT_MASK_
           }
#endif
         So.message_screen("Allocating the values of Y found in the bins of X.");
         // if(this->cwt_used[sua]==I_KNOT)
         //   {
#ifdef _nomas_
#ifdef _USE_REDSHIFT_MASK_
         for(int iz=0;iz<N_REDSHIFT_BINS;++iz)
           {
#endif
#ifdef _USE_MASS_KNOTS_
             for(int mk=0;mk<this->params._n_sknot_massbin();++mk)
               {
#endif
                 vector<vector<real_prec> >DELTAXY;
                 DELTAXY.resize(this->new_nbins_x);

                 // Get the values of Y in each bin of X
                 for(ULONG i=0;i< this->params._NGRID() ;++i)
                   {
#ifdef _USE_BINARY_MASK_
                     if(static_cast<real_prec>(binary_mask[i])>0)
                       {
#endif
#ifdef _USE_MASS_KNOTS_
                             if(this->cwclass.SKNOT_M_info[i]==mk)
                               {
#endif
#ifdef _USE_REDSHIFT_MASK_
                                 if(redshift_mask_bin[i]==iz)
                                   {
#endif
                                    // This last condition is important, specially when log10(0) is involved.
                                     // For Y, we ask to do statistics with cells that have at least ione halo, the > imposes it
                                     {
                                      if((delta_X[i]>=this->Xmin && delta_X[i]<=this->Xmax) && (delta_Y[i]>=this->Ymin && delta_Y[i]<=this->Ymax))
                                      {
                                        int IX=get_bin(delta_X[i],this->Xmin,this->new_nbins_x,DELTAX,this->bin_accumulate_borders);
                                         DELTAXY[IX].push_back(delta_Y[i]);// IX ES EL BIN DE X: allocate all values of Y in the bins of dm.
                                       }
                                      }
#ifdef _USE_REDSHIFT_MASK_
                                   }
#endif
#ifdef _USE_MASS_KNOTS_
                               }
#endif
#ifdef _USE_REDSHIFT_MASK_
                           }
#endif
                   }
                 So.message_screen("Done");
                 this->mean_Y.clear();
                 this->mean_Y.shrink_to_fit();
                 this->mean_Y.resize(new_nbins_x,0);
                 vector<real_prec> var_Y(new_nbins_x,0);
                 vector<real_prec> skew_Y(new_nbins_x,0);
                 vector<real_prec> kurt_Y(new_nbins_x,0);

                 So.message_screen("Getting the PDF  P(Y|X) (mean, rms, higher moments).");

                 for(int i=0;i<new_nbins_x;++i)// loop sobre bines de dark matter
                   {
                     // Ddeefine this vector to perform some statistics
                     vector<real_prec> vEPSILON;
                     for(int j=0;j<DELTAXY[i].size();++j)
                       vEPSILON.push_back(DELTAXY[i][j]);


//                     for(int j=0;j<vEPSILON.size();++j)
 //                       std::cout<<vEPSILON[j]<<endl;

                     // Mean of Y in the current X bin
                     mean_Y[i]= vEPSILON.size()>0 ?  get_mean(vEPSILON): 0 ;

                     // RMS in the current X bin
                     var_Y[i]=vEPSILON.size()> num_1 ? sqrt(get_var(vEPSILON)) : 0.  ;

                     // Skewness
                     //                  skew_Y[i]=gsl_stats_skew(&vEPSILON[0],num_1, vEPSILON.size());

                     // Kurtossis
                     //                  kurt_Y[i]=gsl_stats_kurtosis(&vEPSILON[0],num_1, vEPSILON.size());


                     // THIS IS TO WRITE AND DO HISTOGRAMS, BUT IF WE KEEP TRACK OF THE MOMENTS, WE COULD GET THE DISTRIBUTION
                     if (true==this->params._write_files_for_histograms() && vEPSILON.size()>0)
                       {
                         string file_X_bins=file+"_dist_bin"+to_string(i)+".txt";;
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_MKbin"+to_string(mk)+".txt";
#endif
#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_zbin"+to_string(iz)+".txt";
#endif
#if defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                         file_X_bins=file+"_dist_bin"+to_string(i)+"_zbin"+to_string(iz)+"_MKbin"+to_string(mk)+".txt";
#endif
                         this->sal.open(file_X_bins.c_str());
                         So.message_screen("Writting values of Y in bins of X in file ",file_X_bins);
                         for(int ie=0;ie<vEPSILON.size();++ie)sal<<vEPSILON[ie]<<endl;
                         sal.close();
                       }
                     vEPSILON.clear();
                     vEPSILON.shrink_to_fit();
                   }
                 string out2=file+"_mean"+XNAME+"bins.txt";
#if defined(_USE_MASS_KNOTS_) && !defined(_USE_REDSHIFT_MASK_)
                 out2=file+"_mean"+XNAME+"bins_MKbin"+to_string(mk)+".txt";
#endif
#if defined(_USE_REDSHIFT_MASK_) && !defined(_USE_MASS_KNOTS_)
                 out2=file+"_mean"+XNAME+"bins_zbin"+to_string(iz)+".txt";
#endif
#if  defined(_USE_REDSHIFT_MASK_) && defined(_USE_MASS_KNOTS_)
                 out2=file+"_mean"+XNAME+"bins_zbin"+to_string(iz)+"bins_MKbin"+to_string(mk)+".txt";
#endif
                 this->sal.open(out2.c_str());
                 int lnc=0;
                 for(int i=0;i<X_bins.size();++i)
                   if(mean_Y[i]!=0)
                     {
                       lnc++;
                       this->sal<<X_bins[i]<<"\t"<<mean_Y[i]<<"\t"<<var_Y[i]<<endl;//"\t"<<skew_Y[i]<<"\t"<<kurt_Y[i]<<endl;
                     }
                 So.message_screen("Wrote output in file", out2, ". Number of lines =", lnc);
                 this->sal.close();
#ifdef _USE_MASS_KNOTS_
               }
#endif
#ifdef _USE_REDSHIFT_MASK_
           }
#endif
#ifdef _USE_CWC_
       }
#endif
#endif   //edif del nomas
#ifdef _USE_CWC_
}
#endif   //edif del nomas
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_

struct s_particle_neighbours
{
  vector<ULONG>global_id_candidate;
  vector<ULONG>id_dist_candidate;
  ULONG sort_candidates(){
    ULONG rank=0;
    sort_1d_vectors<ULONG>(this->id_dist_candidate,rank);
    return this->global_id_candidate[rank];
  }

};
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MOCK_MODE
   void BiasMT::collapse_randoms()
   {
     So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
     this->So.message_screen("**Collapsing randoms towards the DM particles**");
     std::cout<<endl;
#endif
     int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
     omp_set_num_threads(NTHREADS);
     s_params_box_mas box_collps;
     box_collps.min1=this->params._xllc();
     box_collps.min2=this->params._yllc();
     box_collps.min3=this->params._zllc();
     box_collps.Lbox=this->params._Lbox();
     box_collps.Nft=this->params._Nft_random_collapse();
     box_collps.d1= box_collps.Lbox/static_cast<real_prec>(box_collps.Nft);		/* grid spacing x-direction */
     box_collps.d2= box_collps.d1;
     box_collps.d3= box_collps.d1;
     box_collps.NGRID=(box_collps.Nft*box_collps.Nft*box_collps.Nft);
     ULONG Ntracers=this->tracer._NOBJS();
     vector<int>dm_count(box_collps.NGRID,0);
     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;
     vector<real_prec> mass_dm;
     vector<ULONG> dm_id;
     vector<ULONG> dm_id_global;
     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<real_prec> mass_random;
     vector<ULONG> ran_id;
     vector<ULONG> ran_id_global;
#ifdef _COLLAPSE_RANDOMS_VELS_
     vector<real_prec> x_random_v;
     vector<real_prec> y_random_v;
     vector<real_prec> z_random_v;
#endif
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Separating DM and random:");
#endif
     ULONG N_dms=0;
     ULONG N_rand=0;
     for(ULONG i=0; i< Ntracers; ++i)
       {
         ULONG id=grid_ID(&box_collps, this->tracer.Halo[i].coord1,this->tracer.Halo[i].coord2,this->tracer.Halo[i].coord3);
         if(this->tracer.Halo[i].identity>0)
           {
             x_dm_pos.push_back(this->tracer.Halo[i].coord1);
             y_dm_pos.push_back(this->tracer.Halo[i].coord2);
             z_dm_pos.push_back(this->tracer.Halo[i].coord3);
             mass_dm.push_back(this->tracer.Halo[i].mass);
             dm_id.push_back(id);// For each dm, keep the ID of the cell in which it s located
             dm_id_global.push_back(i);
             dm_count[id]++;     // Count the numnber of DM particles in this cell
             N_dms++;            // Count the numnber of DM particles
           }
         else if(this->tracer.Halo[i].identity<0)
           {
             x_random_pos.push_back(this->tracer.Halo[i].coord1);
             y_random_pos.push_back(this->tracer.Halo[i].coord2);
             z_random_pos.push_back(this->tracer.Halo[i].coord3);
             mass_random.push_back(this->tracer.Halo[i].mass);
             ran_id.push_back(id);// For each random, keep the ID of the cell in which it s located
             ran_id_global.push_back(i);
             N_rand++;            // Count the numnber of random particles
#ifdef _COLLAPSE_RANDOMS_VELS_
             x_random_v.push_back(this->tracer.Halo[i].vel1);
             y_random_v.push_back(this->tracer.Halo[i].vel2);
             z_random_v.push_back(this->tracer.Halo[i].vel3);
#endif
           }
       }
     this->So.DONE();
     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
#ifdef _FULL_VERBOSE_
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);
#endif
     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry
#define NO_NUM -999
     ULONG Mdm_index_cell=box_collps.NGRID*max_count;
     vector<ULONG> dm_index_cell(Mdm_index_cell, NO_NUM);
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
     vector<ULONG>dm_particle_id(Mdm_index_cell,0);
#endif
#ifdef _FULL_VERBOSE_
     So.message_screen("Getting ids of dm particles in low-res cells");
#endif
     ULONG idm=0;
     dm_count.clear(); // We need to count again in order to assign index in cells
     dm_count.shrink_to_fit();
     dm_count.resize(box_collps.NGRID,0);
     for(ULONG i=0;i<Ntracers;++i) // Loop over the DM particles
     {
       if(this->tracer.Halo[i].identity>0) // Loop over the DM particles
         {
           ULONG idc=dm_id[i];         // Get the ID of the cell where this DM particle is located
           ULONG dmc=dm_count[idc];    // Retrieve the current number of dm particles (local id) in this IDcell, 0, 1, 2, ..., Ndm(cell) in the cell
           dm_index_cell[index_2d(idc, dmc, max_count)]=idm; // Associate the partial id in 0, 1, 2, ...Ndm-1 of the dm particle in the cell
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
           dm_particle_id[index_2d(idc, dmc, max_count)]=i; // Associate the global i in (0,Ntracers) of the dm particle in the cell
#endif
           dm_count[idc]++;              // Count the number of dm particle sin this particular lowres cell
           idm++;                        // Count dm particles, to get the partial id
         }
      }
     dm_id.clear();
     dm_id.shrink_to_fit();
     this->So.DONE();
     vector<ULONG>dm_index_closer_tot(N_rand,0);
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
     vector<s_particle_neighbours> dm_closest_to_ran(N_rand);
#endif
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms IN-cell");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG id_tracer=ran_id[i]; // identify the ID (low-res cell) where this random lives.
         real_prec xra=x_random_pos[i];
         real_prec yra=y_random_pos[i];
         real_prec zra=z_random_pos[i];
         vector<ULONG>i_r_to_dm_dist; // define array to allocate distance in form of index
         vector<ULONG>n_index_tot;    // define array to allocate the label of the dm particle in the cell
         ULONG N_dm_cell=dm_count[id_tracer];// Get the number of dm particles in this cell
         if(N_dm_cell>0)// if there are dm particles in the cell where this random is located, then
           for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm particles in the cell where the current random is
              {
                ULONG jdm=dm_index_cell[index_2d(id_tracer,jc,max_count)]; //indice entre (0,N_dm-1) que tiene cada particula de dm dentro de la celda id_ran
                if(jdm!=NO_NUM)
                   {
                     real_prec xda=xra-x_dm_pos[jdm];// x_random - x_dm
                     real_prec yda=yra-y_dm_pos[jdm];
                     real_prec zda=zra-z_dm_pos[jdm];
                     ULONG i_dist=_get_modulo_squared(xda,yda,zda);  //distance² between rand and each dm particle, converted to floor
                     i_r_to_dm_dist.push_back(i_dist);  // allocate the index for the distance
                     n_index_tot.push_back(jdm);        // allocat the id of the dm particle
                   }
              }
         else
           So.message_screen("No dm particles found in the cell corresponding to the random tracer ", i, ". You might want to increase the cell-size");
            // Now sort and get the identification to the closest dm particle to the current i-random
            // Sort with respect to i_r_to_dm_dist; the n_index_tot vector is sorted correspndingly such that
            // its first element is the id of the nearest dm particle to the random i.
            //The function sort_1dvectors is modified to avoid loopws and returns the ero element of the sorted array,
            // Note that here we allocate jdm in n_index_tot and  dm_particle_id[s_ind] to the vector if structures  dm_closest_to_ran[i]
            // In the loop OUTcell below we allocate instead dm_particle_id[s_ind] directly to n_index_tot in order to sort wrt distances
         ULONG rank=0;
         sort_1d_vectors<ULONG>(i_r_to_dm_dist, rank); // sort the container with distance indices
         ULONG s_ind=n_index_tot[rank];                // get the if dm in (0,.Ndm-1) particle associated to the first element of the sorted distance
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
         ULONG s_val=i_r_to_dm_dist[rank];             // get the distance index
         dm_closest_to_ran[i].id_dist_candidate.push_back(s_val); // This is the distance_id of the closest dm particle in the cell. Keep this value to compare with those from the outcell search
         dm_closest_to_ran[i].global_id_candidate.push_back(s_ind);              // keep the global id of this dm candidate
#else
         dm_index_closer_tot[i]=s_ind;// get the identification (within this cell) of the closest dm particle
#endif
       }
      So.DONE();
// Now look for in neighbour cells for closest dm particles
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying neighbouring cells");
#endif
     vector<s_nearest_cells>nearest_cells(box_collps.NGRID);
     get_neighbour_cells(box_collps.Nft,1,nearest_cells);
     So.DONE();
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms OUT-cell");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG id_tracer=ran_id[i]; // identify the ID (low-res cell) where this random lives.
         real_prec xra=x_random_pos[i];
         real_prec yra=y_random_pos[i];
         real_prec zra=z_random_pos[i];
         vector<ULONG>i_r_to_dm_dist; // define array to allocate distance in form of index
         vector<ULONG>n_index_tot;    // define array to allocate the label of the dm particle in the cell
         int N_neighbour_cells=nearest_cells[id_tracer].close_cell.size(); // Number of neighbour cells oif the current cell
         for (int ic=0; ic < N_neighbour_cells; ++ic)  // Loop over neighbour cells
           {
             ULONG id_neighbour=nearest_cells[id_tracer].close_cell[ic]; // ID of the current neighbout cell
             ULONG N_dm_cell=dm_count[id_neighbour];                     // Get the number of dm particles in this cell
             if(N_dm_cell > 0)// for this "if" I do not write an else, for if there are no dm in neighbourr cells, we do not care
              {
                for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm particles in the neighbour cell id_n
                  {
                    ULONG jdm  = dm_index_cell[index_2d(id_neighbour, jc, max_count)]; //index jdm in(0,N_dm-1) for each dm particle in the cell id_n
                    if(jdm!=NO_NUM)
                     {
                       real_prec xda = xra-x_dm_pos[jdm];
                       real_prec yda = yra-y_dm_pos[jdm];
                       real_prec zda = zra-z_dm_pos[jdm];
                       ULONG i_dist= _get_modulo_squared(xda,yda,zda); //100000 * distance² between rand and each dm particle, converted to floor
                       ULONG id_dm = dm_particle_id[index_2d(id_neighbour, jc, max_count)]; // global index in (0,N_tracer-1) for each dm particle in the cell id_n
                       i_r_to_dm_dist.push_back(i_dist);  // allocate the distance index
                       n_index_tot.push_back(id_dm);// allocate the global id of the dm particle
                     }
                  }
               }
            }
          ULONG rank=0;
          //For all distances in each neighbour cell, sort the distances and get the closest dm
          sort_1d_vectors<ULONG>(i_r_to_dm_dist, rank);
          dm_closest_to_ran[i].id_dist_candidate.push_back(i_r_to_dm_dist[rank]);       // keep this value to compare with those formn the incell search
          dm_closest_to_ran[i].global_id_candidate.push_back(n_index_tot[rank]);     // keep the global id of this dm candidate
       }
      So.DONE();
// Here find the minimum between the two sets in and out cell
#ifdef _FULL_VERBOSE_
      So.message_screen("Identifying id of closest DM particle");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       dm_index_closer_tot[i]=dm_closest_to_ran[i].sort_candidates();// get the global identification
 So.DONE();
 nearest_cells.clear();nearest_cells.shrink_to_fit();
 dm_closest_to_ran.clear();dm_closest_to_ran.shrink_to_fit();
#endif // end of SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
#ifdef _VERBOSE_FREEMEM_
     So.message_screen("Freeing memmory in line", __LINE__);
#endif
   dm_count.clear();
   dm_count.shrink_to_fit();
   ran_id.clear();ran_id.shrink_to_fit();
   So.DONE();
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Collapsing randoms using fraction of distance to closest dm particle = ", this->params._Distance_fraction());
#endif
#if defined _COLLAPSE_RANDOMS_USING_EXCLUSION_ || defined (_APPLY_GLOBAL_EXCLUSION_)
        // Estimate mean density of dark matter halos as rho_mean_background times Delta_vir
     real_prec density=this->Cosmo.mean_matter_density(this->params._redshift(), (void *)&this->s_cosmo_pars)*this->s_cosmo_info.Delta_Vir;
#endif
#ifdef _COLLAPSE_RANDOMS_VELS_
  const gsl_rng_type *  T;
  gsl_rng * rng ;
  gsl_rng_env_setup();
  ULONG seed=55457;//time(NULL);
  gsl_rng_default_seed=seed;
  T = gsl_rng_mt19937;//.gsl_rng_ranlux;
  rng = gsl_rng_alloc (T);
  real_prec Fvcoll=1.0;//pow(1.0-this->params._Distance_fraction(),1.0);   // cuanto varía la magnitud de la velocidad con el colapso:this is a model
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Scaling velocity |v| for randoms using a factor = ", Fvcoll);
#endif
#endif
#ifdef _USE_OMP_
//#pragma omp parallel for
#endif
     for(ULONG counter=0; counter< N_rand; ++counter)
      {
         //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Useful to retrieve coordinates
         // When closesd tm in only searched in the cell (incell), this gives the id of the closest dm tracer in the cell
         // Otherwise, this gives the global id of the closest dm particle
         ULONG index_dm_closer_a = dm_index_closer_tot[counter];
         // Read coordinates of random particles
         real_prec new_x = x_random_pos[counter];
         real_prec new_y = y_random_pos[counter];
         real_prec new_z = z_random_pos[counter];
#ifdef _SEARCH_CLOSEST_DM_IN_NEIGH_CELLS_
         real_prec xdm = this->tracer.Halo[index_dm_closer_a].coord1;
         real_prec ydm = this->tracer.Halo[index_dm_closer_a].coord2;
         real_prec zdm = this->tracer.Halo[index_dm_closer_a].coord3;
#else
         real_prec xdm = x_dm_pos[index_dm_closer_a];
         real_prec ydm = y_dm_pos[index_dm_closer_a];
         real_prec zdm = z_dm_pos[index_dm_closer_a];
#endif
        // redefine ran coords to the ref sistem of its closest dm particle:
         new_x -= xdm;
         new_y -= ydm;
         new_z -= zdm;
         real_prec mass_exclusion=1.0;
         real_prec dm_mass = mass_dm[index_dm_closer_a];              // Mass of dm particle
         real_prec ran_mass = mass_random[counter];                                // Mass of random
#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
         real_prec radius_dm = pow(3.0*dm_mass/(4.*M_PI*density), 1./3.);
         real_prec radius_ran = pow(3.0*ran_mass/(4.*M_PI*density), 1./3.);             // get an estimate of the radius of the
         real_prec RR=radius_dm+radius_ran;              // sum of the radii of the two halos
#endif
         // ARTBIRARY function to make fcol->fcol(M)
         real_prec reduced_mass=(dm_mass*ran_mass)/(dm_mass+ran_mass);
#ifdef _LOW_PROP_TO_RANDOMS_
         mass_exclusion=1.0-tanh(ran_mass/the(this->params._M_exclusion()));
#else
         mass_exclusion=1.0-tanh(reduced_mass/(this->params._M_exclusion()));
#endif
         // get the distance bewteen randoms and their closest DM particle:
         real_prec dist_random_to_dm = _get_modulo(new_x,new_y,new_z);
         //get the angular coordinates:
         real_prec theta = acos(new_z/dist_random_to_dm);
         real_prec phi   = atan2(new_y,new_x);
         real_prec new_sep = this->params._Distance_fraction()*dist_random_to_dm*mass_exclusion;
             // The separation between halos must be greater or equal than the sum of their radii
         real_prec sep_with_exclusion=new_sep;
#ifdef _COLLAPSE_RANDOMS_USING_EXCLUSION_
         sep_with_exclusion=max(new_sep,RR);
#endif
             // transfor to cartesian given a reduced distance and return ot origin of coords:
         this->tracer.Halo[ran_id_global[counter]].coord1=sep_with_exclusion*sin(theta)*cos(phi)+xdm;
         this->tracer.Halo[ran_id_global[counter]].coord2=sep_with_exclusion*sin(theta)*sin(phi)+ydm;
         this->tracer.Halo[ran_id_global[counter]].coord3=sep_with_exclusion*cos(theta)+zdm;
#ifdef _COLLAPSE_RANDOMS_VELS_
         // This section ortates the velocity-vector of the random particle towards
         // Original velocity components in the E system
         real_prec vx=x_random_v[counter];
         real_prec vy=y_random_v[counter];
         real_prec vz=z_random_v[counter];
         real_prec vdotr=vx*new_x+vy*new_y+vz*new_z;          //   Vec V cdot Delta r
         real_prec vv=_get_modulo(vx,vy,vz);   // |v|
         real_prec rxy=sqrt(pow(new_x,2)+pow(new_y,2));   //sqrt(dx²+dy²)
         // Coordinates of the velocity in the refernce frame E' (of the DM with hat y= hat r, r is the distance ran-dm closest)
         real_prec vxp=(new_y*vx-new_x*vy)/static_cast<double>(rxy); //vx'
         real_prec vyp=vdotr/dist_random_to_dm;// vy'=v cdot r/|r|
         real_prec vzp=(new_x*new_z*vx+new_y*new_z*vy)/(rxy*dist_random_to_dm)-(rxy/dist_random_to_dm)*vz; // vzp'
         real_prec vvp=_get_modulo(vxp,vyp,vzp); // modulo, just co cross-check with v
         real_prec alpha=acos(vyp/static_cast<double>(vv));//
         //* Get new velocity components in the E' system *//
         real_prec gamma_col=(1.0-this->params._Distance_fraction())*(alpha<M_PI ? (M_PI-alpha): (alpha -M_PI));// this is a model
//       nice idea: it does not help
/*         real_prec alpha_max=alpha<M_PI ? (M_PI-alpha): alpha;
         real_prec alpha_min=alpha<M_PI ? alpha: (M_PI-alpha);
         real_prec ran_angle= (alpha_min+0.5*(alpha_max-alpha_min))*(1+gsl_ran_gaussian(rng,M_PI/4.0)) ;//gsl_rng_uniform (rng);
         real_prec gamma_col=alpha_min+(alpha_max-alpha_min)*ran_angle; */
         real_prec uyp= alpha<M_PI ? Fvcoll*vv*cos(alpha+gamma_col) : Fvcoll*vv*cos(alpha-gamma_col);
         real_prec sigma=(Fvcoll*pow(vv,2)*cos(gamma_col)-vyp*uyp)/static_cast<double>(vxp);
         real_prec Gam=pow(Fvcoll*vv,2) - pow(uyp,2);
         double det= (dONE+ static_cast<double>((Gam/static_cast<double>(pow(sigma,2))-dONE)*(pow(vxp,2)+pow(vzp,2))/pow(static_cast<double>(vzp),2)));
         real_prec ddet=det <0 ? 0.0: det; // Si det es negativo, la raíz es imaginaria y nos sale una solución compleja, de la cual tomaremos sólo la parte real.
         real_prec uzp=(sigma*vzp*vxp/(pow(vxp,2)+pow(vzp,2))*(dONE+sqrt(ddet)));
         real_prec uxp=(pow(Fvcoll*vv,2)*cos(gamma_col)-vyp*uyp-vzp*uzp)/static_cast<double>(vxp);
         //* New velocity components in the E system *//
         real_prec ux=(new_y/rxy)*uxp+(new_x/dist_random_to_dm)*uyp+(new_x*new_z)*uzp/(rxy*dist_random_to_dm);
         real_prec uy=-(new_x/rxy)*uxp+(new_y/dist_random_to_dm)*uyp+(new_y*new_z)*uzp/(rxy*dist_random_to_dm);
         real_prec uz=(new_z/dist_random_to_dm)*uyp- rxy*uzp/dist_random_to_dm;
        // real_prec uup=_get_modulo(ux,uy,uz); // modulo, just co cross-check with v
         real_prec nux=std::isnan(ux) ? vx: ux;
         real_prec nuy=std::isnan(uy) ? vy: uy;
         real_prec nuz=std::isnan(uz) ? vz: uz;
         this->tracer.Halo[ran_id_global[counter]].vel1=nux;
         this->tracer.Halo[ran_id_global[counter]].vel2=nuy;
         this->tracer.Halo[ran_id_global[counter]].vel3=nuz;
#endif
       }
     this->So.DONE();
#ifdef _APPLY_GLOBAL_EXCLUSION_

#ifdef _FULL_VERBOSE_
     this->So.message_screen("Applying random shift to objects above ",this->params._M_exclusion() );
#endif
#ifdef _USE_OMP_
     vector<ULONG>vseeds(NTHREADS,0);
     for(int i=0;i<vseeds.size();++i)
         vseeds[i]=35+static_cast<ULONG>(i)*565;
      const gsl_rng_type *  rng_t;
      gsl_rng * gBaseRand;
      int jthread=0;
      gsl_rng_env_setup();
#pragma omp parallel private(jthread, rng_t,gBaseRand)
  {
      jthread=omp_get_thread_num();
      gsl_rng_default_seed=vseeds[jthread];
      rng_t = gsl_rng_mt19937;//_default;
      gBaseRand = gsl_rng_alloc (rng_t);
#pragma omp for
#endif
     for(ULONG i=0; i< this->tracer._NOBJS(); ++i)
      {
         if(this->tracer.Halo[i].mass>this->params._M_exclusion())
        {
          real_prec x = this->tracer.Halo[i].coord1;
          real_prec y = this->tracer.Halo[i].coord2;
          real_prec z = this->tracer.Halo[i].coord3;
          real_prec mass = this->tracer.Halo[i].mass;              // Mass of dm particle
          real_prec radius = pow(3.0*mass/(4.*M_PI*density), 1./3.);
//          this->tracer.Halo[i].coord1= x + gsl_ran_gaussian(gBaseRand, radius);
//          this->tracer.Halo[i].coord2= y + gsl_ran_gaussian(gBaseRand, radius);
//          this->tracer.Halo[i].coord3= z + gsl_ran_gaussian(gBaseRand, radius);
          this->tracer.Halo[i].coord1= x + radius*(gsl_rng_uniform(gBaseRand)-0.5);
          this->tracer.Halo[i].coord2= y + radius*(gsl_rng_uniform(gBaseRand)-0.5);
          this->tracer.Halo[i].coord3= z + radius*(gsl_rng_uniform(gBaseRand)-0.5);
         }
      }
#ifdef _USE_OMP_
}
#endif
     this->So.DONE();
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void BiasMT::execute()
  {
#ifdef _FULL_VERBOSE_
    std::cout<<endl;
#endif
    time_t time_bam;
    time(&time_bam);
#ifdef _ONLY_LPT_
    So.message_BiasMT(time_bam);
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
    So.message_screen("Building halo catalogs using Bias Mapping Technique.");
    So.message_screen("Realization =",this->params._IC_index());
#endif
    this->So.initial_time=time_bam;
    So.enter(__PRETTY_FUNCTION__);
    int NTHREADS=_NTHREADS_;   // omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
    this->warnings();
#ifdef MOCK_MODE
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _FULL_VERBOSE_
     So.message_screen("GENERATING NEW TRACER DENSITY FIELD FROM BIAS AND KERNEL");
     std::cout<<endl;
#endif
#else
#ifdef _FULL_VERBOSE_
     std::cout<<endl;
     So.message_screen("CALIBRATING BIAS AND KERNEL FROM REFERNECE SIMULATION");
     std::cout<<endl;
#endif
#endif
#else
       So.message_screen("Statistics of halo bias");
#endif
     // Initialize cosmological functions using the input redshift
     this->tstruct=0;
     this->im=0; //mas
     // Get the cosmological quantities derived from input cosmological parameters.
     this->get_cosmo();
     this->tracer_ref.set_type_of_object("TRACER_REF");
     this->tracer_ref.set_params(this->params);
     this->cwclass_ref.set_params(this->params);
     this->tracer.set_type_of_object("TRACER_MOCK");
     this->tracer.set_params(this->params);
     this->cwclass.set_params(this->params);
     // Define strings for input DF
     string file_X, file_X_ref_pdf, file_Y,file_Y_HR;
     string file_Vx, file_Vy, file_Vz;
     // IF DEFINED, READ INPUT FILE WITH INFO FROM THE TRACER
     // If READ_REF_CAT is undef or GET_BiasMT_RALIZATIONS is def , the code will read these fields and analyze then
#ifdef _READ_REF_CATALOG_
     string file_density_field_tracer=this->params._Output_directory()+"TR_DENS_FIELD";
     string file_Y_mass=this->params._Output_directory()+"TR_MASS_DENS_FIELD";
     string file_Y_sat_frac=this->params._Output_directory()+"TR_SAT_FRACTION_FIELD";
#else
     string file_density_field_tracer=this->params._Output_directory()+"TR_DENS_FIELD.dat";
     string file_Y_mass=this->params._Output_directory()+"TR_MASS_DENS_FIELD.dat";
     string file_Y_sat_frac=this->params._Output_directory()+"TR_SAT_FRACTION_FIELD.dat";
#endif
#ifdef _DO_BiasMT_CALIBRATION_
#ifdef _READ_REF_CATALOG_
     // Here we shall read the reference and produce two density fields,
     // one with number counts and other with mass weighted number counts
     // The one with number counts will be written in the file pointed to
     // by the parameter file_Y, such that it can be read later below
     this->So.message_screen("Reading the Reference Catalog of tracers");
#ifdef _SET_CAT_WITH_MASS_CUT_
     this->tracer.read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),pow(10,params._LOGMASSmin())*params._MASS_units());
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
     this->tracer.read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),params._VMAXmin());
#endif
     this->tracer.get_density_field_grid(_COUNTS_, file_density_field_tracer);
#ifdef _USE_MASS_FIELD_
     this->tracer.get_density_field_grid(_MASS_, file_Y_mass);
     string fname_mass_function_Y = this->params._Output_directory()+"tracer_ref_abundance.txt";
     this->tracer.get_property_function(fname_mass_function_Y);
#endif
#endif
#endif
#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
     this->tracer.read_catalog_bin( this->params._NGRID(),"A","A","A");
     exit(0);
#endif
     // HERE WE CAN USE GET_NW_DM() ALSO
#if defined(_USE_LPT_) || defined (_GET_BiasMT_CAT_)
     const gsl_rng_type *rng_t;
#ifdef _USE_OMP_
     gsl_rng **gBaseRand;
#else
    const gsl_rng_type *  T;
     gsl_rng *gBaseRand;
     gsl_rng_default_seed=seed;
     T = gsl_rng_mt19937;//.gsl_rng_ranlux;
     gBaseRand = gsl_rng_alloc (T);
#endif
     gsl_rng_env_setup();
     rng_t = gsl_rng_mt19937;// gsl_rng_default;
#ifdef _USE_OMP_
     int nt=omp_get_max_threads();
     gBaseRand = (gsl_rng **) malloc(nt * sizeof(gsl_rng *));
#pragma omp parallel for num_threads(nt)
     for(int i=0;i<nt;i++)
       {
         gBaseRand[i] = gsl_rng_alloc(rng_t);
         gsl_rng_set(gBaseRand[i],this->params._seed());
       }
#endif
     // This is used for the calibration *and* the construction of mocks
#ifdef _USE_LPT_
     this->lpt.set_s_cosmo_info(this->s_cosmo_info);
     //share the s_cosmo_info with LPT
     this->lpt.set_params(this->params);
     //LPT reads its parameters from params class
     this->lpt.set_fnames();
#else
     So.message_warning("Check names in line", __LINE__);
#endif    // end use LPT and get_BiasMT_cat line 8905


     // set the output name files for LPTt
#endif    // end use LPT and get_BiasMT_cat line 8905
     // If we need to do the calibration,
#ifdef _DO_BiasMT_CALIBRATION_
     //  we might want to use LPT either to creat the DM, or assign position to partices in mock, or both
#if defined(_USE_LPT_) || defined (_GET_BiasMT_CAT_)
     time_t start_LPT;
     time(&start_LPT);
#endif // end use LPT and get_BiasMT_cat line 8942
#ifdef _USE_LPT_
#ifdef _USE_OMP_
     this->lpt.get_dm_field(gBaseRand);
#else
     this->lpt.get_dm_field(gBaseRand);  // Run LPT!
#endif
#ifdef _FULL_VERBOSE_
     this->So.message_screen("LPT has created DMDF using ALPT");
     So.message_time2(start_LPT);
     std::cout<<endl;
#endif  // end use LPT line 8950
#endif  // end use LPT line 8950
#if defined (_USE_LPT_) & !defined(_READ_VELOCITIES_)
     // Read the file names generated in LPT
     file_X =this->lpt._fnameDM()+".dat";
     file_Vx=this->lpt._fnameVX()+".dat";
     file_Vy=this->lpt._fnameVY()+".dat";
     file_Vz=this->lpt._fnameVZ()+".dat";
#endif
#ifdef _FULL_VERBOSE_
     this->So.message_screen("BAM:");
#endif
#ifdef _READ_VELOCITIES_
     file_Vx=this->params._Input_Directory_X()+this->params._Name_VelFieldx_X();
     file_Vy=this->params._Input_Directory_X()+this->params._Name_VelFieldy_X();
     file_Vz=this->params._Input_Directory_X()+this->params._Name_VelFieldz_X();
#endif
#endif //ifdef do_BiasMT_calibration
#ifdef _USE_LPT_
     file_X=this->lpt._fnameDM()+".dat";
#else
     file_X=this->params._Input_Directory_X()+this->params._Name_Catalog_X();
#endif
     // Input file coptaining the reference tracer
#ifndef _READ_REF_CATALOG_
     file_Y_HR=this->params._Input_Directory_Y()+this->params._Name_Catalog_Y_HR();
     file_Y=this->params._Input_Directory_Y()+this->params._Name_Catalog_Y();
#else
     file_Y = file_density_field_tracer+".dat";
     file_Y_HR = file_density_field_tracer+".dat"; //in thie case we pass the same file
     file_Y_mass += ".dat";
     file_Y_sat_frac += ".dat";
#endif
    this->iteration=0;
     // **************************************************************************************************
#ifdef _DO_BiasMT_CALIBRATION_
#ifdef _USE_VELOCITIES_
     this->read_BiasMT_files(file_X, file_Y, file_Y_HR, file_Y_mass, file_Y_sat_frac, file_X_ref_pdf, file_Vx, file_Vy, file_Vz);
#endif // end ifdef _USE_VELOCITIES_
#endif
#ifdef _ONLY_LPT_  // This region is devoted to mesure the power of the DM generated in LPT when only LPT is used.
         this->delta_X_ini.resize( this->params._NGRID(),0);
         this->File.read_array(file_X,delta_X_ini);
         this->get_power_spectrum("DM_REF");
         this->delta_X_ini.clear(); delta_X.shrink_to_fit();
         exit(0);
#endif
#if defined (_DO_BiasMT_CALIBRATION_) || defined (BIAS_MODE)
#ifndef _test_mass_assign_
         {
#ifndef _USE_VELOCITIES_
           this->read_BiasMT_files(file_X, file_Y, file_Y_HR, file_Y_mass, file_Y_sat_frac, file_X_ref_pdf);
           // Read all input files (density fields in binary )
#endif
         }
#endif
#endif
    // Intialize arrays for power spectrum of input fields
      this->set_Fourier_vectors();
#if defined (_DO_BiasMT_CALIBRATION_)  || defined (BIAS_MODE)
         this->analyze_input(); // analyze the input references, requested also to create mock for limits
#endif
         // Get numbers from the input density fields and their power spectrum.
         // In this function, if requested from parameter file, an iterative process is performed in order to make the approx method
         // DM field match the reference DM power (if N_iterations_dm >0)
         // *******************************************************************************************************
         // load parameters for the cwclass
#if defined (_USE_CWEB_) || defined (_USE_TWEB_) || defined (_USE_IWEB_) || defined (_USE_IKWEB_) || defined (_USE_AWEB_) || defined (_USE_IWEB_V_) || defined (_USE_PWEB_) || defined (_USE_NABLA2DELTA_) || defined (_USE_S2DELTA_) || defined (_USE_S3_)  || defined (_USE_S2_) ||  defined (_USE_DELTA3_) || defined (_USE_TIWEB_)
         this->cwclass.s_cosmo_pars=this->s_cosmo_pars;
#endif
         //  Clean and initialize arrays for power spectrum  and kernels
         this->set_Fourier_vectors();
#ifdef BIAS_MODE
#ifdef _SEVERAL_REAL_BIAS_
         this->get_pdf(id);
#elif !defined _SEVERAL_REAL_BIAS_
         this->get_pdf();
#endif // end _SEVERAL_REAL_
#endif // end ifdef BIAS
#ifdef _TEST_THRESHOLDS_RESIDUALS_
//         this->file_residuals=this->params._Output_directory()+"Resuduals_thresholds_original_vel.txt";
         this->file_residuals=this->params._Output_directory()+"chis_thresholds_original_vel.txt";
         this->output_res.open(this->file_residuals.c_str());
         for(int ilt=0; ilt<100;ilt++)
           for(int ilv=0; ilv<100;ilv++)
             {
               std::cout<<endl;
               this->lambdath_v=-1.0+2.0*(static_cast<real_prec>(ilv))/static_cast<real_prec>(100.0);
               this->params.set_lambdath ((1.0)*(static_cast<real_prec>(ilt))/static_cast<real_prec>(100.0));
               this->cwclass.lambdath=this->params._lambdath();
               this->cwclass.lambdath_v=this->lambdath_v;
#endif
               // DO the V-classification from the velocity field.
               // This is done once, since the Vel field won't change during the iterative process.
               // Otherwise this should have been done within the iterations, in the get_BiasMT_DM method.
#if defined (_USE_CWEB_V_) || defined (_USE_IWEB_V_)
               this->cwclass.get_CWC_V(this->Velx_X, this->Vely_X,this->Velz_X);

#ifdef _USE_VEL_KNOTS_V_
               this->cwclass.get_SigmaVel_collapsing_regions(this->delta_X_ini,  this->Velx_X, this->Vely_X,this->Velz_X, static_cast<real_prec>(this->N_objects_X)/static_cast<real_prec>( this->params._NGRID()));
#endif
#endif
               // ****************************ITERATIVE BAM AND MOCK GENERATION******************************************
#ifdef MOCK_MODE
               // Get the total number of steps adding those of the calibration of the kernal and the number of
               // DM realizations used to get the same number of mock density fields
               int N_steps=this->params._N_iterations_Kernel()+1;
               int init=0;
#ifdef _TEST_THRESHOLDS_RESIDUALS_
               init=0;
               this->params._set_N_iterations_Kernel(0);
#endif
               // If we read the Kernel, we start directly creating the mocks
                this->use_iteration_ini=false;
#ifdef _GET_BiasMT_REALIZATIONS_
               init = this->params._N_iterations_Kernel()+1;
               this->get_min_max_X_Y();
               // in get bias we read the bias and the kernel ,only when we want to get the realizations
               this->get_BIAS(this->params._Name_Property_Y());
#endif
               // i=0, zero order approach
               // from i=1, to i=N_iterations_Kernel, we calibrate the Kernel.
               // From i=N_iterations_Kernel+1 to i=N_steps, we create mocks based on independent realization of approx DM fields
#ifdef _FULL_VERBOSE_
               time_t start_end;
               time_t start_aux;
               time(&start_end);
               start_aux=start_end;
#endif
#ifdef _DO_BiasMT_CALIBRATION_
               // Loop over the number if iterations demanded to perform the calibration of the Kernel and the Bias
#ifdef _USE_GNUPLOT_
              this->gp_kernel<<"set size 1.0, 1.0\n";
              this->gp_kernel<<"set origin 0.0, 0.0\n";
              this->gp_kernel<<"set style function lines\n";
#endif
#ifdef _DISPLACEMENTS_
              this->lpt.Displacement.resize( this->params._NGRID(),0);
              this->Displacement_inicial.resize( this->params._NGRID(),0);
//              File.write_array(this->lpt.fnameICDELTA, this->delta_X_ini); // Lag2Eul reads thois file
              this->lpt.get_displacement(gBaseRand,this->delta_X_ini, 0);// delta_X acts as IC. This function generates the container vector<real_prec>LPT.displacement
              this->Displacement_inicial=this->lpt.Displacement;
              this->File.read_array(this->lpt.fnameICDELTA+".dat", delta_X_ini); // raead
#endif
              for(int i=init; i<=this->params._N_iterations_Kernel() ;++i)
               {
#ifdef _USE_GNUPLOT_
                this->gp_kernel<<"set multiplot \n";
                this->gp_kernel<<"set border linewidth 2.0\n";
                this->gp_kernel<<"set border linecolor '"<<FG_COLOR<<"' \n";
                this->gp_kernel<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
                this->gp_kernel<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
                this->gp_kernel<<"set xrange[*:*] \n";
                this->gp_kernel<<"set yrange[*:*] \n";
#endif
#ifdef _FULL_VERBOSE_
                   time_t start;
                   time(&start);
                   std::cout<<endl;
                   std::cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
                   std::cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
#endif
                   this->iteration=i;
#if defined _USE_CWC_ || defined (_USE_INVARIANT_TIDAL_FIELD_II_) || defined (_USE_INVARIANT_TIDAL_FIELD_III_)
                   this->cwclass.step=i;
#endif
                   // Calculations performed during the iterative process.
#ifdef _FULL_VERBOSE_
                   if(i==0)
                     So.message_screen("BAM raw mapping" , i, start, start_aux);
                   else
                     So.message_screen("Iteration ", i, start,start_aux);
                   std::cout<<CYAN<<"**************************************************************************"<<RESET<<endl;
#endif
#ifdef _DISPLACEMENTS_
                   if(i>0)
                     {
                       this->GetKernel(true,1);// Get the kernel from the ratio of P_disp_new/P_disp_ref
                       this->Konvolve(this->delta_X_ini,this->delta_X);
                     }
                   File.write_array(this->lpt.fnameICDELTA, this->delta_X); // Lag2Eul reads this file
                   this->lpt.get_displacement(gBaseRand,this->delta_X, i);// delta_X acts as IC. This function generates the container vector<real_prec>LPT.displacement
                   for(ULONG i=0;i<  this->params._NGRID();++i)
                     this->delta_Y[i]-=this->Displacement_inicial[i]; //Delta=Psi_TRUE - Psi_APPROX
#ifdef _USE_CWC_
                   this->cwclass.get_CWC(this->delta_X);
#endif
#else
                   // step i) Do the CWC if requested
                   // Get the ratio T, update kernel K, convolve K with original DM field.
                   // The kernel is computed with the Power aspectrum or the mass weighted power spectrum, according to the preproc def _USE_MASS_WEIGHTED_KERNEL_
                   this->get_BiasMT_DM();
#endif
                   // Step ii)
                   // Compute the halo bias from Numnber counts reference and DM from step i)
                   this->get_BIAS(this->params._Name_Property_Y());
                   this->get_mock_grid(this->params._Name_Property_Y());
                   // Step iii)
                   // Generate the Halo number density field by sampling the DM from step i) using the information of the bias from step ii)
                   // argument false indicates that the DM used in the one of the reference (either original Nbody or approximated)
                   // Also gets the power spectrum of the mock
                   // Since the info of the mass distribution and the sat fraction is not used for the clustering analysis to calibrate the bias,
                   // we can compute them in the last step of the iteration, with the DM already transformed with the Bam kernel.
#ifdef _FULL_VERBOSE_
                   start_aux=start;
#endif
#ifdef _USE_GNUPLOT_
                   this->gp_kernel<<"unset multiplot \n";
#endif
                   } // closes iteration
#endif
#ifdef _TEST_THRESHOLDS_RESIDUALS_
               real_prec residuals=0;
               int ncounts=0;
               /*
#pragma omp parallel for reduction(+:residuals,ncounts)
               for(int i=0;i<this->Power_REF.size();++i)
                if(this->Power_NEW[i]>0)
                 {
                   ncounts++;
                  residuals+=fabs(this->Power_REF[i]/this->Power_NEW[i]-1.0);
                 }
               residuals/=static_cast<real_prec>(ncounts);
               So.message_screen("Residuals at this iteration (%) =",100.0*residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<100.0*residuals<<endl;
             }
               So.message_screen("Residuals at this iteration (%) =",100*residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<100.0*residuals<<endl;
 */
               real_prec deltak=this->kvec[1]-this->kvec[0];
               for(int i=0;i<this->Power_REF.size();++i)
                if(this->Power_NEW[i]>0)
                 {
                   ncounts++;
                   real_prec sigma_squared = (1./(2.*M_PI*deltak*pow(this->kvec[i],2)*pow(this->params._Lbox(),3))) *(this->Power_NEW[i]*this->Power_NEW[i]); // ESTE PARA EL KERNEL DE LA DM
                   residuals+=pow(this->Power_REF[i]-this->Power_NEW[i],2)/(sigma_squared);
                 }
               residuals/=static_cast<real_prec>(ncounts);
               /* not rady yet:
               if(this->cwclass_ref.volume_knots>0)
                 residuals+=pow(this->cwclass.volume_knots-this->cwclass_ref.volume_knots,2)/sqrt(this->cwclass_ref.volume_knots);
               if(this->cwclass_ref.volume_filaments>0)
                 residuals+=pow(this->cwclass.volume_filaments-this->cwclass_ref.volume_filaments,2)/sqrt(this->cwclass_ref.volume_filaments);
               if(this->cwclass_ref.volume_sheets>0)
                  residuals+=pow(this->cwclass.volume_sheets-this->cwclass_ref.volume_sheets,2)/sqrt(this->cwclass_ref.volume_sheets);
               if(this->cwclass_ref.volume_voids>0)
                 residuals+=pow(this->cwclass.volume_voids-this->cwclass_ref.volume_voids,2)/sqrt(this->cwclass_ref.volume_voids);
                */
               So.message_screen("chi² =",residuals);
               this->output_res<<this->lambdath<<"\t"<<this->lambdath_v<<"\t"<<residuals<<endl;
           }
         this->output_res.close();
#endif
#if defined(_SEVERAL_REAL_BIAS_) || defined(_SEVERAL_REAL_CAL_)
       }
#endif
#ifdef _GET_BiasMT_REALIZATIONS_
#ifdef _FULL_VERBOSE_
     time(&start_end);
     start_aux=start_end;
#endif
     // Sample the bias measured from the reference into other realzations of the Approximated density field.
     // Loop over the new DM fields
     int i=init;
#ifndef _ONLY_POST_PROC_
     for(i=init; i<=N_steps ;++i)
       {
#endif
#ifdef _FULL_VERBOSE_
         time_t start;
         time(&start);
#endif
         this->iteration=init;
         int i_realization=i-this->params._N_iterations_Kernel()+this->params._N_dm_initial()-1;
         this->iteration=i;
#ifndef _ONLY_POST_PROC_
#ifdef _FULL_VERBOSE_
         So.message_screen("Creating kmock");
#endif
#else
#ifdef _FULL_VERBOSE_
         So.message_screen("Assigning properties to mock");
#endif
#endif
         // Get the new dm Df from ALPT, convolve with the kernel and  and get its properties ({theta}). This will be used for the conditional assignment of properties.
#ifndef _FCOL_TESTS_  // For the test checking different values of fcol, we do not need the dark matter denisty field
         this->get_new_DM_field();
#endif
#ifndef _ONLY_POST_PROC_   // These Lines are used only if the number counts of mocks are to be generated. If not, (i.e, if already generated), this is ognored
#ifdef _USE_TWO_REFS_MOCKS_
     this->get_mock_grid_two(_COUNTS_);
#else
     this->get_mock_grid(_COUNTS_);
#endif
#endif
         // I) Using the information of the bias and kernel, apply these two to a new DM field
         // in order to generate a new halo number denisty field. Measure Power spectrum pof the mock.
    // This sections names the halo number count on which the assignment will be done in makecat()
#ifdef _ASSIGN_TO_CALIBRATION_
         this->fnameMOCK=this->params._Output_directory()+"MOCK_TR_iteration400_MASY"+to_string(this->params._iMAS_Y())+"_Nft"+to_string(this->params._Nft());
#elif defined _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_
         this->fnameMOCK=this->params._Input_dir_cat_new_ref()+this->params._Name_Catalog_Y_new_ref();
#elif defined _ASSIGN_PROPERTIES_TO_REFERENCE_ && !defined (_ASSIGN_PROPERTIES_TO_NEW_REFERENCE_)
         this->fnameMOCK=this->params._Input_Directory_Y()+this->params._Name_Catalog_Y();
#elif defined (_ASSIGN_PROPERTIES_TO_MOCK_)
	 //         this->fnameMOCK=this->params._Output_directory()+"MOCK_TR_Realization"+to_string(this->params._realization())+"_p"+to_string(this->params._unitsim_plabel())+"_TWEB_Nft"+to_string(this->params._Nft());
	 this->fnameMOCK=this->params._Output_directory()+"MOCK_TR_Realization"+to_string(this->params._realization())+"_MASY0_Nft"+to_string(this->params._Nft());
#endif
#ifdef _GET_BiasMT_CAT_
         this->BiasMT_makecat(gBaseRand);
#endif
#ifdef _FULL_VERBOSE_
         start_aux=start;
#endif
#ifndef _ONLY_POST_PROC_
       }
#endif
#endif //end of GET_BiasMT_REALIZATIONS
#endif //end of MOCK_MODE
     this->So.message_time(time_bam);
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _GET_BiasMT_CAT_
#ifdef _USE_OMP_
void BiasMT::BiasMT_makecat(gsl_rng ** gBaseRand)
#else
void BiasMT::BiasMT_makecat(gsl_rng * gBaseRand)
#endif
{
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_LPT_
    string stradd=this->lpt._stradd_bam();
#else
    So.message_warning("CHeck names in line ", __LINE__);
#endif
    int ir=this->params._IC_index();
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
    // In this function we asign coordinates from DM + random particles based on the mock number counts
     // Masses are also assigned and the collapse of randoms towards dm (to correct small scale clusterin) is also performed after that
     // Catalog is then written
#ifdef _FULL_VERBOSE_
     std::cout<<endl;
     So.message_screen("***********************************************************************");
     So.message_screen("********************GENERATING MOCK CATALOG****************************");
     So.message_screen("***********************************************************************");
     std::cout<<endl;
#endif
#if defined (_ASSIGN_PROPERTIES_TO_MOCK_) || defined (_ASSIGN_TO_CALIBRATION_)
#ifdef _FULL_VERBOSE_
     std::cout<<endl;
     this->So.message_screen("*****Assigning position and velocities to tracers using DM particles****");
     std::cout<<endl;
#endif
#ifdef _FCOL_TESTS_
#ifdef _FULL_VERBOSE_
     std::cout<<endl;
     this->So.message_screen("FCOL test is active. No velocities will be loaded.");
     std::cout<<endl;
#endif
#endif
     real_prec redshift=this->s_cosmo_pars.cosmological_redshift;
     s_params_box_mas box;
     box.min1=this->params._xllc();
     box.min2=this->params._yllc();
     box.min3=this->params._zllc();
     box.Lbox=this->params._Lbox();
     box.Nft=this->params._Nft();
     box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
     box.d2= box.d1;
     box.d3= box.d1;
     box.NGRID= this->params._NGRID();
     ULONG N_dms=0;
     ULONG N_rand=0;
#ifdef _COLLAPSE_RANDOMS_AUX_
     s_params_box_mas box_collps;
     box_collps.min1=this->params._xllc();
     box_collps.min2=this->params._yllc();
     box_collps.min3=this->params._zllc();
     box_collps.Lbox=this->params._Lbox();
     box_collps.Nft=this->params._Nft_random_collapse();
     box_collps.d1= box_collps.Lbox/static_cast<real_prec>(box_collps.Nft);		/* grid spacing x-direction */
     box_collps.d2= box_collps.d1;
     box_collps.d3= box_collps.d1;
     box_collps.NGRID=(box_collps.Nft*box_collps.Nft*box_collps.Nft);
     vector<int>dm_count(box_collps.NGRID,0);
     vector<real_prec> x_dm_pos;
     vector<real_prec> y_dm_pos;
     vector<real_prec> z_dm_pos;
     vector<ULONG> dm_id;
     vector<real_prec> x_random_pos;
     vector<real_prec> y_random_pos;
     vector<real_prec> z_random_pos;
     vector<ULONG> ran_id;
#ifdef _COLLAPSE_RANDOMS_VELS_
     vector<real_prec> x_random_v;
     vector<real_prec> y_random_v;
     vector<real_prec> z_random_v;
#endif
#endif  //* END FOR _COLLAPSE_RANDOMS_AUX_
     real_prec d1=this->params._d1();
     real_prec d2=this->params._d2();
     real_prec d3=this->params._d3();
     // ****************************DM POSITIONS *******************************************
     ULONG N_dm= this->params._NGRID();     // Note that here the number of dm particles is that of the Ngrid
     /* Even if we do not use the dm particles coords, we need to read them in order to get the velocity interpolated. But this is
     just because the interpolation is done here, while it could have been made in LPT*/
     vector<real_prec> posx(N_dm,0),posy(N_dm,0),posz(N_dm,0);
     vector<ULONG>index(N_dm,0);
     // Filenames of the bninary files containing the positions of the dark matter particles
     this->File.read_array(this->lpt._fnamePOSX()+".dat",posx);
     this->File.read_array(this->lpt._fnamePOSY()+".dat",posy);
     this->File.read_array(this->lpt._fnamePOSZ()+".dat",posz);
    // I have lost track of the velocities interpolated in the mesh. The trilinear interpolation algorithm
     // expects that the trining set is part of a regular mesh. Hence, if we ise the line3D function, this has the be
     // provided the intepolated v field. We are providing the actual velocities
#ifdef _GET_BiasMT_REALIZATIONS_
     this->lpt.set_fnames();
#endif
#ifndef _FCOL_TESTS_
    //***********************************************************************************************
    //  Cointainers for the velocities of the dm particles evaluated at Lagrangian positions.
     // THese must be somehow used to assign velocities at a Eulerian position
     // Either through intepolation or by searching the closest lagrangian position.
     // This section generates the intepolated vel field which is used to apply the trilinear interpolation
     // HOwever this does not wors as well as the implementation below
     vector<real_prec> dens_field_ngp(N_dm,0);
     getDensity_NGP(box.Nft,box.Nft,box.Nft,this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),box.d1,box.d2,box.d3,box.min1,box.min2,box.min3,posx,posy,posz,posz,dens_field_ngp,false);
     vector<real_prec> vel_aux(N_dm,0);
     this->File.read_array(this->lpt._fnameVXpart()+".dat",vel_aux);
     vector<real_prec> velx(N_dm,0);
     vector<real_prec> vely;
     vector<real_prec> velz;
#ifdef _FULL_VERBOSE_
     So.message_screen("Interpolating Vx");
#endif
     getDensity_NGP(box.Nft,box.Nft,box.Nft,this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),box.d1,box.d2,box.d3,box.min1,box.min2,box.min3,posx,posy,posz,vel_aux,velx,true);
     So.DONE();
#ifdef _USE_NP_VEL_
     vector<ULONG>id_empty_cells;  // COntainer to allocate the index ID of empty cells
#endif
#ifdef _FULL_VERBOSE_
     So.message_screen("Averaging and identifying empty cells");
#endif
     for(ULONG i=0;i<N_dm;i++)
       if(dens_field_ngp[i]>0)
          velx[i]/=static_cast<real_prec>(dens_field_ngp[i]);
#ifdef _USE_NP_VEL_
      else
        id_empty_cells.push_back(i);   // THis is needed only once
#endif#ifdef _USE_LPT_
#ifdef _USE_NP_VEL_
     So.message_screen("Fraction of empty cells =",100.0*id_empty_cells.size()/static_cast<real_prec>(N_dm),"%");
     // Even if these are empty cells (devoid of dm tracers), it might be that we need halos there so we need the correct vel field.
     vector<s_nearest_cells>nearest_cells;
#endif
#ifdef _USE_NP_VEL_
    if(id_empty_cells.size()>0)
    {
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying Neighbouring Cells (this is done once)");
#endif
     nearest_cells.resize(id_empty_cells.size());
     // Obtain the ID of the nearst cells of those empty cells.
     // Nearest cells come ordered from 1 to NUmber of empty cells
     get_neighbour_cells_of_cell(this->params._Nft(), id_empty_cells,1, nearest_cells);
     So.DONE();
#ifdef _FULL_VERBOSE_
     So.message_screen("Assigning Vx to empty cells");
#endif
     get_vel_field_using_neigh_cells(id_empty_cells,nearest_cells,velx);
    }
#endif
    // Read Comp Y of velocities
    vel_aux.resize(N_dm,0);
    this->File.read_array(this->lpt._fnameVYpart()+".dat",vel_aux);
    vely.resize(N_dm,0);
    getDensity_NGP(box.Nft,box.Nft,box.Nft,this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),box.d1,box.d2,box.d3,box.min1,box.min2,box.min3,posx,posy,posz,vel_aux,vely,true);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<N_dm;i++)
      if(dens_field_ngp[i]>0)
        vely[i]/=static_cast<real_prec>(dens_field_ngp[i]);
#ifdef _USE_NP_VEL_
    if(id_empty_cells.size()>0)
     {
#ifdef _FULL_VERBOSE_
       So.message_screen("Assigning Vy to empty cells");
#endif
      get_vel_field_using_neigh_cells(id_empty_cells,nearest_cells,vely);
    }
#endif
    // Read Comp Z of velocities
    vel_aux.resize(N_dm,0);
    this->File.read_array(this->lpt._fnameVZpart()+".dat",vel_aux);
    velz.resize(N_dm,0);
    getDensity_NGP(box.Nft,box.Nft,box.Nft,this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),box.d1,box.d2,box.d3,box.min1,box.min2,box.min3,posx,posy,posz,vel_aux,velz,true);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<N_dm;i++)
          if(dens_field_ngp[i]>0)
            velz[i]/=static_cast<real_prec>(dens_field_ngp[i]);
#ifdef _USE_NP_VEL_
    if(id_empty_cells.size()>0)
    {
#ifdef _FULL_VERBOSE_
     So.message_screen("Assigning Vz to empty cells");
#endif
      get_vel_field_using_neigh_cells(id_empty_cells,nearest_cells,velz);
    }// closes if(id_empty_cells.size()>0)
    id_empty_cells.clear(); id_empty_cells.shrink_to_fit();
#endif
    dens_field_ngp.clear();dens_field_ngp.shrink_to_fit();
    vel_aux.clear();vel_aux.shrink_to_fit();
/*
    // My implementation was: read particle Vels in Lag and interpolate to eulerian
    // This should work better as the positions at which these are interpolated
    // ar ethe Lagrangian, which are on a regular grid.
    vector<real_prec> velx(N_dm,0);
    this->File.read_array(this->lpt.fnameVXpart+".dat",velx);
     vector<real_prec> vely(N_dm,0);
    this->File.read_array(this->lpt.fnameVYpart+".dat",vely);
     vector<real_prec> velz(N_dm,0);
    this->File.read_array(this->lpt.fnameVZpart+".dat",velz);
*/
// COnvolve the velocities with a kernel of the form 1/(1+s² k²). This has to be applin to the intepolated vel-field.
    if(true==this->params._use_vel_kernel())
    {
        bool apply_kernel_to_knots=false;  // defined temporally here.
        if(apply_kernel_to_knots==true)
        {
      vector<real_prec>aux_vx(N_dm,0);
      vector<real_prec>aux_vy(N_dm,0);
      vector<real_prec>aux_vz(N_dm,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<aux_vx.size();++i)
       {
       int I_CWT=this->cwclass.get_Tclassification(i);
       if(this->cwclass.cwt_used[I_CWT]== I_KNOT ||  this->cwclass.cwt_used[I_CWT]== I_FILAMENT)
         {
           aux_vx[i]=velx[i];
           aux_vy[i]=vely[i];
           aux_vz[i]=velz[i];
         }
       }
      this->v_konvolve(aux_vx);
      this->v_konvolve(aux_vy);
      this->v_konvolve(aux_vz);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<aux_vx.size();++i)
       {
       int I_CWT=this->cwclass.get_Tclassification(i);
       if(this->cwclass.cwt_used[I_CWT]== I_KNOT  ||  this->cwclass.cwt_used[I_CWT]== I_FILAMENT)
         {
           velx[i]=aux_vx[i];
           vely[i]=aux_vy[i];
           velz[i]=aux_vz[i];
         }
       }
        aux_vx.clear();aux_vx.shrink_to_fit();
        aux_vy.clear();aux_vy.shrink_to_fit();
        aux_vz.clear();aux_vz.shrink_to_fit();
   }
else{
      this->v_konvolve(velx);
      this->v_konvolve(vely);
      this->v_konvolve(velz);
      }
    }
#endif // endif for #ifndef _FCOL_TESTS_
     real_prec vel_bias_dm=1.+this->params._velbias_dm();// This value is put here to make the vel distribution of the BAM-SLIS mocks be similat to the reference
     real_prec vel_bias_ran=1.+this->params._velbias_random();
     struct s_cell_info{
       vector<real_prec> posx_p;
       vector<real_prec> posy_p;
       vector<real_prec> posz_p;
#ifndef _FCOL_TESTS_
       vector<real_prec> velx_p;
       vector<real_prec> vely_p;
       vector<real_prec> velz_p;
       vector<real_prec> mass_p;
#endif // endif for #ifndef _FCOL_TESTS_
       vector<ULONG> id_p;
     };
#endif   // end of ifdef _ASSIGN_TO_MOCK OR assign_to_calibration
// We make a break here in the ifndef assign to reference in order to open the reference number count file
// When assining to reference or new reference, we di not need the number counts unless we explicitely 
// call the _USE_TRACERS_IN_CELLS_
     vector<real_prec> MOCK_DEN_FIELD;
#if defined (_ASSIGN_PROPERTIES_TO_MOCK_) || defined (_ASSIGN_TO_CALIBRATION_) || defined _USE_TRACERS_IN_CELLS_
     MOCK_DEN_FIELD.resize(this->params._NGRID(),0);
     this->File.read_array(this->fnameMOCK,MOCK_DEN_FIELD);
     ULONG Nobjects_mock=get_nobjects(MOCK_DEN_FIELD);
     ULONG Ntracers = Nobjects_mock;  //reduntant, but still
     this->tracer.set_NOBJS(Nobjects_mock);
#endif

#if defined (_ASSIGN_PROPERTIES_TO_MOCK_) || defined (_ASSIGN_TO_CALIBRATION_)
     // The following section assign coords and vels to mocks. If a test is done over references, we do not need this
     // procedure.
     vector<int> aux_cont( this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel num_threads(4)
     {
         vector<int> aux_cont_par( this->params._NGRID(),0);
#pragma omp for nowait
#endif
     for(ULONG i=0;i<N_dm;++i)
       {
         real_prec x = static_cast<real_prec>(posx[i]);
         real_prec y = static_cast<real_prec>(posy[i]);
         real_prec z = static_cast<real_prec>(posz[i]);
         ULONG ind=grid_ID(&box, x,y,z);
         index[i]=ind;
#ifdef _USE_OMP_
         aux_cont_par[ind]++;  //count the number of dm particle in a cell
#endif
       }
#ifdef _USE_OMP_
#pragma omp critical
     for(ULONG i=0;i< this->params._NGRID();++i)
         aux_cont[i]+=aux_cont_par[i];  //count the number of dm particle in a cell
      }
#endif
     vector<bool>dm_used( this->params._NGRID(),false);
     vector<int>dm_cases( this->params._NGRID(),0);
     //ACA DEBEMOS CONTEMPLAR 4 CASOS:
     // I) EN UNA CELDA HAY MAS TRACERS QUE DM
     // II) EN UNA CELDA HAY DM QUE TRACERS
     // III) En una celda hay Tracers pero no hay DM
     // IV) Celdas vacias, deben permanecer vacias
#ifdef _FULL_VERBOSE_
     So.message_screen("Total number of tracers =", Ntracers);
     So.message_screen("Total number of DM particles = ", posz.size());
     std::cout<<endl;
#endif

     ULONG empty_cells_original=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:empty_cells_original)
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       if(MOCK_DEN_FIELD[id]==0)
         {
           empty_cells_original++;
           dm_cases[id]=4;
           dm_used[id]=false;
         }
     //   So.message_screen("Original number of empty cells ", empty_cells_original);
#ifdef _FULL_VERBOSE_
    So.message_screen("Assigning coordinates and velocities to tracers from coordiantes of ALPT-DM particles");
#endif
     vector<int> aux_cont1( this->params._NGRID(),0);
     vector<int> aux_cont2( this->params._NGRID(),0);
     // define container of structure to allocate coordinates within each cell
     vector<s_cell_info> cell_inf_dm( this->params._NGRID());
     // caso I: mas trcers que dm
     for(ULONG i=0;i<N_dm;++i)
       {
         real_prec x = static_cast<real_prec>(posx[i]);
         real_prec y = static_cast<real_prec>(posy[i]);
         real_prec z = static_cast<real_prec>(posz[i]);
#ifndef _FCOL_TESTS_
         // The dm velocities are interpolated from the vel filed built in eulerian coords
         real_prec vx=vel_bias_dm*this->lpt.linInterp(x,y,z,velx);
         real_prec vy=vel_bias_dm*this->lpt.linInterp(x,y,z,vely);
         real_prec vz=vel_bias_dm*this->lpt.linInterp(x,y,z,velz);
//         real_prec vx= static_cast<real_prec>(vel_bias_dm*velx[i]);
//         real_prec vy= static_cast<real_prec>(vel_bias_dm*vely[i]);
//         real_prec vz= static_cast<real_prec>(vel_bias_dm*velz[i]);
#endif
         ULONG id   =  index[i];    //identify the cell where this particle lives
#ifdef _COLLAPSE_RANDOMS_AUX_
         ULONG id_collapse = grid_ID(&box_collps, x,y,z);
#endif
         if(MOCK_DEN_FIELD[id]>0 &&  aux_cont[id]>0)  // si hay tracers en esta celda y dm tambien
           {
             if(MOCK_DEN_FIELD[id]> aux_cont[id]) //si  hay mas o igual número tracers que dm (y hay dm), tomar todas las dm que hay en cada celda. Al haber mas tracers que dm, vendrá la necesidad de tener randoms
               {
                 cell_inf_dm[id].posx_p.push_back(x);
                 cell_inf_dm[id].posy_p.push_back(y);
                 cell_inf_dm[id].posz_p.push_back(z);
#ifndef _FCOL_TESTS_
                 cell_inf_dm[id].velx_p.push_back(vx);
                 cell_inf_dm[id].vely_p.push_back(vy);
                 cell_inf_dm[id].velz_p.push_back(vz);
#endif
                 cell_inf_dm[id].id_p.push_back(id);
                 dm_used[id]=true;
                 dm_cases[id]=1;
                 aux_cont1[id]++;
#ifdef _COLLAPSE_RANDOMS_AUX_
                 x_dm_pos.push_back(x);
                 y_dm_pos.push_back(y);
                 z_dm_pos.push_back(z);
                 dm_id.push_back(id_collapse);
                 dm_count[id_collapse]++;
                 N_dms++;
#endif
               }
             else
               { // caso II: se necesitan menor numero de tracers que el numero de dm en la celda. Ntr<=Ndm
                 if(aux_cont2[id]< MOCK_DEN_FIELD[id]) // el "if" es para tomar sólo los que necesitamos
                   {
                     cell_inf_dm[id].posx_p.push_back(x);
                     cell_inf_dm[id].posy_p.push_back(y);
                     cell_inf_dm[id].posz_p.push_back(z);
#ifndef _FCOL_TESTS_
                     cell_inf_dm[id].velx_p.push_back(vx);
                     cell_inf_dm[id].vely_p.push_back(vy);
                     cell_inf_dm[id].velz_p.push_back(vz);
#endif
                     cell_inf_dm[id].id_p.push_back(id);
                     dm_used[id]=true;
                     dm_cases[id]=2;
                     aux_cont2[id]++;
#ifdef _COLLAPSE_RANDOMS_AUX_
                     x_dm_pos.push_back(x);
                     y_dm_pos.push_back(y);
                     z_dm_pos.push_back(z);
                     dm_id.push_back(id_collapse);
                     dm_count[id_collapse]++;
                     N_dms++;
#endif
                   }
               }
           }
       }
     index.clear();
     index.shrink_to_fit();
     // Get the number of randoms from case i and ii
     vector<int>Nrandom_tracers( this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       if(MOCK_DEN_FIELD[id]>0)
         if(1==dm_cases[id] || 2==dm_cases[id])
           Nrandom_tracers[id]=MOCK_DEN_FIELD[id]-cell_inf_dm[id].posx_p.size();
     ULONG Ndm_used_I=0;
     ULONG Nrandoms1=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Ndm_used_I, Nrandoms1)
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       if(dm_cases[id]==1)
         {
           Ndm_used_I+=aux_cont1[id];
           Nrandoms1+=Nrandom_tracers[id];
         }
     // So.message_screen("Number of dm used in case I",Ndm_used_I);
     // So.message_screen("Number of randoms demanded from case I ", Nrandoms1);
     //    ULONG Ndm_used_II=0;
     // #pragma omp parallel for reduction(+:Ndm_used_II)
     //    for(ULONG id=0;id< this->params._NGRID();++id)
     //      Ndm_used_II+=aux_cont2[id];
     // So.message_screen("Number of dm used in case II",Ndm_used_II);
     // So.message_screen("(case II demands no randoms)");
     posx.clear();posx.shrink_to_fit();
     posy.clear();posy.shrink_to_fit();
     posz.clear();posz.shrink_to_fit();
     // So.message_screen("Total Number of dm used",Ndm_used_I+Ndm_used_II);
     // caso 3, celdas con tracers pero sin dm particles -> Todas random
     //   So.message_screen("Getting randoms from case 3");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       {
         if(MOCK_DEN_FIELD[id]>0 && 0==aux_cont[id])
           if(dm_cases[id]!=1 && dm_cases[id]!=2)
             {
               dm_used[id]=false;
               Nrandom_tracers[id]=MOCK_DEN_FIELD[id]; //todas randoms
               dm_cases[id]=3;
             }
       }
     ULONG Nrandoms3=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nrandoms3)
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       if(MOCK_DEN_FIELD[id]>0 && aux_cont[id]==0)
         if(3==dm_cases[id])
           Nrandoms3+=Nrandom_tracers[id];

     // So.message_screen("Total number of randoms requested from case 3", Nrandoms3);
     // So.message_screen("Total number of randoms requested", Nrandoms1+Nrandoms3);
     // So.message_screen("Randoms requested + Ndm", Nrandoms3+Nrandoms1+Ndm_used_I+Ndm_used_II);
     // So.message_screen("Number of original tracers",Nobjects_mock);
     Nrandoms1+=Nrandoms3;

     ULONG ncells_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ncells_check)
#endif
     for(ULONG id=0;id< this->params._NGRID();++id)
       {
         if(dm_cases[id]>=1 && dm_cases[id]<=4)
           ncells_check++;
       }
     if(ncells_check!= this->params._NGRID())
       {
         So.message_screen("Cells in cases", ncells_check);
         So.message_screen("Total number of cells",  this->params._NGRID());
         exit(0);
       }
     int jthread=1;
     vector<s_cell_info> cell_inf_ran( this->params._NGRID());
     vector<bool>random_used( this->params._NGRID(),false);
     if(Nrandoms1>0)
       {
#ifdef _FULL_VERBOSE_
         So.message_screen("Generating coordinates and velocities for random tracers");
#endif
         for(ULONG i=0; i<this->params._Nft(); ++i)
           for(ULONG j=0; j<this->params._Nft(); ++j)
             for(ULONG k=0; k<this->params._Nft(); ++k)
               {
                 ULONG id= index_3d(i,j,k,this->params._Nft(),this->params._Nft());
                 if(Nrandom_tracers[id]>0)
                   {
                     cell_inf_ran[id].posx_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].posy_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].posz_p.resize(Nrandom_tracers[id],0);
#ifndef _FCOL_TESTS_
                     cell_inf_ran[id].velx_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].vely_p.resize(Nrandom_tracers[id],0);
                     cell_inf_ran[id].velz_p.resize(Nrandom_tracers[id],0);
#endif
                     random_used[id]=true;
                     for(int ir=0;ir<Nrandom_tracers[id];++ir)
                       {
#ifdef _USE_OMP_
                         real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                         real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
#else
                         real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
#endif
                         real_prec x = d1*static_cast<real_prec>(i) + d1*rx;
                         real_prec y = d2*static_cast<real_prec>(j) + d2*ry;
                         real_prec z = d3*static_cast<real_prec>(k) + d3*rz;
                         cell_inf_ran[id].posx_p[ir]=x;
                         cell_inf_ran[id].posy_p[ir]=y;
                         cell_inf_ran[id].posz_p[ir]=z;
#ifndef _FCOL_TESTS_
                         cell_inf_ran[id].velx_p[ir]=vel_bias_ran*this->lpt.linInterp(x,y,z,velx);
                         cell_inf_ran[id].vely_p[ir]=vel_bias_ran*this->lpt.linInterp(x,y,z,vely);
                         cell_inf_ran[id].velz_p[ir]=vel_bias_ran*this->lpt.linInterp(x,y,z,velz);
#endif
                         cell_inf_ran[id].id_p.push_back(id);
#ifdef _COLLAPSE_RANDOMS_AUX_
                         x_random_pos.push_back(x);
                         y_random_pos.push_back(y);
                         z_random_pos.push_back(z);
                         ULONG id_collapse = grid_ID(&box_collps, x,y,z);
                         ran_id.push_back(id_collapse);
                         N_rand++;
#endif
                       }
                   }
               }
         So.DONE();
       }
#ifndef _FCOL_TESTS_
     velx.clear();velx.shrink_to_fit();
     vely.clear();vely.shrink_to_fit();
     velz.clear();velz.shrink_to_fit();
#endif

#ifdef _COLLAPSE_RANDOMS_AUX_
     int max_count=get_max(dm_count); // get the maximum number of dm particles in one cell.
#ifdef _FULL_VERBOSE_
     So.message_screen("Maximum number of DM particles found in one cell (low resolution) = ", max_count);
#endif
     // Initialize this to NO_NUM, a negative number, for not all elements of the container
     // dm_index_cell will be filled and initializing with 0 makes confusion, for 0 is a valid/used entry
#define NO_NUM -999
     vector<vector<int> > dm_index_cell(box_collps.NGRID, vector<int>(max_count,NO_NUM));
     dm_count.clear();
     dm_count.resize(box_collps.NGRID,0);
#ifdef _FULL_VERBOSE_
     So.message_screen("Getting ids in cells");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG idm=0;idm<N_dms;++idm)
       {
         ULONG id=dm_id[idm];
         dm_index_cell[id][dm_count[id]]=idm;
#ifdef _USE_OMP_
#pragma atomic update
#endif
         dm_count[id]++;
       }
     this->So.DONE();
     dm_id.clear();dm_id.shrink_to_fit();
#ifdef _FULL_VERBOSE_
     So.message_screen("Identifying closest DM particles for randoms at line", __LINE__);
#endif
     vector<int>dm_index_closer_tot(N_rand,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0; i< N_rand; ++i)//loop over randoms
       {
         ULONG count =0;
         int id_ran=ran_id[i];// id of the random
         int N_dm_cell=dm_count[id_ran];//Numbner of dm particles in the low-res cell where the random is located
         vector<ULONG>i_r_to_dm_dist;
         vector<ULONG>n_index_tot;
         if(N_dm_cell>0)// if there are dm particles in the low res cell where this random is located, then
           {
             for(ULONG jc = 0; jc < N_dm_cell; ++jc)//loop over the dm in that cell where the random is
               {
                 ULONG jdm=dm_index_cell[id_ran][jc]; //indice entre (0,Ndm-1) que identifica a cada particula DM dentro de la celda id_ran
                 if(jdm!=NO_NUM)
                   {
                     real_prec xdd=x_random_pos[i]-x_dm_pos[jdm];
                     real_prec ydd=y_random_pos[i]-y_dm_pos[jdm];
                     real_prec zdd=z_random_pos[i]-z_dm_pos[jdm];
                     ULONG i_dist=_get_modulo_squared(xdd,ydd,zdd);  //distance² between rand and each dm particle, converted to floor
                     i_r_to_dm_dist.push_back(i_dist);  // allocate the index for the distance
                     n_index_tot.push_back(jdm);        // allocat the id of the dm particle
                     count++;
                   }
               }
            }
         else
           So.message_screen("No dm particles found in the cell corresponding to random ", i, ". You might want to increase the cell-size");
         ULONG rank=0;  // OJO ACA: 0 es el mas cercano, 1 es el segundo más cercano
         sort_1d_vectors<ULONG>(i_r_to_dm_dist, rank); // sort the container with distance indices
         dm_index_closer_tot[i]=n_index_tot[rank];                // get the if dm in (0,.Ndm-1) particle associated to the first element of the sorted distance
         }
     So.DONE();
     dm_count.clear();
     dm_count.shrink_to_fit();
     ran_id.clear();ran_id.shrink_to_fit();
#endif  // end collapse_random_aux
     // Before write to file we assign masses. This replaces the call of the function in bamrunner
     // Here we pass the coordinates to a prop vctor and pass it as vector& to the mass assignment function,
     // to then pass it to collapse randoms
     ULONG Nobjects_mock_rc=0;
     this->tracer.Halo.resize(this->tracer._NOBJS());
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Filling structure with dm positions:");
#endif

// Problems with this parallelization. LEAVE COMMENTED!!!!!!!!!1
     for(ULONG id=0;id< this->params._NGRID();++id)
       {
         if(MOCK_DEN_FIELD[id]>0)
           if(true==dm_used[id])
             for(int in=0;in < cell_inf_dm[id].posx_p.size(); ++in)
               {
                 this->tracer.Halo[Nobjects_mock_rc].coord1 = cell_inf_dm[id].posx_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].coord2 = cell_inf_dm[id].posy_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].coord3 = cell_inf_dm[id].posz_p[in];
#ifndef _FCOL_TESTS_
                 this->tracer.Halo[Nobjects_mock_rc].vel1 = cell_inf_dm[id].velx_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].vel2 = cell_inf_dm[id].vely_p[in];
                 this->tracer.Halo[Nobjects_mock_rc].vel3 = cell_inf_dm[id].velz_p[in];
#endif
                 this->tracer.Halo[Nobjects_mock_rc].identity = ONE ;
                 this->tracer.Halo[Nobjects_mock_rc].GridID = cell_inf_dm[id].id_p[in];
                 Nobjects_mock_rc++;
               }
       }
     So.DONE();
#ifdef _FULL_VERBOSE_
     So.message_screen("Number of tracers associated to DM particles =", Nobjects_mock_rc);
#endif
     this->tracer.Ntracers_dm=Nobjects_mock_rc;
     this->tracer_ref.Ntracers_dm=Nobjects_mock_rc;
     this->tracer.fraction_tracer_from_dm=static_cast<real_prec>(Nobjects_mock_rc)/(static_cast<real_prec>(Nobjects_mock_rc)+static_cast<real_prec>(Nrandoms1)) ;
     this->tracer_ref.fraction_tracer_from_dm=this->tracer.fraction_tracer_from_dm;

#ifdef _FULL_VERBOSE_
     So.message_screen("(", 100.0*this->tracer.fraction_tracer_from_dm, "%)");
#ifdef _VERBOSE_FREEMEM_
     So.message_screen("Freeing memory");
#endif
#endif
     cell_inf_dm.clear(); cell_inf_dm.shrink_to_fit();
     So.DONE();
     ULONG N_dms_b= Nobjects_mock_rc;
     if(Nrandoms1>0)
       {
#ifdef _FULL_VERBOSE_
         this->So.message_screen("Adding information of random particles:");
#endif
//#ifdef _USE_OMP_
//#pragma omp parallel for reduction (+:Nobjects_mock_rc)  // lio
//#endif
         for(ULONG id=0;id< this->params._NGRID();++id)
           {
             if(MOCK_DEN_FIELD[id]>0)
               if(true==random_used[id])
                 for(ULONG in=0;in<cell_inf_ran[id].posx_p.size();++in)
                   {
                     this->tracer.Halo[Nobjects_mock_rc].coord1=cell_inf_ran[id].posx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord2=cell_inf_ran[id].posy_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].coord3=cell_inf_ran[id].posz_p[in];
#ifndef _FCOL_TESTS_
                     this->tracer.Halo[Nobjects_mock_rc].vel1=cell_inf_ran[id].velx_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel2=cell_inf_ran[id].vely_p[in];
                     this->tracer.Halo[Nobjects_mock_rc].vel3=cell_inf_ran[id].velz_p[in];
#endif
                     this->tracer.Halo[Nobjects_mock_rc].identity=-1 ;
                     this->tracer.Halo[Nobjects_mock_rc].GridID=cell_inf_ran[id].id_p[in];
                     Nobjects_mock_rc++; // This counter keeps on counting starting from the number of dm particles above
                   }
           }
#ifdef _FULL_VERBOSE_
         So.DONE();
         So.message_screen("Number of tracers associated to random particles =", Nrandoms1);
         this->tracer.Ntracers_ran=Nrandoms1;
         this->tracer_ref.Ntracers_ran=Nrandoms1;
         this->tracer.fraction_tracer_from_random=static_cast<real_prec>(Nrandoms1)/(static_cast<real_prec>(N_dms_b)+static_cast<real_prec>(Nrandoms1));
         this->tracer_ref.fraction_tracer_from_random=this->tracer.fraction_tracer_from_random;
         So.message_screen("(", 100.0*this->tracer.fraction_tracer_from_random, "%)");
#ifdef _VERBOSE_FREEMEM_
         So.message_screen("Freeing memory");
#endif
#endif
         cell_inf_ran.clear(); cell_inf_ran.shrink_to_fit();
         random_used.clear();
         random_used.shrink_to_fit();
#ifdef _FULL_VERBOSE_
         So.DONE();
#endif
       }
#ifdef _COLLAPSE_RANDOMS_AUX_
     // THE GrdID kept in the sturcture Tracer is the *ORIGINAL*, not the one computed after the collapse of randoms
#ifdef _FULL_VERBOSE_
     std::cout<<endl;
     this->So.message_screen("Using fraction of distance to closest dm particle = ", this->params._Distance_fraction());
     this->So.message_screen("Computing new position of randoms:");
#endif
       vector<ULONG> ran_id_global;
     ULONG counter=0;
//#ifdef _USE_OMP_
//#pragma omp parallel for reduction(+:counter)  // problem with paralelization and counter, leave commented
//#endif
     for(ULONG i=0; i<Ntracers; ++i)
       {
         if(this->tracer.Halo[i].identity<0)  //This means, use the randoms
           {
             ran_id_global.push_back(i);
             //total index (i.e, from 0 to N_dms) of the closest dm particle to the ramdom i. Useful to retrieve coordinates
             ULONG index_dm_closer_a=dm_index_closer_tot[counter];
             //Cartesian coordinates of the closest DM particle
             real_prec xdm=x_dm_pos[index_dm_closer_a];
             real_prec ydm=y_dm_pos[index_dm_closer_a];
             real_prec zdm=z_dm_pos[index_dm_closer_a];
             // redefine ran coords to the ref sistem of its closest dm particle:

#ifdef _COLLAPSE_RANDOMS_VELS_
             x_random_v.push_back(this->tracer.Halo[i].vel1);
             y_random_v.push_back(this->tracer.Halo[i].vel2);
             z_random_v.push_back(this->tracer.Halo[i].vel3);
#endif

             real_prec new_x=x_random_pos[counter]-xdm;
             real_prec new_y=y_random_pos[counter]-ydm;
             real_prec new_z=z_random_pos[counter]-zdm;
             real_prec dist_random_to_dm=_get_modulo(new_x,new_y,new_z);             // get the distance:
             //get the angular coordinates:
             real_prec theta=acos(new_z/dist_random_to_dm);
             real_prec phi=atan2(new_y,new_x);
             // get the new distance between the random and its closest dm particle
             real_prec new_distance = dist_random_to_dm*this->params._Distance_fraction();
             // transfor to cartesian given a reduced distance and return ot origin of coords:
             // Assign the new cartesian coordinates
             this->tracer.Halo[i].coord1=new_distance*sin(theta)*cos(phi)+xdm;
             this->tracer.Halo[i].coord2=new_distance*sin(theta)*sin(phi)+ydm;
             this->tracer.Halo[i].coord3=new_distance*cos(theta)+zdm;
             counter++;
             }
       }
     real_prec Fvcoll=1.0;
#ifdef _COLLAPSE_RANDOMS_VELS_
     So.message_screen("Collapsing velocities");
     for(ULONG counter=0; counter< N_rand; ++counter)
     {
         ULONG index_dm_closer_a = dm_index_closer_tot[counter];
         // Read coordinates of random particles
         real_prec new_x = x_random_pos[counter];
         real_prec new_y = y_random_pos[counter];
         real_prec new_z = z_random_pos[counter];
         real_prec xdm = x_dm_pos[index_dm_closer_a];
         real_prec ydm = y_dm_pos[index_dm_closer_a];
         real_prec zdm = z_dm_pos[index_dm_closer_a];
         //* Redefine ran coords to the ref sistem of its closest dm particle:*//
         new_x -= xdm;
         new_y -= ydm;
         new_z -= zdm;
         real_prec dist_random_to_dm = _get_modulo(new_x,new_y,new_z);

         // This section ortates the velocity-vector of the random particle towards
         // Original velocity components in the E system
         real_prec vx=x_random_v[counter];
         real_prec vy=y_random_v[counter];
         real_prec vz=z_random_v[counter];
         real_prec vdotr=vx*new_x+vy*new_y+vz*new_z;          //   Vec V cdot Delta r
         real_prec vv=_get_modulo(vx,vy,vz);   // |v|
         real_prec rxy=sqrt(pow(new_x,2)+pow(new_y,2));   //sqrt(dx²+dy²)

         // Coordinates of the velocity in the refernce frame E' (of the DM with hat y= hat r, r is the distance ran-dm closest)
         real_prec vxp=(new_y*vx-new_x*vy)/static_cast<double>(rxy); //vx'
         real_prec vyp=vdotr/dist_random_to_dm;// vy'=v cdot r/|r|
         real_prec vzp=(new_x*new_z*vx+new_y*new_z*vy)/(rxy*dist_random_to_dm)-(rxy/dist_random_to_dm)*vz; // vzp'
         real_prec vvp=_get_modulo(vxp,vyp,vzp); // modulo, just co cross-check with v
         real_prec alpha=acos(vyp/static_cast<double>(vv));//
         //* Get new velocity components in the E' system *//
         real_prec gamma_col=(1.0-this->params._Distance_fraction())*(alpha<M_PI ? (M_PI-alpha): (alpha -M_PI));// this is a model
//       nice idea: it does not help
/*         real_prec alpha_max=alpha<M_PI ? (M_PI-alpha): alpha;
         real_prec alpha_min=alpha<M_PI ? alpha: (M_PI-alpha);
         real_prec ran_angle= (alpha_min+0.5*(alpha_max-alpha_min))*(1+gsl_ran_gaussian(rng,M_PI/4.0)) ;//gsl_rng_uniform (rng);
         real_prec gamma_col=alpha_min+(alpha_max-alpha_min)*ran_angle; */
         real_prec uyp= alpha<M_PI ? Fvcoll*vv*cos(alpha+gamma_col) : Fvcoll*vv*cos(alpha-gamma_col);
         real_prec sigma=(Fvcoll*pow(vv,2)*cos(gamma_col)-vyp*uyp)/static_cast<double>(vxp);
         real_prec Gam=pow(Fvcoll*vv,2) - pow(uyp,2);
         double det= (dONE+ static_cast<double>((Gam/static_cast<double>(pow(sigma,2))-dONE)*(pow(vxp,2)+pow(vzp,2))/pow(static_cast<double>(vzp),2)));
         real_prec ddet=det <0 ? 0.0: det; // Si det es negativo, la raíz es imaginaria y nos sale una solución compleja, de la cual tomaremos sólo la parte real.
         real_prec uzp=(sigma*vzp*vxp/(pow(vxp,2)+pow(vzp,2))*(dONE+sqrt(ddet)));
         real_prec uxp=(pow(Fvcoll*vv,2)*cos(gamma_col)-vyp*uyp-vzp*uzp)/static_cast<double>(vxp);
         //* New velocity components in the E system *//
         real_prec ux=(new_y/rxy)*uxp+(new_x/dist_random_to_dm)*uyp+(new_x*new_z)*uzp/(rxy*dist_random_to_dm);
         real_prec uy=-(new_x/rxy)*uxp+(new_y/dist_random_to_dm)*uyp+(new_y*new_z)*uzp/(rxy*dist_random_to_dm);
         real_prec uz=(new_z/dist_random_to_dm)*uyp- rxy*uzp/dist_random_to_dm;
        // real_prec uup=_get_modulo(ux,uy,uz); // modulo, just co cross-check with v
         real_prec nux=std::isnan(ux) ? vx: ux;
         real_prec nuy=std::isnan(uy) ? vy: uy;
         real_prec nuz=std::isnan(uz) ? vz: uz;
         this->tracer.Halo[ran_id_global[counter]].vel1=nux;
         this->tracer.Halo[ran_id_global[counter]].vel2=nuy;
         this->tracer.Halo[ran_id_global[counter]].vel3=nuz;
     }
#endif
#ifdef _FULL_VERBOSE_
     this->So.DONE();
#endif
     x_random_pos.clear();
     x_random_pos.shrink_to_fit();
     y_random_pos.clear();
     y_random_pos.shrink_to_fit();
     z_random_pos.clear();
     z_random_pos.shrink_to_fit();
     x_dm_pos.clear();
     x_dm_pos.shrink_to_fit();
     y_dm_pos.clear();
     y_dm_pos.shrink_to_fit();
     z_dm_pos.clear();
     z_dm_pos.shrink_to_fit();
#endif // end of _COLLAPSE_RANDOMS_AUX_
     // This space was meant to allocate the vector of struturures tracer_ref.Halo[] with the same elements of the tracer.Halo
     // In order to save time, we are assuming that the mock number count used in the test "assign to reference"
     // is the same obtained from the catalog reference
     // This allocation cannot be done here, for the reference has not been read. It will be read in the function get_scaling_relations_primary_property()
    //  so we make it there internally.
#endif  //end of ifdef _ASSIGN_TO_MOCK or to calibration
     this->params.set_Name_survey("BAM");
#ifdef _WRITE_BINARY_BiasMT_FORMAT_
     string outputFileName=this->params._Output_directory()+"CAT_BiasMT_R"+to_string(this->params._realization())+".dat"; //   this->lpt.stradd_bam+string(".dat");
#else
     string outputFileName=this->params._Output_directory()+"CAT_R"+to_string(this->params._realization())+"_"+this->params._Name_survey()+".txt";
#endif
#ifdef _ASSIGN_PROPERTY_
#ifdef _USE_TWO_REFS_MOCKS_ASSIGNMENT_
     // Get info from the referencet
     this->tracers_multi.resize(this->params._Number_of_MultiLevels());
     for(int i=0;i<this->params._Number_of_MultiLevels();++i)
       this->tracers_multi[i].resize(this->params._Number_of_references(),0);

     for (int ifile=0; ifile< this->params._Number_of_references(); ++ifile)
       this->get_scaling_relations_primary_property_two(ifile); // accumulate the tracers incoming, in bins of THETA
#else
     this->get_scaling_relations_primary_property();
#endif
#ifdef _ASSIGN_TO_CALIBRATION_
     this->tracer.type_of_object="TRACER_MOCK_ONLY_COORDS";
#else
     this->tracer.set_type_of_object("TRACER_MOCK");
#endif
     string primary_prop;
#ifdef  _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     primary_prop = _VMAX_;
#else
     primary_prop = _MASS_;
#endif
      this->assign_tracer_property(true, primary_prop);
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
#ifdef _USE_GNUPLOT_SREL_
      this->plot_scaling_relation_assignment(primary_prop);
#endif
#endif

#ifdef _USE_HYBRID_ASSIGNMENT_NEW_
    this->So.message_warning("Checking hybrid_new");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
     real_prec delta_PROP=log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<real_prec>(pNbins);
#else
     real_prec delta_PROP=(this->params._LOGMASSmax()-this->params._LOGMASSmin())/static_cast<real_prec>(pNbins);
#endif
    ofstream sal;
    ULONG counter_swap=0;
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_  // If this is the case, we just put the assigned int ref_assignd. Otherwise we need to calulate it
#ifdef _USE_SIMD_OMP_
#pragma omp simd
#endif
    for(ULONG i=0; i<this->tracer_ref._NOBJS();++i)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      this->tracer.Halo[i].vmax_parent=this->tracer_ref.Halo[i].vmax;
#else
        this->tracer.Halo[i].mass_parent=this->tracer_ref.Halo[i].mass;
#endif
#endif
    // keep track of the assigned vmax to test method below
#ifdef _USE_SIMD_OMP_
#pragma omp simd
#endif
    for(ULONG i=0; i<this->tracer_ref._NOBJS();++i)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      this->tracer.Halo[i].vmax_assigned=this->tracer.Halo[i].vmax;
#else
      this->tracer.Halo[i].mass_assigned=this->tracer.Halo[i].mass;
#endif
    const gsl_rng_type * rng_s;
    gsl_rng * gBaseRands;
    int jthread=0;
    gsl_rng_env_setup();
    gsl_rng_default_seed=30455;
    rng_s = gsl_rng_mt19937;
    gBaseRands = gsl_rng_alloc (rng_s);
    So.message_screen("Computing some ideal scatter based on a gaussian scaling relation");
    //This must be such that the bias in power is small wrt the reference
    real_prec ideal_corr=0;
#pragma omp parallel for reduction (+:ideal_corr)
    for(ULONG i=0; i<this->tracer._NOBJS();++i)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        ideal_corr+=pow(log10(this->tracer.Halo[i].vmax) - log10(this->tracer.Halo[i].vmax_parent),2);
#else
        ideal_corr+=pow(log10(this->tracer.Halo[i].mass) - log10(this->tracer.Halo[i].mass_parent),2);
#endif
    ideal_corr=sqrt(ideal_corr/static_cast<real_prec>(this->tracer._NOBJS()));
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    So.message_screen("Initial Correlation =",ideal_corr," in log10[km/s]");
#else
    So.message_screen("Initial Correlation =",ideal_corr,"Ms/h");
#endif
    ideal_corr=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction (+:ideal_corr)
#endif
    for(ULONG i=0; i<this->tracer._NOBJS();++i)
        ideal_corr+=pow(gsl_ran_gaussian(gBaseRands,SIGMA_REC),2);
    ideal_corr=sqrt(ideal_corr/static_cast<real_prec>(this->tracer._NOBJS()));
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    So.message_screen("Imposed Correlation =",ideal_corr," in log V[km/s]");
#else
    So.message_screen("Imposed Correlation =",ideal_corr,"in log 10 M [Ms/h]");
#endif

    // It is use just as a diagnosis for plotting
    ULONG counter_rec=0;
    //* This loop aims at obtaining a smoothn trnaistion probability by evaluating a number Nr of times *//
    //* the swapping process changing every time the seed  *//
    ULONG Nreal = 1;
    for (ULONG ir=0;ir<Nreal;++ir)
      {
#ifdef _USE_OMP_re
        vector<ULONG>vseeds(_NTHREADS_,0);
        for(int i=0;i<vseeds.size();++i)
          vseeds[i]=353*ir+static_cast<ULONG>(i)*53215*pow(ir,3);
#endif
        counter_swap=0;
        real_prec corr=1e14;// * This has to be a large number. For M as prim-prop, correlations are of the order of 10e12
        while(corr>ideal_corr)//* Start a while loop until the stvd of the Vprop_prop is that imposed in the Gaussian model*/
          {
#ifdef _USE_OMP_re
#pragma omp parallel private (jthread, gBaseRands, rng_s)
            {
              jthread=omp_get_thread_num();
              gsl_rng_default_seed=vseeds[jthread];
#else
              gsl_rng_default_seed=35;
#endif
              rng_s = gsl_rng_mt19937;//_default;
              gBaseRands = gsl_rng_alloc (rng_s);
#ifdef _USE_OMP_re
#pragma omp for reduction(+:counter_swap)
#endif
              for(ULONG idm=0; idm<this->dm_properties_bins_mock.size();++idm)
                {
                  ULONG Number_of_props_in_Theta_bin=this->dm_properties_bins_mock[idm].tracer_properties.size();
                  if(Number_of_props_in_Theta_bin>0)
                    {
                      vector<ULONG>prop_bin_counter(pNbins);//* for histograms of prop within each theta bin *//
                      vector<vector<ULONG>> id_in_propbin(pNbins*pNbins);//* to allocate gals ID in bins of main propery *//
                      for(ULONG idm_tracer=0;idm_tracer<Number_of_props_in_Theta_bin;++idm_tracer)
                        {
                          ULONG ig=this->dm_properties_bins_mock[idm].GalID_bin_properties[idm_tracer];
                          //* Get the value of vmax of the tracer: *//
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                       //* Get bin of primary_property from reference or parent *//
                          ULONG ibin_prop=get_bin(log10(this->tracer_ref.Halo[ig].vmax),log10(this->params._VMAXmin()),pNbins,delta_PROP,true);
                  //* Get the value of vmax of the tracer_ref:  *//
                  //* Get bin of vmax *//
                         ULONG jbin_prop=get_bin(log10(this->tracer.Halo[ig].vmax),log10(this->params._VMAXmin()),pNbins,delta_PROP,true);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                      //* Get bin of vmax *//
                          ULONG ibin_prop = get_bin(log10(this->tracer_ref.Halo[ig].mass),this->params._LOGMASSmin(),pNbins,delta_PROP,true);
                      //* Get bin of vmax *//
                          ULONG jbin_prop = get_bin(log10(this->tracer.Halo[ig].mass),this->params._LOGMASSmin(),pNbins,delta_PROP,true);
#endif
                          prop_bin_counter[ibin_prop]++;//* Count # of objects in prop bin in each theta bin *//
                          id_in_propbin[index_2d(ibin_prop,jbin_prop,pNbins)].push_back(ig);//*Allocate the GalID in the theta bin of the prop in this prop bin*//
                    }
                   for(ULONG idm_tracer=0;idm_tracer<Number_of_props_in_Theta_bin;++idm_tracer)
                    {
			  //* Read the id of the chosen tracer *//
                      ULONG ig=this->dm_properties_bins_mock[idm].GalID_bin_properties[idm_tracer];
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                      real_prec vmax_tracer     =  this->tracer.Halo[ig].vmax;//* get the value of vmax assigned to the tracer, renewed below swaped from another bin *//
                      real_prec vmax_tracer_ref =  this->tracer_ref.Halo[ig].vmax; //* get the value of the original vmax of the tracer *//
                      ULONG ibin_prop = get_bin(log10(vmax_tracer_ref),log10(this->params._VMAXmin()),pNbins,delta_PROP,true); //* get i bin of vmax *//
                      ULONG jbin_prop = get_bin(log10(vmax_tracer),log10(this->params._VMAXmin()),pNbins,delta_PROP,true); //* get j bin of vmax *//
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                      real_prec mass_tracer     =  this->tracer.Halo[ig].mass;//* get the value of vmax assigned to the tracer, renewed below swaped from another bin *//
                      real_prec mass_tracer_ref =  this->tracer_ref.Halo[ig].mass; //* get the value of the original vmax of the tracer *//
                      ULONG ibin_prop = get_bin(log10(mass_tracer_ref),this->params._LOGMASSmin(),pNbins,delta_PROP,true); //* get i bin of prop in reference *//
                      ULONG jbin_prop = get_bin(log10(mass_tracer),this->params._LOGMASSmin(),pNbins,delta_PROP,true); //* get j bin of prop in  *//
#endif
              //* This defines the diagonal bins as exlusion zone, i.e, we do not swap from those bins *//
                  bool exclusion_zone=false;
                  if(fabs(jbin_prop-ibin_prop)<Exclusion_bins)
                    exclusion_zone=true;
                  if(false==exclusion_zone)// if we are in a bin that needs to be collapsed towards the diagonal
                    {
                      int dir_jump=0; // this gives ths direction of the jump and automatically defines the direction of the swap/jump
                      if(ibin_prop>jbin_prop)// en este caso enviamos la propiedad hacia arriba (para acercarse a la diagonal)
                        dir_jump=1;
                      else if (ibin_prop<jbin_prop)// en este caso enviamos vmax hacia abajo (para acercarse a la diagonal)
                        dir_jump=-1;
                  // elige el número de bins en prop a saltar hacia arriba. Este va desde 1 hasta el i-bin:(ver dibujo)
                      ULONG diff = jbin_prop>ibin_prop ? jbin_prop-ibin_prop : ibin_prop-jbin_prop;
                      ULONG ncells_jumps= diff>1 ? 2 + gsl_rng_uniform_int(gBaseRands, diff)  : 0 ;
                  if(ncells_jumps>0)
                    {
                      long new_jbin_prop  = jbin_prop + dir_jump*ncells_jumps;   //define new j-bin to swap
                      if(new_jbin_prop>=pNbins)//* do not go beyond the numbre of bins *//
                        new_jbin_prop=pNbins-1;
                      if(new_jbin_prop<0)//* do not go beyond the numbre of bins  *//
                        new_jbin_prop=0;
                      ULONG Ntracers_in_new_prop_jbin=prop_bin_counter[new_jbin_prop];//* Get the number of objects in the new j-bin *//
                      if(Ntracers_in_new_prop_jbin>0) //*Be sure that there are properties in the new full j-bin: *//
                        {
                      //* Define the range in which we will find the new i_bin: *//
				      //* this depends on the side of the diagonal in which the selected tracer is (in the vmax-vmax plane) *//
				      //* Starting bin of new i array: *//
                          ULONG baux= (dir_jump> 0 ? 0   : new_jbin_prop-(Exclusion_bins-1));
				      //* Number of bins in the i-dir to use: this is aiming at selecting a new i bin such that when swaping, we do not go above (dir>0) oe below the diagonal*//
                          ULONG naux= (dir_jump> 0 ? ibin_prop-(Exclusion_bins-1) : pNbins-baux); //* asegurarse de qe del lado de la diagonal tenemos celdas con tracers: //
                          vector<ULONG>noeb;
                          for(ULONG ip=baux;ip<naux;++ip)
                            if(id_in_propbin[index_2d(ip,new_jbin_prop,pNbins)].size()>0)
                              noeb.push_back(ip);
                          if(noeb.size()>0)
                            {
                               ULONG new_ibin_prop=noeb[gsl_rng_uniform_int(gBaseRands,noeb.size())];
                               ULONG N_props_in_new_ijbin=id_in_propbin[index_2d(new_ibin_prop,new_jbin_prop,pNbins)].size();
                               if(N_props_in_new_ijbin>0)
                                {
                                  ULONG idprop=gsl_rng_uniform_int(gBaseRands,N_props_in_new_ijbin);
                                  ULONG new_ig = id_in_propbin[index_2d(new_ibin_prop,new_jbin_prop,pNbins)][idprop];
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                                  real_prec vmax_tracer_new= this->tracer.Halo[new_ig].vmax;
                                  this->tracer.Halo[ig].vmax=vmax_tracer_new; //* SWAP *//
                                  this->tracer.Halo[new_ig].vmax=vmax_tracer; //* SWAP *//
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
                                  real_prec mass_tracer_new= this->tracer.Halo[new_ig].mass;
                                  this->tracer.Halo[ig].mass=mass_tracer_new; //* SWAP *//
                                  this->tracer.Halo[new_ig].mass=mass_tracer; //* SWAP *//
#endif
                                  counter_swap++;
                                }//* close      if(N_props_in_new_ijbin>0) *//
                            }//* close  if(noeb.size()>0)  *//
                        }//* close  if(Ntracers_in_new_prop_jbin>0) *//
				}//* close if(ncells_jumps>0) *//
			    }//* close  if(false==exclusion_zone) *//
			}//* close for(ULONG idm_tracer=0 *//
		    }//* close  if(Number_of_props_in_Theta_bin>0) *//
		}//* close  for(ULONG idm=0.. *//
#ifdef _USE_OMP_re
	      gsl_rng_free (gBaseRands);
	    }//* close #pragma omp parallel *//
#endif
        // Now perform some statistics: //
	    real_prec ecorr=0;
#pragma omp parallel for reduction(+:ecorr)
	    for(ULONG i=0;i<this->tracer._NOBJS(); ++i)
	      {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
        ecorr+=pow(log10(this->tracer.Halo[i].vmax)-log10(this->tracer.Halo[i].vmax_parent),2);
#else
         ecorr+=pow(log10(this->tracer.Halo[i].mass)-log10(this->tracer.Halo[i].mass_parent),2);
#endif
        }
        corr=sqrt(ecorr/static_cast<real_prec>(this->tracer._NOBJS()));
        real_prec ebias=0;
#pragma omp parallel for reduction(+:ebias)
	    for(ULONG i=0;i<this->tracer._NOBJS(); ++i)
	      {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
            ebias+=this->tracer.Halo[i].vmax/this->tracer.Halo[i].vmax_parent;
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
            ebias+=this->tracer.Halo[i].mass/this->tracer.Halo[i].mass_parent;
#endif
         }
        ebias=sqrt(ebias/static_cast<real_prec>(this->tracer._NOBJS()));
        real_prec sbias=0;
#pragma omp parallel for reduction(+:ebias)
	    for(ULONG i=0;i<this->tracer._NOBJS(); ++i)
	      {
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
            sbias+=pow(this->tracer.Halo[i].vmax/this->tracer.Halo[i].vmax_parent-1,2);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
            sbias+=pow(this->tracer.Halo[i].mass/this->tracer.Halo[i].mass_parent-1,2);
#endif
	      }
        sbias=sqrt(sbias/static_cast<real_prec>(this->tracer._NOBJS()));
        So.message_screen_flush("Sigma=", corr, "Bias=", ebias, "sigma_bias=", sbias);
#ifdef _USE_GNUPLOT_SREL_
	    this->plot_scaling_relation_assignment();
        this->get_mean_scaling_relation_assignment_bias(_VMAX_);
#endif
	    counter_rec++;
          }//* close while *//
    //* The next lines aim to compute the transition probability by checking the initial (assigned) and final bin of Vmax in which the *//
	//* tracer ends after the "relaxation" *//
#ifdef _CHECK_HYBRID_V2_
#ifdef _USE_OMP_re
#pragma omp parallel for
#endif
	for(ULONG idm=0; idm<this->dm_properties_bins_mock.size();++idm)
	  {
	    ULONG Number_of_props_in_Theta_bin=this->dm_properties_bins_mock[idm].tracer_properties.size();
            this->dm_properties_bins_mock[idm].Pfj_upward.resize(pNbins*pNbins,0);
            this->dm_properties_bins_mock[idm].Pfj_downward.resize(pNbins*pNbins,0);
            this->dm_properties_bins_mock[idm].Pjj.resize(pNbins*pNbins,0);
            vector<ULONG> prop_bin_counter(pNbins,0);
	    if(Number_of_props_in_Theta_bin>0)
	      {
		for(ULONG idm_tracer=0;idm_tracer<Number_of_props_in_Theta_bin;++idm_tracer)
		  {
		    ULONG ig=this->dm_properties_bins_mock[idm].GalID_bin_properties[idm_tracer];
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                    prop_bin_counter[get_bin(log10(this->tracer.Halo[ig].vmax_assigned),log10(this->params._VMAXmin()),pNbins,delta_PROP,true)]++;
#else
                    prop_bin_counter[get_bin(log10(this->tracer.Halo[ig].mass_assigned),this->params._LOGMASSmin(),pNbins,delta_PROP,true)]++;
#endif
                }
		for(ULONG idm_tracer=0;idm_tracer<Number_of_props_in_Theta_bin;++idm_tracer)
		  {
		    ULONG ig=this->dm_properties_bins_mock[idm].GalID_bin_properties[idm_tracer];


#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
                    ULONG jbin_prop_assigned=get_bin(log10(this->tracer.Halo[ig].vmax_assigned),log10(this->params._VMAXmin()),pNbins,delta_PROP,true);
                    ULONG jbin_prop_new=get_bin(log10(this->tracer.Halo[ig].vmax),log10(this->params._VMAXmin()),pNbins,delta_PROP,true);
#else
                    ULONG jbin_prop_assigned=get_bin(log10(this->tracer.Halo[ig].mass_assigned),this->params._LOGMASSmin(),pNbins,delta_PROP,true);
                    ULONG jbin_prop_new=get_bin(log10(this->tracer.Halo[ig].mass),this->params._LOGMASSmin(),pNbins,delta_PROP,true);
#endif
                    this->dm_properties_bins_mock[idm].Pjj[index_2d(jbin_prop_assigned,jbin_prop_new,pNbins)]+=1./static_cast<real_prec>(Nreal);

                    //* Exclude diagonals*//
                    if(jbin_prop_assigned<jbin_prop_new) //* These are those who went up  *//
		      {
#ifdef _USE_OMP_re
#pragma omp atomic
#endif
                        //* Counting the fraction of tracers in this Vmax bin going up from a bin jbin_prop_assigned to a bin jbin_prop_new  *//
                       this->dm_properties_bins_mock[idm].Pfj_upward[index_2d(jbin_prop_assigned,jbin_prop_new,pNbins)]+=1./static_cast<real_prec>(Nreal);
			
		      } //* Cloese if() *//
                    else if(jbin_prop_assigned>jbin_prop_new) // * These are now those that went down*//)
		      {
#ifdef _USE_OMP_re
#pragma omp atomic
#endif
                        this->dm_properties_bins_mock[idm].Pfj_downward[index_2d(jbin_prop_assigned,jbin_prop_new,pNbins)]+=1./static_cast<real_prec>(Nreal); // Transition probability between original bin j and new bin j
		      } //* Closes else*//
		  }//Closes loop*//
	      }//* closes if *//
	  }//Closes loop*//

    //* Go back to the assignmment to start the process again*//
#pragma omp parallel for
	for(ULONG i=0; i<this->tracer._NOBJS();++i)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
          this->tracer.Halo[i].vmax=this->tracer.Halo[i].vmax_assigned;
#else
            this->tracer.Halo[i].mass=this->tracer.Halo[i].mass_assigned;
#endif


#endif   //#if defined _CHECK_HYBRID_  || defined _CHECK_HYBRID_V2_
    } //*  Closes loop over Nr *//
    So.message_screen("Swap", counter_swap, "props");
    //* -----------------------------------------------------------------------------------------------------*//
#ifdef _CHECK_HYBRID_V2_
    So.message_screen("Testing v2");
    counter_swap=0;
    for(int it=0; it <1; ++it)
      {
#ifdef _USE_OMP_rea
#pragma omp parallel private (jthread, gBaseRands, rng_s)
	{
	  jthread=omp_get_thread_num();
	  gsl_rng_default_seed=vseeds[jthread];
#else
	  gsl_rng_default_seed=35845;
#endif
	  rng_s = gsl_rng_mt19937;//_default;
	  gBaseRands = gsl_rng_alloc (rng_s);
      vector<bool> swap(this->tracer._NOBJS(), false);
#ifdef _USE_OMP_rea
#pragma omp for
#endif
	  for(ULONG idm=0; idm<this->dm_properties_bins_mock.size();++idm)
	    {
	      ULONG Number_of_props_in_Theta_bin=this->dm_properties_bins_mock[idm].tracer_properties.size();
         if(Number_of_props_in_Theta_bin>0)
		{
          vector<vector<ULONG>> id_in_propbin(pNbins);
           vector<ULONG> prop_bin_counter(pNbins,0);
          for(ULONG idm_tracer=0;idm_tracer<Number_of_props_in_Theta_bin;++idm_tracer)
		    {
		      ULONG ig=this->dm_properties_bins_mock[idm].GalID_bin_properties[idm_tracer];
		      real_prec vmax_tracer = this->tracer.Halo[ig].vmax_assigned;//get the value of vmax assigned to  the tracer, renewed below swaped from another bin
                      ULONG jbin_prop = get_bin(log10(vmax_tracer),log10(this->params._VMAXmin()),pNbins,delta_PROP,true); // get j bin of vmax
                      id_in_propbin[jbin_prop].push_back(ig);
                      prop_bin_counter[get_bin(log10(vmax_tracer),log10(this->params._VMAXmin()),pNbins,delta_PROP,true)]++;
		    }
                  //* Container to allocate ID of galaxies which in each vmax bin went up*/
		  vector<vector<ULONG>>ig_up(pNbins*pNbins);
                  //* Container to allocate ID os galaxies which in each vmax bin went down*/
		  vector<vector<ULONG>>ig_down(pNbins*pNbins);
          //* Do loop over bins of vmax *//
		  for(ULONG iv=0;iv<pNbins;++iv)
		    {
		      //* In this bin: *//
		      //* Identify the number of tracers *//
                      ULONG Ntracers_vbin=id_in_propbin[iv].size();
		      if(Ntracers_vbin>0)
			{
              //* Get from the learning phase the fraction of tracers moved upwards from the current bin to the other bins up *//
           //* Choose the tracers that in the bin iv went up to different bins .//
			  vector<bool>sel(Ntracers_vbin, false);
              ULONG Nup_m=0;
			  for(ULONG iav=0;iav<pNbins;++iav)
			    Nup_m+=this->dm_properties_bins_mock[idm].Pfj_upward[index_2d(iv,iav,pNbins)];	  
              if(Nup_m>0)
			    {
			      for(ULONG iv_up=iv+1;iv_up<pNbins;++iv_up)
				{
				  ULONG index=index_2d(iv,iv_up,pNbins);
                                  // Number of tracers going up
				  ULONG Ntracers_bin_bin_up=static_cast<ULONG>(this->dm_properties_bins_mock[idm].Pfj_upward[index]);
                                  //Number of tracers in th bin j j'
                 ULONG Ntracers_bin_bin_all=static_cast<ULONG>(this->dm_properties_bins_mock[idm].Pjj[index]);
                  if(Ntracers_bin_bin_up>0)
				    {
				      ULONG counter_du=0;				
				      while(counter_du<Ntracers_bin_bin_up)
					{
					  //* This is the key point: here we decide who is going up or who is going down*//
					  //* And here is where we can mix things.*//
                                          ULONG igr = gsl_rng_uniform_int(gBaseRands,Ntracers_vbin);
                                          ULONG gal_id=id_in_propbin[iv][igr];
                                          ig_up[index].push_back(gal_id);
                                          sel[igr]=true;
					  counter_du++;
					}
				    }
				}
			    }
                          ULONG counter_partial=0;
                          for(ULONG ip=0;ip<sel.size();++ip)
                            if(false==sel[ip])
                              counter_partial++;
                          //* Select the Ndown tracers. We use another while-loop as there might be tracers which did not move, hence fup!=1-fdown *//
                          ULONG Ndown_m=0;
                          for(ULONG iav=0;iav<pNbins;++iav)
                            Ndown_m+=this->dm_properties_bins_mock[idm].Pfj_downward[index_2d(iv,iav,pNbins)];
                          if(Ndown_m>0)
			    {
                              for(ULONG iv_down=0;iv_down<iv ;++iv_down)
				{
				  ULONG index=index_2d(iv,iv_down,pNbins);
				  ULONG Ntracers_bin_bin_down=static_cast<ULONG>(this->dm_properties_bins_mock[idm].Pfj_downward[index]);
                                  ULONG Ntracers_bin_bin_all=static_cast<ULONG>(this->dm_properties_bins_mock[idm].Pjj[index]);
				  if(Ntracers_bin_bin_down>0)
				    {
				      ULONG counter_du=0;
                                      while(counter_du< Ntracers_bin_bin_down)
					{
                                          ULONG igr = gsl_rng_uniform_int(gBaseRands,Ntracers_vbin);//Ntracers_bin_bin_all
                                          ULONG gal_id=id_in_propbin[iv][igr];
//                                          if(false==sel[gal_id])
					    {
                                              ig_down[index].push_back(gal_id);
                                              sel[igr]=true;
					      counter_du++;
					    }
					}
				    }
				}
			    }
                        }
		    }
		  
		  //*  At this point we have, for each Vmax bin, identified (randomly though) the number of tracers going up and down *//
		  //* So now we go bin by bin again:*//
		  // loop over parent bin dedicated to push objects up (and bring down those at thee arrival bin):*//
		  for(ULONG iv_parent=0;iv_parent<pNbins-1;++iv_parent)
		    {
                      ULONG Ntracers_vbin=id_in_propbin[iv_parent].size();
		      if(Ntracers_vbin>0)
			{
			  for(ULONG iv_final=0;iv_final<pNbins-1;++iv_final)
			    {
			      ULONG index_parent_final=index_2d(iv_parent, iv_final, pNbins);
			      //* Number of tracers going up in from bin iv_parent to bin iv_final*//
			      ULONG Ntracers_up = ig_up[index_parent_final].size();
			      //* Loop over the tracers ready to go-up from the iv_parent to the bin iv_final: *//
			      for(ULONG i_up = 0; i_up< Ntracers_up;++i_up)
				{
				  //* For the current tracer to be pushed-up: *//
				  ULONG Ig_parent=ig_up[index_parent_final][i_up];
				  //if(false==swap[Ig_parent])
				    {
				      //* Read its vmax *//
				      real_prec vmax_parent = this->tracer.Halo[Ig_parent].vmax_assigned;
				      //* If that bin has trcers willing to go down, *//
                                     //                                   cout<<ig_down[index_parent_final].size()<<endl;
				      if(ig_down[index_parent_final].size()>0)
					{
					  //* Pick up the id one tracer in that bin willing to go down :*//
					  ULONG Ig_new = ig_down[index_parent_final][gsl_rng_uniform_int(gBaseRands,ig_down[index_parent_final].size())];
					  //* Read its vmax *//
					  real_prec vmax_new=this->tracer.Halo[Ig_new].vmax_assigned;
					  //swap only if these have not been swaped before
					  //if(false==swap[Ig_new])
					    {
					      //					      cout<<iv_parent<<"  "<<iv_final<<endl;
					      this->tracer.Halo[Ig_new].vmax_assigned=vmax_parent;
					      this->tracer.Halo[Ig_parent].vmax_assigned=vmax_new;
					      swap[Ig_new]=true;
					      swap[Ig_parent]=true;
					      counter_swap++;
					    }//* closes if*/
					}//* closes if *//
				    }//* closes for*//
				}//* closes if*//
			    }//* closes for*//
			}//* closes if*//
		    }//* closes for loop*//

		  // loop over parent bin dedicated to push objects down (and bring up those at thee arrival bin):*//
		  for(ULONG iv_parent=pNbins-1;iv_parent--> 1;)
		    {
                      ULONG Ntracers_vbin=id_in_propbin[iv_parent].size();
		      if(Ntracers_vbin>0)
			{
			  for(ULONG iv_final=0;iv_final<iv_parent;++iv_final)
			    {
			      ULONG index_parent_final=index_2d(iv_parent, iv_final, pNbins);
			      
			      ULONG Ntracers_down=ig_down[index_parent_final].size();
			      //* Loop over the tracers ready to go-up
			      for(ULONG i_down=0; i_down< Ntracers_down;++i_down)
				{
				  //* Select one tracer to be pushed-up:
				  ULONG Ig_parent=ig_down[index_parent_final][i_down];
				  if(false==swap[Ig_parent])
				    {
				      real_prec vmax_parent=this->tracer.Halo[Ig_parent].vmax_assigned;
				      //* New bin down: *//
				      //* Pick up one tracer in that bin willing to go up:
				      if(ig_up[index_parent_final].size()>0)
					{
					  ULONG Ig_new = ig_up[index_parent_final][gsl_rng_uniform_int(gBaseRands,ig_up[index_parent_final].size())];
					  real_prec vmax_new = this->tracer.Halo[Ig_new].vmax_assigned;
					  //swap:
					  if(false==swap[Ig_new])
					    {
					      this->tracer.Halo[Ig_new].vmax_assigned=vmax_parent;
					      this->tracer.Halo[Ig_parent].vmax_assigned=vmax_new;
					      swap[Ig_new]=true;
					      swap[Ig_parent]=true;
					      counter_swap++;
					    }//* closes if *//
					} //* closes if*//
				    }//* closes if*//
				}//* closes if*//
			    }//*closes for loop*//
			}//* closes if*//
		    }//* closes for loop*//
		}
	    }
#ifdef _USE_OMP_rea
	}//* closes parallel regions *//
#endif
	//
#pragma omp parallel for
    for(ULONG i=0; i<this->tracer._NOBJS();++i)
      this->tracer.Halo[i].vmax=this->tracer.Halo[i].vmax_assigned;
    real_prec ecorr=0;
#pragma omp parallel for reduction(+:ecorr)
    for(ULONG i=0;i<this->tracer._NOBJS(); ++i)
      ecorr+=pow(this->tracer.Halo[i].vmax-this->tracer.Halo[i].vmax_parent,2);
    real_prec corr=sqrt(ecorr/static_cast<real_prec>(this->tracer._NOBJS()));
    So.message_screen_flush("Sigma ",corr);
#ifdef _USE_GNUPLOT_SREL_
    this->plot_scaling_relation_assignment();
#endif
    }//*Closes loop over trials with label it **//
    So.DONE();
    cout<<endl;
    So.message_screen("NUmber of swap pairs:",counter_swap);
    
#endif // end for _CHECK_HYBRID_V2_
    this->dm_properties_bins_mock.clear();
    this->dm_properties_bins_mock.shrink_to_fit();
    So.DONE();
#ifdef _USE_GNUPLOT_
    this->gp_abundance_v<<"set log x \n";
    this->tracer.params.set_i_spin_g(NEGATIVE_INT); //no spin information here
    this->tracer.params.set_i_rs_g(NEGATIVE_INT); //no rs information here
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    this->tracer.params.set_i_vmax_g(POSITIVE_INT); //this allows the tracer to get the vmax function
    this->tracer.params.set_i_mass_g(NEGATIVE_INT); //no mass information here
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    this->tracer.params.set_i_vmax_g(NEGATIVE_INT); //no vmax  information here
    this->tracer.params.set_i_mass_g(POSITIVE_INT); //this allows the tracer to get the mass function
#endif
    this->tracer.get_property_function(this->params._Output_directory()+"test");
    this->tracer_ref.get_property_function(this->params._Output_directory()+"test");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    this->gp_abundance_v<<"set xlabel 'Vrms [km/s]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set ylabel 'log n(Vrms)' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    vector<pair<real_prec, real_prec> > xy_pts_v;
    for(ULONG i=0; i<this->tracer.vmax_function.size(); ++i)
      xy_pts_v.push_back(std::make_pair(this->tracer.VMAXBmin[i]+(i+0.5)*(this->tracer.VMAXBmax[i]-this->tracer.VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.vmax_function[i])));
    this->gp_abundance_v<<"plot [120:1000]"<<this->gp_abundance_v.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Rec' "<<endl;
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    this->gp_abundance_v<<"set xlabel 'Mvir [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_abundance_v<<"set ylabel 'log n(Mvir)' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    vector<pair<real_prec, real_prec> > xy_pts_v;
    for(ULONG i=0; i<this->tracer.mass_function.size(); ++i)
      xy_pts_v.push_back(std::make_pair(this->tracer.MBmin[i]+(i+0.5)*(this->tracer.MBmax[i]-this->tracer.MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->tracer.mass_function[i])));
    this->gp_abundance_v<<"plot [11:15]"<<this->gp_abundance_v.file1d(xy_pts_v) << "w l lw 3 lt 3 title 'Rec' "<<endl;
#endif
#endif
//    exit(0);
#endif   //end for _USE_HYBRID_ASSIGNMENT_NEW_
#if defined _ASSIGN_MASS_POST_ || defined _ASSIGN_VMAX_POST_
    // Here we can now complement the property assignment, using the information of vmax already assigned
#ifdef _FULL_VERBOSE_
    So.message_screen("******************************************************");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    So.message_screen("Assigning Mvir using VMAX information::");
#else
    So.message_screen("Assigning Vmax using Mvir information::");
#endif
    So.message_screen("******************************************************");
    std::cout<<endl;
#endif
    this->params.set_i_vmax_g(POSITIVE_INT);// this number only needs to be positive
    this->tracer_ref.set_type_of_object("TRACER_REF"); // THis name is important.
    this->tracer.set_type_of_object("TRACER_MOCK"); // THis name is important.
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
    if(this->counter_assigned_secondary>0)
      {
#endif
    this->get_scaling_relations_secondary_property(_MASS_); // Learn P(M|V,delta) from reference if Vmax has been assigned
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
      }
#endif
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    this->get_scaling_relations_secondary_property(_VMAX_);   // Learn P(sigma|M,delta) from reference if M has been assigned
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
     if(this->counter_assigned_secondary>0)
#endif
       this->assign_tracer_property(false, _MASS_);  // Drawn M from P(M|V_mock=Vref, delta_mock=delta_ref)
#ifdef _GET_DIST_MIN_SEP_MOCK_
    this->tracer.get_distribution_min_separations(this->ncells_info);
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
#ifdef _USE_GNUPLOT_SREL_
      this->plot_scaling_relation_assignment(_MASS_);
#endif
      this->get_mean_scaling_relation_assignment_bias(_VMAX_);
      this->get_mean_scaling_relation_assignment_bias(_MASS_);
#endif
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
     this->assign_tracer_property(false, _VMAX_);  // Drawn Vmax from P(Vmax|M_mock=Mref, delta_mock=delta_ref)
#ifdef _USE_GNUPLOT_SREL_
     this->plot_scaling_relation_assignment(_VMAX_);
#endif
#ifdef _ASSIGN_PROPERTIES_TO_REFERENCE_
    this->plot_mean_scaling_relation_assignment_bias(_MASS_);
#endif
#endif
#endif
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// Now we shall assign propertie sin the order Rs, Cvir, Spin. Spin is the last to be assigned
// But if spin is not used, we freee memmory at Cvir
//------------RS-----------------------------
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _FULL_VERBOSE_
     So.message_screen("******************************************************");
     So.message_screen("Assigning Rs using P(Rs|VMAX, Mvir):");
     So.message_screen("******************************************************");
#endif
    this->params.set_i_mass_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_vmax_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_spin_bullock_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_rs_g(POSITIVE_INT);// this number only needs to be positive
    this->tracer.set_type_of_object("TRACER_MOCK"); // THis name is important.
    this->get_scaling_relations_secondary_property(_RS_); // Learn P(RS|M,V) from reference
    this->assign_tracer_property(false, _RS_);            // Assign property
    this->get_mean_scaling_relation_assignment_bias(_RS_);
    this->get_mean_scaling_relation_assignment_secondary_bias(_VMAX_,_RS_); // uses methods from Catalog to get mean primary bias
    this->get_mean_scaling_relation_assignment_secondary_bias(_MASS_,_RS_); // uses methods from Catalog to get mean primary bias
#ifndef _USE_SPIN_AS_DERIVED_OBSERVABLE_ && !defined _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
   if(this->tracer_ref.Halo.size()>0)// We kept tracer_ref in memory for comparisons of bias. Here we clar it and read  reference catalog below if we need it again
   {
     So.message_screen("**Freeing memmory from tracer_ref, line", __LINE__);
     this->tracer_ref.Halo.clear();
     this->tracer_ref.Halo.shrink_to_fit();
   }
#endif
#endif
//------------Cvir-----------------------------
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
#ifdef _FULL_VERBOSE_
     So.message_screen("******************************************************");
     So.message_screen("Assigning Cvir using P(Cvir|VMAX, Mvir):");
     So.message_screen("******************************************************");
#endif
    this->params.set_i_mass_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_vmax_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_spin_bullock_g(NEGATIVE_INT);// this number only needs to be positive
    this->params.set_i_rs_g(NEGATIVE_INT);// this number only needs to be positive
    this->tracer.set_type_of_object("TRACER_MOCK"); // THis name is important.
    this->get_scaling_relations_secondary_property(_CONCENTRATION_); // Learn P(RS|M,V) from reference
    this->assign_tracer_property(false, _CONCENTRATION_);            // Assign property
    this->get_mean_scaling_relation_assignment_bias(_CONCENTRATION_);
    this->get_mean_scaling_relation_assignment_secondary_bias(_VMAX_,_CONCENTRATION_); // uses methods from Catalog to get mean primary bias
    this->get_mean_scaling_relation_assignment_secondary_bias(_MASS_,_CONCENTRATION_); // uses methods from Catalog to get mean primary bias
#ifndef _USE_SPIN_AS_DERIVED_OBSERVABLE_
   if(this->tracer_ref.Halo.size()>0)// We kept tracer_ref in memory for comparisons of bias. Here we clar it and read  reference catalog below if we need it again
   {
     So.message_screen("**Freeing memmory from tracer_ref, line", __LINE__);
     this->tracer_ref.Halo.clear();
     this->tracer_ref.Halo.shrink_to_fit();
   }
#endif
#endif
//------------SPIN-----------------------------
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
#ifdef _FULL_VERBOSE_
  So.message_screen("******************************************************");
  So.message_screen("Assigning Spin using P(Spin|VMAX, Mvir):");
  So.message_screen("******************************************************");
#endif
  this->params.set_i_mass_g(NEGATIVE_INT);// this number only needs to be positive
  this->params.set_i_vmax_g(NEGATIVE_INT);// this number only needs to be positive
  this->params.set_i_rs_g(NEGATIVE_INT);// this number only needs to be positive
  this->params.set_i_spin_bullock_g(POSITIVE_INT);// this number only needs to be positive
  this->tracer.set_type_of_object("TRACER_MOCK"); // THis name is important.
  this->get_scaling_relations_secondary_property(_SPIN_BULLOCK_); // Learn P(SPIN|M,V)// Here we need to enlarge dependencies to environmental.
  this->assign_tracer_property(false, _SPIN_BULLOCK_);            // Assign property
  this->get_mean_scaling_relation_assignment_bias(_SPIN_BULLOCK_); // uses methods from Catalog to get mean primary bias
  this->get_mean_scaling_relation_assignment_secondary_bias(_VMAX_,_SPIN_BULLOCK_); // uses methods from Catalog to get mean primary bias
  this->get_mean_scaling_relation_assignment_secondary_bias(_MASS_,_SPIN_BULLOCK_); // uses methods from Catalog to get mean primary bias
  if(this->tracer_ref.Halo.size()>0)// We kept tracer_ref in memory for comparisons of bias. Here we clar it and read  reference catalog below if we need it again
   {
     So.message_screen("**Freeing memmory from tracer_ref, line", __LINE__);
     this->tracer_ref.Halo.clear();
     this->tracer_ref.Halo.shrink_to_fit();
   }
#endif
#endif // end if ASSIGN_PROPERTY
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------
// Collapse randoms:  this applies if collapse_randoms_aux is *undef*
#ifdef _COLLAPSE_RANDOMS_
     this->collapse_randoms();
#endif
#ifdef _USE_LPT_
     this->lpt.set_fnameTRACERCAT(outputFileName);
#else
    So.message_warning("Check file name in line ", __LINE__);
#endif
// ******** here we correct for mean vels and boundary conditions
#if !defined _ASSIGN_PROPERTIES_TO_REFERENCE_ || !defined _ASSIGN_PROPERTIES_TO_NEW_REFERENCE_  // THESE TWO CASES DO NOT NEED bc, 
#ifdef _APPLY_PERIODIC_BC_
#ifdef _FULL_VERBOSE_
    this->So.message_screen("Reinforcing boundary conditions");
#endif
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
    {
      if(this->tracer.Halo[i].coord1>this->params._Lbox())
        this->tracer.Halo[i].coord1-=this->params._Lbox();
      if(this->tracer.Halo[i].coord2>this->params._Lbox())
        this->tracer.Halo[i].coord2-=this->params._Lbox();
      if(this->tracer.Halo[i].coord3>this->params._Lbox())
        this->tracer.Halo[i].coord3-=this->params._Lbox();
      if(this->tracer.Halo[i].coord1<0)
        this->tracer.Halo[i].coord1+=this->params._Lbox();
      if(this->tracer.Halo[i].coord2<0)
        this->tracer.Halo[i].coord2+=this->params._Lbox();
      if(this->tracer.Halo[i].coord3<0)
        this->tracer.Halo[i].coord3+=this->params._Lbox();
    }
    this->So.DONE();
#endif
#endif
#ifndef _FCOL_TESTS_
#ifdef _CORRECT_MEAN_VELOCITIES_
#ifdef _FULL_VERBOSE_
    this->So.message_screen("Correcting velocities to zero mean");
#endif
    real_prec meanvel=0;
#pragma omp parallel for reduction(+:meanvel)
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        meanvel+=this->tracer.Halo[i].vel1;
    meanvel/=static_cast<real_prec>(this->tracer._NOBJS());
#pragma omp parallel for
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        this->tracer.Halo[i].vel1-=meanvel;
    meanvel=0;
#pragma omp parallel for reduction(+:meanvel)
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        meanvel+=this->tracer.Halo[i].vel2;
    meanvel/=static_cast<real_prec>(this->tracer._NOBJS());
#pragma omp parallel for
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        this->tracer.Halo[i].vel2-=meanvel;
    meanvel=0;
#pragma omp parallel for reduction(+:meanvel)
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        meanvel+=this->tracer.Halo[i].vel3;
    meanvel/=static_cast<real_prec>(this->tracer._NOBJS());
#pragma omp parallel for
    for(ULONG i = 0; i< this->tracer._NOBJS(); ++i)
        this->tracer.Halo[i].vel3-=meanvel;
    this->So.message_screen("Done");
#endif // end for fcol_tests_
#endif // end for correct_mean_velocities
 // ------------------------------------------------------------------------
// Herew we get the power spectrum from the catalog, before actually writting it as an output.
#ifdef _GET_POWER_FROM_CATS_
    Params params_aux=this->params;
    this->tracer.set_type_of_object("TRACER_MOCK"); // THis name is important.
#ifdef _FCOL_TESTS_
    ULONG ifcol=static_cast<ULONG>(floor(1000.*this->params._Distance_fraction()));
    params_aux.set_Name_survey("BAM_HALOS_fcol"+to_string(ifcol));
#else
#ifdef _MULTISCALE_
    //      params_aux.set_Name_survey("BAM_HALOS_ML"+to_string(this->params._Number_of_MultiLevels()));
    params_aux.set_Name_survey("BAM_HALOS");
#else
    params_aux.set_Name_survey("BAM_HALOS");
#endif
#endif // end _FCOL_TESTS_

#ifdef _CORRECT_VEL_DENSITY_  // This contributes to get teh right velocity distribution and Large-scale multipole spectra in redshift space
 vector<real_prec>corr_den;
//#endif
    corr_den.push_back(_VEL_CORR_EXP_KNOTS_);
    corr_den.push_back(_VEL_CORR_EXP_FILAMENTS_);
    corr_den.push_back(_VEL_CORR_EXP_SHEETS_);
    corr_den.push_back(_VEL_CORR_EXP_VOIDS_);
    real_prec ibvel= 1.0;// 1.03;
    
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->tracer._NOBJS();++i)
      {
#ifdef _CORRECT_VEL_HALOMASS_
	if (this->tracer.Halo[i].mass >1e12  && this->tracer.Halo[i].mass <1e13)
	  ibvel=0.9;
	else if (this->tracer.Halo[i].mass >=1e13)
	  ibvel=0.83;
#endif //* end of _CORRECT_VEL_HALOMASS_
	int cwc_exponent=this->tracer.Halo[i].gal_cwt;
	real_prec density=pow(abs(1.00+this->tracer.Halo[i].local_dm),corr_den[cwc_exponent]);
	this->tracer.Halo[i].vel1*=density*ibvel;
	this->tracer.Halo[i].vel2*=density*ibvel;
	this->tracer.Halo[i].vel3*=density*ibvel;
      }
#endif //* end of _CORRECT_VEL_DENSITY_
    //* If no properties are assigned, we might still want to use the cwt of each tracer. THis is done in get_scaling_relations_primary_property and _assign_property_new_new *//
    //* but these functions are called only if assignment is to be done. So, if that is not the case, we do it here*//
#ifndef _ASSIGN_PROPERTY_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->tracer._NOBJS();++i)
        this->tracer.Halo[i].gal_cwt=this->cwclass.CWClass[this->tracer.Halo[i].GridID];
#endif
#ifdef _ADD_RANDOM_VEL_
    const gsl_rng_type * rng_p;
    gsl_rng * gBaseRand_p;
    rng_p = gsl_rng_mt19937;//_default;
    gsl_rng_env_setup();
    gsl_rng_default_seed=125*this->params._realization();
    gBaseRand_p = gsl_rng_alloc (rng_p);
    for(ULONG i=0;i<this->tracer._NOBJS();++i)
      {
	if(this->tracer.Halo[i].gal_cwt==I_KNOT)
	  {
	    //	    real_prec epsilon=
	    // This otion adds scaltter to individual vel components
	    this->tracer.Halo[i].vel1+=gsl_ran_gaussian(gBaseRand_p,this->params._slengthv());
	    this->tracer.Halo[i].vel2+=gsl_ran_gaussian(gBaseRand_p,this->params._slengthv());
	    this->tracer.Halo[i].vel3+=gsl_ran_gaussian(gBaseRand_p,this->params._slengthv());
	    // this option adds scatte to the mvelocity without changing the direction
	    /*
	      real_prec vmod=sqrt(pow(this->tracer.Halo[i].vel1,2)+pow(this->tracer.Halo[i].vel2,2)+pow(this->tracer.Halo[i].vel3,2));
	      real_prec vmod_ep=(vmod+epsilon)/static_cast<real_prec>(vmod);
	      this->tracer.Halo[i].vel1*=vmod_ep;
	      this->tracer.Halo[i].vel2*=vmod_ep;
	      this->tracer.Halo[i].vel3*=vmod_ep;
	    */
	  }
      }
#endif
    // Now get power:
#ifdef _GET_POWER_BMT_
    params_aux.set_Name_survey("BTM_CAL");
    params_aux.set_mass_assignment_scheme("TSC");
    params_aux.set_MAS_correction(true);
    params_aux.set_SN_correction(true);
    params_aux.set_vel_units_g("kmps");
    params_aux.set_Nft(400);
    params_aux.set_input_type("catalog");
    params_aux.derived_pars(); // Since we have new parmeter we ned to compute the derived params again.
#ifdef _COMPUTE_REF_POWER_
    this->params_original.set_Name_survey("REF_GALS");
    this->params_original.set_mass_assignment_scheme("TSC");
    this->params_original.set_MAS_correction(true);
    this->params_original.set_SN_correction(true);
    this->params_original.set_Nft(400);
    this->params_original.set_input_type("catalog");
    this->params_original.derived_pars();
#endif //* end for _COMPUTE_REF_POWER_ *//
    PowerSpectrumF Power_BiasMT(params_aux);
    PowerSpectrumF Power_BiasMT_ref(this->params_original);
#ifdef _COMPUTE_REF_POWER_
    Catalog atracer_ref;
#ifdef _SLICS_
      if(this->params._IC_index()<1000)
        {
#endif
        atracer_ref.set_params(this->params_original);
        atracer_ref.read_catalog(this->params_original._Input_dir_cat()+this->params_original._file_catalogue(),pow(10,this->params._LOGMASSmin())*this->params._MASS_units());
#ifdef _CAT_POWER_ONLY_COORDS_
        Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"real","_NONE_"); //*Get power*//
#else
        Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"real",_MASS_); //* Get the reference power in mass bins*//
#endif
#ifdef _SLICS_
	}
#endif
#endif //* end for _COMPUTE_REF_POWER_ *//
#ifdef _CAT_POWER_ONLY_COORDS_
    Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"real","_NONE_"); //*Get power*//
#else
    Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"real",_MASS_); //* Get the REFERENCE power in mass bins*//
#endif
#ifdef _USE_GNUPLOT_
#ifdef _COMPUTE_REF_POWER_
#ifdef _SLICS_
      if(this->params._IC_index()<1000)
#endif  //* closes _SLICS_
        this->plot_power(Power_BiasMT_ref.power_in_bins,Power_BiasMT.power_in_bins,_MASS_,"real");
#else
      this->plot_power(Power_BiasMT.power_in_bins,"real","_NONE_");
#endif //* closes _COMPUTE_REF_POWER_  *//
#endif //* Closes  _USE_GNUPLOT_  *//
#ifndef _CAT_POWER_ONLY_COORDS_
#ifdef _COMPUTE_REF_POWER_
      Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"real",_VMAX_); //* Get the reference power in vmax bins*//
#endif
      Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"real",_VMAX_); //* Get the bam power in vmax bins*//
#ifdef _USE_GNUPLOT_
      this->plot_power(Power_BiasMT_ref.power_in_bins,Power_BiasMT.power_in_bins,_VMAX_,"real");
#endif
#endif
#ifdef _REDSHIFT_SPACE_
#ifndef _CAT_POWER_ONLY_COORDS_
      // *******************
      Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"redshift",_MASS_);
#ifdef _COMPUTE_REF_POWER_
      Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"redshift",_MASS_);
#endif
#ifdef _USE_GNUPLOT_
      this->plot_power(Power_BiasMT_ref.power_in_bins,Power_BiasMT.power_in_bins, "redshift",_MASS_);
#endif

      Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"redshift",_VMAX_);
#ifdef _COMPUTE_REF_POWER_
      Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"redshift",_VMAX_);
#endif
#ifdef _USE_GNUPLOT_
      this->plot_power(Power_BiasMT_ref.power_in_bins,Power_BiasMT.power_in_bins, "redshift",_VMAX_);
#endif

#else  // else for ifndef _CAT_POWER_ONLY_COORDS_
      Power_BiasMT.compute_power_spectrum(this->tracer.Halo,"redshift","_NONE_");
#ifdef _COMPUTE_REF_POWER_
      Power_BiasMT_ref.compute_power_spectrum(atracer_ref.Halo,"redshift","_NONE_");
#endif

#ifdef _USE_GNUPLOT_
      this->plot_power(Power_BiasMT_ref.power_in_bins,Power_BiasMT.power_in_bins,"Redshift space","_NONE_");
#endif


#endif //* end for _CAT_POWER_ONLY_COORDS_*//

#endif //* end for _REDSHIFT_SPACE_ *//
     So.message_screen("PowerSpectrum has been measured.");
#endif //* End for _GET_POWER_FROM_CATS_
#endif //* End for _GET_POWER_BMT
#if defined _ASSIGN_PROPERTIES_TO_MOCK_  || defined (_WRITE_COORDINATES_)
#ifdef _WRITE_BiasMT_CATALOGS_
  this->tracer.write_catalog(outputFileName.c_str());
#ifdef _VERBOSE_FREEMEM_
         So.message_screen("Freeing memmory from tracer");
#endif
     this->tracer.Halo.clear();
     this->tracer.Halo.shrink_to_fit();
     So.DONE();
#endif //* end for _WRITE_BiasMT_CATALOGS_
#endif //* end for defined _ASSIGN_PROPERTIES_TO_MOCK_  || defined (_WRITE_COORDINATES_)
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef MOCK_MODE
#ifdef _CORRECT_EXCLUSION_
void BiasMT::correct_for_exclusion(ULONG LENGHTdm){
  So.enter(__PRETTY_FUNCTION__);
  //A bin in DM properties have contributions from different sections of the volume (i.e, differnet IDGrids).
  // We study here the distribution of separations within objects in a THETA bin.
  // The idea is to move the masses in that particular theta_bin such that the distribution of the mock follow that of the reference, computed here.
  // An alternative idea (used in Hadron) is to  define a Mass threshold such that masses in a theta-bin above that mass threshold
  //are to be associated to distant cells.
  const gsl_rng_type *Tn;
  gsl_rng *rn ;
  Tn = gsl_rng_default;
  rn = gsl_rng_alloc (Tn);
  this->tracer.masses_in_cells_min_sep.clear();
  this->tracer.masses_in_cells_min_sep.shrink_to_fit();
  this->tracer.masses_in_cells_min_sep.resize(LENGHTdm);
  vector<int>number_in_theta_mock(LENGHTdm,0);
#pragma omp parallel for
  for(ULONG i=0;i<LENGHTdm; ++i)
    number_in_theta_mock[i]=this->dm_properties_bins_mock[i].tracer_properties.size();
  So.DONE();
  // Obtain the list of pairs of masses allocated in a bin of theta that are within a distance min_halo_sep, EXCLUSION_SCALE. output in  this->tracer.masses_in_cells_min_sep
  this->tracer.get_masses_of_pairs_in_min_separation_bin_in_theta_bin(this->min_halo_separation, this->dm_properties_bins_mock);
  for(int ih=0;ih< LENGHTdm;++ih) //Loop over the bins in Theta
    {
      real_prec M1_max_ref; // min and max of reference masses in theta_bin in the min_separation bin
      real_prec M2_max_ref;
      real_prec M1_max_aux=-1e5; // min and max of reference masses in theta_bin in the min_separation bin
      real_prec M2_max_aux=-1e5;
      for(int j=0; j< this->tracer_ref.masses_in_cells_min_sep[ih].M1.size();++j)
        {
          M1_max_ref=max(this->tracer_ref.masses_in_cells_min_sep[ih].M1[j], M1_max_aux);
          M1_max_aux=M1_max_ref;
          M2_max_ref=max(this->tracer_ref.masses_in_cells_min_sep[ih].M2[j], M1_max_aux);
          M2_max_aux=M2_max_ref;
        }
      // search for the masses M1 above the max in the bin. It is enough to look at one particle. These are the masses to be swaped
      for(int j=0; j< this->tracer.masses_in_cells_min_sep[ih].M1.size();++j)
        {
          real_prec M1=this->tracer.masses_in_cells_min_sep[ih].M1[j];
          if(M1>M1_max_ref) //allocate the galID of this mass
            this->tracer.masses_in_cells_min_sep[ih].mass_to_swap.push_back(this->dm_properties_bins_mock[ih].GalID_bin_properties[j]);
        }
    }
  So.message_screen("Swaping masses...");
  for(int ih=0;ih< LENGHTdm;++ih)//Loop over the bins in Theta
    {
      int N_props_in_bin = number_in_theta_mock[ih]; //Available Masses in bin of THETA, regardless of the separation between objects
      for(int j=0; j< this->tracer.masses_in_cells_min_sep[ih].mass_to_swap.size();++j)//loopover the masses to be swaped
        {
          ULONG old_ig=this->tracer.masses_in_cells_min_sep[ih].mass_to_swap[j]; // galaxy index IDg [0,NOBJS) of the mass requested to be swaped
          real_prec old_mass=this->tracer.Halo[old_ig].mass; //mass requested to be swaped
          bool baux=false;
          while(baux==false)
            {
              int i_mass_halo_label= gsl_rng_uniform_int(rn,N_props_in_bin); // pick-up randomly some other mass in the same bin-theta
              ULONG new_ig=this->dm_properties_bins_mock[ih].
                      [i_mass_halo_label];  // galaxy index IDg [0,NOBJS) of an randomly selected object
              real_prec new_mass=this->tracer.Halo[new_ig].mass; //mass of the randomly selected object
              //swap masses:
              if(new_ig!=this->tracer.masses_in_cells_min_sep[ih].mass_to_swap[j])
                {
                  this->tracer.Halo[old_ig].mass=new_mass; // assign mass
                  this->tracer.Halo[new_ig].mass=old_mass;
                  number_in_theta_mock[ih]--;
                  baux=true;
                }
            }
        }
    }
    So.DONE();
    number_in_theta_mock.clear();number_in_theta_mock.shrink_to_fit();
    this->dm_properties_bins_mock.clear();
    this->dm_properties_bins_mock.shrink_to_fit();
}
#endif
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::v_konvolve(vector<real_prec> &vel)
{
   // theta = nabla V. Getting V given theta is the main goal.
    // For norm=false, this function solves for the Displacement(vei) from the Divergencve of the Displacement (nabla delta) in order to obtain the component comp of the displacement.
    // For norm=true, this function solves for the Displacement from the Divergencve of the Displacement and multiply for
    // the factor fHa in order to get the component comp of te velocity filed.
    //
   So.enter(__PRETTY_FUNCTION__ );
   ULONG N= this->params._NGRID();
   real_prec Lzp1=this->params._Lbox();
   real_prec Lzp2=this->params._Lbox();
   real_prec Lzp3=this->params._Lbox();
   ULONG Nzp1=this->params._Nft();
   ULONG Nzp2=this->params._Nft();
   ULONG Nzp3=this->params._Nft();
   ULONG Nhalf=this->NTT;
#ifdef DOUBLE_PREC
   complex_prec *vel_f= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
   complex_prec *vel_f= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
   // Fourier transform of the nabla X
   do_fftw_r2c(this->params._Nft(),vel,vel_f);
  vector<real_prec> coords(this->params._Nft(),0);
  real_prec deltak=2.*M_PI/this->params._Lbox();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->params._Nft() ;++i)
    coords[i]=deltak*(i<=this->params._Nft()/2? static_cast<real_prec>(i): -static_cast<real_prec>(this->params._Nft()-i));
  real_prec kstar2=pow(this->params._slengthv(),2.);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
   for (ULONG i=0;i<this->params._Nft();i++)
     for (ULONG j=0;j<this->params._Nft();j++)
       for (ULONG k=0;k<this->params._Nft()/2+1;k++)
         {
           real_prec kmod2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);  // Get k**2
           ULONG index=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
           real_prec vel_kernel = pow(1+2.*kmod2*kstar2,this->params._vkernel_exponent());
//           real_prec vel_kernel = exp(-kmod2*kstar2);
           real_prec vreal= vel_f[index][REAL];
           real_prec vimag= vel_f[index][IMAG];
           vel_f[index][REAL]= vreal*vel_kernel;
           vel_f[index][IMAG]= vimag*vel_kernel;
         }
   do_fftw_c2r(this->params._Nft(),vel_f,vel);
#ifdef DOUBLE_PREC
   fftw_free(vel_f);
#else
   fftwf_free(vel_f);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_minimum_multiscale_property(){
  this->minimum_multiscale_property  =0.0;
  vector<real_prec>aux;
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
  for(int i=0;i<this->params._Number_of_MultiLevels();++i)
    aux.push_back(this->params.get_PropThreshold_MultiLevels(i)) ;
#endif
  this->minimum_multiscale_property=static_cast<real_prec>(get_min(aux));
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Minimum Multiscale property =", this->minimum_multiscale_property);
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::randomize_properties()
{
    gsl_rng_env_setup();
    const gsl_rng_type *Tn;
    gsl_rng_default_seed=1015;
    Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
    gsl_rng *rn = gsl_rng_alloc (Tn);
    So.message_screen("Randomizing container of properties in bins of theta");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG it=0; it <this->dm_properties_bins.size(); ++it)
    {
      int ll=this->dm_properties_bins[it].tracer_properties.size();
      if(ll>0)
       {
        vector<ULONG>aux_l;
        vector<real_prec>aux_p;
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
        vector<real_prec>aux_p_s;
#endif
        vector<bool>aux_u;
        for(int im=0; im< ll; ++im )
         {
            aux_l.push_back(im);
            aux_p.push_back(this->dm_properties_bins[it].tracer_properties[im]);
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
            aux_p_s.push_back(this->dm_properties_bins[it].tracer_properties_secondary[im]);
#endif
            aux_u.push_back(this->dm_properties_bins[it].used_property[im]);
         }
         if(aux_l.size()>0)
           {
              gsl_ran_shuffle(rn,&aux_l[0], aux_l.size() ,sizeof(ULONG));  // randomize within the level 4
              for(int im=0; im< aux_l.size(); ++im )
                {
                  this->dm_properties_bins[it].tracer_properties[im]=aux_p[aux_l[im]];
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
                  this->dm_properties_bins[it].tracer_properties_secondary[im]=aux_p_s[aux_l[im]];
#endif
                  this->dm_properties_bins[it].used_property[im]=aux_u[aux_l[im]];
                }
            }
        }
    }
So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::sort_properties(){
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   So.message_screen("Sorting masses in bins of theta");
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
   So.message_screen("Sorting (top-bottom) the values of Vmax in each bin of theta ");
#endif
   So.message_screen("used in the multiscale approach: ");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG it=0; it < this->dm_properties_bins.size(); ++it)  // We now sort all properties in the different theta bins from top-to-bottom
     {
       ULONG ll=this->dm_properties_bins[it].tracer_properties.size();
       if(ll>0)
        {
          gsl_vector *aux_properties=gsl_vector_alloc(ll);
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
           gsl_vector *aux_properties_secondary=gsl_vector_alloc(ll);
#endif
#if defined (_ASSIGN_MASS_LINKED_TO_V_FROM_REF_) || defined _NEW_APPROACH_ASS_
           gsl_vector *aux_properties_index=gsl_vector_alloc(ll);
#endif
           for(int im=0; im< ll; ++im )
           {
             gsl_vector_set(aux_properties,im,this->dm_properties_bins[it].tracer_properties[im]);
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
             gsl_vector_set(aux_properties_secondary,im,this->dm_properties_bins[it].tracer_properties_secondary[im]);
#endif
#if defined (_ASSIGN_MASS_LINKED_TO_V_FROM_REF_) || defined _NEW_APPROACH_ASS_
             gsl_vector_set(aux_properties_index,im,this->dm_properties_bins[it].index_reference[im]);
#endif
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
             gsl_vector_set(aux_properties_secondary,im,this->dm_properties_bins[it].tracer_properties_secondary[im]);
#endif
           }
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
           gsl_sort_vector2(aux_properties,aux_properties_secondary);
#else
           gsl_sort_vector(aux_properties);
#endif
           for(int im_new=0; im_new< ll; ++im_new )
            {
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
               this->dm_properties_bins[it].tracer_properties[im_new]=gsl_vector_get(aux_properties,ll-1-im_new);//top-bottom
               this->dm_properties_bins[it].tracer_properties_secondary[im_new]=gsl_vector_get(aux_properties_secondary,ll-1-im_new);
#else   // else for if def ASSIGN_MASS_LINKED_TO_V_FROM_REF
               this->dm_properties_bins[it].tracer_properties[im_new]=gsl_vector_get(aux_properties,ll-1-im_new);//top-bottom
#endif
#if defined (_ASSIGN_MASS_LINKED_TO_V_FROM_REF_)  || defined _NEW_APPROACH_ASS_
               this->dm_properties_bins[it].index_reference[im_new]=gsl_vector_get(aux_properties_index,ll-1-im_new);
#endif
           }
           gsl_vector_free(aux_properties);
#ifdef _ASSIGN_MASS_LINKED_TO_V_FROM_REF_
            gsl_vector_free(aux_properties_secondary);
#endif
#ifdef  _NEW_APPROACH_ASS_
            gsl_vector_free(aux_properties_index);
#endif
       }
     }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::sort_properties(int trial){
#ifdef _FULL_VERBOSE_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   So.message_screen("Sorting (top-bottom) masses in bins of {theta}_dm in dm container");
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
   So.message_screen("Sorting (top-bottom) the values of Vmax in each bin of {theta}_dm in dm container");
#endif
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG it=0; it < this->dm_properties_bins_mock.size(); ++it)  // We now sort all properties in the different theta bins from top-to-bottom
     {
       ULONG ll=this->dm_properties_bins_mock[it].tracer_properties.size();
       if(ll>0)
        {
          gsl_vector *aux_properties=gsl_vector_alloc(ll);
          gsl_vector *aux_properties_id=gsl_vector_alloc(ll);
           for(ULONG im=0; im< ll; ++im )
           {
             gsl_vector_set(aux_properties,im,this->dm_properties_bins_mock[it].tracer_properties[im]);
             gsl_vector_set(aux_properties_id,im,this->dm_properties_bins_mock[it].GalID_bin_properties[im]);
           }
           gsl_sort_vector2(aux_properties,aux_properties_id);
           for(ULONG im_new=0; im_new< ll; ++im_new )
            {
              this->dm_properties_bins_mock[it].tracer_properties[im_new]=gsl_vector_get(aux_properties,ll-1-im_new);//top-bottom
              this->dm_properties_bins_mock[it].GalID_bin_properties[im_new]=gsl_vector_get(aux_properties_id,ll-1-im_new);
           }
           gsl_vector_free(aux_properties);
           gsl_vector_free(aux_properties_id);
       }
     }
   So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::swap_mock_properties_in_theta_bins(int asy){
    // sort properties in the mock dm_properties structure
#ifdef _FULL_VERBOSE_
  So.message_screen("Swaping properties in mock:");
#endif
  this->sort_properties(asy);
  for(ULONG ih=0;ih<this->dm_properties_bins_mock.size();++ih)
       {
         ULONG size_dmp=this->dm_properties_bins_mock[ih].tracer_properties.size();
         if(size_dmp>1){
         int extra_n=0;
         if(size_dmp%2 !=0)extra_n=1;
         ULONG nswaps= static_cast<ULONG>(floor(size_dmp/2));
         for (ULONG is=0;is<nswaps;++is){ // loop over the possible swaps in this theta bin
             real_prec prop_max=this->dm_properties_bins_mock[ih].tracer_properties[is]; // this is thye largest vmax value in this bin
             ULONG ig_max=this->dm_properties_bins_mock[ih].GalID_bin_properties[is]; // this is the ID of the tracer with largest vmax value in this bin
             ULONG is_low=nswaps*2-is+extra_n-1;
             real_prec prop_min=this->dm_properties_bins_mock[ih].tracer_properties[is_low]; // this is the lowest vmax value in this bin
             ULONG ig_min=this->dm_properties_bins_mock[ih].GalID_bin_properties[is_low]; // this is the ID of the tracer with lowest vmax value in this bin
             //swap:
             this->tracer.Halo[ig_min].vmax=prop_max;
             this->tracer.Halo[ig_max].vmax=prop_min;
          }
        }
    }
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::extrapolate_kernel(vector<real_prec>&kernel_old, vector<real_prec>&kernel_new){
#ifdef _FULL_VERBOSE_
  So.message_screen("Extrapolating kernel");
#endif
  ULONG Nft=this->params._Nft();
  ULONG Nft_low=this->params._Nft_low();
  vector<real_prec> coords_low(Nft_low,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<coords_low.size() ;++i)
    coords_low[i]= (i<=Nft_low/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft_low-i));
  vector<int>nmodes_kk(Nft_low/2/this->params._ndel_data(), 0);
  vector<real_prec>kernel_1d(Nft_low/2/this->params._ndel_data(), 0);
  for(ULONG i=0;i< Nft_low; ++i)
    {
      for(ULONG j=0;j< Nft_low; ++j)
        {
          for(ULONG k=0;k< Nft_low/2+1; ++k)
            {
              ULONG ind=index_3d(i,j,k,Nft_low, Nft_low/2+1);
              real_prec kv=_get_modulo(coords_low[i]*this->params._d_deltak_x_low(),coords_low[j]*this->params._d_deltak_y_low(),coords_low[k]*this->params._d_deltak_z_low());
              ULONG kmod=static_cast<int>(floor((kv-this->params._d_kmin_low())/this->params._d_DeltaK_data_low()));
              if(kmod<nmodes_kk.size())
                {
                  kernel_1d[kmod]+=kernel_old[ind]; //weight is 1 of no improvement; the new kernel if it gets closer
                  nmodes_kk[kmod]++;
                }
            }
        }
    }
  for(ULONG i=0;i<kernel_1d.size(); ++i)
    kernel_1d[i]/=static_cast<real_prec>(nmodes_kk[i]);
  //extrapolate to a new kerneĺ
  vector<real_prec>kernel_1d_new(Nft/2/this->params._ndel_data(), 0);
  real_prec av_kernel=kernel_1d[0];
  for(ULONG i=0;i< kernel_1d_new.size(); ++i)
    {
      ULONG jlim=Nft/2-Nft_low/2;
      kernel_1d_new[i]=  i< jlim? av_kernel: kernel_1d[i-jlim];
    }
  kernel_1d.clear();kernel_1d.shrink_to_fit();
  coords_low.clear(); coords_low.shrink_to_fit();
  ofstream kal; kal.open("kernel_extra1500.txt");
  for(ULONG i=0;i< kernel_1d_new.size(); ++i)
    kal<<this->kvec[i]<<"   "<<kernel_1d_new[i]<<endl;
  kal.close();
  vector<real_prec> coords(Nft,0);
  for(ULONG i=0;i<coords.size() ;++i)
    coords[i]= (i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));
  for(ULONG i=0;i< Nft; ++i)
    for(ULONG j=0;j< Nft; ++j)
      for(ULONG k=0;k< Nft/2+1; ++k)
        {
          ULONG ind=index_3d(i,j,k, Nft, Nft/2+1);
          real_prec kv=_get_modulo(coords[i]*this->params._d_deltak_x(),coords[j]*this->params._d_deltak_y(),coords[k]*this->params._d_deltak_z());
          ULONG kmod=static_cast<ULONG>(floor( (kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
          if(kmod<kernel_1d_new.size())
            kernel_new[ind]=kernel_1d_new[kmod];
        }
  kernel_1d_new.clear();kernel_1d_new.shrink_to_fit();
  coords.clear(); coords.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_IC_from_particles()
{
  ULONG n=this->params._N_lines_binary();
#ifdef _FULL_VERBOSE_
  So.message_screen("Reading binary files for pos of DM, and interpolating into a grid.");
#endif
  cout<<this->params._Nft_HR()<<endl;
  cout<<this->params._Nft()<<endl;
  cout<<this->params._Lbox()<<endl;
  cout<<this->params._d2_HR()<<endl;
  cout<<this->params._masskernel()<<endl;
  vector<real_prec>xg(this->params._N_lines_binary(),0);
  this->File.read_array(this->params._file_bin_x_coord(), xg);
  vector<real_prec>yg(this->params._N_lines_binary(),0);
  this->File.read_array(this->params._file_bin_y_coord(), yg);
  vector<real_prec>zg(this->params._N_lines_binary(),0);
  this->File.read_array(this->params._file_bin_z_coord(), zg);
  string  file=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft())+"_MAS"+to_string(this->params._masskernel())+"_"+this->params._Name_survey();
  if(true==this->params._use_low_pass_filter())
  {
    string  file=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft())+"_MAS"+to_string(this->params._masskernel())+"_"+this->params._Name_survey();
    string  file_hr=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft_HR())+"_MAS"+to_string(this->params._masskernel())+"_"+this->params._Name_survey();
    vector<real_prec>field_hr(this->params._NGRID_HR() ,0);
    string mass;
    if(this->params._masskernel()==0)
    {
        getDensity_NGP(this->params._Nft_HR(),this->params._Nft_HR(),this->params._Nft_HR(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(), this->params._d1_HR(), this->params._d2_HR(), this->params._d3_HR(),0,0,0,xg,yg,zg,zg,field_hr,false);
        mass="NGP";
    }
    else if(this->params._masskernel()==1){
        getDensity_CIC(this->params._Nft_HR(),this->params._Nft_HR(),this->params._Nft_HR(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(), this->params._d1_HR(), this->params._d2_HR(), this->params._d3_HR(),0,0,0,xg,yg,zg,zg,field_hr,false);
        mass="CIC";
    }
    if(this->params._masskernel()==2){
       getDensity_TSC(this->params._Nft_HR(),this->params._Nft_HR(),this->params._Nft_HR(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(), this->params._d1_HR(), this->params._d2_HR(), this->params._d3_HR(),0,0,0,xg,yg,zg,zg,field_hr,false);
        mass="TSC";
   }
    xg.clear(); xg.shrink_to_fit();
    yg.clear(); yg.shrink_to_fit();
    zg.clear(); zg.shrink_to_fit();
    real_prec mean_field_hr=get_mean(field_hr);
    this->File.write_array(file_hr, field_hr);
    vector<real_prec>field(this->params._NGRID() ,0);
    So.message_screen("Appliying low-pass filter");
    low_pass_filter(this->params._Nft_HR(), this->params._Nft(),this->params._masskernel(),false, field_hr, field,this->params._Lbox());
    real_prec mean_field_lr=get_mean(field);
    So.message_screen("Mean from HR =", mean_field_hr);
    So.message_screen("Mean from LR =", mean_field_lr);
#ifdef _USE_GNUPLOT_
        Params params_aux=this->params; // copy the params
        this->params.set_Nft(this->params._Nft_HR()); //change Nft by Nft_hr in params for PowerP measurement
        this->params.derived_pars();
        this->params.set_mass_assignment_scheme(mass); //change Nft by Nft_hr in params for PowerP measurement
        this->params.set_MAS_correction(true); //change Nft by Nft_hr in params for PowerP measurement
        this->params.set_SN_correction(false); //change Nft by Nft_hr in params for PowerP measurement
        this->params.set_input_type("density_grid"); //change Nft by Nft_hr in params for PowerP measurement
        PowerSpectrumF cPSFa(this->params);
        cPSFa.compute_power_spectrum_grid(field_hr,true);
        cPSFa.write_power_and_modes();
        std::vector<std::pair<double, double> > xy_ptr;
        for(ULONG i=1; i<cPSFa._kvector_data_size(); ++i)
            xy_ptr.push_back(std::make_pair(cPSFa._kvector_data(i), cPSFa._pk0(i)));
        params_aux.set_Name_survey(this->params._Name_survey()+"_HR");
        params_aux.set_mass_assignment_scheme(mass); //change Nft by Nft_hr in params for PowerP measurement
        params_aux.set_MAS_correction(false); //change Nft by Nft_hr in params for PowerP measurement
        params_aux.set_SN_correction(false); //change Nft by Nft_hr in params for PowerP measurement
        params_aux.set_input_type("density_grid"); //change Nft by Nft_hr in params for PowerP measurement
        PowerSpectrumF cPSFb(params_aux);
        cPSFb.compute_power_spectrum_grid(field,true);
        cPSFb.write_power_and_modes();
        std::vector<std::pair<double, double> > xy_pts;
        for(ULONG i=1; i<cPSFb._kvector_data_size(); ++i)
            xy_pts.push_back(std::make_pair(cPSFb._kvector_data(i), cPSFb._pk0(i)));
        cout<<cPSFa._kvector_data_size()<<"   "<<cPSFb._kvector_data_size()<<endl;
        this->gp_power<<"set border linewidth 2.\n";
        this->gp_power<<"set size square \n";
        this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
        this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
        this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
        this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
        this->gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
        this->gp_power<<"set log x \n";
        this->gp_power<<"set log y \n";
        this->gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,14'\n";
        this->gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,14'\n";
        this->gp_power<<"plot [0.001:2.0]"<<this->gp_power.file1d(xy_ptr) << "w l lw 3.3 lt 2 title 'HighRes' \n";
        this->gp_power<<"replot "<<this->gp_power.file1d(xy_pts) << " w l lw 3.3 lt 6 title 'LowRes"<<endl;
#endif
        get_overdens(field,mean_field_hr, field); // transform to overdensities
        this->File.write_array(file, field);
        field.clear(); field.shrink_to_fit();
        field_hr.clear(); field_hr.shrink_to_fit();
    }
    else
    {
     vector<real_prec>field(this->params._NGRID() ,0);
     string  file=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft())+"_MAS"+to_string(this->params._masskernel())+"_nc_"+this->params._Name_survey();
     getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(), this->params._d1(), this->params._d2(), this->params._d3(),0,0,0,xg,yg,zg,zg,field,false);
     xg.clear(); xg.shrink_to_fit();
     yg.clear(); yg.shrink_to_fit();
     zg.clear(); zg.shrink_to_fit();
     this->File.write_array(file, field);
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_GNUPLOT_
void BiasMT::plot_power(vector<s_info_in_bins> power_in_bins_ref, vector<s_info_in_bins> power_in_bins, string prop, string space)
{
  this->So.enter(__PRETTY_FUNCTION__);
  Gnuplot gp_power;
  gp_power<<"set style function lines\n";
  gp_power<<"set border linewidth 2.0\n";
  gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
  gp_power<<"set title 'Power Spectrum in "<<space<<"' \n";
  gp_power<<"set size square\n";
  gp_power<<"set log x \n";
  gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,15'  textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,15'  textcolor rgb '"<<FG_COLOR<<"' \n";
  int nbins=1;
  if("_NONE_"!=prop)
    nbins= prop == _MASS_ ? this->params._NMASSbins_power() : this->params._NVMAXbins_power();
  if("real"==space)
    {
      for(int im_bins=0; im_bins <nbins ;++im_bins) // start from 1 to avoid weird mode
	{
	  std::vector<std::pair<double, double> > xy_pts;
	  for(ULONG i=1; i< power_in_bins[im_bins].vbin.size(); ++i)
	    xy_pts.push_back(std::make_pair(power_in_bins[im_bins].vbin[i], log10(power_in_bins[im_bins].vq1[i])));
	  std::vector<std::pair<double, double> > xy_pts_ref;
	  for(ULONG i=1; i< power_in_bins_ref[im_bins].vbin.size(); ++i)
	    xy_pts_ref.push_back(std::make_pair(power_in_bins_ref[im_bins].vbin[i], log10(power_in_bins_ref[im_bins].vq1[i])));
	  if("_NONE_"==prop)
        gp_power<<"plot[:1.0][1.5:6.0]"<<gp_power.file1d(xy_pts) << " w l lw 3 lt 1 title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt 1 title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
	  else
	    {
	      if(0==im_bins)
        gp_power<<"plot[:1.0][1.5:6.0]"<<gp_power.file1d(xy_pts) <<" w l lw 3 lt 2 title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt 2 title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
	      else
        gp_power<<"replot"<<gp_power.file1d(xy_pts) << "w l lw 3 lt 3 title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt 2 title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
	    }
	  xy_pts.clear(); xy_pts.shrink_to_fit();
	  xy_pts_ref.clear(); xy_pts_ref.shrink_to_fit();
	}
    }
  else
    {
      for(int im_bins=0; im_bins <nbins ;++im_bins) // start from 1 to avoid weird mode
	{
	  std::vector<std::pair<double, double> > xy_pts;
	  for(ULONG i=1; i< power_in_bins[im_bins].vbin.size(); ++i)
	    xy_pts.push_back(std::make_pair(power_in_bins[im_bins].vbin[i], log10(power_in_bins[im_bins].vq1[i])));
	  std::vector<std::pair<double, double> > xy_pts2;
	  for(ULONG i=1; i< power_in_bins[im_bins].vbin.size(); ++i)
	    xy_pts2.push_back(std::make_pair(power_in_bins[im_bins].vbin[i], log10(power_in_bins[im_bins].vq2[i])));
	  std::vector<std::pair<double, double> > xy_pts4;
	  for(ULONG i=1; i< power_in_bins[im_bins].vbin.size(); ++i)
	    xy_pts4.push_back(std::make_pair(power_in_bins[im_bins].vbin[i], log10(fabs(power_in_bins[im_bins].vq3[i]))));
	  std::vector<std::pair<double, double> > xy_pts_ref;
	  for(ULONG i=1; i< power_in_bins_ref[im_bins].vbin.size(); ++i)
	    xy_pts_ref.push_back(std::make_pair(power_in_bins_ref[im_bins].vbin[i], log10(power_in_bins_ref[im_bins].vq1[i])));
      std::vector<std::pair<double, double> > xy_pts_ref2;
	  for(ULONG i=1; i< power_in_bins_ref[im_bins].vbin.size(); ++i)
	    xy_pts_ref2.push_back(std::make_pair(power_in_bins_ref[im_bins].vbin[i], log10(power_in_bins_ref[im_bins].vq2[i])));
	  std::vector<std::pair<double, double> > xy_pts_ref4;
	  for(ULONG i=1; i< power_in_bins_ref[im_bins].vbin.size(); ++i)
	    xy_pts_ref4.push_back(std::make_pair(power_in_bins_ref[im_bins].vbin[i], log10(fabs(power_in_bins_ref[im_bins].vq3[i]))));
	  if("_NONE_"!=prop)
	    {
	      if(im_bins==0)
		{
          gp_power<<"plot[:1.1][1.5:5.0]"<<gp_power.file1d(xy_pts) << " w l lw 3 lt 1 title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt 1 title ' Ref, "<<prop<<" bin "<<im_bins<<"' \n";
          gp_power<<"replot"<<gp_power.file1d(xy_pts2) <<" w l lw 3 lt 1 title 'Rec', "<<gp_power.file1d(xy_pts_ref2)<<" w l lw 1 lt 1 title ' Ref, "<<prop<<" bin "<<im_bins<<"' \n";
          gp_power<<"replot"<<gp_power.file1d(xy_pts4) <<" w l lw 3 lt 1 title 'Rec', "<<gp_power.file1d(xy_pts_ref4)<<" w l lw 1 lt 1 title ' Ref, "<<prop<<" bin "<<im_bins<<"' \n";
		}
	      else
		{
		  gp_power<<"plot[:1.1][1.5:5.0]"<<gp_power.file1d(xy_pts) << "w l lw 3 lt "<<im_bins<<" title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt "<<im_bins<<" title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
		  gp_power<<"replot"<<gp_power.file1d(xy_pts2) << "w l lw 3 lt "<<im_bins<<" title 'Rec', "<<gp_power.file1d(xy_pts_ref2)<<" w l lw 1 lt "<<im_bins<<" title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
		  gp_power<<"replot"<<gp_power.file1d(xy_pts4) << "w l lw 3 lt "<<im_bins<<" title 'Rec', "<<gp_power.file1d(xy_pts_ref4)<<" w l lw 1 lt "<<im_bins<<" title 'Ref, "<<prop<<" bin "<<im_bins<<"' \n";
		}
	    }
	  else
	    {
	      gp_power<<"plot[:1.1][1.5:5.0]"<<gp_power.file1d(xy_pts) << "w l lw 3 lt 1 title 'Rec', "<<gp_power.file1d(xy_pts_ref)<<" w l lw 1 lt 1 title ' Ref l=0' \n";
	      gp_power<<"replot"<<gp_power.file1d(xy_pts2) << "w l lw 3 lt 2 title 'Rec', "<<gp_power.file1d(xy_pts_ref2)<<" w l lw 1 lt 2 title ' Ref l=2' \n";
	      gp_power<<"replot"<<gp_power.file1d(xy_pts4) << "w l lw 3 lt 3 title 'Rec', "<<gp_power.file1d(xy_pts_ref4)<<" w l lw 1 lt 3 title ' Ref l=4' \n";
	    }
	  xy_pts.clear(); xy_pts.shrink_to_fit();
	  xy_pts_ref.clear(); xy_pts_ref.shrink_to_fit();
  	}
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::plot_power(vector<s_info_in_bins> power_in_bins, string prop, string space)
{
  this->So.enter(__PRETTY_FUNCTION__);
  Gnuplot gp_power;
  gp_power<<"set style function lines\n";
  gp_power<<"set border linewidth 2.0\n";
  gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    gp_power<<"set title 'Power Spectrum in "<<space<<"' \n";
  gp_power<<"set size square\n";
  gp_power<<"set log x \n";
  //gp_power<<"set log y \n";
  gp_power<<"set xlabel 'k [h / Mpc]' font 'Times-Roman,15'  textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_power<<"set ylabel 'log P(k) [(Mpc / h)³]'  font 'Times-Roman,15'  textcolor rgb '"<<FG_COLOR<<"' \n";
    int nbins=1;
  if("_NONE_"!=prop)
    nbins= prop == _MASS_ ? this->params._NMASSbins_power() : this->params._NVMAXbins_power();
  for(int im_bins=0; im_bins <nbins ;++im_bins)
    {
      std::vector<std::pair<double, double> > xy_pts;
      for(ULONG i=1; i< power_in_bins[im_bins].vbin.size(); ++i)
	xy_pts.push_back(std::make_pair(power_in_bins[im_bins].vbin[i], log10(power_in_bins[im_bins].vq1[i])));
      
      gp_power<<"plot[][1.5:5.0]"<<gp_power.file1d(xy_pts) << "w l lw 3 lt 3 title 'Power'"<<"' \n";
      xy_pts.clear(); xy_pts.shrink_to_fit();
    }
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_GNUPLOT_
void BiasMT::plot_scaling_relation_assignment(string prop){
  this->So.enter(__PRETTY_FUNCTION__);
  const gsl_rng_type * rng_t;
  gsl_rng * gBaseRando;
  gsl_rng_env_setup();
  gsl_rng_default_seed=35;
  rng_t = gsl_rng_mt19937;//_default;
  gBaseRando = gsl_rng_alloc (rng_t);
  real_prec fraction_tracers_to_file=0.06;
  ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction_tracers_to_file*this->tracer._NOBJS()));
  vector<pair<real_prec, real_prec> > xy_pts_ref;
  vector<pair<real_prec, real_prec> > xy_pts_ref_k;
  vector<pair<real_prec, real_prec> > xy_pts_ref_f;
  vector<pair<real_prec, real_prec> > xy_pts_ref_s;
  vector<pair<real_prec, real_prec> > xy_sr_ref;
  vector<pair<real_prec, real_prec> > xy_sr_rec;
  ULONG counter=0;
  Gnuplot gp_sr;
  gp_sr<<"set border linecolor '"<<FG_COLOR<<"' \n";
  gp_sr<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  gp_sr<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  gp_sr<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_sr<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
  gp_sr<<"set border linewidth 2.0\n";
  gp_sr<<"set size square \n";
  gp_sr<<"set log  \n";
  Gnuplot gp_sr_ref;
  gp_sr_ref<<"set border linecolor '"<<FG_COLOR<<"' \n";
  gp_sr_ref<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  gp_sr_ref<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  gp_sr_ref<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_sr_ref<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
  gp_sr_ref<<"set border linewidth 2.0\n";
  gp_sr_ref<<"set size square \n";
  gp_sr_ref<<"set log  \n";
  gp_sr_ref<<"set xlabel 'Reference' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  gp_sr_ref<<"set ylabel 'Assigned'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  if(prop==_MASS_)
    {
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
	  xy_sr_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer_ref.Halo[i].vmax));
	  xy_sr_rec.push_back(std::make_pair(this->tracer.Halo[i].mass, this->tracer.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_KNOT)
	    xy_pts_ref_k.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
	  counter++;
	}while(counter<Nobjs_fraction);

      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_SHEET)
	    xy_pts_ref_s.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
	  counter++;
	}while(counter<Nobjs_fraction);

      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_FILAMENT)
	    xy_pts_ref_f.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
	  counter++;
	}while(counter<Nobjs_fraction);    
    }
  else if(_VMAX_==prop)
    {
      /*               for(ULONG i=0;i<10000;++i)
		       xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
      */
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  // if(this->tracer.Halo[i].GridID==393708){
	  xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
	  //  }
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_KNOT)
	    xy_pts_ref_k.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_SHEET)
	    xy_pts_ref_s.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_FILAMENT)
	    xy_pts_ref_f.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
    }
  else if(_RS_==prop)
    {
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].rs, this->tracer.Halo[i].rs));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_KNOT)
	    xy_pts_ref_k.push_back(std::make_pair(this->tracer_ref.Halo[i].rs, this->tracer.Halo[i].rs));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_SHEET)
	    xy_pts_ref_s.push_back(std::make_pair(this->tracer_ref.Halo[i].rs, this->tracer.Halo[i].rs));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_FILAMENT)
	    xy_pts_ref_f.push_back(std::make_pair(this->tracer_ref.Halo[i].rs, this->tracer.Halo[i].rs));
	  counter++;
	}while(counter<Nobjs_fraction);

    }
  else if(_SPIN_==prop  )
    {
      do
    {
      ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
      xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].spin, this->tracer.Halo[i].spin));
      counter++;
    }while(counter<Nobjs_fraction);
  }  else if(_SPIN_BULLOCK_==prop  )
   {
      do
      {
      ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
      xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].spin_bullock, this->tracer.Halo[i].spin_bullock ));
      counter++;
    }while(counter<Nobjs_fraction);
   }
  else if("MASS_DENS"==prop)
    {
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].local_dm, this->tracer_ref.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_KNOT)
	    xy_pts_ref_k.push_back(std::make_pair(this->tracer_ref.Halo[i].local_dm, this->tracer_ref.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_SHEET)
	    xy_pts_ref_s.push_back(std::make_pair(this->tracer_ref.Halo[i].local_dm, this->tracer_ref.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
      counter=0;
      do
	{
	  ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
	  if(this->tracer_ref.Halo[i].gal_cwt==I_FILAMENT)
	    xy_pts_ref_f.push_back(std::make_pair(this->tracer_ref.Halo[i].local_dm, this->tracer_ref.Halo[i].vmax));
	  counter++;
	}while(counter<Nobjs_fraction);
    }
  So.message_screen("Number of objects shown in scatter plot =", counter,"objects");
  if(_VMAX_==prop)
    {
      gp_sr<<"set log"<<"\n";
      gp_sr<<"set xlabel 'Vmax reference [km/s]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gp_sr<<"set ylabel 'Vmax assigned [km/s]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gp_sr<<"plot [140:1000][140:1000] x, "<<gp_sr.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction "<<prop<<"'"<<endl;
      // gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_k)<<" w p ps 0.1 pt 7 title 'Knots'"<<endl;
      // gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_f)<<" w p ps 0.1 pt 5 title 'Filaments'"<<endl;
      // gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_s)<<" w p ps 0.1 pt 9 title 'Sheets'"<<endl;
    }
  if(_MASS_==prop)
    {
      gp_sr<<"set xlabel 'Halo Mass reference [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gp_sr<<"set ylabel 'Halo Mass assigned [Ms/h]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gp_sr<<"plot [1e11:][1e11:]"<<gp_sr.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction'"<<endl;
      //      gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_k)<<" w p ps 0.1 pt 7 title 'Knots'"<<endl;
      //gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_f)<<" w p ps 0.1 pt 5 title 'Filaments'"<<endl;
      // gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_s)<<" w p ps 0.1 pt 9 title 'Sheets'"<<endl;
      //      gp_sr_ref<<"set xlabel 'Halo Mass [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      // gp_sr_ref<<"set ylabel 'Halo Mass [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      //gp_sr_ref<<"plot "<<gp_sr.file1d(xy_sr_ref)<<" w p ps 0.1 pt 7 title 'Reference"<<prop<<"'"<<endl;
      //gp_sr_ref<<"replot "<<gp_sr.file1d(xy_sr_rec)<<" w p ps 0.1 pt 5 title 'Reconstruction "<<prop<<"'"<<endl;
    }
  if(_RS_==prop)
    {
      gp_sr<<"set xlabel 'Rs reference [Kpc/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"'\n";
      gp_sr<<"set ylabel 'Rs assigned [Kpc/h]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"'\n";
      gp_sr<<"plot "<<gp_sr.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction "<<prop<<"'"<<endl;
      gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_k)<<" w p ps 0.1 pt 7 title 'Knots '"<<endl;
      gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_f)<<" w p ps 0.1 pt 5 title 'Filaments '"<<endl;
      gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_s)<<" w p ps 0.1 pt 9 title 'Sheets'"<<endl;
    }
  else if ("MASS_DENS"==prop){
    gp_sr<<"set xlabel 'log(1+{/Symbol d})' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    gp_sr<<"set ylabel 'Vmax [km/s]'  font 'Times-Roman,15' textcolor rgb' "<<FG_COLOR<<"' \n";
    gp_sr<<"plot "<<gp_sr.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction "<<prop<<"'"<<endl;
    gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_k)<<" w p ps 0.1 pt 7 title 'Knots "<<prop<<"'"<<endl;
    gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_f)<<" w p ps 0.1 pt 5 title 'Filaments "<<prop<<"'"<<endl;
    gp_sr<<"replot "<<gp_sr.file1d(xy_pts_ref_s)<<" w p ps 0.1 pt 9 title 'Sheets"<<prop<<"'"<<endl;
  }
  else if (_SPIN_==prop || _SPIN_BULLOCK_==prop)
    gp_sr<<"plot "<<gp_sr.file1d(xy_pts_ref)<<" w p ps 0.2 pt 7 title 'Reconstruction "<<prop<<"', x"<<endl;
   So.DONE();
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_GNUPLOT_
void BiasMT::plot_scaling_relation_assignment(){
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35454;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
    real_prec fraction_tracers_to_file=0.06;
    ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction_tracers_to_file*this->tracer._NOBJS()));
    vector<pair<real_prec, real_prec> > xy_pts_ref;
    ULONG counter=0;
    this->gp_swap<<"set border linecolor '"<<FG_COLOR<<"' \n";
    this->gp_swap<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_swap<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_swap<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    this->gp_swap<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    this->gp_swap<<"set border linewidth 2.0\n";
    this->gp_swap<<"set size square \n";
    do{
      ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].vmax, this->tracer.Halo[i].vmax));
#else
      xy_pts_ref.push_back(std::make_pair(this->tracer_ref.Halo[i].mass, this->tracer.Halo[i].mass));
#endif
      counter++;
  }while(counter<Nobjs_fraction);
  this->gp_swap<<"set log"<<"\n";
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    this->gp_swap<<"set xlabel 'Vmax reference [km/s]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap<<"set ylabel 'Vmax assigned [km/s]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap<<"set grid \n";
  this->gp_swap<<"plot [140:1000][140:1000] x, "<<gp_swap.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction'"<<endl;
#else
  this->gp_swap<<"set xlabel 'Halo Mass reference [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap<<"set ylabel 'Halo Mass assigned [Ms/h]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap<<"set grid \n";
  this->gp_swap<<"plot [1e11:][1e11:] x, "<<gp_swap.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction'"<<endl;
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_GNUPLOT_
void BiasMT::plot_scaling_relation_assignment_bias(){
   const gsl_rng_type * rng_t;
   gsl_rng * gBaseRando;
   gsl_rng_env_setup();
   gsl_rng_default_seed=354;
   rng_t = gsl_rng_mt19937;//_default;
   gBaseRando = gsl_rng_alloc (rng_t);
   real_prec fraction_tracers_to_file=0.06;
   ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction_tracers_to_file*this->tracer._NOBJS()));
   vector<pair<real_prec, real_prec> > xy_pts_ref;
   ULONG counter=0;
   this->gp_swap_b<<"set border linecolor '"<<FG_COLOR<<"' \n";
   this->gp_swap_b<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
   this->gp_swap_b<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
   this->gp_swap_b<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
   this->gp_swap_b<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
   this->gp_swap_b<<"set border linewidth 2.0\n";
   this->gp_swap_b<<"set size square \n";
   do{
     ULONG i= gsl_rng_uniform_int(gBaseRando,this->tracer._NOBJS());
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      xy_pts_ref.push_back(std::make_pair(log10(this->tracer.Halo[i].vmax), this->tracer.Halo[i].bias));
#else
     xy_pts_ref.push_back(std::make_pair(log10(this->tracer.Halo[i].mass), this->tracer.Halo[i].bias));
#endif
     counter++;
  }while(counter<Nobjs_fraction);
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  this->gp_swap_b<<"set xlabel 'log Vmax reference ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap_b<<"set ylabel 'Effective bias '  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap_b<<"set grid \n";
  this->gp_swap_b<<"plot [1.5:3.5][-80:80] "<<this->gp_swap_b.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction'"<<endl;
#else
  this->gp_swap_b<<"set xlabel 'Halo Mass reference [Ms/h]' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap_b<<"set ylabel 'Halo Mass assigned [Ms/h]'  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_swap_b<<"set grid \n";
  this->gp_swap_b<<"plot [11:][-50:50] x, "<<gp_swap.file1d(xy_pts_ref)<<" w p ps 0.1 pt 7 title 'Reconstruction'"<<endl;
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_mean_scaling_relation_assignment_bias(string prop){  //esto hay que pedirselo a Catalog or split loops for ref and mock separately
  So.enter(__PRETTY_FUNCTION__);
  int Nbins=N_BINS_BIAS;
  s_info_in_bins bias_info;
  bias_info.allocate_all(Nbins);
  bias_info.name_info=prop;
  this->tracer.get_mean_bias_relation(bias_info);                // Get bias_prop relation from mock
  s_info_in_bins bias_info_ref;
  bias_info_ref.allocate_all(Nbins);
  bias_info_ref.name_info=prop;
  this->tracer_ref.get_mean_bias_relation(bias_info_ref); // Get bias_prop relation from reference
  string file_o=this->params._Output_directory()+"halo_bias";
  string bfile=file_o+prop+"_with_dm_";
#ifdef  _USE_CWC_
  bfile+="cwc_";
#endif
#ifdef  _USE_MKNOTS_
  bfile+="mk_";
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
  bfile+="indiv_bias_";
#endif
#ifdef _USE_MACH_NUMBER_
  bfile+="mach5_";
#endif
#ifdef  _USE_LOCAL_OVERDENSITY_
  bfile+="delta5_";
#endif
#ifdef  _USE_TIDAL_ANISOTROPY_SEC_PROP_
  bfile+="ta_";
#endif
#ifdef _MULTISCALE_
  bfile+="MS_";
#endif
  bfile+=".txt";
  ofstream sbal;
  So.message_screen("Writing bias in file", bfile);
  sbal.open(bfile.c_str());
  for(int i=0;i<Nbins;++i)
    if(bias_info_ref.vq1[i]>0)
       sbal<<bias_info.vbin[i]<<"\t"<<bias_info_ref.vq1[i]<<"\t"<<bias_info_ref.vq3[i]<<"\t"<<bias_info.vq1[i]<<"\t"<<bias_info.vq3[i]<<endl;
  sbal.close();
#ifdef _USE_GNUPLOT_
  vector<pair<real_prec, real_prec> > xy_pts;
  for(int i=0;i<Nbins;++i)
   xy_pts.push_back(std::make_pair(bias_info.vbin[i], bias_info.vq1[i]));
  vector<pair<real_prec, real_prec> > xy_pts_ref;
  for(int i=0;i<Nbins;++i)
   xy_pts_ref.push_back(std::make_pair(bias_info_ref.vbin[i], bias_info_ref.vq1[i]));
  Gnuplot gplot;
  gplot<<"set border linecolor '"<<FG_COLOR<<"' \n";
  gplot<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
  gplot<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
  this->gp_swap_b<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
  gplot<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
  gplot<<"set border linewidth 2.0\n";
  gplot<<"set size square \n";
  gplot<<"set grid \n";
  gplot<<"set ylabel 'Effective bias '  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
  if(prop==_VMAX_)
   {
      gplot<<"set xlabel 'log Vmax reference ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gplot<<"plot [1.5:3.5][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
    }
   else if (prop==_MASS_)
    {
      gplot<<"set xlabel 'log Mvir' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gplot<<"plot [11:15][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
   }
   else if (prop==_SPIN_BULLOCK_)
    {
      gplot<<"set xlabel 'log Spin' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gplot<<"plot [-4:0.5][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
   }
   else if (prop==_RS_)
    {
      gplot<<"set xlabel 'log R_s' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
      gplot<<"plot [0:3][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
   }
   gplot<<"replot  "<<gplot.file1d(xy_pts_ref)<<" w l lw 2 lt 2  title 'Reference'"<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void BiasMT::get_mean_scaling_relation_assignment_secondary_bias(string primary_prop, string secondary_prop){  //esto hay que pedirselo a Catalog or split loops for ref and mock separately
  So.enter(__PRETTY_FUNCTION__);
  int Nbins=N_BINS_BIAS;
  int Nqrts=4+1; // NUmber of qaurtiles plus the first
// -------------------------------------------------------
  vector<s_info_in_bins>secondary_bias_info_ref(Nqrts);  //container for seconday bias: 2  the number of quartiles to use 
  for(int i=0;i<secondary_bias_info_ref.size();++i)
    secondary_bias_info_ref[i].allocate_all(Nbins);  
  for(int i=0;i<secondary_bias_info_ref.size();++i)
    secondary_bias_info_ref[i].name_info=primary_prop;  
  for(int i=0;i<secondary_bias_info_ref.size();++i)
    secondary_bias_info_ref[i].name_info_sec =secondary_prop;  
  this->tracer_ref.get_mean_secondary_bias_relation(secondary_bias_info_ref); // Get bias_prop relation from mock
// -------------------------------------------------------
  vector<s_info_in_bins>secondary_bias_info(Nqrts);  //container for seconday bias: 2  the number of quartiles to use 
  for(int i=0;i<secondary_bias_info.size();++i)
    secondary_bias_info[i].allocate_all(Nbins);  
  for(int i=0;i<secondary_bias_info.size();++i)
    secondary_bias_info[i].name_info=primary_prop;  
  for(int i=0;i<secondary_bias_info.size();++i)
    secondary_bias_info[i].name_info_sec=secondary_prop;  
  this->tracer.get_mean_secondary_bias_relation(secondary_bias_info); // Get bias_prop relation from reference
// -------------------------------------------------------
  string file_o=this->params._Output_directory()+"secondary_halo_bias";
  string bfile=file_o+primary_prop+secondary_prop+"_with_dm_";
#ifdef  _USE_CWC_
  bfile+="cwc_";
#endif
#ifdef  _USE_MKNOTS_
  bfile+="mk_";
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
  bfile+="indiv_bias_";
#endif
#ifdef _USE_MACH_NUMBER_
  bfile+="mach5_";
#endif
#ifdef  _USE_LOCAL_OVERDENSITY_
  bfile+="delta5_";
#endif
#ifdef  _USE_TIDAL_ANISOTROPY_SEC_PROP_
  bfile+="ta_";
#endif
#ifdef _MULTISCALE_
  bfile+="MS_";
#endif
  bfile+=".txt";
 for(int iq=1;iq<Nqrts;iq++)//loop over quartiles
  {
    string bbfile=bfile+"_quartile"+to_string(iq);
    ofstream sbal;
    So.message_screen("Writing bias in file", bbfile);
    sbal.open(bbfile.c_str());
    for(int i=0;i<Nbins;++i)
      if(secondary_bias_info_ref[iq].vq1[i]>0)
         sbal<<secondary_bias_info[0].vbin[i]<<"\t"<<secondary_bias_info_ref[iq].vq1[i]<<"\t"<<secondary_bias_info_ref[iq].vq3[i]<<"\t"<<secondary_bias_info[iq].vq1[i]<<"\t"<<secondary_bias_info[iq].vq3[i]<<endl;
    sbal.close();
#ifdef _USE_GNUPLOT_
    vector<pair<real_prec, real_prec> > xy_pts_ref;
    vector<pair<real_prec, real_prec> > xy_pts;
    for(int i=0;i<Nbins;++i)
      xy_pts.push_back(std::make_pair(secondary_bias_info[0].vbin[i], secondary_bias_info[iq].vq1[i]));
    for(int i=0;i<Nbins;++i)
      xy_pts_ref.push_back(std::make_pair(secondary_bias_info_ref[0].vbin[i], secondary_bias_info_ref[iq].vq1[i]));
    Gnuplot gplot;
    gplot<<"set border linecolor '"<<FG_COLOR<<"' \n";
    gplot<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
    gplot<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
    this->gp_swap_b<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
    gplot<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
    gplot<<"set border linewidth 2.0\n";
    gplot<<"set size square \n";
    gplot<<"set grid \n";
    gplot<<"set ylabel 'Effective bias '  font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
    if(primary_prop==_VMAX_)
    {
        gplot<<"set xlabel 'log Vmax reference ' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
        gplot<<"plot [1.5:3.5][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
      }
    else if (primary_prop==_MASS_)
      {
        gplot<<"set xlabel 'log Mvir' font 'Times-Roman,15' textcolor rgb '"<<FG_COLOR<<"' \n";
        gplot<<"plot [11:15][-2:25] "<<gplot.file1d(xy_pts)<<" w p ps 2 pt 7 title 'Reconstruction'"<<endl;
    }
    gplot<<"replot  "<<gplot.file1d(xy_pts_ref)<<" w l lw 2 lt 2  title 'Reference, quartile "<<iq<<"'"<<endl;
#endif
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

