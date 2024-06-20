///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @file PowerSpectrumF.cpp
 * @brief This file contains a public PowerSpectrum-class member function
 * @details Generates thee estimates of power spectrum
 * @author Andres Balaguera Antolinez
 * @date 2008-2023
 */
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../headers/PowerSpectrumF.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::set_output_filenames ()
{
#ifdef _USE_REDSHIFT_BINS_
  string aux_name=this->params._Name_survey()+"_zmin_"+to_string(this->params._redshift_min_sample()) +"_zmax_"+to_string(this->params._redshift_max_sample());
  this->params.set_Name_survey(aux_name);
#endif
#ifdef _REDSHIFT_SPACE_
  this->file_MCF = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_mark"+this->params._mark()+"_Real"+to_string(this->params._realization())+"_RSS.txt";
  this->file_power  = this->params._Output_directory()+this->params._statistics()+"_"+params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+"_RSS.txt";
  this->file_power_marked  = this->params._Output_directory()+this->params._statistics()+"_marked_"+params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+"_RSS.txt";
  this->file_power_log      = this->params._Output_directory()+this->params._statistics()+"_"+params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power_log()+"_RSS.log";
#else
  this->file_power  = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+this->params._file_power()+".txt";
  this->file_power_marked  = this->params._Output_directory()+this->params._statistics()+"_marked_"+params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+".txt";
  this->file_power_log = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power_log()+".log";
  this->file_MCF = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_mark"+this->params._mark()+"_Real"+to_string(this->params._realization())+".txt";
#endif
  this->file_power  = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+this->params._file_power()+".txt";
  this->file_power_real_space  = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+".txt";
  this->file_power_redshift_space  = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+"_RSS.txt";
  this->file_power_marked_real_space  = this->params._Output_directory()+this->params._statistics()+"_marked_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+".txt";
  this->file_power_marked_redshift_space  = this->params._Output_directory()+this->params._statistics()+"_marked_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+"_RSS.txt";
  this->file_power_cross  = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power()+".txt_CROSS";
  this->file_dndz           = this->params._Output_directory()+"dndz_"+params._Name_survey()+this->params._file_dndz()+".txt";
  this->file_power2d        = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power2d()+".txt";
  this->file_power2d_mk     = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_power2d_mk()+".txt";
  this->file_window         = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_window()+".txt";
  this->file_bispectrum     = this->params._Output_directory()+this->params._statistics()+"_"+this->params._Name_survey()+"_Nft"+to_string(this->params._Nft())+"_"+this->params._mass_assignment_scheme()+"_"+params._file_bispectrum()+".txt";
  //file_power_fb       = file_power_cl; // ???
  this->file_data = this->params._file_catalogue();
  this->file_random = this->params._file_random();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_fftw_parameters()
{
  ofstream out;
  out.open(this->file_power_log.c_str() , ios::out);
  out.precision(12);
  out.setf(ios::showpoint);
  out.width(12);
  time_t rawtime;
  time ( &rawtime );
  out<<"******************************************************************************"<<endl;
  out<<"******************************************************************************"<<endl;
  out<<endl;
  out<<BLUE<<"Log file generated on"<<RESET<<endl;
  out<<ctime(&rawtime)<<endl;
  out<<"******************************************************************************"<<endl;
  out<<BLUE<<"INPUT FILES"<<RESET<<endl;
  out<<"Input parameter file        "<<this->params._par_file()<<endl;
  out<<"Name Survey                 "<<this->params._Name_survey()<<endl;
  out<<"Tracer Catalog              "<<this->params._file_catalogue()<<endl;
  out<<"Number of columns           "<<this->tracer_cat._NCOLS()<<endl;
  if(this->params._use_random_catalog())
  {
      out<<"Random Catalog              "<<this->params._file_random()<<endl;
      out<<"Number of columns           "<<this->random_cat._NCOLS()<<endl;
      out<<"Smoothed dNdz               "<<this->params._file_nbar()<<endl;
  }
  else{
      out<<"No random catalog requested"<<endl;
  }
  out<<"******************************************************************************"<<endl;
  out<<BLUE<<"SELECTED OPTIONS"<<RESET<<endl;
  out<<"Statistics                  "<<this->params._statistics()<<endl;
  out<<"Mass Assignment Scheme      "<<this->params._mass_assignment_scheme()<<endl;
  string enabbled = true == this->params._MAS_correction() ? "true" : "false";
  out<<"MAS correction              "<<enabbled<<endl;
  out<<"Lbox                        "<<this->params._Lbox()<<" Mpc/h"<<endl;
  out<<"Nft                         "<<this->params._Nft()<<endl;
  enabbled = true == this->params._FKP_weight() ? "true" : "false";
  out<<"FKP weight                  "<<enabbled<<endl;
  out<<"P_0                         "<<this->params._Pest()<<" (Mpc/h)^3"<<endl;
  enabbled = true == this->params._SN_correction()  ? "true" : "false";
  out<<"Shot-noise correction       "<<enabbled<<endl;
  out<<"******************************************************************************"<<endl;
  out<<BLUE<<"DERIVED QUANTITIES"<<RESET<<endl;
  out<<"Number of objects           "<<this->fftw_functions._n_gal()<<endl;
  out<<"Weighted number of objects  "<<this->fftw_functions._w_g()<<endl;
  if(this->params._use_random_catalog())
      out<<"Number of random objects    "<<this->fftw_functions._n_ran()<<endl;
  if(this->params._use_random_catalog())
      out<<"Weighted number of randoms  "<<this->fftw_functions._w_r()<<endl;
  out<<"alpha                       "<<this->fftw_functions._alpha()<<endl;
  out<<"Shot_noise                  "<<this->fftw_functions._shot_noise()<<" (Mpc/h)^3"<<endl;
  if(this->params._use_random_catalog())
      out<<"Shot_noise (window)         "<<this->fftw_functions._shot_noise_window()<<endl;
  out<<"Normalization               "<<this->fftw_functions._normal_power()<<endl;
  if(this->params._use_random_catalog())
    out<<"Mean nbar(from weights)     "<<(this->fftw_functions._n_ran()-this->fftw_functions._w_r())/(this->params._Pest()*this->fftw_functions._w_r())<<" (Mpc h^-1)^(-3)"<<endl;
  out<<"Sum nw² galaxies            "<<this->fftw_functions._s_g()<<endl;
  if(this->params._use_random_catalog())
      out<<"Sum nw² randoms             "<<this->fftw_functions._s_r()<<endl;
  out<<"******************************************************************************"<<endl;
  out<<BLUE<<"FOURIER INFORMATION"<<RESET<<endl;
  out<<"Fundamental mode            "<<this->params._d_deltak_0()<<" h/Mpc"<<endl;
  out<<"Bin size for P(k)           "<<this->params._d_DeltaK_data()<<" h/Mpc"<<endl;
  if(this->params._use_random_catalog())
      out<<"Bin size for W(k)           "<<this->params._d_DeltaK_window()<<" h/Mpc"<<endl;
  out<<"Nyquist Frequency           "<<0.5*this->params._Nft()*(this->params._d_deltak_0())<<" h/Mpc"<<endl;
  out<<"******************************************************************************"<<endl;
  if(this->params._sys_of_coord_g()>0)
  {
    out<<BLUE<<"COSMOLOGICAL PARAMETERS"<<RESET<<endl;
    out<<"Omega Matter                "<<this->params.s_cosmo_pars.Om_matter<<endl;
    out<<"Omega Vac                   "<<this->params.s_cosmo_pars.Om_vac<<endl;
    out<<"Omega baryons               "<<this->params.s_cosmo_pars.Om_baryons<<endl;
    out<<"Omega curvature             "<<this->params.s_cosmo_pars.Om_k<<endl;
    out<<"wde_eof                       "<<this->params.s_cosmo_pars.wde_eos<<endl;
    out<<"Hubble parameter            "<<this->params.s_cosmo_pars.Hubble<<endl;
    out<<"Sigma_8                     "<<this->params.s_cosmo_pars.sigma8<<endl;
    out<<"Minimim redshift            "<<this->params._redshift_min_sample()<<endl;
    out<<"Maximum redshift            "<<this->params._redshift_max_sample()<<endl;
  }
  else
    out<<"No cosmologial parameters have been used. "<<endl;
  out.close();
  out<<"******************************************************************************"<<endl;
  out<<"******************************************************************************"<<endl;
  So.message_screen("Log-file written in ",this->file_power_log);
}
////////////////////////////////////////////////////////////////////////////
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void PowerSpectrumF::add_catalogues(real_prec mcut)
{
#elif defined (_USE_MASS_BINS_PK_)
  void PowerSpectrumF::add_catalogues(real_prec m_min, real_prec m_max)
  {
#endif
  this->So.enter(__PRETTY_FUNCTION__);
  this->tracer_cat.set_params(this->params);
  this->tracer_cat.set_type_of_object(this->params._type_of_object());
  this->tracer_cat.read_catalog(this->params._Input_dir_cat()+ this->file_data, mcut);
  this->N_galaxy=this->tracer_cat._NOBJS();
  if(true==this->params._use_random_catalog())
    {
        this->random_cat.set_params(this->params);// this must come first than set_type_of_object
        this->random_cat.set_type_of_object("RANDOM");
        this->random_cat.read_catalog(this->params._Input_dir_cat()+this->file_random, mcut);
        this->N_random=this->random_cat._NOBJS();
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::add_catalogues(int b1, int b2, int b3) // three bins
{
   this->So.enter(__PRETTY_FUNCTION__);
   this->tracer_cat.set_params(this->params);
   this->tracer_cat.set_type_of_object(this->params._type_of_object());
   this->tracer_cat.read_catalog(this->params._Input_dir_cat()+ this->file_data, b1, b2, b3);
   this->N_galaxy=this->tracer_cat._NOBJS();
   if(true==this->params._use_random_catalog())
     {
       this->random_cat.set_params(this->params);// this must come first than set_type_of_object
       this->random_cat.set_type_of_object("RANDOM");
       this->random_cat.read_catalog(this->params._Input_dir_cat()+this->file_random, b1, b2, b3);
       this->N_random=this->random_cat._NOBJS();
   }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_SEVERAL_RANDOM_FILES_
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void PowerSpectrumF::add_catalogues_ran(real_prec mcut, int ir)
{
#elif defined (_USE_MASS_BINS_PK_)
  void PowerSpectrumF::add_catalogues_ran(real_prec m_min, real_prec m_max,int it)
  {
#endif
    this->So.enter(__PRETTY_FUNCTION__);
    this->random_cat.set_params(this->params);// this must come first than set_type_of_object
    this->random_cat.set_type_of_object("RANDOM");
#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->params._RANDOMfile(ir) , mcut);
#elif defined (_USE_MASS_BINS_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->params._RANDOMfile(ir),m_min,m_max);
#endif
    this->N_random=this->random_cat._NOBJS();
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::add_catalogues_ran(int ir, int b1, int b2, int b3)
{
   this->So.enter(__PRETTY_FUNCTION__);
    this->random_cat.set_params(this->params);// this must come first than set_type_of_object
    this->random_cat.set_type_of_object("RANDOM");
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->params._RANDOMfile(ir) , b1, b2, b3);
    this->N_random=this->random_cat._NOBJS();
  }
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::add_catalogues_ran(int b1, int b2, int b3)
{
    this->So.enter(__PRETTY_FUNCTION__);
    this->random_cat.set_params(this->params);// this must come first than set_type_of_object
    this->random_cat.set_type_of_object("RANDOM");
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->params._file_random(), b1, b2, b3);
    this->N_random=this->random_cat._NOBJS();
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void PowerSpectrumF::add_catalogues_ran(real_prec mcut)
{
#elif defined (_USE_MASS_BINS_PK_)
  void PowerSpectrumF::add_catalogues_ran(real_prec m_min, real_prec m_max)
  {
#endif
    this->So.enter(__PRETTY_FUNCTION__);
    this->random_cat.set_params(this->params);// this must come first than set_type_of_object
    this->random_cat.set_type_of_object("RANDOM");
#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->file_random, mcut);
#elif defined (_USE_MASS_BINS_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+file_random,m_min,m_max);
#endif
    this->N_random=this->random_cat._NOBJS();
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void PowerSpectrumF::add_catalogues_gal(real_prec mcut)
{
#elif defined (_USE_MASS_BINS_PK_)
  void PowerSpectrumF::add_catalogues_gal(real_prec m_min, real_prec m_max)
  {
#endif
    this->So.enter(__PRETTY_FUNCTION__);
    this->tracer_cat.set_params(this->params);
    this->tracer_cat.set_type_of_object(this->params._type_of_object());

#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    this->tracer_cat.read_catalog(this->params._Input_dir_cat()+ this->file_data, mcut);
#elif defined (_USE_MASS_BINS_PK_)
    this->tracer_cat.read_catalog(this->params._Input_dir_cat()+file_data,m_min,m_max);
#endif
    this->N_galaxy=this->tracer_cat._NOBJS();
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::add_catalogues_gal(int i1, int i2, int i3)
{
    this->So.enter(__PRETTY_FUNCTION__);
    this->tracer_cat.set_params(this->params);
    this->tracer_cat.set_type_of_object(this->params._type_of_object());
    this->tracer_cat.read_catalog(this->params._Input_dir_cat()+ this->file_data, i1, i2, i3);
    this->N_galaxy=this->tracer_cat._NOBJS();
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::add_random_catalogue()
  {
    this->So.enter(__PRETTY_FUNCTION__);
    this->random_cat.set_params(this->params);// this must come first than set_type_of_object
    this->random_cat.set_type_of_object("RANDOM");
#if defined(_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+this->file_random, 0.);
#elif defined (_USE_MASS_BINS_PK_)
    this->random_cat.read_catalog(this->params._Input_dir_cat()+file_random,0,LARGE_NUMBER);
#endif
    this->N_random=this->random_cat._NOBJS();
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_power_spectrum()
  {
    this->So.enter(__PRETTY_FUNCTION__);
    if(this->params._statistics()=="Pk_fkp")
      {
#ifdef _WRITE_MULTIPOLES_
        if(true==this->params._FKP_error_bars())
            File.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
        else
            File.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->modes_g);
#else
    File.write_to_file(this->file_power,this->kvector_data,this->pk0,this->modes_g);
#endif
    // Write P(kperp, kpar) to file:
#ifdef _WRITE_2DPOWER_
    File.write_to_file(this->file_power2d,this->kvector_data2d,this->kvector_data2d,this->pkk);
    //Write P(k, mu) to file:
    File.write_to_file(file_power2d_mk,this->muvector,this->kvector_data2d,this->pmk);
#endif
    // Write W(k):
    if(true==this->params._use_random_catalog())
      File.write_to_file(this->file_window,this->kvector_window,this->pk_w);
      }
    else if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_yb"  || this->params._statistics()=="Pk_ybc" || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_y_ds")
      {
        if(true==this->params._FKP_error_bars())
            File.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
        else
            File.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->modes_g);
    }
    else if(this->params._statistics()=="Bk_fkp")
      File.write_to_file(file_bispectrum,kvector_data_b,bispectrum,sn_bispectrum, modes_tri);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_power_spectrum(bool write_sigma) ///PLeas unify the name of this function with write_power and modes
  {
#ifdef _VERBOSE_POWER_
    this->So.enter(__PRETTY_FUNCTION__);
#endif
    if(this->params._statistics()=="Pk_fkp")
      {
#ifdef _WRITE_MULTIPOLES_
        if(true==write_sigma)
      File.write_to_file(this->file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
        else
      File.write_to_file(this->file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->modes_g);
#else
        File.write_to_file2(this->file_power,this->kvector_data,this->pk0,this->modes_g,write_sigma);
#endif
    // Write P(kperp, kpar) to file:
#ifdef _WRITE_2DPOWER_
    File.write_to_file(this->file_power2d,this->kvector_data2d,this->kvector_data2d,this->pkk);
    //Write P(k, mu) to file:
    File.write_to_file(file_power2d_mk,this->muvector,this->kvector_data2d,this->pmk);
#endif
    // Write W(k):
    if(true==this->params._use_random_catalog())
      File.write_to_file(this->file_window,this->kvector_window,this->pk_w);
      }
    else if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_yb"  || this->params._statistics()=="Pk_ybc" || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_y_ds")
      {
    File.write_to_file(file_power,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
      }
    else if(this->params._statistics()=="Bk_fkp")
      File.write_to_file(file_bispectrum,kvector_data_b,bispectrum,sn_bispectrum, modes_tri);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_power_and_modes()
  {
#ifdef _VERBOSE_POWER_
    this->So.enter(__PRETTY_FUNCTION__);
#endif
    File.write_to_file2(this->file_power,this->kvector_data,this->pk0,this->modes_g,true);
    if(this->params._use_random_catalog())
      File.write_to_file2(this->file_window,this->kvector_data,this->pk_w,this->modes_g,true);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_power_and_modes(string file)
  {
    this->So.enter(__PRETTY_FUNCTION__);
    File.write_to_file2(file,this->kvector_data,this->pk0,this->modes_g, true);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::write_power_spectrum_grid(string output_file)
  {
#ifdef _VERBOSE_POWER_
    this->So.enter(__PRETTY_FUNCTION__);
    So.message_screen("Writing outputs");
#endif
    File.write_to_file(output_file,this->kvector_data,this->pk0,this->pk2,this->pk4,this->sigma_fkp,this->modes_g);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_marked_correlation_function()
  {
    this->So.enter(__PRETTY_FUNCTION__);
    this->params.set_statistics("MCF");
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    this->add_catalogues(0);
#else
    this->add_catalogues(0,1e20);
#endif
    vector<ULONG>count(this->params._Nbins_cf(),0);
    vector<real_prec>mcount(this->params._Nbins_cf(),0);
    vector<real_prec>vcount(this->params._Nbins_cf(),0);
    vector<real_prec>mvcount(this->params._Nbins_cf(),0);
    vector<real_prec> rbin(this->params._Nbins_cf(),0);
    real_prec Deltar=this->params._rbin_type() == "linear"? (this->params._rmax_cf()-this->params._rmin_cf())/(static_cast<real_prec>(this->params._Nbins_cf())) :  (log10(this->params._rmax_cf()/this->params._rmin_cf()))/static_cast<real_prec>(this->params._Nbins_cf());
    int NTHREADS = 1;
#ifdef _USE_OMP_
    NTHREADS=_NTHREADS_;
#endif
    So.message_screen("Measuring Marked Correlation function using",NTHREADS," threads");
    if(this->params._rbin_type() =="linear")
      for (int i=0;i<rbin.size();++i)
        rbin[i]=this->params._rmin_cf()+(i+0.5)*Deltar;
    else
      for (int i=0;i<rbin.size();++i)
    rbin[i]=pow(10, log10(this->params._rmin_cf())+(i+0.5)*Deltar);
    real_prec mean_mass=0;
    if(this->params._i_mass_g()>0){
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_mass)
#endif
      for (ULONG i=0;i<this->N_galaxy; ++i)
    mean_mass+=this->tracer_cat.Halo[i].mass;
      mean_mass=(static_cast<double>(mean_mass))/(static_cast<double>(this->N_galaxy));
#ifdef _FULL_VERBOSE_
      So.message_screen("Mean mass = ", mean_mass);
#endif
    }
    real_prec mean_vmax=0;
    if(this->params._i_vmax_g()>0){
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_vmax)
#endif
      for (ULONG i=0;i<this->N_galaxy; ++i)
    mean_vmax+=this->tracer_cat.Halo[i].vmax;
      mean_vmax=(static_cast<double>(mean_vmax))/(static_cast<double>(this->N_galaxy));
#ifdef _FULL_VERBOSE_
      So.message_screen("Mean vmax = ", mean_vmax);
#endif
    }
#ifdef _REDSHIFT_SPACE_
    real_prec rsd_x=1.0;
    real_prec rsd_y=1.0;
    real_prec rsd_z=1.0;
    if(false==this->params._redshift_space_coords_g())
      switch(LOS)
    {
    case(I_LOSX):rsd_x*=1.0;rsd_y*=0.0;rsd_z*=0.0; break;
    case(I_LOSY):rsd_x*=0.0;rsd_y*=1.0;rsd_z*=0.0; break;
    case(I_LOSZ):rsd_x*=0.0;rsd_y*=0.0;rsd_z*=1.0; break;
    }
    else // of coordinates have already the rsd included (as in a real catalog) and we still use _REDSHIFT_SPACE_ we have to set all rsd to 0
      {
    rsd_x=0.0;
    rsd_y=0.0;
    rsd_z=0.0;
      }
    // Conversion from km/s to Mpc/h for RSD
    real_prec conversion_factor=1.;
    if(this->params._redshift_space_coords_g() == false)
      {
        if("kmps"==this->params._vel_units_g())
          conversion_factor=(1.+this->params._redshift())/(this->cosmology.Hubble_function(this->params._redshift()));
        else if("alpt"==this->params._vel_units_g())
          conversion_factor= cgs_Mpc/(this->cosmology.Hubble_function(this->params._redshift()));
        else if("Mpcph"==this->params._vel_units_g())
          conversion_factor=1;
    #ifdef _FULL_VERBOSE_
        So.message_screen("Current redshift =",this->params._redshift());
        So.message_screen("Hubble function at current redshift =",this->cosmology.Hubble_function(this->params._redshift()));
        So.message_screen("Conversion factor =",conversion_factor);
    #endif
        So.message_screen("Hubble function at current redshift =",this->cosmology.Hubble_function(this->params._redshift()));
          }
#endif
    real_prec lrmin= this->params._rbin_type() == "linear" ? this->params._rmin_cf():  log10(this->params._rmin_cf());
    vector<ULONG> count_priv(this->params._Nbins_cf()*NTHREADS,0);
    vector<real_prec> mcount_priv(this->params._Nbins_cf()*NTHREADS,0);
    vector<real_prec> vcount_priv(this->params._Nbins_cf()*NTHREADS,0);
    vector<real_prec> mvcount_priv(this->params._Nbins_cf()*NTHREADS,0);
#pragma omp parallel num_threads(NTHREADS)
    {
      const int ithread = omp_get_thread_num();
      for(ULONG i=0;i<this->params._Nbins_cf();++i)count_priv[this->params._Nbins_cf()*ithread+i]=0;
      for(ULONG i=0;i<this->params._Nbins_cf();++i)count_priv[this->params._Nbins_cf()*ithread+i]=0;
#pragma omp for
      for (ULONG i=0;i<this->N_galaxy;++i)
    {
#ifdef _REDSHIFT_SPACE_
      real_prec vx=this->tracer_cat.Halo[i].vel1;
      real_prec vy=this->tracer_cat.Halo[i].vel2;
      real_prec vz=this->tracer_cat.Halo[i].vel3;
      real_prec x=this->tracer_cat.Halo[i].coord1+rsd_x*conversion_factor*vx;
      real_prec y=this->tracer_cat.Halo[i].coord2+rsd_y*conversion_factor*vy;
      real_prec z=this->tracer_cat.Halo[i].coord3+rsd_z*conversion_factor*vz;
#else
      real_prec x=this->tracer_cat.Halo[i].coord1;
      real_prec y=this->tracer_cat.Halo[i].coord2;
      real_prec z=this->tracer_cat.Halo[i].coord3;
#endif
      real_prec mass=this->tracer_cat.Halo[i].mass;
      real_prec sigma=this->tracer_cat.Halo[i].vmax;
      for (ULONG j=i+1;j<this->N_galaxy;++j)
        {
#ifdef _REDSHIFT_SPACE_
          real_prec vxp=this->tracer_cat.Halo[j].vel1;
          real_prec vyp=this->tracer_cat.Halo[j].vel2;
          real_prec vzp=this->tracer_cat.Halo[j].vel3;
          real_prec dx=x-(this->tracer_cat.Halo[j].coord1+rsd_x*conversion_factor*vxp);
          real_prec dy=y-(this->tracer_cat.Halo[j].coord2+rsd_y*conversion_factor*vyp);
          real_prec dz=z-(this->tracer_cat.Halo[j].coord3+rsd_z*conversion_factor*vzp);
#else
          real_prec dx=x-this->tracer_cat.Halo[j].coord1;
          real_prec dy=y-this->tracer_cat.Halo[j].coord2;
          real_prec dz=z-this->tracer_cat.Halo[j].coord3;
#endif
          real_prec pmass=this->tracer_cat.Halo[j].mass;
          real_prec psigma=this->tracer_cat.Halo[j].vmax;
          real_prec idr=dx*dx+dy*dy+dz*dz;
          real_prec rd=  this->params._rbin_type() == "linear" ? sqrt(idr): 0.5*log10(idr);
          ULONG ind = get_bin(rd,lrmin,this->params._Nbins_cf(),Deltar,false);
          if(ind < this->params._Nbins_cf())
        {
          count_priv[ind+this->params._Nbins_cf()*ithread]++;
          mcount_priv[ind+this->params._Nbins_cf()*ithread]+=0.5*(mass+pmass);
          vcount_priv[ind+this->params._Nbins_cf()*ithread]+=0.5*(sigma+psigma);
          mvcount_priv[ind+this->params._Nbins_cf()*ithread]+=0.25*((mass+pmass)*(psigma+sigma));
        }
    }
}
#pragma omp for
   for(ULONG i=0;i<this->params._Nbins_cf();++i)
    for(ULONG t=0;t<NTHREADS;++t){
      mcount[i]+= mcount_priv[i+this->params._Nbins_cf()*t];
      vcount[i]+= vcount_priv[i+this->params._Nbins_cf()*t];
      mvcount[i]+= mvcount_priv[i+this->params._Nbins_cf()*t];
    }
#pragma omp for
      for(ULONG i=0;i<this->params._Nbins_cf() ;++i)
    for(ULONG t=0;t<NTHREADS;++t)
      count[i]+= count_priv[i+this->params._Nbins_cf()*t];
    }
    this->tracer_cat.Halo.clear();
    this->tracer_cat.Halo.shrink_to_fit();
    So.DONE();
    cout<<CYAN<<" Done"<<RESET<<endl;
    for(int i=0;i<mcount.size();++i)
      mcount[i]/=(pow(mean_mass,1)*static_cast<real_prec>(count[i]));
    for(int i=0;i<mcount.size();++i)
      vcount[i]/=(pow(mean_vmax,1)*static_cast<real_prec>(count[i]));
    for(int i=0;i<mcount.size();++i)
      mvcount[i]/=(mean_vmax*mean_mass*static_cast<real_prec>(count[i]));
    ofstream mcor; mcor.open(this->file_MCF.c_str());
    mcor.precision(8);
    for(int i=0;i<mcount.size();++i)mcor<<rbin[i]<<"\t"<<mcount[i]<<"\t"<<vcount[i]<<"\t"<<mvcount[i]<<"\t"<<count[i]<<endl;
    mcor.close();
    cout<<CYAN<<"Output in file "<<GREEN<<this->file_MCF.c_str()<<RESET<<endl;
    cout<<CYAN<<"Done"<<endl;
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_cross_power_spectrum_grid(bool dm, vector<real_prec>&X,vector<real_prec>&Y, bool  get_cross_coeff)
  {
    this->So.enter(__PRETTY_FUNCTION__);
    real_prec ngal_new=0;
    this->fftw_functions.data_g.clear();
    this->fftw_functions.data_g.shrink_to_fit();
    this->fftw_functions.data_g.resize(this->params._NGRID(),0);
    if (this->params._input_type_two()=="delta_grid")
      {
#pragma omp parallel for
        for(ULONG i=0;i<this->params._NGRID();++i)
          this->fftw_functions.set_data_g(i,static_cast<real_prec>(Y[i]));
        ngal_new=this->params._ngal_delta();
      }
    else if (this->params._input_type_two()=="density_grid")
      {
#pragma omp parallel for reduction(+:ngal_new)
    for(ULONG i=0;i<this->params._NGRID();++i)
      ngal_new+=static_cast<real_prec>(Y[i]);
    real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
    So.message_screen("Y-type density_grid with mean",nmean);
#pragma omp parallel for
        for(ULONG i=0;i<this->params._NGRID();++i)
          this->fftw_functions.set_data_g(i,(static_cast<real_prec>(Y[i])/static_cast<real_prec>(nmean))-1.);
    this->params.set_ngal_delta(ngal_new);
      }
    this->fftw_functions.set_shot_noise(0);
    real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
    this->fftw_functions.set_normal_power(pow(factor,-2));
    this->fftw_functions.data_gp.clear();
    this->fftw_functions.data_gp.shrink_to_fit();
    this->fftw_functions.data_gp.resize(this->params._NGRID(),0);
    ngal_new=0;
    if (this->params._input_type()=="delta_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
          this->fftw_functions.set_data_gp(i,static_cast<real_prec>(X[i]));
    So.message_screen("X-type delta_grid");
      }
    else if (this->params._input_type()=="density_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
      ngal_new+=static_cast<real_prec>(X[i]);
    real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
    So.message_screen("X-type density_grid with mean", nmean);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
          this->fftw_functions.set_data_gp(i,(static_cast<real_prec>(X[i])/static_cast<real_prec>(nmean))-1.);
    }
    this->fftw_functions.set_shot_noise2(0);
    factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
    this->fftw_functions.set_normal_power_two(pow(factor,-2));
    this->fftw_functions.resize_fftw_vectors();
    this->kvector_data.clear();
    this->kvector_data.shrink_to_fit();
    for(ULONG i=0;i<this->params._d_Nnp_data();i++)
      this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
    this->pk0.clear();
    this->pk0.shrink_to_fit();
    this->pk0.resize(this->params._d_Nnp_data(),0);
    this->modes_g.clear();
    this->modes_g.resize(this->params._d_Nnp_data(),0);
    //#if !defined (_USE_BIAS_OBJECT_TO_OBJECT_) || !defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
#if !defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
    fftw_functions.cross_power_spectrum_fkp(this->pk0,this->modes_g,  get_cross_coeff);
#else
    vector<real_prec>corr(this->params._NGRID(),0);
    fftw_functions.cross_power_spectrum_fkp(this->pk0,this->modes_g,corr);
#endif
    this->write_power_and_modes(this->file_power_cross);
    fftw_functions.data_g.clear();
    fftw_functions.data_g.shrink_to_fit();
    fftw_functions.data_gp.clear();
    fftw_functions.data_gp.shrink_to_fit();
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_cross_power_spectrum_grid(bool dm,string file_X, string file_Y, bool get_cross_coeff)
  {
#ifdef _VERBOSE_POWER_
    this->So.enter(__PRETTY_FUNCTION__);
#endif
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
    FileOutput File;
    real_prec ngal_new=0;
    fftw_functions.data_g.clear();
    fftw_functions.data_g.shrink_to_fit();
    fftw_functions.data_g.resize(this->params._NGRID(),0);
    vector<real_prec>Y(this->params._NGRID(),0);
    File.read_array_t<PrecType_Y>(file_Y, Y);
    if (this->params._input_type_two()=="delta_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<this->params._NGRID();++i)
          fftw_functions.data_g[i]=static_cast<real_prec>(Y[i]);
        ngal_new=this->params._ngal_delta();
        So.message_screen("Y-type delta_grid");
      }
    else if (this->params._input_type_two()=="density_grid")
      {
#pragma omp parallel for reduction(+:ngal_new)
    for(ULONG i=0;i<this->params._NGRID();++i)
          ngal_new+=static_cast<real_prec>(Y[i]);
    real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
        So.message_screen("Mean = ", nmean);
        So.message_screen("Y-type density_grid");
#pragma omp parallel for
    for(ULONG i=0;i<this->params._NGRID();++i)
          fftw_functions.data_g[i]=(static_cast<real_prec>(Y[i])/static_cast<real_prec>(nmean))-1.;
    this->params.set_ngal_delta(ngal_new);
      }
    if(true==this->params._SN_correction())
      fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
    real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
    fftw_functions.set_normal_power(pow(factor,-2));
    Y.clear(); Y.shrink_to_fit();
    fftw_functions.data_gp.clear();
    fftw_functions.data_gp.shrink_to_fit();
    fftw_functions.data_gp.resize(this->params._NGRID(),0);
    ngal_new=0;
    vector<float>X(this->params._NGRID(),0);
    File.read_array_t<float>(file_X, X);
    if (this->params._input_type()=="delta_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
          fftw_functions.data_gp[i]=static_cast<real_prec>(X[i]);
    ngal_new=this->params._ngal_delta();
    So.message_screen("X-type delta_grid");
      }
    else if (this->params._input_type()=="density_grid")
      {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
        ngal_new+=static_cast<real_prec>(X[i]);
    real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
    So.message_screen("Mean = ", nmean);
    So.message_screen("X-type = density_grid");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->params._NGRID();++i)
          fftw_functions.data_gp[i]=(static_cast<real_prec>(X[i])/static_cast<real_prec>(nmean))-1.;
   }
    fftw_functions.shot_noise2=0;
    factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
    fftw_functions.set_normal_power_two(pow(factor,-2));
    X.clear(); X.shrink_to_fit();
    if(false==dm && true==this->params._SN_correction())
      fftw_functions.shot_noise2=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
    fftw_functions.write_fftw_parameters();
    fftw_functions.resize_fftw_vectors();
    this->kvector_data.clear();
    this->kvector_data.shrink_to_fit();
    for(int i=0;i<this->params._d_Nnp_data();i++)
      this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
    this->pk0.clear();
    this->pk0.shrink_to_fit();
    this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
    this->modes_g.clear();
    this->modes_g.resize(this->params._d_Nnp_data(),0); //Monopole
#if !defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
    fftw_functions.cross_power_spectrum_fkp(this->pk0,this->modes_g, get_cross_coeff);
#else
    vector<real_prec>corr(this->params._NGRID(),0);
    fftw_functions.cross_power_spectrum_fkp(this->pk0,this->modes_g,corr);
#endif
    this->write_power_and_modes();
    fftw_functions.data_g.clear();
    fftw_functions.data_g.shrink_to_fit();
    fftw_functions.data_g.resize(this->params._NGRID(),0);
    fftw_functions.data_gp.clear();
    fftw_functions.data_gp.shrink_to_fit();
    fftw_functions.data_gp.resize(this->params._NGRID(),0);
  }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_power_spectrum(bool verbose, bool mcut){
  So.enter(__PRETTY_FUNCTION__);
  this->So.welcome_message();
#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
    // Output file, keep roor in case intervals are used
  string file_pow=this->file_power;
    // Initizalize number of intervals
  int N_intervals=1;
#ifdef  _JPAS_BINS_POWER_
#ifdef _USE_MASS_AS_OBSERVABLE_POWER_ //deprec
  N_intervals=this->params._NMASSbins_power(); //deprec
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_ //deprec
  N_intervals=this->params._NVMAXbins_power(); //deprec
#endif
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_ //deprec
  for(int i=0;i<N_intervals;++i)cout<<"Bin "<<i<<"  :"<<this->params._VMAXbins_min(i)<<"  "<<this->params._VMAXbins_max(i)<<endl;
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_  //deprec
  for(int i=0;i<N_intervals;++i)cout<<"Bin "<<i<<"  :"<<pow(10,this->params._MASSbins_min(i))<<"  "<<pow(10,this->params._MASSbins_max(i))<<endl;
#endif
  if(I_EQZ==this->params._sys_of_coord_g())
    So.write_cosmo_parameters((void *)&this->params.s_cosmo_pars);
// ***********************************************************
// Define some structures here
  s_data_structure s_data_struct_r; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
  s_data_structure_direct_sum s_data_struct_r_ds; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
// ***********************************************************
  this->File.input_type=this->params._input_type();
#if defined (_USE_MASS_CUTS_PK_) ||  defined (_USE_ALL_PK_)
  if(this->params._input_type()!="catalog")
    So.message_warning("Warning: expecting catalog to perform mass cuts when density field is provided");
  real_prec mm=0;
  if("catalog"==this->params._input_type())
#ifdef _USE_ALL_PK_
    if(this->params._i_mass_g()>0)
      mm=pow(10,this->params._LOGMASSmin());
    else
      mm=0;
#else
    mm=log10(this->params._MASScuts(i));
#endif
  if(this->params._input_type()=="catalog")// This is a long if
    {
           // Tabulate r-z relation if radial coordinate is redshift
      vector<gsl_real>vzz;
      vector<gsl_real>vrc;
      if(I_EQR==this->params._sys_of_coord_g() || I_EQZ==this->params._sys_of_coord_g())
        {
          So.message_screen("Computing comoving distances in the range",this->params._redshift_min_sample(),"<z<", this->params._redshift_max_sample());
          vzz.resize(params._Nbins_redshift(),0);
          vrc.resize(params._Nbins_redshift(),0);
          this->cosmology.Comoving_distance_tabulated(this->params._redshift_min_sample(),this->params._redshift_max_sample(), vzz,vrc);
          So.DONE();
        }
        // Now read the galaxy and random catalog
#elif defined (_USE_MASS_BINS_PK_)
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_
                this->add_catalogues(this->params._VMAXbins_min(im),this->params._VMAXbins_max(im));
#elif defined _USE_MASS_AS_OBSERVABLE_POWER_
                this->add_catalogues(pow(10,this->params._MASSbins_min(im)),pow(10,this->params._MASSbins_max(im)));
#endif
#endif
    // If applies, give a first estimate of mean number density from a box
    // in case we do not use random catalog                                                    *
    //If a random catalog is used, set this numer to 1
    //and use the nmbar tabulated in the catalogs or computed in this code
    mean_density=1.0;
    if(false==this->params._use_random_catalog())
        So.message_screen("Using particles in a box.");
    s_data_structure s_data_struct_prop;
    vector<gsl_real> z_v;
    vector<gsl_real> dndz_v;
    bool compute_dndz=false;
    if(true==this->params._use_random_catalog())
      {
#ifdef _USE_SEVERAL_RANDOM_FILES_
        for(int ir=0;ir<this->params._NRANDOMfiles();++ir)
        {
          this->add_catalogues_ran(mm, ir);
#else
        this->add_catalogues_ran(mm);
#endif
// 2
        if(this->params._sys_of_coord_r()>0)
          {
            // The new Lbox is updated in this method into the instance params belonging to the class random_cat
            if(false==this->params._nbar_tabulated() && false==this->params._use_file_nbar())// extra security check
              {
                 compute_dndz=true;
#ifdef _USE_SEVERAL_RANDOM_FILES_
                if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
                 {
                   z_v.resize(params._new_N_dndz_bins(),0);
                   dndz_v.resize(z_v.size(),0);
                  }
#else
                 z_v.resize(params._new_N_dndz_bins(),0);
                 dndz_v.resize(z_v.size(),0);
#endif
                 // First argument of the following method is alpha. Below the dndz_v will be rescaled by the ratio Ngal/Nran once Ngal is read
                 this->random_cat.get_mean_number_density(1.0,vzz,vrc,z_v,dndz_v);
                 this->So.message_warning("In PowerSpectrumF.cpp, variable alpha has been set to one. Must be Ngal/Nran. Line ",__LINE__);
              }
            else if(false==this->params._nbar_tabulated() && true==this->params._use_file_nbar())
              {
#ifdef _USE_SEVERAL_RANDOM_FILES_
                if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
                 {
#endif
          // Read the nbar file provided from the parameter file, done only for the first "ir"
               vector<real_prec> nb;
               ULONG nzl=this->File.read_file(this->params._file_nbar(),nb, 1);
               ULONG ncols_nb=(static_cast<ULONG>(nb.size()/nzl));
               z_v.resize(nzl,0);
               dndz_v.resize(nzl,0);
               real_prec fsky=this->params._area_survey()*pow(M_PI/180.0,2)/(4.0*M_PI);
               So.message_screen("Computed f_sky", fsky);
               for(ULONG i=0;i<nzl;++i)
                {
                  z_v[i]=nb[i*ncols_nb];  // central redshift, assumed to be in the first column of the dNdz file
                  dndz_v[i]=27.62*nb[1+i*ncols_nb];///fsky;//  /(fsky*nb[2+i*ncols_nb]); // mean number density being read here, column 2 is bin volume (area included?)
                // 27.3 es lo que neceisto par allegar al dndz anterior, puede ser el area?, 
                }
               nb.clear(); nb.shrink_to_fit();
#ifdef _USE_SEVERAL_RANDOM_FILES_
                }
#endif
              }
            else if(true==params._nbar_tabulated())
            {
#ifdef _USE_SEVERAL_RANDOM_FILES_
                if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
                 {
#endif
                   z_v.resize(params._new_N_dndz_bins(),0);
#ifdef _USE_SEVERAL_RANDOM_FILES_
                }
#endif
            }
           s_data_struct_prop={
              mean_density,
              compute_dndz,
              z_v,
              dndz_v,
              vzz,
              vrc
            };
            // 3
#ifdef _USE_SEVERAL_RANDOM_FILES_
            this->random_cat.ang_to_cart_coordinates(&s_data_struct_prop,ir);
#else
            this->random_cat.ang_to_cart_coordinates(&s_data_struct_prop);
#endif
            // 4: Interpolate on the mesh
            if(this->params._statistics()=="Pk_y_ds")
              {
                s_data_struct_r_ds={
                  this->random_cat.Halo,
                  this->random_cat._NCOLS(),
                  this->random_cat._type_of_object(),
                  mean_density,
                  compute_dndz,
                  z_v,
                  dndz_v,
                  vzz,
                  vrc
                };
                fftw_functions.get_power_moments_fourier_grid_ds_yam(&s_data_struct_r_ds);
              }
            else
#ifdef _USE_SEVERAL_RANDOM_FILES_
              this->random_cat.get_interpolated_density_field(false,"any",ir);  // Interpolate random catalog on a mesh
#else
             this->random_cat.get_interpolated_density_field(false,"any");  // Interpolate random catalog on a mesh
#endif
            // 5 Update
#ifdef _USE_SEVERAL_RANDOM_FILES_
            if(ir==0)
            {
                //Here we are assuming that *all random files are filling the same volume*
                // in other words, that the randoms are not sorted.
                //so that we can get Lbox from the first one and define the
                //mesh onto whcih the rest will be interpolated! This is important
#endif
                // ############################# This is important ################################################################################################
                // update even if no new_Lbox is requested, as in any case the offsets are to be computed and used by the fftwfunctions' methods (multipoles mainly)
               this->set_params(this->random_cat.params);                // Update the instance parms in this class from the instance params in random_cat
               this->fftw_functions.set_params(this->random_cat.params); //Update the instance parms  in the class fftw_functions from the instance params in random_cat
                // ############################# This is important ################################################################################################
#ifdef _USE_SEVERAL_RANDOM_FILES_
           }
#endif
#ifndef _USE_SEVERAL_RANDOM_FILES_
            // 6 Compute window matrix with the randoms, only if read from one random file
            if(true==this->params._get_window_matrix())
              {
                So.message_screen("Computing mixing matrix");
                this->fftw_functions.set_alpha(1.0);// We set alpha=1 here, though below we update its value
                this->get_window_matrix_multipole();
              }
#endif
            // 7 Release memmory
            this->random_cat.Halo.clear();
            this->random_cat.Halo.shrink_to_fit();
           }// closes if(true==this->params._use_random_catalog() && this->params._sys_of_coord_r()>0)
#ifdef _USE_SEVERAL_RANDOM_FILES_
           this->So.message_screen("\t\tPartial number of randoms = ", get_nobjects(this->random_cat.field_external));
            }
#endif
         } // closes if(true==this->params._use_random_catalog())
#ifdef _USE_SEVERAL_RANDOM_FILES_
           this->So.message_screen("\tTotal number of randoms = ", get_nobjects(this->random_cat.field_external));
#endif
        // *****************************************************************************************
        // Dark matter tracers: Steps are
        // 1 Read catalog
        // 3 Transform to cartessian coord (assigning nbar if needed and searching for box side lenght, offsets and mins)
        // 4 Interpolate on the mesh
        // 5 Update params
        // 6 Release memmory from the catalog if relative bias is not to be computed.
        this->add_catalogues_gal(mm);
        // Read
        // Convert to cartessian coords
#ifdef _USE_SEVERAL_RANDOM_FILES_
        this->tracer_cat.ang_to_cart_coordinates(&s_data_struct_prop,0);
#else
        this->tracer_cat.ang_to_cart_coordinates(&s_data_struct_prop);
#endif
          // Do direct sum or interpolate on a mesh to do fftw
        if(this->params._statistics()=="Pk_y_ds")//will be deprecated
          {
            s_data_structure_direct_sum s_data_struct_g_ds={
              this->tracer_cat.Halo,
              this->tracer_cat._NCOLS(),
              this->tracer_cat._type_of_object(),
              mean_density,
              compute_dndz,
              z_v,
              dndz_v,
              vzz,
              vrc
            };
            fftw_functions.get_power_moments_fourier_grid_ds_yam(&s_data_struct_g_ds);
          }
        else
#ifdef _USE_SEVERAL_RANDOM_FILES_ // we do not need it explcitely. but as we did not create a method only for randoms, we pass it as needed filling with 0 if several files are used
          this->tracer_cat.get_interpolated_density_field(false,"any",0); // Interpolate tracer catalog on a mesh
#else
          this->tracer_cat.get_interpolated_density_field(false,"any"); // Interpolate tracer catalog on a mesh
#endif
        //Update the instance params in  in the class tracer_cat from the instance params in random_cat
        // such that the offsets , and derived parameters are applying the same to the tracer cat
        if(true==this->params._use_random_catalog())
          this->tracer_cat.set_params(this->random_cat.params);
        // Free memmory if necessary. If not, released below
        if(false==this->params._Get_tracer_relative_bias())
          {
            this->tracer_cat.Halo.clear();
            this->tracer_cat.Halo.shrink_to_fit();
          }
        this->fftw_functions.resize_fftw_vectors();
        if(true==verbose)
          fftw_functions.write_fftw_parameters();
        // PASS parameter computed in the class Catalog when interpolating to the fftwFunction class
        this->fftw_functions.set_n_gal(this->tracer_cat._n_gal());
        this->fftw_functions.set_w_g(this->tracer_cat._w_g());
        this->fftw_functions.set_s_g(this->tracer_cat._s_g());
        this->fftw_functions.set_sg1(this->tracer_cat._sg1());
        this->fftw_functions.set_sg2(this->tracer_cat._sg2());
        if(true==this->params._use_random_catalog())
          {
            this->fftw_functions.set_n_ran(this->random_cat._n_gal());
            this->fftw_functions.set_w_r(this->random_cat._w_g());
            this->fftw_functions.set_s_r(this->random_cat._s_g());
            this->fftw_functions.set_sr1(this->random_cat._sg1());
            this->fftw_functions.set_sr2(this->random_cat._sg2());
            this->fftw_functions.set_normal_p(this->random_cat._normal_p());
            this->fftw_functions.set_normal_b(this->random_cat._normal_b());
          }
        this->fftw_functions.get_parameters_estimator(verbose);
#ifdef _USE_OMP_
#pragma omp parallel for 
#endif
        for(ULONG i=0;i<this->params._NGRID();++i)
           this->fftw_functions.set_data_g(i,this->tracer_cat.field_external[i]);// Set rho_r to its container in fftw_functions

        if(false==this->params._use_random_catalog())
          this->fftw_functions.raw_sampling(pow(this->params._Lbox(),3));//This set normalization and SN for box
        else  
          {
        // Here we need to pass the fields from catalog to fftw:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<this->params._NGRID();++i)
            this->fftw_functions.set_data_r(i,this->random_cat.field_external[i]); // Set rho_r to its container in fftw_functions
          }
       // If relative bias is to be computed, do it  and release memmory
        if(true==this->params._Get_tracer_relative_bias())
          {
                // Vamos a generar uan muestra pequeña aleatoria y a esa muestra le asignamos el bias
            this->tracer_cat.select_random_subsample_v(0.001); // Esto generara un .observed=1 si lo tomamos o 0 si no
            this->object_by_object_rbias(); // Assign and print downsampled catalog
            this->tracer_cat.Halo.clear();
            this->tracer_cat.Halo.shrink_to_fit();
          }
          // Free the memmory used for the interpolation of tracer and randoms
        this->tracer_cat.field_external.clear();
        this->tracer_cat.field_external.shrink_to_fit();
        if(true==this->params._use_random_catalog())
          {
            this->random_cat.field_external.clear();
            this->random_cat.field_external.shrink_to_fit();
          }
        this->fftw_functions.get_fluctuation(); // The container fftw_functions.delta_h is generated
      } // CLOSES if this is catalog or field
    else if (this->params._input_type()=="delta_grid" || this->params._input_type()=="density_grid"  ) // Else, we read the delta from this input file.
         {
           So.message_screen("Starting with density field on a grid");
           real_prec ngal_new;
           bool measure_diff=false;
           if(false==this->params._measure_cross() && false==measure_diff) // If no crossed power, read per default theparams._delta_grid_file()
            {
             vector<real_prec> field(this->params._NGRID());
             switch(this->params._measure_cross_from_1())
              { //choose the file to get the auto power from
                case(1):
                this->File.read_array_fast(this->params._delta_grid_file(), field);
                break;
                case(2):
                this->File.read_array_fast(this->params._delta_grid_file2(), field);
                break;
                case(3):
                this->File.read_array_fast(this->params._delta_grid_file3(), field);
                break;
                case(4):
                this->File.read_array_fast(this->params._delta_grid_file4(), field);
                break;
             }
             fftw_functions.data_g.clear();
             fftw_functions.data_g.shrink_to_fit();
             fftw_functions.data_g.resize(this->params._NGRID(),0);
             if (this->params._input_type()=="delta_grid")
              {
#ifdef _USE_OMP_
#pragma omp parallel
#endif
                for(ULONG i=0;i<this->params._NGRID();++i)
                  fftw_functions.data_g[i]=static_cast<real_prec>(field[i]);
                ngal_new=this->params._ngal_delta();
              }
            else if (this->params._input_type()=="density_grid")
              {
                ngal_new=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
                for(ULONG i=0;i<this->params._NGRID();++i)
                 ngal_new+=static_cast<real_prec>(field[i]);
                real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                for(ULONG i=0;i<this->params._NGRID();++i)
                  fftw_functions.data_g[i]=(static_cast<real_prec>(field[i])/static_cast<real_prec>(nmean))-1.;
                this->params.set_ngal_delta(ngal_new);
                fftw_functions.set_n_gal(ngal_new);
             }
            real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
            fftw_functions.set_normal_power(pow(factor,-2));
            fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
            if(true==verbose)
              {
                fftw_functions.write_fftw_parameters();
#ifdef _FULL_VERBOSE_
                So.message_screen("Shot Noise =",fftw_functions.shot_noise);
#endif
              }
              fftw_functions.resize_fftw_vectors();
          }
        else if(true==measure_diff)
          {
            vector<float> field(this->params._NGRID());
            vector<float> field2(this->params._NGRID());
            this->File.read_array_fast(this->params._delta_grid_file(), field);
            this->File.read_array_fast(this->params._delta_grid_file2(), field2);
            this->fftw_functions.data_g.resize(this->params._NGRID(),0);
            ULONG ngal_new2=0;
            ULONG ngal_new1=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new1)
#endif
            for(ULONG i=0;i<this->params._NGRID();++i)
              ngal_new1+=static_cast<real_prec>(field[i]);
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new2)
#endif
            for(ULONG i=0;i<this->params._NGRID();++i)
              ngal_new2+=static_cast<real_prec>(field2[i]);
            real_prec nmean1=static_cast<real_prec>(ngal_new1)/static_cast<real_prec>(this->params._NGRID());
            real_prec nmean2=static_cast<real_prec>(ngal_new2)/static_cast<real_prec>(this->params._NGRID());
                    fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new1);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
            for(ULONG i=0;i<this->params._NGRID();++i)
                      //	    fftw_functions.data_g[i]=(static_cast<real_prec>(field[i])/static_cast<real_prec>(nmean1) - static_cast<real_prec>(field2[i])/static_cast<real_prec>(nmean2));///fftw_functions.shot_noise;
                      fftw_functions.data_g[i]=(static_cast<real_prec>(field[i]) - static_cast<real_prec>(field2[i])/static_cast<real_prec>(nmean2));///fftw_functions.shot_noise;
            this->params.set_ngal_delta(ngal_new1);
                    fftw_functions.set_n_gal(ngal_new1);
                    real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
                    fftw_functions.set_normal_power(pow(factor,-2));
                    fftw_functions.shot_noise=0;
            if(true==verbose)
                      fftw_functions.write_fftw_parameters();
                    fftw_functions.resize_fftw_vectors();
          }
        else  if(true==this->params._measure_cross()) // If no crossed power, read per default theparams._delta_grid_file()  // if we measure the cross, then
          {
            fftw_functions.data_g.resize(this->params._NGRID(),0);
            vector<float> field(this->params._NGRID());
            switch(this->params._measure_cross_from_1())
            { //choose the file to get the auto power from
              case(1):
                this->File.read_array_fast(this->params._delta_grid_file(), field);
                break;
              case(2):
                this->File.read_array_fast(this->params._delta_grid_file2(), field);
                break;
              case(3):
                this->File.read_array_fast(this->params._delta_grid_file3(), field);
                break;
              case(4):
                this->File.read_array_fast(this->params._delta_grid_file4(), field);
                break;
            }
          if(this->params._input_type()=="delta_grid")
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                  for(ULONG i=0;i<this->params._NGRID();++i)
                    fftw_functions.data_g[i]=static_cast<real_prec>(field[i]);
                  ngal_new=this->params._ngal_delta();
           }
          else if (this->params._input_type()=="density_grid")
            {
              ngal_new=0;
              for(ULONG i=0;i<this->params._NGRID();++i)ngal_new+=static_cast<real_prec>(field[i]);
              real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
    #ifdef _USE_OMP_
#pragma omp parallel for
#endif
              for(ULONG i=0;i<this->params._NGRID();++i)
                fftw_functions.data_g[i]=(static_cast<real_prec>(field[i])/static_cast<real_prec>(nmean))-1.0;
            }
          fftw_functions.shot_noise=0;
          if(true==this->params._SN_correction())
            fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
          So.message_screen("Shot Noise =",fftw_functions.shot_noise);
          switch(this->params._measure_cross_from_2())
            { //choose the file to get the auto power from
              case(1):
                this->File.read_array_fast(this->params._delta_grid_file(), field);
              break;
              case(2):
                this->File.read_array_fast(this->params._delta_grid_file2(), field);
              break;
              case(3):
                this->File.read_array_fast(this->params._delta_grid_file3(), field);
              break;
              case(4):
                this->File.read_array_fast(this->params._delta_grid_file4(), field);
              break;
            }
          fftw_functions.data_gp.resize(this->params._NGRID(),0);
          if (this->params._input_type()=="delta_grid")
            {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                for(ULONG i=0;i<this->params._NGRID();++i)
                  fftw_functions.data_gp[i]=static_cast<real_prec>(field[i]);
                ngal_new=this->params._ngal_delta();
            }
            else if (this->params._input_type()=="density_grid")
              {
                ngal_new=0;
                            //                  for(ULONG i=0;i<this->params._NGRID();++i)fftw_functions.data_gp[i]=static_cast<real_prec>(field[i]);
                this->So.message_warning("In PowerSpectrumF.cpp I have manually set that the field2 is already a delta. This can be solved if I define a new variable input_type for the second field in the cross power");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                 for(ULONG i=0;i<this->params._NGRID();++i)
                  fftw_functions.data_gp[i]=static_cast<real_prec>(field[i]);
              }
            // Here we find 1+delta in the grid
            real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
            fftw_functions.set_normal_power(pow(factor, -2));
            //IF IT IS DM
            fftw_functions.shot_noise2=0;//static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
            if(true==this->params._SN_correction())
             {
                 So.message_screen("Shot noise 1 =", fftw_functions.shot_noise);
                 So.message_screen("Shot noise 2 =", fftw_functions.shot_noise2);
              }
            else
              So.message_screen("No SN correction");

            if(true==verbose)
               fftw_functions.write_fftw_parameters();
            fftw_functions.resize_fftw_vectors();
          }
        } // closes else if (this->params._input_type()=="delta_grid" || this->params._input_type()=="density_grid"  ) 
        // *****************************************************************************************
        // WELCOME TO FOURIER SPACE
        if(this->params._statistics()=="Pk_fkp" || this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_yb"  || this->params._statistics()=="Pk_ybc"   || this->params._statistics()=="Pk_ysc" ||  this->params._statistics()=="Pk_y_ds")
          {
            kvector_data.clear();
            kvector_data.shrink_to_fit();
            kvector_window.clear();
            kvector_window.shrink_to_fit();
            if("linear" == this->params._type_of_binning())
            {
              for(int i=0;i<this->params._d_Nnp_data();i++)
                this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
              for(int i=0;i<this->params._d_Nnp_window();i++)
                this->kvector_window.push_back(this->params._d_kmin()+this->params._d_DeltaK_window()*(i+0.5));
            }
            else if("log"==this->params._type_of_binning())
              {
                for(int i=0;i<kvector_data.size();i++)
                  this->kvector_data.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
                for(int i=0;i<kvector_window.size();i++)
                this->kvector_window.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
              }
            for(int i=0;i<this->params._d_Nnp_data();i++)
              kvector_data2d.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
            for(int i=0;i<params._N_mu_bins();i++)
              muvector.push_back(-1.0+this->params._d_Deltamu()*(i+0.5));
            // *****************************************************************************
            // Resize arrays for P(k), and 2d P(k). Compute and write to file
            this->pk0.clear();
            this->pk0.shrink_to_fit();
            this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
            //#ifdef _WRITE_MULTIPOLES_
            this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
            this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole
            //#endif
            this->pk_w.clear();
            this->pk_w.shrink_to_fit();
            this->pk_w.resize(this->params._d_Nnp_window(),0); //W(k)
            this->modes_g.clear();
            this->modes_g.shrink_to_fit();
            this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
            //#ifdef _WRITE_2DPOWER_
            this->pkk.resize(this->params._d_Nnp_data());
            this->pmk.resize(params._N_mu_bins());
            for(int i=0;i<this->params._d_Nnp_data();i++)this->pkk[i].resize(this->params._d_Nnp_data(),0);
            for(int i=0;i<params._N_mu_bins();i++)this->pmk[i].resize(this->params._d_Nnp_data(),0);
            //#endif
            this->sigma_fkp.clear();
            this->sigma_fkp.shrink_to_fit();
            this->sigma_fkp.resize(this->params._d_Nnp_data(),0);
            // ****************************************************************************
            // Get power spectrum and more
            if(this->params._statistics()=="Pk_fkp")
              {
                fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
                sigma_y_l2.resize(this->params._d_Nnp_data(),0);
                sigma_y_l4.resize(this->params._d_Nnp_data(),0);
                if(true==params._FKP_error_bars())
                  {
                    So.message("Computing FKP error bars");
                    this->fftw_functions.get_fkp_error_bars(&s_data_struct_r_ds, kvector_data, this->pk0, this->modes_g, this->sigma_fkp);
                  }
              }
            else if(this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc" || this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_y_ds" || this->params._statistics()=="Pk_ysc" )
              this->fftw_functions.get_power_spectrum_yamamoto(this->pk0,this->pk2,this->pk4,this->modes_g);
            //MISSINGN ERROR BARS FROM YAMAMOTO HERE.
    #ifdef _USE_GNUPLOT_POWER_
              this->gp.plot_power_spectrum_redshift_space(this->kvector_data,this->pk0,this->pk2);
    #endif
            }
            // Estimates of Bispectrum. Using the DFT already done for P(k)
            else if(this->params._statistics()=="Bk_fkp")
              {
            if(this->params._type_of_binning()=="linear")
              for(int i=0;i<this->params._d_Nnp_data();i++)
                this->kvector_data_b.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5)); //Oficcial binning
            else
              if(this->params._type_of_binning()=="log"){
                for(int i=0;i<this->params._d_Nnp_data();i++)
                  this->kvector_data_b.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
              }
            bispectrum.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
            sn_bispectrum.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
            modes_tri.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
            this->fftw_functions.get_bispectrum_fkp('d', bispectrum, sn_bispectrum, modes_tri);
            File.write_to_file(file_bispectrum,kvector_data_b,bispectrum,modes_tri);
              }
            // Estimates of Bispectrum for FKP using fast version
            else if(this->params._statistics()=="Bk_fkp_fast")
              {
            this->pk0.resize(this->params._d_Nnp_data(),0);
            for(int i=0;i<fftw_functions.Nshells_bk;i++)
              kvector_data_b.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5)); //Jennifer's binning
            bispectrum.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
            sn_bispectrum.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
            modes_tri.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
            fftw_functions.get_power_spectrum_for_bispectrum(this->pk0);
            fftw_functions.get_bispectrum_fkp_fast(this->pk0,bispectrum,modes_tri,file_bispectrum);
          }
#ifndef _WRITE_MULTIPOLES_
        write_power_and_modes();
#else
        write_power_spectrum();
#endif
        this->write_fftw_parameters();
   }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_power_spectrum(){
    string file_p=this->file_power;
    string file_lg=this->file_power_log;
    for(int bin1=0; bin1<3; bin1++) // Bin in redshift
    for(int bin2=0; bin2<3; bin2++) // Bin in color
    for(int bin3=0; bin3<3; bin3++) // Bin in log stellar mass
	  {
	    // Just in case, will be deprecated
	    real_prec mean_density=1.0;
	    // Initizalize number of intervals
	    int N_intervals=1;
       if(I_EQZ==this->params._sys_of_coord_g())
	      So.write_cosmo_parameters((void *)&this->params.s_cosmo_pars);
	    s_data_structure s_data_struct_r; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
	    s_data_structure_direct_sum s_data_struct_r_ds; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
	    this->File.input_type=this->params._input_type();
	    this->file_power=file_p+"_zbin"+to_string(bin1)+"_cbin"+to_string(bin2)+"_msbin"+to_string(bin3);
	    this->file_power_log=file_lg+"_zbin"+to_string(bin1)+"_cbin"+to_string(bin2)+"_msbin"+to_string(bin3);
	    // **************************************************************************
	    vector<gsl_real>vzz;
	    vector<gsl_real>vrc;
	    So.message_screen("Computing comoving distances in the range",this->params._redshift_min_sample(),"<z<", this->params._redshift_max_sample());
	    vzz.resize(params._Nbins_redshift(),0);
	    vrc.resize(params._Nbins_redshift(),0);
	    cosmology.Comoving_distance_tabulated(this->params._redshift_min_sample(),this->params._redshift_max_sample(), vzz,vrc);
	    So.DONE();
	    // If applies, give a first estimate of mean number density from a box
	    // in case we do not use random catalog                                                    *
	    //If a random catalog is used, set this numer to 1
	    //and use the nmbar tabulated in the catalogs or computed in this code
	    // *****************************************************************************************
	    // *****************************************************************************************
	    // *****************************************************************************************
	    // Randoms: Steps are
	    // 1 Read catalog
	    // 2 Get ready for nbar
	    // 3 Transform to cartessian coord (assigning nbar if needed and searching for box side lenght, offsets and mins)
	    // 4 Interpolate on the mesh
	    // 5 Update params
	    // 6 Get window matrix if requested
	    // 7 Release memmory from the catalog
	    // If random file is too large, a loop over the chunks is done,
	    s_data_structure s_data_struct_prop;
	    vector<gsl_real> z_v;
	    vector<gsl_real> dndz_v;
	    bool compute_dndz=false;
#ifdef _USE_SEVERAL_RANDOM_FILES_
	    for(int ir=0;ir<this->params._NRANDOMfiles();++ir)
	      {
		this->add_catalogues_ran(ir, bin1, bin2, bin3);
#else
		this->add_catalogues_ran(bin1, bin2, bin3);
#endif
		// 2
		// The new Lbox is updated in this method into the instance params belonging to the class random_cat
		if(false==this->params._nbar_tabulated() && false==this->params._use_file_nbar())// extra security check
		  {
		    compute_dndz=true;
#ifdef _USE_SEVERAL_RANDOM_FILES_
		    if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
		      {
			z_v.resize(params._new_N_dndz_bins(),0);
			dndz_v.resize(z_v.size(),0);
		      }
#else
		    z_v.resize(params._new_N_dndz_bins(),0);
		    dndz_v.resize(z_v.size(),0);
#endif
		    // First argument of the following method is alpha. Below the dndz_v will be rescaled by the ratio Ngal/Nran once Ngal is read
		    this->random_cat.get_mean_number_density(1.0,vzz,vrc,z_v,dndz_v);
		    this->So.message_warning("In PowerSpectrumF.cpp, variable alpha has been set to one. Must be Ngal/Nran. Line ",__LINE__);
		  }
		else if(false==this->params._nbar_tabulated() && true==this->params._use_file_nbar())
		  {
#ifdef _USE_SEVERAL_RANDOM_FILES_
		    if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
		      {
#endif
			// Read the nbar file provided from the parameter file, done only for the first "ir"
			vector<real_prec> nb;
			ULONG nzl=this->File.read_file(this->params._file_nbar(),nb, 1);
			ULONG ncols_nb=(static_cast<ULONG>(nb.size()/nzl));
			z_v.resize(nzl,0);
			dndz_v.resize(nzl,0);
			for(ULONG i=0;i<nzl;++i)
			  {
			    z_v[i]=nb[i*ncols_nb];  // redshift, assumed to be in the first column of the dNdz file
			    dndz_v[i]=nb[1+i*ncols_nb]; // mean number density
			  }
			nb.clear(); nb.shrink_to_fit();
#ifdef _USE_SEVERAL_RANDOM_FILES_
		      }
#endif
		  }
		else if(true==params._nbar_tabulated())
		  {
#ifdef _USE_SEVERAL_RANDOM_FILES_
		    if(0==ir)// resize only when reading the first file, For the others, we just keep on filling it.
		      {
#endif
			z_v.resize(params._new_N_dndz_bins(),0);
#ifdef _USE_SEVERAL_RANDOM_FILES_
		      }
#endif
		  }
		s_data_struct_prop={
		  mean_density,
		  compute_dndz,
		  z_v,
		  dndz_v,
		  vzz,
		  vrc
		};
	// 3
#ifdef _USE_SEVERAL_RANDOM_FILES_
		this->random_cat.ang_to_cart_coordinates(&s_data_struct_prop,ir);
#else
		this->random_cat.ang_to_cart_coordinates(&s_data_struct_prop);
#endif
		
		// 4: Interpolate on the mesh
#ifdef _USE_SEVERAL_RANDOM_FILES_
		this->random_cat.get_interpolated_density_field(false,"any",ir);  // Interpolate random catalog on a mesh
#else
		this->random_cat.get_interpolated_density_field(false,"any");  // Interpolate random catalog on a mesh
#endif
		
		// 5 Update
#ifdef _USE_SEVERAL_RANDOM_FILES_
		if(ir==0)
		  {
		    //Here we are assuming that *all random files are filling the same volume*
		    // in other words, that the randoms are not sorted.
		    //so that we can get Lbox from the first one and define the
		    //mesh onto whcih the rest will be interpolated! This is important
#endif
		    // *****************This is important**********************
		    // update even if no new_Lbox is requested, as in any case the offsets are to be computed and used by the fftwfunctions' methods (multipoles mainly)
		    this->set_params(this->random_cat.params);                // Update the instance parms in this class from the instance params in random_cat
		    this->fftw_functions.set_params(this->random_cat.params); //Update the instance parms  in the class fftw_functions from the instance params in random_cat
#ifdef _USE_SEVERAL_RANDOM_FILES_
		  }
#endif
		// 7 Release memmory
		this->random_cat.Halo.clear();
		this->random_cat.Halo.shrink_to_fit();
		#ifdef _USE_SEVERAL_RANDOM_FILES_
		this->So.message_screen("\t\tPartial number of randoms = ", get_nobjects(this->random_cat.field_external));
          }// closes loop over random files
#endif 
#ifdef _USE_SEVERAL_RANDOM_FILES_
	    this->So.message_screen("\tTotal number of randoms = ", get_nobjects(this->random_cat.field_external));
#endif
	    // *****************************************************************************************
	    // Dark matter tracers: Steps are
	    // 1 Read catalog
	    // 3 Transform to cartessian coord (assigning nbar if needed and searching for box side lenght, offsets and mins)
	    // 4 Interpolate on the mesh
	    // 5 Update params
	    // 6 Release memmory from the catalog if relative bias is not to be computed.
	    this->add_catalogues_gal(bin1, bin2, bin3);
	    // Read
	    // Convert to cartessian coords
#ifdef _USE_SEVERAL_RANDOM_FILES_
	    this->tracer_cat.ang_to_cart_coordinates(&s_data_struct_prop,0);
#else
	    this->tracer_cat.ang_to_cart_coordinates(&s_data_struct_prop);
#endif
    // Do direct sum or interpolate on a mesh to do fftw
 #ifdef _USE_SEVERAL_RANDOM_FILES_ // we do not need it explcitely. but as we did not create a method only for randoms, we pass it as needed filling with 0 if several files are used
            this->tracer_cat.get_interpolated_density_field(false,"any",0); // Interpolate tracer catalog on a mesh
#else
            this->tracer_cat.get_interpolated_density_field(false,"any"); // Interpolate tracer catalog on a mesh
#endif
	    //Update the instance params in  in the class tracer_cat from the instance params in random_cat
	    // such that the offsets , and derived parameters are applying the same to the tracer cat
	    this->tracer_cat.set_params(this->random_cat.params);
	    // Free memmory if necessary. If not, released below
	    if(false==this->params._Get_tracer_relative_bias())
	      {
            this->tracer_cat.Halo.clear();
            this->tracer_cat.Halo.shrink_to_fit();
	      }
	    // ***********************************************************************************
	    this->fftw_functions.resize_fftw_vectors();
	    this->fftw_functions.write_fftw_parameters();
	    // PASS parameter computed in the class Catalog when interpolating to the fftwFunction class
        this->fftw_functions.set_n_gal(this->tracer_cat._n_gal());
        this->fftw_functions.set_w_g(this->tracer_cat._w_g());
        this->fftw_functions.set_s_g(this->tracer_cat._s_g());
        this->fftw_functions.set_sg1(this->tracer_cat._sg1());
        this->fftw_functions.set_sg2(this->tracer_cat._sg2());
        this->fftw_functions.set_n_ran(this->random_cat._n_gal());
        this->fftw_functions.set_w_r(this->random_cat._w_g());
        this->fftw_functions.set_s_r(this->random_cat._s_g());
        this->fftw_functions.set_sr1(this->random_cat._sg1());
        this->fftw_functions.set_sr2(this->random_cat._sg2());
        this->fftw_functions.set_normal_p(this->random_cat._normal_p());
        this->fftw_functions.set_normal_b(this->random_cat._normal_b());
        this->fftw_functions.get_parameters_estimator(true);
	    // Here we need to pass the fields from catalog to fftw:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
	    for(ULONG i=0;i<this->params._NGRID();++i)
	      this->fftw_functions.set_data_g(i,this->tracer_cat.field_external[i]);// Set rho_r to its container in fftw_functions

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
	    for(ULONG i=0;i<this->params._NGRID();++i)
	      this->fftw_functions.set_data_r(i,this->random_cat.field_external[i]); // Set rho_r to its container in fftw_functions
	    // If relative bias is to be computed, do it  and release memmory
	    if(true==this->params._Get_tracer_relative_bias())
	      {
		// Vamos a generar uan muestra pequeña aleatoria y a esa muestra le asignamos el bias
            this->tracer_cat.select_random_subsample_v(0.001); // Esto generara un .observed=1 si lo tomamos o 0 si no
            this->object_by_object_rbias(); // Assign and print downsampled catalog
            this->tracer_cat.Halo.clear();
            this->tracer_cat.Halo.shrink_to_fit();
	      }
	    // Free the memmory used for the interpolation of tracer and randoms
	    this->tracer_cat.field_external.clear();
	    this->tracer_cat.field_external.shrink_to_fit();
	    this->random_cat.field_external.clear();
	    this->random_cat.field_external.shrink_to_fit();
	    // WELCOME TO FOURIER SPACE
        //Get fluctuation
	    this->fftw_functions.get_fluctuation(); // The container fftw_functions.delta_h is generated
        this->fftw_functions.resize_fftw_vectors();
        this->kvector_data.clear();
	    this->kvector_data.shrink_to_fit();
	    this->kvector_window.clear();
	    this->kvector_window.shrink_to_fit();
	    this->pk0.clear();
	    this->pk0.shrink_to_fit();
	    this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
	    this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
	    this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole
	    this->pk_w.clear();
	    this->pk_w.shrink_to_fit();
	    this->pk_w.resize(this->params._d_Nnp_window(),0); //W(k)
	    this->modes_g.clear();
	    this->modes_g.shrink_to_fit();
	    this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
	    this->pkk.resize(this->params._d_Nnp_data());
	    this->pmk.resize(params._N_mu_bins());
	    if("linear" == this->params._type_of_binning())
	      {
		for(int i=0;i<this->params._d_Nnp_data();i++)
		  this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
		for(int i=0;i<this->params._d_Nnp_window();i++)
		  this->kvector_window.push_back(this->params._d_kmin()+this->params._d_DeltaK_window()*(i+0.5));
	      }
	    else if("log"==this->params._type_of_binning())
	      {
		for(int i=0;i<kvector_data.size();i++)
		  this->kvector_data.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
		for(int i=0;i<kvector_window.size();i++)
		  this->kvector_window.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
	      }
	    for(int i=0;i<this->params._d_Nnp_data();i++)
	      kvector_data2d.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
	    for(int i=0;i<params._N_mu_bins();i++)
	      muvector.push_back(-1.0+this->params._d_Deltamu()*(i+0.5));
	    for(int i=0;i<this->params._d_Nnp_data();i++)
	      this->pkk[i].resize(this->params._d_Nnp_data(),0);
	    for(int i=0;i<params._N_mu_bins();i++)
	      this->pmk[i].resize(this->params._d_Nnp_data(),0);
	    if(this->params._statistics()=="Pk_fkp")
	      this->fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
	    else
	      this->fftw_functions.get_power_spectrum_yamamoto(this->pk0,this->pk2,this->pk4,this->modes_g);
#ifdef _USE_GNUPLOT_POWER_
	    this->gp.plot_power_spectrum_redshift_space(this->kvector_data,this->pk0,this->pk2);
#endif
	    this->write_power_spectrum();
	    this->write_fftw_parameters();
	  } // close loop over bins
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_window_function(){
    time_t start;
    time (&start);
#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
    string file_pow=this->file_power;
#ifdef _FULL_VERBOSE_
      if(this->params._statistics()=="Pk_fkp")So.welcome_message();
       if(this->params._statistics()=="Bk_fkp")So.welcome_message_bispectrum();
       if(this->params._statistics()=="Bk_fkp_fast")So.welcome_message_bispectrum_fast();
      if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc" || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_y_ds")So.welcome_message_yama();
#endif
        if(2==this->params._sys_of_coord_g())
          So.write_cosmo_parameters((void *)&this->params.s_cosmo_pars);
        s_parameters_box s_p_box; // This will be deprecated as we pas the params object to FFTFunctions directly.
        this->File.input_type=this->params._input_type();
        this->add_random_catalogue();
        vector<gsl_real>vzz(params._Nbins_redshift(),0);
        vector<gsl_real>vrc(params._Nbins_redshift(),0);
        So.message_screen("Computing comoving distances");
        if(this->params._sys_of_coord_g()==1 ||this->params._sys_of_coord_g()==2 )
          cosmology.Comoving_distance_tabulated(this->params._redshift_min_sample(),this->params._redshift_max_sample(), vzz,vrc);
        So.DONE();
        // Give a first estimate of the alpha-parameter                                            *
        real_prec alpha_0=1.0;
        if(true==this->params._use_random_catalog())
          alpha_0 = ((real_prec)this->N_galaxy)/((real_prec)this->N_random);
        // This should only be defined if we do not use randoms with nbar tabulated.
        vector<gsl_real> z_v;
        vector<gsl_real> dndz_v;
        vector< vector<gsl_real> > dndz_matrix;
        // If nbar is not tabulated, Compute a smoothed version of dN/dz from randoms to get nbar  *
        bool compute_dndz=false;
        if(true==this->params._use_random_catalog() && false==params._nbar_tabulated())
          {
            compute_dndz=true;
            vector<gsl_real> z_v;
            z_v.resize(params._new_N_dndz_bins(),0);
            vector<gsl_real> dndz_v;
            dndz_v.resize(z_v.size(),0);
           this->random_cat.get_mean_number_density(alpha_0,vzz,vrc,z_v,dndz_v);
          }
        // Allocate structure with random catalogue and information of dN/dz                       *
        So.message_screen("Setting data in structure for RANDOMS");
        s_data_structure s_data_struct_r={
          mean_density,
          compute_dndz,
          z_v,
          dndz_v,
          vzz,
         vrc
        };
        So.DONE();
        s_p_box.npixels=1;
        s_p_box.nside=1;
        s_p_box.file_dndz=file_dndz;
        // Estimate of the mean number density
        real_prec Lside=this->params._Lbox();
        mean_density=1.0;
    // Transforming to cartessian coord.- and searching for box side lenght.
        if( this->params._sys_of_coord_r()!=0)
            {
            So.message_screen("Transform to cartesian coordinates in random catalogue ");
#ifdef _USE_SEVERAL_RANDOM_FILES_  // this is accidental, we only use one for the window
            this->random_cat.ang_to_cart_coordinates(&s_data_struct_r,0);
#else
            this->random_cat.ang_to_cart_coordinates(&s_data_struct_r);
#endif
            if(true==this->params._new_Lbox())//if this is true, the params need to be updated
            {
              this->set_params(this->random_cat.params); // Update the instance parms in this class from the instance params in random_cat
              this->fftw_functions.set_params(this->random_cat.params); //Update the instance parms  in the class fftw_functions from the instance params in random_cat
              this->tracer_cat.set_params(this->random_cat.params); //Update the instance params in  in the class tracer_cat from the instance params in random_cat to the one
            }
      }
     this->fftw_functions.resize_fftw_vectors();
#ifdef _FULL_VERBOSE_
        fftw_functions.write_fftw_parameters();
#endif
    if(true==this->params._nbar_tabulated())
      {
#ifdef _FULL_VERBOSE_
            So.message_screen("Computing mixing matrix");
#endif
        this->get_window_matrix_multipole();
     }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::compute_power_spectrum_grid(const vector<real_prec> &data_in, bool write_to_file)
 {
   this->So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   this->fftw_functions.resize_fftw_vectors();
   real_prec ngal_new,nmean;
   real_prec vol = static_cast<real_prec>(pow(this->params._Lbox(),3));
   // if delta  = rho - mean, (with Ngal given from par file) use normal=Ngal²/Vol (i.e, constructing delta from a catalog)
   // if delta = rho/mean -1, use normal = Ncells²/Vol (i.e, givng delta from outside.
   if(this->params._input_type()=="density_grid")
     {
       ngal_new=get_nobjects(data_in);
       this->params.set_ngal_delta(ngal_new);
       this->fftw_functions.set_n_gal(ngal_new);
       nmean=static_cast<double>(ngal_new)/static_cast<double>(this->params._NGRID());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<this->params._NGRID();++i)
         this->fftw_functions.data_g[i]=static_cast<real_prec>(data_in[i])-static_cast<real_prec>(nmean);
       this->fftw_functions.set_normal_power(static_cast<double>((static_cast<double>(ngal_new)/static_cast<double>(vol))*static_cast<double>(ngal_new)));
     }
   else if(this->params._input_type()=="delta_grid")
     {
       ngal_new=1.0;
       this->params.set_ngal_delta(ngal_new);
       this->fftw_functions.set_n_gal(ngal_new);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<this->params._NGRID();++i)
         this->fftw_functions.data_g[i]=data_in[i];
       this->fftw_functions.set_normal_power(static_cast<real_prec>(this->params._NGRID())*static_cast<real_prec>(this->params._NGRID())/static_cast<real_prec>(vol));
       So.message_screen("Mean of inout field", get_mean(data_in));
     }
#ifdef _FULL_VERBOSE_POWER_
   if(ngal_new>1)
     So.message_screen("Number of tracers in input density field = ",ngal_new);
   So.message_screen("Normalization of power = ",fftw_functions._normal_power());
#endif
   this->fftw_functions.set_shot_noise(vol/static_cast<real_prec>(ngal_new));
#ifdef _FULL_VERBOSE_POWER_
   if(true==this->params._SN_correction())
     So.message_screen("Poisson Shot noise = ",fftw_functions.shot_noise);
   else
     So.message_screen("Poisson Shot noise (not applied) = ",fftw_functions.shot_noise);
#endif
   kvector_data.resize(this->params._d_Nnp_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->params._d_Nnp_data();i++)
     this->kvector_data[i]=this->params._d_kmin()+(i+0.5)*this->params._d_DeltaK_data();
   this->pk0.resize(this->params._d_Nnp_data(),0);
   this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
   fftw_functions.power_spectrum(this->pk0,this->modes_g);
   if(true==write_to_file)
     this->write_power_and_modes();
#ifdef _USE_GNUPLOT_POWER_
   this->gp.plot_power_spectrum(this->kvector_data,this->pk0);
#endif
 }//closes memmber function
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::compute_power_spectrum(bool verbose, vector<s_Halo>& tracer_cat){
      this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
      int NTHREADS=_NTHREADS_;
      omp_set_num_threads(NTHREADS);
#endif
#ifdef _MASS_WEIGHT_POWER_
      if(true==this->params._weight_with_mass())
    So.message_screen("Measuring the Mass-power spectrum");
#elif defined(_USE_MASS_CUTS_PK_)
      So.message_screen("Measuring power spectrum in mass cuts");
#endif
      string file_pow=this->file_power;
#ifdef _USE_MASS_CUTS_PK_
      vector<real_prec>mass_cuts;

#ifdef _SET_GLOBAL_MASS_CUT_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      So.message_screen("Using one global Vmax-cut at Vmax = ",MINIMUM_PROP_CUT," km/s");
#elif defined _USE_MASS_AS_OBSERVABLE_
      So.message_screen("Using one global Halo Mass-cut at M = ",MINIMUM_PROP_CUT, "Ms/h");
#endif
#else
      So.message_warning("Masss-cuts are defined in PowerSpectrumF::compute_power_spectrum");
#endif
      mass_cuts.clear();
      bool mcut=true;
      if(true==this->params._weight_with_mass())
    mcut=false;
   if(true==mcut)
    {
#ifdef _SET_GLOBAL_MASS_CUT_
      mass_cuts.push_back(MINIMUM_PROP_CUT);
#endif
    }
      else
    {
      mass_cuts.push_back(pow(10,this->params._LOGMASSmin()*this->params._MASS_units()));
    }
      for(int im=0; im< mass_cuts.size();++im)
    {
#endif
#ifdef _USE_MASS_CUTS_PK_
      this->file_power=file_pow+"_masscut"+to_string(im);
#ifdef _SET_GLOBAL_MASS_CUT_
      this->file_power=file_pow+"_global_cut";
#endif
#endif
      if(true==verbose)
        {
          if(this->params._statistics()=="Pk_fkp")So.welcome_message();
          if(this->params._statistics()=="Bk_fkp")So.welcome_message_bispectrum();
          if(this->params._statistics()=="Bk_fkp_fast")So.welcome_message_bispectrum_fast();
          if(this->params._statistics()=="Pk_ys" || this->params._statistics()=="Pk_yb" || this->params._statistics()=="Pk_ybc" || this->params._statistics()=="Pk_ysc" || this->params._statistics()=="Pk_y_ds")So.welcome_message_yama();
        }
      if( this->params._sys_of_coord_g()==2)
        So.write_cosmo_parameters((void *)&this->params.s_cosmo_pars);
      // Define some structures here
      s_parameters_box s_p_box;
      s_data_structure_direct_sum s_data_struct_r_ds; // this is defiend here for the fkp error bars need the info from the randoms. So far at this point it is only defiend but not filled
      // Add random and galaxy catalogues
      this->File.input_type=this->params._input_type();
      this->N_galaxy=tracer_cat.size();
      mean_density=static_cast<real_prec>(this->N_galaxy)/pow(this->params._Lbox(),3);
      So.message_screen("Attention: Using particles in a box.");
      So.message_screen("Mean Number density = ",mean_density," (Mpc/h)^(-3)");
      // Give a first estimate of the alpha-parameter                                            *
      real_prec alpha_0=1.0;
      fftw_functions.resize_fftw_vectors();
      fftw_functions.write_fftw_parameters();
      // Build interpolated galaxy density field.
      // The direct sum is spetial , for we alrady mix there the interpolation and wavenumber modes, so
      //I propose to build a structure lile s_data_struct to be filled only if the direct sum is to be used.
#ifdef _USE_SEVERAL_RANDOM_FILES_
      this->tracer_cat.get_interpolated_density_field(false, "any",0);
#else
      this->tracer_cat.get_interpolated_density_field(false, "any");
#endif
     // Here we need to pass the fields from catalog to fftw:
       for(ULONG i=0;i<this->params._NGRID();++i)
          this->fftw_functions.set_data_g(i,this->tracer_cat.field_external[i]);
       this->tracer_cat.field_external.clear();this->tracer_cat.field_external.shrink_to_fit();
      // Build interpolated random density field
      real_prec vol=pow(this->params._Lbox(),3);
      fftw_functions.raw_sampling(vol);
      fftw_functions.get_parameters_estimator( verbose);
      // Build fluctuation                                                                       *
      fftw_functions.get_fluctuation();
      // WELCOME TO FOURIER SPACE
       this->kvector_data.clear();
       this->kvector_data.shrink_to_fit();
       this->kvector_window.clear();
       this->kvector_window.shrink_to_fit();
       if(this->params._statistics()=="Pk_fkp")
         {
         if(this->params._type_of_binning()=="linear")
           {
            for(int i=0;i<this->params._d_Nnp_data();i++)
            kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
            for(int i=0;i<this->params._d_Nnp_window();i++)
             kvector_window.push_back(this->params._d_kmin()+fftw_functions.DeltaK_window*(i+0.5));
          }
          else  if(this->params._type_of_binning()=="log")
          {
            for(int i=0;i<kvector_data.size();i++)
              kvector_data.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
            for(int i=0;i<kvector_window.size();i++)
             kvector_window.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
          }
        for(int i=0;i<this->params._d_Nnp_data();i++)
        kvector_data2d.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
          for(int i=0;i<params._N_mu_bins();i++)
         muvector.push_back(-1.0+this->params._d_Deltamu()*(i+0.5));
         this->pk0.clear();
         this->pk0.shrink_to_fit();
         this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
         this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
         this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole
          this->pk_w.clear();
          this->pk_w.shrink_to_fit();
          this->pk_w.resize(this->params._d_Nnp_window(),0); //W(k)
          this->modes_g.clear();
          this->modes_g.shrink_to_fit();
          this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
          this->pkk.resize(this->params._d_Nnp_data());
          this->pmk.resize(params._N_mu_bins());
          for(int i=0;i<this->params._d_Nnp_data();i++)this->pkk[i].resize(this->params._d_Nnp_data(),0);
          for(int i=0;i<params._N_mu_bins();i++)this->pmk[i].resize(this->params._d_Nnp_data(),0);
           this->sigma_fkp.clear();
          this->sigma_fkp.shrink_to_fit();
          this->sigma_fkp.resize(this->params._d_Nnp_data(),0);
         fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
          sigma_y_l2.resize(this->params._d_Nnp_data(),0);
          sigma_y_l4.resize(this->params._d_Nnp_data(),0);
          if(true==params._FKP_error_bars())
            {
              So.message("Computing FKP error bars");
                      fftw_functions.get_fkp_error_bars(&s_data_struct_r_ds, kvector_data, this->pk0, this->modes_g, this->sigma_fkp);
            }
     }
      // Estimates of Bispectrum. Using the DFT already done for P(k)
      else if(this->params._statistics()=="Bk_fkp")
        {
          if(this->params._type_of_binning()=="linear")
        for(int i=0;i<this->params._d_Nnp_data();i++)
          kvector_data_b.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5)); //Oficcial binning
          else
        if(this->params._type_of_binning()=="log"){
          for(int i=0;i<this->params._d_Nnp_data();i++)
            kvector_data_b.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
        }
          bispectrum.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
          sn_bispectrum.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
          modes_tri.resize(this->params._Nft()*this->params._Nft()*this->params._Nft());
          fftw_functions.get_bispectrum_fkp('d',bispectrum, sn_bispectrum, modes_tri);
          File.write_to_file(file_bispectrum,kvector_data_b,bispectrum,modes_tri);
        }
      // Estimates of Bispectrum for FKP using fast version
      else if(this->params._statistics()=="Bk_fkp_fast")
        {
           this->pk0.resize(this->params._d_Nnp_data(),0);
          for(int i=0;i<fftw_functions.Nshells_bk;i++)
        kvector_data_b.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5)); //Jennifer's binning
        bispectrum.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
           sn_bispectrum.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
              modes_tri.resize(fftw_functions.Nshells_bk*fftw_functions.Nshells_bk*fftw_functions.Nshells_bk,0);
              fftw_functions.get_power_spectrum_for_bispectrum(this->pk0);
              fftw_functions.get_bispectrum_fkp_fast(this->pk0,bispectrum,modes_tri,file_bispectrum);
        }
      So.DONE();
#ifndef _WRITE_MULTIPOLES_
      write_power_and_modes();
#else
      write_power_spectrum();
#endif
#ifdef _USE_MASS_CUTS_PK_
    }
#endif
    }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::compute_power_spectrum(vector<s_Halo>& tracer_cat, string space_p,string property){
   this->So.enter(__PRETTY_FUNCTION__);
   if(tracer_cat.size()==0)
   {
     So.message_warning("Empty halo catalog in ", __LINE__);
     exit(0);
   }
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   string file_pow;
   if(space_p=="real")
     file_pow=this->file_power_real_space;
   else if (space_p=="redshift")
     file_pow=this->file_power_redshift_space;
   int N_intervals=1;
   if(property==_MASS_)
     N_intervals=this->params._NMASSbins_power();
   else if (property==_VMAX_)
     N_intervals=this->params._NVMAXbins_power();
   So.message_screen("Selected propertyy ", property);
   So.message_screen("Number of intervals to measure power in, ", N_intervals);
   this->power_in_bins.resize(N_intervals);
   for(int im=0; im< N_intervals;++im)
     {
     this->file_power=file_pow;
     if (property==_VMAX_)
     {
       So.message_screen("Minimum value of vmax selected = ",pow(10,this->params._VMAXbins_min(im)));
       So.message_screen("Maximum value of vmax selected = ",pow(10,this->params._VMAXbins_max(im)));
       this->file_power=file_pow+"_massbin"+to_string(im);
     }
     else if(property==_MASS_)
     {
       So.message_screen("Minimum value of mass selected = ",pow(10,this->params._MASSbins_min(im)));
       So.message_screen("Maximum value of mass selected = ",pow(10,this->params._MASSbins_max(im)));
       this->file_power=file_pow+"_vmaxbin"+to_string(im);
     }
     this->File.input_type=this->params._input_type();
     ULONG in_new=0;
     vector<s_Halo> final_cat; //* Define a s_Halo container
     if(property==_MASS_)
       {
        in_new=0;
        if(space_p=="redshift")
         {
          for(ULONG i=0;i<tracer_cat.size();++i)
           if(tracer_cat[i].mass>=pow(10,this->params._MASSbins_min(im)) && tracer_cat[i].mass<pow(10,this->params._MASSbins_max(im)))
           {
           final_cat.push_back(s_Halo());
           final_cat[in_new].mass=tracer_cat[i].mass;
           final_cat[in_new].coord1=tracer_cat[i].coord1;
           final_cat[in_new].coord2=tracer_cat[i].coord2;
           final_cat[in_new].coord3=tracer_cat[i].coord3;
           final_cat[in_new].vel1=tracer_cat[i].vel1;
           final_cat[in_new].vel2=tracer_cat[i].vel2;
           final_cat[in_new].vel3=tracer_cat[i].vel3;
           in_new++;
         }
       }
       else if(space_p=="real")
         {
         for(ULONG i=0;i<tracer_cat.size();++i)
          if(tracer_cat[i].mass>=pow(10,this->params._MASSbins_min(im)) && tracer_cat[i].mass<pow(10,this->params._MASSbins_max(im)))
           {
             final_cat.push_back(s_Halo());
//             cout<<tracer_cat[i].mass<<"  "<<in_new<<endl;
             final_cat[in_new].mass=tracer_cat[i].mass;
             final_cat[in_new].coord1=tracer_cat[i].coord1;
             final_cat[in_new].coord2=tracer_cat[i].coord2;
             final_cat[in_new].coord3=tracer_cat[i].coord3;
             in_new++;
           }
         }
     }
       else if (property==_VMAX_)
     {
       in_new=0;
       if(space_p=="redshift")
         {
           for(ULONG i=0;i<tracer_cat.size();++i)
         if(log10(tracer_cat[i].vmax)>=this->params._VMAXbins_min(im) && log10(tracer_cat[i].vmax)< this->params._VMAXbins_max(im))
           {
             final_cat.push_back(s_Halo());
             final_cat[in_new].vmax=tracer_cat[i].vmax;
             final_cat[in_new].coord1=tracer_cat[i].coord1;
             final_cat[in_new].coord2=tracer_cat[i].coord2;
             final_cat[in_new].coord3=tracer_cat[i].coord3;
             final_cat[in_new].vel1=tracer_cat[i].vel1;
             final_cat[in_new].vel2=tracer_cat[i].vel2;
             final_cat[in_new].vel3=tracer_cat[i].vel3;
             in_new++;
           }
         }
       else
        {
         for(ULONG i=0;i<tracer_cat.size();++i)
          if(log10(tracer_cat[i].vmax)>= this->params._VMAXbins_min(im) && log10(tracer_cat[i].vmax)< this->params._VMAXbins_max(im))
           {
             final_cat.push_back(s_Halo());
             final_cat[in_new].vmax=tracer_cat[i].vmax;
             final_cat[in_new].coord1=tracer_cat[i].coord1;
             final_cat[in_new].coord2=tracer_cat[i].coord2;
             final_cat[in_new].coord3=tracer_cat[i].coord3;
             in_new++;
           }
         }
      }
       else if (property=="_NONE_")
      {
       in_new=tracer_cat.size();
       final_cat.resize(tracer_cat.size());
       if(space_p=="redshift")
         {
           for(ULONG i=0;i<tracer_cat.size();++i)
          {
           final_cat[i].coord1=tracer_cat[i].coord1;
           final_cat[i].coord2=tracer_cat[i].coord2;
           final_cat[i].coord3=tracer_cat[i].coord3;
           final_cat[i].vel1=tracer_cat[i].vel1;
           final_cat[i].vel2=tracer_cat[i].vel2;
           final_cat[i].vel3=tracer_cat[i].vel3;
           }
        }
       else
        {
          for(ULONG i=0;i<tracer_cat.size();++i)
           {
            final_cat[i].coord1=tracer_cat[i].coord1;
            final_cat[i].coord2=tracer_cat[i].coord2;
            final_cat[i].coord3=tracer_cat[i].coord3;
           }
        }
     }
     So.message_screen("Number of tracers selected (in power func)= ",in_new);
     Catalog new_catalog(this->params); // Make this sub-sample a member of an instnace of the class Catalog
     new_catalog.set_tracer_catalog(final_cat);//bring final_cat into new_catalog.Halo
     final_cat.clear();final_cat.shrink_to_fit();//release memmory
     this->fftw_functions.write_fftw_parameters();//get parameter for estimator
     So.message_screen("Interpolating galaxy density field on a grid");
     if(space_p=="real")
       new_catalog.get_interpolated_density_field_real_space(false,"any");
     else if (space_p=="redshift")
       new_catalog.get_interpolated_density_field(false,"any");
     So.message_screen("Number of tracers interpolated",new_catalog._n_gal());
     this->fftw_functions.resize_fftw_vectors();
     this->fftw_functions.write_fftw_parameters();
     this->fftw_functions.set_n_gal(new_catalog._n_gal());
     this->fftw_functions.set_w_g(new_catalog._w_g());
     this->fftw_functions.set_s_g(new_catalog._s_g());
     this->fftw_functions.set_sg1(new_catalog._sg1());
     this->fftw_functions.set_sg2(new_catalog._sg2());
     this->fftw_functions.get_parameters_estimator(true);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<this->params._NGRID();++i)// Here we need to pass the fields from catalog to fftw:
       this->fftw_functions.set_data_g(i,new_catalog.field_external[i]);// Set rho_r to its container in fftw_functions
       // Free memmory
     new_catalog.field_external.clear();
     new_catalog.field_external.shrink_to_fit();
     fftw_functions.raw_sampling(pow(this->params._Lbox(),3));
     fftw_functions.get_fluctuation();       //Get fluctuation
     kvector_data.clear();
     kvector_data.shrink_to_fit();
     kvector_window.clear();
     kvector_window.shrink_to_fit();
     if(this->params._type_of_binning()=="linear")
     {
       for(int i=0;i<this->params._d_Nnp_data();i++)
         kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
       for(int i=0;i<this->params._d_Nnp_window();i++)
         kvector_window.push_back(this->params._d_kmin()+this->params._d_DeltaK_window()*(i+0.5));
     }
       else  if(this->params._type_of_binning()=="log")
     {
       for(int i=0;i<kvector_data.size();i++)
         kvector_data.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
       for(int i=0;i<kvector_window.size();i++)
         kvector_window.push_back(this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal()));
     }
       for(int i=0;i<this->params._d_Nnp_data();i++)
         kvector_data2d.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
       for(int i=0;i<params._N_mu_bins();i++)
         muvector.push_back(-1.0+this->params._d_Deltamu()*(i+0.5));
       this->pk0.clear();
       this->pk0.shrink_to_fit();
       this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
       this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
       this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole
       this->pk_w.clear();
       this->pk_w.shrink_to_fit();
       this->pk_w.resize(this->params._d_Nnp_window(),0); //W(k)
       this->modes_g.clear();
       this->modes_g.shrink_to_fit();
       this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
       this->pkk.resize(this->params._d_Nnp_data());
       this->pmk.resize(params._N_mu_bins());
       fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
       this->power_in_bins[im].i_ncbin=this->modes_g;
       this->power_in_bins[im].vbin=this->kvector_data;
       this->power_in_bins[im].vq1=this->pk0;
       this->power_in_bins[im].vq2=this->pk2;
       this->power_in_bins[im].vq3=this->pk4;
       this->write_power_and_modes();
       this->write_fftw_parameters();
     }// closes loop of bins or cuts
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define _use_two_quartiles_
 void PowerSpectrumF::halo_bias_analysis(string space_p)
 {
     time_t time_POWER;
     time(&time_POWER);
     bool get_marked=this->params._Get_marked_power_spectrum();
#ifdef _VERBOSE_POWER_
   this->So.enter(__PRETTY_FUNCTION__);
#endif
   time_t start;
   time (&start);
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   // **************************************************************************************************
   Catalog tracer(params);
   tracer.read_catalog(params._Input_dir_cat()+params._file_catalogue(),0);
   this->set_params(tracer.params); // params has been updated inside the catalog class.
   tracer.Get_SO_tracer();
   // **************************************************************************************************
   this->fftw_functions.set_params(this->params);
   // **************************************************************************************************
   ULONG Nbins=this->params._Nbins_hist();  // Number of bins for 1d or 2d histograms
   // **************************************************************************************************
   string file_pow;   //Output file name
   if (space_p=="redshift_space")
     file_pow=this->file_power_redshift_space;
   else
     file_pow=this->file_power_real_space;

   string file_pow_marked;
   if (space_p=="redshift_space")
     file_pow_marked=this->file_power_marked_redshift_space;
   else
     file_pow_marked=this->file_power_marked_real_space;
   string file_pow_cross=this->file_power_cross;
   // **************************************************************************************************
   vector<string>prop_name;
   vector<int>prop_Nbins;
   vector<real_prec>variance_property;
   vector<real_prec>mean_prop;
   vector<bool>used_prop;
   vector<bool>secondary_prop;
   // *Get the total number of propertye bis* //
   // * The ordeing of alocation must coincide with teh ordering of the loops */
  // **************************************************************************************************
    // Counter on properties
   int Number_of_properties=0;
   // **************************************************************************************************
   int label_mass=Number_of_properties;
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
   if(this->params._i_mass_g()<0)
   {
     this->So.message_screen("Warning, Mass used as primary property but not specified in input par file. Check");
    exit(0);
   }
   int Nbins_mass=this->params._NMASSbins_power();
   prop_name.push_back("_MASS_");
   used_prop.push_back(true);
   prop_Nbins.push_back(Nbins_mass);
   secondary_prop.push_back(false);
   Number_of_properties++;
#else
   bool bmass=false;
   if(this->params._i_mass_g()>0){
     used_prop.push_back(true);
     bmass=true;
   }
   else
     used_prop.push_back(false);
   int Nbins_mass=this->params._NVMAXbins_power();
   prop_Nbins.push_back(Nbins_mass);
   prop_name.push_back("_MASS_");
   Number_of_properties++;
   secondary_prop.push_back(true);
#endif
   // ---------------------------------------
   int label_vmax=Number_of_properties;
#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
   if(this->params._i_vmax_g()<0)
   {
     this->So.message_screen("Warning, Vmax used as primary property but not specified in input par file. Check");
    exit(0);
   }
   int Nbins_vmax=this->params._NVMAXbins_power();
   prop_name.push_back("_VMAX_");
   used_prop.push_back(true);
   prop_Nbins.push_back(Nbins_vmax);
   secondary_prop.push_back(false);
   Number_of_properties++;
#else
   bool bvmax=false;
   if(this->params._i_vmax_g()>0)
   {
     used_prop.push_back(true);
    bvmax=true;
   }
   else
     used_prop.push_back(false);
   int Nbins_vmax=this->params._NVMAXbins_power();
   prop_Nbins.push_back(Nbins_vmax);
   prop_name.push_back("_VMAX_");
   secondary_prop.push_back(true);

   Number_of_properties++;
#endif

   // ---------------------------------------
   int label_rs=Number_of_properties;
   if(this->params._i_rs_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_rs=this->params._NRSbins_power();
   prop_Nbins.push_back(Nbins_rs);
   prop_name.push_back("_RS_");
   Number_of_properties++;
   // ---------------------------------------
   int label_rvir=Number_of_properties;
   if(this->params._i_rvir_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_rvir=this->params._NRVIRbins_power();
   prop_Nbins.push_back(Nbins_rvir);
   prop_name.push_back("_RVIR_");
   Number_of_properties++;
   // ---------------------------------------
   int label_concentration=Number_of_properties;
#ifdef _USE_CONCENTRATION_
   if(this->params._i_rvir_g()<0 || this->params._i_rs_g()<0)
    {
       So.message_warning("Rvir or Rs not used. Concentration cannot be defined");
       used_prop.push_back(false);
    }
    else
       used_prop.push_back(true);
#endif
   secondary_prop.push_back(true);
   int Nbins_concentration=this->params._NCONCENTRATIONbins_power();
   prop_Nbins.push_back(Nbins_concentration);
   prop_name.push_back("_CONCENTRATION_");
   Number_of_properties++;
   // ---------------------------------------
   int label_spin=Number_of_properties;
   if(this->params._i_spin_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_spin=this->params._NSPINbins_power();
   prop_Nbins.push_back(Nbins_spin);
   prop_name.push_back("_SPIN_");
   Number_of_properties++;
   // ---------------------------------------
   int label_spin_bullock=Number_of_properties;
   if(this->params._i_spin_bullock_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_spin_bullock=this->params._NSPINBULLOCKbins_power();
   prop_Nbins.push_back(Nbins_spin_bullock);
   prop_name.push_back("_SPIN_BULLOCK_");
   Number_of_properties++;
   // ---------------------------------------
   int label_vrms=Number_of_properties;
   if(this->params._i_vrms_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_vrms=this->params._NVRMSbins_power();
   prop_Nbins.push_back(Nbins_vrms);
   prop_name.push_back("_VRMS_");
   Number_of_properties++;
   // ---------------------------------------
   int label_virial=Number_of_properties;
   if(this->params._i_virial_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_virial=this->params._NVIRIALbins_power();
   prop_Nbins.push_back(Nbins_virial);
   prop_name.push_back("_VIRIAL_");
   Number_of_properties++;
   // ---------------------------------------
   int label_btoa=Number_of_properties;
   if(this->params._i_b_to_a_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_btoa=this->params._NBTOAbins_power();
   prop_Nbins.push_back(Nbins_btoa);
   prop_name.push_back("_BTOA_");
   Number_of_properties++;
   // ---------------------------------------
   int label_ctoa=Number_of_properties;
   if(this->params._i_c_to_a_g()>0)
     used_prop.push_back(true);
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_ctoa=this->params._NCTOAbins_power();
   prop_Nbins.push_back(Nbins_ctoa);
   prop_name.push_back("_CTOA_");
   Number_of_properties++;
   // ---------------------------------------
   int label_mach=Number_of_properties;
   int Nbins_mach=this->params._NMACHbins_power(); // ya viene seleccionad a 1 enc aso que no queramos mach
   prop_Nbins.push_back(Nbins_mach);
   bool mach=false;
   prop_name.push_back("_MACH_");
   secondary_prop.push_back(true);
   Number_of_properties++;
   if(this->params._Get_tracer_local_mach_number())
     {
       used_prop.push_back(true);
       mach=true;
     }
   else
   used_prop.push_back(false);
   // ---------------------------------------
   int label_bias=Number_of_properties;
   bool tbias=false;
   int Nbins_bias=this->params._NBIASbins_power();// ya viene seleccionado a 1 en caso que no queramos bias
   prop_Nbins.push_back(Nbins_bias);
   prop_name.push_back("_BIAS_");
   Number_of_properties++;
   if(this->params._Get_tracer_bias()){
     used_prop.push_back(true);
     tbias=true;
   }
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   // ---------------------------------------
   int label_lc=Number_of_properties;
   bool tlc=false;
   int Nbins_lc=this->params._NLCbins_power();// ya viene seleccionado a 1 en caso que no queramos bias
   prop_Nbins.push_back(Nbins_lc);
   prop_name.push_back("_LOCAL_OVERDENSITY_");
   Number_of_properties++;
   if(this->params._Get_local_overdensity()){
     tlc=true;
     used_prop.push_back(true);
   }
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   // ---------------------------------------
   int label_ta=Number_of_properties;
   prop_name.push_back("_TIDAL_ANISOTROPY_");
   Number_of_properties++;
   int Nbins_ta=this->params._NTAbins_power();// ya viene seleccionado a 1 en caso que no queramos bias
   prop_Nbins.push_back(Nbins_ta);
   bool bta=false;
   if(this->params._Get_tidal_anisotropy_at_halo()) {
     bta=true;
     used_prop.push_back(true);
   }
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   // ---------------------------------------
   int label_ph=Number_of_properties;
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
   int Nbins_ph=this->params._NPHbins_power();
   prop_name.push_back("_PEAK_HEIGHT_");
   used_prop.push_back(true);
   prop_Nbins.push_back(Nbins_ph);
   secondary_prop.push_back(false);
   Number_of_properties++;
#else
   prop_name.push_back("_PEAK_HEIGHT_");
   Number_of_properties++;
   bool bph=false;
   if(this->params._Get_peak_height_at_halo()>0){
     used_prop.push_back(true);
    bph=true;
   }
   else
     used_prop.push_back(false);
   secondary_prop.push_back(true);
   int Nbins_ph=this->params._NPHbins_power();
   prop_Nbins.push_back(Nbins_ph);
#endif
  // ---------------------------------------
   int label_da=Number_of_properties;
   prop_name.push_back("_DACH_");
   Number_of_properties++;
   int Nbins_dach=this->params._NDACHbins_power();// ya viene seleccionado a 1 en caso que no queramos bias
   prop_Nbins.push_back(Nbins_dach);
   bool dach=false;
   secondary_prop.push_back(true);
   if(true==this->params._Get_tracer_local_dach_number())
   {
     dach=true;
     used_prop.push_back(true);
   }
   else
     used_prop.push_back(false);
   // ---------------------------------------
   int label_dm_local=Number_of_properties;
   prop_name.push_back("_LOCALDM_");
   Number_of_properties++;
   int Nbins_dm_local=this->params._NLOCALDMbins_power();// ya viene seleccionado a 1 en caso que no queramos bias
   prop_Nbins.push_back(Nbins_dm_local);
   bool dm_local=false;
   secondary_prop.push_back(true);
   if(true==this->params._Get_tracer_local_dm_density())
   {
     dm_local=true;
     used_prop.push_back(true);
   }
   else
     used_prop.push_back(false);

    So.message_screen("Total number of properties:", Number_of_properties);
   // **************************************************************************************************
   int Ncwt=0;
   vector<real_prec> DM_DEN_FIELD( this->params._NGRID(),0);
   this->File.read_array(this->params._Input_Directory_X()+this->params._Name_Catalog_X(),DM_DEN_FIELD);
   for(ULONG i=0;i<DM_DEN_FIELD.size();++i)
    if(DM_DEN_FIELD[i]<0)
        DM_DEN_FIELD[i]=0;
    So.message_screen("\tMin DM field:",get_min_nm(DM_DEN_FIELD));
    So.message_screen("\tMax DM field:",get_max_nm(DM_DEN_FIELD));
    So.message_screen("\tMean DM field:",get_mean(DM_DEN_FIELD));
#ifdef _USE_CWC_HALO_ANALYSIS_
   vector<real_prec>DM_OVERDENSITY_FIELD(DM_DEN_FIELD.size(),0 );
   get_overdens(DM_DEN_FIELD,DM_OVERDENSITY_FIELD);
#endif
    // **************************************************************************************************
   // Get DM power spectrum
   this->pk0.clear();
   this->pk0.shrink_to_fit();
   this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
   this->params.set_input_type("density_grid");
   string original_name=this->params._Name_survey();
   this->params.set_Name_survey("UNITSIM_DM");
   this->params.set_mass_assignment_scheme("CIC");
   this->params.set_SN_correction(false);
   this->set_output_filenames();
   this->compute_power_spectrum_grid(DM_DEN_FIELD, true);
   vector<real_prec>power_dm(this->params._d_Nnp_data(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<power_dm.size();++i)
     power_dm[i]=this->pk0[i];      // The container pk0 will be cleared and reused. power_dm iwll be used to get bias object by object
   // -------------------------------------------------------------------------------------------------------
   if(true==this->params._Get_tracer_bias())
     {
       this->So.message_screen("Getting object-to-object bias");
       if(true==this->params._Get_tracer_relative_bias())
       {
         vector<real_prec> TR_DEN_FIELD( this->params._NGRID(),0);
         this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),TR_DEN_FIELD);
         this->object_by_object_bias(tracer.Halo,DM_DEN_FIELD,TR_DEN_FIELD);
       }
       else
        this->object_by_object_bias(tracer.Halo,DM_DEN_FIELD);
    this->So.message_screen("\tMin Bias  =", tracer.get_min("_BIAS_"));
    this->So.message_screen("\tMax Bias  =", tracer.get_max("_BIAS_"));
    this->So.DONE();
    cout<<""<<endl;
   if(true==this->params._Get_tracer_quadratic_bias())
       {
         this->So.message_screen("Getting object-to-object quadratic bias");
         this->object_by_object_qbias(tracer.Halo,DM_DEN_FIELD);
       }
    }
   // **************************************************************************************************
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
   // If this preproc is disableD, we can still want to use peak height as secondary property so we follow whatever the parameter file said
   tracer.get_peak_height_at_tracer();
#else
   if(true==this->params._Get_peak_height_at_halo())
     {
       tracer.get_peak_height_at_tracer();
       this->So.message_screen("   Min Peak Height  =", tracer.get_min("_PEAK_HEIGHT_"));
       this->So.message_screen("   Max Peak Height =", tracer.get_max("_PEAK_HEIGHT_"));
       cout<<""<<endl;
     }
#endif

   // **************************************************************************************************
   // Get CW classification
#ifdef _USE_CWC_HALO_ANALYSIS_
   // We can here read the dark matter density field and determine the CWC
   this->cwclass.set_params(this->params);
   this->cwclass.get_CWC(DM_OVERDENSITY_FIELD);
   string file_cwc = this->params._Output_directory()+"CWC_vf_lambdath"+to_string(this->params._lambdath())+"_p"+to_string(this->params._unitsim_plabel())+"_Nft"+to_string(this->params._Nft())+".txt";
   ofstream cwcf; cwcf.open(file_cwc.c_str());
   cwcf<<this->params._lambdath()<<"\t"<<this->cwclass._knots_fraction()<<"\t"<<this->cwclass._filaments_fraction()<<"\t"<<this->cwclass._sheets_fraction()<<"\t"<<this->cwclass._voids_fraction()<<endl;
   cwcf.close();

    cout<<endl;
   this->So.message_screen("Assigning CWT to tracers=");
#pragma omp parallel for
   for(ULONG i=0;i<tracer._NOBJS();++i)
      tracer.Halo[i].gal_cwt=this->cwclass.CWClass[tracer.Halo[i].GridID]; //esto asigna 1, 2, 3, 4

   So.DONE();

   Ncwt=this->cwclass.cwt_used.size();
   So.message_screen("Number cosmic-web types requested: ",Ncwt);
   if(true==this->params._Get_tidal_anisotropy_at_halo())
     {
       vector<real_prec>tidal(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<this->params._NGRID() ;++i)
         tidal[i]=tidal_anisotropy(this->cwclass.lambda1[i], this->cwclass.lambda2[i], this->cwclass.lambda3[i]);
       tracer.get_tracer_tidal_anisotropy(tidal);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<tracer._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
       tracer.Halo[i].lambda1 =this->cwclass.lambda1[tracer.Halo[i].GridID];
       tracer.Halo[i].lambda2 =this->cwclass.lambda2[tracer.Halo[i].GridID];
       tracer.Halo[i].lambda3 =this->cwclass.lambda3[tracer.Halo[i].GridID];
    }
     this->So.message_screen("    Min tidal_anisotropy  =", tracer.get_min("_TIDAL_ANISOTROPY_"));
     this->So.message_screen("    Max tidal_anisotropy  =", tracer.get_max("_TIDAL_ANISOTROPY_"));
     cout<<""<<endl;
     tidal.clear();tidal.shrink_to_fit();
     }
#endif

    // **************************************************************************************************
    if(true==this->params._Get_tracer_local_dm_density())
    {
    this->So.message_screen("Assigning log(delta_dm+1) from cells to halos");
 #pragma omp parallel for
    for(ULONG i=0;i<tracer._NOBJS();++i)
       tracer.Halo[i].local_dm=log10(1+DM_OVERDENSITY_FIELD[tracer.Halo[i].GridID]);
    this->So.message_screen("    Min localdm  =", tracer.get_min("_LOCALDM_"));
    this->So.message_screen("    Max localdm  =", tracer.get_max("_LOCALDM_"));
    }

#ifdef _USE_CWC_HALO_ANALYSIS_
    DM_OVERDENSITY_FIELD.clear();
    DM_OVERDENSITY_FIELD.shrink_to_fit();
#endif
   this->params.set_Name_survey(original_name);
   this->params.set_SN_correction(true);
   this->params.set_mass_assignment_scheme("TSC");
   this->set_output_filenames();
   // **************************************************************************************************
   // **************************************************************************************************
   // **************************************************************************************************
   // Define number of bins in each property
   int Nmass_ind=Nbins_mass;
   int Nvmax_ind=Nbins_vmax;
   int Nrs_ind=Nbins_rs;
   int Nrvir_ind=Nbins_rvir;
   int Nrconcentration_ind=Nbins_concentration;
   int Nspin_ind=Nbins_spin;
   int Nspinbullock_ind=Nbins_spin_bullock;
   int Nvrms_ind=Nbins_vrms;
   int Nvirial_ind=Nbins_virial;
   int Nbtoa_ind=Nbins_btoa;
   int Nctoa_ind=Nbins_ctoa;
   int Nmach_ind=Nbins_mach;
   int Nbias_ind=Nbins_bias;
   int Nlc_ind=Nbins_lc;
   int Nta_ind=Nbins_ta;
   int Nph_ind=Nbins_ph;
   int Ndach_ind=Nbins_dach;
   int Ndmlocal_ind=Nbins_dm_local;
   // **************************************************************************************************
   // **************************************************************************************************
   // **************************************************************************************************
   int Nprimary_bins=0;
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
    Nprimary_bins=Nbins_mass;
#endif

#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
    Nprimary_bins=Nbins_vmax;
#endif
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
    Nprimary_bins=Nbins_ph;
#endif

    string rcatalog=this->params._Output_directory()+"catalog_reduced.txt";
    tracer.select_random_subsample(this->params._fraction_dilute(), rcatalog); // DONE IN PCA NOW

    // **************************************************************************************************
    // WE can write a random sub-sample of the total catalogs to be read in python and do 2d historgrams.
//    string rcatalog=this->params._Output_directory()+"catalog_reduced.txt";
//    tracer.select_random_subsample(0.2, rcatalog); // DONE IN PCA NOW
/*
    // get randomized versions
    Catalog tracer_randomize(params);
    tracer_randomize.define_property_bins(this->params._NPROPbins_bam(), "_MASS_");
    string rcatalog;
    tracer_randomize.set_NOBJS(tracer._NOBJS());

    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_MACH_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_", "_MACH_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_MACH_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_mach_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);


    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_DACH_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_", "_DACH_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_DACH_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_dach_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);


    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_VIRIAL_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_", "_VIRIAL_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_VIRIAL_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_virial.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_CTOA_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_", "_CTOA_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_CTOA_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_ctoa.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCALDM_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_", "_LOCALDM_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCALDM_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_localdm_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_TIDAL_ANISOTROPY_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_TIDAL_ANISOTROPY_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_TIDAL_ANISOTROPY_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_ta_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCAL_OVERDENSITY_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_LOCAL_OVERDENSITY_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCAL_OVERDENSITY_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_twopr_lc_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

   tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCAL_OVERDENSITY_", "_MACH_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_LOCAL_OVERDENSITY_", "_MACH_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCAL_OVERDENSITY_", "_MACH_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_mach_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCAL_OVERDENSITY_", "_LOCALDM_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_LOCAL_OVERDENSITY_", "_LOCALDM_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCAL_OVERDENSITY_", "_LOCALDM_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_dm_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCAL_OVERDENSITY_", "_CTOA_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_LOCAL_OVERDENSITY_", "_CTOA_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCAL_OVERDENSITY_", "_CTOA_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_ctoa_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);


    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_VMAX_", "_LOCAL_OVERDENSITY_", "_VIRIAL_");
    tracer_randomize.get_scaling_relation_primary_property("_CONCENTRATION_","_LOCAL_OVERDENSITY_", "_VIRIAL_");
    tracer_randomize.get_scaling_relation_primary_property("_SPIN_BULLOCK_", "_LOCAL_OVERDENSITY_", "_VIRIAL_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_virial_cwt.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);


    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_BIAS_", "_LOCAL_OVERDENSITY_", "_LOCALDM_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_dm_cwt_newbias.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_BIAS_", "_LOCAL_OVERDENSITY_", "_MACH_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_mach_cwt_newbias.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_BIAS_", "_LOCAL_OVERDENSITY_", "_TIDAL_ANISOTROPY_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_lc_ta_cwt_newbias.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_primary_property("_BIAS_", "_LOCALDM_", "_MACH_");
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_threepr_dm_mach_cwt_newbias.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);

    tracer_randomize.Halo.clear();
    tracer_randomize.Halo.shrink_to_fit() ;
    tracer_randomize.Halo.resize(tracer._NOBJS());
    tracer_randomize.Halo=tracer.Halo;
    tracer_randomize.get_scaling_relation_bias();
    rcatalog=this->params._Output_directory()+"catalog_reduced_randomized_newbias_newprops.txt";
    tracer_randomize.select_random_subsample(0.3,rcatalog);
*/
    // **************************************************************************************************
    // PRINCIPAL COMPONENT ANALYSIS. Reduced catalogs are written in the PCA method
    if(this->params._Get_PCA())
    {
        So.message_screen("Doing PCA");
         tracer.PCA(prop_name, used_prop, "", true);
         So.DONE();
         exit(0);
    }
    // **************************************************************************************************
    // **************************************************************************************************
    // SECTION TO DO 2D HISTOGRAMS-CORRELATION RELATIONS BETWEEN HALOS
/*      {
        vector<real_prec>Vaux(Nbins*Nbins,0);
        get_2d_histogram(this->params._ldelta_X_min(),this->params._ldelta_X_max(),this->params._ldelta_Y_min(),this->params._ldelta_Y_max(),Nbins, DM_DEN_FIELD, aux_field,Vaux,true);
        this->File.write_array(file_hist,Vaux);
        this->mcmc.get_contour_levels(file_contours+"_contour_levels",Nbins, Vaux);
        Vaux.clear();Vaux.shrink_to_fit();
      }*/
    // **************************************************************************************************
#ifdef _USE_CWC_HALO_ANALYSIS_
   for(int icw=0;icw<=Ncwt;icw++)  // loop over nuber of cosmic-web types. ALL=0 knots. 1 = filaments: 2=sheets; 3 =voids
     {
       int ict=icw;// Ojo, esto asume que Halo[].gal_cwt es 1, 2, 3, 4 para k, f, s, v respectivamente
       So.message_screen("CWT ",icw);
#endif
       // ************************************LOOP OVER THE BINS OF THE PRIMARY PROPERTY**********************************************+
       for(int im_primary=0; im_primary< Nprimary_bins;++im_primary)
     {
        vector<real_prec> primary_field(this->params._NGRID(), 0);
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
       So.message_screen("Mass bin ",im_primary);
#elif defined _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
       So.message_screen("Vmax bin ",im_primary);
#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
       So.message_screen("Nu-bin ",im_primary);
#endif
       Catalog tracer_aux(this->params);
       // if the bins are to be defined with same number of objects, here is the place to redefine those bins for each mass bin
       // using inside this loop the  members f the params class setbins*()
       // *Selecting the trabcers and its propertis in each mass bin*//
       ULONG counter_m=0;
       So.message_screen("Selecting tracers in bin of primary property");
#ifdef _USE_CWC_HALO_ANALYSIS_
       if(icw==0)// THis selects all objects regardless their CW clasiffication.
         {
#endif
           for(ULONG i=0;i<tracer._NOBJS();++i)
         {
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
         if(tracer.Halo[i].peak_height>= this->params._PHbins_min(im_primary) && tracer.Halo[i].peak_height< this->params._PHbins_max(im_primary))
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
           if(tracer.Halo[i].mass>= pow(10,this->params._MASSbins_min(im_primary)) && tracer.Halo[i].mass< pow(10,this->params._MASSbins_max(im_primary)))
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
             if(tracer.Halo[i].vmax>= pow(10,this->params._VMAXbins_min(im_primary)) && tracer.Halo[i].vmax< pow(10,this->params._VMAXbins_max(im_primary)))
#endif
               {
             tracer_aux.Halo.push_back(s_Halo());
             tracer_aux.Halo[counter_m].GridID = tracer.Halo[i].GridID;
             tracer_aux.Halo[counter_m].coord1 = tracer.Halo[i].coord1;
             tracer_aux.Halo[counter_m].coord2 = tracer.Halo[i].coord2;
             tracer_aux.Halo[counter_m].coord3 = tracer.Halo[i].coord3;
             tracer_aux.Halo[counter_m].vel1 = tracer.Halo[i].vel1;
             tracer_aux.Halo[counter_m].vel2 = tracer.Halo[i].vel2;
             tracer_aux.Halo[counter_m].vel3 = tracer.Halo[i].vel3;
             tracer_aux.Halo[counter_m].mass = tracer.Halo[i].mass;
             tracer_aux.Halo[counter_m].vmax = tracer.Halo[i].vmax;
             tracer_aux.Halo[counter_m].rs = tracer.Halo[i].rs;
             tracer_aux.Halo[counter_m].rvir = tracer.Halo[i].rvir;
             tracer_aux.Halo[counter_m].concentration = tracer.Halo[i].concentration;
             tracer_aux.Halo[counter_m].spin = tracer.Halo[i].spin;
             tracer_aux.Halo[counter_m].spin_bullock = tracer.Halo[i].spin_bullock;
             tracer_aux.Halo[counter_m].vrms = tracer.Halo[i].vrms;
             tracer_aux.Halo[counter_m].virial = tracer.Halo[i].virial;
             tracer_aux.Halo[counter_m].b_to_a = tracer.Halo[i].b_to_a;
             tracer_aux.Halo[counter_m].c_to_a = tracer.Halo[i].c_to_a;
             tracer_aux.Halo[counter_m].mach_number = tracer.Halo[i].mach_number;
             tracer_aux.Halo[counter_m].local_overdensity = tracer.Halo[i].local_overdensity;
             tracer_aux.Halo[counter_m].bias = tracer.Halo[i].bias;
             tracer_aux.Halo[counter_m].tidal_anisotropy = tracer.Halo[i].tidal_anisotropy;
             tracer_aux.Halo[counter_m].peak_height = tracer.Halo[i].peak_height;
             tracer_aux.Halo[counter_m].dach_number= tracer.Halo[i].dach_number;
             tracer_aux.Halo[counter_m].local_dm= tracer.Halo[i].local_dm;
             counter_m++;
               }//closes if
         } // closes loop
         }//closes if cwc
#ifdef _USE_CWC_HALO_ANALYSIS_
       else if(icw>0)
         {
           for(ULONG i=0;i<tracer._NOBJS();++i)
         {
           if(tracer.Halo[i].gal_cwt==ict)
             {
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
               if(tracer.Halo[i].peak_height>= this->params._PHbins_min(im_primary) && tracer.Halo[i].peak_height< this->params._PHbins_max(im_primary))
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
             if(tracer.Halo[i].mass>= pow(10,this->params._MASSbins_min(im_primary)) && tracer.Halo[i].mass< pow(10,this->params._MASSbins_max(im_primary)))
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
               if(tracer.Halo[i].vmax>= pow(10,this->params._VMAXbins_min(im_primary)) && tracer.Halo[i].vmax< pow(10,this->params._VMAXbins_max(im_primary)))
#endif
                 {
                   tracer_aux.Halo.push_back(s_Halo());
                   tracer_aux.Halo[counter_m].GridID = tracer.Halo[i].GridID;
                   tracer_aux.Halo[counter_m].coord1 = tracer.Halo[i].coord1;
                   tracer_aux.Halo[counter_m].coord2 = tracer.Halo[i].coord2;
                   tracer_aux.Halo[counter_m].coord3 = tracer.Halo[i].coord3;
                   tracer_aux.Halo[counter_m].vel1 = tracer.Halo[i].vel1;
                   tracer_aux.Halo[counter_m].vel2 = tracer.Halo[i].vel2;
                   tracer_aux.Halo[counter_m].vel3 = tracer.Halo[i].vel3;
                   tracer_aux.Halo[counter_m].mass = tracer.Halo[i].mass;
                   tracer_aux.Halo[counter_m].vmax = tracer.Halo[i].vmax;
                   tracer_aux.Halo[counter_m].rs = tracer.Halo[i].rs;
                   tracer_aux.Halo[counter_m].rvir = tracer.Halo[i].rvir;
                   tracer_aux.Halo[counter_m].concentration = tracer.Halo[i].concentration;
                   tracer_aux.Halo[counter_m].spin = tracer.Halo[i].spin;
                   tracer_aux.Halo[counter_m].spin_bullock = tracer.Halo[i].spin_bullock;
                   tracer_aux.Halo[counter_m].vrms = tracer.Halo[i].vrms;
                   tracer_aux.Halo[counter_m].virial = tracer.Halo[i].virial;
                   tracer_aux.Halo[counter_m].b_to_a = tracer.Halo[i].b_to_a;
                   tracer_aux.Halo[counter_m].c_to_a = tracer.Halo[i].c_to_a;
                   tracer_aux.Halo[counter_m].mach_number = tracer.Halo[i].mach_number;
                   tracer_aux.Halo[counter_m].local_overdensity = tracer.Halo[i].local_overdensity;
                   tracer_aux.Halo[counter_m].bias = tracer.Halo[i].bias;
                   tracer_aux.Halo[counter_m].tidal_anisotropy = tracer.Halo[i].tidal_anisotropy;
                   tracer_aux.Halo[counter_m].peak_height = tracer.Halo[i].peak_height;
                   tracer_aux.Halo[counter_m].dach_number= tracer.Halo[i].dach_number;
                   tracer_aux.Halo[counter_m].local_dm= tracer.Halo[i].local_dm;
                   counter_m++;
                 } // closes if
             } // closes if
         } // closes loop
         }// closes else if
#endif
       So.DONE();
       if(counter_m==0)
         So.message_warning("No tracers found in the interval", im_primary);
       else
         So.message_screen("Found ",counter_m);
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
       tracer.set_Number_of_tracers_in_ph_bins(im_primary,counter_m);
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
       tracer.set_Number_of_tracers_in_mass_bins(im_primary,counter_m);
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
       tracer.set_Number_of_tracers_in_vmax_bins(im_primary,counter_m);
#endif
       tracer_aux.set_NOBJS(counter_m);
       // ************************************************************************************************************
       // Here we masure marked power before starting the bins ins ec properties
       cout<<endl;
      if(true==this->params._Get_marked_power_spectrum() && counter_m>0)
      {
      So.message_screen("Measuring marked power spectrum:");
       s_data_structure s_data_struct_g_marked;
       this->N_galaxy=tracer_aux.Halo.size();
       real_prec mean_density=static_cast<real_prec>(tracer_aux.Halo.size())/pow(this->params._Lbox(),3);
       s_data_struct_g_marked.mean_density=mean_density;
       fftw_functions.set_n_gal(tracer_aux.Halo.size());
       int first_sec_property_used=0;// this variable us used to computu/allocate/free memmory
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
       first_sec_property_used=1;//the first secondary property used the vmax, labled 1
#elif defined _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
       first_sec_property_used=0;//the first secondary property used the mass, labled 0
#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
           first_sec_property_used=0;//the first secondary property used the mass, labled 0
#endif
        for(int Nbb=0;Nbb<Number_of_properties;++Nbb) // Loop over the  properties to be used as marks, in bins of the primary poroperty
         {
         if(true==used_prop[Nbb] && true==secondary_prop[Nbb])//THis If is importante. Continues only if these properteis are asked to be used and are secondary
         {
           pair<real_prec,real_prec>stats=tracer_aux.get_variance(prop_name[Nbb],false);// false means that it gives <prop²> istead of <(prop-mean)²>
           real_prec mean_prop=stats.first;
           real_prec variance_property=stats.second;
           So.message_screen("Property ", prop_name[Nbb]);
           So.message_screen("Mean     ", mean_prop);
           So.message_screen("Variance ", variance_property);
           fftw_functions.set_mean_property_value(mean_prop);// Useful to divide by the mean each prop, to get prop/<mean> and then do the interpolation
#ifdef _USE_CWC_HALO_ANALYSIS_
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._Output_directory()+"lss_bias_marked_nubin"+to_string(im_primary)+"_cwc"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._Output_directory()+"lss_bias_marked_massbin"+to_string(im_primary)+"_cwc"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._Output_directory()+"lss_bias_marked_vmaxbin"+to_string(im_primary)+"_cwc"+to_string(icw)+prop_name[Nbb];
#endif
#else
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._Output_directory()+"lss_bias_marked_nubin"+to_string(im_primary)+"_"+prop_name[Nbb];
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._Output_directory()+"lss_bias_marked_massbin"+to_string(im_primary)+"_"+prop_name[Nbb];
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
           string file_lss_bias_marked=this->params._/home/balaguera/data/Numerics/UNITSIM/REFERENCE/001/HALOS/SNAPSHOT_p65/NEW/BINNING_massbin2_cwt0_p65.txtOutput_directory()+"lss_bias_marked_vmaxbin"+to_string(im_primary)+"_"+prop_name[Nbb];
#endif
#endif
           ofstream blsm; blsm.open(file_lss_bias_marked.c_str());
           this->File.input_type=this->params._input_type();
           if(Nbb==first_sec_property_used)// Compute the raw number counts, only needed once, so we do it for Nbb=first index used
             {
               this->tracer_cat.field_external.clear();
               this->tracer_cat.field_external_s.clear();
             if(false==this->params._use_real_and_redshift_space())
             {
               if (space_p=="redshift_space")
                 {
                   So.message_screen("Interpolating galaxy density field on a grid");
#ifdef _USE_SEVERAL_RANDOM_FILES_
                   this->tracer_cat.get_interpolated_density_field(false, "any",0); // this gives  fftw_functions.field_external
#else
                   this->tracer_cat.get_interpolated_density_field(false, "any"); // this gives  fftw_functions.field_external
#endif
               }
               else if (space_p=="real_space")
                 {
                   So.message_screen("Interpolating galaxy density field on a grid");
                   this->tracer_cat.get_interpolated_density_field_real_space(false,"any"); // this gives  fftw_functions.field_external
                 }
             }
               else if(true==this->params._use_real_and_redshift_space())
                 this->tracer_cat.get_interpolated_density_field_real_and_redshift_space(false, "any");
           }
           if(false==this->params._use_real_and_redshift_space())
            {
               if (space_p=="redshift_space")
             {
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked_redshift_space+"_nubin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked_redshift_space+"_massbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked_redshift_space+"_vmaxbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#endif
               So.message_screen("Interpolating weighted galaxy density field on a grid for", prop_name[Nbb]);
               this->tracer_cat.field_external_marked.clear();
#ifdef _USE_SEVERAL_RANDOM_FILES_
               this->tracer_cat.get_interpolated_density_field(true, prop_name[Nbb],0);// this gives  fftw_functions.field_external_marked
#else
               this->tracer_cat.get_interpolated_density_field(true, prop_name[Nbb]);// this gives  fftw_functions.field_external_marked
#endif
               }//closes if
            else if(space_p=="real_space")
             {
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_nubin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_massbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_vmaxbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#endif
               So.message_screen("Interpolating weighted galaxy density field on a grid for", prop_name[Nbb]);
               this->tracer_cat.field_external_marked.clear();
               this->tracer_cat.get_interpolated_density_field_real_space(true, prop_name[Nbb]); // this gives  fftw_functions.field_external_marked
             }
           }
           else if(true==this->params._use_real_and_redshift_space())
             {
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_nubin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_massbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
               file_pow_marked=this->file_power_marked+"_vmaxbin"+to_string(im_primary)+"_cwt"+to_string(icw)+prop_name[Nbb];
#endif
               So.message_screen("Interpolating weighted galaxy density field on a grid for", prop_name[Nbb]);
               this->tracer_cat.field_external_marked.clear();
               this->tracer_cat.field_external_marked_s.clear();
               this->tracer_cat.get_interpolated_density_field_real_and_redshift_space(true, prop_name[Nbb]); // this gives  fftw_functions.field_external_marked

               // Write the weighted fields to binary:
               if(im_primary==0 && icw==0)
               {
                string outf=this->params._Output_directory()+"property_field"+prop_name[Nbb];
                this->File.write_array(outf,this->tracer_cat.field_external_marked);
               }
              So.DONE();
           }
            // NOTE:
           // We avoid get_overdens as (delta_w - delta) = n_w/barn - n/ barn. Below we assign to data_g the fluctuation(n_w/barn - n/ barn)
           // get_overdens(fftw_functions.field_external_marked,w_delta);
           //  get_overdens(fftw_functions.field_external,delta);
           fftw_functions.resize_fftw_vectors();
           fftw_functions.raw_sampling(pow(this->params._Lbox(),3));
           fftw_functions.get_parameters_estimator(true);
           real_prec nmean=static_cast<real_prec>(tracer_aux.Halo.size())/static_cast<real_prec>(this->params._NGRID()); // mean number of tracers in cells
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
           for(ULONG i=0;i<this->params._NGRID();++i) // Build the property-fluctuation as n_w/barn - n/ barn
                     fftw_functions.set_data_g(i, (static_cast<real_prec>(this->tracer_cat.field_external_marked[i])-static_cast<real_prec>(this->tracer_cat.field_external[i]))/static_cast<real_prec>(nmean));
           if (true==this->params._use_real_and_redshift_space())
             {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
               for(ULONG i=0;i<this->params._NGRID();++i) // Build the property-fluctuation as n_w/barn - n/ barn
             fftw_functions.set_data_g_rss(i, (static_cast<real_prec>(this->tracer_cat.field_external_marked_s[i])-static_cast<real_prec>(this->tracer_cat.field_external_s[i]))/static_cast<real_prec>(nmean));
             }
           this->tracer_cat.field_external_marked.clear(); this->tracer_cat.field_external_marked.shrink_to_fit();
           if(Nbb==Number_of_properties)
             this->tracer_cat.field_external.clear();this->tracer_cat.field_external.shrink_to_fit();
           if (true==this->params._use_real_and_redshift_space())
               if(Nbb==Number_of_properties)
                 this->tracer_cat.field_external_s.clear();this->tracer_cat.field_external_s.shrink_to_fit();
           real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
           fftw_functions.set_normal_power(pow(factor,-2));
           fftw_functions.set_shot_noise((variance_property -1.0)/static_cast<real_prec>(nmean));
           kvector_data.clear();
           kvector_data.shrink_to_fit();
           for(ULONG i=0;i<this->params._d_Nnp_data();i++)
             kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
           if (true==this->params._use_real_and_redshift_space())
             {
               this->pk.clear();
               this->pk.shrink_to_fit();
               this->pk.resize(this->params._d_Nnp_data(),0); // Real space
             }
           this->pk0.clear();
           this->pk0.shrink_to_fit();
           this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
           this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
           this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole
           this->modes_g.clear();
           this->modes_g.shrink_to_fit();
           this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
           this->pk_w.clear();
           this->pk_w.shrink_to_fit();
           this->pk_w.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
           if (false==this->params._use_real_and_redshift_space())
             {
               fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
               this->write_power_spectrum(false); // argument asks to write sigma or not
             }
           else
             {
               fftw_functions.get_power_spectrum_fkp(this->pk,this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
               this->File.write_to_file(file_pow_marked,this->kvector_data,this->pk, this->pk0,this->pk2,this->pk4,this->modes_g);
             }
           fftw_functions.free_fftw_vectors();
           real_prec lss_bias=0;
           ULONG kcount=0;
           int ik=1;//to avoid zeros
           while(this->kvector_data[ik]<this->params._kmax_tracer_bias())
             {
               lss_bias+=sqrt(this->pk0[ik]/power_dm[ik]);
               ++ik;
               kcount++;
             }
           lss_bias/=kcount;
           So.message_screen("Large scale marked bias =", lss_bias);
           ik=1;//to avoid zeros
           real_prec var_lss_bias=0;
           while(this->kvector_data[ik]<this->params._kmax_tracer_bias())
             {
               var_lss_bias+=pow(this->pk0[ik]/power_dm[ik]-lss_bias*lss_bias,2);
               ++ik;
               kcount++;
             }
           var_lss_bias=sqrt(var_lss_bias/kcount);
           blsm<<lss_bias<<"\t"<<var_lss_bias<<endl;
           So.DONE();
           blsm.close();
         }//closes if  if(true==used_prop[Nbb] && true==secondary_prop[Nbb])
         } //closes loop over properties
       }// closes if(true==get_marked)
       // ************************************************************************************************************
       // END OF MARKED POWER SECTION
       //if marked is false, we still would like to have the primary field to have cross power computed below:
       else{
           if(true==this->params._Get_cross_power_spectrum())
           {
           s_data_structure s_data_struct_g_marked;
           this->N_galaxy=tracer_aux.Halo.size();
           real_prec mean_density=static_cast<real_prec>(tracer_aux.Halo.size())/pow(this->params._Lbox(),3);
           s_data_struct_g_marked.mean_density=mean_density;
           fftw_functions.resize_fftw_vectors();
           if(false==this->params._use_real_and_redshift_space())
             {
               if (space_p=="redshift_space")
                 {
                   So.message_screen("Interpolating galaxy density field on a grid");
#ifdef _USE_SEVERAL_RANDOM_FILES_
                  this->tracer_cat.get_interpolated_density_field(false, "any",0); // this gives  fftw_functions.field_external
#else
                   this->tracer_cat.get_interpolated_density_field(false, "any"); // this gives  fftw_functions.field_external
#endif
               }
                else if (space_p=="real_space")
                {
                    So.message_screen("Interpolating galaxy density field on a grid");
                     this->tracer_cat.get_interpolated_density_field_real_space(false,"any"); // this gives  fftw_functions.field_external
                }
             }
           else if(true==this->params._use_real_and_redshift_space())
              this->tracer_cat.get_interpolated_density_field_real_and_redshift_space(false, "any");

            primary_field=this->tracer_cat.field_external;
           fftw_functions.free_fftw_vectors();
            }
         }// end else
       // ************************************************************************************************************
       // ************************************************************************************************************
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
         if(true==this->params._set_bins_equal_number_tracers() && tracer._Number_of_tracers_in_mass_bins(im_primary)>0)// if not, read those from the para files as have been provided
#elif defined _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
         if(true==this->params._set_bins_equal_number_tracers() && tracer._Number_of_tracers_in_vmax_bins(im_primary)>0)// if not, read those from the para files as have been provided
#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
       if(true==this->params._set_bins_equal_number_tracers() && tracer._Number_of_tracers_in_ph_bins(im_primary)>0)// if not, read those from the para files as have been provided
#endif
         {
         /*
           if(abs(pow(counter_m  - tracer._Number_of_tracers_in_mass_bins(im),2))>0)
           {
           cerr<<RED<<"WARNING: missmatch between number in mass bins at "<<__PRETTY_FUNCTION__<<", line "<<__LINE__<<RESET<<endl;
           cerr<<im<<"  "<< pow(10,this->params._MASSbins_min(im))<<"   "<< pow(10,this->params._MASSbins_max(im))<<"   Current: "<<counter_m<<"  Expected: "<<tracer._Number_of_tracers_in_mass_bins(im)<<endl;
           }
         */
           string nname;
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
        nname="massbin";
#endif
#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
         nname="vmaxbin";
#endif
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
         nname="nubin";
#endif
        So.message_screen("Identifying bins of secondary properties");
         string file_bins_log = this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".log";
         So.message_screen("Writing bining scheme in file ", file_bins_log );
         ofstream binn_l;binn_l.open(file_bins_log.c_str());
         ofstream binn;
         string file_bins;
         //***********************************************************************
         // Headers
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
         binn_l<<"#MASS bin (log10): "<<im_primary<<endl;
         binn_l<<"#logMass min = "<<params._MASSbins_min(im_primary)<<"\tlog Mass max = "<<params._MASSbins_max(im_primary)<<endl;
         binn_l<<"#Number of tracers in current mass bin: "<< tracer._Number_of_tracers_in_mass_bins(im_primary)<<endl;
         binn_l<<"#Number of bins of secondary properties: "<< Nvmax_ind<<"  (0-th: full sample)"<<endl;  // I choose Nvmax_ind as this will be always used.
         binn_l<<"#Number of tracers in bins of secondary properties: "<< tracer._Number_of_tracers_in_mass_bins(im_primary)/(Nvmax_ind-1)<<endl;
#endif
         // Headers
#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
         binn_l<<"#VMAX bin (log10): "<<im_primary<<endl;
         binn_l<<"#logVmax min = "<<params._VMAXbins_min(im_primary)<<"\tlog VMax max = "<<params._VMAXbins_max(im_primary)<<endl;
         binn_l<<"#Number of tracers in current vmax bin: "<< tracer._Number_of_tracers_in_vmax_bins(im_primary)<<endl;
         binn_l<<"#Number of bins of secondary properties: "<< Nmass_ind<<endl;  // I choose Nvmax_ind as this will be always used.
         binn_l<<"#Number of tracers in bins of secondary properties: "<< tracer._Number_of_tracers_in_vmax_bins(im_primary)/(Nvmax_ind-1)<<endl;
#endif
         // Headers
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
         binn_l<<"#NU bin (log10): "<<im_primary<<endl;
         binn_l<<"#logNu min = "<<params._PHbins_min(im_primary)<<"\tlog Nu max = "<<params._PHbins_max(im_primary)<<endl;
         binn_l<<"#Number of tracers in current nu bin: "<< tracer._Number_of_tracers_in_ph_bins(im_primary)<<endl;
         binn_l<<"#Number of bins of secondary properties: "<< Nvmax_ind<<endl;  // I choose Nvmax_ind as this will be always used.
         binn_l<<"#Number of tracers in bins of secondary properties: "<< tracer._Number_of_tracers_in_ph_bins(im_primary)/(Nvmax_ind-1)<<endl;
#endif
        //***********************************************************************
        // Let us try that there is at least one object per quartile
        if(tracer._Number_of_tracers_in_mass_bins(im_primary)/(this->params._Number_of_bins_equal_number_tracers()) >= 1 )
        {
#ifndef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
         if(true==bmass)
           {
             binn_l<<"#Mass:"<<endl;
             tracer_aux.get_intervals_equal_number_aux("_MASS_");
             So.message_screen("Writing bining scheme in file ", file_bins );
             for(int i=0;i<tracer_aux.params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._MASSbins_min(i)<<"  "<<tracer_aux.params._MASSbins_max(i)<<"\t"<< endl;
           }
#endif
#ifndef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
         if(true==bvmax)
           {

             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[1];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             binn_l<<"#Vmax:"<<endl;
             tracer_aux.get_intervals_equal_number_aux("_VMAX_");
             for(int i=0;i<tracer_aux.params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._VMAXbins_min(i)<<"  "<<tracer_aux.params._VMAXbins_max(i)<<"\t"<< endl;
             for(int i=0;i<tracer_aux.params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._VMAXbins_min(i)<<"  "<<tracer_aux.params._VMAXbins_max(i)<<"\t"<< endl;
            binn.close();
         }
#endif
         if(true==used_prop[label_rs])
           {
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[2];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
         binn_l<<"#RS:"<<endl;
         tracer_aux.get_intervals_equal_number_aux("_RS_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._RSbins_min(i)<<"  "<<tracer_aux.params._RSbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._RSbins_min(i)<<"  "<<tracer_aux.params._RSbins_max(i)<<endl;
         binn.close();
        }
         if(true==used_prop[label_rvir])
           {
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[3];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
         binn_l<<"#RVIR:"<<endl;
         tracer_aux.get_intervals_equal_number_aux("_RVIR_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._RVIRbins_min(i)<<"  "<<tracer_aux.params._RVIRbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._RVIRbins_min(i)<<"  "<<tracer_aux.params._RVIRbins_max(i)<<endl;
            binn.close();
         }
         if(true==used_prop[label_concentration])
           {
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[4];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
         binn_l<<"#CONCENTRATION:"<<endl;
         tracer_aux.get_intervals_equal_number_aux("_CONCENTRATION_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._CONCENTRATIONbins_min(i)<<"  "<<tracer_aux.params._CONCENTRATIONbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._CONCENTRATIONbins_min(i)<<"  "<<tracer_aux.params._CONCENTRATIONbins_max(i)<<endl;
         binn.close();
        }
         if(true==used_prop[label_spin])
           {
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[5];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
         binn_l<<"#SPIN:"<<endl;
         tracer_aux.get_intervals_equal_number_aux("_SPIN_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._SPINbins_min(i)<<"  "<<tracer_aux.params._SPINbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._SPINbins_min(i)<<"  "<<tracer_aux.params._SPINbins_max(i)<<endl;
         binn.close();
    }
         if(true==used_prop[label_spin_bullock])
           {
             binn_l<<"#SPIN_BULLOCK:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[6];
             binn.open(file_bins.c_str());
             So.message_screen("Writing bining scheme in file ", file_bins );
          tracer_aux.get_intervals_equal_number_aux("_SPIN_BULLOCK_");
          for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
            binn_l<<i<<"  "<<tracer_aux.params._SPINBULLOCKbins_min(i)<<"  "<<tracer_aux.params._SPINBULLOCKbins_max(i)<<endl;
          for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
            binn<<i<<"  "<<tracer_aux.params._SPINBULLOCKbins_min(i)<<"  "<<tracer_aux.params._SPINBULLOCKbins_max(i)<<endl;
          binn.close();
         }
         if(true==used_prop[label_vrms])
           {
         binn_l<<"#VRMS:"<<endl;
         file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[7];
         binn.open(file_bins.c_str());
         So.message_screen("Writing bining scheme in file ", file_bins );
         tracer_aux.get_intervals_equal_number_aux("_VRMS_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._VRMSbins_min(i)<<"  "<<tracer_aux.params._VRMSbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._VRMSbins_min(i)<<"  "<<tracer_aux.params._VRMSbins_max(i)<<endl;
         binn.close();
        }
         if(true==used_prop[label_virial])
           {
         binn_l<<"#VIRIAL:"<<endl;
         file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[8];
         binn.open(file_bins.c_str());
         So.message_screen("Writing bining scheme in file ", file_bins );
         tracer_aux.get_intervals_equal_number_aux("_VIRIAL_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._VIRIALbins_min(i)<<"  "<<tracer_aux.params._VIRIALbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._VIRIALbins_min(i)<<"  "<<tracer_aux.params._VIRIALbins_max(i)<<endl;
         binn.close();
        }
         if(true==used_prop[label_btoa])
           {
         binn_l<<"#B_TO_A:"<<endl;
         file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[9];
         binn.open(file_bins.c_str());
         So.message_screen("Writing bining scheme in file ", file_bins );
         tracer_aux.get_intervals_equal_number_aux("_BTOA_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._BTOAbins_min(i)<<"  "<<tracer_aux.params._BTOAbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._BTOAbins_min(i)<<"  "<<tracer_aux.params._BTOAbins_max(i)<<endl;
         binn.close();
         }
         if(true==used_prop[label_ctoa])
           {
         binn_l<<"#C_TO_A:"<<endl;
         file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[10];
         binn.open(file_bins.c_str());
         So.message_screen("Writing bining scheme in file ", file_bins );
         tracer_aux.get_intervals_equal_number_aux("_CTOA_");
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn_l<<i<<"  "<<tracer_aux.params._CTOAbins_min(i)<<"  "<<tracer_aux.params._CTOAbins_max(i)<<endl;
         for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
           binn<<i<<"  "<<tracer_aux.params._CTOAbins_min(i)<<"  "<<tracer_aux.params._CTOAbins_max(i)<<endl;
         binn.close();
            }
         if(true==mach)
           {
             binn_l<<"#MACH:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[11];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_MACH_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._MACHbins_min(i)<<"  "<<tracer_aux.params._MACHbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._MACHbins_min(i)<<"  "<<tracer_aux.params._MACHbins_max(i)<<endl;
             binn.close();
           }
         if(true==tbias)
           {
             binn_l<<"#BIAS:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[12];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_BIAS_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._BIASbins_min(i)<<"  "<<tracer_aux.params._BIASbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._BIASbins_min(i)<<"  "<<tracer_aux.params._BIASbins_max(i)<<endl;
             binn.close();
           }
         if(true==tlc)
           {
             binn_l<<"#LC:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[13];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_LOCAL_OVERDENSITY_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._LCbins_min(i)<<"  "<<tracer_aux.params._LCbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._LCbins_min(i)<<"  "<<tracer_aux.params._LCbins_max(i)<<endl;
             binn.close();
           }
         if(true==bta)
           {
             binn_l<<"#TA:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[14];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_TIDAL_ANISOTROPY_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._TAbins_min(i)<<"  "<<tracer_aux.params._TAbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._TAbins_min(i)<<"  "<<tracer_aux.params._TAbins_max(i)<<endl;
             binn.close();
           }
#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
         if(true==bph)
           {
             binn_l<<"#PH:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[15];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_PEAK_HEIGHT_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._PHbins_min(i)<<"  "<<tracer_aux.params._PHbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._PHbins_min(i)<<"  "<<tracer_aux.params._PHbins_max(i)<<endl;
             binn.close();
           }
#endif

         if(true==dach)
           {
             binn_l<<"#DACH:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[16];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_DACH_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._DACHbins_min(i)<<"  "<<tracer_aux.params._DACHbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._DACHbins_min(i)<<"  "<<tracer_aux.params._DACHbins_max(i)<<endl;
             binn.close();
           }
         if(true==dm_local)
           {
             binn_l<<"#lDM:"<<endl;
             file_bins=this->params._Output_directory()+"BINNING_"+nname+to_string(im_primary)+"_cwt"+to_string(icw)+"_p"+to_string(this->params._unitsim_plabel())+".txt"+prop_name[17];
             So.message_screen("Writing bining scheme in file ", file_bins );
             binn.open(file_bins.c_str());
             tracer_aux.get_intervals_equal_number_aux("_LOCALDM_");
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn_l<<i<<"  "<<tracer_aux.params._LOCALDMbins_min(i)<<"  "<<tracer_aux.params._LOCALDMbins_max(i)<<endl;
             for(int i=0;i<params._Number_of_bins_equal_number_tracers()+1;++i)
               binn<<i<<"  "<<tracer_aux.params._LOCALDMbins_min(i)<<"  "<<tracer_aux.params._LOCALDMbins_max(i)<<endl;
             binn.close();
           }
         binn_l.close();
         So.DONE();
        }//closes if Ntracers/Nbins > Nbins so that at least there is one per quartile
         cout<<endl;
          } //closes  if(true==this->params._set_bins_equal_number_tracers())
       // ----------------------------------------------------- INICIA LOOP SOBRE SECONDARY PROPERTIES
#ifdef  _use_two_quartiles_
#ifdef  _ONLY_CWT_AND_PRIMARY_HBIAS_  // if this is defined, we only go through bin0 for all secondary properties
       vector<int>used_bins{0};
#else
       vector<int>used_bins{0,1,4};
#endif
       vector<int> initial_bin(Number_of_properties,0);//where to start. Initialized to 0: it is updated below if the property is used and if its secondary
       vector<int> N_bins_s(Number_of_properties, 1); // number pf bins. Initialized to 1: it is updated below if the property is used and if its secondary
#endif


    int Nbb=0;
#ifndef _ONLY_CWT_AND_PRIMARY_HBIAS_  // if this is defined, we only go through bin0 for all secondary properties


#ifdef _ASSEMBLY_BIAS_MASS_ONLY_
       int counter_sec=0;// this counter helps to do all 0 bins only once
       // Loop over the total number of properties, regardless whethere these are prim or sec (it will be specified below with "if" statements).
       for(int Nbb=0;Nbb<Number_of_properties;++Nbb)
         {
           if(true==secondary_prop[Nbb] && true==used_prop[Nbb])//only valid if these are asked to be used and secondary
         {
            counter_sec++;
#ifndef _use_two_quartiles_
           switch(Nbb)
             {
           Será cuestión de tiempo para verlo.
               //case 0 missing here
             case(1): Nvmax_ind=Nbins_vmax;Nspin_ind=1;Nrs_ind=1;Nvirial_ind=1; Nvrms_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nmach_ind=1; Nbias_ind=1; Nlc_ind=1;Nta_ind=1;Nph_ind=1; break;
             case(2): Nrs_ind=Nbins_rs;Nvmax_ind=1;Nspin_ind=1;Nvirial_ind=1; Nvrms_ind=1; Nbtoa_ind=1;Nctoa_ind=1;Nmach_ind=1; Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(3): Nspin_ind=Nbins_spin;Nvmax_ind=1;Nrs_ind=1;Nvirial_ind=1 ; Nvrms_ind=1;Nbtoa_ind=1;Nctoa_ind=1;Nmach_ind=1; Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(4): Nvrms_ind=Nbins_vrms;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1;Nmach_ind=1;Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(5): Nvirial_ind=Nbins_virial;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1; Nbtoa_ind=1;Nctoa_ind=1;Nmach_ind=1;Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(6): Nbtoa_ind=Nbins_btoa;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nctoa_ind=1;Nmach_ind=1;Nbias_ind=1; Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(7): Nctoa_ind=Nbins_ctoa;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nmach_ind=1; Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(8): Nmach_ind=Nbins_mach;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nbias_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(9): Nbias_ind=Nbins_bias;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nmach_ind=1;Nlc_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(10): Nlc_ind=Nbins_lc;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nmach_ind=1;Nbias_ind=1;Nta_ind=1;Nph_ind=1;break;
             case(11): Nta_ind=Nbins_ta;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nmach_ind=1;Nbias_ind=1;Nlc_ind=1;Nph_ind=1;break;
             case(12): Nph_ind=Nbins_ph;Nvmax_ind=1;Nrs_ind=1;Nspin_ind=1; Nvrms_ind=1;Nvirial_ind=1; Nbtoa_ind=1;Nctoa_ind=1; Nmach_ind=1;Nbias_ind=1;Nlc_ind=1;Nta_ind=1;break;
             case(13):aca falta arreglar
             }
#else
           // These lines are meant to do the loops in the secondary propperties covering only the 0 (full) ,1 and 4 quartiles.
           // This first if helps to do all bins 0 only for the first used secondary property.
#ifndef  _ONLY_CWT_AND_PRIMARY_HBIAS_
            if(counter_sec==1)
            {
                     initial_bin[Nbb] = 0;    // Where to start. If the property is not used, we wull never pass through its label so no need for another if
                     N_bins_s[Nbb] = used_bins.size(); // If in the loop over secondary properties you are varying a particular prop, use 0,1,2 bins. Otherwise, the number of bins is 1
            }
           else
               for(int Nab=0;Nab<Number_of_properties;++Nab)
            {
                    initial_bin[Nab] = (Nab==Nbb ? 1 : 0);    // Where to start. If the property is not used, we wull never pass through its label so no need for another if
                    N_bins_s[Nab] = (Nab==Nbb ? used_bins.size(): 1); // If in the loop over secondary properties you are varying a particular prop, use 0,1,2 bins. Otherwise, the number of bins is 1
            }
#endif
#endif // end for _use_two_quartiles_
#endif // end for _ASSEMBLY_BIAS_MASS_ONLY_
#endif // for _ONLY_CWT_AND_PRIMARY_HBIAS_
#ifndef _ONLY_CWT_AND_PRIMARY_HBIAS_  // if this is defined, we only go through bin0 for all secondary properties
#ifdef _use_two_quartiles_
#ifndef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
           for(int qimass=initial_bin[label_mass]; qimass< N_bins_s[label_mass];++qimass) // loop over bins of mass *only if it is secondary*
#endif
#ifndef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
             for(int qiv=initial_bin[label_vmax]; qiv< N_bins_s[label_vmax];++qiv) // // loop over bins of vmax only if it is secondary
#endif
               for(int qir=initial_bin[label_rs]; qir< N_bins_s[label_rs];++qir) //loop over bins of Rs
                   for(int qirv=initial_bin[label_rvir]; qirv< N_bins_s[label_rvir];++qirv) //loop over bins of Rvir
                       for(int qicon=initial_bin[label_concentration]; qicon< N_bins_s[label_concentration];++qicon) //loop over bins of concentration
             for(int qis=initial_bin[label_spin]; qis< N_bins_s[label_spin];++qis) // loop over bins of spin
                 for(int qisb=initial_bin[label_spin_bullock]; qisb< N_bins_s[label_spin_bullock];++qisb) // loop over bins of spin_bullock
               for(int qirm=initial_bin[label_vrms]; qirm< N_bins_s[label_vrms];++qirm) // loop over bins of VRms
                 for(int qivir=initial_bin[label_virial]; qivir< N_bins_s[label_virial];++qivir) // loop over bins of Virial
                   for(int qiba=initial_bin[label_btoa]; qiba< N_bins_s[label_btoa];++qiba) // loop over bins of b_to_a
                 for(int qica=initial_bin[label_ctoa]; qica< N_bins_s[label_ctoa];++qica) // loop over bins of c_to_a
                   for(int qima=initial_bin[label_mach]; qima< N_bins_s[label_mach];++qima) // loop over bins of mach
                     for(int qibi=initial_bin[label_bias]; qibi< N_bins_s[label_bias];++qibi) // loop over bins of bias
                       for(int qilc=initial_bin[label_lc]; qilc< N_bins_s[label_lc];++qilc) // loop over bins of bias
                     for(int qita=initial_bin[label_ta]; qita< N_bins_s[label_ta];++qita) // loop over bins of tidal_anisotropy
#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
                       for(int qiph=initial_bin[label_ph]; qiph< N_bins_s[label_ph];++qiph) // loop over bins of nu only if it is not primary
#endif
                         for(int qida=initial_bin[label_da]; qida< N_bins_s[label_da];++qida) // loop over bins of dach number
                             for(int qildm=initial_bin[label_dm_local]; qildm< N_bins_s[label_dm_local];++qildm) // loop over bins of local dm overdenisty (log of)
                           {
#ifndef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
                         int imass=used_bins[qimass];
#endif
#ifndef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
                         int iv=used_bins[qiv];
#endif
                         int ir= used_bins[qir];
                         int irv= used_bins[qirv];
                         int icon= used_bins[qicon];
                         int is= used_bins[qis];
                         int isb= used_bins[qisb];
                         int irm=used_bins[qirm];
                         int ivir=used_bins[qivir];
                         int iba=used_bins[qiba];
                         int ica=used_bins[qica];
                         int ima=used_bins[qima];
                         int ibi=used_bins[qibi];
                         int ilc=used_bins[qilc];
                         int ita=used_bins[qita];
#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
                         int iph=used_bins[qiph];
#endif
                         int ida=used_bins[qida];
                         int ildm=used_bins[qildm];

#else //else for _use_two_quartiles_
                         for(int iv=0; iv< Nvmax_ind;++iv) // loop over bins of vmax
                           for(int ir=0; ir< Nrs_ind;++ir) //loop over bins of Rs
                             for(int is=0; is< Nspin_ind;++is) // loop over bins of spin
                               for(int irm=0; irm< Nvrms_ind;++irm) // loop over bins of VRms
                             for(int ivir=0; ivir< Nvirial_ind;++ivir) // loop over bins of Virial
                               for(int iba=0; iba< Nbtoa_ind;++iba) // loop over bins of b_to_a
                                 for(int ica=0; ica< Nctoa_ind;++ica) // loop over bins of c_to_a
                                   for(int ima=0; ima< Nmach_ind;++ima) // loop over bins of mach
                                 for(int ibi=0; ibi< Nbias_ind;++ibi) // loop over bins of bias
                                   for(int ilc=0; ilc< Nlc_ind;++ilc) // loop over bins of bias
                                     for(int ita=0; ita< Nta_ind;++ita) // loop over bins of tidal_anisotropy
                                       for(int iph=0; iph< Nph_ind;++iph) // loop over bins of peak height
                                     for(int ida=0; iph< Nda_ind;++ida) // loop over bins of dach number
                                         for(int ildm=0; ildm< Ndmlocal_ind;++ildm) // loop over bins of dach number
                                       {
#endif // endif for _use_two_quartiles_
#else
    int imass=0;int ir=0;int iv=0;int irv=0;int icon=0; int is=0; int isb=0;int irm=0;int ivir=0; int iba=0; int ica=0; int ima=0; int ibi=0; int ilc=0; int ita=0; int ildm=0; int ida=0;int iph=0;

#endif  // endif for  _ONLY_CWT_AND_PRIMARY_HBIAS_  //
                                         So.message_screen("\tCWT",icw);
                                         So.message_screen("\tProperty label ",Nbb);
                                         So.message_screen("\tProperty name ",prop_name[Nbb]);
                                         cout<<endl;
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
                                         So.message_screen("\tMvir bin (PRIMARY)",im_primary);
#else
                                         So.message_screen("\t\tMvir bin",imass);
#endif
#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
                                         So.message_screen("\tVMAX (PRIMARY) bin",im_primary);
#else
                                         So.message_screen("\t\tVMAX bin",iv);
#endif
                                         if(this->params._i_rs_g()>0)
                                             So.message_screen("\t\tRs bin",ir);
                                         if(this->params._i_rvir_g()>0)
                                             So.message_screen("\t\tRvir bin",irv);
#ifdef _USE_CONCENTRATION_
                                         if(this->params._i_rs_g()>0 && this->params._i_rvir_g()>0)
                                             So.message_screen("\t\tConcentration bin",icon);
#endif
                                         if(this->params._i_spin_g()>0)
                                             So.message_screen("\t\tSpin bin ",is);
                                         if(this->params._i_spin_bullock_g()>0)
                                             So.message_screen("\t\tSpinB bin ",isb);
                                         if(this->params._i_vrms_g()>0)
                                             So.message_screen("\t\tVrms bin ",irm);
                                         if(this->params._i_virial_g()>0)
                                             So.message_screen("\t\tVirial bin ",ivir);
                                         if(this->params._i_b_to_a_g()>0)
                                             So.message_screen("\t\tb2a bin ",iba);
                                         if(this->params._i_c_to_a_g()>0)
                                             So.message_screen("\t\tc2a bin ",ica);
                                         if(this->params._Get_tracer_local_mach_number() ==true)
                                             So.message_screen("\t\tmach bin ",ima);
                                         if(this->params._Get_tracer_bias() ==true)
                                             So.message_screen("\t\tbias bin ",ibi);
                                         if(this->params._Get_local_overdensity() ==true)
                                             So.message_screen("\t\tLC bin ",ilc);
                                         if(this->params._Get_tidal_anisotropy_at_halo() ==true)
                                            So.message_screen("\t\tTA bin ",ita);
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
                                         So.message_screen("\tPH (PRIMARY) bin ",im_primary);
#else
                                         if(this->params._Get_peak_height_at_halo() ==true)
                                             So.message_screen("\t\tPH bin ",iph);
#endif
                                         if(this->params._Get_tracer_local_dach_number() ==true)
                                             So.message_screen("\t\tdach bin ",ida);
                                         if(this->params._Get_tracer_local_dm_density() ==true)
                                             So.message_screen("\t\tlog (delta dm +1) bin ",ildm);

#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
                                         int ph_index=iph;
                                         int vmax_index=iv;
                                         int mass_index=im_primary;
#elif defined _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
                                         int ph_index=iph;
                                         int vmax_index=im_primary;
                                         int mass_index=imass;
#elif defined _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
                                         int ph_index=im_primary;
                                         int vmax_index=iv;
                                         int mass_index=imass;
#endif


                                          string extra_info="_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);

                                          string extra_info_dm="_DM_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
#ifdef _USE_CWC_HALO_ANALYSIS_
                                         this->file_power=file_pow+extra_info;
                                         this->file_power_cross=file_pow_cross+extra_info;

                                         string file_field=this->params._Output_directory()+"halo_field_Nres_"+to_string(this->params._Nft())+extra_info;
#else
                                         if(false==this->params._set_bins_equal_number_tracers())
                                           {
                                         this->file_power=file_pow+"_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_massbin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_fixedbins";
                                         string file_field=this->params._Output_directory()+"halo_field_Nres_"+to_string(this->params._Nft())+"_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_fixedbins";
                                           }
                                         else
                                           {
                                         this->file_power=file_pow+"_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_massbin"+to_string(iph)+"_dabin"+"_ldmbin"+to_string(ildm)+to_string(ida);
                                         string file_field=this->params._Output_directory()+"halo_field_Nres_"+to_string(this->params._Nft())+"_massbin"+to_string(im)+"_vmaxbin"+to_string(iv)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(iph)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm);
                                           }
#endif
                                         string file_contours =this->params._Output_directory()+"halo_field_Nres_"+to_string(this->params._Nft())+extra_info;
                                         string file_hist = this->params._Output_directory()+"halo_field_Nres_"+to_string(this->params._Nft())+extra_info;
                                         ofstream blss;
                                         ofstream blssc;
#ifdef _USE_CWC_HALO_ANALYSIS_
                                         string file_lss_bias=this->params._Output_directory()+"lss_bias"+extra_info;
                                         string file_lss_bias_cross=this->params._Output_directory()+"lss_bias_cross"+extra_info;
                                         string file_lss_bias_object=this->params._Output_directory()+"lss_bias_object"+extra_info;
#else
                                         string file_lss_bias=this->params._Output_directory()+"lss_bias_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_phbin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm);
#endif
                                         blss.open(file_lss_bias.c_str());
                                         blssc.open(file_lss_bias_cross.c_str());
                                         Catalog tracer_final(this->params);
                                         ULONG in_new=0;  // counter for sub sample
                                         real_prec bias_from_objects=0;
                                         vector<real_prec>bias_ob;
                                         // We select now properties in bins of secondary properties. Selection in CWT and Mass has been already imposed in the tracer_aux catalog.
                                         for(ULONG i=0;i<tracer_aux.Halo.size();++i)
                                         {
#ifndef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
                                           if(log10(tracer_aux.Halo[i].mass)>= tracer_aux.params._MASSbins_min(mass_index) && log10(tracer_aux.Halo[i].mass) < tracer_aux.params._MASSbins_max(mass_index))
#endif
#ifndef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
                                         if(log10(tracer_aux.Halo[i].vmax)>= tracer_aux.params._VMAXbins_min(vmax_index) && log10(tracer_aux.Halo[i].vmax)< tracer_aux.params._VMAXbins_max(vmax_index))
#endif
                                           if(tracer_aux.Halo[i].rs>= tracer_aux.params._RSbins_min(ir) && tracer_aux.Halo[i].rs< tracer_aux.params._RSbins_max(ir))
                                               if(tracer_aux.Halo[i].rvir>= tracer_aux.params._RVIRbins_min(irv) && tracer_aux.Halo[i].rvir< tracer_aux.params._RVIRbins_max(irv))
                                                   if(tracer_aux.Halo[i].concentration>= tracer_aux.params._CONCENTRATIONbins_min(icon) && tracer_aux.Halo[i].concentration< tracer_aux.params._CONCENTRATIONbins_max(icon))
                                               if(tracer_aux.Halo[i].spin>= tracer_aux.params._SPINbins_min(is) && tracer_aux.Halo[i].spin< tracer_aux.params._SPINbins_max(is))
                                                   if(tracer_aux.Halo[i].spin_bullock>= tracer_aux.params._SPINBULLOCKbins_min(isb) && tracer_aux.Halo[i].spin_bullock< tracer_aux.params._SPINBULLOCKbins_max(isb))
                                               if(tracer_aux.Halo[i].vrms>= tracer_aux.params._VRMSbins_min(irm) && tracer_aux.Halo[i].vrms< tracer_aux.params._VRMSbins_max(irm))
                                             if(tracer_aux.Halo[i].virial>= tracer_aux.params._VIRIALbins_min(ivir) && tracer_aux.Halo[i].virial< tracer_aux.params._VIRIALbins_max(ivir))
                                               if(tracer_aux.Halo[i].b_to_a>= tracer_aux.params._BTOAbins_min(iba) && tracer_aux.Halo[i].b_to_a< tracer_aux.params._BTOAbins_max(iba))
                                                 if(tracer_aux.Halo[i].c_to_a>= tracer_aux.params._CTOAbins_min(ica) && tracer_aux.Halo[i].c_to_a< tracer_aux.params._CTOAbins_max(ica))
                                                   if(tracer_aux.Halo[i].mach_number>= tracer_aux.params._MACHbins_min(ima) && tracer_aux.Halo[i].mach_number< tracer_aux.params._MACHbins_max(ima))
                                                 if(tracer_aux.Halo[i].bias>= tracer_aux.params._BIASbins_min(ibi) && tracer_aux.Halo[i].bias< tracer_aux.params._BIASbins_max(ibi))
                                                   if(tracer_aux.Halo[i].local_overdensity>= tracer_aux.params._LCbins_min(ilc) && tracer_aux.Halo[i].local_overdensity< tracer_aux.params._LCbins_max(ilc))
                                                     if(tracer_aux.Halo[i].tidal_anisotropy>= tracer_aux.params._TAbins_min(ita) && tracer_aux.Halo[i].tidal_anisotropy< tracer_aux.params._TAbins_max(ita))
#ifndef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
                                                       if(tracer_aux.Halo[i].peak_height>= tracer_aux.params._PHbins_min(ph_index) && tracer_aux.Halo[i].peak_height< tracer_aux.params._PHbins_max(ph_index))
#endif
                                                     if(tracer_aux.Halo[i].dach_number>= tracer_aux.params._DACHbins_min(ida) && tracer_aux.Halo[i].dach_number< tracer_aux.params._DACHbins_max(ida))
                                                         if(tracer_aux.Halo[i].local_dm>= tracer_aux.params._LOCALDMbins_min(ida) && tracer_aux.Halo[i].local_dm< tracer_aux.params._LOCALDMbins_max(ida))
                                                       {
                                                         tracer_final.Halo.push_back(s_Halo());//Push back new subject created with default constructor.
                                                         tracer_final.Halo[in_new].coord1=tracer_aux.Halo[i].coord1;
                                                         tracer_final.Halo[in_new].coord2=tracer_aux.Halo[i].coord2;
                                                         tracer_final.Halo[in_new].coord3=tracer_aux.Halo[i].coord3;
                                                         tracer_final.Halo[in_new].mass=tracer_aux.Halo[i].mass;
                                                         tracer_final.Halo[in_new].vmax=tracer_aux.Halo[i].vmax;
                                                         tracer_final.Halo[in_new].rs=tracer_aux.Halo[i].rs;
                                                         tracer_final.Halo[in_new].rvir=tracer_aux.Halo[i].rvir;
                                                         tracer_final.Halo[in_new].concentration=tracer_aux.Halo[i].concentration;
                                                         tracer_final.Halo[in_new].spin=tracer_aux.Halo[i].spin;
                                                         tracer_final.Halo[in_new].spin_bullock=tracer_aux.Halo[i].spin_bullock;
                                                         tracer_final.Halo[in_new].vrms=tracer_aux.Halo[i].vrms;
                                                         tracer_final.Halo[in_new].virial=tracer_aux.Halo[i].virial;
                                                         tracer_final.Halo[in_new].b_to_a=tracer_aux.Halo[i].b_to_a;
                                                         tracer_final.Halo[in_new].c_to_a=tracer_aux.Halo[i].c_to_a;
                                                         tracer_final.Halo[in_new].mach_number=tracer_aux.Halo[i].mach_number;
                                                         tracer_final.Halo[in_new].bias=tracer_aux.Halo[i].bias;
                                                         tracer_final.Halo[in_new].local_overdensity=tracer_aux.Halo[i].local_overdensity;
                                                         tracer_final.Halo[in_new].tidal_anisotropy=tracer_aux.Halo[i].tidal_anisotropy;
                                                         tracer_final.Halo[in_new].peak_height=tracer_aux.Halo[i].peak_height;
                                                         tracer_final.Halo[in_new].dach_number=tracer_aux.Halo[i].dach_number;
                                                         bias_from_objects+=tracer_aux.Halo[i].bias;
                                                         bias_ob.push_back(tracer_aux.Halo[i].bias);
                                                         tracer_final.Halo[in_new].vel1=tracer_aux.Halo[i].vel1;
                                                         tracer_final.Halo[in_new].vel2=tracer_aux.Halo[i].vel2;
                                                         tracer_final.Halo[in_new].vel3=tracer_aux.Halo[i].vel3;
                                                         tracer_final.Halo[in_new].local_dm=tracer_aux.Halo[i].local_dm;
                                                         in_new++;
                                                       }
                                             }
                                         cout<<endl;
                                         So.message_screen("\tNumber of tracers selected in bin of secondary property =",in_new);

                             //            string file_ntracers = this->params._Output_directory()+"Ntracers_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_phbin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
                             //            So.message_screen("Saving number of tracers in file", file_ntracers);
                             //            ofstream ntr; ntr.open(file_ntracers.c_str());
                             //            ntr<<in_new<<endl;
                             //            ntr.close();
                             //            So.DONE();
                                         tracer_final.set_NOBJS(in_new);// IMPORTANT. OTHERWISE METHODS OF tracer_aux WILL SEE DEFAULT VALUE NOBS=0
                                         cout<<endl;
                                         if(in_new>0)
                                           {
                                         // ************************************************************************************************ //
                                             if(true==this->params._Get_tracer_bias())
                                             {
                                                 So.message_screen("\tComputing population bias and writting in file =", file_lss_bias_object);
                                                  if(in_new>0)
                                                     bias_from_objects/=static_cast<real_prec>(in_new);
                                                  real_prec e_bias=sqrt(get_var(bias_from_objects,bias_ob));
                                                 ofstream nb; nb.open(file_lss_bias_object.c_str());
                                                 nb<<bias_from_objects<<"\t"<<e_bias<<endl;
                                                 nb.close();
                                                 So.DONE();
                                                 cout<<endl;
                                             }
                                             if(true==this->params._Get_pearson_coefficient())
                                            {
                                             // *Compute the pearson coefficients among halo proeprties. This is the intrinsic correlation *//
                                         ofstream pearsout;
                                         ofstream pearsout_b;
                                         string file_pearson_intrinsic;
                                         string file_pearson_bias;
                                         if(false==this->params._set_bins_equal_number_tracers())
                                           {
                                             file_pearson_intrinsic=this->params._Output_directory()+"pearson_coef_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_ldmbin"+to_string(ildm)+"_fixedbins";

                                           }
                                         else
                                         {
#ifdef _USE_CWC_HALO_ANALYSIS_
                                           file_pearson_intrinsic=this->params._Output_directory()+"pearson_coef_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
                                           file_pearson_bias=this->params._Output_directory()+"bias_pearson_coef_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
#else
                                         file_pearson_intrinsic=this->params._Output_directory()+"pearson_coef_massbin"+to_string(im)+"_vmaxbin"+to_string(iv)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(iph)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm);
#endif
                                         }
                                         So.message_screen("\tComputing Pearson coefficient in file", file_pearson_intrinsic );
                                         vector<real_prec>b_p;
                                         pearsout.open(file_pearson_intrinsic.c_str());
                                         for(int ka=0;ka<prop_name.size();++ka)
                                           {
                                             for(int kb=ka+1;kb<prop_name.size();++kb)
                                               {
                                                 real_prec pearson_c=tracer_final.pearson_correlation(prop_name[ka],prop_name[kb]);
                                                 pearsout<<prop_name[ka]<<"\t"<<prop_name[kb]<<"\t"<<pearson_c<<endl;
                                                 if(prop_name[kb]=="_BIAS_" || prop_name[ka]=="_BIAS_")
                                                    b_p.push_back(pearson_c);
                                             }
                                           }
                                         pearsout.close();

                                         pearsout_b.open(file_pearson_bias.c_str());
                                         for(int ka=0;ka<b_p.size();++ka)
                                             pearsout_b<<ka<<"  "<<b_p[ka]<<endl;
                                         pearsout_b.close();
                                         b_p.clear();b_p.shrink_to_fit();
                                         So.DONE();
                                         }// close if get pearson
                                          if(true==this->params._Get_spearman_coefficient())
                                           {
                                              ofstream spearman;
                                              ofstream spearman_b;
                                              string file_spearman_intrinsic;
                                              string file_spearman_bias;
                                              file_spearman_intrinsic=this->params._Output_directory()+"spearman_coef_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
                                              file_spearman_bias=this->params._Output_directory()+"bias_spearman_coef_massbin"+to_string(mass_index)+"_vmaxbin"+to_string(vmax_index)+"_rsbin"+to_string(ir)+"_rvirbin"+to_string(irv)+"_concentrationbin"+to_string(icon)+"_spinbin"+to_string(is)+"_spinbullockbin"+to_string(isb)+"_vrmsbin"+to_string(irm)+"_virialbin"+to_string(ivir)+"_btoabin"+to_string(iba)+"_ctoabin"+to_string(ica)+"_machbin"+to_string(ima)+"_biasbin"+to_string(ibi)+"_lcbin"+to_string(ilc)+"_tabin"+to_string(ita)+"_nubin"+to_string(ph_index)+"_dabin"+to_string(ida)+"_ldmbin"+to_string(ildm)+"_cwc"+to_string(icw);
                                              vector<real_prec>b_p;
                                              So.message_screen("\tComputing Spearman-rank coefficient in file", file_spearman_intrinsic );
                                              for(int ka=0;ka<prop_name.size();++ka)
                                                 tracer_final.Get_Ranked_Props(prop_name[ka]);
                                              spearman.open(file_spearman_intrinsic.c_str());
                                              for(int ka=0;ka<prop_name.size();++ka)
                                                {
                                                  for(int kb=ka+1;kb<prop_name.size();++kb)
                                                    {
                                                      real_prec pearson_c=tracer_final.spearman_correlation(prop_name[ka],prop_name[kb]);
                                                      spearman<<prop_name[ka]<<"\t"<<prop_name[kb]<<"\t"<<pearson_c<<endl;
                                                      if(prop_name[kb]=="_BIAS_" || prop_name[ka]=="_BIAS_")
                                                        b_p.push_back(pearson_c);
                                                    }
                                                 }
                                               spearman.close();
                                               spearman_b.open(file_spearman_bias.c_str());
                                               for(int ka=0;ka<b_p.size();++ka)
                                                  spearman_b<<ka<<"  "<<b_p[ka]<<endl;
                                               spearman_b.close();
                                               b_p.clear();b_p.shrink_to_fit();
                                               tracer_final.Halo_ranked.clear();
                                               tracer_final.Halo_ranked.shrink_to_fit();
                                          }
                                          // ************************************************************************************************ //
                                        if(this->params._Get_PCA())
                                             {
					       So.message_screen("Doing PCA");
					       tracer_final.PCA(prop_name, used_prop,extra_info, false);
					       So.DONE();
					     }
			     
					     if(true==this->params._Get_power_spectrum())
					       {			      
					      this->File.input_type=this->params._input_type();
					      this->N_galaxy=in_new;
					      real_prec mean_density=static_cast<real_prec>(in_new)/pow(this->params._Lbox(),3);
					      s_data_structure s_data_struct_g;
					      s_data_struct_g.mean_density=mean_density;
					      tracer_final.Halo.clear();tracer_final.Halo.shrink_to_fit();
					      this->fftw_functions.set_n_gal(in_new);
					      this->fftw_functions.resize_fftw_vectors();
#ifdef _VERBOSE_POWER_
					      this->fftw_functions.write_fftw_parameters();
					      So.message_screen("Interpolating galaxy density field on a grid");
#endif
					      if (false==this->params._use_real_and_redshift_space())
						{
						  if (space_p=="redshift_space")
#ifdef _USE_SEVERAL_RANDOM_FILES_
                            tracer_final.get_interpolated_density_field(false, "any",0);
#else
                              tracer_final.get_interpolated_density_field(false, "any");
#endif

                          else if (space_p=="real_space")
                            tracer_final.get_interpolated_density_field_real_space(false, "any");
						}
                        else if (true==this->params._use_real_and_redshift_space())
                          tracer_final.get_interpolated_density_field_real_and_redshift_space(false, "any");

					      // *****************************************************************************************
					      if(true==this->params._write_files_for_histograms())
						{
                          this->File.write_array(file_field,this->tracer_cat.field_external);
                          vector<real_prec>aux_field (this->tracer_cat.field_external.size(),0);
                          get_overdens(this->tracer_cat.field_external,aux_field);
						  vector<real_prec>Vaux(Nbins*Nbins,0);
						  get_2d_histogram(this->params._ldelta_X_min(),this->params._ldelta_X_max(),this->params._ldelta_Y_min(),this->params._ldelta_Y_max(),Nbins, DM_DEN_FIELD, aux_field,Vaux,true);
						  this->File.write_array(file_hist,Vaux);
						  this->mcmc.get_contour_levels(file_contours+"_contour_levels",Nbins, Vaux);
						  Vaux.clear();Vaux.shrink_to_fit();
						}
					      // *****************************************************************************************
					      So.message_screen("\tComputing Prop-power");
					      this->fftw_functions.raw_sampling(pow(this->params._Lbox(),3));
					      this->fftw_functions.get_parameters_estimator(true);
					      if(true==this->params._use_real_and_redshift_space())
                            this->fftw_functions.get_fluctuation(true);
                         else
                            this->fftw_functions.get_fluctuation();
                         // *****************************************************************************************
                         kvector_data.clear(); // Bins here are linear by default
                         kvector_data.shrink_to_fit();
                         for(ULONG i=0;i<this->params._d_Nnp_data();i++)
                           kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
                         if (true==this->params._use_real_and_redshift_space())
                         {
                            this->pk.clear();
                            this->pk.shrink_to_fit();
                            this->pk.resize(this->params._d_Nnp_data(),0); // Real space
                         }

                         this->pk0.clear();
                         this->pk0.shrink_to_fit();
                         this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
                         this->pk2.resize(this->params._d_Nnp_data(),0); //Quadrupole
                         this->pk4.resize(this->params._d_Nnp_data(),0); //Hexadecapole

                         this->modes_g.clear();
                         this->modes_g.shrink_to_fit();
                         this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
                         this->pk_w.clear();
                         this->pk_w.shrink_to_fit();
                         this->pk_w.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
                         if (false==this->params._use_real_and_redshift_space())
                           {
                             fftw_functions.get_power_spectrum_fkp(this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
                             this->write_power_spectrum(false); // argument asks to write sigma or not
                           }
                         else
                        {
                          fftw_functions.get_power_spectrum_fkp(this->pk,this->pk0,this->pk2,this->pk4,this->pk_w,this->pkk,this->pmk,this->modes_g);
                          this->File.write_to_file(this->file_power,this->kvector_data,this->pk, this->pk0,this->pk2,this->pk4,this->modes_g);
                        }
					 
                     if(true==this->params._Get_cross_power_spectrum())
					   {
					     // Measure cross correlation between halo full population and sub-population
					     So.message_screen("\tComputing cross power spectrum between primary and secondary selected tracers");
                    //     this->compute_cross_power_spectrum_grid(false, this->fftw_functions.field_external, primary_field, true);
					     // Measure cross correlation between DM field and sub-population
					     this->file_power_cross=file_pow_cross+extra_info_dm;
                         this->compute_cross_power_spectrum_grid(true, DM_DEN_FIELD, tracer_final.field_external,true);
					     
                         pair<real_prec, real_prec>lss_bias_cpower=this->get_lss_cross_bias(this->pk0, power_dm, this->kvector_data, this->modes_g, 2, this->params._kmax_tracer_bias()); // 0.08
                         pair<real_prec, real_prec>lss_bias_cpower2=this->get_lss_cross_bias(this->pk0, power_dm, this->kvector_data, this->modes_g, 2, 0.06);
                         pair<real_prec, real_prec>lss_bias_cpower3=this->get_lss_cross_bias(this->pk0, power_dm, this->kvector_data, this->modes_g, 2, 0.04);

                         So.message_screen("\tLarge scale bias from cross power =", lss_bias_cpower.first);
                         blssc<<lss_bias_cpower.first<<"\t"<<lss_bias_cpower.second<<"\t"<<lss_bias_cpower2.first<<"\t"<<lss_bias_cpower2.second<<"\t"<<lss_bias_cpower3.first<<"\t"<<lss_bias_cpower3.second<<endl;
					     So.DONE();
					     blssc.close();
					   }
					 // *****************************************************************************************
                                         if(Nbb !=Number_of_properties) //this means, free memory except the last pass, for it will be done directly by the destructor and it compaints if it finds nothing to free
                                           fftw_functions.free_fftw_vectors();
					 // *****************************************************************************************
            // *****************************************************************************************
                                         // here we need to cpmpute kmax as Min(this->params._kmax_tracer_bias() and k such that P=0.1Pshot)
					 
                                         real_prec shot_noise=pow(this->params._Lbox(),3)/static_cast<real_prec>(in_new); // Poisson Shot noise
					 
                                         int ik=1;//to avoid zero mode
                                         while(this->pk[ik]>=SN_TOLERANCE_HBIAS*shot_noise)
                                           ++ik;
					 
                                         int initial_kmode=2;
                                         real_prec kmax_sn=min(this->params._kmax_tracer_bias(),this->kvector_data[ik]);
                                         real_prec kmax_sn2=min(static_cast<real_prec>(0.06),this->kvector_data[ik]);
                                         real_prec kmax_sn3=min(static_cast<real_prec>(0.04),this->kvector_data[ik]);
                                         if(kmax_sn>=this->kvector_data[initial_kmode])// this prevents that bias is computed from a shot-noise dominated signal
                                           {
					     
                                             this->So.message_screen("\tMaximum wave-number before shot-noise dominance:",this->kvector_data[ik]);
                                             this->So.message_screen("\tShot-noise:", shot_noise);
                                             this->So.message_screen("\tTolerance power ", SN_TOLERANCE_HBIAS*shot_noise);
                                             this->So.message_screen("\tMaximum wave-number for bias calculation ", kmax_sn);
					     
                                             pair<real_prec, real_prec>lss_bias_power=this->get_lss_bias(this->pk, power_dm, this->kvector_data, this->modes_g, initial_kmode, kmax_sn);
                                             pair<real_prec, real_prec>lss_bias_power2=this->get_lss_bias(this->pk, power_dm, this->kvector_data, this->modes_g, initial_kmode, kmax_sn2);
                                             pair<real_prec, real_prec>lss_bias_power3=this->get_lss_bias(this->pk, power_dm, this->kvector_data, this->modes_g, initial_kmode, kmax_sn3);
                                             So.message_screen("\tLarge scale bias =", lss_bias_power.first);
                                             blss<<lss_bias_power.first<<"\t"<<lss_bias_power.second<<"\t"<<lss_bias_power2.first<<"\t"<<lss_bias_power2.second<<"\t"<<lss_bias_power3.first<<"\t"<<lss_bias_power3.second<<endl;
                                             So.DONE();
                                             blss.close();
                                           }
					 else
                                           this->So.message_screen("\t\tP(k) dominated by shot-noise. LSS bias is not computed");
					 cout<<endl;
                                         }//cloess do power, andres
                                           }// closes if number of tracers>0

#ifndef _ONLY_CWT_AND_PRIMARY_HBIAS_  // if this is defined, we only go through bin0 for all secondary properties
                                       }// closes loop of bins or cuts, in particualr, loop over bin of  c_to_a
#ifdef _ASSEMBLY_BIAS_MASS_ONLY_
                           }//closes if properties are seconday and used
         }//closes loop over N bins of properties
#endif
#endif
         }//closes loop over primary properties
#ifdef _USE_CWC_HALO_ANALYSIS_
       }//closes loop over CTW
#endif
    this->So.message_time(time_POWER);
 }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::compute_marked_power_spectrum_grid(const vector<real_prec> &data_in,const vector<real_prec> &data_in_MW)
 {

#ifdef _VERBOSE_POWER_
   this->So.enter(__PRETTY_FUNCTION__);
#endif

#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif

   fftw_functions.resize_fftw_vectors();

   real_prec ngal_new=get_nobjects(data_in);

   real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     fftw_functions.data_g[i]=(static_cast<real_prec>(data_in_MW[i])-static_cast<real_prec>(data_in[i]))/static_cast<real_prec>(nmean);

   this->params.set_ngal_delta(ngal_new);
   fftw_functions.set_n_gal(ngal_new);

   real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
   fftw_functions.set_normal_power(pow(factor,-2));
   fftw_functions.shot_noise=(this->var_prop-1.0)*static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);

   kvector_data.resize(this->params._d_Nnp_data(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<this->params._d_Nnp_data();i++)
     kvector_data[i]=this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5);

   this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
   this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance

#if !defined _USE_BIAS_OBJECT_TO_OBJECT_ || !defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
   fftw_functions.power_spectrum(this->pk0,this->modes_g);
#else
   //we're done, no pwoer haja
#endif

#ifdef _USE_GNUPLOT_POWER_
   this->gp.plot_power_spectrum(this->kvector_data,this->pk0);
#endif
 }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::compute_power_spectrum_grid()
 {
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   vector<real_prec> data_in(this->params._NGRID(),0);
   So.message_screen("Reading input file on a mesh, with precision set by PrecType_Y directive");
   this->File.read_array_t<PrecType_Y>(this->params._delta_grid_file(), data_in);
   fftw_functions.resize_fftw_vectors();
#ifdef _FULL_VERBOSE_
   fftw_functions.write_fftw_parameters();
#endif
#ifdef _NCUTS_POWER_
   string ofile=this->file_power;
   for(int Ni=0;Ni<N_MAX_OCCUPATION;Ni+=2)
     {
#endif
#ifdef _NCUTS_POWER_
       if(this->params._input_type()=="density_grid")
     for(ULONG i=0;i<this->params._NGRID();++i)
       if(data_in[i]<Ni)
         data_in[i]=0;
       else
         if(this->params._input_type()=="dela_grid")
           for(ULONG i=0;i<this->params._NGRID();++i)
         if(data_in[i]<Ni)
           data_in[i]=-1;
#endif
      real_prec ngal_new=get_nobjects(data_in);
       if(this->params._input_type()=="density_grid")
#ifdef _VERBOSE_POWER_
     So.message_screen("Number of objects =",ngal_new);
#endif
       if(this->params._input_type()=="density_grid")
     {
       real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
       So.message_screen("Mean of input field", ngal_new);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<this->params._NGRID();++i)
             fftw_functions.data_g[i]=(static_cast<real_prec>(data_in[i])/static_cast<real_prec>(nmean))-1.;
     }
       else
     {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<this->params._NGRID();++i)
             fftw_functions.data_g[i]=data_in[i];
       So.message_screen("Mean of input field", get_mean(data_in));


     }
       this->params.set_ngal_delta(ngal_new);
       fftw_functions.set_n_gal(ngal_new);
       real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
       fftw_functions.set_normal_power(pow(factor,-2));
       fftw_functions.shot_noise=0;
       if(true==this->params._SN_correction())
         fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
#ifdef _VERBOSE_POWER_
       So.message_screen("Shot Noise =",fftw_functions.shot_noise);
       So.message_screen("Normalization =",fftw_functions._normal_power());
#endif
       this->kvector_data.resize(this->params._d_Nnp_data(), 0);
       if("linear" == this->params._type_of_binning())
     {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(int i=0;i<this->params._d_Nnp_data();i++)
         this->kvector_data[i]=this->params._d_kmin()+ this->params._d_DeltaK_data()*(i+0.5);

     }
       else if("log"==this->params._type_of_binning())
     {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(int i=0;i<this->params._d_Nnp_data();i++)
         this->kvector_data[i]=this->params._d_kmin()*pow(10,(i-0.5)*this->params._d_Deltal());
     }
       this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
       this->modes_g.resize(this->params._d_Nnp_data(),0); //Needed in case we use the Veff for the variance
#if  !defined (_USE_CROSS_CORRELATION_CONF_SPACE_)
       fftw_functions.power_spectrum(this->pk0,this->modes_g);
#endif
#ifdef _NCUTS_POWER_
       this->file_power=ofile+"_Ncuts"+to_string(Ni);
#endif
      write_power_and_modes();
#ifdef _USE_GNUPLOT_POWER_
       this->gp.plot_power_spectrum(this->kvector_data,this->pk0);
#endif
#ifdef _NCUTS_POWER_
     }
#endif
 }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::get_window_matrix_multipole()
{
    this->So.enter(__PRETTY_FUNCTION__);
    if(true==this->params._FKP_weight() && true==this->params._nbar_tabulated())
#ifdef _USE_OMP_
#pragma omp parlalel for
#endif
    for (ULONG i=0; i<this->N_random;++i)
      this->random_cat.Halo[i].weight1=1./(1+this->params._Pest()*this->random_cat.Halo[i].mean_density);
  else
#ifdef _USE_OMP_
#pragma omp parlalel for
#endif
    for (ULONG i=0; i<this->N_random;++i)
      this->random_cat.Halo[i].weight1=1.;
  this->alpha=fftw_functions._alpha();
  this->normal_power=fftw_functions._normal_power();
  So.DONE();
  // Compute here the squared of the distance to the origin d² for each random tracer
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tComputing r**2 for randoms");
#endif
  vector<real_prec>distances_cat(this->N_random,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0; i<this->N_random;++i)
    distances_cat[i]=pow(this->random_cat.Halo[i].coord1,2)+pow(this->random_cat.Halo[i].coord2,2)+pow(this->random_cat.Halo[i].coord3,2);
  So.DONE();
  this->kvector_data.clear();
  this->kvector_data.shrink_to_fit();
  for(int i=0;i<this->params._d_Nnp_data();i++)
    this->kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
  //GAUSS LEGENDRE FOR INTERGRATION WRT k
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tComputing GL-weights");
#endif
  vector<real_prec>kvector_data_GL(GAUSS_LEGANDERE_NODES_INTEGRATION,0);
  vector<real_prec>weight_GL(GAUSS_LEGANDERE_NODES_INTEGRATION,0);
  gsl_get_GL_weights(this->params._d_kmin(),this->params._d_kmax(), kvector_data_GL, weight_GL);
  So.DONE();
  real_prec pre_factor=2.0*pow(this->alpha,2)/pow(this->normal_power,2);
  ULONG Npairs=this->N_random*(this->N_random-1)/2;
  // Allocate structure for pre-computed pair_info:
  struct pair_info{
    ULONG ipair; // index of the i tracer in the pair
    ULONG jpair; // index of the j tracer in the pair
    real_prec weight_pair;  // peoduct of the weights w_i * w_j
    real_prec Rminus;   // rmin, modulus of the separation vector between tracers
    real_prec DeltaR;   // (r²_i - r²_j)/(rmax*rmin)
    vector<real_prec> Legpol; // Legandre_{l}(DeltaR)
  };
  vector<pair_info> pairs(Npairs);
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tAllocating pair indices");
#endif
  ULONG counter=0;
  for (ULONG i=0; i<this->N_random;++i)
    for (ULONG j=i+1; j<this->N_random;++j)
      {
        pairs[counter].ipair=i;
        pairs[counter].jpair=j;
        counter++;
      }
  So.DONE();
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tAllocating weights for pairs");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG ip=0; ip<counter;++ip)
  {
      ULONG i=pairs[ip].ipair;
      ULONG j=pairs[ip].jpair;
      real_prec rp=sqrt(pow(this->random_cat.Halo[i].coord1+this->random_cat.Halo[j].coord1,2)+pow(this->random_cat.Halo[i].coord2+this->random_cat.Halo[j].coord2,2)+pow(this->random_cat.Halo[i].coord3+this->random_cat.Halo[j].coord3,2));
   	  real_prec rm=sqrt(pow(this->random_cat.Halo[i].coord1-this->random_cat.Halo[j].coord1,2)+pow(this->random_cat.Halo[i].coord2-this->random_cat.Halo[j].coord2,2)+pow(this->random_cat.Halo[i].coord3-this->random_cat.Halo[j].coord3,2));
	    real_prec deltaR=(distances_cat[i]-distances_cat[j])/(rp*rm);
      pairs[ip].Rminus=rm;
      if(abs(deltaR)<=1)
        pairs[ip].DeltaR=deltaR;
      else
        pairs[ip].DeltaR=0;
      pairs[ip].weight_pair=this->random_cat.Halo[i].weight1*this->random_cat.Halo[j].weight1;
    }
  So.DONE();
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tComputing Legendre Pols");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG ip=0; ip<Npairs;++ip)                   // Loop over the set of random pairs
     {
       pairs[ip].Legpol.resize(MAXIMUM_MULTIPOLE_WINDOW_MATRIX,0);
       for (int lp=0;lp<MAXIMUM_MULTIPOLE_WINDOW_MATRIX;++lp)              //Loop over multipole lp
          pairs[ip].Legpol[lp]=gsl_sf_legendre_Pl(lp,pairs[ip].DeltaR);
     }
   So.DONE();
#ifdef _FULL_VERBOSE_
  this->So.message_screen("\tCounted ", counter, "pairs");
#endif
  string wfile=this->params._Output_directory()+"window_matrix";
  //this->window_matrix.resize(kvector_data_GL.size()*kvector_data.size()*MAXIMUM_MULTIPOLE_POWER,0);
  for (int l=0;l<MAXIMUM_MULTIPOLE_POWER;l+=2)              // Loop over multipole l
    {
      this->So.message_screen("\tComputing for multipole l =", l);
      for (ULONG ik=1;ik<this->kvector_data.size();++ik)    // Loop over wavenumbers of measured power. Starts from the first mode
      {
       for (int lp=0;lp<MAXIMUM_MULTIPOLE_WINDOW_MATRIX;lp+=2)              //Loop over multipole lp
	    {
          real_prec ii=lp==0? 1: -1;
          this->So.message_screen("\t \t Computing for Fourier kmode =", this->kvector_data[ik]);
	      ofstream wind_out;
          string  wfile_kmode=wfile+"_kbin"+to_string(ik)+"_l"+to_string(l)+"_lp"+to_string(lp)+".txt";
          So.message_screen("\t\tWriting mixing matrix in file ", wfile_kmode);
	      wind_out.open(wfile_kmode.c_str());
          for (ULONG jk=0;jk<kvector_data_GL.size() ;++jk)     //Loop over wanvenumbers GS integration
          {
            real_prec factor_n=ii*(2*l+1)*pre_factor*kvector_data_GL[jk]*kvector_data_GL[jk]*weight_GL[jk];
            real_prec window=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:window)
#endif
           for (ULONG ip=0; ip<Npairs;++ip)                   // Loop over the set of random pairs
             {
              real_prec bess_l= gsl_sf_bessel_jl(l, this->kvector_data[ik]*pairs[ip].Rminus); //j_l(kdeltaR)
		      real_prec leg_l = pairs[ip].Legpol[l];  // Leg_l()
              real_prec bess_lp=gsl_sf_bessel_jl(lp, kvector_data_GL[jk]*pairs[ip].Rminus);  //j_lp(kdeltaR)
		      real_prec leg_lp= pairs[ip].Legpol[lp]; //Leg_lp()
              window += factor_n*pairs[ip].weight_pair*bess_l*bess_lp*leg_l*leg_lp;
           }
           wind_out<<kvector_data_GL[jk]<<"  "<<window<<endl;
		  //    this->window_matrix[index_3d(l,nk,nkp,kvector_data.size(),kvector_data_GL.size())]=factor_n*window; //allocate in a container
          }
	      wind_out.close();
	    }
	}
    }
 this->So.DONE();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_CROSS_CORRELATION_CONF_SPACE_
 void PowerSpectrumF::get_cross_correlation_config_space(vector<real_prec>& X,vector<real_prec>& Y,vector<real_prec>& corr){
   this->So.enter(__PRETTY_FUNCTION__);
   this->params.set_measure_cross(true);
   this->params.set_input_type("density_grid");
   this->params.set_input_type_two("density_grid");
   this->params.set_SN_correction(false);
   this->fftw_functions.set_params(this->params);
   s_parameters_box s_p_box;
   real_prec ngal_new=0;
   fftw_functions.data_g.clear();
   fftw_functions.data_g.shrink_to_fit();
   fftw_functions.data_g.resize(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     ngal_new+=static_cast<real_prec>(Y[i]);
   real_prec nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     fftw_functions.data_g[i]=(static_cast<real_prec>(Y[i])/static_cast<real_prec>(nmean))-1.;
   this->params.set_ngal_delta(ngal_new);
   if(true==this->params._SN_correction())
     fftw_functions.shot_noise=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
   real_prec factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
   fftw_functions.set_normal_power(pow(factor,-2));
   fftw_functions.data_gp.clear();
   fftw_functions.data_gp.shrink_to_fit();
   fftw_functions.data_gp.resize(this->params._NGRID(),0);
   ngal_new=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     fftw_functions.data_gp[i]=static_cast<real_prec>(X[i]);
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_new)
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     ngal_new+=static_cast<real_prec>(X[i]);
   nmean=static_cast<real_prec>(ngal_new)/static_cast<real_prec>(this->params._NGRID());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->params._NGRID();++i)
     fftw_functions.data_gp[i]=(static_cast<real_prec>(X[i])/static_cast<real_prec>(nmean))-1.;
   fftw_functions.shot_noise2=0;
   factor=pow(static_cast<real_prec>(this->params._Lbox()),1.5)/static_cast<real_prec>(fftw_functions.data_g.size());
   fftw_functions.set_normal_power_two(pow(factor,-2));
   bool dm=true;
   if(false==dm && true==this->params._SN_correction())
     fftw_functions.shot_noise2=static_cast<real_prec>(pow(this->params._Lbox(),3))/static_cast<real_prec>(ngal_new);
   fftw_functions.resize_fftw_vectors();
   this->pk0.clear();
   this->pk0.shrink_to_fit();
   this->pk0.resize(this->params._d_Nnp_data(),0); //Monopole
   this->modes_g.clear();
   this->modes_g.resize(this->params._d_Nnp_data(),0); //Monopole
   fftw_functions.cross_power_spectrum_fkp(this->pk0,this->modes_g,corr);
 }
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // Assignment of individual bias. Halo catalog is in the vector tracer_cat. Dark matter field in the dm_field
 void PowerSpectrumF::object_by_object_bias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field){
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
   this->So.message_screen("Getting object-to-object bias");
#endif
     kvector_data.clear();
#if defined _TNG_ || defined _TNG_GAL_ || defined _UNITSIM_
     kvector_data.shrink_to_fit();
     for(int i=0;i<this->params._d_Nnp_data();i++)
       kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
#endif
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   // Given kmax, determine the maximum number of bins requested
   ULONG Kmin_bin=static_cast<ULONG>(floor((this->params._kmin_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));
   ULONG Kmax_bin=static_cast<ULONG>(floor((this->params._kmax_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   ULONG new_Nft=2*Kmax_bin;  // Define the new Nft=2*Kmax_bin
   ULONG new_ngrid_h =new_Nft*new_Nft*(new_Nft/2+1);
#ifdef DOUBLE_PREC
   complex_prec * Delta_dm = (complex_prec *)fftw_malloc(2*new_ngrid_h*sizeof(real_prec));
#else
//    complex_prec * Delta_dm =(complex_prec *)fftwf_malloc(2*new_ngrid_h*sizeof(real_prec));
    complex_prec * Delta_dm =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
    real_prec mean_field=get_mean(dm_field);
    So.message_screen("\tMean of field:", mean_field);
    vector<real_prec> over_dm_field(dm_field.size(),0);
    get_overdens(dm_field,mean_field,over_dm_field);
    So.message_screen("\tFourier transforming");
    do_fftw_r2c(this->params._Nft(),over_dm_field,Delta_dm);
    over_dm_field.clear(); over_dm_field.shrink_to_fit();
   vector<real_prec> kcoords(new_Nft,0);// Build k-coordinates
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<kcoords.size() ;++i)
    kcoords[i]=(i<=new_Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(new_Nft-i));
#ifdef _FULL_VERBOSE_
   this->So.message_screen("\tNew Nft ",new_Nft);
   this->So.message_screen("\tComputing from k =", (Kmin_bin+0.5)*this->params._d_DeltaK_data());
   this->So.message_screen("\t            to k =", this->params._kmax_tracer_bias());
#endif
   real_prec use_imag=0;
  // This is the normalization of the halo power spectrum, which being computed object-by object, is just the volume
   // FOr the full sample it would be the volume/Ntracer, = nbar. But Ntracer =1 for each object!
   real_prec dk_x=this->params._d_deltak_x();
   real_prec dk_y=this->params._d_deltak_y();
   real_prec dk_z=this->params._d_deltak_z();
   this->So.message_screen("\tGetting power dark matter");
  // This is done once for all tracers
   vector<real_prec>power_dmat(Kmax_bin,0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
   for(ULONG i=Kmin_bin; i< new_Nft/2;++i)
     for(ULONG j=Kmin_bin; j< new_Nft/2;++j)
         for(ULONG k=Kmin_bin; k< new_Nft/2+1;++k)
         {
           ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
           real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
           ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
           if(lp!=0)
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
           if(j>0  && k>0)
             {
               lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
            }
           if(i>0  && (j>0 || k>0))
             {
               lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_object_b
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
             }
             if(i>0  && j>0  && k>0)
             {
               lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
             }
         }
    So.DONE();
    real_prec conversion_factor=(1.+this->params._redshift())/(this->cosmology.Hubble_function(this->params._redshift()));
#ifndef _TNG_GAL
    real_prec kmax_b=0.08;
    real_prec kmax_c=0.04;
#elif defined _TNG_GAL_
     real_prec kmax_b=0.25;
     real_prec kmax_c=0.2;
#endif
     real_prec lss_bias_halo=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:lss_bias_halo)
#endif
   for(ULONG itr=0;itr<tracer_cat.size();++itr)
     {
        real_prec xtracer=tracer_cat[itr].coord1;
        real_prec ytracer=tracer_cat[itr].coord2;
        real_prec ztracer=tracer_cat[itr].coord3;
#ifndef _TNG_GAL_
        real_prec vx=tracer_cat[itr].vel1*conversion_factor;
#endif
        vector<real_prec>power_cross(Kmax_bin,0);
#ifndef _TNG_GAL_
       vector<real_prec>power_cross_s(Kmax_bin,0);
       vector<real_prec>Gamma_num(Kmax_bin,0);
#endif
       //*******   The cross power spectrum has real and imaginary parts. C=A+iB.If the fieidsl are not properly correlated, the imaginary parts will be non-zero
       //*******   I.e, for perfectly correlated fields, B== and the cross power is C=A.
       //*******   Now, for bias. we compare the bias b_1 fro the one measured with power spectrum b²=Phh/Pmm.
       //*******   while here, with the cross, we compute b = Phm/Pmm. This is chosen such that Phm is decomposed in delta_dm * exp(-ikr)
       // Loop over the Fourier box up to the maximum k (bin) used to get bias
       for(ULONG i=Kmin_bin; i< new_Nft/2;++i)
         for(ULONG j=Kmin_bin; j< new_Nft/2;++j)
             for(ULONG k=Kmin_bin; k< new_Nft/2+1;++k)
             {
               ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               real_prec k_dot_r=0;
               real_prec k_dot_s=0;
               real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
               ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
               /****************************************************************/
               if(lp!=0)
               {
                 k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;  // vec k dot vec r
#ifndef _TNG_GAL_
                 k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                 if(kbin<Kmax_bin){
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
#ifndef _TNG_GAL_
                     power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                     Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                 }
               }
               /****************************************************************/
               if(j>0  && k>0)
                 {
                   lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin){
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
#ifndef _TNG_GAL_
                     power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                     Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                   }
                   }
               if(i>0  && (j>0 || k>0))
                 {
                   lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin){
                       power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
#ifndef _TNG_GAL_
                       power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                       Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                   }
               }
                 if(i>0  && j>0  && k>0)
                 {
                   lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer+dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin)
                   {
                       power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
#ifndef _TNG_GAL_
                       power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                       Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                  }
               }
          }
        real_prec power_hm=0;
        real_prec power_hm2=0;
        real_prec power_hm3=0;
        real_prec p_dm=0;
        real_prec p_dm2=0;
        real_prec p_dm3=0;
        real_prec power_hm_s=0;
        real_prec gama=0;
        for(ULONG i=Kmin_bin; i< power_cross.size();++i)
          {
            power_hm+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
#ifndef _TNG_GAL_
            power_hm_s+=power_cross_s[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
#endif
            p_dm+=power_dmat[i];
#ifndef _TNG_GAL_
            gama+=Gamma_num[i];
#endif
            if(this->kvector_data[i]<kmax_b)
            {
                power_hm2+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
                p_dm2+=power_dmat[i];
             }
            if(this->kvector_data[i]<kmax_c)
            {
                power_hm3+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
                p_dm3+=power_dmat[i];
             }
        }
        real_prec hb=(power_hm/p_dm)*this->params._NGRID();
        real_prec hb2=(power_hm2/p_dm2)*this->params._NGRID();
        real_prec hb3=(power_hm3/p_dm3)*this->params._NGRID();
#ifndef _TNG_GAL_
        real_prec hbs=power_hm_s*this->params._NGRID();
        gama= (gama/p_dm)*this->params._NGRID();
#endif
        tracer_cat[itr].bias=hb;
        tracer_cat[itr].bias2=hb2;
        tracer_cat[itr].bias3=hb3;
#ifndef _TNG_GAL_
        tracer_cat[itr].bias_rs=hbs;
        tracer_cat[itr].rs_factor=gama;
#endif
        lss_bias_halo+=hb;
   }
   lss_bias_halo/=static_cast<real_prec>(tracer_cat.size());
   So.message_screen("\tMean large-scale bias from individual bias =", lss_bias_halo);
#ifdef _TNG_GAL_
   vector<real_prec>baux(tracer_cat.size(),0);
   vector<real_prec>xaux(tracer_cat.size(),0);
   vector<real_prec>yaux(tracer_cat.size(),0);
   vector<real_prec>zaux(tracer_cat.size(),0);

   for(ULONG i=0;i<xaux.size(); ++i)
    {
        baux[i]=tracer_cat[i].bias;
        xaux[i]=tracer_cat[i].coord1;
        yaux[i]=tracer_cat[i].coord2;
        zaux[i]=tracer_cat[i].coord3;
    }
    string file_bias=this->params._Output_directory()+"Bias_gal1";
    this->File.write_array(file_bias, baux);

    vector<real_prec>bfaux(this->params._NGRID(),0);

    string file_bias_f=this->params._Output_directory()+"Bias_gal1_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);


/*
    file_bias=this->params._Output_directory()+"Bias_gal2";
    for(ULONG i=0;i<baux.size(); ++i)
       baux[i]=tracer_cat[i].bias2;
     this->File.write_array(file_bias, baux);
    file_bias_f=this->params._Output_directory()+"Bias_gal2_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);



    file_bias=this->params._Output_directory()+"Bias_gal3";
    for(ULONG i=0;i<baux.size(); ++i)
       baux[i]=tracer_cat[i].bias3;
     this->File.write_array(file_bias, baux);
    file_bias_f=this->params._Output_directory()+"Bias_gal3_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);
*/

#endif
   So.DONE();
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::object_by_object_rbias(){
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
   this->So.message_screen("Getting object-to-object relative bias");
#endif
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   // Given kmax, determine the maximum number of bins requested
   ULONG initial_mode=static_cast<ULONG>(floor((this->params._kmin_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   ULONG Kmax_bin=static_cast<ULONG>(floor((this->params._kmax_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   // Define the new Nft=2*Kmax_bin
   ULONG new_Nft=2*Kmax_bin;
   ULONG new_ngrid_h =new_Nft*new_Nft*(new_Nft/2+1);
#ifdef DOUBLE_PREC
   complex_prec * Delta_tr = (complex_prec *)fftw_malloc(2*new_ngrid_h*sizeof(real_prec));
#else
    complex_prec * Delta_tr =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
    complex_prec * Delta_W =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
   // ---------------------------------------------------
   // Get the Fourier transform of the dm filed alreay filtered
    vector<real_prec> over_field(this->params._NGRID(),0);
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
     for(ULONG i=0;i<this->params._NGRID();++i)
        over_field[i]=this->tracer_cat.field_external[i]-this->fftw_functions._alpha()*this->random_cat.field_external[i];
    So.message_screen("\tFourier transforming in ", __PRETTY_FUNCTION__);
    do_fftw_r2c(this->params._Nft(),over_field,Delta_tr);
#ifdef _USE_OMP_
#pragma omp for nowait
#endif
     for(ULONG i=0;i<this->params._NGRID();++i)
        over_field[i]=this->fftw_functions._alpha()*this->random_cat.field_external[i];
     do_fftw_r2c(this->params._Nft(),over_field,Delta_W);// Window funtion in Fourier
     over_field.clear(); over_field.shrink_to_fit();
    // ---------------------------------------------------
   // Build k-coordinates
    real_prec dk_x=this->params._d_deltak_x();
    real_prec dk_y=this->params._d_deltak_y();
    real_prec dk_z=this->params._d_deltak_z();
   vector<real_prec> kcoords(new_Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<kcoords.size() ;++i)
    kcoords[i]=(i<=new_Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(new_Nft-i));
   // ---------------------------------------------------
#ifdef _FULL_VERBOSE_
   this->So.message_screen("\tNew Nft ",new_Nft);
   this->So.message_screen("\tComputing from k =", (initial_mode+0.5)*this->params._d_DeltaK_data());
   this->So.message_screen("\t            to k =", this->params._kmax_tracer_bias());
#endif
   // ---------------------------------------------------
  // This is the normalization of the halo power spectrum, which being computed object-by object, is just the volume
   // FOr the full sample it would be the volume/Ntracer, = nbar. But Ntracer =1 for each object!
   this->So.message_screen("\tGetting power dark matter:");
  // This is done once for all tracers
   vector<real_prec>power_tracer(Kmax_bin,0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
   for(ULONG i=initial_mode; i< new_Nft/2;++i)
     for(ULONG j=initial_mode; j< new_Nft/2;++j)
         for(ULONG k=initial_mode; k< new_Nft/2+1;++k)
         {
           ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
           real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
           ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
           if(lp!=0)
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
           /****************************************************************/
           if(j>0  && k>0)
             {
               lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                 if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
            }
           /****************************************************************/
           if(i>0  && (j>0 || k>0))
             {
               lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
             }
           /****************************************************************/
             if(i>0  && j>0  && k>0)
             {
               lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
             }
       }

   So.DONE();
    this->So.message_screen("\tAssigining relative bias to downsampled catalog:");
   ULONG Ntr=0;
   real_prec w_g_reduced=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Ntr, w_g_reduced)
#endif
   for(ULONG itr=0;itr<tracer_cat.Halo.size();++itr)
    if(this->tracer_cat.Halo[itr].observed==1)
    {
        real_prec we_fkp=1.0;
        if(true==this->params._FKP_weight())
          we_fkp=1.0/(1.0+this->params._Pest()*this->tracer_cat.Halo[itr].mean_density);
        w_g_reduced+=we_fkp;
        Ntr++;
    }
   real_prec lss_rbias_halo=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:lss_rbias_halo)
#endif
   for(ULONG itr=0;itr<tracer_cat.Halo.size();++itr)
   {
     if(tracer_cat.Halo[itr].observed==1)
       {
        real_prec xtracer=this->tracer_cat.Halo[itr].coord1;
        real_prec ytracer=this->tracer_cat.Halo[itr].coord2;
        real_prec ztracer=this->tracer_cat.Halo[itr].coord3;
        real_prec we_fkp=1.0;
        if(true==this->params._FKP_weight())
          we_fkp=1.0/(1.0+this->params._Pest()*this->tracer_cat.Halo[itr].mean_density);
        vector<real_prec>power_cross_tr(Kmax_bin,0);
       //*******   The cross power spectrum has real and imaginary parts. C=A+iB.If the fieidsl are not properly correlated, the imaginary parts will be non-zero
       //*******   I.e, for perfectly correlated fields, B== and the cross pwoer is C=A.
       //*******   Now, for bias. we compare the bias b_1 fro the one measured with power spectrum b²=Phh/Pmm.
       //*******   while here, with the cross, we compute b = Phm/Pmm. This is chosen as Phm is decomposed in delta_dm * exp(-ikr)
      // Loop over the Fourier box up to the maximum k (bin) used to get bias
       for(ULONG i=initial_mode; i< new_Nft/2;++i)
         for(ULONG j=initial_mode; j< new_Nft/2;++j)
             for(ULONG k=initial_mode; k< new_Nft/2+1;++k)
             {
               ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               real_prec k_dot_r=0;
               real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
               ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
               /****************************************************************/
               if(lp!=0)
               {
                 k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
                 if(kbin<Kmax_bin){
                     real_prec term_a=we_fkp*(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
                     real_prec term_b=(cos(k_dot_r)*Delta_W[lp][REAL]+sin(k_dot_r)*Delta_W[lp][IMAG]);
                     power_cross_tr[kbin]+=term_a-term_b;
                 }
                }
               if(j>0  && k>0)
                 {
                   lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin){
                       real_prec term_a=we_fkp*(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
                       real_prec term_b=(cos(k_dot_r)*Delta_W[lp][REAL]+sin(k_dot_r)*Delta_W[lp][IMAG]);
                       power_cross_tr[kbin]+=term_a-term_b;
                   }
                 }
               if(i>0  && (j>0 || k>0))
                 {
                   lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin){
                       real_prec term_a=we_fkp*(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
                       real_prec term_b=(cos(k_dot_r)*Delta_W[lp][REAL]+sin(k_dot_r)*Delta_W[lp][IMAG]);
                       power_cross_tr[kbin]+=term_a-term_b;
                   }
               }
                 if(i>0  && j>0  && k>0)
                 {
                   lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer+dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin)
                   {
                       real_prec term_a=we_fkp*(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
                       real_prec term_b=(cos(k_dot_r)*Delta_W[lp][REAL]+sin(k_dot_r)*Delta_W[lp][IMAG]);
                       power_cross_tr[kbin]+=term_a-term_b;
                   }
                  }
             }
        real_prec power_hmt=0;
        real_prec p_t=0;
        for(ULONG i=initial_mode; i< power_cross_tr.size();++i)
          {
            power_hmt+=power_cross_tr[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
            p_t+=power_tracer[i];
        }
        // Here the fist vol comes from the discrete FT. The term vol/Ngrid is the normzalization of the term 1/dm, i.e the DM power is normalized with Ngrid/vol
//         real_prec hb=(vol)*(power_hm/p_dm)/(vol/this->params._NGRID()); THIS WAS THE ORIGINAL SHAPE, below I just rearrange terms
//        real_prec hbtr=(power_hmt/p_t)*(this->params._NGRID())*(static_cast<real_prec>(this->tracer_cat.Halo.size())/static_cast<real_prec>(Ntr));
        real_prec hbtr=(power_hmt/p_t)*(this->params._NGRID())*(this->fftw_functions._w_g()/w_g_reduced);
        this->tracer_cat.Halo[itr].relative_bias=hbtr;
        lss_rbias_halo+=hbtr;
       } // closes if
     }// Closes loop over randoms
   lss_rbias_halo/=static_cast<real_prec>(Ntr);
   So.message_screen("\tMean large-scale relative bias from individual bias =", lss_rbias_halo);
   this->tracer_cat.write_catalog(this->params._Output_directory()+"reduced_gal_cat.txt");
   So.DONE();
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void PowerSpectrumF::object_by_object_bias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field, vector<real_prec>&tracer_field){
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
   this->So.message_screen("Getting object-to-object bias");
#endif
#if defined _TNG_ || defined _TNG_GAL_
     kvector_data.clear();
     kvector_data.shrink_to_fit();
     for(int i=0;i<this->params._d_Nnp_data();i++)
       kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
#endif
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif

   // Given kmax, determine the maximum number of bins requested
   ULONG initial_mode=static_cast<ULONG>(floor((this->params._kmin_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   ULONG Kmax_bin=static_cast<ULONG>(floor((this->params._kmax_tracer_bias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   // ---------------------------------------------------
   // Define the new Nft=2*Kmax_bin
   ULONG new_Nft=2*Kmax_bin;
   // ---------------------------------------------------
   ULONG new_ngrid_h =new_Nft*new_Nft*(new_Nft/2+1);
#ifdef DOUBLE_PREC
   complex_prec * Delta_dm = (complex_prec *)fftw_malloc(2*new_ngrid_h*sizeof(real_prec));
#else
    complex_prec * Delta_dm =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
#ifdef DOUBLE_PREC
   complex_prec * Delta_tr = (complex_prec *)fftw_malloc(2*new_ngrid_h*sizeof(real_prec));
#else
    complex_prec * Delta_tr =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
    real_prec mean_field=get_mean(dm_field);
    So.message_screen("\tMean of field:", mean_field);
    vector<real_prec> over_dm_field(dm_field.size(),0);
    get_overdens(dm_field,mean_field,over_dm_field);
    So.message_screen("\tFourier transforming");
    do_fftw_r2c(this->params._Nft(),over_dm_field,Delta_dm);
    over_dm_field.clear(); over_dm_field.shrink_to_fit();
    mean_field=get_mean(tracer_field);
    So.message_screen("\tMean of field:", mean_field);
    get_overdens(tracer_field,mean_field,tracer_field);
    So.message_screen("\tFourier transforming");
    do_fftw_r2c(this->params._Nft(),tracer_field,Delta_tr);
   // Build k-coordinates
   vector<real_prec> kcoords(new_Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<kcoords.size() ;++i)
    kcoords[i]=(i<=new_Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(new_Nft-i));
#ifdef _FULL_VERBOSE_
   this->So.message_screen("\tNew Nft ",new_Nft);
   this->So.message_screen("\tComputing from k =", (initial_mode+0.5)*this->params._d_DeltaK_data());
   this->So.message_screen("\t            to k =", this->params._kmax_tracer_bias());
#endif
   real_prec use_imag=0;
   real_prec vol=pow(this->params._Lbox(),3);
   real_prec normal_dm= sqrt(vol)/static_cast<real_prec>(this->params._NGRID());
   real_prec normal_h =static_cast<real_prec>(1./vol);
   real_prec normal_cross= 1; //normal_dm * normal_h;
   real_prec dk_x=this->params._d_deltak_x();
   real_prec dk_y=this->params._d_deltak_y();
   real_prec dk_z=this->params._d_deltak_z();
   this->So.message_screen("\tGetting power dark matter:");
  // This is done once for all tracers
   vector<real_prec>power_dmat(Kmax_bin,0);
   vector<real_prec>power_tracer(Kmax_bin,0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
   for(ULONG i=initial_mode; i< new_Nft/2;++i)
     for(ULONG j=initial_mode; j< new_Nft/2;++j)
         for(ULONG k=initial_mode; k< new_Nft/2+1;++k)
         {
           ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
           real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
           ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
           if(lp!=0)
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);

             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
          if(j>0  && k>0)
             {
               lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
            }
           if(i>0  && (j>0 || k>0))
             {
               lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_object_b
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
             }
             if(i>0  && j>0  && k>0)
             {
               lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
             if(kbin<Kmax_bin)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_tracer[kbin]+=(Delta_tr[lp][REAL]*Delta_tr[lp][REAL]+Delta_tr[lp][IMAG]*Delta_tr[lp][IMAG]);
             }
         }
    So.DONE();
    real_prec conversion_factor=(1.+this->params._redshift())/(this->cosmology.Hubble_function(this->params._redshift()));
#ifndef _TNG_GAL_
    real_prec kmax_b=0.06;
    real_prec kmax_c=0.04;
#elif defined _TNG_GAL_
     real_prec kmax_b=0.25;
     real_prec kmax_c=0.2;
#endif
     real_prec lss_bias_halo=0;
     real_prec lss_rbias_halo=0;
     this->So.message_screen("\tAssigining bias:");
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:lss_bias_halo)
#endif
   for(ULONG itr=0;itr<tracer_cat.size();++itr)
     {
        real_prec xtracer=tracer_cat[itr].coord1;
        real_prec ytracer=tracer_cat[itr].coord2;
        real_prec ztracer=tracer_cat[itr].coord3;
#ifndef _TNG_GAL_
        real_prec vx=tracer_cat[itr].vel1*conversion_factor;
#endif
        vector<real_prec>power_cross(Kmax_bin,0);
        vector<real_prec>power_cross_tr(Kmax_bin,0);
#ifndef _TNG_GAL_
       vector<real_prec>power_cross_s(Kmax_bin,0);
       vector<real_prec>Gamma_num(Kmax_bin,0);
#endif
        // Loop over the Fourier box up to the maximum k (bin) used to get bias
       for(ULONG i=initial_mode; i< new_Nft/2;++i)
         for(ULONG j=initial_mode; j< new_Nft/2;++j)
             for(ULONG k=initial_mode; k< new_Nft/2+1;++k)
             {
               ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               real_prec k_dot_r=0;
               real_prec k_dot_s=0;
               real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
               ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
               /****************************************************************/
               if(lp!=0)
               {
                 k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                 k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                 if(kbin<Kmax_bin){
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                     power_cross_tr[kbin]+=(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
#ifndef _TNG_GAL_
                     power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                     Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                 }

                 }
               if(j>0  && k>0)
                 {
                   lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin){
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                     power_cross_tr[kbin]+=(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
#ifndef _TNG_GAL_
                     power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                     Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                   }
                   }
               if(i>0  && (j>0 || k>0))
                 {
                   lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin){
                       power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                      power_cross_tr[kbin]+=(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
#ifndef _TNG_GAL_
                       power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                       Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                   }
               }
                 if(i>0  && j>0  && k>0)
                 {
                   lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer+dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
#ifndef _TNG_GAL_
                   k_dot_s=dk_x*kcoords[i]*(xtracer+vx) + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
#endif
                   if(kbin<Kmax_bin)
                   {
                       power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                       power_cross_tr[kbin]+=(cos(k_dot_r)*Delta_tr[lp][REAL]-sin(k_dot_r)*Delta_tr[lp][IMAG]);
#ifndef _TNG_GAL_
                       power_cross_s[kbin]+=(cos(k_dot_s)*Delta_dm[lp][REAL]-sin(k_dot_s)*Delta_dm[lp][IMAG]);
                       Gamma_num[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG])*pow(cos(dk_x*kcoords[i]/kv),2);
#endif
                 }
               }
           }
        real_prec power_hm=0;
        real_prec power_hmt=0;
        real_prec power_hm2=0;
        real_prec power_hm3=0;
        real_prec p_dm=0;
        real_prec p_t=0;
        real_prec p_dm2=0;
        real_prec p_dm3=0;
        real_prec power_hm_s=0;
        real_prec gama=0;
        for(ULONG i=initial_mode; i< power_cross.size();++i)
          {
            power_hm+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
            power_hmt+=power_cross_tr[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
#ifndef _TNG_GAL_
            power_hm_s+=power_cross_s[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
#endif
            p_dm+=power_dmat[i];
            p_t+=power_tracer[i];
#ifndef _TNG_GAL_
            gama+=Gamma_num[i];
#endif
            if(this->kvector_data[i]<kmax_b)
            {
                power_hm2+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
                p_dm2+=power_dmat[i];
             }
            if(this->kvector_data[i]<kmax_c)
            {
                power_hm3+=power_cross[i]; // *** here it must be nmodes*(<power>_av)= nmodes*(power/nmodes)=power. That's why I do not need nmodes
                p_dm3+=power_dmat[i];
             }

        }
        real_prec hb=(power_hm/p_dm)*(this->params._NGRID());
        real_prec hbtr=(power_hmt/p_t)*(this->params._NGRID());
        real_prec hb2=(power_hm2/p_dm2)*(this->params._NGRID());
        real_prec hb3=(power_hm3/p_dm3)*(this->params._NGRID());
#ifndef _TNG_GAL_
        real_prec hbs=(power_hm_s/p_dm)*(this->params._NGRID());
        gama= (vol)*(gama/p_dm)/(vol/this->params._NGRID());
#endif
        tracer_cat[itr].bias=hb;
        tracer_cat[itr].relative_bias=hbtr;
        tracer_cat[itr].bias2=hb2;
        tracer_cat[itr].bias3=hb3;
#ifndef _TNG_GAL_
        tracer_cat[itr].bias_rs=hbs;
        tracer_cat[itr].rs_factor=gama;
#endif
        lss_bias_halo+=hb;
        lss_rbias_halo+=hbtr;
   }
   lss_bias_halo/=static_cast<real_prec>(tracer_cat.size());
   lss_rbias_halo/=static_cast<real_prec>(tracer_cat.size());
   So.message_screen("\tMean large-scale bias from individual bias =", lss_bias_halo);
   So.message_screen("\tMean large-scale relative bias from individual bias =", lss_rbias_halo);
#ifdef _TNG_GAL_
   vector<real_prec>baux(tracer_cat.size(),0);
   vector<real_prec>xaux(tracer_cat.size(),0);
   vector<real_prec>yaux(tracer_cat.size(),0);
   vector<real_prec>zaux(tracer_cat.size(),0);
    for(ULONG i=0;i<xaux.size(); ++i)
    {
        baux[i]=tracer_cat[i].bias;
        xaux[i]=tracer_cat[i].coord1;
        yaux[i]=tracer_cat[i].coord2;
        zaux[i]=tracer_cat[i].coord3;
    }
    string file_bias=this->params._Output_directory()+"Bias_gal1";
    this->File.write_array(file_bias, baux);
    vector<real_prec>bfaux(this->params._NGRID(),0);
    string file_bias_f=this->params._Output_directory()+"Bias_gal1_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);
    file_bias=this->params._Output_directory()+"Bias_gal2";
    for(ULONG i=0;i<baux.size(); ++i)
       baux[i]=tracer_cat[i].bias2;
     this->File.write_array(file_bias, baux);
    file_bias_f=this->params._Output_directory()+"Bias_gal2_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);
    file_bias=this->params._Output_directory()+"Bias_gal3";
    for(ULONG i=0;i<baux.size(); ++i)
       baux[i]=tracer_cat[i].bias3;
     this->File.write_array(file_bias, baux);
    file_bias_f=this->params._Output_directory()+"Bias_gal3_field";
    getDensity_CIC(this->params._Nft(),this->params._Nft(),this->params._Nft(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(),this->params._d_delta_x(),this->params._d_delta_x(),this->params._d_delta_x(),0,0,0,xaux,yaux,zaux,baux,bfaux,true);
    this->File.write_array(file_bias_f, bfaux);

#endif
   So.DONE();
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::object_by_object_qbias(vector<s_Halo>& tracer_cat, vector<real_prec>& dm_field){
   this->So.enter(__PRETTY_FUNCTION__);
#ifdef _FULL_VERBOSE_
   this->So.message_screen("Getting object-to-object quadratic bias");
#endif
    real_prec kmin_int=0.001;
    real_prec kmax_int=1e3;
    PowerSpectrum PS(this->params.s_cosmo_pars,kmin_int, kmax_int, 1000, 30);
    this->params.s_cosmo_pars.pk_normalization=PS.normalization();
    PS.set_cosmo_pars(this->params.s_cosmo_pars); // update
    kvector_data.clear();
    kvector_data.shrink_to_fit();
    for(int i=0;i<this->params._d_Nnp_data();i++)
      kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
   ULONG Kmax_bin=static_cast<ULONG>(floor((this->params._kmax_tracer_qbias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   ULONG Kmin_bin=static_cast<ULONG>(floor((this->params._kmin_tracer_qbias()-this->params._d_kmin())/this->params._d_DeltaK_data()));//get the bin to go only up to Kmax in the loops
   ULONG new_Nft=2*Kmax_bin;
   ULONG new_ngrid_h =new_Nft*new_Nft*(new_Nft/2+1);
#ifdef DOUBLE_PREC
   complex_prec * Delta_dm = (complex_prec *)fftw_malloc(2*new_ngrid_h*sizeof(real_prec));
#else
    complex_prec * Delta_dm =(complex_prec *)fftwf_malloc(2*this->params._NGRID_h()*sizeof(real_prec));
#endif
    real_prec mean_field=get_mean(dm_field);
    So.message_screen("\tMean of field:", mean_field);
    vector<real_prec> over_dm_field(dm_field.size(),0);
    get_overdens(dm_field,mean_field,over_dm_field);
    So.message_screen("\tFourier transforming");
    do_fftw_r2c(this->params._Nft(),over_dm_field,Delta_dm);
    over_dm_field.clear(); over_dm_field.shrink_to_fit();
   vector<real_prec> kcoords(new_Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<kcoords.size() ;++i)
    kcoords[i]=(i<=new_Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(new_Nft-i));
#ifdef _FULL_VERBOSE_
   this->So.message_screen("\tNew Nft ",new_Nft);
   this->So.message_screen("\tComputing from k =", this->params._kmin_tracer_qbias());
   this->So.message_screen("\t            to k =", this->params._kmax_tracer_qbias());
#endif
   real_prec vol=pow(this->params._Lbox(),3);
   real_prec normal_dm= sqrt(vol)/static_cast<real_prec>(this->params._NGRID());
   real_prec normal_h =static_cast<real_prec>(1./vol);
   real_prec normal_cross= 1; //normal_dm * normal_h;
   real_prec dk_x=this->params._d_deltak_x();
   real_prec dk_y=this->params._d_deltak_y();
   real_prec dk_z=this->params._d_deltak_z();
   this->So.message_screen("\tGetting power dark matter P11 and P22:");
   vector<real_prec>power_dmat(Kmax_bin,0); // container for the dm power
   vector<real_prec>Apower(Kmax_bin,0);    // container for the function A(k)
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
   for(ULONG i=Kmin_bin; i< new_Nft/2;++i)
     for(ULONG j=Kmin_bin; j< new_Nft/2;++j)
         for(ULONG k=Kmin_bin; k< new_Nft/2+1;++k)
         {
           ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
           real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
           ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
           if(lp!=0)
             if(kbin<Kmax_bin){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 Apower[kbin]+=PS.P1loop(kv);
             }
           if(j>0  && k>0)
             {
               lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin)
               {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 Apower[kbin]+=PS.P1loop(kv);
             }
            }
           if(i>0  && (j>0 || k>0))
             {
               lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                 Apower[kbin]+=PS.P1loop(kv);
             }
             }
             if(i>0  && j>0  && k>0)
             {
               lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
               if(kbin<Kmax_bin){
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                   power_dmat[kbin]+=(Delta_dm[lp][REAL]*Delta_dm[lp][REAL]+Delta_dm[lp][IMAG]*Delta_dm[lp][IMAG]);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                  Apower[kbin]+=PS.P1loop(kv);
             }
          }
       }
    So.DONE();
    this->So.message_screen("\tAssigining quadratic bias:");
    real_prec lss_bias_halo=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:lss_bias_halo)
#endif
   for(ULONG itr=0;itr<tracer_cat.size();++itr)
     {
        real_prec xtracer=tracer_cat[itr].coord1;
        real_prec ytracer=tracer_cat[itr].coord2;
        real_prec ztracer=tracer_cat[itr].coord3;
        vector<real_prec>power_cross(Kmax_bin,0);
       for(ULONG i=Kmin_bin; i< new_Nft/2;++i)
         for(ULONG j=Kmin_bin; j< new_Nft/2;++j)
             for(ULONG k=Kmin_bin; k< new_Nft/2+1;++k)
             {
               ULONG lp=index_3d(i,j,k,this->params._Nft(),this->params._Nft()/2+1);
               real_prec k_dot_r=0;
               real_prec  kv=sqrt(pow(dk_x*kcoords[i],2)+pow(dk_y*kcoords[j],2)+pow(dk_z*kcoords[k],2));
               ULONG kbin=static_cast<ULONG>(floor((kv-this->params._d_kmin())/this->params._d_DeltaK_data()));
               if(lp!=0)
               {
                 k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
                 if(kbin<Kmax_bin)
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                 }
               if(j>0  && k>0)
                 {
                   lp=index_3d(i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[i]*xtracer + dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin)
                     power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                   }
               if(i>0  && (j>0 || k>0))
                 {
                   lp=index_3d(this->params._Nft()-i,j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer + dk_y*kcoords[j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin)
                       power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
               }
                 if(i>0  && j>0  && k>0)
                 {
                   lp=index_3d(this->params._Nft()-i,this->params._Nft()-j,k,this->params._Nft(),this->params._Nft()/2+1);
                   k_dot_r=dk_x*kcoords[new_Nft-i]*xtracer+dk_y*kcoords[new_Nft-j]*ytracer + dk_z*kcoords[k]*ztracer;
                   if(kbin<Kmax_bin)
                    power_cross[kbin]+=(cos(k_dot_r)*Delta_dm[lp][REAL]-sin(k_dot_r)*Delta_dm[lp][IMAG]);
                }
             }
        real_prec power_hm=0;
        real_prec p_dm=0;
        for(ULONG i=Kmin_bin; i< power_cross.size();++i)
          {
              power_hm+=(power_cross[i]*this->params._NGRID() - tracer_cat[itr].bias*power_dmat[i]);
              p_dm+=Apower[i];
            }
        real_prec hb=(power_hm/p_dm)/this->params._NGRID();
        tracer_cat[itr].qbias=hb;
        lss_bias_halo+=hb;
   }
   lss_bias_halo/=static_cast<real_prec>(tracer_cat.size());
   So.DONE();
   So.message_screen("\tMean quadratic bias from individual bias =", lss_bias_halo);
 }
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //#define paired
// #define parallel_grf
#define white_noise
 void PowerSpectrumF::get_GaussianRandomField(){
#ifdef _USE_OMP_
   omp_set_num_threads(_NTHREADS_);
#endif
   ULONG ir=0;
   int ndel=1;
   real_prec Lside=this->params._Lbox();
   real_prec kmax=this->params._Kmax_FA();
   int Nft=static_cast<int>(this->params._Nft());
   int Nft_LR=Nft;
   ULONG Nft_HR=this->params._Nft_HR();
   ULONG NTT=Nft*Nft*(Nft/2+1);
   ULONG Ngrid_LR=Nft_LR*Nft_LR*Nft_LR;
   ULONG Ngrid_HR=Nft_HR*Nft_HR*Nft_HR;
   ULONG Ngrid=this->params._NGRID();
   if(true==this->params._use_low_pass_filter())
     {
       this->params.set_Nft(this->params._Nft_HR());// assign Nft_HR to Nft if we want to apply the low-pass filter
       Nft=this->params._Nft();
       Ngrid=this->params.d_NGRID();
       NTT=Nft*Nft*(Nft/2+1);
     }
   // ******************************************************************************
   // Computation of cosmological growth
   real_prec g1=this->cosmology.growth_factor(this->params._Initial_Redshift_DELTA());
   real_prec g2=this->cosmology.growth_factor(this->params._Initial_Redshift_TH_power_file());
   real_prec ic_factor=1.;
   if(this->params._Normalize_IC_to_initial_redshift()==true)
     ic_factor= pow(g1/g2,2);
#ifdef _ABACUS_
   real_prec abacus_factor=47.304805056;
   ic_factor=1./(pow(abacus_factor,2));
#endif
   So.message_screen("Growth at redshift of theoretical power", g1);
   So.message_screen("Growth at redshift of initial conditions", g2);
   So.message_screen("IC Factor", ic_factor);
   // ******************************************************************************
   real_prec norm=(Ngrid*Ngrid) / pow(Lside,3);  // Normalization of power spectrum. Multiply by this factor to un-normalize the theoretical power
   vector<real_prec>prop;
   string power_file=this->params._dir()+this->params._ic_power_file();
   ULONG nlines_p=this->File.read_file(power_file,prop,8);
   ULONG ncols=(static_cast<ULONG>(prop.size()/nlines_p));
   vector<double>kv(nlines_p,0);
   vector<double>pv(nlines_p,0);
   for(ULONG i=0;i<nlines_p;++i)
     {
       kv[i]=static_cast<double>(prop[ncols*i]);
#ifdef _USE_IC_INPUT_POWER_DELTA_
       pv[i]=static_cast<double>(prop[1+ncols*i])*(norm*ic_factor)*(2.*M_PI*M_PI)/pow(kv[i],3); //con factor dejo el P(k) a z=z_initial_simulation
#else
#ifdef _correct_shape_theoretica_power_
       pv[i]=static_cast<double>(prop[1+ncols*i])*(norm*ic_factor)/(1.0+0.05*kv[i]*kv[i]); //con factor dejo el P(k) a z=z_initial_simulation
#else
       pv[i]=static_cast<double>(prop[1+ncols*i])*(norm*ic_factor); //con factor dejo el P(k) a z=z_initial_simulation
#endif
#endif
     }
   So.DONE();
   prop.clear();prop.shrink_to_fit();
   gsl_interp_accel *acc = gsl_interp_accel_alloc();
   gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,kv.size());
   gsl_spline_init(spline,&kv[0],&pv[0],pv.size());
   real_prec deltak=this->params._d_DeltaK_data();
   real_prec ideltak=1.0/static_cast<double>(deltak);
   vector<real_prec> coords(Nft,0);// Coordinates of wavevectors in a regular grid
   for(ULONG i=0;i<Nft ;++i)
     coords[i]=deltak*(i<=Nft/2  ? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));
   vector<real_prec> mean_power(Nft/2/ndel,0);// Mean power to be updated
   vector<real_prec>IC_field_grf(Ngrid,0);
   for(ir=0;ir<this->params._Number_of_GRF();++ir)
     {
       So.message_screen("Realization ", ir);
#ifdef DOUBLE_PREC
       complex_prec *IC_FOURIER_GRF= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#else
       complex_prec *IC_FOURIER_GRF= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
#endif

#ifdef paired
       complex_prec *IC_FOURIER_GRF_PAIRED= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
       complex_prec *IC_FOURIER_GRF_PAIRED_2= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
       complex_prec *IC_FOURIER_GRF_PAIRED_4= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#endif
       complex_prec *IC_FOURIER_MIXED;
       complex_prec *IC_FOURIER_FA;

       vector<ULONG> nmodes(Nft/2/ndel,0);
       vector<real_prec> IC_GRF(Nft/2/ndel,0);
       vector<real_prec> IC_TH(Nft/2/ndel,0);
#ifdef paired
       vector<real_prec> IC_GRF_PAIRED(Nft/2/ndel,0);
       vector<real_prec> IC_GRF_PAIRED_2(Nft/2/ndel,0);
       vector<real_prec> IC_GRF_PAIRED_4(Nft/2/ndel,0);
#endif
       vector<real_prec> IC_MIXED;
       vector<real_prec> IC_FA;
       vector<real_prec> IC_GRF_lr;
       vector<real_prec> IC_MIXED_lr;
       vector<real_prec> IC_FA_lr;

       if(true==this->params._use_low_pass_filter())
     {
       IC_GRF_lr.resize(Nft_LR/2/ndel,0);
       if(true==this->params._Generate_FA())
         {
           IC_MIXED_lr.resize(Nft_LR/2/ndel,0);
           IC_FA_lr.resize(Nft_LR/2/ndel,0);
         }
     }

       if(true==this->params._Generate_FA())
     {
       IC_FA.resize(Nft/2/ndel,0);
       IC_MIXED.resize(Nft/2/ndel,0);
     }
       //************************************ parallel stuff****************
       // This is the right way to implement parallelization with random number generators //
#ifdef parallel_grf
       gsl_rng **gBaseRand = new gsl_rng *[_NTHREADS_];
       for (int b = 0; b < _NTHREADS_; b++)
     {
       gBaseRand[b] = gsl_rng_alloc(gsl_rng_mt19937);
       gsl_rng_set(gBaseRand[b], b * 101*(ir+1));
     }
#else
       gsl_rng * gBaseRand;
       const gsl_rng_type * rng_t;
#endif
       gsl_rng_env_setup();
       if(true==this->params._Generate_FA())
     {
       IC_FOURIER_MIXED= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
       IC_FOURIER_FA= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
       So.message_screen("Going for fixed amplitude");
#ifdef parallel_grf
#pragma omp parallel
       {
         int jthread=omp_get_thread_num();
#pragma omp for nowait
#else
         gsl_rng_default_seed= 5658*(ir+1);
         rng_t = gsl_rng_mt19937;
         gBaseRand = gsl_rng_alloc (rng_t);
#endif

         for(ULONG i=0;i<Nft;++i)
           {
         for(ULONG j=0;j< Nft;++j)
           {
             for(ULONG k=0;k<Nft/2+1;++k)
               {
             ULONG lp=index_3d(i,j,k,Nft,Nft/2+1);                                      // 3D C-like index
             real_prec kmod=sqrt(pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2)); // k=|k|= Magnitude of wave-vector
             ULONG ik=static_cast<ULONG>(floor(kmod/deltak));                                               // Bin in k
             real_prec Power_th = (kmod<kv[0] ? 0 : gsl_spline_eval(spline,kmod, acc));  // Theoretial Power interpolated at k
#ifdef parallel_grf
             real_prec phi_ic = 2.*M_PI*gsl_rng_uniform(gBaseRand[jthread]);                     // GR: Random phases. Add pi to get the inverted-phase pair
             real_prec Amp_ic = gsl_ran_rayleigh(gBaseRand[jthread],sqrt(0.5*Power_th));         // GR: Random Amplitude->Reyleigh distribution
#else
             real_prec phi_ic = 2.*M_PI*gsl_rng_uniform(gBaseRand);                     // GR: Random phases. Add pi to get the inverted-phase pair
             real_prec Amp_ic = gsl_ran_rayleigh(gBaseRand,sqrt(0.5*Power_th));         // GR: Random Amplitude->Reyleigh distribution
#endif
             real_prec cphi = cos(phi_ic);                                              // Cosinos of phase
             real_prec sphi = sin(phi_ic);                                              // Sin of phase
             IC_FOURIER_GRF[lp][REAL] = Amp_ic*cphi;                             // Real part of GR field
             IC_FOURIER_GRF[lp][IMAG] = Amp_ic*sphi;                             // Imag part of Mixed field
#ifdef paired
             real_prec cphi_paired = cos(phi_ic+M_PI);                                              // Cosinos of phase
             real_prec sphi_paired = sin(phi_ic+M_PI);                                              // Sin of phase
             IC_FOURIER_GRF_PAIRED[lp][REAL] = Amp_ic*cphi_paired;                             // Real part of GR field
             IC_FOURIER_GRF_PAIRED[lp][IMAG] = Amp_ic*sphi_paired;                             // Imag part of Mixed field
#endif
             real_prec Amp_faf = sqrt(Power_th);                                        // FA: Amplitude fixed to the power at k
             real_prec Amp_ic_final = kmod<=kmax ? Amp_faf: Amp_ic;                    // Final amplitude: below kmax, let it FA: above, let it GR
             IC_FOURIER_MIXED[lp][REAL] = Amp_ic_final*cphi;                                // Real part of Mixed field
             IC_FOURIER_MIXED[lp][IMAG] = Amp_ic_final*sphi;                                // Imag part of Mixed field
             IC_FOURIER_FA[lp][REAL] = Amp_faf*cphi;                                   // Real part of FA field
             IC_FOURIER_FA[lp][IMAG] = Amp_faf*sphi;                                   // Imag part of FA field
             // Do shell average to measure power of the three fields
             if(ik<Nft/2/ndel)
               {
#ifdef parallel_grf
#pragma omp atomic
#endif
                 nmodes[ik]++;    // Number of modes in k-shells
#ifdef parallel_grf
#pragma omp atomic
#endif
                 IC_MIXED[ik]+= pow(Amp_ic_final,2);  // Mixed
#ifdef parallel_grf
#pragma omp atomic
#endif
                 IC_FA[ik]+= pow(Amp_faf,2);      // FA
#ifdef parallel_grf
#pragma omp atomic
#endif
                 IC_TH[ik]+= pow(sqrt(Power_th),2); //Theoretical power evaluated at the Fourier mesh
#ifdef parallel_grf
#pragma omp atomic
#endif
                 IC_GRF[ik]+= pow(Amp_ic,2);  // GRF
               }
               }
           }
           }
#ifdef parallel_grf
       } // closes parallel region
#endif

     }// closes if

       else if (false==this->params._Generate_FA()) // this is only gaussian
     {

#ifdef parallel_grf
#pragma omp parallel
       {
         int  jthread=omp_get_thread_num();
#pragma omp for nowait
#else
             gsl_rng_default_seed= 5658*(ir+1587);
         rng_t = gsl_rng_mt19937;
         gBaseRand = gsl_rng_alloc (rng_t);
#endif

#ifdef white_noise
         // get a white noise
         So.message_screen("Getting white noise");
         for(ULONG i=0;i< IC_field_grf.size();++i)
           IC_field_grf[i]=gsl_rng_uniform(gBaseRand);
         // Get mean of WN
             real_prec meanWN=get_mean(IC_field_grf);
         // Ensure that WN has zero mean
         for(ULONG i=0;i< IC_field_grf.size();++i)
           IC_field_grf[i]-=meanWN;
         meanWN=0;
         for(ULONG i=0;i< IC_field_grf.size();++i)
           meanWN+=IC_field_grf[i];
         meanWN/=static_cast<real_prec>(IC_field_grf.size());
         this->So.message_screen("Mean WN =",meanWN);
         //Compute variance of the field sigma²
         real_prec sigmaWN=0;
         for(ULONG i=0;i< IC_field_grf.size();++i)
           sigmaWN+=(IC_field_grf[i]-meanWN)*(IC_field_grf[i]-meanWN);
         sigmaWN/=static_cast<real_prec>(IC_field_grf.size()-1);
         this->So.message_screen("Sigma² WN =",sigmaWN);

             // Nrmalize field to get sigma² = 1
         real_prec factor=static_cast<double>(sqrt(sigmaWN));
         for(ULONG i=0;i< IC_field_grf.size();++i)
           IC_field_grf[i]/=factor;

             for(int ui=0;ui<2;ui++)
           {
         meanWN=0;
         for(ULONG i=0;i< IC_field_grf.size();++i)
           meanWN+=IC_field_grf[i];
         meanWN/=static_cast<real_prec>(IC_field_grf.size());
         for(ULONG i=0;i< IC_field_grf.size();++i)
           IC_field_grf[i]-=meanWN;
         sigmaWN=0;
         for(ULONG i=0;i< IC_field_grf.size();++i)
           sigmaWN+=(IC_field_grf[i]-meanWN)*(IC_field_grf[i]-meanWN);
         sigmaWN/=static_cast<real_prec>(IC_field_grf.size()-1);
         factor=static_cast<double>(sqrt(sigmaWN));
         for(ULONG i=0;i< IC_field_grf.size();++i)
           IC_field_grf[i]/=factor;
         sigmaWN=0;
         for(ULONG i=0;i< IC_field_grf.size();++i)
           sigmaWN+=(IC_field_grf[i]-meanWN)*(IC_field_grf[i]-meanWN);
         sigmaWN/=static_cast<real_prec>(IC_field_grf.size()-1);
           }
         this->So.message_screen("New Sigma² WN =",sigmaWN);


         do_fftw_r2c(Nft ,IC_field_grf,IC_FOURIER_GRF);

         real_prec factorN=static_cast<real_prec>(Ngrid);
#endif
         So.message_screen("Going for GRF");
         // Loop over half of the Fourier box
             for(int i=0; i < Nft;++i)
               {
                 int i_m= (i <= Nft/2 ? i : Nft-i);
                 real_prec i_deltak_x =  static_cast<real_prec>(i_m *i_m) *this->params._d_deltak_x() * this->params._d_deltak_x();
                 for(int j=0; j < Nft; ++j)
           {
                     int j_m= (j <= Nft/2 ? j : Nft-j);
                     real_prec i_deltak_y = static_cast<real_prec>(j_m * j_m) * this->params._d_deltak_y() * this->params._d_deltak_y();
                     for(int k=0; k<=Nft/2 ;++k)
               {
                         real_prec i_deltak_z = static_cast<real_prec>(k * k) * this->params._d_deltak_z() * this->params._d_deltak_z();
             real_prec kmod = sqrt(i_deltak_x+i_deltak_y+i_deltak_z); // k=|k|= Magnitude of wave-vector
                         real_prec Power_th = 0;
                         if(kmod>=kv[0] && kmod<kv[kv.size()-1])
                           Power_th = gsl_spline_eval(spline,kmod, acc); // Theoretial Power interpolated at k

#ifndef white_noise
#ifdef parallel_grf
             real_prec Amp_ic = gsl_ran_rayleigh(gBaseRand[jthread],sqrt(0.5*Power_th));
                         real_prec phi_ic = 2.*M_PI*gsl_rng_uniform(gBaseRand[jthread]);
#else
             real_prec Amp_ic = static_cast<real_prec>(gsl_ran_rayleigh(gBaseRand,sqrt(0.5*Power_th)));
             real_prec phi_ic = 2.*M_PI*gsl_rng_uniform(gBaseRand);
#endif

             IC_FOURIER_GRF[lp][REAL] = Amp_ic*cos(phi_ic);
                         IC_FOURIER_GRF[lp][IMAG] = -Amp_ic*sin(phi_ic);
#else
             real_prec Amp_ic = sqrt(Power_th/factorN);
             ULONG lp = index_3d(i,j,k,Nft,Nft/2+1);
                         IC_FOURIER_GRF[lp][REAL] *= Amp_ic;
                         IC_FOURIER_GRF[lp][IMAG] *= Amp_ic;
#endif
                         int ik = static_cast<int>(floor((kmod-this->params._d_kmin())*ideltak));
                         if(ik<nmodes.size())
               {
#ifdef parallel_grf
#pragma omp atomic
#endif
                 nmodes[ik]++;    // Number of modes in k-shells
#ifdef parallel_grf
#pragma omp atomic
#endif
                 IC_TH[ik]+=Power_th; //Theoretical power evaluated at the Fourier mesh
#ifdef parallel_grf
#pragma omp atomic
#endif
                             IC_GRF[ik]+= IC_FOURIER_GRF[lp][REAL]*IC_FOURIER_GRF[lp][REAL]+ IC_FOURIER_GRF[lp][IMAG]*IC_FOURIER_GRF[lp][IMAG];
                           }
               }
           }
           }
#ifdef parallel_grf
       } // closes parallel region
#endif
     }//closes else
       So.DONE();
       So.message_screen("Getting power");
       if(true==this->params._Generate_FA())
     {
       for(ULONG i=0;i<IC_MIXED.size();++i)
         if(nmodes[i]>0)
           IC_MIXED[i]/=static_cast<real_prec>(nmodes[i]*norm);


       for(ULONG i=0;i<IC_FA.size();++i)
         if(nmodes[i]>0)
           IC_FA[i]/=static_cast<real_prec>(nmodes[i]*norm);
     }
       else if(false==this->params._Generate_FA())
     {
       for(ULONG i=0;i<IC_TH.size();++i)
         if(nmodes[i]>0)
           IC_TH[i]/=static_cast<real_prec>(nmodes[i]*norm);
       for(ULONG i=0;i<IC_GRF.size();++i)
         if(nmodes[i]>0)
               IC_GRF[i]/=(static_cast<real_prec>(nmodes[i])*norm);
     }
       So.DONE();
       // Vectors to allocate overdensities
       vector<real_prec>IC_field_fa;
       vector<real_prec>IC_field_mixed;
       do_fftw_c2r(Nft ,IC_FOURIER_GRF,IC_field_grf);
#ifdef DOUBLE_PREC
       fftw_free(IC_FOURIER_GRF);
#else
       fftwf_free(IC_FOURIER_GRF);
#endif
       if(true==this->params._Generate_FA())
     {
       IC_field_fa.resize(Ngrid,0);
       do_fftw_c2r(Nft ,IC_FOURIER_FA,IC_field_fa);
       IC_field_mixed.resize(Ngrid,0);
       do_fftw_c2r(Nft ,IC_FOURIER_MIXED,IC_field_mixed);
       fftw_free(IC_FOURIER_MIXED);
       fftw_free(IC_FOURIER_FA);
     }
       So.DONE();
       vector<real_prec>IC_field_grf_lr(Ngrid_LR,0);
       vector<real_prec>IC_field_grf_paired_lr(Ngrid_LR,0);
       vector<real_prec>IC_field_fa_lr;
       vector<real_prec>IC_field_mixed_lr;
       if(true==this->params._use_low_pass_filter())
     {
       So.message_screen("Applying low pass filter: ");

       real_prec mean_HR=static_cast<real_prec>(Ngrid_HR)/pow(params._Lbox(),3);
       for(ULONG i=0; i<IC_field_grf.size();++i)
         IC_field_grf[i]=mean_HR*(1.0+IC_field_grf[i]);

           low_pass_filter(Nft_HR,Nft_LR,1,false, IC_field_grf,IC_field_grf_lr, params._Lbox());
       get_overdens(IC_field_grf_lr,IC_field_grf_lr);

       this->File.write_array(this->params._Output_directory()+"IC_GRF_lpass_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_grf_lr);
       this->File.write_array(this->params._Output_directory()+"IC_GRF_lpass_paired_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_grf_paired_lr);
       this->params.set_input_type("delta_grid");
       this->params.set_Nft(Nft_LR);
       this->params.set_NGRID(Ngrid_LR);
       this->params.set_SN_correction(false);
       this->params.set_MAS_correction(false);
           this->compute_power_spectrum_grid(IC_field_grf_lr, false);
       for(ULONG i=0;i<IC_GRF_lr.size();++i)
         IC_GRF_lr[i]=this->pk0[i];

       if(true==this->params._Generate_FA())
         {
           IC_field_fa_lr.resize(Ngrid_LR,0);
           IC_field_mixed_lr.resize(Ngrid_LR,0);
               low_pass_filter(Nft_HR,Nft_LR,1,false, IC_field_fa,IC_field_fa_lr, params._Lbox());
           this->File.write_array(this->params._Output_directory()+"IC_FA_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_fa_lr) ;
               low_pass_filter(Nft_HR,Nft_LR,1,false, IC_field_mixed,IC_field_mixed_lr, params._Lbox());
           this->File.write_array(this->params._Output_directory()+"IC_mixed_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_mixed_lr);
         }
       So.DONE();
     }
       else
     {
       this->File.write_array(this->params._Output_directory()+"IC_GRF_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_grf);
       if(true==this->params._Generate_FA())
         {
           this->File.write_array(this->params._Output_directory()+"IC_FA_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_fa) ;
           this->File.write_array(this->params._Output_directory()+"IC_mixed_realization"+to_string(ir)+"_Nres"+to_string(Nft_LR), IC_field_mixed);
         }
     }
#ifdef _USE_GNUPLOT_
       vector<pair<real_prec, real_prec> > xy_pts_fa;
       vector<pair<real_prec, real_prec> > xy_pts_mixed;
       vector<pair<real_prec, real_prec> > xy_pts_grf;
       vector<pair<real_prec, real_prec> > xy_pts_th;
       vector<pair<real_prec, real_prec> > xy_pts_grf_lr;
       vector<pair<real_prec, real_prec> > xy_pts_mean_power;
       for(ULONG i=0;i<IC_GRF.size();++i) // Update mean power
          mean_power[i]+=(IC_GRF[i]-mean_power[i])/static_cast<real_prec>(ir+1);
       for(ULONG i=1;i<mean_power.size();++i)
          xy_pts_mean_power.push_back(std::make_pair(deltak*(i+0.5), (mean_power[i])));
       for(ULONG i=1;i<IC_GRF.size();++i)
         xy_pts_grf.push_back(std::make_pair(deltak*(i+0.5), (IC_GRF[i])));
      if(true==this->params._use_low_pass_filter())
         for(ULONG i=1;i<IC_GRF_lr.size();++i)
           xy_pts_grf_lr.push_back(std::make_pair(deltak*(i+0.5), (IC_GRF_lr[i])));
       for(ULONG i=1;i<IC_TH.size();++i)
     xy_pts_th.push_back(std::make_pair(deltak*(i+0.5), (IC_TH[i])));
       if(true==this->params._Generate_FA())
     {
       for(ULONG i=1;i< IC_FA.size();++i)
         xy_pts_fa.push_back(std::make_pair(deltak*(i+0.5), IC_FA[i]));

       for(ULONG i=1;i<IC_MIXED.size();++i)
         xy_pts_mixed.push_back(std::make_pair(deltak*(i+0.5), (IC_MIXED[i])));
     }
       this->gp_pdf<<"set size 1,1\n";
       this->gp_pdf<<"set origin 0., 0.\n";
       this->gp_pdf<<"set grid\n";
       this->gp_pdf << "set log x \n";
       this->gp_pdf << "set log y \n";
       this->gp_pdf << "set ylabel 'Power Spectrum P(k)'\n";
       this->gp_pdf << "set xlabel 'k [Mpc/h]'\n";
       if(ir==0)
     {
       if(true==this->params._Generate_FA())
         {
#ifdef _USE_IC_INPUT_POWER_DELTA_
           this->gp_pdf << "plot[0.003:3][0.001:50]" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 lt 11  title 'GRF', "<<this->gp_pdf.file1d(xy_pts_fa)<<"w l lw 3 lt 2 title 'FA',"<<this->gp_pdf.file1d(xy_pts_mixed)<<"w l lw 3 lt 4 title 'MIXED', "<<gp_pdf.file1d(xy_pts_th)<<" w l lw 1 title 'CAMB_g','"<<power_file<<"' u 1:($2*(2.0*acos(-1.0)*acos(-1.0))/($1*$1*$1)) w l lt 2 ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 1 title 'GRF_lr'"<<endl;
#else
#ifdef _ABACUS_
               this->gp_pdf<<"factor=(1./47.304805056)**2 \n";
#elif defined _UNITSIM_
               this->gp_pdf<<"factor=0.000158091416709513\n";
#endif
               this->gp_pdf<<"plot[0.003:3][0.001:50]" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 lt 11  title 'GRF', "<<this->gp_pdf.file1d(xy_pts_fa)<<"w l lw 3 lt 2 title 'FA',"<<this->gp_pdf.file1d(xy_pts_mixed)<<"w l lw 3 lt 4 title 'MIXED', "<<gp_pdf.file1d(xy_pts_th)<<" w l lw 1 title 'CAMB_g','"<<power_file<<"' u 1:($2*factor) w l lt 2 ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 1 title 'GRF_lr'"<<endl;
#endif
         }
       else if (false==this->params._Generate_FA())
         {
#ifdef _USE_IC_INPUT_POWER_DELTA_
           this->gp_pdf << "plot[0.003:3][0.001:50]" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 lt 11  title 'GRF', "<<gp_pdf.file1d(xy_pts_th)<<"w l lw 3 title 'CAMB_g', '"<<power_file<<"' u 1:($2*(2.0*acos(-1.0)*acos(-1.0))/($1*$1*$1)) w l lt 2 ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 4 title 'GRF_lr'"<<endl;
#else
#ifdef _ABACUS_
               this->gp_pdf<<"factor=(1./47.304805056)**2 \n";
#elif defined _UNITSIM_
               this->gp_pdf<<"factor=0.000158091416709513\n";
#endif
           this->gp_pdf << "plot[0.003:3][0.001:50]" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 lt 11  title 'GRF', "<<gp_pdf.file1d(xy_pts_th)<<"w l lw 3 title 'CAMB_g', '"<<power_file<<"'  u 1:($2*factor) w l lt 2 ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 4 title 'GRF_lr'"<<endl;
#endif
         }
     }
       else if(ir>0)
     {
       if(true==this->params._Generate_FA())
         this->gp_pdf << "replot" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 lt 11 notitle , "<<this->gp_pdf.file1d(xy_pts_fa)<<"w l lw 1 lt 2 notitle ,"<<this->gp_pdf.file1d(xy_pts_mixed)<<"w l lw 3 lt 4 notitle ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 1 notitle"<<endl;
       else
             this->gp_pdf << "replot" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 notitle ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 1 notitle"<<endl;
//           this->gp_pdf << "replot" << this->gp_pdf.file1d(xy_pts_grf) << "w l lw 1 notitle ,"<<gp_pdf.file1d(xy_pts_grf_lr)<<" w l lw 1 notitle, "<<gp_pdf.file1d(xy_pts_mean_power)<<" w l lw 3 lt 3 notitle "<<endl;
      // this->gp_pdf << "replot" <<gp_pdf.file1d(xy_pts_mean_power)<<" w l lw 3 lt 3 notitle "<<endl;
     }
#endif
       if(ir==0)
        {
       ofstream atea; atea.open(this->params._Output_directory()+"mean_power_th.txt");
       atea.precision(6);
       atea.setf(ios::showpoint);
       atea.setf(ios::scientific);
       for(int i=1;i<IC_TH.size();++i)
        atea<<deltak*(i+0.5)<<"  "<<IC_TH[i]<<endl;
       atea.close();
        }
       ofstream btea; btea.open(this->params._Output_directory()+"power_grf_real"+to_string(ir)+".txt");
       btea.precision(6);
       btea.setf(ios::showpoint);
       btea.setf(ios::scientific);

       for(int i=1;i<IC_GRF.size();++i)
        btea<<deltak*(i+0.5)<<"\t"<<IC_GRF[i]<<"\t"<<nmodes[i]<<endl;
       btea.close();
       }// closes loop over ir
   ofstream tea; tea.open(this->params._Output_directory()+"mean_power.txt");
   for(ULONG i=0;i<mean_power.size();++i)
    tea<<deltak*(i+0.5)<<"  "<<mean_power[i]<<endl;
   tea.close();
   gsl_spline_free(spline);
   gsl_interp_accel_free(acc);
   So.DONE();
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 pair<real_prec, real_prec> PowerSpectrumF::get_lss_cross_bias(vector<real_prec>&p1, vector<real_prec>&p2,  vector<real_prec>&kv,  vector<int>&nmodes, int init_mode, real_prec kmax)
 {
   real_prec lss_bias_k=0;
   real_prec lss_bias=0;
   int ik=init_mode;
   while(kv[ik]<kmax)
     {
       lss_bias+=(p1[ik]/p2[ik])*static_cast<real_prec>(nmodes[ik]);
       lss_bias_k+=static_cast<real_prec>(nmodes[ik]);
       ++ik;
     }
   lss_bias/=lss_bias_k;
   real_prec var_lss_bias=0;
   ik=init_mode;
   while(kv[ik]<kmax)
     {
       var_lss_bias+=pow(p1[ik]/p2[ik]-lss_bias,2)*static_cast<real_prec>(nmodes[ik]);
       ++ik;
     }
   var_lss_bias=sqrt(var_lss_bias/lss_bias_k);
   return std::make_pair(lss_bias, var_lss_bias);
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
pair<real_prec, real_prec> PowerSpectrumF::get_lss_bias(vector<real_prec>&p1, vector<real_prec>&p2,  vector<real_prec>&kv,  vector<int>&nmodes, int init_mode, real_prec kmax)
 {
   real_prec lss_bias_k=0;
   real_prec lss_bias=0;
   int ik=init_mode;
   while(kv[ik]<kmax)
     {
       lss_bias+=sqrt(p1[ik]/p2[ik])*nmodes[ik];
       lss_bias_k+=static_cast<real_prec>(nmodes[ik]);
       ++ik;
     }
   lss_bias/=lss_bias_k;
   real_prec var_lss_bias=0;
   ik=init_mode;
   while(kv[ik]<kmax)
     {
       var_lss_bias+=pow(sqrt(p1[ik]/p2[ik])-lss_bias,2)*nmodes[ik];
       ++ik;
     }
   var_lss_bias=sqrt(var_lss_bias/lss_bias_k);
  return std::make_pair(lss_bias, var_lss_bias);
 }
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::GetSuperClusters(string realm, string tbias){
    int NTHREADS = 1;
#ifdef _USE_OMP_
    NTHREADS=_NTHREADS_;
#endif
    string namef;
    if(tbias=="_BIAS_")
        namef= (realm == "over"? "super": "usuper");
    else if(tbias=="_RELATIVE_BIAS_")
        namef= (realm == "over"? "rsuper": "rusuper");
   vector<real_prec> DM_DEN_FIELD(this->params._NGRID(),0);
   vector<real_prec> TR_DEN_FIELD(this->params._NGRID(),0);   
   this->File.read_array(this->params._Input_Directory_X()+this->params._Name_Catalog_X(),DM_DEN_FIELD);
   this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),TR_DEN_FIELD);
   this->params.set_input_type("density_grid");
   this->set_output_filenames();
   this->So.message_screen("Getting object-to-object bias");
   Catalog cat(this->params);
   cat.read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),0);
    this->set_params(cat.params); // params has been updated inside the catalog class.
   kvector_data.clear();
   kvector_data.shrink_to_fit();
   for(int i=0;i<this->params._d_Nnp_data();i++)
     kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
   this->object_by_object_bias(cat.Halo,DM_DEN_FIELD,TR_DEN_FIELD);
   this->So.message_screen("\tMin Bias  =", cat.get_min("_BIAS_"));
   this->So.message_screen("\tMax Bias  =", cat.get_max("_BIAS_"));
   this->So.message_screen("\tMin rBias  =", cat.get_min("_RELATIVE_BIAS_"));
   this->So.message_screen("\tMax rBias  =", cat.get_max("_RELATIVE_BIAS_"));
   cat.get_superclusters(realm, tbias); // Here we selec tracers with bias above the fourth quartile. This gives tracer_aux
   ULONG size_cat=cat.tracer_aux.size();
   int Nrbins=this->params._Nbins_cf();
   real_prec rmax=this->params._rmax_cf();
   real_prec rmin=this->params._rmin_cf();
   real_prec rmin_squared = rmin*rmin;
   real_prec rmax_squared = rmax*rmax;
   real_prec lrmin=log10(rmin);
   real_prec deltar=log10(rmax/rmin)/static_cast<real_prec>(Nrbins);
   // here we shuffle properties
   So.message_screen("  Ranking properties");
   cat.Get_Ranked_Props(cat.tracer_aux,"_MASS_");
   cat.Get_Ranked_Props(cat.tracer_aux,tbias);
   cat.Get_Ranked_Props(cat.tracer_aux,"_SPIN_");
   cat.Get_Ranked_Props(cat.tracer_aux,"_CONCENTRATION_");
   cat.Get_Ranked_Props(cat.tracer_aux,"_MACH_"); // this returns Halo_ranked
   cat.Get_Ranked_Props(cat.tracer_aux,"_BTOA_");
   cat.Get_Ranked_Props(cat.tracer_aux,"_CTOA_");
   cat.Get_Ranked_Props(cat.tracer_aux,"_LOCALDM_");
   cat.Get_Ranked_Props(cat.tracer_aux,"_TIDAL_ANISOTROPY_");
   So.DONE();
   //******************************************************
   //   Halo_ranked and tracer_aux have the same dimension, 
   // same coordinates in the i-label, only properties have been ranked
   //******************************************************
   //   jacknife
     So.message_screen("Getting Jackknife estimates");
   So.message_screen("using ",this->params.d_NGRID_JK()," 1time deleted sub-samples");
#pragma omp parallel for
   for(ULONG i=0;i<size_cat;++i) 
        cat.tracer_aux[i].GridID_n = grid_ID(&cat.box_JK, cat.tracer_aux[i].coord1,cat.tracer_aux[i].coord2,cat.tracer_aux[i].coord3);
   for(ULONG rr=0;rr<=this->params.d_NGRID_JK();++rr) // Loop over each of the sub-samples JK. The 0th item is the full sample
     {

       vector<real_prec>xc;
       vector<real_prec>yc;
       vector<real_prec>zc;
       vector<real_prec>mc;
       vector<real_prec>sc;
       vector<real_prec>cc;
       vector<real_prec>bc;
       vector<real_prec>mach;
       vector<real_prec>tidalc;
       vector<real_prec>localdmc;
       vector<real_prec>btoac;
       vector<real_prec>ctoac;
       real_prec mean_mark_m=0;
       real_prec mean_mark_b=0;
       real_prec mean_mark_s=0;
       real_prec mean_mark_mach=0;
       real_prec mean_mark_concentration=0;
       real_prec mean_mark_btoa=0;
       real_prec mean_mark_ctoa=0;
       real_prec mean_mark_lcdm=0;
       real_prec mean_mark_ta=0;
    if(rr==0)
    {
       for(ULONG i=0;i<size_cat;++i)// Loop over the catalog
       {
         xc.push_back(cat.tracer_aux[i].coord1);
         yc.push_back(cat.tracer_aux[i].coord2);
         zc.push_back(cat.tracer_aux[i].coord3);
         mc.push_back(cat.Halo_ranked[i].mass);
         if(tbias=="_BIAS_")
             bc.push_back(cat.Halo_ranked[i].bias);
         if(tbias=="_RELATIVE_BIAS_")
             bc.push_back(cat.Halo_ranked[i].relative_bias);
         sc.push_back(cat.Halo_ranked[i].spin);
         mach.push_back(cat.Halo_ranked[i].mach_number);
         cc.push_back(cat.Halo_ranked[i].concentration);
         tidalc.push_back(cat.Halo_ranked[i].tidal_anisotropy);
         btoac.push_back(cat.Halo_ranked[i].b_to_a);
         ctoac.push_back(cat.Halo_ranked[i].c_to_a);
         localdmc.push_back(cat.Halo_ranked[i].local_dm);
         mean_mark_m+=cat.Halo_ranked[i].mass;
         if(tbias=="_BIAS_")
            mean_mark_b+=cat.Halo_ranked[i].bias;
         if(tbias=="_RELATIVE_BIAS_")
            mean_mark_b+=cat.Halo_ranked[i].relative_bias;
         mean_mark_s+=cat.Halo_ranked[i].spin;
         mean_mark_mach+=cat.Halo_ranked[i].mach_number;
         mean_mark_concentration+=cat.Halo_ranked[i].concentration;
         mean_mark_btoa+=cat.Halo_ranked[i].b_to_a;
         mean_mark_ctoa+=cat.Halo_ranked[i].c_to_a;
         mean_mark_lcdm+=cat.Halo_ranked[i].local_dm;
         mean_mark_ta+=cat.Halo_ranked[i].tidal_anisotropy;
       }
    }
      else{
       for(ULONG i=0;i<size_cat;++i)// Loop over the catalog
    	 if(cat.tracer_aux[i].GridID_n != rr)
    	   {
    	     xc.push_back(cat.tracer_aux[i].coord1);
    	     yc.push_back(cat.tracer_aux[i].coord2);
    	     zc.push_back(cat.tracer_aux[i].coord3);
    	     mc.push_back(cat.Halo_ranked[i].mass);
             if(tbias=="_BIAS_")
        	     bc.push_back(cat.Halo_ranked[i].bias);
    	     if(tbias=="_RELATIVE_BIAS_")
                 bc.push_back(cat.Halo_ranked[i].relative_bias);
             sc.push_back(cat.Halo_ranked[i].spin);
             mach.push_back(cat.Halo_ranked[i].mach_number);
             cc.push_back(cat.Halo_ranked[i].concentration);
             tidalc.push_back(cat.Halo_ranked[i].tidal_anisotropy);
             btoac.push_back(cat.Halo_ranked[i].b_to_a);
             ctoac.push_back(cat.Halo_ranked[i].c_to_a);
             localdmc.push_back(cat.Halo_ranked[i].local_dm);
             mean_mark_m+=cat.Halo_ranked[i].mass;
             if(tbias=="_BIAS_")
                 mean_mark_b+=cat.Halo_ranked[i].bias;
             if(tbias=="_RELATIVE_BIAS_")
                 mean_mark_b+=cat.Halo_ranked[i].relative_bias;
             mean_mark_b+=cat.Halo_ranked[i].bias;
    	     mean_mark_s+=cat.Halo_ranked[i].spin;
    	     mean_mark_mach+=cat.Halo_ranked[i].mach_number;
             mean_mark_concentration+=cat.Halo_ranked[i].concentration;
             mean_mark_btoa+=cat.Halo_ranked[i].b_to_a;
             mean_mark_ctoa+=cat.Halo_ranked[i].c_to_a;
             mean_mark_lcdm+=cat.Halo_ranked[i].local_dm;
             mean_mark_ta+=cat.Halo_ranked[i].tidal_anisotropy;
	      }
       }
       ULONG nsize_cat=xc.size();    
       So.message_screen("Number of tracers before removing subsample:", size_cat);
       So.message_screen("Number of tracers after removing subsample:", nsize_cat);
       mean_mark_m/=static_cast<real_prec>(nsize_cat);
       mean_mark_b/=static_cast<real_prec>(nsize_cat);
       mean_mark_s/=static_cast<real_prec>(nsize_cat);
       mean_mark_mach/=static_cast<real_prec>(nsize_cat);
       mean_mark_concentration/=static_cast<real_prec>(nsize_cat);
       mean_mark_ctoa/=static_cast<real_prec>(nsize_cat);
       mean_mark_btoa/=static_cast<real_prec>(nsize_cat);
       mean_mark_ta/=static_cast<real_prec>(nsize_cat);
       mean_mark_lcdm/=static_cast<real_prec>(nsize_cat);
       So.message_screen("Computing marked statistics for JK=", rr);
       time_t start;
       time_t end;
       time (&start);
       vector<real_prec>gcount(Nrbins,0);
       vector<real_prec>mcount(Nrbins,0);
       vector<real_prec>bcount(Nrbins,0);
       vector<real_prec>scount(Nrbins,0);
       vector<real_prec>machcount(Nrbins,0);
       vector<real_prec>ccount(Nrbins,0);
       vector<real_prec>btoacount(Nrbins,0);
       vector<real_prec>ctoacount(Nrbins,0);
       vector<real_prec>lcdmcount(Nrbins,0);
       vector<real_prec>tacount(Nrbins,0);
#pragma omp parallel 
       {
     vector<real_prec>gcount_priv(Nrbins,0);
     vector<real_prec>mcount_priv(Nrbins,0);
	 vector<real_prec>bcount_priv(Nrbins,0);
     vector<real_prec>scount_priv(Nrbins,0);
	 vector<real_prec>machcount_priv(Nrbins,0);
     vector<real_prec>ccount_priv(Nrbins,0);
     vector<real_prec>btoacount_priv(Nrbins,0);
     vector<real_prec>ctoacount_priv(Nrbins,0);
     vector<real_prec>lcdmcount_priv(Nrbins,0);
     vector<real_prec>tacount_priv(Nrbins,0);
#pragma omp for collapse(2)
	 for(ULONG i=0;i<nsize_cat;++i)
	   for(ULONG j=i+1;j<nsize_cat;++j)
	     {
	       real_prec dx=xc[i]-xc[j];
	       real_prec dy=yc[i]-yc[j];
	       real_prec dz=zc[i]-zc[j];
	       real_prec dist_squared=dx*dx+dy*dy+dz*dz;
	       if (dist_squared>=rmin_squared && dist_squared<rmax_squared)
		 {
		   real_prec ldist=0.5*log10(dist_squared); // this is log (sqrt(dist_squared))
		   ULONG ibin=get_bin(ldist,lrmin,Nrbins,deltar,false);
		   if(ibin<Nrbins) 
		     {
		       gcount_priv[ibin]++;  // marking with property=mass
               mcount_priv[ibin]+=mc[i]*mc[j];
		       bcount_priv[ibin]+=bc[i]*bc[j]; // marking with bias
		       scount_priv[ibin]+=sc[i]*sc[j];
               machcount_priv[ibin]+=mach[i]*mach[j];
               ccount_priv[ibin]+=cc[i]*cc[j];
               btoacount_priv[ibin]+=btoac[i]*btoac[j];
               ctoacount_priv[ibin]+=ctoac[i]*ctoac[j];
               lcdmcount_priv[ibin]+=localdmc[i]*localdmc[j];
               ctoacount_priv[ibin]+=tidalc[i]*tidalc[j];
		     }
		 }
	     }
#pragma omp critical
	 {
	   for(ULONG i=0;i<Nrbins;++i)
	     {
	       gcount[i]+= gcount_priv[i];
           mcount[i]+= mcount_priv[i]/pow(mean_mark_m,2);
	       bcount[i]+= bcount_priv[i]/pow(mean_mark_b,2);
           scount[i]+= scount_priv[i]/pow(mean_mark_s,2);
	       machcount[i]+= machcount_priv[i]/pow(mean_mark_mach,2);
           ccount[i]+= ccount_priv[i]/pow(mean_mark_concentration,2);
           btoacount[i]+= btoacount_priv[i]/pow(mean_mark_btoa,2);
           ctoacount[i]+= ctoacount_priv[i]/pow(mean_mark_ctoa,2);
           lcdmcount[i]+= lcdmcount_priv[i]/pow(mean_mark_lcdm,2);
           tacount[i]+= tacount_priv[i]/pow(mean_mark_lcdm,2);
	     }
	 }
	 scount_priv.clear();scount_priv.shrink_to_fit();
     ccount_priv.clear();ccount_priv.shrink_to_fit();
	 gcount_priv.clear();gcount_priv.shrink_to_fit();
	 mcount_priv.clear();mcount_priv.shrink_to_fit();
	 machcount_priv.clear();machcount_priv.shrink_to_fit();
	 bcount_priv.clear();bcount_priv.shrink_to_fit();
     ctoacount_priv.clear();ctoacount_priv.shrink_to_fit();
     btoacount_priv.clear();btoacount_priv.shrink_to_fit();
     lcdmcount_priv.clear();lcdmcount_priv.shrink_to_fit();
     tacount_priv.clear();tacount_priv.shrink_to_fit();

       }//Closes parallel region
       
       So.DONE();
       time (&end);
       So.message_screen("CT[secs]:", difftime(end,start));
       // We now write rbin,      
       ofstream profd;
       string file;

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_mass_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
	 profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<mcount[j]<<"\t"<<mcount[j]/gcount[j]<<endl;
       profd.close();
       
       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_bias_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
	 profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<bcount[j]<<"\t"<<bcount[j]/gcount[j]<<endl;
       profd.close();
       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_spin_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
	 profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<scount[j]<<"\t"<<scount[j]/gcount[j]<<endl;
       profd.close();
       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_mach_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
	 profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<machcount[j]<<"\t"<<machcount[j]/gcount[j]<<endl;
       profd.close();

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_concentration_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
        profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<ccount[j]<<"\t"<<ccount[j]/gcount[j]<<endl;
       profd.close();

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_btoa_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
        profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<btoacount[j]<<"\t"<<btoacount[j]/gcount[j]<<endl;
       profd.close();

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_ctoa_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
        profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<ctoacount[j]<<"\t"<<ctoacount[j]/gcount[j]<<endl;
       profd.close();

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_lcdm_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
        profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<lcdmcount[j]<<"\t"<<lcdmcount[j]/gcount[j]<<endl;
       profd.close();

       file=this->params._Output_directory()+"separation_distribtion_"+namef+"_ta_jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
        profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<tacount[j]<<"\t"<<tacount[j]/gcount[j]<<endl;
       profd.close();

     }// End loop over the JK regions


   //******************************************************
   // RANDOMIZAITION OF MARks

  int Nruns=100;// Numbner of shuffling marks to get errors

  gsl_rng_env_setup();
  const gsl_rng_type *Tn;
  gsl_rng_default_seed=1015;
  Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (Tn);

  vector<ULONG>aux_p;
  for(ULONG i=0;i<size_cat;++i)
    aux_p.push_back(i);

  for(int rr=0;rr<=Nruns;++rr)
    {

   if(rr>0)
   {
   So.message_screen("  Shuffling properties for run ", rr);
      gsl_ran_shuffle(rn,&aux_p[0], aux_p.size() ,sizeof(ULONG));
      for(ULONG i=0;i<size_cat;++i)
       {
        cat.Halo_ranked[i].mass=cat.Halo_ranked[aux_p[i]].mass;
        if(tbias=="_BIAS_")
          cat.Halo_ranked[i].bias=cat.Halo_ranked[aux_p[i]].bias;
        if(tbias=="_RELATIVE_BIAS_")
          cat.Halo_ranked[i].relative_bias=cat.Halo_ranked[aux_p[i]].relative_bias;
        cat.Halo_ranked[i].spin=cat.Halo_ranked[aux_p[i]].spin;
        cat.Halo_ranked[i].mach_number=cat.Halo_ranked[aux_p[i]].mach_number;
        cat.Halo_ranked[i].concentration=cat.Halo_ranked[aux_p[i]].concentration;
        cat.Halo_ranked[i].c_to_a=cat.Halo_ranked[aux_p[i]].c_to_a;
        cat.Halo_ranked[i].b_to_a=cat.Halo_ranked[aux_p[i]].b_to_a;
        cat.Halo_ranked[i].local_dm=cat.Halo_ranked[aux_p[i]].local_dm;
        cat.Halo_ranked[i].tidal_anisotropy=cat.Halo_ranked[aux_p[i]].tidal_anisotropy;
       }
       So.DONE();
    }
   real_prec mean_mark_m=0;
   real_prec mean_mark_b=0;
   real_prec mean_mark_s=0;
   real_prec mean_mark_concentration=0;
   real_prec mean_mark_mach=0;
   real_prec mean_mark_btoa=0;
   real_prec mean_mark_ctoa=0;
   real_prec mean_mark_lcdm=0;
   real_prec mean_mark_ta=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_mark_b,mean_mark_m,mean_mark_s,mean_mark_concentration, mean_mark_mach, mean_mark_lcdm, mean_mark_btoa, mean_mark_ctoa, mean_mark_ta)
#endif  
   for(ULONG i=0;i<size_cat;++i)
     {
       mean_mark_m+=cat.Halo_ranked[i].mass;
      if(tbias=="_BIAS_")
           mean_mark_b+=cat.Halo_ranked[i].bias;
      if(tbias=="_RELATIVE_BIAS_")
           mean_mark_b+=cat.Halo_ranked[i].relative_bias;
       mean_mark_s+=cat.Halo_ranked[i].spin;
       mean_mark_mach+=cat.Halo_ranked[i].mach_number;
       mean_mark_concentration+=cat.Halo_ranked[i].concentration;
       mean_mark_btoa+=cat.Halo_ranked[i].b_to_a;
       mean_mark_ctoa+=cat.Halo_ranked[i].c_to_a;
       mean_mark_lcdm+=cat.Halo_ranked[i].local_dm;
       mean_mark_ta+=cat.Halo_ranked[i].tidal_anisotropy;
     }
   mean_mark_m/=static_cast<real_prec>(size_cat);
   mean_mark_b/=static_cast<real_prec>(size_cat);
   mean_mark_s/=static_cast<real_prec>(size_cat);
   mean_mark_concentration/=static_cast<real_prec>(size_cat);
   mean_mark_mach/=static_cast<real_prec>(size_cat);
   mean_mark_ctoa/=static_cast<real_prec>(size_cat);
   mean_mark_btoa/=static_cast<real_prec>(size_cat);
   mean_mark_ta/=static_cast<real_prec>(size_cat);
   mean_mark_lcdm/=static_cast<real_prec>(size_cat);

   So.message_screen("Computing marked statistics");
   time_t start;
   time_t end;
   time (&start);


   vector<real_prec>gcount(Nrbins,0);
   vector<real_prec>mcount(Nrbins,0);
   vector<real_prec>bcount(Nrbins,0);
   vector<real_prec>scount(Nrbins,0);
   vector<real_prec>machcount(Nrbins,0);
   vector<real_prec>ccount(Nrbins,0);
   vector<real_prec>btoacount(Nrbins,0);
   vector<real_prec>ctoacount(Nrbins,0);
   vector<real_prec>lcdmcount(Nrbins,0);
   vector<real_prec>tacount(Nrbins,0);
#pragma omp parallel 
   {
     vector<real_prec> gcount_priv(Nrbins,0);
     vector<real_prec> mcount_priv(Nrbins,0);
     vector<real_prec> bcount_priv(Nrbins,0);
     vector<real_prec> scount_priv(Nrbins,0);
     vector<real_prec> machcount_priv(Nrbins,0);
     vector<real_prec> ccount_priv(Nrbins,0);
      vector<real_prec>btoacount_priv(Nrbins,0);
     vector<real_prec>ctoacount_priv(Nrbins,0);
     vector<real_prec>lcdmcount_priv(Nrbins,0);
     vector<real_prec>tacount_priv(Nrbins,0);
#pragma omp for collapse(2)
     for(ULONG i=0;i<size_cat;++i)
       for(ULONG j=i+1;j<size_cat;++j)
    	 {
    	   real_prec dx=cat.tracer_aux[i].coord1-cat.tracer_aux[j].coord1;
    	   real_prec dy=cat.tracer_aux[i].coord2-cat.tracer_aux[j].coord2;
    	   real_prec dz=cat.tracer_aux[i].coord3-cat.tracer_aux[j].coord3;
    	   real_prec dist_squared=dx*dx+dy*dy+dz*dz;
           if (dist_squared>=rmin_squared && dist_squared<rmax_squared)
           {
        	   real_prec ldist=0.5*log10(dist_squared); // this is log (sqrt(dist_squared))
    	       ULONG ibin=get_bin(ldist,lrmin,Nrbins,deltar,false);
    	       if(ibin<Nrbins) 
    	       {
        	       gcount_priv[ibin]++;  // marking with property=mass
                   mcount_priv[ibin]+=cat.Halo_ranked[i].mass*cat.Halo_ranked[j].mass;
                   if(tbias=="_BIAS_")
            	       bcount_priv[ibin]+=cat.Halo_ranked[i].bias*cat.Halo_ranked[j].bias; // marking with bias
                   if(tbias=="_RELATIVE_BIAS_")
                       bcount_priv[ibin]+=cat.Halo_ranked[i].relative_bias*cat.Halo_ranked[j].relative_bias; // marking with bias
        	       scount_priv[ibin]+=cat.Halo_ranked[i].spin*cat.Halo_ranked[j].spin;
                   machcount_priv[ibin]+=cat.Halo_ranked[i].mach_number*cat.Halo_ranked[j].mach_number; // marking with bias       
                   ccount_priv[ibin]+=cat.Halo_ranked[i].concentration*cat.Halo_ranked[j].concentration;
                   btoacount_priv[ibin]+=cat.Halo_ranked[i].b_to_a*cat.Halo_ranked[j].b_to_a;
                   ctoacount_priv[ibin]+=cat.Halo_ranked[i].c_to_a*cat.Halo_ranked[j].c_to_a;
                   lcdmcount_priv[ibin]+=cat.Halo_ranked[i].local_dm*cat.Halo_ranked[j].local_dm;
                   ctoacount_priv[ibin]+=cat.Halo_ranked[i].tidal_anisotropy*cat.Halo_ranked[j].tidal_anisotropy;
    	       }
    	    }
         }
#pragma omp critical
     {
       for(ULONG i=0;i<Nrbins;++i)
	    {
	       gcount[i]+= gcount_priv[i];
           mcount[i]+= mcount_priv[i]/pow(mean_mark_m,2);
    	   bcount[i]+= bcount_priv[i]/pow(mean_mark_b,2);
           scount[i]+= scount_priv[i]/pow(mean_mark_s,2);
    	   machcount[i]+= machcount_priv[i]/pow(mean_mark_mach,2);
           ccount[i]+= ccount_priv[i]/pow(mean_mark_concentration,2);
           btoacount[i]+= btoacount_priv[i]/pow(mean_mark_btoa,2);
           ctoacount[i]+= ctoacount_priv[i]/pow(mean_mark_ctoa,2);
           lcdmcount[i]+= lcdmcount_priv[i]/pow(mean_mark_lcdm,2);
           tacount[i]+= tacount_priv[i]/pow(mean_mark_lcdm,2);
    	}
     }
     scount_priv.clear();scount_priv.shrink_to_fit();
     ccount_priv.clear();ccount_priv.shrink_to_fit();
     gcount_priv.clear();gcount_priv.shrink_to_fit();
     mcount_priv.clear();mcount_priv.shrink_to_fit();
     machcount_priv.clear();machcount_priv.shrink_to_fit();
     bcount_priv.clear();bcount_priv.shrink_to_fit();
     ctoacount_priv.clear();ctoacount_priv.shrink_to_fit();
     btoacount_priv.clear();btoacount_priv.shrink_to_fit();
     lcdmcount_priv.clear();lcdmcount_priv.shrink_to_fit();
     tacount_priv.clear();tacount_priv.shrink_to_fit();
   }
   So.DONE();
   time (&end);
   So.message_screen("CT[secs]:", difftime(end,start));
   // We now write rbin,      
   ofstream profd;
   string file;
   file=this->params._Output_directory()+"separation_distribtion_"+namef+"_mass_run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<mcount[j]<<"\t"<<mcount[j]/gcount[j]<<endl;
   profd.close();
   file=this->params._Output_directory()+"separation_distribtion_"+namef+"_bias_run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<bcount[j]<<"\t"<<bcount[j]/gcount[j]<<endl;
   profd.close();
   file=this->params._Output_directory()+"separation_distribtion_"+namef+"_spin_run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<scount[j]<<"\t"<<scount[j]/gcount[j]<<endl;
   profd.close();
   file=this->params._Output_directory()+"separation_distribtion_"+namef+"_mach_run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<machcount[j]<<"\t"<<machcount[j]/gcount[j]<<endl;
   profd.close();
   file=this->params._Output_directory()+"separation_distribtion_"+namef+"_concentration_run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<ccount[j]<<"\t"<<ccount[j]/gcount[j]<<endl;
   profd.close();
}
   So.DONE();
}
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void PowerSpectrumF::GetSuperClusters(string prop){

    int NTHREADS = 1;
#ifdef _USE_OMP_
    NTHREADS=_NTHREADS_;
#endif
   vector<real_prec> DM_DEN_FIELD;
   vector<real_prec> TR_DEN_FIELD;
   if(prop=="_BIAS_"|| prop=="_RELATIVE_BIAS_")
    {
      DM_DEN_FIELD.resize(this->params._NGRID(),0);   
       this->File.read_array(this->params._Input_Directory_X()+this->params._Name_Catalog_X(),DM_DEN_FIELD);
      TR_DEN_FIELD.resize(this->params._NGRID(),0);   
      this->File.read_array(this->params._Input_Directory_Y()+this->params._Name_Catalog_Y(),TR_DEN_FIELD);
    }    
   this->params.set_input_type("density_grid");
   this->set_output_filenames();
   Catalog cat(this->params);
   cat.read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),0);
   this->set_params(cat.params); // params has been updated inside the catalog class.
   if(prop=="_BIAS_"|| prop=="_RELATIVE_BIAS_")
    {
       kvector_data.clear();
       kvector_data.shrink_to_fit();
       for(int i=0;i<this->params._d_Nnp_data();i++)
        kvector_data.push_back(this->params._d_kmin()+this->params._d_DeltaK_data()*(i+0.5));
       this->So.message_screen("Getting object-to-object bias");
       this->object_by_object_bias(cat.Halo,DM_DEN_FIELD,TR_DEN_FIELD);
       this->So.message_screen("\tMin Bias  =", cat.get_min("_BIAS_"));
       this->So.message_screen("\tMax Bias  =", cat.get_max("_BIAS_"));
       this->So.message_screen("\tMin RBias  =", cat.get_min("_RELATIVE_BIAS_"));
       this->So.message_screen("\tMax RBias  =", cat.get_max("_RELATIVE_BIAS_"));
    }
   cat.get_superclusters("over", prop); // Here we selec tracers with bias above the fourth quartile. This gives tracer_aux
   ULONG size_cat=cat.tracer_aux.size();
   int Nrbins=this->params._Nbins_cf();
   real_prec rmax=this->params._rmax_cf();
   real_prec rmin=this->params._rmin_cf();
   real_prec rmin_squared = rmin*rmin;
   real_prec rmax_squared = rmax*rmax;
   real_prec lrmin=log10(rmin);
   real_prec deltar=log10(rmax/rmin)/static_cast<real_prec>(Nrbins);
   // here we shuffle properties
   So.message_screen("  Ranking properties");
//   cat.Get_Ranked_Props(cat.tracer_aux,prop);
   So.DONE();
   //******************************************************
   //   Halo_ranked and tracer_aux have the same dimension, 
   // same coordinates in the i-label, only properties have been ranked
   //******************************************************
   //******************************************************
   //******************************************************
   //******************************************************
   // RANDOMIZAITION OF MARks
  int Nruns=100;// Numbner of shuffling marks to get errors
  gsl_rng_env_setup();
  const gsl_rng_type *Tn;
  gsl_rng_default_seed=1015;
  Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
  gsl_rng *rn = gsl_rng_alloc (Tn);

  vector<real_prec>rbias_tr(size_cat,0);
  for(ULONG i=0;i<size_cat;++i)
    rbias_tr[i]=cat.Halo_ranked[i].relative_bias; // keep this for JK
  vector<ULONG>aux_p;
  for(ULONG i=0;i<size_cat;++i)
    aux_p.push_back(i);
  for(int rr=0;rr<=Nruns;++rr)
    {
      vector<real_prec>gcount(Nrbins,0);
      vector<real_prec>mcount(Nrbins,0);
   if(rr>0)
   {
      So.message_screen("  Shuffling properties for run ", rr);
      gsl_ran_shuffle(rn,&aux_p[0], aux_p.size() ,sizeof(ULONG));

      for(ULONG i=0;i<size_cat;++i)
          if(prop=="_RELATIVE_BIAS_")
//            cat.Halo_ranked[i].relative_bias=cat.Halo_ranked[aux_p[i]].relative_bias;
            cat.Halo_ranked[i].relative_bias=rbias_tr[aux_p[i]];
         So.DONE();

    }

   real_prec mean_mark_m=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_mark_m)
#endif
   for(ULONG i=0;i<size_cat;++i)
     {
      if(prop=="_RELATIVE_BIAS_")
           mean_mark_m+=cat.Halo_ranked[i].relative_bias;
     }

   mean_mark_m/=static_cast<real_prec>(size_cat);
   So.message_screen("Computing marked statistics");
   time_t start;
   time_t end;
   time (&start);
#pragma omp parallel
   {
     vector<real_prec> gcount_priv(Nrbins,0);
     vector<real_prec> mcount_priv(Nrbins,0);
#pragma omp for collapse(2)
     for(ULONG i=0;i<size_cat;++i)
       for(ULONG j=i+1;j<size_cat;++j)
         {
           real_prec dx=cat.tracer_aux[i].coord1-cat.tracer_aux[j].coord1;
           real_prec dy=cat.tracer_aux[i].coord2-cat.tracer_aux[j].coord2;
           real_prec dz=cat.tracer_aux[i].coord3-cat.tracer_aux[j].coord3;
           real_prec dist_squared=dx*dx+dy*dy+dz*dz;
           if (dist_squared>=rmin_squared && dist_squared<rmax_squared)
           {
               real_prec ldist=0.5*log10(dist_squared); // this is log (sqrt(dist_squared))
               ULONG ibin=get_bin(ldist,lrmin,Nrbins,deltar,false);
               if(ibin<Nrbins)
               {
                   gcount_priv[ibin]++;  // marking with property=mass
                   if(prop=="_RELATIVE_BIAS_")
                     mcount_priv[ibin]+=cat.Halo_ranked[i].relative_bias*cat.Halo_ranked[j].relative_bias;
               }
            }
         }
#pragma omp critical
     {
       for(ULONG i=0;i<Nrbins;++i)
        {
           gcount[i]+= gcount_priv[i];
           mcount[i]+= mcount_priv[i]/pow(mean_mark_m,2);
        }
     }
     mcount_priv.clear();mcount_priv.shrink_to_fit();
   }
   So.DONE();
   time (&end);
   So.message_screen("CT[secs]:", difftime(end,start));
   // We now write rbin,
   ofstream profd;
   string file;
   file=this->params._Output_directory()+"separation_distribtion_superk"+prop+"run"+to_string(rr)+".txt";
   profd.open(file);
   for(ULONG j=0;j<Nrbins;j++)
     profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<mcount[j]<<"\t"<<mcount[j]/gcount[j]<<endl;
   profd.close();
}
   //******************************************************
   //******************************************************
   //******************************************************
   //******************************************************
   //   jacknife
      So.message_screen("Getting Jackknife estimates");
#pragma omp parallel for
   for(ULONG i=0;i<size_cat;++i) 
        cat.tracer_aux[i].GridID_n = grid_ID(&cat.box_JK, cat.tracer_aux[i].coord1,cat.tracer_aux[i].coord2,cat.tracer_aux[i].coord3);

   for(ULONG rr=0;rr<=this->params.d_NGRID_JK();++rr) // Loop over each of the sub-samples JK. The 0th item is the full sample
     {

       vector<real_prec>gcount(Nrbins,0);
       vector<real_prec>mcount(Nrbins,0);

       
       vector<real_prec>xc;
       vector<real_prec>yc;
       vector<real_prec>zc;
       vector<real_prec>mc;
       real_prec mean_mark_m=0;

    if(rr==0)
    {
       for(ULONG i=0;i<size_cat;++i)// Loop over the catalog
       {
         xc.push_back(cat.tracer_aux[i].coord1);
         yc.push_back(cat.tracer_aux[i].coord2);
         zc.push_back(cat.tracer_aux[i].coord3);
         if(prop=="_RELATIVE_BIAS_")
         {
//             mc.push_back(cat.Halo_ranked[i].relative_bias);
 //           mean_mark_m+=cat.Halo_ranked[i].relative_bias;
            mc.push_back(rbias_tr[i]);
           mean_mark_m+=rbias_tr[i];
         }
       }
    }
      else{
       for(ULONG i=0;i<size_cat;++i)// Loop over the catalog
         if(cat.tracer_aux[i].GridID_n != rr)
           {
             xc.push_back(cat.tracer_aux[i].coord1);
             yc.push_back(cat.tracer_aux[i].coord2);
             zc.push_back(cat.tracer_aux[i].coord3);
             if(prop=="_RELATIVE_BIAS_")
             {
//                mc.push_back(cat.Halo_ranked[i].relative_bias);
//                mean_mark_m+=cat.Halo_ranked[i].relative_bias;
                mc.push_back(rbias_tr[i]);
                mean_mark_m+=rbias_tr[i];
             }
          }
       }

       ULONG nsize_cat=xc.size();    
       So.message_screen("Number of tracers before removing subsample:", size_cat);
       So.message_screen("Number of tracers after removing subsample:", nsize_cat);
       mean_mark_m/=static_cast<real_prec>(nsize_cat);

       So.message_screen("Computing marked statistics for JK=", rr);
       time_t start;
       time_t end;
       time (&start);
#pragma omp parallel 
       {
     vector<real_prec> gcount_priv(Nrbins,0);
     vector<real_prec> mcount_priv(Nrbins,0);
#pragma omp for collapse(2)
     for(ULONG i=0;i<nsize_cat;++i)
       for(ULONG j=i+1;j<nsize_cat;++j)
         {
           real_prec dx=xc[i]-xc[j];
           real_prec dy=yc[i]-yc[j];
           real_prec dz=zc[i]-zc[j];
           real_prec dist_squared=dx*dx+dy*dy+dz*dz;
           if (dist_squared>=rmin_squared && dist_squared<rmax_squared)
         {
           real_prec ldist=0.5*log10(dist_squared); // this is log (sqrt(dist_squared))
           ULONG ibin=get_bin(ldist,lrmin,Nrbins,deltar,false);
           if(ibin<Nrbins) 
             {
               gcount_priv[ibin]++;  // marking with property=mass
               mcount_priv[ibin]+=mc[i]*mc[j];
             }
         }
         }
#pragma omp critical
     {
       for(ULONG i=0;i<Nrbins;++i)
         {
           gcount[i]+= gcount_priv[i];
           mcount[i]+= mcount_priv[i]/pow(mean_mark_m,2);
         }
     }
     gcount_priv.clear();gcount_priv.shrink_to_fit();
     mcount_priv.clear();mcount_priv.shrink_to_fit();
       }//Closes parallel region
       
       So.DONE();
       time (&end);
       So.message_screen("CT[secs]:", difftime(end,start));
       
       // We now write rbin,      
       ofstream profd;
       string file;
       
       file=this->params._Output_directory()+"separation_distribtion_superk"+prop+"jkrun"+to_string(rr)+".txt";
       profd.open(file);
       for(ULONG j=0;j<Nrbins;j++)
         profd<<rmin*pow(10,(j+0.5)*deltar)<<"\t"<<gcount[j]<<"\t"<<mcount[j]<<"\t"<<mcount[j]/gcount[j]<<endl;
       profd.close();

     }// End loop over the JK regions
  So.DONE();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
