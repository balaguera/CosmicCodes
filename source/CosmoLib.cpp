////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @file CosmoLib.cpp
 *
 * @brief This file contains methods of the class CosmoLib
 * @details The class Cosmolib generates estaimates of a number of cosmlogical observables
 * @author Andres Balaguera Antolinez
 * @date 2007-2023
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../headers/CosmoLib.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void message(string mess){
  std::cout<<BOLDRED<<mess<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void welcome_message(){
 time_t rawtime;
 time (&rawtime);
 cout<<CYAN<<"***********************************************************************************"<<endl;
 cout<<"COSMOLIB: Some cosmology-related numbers "<<endl;
 cout<<"VERSION 1.2"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void welcome_read_chains(){
  time_t rawtime;
 time (&rawtime);
 cout<<"***********************************************************************************"<<endl;
 cout<<"Reading MCMChains"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void message_screen(string s, real_prec v, string units)
{
  cout<<"Derived quantity"<<endl;
  cout<<CYAN;
  cout<<s<<" = "<<v<<" "<<units<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void message_screen(string s, real_prec v)
{
  cout<<"Derived quantity"<<endl;
  cout<<CYAN<<s<<" = "<<v<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void message_screen_ini(string s, real_prec v)
{
  cout<<"Input parameter"<<endl;
  cout<<CYAN<<s<<" = "<<v<<RESET<<endl;
}
///////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void welcome_chains(){
 time_t rawtime;
 time (&rawtime);
 cout<<"***********************************************************************************"<<endl;
 cout<<"MCMChains"<<endl;
 cout<<"\t"<<endl;
 cout<<"Starting time "<<ctime (&rawtime);
 cout<<"***********************************************************************************"<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void enter(string text){
  cout<<"\t"<<endl;
  cout<<"Going to:  "<<text<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void ending_message(){
  cout<<"    "<<endl;
  time_t rawtime;
  time ( &rawtime );
  cout<<"Ending time:"<<ctime (&rawtime)<<endl;
  cout<<"***********************************************************************************"<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void CosmoLib::set_par_file(string par_file)
{
  Params param (par_file);
  this->params=param;
  Statistics Csa(this->params._M_min_mf(), params._M_max_mf(), params._npoints_mf());
  this->Cs=Csa;
  PowerSpectrum Psa(s_cosmo_par,params._kmin_integration() ,params._kmax_integration() ,params._npoints_dp_k(),10, params._M_min_mf(),params._M_max_mf(),params._npoints_mf());
  this->Ps=Psa;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void comp_time(time_t start, long full, long step){
  real_prec fraction=100.0*((real_prec)(step+1))/((real_prec)full);
  time_t end;
  time(&end);
  real_prec lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs \r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60.<<"  mins \r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600.<<"  hrs \r";cout.flush();
  if (fraction==25) cout <<"\r  ..25% completed \r";cout.flush();
  if (fraction==50) cout <<"\r  ..50% completed \r";cout.flush();
  if (fraction==75) cout <<"\r  ..75% completed \r";cout.flush();
  if (fraction==100) cout<<"\r ..100% completed \r";cout.flush();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void comp_time_MH(time_t start, long full, long step, real_prec H){
  std::cout<<RED;
  real_prec fraction=100.0*((real_prec)(step+1))/((real_prec)full);
  time_t end;
  time(&end);
  real_prec lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600<<"  secs.  H(q,p) = "<<H<<"\r";cout.flush();
  std::cout<<RESET;
  if (fraction==25) cout <<"\r  ..25% completed \r";cout.flush();
  if (fraction==50) cout <<"\r  ..50% completed \r";cout.flush();
  if (fraction==75) cout <<"\r  ..75% completed \r";cout.flush();
  if (fraction==100) cout<<"\r ..100% completed \r";cout.flush();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void CosmoLib::show_cosmo_struct()
{
#ifdef _USE_UNITSIM_COSMOLOGY_
  So.message_screen("Using cosmological parameters from UNITsim (Planck 16)");
#endif
#ifdef _USE_SLICS_COSMOLOGY_
  So.message_screen("Using cosmological parameters from SLICS ");
#endif
#ifdef _USE_ABACUS_COSMOLOGY_
  So.message_screen("Using cosmological parameters from ABACUS simulation ");
#endif
  Cf.check_cosmo_pars();
  So.message_screen("Omega Matter =", this->s_cosmo_par.Om_matter);
  So.message_screen("Omega CDM =", this->s_cosmo_par.Om_cdm);
  So.message_screen("Omega Vacuum =", this->s_cosmo_par.Om_vac);
  So.message_screen("Omega Baryons =", this->s_cosmo_par.Om_baryons);
  So.message_screen("Omega Radiation =", this->s_cosmo_par.Om_radiation);
  So.message_screen("Omega Curvature =", this->s_cosmo_par.Om_k);
  So.message_screen("DE eos =", this->s_cosmo_par.wde_eos);
  So.message_screen("Hubble parameter h= ", this->s_cosmo_par.hubble);
  So.message_screen("Spectral index =", this->s_cosmo_par.spectral_index);
  So.message_screen("Sigma 8 =", this->s_cosmo_par.sigma8);
  So.message_screen("Primordial amplitude =", this->s_cosmo_par.A_s);
  So.message_screen("Running index =", this->s_cosmo_par.alpha_s);
  So.message_screen("Use wiggles =", this->s_cosmo_par.use_wiggles);
  So.message_screen("Mean CMB temperature =", this->s_cosmo_par.Tcmb);
  So.message_screen("Delta Spherical overdensity =", this->s_cosmo_par.Delta_SO);
  So.message_screen("f_baryon =", this->s_cosmo_par.f_baryon);
  this->Ps.set_cosmo_pars(this->s_cosmo_par);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Derive the abundance and dNdz as a function of redshift given some HOD and halo mass function.
// Also, we must give the minimum mass as a function of redhsift.
void CosmoLib::get_dndz_gal(){
    this->show_cosmo_struct();
    int NTHREADS=1;
#ifdef _USE_OMP_
    NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
    // ----------------------------------------------
    // Open file with minimum halo mass
    string mmin_f="/home/andres/data/Numerics/JPAS/Mmin_Watson_comp.dat";
    vector<real_prec> nb;
    vector<gsl_real> z_v;
    vector<gsl_real> Mmin_v;
    ULONG nzl=this->File.read_file(mmin_f,nb, 1);
    ULONG ncols_nb=(static_cast<ULONG>(nb.size()/nzl));
    z_v.resize(nzl,0);
    Mmin_v.resize(nzl,0);
    for(ULONG i=0;i<nzl;++i)
     {
       z_v[i]=nb[i*ncols_nb];  // redshift, assumed to be in the first column of the file
       Mmin_v[i]=2.0*nb[1+i*ncols_nb]; // mean number density
     }
    nb.clear(); nb.shrink_to_fit();
    So.message_screen("Get ready for interpolation of MMin(z)");
    gsl_spline *spline_mm = gsl_spline_alloc (gsl_interp_linear,nzl);
    gsl_spline_init (spline_mm, &z_v[0], &Mmin_v[0], nzl);
    gsl_interp_accel *spline_acc_mm = gsl_interp_accel_alloc ();
    So.DONE();
    //This has to be modified. We need to find:
        // 1 ) Halo - Central stellar mass relation
       //  2 ) For each redshift and minium app magnitude, the corresponding Lmin=Lmin(z)
       //  3 ) The corresponding stellar mass for such minimum luminosity
       //  4 ) The correspoinding halo mass for that minimum stellar mass
    // This is the mass we must use.
    // ----------------------------------------------
    // Feed Cosmology with cosmological parameters
    this->Cf.set_cosmo_pars(s_cosmo_par);
    // Compute normalization of power spectrum (at z=0)
    real_prec pk_normalization;
    this->Ps.normalization(pk_normalization);
    this->s_cosmo_par.pk_normalization=pk_normalization;
    Ps.set_cosmo_pars(s_cosmo_par);// update parameters
    So.message_screen("Get ready for integration in k-space");
    this->Cs.compute_int_table_wavenumber(this->params._kmin_ps(),this->params._kmax_ps());
    So.DONE();
    // Get gwoth factor
    real_prec growth_factor_0=Cf.growth_factor(0);
    // ----------------------------------------------
    // BUild redshift vector
    real_prec delta_redshift=(this->params._redshift_max()-this->params._redshift_min())/static_cast<real_prec>(this->params._nbins_redshift());
    // We open a loop over redshift and compute cosmological functions at each step.
    vector<real_prec> redb(params._nbins_redshift(),0);
    vector<real_prec> Vol(params._nbins_redshift(),0);
    vector<real_prec> nbar(params._nbins_redshift(),0);
    vector<real_prec> dndz(params._nbins_redshift(),0);
    // ----------------------------------------------
    this->v_mass_function.resize(params._npoints_mf(),0);
    this->v_mass.resize(params._npoints_mf(),0);
    for(int i=0;i<this->params._npoints_mf();++i)
        this->v_mass[i]=static_cast<gsl_real>(log10(Mmin_v[0])+i*(log10(params._M_max_effective())-log10(Mmin_v[0]))/static_cast<double>(v_mass.size()));
    this->s_cosmo_par.v_mass=this->v_mass;// Update
    // ----------------------------------------------
    for(int iz=0;iz<this->params._nbins_redshift();++iz)
    {
        redb[iz]=this->params._redshift_min()+iz*delta_redshift;
        real_prec cdist=this->Cf.comoving_distance(redb[iz]);
        Vol[iz]=this->params._area_survey()*cdist*cdist*cdist;
    }
    // ----------------------------------------------
    So.message_screen("Starting loop over redshift");
    // This sets the global limits of integration. Inside the redshift loop
    // we use Mmin in the calculation of nbar gal
    this->Cs.compute_int_table_mass(params._M_min_effective(), params._M_max_effective());
    ofstream fo;
    string fof=this->params._Output_directory()+"abundance_jpas_theoretical.txt";
    fo.open(fof.c_str());
#ifdef _USE_OMP_
//#pragma omp parallel for
#endif
    for(int iz=0;iz<this->params._nbins_redshift();++iz)
    {
       s_cosmo_par.growth_factor=Cf.growth_factor(redb[iz])/growth_factor_0;
       real_prec critical_density=Cf.critical_overdensity(redb[iz]);
       s_cosmo_par.critical_density=critical_density;
       real_prec density_contrast_top_hat=Cf.density_contrast_top_hat(redb[iz]);
       s_cosmo_par.density_contrast_top_hat=density_contrast_top_hat;
       real_prec mean_matter_density=Cf.mean_matter_density(redb[iz]);
       s_cosmo_par.mean_matter_density=mean_matter_density;
       this->Cs.set_cosmo_pars(s_cosmo_par); // Update params in Statistics (and in its Cosmology instance)
       //Define minimum mass as a function of redshift
       // Get mass function
       for(int im=0;im<this->params._npoints_mf();++im)
           v_mass_function[im]=static_cast<gsl_real>(this->Cs.mass_function_D(static_cast<real_prec>(this->v_mass[im]),redb[iz]));
       this->s_cosmo_par.v_mass_function=v_mass_function;
       this->Cs.set_cosmo_pars(s_cosmo_par); // Update params in Statistics
       real_prec min_mass=gsl_spline_eval(spline_mm, redb[iz], spline_acc_mm); // Minimumn mass at this redshift
       nbar[iz]=this->Cs.mean_galaxy_number_density(redb[iz],min_mass);// get nbar_gal
       dndz[iz]=nbar[iz]*Vol[iz];
       fo<<redb[iz]<<"  "<<nbar[iz]<<"  "<<dndz[iz]<<endl;
       cout<<redb[iz]<<"  "<<nbar[iz]<<"  "<<dndz[iz]<<endl;
   }
    fo.close();
    So.DONE();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CosmoLib::get_cosmolib(){
  welcome_message();
  int NTHREADS=1;
  this->show_cosmo_struct();
#ifdef _USE_OMP_
  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  real_prec redshift=0;
#ifdef _GET_EFF_BIAS_REDSHIFT_
`  vector<real_prec> redb(params._nbins_redshift(),0);
  vector<real_prec> effb(params._nbins_redshift(),0);
//#pragma omp parallel for
  for(int iz=0;iz<params._nbins_redshift();++iz)
      {
         redshift=params._redshift_min() + iz*(params._redshift_max()-params._redshift_min())/(static_cast<real_prec>(params._nbins_redshift()));
         redb[iz]=redshift;
#else
  redshift=this->params._redshift();
#endif
  ofstream hout;
  hout.open("hm_check.log");
  // *********************************************************
  this->Cf.set_cosmo_pars(this->s_cosmo_par);
  real_prec critical_density=Cf.critical_overdensity(redshift);
  real_prec density_contrast_top_hat=Cf.density_contrast_top_hat(redshift);
  if(this->params._Get_SO_from_BN()==true)
    {
      this->s_cosmo_par.Delta_SO=density_contrast_top_hat;
      So.message_screen("Using density contrast top hat as SO");
  }
  this->Cf.set_cosmo_pars(this->s_cosmo_par);// update
  this->Ps.set_cosmo_pars(this->s_cosmo_par);// update
  real_prec Hubble_function=Cf.Hubble_function(redshift);
  real_prec comoving_distance=this->Cf.comoving_distance(redshift);
  real_prec comoving_angular_diameter_distance=Cf.comoving_angular_diameter_distance(redshift);
  real_prec mean_matter_density=Cf.mean_matter_density(redshift);
  real_prec age_universe=Cf.age_universe(redshift);
  real_prec comoving_sound_horizon=Cf.comoving_sound_horizon(this->Cf.drag_redshift());
  real_prec growth_factor=Cf.growth_factor(redshift);
  real_prec growth_factor_z0=Cf.growth_factor(0.0);
  real_prec growth_index=Cf.growth_index(redshift);
  real_prec halo_dynamical_time=Cf.halo_dynamical_time(redshift);
  real_prec omega_matter=Cf.omega_matter(redshift);
  real_prec Distance_Modulus=Cf.Distance_Modulus(redshift);
  real_prec pk_normalization=1.;
  // Normalize power spectrum
  if(false==this->params._use_file_power())
    pk_normalization=Ps.normalization();
  this->s_cosmo_par.pk_normalization=pk_normalization;
  this->Ps.set_cosmo_pars(this->s_cosmo_par);// update
  this->Cf.set_cosmo_pars(this->s_cosmo_par);
  So.message_screen("Redshift z =", this->params._redshift());
  So.message_screen("Derived quantities");
  So.message_screen("Omega matter at this redshift =", omega_matter);
  So.message_screen("Hubble parameter at this redshift =", Hubble_function);
  So.message_screen("Mean matter density =", mean_matter_density, "(Ms/h)/(Mpc/h)^(-3)");
  So.message_screen("Comoving distance to current redshift =", comoving_distance, "Mpc/h");
  So.message_screen("Comoving angular diameter distance to current resdhift =", comoving_angular_diameter_distance, "Mpc/h");
  So.message_screen("Drag reedshift =", this->Cf.drag_redshift());
  So.message_screen("Comoving sound horizon @ z_drag =", comoving_sound_horizon, "Mpc/h");
  So.message_screen("Comoving sound horizon (fit) =", Cf.sound_horizon(), "Mpc");
  So.message_screen("Age of the Universe at current redshift =", age_universe/1e9, "Gyr/h");
  So.message_screen("Distance Modulus =", Distance_Modulus);
  So.message_screen("Growth D(z) at redshift z =", growth_factor);
  So.message_screen("D(z=0) =", growth_factor_z0);
  So.message_screen("D(z=0)/D(z) =", pow(growth_factor_z0/growth_factor,1));
  So.message_screen("g(z)=D(z)/D(z=0) =", growth_factor/growth_factor_z0);
  So.message_screen("Growth index f(z) = Om(z)⁰.⁵⁵ =", growth_index);
  So.message_screen("Halo-dynamical time =", halo_dynamical_time/1e9, "Gyr/h");
  So.message_screen("Critical overdensity linearly extrapolated =", critical_density);
  So.message_screen("Top-hat density contrast at virial =",density_contrast_top_hat);
  // ***********************************************************************************************
  // Compute tables to perform integrals using GL weights
  this->Cs.compute_int_table_wavenumber(this->params._kmin_ps(),this->params._kmax_ps());
  // ***********************************************************************************************
  this->Cs.compute_int_table_mass(this->params._M_min_effective(), this->params._M_max_effective());
  if(true==this->params._use_file_power())
  {
    vector<real_prec>prop;
    ULONG NLINES= File.read_file(this->params._file_power(),prop,1);
    ULONG NCOLS=(static_cast<ULONG>(prop.size()/NLINES));
    s_cosmo_par.kvector_external.resize(NLINES,0);
    s_cosmo_par.power_external.resize(NLINES,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<NLINES;++i)
    {
      s_cosmo_par.kvector_external[i]=static_cast<gsl_real>(prop[NCOLS*i]);
      s_cosmo_par.power_external[i]=static_cast<gsl_real>(prop[1+NCOLS*i]);
    }
    prop.clear(); prop.shrink_to_fit();
    Cs.use_external_power=this->params._use_file_power();
    s_cosmo_par.use_external_power=this->params._use_file_power();
    this->params.set_npoints_ps(NLINES);
    this->s_cosmo_par.kmax_int=s_cosmo_par.kvector_external[NLINES-1];
    this->s_cosmo_par.kmin_int=s_cosmo_par.kvector_external[0];
    this->params.set_kmax_ps(this->s_cosmo_par.kmax_int);
    this->params.set_kmin_ps(this->s_cosmo_par.kmin_int);
    So.message_screen("Minimum k-in external file = ", this->params._kmin_ps());
    So.message_screen("Maximum k-in external file = ", this->params._kmax_ps());
  }
  else
    So.message_screen("Normalization of matter P(k,z) =", pk_normalization);
  growth_factor/=Cf.growth_factor(0);
  // Normalize the growth factor to compute the processed linear matter power spectrum
  // Aca reacomodo algunos de estos factores en la estructura grande
  this->s_cosmo_par.critical_density=critical_density;
  this->s_cosmo_par.density_contrast_top_hat=density_contrast_top_hat;
  this->s_cosmo_par.mean_matter_density=mean_matter_density;
  this->s_cosmo_par.growth_factor=growth_factor;
  // ***********************************************************************************************
  this->Cf.set_cosmo_pars(this->s_cosmo_par);
  this->Ps.set_cosmo_pars(this->s_cosmo_par);// update
  this->Cs.set_cosmo_pars(this->s_cosmo_par);
  // ***********************************************************************************************
//  real_prec sigma_from_As=Cs.As2sigma8(&s_cosmo_par);
  // ***********************************************************************************************
  // Generating mass function
  this->v_mass.resize(params._npoints_mf(),0);
  this->v_sigma_mass.resize(params._npoints_mf(),0);
  this->v_mass_function.resize(params._npoints_mf(),0);
  this->v_halo_mass_bias.resize(params._npoints_mf(),0);
  time_t start;
  time (&start);
 #ifdef _GET_ONLY_MASS_FUNCTION_
  So.message_screen("Computing Sigma(M)");
  if(params._scale_mf()=="linear")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i<v_mass.size();i++)
     this->v_mass[i]=static_cast<gsl_real>(log10(params._M_min_mf()+i*(params._M_max_mf()-params._M_min_mf())/static_cast<double>(v_mass.size())));
  else if(params._scale_mf()=="log")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i<v_mass.size();i++)
      this->v_mass[i]=static_cast<gsl_real>(log10(params._M_min_mf())+i*(log10(params._M_max_mf())-log10(params._M_min_mf()))/static_cast<double>(v_mass.size()));
  for(int i=0; i<v_mass.size();i++) // DO not parallelize this loop.
    this->v_sigma_mass[i]=static_cast<gsl_real>(Cs.sigma_masa(v_mass[i],redshift));
  So.DONE();
  s_cosmo_par.v_mass=this->v_mass;
  s_cosmo_par.v_sigma_mass=this->v_sigma_mass;
  this->Cs.set_cosmo_pars(s_cosmo_par);
  this->Cs.set_spline_vars_sigma();
  So.message_screen("Computing halo mass-function");
  for(int i=0; i<v_mass.size();i++)
   {
#ifdef TIME
    comp_time(start, v_mass.size(), i);
#endif    
    this->v_mass_function[i]=static_cast<gsl_real>(Cs.mass_function_D(static_cast<real_prec>(this->v_mass[i]), redshift));
    this->v_halo_mass_bias[i]=static_cast<gsl_real>(Cs.bias(static_cast<real_prec>(this->v_mass[i]),redshift,&s_cosmo_par));
  }
  So.DONE();
#ifndef _GET_EFF_BIAS_REDSHIFT_
  File.write_to_file(this->params._mass_function_output_file(),v_mass,v_mass_function,v_sigma_mass);
  File.write_to_file(this->params._halo_mass_bias_output_file(),v_mass,v_halo_mass_bias);
#endif
  // ALlcate new comp[uted vectors in a structure to be interpolated later
  s_cosmo_par.M_max_mf=params._M_max_mf();
  s_cosmo_par.M_min_mf=params._M_min_mf();
  s_cosmo_par.v_mass=this->v_mass;
  s_cosmo_par.v_mass_function=this->v_mass_function;
  s_cosmo_par.v_halo_mass_bias=this->v_halo_mass_bias;
  //this->Cs.set_spline_vars_mfunction();// Ready for interpolations, buut with problems, as it uses this in static memebr functiuons, impossible
  // Once this has been created, we can compute integrals with respect to the mass.
  // Let us go for effective halo_mass bias as a function of a minimum mass
  v_effective_halo_mass_bias.resize(this->params._npoints_mf(),0);
  v_effective_halo_mean_number_density.resize(this->params._npoints_mf(),0);
  s_cosmo_par.M_min_effective=params._M_min_effective();
  s_cosmo_par.M_max_effective=params._M_max_effective();
  So.message_screen("Computing effective halo mass-bias");
//#pragma omp parallel for
  for(int i=0; i<v_mass.size();i++){
#ifdef TIME
    comp_time(start, v_mass.size(), i);
#endif
    this->v_effective_halo_mass_bias[i]=Cs.effective_halo_mass_bias(v_mass[i],redshift,&s_cosmo_par);
    this->v_effective_halo_mean_number_density[i]=Cs.effective_halo_mean_number_density(v_mass[i],redshift,&s_cosmo_par);
  }  
  So.DONE();
#ifndef _GET_EFF_BIAS_REDSHIFT_
  File.write_to_file(params._effective_halo_mass_bias_output_file(),v_mass, v_effective_halo_mass_bias);
  File.write_to_file(params._effective_halo_mean_number_density_output_file(),v_mass,v_effective_halo_mean_number_density);
#endif
  // ***********************************************************************************************
  // ***********************************************************************************************
  real_prec eff_bias=Cs.effective_halo_mass_bias( log10(params._M_min_effective()),redshift,&s_cosmo_par);
  So.message_screen("Effective halo-mass bias at the resolution (min) mhalo mass", eff_bias);
  real_prec kaiser = 1.+(2./3.)*(growth_index/eff_bias)+(1./5.)*pow(growth_index/eff_bias,2);
  So.message_screen("Kaiser Factor for RSD ", kaiser);
#endif  // end of only mass function
#ifdef _GET_EFF_BIAS_REDSHIFT_
  effb[iz]=eff_bias;
}
ofstream bout; bout.open("../MOCKS_JPAS/effective_bias_redshift.txt");
for(int iz=0;iz<params._nbins_redshift();++iz)
        bout<<redb[iz]<<"  "<<effb[iz]<<endl;
bout.close();
#endif
#ifdef _GET_POWER_SPECTRUM_
  // ***********************************************************************************************
  // Calculamos las escalas que definen que es linal approx, usando sigma = 1
  real_prec Mnl, rnl, sigman, knl;
 #ifdef _GET_NL_POWER_SPECTRUM_
  s_cosmo_par.aux_var1=params._redshift(); //ESTO LO PUEDO HACER MEJOR CON LAS NON SCALES DEL HALO FIT, VER CODIO cl_functions
  this->Cs.set_cosmo_pars(s_cosmo_par);
  this->Cs.non_linear_scales(knl,Mnl,rnl,sigman);
  s_cosmo_par.Mnl=Mnl;
  s_cosmo_par.knl=knl;
  s_cosmo_par.rnl=rnl;
  So.message_screen("Non linear mass scale Mnl",Mnl,"Ms/h");
  So.message_screen("Non linear wave number knl",knl,"h/Mpc");
  So.message_screen("Non linear scales rnl",rnl,"Mpc/h");
  So.message_screen("\tSanity check: sigma(Mnl)",sigman);
  hout<<"Non linear scales : Mass = "<<Mnl<<" Ms/h "<<endl;
  hout<<"Non linear scales : k    = "<<knl<<" h/Mpc "<<endl;
  hout<<"Non linear scales : r    = "<<rnl<<" Mpc/h "<<endl;
  hout<<"Sanity check: sigma(Mnl) = "<<sigman<<endl;
#endif
  // ***********************************************************************************************
  // NOW WE CAN COMPUTE HALO FIT
  // ***********************************************************************************************
  v_k_ps.resize(params._npoints_ps(),0);
  if(params._scale_ps()=="linear")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i<v_k_ps.size();i++)
        this->v_k_ps[i]=(params._kmin_ps()+i*(params._kmax_ps()-params._kmin_ps())/v_k_ps.size());
  else if(params._scale_ps()=="log")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0; i<v_k_ps.size();i++)
       this->v_k_ps[i]=pow(10,(log10(params._kmin_ps())+i*(log10(params._kmax_ps())-log10(params._kmin_ps()))/static_cast<real_prec>(v_k_ps.size())));
  this->v_l_power_spectrum.resize(params._npoints_ps(),0);
#ifdef _GET_NL_POWER_SPECTRUM_
  this->v_nl_power_spectrum.resize(params._npoints_ps(),0);
#ifdef _GET_NL_PT_POWER_SPECTRUM_
  v_nl_power_spectrum_pt.resize(params._npoints_ps(),0);
#endif
  if(params._scale_ps()=="linear")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<v_k_ps.size();++i)
     this->v_k_ps[i]=(0.5*params._kmin_ps()+i*(2.*params._kmax_ps()-0.5*params._kmin_ps())/static_cast<double>(v_k_ps.size()));
  else if(params._scale_ps()=="log")
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<v_k_ps.size();++i)
      this->v_k_ps[i]=pow(10,(log10(0.5*params._kmin_ps())+i*(log10(2.*params._kmax_ps())-log10(0.5*params._kmin_ps()))/static_cast<double>(v_k_ps.size())));
  this->s_cosmo_par.v_k_ps=v_k_ps;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
  for(int i=0; i<v_k_ps.size();i++)
      v_l_power_spectrum[i]=Ps.Linear_Matter_Power_Spectrum(this->v_k_ps[i]);
  File.write_to_file(params._linear_matter_ps_output_file(),this->v_k_ps,this->v_l_power_spectrum);
  s_cosmo_par.v_lin_power_spectrum=v_l_power_spectrum;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
  vector<real_prec> rrs(100,0);
  vector<real_prec> sss(100,0);
  real_prec knl_hf, rnl_hf;
  Ps.nl_scales_halo_fit( knl_hf, rnl_hf, rrs, sss,true);
  s_cosmo_par.knl_hf=knl_hf;
  s_cosmo_par.rnl_hf=rnl_hf;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
  real_prec kstar;
  Ps.kstar_integral(kstar);
  s_cosmo_par.kstar=kstar;
  Ps.set_cosmo_pars(s_cosmo_par); // update
  hout<<"Non linear scales halo fit: k       = "<<knl_hf<<" h/Mpc "<<endl;
  hout<<"Non linear scales halo fit: r       = "<<rnl_hf<<" Mpc/h "<<endl;
  So.message_screen("Non linear scales halo fit k_nl =",knl_hf," h/Mpc");
  So.message_screen("Non linear scales halo fit: r_nl =", rnl_hf," Mpc/h");
  //  hout<<"Sigma8 from As                      = "<<sigma_from_As<<endl;
  So.message_screen("k* = ",kstar," h/Mpc");
  time (&start);
  So.message_screen("Defining bins in k-space");
  So.message_screen("Computing Non linear matter power spectrum");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0; i<v_k_ps.size();i++){
      v_nl_power_spectrum[i]=Ps.Non_Linear_Matter_Power_Spectrum_Halo_Fit(v_k_ps[i]);
#ifdef _GET_NL_PT_POWER_SPECTRUM_
      v_nl_power_spectrum_pt[i]=Ps.Non_Linear_Matter_Power_Spectrum_PT(v_k_ps[i]);
#endif
  #ifdef TIME
    comp_time(start, v_k_ps.size(), i);
#endif
       }
#endif
#ifdef _GET_NL_POWER_SPECTRUM_
  File.write_to_file(params._non_linear_matter_ps_halo_fit_output_file(),this->v_k_ps,this->v_nl_power_spectrum);
#ifdef _GET_NL_PT_POWER_SPECTRUM_
  File.write_to_file(params._non_linear_matter_ps_pt_output_file(),v_k_ps,v_nl_power_spectrum_pt);
#endif
  s_cosmo_par.v_nl_power_spectrum=v_nl_power_spectrum;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
#endif
  s_cosmo_par.v_k_ps=this->v_k_ps;
  Ps.set_cosmo_pars(s_cosmo_par); // updates
#endif
  // ***********************************************************************************************
  // ***********************************************************************************************
#ifdef _GET_CORRELATION_FUNCTION_
  // ***********************************************************************************************
  // CORRELATION FUNCTION
  CorrelationFunctionTH SCf(1000, 1.);  // Check NumericalMethods to understand these inputs
  
  v_r_cf.resize(params._npoints_cf(),0);
  if(params._compute_output_linear_correlation_function())
    {
      So.message_screen("Computing matter correlation function");
      v_nl_correlation_function.resize(params._npoints_cf(),0);
      v_l_correlation_function.resize(params._npoints_cf(),0);
      time (&start);
      
      for(int i=0;i<v_r_cf.size();i++){
	//      comp_time(start, v_r_cf.size(), i);
	if(params._scale_cf()=="linear")v_r_cf[i]=(params._rmin_cf()+i*(params._rmax_cf()-params._rmin_cf())/v_r_cf.size());
	else if(params._scale_cf()=="log")v_r_cf[i]=pow(10,(log10(params._rmin_cf())+i*(log10(params._rmax_cf())-log10(params._rmin_cf()))/v_r_cf.size()));
	v_l_correlation_function[i]=SCf.Linear_Matter_Correlation_Function(&s_cosmo_par,v_r_cf[i]);
	v_nl_correlation_function[i]=SCf.Non_Linear_Matter_Correlation_Function_Halo_Fit(&s_cosmo_par,v_r_cf[i]);
      }
      File.write_to_file(params._linear_matter_cf_output_file(),v_r_cf,v_l_correlation_function);
      File.write_to_file(params._non_linear_matter_cf_halo_fit_output_file(),v_r_cf,v_nl_correlation_function);
      So.DONE();
  }
#endif
  // ***********************************************************************************************
  // ***********************************************************************************************
  // ***********************************************************************************************
  // *********************************************************************************************** 
#ifdef _GET_DENSITY_PROFILES_
  So.message_screen("Density profiles in configuration space");
  DensityProfiles Dp;
  v_r_dp.resize(params._npoints_dp_r(),0);
  v_density_profile_r.resize(params._npoints_dp_r(),0);
  for(int i=0;i<v_r_dp.size();i++){
#ifdef TIME
    comp_time(start, v_r_dp.size(), i);
#endif
    if(params._scale_dp_r()=="linear")v_r_dp[i]=(params._rmin_dp()+i*(params._rmax_dp()-params._rmin_dp())/v_r_dp.size());
    else if(params._scale_dp_r()=="log")v_r_dp[i]=pow(10,(log10(params._rmin_dp())+i*(log10(params._rmax_dp())-log10(params._rmin_dp()))/v_r_dp.size()));
    v_density_profile_r[i]=Dp.density_r(v_r_dp[i], log10(Mnl/M_reference), redshift, (void *)&s_cosmo_par);
  }
  So.DONE();
  File.write_to_file(params._density_profile_r_output_file(),v_r_dp,v_density_profile_r);
  So.message_screen("Density profiles in Fourier space");
  v_k_dp.resize(params._npoints_dp_k(),0);
  v_density_profile_k.resize(params._npoints_dp_k(),0);
  for(int i=0;i<v_k_dp.size();i++){
#ifdef TIME
    comp_time(start, v_k_dp.size(), i);
#endif
    if(params._scale_dp_k()=="linear")v_k_dp[i]=(params._k_min_dp()+i*(params._k_max_dp()-params._k_min_dp())/v_r_dp.size());
    else if(params._scale_dp_k()=="log")v_k_dp[i]=pow(10,(log10(params._k_min_dp())+i*(log10(params._k_max_dp())-log10(params._k_min_dp()))/v_k_dp.size()));
    v_density_profile_k[i]=Dp.density_k(v_k_dp[i], log10(Mnl/M_reference), redshift, &s_cosmo_par);
  }
  File.write_to_file(params._density_profile_k_output_file(),v_k_dp,v_density_profile_k);
  So.DONE();
  s_cosmo_par.v_density_profile_k=v_density_profile_k;
  s_cosmo_par.v_density_profile_r=v_density_profile_r;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
 #endif
  // ***********************************************************************************************
#ifdef  _GET_HM_POWER_SPECTRUM_
  // ***********************************************************************************************
  // POWER SPECTRUM IN THE HALO MODEL.
  // ***********************************************************************************************
  s_cosmo_par.hod_model=params._hod_model();
  s_cosmo_par.mmin_hod=params._mmin_hod();
  s_cosmo_par.alpha_hod=params._alpha_hod();
  s_cosmo_par.scatter_hod=params._scatter_hod();
  s_cosmo_par.muno_hod=params._muno_hod();
  // ***********************************************************************************************
  // 1 HALO TERM: SATELLITE-SATELLIT AND CENTRAL-SATELLITE.
  // 2 HALO TERM: SATELLITE1-SATELLIT2, CENTRAL 1-SATELLITE 2, CENTRAL-CENTRAL:
  // condensed in a single scale dependnet bias and the linear power spectrum
    // ***********************************************************************************************
  v_galaxy_power_1h_ss.resize(params._npoints_ps(),0);
  v_galaxy_power_1h_sc.resize(params._npoints_ps(),0);
  v_galaxy_power_2h.resize(params._npoints_ps(),0);
  v_galaxy_matter_bias.resize(params._npoints_ps(),0);
  v_galaxy_power_spectrum.resize(params._npoints_ps(),0);
  real_prec mean_gal_density=Cs.mean_galaxy_number_density(redshift);
  So.message_screen("Mean galaxy number density",mean_gal_density,"(h/Mpc)^3 ");
  real_prec mean_halo_density=Cs.mean_halo_number_density(&s_cosmo_par);
  So.message_screen("Mean halo number density",mean_halo_density,"(h/Mpc)^3");
  Ps.v_mass=v_mass;
  Ps.v_mass_function=v_mass_function;
  Ps.v_halo_mass_bias=v_halo_mass_bias;
  // ***********************************************************************************************
  So.message_screen("Computing galaxy power spectrum");
  for(int i=0;i<v_k_ps.size();i++)
    {
#ifdef TIME
      comp_time(start, v_k_ps.size(), i);
#endif
      v_galaxy_power_1h_ss[i]=Ps.Galaxy_power_spectrum_h1_ss(v_k_ps[i], redshift)/pow(mean_gal_density,2);
      v_galaxy_power_1h_sc[i]=Ps.Galaxy_power_spectrum_h1_sc( v_k_ps[i], redshift)/pow(mean_gal_density,2);
      v_galaxy_matter_bias[i]=Ps.Galaxy_matter_bias( v_k_ps[i], redshift)/mean_gal_density;
      v_galaxy_power_2h[i]=v_l_power_spectrum[i]*(v_galaxy_matter_bias[i],2);
      v_galaxy_power_spectrum[i]=v_galaxy_power_1h_ss[i]+v_galaxy_power_1h_sc[i]+v_galaxy_power_2h[i];
    }
  File.write_to_file(params._galaxy_power_spectrum_halo_model_output_file(),v_k_ps,v_galaxy_power_1h_ss,v_galaxy_power_1h_sc,v_galaxy_matter_bias,v_galaxy_power_2h,v_galaxy_power_spectrum);
  s_cosmo_par.v_galaxy_power_spectrum_1h_ss=this->v_galaxy_power_1h_ss;
  s_cosmo_par.v_galaxy_power_spectrum_1h_sc=this->v_galaxy_power_1h_sc;
  s_cosmo_par.v_galaxy_power_spectrum_2h=this->v_galaxy_power_2h;
  Ps.set_cosmo_pars(this->s_cosmo_par); // update
// ***********************************************************************************************
// GALAXY_CORRELATION FUNCTION HALO MODEL
// Para calcular esto es obligatorio calcular el espectro
#ifdef _GET_HM_CORRELATION_FUNCTIONS_
  // ***********************************************************************************************
  if(params._compute_output_non_linear_correlation_function())
    {
      v_galaxy_correlation_1h_ss.resize(params._npoints_cf(),0);
      v_galaxy_correlation_1h_sc.resize(params._npoints_cf(),0);
      v_galaxy_correlation_2h.resize(params._npoints_cf(),0);
      v_galaxy_correlation.resize(params._npoints_cf(),0);
      So.message_screen("Computing galaxy correlation function");
      for(int i=0;i<v_r_cf.size();i++){
#ifdef TIME
	comp_time(start, v_r_cf.size(), i);
#endif
	v_galaxy_correlation_1h_ss[i]=SCf.Galaxy_Correlation_Function_1h_ss(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation_1h_sc[i]=SCf.Galaxy_Correlation_Function_1h_sc(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation_2h[i]=SCf.Galaxy_Correlation_Function_2h(&s_cosmo_par, v_r_cf[i]);
	v_galaxy_correlation[i]=v_galaxy_correlation_2h[i]+v_galaxy_correlation_1h_ss[i]+v_galaxy_correlation_1h_sc[i];
      }
      So.DONE();
      File.write_to_file(params._galaxy_correlation_function_halo_model_output_file(),
		       v_r_cf,
		       v_galaxy_correlation_1h_ss,
		       v_galaxy_correlation_1h_sc,
		       v_galaxy_correlation_2h,
		       v_galaxy_correlation
		       );
      
    }
#endif // endif for _GET_HM_CORRELATION_FUNCTIONS_
#endif // endif for _GET_HM_POWER_SPECTRUM_
 ending_message();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
