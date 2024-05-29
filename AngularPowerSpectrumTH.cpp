/** @file AngularPowerSpectrumTH.cpp
 *
 * @brief This file contains the class AngularPowerTH 
 * @details Generates theoretical predictions of angular power spectrum
 * @author Andres Balaguera Antolinez
 */

# include "../Headers/AngularPowerSpectrumTH.h"

#define BESS_LIMIT_1 2E-3
#define BESS_LIMIT_2 12.0
#define nbess 100

double top_hat(double x){
  double ans;
  if(fabs(x)>0.5)ans=0.0;
  else if(fabs(x)==0.5)ans=0.5;
  else if(fabs(x)<0.5)ans=1.0;
  return ans;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::read_pars(string &file){
  ifstream fin_parameters (file.c_str());
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<file<<"!"<<endl; exit(1); }
  // a line in parameter file has a form: par_name = par_value

  string line_in_file;
  string par_name;
  string equality;
  string par_value;

  while (getline(fin_parameters,line_in_file)) {
   
    // ignore lines starting with hashtag
    if (line_in_file[0] != '#' && line_in_file.empty()==0) {

      // read parameter name 
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;

      // check that second word is "="
      line_string >> equality;

      if (equality != "=") {
  cerr << "Errorr in parameters file at " << par_name << " = " << "???" << endl;
  cerr << "Using a default value for " << par_name << endl; exit(1);
      }

      // read parameter value 
      line_string >> par_value;
      
      if (par_value.empty()) {
        cout << "Value of " << par_name << " not specified in " <<file << endl;
        cout << "Assuming a default value for " << par_name << endl;
        continue;
      }
      
      if (par_name == "Lmax")Lmax = atoi(par_value.c_str());
      else if (par_name == "Lmin")Lmin = atoi(par_value.c_str());
      else if (par_name == "N_L_bins")N_L_bins = atoi(par_value.c_str());
      else if (par_name == "nside")nside = atoi(par_value.c_str());

      else if (par_name == "code")code = par_value;
     
      else if (par_name == "use_non_linear_pk"){
        if(par_value=="true")this->use_non_linear_pk=true;
        else this->use_non_linear_pk=false;
      }
      else if (par_name == "type_of_nl_power") this->type_of_nl_power = par_value;

      else if (par_name == "use_zmin_zmax_from_dndz_file"){
        if(par_value=="true")this->use_zmin_zmax_from_dndz_file=true;
        else this->use_zmin_zmax_from_dndz_file=false;
      }
      else if (par_name == "pdf_zerrors") this->pdf_zerrors = par_value;
      else if (par_name == "redshift")this->redshift = par_value;
      else if (par_name == "sigma_p_errors")this->sigma_p_errors = (real_prec)(atof(par_value.c_str()));
      else if (par_name == "window_type")this->window_type = par_value;
      else if (par_name == "z_min_clth")this->z_min_clth= (real_prec)(atof(par_value.c_str()));
      else if (par_name == "z_max_clth")this->z_max_clth= (real_prec)(atof(par_value.c_str()));
      else if (par_name == "zmin_bin")this->zmin_bin= (real_prec)(atof(par_value.c_str()));
      else if (par_name == "zmax_bin")this->zmax_bin= (real_prec)(atof(par_value.c_str()));
      else if (par_name == "L_break")this->L_break= atoi(par_value.c_str());
      else if (par_name == "input_cl_file")this->input_cl_file = par_value;
      else if (par_name == "input_mixing_matrix_lbins")this->input_mixing_matrix_lbins = par_value;
      else if (par_name == "input_mixing_matrix_l")this->input_mixing_matrix_l = par_value;
      else if (par_name == "input_dndz_file")this->input_dndz_file = par_value;
      else if (par_name == "input_nbar_file")this->input_nbar_file = par_value;
      else if (par_name == "input_pk_file")this->input_pk_file = par_value;
      else if (par_name == "output_clth_file")this->output_clth_file = par_value;
      else if (par_name == "number_of_kmodes")this->number_of_kmodes = atoi(par_value.c_str());
      else if (par_name == "k_max_integration_cl")this->k_max_integration_cl = (real_prec)atof(par_value.c_str());
      else if (par_name == "k_min_integration_cl")this->k_min_integration_cl = (real_prec)atof(par_value.c_str());
      // Get cosmological parameters from parameter file and put them in a structure
      // of type s_CosmologicalParameters
      else if (par_name == "Om_matter")this->gal.scp.Om_matter = (real_prec)atof(par_value.c_str());  //anyway will be derived
      else if (par_name == "Om_cdm") this->gal.scp.Om_cdm = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_radiation")this->gal.scp.Om_radiation = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_baryons")this->gal.scp.Om_baryons = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_vac")this->gal.scp.Om_vac = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_k")this->gal.scp.Om_k = (real_prec)atof(par_value.c_str());
      else if (par_name == "Hubble")this->gal.scp.Hubble = (real_prec)atof(par_value.c_str());
      else if (par_name == "hubble")this->gal.scp.hubble = (real_prec)atof(par_value.c_str());
      else if (par_name == "n_s")this->gal.scp.n_s = (real_prec)atof(par_value.c_str());
      else if (par_name == "alpha_s")this->gal.scp.alpha_s = (real_prec)atof(par_value.c_str());
      else if (par_name == "w_eos")this->gal.scp.w_eos = (real_prec)atof(par_value.c_str());
      else if (par_name == "N_eff")this->gal.scp.N_eff = (real_prec)atof(par_value.c_str());
      else if (par_name == "sigma8")this->gal.scp.sigma8 = (real_prec)atof(par_value.c_str());
      else if (par_name == "Tcmb")this->gal.scp.Tcmb = (real_prec)atof(par_value.c_str());
      else if (par_name == "GAL_BIAS")this->gal.scp.GAL_BIAS = (real_prec)atof(par_value.c_str());
      else if (par_name == "alpha_BIAS")this->gal.scp.alpha_BIAS = (real_prec)atof(par_value.c_str());
      else if (par_name == "beta_rsd")this->gal.scp.beta_rsd = (real_prec)atof(par_value.c_str());
      else if (par_name == "Amc")this->gal.scp.Amc = (real_prec)atof(par_value.c_str());

      else if (par_name == "kstar")this->gal.scp.kstar = (real_prec)atof(par_value.c_str());

      else if (par_name == "n_points_mass")this->gal.scp.n_points_mass= (int)atof(par_value.c_str());
      else if (par_name == "n_points_mass_integration")this->gal.scp.n_points_mass_integration= (int)atof(par_value.c_str());


      else if (par_name == "mass_function_fit")this->gal.scp.mass_function_fit= par_value.c_str();
      else if (par_name == "halo_mass_bias_fit")this->gal.scp.halo_mass_bias_fit= par_value.c_str();
      else if (par_name == "density_profile")this->gal.scp.density_profile= par_value.c_str();
      else if (par_name == "hod_model")this->gal.scp.hod_model= (int)atof(par_value.c_str());

      else if (par_name == "scatter_hod")this->gal.scp.scatter_hod= (real_prec)atof(par_value.c_str());

      else if (par_name == "alpha_hod")this->gal.scp.alpha_hod= (real_prec)atof(par_value.c_str());
      else if (par_name == "mmin_hod")this->gal.scp.mmin_hod= (real_prec)atof(par_value.c_str());
      else if (par_name == "muno_hod")this->gal.scp.muno_hod= (real_prec)atof(par_value.c_str());
      else if (par_name == "M_min_effective")this->gal.scp.M_min_effective= (real_prec)atof(par_value.c_str());
      else if (par_name == "M_max_effective")this->gal.scp.M_max_effective= (real_prec)atof(par_value.c_str());

      else if (par_name == "Delta_SO")this->gal.scp.Delta_SO= (real_prec)atof(par_value.c_str());


      if(this->gal.scp.n_points_mass< this->gal.scp.n_points_mass_integration){
          cout<<RED<<"Warning. You might have problems with integration and interpolation in mass"<<RESET<<endl;exit(1);
      }


//      this->gal.scp.hod_model=3;
//      this->gal.scp.scatter_hod=0.15;
//      this->gal.scp.mmin_hod=pow(10, 11.84);
//      this->gal.scp.alpha_hod=0.85;
//      this->gal.scp.muno_hod=pow(10,11.98);

      this->gal.scp.RR=8.0;
      this->gal.scp.use_wiggles=1;
      this->gal.scp.M_min_mf=this->gal.scp.M_min_effective;
      this->gal.scp.M_max_mf=this->gal.scp.M_max_effective;

   //PASAR DIRECTAMENTE DE LEERLOS DEL PARAMETER FILE A LA ESTRUCTURA SIN DEFINIR VARIABLES EN CL QUE NO VAMOS A USAR
/*
      this->gal.scp.Om_baryons=this->Om_baryons;
//      this->gal.scp.Om_cdm=this->Om_cdm;
      this->gal.scp.beta_rsd=this->beta_rsd; // useful when fitting RSD
  //    this->gal.scp.Om_matter=this->Om_cdm+this->Om_baryons;
      this->gal.scp.f_baryon=this->Om_baryons/this->Om_matter;
      this->gal.scp.Om_radiation=this->Om_radiation;
      this->gal.scp.Om_k=this->Om_k;
      this->gal.scp.Hubble=this->Hubble;
      this->gal.scp.hubble=this->hubble;
      this->gal.scp.n_s=this->n_s;
      this->gal.scp.alpha_s=this->alpha_s;
      this->gal.scp.w_eos=this->w_eos;
      this->gal.scp.N_eff=this->N_eff;
      this->gal.scp.sigma8=this->sigma8;
      this->gal.scp.Tcmb=this->Tcmb ;
      this->gal.scp.mlim=this->mKlim;
      this->gal.scp.use_wiggles=1;
      this->gal.scp.RR=8.0;
      this->gal.scp.GAL_BIAS=this->GAL_BIAS;
      this->gal.scp.alpha_BIAS=this->alpha_BIAS;
      this->gal.scp.kstar=this->kstar;
      this->gal.scp.Amc=this->Amc;
*/




//      this->gal.scp.mass_function_fit="Tinker";
    //  this->gal.scp.halo_mass_bias_fit="Sheth_Tormen";
   //   this->gal.scp.density_profile="nfw";



      this->gal.scp.kmin_int=this->k_min_integration_cl;
      this->gal.scp.kmax_int=this->k_max_integration_cl;


    }
  }



  //---------------------------------------------------------
  /*
  // Check for the correctness of the input parameters
  // e.g
  if(this->window_type != "tophat" || this->window_type!= "gaussian" || this->window_type != "tophat_c_gaussian"){
      cout<<this->window_type<<endl;

      cout<<RED<<"ALERT. Parameter "<<this->window_type<<" does not match any options assigned to parameter 'window_type'. Please set either 'tophat, gaussian or tophat_c_gaussian. "<<RESET<<endl;
    exit(1);
  }
*/

  //---------------------------------------------------------
  this->gal.Cf.check_cosmo_pars(&this->gal.scp);
  //---------------------------------------------------------

  output_clth_file=output_clth_file+"_zerrors_"+this->pdf_zerrors;
  if(!this->use_external_pk_file && this->use_non_linear_pk){
   output_clth_file_bin = output_clth_file+"_nl_"+this->type_of_nl_power+"_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_lbin";
   output_clth_file = output_clth_file+"_nl_"+this->type_of_nl_power+"_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_raw";
  }
  if(!this->use_external_pk_file && !this->use_non_linear_pk){
     output_clth_file_bin = output_clth_file+"_l_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_lbin";
     output_clth_file = output_clth_file+"_l_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_raw";
  }
  else if(this->use_external_pk_file){
    output_clth_file_bin = output_clth_file+"external_pk_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_lbin";
    output_clth_file = output_clth_file+"external_pk_Lbreak_"+to_string(L_break)+"_zbin_"+to_string(this->Bin_index)+"_raw";
   }

  output_clth_file_bin =output_clth_file_bin +".dat";
  output_clth_file =output_clth_file +".dat";

  cout<<"Done"<<endl;
  std::cout<<"*******************************************************"<<endl;
  cout<<RESET<<endl;
 
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::write_cosmo_pars(){
  std::cout<<CYAN;
  std::cout<<"Cosmological parameters:"<<endl;
  std::cout<<"Omega_matter = "<<this->gal.scp.Om_matter<<endl;
  std::cout<<"Omega_vac = "<<this->gal.scp.Om_vac<<endl;
  std::cout<<"Omega_baryons = "<<this->gal.scp.Om_baryons<<endl;
  std::cout<<"Omega_k = "<<this->gal.scp.Om_k<<endl;
  std::cout<<"w_eos = "<<this->gal.scp.w_eos <<endl;
  std::cout<<"h = "<<this->gal.scp.hubble <<endl;
  std::cout<<"sigma_8 = "<<this->gal.scp.sigma8<<endl;
  std::cout<<"n_s = "<<this->gal.scp.n_s <<endl;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::set_vectors(){
    this->lvec.resize(this->Lmax+1);
    this->Clvec.resize(this->Lmax+1);
    this->lbin.resize(this->N_L_bins,0);
    this->lbin_min.resize(this->N_L_bins,0);
    this->lbin_max.resize(this->N_L_bins,0);
    this->Clbin.resize(this->N_L_bins,0);
    this->Rll_bin.resize(this->N_L_bins);
    for(int i=0;i<this->Rll_bin.size();i++)
      this->Rll_bin[i].resize(this->Lmax+1,0);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::set_Lbins(){

   if(this->bin_type=="linear"){
    real_prec deltal=((real_prec)(this->Lmax-this->Lmin))/((real_prec)lbin.size()); 
    //here and just here I define N_L_bins as a real_prec instead of an integer
    cout<<CYAN<<endl;
    std::cout<<"*********************************************************"<<endl;
    cout<<BLUE<<"L-bins INFO "<<endl;
    cout<<"Bin-averaged power spectrum"<<endl;
    cout<<"Number of bins = "<<this->lbin.size()<<endl;
    cout<<"l-Bin width = "<<deltal<<RESET<<endl;
    std::cout<<"*********************************************************"<<endl;
    cout<<RESET<<endl;
    for(int i=0;i<lbin.size();i++)this->lbin[i]=this->Lmin+(i+0.5)*deltal;
    for(int i=0;i<lbin.size();i++)this->lbin_min[i]=this->Lmin+(i)*deltal;
    for(int i=0;i<lbin.size();i++)this->lbin_max[i]=this->Lmin+(i+1)*deltal;
  }
  else{
    if(this->bin_type=="log"){
      real_prec deltal=log10(Lmax/1.0)/((real_prec)lbin.size());
      for(int i=0;i<lbin.size();i++)lbin[i]=pow(10,log10(1)+(i+0.5)*deltal);
      for(int i=0;i<lbin.size();i++)lbin_min[i]=pow(10,log10(1)+i*deltal);
      for(int i=0;i<lbin.size();i++)lbin_max[i]=pow(10,log10(1)+(i+1)*deltal);
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerTH::get_lg_map(Healpix_Map<real_prec>map, Healpix_Map<real_prec>&lmap){

  // THis returns a delta that is now log-normal distributed
  real_prec mean=0.0;
  for(int i=0;i<map.Npix();i++)mean+=map[i];

  mean/=(real_prec(map.Npix()));
  for(int i=0;i<map.Npix();i++)map[i]-=mean;

  real_prec var=0.0;
  for(int i=0;i<map.Npix();i++)var+=pow(map[i]-mean,2);
  var/=(map.Npix());
  real_prec m=-0.5*log(1+var);
  for(int i=0;i<map.Npix();i++)lmap[i]=exp(map[i]+m)-1.;

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_pars_mask(){
  cout<<RED;
  this->area_pixel=4.*M_PI/(real_prec)n_total_pixels;
  this->area_survey=(real_prec)n_observed_pixels*area_pixel;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_Cl_gaussian(int seed, vector<real_prec>&cl_th, vector<real_prec>&cl_fs){
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  gsl_rng_default_seed=seed;
  T = gsl_rng_ranlux;
  r = gsl_rng_alloc (T);
  for(int l=0;l<=Lmax;l++)cl_fs[l]=cl_th[l]+gsl_ran_gaussian(r, sqrt(1/(2*l+1))*cl_th[l]);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::get_bessel(){


 cout<<BLUE<<"Computing Bessel functions using a mixture of techniques"<<RESET<<endl;

 this->xBessel.resize(nbess+1);
 this->Bessel.resize(this->L_break+1);
 for(int i=0;i<Bessel.size();i++)Bessel[i].resize(nbess+1);
 //#pragma omp parallel for
  for(int i=0;i<this->xBessel.size();i++)xBessel[i]=pow(10,-3+i*(1+3)/((real_prec)nbess));
  //Using the gsl_sf_bessel_jl_steed_arrayI have problems for x<12
  // For small values of x, say, x<1e-3, I use the asymptotic limit
  for(int i=0;i<this->xBessel.size();i++)
  {
    if(xBessel[i]<BESS_LIMIT_1)
     for(int li=0;li<=this->L_break;++li)
      Bessel[li][i]= pow(2, li)*factorial(li)*pow(xBessel[i],li)/factorial(2.0*li+1.0);
    else if (xBessel[i]>=BESS_LIMIT_1 && xBessel[i]<=BESS_LIMIT_2)
      for(int li=0;li<=this->L_break;++li)Bessel[li][i]=gsl_sf_bessel_jl(li,xBessel[i]);
    else if (xBessel[i]>BESS_LIMIT_2)
    {
      real_prec *jl_array= new real_prec[this->L_break+1];
      gsl_sf_bessel_jl_steed_array(this->L_break,xBessel[i],jl_array);
      for(int li=0;li<=this->L_break;++li)Bessel[li][i]=jl_array[li];
      delete[] jl_array;
    }
  }
  cout<<BLUE<<"Done "<<RESET<<endl;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::set_interpolation_tools_dndz(){
// This funciton has to be used after having built the arrays this->dz and this->dndz

  So.enter(__PRETTY_FUNCTION__);
  if(this->dz.size()==0  || this->dn_photo.size() ==0)
   {
     So.message_screen("Requested arrays for interpolation in dndz are not yet resized. Code ends here.");
     So.message_screen("Size of dz=",dz.size());
     So.message_screen("Size of dN=",dn_photo.size());
     exit(0);
   }
/*  else
   {
     this->acc_dndz = gsl_interp_accel_alloc ();
     this->spline_dndz= gsl_spline_alloc (gsl_interp_linear, this->dz.size());
     gsl_spline_init (this->spline_dndz, &(this->dz[0]), &(this->dn_photo[0]), this->dz.size());
   }*/

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



void AngularPowerTH::compute_int_table(){

  So.enter(__PRETTY_FUNCTION__);
  
  int nss_z=100;
  wf =  gsl_integration_glfixed_table_alloc (nss_z);
  WW_z.resize(nss_z);
  XX_z.resize(nss_z);
//  gsl_get_GL_weights(1.02*this->z_min_clth,0.98*z_max_clth,wf,this->XX_z,this->WW_z);
  gsl_get_GL_weights(1.05*this->z_min_clth,0.95*this->z_max_clth,wf,this->XX_z,this->WW_z);

//  if(this->use_non_linear_pk && this->type_of_nl_power=="pt")
 //   this->Ps.compute_int_table_k_mu(this->k_min_integration_cl, this->k_max_integration_cl, 20, 20);

  if(this->L_break>0){
    int nss_z_cl=100;
    wfc =  gsl_integration_glfixed_table_alloc (nss_z_cl);
    WW_z_cl.resize(nss_z_cl);
    XX_z_cl.resize(nss_z_cl);
    gsl_get_GL_weights(1.02*this->z_min_clth,0.98*z_max_clth,wfc,this->XX_z_cl,this->WW_z_cl);

    int nss_k=900;
    if(nss_k<this->number_of_kmodes){
      cout<<RED<<"Please, increase the number of steps for LG integration, at least equal or greater than the parameter number_of_kmodes."<<endl;
      cout<<"Exit"<<RESET<<std::endl;
      exit(1);
    }
    wfd =  gsl_integration_glfixed_table_alloc (nss_k);
    WW_k.resize(nss_k);
    XX_k.resize(nss_k);
    gsl_get_GL_weights(log10(1.2*this->k_min_integration_cl),log10(0.8*this->k_max_integration_cl),wfd,XX_k,WW_k);
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

real_prec AngularPowerTH::window(real_prec z, void *p)
{
  struct params_clth * s_p= (struct params_clth *)p;

  real_prec width=0.5*(s_p->zmax_bin - s_p->zmin_bin);
  real_prec zmean = 0.5*(s_p->zmax_bin + s_p->zmin_bin);
  real_prec x = (z-zmean)/(2.*width);
  real_prec ans=0;
  if(s_p->wtype=="tophat"){
    ans= top_hat(x);
  }
  else if(s_p->wtype=="gaussian"){
    ans= exp(-2.0*pow(x,2));
  }
  else if (s_p->wtype=="tophat_c_gaussian"){
    real_prec sigma_p=s_p->sigma_p_errors;
    ans= gsl_sf_erf((z-s_p->zmin_bin)/sigma_p)-gsl_sf_erf((z-s_p->zmax_bin)/sigma_p);
  }
  return ans;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Spectroscopic redshift, as a convolution of the
// photometric with a pdf

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::read_dndz_photo(){


    if(true==this->read_dndz_from_file)
    {
    vector< vector<real_prec> > dndz;
    this->File.read_file(this->input_dndz_file,dndz);
    for(int i=0;i<dndz.size();++i)this->dz.push_back(dndz[i][0]);
    for(int i=0;i<dndz.size();++i)this->dn_photo.push_back(dndz[i][1]);
    if(this->use_zmin_zmax_from_dndz_file){
      this->z_min_clth= dndz[0][0];
      this->z_max_clth= dndz[this->  dz.size()-1][0];
    }
    else {
      if(this->z_max_clth > dndz[this->  dz.size()-1][0]){
        std::cout<<RED<<"Warning: Maximum z in parameter file ("<<zmax<<") greater than maximum value found in the dNdz ("<<this->dz[this->dz.size()-1] <<"). Setting value from catalog"<<RESET<<endl;
        this->z_max_clth= dndz[this->  dz.size()-1][0];
      }
      else if(this->zmin < dndz[0][0]){
      std::cout<<RED<<"Warning: Minimum z in parameter file ("<<zmin<<") smaller than minimum value found in the dNdz ("<<dndz[0][2]<<"). Setting value from catalog"<<RESET<<endl;
      this->z_min_clth= dndz[0][0];
      }
    }
    dndz.clear();
   }

   cout<<CYAN<<"Minimum z for integration = "<<this->z_min_clth<<RESET<<endl;
   cout<<CYAN<<"Maximum z for integration = "<<this->z_max_clth<<RESET<<endl;
   cout<<CYAN<<"Minimum z for this bin = "<<this->zmin_bin <<RESET<<endl;
   cout<<CYAN<<"Maximum z for this bin = "<<this->zmax_bin<<RESET<<endl;

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::read_nbar_photo(){
    vector< vector<real_prec> > nbar;
    this->File.read_file(this->input_nbar_file,nbar);
    for(int i=0;i<nbar.size();++i)this->dz.push_back(nbar[i][0]);
    for(int i=0;i<nbar.size();++i)this->nbar_photo.push_back(nbar[i][1]);
    if(this->use_zmin_zmax_from_dndz_file){
      this->z_min_clth= nbar[0][0];
      this->z_max_clth= nbar[this->  dz.size()-1][0];
    }
    else {
      if(this->z_max_clth > nbar[this->  dz.size()-1][0]){
        std::cout<<RED<<"Warning: Maximum z in parameter file ("<<zmax<<") greater than maximum value found in the dNdz ("<<this->dz[this->dz.size()-1] <<"). Setting value from catalog"<<RESET<<endl;
        this->z_max_clth= nbar[this->  dz.size()-1][0];
      }
      else if(this->zmin < nbar[0][0]){
      std::cout<<RED<<"Warning: Minimum z in parameter file ("<<zmin<<") smaller than minimum value found in the dNdz ("<<nbar[0][2]<<"). Setting value from catalog"<<RESET<<endl;
      this->z_min_clth= nbar[0][0];
      }
    }
    nbar.clear();
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

struct sL{
  real_prec zs;
  real_prec s;
};

real_prec iLorentz(real_prec zp, void *p){
   struct sL * sL = (struct sL *)p;
   return pow(pow((sL->zs-zp)/(2.*sL->s),2)+1, -3);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::PDF_errors(string te, real_prec zp, real_prec zs, real_prec s){

  s=0.03*tanh(-20.78*zp*zp+7.76*zp+0.05)/(1.+zp);
  real_prec ans=0;
  if(te=="gaussian"){
    ans=exp(-0.5*pow(zs-zp,2)/pow(s,2))/(sqrt(2*M_PI)*s);
  }
  else if(te=="lorentz"){
//    return  (s/(pow(zs-zp,2)+pow(s,2)))/M_PI;
    struct sL sL;
    sL.s=s; sL.zs=zs;
    real_prec Normal_l=gsl_integration3(100, iLorentz, (void *)&sL, 0.0, 0.5);
    ans=  pow(pow((zs-zp)/(2.*s),2)+1, -3)/Normal_l;
 }
 return ans;

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// This integrand is used in the function get_dndz_spect
gsl_real AngularPowerTH::idndz_spect(real_prec zp, void *p){
  AngularPowerTH cl;
  struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum> *)p;
  params_clth * s_p = s_aa->s_clth;
  return gsl_inter_new(s_p->dz,s_p->dn_photo,zp)*cl.window(zp,(void *)s_p)*cl.PDF_errors(s_p->pdf_zerrors, zp,s_aa->zaux,s_p->sigma_p_errors);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// This is a simple interpolation of the results obtained by the function get_dndz_spect
real_prec AngularPowerTH::dndz_spect(real_prec z, void *p){
  struct params_clth * s_p= (struct params_clth *)p;
  return gsl_inter_new(s_p->dz,s_p->dn_spect,z);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
gsl_real AngularPowerTH::idndz_spect_norm(real_prec zs, void *p)
{
  AngularPowerTH cl;
  return cl.dndz_spect(zs,p);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// DNDZ PHTOMETRic, this is an interpolation of the input file with the
// photo.z multiplied by top hat window function defining the photometric z-bin
real_prec AngularPowerTH::dndz_photo(real_prec z, void *p){
  struct params_clth * s_p= (struct params_clth *)p;
//  real_prec dndz= gsl_spline_eval (this->spline_dndz,z, this->acc_dndz);
  real_prec dndz=gsl_inter_new(s_p->dz,s_p->dn_photo,z);
  real_prec _window=window(z,p);
  return dndz*_window;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Integrand of int dndz
gsl_real AngularPowerTH::idndz(real_prec zp, void *p){
  AngularPowerTH cl;
  return cl.dndz_photo(zp,p);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

gsl_real AngularPowerTH::idndz_spect_mean_redshift(real_prec z, void *p){
  AngularPowerTH cl;
  return z*cl.dndz_spect(z,p);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::get_dndz_spect(params_clth sp, int in){
// The variable in works for outputs

    std::cout<<CYAN<<"Obtaining dN/dz spectroscopic...";
    time_t start_all;
    time(&start_all);
    time_t end;
   dn_spect.resize(dn_photo.size(),0);
    real_prec res=this->dz[this->dz.size()-1]/((real_prec)this->dz.size());
   if(res >this->sigma_p_errors){
    std::cout<<RED<<"ALERT:"<<endl;
    std::cout<<"Value of parameter 'sigma_p_errors' ("<<this->sigma_p_errors<<") smaller than the resolution found in the dNdz file ("<<res<<")"<<std::endl;
    std::cout<<"This will generate an unstable convolution with the PDF P(zs|zp)."<<std::endl;
    std::cout<<"Please increase value of  'sigma_p_errors' or use a dNdz input file with more lines"<<std::endl;
    std::cout<<"Stoping here. Aufwiedersehen"<<RESET<<std::endl;
    exit(1);
   }
#pragma omp parallel for
    for(int i=0;i<this->dn_photo.size();++i){
      s_aux<PowerSpectrum>  ssb;
      ssb.s_clth=&sp;
      ssb.zaux=this->dz[i];
      this->dn_spect[i]=gsl_integration2(idndz_spect,(void *)&ssb, this->XX_z, this->WW_z);
    }

    //string dbin="2MPZ_dndz_zbin"+to_string(this->Bin_index)+"_"+this->pdf_zerrors+"_sigmaz_"+to_string(in)+".dat";
    string dbin="2MPZ_dndz_zbin"+to_string(this->Bin_index)+"_"+this->pdf_zerrors+"_sigmaz_var"+".dat";

    this->File.write_to_file(dbin, dz,dn_photo,dn_spect);

    time(&end);
    real_prec lapse=difftime(end,start_all);
    cout<<"Done"<<RESET<<endl;
    if(this->code!="cl_mcmc")std::cout<<BLUE<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_normal(void *p){
   So.enter(__PRETTY_FUNCTION__);
    real_prec aa=gsl_integration2(idndz, p, this->XX_z,this->WW_z );
    this->normal_cl_power_photo=aa;
    if(this->pdf_zerrors=="gaussian" || this->pdf_zerrors=="lorentz"){
        this->normal_cl_power=gsl_integration2(idndz_spect_norm, p, this->XX_z,this->WW_z);
    }
    else if (this->pdf_zerrors=="delta"){
        this->normal_cl_power=aa;
    }
}
// *******************************************************************************************************************************************************
void AngularPowerTH::get_normal(){
    So.enter(__PRETTY_FUNCTION__);
    real_prec aa=gsl_integration2(idndz, (void *)&this->pclth, this->XX_z,this->WW_z );
    this->normal_cl_power_photo=aa;
    this->normal_cl_power=aa;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_mean_redshift(void *p){
  this->mean_redshift= gsl_integration2(idndz_spect_mean_redshift, p, this->XX_z,this->WW_z)/this->normal_cl_power;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

real_prec AngularPowerTH::iKernel(real_prec z, void *p)
{
  AngularPowerTH cl;
  PowerSpectrum Ps;
  struct s_aux<PowerSpectrum> * s_aa = (struct s_aux <PowerSpectrum>*)p;
  params_clth * s_p = s_aa->s_clth;
  s_CosmologicalParameters * s_cp = s_aa->scp_a;

  real_prec Bessel, Bessel_l=0,Bessel_la=0, Bessel_lb=0;
  real_prec r= gsl_inter_new(s_p->zv,s_p->rv,z);
  real_prec bb= gsl_inter_new(s_p->zv,s_p->bias_zv,z);
 // real_prec gg= gsl_inter_new(s_p->zv,s_p->gv,z);
 // real_prec knl= gsl_inter_new(s_p->zv, s_p->klnv,z);
 // real_prec h4= gsl_inter_new(s_p->zv,s_p->HF_i4,z);
 // real_prec h2= gsl_inter_new(s_p->zv,s_p->HF_i2,z);

  real_prec x= (s_aa->kaux)*r;
  //real_prec beta= gsl_inter_new(s_p->zv,s_p->gfv,z)/bb;
  real_prec beta=s_cp->beta_rsd;
  int l = s_aa->laux;

  real_prec power=1.0; // P ouTSIDE

 // real_prec power=(s_p->use_non_linear_pk? Ps.Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(s_cp,s_aa->kaux,z,gg, h4,h2,knl):Ps.Linear_Matter_Power_Spectrum_z(s_cp,z,gg)) ;
  //power*= (bb*bb*exp(-pow(s_aa->kaux/s_cp->kstar,2)));  // P inside


  real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);

  Bessel= (l<=1? 0 : gsl_sf_bessel_jl(l,x));

  // We start RSD from l=2, not t take the dipole into account

  real_prec Al= (2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1.));
  real_prec Bl= -l*(l-1.)/((2*l-1.)*(2.*l+1.));
  real_prec Cl= -(l+2.)*(l+1.)/((2.*l+1.)*(2.*l+3.));

  Bessel_l= (l<=1? 0 : Al*Bessel);

  Bessel_la=(l<=1? 0 : Bl*gsl_sf_bessel_jl(l-2,x));

  Bessel_lb=(l<=1? 0 : Cl*gsl_sf_bessel_jl(l+2,x));

  return dndz*sqrt(power)*(Bessel + beta*(Bessel_l+Bessel_la+Bessel_lb));
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

real_prec AngularPowerTH::get_Kernel(void *p){
 return gsl_integration2(iKernel, p, this->XX_z_cl,this->WW_z_cl); //esta da problema
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

real_prec AngularPowerTH::iCl(real_prec kl, void *p){
  real_prec k=pow(10,kl);
  struct params_clth * s_p= (struct params_clth *)p;
  real_prec Fkernel = gsl_inter_new(s_p->kv,s_p->Fkernel,k);
  real_prec power=gsl_inter_new(s_p->kv, s_p->pk, k);
//  real_prec power=1; // Pinside
  return (log(10.0)*k)*pow(k*Fkernel,2.0)*power;
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::get_cl_exact(void *p){
 return (2./M_PI)*gsl_integration2(iCl, p, this->XX_k, this->WW_k);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::iget_ngal_halo_model(real_prec z, void *p)
{
   AngularPowerTH cl;
   struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum> *)p;
   params_clth * s_p = s_aa->s_clth;
   s_CosmologicalParameters * s_cp = s_aa->scp_a;
   PowerSpectrum Pst_hm =  s_aa->Ps;
   real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
   real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
   //real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
   real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);

   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass[i]=(log10(s_cp->M_min_effective)+i*(log10(s_cp->M_max_effective)-log10(s_cp->M_min_effective))/Pst_hm.v_mass.size());
     Pst_hm.v_sigma_mass[i]=gg*s_p->sigma_mass[i];
   }
   // This we do in other loop for use interpolation over the vectors created in the previous one
   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass_function[i]= Pst_hm.mass_function_D(Pst_hm.v_mass[i],z,s_cp);
   }
   real_prec mean_gal_density=Pst_hm.mean_galaxy_number_density(s_cp);
   return mean_gal_density*(4.*M_PI*0.693)*Constants::speed_light*(r*r)/HH;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// ******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// test para hcer el hm  mas rapido,  ESTA YA NO SERA ESTATICA, PUEDE VER TODO!!
real_prec AngularPowerTH::iCl_limber_non_linear_hm_FAST(int iz, void *p)
{
   AngularPowerTH cl;
   struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum> *)p;
   params_clth * s_p = s_aa->s_clth;
   s_CosmologicalParameters * s_cp = s_aa->scp_a;
   PowerSpectrum Pst_hm =  s_aa->Ps;


   real_prec z = this->XX_z[iz];

   real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
   real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
   real_prec ks=(s_aa->laux+0.5)/r;
   real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
   real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
   real_prec power;

   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass_function[i]= Pst_hm.MASS_FUNCTION_M_Z[iz][i];
     Pst_hm.v_halo_mass_bias[i]= Pst_hm.MASS_BIAS_M_Z[iz][i];
     //cout<< Pst_hm.v_mass[i]<<"  "<<Pst_hm.v_halo_mass_bias[i]<<endl;
   }

   s_cp->Mnl=gsl_inter_new(s_p->zv,s_p->M_nl,z);     // Non linear masses, needed for the density profile omputed at z=0 in get-HM_scales
    // el z en mean_gal_density aca no juega ningun papel, por eso mando un z=0
   real_prec mean_gal_density=Pst_hm.mean_galaxy_number_density(s_cp);
   real_prec b0=Pst_hm.Galaxy_matter_bias(s_cp,ks,z);
   real_prec bb=b0/mean_gal_density;
   real_prec power_ss = Pst_hm.Galaxy_power_spectrum_h1_ss(s_cp,ks,z)/pow(mean_gal_density,2);
   real_prec power_sc = Pst_hm.Galaxy_power_spectrum_h1_sc(s_cp,ks,z)/pow(mean_gal_density,2);
   real_prec power1H=  power_ss+power_sc;
   real_prec power2H= bb*bb*Pst_hm.Linear_Matter_Power_Spectrum_z(s_cp,ks,gg);
   power= power2H + power1H;
   return dndz*dndz*power*HH/(r*r);
}
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::iCl_limber_non_linear_hm(real_prec z, void *p)
{
   AngularPowerTH cl;
   struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum> *)p;
   params_clth * s_p = s_aa->s_clth;
   s_CosmologicalParameters * s_cp = s_aa->scp_a;
   PowerSpectrum Pst_hm =  s_aa->Ps;
   real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
   real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
   real_prec ks=(s_aa->laux+0.5)/r;
   real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
   real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
   real_prec power;



   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass[i]=(log10(s_cp->M_min_effective)+i*(log10(s_cp->M_max_effective)-log10(s_cp->M_min_effective))/Pst_hm.v_mass.size());
     Pst_hm.v_sigma_mass[i]=gg*s_p->sigma_mass[i];
   }
   // This we do in other loop for use interpolation over the vectors created in the previous one
   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass_function[i]= Pst_hm.mass_function_D(Pst_hm.v_mass[i],z,s_cp);
     Pst_hm.v_halo_mass_bias[i]=Pst_hm.bias(Pst_hm.v_mass[i],z,s_cp);

   }


   // **************************************************************
/*   int IZZ;
   for(int IZ=0;IZ<s_aa->XX_z.size();IZ++)if(s_aa->XX_z[IZ]=z)IZZ=IZ;
   cout<<IZZ<<"  "<<s_aa->XX_z[IZZ]<<endl;
   for(int i=0; i<Pst_hm.v_mass.size();i++){
     Pst_hm.v_mass_function[i]= 1;//Pst_hm.MASS_FUNCTION_M_Z[IZZ][i];
     Pst_hm.v_halo_mass_bias[i]=1;//Pst_hm.bias(Pst_hm.v_mass[i],z,s_cp);
   }
   // **************************************************************

*/


   s_cp->Mnl=gsl_inter_new(s_p->zv,s_p->M_nl,z);     // Non linear masses, needed for the density profile
    // el z en mean_gal_density aca no juega ningun papel, por eso mando un z=0
   real_prec mean_gal_density=Pst_hm.mean_galaxy_number_density(s_cp);
   real_prec b0=Pst_hm.Galaxy_matter_bias(s_cp,ks,z);
   real_prec bb=b0/mean_gal_density;
   real_prec power_ss = Pst_hm.Galaxy_power_spectrum_h1_ss(s_cp,ks,z)/pow(mean_gal_density,2);
   real_prec power_sc = Pst_hm.Galaxy_power_spectrum_h1_sc(s_cp,ks,z)/pow(mean_gal_density,2);
   real_prec power1H=  power_ss+power_sc;
   real_prec power2H= bb*bb*Pst_hm.Linear_Matter_Power_Spectrum_z(s_cp,ks,gg);
   power= power2H + power1H;
   return dndz*dndz*power*HH/(r*r);
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::iCl_limber_non_linear_hf(real_prec z, void *p){
  AngularPowerTH cl;
  struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum>  *)p;
  params_clth * s_p = s_aa->s_clth;
  s_CosmologicalParameters * s_cp = s_aa->scp_a;
  PowerSpectrum Pst =  s_aa->Ps;

  real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
  real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
  real_prec bb = gsl_inter_new(s_p->zv,s_p->bias_zv,z);
  real_prec ks=(s_aa->laux+0.5)/r;
  real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
  real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
  real_prec power;
  real_prec knl = gsl_inter_new(s_p->zv, s_p->klnv,z);
  real_prec h4= gsl_inter_new(s_p->zv,s_p->HF_i4,z);
  real_prec h2= gsl_inter_new(s_p->zv,s_p->HF_i2,z);
  power=Pst.Non_Linear_Matter_Power_Spectrum_Halo_Fit_z(s_cp,ks,z,gg, h4,h2,knl);
 return dndz*dndz*bb*bb*power*HH/(r*r);

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::iCl_limber_non_linear_pt(real_prec z, void *p){
  AngularPowerTH cl;
  struct s_aux<PowerSpectrum>  * s_aa = (struct s_aux<PowerSpectrum>  *)p;
  params_clth * s_p = s_aa->s_clth;
  s_CosmologicalParameters * s_cp = s_aa->scp_a;
  PowerSpectrum Pst =  s_aa->Ps;

  real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
  real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
  real_prec bb = gsl_inter_new(s_p->zv,s_p->bias_zv,z);
  real_prec ks=(s_aa->laux+0.5)/r;
  real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
  real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
  real_prec power;
  s_cp->growth_factor=gg;
  power=Pst.Non_Linear_Matter_Power_Spectrum_PT(s_cp,ks);
  return dndz*dndz*bb*bb*power*HH/(r*r);

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::iCl_limber_non_linear_pt_inter(real_prec z, void *p){
    /*This creates the Cl using the RPT, in case we want to fit the parameters b, kstar, Amc
     * Here we read the once created P1loop and interpolate, for it only depens
     * on cosmological paramerters and for mcmc will be useless to compute at every step,
     * if we fic the cosmo
     */

  AngularPowerTH cl;
  struct s_aux<PowerSpectrum>  * s_aa = (struct s_aux<PowerSpectrum>  *)p;
  params_clth * s_p = s_aa->s_clth;
  s_CosmologicalParameters * s_cp = s_aa->scp_a;
  PowerSpectrum Pst =  s_aa->Ps;

  real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
  real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
  real_prec bb = gsl_inter_new(s_p->zv,s_p->bias_zv,z);
  real_prec ks=(s_aa->laux+0.5)/r;
  real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
  real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
  s_cp->growth_factor=1.0;
  real_prec P1loop=gsl_inter_new(s_p->kv,s_p->pk,ks);
  real_prec power= gg*gg*exp(-0.5*pow(ks/s_cp->kstar,2))*Pst.Linear_Matter_Power_Spectrum_z(s_cp,ks,1.0)+pow(gg,4)*s_cp->Amc*P1loop;

  return dndz*dndz*bb*bb*power*HH/(r*r);

}
// *******************************************************************************************************************************************************

real_prec AngularPowerTH::iCl_limber_linear(real_prec z, void *p){
  AngularPowerTH cl;
  PowerSpectrum Ps;
  struct s_aux<PowerSpectrum>  * s_aa = (struct s_aux<PowerSpectrum>  *)p;
  params_clth * s_p = s_aa->s_clth;
  s_CosmologicalParameters * s_cp = s_aa->scp_a;
  real_prec gg = gsl_inter_new(s_p->zv,s_p->gv,z);
  real_prec r=gsl_inter_new(s_p->zv,s_p->rv,z);
  real_prec ks=(s_aa->laux+0.5)/r;
  real_prec HH=gsl_inter_new(s_p->zv,s_p->Hv,z);
  real_prec dndz = (s_p->pdf_zerrors=="gaussian" || s_p->pdf_zerrors=="lorentz") ? cl.dndz_spect(z,(void *)s_p) : cl.dndz_photo(z,(void *)s_p);
  real_prec bb = gsl_inter_new(s_p->zv,s_p->bias_zv,z);
  real_prec power= Ps.Linear_Matter_Power_Spectrum_z(s_cp,ks,gg);
  real_prec iCl=dndz*dndz*power*HH/(r*r);
  return iCl;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
real_prec AngularPowerTH::get_cl_limber(void *p){
/* Here I have acces to the definition of the type of NL.
 * then I will select  pointing to differnet functions for HF, PT, HM. This will allow me to pass an
 * object type of PowerSpectrum with the tables of integration with respect to k and masses already computed,
 * saving time.*/
    struct s_aux<PowerSpectrum> * s_aa = (struct s_aux<PowerSpectrum> *)p;
    params_clth * s_p = s_aa->s_clth;

    real_prec ans=0;
    if(false==this->use_non_linear_pk)
      return gsl_integration2(iCl_limber_linear, p, this->XX_z,this->WW_z);
    else
      if(this->type_of_nl_power=="hf"){
        ans=gsl_integration2(iCl_limber_non_linear_hf,p, this->XX_z,this->WW_z);
      }
      else  if(this->type_of_nl_power=="pt"){
        ans=gsl_integration2(iCl_limber_non_linear_pt,p, this->XX_z,this->WW_z);
      }
      else  if(this->type_of_nl_power=="hm")
      {
        real_prec result=0;
        for(int i=0;i<this->WW_z.size();++i)
          result+=this->WW_z[i]*iCl_limber_non_linear_hm_FAST(i,p);
         ans=result;
      }
   
   return ans;


}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


real_prec AngularPowerTH::get_cl_limber_pt_inter(void *p){
/* Here I have acces to the definition of the type of NL.
 * then I will select  pointing to differnet functions for HF, PT, HM. This will allow me to pass an
 * object type of PowerSpectrum with the tables of integration with respect to k and masses already computed,
 * saving time.*/
   return gsl_integration2(iCl_limber_non_linear_pt_inter,p, this->XX_z,this->WW_z);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



void AngularPowerTH::get_pk(real_prec redshift){
    // Getting Pk evaluated at a redshift z. The linear power spectrum is normalized at z=0 using the value of sigma8
    // taken from CMB (input parameter file). Then, the non-linear power spectrum is computed at the redshift z (input)
   // which is meant to be the mean redshift of the spectroscopic bin. With this, the power spectrum takes the form
    // Cl = int dk k**2 P(k,<z>_i)*F_l(k)**2, which should be good if the redshift bins are narrow.
    // If not, then the Pk should go inside the integral definig the Fl


  if(this->use_external_pk_file){
    vector< vector<real_prec> > pk_input;
    this->File.read_file(this->input_pk_file, pk_input);
    for(int i=0;i<pk_input.size();++i)this->k_th.push_back(pk_input[i][0]);
    for(int i=0;i<pk_input.size();++i)this->pk_th.push_back(pk_input[i][1]);
    this->k_min_integration_cl = this->k_th[0];
    this->k_max_integration_cl = k_th[k_th.size()-1];
    this->number_of_kmodes = k_th.size();
    for(int i=0;i<pk_input.size();++i)this->pk_th[i]*=pow(this->gal.scp.GAL_BIAS,2);
    pk_input.clear();
  }
  else{

    this->pk_th.resize(this->number_of_kmodes,0);
    this->k_th.resize(this->number_of_kmodes,0);

    // Growth factor at redshift zero. Recall we have P(k)= T(k)**2 k**n (D(z)/D(z=0))
    // where D(z) is normalized such that D(a)=a for high redshift ,
    // and we simply pass a factor 1.0 to the cosmological functions.

    this->gal.scp.growth_factor=1.0;
    // Normalization of linear matter power spectrum using the value of omega8 from CMB,
    // at redshift zero.
    real_prec pk_normalization;
    this->Ps.normalization((void *)&this->gal.scp,pk_normalization);
    this->gal.scp.pk_normalization=pk_normalization;
    if(this->code!="cl_mcmc")std::cout<<GREEN<<"Normalization of linear P(k) at z = 0  ( sigma 8 = "<<this->gal.scp.sigma8<<" ) : "<<pk_normalization<<RESET<<std::endl;
    vector<real_prec>pk_aux(this->k_th.size(),0);
    vector<real_prec>kk_aux(this->k_th.size(),0);

    omp_set_num_threads(omp_get_max_threads());

   // Now recompute the growth factor at the mean rdshift of the bin,
    this->gal.scp.cosmological_redshift=redshift;
    real_prec gz=gsl_inter_new(this->gal.zv, this->gal.gv,redshift);
    this->gal.scp.growth_factor=gz; // to be passed now to the power at the right z


#pragma omp parallel for
    for(int i=0;i<k_th.size();++i)this->k_th[i]=pow(10, log10(this->k_min_integration_cl)+i*log10(this->k_max_integration_cl/this->k_min_integration_cl)/((real_prec)this->number_of_kmodes));


    if(this->use_non_linear_pk){


      if(this->type_of_nl_power=="hf"){  // Get P(k) from the Halo-Fit

        vector<real_prec>rr(1000,0); // The size of these vectors are related to the precision to which we get Rnl
        vector<real_prec>sum(1000,0);

        real_prec kln, rnl;

        this->gal.scp.kmin_int=this->k_min_integration_cl;
        this->gal.scp.kmax_int=this->k_max_integration_cl;

        // Get the non linear scales. This returns kln and rnl at z=0,
        /*
        // **************************************************************************
        //THIS IS A TEST, PLEASE DELETE AFTER DONE:
        cout<<RED<<"THIS IS A TEST, PLEASE DELETE AFTER DONE:"<<RESET<<endl;
        this->gal.scp.cosmological_redshift=0.;
        this->gal.scp.growth_factor=1.0; // to be passed now to the power at the right z
        real_prec gz=1.0;
        // **************************************************************************
        */

        this->Ps.nl_scales_halo_fit((void *)&this->gal.scp,&kln,&rnl,rr,sum,true);
        real_prec igz=1./pow(gz,2);
        // so we recale-thm at the input redshift
        kln=1./gsl_inter_new(sum,rr,igz);
        rnl=1./kln;
        this->gal.scp.knl_hf =kln;
        this->gal.scp.rnl_hf =rnl;
        #pragma omp parallel for
        for(int i=0;i<k_th.size();++i)pk_th[i]=pow(this->gal.scp.GAL_BIAS,2)*this->Ps.Non_Linear_Matter_Power_Spectrum_Halo_Fit(&this->gal.scp,this->k_th[i]);
      }
      else if(this->type_of_nl_power=="pt"){  // Get P(k) from the PT
          PowerSpectrum Pst(this->k_min_integration_cl, this->k_max_integration_cl ,20,10);
           this->gal.scp.growth_factor=1.0;
           #pragma omp parallel for
       for(int i=0;i<k_th.size();++i)pk_th[i]=Pst.P1loop(&this->gal.scp,this->k_th[i]);
       }
    }
    else if(!this->use_non_linear_pk){  // We use the Liner Power and multiply here for the damping factor
    #pragma omp parallel for
      for(int i=0;i<this->k_th.size();++i)this->pk_th[i]=exp(-0.5*pow(this->k_th[i]/this->gal.scp.kstar, 2))*pow(this->gal.scp.GAL_BIAS,2)*this->Ps.Linear_Matter_Power_Spectrum(&this->gal.scp, this->k_th[i]);
    }
     cout<<"*******************************************************"<<endl;
   }
   string pf = (this->use_non_linear_pk? "pk_nl.dat":"pk_l.dat");
   ofstream po; po.open(pf.c_str());
   for(int i=0;i<this->k_th.size();++i)po<<this->k_th[i]<<" "<<this->pk_th[i]<<endl;
   po.close();

}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_HALO_FIT_scales(){

    // Here we get some quantities used in the HF
    // specially integrals thatn can be made at z=0 and then scaled
    // by the growth factor

   this->gal.scp.kmin_int=this->k_min_integration_cl;
   this->gal.scp.kmax_int=this->k_max_integration_cl;

    // --------------------------------------------------
    // Normalization of the linear matter power spectrum at z=0
    // using the value of omega8 taken from CMB
    // --------------------------------------------------
   this->gal.scp.growth_factor=1.0;
   real_prec pk_normalization;
   this->Ps.normalization((void *)&this->gal.scp,pk_normalization);
   this->gal.scp.pk_normalization=pk_normalization;

   if(this->code!="cl_mcmc")std::cout<<GREEN<<"Normalization of linear P(k) at z = 0  ( sigma 8 = "<<this->gal.scp.sigma8<<" ) : "<<pk_normalization<<RESET<<std::endl;

   this->k_th.resize(this->number_of_kmodes,0);
   for(int i=0;i<k_th.size();++i)this->k_th[i]=pow(10, log10(this->k_min_integration_cl)+i*log10(this->k_max_integration_cl/this->k_min_integration_cl)/((real_prec)this->number_of_kmodes));

   // --------------------------------------------------
   // Escalas no lineales para el HF,
   // Encuentro las escalas a z=0
   vector<real_prec>rr(1000,0);
   vector<real_prec>sum(1000,0);
   this->klnv.resize(this->gal.zv.size(), 0);
   this->HF_i4.resize(this->gal.zv.size(), 0);
   this->HF_i2.resize(this->gal.zv.size(), 0);

   if(this->use_non_linear_pk){

     // Important, since we will evaluate integrals at z=0 and then re-scale
     this->gal.scp.growth_factor=1.0;  // Set g(z) = 1, i.e, z=0
     //  this->gal.scp.cosmological_redshift=0.0;

     real_prec kln,rnl;

     // Get array with rr(scale) and the sigma**2 at z=0
     this->Ps.nl_scales_halo_fit((void *)&this->gal.scp,&kln,&rnl, rr, sum, true);

      // A luego reescalo por el gg, ya que esto se hace con el P lineal
//      ofstream sa; sa.open("kln.dat");

     if(this->code=="cl_theo")cout<<CYAN<<"Evaluating non linear scales at different redshifts ...";
    //      Loop over redshift to get non linear scales at diff z

     ofstream sa;sa.open("cosmo_func.dat");
     for(int i=0;i<this->gal.zv.size();++i){

        //Get 1/growth factor**2 at this redshift
       real_prec gi=1.0/pow(this->gal.gv[i],2);
//        cout<<this->gal.gv[i]<<endl;


        // Find scale k at which sigma**2 = 1/g**2(z)
//       cout<<RED<<pow(this->gal.gv[i],2)<<"  "<<gi<<RESET<<endl;
       real_prec kkln=1.0/gsl_inter_new(sum,rr,gi);

        // Fill vector with these non linear scales kln(z)
       this->klnv[i]=kkln;

        // Get the scale r non linear and put it in structure to compute redshift dependent integrals in HF
       this->gal.scp.aux_var3=1./kkln;

        // Put the growth-factor(z) in the structure to compute reshift dependent integrals in HF
       this->gal.scp.growth_factor=this->gal.gv[i];

        // Compute the integrals for HF whcih already contain the redshift dependence in the
        // linear power spectrum, for we have passed the g(z) to the structure
       real_prec ii4=0, ii2=0;
       this->Ps.halo_fit_integrals((void *)&this->gal.scp, &ii4, &ii2);
        // Fill vectors with these intergrals. as a function of Rnl, -> z
       this->HF_i4[i]=ii4;//integral for diff kln ie diff z
       this->HF_i2[i]=ii2;
 //        sa<<gal.zv[i]<<"  "<<this->gal.gv[i]<<"  "<<this->gal.rv[i]<<"  "<<this->gal.Cf.omega_matter(gal.zv[i], (void *)&this->gal.scp)<<"  "<<ii4<<"  "<<ii2<<"  "<<kkln<<endl;
     }
     sa.close();
     std::cout<<CYAN<<"Done"<<RESET<<endl;
   }
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_HM_scales(){

    // Here we get some quantities used in the HF
    // specially integrals thatn can be made at z=0 and then scaled
    // by the growth factor

    Cosmology Cf;
    PowerSpectrum Psk(this->gal.scp.kmin_int, this->gal.scp.kmax_int,40,10);//40 ponts in k, 2 in mu, low in mu, we do not need it here

    // --------------------------------------------------
    // Normalization of the linear matter power spectrum at z=0
    // using the value of omega8 taken from CMB
    // --------------------------------------------------
   this->gal.scp.growth_factor=1.0;
   real_prec pk_normalization;
   this->Ps.normalization((void *)&this->gal.scp,pk_normalization);
   this->gal.scp.pk_normalization=pk_normalization;

   if(this->code!="cl_mcmc")std::cout<<GREEN<<"Normalization of linear P(k) at z = 0  ( sigma 8 = "<<this->gal.scp.sigma8<<" ) : "<<pk_normalization<<RESET<<std::endl;

   this->k_th.resize(this->number_of_kmodes,0);
   for(int i=0;i<k_th.size();++i)this->k_th[i]=pow(10, log10(this->k_min_integration_cl)+i*log10(this->k_max_integration_cl/this->k_min_integration_cl)/((real_prec)this->number_of_kmodes));

   // --------------------------------------------------
   // Escalas no lineales para el HF,
   // Encuentro las escalas a z=0
   vector<real_prec>rr(1000,0);
   vector<real_prec>sum(1000,0);

   this->Mass_nl.resize(this->gal.zv.size(), 0);

   if(this->use_non_linear_pk){

     // Important, since we will evaluate integrals at z=0 and then re-scale
     this->gal.scp.growth_factor=1.0;  // Set g(z) = 1, i.e, z=0

     real_prec kln,rnl;
     // Get array with rr(scale) and the sigma**2 at z=0
     this->Ps.nl_scales_halo_fit((void *)&this->gal.scp,&kln,&rnl, rr, sum, true);

      // A luego reescalo por el gg, ya que esto se hace con el P lineal

     omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for
     for(int i=0;i<this->gal.zv.size();++i){
        //Get 1/growth factor**2 at this redshift
       real_prec gi=1.0/pow(this->gal.gv[i],2);
       real_prec rrln=gsl_inter_new(sum,rr,gi);
       this->Mass_nl[i]=(4./3.)*M_PI*pow(rrln,3)*Cf.mean_matter_density(this->gal.zv[i],(void *)&this->gal.scp);
     }

    this->sigma_mass.resize(this->gal.scp.n_points_mass,0);
     //Loop over masses
//#pragma omp parallel for
     for(int i=0; i<this->gal.scp.n_points_mass;i++){
       real_prec mass=log10(this->gal.scp.M_min_effective)+i* log10(this->gal.scp.M_max_effective/this->gal.scp.M_min_effective)/(real_prec)this->gal.scp.n_points_mass;
       this->sigma_mass[i]=Psk.sigma_masa(mass,0.0,&this->gal.scp);
     }




   }
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::read_mixing_matrix_lbins(){
   vector< vector<real_prec> > pR;
   File.read_file(this->input_mixing_matrix_lbins,pR);
   int ik=-1;
   for(int i=0;i<this->N_L_bins;i++){
     for(int j=this->Lmin;j<this->Lmax+1;++j){
       ++ik;
       this->Rll_bin[i][j]=pR[ik][2];
     }
   }
  pR.clear();
  // Normalize mixing matrix
  for(int l=0;l<this->lbin.size();l++){
      real_prec ch=0;
      for(int lp=this->Lmin;lp<this->Lmax+1;lp++)ch+=this->Rll_bin[l][lp];
      for(int lp=this->Lmin;lp<this->Lmax+1;lp++)this->Rll_bin[l][lp]/=ch;
  }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::read_mixing_matrix_l(){
   vector< vector<real_prec> > pR;
   File.read_file(this->input_mixing_matrix_l,pR);

   int ik=-1;
   for(int i=0;i<=this->Lmax;i++){
     for(int j=0;j<=this->Lmax;++j){
       ++ik;
       this->R[i][j]=pR[ik][2];
     }
   }
  pR.clear();

  // Normalize mixing matrix

  cout<<"Done"<<endl;


}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::get_Cl_bin_convolved(){
  for(int l=0;l<this->lbin.size();l++)this->Clbin[l]=0;
  for(int l=0;l<this->lbin.size();l++)
    for(int lp=this->Lmin;lp<this->Lmax+1;lp++)
      this->Clbin[l]+=this->Rll_bin[l][lp]*this->Clvec[lp];
#ifdef _WRTITE_OUTPUTS_APSTH_
  File.write_to_file(output_clth_file_bin,lbin,Clbin);
#endif
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

// This is to write models varying two parameters
// specifically omega matter and bias
void AngularPowerTH::get_Cl_bin_convolved(int i, int j){
  for(int l=0;l<this->lbin.size();l++)this->Clbin[l]=0;

  for(int l=0;l<this->lbin.size();l++)
    for(int lp=this->Lmin;lp<this->Lmax+1;lp++)
    this->Clbin[l]+=this->Rll_bin[l][lp]*this->Clvec[lp];
#ifdef _WRTITE_OUTPUTS_APSTH_
  string of=output_clth_file_bin+"_om_"+to_string(i)+"_sbias_"+to_string(j);
  File.write_to_file(of,lbin,Clbin);
#endif
}

// This is to write models varying one parameters,
// specifically the width of the pdf for photoz
void AngularPowerTH::get_Cl_bin_convolved(int i){
  for(int l=0;l<this->lbin.size();l++)this->Clbin[l]=0;
  for(int l=0;l<this->lbin.size();l++)for(int lp=this->Lmin;lp<this->Lmax+1;lp++)this->Clbin[l]+=this->Rll_bin[l][lp]*this->Clvec[lp];
#ifdef _WRTITE_OUTPUTS_APSTH_
  string of=output_clth_file_bin+"_sigmaz_"+to_string(i);
  File.write_to_file(of,lbin,Clbin);
#endif
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerTH::free_gsl_table(){
  gsl_integration_glfixed_table_free(this->wf);
  gsl_integration_glfixed_table_free(this->wfd);
}
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_Cl_theo(params_clth sp, string fit_type){

    /* Fit type is a parameter used for mcmc purposes. If we want to fit
     * to cosmolgical parameters, then we need to cmpute every time the P(k) in the chains
     * and hence this function provides the access to cosmo-par dependent function.
     * Inwe want to fit to shape parameters, then we need to compute the P(k) only once and then
     * interpolate. So basically, this function access the same P(k) and in one case generates the Cl with
     * the full function or with the interpolation, which HAS TO BE DONE BEFORE THE CALL
     * If this function. This is basially suitable for RPT power spectrum, for there the
     * factorization with respect to time and k is possible. */

    omp_set_num_threads(_NTHREADS_);
    time_t start_all;
    time(&start_all);
    time_t end;
    real_prec lapse;


    if(this->code!="cl_mcmc")
      std::cout<<BLUE<<"Theoretical Cl using "<<omp_get_max_threads()<<" threads"<<RESET<<std::endl;
    
    // ------------------------------------------------------------------------------------------
    real_prec norm=(this->pdf_zerrors=="gaussian" || this->pdf_zerrors=="lorentz")? pow(this->normal_cl_power,2): pow(this->normal_cl_power_photo,2);
    norm*=Constants::speed_light;

    // ------------------------------------------------------------------------------------------
    // We calculate the Cl every 'step', then interpolate.
    int step=2;
    int Nll = (int)((this->Lmax)/step) + 1;
    int Nll_break = (int)floor(this->L_break/step);
    vector<gsl_real>Cl_two(Nll,0);
    vector<gsl_real>vl_two(Nll,0);
    for(int l=0;l<vl_two.size();++l)vl_two[l]=step*l;
    std::vector<real_prec>F(this->number_of_kmodes,0);
    // ------------------------------------------------------------------------------------------
    if(fit_type=="cosmo")
      {
        time(&start_all);
      // -------------------------------------------
      /* Creating objects type PowerSpectrum that initates integration material only once
       *  This is donde only if we are to use the fit_type comso */
        PowerSpectrum PS;
        PowerSpectrum Pst_hm(this->gal.scp.kmin_int, this->gal.scp.kmax_int,40,10,this->gal.scp.M_min_effective, this->gal.scp.M_max_effective, this->gal.scp.n_points_mass_integration);
        PowerSpectrum Pst_pt(this->gal.scp.kmin_int, this->gal.scp.kmax_int,20,10);
        Pst_hm.v_mass.resize(this->gal.scp.n_points_mass,0);
        Pst_hm.v_sigma_mass.resize(this->gal.scp.n_points_mass,0);
        Pst_hm.v_mass_function.resize(this->gal.scp.n_points_mass,0);
        Pst_hm.v_halo_mass_bias.resize(this->gal.scp.n_points_mass,0);
        if(!this->use_non_linear_pk)
          PS=this->Ps;
        if(this->use_non_linear_pk && this->type_of_nl_power=="hf")
          PS=this->Ps;
        if(this->use_non_linear_pk && this->type_of_nl_power=="hm")
          PS=Pst_hm;
        if(this->use_non_linear_pk && this->type_of_nl_power=="pt")
          PS=Pst_pt;
        if(this->code!="cl_mcmc" && this->L_break>0)
          cout<<BLUE<<"Cl using exact expression up to l ="<<this->L_break<<RESET<<std::endl;

     // let us create MASS_function[M][z] Aand MASS_BIAS_MZ
      // as an object ot type PPOWER SPECTRUM, using the same XX_MASS and XX_z
      // such that w can mass it to the Cl already to be used,
      this->gal.scp.gv=sp.gv;
      Pst_hm.mass_function_M_Z(this->XX_z,&this->gal.scp);
      PS=Pst_hm;

      int ik;
#ifdef _USE_OMP_
#pragma omp parallel for private(ik)
#endif
      for(int l=0;l<Nll_break;l++)
      {
        s_aux <PowerSpectrum> ssa;
        ssa.Ps=PS;
        ssa.s_clth=&sp;
        ssa.scp_a=&this->gal.scp;
        ssa.XX_z=this->XX_z;
        ssa.laux=step*l;
        fill(F.begin(), F.end(),0);
        for(ik=0;ik<this->k_th.size();++ik)
         {
           ssa.kaux=this->k_th[ik];
           F[ik]=this->get_Kernel((void*)&ssa);
         }
        sp.Fkernel=F;
        sp.l=step*l;
        Cl_two[l]= pow(this->pixel_window[step*l],2)*this->get_cl_exact((void *)&sp)*Constants::speed_light/norm;
      }
      time(&end);
      lapse=difftime(end,start_all);
      if(this->code!="cl_mcmc")std::cout<<BLUE<<"Time elapsed for Cl-exact: "<<lapse<<"  secs \r"<<RESET<<std::endl;

      time(&start_all);
      if(this->code!="cl_mcmc")std::cout<<BLUE<<"Cl using Limber approximation for l>= "<<this->L_break<<RESET<<std::endl;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int l=Nll_break;l<Nll;l++){
          s_aux<PowerSpectrum>  ssa;
          ssa.Ps=PS;
          ssa.s_clth=&sp;
          ssa.XX_z=this->XX_z;
          ssa.scp_a=&this->gal.scp;
          int l_i= step*l;
          ssa.laux=l_i;
          Cl_two[l]=pow(this->pixel_window[l_i],2)*this->get_cl_limber((void*)&ssa)/norm;
      }

      time(&end);
      lapse=difftime(end,start_all);
      if(this->code!="cl_mcmc")std::cout<<BLUE<<"Time elapsed for Cl-Limber: "<<lapse<<"  secs \r"<<RESET<<std::endl;

     // Now interpolate the values obtained every steps
      Clvec[0]=Cl_two[0];
      for(int l=this->Lmin+1;l<=this->Lmax;++l)Clvec[l]=gsl_inter_new(vl_two, Cl_two, l);
      for(int l=0;l<=this->Lmax;++l)this->lvec[l]=l;
#ifdef _WRTITE_OUTPUTS_APSTH_
      if(this->code!="cl_mcmc")
        this->File.write_to_file(this->output_clth_file, this->lvec, this->Clvec);
#endif
  }

        // WARNING, SO FAR ONLY FOR LIMBER THIS PART. THIS USES INTERPOLATED FUNCTIONS COMPUTED IN get-pk
       // THAT ARE USEFUL ONLY FOR THE RPT POWER SPECTRUM WHEN CONSTRAINGINH b, AMC, kstar.
    else if(fit_type=="shape")
    {
        time(&start_all);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int l=Nll_break;l<Nll;l++){
          s_aux<PowerSpectrum>  ssa;
          ssa.Ps=this->Ps;
          ssa.s_clth=&sp;
          ssa.scp_a=&this->gal.scp;
          int l_i= step*l;
          ssa.laux=l_i;
          Cl_two[l]=pow(this->pixel_window[l_i],2)*this->get_cl_limber_pt_inter((void*)&ssa)/norm;
      }

      time(&end);
      lapse=difftime(end,start_all);
      if(this->code!="cl_mcmc")std::cout<<BLUE<<"Time elapsed for Cl-Limber: "<<lapse<<"  secs \r"<<RESET<<std::endl;

     // Now interpolate the values obtained every steps
      Clvec[0]=Cl_two[0];
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int l=this->Lmin+1;l<=this->Lmax;++l)
        Clvec[l]=gsl_inter_new(vl_two, Cl_two, l);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int l=0;l<=this->Lmax;++l)
        this->lvec[l]=l;

#ifdef _WRTITE_OUTPUTS_APSTH_
      if(this->code!="cl_mcmc")
        this->File.write_to_file(this->output_clth_file, this->lvec, this->Clvec);
#endif
    }



}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerTH::get_Cl_theo(){

  So.enter(__PRETTY_FUNCTION__);


  omp_set_num_threads(_NTHREADS_);
  time_t start_all;
  time(&start_all);
  time_t end;
  real_prec lapse;



    this->get_normal((void *)&this->pclth);
    // ------------------------------------------------------------------------------------------
    real_prec norm=(this->pdf_zerrors=="gaussian" || this->pdf_zerrors=="lorentz")? pow(this->normal_cl_power,2): pow(this->normal_cl_power_photo,2);
    norm*=Constants::speed_light;
    vector<real_prec>F(this->number_of_kmodes,0);
     // -------------------------------------------
    /* Creating objects type PowerSpectrum that initates integration material only once
     *  This is donde only if we are to use the fit_type cosmo */

    PowerSpectrum PS;
    if(true==use_non_linear_pk)
      {
     PowerSpectrum Pst_hm(this->gal.scp.kmin_int, this->gal.scp.kmax_int,40,10,this->gal.scp.M_min_effective, this->gal.scp.M_max_effective, this->gal.scp.n_points_mass_integration);
     PowerSpectrum Pst_pt(this->gal.scp.kmin_int, this->gal.scp.kmax_int,20,10);
     Pst_hm.v_mass.resize(this->gal.scp.n_points_mass,0);
     Pst_hm.v_sigma_mass.resize(this->gal.scp.n_points_mass,0);
     Pst_hm.v_mass_function.resize(this->gal.scp.n_points_mass,0);
     Pst_hm.v_halo_mass_bias.resize(this->gal.scp.n_points_mass,0);
     if(!this->use_non_linear_pk)
       PS=this->Ps;
     if(this->use_non_linear_pk && this->type_of_nl_power=="hf")
       PS=this->Ps;
     if(this->use_non_linear_pk && this->type_of_nl_power=="hm")
       PS=Pst_hm;
     if(this->use_non_linear_pk && this->type_of_nl_power=="pt")
       PS=Pst_pt;
     if(this->code!="cl_mcmc" && this->L_break>0)
       cout<<BLUE<<"Cl using exact expression up to l ="<<this->L_break<<RESET<<std::endl;

     // let us create MASS_function[M][z] Aand MASS_BIAS_MZ
      // as an object ot type PPOWER SPECTRUM, using the same XX_MASS and XX_z
      // such that w can mass it to the Cl already to be used,
      this->gal.scp.gv=this->pclth.gv;
      Pst_hm.mass_function_M_Z(this->XX_z,&this->gal.scp);
      PS=Pst_hm;

      int ik;
      if(this->L_break>0)
      {

#ifdef _USE_OMP_
#pragma omp parallel for private(ik)
#endif
      for(int l=0;l<this->L_break;l++)
      {
        s_aux <PowerSpectrum> ssa;
        ssa.Ps=PS;
        ssa.s_clth=&pclth;
        ssa.scp_a=&this->gal.scp;
        ssa.XX_z=this->XX_z;
        ssa.laux=l;
        fill(F.begin(), F.end(),0);
        for(ik=0;ik<this->k_th.size();++ik)
         {
           ssa.kaux=this->k_th[ik];
           F[ik]=this->get_Kernel((void*)&ssa);
         }
        this->pclth.Fkernel=F;
        this->pclth.l=l;
        this->Clvec[l]= this->get_cl_exact((void *)&pclth)*Constants::speed_light/norm;
      }
    }

    for(int l=this->L_break;l<this->Lmax;l++)
      {
          s_aux<PowerSpectrum>  ssa;
          ssa.Ps=PS;
          ssa.s_clth=&pclth;
          ssa.XX_z=this->XX_z;
          ssa.scp_a=&this->gal.scp;
          int l_i= l;
          ssa.laux=l_i;
          this->Clvec[l]=this->get_cl_limber((void*)&ssa)/norm;
          this->lvec[l]=l; 
       }
  }
  else
  {

    for(int l=this->Lmin;l<this->Lmax;l++)
      {
          s_aux<PowerSpectrum>  ssa;
          ssa.Ps=this->Ps;
          ssa.s_clth=&this->pclth;
          ssa.XX_z=this->XX_z;
          ssa.scp_a=&this->gal.scp;
          ssa.laux=l;
          this->Clvec[l]=this->get_cl_limber((void*)&ssa)/norm;
          this->lvec[l]=l; 
       }
  }



Gnuplot gp_pdf;
 vector<pair<real_prec, real_prec> > xy_pts;
 for(int i=0;i<Clvec.size();++i)
   xy_pts.push_back(std::make_pair(lvec[i], Clvec[i])); 

 gp_pdf<<"set border linewidth 1.5\n";
 gp_pdf<<"set grid\n";
 gp_pdf << "set log \n";
 gp_pdf << "set xlabel 'l'\n";
 gp_pdf << "set ylabel 'C_l'\n";
 gp_pdf << "plot[][1e-8: 1e-4]" << gp_pdf.file1d(xy_pts) << "w l lw 2 title 'Cl'"<<endl;


#ifdef _WRTITE_OUTPUTS_APSTH_
        this->File.write_to_file(this->output_clth_file, this->lvec, this->Clvec);
#endif

  

}



