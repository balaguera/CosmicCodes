////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** @class<Galaxy>
 *  @file Galaxy.cpp
 *  @brief Methods of the class Galaxy
 *  @details Reads and administrates input parameters
 *  @author Andrés Balaguera-Antolínez,
 *  @date 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../headers/Galaxy.h"
#define EPSILON_MK 1.0

//#define _USE_PDF_PHOT2SPEC_   //tjis is to speed up the code, ut has to be given in accordance to parameter file
#undef _USE_PDF_PHOT2SPEC_   //this is to speed up the code, ut has to be given in accordance to parameter file

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::read_pars(string &file){

  cout<<CYAN<<endl;
  std::cout<<"*******************************************************"<<endl;
  cout<<CYAN<<"Reading parameter file"<<RESET<<endl;
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
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
      }

      // read parameter value
      line_string >> par_value;

      if (par_value.empty()) {
        cout << "Value of " << par_name << " not specified in " <<file << endl;
        cout << "Assuming a default value for " << par_name << endl;
        continue;
      }
      if (par_name == "name_cat")this->name_cat = par_value;
      if (par_name == "name_input_file_cat")this->name_input_file_cat = par_value;
      else if (par_name == "name_file_mask")this->name_file_mask = par_value;
      else if (par_name == "name_file_mask_north_gal")this->name_file_mask_north_gal = par_value;
      else if (par_name == "name_file_mask_south_gal")this->name_file_mask_south_gal = par_value;
      else if (par_name == "name_file_Kcorr")this->name_file_Kcorr = par_value;
      else if (par_name == "i_z_Kcorr")this->i_z_Kcorr = atoi(par_value.c_str());
      else if (par_name == "i_flux_factor_Kcorr")this->i_flux_factor_Kcorr = atoi(par_value.c_str());
      else if (par_name == "i_color_correction Kcorr")this->i_color_correction_Kcorr = atoi(par_value.c_str());
      else if (par_name == "i_ra")this->i_ra = atoi(par_value.c_str());
      else if (par_name == "i_dec")this->i_dec = atoi(par_value.c_str());
      else if (par_name == "i_mass")this->i_mass = atoi(par_value.c_str());
      else if (par_name == "i_lgal")this->i_lgal = atoi(par_value.c_str());
      else if (par_name == "i_bgal")this->i_bgal = atoi(par_value.c_str());
      else if (par_name == "i_mJ")this->i_mJ = atoi(par_value.c_str());
      else if (par_name == "i_mH")this->i_mH = atoi(par_value.c_str());
      else if (par_name == "i_mK")this->i_mK = atoi(par_value.c_str());
      else if (par_name == "i_EBV")this->i_EBV = atoi(par_value.c_str());
      else if (par_name == "i_lsd")this->i_lsd = atoi(par_value.c_str());
      else if (par_name == "i_zs")this->i_zs = atoi(par_value.c_str());
      else if (par_name == "i_zp")this->i_zp = atoi(par_value.c_str());
      else if (par_name == "i_Kcsb")this->i_Kcsb = atoi(par_value.c_str());
      else if (par_name == "redshift_type")this->redshift_type = par_value;
      else if (par_name == "hemisphere")this->hemisphere = par_value;

      else if (par_name == "type_of_K_correction")this->type_of_K_correction = par_value;

      
      else if (par_name == "Redefine_MKminmax"){
	if(par_value=="true")this->Redefine_MKminmax=true;
        else if(par_value=="false") this->Redefine_MKminmax=false;
      }


      

      else if (par_name == "Get_new_catalog"){
        if(par_value=="true")this->Get_new_catalog=true;
        else if(par_value=="false") this->Get_new_catalog=false;
      }
      
      else if (par_name == "Measure_LF"){
        if(par_value=="true")this->Measure_LF =true;
        else if(par_value=="false") this->Measure_LF=false;
      }

      else if (par_name == "use_K_correction"){
        if(par_value=="true")this->use_K_correction =true;
        else if(par_value=="false") this->use_K_correction=false;
      }


      else if (par_name == "use_e_correction"){
        if(par_value=="true")this->use_e_correction =true;
        else if(par_value=="false") this->use_e_correction=false;
      }
      else if (par_name == "Measure_M_pdf"){
        if(par_value=="true")this->Measure_M_pdf =true;
        else if(par_value=="false")this->Measure_M_pdf=false;
      }

      else if (par_name == "observed_pixels_in_mask")this->observed_pixels_in_mask  = atoi(par_value.c_str());
      
      else if (par_name == "i_mask_pixel")i_mask_pixel = atoi(par_value.c_str());
      else if (par_name == "i_lpix")i_lpix = atoi(par_value.c_str());
      else if (par_name == "i_bpix")i_bpix = atoi(par_value.c_str());
      else if (par_name == "i_mask_flag")i_mask_flag = atoi(par_value.c_str());
      else if (par_name == "Healpix_resolution")Healpix_resolution = atoi(par_value.c_str());
      else if (par_name == "N_iterations")N_iterations = atoi(par_value.c_str());

      else if (par_name == "N_bin_color")N_bin_color = atoi(par_value.c_str());
      else if (par_name == "N_bin_lmass")N_bin_lmass = atoi(par_value.c_str());
      else if (par_name == "N_bin_lmass_lowres")N_bin_lmass_lowres = atoi(par_value.c_str());
      else if (par_name == "lmass_min")lmass_min = atoi(par_value.c_str());
      else if (par_name == "lmass_max")lmass_max = atoi(par_value.c_str());

      else if (par_name == "N_bin_Mag")N_bin_Mag = atoi(par_value.c_str());
      else if (par_name == "N_bin_Mag_low_res")N_bin_Mag_low_res = atoi(par_value.c_str());

      else if (par_name == "N_bin_mag")N_bin_mag = atoi(par_value.c_str());

      else if (par_name == "LF_estimator")this->LF_estimator= par_value;


      else if (par_name == "N_bin_z")N_bin_z = atoi(par_value.c_str());
      else if (par_name == "N_bin_z_low_res")N_bin_z_low_res = atoi(par_value.c_str());
      else if (par_name == "z_min")z_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "z_max")z_max = static_cast<real_prec>(atof(par_value.c_str()));
      
      else if (par_name == "z_min_low_res")z_min_low_res = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "z_max_low_res")z_max_low_res = static_cast<real_prec>(atof(par_value.c_str()));

      else if (par_name == "kcorr_index")this->kcorr_index= static_cast<real_prec>(atof(par_value.c_str()));
      
      
      else if (par_name == "mK_min")mK_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "mK_max")mK_max = static_cast<real_prec>(atof(par_value.c_str()));

      else if (par_name == "MK_min")MK_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "MK_max")MK_max = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "color_min")color_min = static_cast<real_prec>(atof(par_value.c_str()));
      else if (par_name == "color_max")color_max = static_cast<real_prec>(atof(par_value.c_str()));

      else if (par_name == "N_z_cosmo")N_z_cosmo = atoi(par_value.c_str());

      // Get cosmological parameters from parameter file and put them in a structure
      // of type s_cosmological_parameters




      else if (par_name == "Om_matter")this->Om_matter = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_cdm")this->Om_cdm = (real_prec)atof(par_value.c_str());

      else if (par_name == "Om_radiation")this->Om_radiation = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_baryons")this->Om_baryons = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_vac")this->Om_vac = (real_prec)atof(par_value.c_str());
      else if (par_name == "Om_k")this->Om_k = (real_prec)atof(par_value.c_str());
      else if (par_name == "Hubble")this->Hubble = (real_prec)atof(par_value.c_str());
      else if (par_name == "hubble")this->hubble = (real_prec)atof(par_value.c_str());
      else if (par_name == "n_s")this->n_s = (real_prec)atof(par_value.c_str());
      else if (par_name == "alpha_s")this->alpha_s = (real_prec)atof(par_value.c_str());
      else if (par_name == "w_eos")this->w_eos = (real_prec)atof(par_value.c_str());
      else if (par_name == "N_eff")this->N_eff = (real_prec)atof(par_value.c_str());
      else if (par_name == "sigma8")this->sigma8 = (real_prec)atof(par_value.c_str());
      else if (par_name == "Tcmb")this->Tcmb = (real_prec)atof(par_value.c_str());
      else if (par_name == "GAL_BIAS")this->GAL_BIAS = (real_prec)atof(par_value.c_str());
      else if (par_name == "alpha_BIAS")this->alpha_BIAS = (real_prec)atof(par_value.c_str());
      else if (par_name == "kstar")this->kstar = (real_prec)atof(par_value.c_str());


      // this is to read magnitude limits

      else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'x'))
	{
	  stringstream line_string (line_in_file);
	  line_string << line_in_file;
	  line_string >> par_name;   // read the first character, the name of the parameter
	  line_string >> equality; 	  // check that second character is "="
	  if (equality != "=") {
	    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	    cerr << "Using a default value for " << par_name << endl; exit(1);
	  }

	  
	  string par_value0;
	  int ending=0;
	  int n_m_lims=0;
	  do{
	    line_string >> par_value0;  //read value
	    n_m_lims++;
	    if(par_value0=="END")ending=1;
	    this->mK_limits.push_back(atof(par_value0.c_str()));
	  }while(ending!=1);
	  this->mK_limits.pop_back(); //delete the last element, which is END
	}
      
      
      else  if ((line_in_file[0] != '#' && line_in_file.empty()==0) && (line_in_file[0] == 'Z'))
	{
	  stringstream line_string (line_in_file);
	  line_string << line_in_file;
	  line_string >> par_name;   // read the first character, the name of the parameter
	  line_string >> equality; 	  // check that second character is "="
	  if (equality != "=") {
	    cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	    cerr << "Using a default value for " << par_name << endl; exit(1);
	  }
	  
	  string par_value0;
	  int ending=0;
	  int n_m_lims=0;
      vector<real_prec>z_limits_aux;
	  do{
	    line_string >> par_value0;  //read value
	    n_m_lims++;
	    if(par_value0=="END")ending=1;
	    z_limits_aux.push_back(atof(par_value0.c_str()));
	  }while(ending!=1);
	  z_limits_aux.pop_back(); //delete the last element, which is END

	  this->z_limits.resize(z_limits_aux.size(),0);
	  for (int iz=0;iz<z_limits_aux.size();++iz)
	    this->z_limits[iz]=z_limits_aux[iz];
	  

	}
      
      
    }
  }


  // Cosmological parameters
#ifdef _USE_COSMO_PARS_
  this->Om_matter = COSMOPARS::Om_matter;
  this->Om_radiation = COSMOPARS::Om_radiation;
  this->Om_baryons = COSMOPARS::Om_baryons;
  this->Om_cdm = this->Om_matter-this->Om_baryons;
  this->Om_vac = COSMOPARS::Om_vac;
  this->Om_k   = COSMOPARS::Om_k;
  this->Hubble = COSMOPARS::Hubble;
  this->hubble = COSMOPARS::hubble;
  this->n_s = COSMOPARS::n_s;
  this->w_eos = COSMOPARS::w_eos;
  this->N_eff =  COSMOPARS::N_eff;
  this->sigma8 = COSMOPARS::sigma8;
  this->Tcmb = COSMOPARS::Tcmb;
  this->use_wiggles = true;
  this->RR = COSMOPARS::RR;
  this->alpha_s=COSMOPARS::alpha_s;
#endif
  this->scp.Om_matter=this->Om_matter;
  this->scp.Om_cdm=this->Om_cdm;
  this->scp.Om_radiation=this->Om_radiation;
  this->scp.Om_baryons=this->Om_baryons;
  this->scp.Om_k=this->Om_k;
  this->scp.Om_vac=1.- this->scp.Om_matter-this->scp.Om_radiation- this->scp.Om_k;
  this->scp.Hubble=this->Hubble;
  this->scp.hubble=this->hubble;
  this->scp.n_s=this->n_s;
  this->scp.alpha_s=this->alpha_s;
  this->scp.w_eos=this->w_eos;
  this->scp.N_eff=this->N_eff;
  this->scp.sigma8=this->sigma8;
  this->scp.Tcmb=this->Tcmb;
  this->scp.mlim=this->mK_max;
  this->scp.f_baryon = this->Om_baryons/(this->Om_matter-this->Om_baryons);  // This is the one required by the EH Transfer function
  this->scp.use_wiggles=1;
  this->scp.RR=8.0;
  this->scp.GAL_BIAS=this->GAL_BIAS;
  this->scp.alpha_BIAS=this->alpha_BIAS;
  this->scp.kstar=this->kstar;
  this->scp.use_K_correction=this->use_K_correction;
  this->scp.use_e_correction=this->use_e_correction;

  this->cosmology.set_cosmo_pars(this->scp);
  
  if(this->redshift_type=="ps" || this->redshift_type=="s" || this->redshift_type=="ps_s")
    this->zzns=this->z_min;
  else
    if(this->redshift_type=="p" || this->redshift_type=="p_s")
      this->zzns=-2000.0;
  
  So.DONE();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::set_file_names(){

    So.message_screen("Setting file names");
// tHIS IS TO SELECT WHETHER WE WORK WITH S, P or objects WITH BOTH S AND P redshifts
  // en este caso no importa si hay zs o no
  //this->zzns=this->z_min;

  string nmask="";
#ifdef _USE_JK_
  nmask="_JK"+to_string(this->index_mask);
#endif
  // ********************************

  // DNDZ MEASURED

  string output_dir = this->output_dir;

  this->output_file_dndz=output_dir+this->name_cat+"_dndz_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_"+this->hemisphere+nmask+".txt";
  
  // ********************************
  // DNDZ FROM LF
  this->output_file_dndz_lf=output_dir+this->name_cat+"_dndz_LF_"+this->LF_estimator+"_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // DNDZ FROM dN/dmk
  this->output_file_dndmk=output_dir+this->name_cat+"_dndmk_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // SMOOTHED DNDZ
  this->output_file_dndz_smooth =output_dir+this->name_cat+"_dndz_smooth_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // SMOOTHED NBAR
  
  this->output_file_nbar_smooth =output_dir+this->name_cat+"_nbar_smooth_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // LF
  this->output_file_LF=output_dir+this->name_cat+"_LF_"+this->LF_estimator+"_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zbin"+to_string(this->i_z_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // mAGNITUD DISTRIBUTION
  this->output_file_NM=output_dir+this->name_cat+"_NM_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zbin"+to_string(this->i_z_limit)+"_"+this->hemisphere+nmask+".txt";
  // ********************************
  // NEW CAT
  this->output_file_new_cat=output_dir+this->name_cat+"_NewCat_"+this->redshift_type+"_mk"+to_string(this->i_magnitude_limit)+"_zbin"+to_string(this->i_z_limit)+nmask+".txt";
  
  So.DONE();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_cosmo(int nz, bool silence){

  time_t start_all;
  time(&start_all);
  time_t end;
  if(silence)So.message_screen("Getting cosmological functions...");

  omp_set_num_threads(_NTHREADS_);

  // Default values for E-correction  e(z)=a*z*z + b*z+c
  this->scp.K_index_a= 0;
  this->scp.K_index_b= 3.02;
  this->scp.K_index_c= 0;
  this->scp.K_index_d= 0;
  this->scp.K_index_e= 0;
  this->scp.e_index_zstar= 0.00;
  this->scp.k_index=this->kcorr_index;
  
  //Vector for comovong distance
  this->rv.resize(nz,0);
  //Vector for the transverse comoving distance
  this->trv.resize(nz,0);
  //Vector for the redshifts
  this->zv.resize(nz,0);
  //Vector for the growth factor D(z) normalized to 1 at z=0, satisfying D(a)=a for a matter dom. time.
  this->gv.resize(nz,0);
  // Vector for the redshift dependent bias
  this->bias_zv.resize(nz,0);
  //Vector for the growth index f(z)
  this->gfv.resize(nz,0);
  //Vector for the Hubble parameter
  this->Hv.resize(nz,0);
  //Vector for the Distance modulus
  this->Dm.resize(nz,0);

  real_prec  g_today = cosmology.growth_factor(0.0);


  real_prec  zmin_inter=0.001;
  real_prec  zmax_inter=1.1;
 #pragma omp parallel for
  for(int i=0;i<nz;++i)
      this->zv[i]=zmin_inter+i*(zmax_inter-zmin_inter)/static_cast<real_prec>(nz-1);
  


// #pragma omp parallel for
  for(int i=0;i<nz;++i)
    {
      this->Hv[i]=cosmology.Hubble_function(this->zv[i]);
      this->rv[i]=cosmology.comoving_distance(this->zv[i]);
//      this->gv[i]=cosmology.growth_factor(this->zv[i])/g_today;
//      this->bias_zv[i]=this->scp.GAL_BIAS*pow(1+this->zv[i], this->scp.alpha_BIAS);
//      this->gfv[i]=cosmology.growth_index(this->zv[i]);
//      this->Dm[i]=cosmology.Distance_Modulus(this->zv[i]);
//       cout<<rv[i]<<"  "<<zv[i]<<endl;
      }
  
  // ALLOCATE VECTORS TO STRUCTURE cosmological parameters TO USE THEM IN OTHER OPERATIONS
  this->scp.zv=this->zv;
  this->scp.rv=this->rv;
  
  // HERE THESE ARE COMPUTED FROM ALREADY TABULATED QUANTITIES

/*
#pragma omp parallel for
  for(int i=0;i<nz;++i)
    this->trv[i]=cosmology.inter_transverse_comoving_distance(this->zv[i]);
  
  this->scp.trv=this->trv;
*/

  if(silence){
    time(&end);
    real_prec  lapse=difftime(end,start_all);
    std::cout<<BLUE<<"Done in "<<lapse<<" secs "<<RESET<<std::endl;
    time(&start_all);
  }


}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_mask(){

  this->mask.clear();
/* ANDRES WARNING
    if(this->hemisphere=="all")
      this->NMASK=this->File.read_file(name_file_mask,this->mask, omp_get_max_threads());
  else if(this->hemisphere=="south")
      this->NMASK=this->File.read_file(name_file_mask_south_gal,this->mask,omp_get_max_threads());
  else if(this->hemisphere=="north")
      this->NMASK=this->File.read_file(name_file_mask_north_gal,this->mask,omp_get_max_threads());
*/

  this->NCOLS_MASK=static_cast<unsigned long>(this->mask.size()/NMASK);
  
  unsigned int pix_count=0;

#pragma omp parallel for reduction(+:pix_count)
  for(int i=0;i<this->NMASK;++i)
    if(this->mask[this->i_mask_flag+i*NCOLS_MASK]==this->observed_pixels_in_mask)
      pix_count++;
  
  this->n_pixels=this->NMASK;
  this->nside=static_cast<int>(sqrt((n_pixels/12)));
  this->area_pixel=4.*M_PI/n_pixels;

  //  With the mask we count the number of pixels
  //  taken into the sample and multiply them
  //  by the area of each pixesl This will give us the
  //  total area of the survey, to be used in the Vmax.
  //  We need to use the mask,since the catalogue read
  //  in this code is not yet masked.
  //  The pixels of the mask are defined in Galactic coordinates.
  this->total_area=static_cast<real_prec>(pix_count*area_pixel);
  So.message_screen("Total area covered by the mask = ", this->total_area/(3.0462e-4)," squared degrees");
  So.message_screen("fsky = ", this->total_area/(3.0462e-4)/(129600.0/M_PI));

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_mask_fits(string file){
  cout<<CYAN<<"Reading mask from "<<file<<endl;
  cout<<RESET<<endl;
  Healpix_Map<real_prec>mask_aux(log2(this->nside), RING);
  //  read_Healpix_map_from_fits(file,mask_aux, 1,2);
  this->nside=mask_aux.Nside();
  n_pixels=12*this->nside*this->nside;
  this->mask.clear();
  this->NCOLS_MASK=4;
  mask.resize(n_pixels*this->NCOLS_MASK,0);

#pragma omp parallel for
  for(int i=0;i<n_pixels;i++)
    mask[this->i_mask_flag+i*this->NCOLS_MASK]=mask_aux[i];

  /*
    int ip=0;
    for(int ir=0;ir<nr;ir++){
    int Npixels_ring=npix_ring(nside, ir);
    for(int ipr=0;ipr<Npixels_ring;++ipr){
    point_rev=mask_aux.pix2ang(ip);
    this->theta_new[ir]=point_rev.theta;
    this->phi_new[ip]=point_rev.phi;
    }
    ip++;
    }

  */
  pointing point_rev;
  ofstream mas;
  mas.open("../../NVSS/DATA/MASK_NVSS.dat");
  for(int i=0;i<12*nside*nside;i++){
    point_rev=mask_aux.pix2ang(i);
    mas<<i<<"\t"<<point_rev.phi*180./M_PI<<"\t"<<90.-point_rev.theta*180./M_PI<<"\t"<<mask[this->i_mask_flag+i*this->NCOLS_MASK]<<endl;
  }
  mas.close();
  So.DONE();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::read_input_cats(string cat)
{

  
  this->prop.clear();
  this->prop.shrink_to_fit();

  if(cat=="g")
    {
      this->NGAL=this->File.read_file(this->name_input_file_cat, this->prop,_NTHREADS_);
      this->NCOLS=this->prop.size()/this->NGAL;
    }
  else
    if(cat=="r"){
      this->prop_r.clear();
      this->NRAN=this->File.read_file(this->name_input_file_cat, this->prop_r,_NTHREADS_);
      this->NCOLS_R=this->prop_r.size()/this->NRAN;
    }

#ifdef _USE_HEALPIX_
  this->set_catalog(cat);
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::set_catalog(string cat){
  if(cat=="g")
    {
      So.message_screen("Adding information of the mask to each galaxy ...");
      Healpix_Map<real_prec>map_aux(log2(this->nside), RING);

      gal_mask.resize(this->NGAL);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->NGAL;++i)
	{
	  this->point.phi=fac*this->prop[this->i_lgal+i*this->NCOLS];
	  this->point.theta=0.5*M_PI-fac*this->prop[this->i_bgal+i*this->NCOLS];
	  unsigned long ipix=map_aux.ang2pix(this->point);
	  //gal_mask[i]=  this->prop[this->i_bgal+i*this->NCOLS]<=60.0 ? 0 :   this->mask[this->i_mask_flag+ipix*this->NCOLS_MASK];  //WARNING, LOOK UP MACIEJ'S CONVENTION
	  gal_mask[i]=  this->mask[this->i_mask_flag+ipix*this->NCOLS_MASK];  //WARNING, LOOK UP MACIEJ'S CONVENTION
	}
    }
  
  else if(cat=="r")
    {
      Healpix_Map<real_prec>map_aux(log2(this->nside), RING);
      So.message_screen("Adding information of the mask to each galaxy ...");
      gal_mask.resize(this->NRAN);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->NRAN;i++)
	{
	  this->point.phi=fac*this->prop_r[this->i_lgal+i*this->NCOLS_R];
	  this->point.theta=0.5*M_PI-fac*this->prop_r[this->i_bgal+i*this->NCOLS_R];
	  long ipix=map_aux.ang2pix(this->point);
	  ran_mask[i]=this->mask[this->i_mask_flag+ipix*NCOLS_MASK];
	}
    }
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::set_pars(int izvls){
  this->i_z = (this->redshift_type == "s" ? this->i_zs : this->i_zp);

  this->iz_vls=izvls;
  switch(iz_vls)
    {
    case(0):this->z_min_vls = 0.0001;break;
    case(1):this->z_min_vls = 0.01;break;
    case(2):this->z_min_vls = 0.015;break;
    case(3):this->z_min_vls = 0.02;break;
    case(4):this->z_min_vls = 0.04;break;
    }

  this->mK_max=this->mK_limits[this->i_magnitude_limit];
  
  this->deltaz=(this->z_max-this->z_min)/(real_prec)this->N_bin_z;
  this->deltaz_low_res=(this->z_max_low_res-this->z_min_low_res)/(real_prec)this->N_bin_z_low_res;

//  this->deltaz=(0.4-0.001)/(real_prec)this->N_bin_z;
//  this->deltaz_low_res=(0.4-0.001)/(real_prec)this->N_bin_z_low_res;

  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::set_vec(){
  So.message_screen("Initializing vectors");
  real_prec  zz_min=this->z_min;
  real_prec  zz_min_low_res=this->z_min;
  this->deltaz=(this->z_max-this->z_min)/(real_prec)this->N_bin_z;
  zn.clear();
  zn.resize(this->N_bin_z,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z;++i)
    this->zn[i]=zz_min+(i+0.5)*deltaz;

  zn_min.clear();
  zn_min.resize(this->N_bin_z,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z;++i)
    this->zn_min[i]=zz_min+(i)*deltaz;

  zn_max.clear();
  zn_max.resize(this->N_bin_z,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z;++i)
    this->zn_max[i]=zz_min+static_cast<real_prec>(i+1)*deltaz;

  zn_low_res.clear();
  zn_low_res.resize(this->N_bin_z_low_res,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z_low_res;i++)
    this->zn_low_res[i]=zz_min_low_res+(i+0.5)*deltaz_low_res;

  Vshell.clear();
  Vshell.resize(this->N_bin_z,0);

#pragma omp parallel for
  for(int i=0;i<this->N_bin_z;i++)
    {
      real_prec  dr_dz=Constants::speed_light/cosmology.Hubble_function(zn[i]);
      real_prec  distance=cosmology.comoving_distance(zn[i]); // gsl_inter_new(this->zv,this->rv,zn[i]);
      this->Vshell[i]=this->total_area*deltaz*pow(distance,2)*dr_dz;
  }


  Vshell_lowres.clear();
  Vshell_lowres.resize(this->N_bin_z_low_res,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z_low_res;i++)
    {
      real_prec  dr_dz=Constants::speed_light/cosmology.Hubble_function(zn_low_res[i]);
      real_prec  distance=gsl_inter_new(this->zv,this->rv,zn_low_res[i]);
      this->Vshell_lowres[i]=this->total_area*deltaz*pow(distance,2)*dr_dz;
    }

  //  this->Pc.resize(this->N_bin_color,0);  
  //  this->number_c.resize(this->N_bin_color,0);
  //this->number_cm.resize(this->N_bin_color);
  //for(int i=0;i<this->N_bin_color;++i)number_cm[i].resize(this->N_bin_Mag, 0);

  v_color.clear();
  v_color.resize(this->N_bin_color,0);
#pragma omp parallel for
  for(int i=0;i<v_color.size();++i)
    v_color[i]=this->color_min+(i+0.5)*(this->color_max-this->color_min)/(static_cast<real_prec>(this->N_bin_color));


  So.DONE();
}  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::set_nbar(vector<real_prec> zz,vector<real_prec>nnb,vector<real_prec>pp){
  zn_new.resize(zz.size());
  nbar_new.resize(nnb.size());
  prob.resize(pp.size());
  for(int i=0;i<zz.size();++i)zn_new[i]=zz[i];
  for(int i=0;i<nnb.size();++i)nbar_new[i]=nnb[i];
  for(int i=0;i<pp.size();++i)prob[i]=pp[i];
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::read_get_Kcorr()
{
    /* The K-correction, obtained from the File provided by Jarret 2014
   * in which we have the rest frame flux = F * observed flux, and
   * here we read the factor F. Converting to magnidutes,
   * we get m_rest  = -2.5log10(F) + m_observed
   * such that in terms of the absolute magnitude M
   * Mrest = mobs +DM(z) + K, where K = -2.5log10(F).
   * Hence, we provide here the K correction to be applied
   * to the absolute magnitude. We here not not have yet the
   * e-correction.
   * Color correction factor is Color_rest = Color_obs - Color correction
   * We therefore allocate two vectors v_z_Kcorr and v_Kcorr
   */
  
  vector<real_prec> propK;
  int nlinesK=this->File.read_file(this->name_file_Kcorr, propK,omp_get_max_threads());
  int ncolsK=propK.size()/nlinesK;
  this->v_z_Kcorr.resize(nlinesK,0);
  this->v_color_correction_Kcorr.resize(nlinesK,0);
  this->v_Kcorr.resize(nlinesK,0);

#pragma omp parallel for
  for(int i=0;i<nlinesK;++i){
    this->v_z_Kcorr[i]=propK[this->i_z_Kcorr+i*ncolsK];
    this->v_Kcorr[i]=-2.5*log10(propK[this->i_flux_factor_Kcorr+i*ncolsK]);
    this->v_color_correction_Kcorr[i]=propK[this->i_color_correction_Kcorr+i*ncolsK];
  }
  
  propK.clear();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec  GALAXY::get_mina(int ii){
  real_prec  ZX=10000;
  real_prec  zmx;
  for(int i=0;i<this->NGAL;i++){
    zmx=min(ZX,this->prop[ii+i*this->NCOLS]);
    ZX=zmx;
  }
  return zmx;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec  GALAXY::get_min_r(int ii){
  real_prec  ZX=10000;
  real_prec  zmx;
  for(int i=0;i<this->NRAN;i++){
    zmx=min(ZX,this->prop_r[ii+i*this->NCOLS_R]);
    ZX=zmx;
  }
  return zmx;
}////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec  GALAXY::get_maxa(int ii ){
  real_prec  ZX=-100;
  real_prec  zmx;
  for(int i=0;i<this->NGAL;i++){
    zmx=max(ZX,this->prop[ii+i*NCOLS]);
    ZX=zmx;
  }
  return zmx;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec  GALAXY::get_max_r(int ii ){
  real_prec  ZX=-100;
  real_prec  zmx;
  for(int i=0;i<this->NRAN;i++){
    zmx=max(ZX,this->prop_r[ii+i*this->NCOLS_R]);
    ZX=zmx;
  }
  return zmx;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GALAXY::get_zbins_same_ngal(int nbins, int iiz, real_prec  zmn, real_prec  zmx, vector<vector<real_prec> >&zzv){

  // ************************************
  // Construct the dNdz with many bins
  real_prec  ZX=-100;
  int na=20000; //This value is critical to avoid seg foults. The higher, the best.
  real_prec  delta=(zmx-zmn)/((real_prec)na);
  vector<int>dn(na,0);
  for(int i=0;i<this->NGAL;++i){
    if(this->prop[iiz+i*this->NCOLS]<zmx && this->prop[iiz+i*this->NCOLS]>=zmn){
      int iza=floor((this->prop[iiz+i*this->NCOLS]-zmn)/delta);
      dn[iza]++;
    }
  }
  // ************************************
  
  int Nca=(int)(floor)(prop.size()/((real_prec)nbins)); //Desired number of galaxies per redshift bin:
  cout<<"Number of galaxies per redshift bin = "<<Nca<<endl;
  vector<real_prec>zan(na,0);
  for(int i=0;i<dn.size();++i)zan[i]=zmn+(i+0.5)*delta;

  // Set the full z-interval
  zzv[0][0]=zmn;
  zzv[0][1]=zmx;
  if(nbins==1){
    zzv[1][0]=zmn;
    zzv[1][1]=zmx;
  }

  if(nbins>1){
    for(int ib=1;ib<=nbins;++ib){
      int caa=0;
      vector<real_prec>zaux;
      for(int i=0;i<dn.size();++i){
        caa+=dn[i];  //Cumulative number of galaxies
        if((caa>=Nca*(ib-1)) &&  (caa<Nca*ib))
            zaux.push_back(zan[i]);
      }
      zzv[ib][0]=zaux[0]-0.5*delta;  //Allocate the zmin of the ib zbin
      zzv[ib][1]=zaux[zaux.size()-1]+0.5*delta; //Allocate the zmax of the ib zbin
      zaux.clear();
    }
  }
  dn.clear();
  return ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_dndz(){
  this->dNdz.clear();
  this->dNdz.resize(this->N_bin_z,0);
  this->dNdz.shrink_to_fit();
  this->edNdz.clear();
  this->edNdz.shrink_to_fit();
  this->edNdz.resize(this->N_bin_z,0);
  
  this->dNdz_low_res.clear();
  this->dNdz_low_res.shrink_to_fit();

  this->dNdz_low_res.resize(this->N_bin_z_low_res,0);
  this->edNdz_low_res.clear();
  this->edNdz_low_res.shrink_to_fit();
  this->edNdz_low_res.resize(this->N_bin_z_low_res,0);
  
  // ***************************************************************
  // select objects in pixels observed following the mask. 
  // The mask cuts includes the cts in EBV and log_stellar_den, 
  // so no need to make them here again.
  // ***************************************************************
#ifdef _USE_OMP_
    omp_set_num_threads(_NTHREADS_);
#endif

  time_t start;
  time(&start);
  
  So.message_screen("Measuring P(z) (Number of objects in z bins)");

#ifdef _USE_COLOR_
  this->Gal_inside_intervals_z_m_c.resize(this->NGAL,0);
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<this->NGAL;i++)
    {
#ifdef _USE_MASK_
      if(this->gal_mask[i]==this->observed_pixels_in_mask)
	{
#endif
#ifdef _USE_MK_
              if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	    {
#endif
              if(this->prop[i_zs+i*this->NCOLS]>=zzns)
		{
#ifdef _USE_COLOR_
                   real_prec  JKcolor = this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS];
                   if(JKcolor>=this->color_min && JKcolor<this->color_max)
                     {
#endif
                       real_prec  zgal=this->prop[i_z+i*this->NCOLS];
                       if(zgal>=this->z_min && zgal<this->z_max)
                         {
#ifdef _USE_COLOR_
                          this->Gal_inside_intervals_z_m_c[i]=num_1;
#endif

                               int iz=get_bin(this->prop[i_z+i*this->NCOLS],this->z_min,this->N_bin_z,deltaz,true);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif

                          dNdz[iz]++;
                         }

                      if((this->prop[i_z+i*this->NCOLS]>=this->z_min_low_res && prop[this->i_z+i*this->NCOLS]<this->z_max_low_res))
                         {
                          int iz_low_res=get_bin(this->prop[i_z+i*this->NCOLS],this->z_min,this->N_bin_z,deltaz_low_res,true);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                          dNdz_low_res[iz_low_res]++;
                         }
                     }
#ifdef _USE_COLOR_
               }
#endif
#ifdef _USE_MK_
          }
#endif
#ifdef _USE_MASK_
              }
#endif
      }
 So.DONE();
  this->count=0;
  for(auto i=0;i<dNdz.size();++i)
    this->count+=dNdz[i];

  this->Ngal_dndz=count;
  
  // Here we divide by the width of the redshift intercal such that
  // When we plot the dNdz, hte heigh will be independent of the chosen Deltaz
  // Hence, at those plots, we will be seeing
  // dN/dz, and if we want to know the true number of objects per bin, then
  // we would need to multiply by the value of Dz  

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<dNdz.size();i++)
    this->dNdz[i]/=this->deltaz;


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<dNdz_low_res.size();i++)
      this->dNdz_low_res[i]/=(this->deltaz_low_res*static_cast<real_prec>(count));


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<dNdz.size();i++){
    auto Nk=(real_prec)this->dNdz[i];
    this->edNdz[i]= Nk==0? 0 : sqrt(Nk);
  }

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<dNdz_low_res.size();i++){
    auto Nk=(real_prec)this->dNdz_low_res[i];
    this->edNdz_low_res[i]= Nk==0? 0 : sqrt(Nk);
  }

  //this->File.write_to_file(this->output_dir+"dndz.txt",zn_low_res, dNdz_low_res, edNdz_low_res);
  So.message_screen("Number of objects in requested (redshift, mK,color) interval =",this->count);
  cout<<endl;
 // bring dNdz low res to dN/dz
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<dNdz_low_res.size();i++)
      this->dNdz_low_res[i]*=static_cast<real_prec>(count);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_dndmk(){
    this->dNdmk.resize(this->N_bin_mag ,0);
  real_prec  delta_mk=(this->mK_max-this->mK_min)/((real_prec)this->N_bin_mag);
  v_mk.clear();
  for(int i=0;i<this->N_bin_mag;++i)v_mk.push_back(this->mK_min+((real_prec)i)*delta_mk);
  // ***************************************************************
  // select objects in pixels observed following the mask.
  // The mask cuts includes the cts in EBV and log_stellar_den,
  // so no need to make them here again.
  // ***************************************************************
  time_t start;
  time(&start);
  // omp_set_num_threads(omp_get_max_threads());
  cout<<this->mK_max<<endl;
  std::cout<<CYAN<<"Measuring dN/dK";
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<this->NGAL;i++)
    {
      if(this->gal_mask[i]==this->observed_pixels_in_mask)
	{
	  if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	    {
	      if(this->prop[i_zs+i*this->NCOLS]>=zzns)
		{
		  if((this->prop[i_z+i*this->NCOLS]>=this->z_min && prop[this->i_z+i*this->NCOLS]<=this->z_max))
		    {
                      int im=get_bin(this->prop[i_mK+i*this->NCOLS],this->mK_min,this->N_bin_mag,delta_mk,true);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
                      this->dNdmk[im]++;
		    }
		}
	    }
	}
    }
  
  int ngal_mk=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ngal_mk)
#endif
  for(int i=0;i<dNdmk.size();++i)
    {
      real_prec  Nk=static_cast<real_prec>(this->dNdmk[i]);
      this->edNdz[i]= Nk==0? 0 : sqrt(Nk);
      ngal_mk+=this->dNdmk[i];
    }
  
  // Normalize the counts in app magnitude

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<dNdmk.size();++i)
    this->dNdmk[i]/=static_cast<real_prec>(ngal_mk);
  
  
  //this->File.write_to_file(this->output_file_dndmk,v_mk,dNdmk, edNdz);
  std::cout<<BLUE<<count<<" objects in N(z), in file"<<this->output_file_dndz<<RESET<<endl ;
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_z_pdf(){
  time_t start;
  time(&start);
  vector<vector<real_prec  > > z_ps_grid(this->N_bin_z,vector<real_prec>(this->N_bin_z,0));
  vector<vector<real_prec  > > z_ps_grid2(this->N_bin_z,vector<real_prec>(this->N_bin_z,0));

  auto Nbins_z_new=16; //We'll show 16 plots in SM
  vector<vector<real_prec  > > z_ps;
  z_ps.resize(Nbins_z_new);
  auto deltaz_new=  (this->z_max-this->z_min)/((real_prec)Nbins_z_new);
  vector<real_prec> ZP(this->N_bin_z,0);
  vector<real_prec> ZS(this->N_bin_z,0);
  omp_set_num_threads(omp_get_max_threads());
  int count=0;
#pragma omp parallel for reduction(+:count)
  for(int i=0;i<this->NGAL;++i)
    {
      if(this->gal_mask[i]==this->observed_pixels_in_mask)
	{
	  if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	    {
	      if(this->prop[this->i_zs+i*this->NCOLS]>=0.)
		{
		  if((this->prop[this->i_zp+i*this->NCOLS]>=this->z_min && prop[this->i_zp+i*this->NCOLS]<this->z_max) && (this->prop[this->i_zs+i*this->NCOLS]>=this->z_min && prop[this->i_zs+i*this->NCOLS]<this->z_max))
		    {
		      
		      int izp=floor((this->prop[this->i_zp+i*this->NCOLS]-this->z_min)/this->deltaz);
		      int izs=floor((this->prop[this->i_zs+i*this->NCOLS]-this->z_min)/this->deltaz);
		      z_ps_grid[izs][izp]++;
#pragma omp atomic update
		      ZP[izp]++;
#pragma omp atomic update
		      ZS[izs]++;
		      int izpnew=floor((this->prop[this->i_zp+i*this->NCOLS]-this->z_min)/deltaz_new);
		      //z_ps[izpnew].push_back(this->prop[this->i_zs+i*this->NCOLS]);
		      count++;
		    }
		}
	    }
	}
    }


  vector<real_prec> z_photi(count,0);
  vector<real_prec> z_spec(count,0);
  
  count=0;
  for(int i=0;i<this->NGAL;++i)
    {
      if(this->gal_mask[i]==this->observed_pixels_in_mask)
	{
	  if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	    {
	      if(this->prop[this->i_zs+i*this->NCOLS]>=0.)
		{
		  if((this->prop[this->i_zp+i*this->NCOLS]>=this->z_min && prop[this->i_zp+i*this->NCOLS]<this->z_max) && (this->prop[this->i_zs+i*this->NCOLS]>=this->z_min && prop[this->i_zs+i*this->NCOLS]<this->z_max))
		    {
		      z_photi[count]=this->prop[this->i_zp+i*this->NCOLS];
		      z_spec[count]=(this->prop[this->i_zs+i*this->NCOLS]);
		      count++;
		    }
		}
	    }
	}
    }
#pragma omp parallel for
  for(int i=0;i< z_photi.size();++i)
    z_photi[i]-=z_spec[i];
  real_prec  mean=0;
#pragma omp parallel for reduction(+:mean)
  for(int i=0;i< z_photi.size();++i)
    mean+=z_photi[i];
  mean/=static_cast<real_prec>(z_photi.size());
    real_prec  var=0;
#pragma omp parallel for reduction(+:var)
  for(int i=0;i< z_photi.size();++i)
    var+=pow(z_photi[i]-mean,2);
  var/=static_cast<real_prec>(z_photi.size()-1);
  cout<<"Mean "<<mean<<"   Var = "<<sqrt(var)<<endl;
  
  cout<<CYAN<<count<<" counted objects with two types of redshifts"<<RESET<<endl;

  // Her we do the zp_mk histograms, taking all zp available
  for(int i=0;i<this->N_bin_z;++i)for(int j=0;j<this->N_bin_z;++j)z_ps_grid2[i][j]= (ZS[i]==0? 0 : ((real_prec)z_ps_grid[i][j])/((real_prec)ZS[i]));

  string filezs="2d_grid_joint_zp_zs_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+".txt";
  string filezs2="2d_grid_zp_zs_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+".txt";

  ofstream sal, sal2;
  sal.open(filezs);
  sal2.open(filezs2);
  real_prec  maxx=-1000.,maxi;
  for(int i=0;i<this->N_bin_z;++i)
    {
      for(int j=0;j<this->N_bin_z;++j)
	{
	  if(z_ps_grid2[i][j]!=0)
	    {
	      maxi=max(maxx,z_ps_grid[i][j]);
	      maxx=maxi;
	    }
	}
    }
  for(int i=0;i<this->N_bin_z;++i)
    for(int j=0;j<this->N_bin_z;++j)
      sal<<zn[i]<<"\t"<<zn[j]<<"\t"<<log10(z_ps_grid[i][j])<<endl;

  for(int i=0;i<this->N_bin_z;++i)
    for(int j=0;j<this->N_bin_z;++j)
      sal2<<zn[i]<<"\t"<<zn[j]<<"\t"<<log10(z_ps_grid2[i][j]/maxi)<<endl;
  
  sal.close();
  sal2.close();

  vector<real_prec> meanz(Nbins_z_new,0);
  for(int i=0;i<Nbins_z_new;++i)
      for(int j=0;j<z_ps[i].size();++j)
          meanz[i]+=z_ps[i][j]/(z_ps[i].size());

  vector<real_prec> varz(Nbins_z_new,0);
  for(int i=0;i<Nbins_z_new;++i)
      for(int j=0;j<z_ps[i].size();++j)
          varz[i]+=pow(meanz[i]-z_ps[i][j],2)/(z_ps[i].size()-1);

  filezs="zs_zp_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+".txt";
  sal.open(filezs);
  for(int i=0;i<Nbins_z_new;++i)
      sal<<this->z_min+(i+0.5)*deltaz_new<<"\t"<<0.5*deltaz_new<<"\t"<<"\t"<<meanz[i]<<"\t"<<sqrt(varz[i])<<"  "<<varz[i]*sqrt((20.0/(z_ps[i].size()-1)))<<endl;
  sal.close();


if(this->redshift_type=="p")
    {
  for(int i=0;i<Nbins_z_new;++i)
    {
      filezs="zs_zpbin"+to_string(i)+".txt";
      sal.open(filezs.c_str());
      for(int j=0;j<z_ps[i].size();++j)
	sal<<z_ps[i][j]-meanz[i]<<endl;
      sal.close();
    }
    }
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_M_pdf(){
  time_t start;
  time(&start);

  int Nbins_new=16; //We'll show 16 plots in SM
  vector<vector<real_prec  > > z_ps;
  z_ps.resize(Nbins_new);
  real_prec  delta_new=  (this->MK_max-this->MK_min)/((real_prec)Nbins_new);
  count=0;
  for(int i=0;i<this->NGAL;++i){
    if(this->gal_mask[i]==this->observed_pixels_in_mask)
      {
	if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	  {
	    if(this->prop[this->i_zs+i*this->NCOLS]>=0.)
	      {
		if((this->prop[this->i_zp+i*this->NCOLS]>=this->z_min && prop[this->i_zp+i*this->NCOLS]<this->z_max) && (this->prop[this->i_zs+i*this->NCOLS]>=this->z_min && prop[this->i_zs+i*this->NCOLS]<this->z_max))
		  {
            real_prec  MKs;
		    this->i_z=this->i_zs;
		    get_MK(prop[this->i_zs+i*this->NCOLS], prop[this->i_mK+i*this->NCOLS], MKs); //Compute ABsolute magnitude spectroscopic
            real_prec  MKp;
		    this->i_z=this->i_zp;
		    get_MK(prop[this->i_z+i*this->NCOLS], prop[this->i_mK+i*this->NCOLS], MKp); //Compute ABsolute magnitude photometric
		    //            cout<<MKp<<"  "<<MKs<<endl;
		    if((MKs>=this->MK_min && MKs<this->MK_max) && (MKp>=this->MK_min && MKp<this->MK_max)){ 	  // Check boundaries in abs Mag
		      int iMpnew=floor((MKp-this->MK_min)/delta_new);
		      z_ps[iMpnew].push_back(MKs);
		      count++;
		    }
		  }
	      }
	  }
      }
  }
  cout<<CYAN<<count<<" counted objects with two types of redshifts"<<RESET<<endl;
  vector<real_prec> meanz(Nbins_new,0);
  for(int i=0;i<Nbins_new;++i)
    for(int j=0;j<z_ps[i].size();++j)
      meanz[i]+=z_ps[i][j]/(z_ps[i].size());

  vector<real_prec> varz(Nbins_new,0);
  for(int i=0;i<Nbins_new;++i)
    for(int j=0;j<z_ps[i].size();++j)
      varz[i]+=pow(meanz[i]-z_ps[i][j],2)/(z_ps[i].size()-1);
  
  ofstream sal;
  string filezs="Magnitudes/"+this->name_cat+"_mean_Ms_Mpbins_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+".txt";
  sal.open(filezs);
  for(int i=0;i<Nbins_new;++i)
    sal<<this->MK_min+(i+0.5)*delta_new<<"\t"<<0.5*delta_new<<"\t"<<"\t"<<meanz[i]<<"\t"<<sqrt(varz[i])<<"  "<<varz[i]*sqrt((20.0/(z_ps[i].size()-1)))<<endl;
  sal.close();
  
  
  // Hwew we write for each Mp bin the difference between Ms and <Ms|Mp> to be shown in supermongo, plots.sm, dms
  for(int i=0;i<Nbins_new;++i)
    {
      filezs="Magnitudes/"+this->name_cat+"_Ms_Mpbin"+to_string(i)+".txt";
      sal.open(filezs);
      for(int j=0;j<z_ps[i].size();++j)
	sal<<z_ps[i][j]-meanz[i]<<endl;
      sal.close();
    }  
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_ipix_from_mask(int i, Healpix_Map<real_prec>map_aux, long & ipix){
  this->point.phi=fac*this->prop[this->i_lgal+i*this->NCOLS];
  this->point.theta=0.5*M_PI-fac*this->prop[this->i_bgal+i*this->NCOLS];
  ipix=map_aux.ang2pix(this->point);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void GALAXY::get_Vmax(){
  // Here we obtain the Vmax and the Sij matrix , used for the LF as well as
  // for the joint ptobability of color-MK

  omp_set_num_threads(omp_get_max_threads());

  this->Vmax_v.resize(this->N_bin_Mag,0);
  // -------------------------------------------------------
  // Here we redefine the limits in Absolute magnitude
  // in order to include all the objects that were included when doing the dNdz
  // when no MK were already assigned.
  real_prec  MK_MAX_NEW=-1000.0;
  real_prec  MK_MIN_NEW=1000.0;
  real_prec  MK_MAX_NEW_AUX, MK_MIN_NEW_AUX;
    
  int i;
    
  if(this->Redefine_MKminmax)
    {
      for(i=0;i<this->NGAL;++i)
	{
	  if(this->gal_mask[i]==this->observed_pixels_in_mask)
	    {
	      if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
		{
		  if(this->prop[i_zs+i*this->NCOLS]>=zzns)
		    {
		      if((this->prop[i_z+i*this->NCOLS]>=this->z_min && prop[this->i_z+i*this->NCOLS]<this->z_max))
			{

              real_prec  MKs;
              real_prec  color1=this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mK+i*this->NCOLS];
              real_prec  color2=this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mH+i*this->NCOLS];

			  this->get_MK(this->prop[i_z+i*this->NCOLS],this->prop[i_mK+i*this->NCOLS], color1, color2,MKs); //Compute ABsolute magnitude
			  MK_MAX_NEW_AUX=max(MK_MAX_NEW,MKs);
			  MK_MIN_NEW_AUX=min(MK_MIN_NEW,MKs);
			  MK_MIN_NEW=MK_MIN_NEW_AUX;
			  MK_MAX_NEW=MK_MAX_NEW_AUX;
			}
		    }
		}
	    }
	}
      this->MK_max=MK_MAX_NEW+EPSILON_MK;
      this->MK_min=MK_MIN_NEW-EPSILON_MK;
    }
        
  real_prec  deltaMK_h=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag;
  real_prec  deltaMK_l=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag_low_res;
  this->v_Magnitude.resize(this->N_bin_Mag,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_Mag;++i)
    this->v_Magnitude[i]=this->MK_min+(i+0.5)*deltaMK_h;

  this->v_Magnitude_low_res.resize(this->N_bin_Mag_low_res,0);
  for(int i=0;i<this->N_bin_Mag_low_res;++i)
    this->v_Magnitude_low_res[i]=this->MK_min+(i+0.5)*deltaMK_l;

  if(this->Redefine_MKminmax)
    {
      cout<<GREEN<<"New MK min = "<<this->MK_max<<endl;
      cout<<GREEN<<"New MK max = "<<this->MK_min<<RESET<<endl;
      cout<<CYAN<<"*******************************************************"<<RESET<<endl;
    }
  // -------------------------------------------------------
  //    omp_set_num_threads(omp_get_max_threads());
      
  if(this->LF_estimator=="Vmax_dc" || this->LF_estimator=="Vmax")
    {
      So.message_screen("Computing Sij ...");
      this->S.resize(this->N_bin_z);
      for(int i=0;i<this->S.size();i++)
	this->S[i].resize(this->N_bin_Mag,0);
	
      bool cw=false;  //true if the color is used in the K-correction. False if not


      if(true==cw)
	get_Pzc(cw); //get the total color in z bins

#pragma omp parallel for
      for(int i=0;i<this->N_bin_z;++i)
	{
	  // Get the limiting absolute magnitude given the mK_max and the minimum redshift of the z-bin
      real_prec  mean_color_zbin=0;
	  // Get the mean color in the redhsift bin
	  if(true==cw)
            {
	      for(int j=0;j<this->N_bin_color;++j)
        mean_color_zbin+=this->Pzc[i][j]/static_cast<real_prec>(this->dNdz[i]);
            }
	  // Get MABS for m_Kmax in the lower limit of the z-bin
      real_prec  Dma1=gsl_inter_new(this->zv, this->Dm, this->zn_min[i]);
      real_prec  K_e1 =cosmology.K_correction(this->zn_min[i],mean_color_zbin, 0.0)+cosmology.e_correction(this->zn_min[i]);
      real_prec  M_lim_zmin1=this->mK_max - Dma1 - K_e1;
	  // Get MABS for m_Kmin in the low limit of the z-bin
      real_prec  M_lim_zmax1=this->mK_min - Dma1 - K_e1;
          // Get MABS for m_Kmax in the upper limit of the z-bin
      real_prec  Dma2=gsl_inter_new(this->zv, this->Dm, this->zn_max[i]);
      real_prec  K_e2 =cosmology.K_correction(this->zn_max[i],mean_color_zbin,0)+cosmology.e_correction(this->zn_max[i]);
      real_prec  M_lim_zmin2=this->mK_max - Dma2 - K_e2;
	  // Get MABS for m_Kmin in the low limit of the z-bin
      real_prec  M_lim_zmax2=this->mK_min - Dma2 - K_e2;
      real_prec  M_lim_zmin=M_lim_zmin1+0.5*(M_lim_zmin2-M_lim_zmin1);
      real_prec  M_lim_zmax=M_lim_zmax1+0.5*(M_lim_zmax2-M_lim_zmax1);
      real_prec  M_zmin=min(M_lim_zmin,this->MK_max);
      real_prec  M_zmax=min(M_lim_zmax,this->MK_min);
      for(int j=0;j<this->N_bin_Mag;++j)
	    {
	      this->S[i][j]=1.0;
	      if(this->v_Magnitude[j]> M_zmin  || this->v_Magnitude[j]< M_zmax )
		this->S[i][j]=0.0;
	    }
	}	
      So.DONE();
      So.message_screen("Computing Vmax(M) ...");
#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	this->Vmax_v[j]=0;

#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	for(int i=0;i<this->N_bin_z;++i)
	  this->Vmax_v[j]+=this->S[i][j]*this->Vshell[i];//*pow(1+zn[i],this->scp.e_index_zstar);
      So.DONE();
	    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_PzK(){
    this->Pzm.resize(this->N_bin_z_low_res);
  for(int i=0;i<Pzm.size();++i)
      Pzm[i].resize(this->N_bin_mag, 0);

  real_prec  deltaz=(this->z_max-this->z_min)/(real_prec)this->N_bin_z_low_res;
  real_prec  deltamK=(this->mK_max-this->mK_min)/(real_prec)this->N_bin_mag;
  zn.clear();
  for(int i=0;i<this->N_bin_z_low_res;++i)zn.push_back(this->z_min+(i+0.5)*deltaz);
  v_mk.clear();
  for(int i=0;i<this->N_bin_mag;++i)v_mk.push_back(this->mK_min+(i+0.5)*deltamK);
  So.message_screen("Measuring  N(z,m)");
  for(int i=0;i<this->NGAL;++i){
    if(this->gal_mask[i]==this->observed_pixels_in_mask){
      if(this->prop[i_z+i*this->NCOLS]<this->z_max && this->prop[i_z+i*this->NCOLS]>=this->z_min){
        if(this->prop[this->i_mK+i*this->NCOLS]<this->mK_max && this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min)
	  {
	    int iz=get_bin(this->prop[this->i_z+i*this->NCOLS],this->z_min,this->N_bin_z,deltaz,true);
	    int im=get_bin(this->prop[this->i_mK+i*this->NCOLS],this->mK_min,this->N_bin_mag,deltamK,true);
	    this->Pzm[iz][im]++;
	  }
      }
    }
  }

  So.DONE();

/*
  string of="2d_grid_z_K.txt";
  ofstream te;
  te.open(of);
  for(int i=0;i<this->N_bin_z_low_res;++i)
     for(int j=0;j<this->N_bin_mag;j++)
       te<<this->zn[i]<<"\t"<<this->v_mk[j]<<"\t"<<log10(this->Pzm[i][j])<<endl;
  te.close();

*/
real_prec  acount =0;
#pragma omp parallel for reduction (+:acount) collapse (2)
  for(int i=0;i<this->N_bin_z_low_res;++i)
    for(int j=0;j<this->N_bin_mag;j++)
      acount+=this->Pzm[i][j];
  So.message_screen("Number of ojects from P(z,M) = ", acount);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_PcMk(){
  this->PcMk.clear();
  this->PcMk.resize(this->N_bin_color);
  for(int i=0;i<PcMk[0].size();++i)
      PcMk[i].resize(this->N_bin_Mag_low_res, 0);
  So.message_screen("Measuring  N(c,M)");
  real_prec  deltac=(this->color_max-this->color_min)/(real_prec)this->N_bin_color;
  real_prec  deltaMK_low_res=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag_low_res;
  real_prec  deltaMK=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag;
  this->v_Magnitude_low_res.clear();
  for(int i=0;i<this->N_bin_Mag_low_res;++i)
      this->v_Magnitude_low_res.push_back(this->MK_min+(i+0.5)*deltaMK_low_res);
  v_color.clear();
  for(int i=0;i<this->N_bin_color;++i)
      v_color.push_back(this->color_min+(i+0.5)*deltac);
  for(int i=0;i<this->NGAL;++i){
    if(this->gal_mask[i]==this->observed_pixels_in_mask){
      if(this->prop[i_zs+i*this->NCOLS]>=this->zzns){
        if(this->prop[i_z+i*this->NCOLS]<this->z_max && this->prop[i_z+i*this->NCOLS]>=this->z_min){
          if(this->prop[this->i_mK+i*this->NCOLS]<this->mK_max && this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min){
            real_prec  JKcolor=(this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS])-gsl_inter_new(this->v_z_Kcorr,this->v_color_correction_Kcorr,prop[i_zp+i*this->NCOLS]);
	    if(JKcolor>=this->color_min && JKcolor<this->color_max){
              real_prec  MKs;
              real_prec  color1=this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mK+i*this->NCOLS];
              this->get_MK(this->prop[i_z+i*this->NCOLS], this->prop[i_mK+i*this->NCOLS],color1,0.0, MKs); //Compute ABsolute magnitude

	      if(MKs>=this->MK_min && MKs<this->MK_max){ 	  // Check boundaries in abs Mag
		int ic=(int)floor((JKcolor- this->color_min)/deltac);
		int iM=(int)floor((MKs-this->MK_min)/deltaMK_low_res);
                this->PcMk[ic][iM]++;
              }
	    }
	  }
        }
      }
    }
  }
  string of="2d_grid_c_MK.txt";
  ofstream te;
  te.open(of);
  for(int i=0;i<this->N_bin_Mag_low_res;++i)for(int j=0;j<this->N_bin_color;j++)te<<this->v_Magnitude_low_res[i]<<"\t"<<this->v_color[j]<<"\t"<<log10(this->PcMk[i][j])<<endl;
  te.close();
  count =0;
  for(int j=0;j<this->N_bin_color;j++)
     for(int i=0;i<this->N_bin_Mag_low_res;++i)
  count+=this->PcMk[j][i];
  cout<<BLUE<<count<<" objects in PcM, in file "<<of<<RESET<<endl;
  So.DONE();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_PcMk_Vmax(){

  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);

  So.message_screen("Measuring  P(c,M) using 1/Vmax");

  real_prec  deltac=(this->color_max-this->color_min)/static_cast<real_prec>(this->N_bin_color);
  real_prec  deltaMK=(this->MK_max-this->MK_min)/static_cast<real_prec>(this->N_bin_Mag);

  this->PcMk.clear();
  this->PcMk.resize(this->N_bin_color);
  for(int i=0;i<this->N_bin_color;++i)
      PcMk[i].resize(this->N_bin_Mag, 0);


#pragma omp parallel for
  for(int i=0;i<this->NGAL;++i)
    {
    if(this->Gal_inside_intervals_z_m_c[i]>0)
      {
         real_prec  MKs;
         real_prec  color1=this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mK+i*this->NCOLS];
         this->get_MK(this->prop[i_z+i*this->NCOLS], this->prop[i_mK+i*this->NCOLS],color1,0.0, MKs); //Compute ABsolute magnitude
         if(MKs>=this->MK_min && MKs<this->MK_max)
          { 	  // Check boundaries in abs Mag
            int ic=get_bin(color1,this->color_min,this->N_bin_color,deltac,true);
            int iM=get_bin(MKs,this->MK_min,this->N_bin_Mag,deltaMK,true);
            real_prec  vmax=this->Vmax_v[iM];
#pragma omp atomic update
            this->PcMk[ic][iM]+=1./static_cast<real_prec>(vmax*deltaMK);
          }
       }
    }
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_PzMK(){
  So.message_screen("Measuring  P(z,M)");
  this->PzM.clear();
  this->PzM.resize(this->N_bin_z_low_res);
  for(int i=0;i<PzM.size();++i)PzM[i].resize(this->N_bin_Mag_low_res, 0);
  real_prec  deltaz=(this->z_max-this->z_min)/(real_prec)this->N_bin_z_low_res;
  real_prec  deltaMK=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag_low_res;
  v_Magnitude.clear();
  v_Magnitude.resize(this->N_bin_Mag,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_Mag;++i)
    v_Magnitude[i]=this->MK_min+(static_cast<real_prec>(i))*deltaMK;
  zn.clear();
  zn.resize(this->N_bin_z_low_res,0);
#pragma omp parallel for
  for(int i=0;i<this->N_bin_z_low_res;++i)
    zn[i]=(this->z_min+(i+0.5)*deltaz);
  for(int i=0;i<this->NGAL;++i)
    {
     if(this->Gal_inside_intervals_z_m_c[i]>0)
       {
         real_prec  MKs;
	  get_MK(this->prop[i_z+i*this->NCOLS],this->prop[i_mK+i*this->NCOLS], MKs); //Compute ABsolute magnitude
          if(MKs>=this->MK_min && MKs<this->MK_max)
            { 	  // Check boundaries in abs Mag
	    int iz=(int)floor((this->prop[i_z+i*this->NCOLS]- this->z_min)/deltaz);
	    int im=(int)floor((MKs-this->MK_min)/deltaMK);
	    this->PzM[iz][im]++;
          }
        }
  }
  So.DONE();
  string of="2d_grid_z_MK.txt";
  ofstream te;
  te.open(of);
  for(int i=0;i<this->N_bin_z_low_res;++i)for(int j=0;j<this->N_bin_Mag_low_res;j++)
      te<<this->zn[i]<<"\t"<<this->v_Magnitude[j]<<"\t"<<log10(this->PzM[i][j])<<endl;
  te.close();
  real_prec  acount =0;
#pragma omp parallel for reduction (+:acount) collapse (2)
  for(int i=0;i<this->N_bin_z;++i)
      for(int j=0;j<this->N_bin_Mag;j++)
          acount+=this->PzM[i][j];
  So.message_screen("Number of ojects from P(z,M) = ", acount);
  So.DONE();
  zn.clear(); zn.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_P_X_Y(s_bins_info *bx, s_bins_info *by, string type){
  So.message_screen("Measuring  P(z,M)");
  this->PXY.clear(); // let us recyclye this
  this->PXY.resize(bx->Nbins*by->Nbins,0);
  this->NXY.clear(); // let us recyclye this
  this->NXY.resize(bx->Nbins*by->Nbins,0);
  this->nbarXY.clear();
  this->nbarXY.resize(bx->Nbins*by->Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<this->NGAL;++i)
    {
      int lix,liy;
      if(bx->type=="log")
        lix=get_bin(log10(this->prop[bx->index+i*this->NCOLS]),bx->min,bx->Nbins, bx->delta,true);
      else
        lix=get_bin(this->prop[bx->index+i*this->NCOLS],bx->min,bx->Nbins, bx->delta,true);
      if(by->type=="log")
        liy=get_bin(log10(this->prop[by->index+i*this->NCOLS]),by->min,by->Nbins, by->delta,true);
      else
        liy=get_bin(this->prop[by->index+i*this->NCOLS],by->min,by->Nbins, by->delta,true);

#ifdef _USE_OMP_
#pragma omp atomic
      this->NXY[index_2d(liy,lix,bx->Nbins)]++;
#endif
      }

  ofstream te;
  te.open(this->output_dir+"PzMass_"+type+".txt");
  for(int i=0;i<bx->Nbins;++i)
    for(int j=0;j<by->Nbins;++j)
       te<< bx->min+(i+0.5)*bx->delta <<"\t"<<by->min+(j+0.5)*by->delta<<"\t"<<NXY[index_2d(j,i,bx->Nbins)]<<endl;
  te.close();
 // Divide by the volume of the z shell to get nbar
  for(int i=0;i<bx->Nbins;++i)
    for(int j=0;j<by->Nbins;++j)
        this->nbarXY[index_2d(j,i,bx->Nbins)]/=this->Vshell_lowres[i];


  // normalization for every bin in z
  for(int i=0;i<bx->Nbins;++i)
    {
    double  aux_a;
    double  aux_b=-1000.0;
    for(int j=0;j<by->Nbins;++j)
      {
        ULONG ind=index_2d(j,i,bx->Nbins);
        aux_a=max(aux_b,nbar[ind]);
        aux_b=aux_a;
     }

    for(int j=0;j<by->Nbins;++j)
        {
        ULONG ind=index_2d(j,i,bx->Nbins);
        PXY[ind]= aux_a==0 ? 0 : nbar[ind]/static_cast<real_prec>(aux_a);
        }
      }
  //PXY must be smoothed
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_Pzc_from_PcMk(){
    int NTHREADS = omp_get_max_threads();
    omp_set_num_threads(NTHREADS);
    So.message_screen("Computing P(z,c) from P(c,Mk(z))");
    real_prec  deltaMK=(this->MK_max-this->MK_min)/static_cast<real_prec>(this->N_bin_Mag);
    this->v_Magnitude.resize(this->N_bin_Mag,0);
#pragma omp parallel for
    for(int i=0;i<this->N_bin_Mag;++i)
      this->v_Magnitude[i]=this->MK_min+(i+0.5)*deltaMK;
    this->Pzc_from_PzMk.clear();
    this->Pzc_from_PzMk.resize(this->N_bin_z);
#pragma omp parallel for
    for(int i=0;i<Pzc_from_PzMk.size();++i)
       Pzc_from_PzMk[i].resize(this->N_bin_Mag, 0);
#pragma omp parallel for collapse(3)
    for(int i=0 ; i< this->N_bin_z; ++i)
       for(int j=0 ; j < this->N_bin_color; ++j)
           for(int k=0; k < this->N_bin_Mag; ++k)
#pragma omp atomic update
              this->Pzc_from_PzMk[i][j]+=this->PcMk[j][k]*this->S[i][k]*this->Vshell[i]*deltaMK;
    So.DONE();

    string of="2d_grid_z_c_mcmc.txt";
    ofstream te;
    te.open(of);
    for(int i=0;i<this->N_bin_z;++i)
      for(int j=0;j<this->N_bin_color;j++)
        te<<this->zn[i]<<"\t"<<this->v_color[j]<<"\t"<<(this->Pzc_from_PzMk[i][j])<<endl;
    te.close();
    So.message_screen("Computing Phi(Mk) from P(c,Mk)");
    this->Pm.clear();
    this->Pm.resize(this->N_bin_Mag,0);
#pragma omp parallel for collapse(2)
     for(int k=0; k < this->N_bin_Mag; ++k)
       for(int j=0 ; j < this->N_bin_color; ++j)
         this->Pm[k]+=this->PcMk[j][k];
   So.DONE();
   //this->File.write_to_file(this->output_file_LF,this->v_Magnitude,this->Pm);
  So.message_screen("Computing N(z) from P(z,c)");
  this->dNdz_lf.clear();
  this->dNdz_lf.resize(this->N_bin_z,0);
#pragma omp parallel for
    for(int i=0;i<this->N_bin_z;++i)
      for(int j=0; j<this->N_bin_Mag ; ++j)
#pragma omp atomic update
        this->dNdz_lf[i]+=this->Pzc_from_PzMk[i][j];
    So.DONE();
    real_prec  acount=0;
#pragma omp parallel for reduction(+:acount)
    for(int i=0;i<this->N_bin_z;++i)
      acount+=this->dNdz_lf[i];
    So.message_screen("Check: Predicted number of objects from P(z,c) = ", static_cast<ULONG>(acount));
    this->Ngal_expected=acount;
      // Divide to obtain the number of objects in the bin
#pragma omp parallel for
  for(int i=0;i<this->dNdz_lf.size();i++)
    this->dNdz_lf[i]/=static_cast<real_prec>(deltaz);
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_Pzc(bool c_weight){
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  So.message_screen("Measuring  P(z,c)");
  real_prec  deltac=(this->color_max-this->color_min)/(real_prec)this->N_bin_color;
  real_prec  deltaz=(this->z_max-this->z_min)/(real_prec)this->N_bin_z;
  this->Pzc.resize(this->N_bin_z);
#pragma omp parallel for
  for(int i=0;i<Pzc.size();++i)
     this->Pzc[i].resize(this->N_bin_color, 0);
  if(false==c_weight)
    {
      v_color.resize(this->N_bin_color,0);
#pragma omp parallel for
      for(int i=0;i<this->N_bin_color;++i)
        v_color[i]=this->color_min+(i+0.5)*deltac;
    }

  if(false==c_weight)
  {
#pragma omp parallel for
  for(int i=0;i<this->NGAL;++i)
    if(this->Gal_inside_intervals_z_m_c[i]>0)
      {
         real_prec  zgal=this->prop[i_z+i*this->NCOLS];
         real_prec  JKcolor=this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS];
         int iz=get_bin(zgal,this->z_min,this->N_bin_z,deltaz,true);
         int ic=get_bin(JKcolor,this->color_min,this->N_bin_color,deltac,true);
#pragma omp atomic update
         this->Pzc[iz][ic]++;  //the histogram
      }
  }
  else
   {
#pragma omp parallel for
  for(int i=0;i<this->NGAL;++i)
    {
      if(this->Gal_inside_intervals_z_m_c[i]>0)
        {
          real_prec  JKcolor=this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS];
          int iz=get_bin(this->prop[i_z+i*this->NCOLS],this->z_min,this->N_bin_z,deltaz,true);
          int ic=get_bin(JKcolor,this->color_min,this->N_bin_color,deltac,true);
#pragma omp atomic update
               this->Pzc[iz][ic]+=JKcolor;  //the histogram
        }
    }
  }
  So.DONE();
  if(false==c_weight)
    {
      string of="2d_grid_z_c.txt";
      ofstream te;
      te.open(of);
      for(int i=0;i<this->N_bin_z;++i)
	for(int j=0;j<this->N_bin_color;j++)
          te<<this->zn[i]<<"\t"<<this->v_color[j]<<"\t"<<(this->Pzc[i][j])<<endl;
      te.close();

      real_prec  acount =0;
#pragma omp parallel for reduction (+:acount) collapse(2)
      for(int i=0;i<this->Pzc.size();++i)
          for(int j=0;j<this->Pzc[0].size();j++)
              acount+=this->Pzc[i][j];

      So.message_screen("Check: Number of objects from P(z,c) = ", acount);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_PKc(){
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  real_prec  deltaK=(this->mK_max-this->mK_min)/(real_prec)this->N_bin_mag;
  real_prec  deltac=(this->color_max-this->color_min)/(real_prec)this->N_bin_color;
  this->PKc.resize(this->N_bin_mag);
  for(int i=0;i<PKc.size();++i)PKc[i].resize(this->N_bin_color, 0);
  v_color.clear();
  for(int i=0;i<this->N_bin_color;++i)v_color.push_back(this->color_min+(i+0.5)*deltac);
  v_mk.clear();
  for(int i=0;i<this->N_bin_mag;++i)v_mk.push_back(this->mK_min+((real_prec)i)*deltaK);
  std::cout<<CYAN<<"Measuring  N(mK,c)"<<RESET<<std::endl;
  for(int i=0;i<this->NGAL;++i){
    if(this->gal_mask[i]==this->observed_pixels_in_mask){
      if(this->prop[i_z+i*this->NCOLS]>=zzns){
        if(this->prop[this->i_mK+i*this->NCOLS]<this->mK_max && this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min){
      real_prec  JKcolor=this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS]-gsl_inter_new(this->v_z_Kcorr,this->v_color_correction_Kcorr,prop[i_zp+i*this->NCOLS]);
	  if(JKcolor>=this->color_min && JKcolor<this->color_max){
	    int ik=(int)floor((this->prop[i_mK+i*this->NCOLS]- this->mK_min)/deltaK);
	    int ic=(int)floor((JKcolor-this->color_min)/deltac);
	    this->PKc[ik][ic]++;
          }
        }
      }
    }
  }
  string of="2d_grid_m_c.txt";
  ofstream te;
  te.open(of);
  for(int i=0;i<this->N_bin_mag;++i)for(int j=0;j<this->N_bin_color;j++)te<<this->v_mk[i]<<"\t"<<this->v_color[j]<<"\t"<<log10(this->PKc[i][j])<<endl;
  te.close();
  count =0;
  for(int i=0;i<this->N_bin_mag;++i)for(int j=0;j<this->N_bin_color;j++)count+=this->PKc[i][j];
  cout<<BLUE<<count<<" objects in Pkc, in file "<<of<<RESET<<endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_NS(){
  std::cout<<CYAN<<"Using "<< omp_get_max_threads()<<" threads"<<RESET<<std::endl;
  time_t start_all;
  time(&start_all);
  time_t end;
  real_prec  delta_mk=(this->mK_max-this->mK_min)/((real_prec)this->N_bin_mag);
  this->Amax_v.resize(this->N_bin_mag,1.0);
  //  this->Amax_v.resize(this->N_bin_mag,this->total_area);
  this->v_mk.clear();
  for(int i=0;i<this->N_bin_mag;++i)v_mk.push_back(this->mK_min+((real_prec)(i+0.5))*delta_mk);
  omp_set_num_threads(omp_get_max_threads());
  // Get here Amax
  this->dNdmk.clear();
  this->dNdmk.resize(this->N_bin_mag,0);
  /*
    for(int k=0;k<this->N_bin_mag;++k){
    for(int i=0;i<this->prop.size();++i){
    if(this->prop[i][this->i_gal_mask]==this->observed_pixels_in_mask){
    if(this->prop[i][i_zs]>=this->zzns){
    real_prec  aux_var=(prop[i][i_mK]<v_mk[k]? 1./this->Amax_v[k]:0);
    this->dNdmk[k]+=aux_var;
    }
    }
    }
    }
  */
  for(int i=0;i<this->NGAL;++i){
    if(this->gal_mask[i]==this->observed_pixels_in_mask){
      if(this->prop[i_zs+i*this->NCOLS]>=this->zzns){
        if(this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max && this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min){
	  int K=(int)floor((this->prop[i_mK+i*this->NCOLS]- this->mK_min)/delta_mk);
	  this->dNdmk[K]+=1/this->Amax_v[K];
        }
      }
    }
  }
  So.DONE();
  for(auto i=0;i<dNdmk.size();i++){
    auto Nk=(real_prec)this->dNdmk[i]*this->total_area;
    this->edNdz[i]= Nk==0? 0 : dNdmk[i]/sqrt(Nk);
  }
  //this->File.write_to_file(this->output_file_dndmk,v_mk,dNdmk, edNdz);
  time(&end);
  real_prec  lapse=difftime(end,start_all);
  std::cout<<BLUE<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
  time(&start_all);
  std::cout<<CYAN<<"*******************************************************"<<RESET<<endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_LF_Vmax(ULONG number_gal){
  const gsl_rng_type * T;
  gsl_rng * r;
   So.message_screen("Measuring LF (Vmax standard estimator)");
  time_t start_all;
  time(&start_all);
  time_t end;
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  real_prec  deltaMK_h=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag;
  real_prec  deltaMK_l=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag_low_res;
  this->number_z.clear();
  this->number_z.resize(this->N_bin_z,0);
  this->Pm.clear();
  this->Pm.resize(this->N_bin_Mag,0);
  this->Pm_low_res.clear();
  this->Pm_low_res.resize(this->N_bin_Mag_low_res,0);
  this->number_m.clear();
  this->number_m.resize(this->N_bin_Mag,0);
  this->number_m_low_res.clear();
  this->number_m_low_res.resize(this->N_bin_Mag_low_res,0);
  vector<int>vseeds(NTHREADS,0);
  for(int i=0;i<vseeds.size();++i)vseeds[i]=12+105*i;   
#ifdef _USE_LF_COLE_
  this->index_z_Gal.resize(this->NGAL,0); //This iwll be use only for COle estimator: allocate mags and read them in the iterations
  this->index_Mk_Gal.resize(this->NGAL,0); //This iwll be use only for COle estimator: allocate mags and read them in the iterations
  this->index_Mk_low_Gal.resize(this->NGAL,0); //This iwll be use only for COle estimator: allocate mags and read them in the iterations
  this->Gal_inside_intervals.resize(this->NGAL,0); //This iwll be use only for COle estimator: allocate mags and read them in the iterations
#endif
#ifdef _USE_PDF_PHOT2SPEC_   //tjis is to speed up the code, ut has to be given in accordance to parameter file
#pragma omp parallel private(T, r)
  {

    gsl_rng_env_setup();
    gsl_rng_default_seed=vseeds[omp_get_thread_num()];
    T = gsl_rng_ranlux;
    r = gsl_rng_alloc (T);
    
#pragma omp for
#else
#pragma omp parallel for
#endif
    for(int i=0;i<this->NGAL;++i)
      {
        if(this->Gal_inside_intervals_z_m_c[i]>0) //This was filled when measuring the dNdz
          {

            real_prec  zgal=this->prop[i_z+i*this->NCOLS];
		    
	    // If p_s or ps_s, assign spectroscopic redshift to this photometric ones
	    // Using the PDF calibrated in paper I
#ifdef _USE_PDF_PHOT2SPEC_   //tjis is to speed up the code, ut has to be given in accordance to parameter file

	    if(this->redshift_type=="ps_s" || this->redshift_type=="p_s")
	      {
		// Width of the Gaussian
        real_prec  zvpdf=0.03*tanh(-20.78*zgal*zgal+7.76*zgal+0.05)/(1.0+zgal);
        real_prec  dist=gsl_ran_gaussian(r,zvpdf);
		zgal+=dist;
	      }
#endif
        real_prec  MKs;
        real_prec  color1=this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mK+i*this->NCOLS];
        real_prec  color2=0.0;//this->prop[i_mJ+i*this->NCOLS]-this->prop[i_mH+i*this->NCOLS];
	    this->get_MK(zgal, this->prop[i_mK+i*this->NCOLS],color1,color2, MKs); //Compute ABsolute magnitude
	    if(MKs>=this->MK_min && MKs<this->MK_max)
	      { 	  // Check boundaries in abs Mag
			    
		int im_h=get_bin(MKs,this->MK_min,this->N_bin_Mag,deltaMK_h,true);
		int im_l=get_bin(MKs,this->MK_min,this->N_bin_Mag_low_res,deltaMK_l,true);
		int iz=get_bin(zgal,this->z_min,this->N_bin_z,deltaz,true);

#ifdef _USE_LF_COLE_
		this->index_Mk_Gal[i]=im_h;
		this->index_Mk_low_Gal[i]=im_l;
		this->index_z_Gal[i]=iz;
#endif
        real_prec  vmax=1;
		if(this->LF_estimator=="Vmax" || this->LF_estimator=="Vmax_dc")
		  vmax=this->Vmax_v[im_h];
		else
		  if(this->LF_estimator=="Vmax_o" )
		    {
		      s_CosmologicalParameters scp_aux;
		      scp_aux=scp;
		      scp_aux.Mabs=MKs;
		      scp_aux.mlim=this->mK_max;
              real_prec  zmax=this->cosmology.zmax();
              real_prec  rmax_v=gsl_inter_new(scp_aux.zv, scp_aux.rv, zmax);
		      vmax=this->total_area*pow(rmax_v,3)/3.;
		    }
#ifdef _USE_LF_COLE_
		this->Gal_inside_intervals[i]=1;
#endif
#pragma omp atomic update
		this->number_z[iz]++;
#pragma omp atomic update
		this->number_m[im_h]++;
#pragma omp atomic update
		this->number_m_low_res[im_l]++;
#pragma omp atomic update
		this->Pm[im_h]+=1./(vmax*deltaMK_h);
#pragma omp atomic update
		this->Pm_low_res[im_l]+=1./(vmax*deltaMK_l);
	      }
	  }	
      }
#ifdef _USE_PDF_PHOT2SPEC_   //tjis is to speed up the code, ut has to be given in accordance to parameter file
  }
#endif
  So.DONE();

  number_gal  =0;
#pragma omp parallel for reduction(+:number_gal)
  for(int i=0;i<this->number_z.size();++i)
    number_gal+=this->number_z[i];

  if(this->LF_estimator=="Vmax" || this->LF_estimator=="Vmax_o" || this->LF_estimator=="Vmax_dc")
    {
      ofstream te;
      te.open(this->output_file_LF);
      
      for(int i=0;i<this->N_bin_Mag_low_res;++i)
	if(this->number_m_low_res[i]!=0)
	  te<<this->v_Magnitude_low_res[i]<<"\t"<<this->Pm_low_res[i]<<"  "<<this->Pm_low_res[i]/sqrt(this->number_m_low_res[i])<<endl;
      te.close();
      So.message_screen("Writting output file ",this->output_file_LF);
      So.DONE();
      //this->File.write_to_file(this->output_file_NM,this->v_Magnitude_low_res,this->number_m_low_res);

    }
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_LF_Cole2(){
  So.message_screen("Using iteration method by Cole (2012)");
  omp_set_num_threads(omp_get_max_threads());
#ifdef  _USE_PDF_PHOT2SPEC_
  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  gsl_rng_default_seed=785488;
  T = gsl_rng_ranlux;
  r = gsl_rng_alloc (T);
#endif
  
  time_t start_all;
  time(&start_all);
  time_t end;


  real_prec  deltaMK_h=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag;
  real_prec  deltaMK_l=(this->MK_max-this->MK_min)/(real_prec)this->N_bin_Mag_low_res;


  real_prec  factor=1.0;
  
  this->Delta_it.clear();
  this->Delta_it.resize(this->N_bin_z,1.0);
  ULONG ngal_lf=0;
     
  for(int IT=0;IT<this->N_iterations;++IT)
    {

      if(IT==0)
	{
	  So.message_screen("Iteration ",IT," uses standard Vmax as initial guess");
	  this->get_LF_Vmax(ngal_lf);
	}
      else  // The first iteration is done with the standard form
	{
          So.message_screen("Iteration ",IT);

	  this->number_z.clear();
	  this->number_z.resize(this->N_bin_z,0);
	  
	  this->Pm.clear();
	  this->Pm.resize(this->N_bin_Mag,0);
	  
	  this->Pm_low_res.clear();
	  this->Pm_low_res.resize(this->N_bin_Mag_low_res,0);
	  
	  this->number_m.clear();
	  this->number_m.resize(this->N_bin_Mag,0);
	  
	  this->number_m_low_res.clear();
	  this->number_m_low_res.resize(this->N_bin_Mag_low_res,0);

  
	  int NTHREADS = omp_get_max_threads();
	  omp_set_num_threads(NTHREADS);
	  vector<int>vseeds(NTHREADS,0);   
	  for(int i=0;i<vseeds.size();++i)vseeds[i]=12+105*i+IT;   
	  
	  
          //For this loop in the galaxies, we use the quantity Gal_inside_intervals whiwh tells
          //whether a galaxy is inside the itervals of Mk, z , mk and color
#ifdef  _USE_PDF_PHOT2SPEC_
#pragma omp parallel private(T, r)
	  {
	    
	    gsl_rng_env_setup();
	    gsl_rng_default_seed=vseeds[omp_get_thread_num()];
	    T = gsl_rng_ranlux;
	    r = gsl_rng_alloc (T);

	    
#pragma omp for
#else
#pragma omp parallel for
#endif
          for(int i=0;i<this->NGAL;++i)
	      {
                if(this->Gal_inside_intervals[i]>0) //This was assigned in Vmax_LF, where Mk is alrady computed and intervals asked
		  {


#ifdef  _USE_PDF_PHOT2SPEC_
                    real_prec  zgal=this->prop[i_z+i*this->NCOLS];
                    // If p_s or ps_s, assign spectroscopic redshift to this photometric ones
                    // Using the PDF calibrated in paper I
                    if(this->redshift_type=="ps_s" || this->redshift_type=="p_s")
                      {
                        // Width of the Gaussian
                        real_prec  zvpdf=0.03*tanh(-20.78*zgal*zgal+7.76*zgal+0.05)/(1.0+zgal);
                        real_prec  dist=gsl_ran_gaussian(r,zvpdf);
                        zgal+=dist;
                      }
#endif
		    // Here we compute the density weighted Vmax.
		    // The function Sij already contains the information on the limits in redshift
		    // for the given magnitde and limiting flux, so we do the intergral up to the maximum index
            // instead of computing zmax from cosmology.zmax() and get the index up to which izh must go.
                    // This also generates a factor of 1 below, so the LF so computed predicts the right number of galaxies.

                    int im_h=this->index_Mk_Gal[i];

                    real_prec  vmax=0;
#pragma omp parallel for reduction (+:vmax)
                    for(int izh=0;izh< this->Vshell.size() ;++izh)
                      vmax+=this->Vshell[izh]*this->Delta_it[izh]*this->S[izh][im_h];

            real_prec  control= isnan(vmax) || vmax==0 ? 0.0: factor;
            real_prec  vmax_mm = vmax ==0? 1.0 : vmax;

                    int iz=this->index_z_Gal[i];
#pragma omp atomic update
		    this->number_z[iz]++;
#pragma omp atomic update
		    this->number_m[im_h]++;

                    int im_l=this->index_Mk_low_Gal[i];
#pragma omp atomic update
		    this->number_m_low_res[im_l]++;

		    //				    cout<<control<<" "<<vmax<<"  "<<deltaMK_h<<endl;
#pragma omp atomic update
		    this->Pm[im_h]+=control/(vmax_mm*deltaMK_h);
#pragma omp atomic update
		    this->Pm_low_res[im_l]+=control/(vmax_mm*deltaMK_l);

	      }	  
          }
	  
          for(int j=0;j<this->N_bin_Mag;++j)
	    if(isnan(Pm[j]))Pm[j]=0;
	  
	  //GET THE NUMBER OF OBSERVED OBJECTS AS GIVEN BY THE LUM FUNCTION
          real_prec  nobs=0.;
#pragma omp parallel for collapse(2) reduction(+:nobs)
	  for(int i=0;i<this->N_bin_z;++i)
	    for(int j=0;j<this->N_bin_Mag;++j)
	      nobs+=this->Vshell[i]*this->Delta_it[i]*this->Pm[j]*this->S[i][j]*deltaMK_h;
	  
          ULONG nobs0=0;
#pragma omp parallel for reduction(+:nobs0)
	  for(int i=0;i<this->N_bin_z;++i)
	    nobs0+= this->number_z[i];
	  
          So.message_screen(" Observed  Number of observed objects =",nobs0);
          So.message_screen(" Predicted Number of observed objects =",nobs);
          factor=1;//(static_cast<real_prec>(nobs))/(static_cast<real_prec>(nobs0));

	} //end if IT >0
      // ****************************************
      
      // get the nbar from LF

      vector<real_prec>v_nbar(this->N_bin_z,0);
#pragma omp parallel for
      for(int i=0;i<this->N_bin_z;++i)
	for(int j=0;j<this->N_bin_Mag;++j)
	  v_nbar[i]+=this->Pm[j]*this->S[i][j]*deltaMK_h;
      
#pragma omp parallel for
      for(int i=0;i<this->N_bin_z;++i)
	this->Delta_it[i]= (this->number_z[i]==0 || v_nbar[i]== 0 || this->Vshell[i]==0) ? 1.0:  factor*(this->number_z[i]/this->Vshell[i])/v_nbar[i];

      // for(int i=0;i<this->N_bin_z;++i)
      // 	cout<<YELLOW<<i_z_limit<<"   "<<IT<<"  "<<this->Delta_it[i]<<"  "<<this->number_z[i]<<"  "<<this->Vshell[i]<<"   "<<v_nbar[i]<<RESET<<endl;
      
      
      v_nbar.clear();
      
      // **************** write to file***********
      if(IT>0)
	{
	  ofstream te, ta;
	  string file= IT==N_iterations-1? this->output_file_LF: this->output_file_LF+"it_"+to_string(IT);
	  te.open(file.c_str());
	  
	  cout<<BLUE<<"Writing file  "<<CYAN<<file<<RESET<<endl;
	  
	  for(int i=0;i<this->N_bin_Mag_low_res;++i)
	    if(this->number_m_low_res[i]!=0)
	      te<<this->v_Magnitude_low_res[i]<<"\t"<< this->Pm_low_res[i]<<"  "<<this->Pm_low_res[i]/sqrt(this->number_m_low_res[i])<<endl;
	  
	  te.close();
	  
	  
	  
	  string file2="Luminosity_Function/delta"+this->redshift_type+"_mk" +to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+"_"+this->hemisphere+"_it"+to_string(IT)+".txt";
	  ta.open(file2.c_str());
	  cout<<BLUE<<"Writing file  "<<CYAN<<file2<<RESET<<endl;

      vector<double>ddelta(this->Delta_it.size(), 0);
      for(int i =0;i<ddelta.size();++i)ddelta[i]=static_cast<double>(this->Delta_it[i]);
      for(int i=0;i<this->N_bin_z_low_res;i++)
	    {
          real_prec  dmm=gsl_inter_new(this->zn,ddelta, this->zn_low_res[i]);
	      ta<<this->zn_low_res[i]<<"\t"<<dmm<<endl;
	    }
	  ta.close();
	}
    }
#ifdef  _USE_PDF_PHOT2SPEC_
}
#endif
  
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_MK(real_prec  z, real_prec  mk, real_prec  &MKS)
{
  real_prec  DMk=gsl_inter_new(this->zv, this->Dm, z);
  real_prec  Kcorr;
  
  if(true==this->use_K_correction)
    {
      Kcorr= type_of_K_correction=="analytical" ?  this->cosmology.K_correction(z): gsl_inter_new(this->v_z_Kcorr,this->v_Kcorr,z);
    }
  else
    Kcorr=0;
  
  real_prec  e_corr= this->use_e_correction ? this->cosmology.e_correction(z): 0.0 ;
  DMk+=Kcorr+e_corr;
  MKS=mk-DMk;  //Compute M=m-(DM+Kcorreection)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_MK(real_prec  z, real_prec  mk, real_prec  color1, real_prec  color2,real_prec  &MKS)
{
  real_prec  DMk=gsl_inter_new(this->zv, this->Dm, z);
  real_prec  Kcorr;

  if(true==this->use_K_correction)
    {
      Kcorr= type_of_K_correction=="analytical" ?  this->cosmology.K_correction(z,color1, color2): gsl_inter_new(this->v_z_Kcorr,this->v_Kcorr,z);
    }
  else
    Kcorr=0;

  real_prec  e_corr= this->use_e_correction ? this->cosmology.e_correction(z): 0.0 ;
  DMk+=Kcorr+e_corr;
  MKS=mk-DMk;  //Compute M=m-(DM+Kcorreection)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_LF_Cole()
{
  
  So.message_screen("Measuring luminosity function using the Vmax_density method");

  time_t start_all;
  time(&start_all);
  time_t end;

  omp_set_num_threads(omp_get_max_threads());

  this->Delta_it.resize(this->N_bin_z,1);
  this->nbarq.resize(this->N_bin_z,0);

  real_prec  deltaMK_h=(this->MK_max-this->MK_min)/static_cast<real_prec>(this->N_bin_Mag);

  // -----------------------------------------------------
  // Get the LF from the /Vmax. (this->Pm)
  // This function also gives the Sij matrix, the number_z, and the number_m

  ULONG ngal_lf;
  this->get_LF_Vmax(ngal_lf);
  // -----------------------------------------------------

  real_prec  factor=1.00;
  

  vector<real_prec>Vm(Vmax_v.size(),9);
  for(int i =0;i<Vm.size();++i)Vm[i]=Vmax_v[i];
  for(int IT=0;IT<this->N_iterations;++IT)
    {
     // Get estimate of mean number density from luminosity function
      for(int i=0;i<this->N_bin_z;++i)
	nbarq[i]=0;
      
#pragma omp parallel for
      for(int i=0;i<this->N_bin_z;++i)
	for(int j=0;j<this->N_bin_Mag;++j)
	  nbarq[i]+= this->Pm[j]*this->S[i][j]*deltaMK_h;
      
      // Get new values of delta
#pragma omp parallel for
      for(int i=0;i<dNdz.size();++i)
    Delta_it[i]= factor*((static_cast<real_prec>(this->number_z[i]))/(this->Vshell[i]))/this->nbarq[i];
      
      
      // Get Vmax in Magnitude bins from values of delta
      // and luminosity function updated
#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	Vmax_v[j]=0;

#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	for(int i=0;i<this->N_bin_z;++i)
	  this->Vmax_v[j]+=Delta_it[i]*this->Vshell[i]*this->S[i][j]*deltaMK_h;


      // vector<real_prec>Mu;
      // vector<real_prec>Res;
      // for(int im=100;im>=0;--im)
      // 	{
      // 	  real_prec  mu=0.0+(20.0)*im/100.;
      // 	  real_prec  res=0;
      // 	  for(int j=0;j<this->N_bin_Mag;++j)
      // 	    res+=this->Vmax_v[j]/(this->Vmax_v[j]+mu*Vm[j]);
      // 	  Res.push_back(res);
      // 	  Mu.push_back(mu);
      // 	  //	  cout<<mu<<"  "<<res<<endl;
      // 	} 

      real_prec  Muu=0;//gsl_inter_new(Res,Mu,1.0);

#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	this->Pm[j]=0;

#pragma omp parallel for
      for(int j=0;j<this->N_bin_Mag;++j)
	this->Pm[j]=(Vmax_v[j]==0 ? 0.0 : factor*this->number_m[j]/(this->Vmax_v[j] +Muu*Vm[j]));
      
      //      real_prec  nobs=0.;
      // for(int i=0;i<this->N_bin_z;++i)
      // 	for(int j=0;j<this->N_bin_Mag;++j)
      // 	  nobs+=this->Vshell[i]*this->Delta_it[i]*this->Pm[j]*this->S[i][j]*deltaMK_h;
      
      //      factor=(static_cast<real_prec>(nobs))/(static_cast<real_prec>(this->Ngal_dndz));
      
      if(IT==N_iterations-1)
	{
	  ofstream te, ta;
	  string file=this->output_file_LF;
	  te.open(file.c_str());
	  
	  cout<<BLUE<<"Writing file  "<<CYAN<<file<<RESET<<endl;

#ifdef SINGLE_PREC
      vector<double>dPm(Pm.size(),0);
      for(int i =0;i<Pm.size();++i)dPm[i]=static_cast<double>(Pm[i]);
      vector<double>dMag(v_Magnitude.size(),0);
      for(int i =0;i<dMag.size();++i)dMag[i]=static_cast<double>(this->v_Magnitude[i]);
#endif

      for(int i=0;i<this->N_bin_Mag_low_res;++i)
	    {
	      if(this->number_m_low_res[i]!=0)
		{
#ifdef SINGLE_PREC
             real_prec  pmm=gsl_inter_new(dMag,dPm,this->v_Magnitude_low_res[i]);
#elif
              real_prec  pmm=gsl_inter_new(this->v_Magnitude,Pm,this->v_Magnitude_low_res[i]);
#endif
              te<<this->v_Magnitude_low_res[i]<<"\t"<<pmm<<"  "<<pmm/sqrt(this->number_m_low_res[i])<<endl;
		}
	    }
	  te.close();
	  
	  string file2="Luminosity_Function/delta_it"+to_string(IT)+".txt";
	  
	  ta.open(file2.c_str());
      vector<double>ddelta(this->Delta_it.size(), 0);
      for(int i =0;i<ddelta.size();++i)ddelta[i]=static_cast<double>(this->Delta_it[i]);

	  cout<<BLUE<<"Writing file  "<<CYAN<<file2<<RESET<<endl;
	  for(int i=0;i<this->N_bin_z_low_res;i++){
        real_prec  dmm=gsl_inter_new(this->zn,ddelta, this->zn_low_res[i]);
	    ta<<this->zn_low_res[i]<<"\t"<<dmm<<endl;
	    //	    cout<<this->zn_low_res[i]<<"\t"<<dmm<<endl;
	    
	  }
	  ta.close();
	}
   }
  time(&end);
  real_prec  lapse=difftime(end,start_all);
  std::cout<<BLUE<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
  time(&start_all);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GALAXY::get_dndz_LF()
{
  omp_set_num_threads(omp_get_max_threads());
  So.message_screen("Measuring dNdz from LF");
  this->dNdz_lf.resize(this->N_bin_z,0);
  real_prec  deltaMK=(this->MK_max-this->MK_min)/static_cast<real_prec>(this->N_bin_Mag);
#pragma omp parallel for
  for(int i=0;i<this->dNdz_lf.size();++i)
    this->dNdz_lf[i]=0;

#pragma omp parallel for
  for(int i=0;i<this->N_bin_z;i++)
    for(int j=0;j<=this->N_bin_Mag;++j)
      this->dNdz_lf[i]+=this->Pm[j]*this->S[i][j]*this->Vshell[i]*deltaMK;


  So.DONE();

  real_prec  acount=0.;
#pragma omp parallel for reduction(+:count)
  for(int i=0;i<this->dNdz_lf.size();i++)
    acount+=static_cast<real_prec>(this->dNdz_lf[i]);

  real_prec  Nfact=(static_cast<real_prec>(this->Ngal_dndz))/(static_cast<real_prec>(acount));

  //  for(int i=0;i<this->dNdz_lf.size();i++)this->dNdz_lf[i]*=Nfact;

  acount=0;
#pragma omp parallel for reduction(+:acount)
  for(int i=0;i<this->dNdz_lf.size();i++)
    acount+=this->dNdz_lf[i];

  So.message_screen("Number of objects LF = ", static_cast<ULONG>(acount));
  So.message_screen("Observed number of objects = ", this->Ngal_dndz);

#pragma omp parallel for
  for(int i=0;i<this->dNdz_lf.size();i++)
    this->dNdz_lf[i]/=static_cast<real_prec>(deltaz);
  So.DONE();

}


// ######################################################################
// ######################################################################

struct aux_s{
  vector<real_prec>znv;
  vector<real_prec>dndzv;

};

// ######################################################################
// ######################################################################

gsl_real idndz_lf(gsl_real z, void *p){
  aux_s *as = (aux_s *)p;
  if(z<as->znv[0])
    z=as->znv[0];

  if(z>as->znv[as->znv.size()-1])
    z=as->znv[as->znv.size()-1];

  return gsl_inter_new(as->znv, as->dndzv, z);
}

// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################

void GALAXY::get_dndz_LF_bins(){

  this->dNdz_lf_low_res.resize(this->N_bin_z_low_res,0);

  aux_s aux_s;
  for(int i =0;i<zn.size();++i)aux_s.znv[i]=static_cast<double>(zn[i]);
  for(int i =0;i<dNdz_lf.size();++i)aux_s.dndzv[i]=static_cast<double>(dNdz_lf[i]);


  int nz=(int)((real_prec)this->N_bin_z/(real_prec)this->N_bin_z_low_res);
  for(int i=0;i<this->dNdz_lf_low_res.size();i++)
    this->dNdz_lf_low_res[i]=gsl_integration(nz,idndz_lf,(void *)&aux_s,zn_low_res[i]-0.5*deltaz_low_res, zn_low_res[i]+0.5*deltaz_low_res)/deltaz_low_res;

  // Divide to be written in a file and plotted
  //  for(int i=0;i<this->dNdz_lf_low_res.size();i++)this->dNdz_lf_low_res[i]/=deltaz_low_res;

  //this->File.write_to_file(this->output_file_dndz_lf,this->zn_low_res,this->dNdz_lf_low_res);


}



// ######################################################################
// ######################################################################
// ######################################################################v
void GALAXY::get_new_mask(string coord){


  vector<vector<real_prec> > mask_temp;
  this->File.read_file(this->name_file_mask, mask_temp);
  
  ofstream mnorth, msouth;
  
  //  Here we careate a mask for south and north equatorial hemispheres
  if(coord=="equatorial"){
    mnorth.open("../MASK/new_mask_north_equ.dat");
    msouth.open("../MASK/new_mask_south_equ.dat");
    int tnorth, tsouth;

    for(int i=0;i<mask_temp.size();++i){
      real_prec  b= mask_temp[i][i_bpix];  //convert to galactic latitude
      real_prec  la= mask_temp[i][i_lpix];
      real_prec  alpha,  delta;
      galactic_to_equatorial(b,la,alpha,delta);
      if(delta>0){
        tsouth=0;
        if(mask_temp[i][this->i_mask_flag]==1)tnorth=1;
        else if(mask_temp[i][this->i_mask_flag]==0)tnorth=0;
      }
      else if(delta<0){
        tnorth=0;
        if(mask_temp[i][this->i_mask_flag]==1)tsouth=1;
        else if(mask_temp[i][this->i_mask_flag]==0) tsouth=0;
      }
      mnorth<<i<<"\t"<<la<<"\t"<<b<<"\t"<<tnorth<<endl;
      msouth<<i<<"\t"<<la<<"\t"<<b<<"\t"<<tsouth<<endl;
    }
    mnorth.close();
    msouth.close();
  }


  //  Here we careate a mask for south and north galactic hemispheres
  else if (coord=="galactic"){
    int tnorth, tsouth;

    mnorth.open("../MASK/new_mask_north_galactic.dat");
    for(int i=0;i<mask_temp.size();++i){
      real_prec  la= mask_temp[i][i_lpix];
      real_prec  b= mask_temp[i][i_bpix];
      if(b>0)tnorth=mask_temp[i][this->i_mask_flag];
      else tnorth=0;
      mnorth<<i<<"\t"<<la<<"\t"<<b<<"\t"<<tnorth<<endl;
    }
    mnorth.close();


    msouth.open("../MASK/new_mask_south_galactic.dat");
    for(int i=0;i<mask_temp.size();++i){
      real_prec  b= mask_temp[i][i_bpix];
      real_prec  la= mask_temp[i][i_lpix];
      if(b<0)tsouth=mask_temp[i][this->i_mask_flag];
      else tsouth=0;
      msouth<<i<<"\t"<<la<<"\t"<<b<<"\t"<<tsouth<<endl;
    }


    msouth.close();
  }



  mask_temp.clear();



}

// ######################################################################
// ######################################################################
// ######################################################################



void GALAXY::get_new_cat(){
  time_t start;
  time(&start);   

  //  this->output_file_new_cat="2mpz_NewCat_p_AND_s_mk0_SLICE.dat";

  So.message_screen("New cat with tabulated nbar");
  ofstream vdina;
  vdina.precision(12);
  vdina.setf(ios::showpoint);
  vdina.setf(ios::scientific);
  vdina.width(6);
  vdina.open(this->output_dir+"pinocchio.Miriam_040.plc_new.txt");
  int count=0;


/*
  
  // This uses the K-correction used by Enzo
  //  real_prec  DMl_kc=gsl_inter_new(this->zv, this->Dm, z_min_vls)+cosmology.K_correction(z_min_vls,(void *)&(this->scp))+cosmology.e_correction(z_min_vls,(void *)&(this->scp));

  // This uses the K-correction privided by J Thomas
  real_prec  DMl_kc=gsl_inter_new(this->zv, this->Dm, z_min_vls)+gsl_inter_new(this->v_z_Kcorr, this->v_Kcorr,z_min_vls)+cosmology.e_correction(z_min_vls,(void *)&(this->scp));

  int izz=(this->i_magnitude_limit==6? i_zs:i_z);
  Healpix_Map<real_prec>mask_aux(log2(nside), RING);
  for(int i=0;i<n_pixels;i++)mask_aux[i]=this->mask[this->i_mask_flag+i*this->NCOLS_MASK];



  for(int i=0;i<this->NGAL;i++){
    if(this->gal_mask[i]==this->observed_pixels_in_mask)
      {
	if(this->prop[i_z+i*this->NCOLS]>=this->z_min  && this->prop[i_z+i*this->NCOLS]<=this->z_max)
	  {
	    if(this->prop[this->i_mK+i*this->NCOLS]>=this->mK_min && this->prop[this->i_mK+i*this->NCOLS]<=this->mK_max)
	      {
        real_prec  JKcolor=(this->prop[this->i_mJ+i*this->NCOLS]-this->prop[this->i_mK+i*this->NCOLS])-gsl_inter_new(this->v_z_Kcorr,this->v_color_correction_Kcorr,prop[i_z+i*this->NCOLS]);
		//if(this->prop[i][i_zs]>=this->zzns){
        real_prec  MKs;
		get_MK(this->prop[i_z+i*this->NCOLS], this->prop[i_mK+i*this->NCOLS], MKs); //Compute ABsolute magnitude
		if(iz_vls==0)
		  { //full sample
		    count++;
            //              real_prec  alpha=this->prop[i][this->i_lgal]+180.0;real_prec  alpha2;
		    //             if(alpha>=360)alpha2=alpha-360.0; else alpha2=alpha;
		    //            vdina<<alpha2<<"\t"<<this->prop[i][this->i_bgal]<<"\t"<<prop[i][i_z]<<"\t"<<this->prop[i][this->i_mK]<<"\t"<<MKs<<"\t"<<JKcolor<<"\t"<<endl;
            real_prec  zzm = this->prop[i_z+i*this->NCOLS]<zn[0] ? zn[0]: this->prop[i_z+i*this->NCOLS];
            real_prec  nb=gsl_inter_new(zn,nbar,zzm);
		    

                    //  if((this->prop[i][this->i_dec]<=10 && this->prop[i][this->i_dec]>=-10) && (this->prop[i][this->i_ra]<=225 && this->prop[i][this->i_ra]>=150)){
                    //  if(prop[i][i_zs]>=0){
                    //  real_prec  x,y,z,xp,yp,zp;
                    //  GO.equatorial_to_cartesian(this->prop[i][this->i_ra],0,prop[i][i_zp],xp,yp,zp);
                    //  GO.equatorial_to_cartesian(this->prop[i][this->i_ra],0,prop[i][i_zs],x,y,z);
                    //  if(x*x+y*y<=pow(0.2,2))vdina<<x<<"\t"<<y<<"\t"<<xp<<"\t"<<yp<<"\t"<<prop[i][i_zs]<<"\t"<<prop[i][i_zp]<<"\t"<<nb<<endl;
                    //  }
                    //  }

		    
		    vdina<<this->prop[this->i_lgal+i*this->NCOLS]<<"\t"<<this->prop[this->i_bgal+i*this->NCOLS]<<"\t"<<prop[i_zp+i*this->NCOLS]<<"\t"<<prop[i_zs+i*this->NCOLS]<<"\t"<<prop[this->i_mK+i*this->NCOLS]<<"\t"<<MKs<<"\t"<<JKcolor<<"\t"<<nb<<endl;
		    
		  }
		else{//semi VLS
		  if(iz_vls>0){
		    if(this->prop[i_z+i*this->NCOLS]>z_min_vls || (this->prop[i_z+i*this->NCOLS]<=z_min_vls && MKs<=this->mK_max-DMl_kc)  ){
		      count++;
              real_prec  nb=gsl_inter_new(zn,nbar,this->prop[i_z+i*this->NCOLS]);
		      vdina<<this->prop[this->i_lgal+i*this->NCOLS]<<"\t"<<this->prop[this->i_bgal+i*this->NCOLS]<<"\t"<<prop[i_z+i*this->NCOLS]<<"\t"<<nb<<"\t"<<this->prop[this->i_mK+i*this->NCOLS]<<"\t"<<MKs<<"\t"<<JKcolor<<"\t"<<endl;
		    }
		  }
		}
		//}
	      } 
	  }
      }
    
  }  

*/


 gsl_interp_accel *acc;
 acc = gsl_interp_accel_alloc ();
 gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, nbar_new.size());
 gsl_spline_init (spline, &(zn_new[0]), &(nbar_new[0]), nbar_new.size());


 this->dNdM.clear();
 this->dNdM.shrink_to_fit();
 this->dNdM.resize(this->N_bin_lmass,0);

 real_prec  deltam=(lmass_max-lmass_min)/static_cast<real_prec>(N_bin_lmass);

  for(int i=0;i<this->NGAL;i++)
      {
        real_prec  zr= this->prop[this->i_z+i*this->NCOLS];
        real_prec  nb=0;
        int iz_low=get_bin(zr,z_min_low_res,N_bin_z_low_res,deltaz_low_res,true);
        if(zr<=this->z_max && zr>=this->z_min)
            {
                if(zr>this->zn_new[zn_new.size()-1])
                     nb = this->nbar_new[zn_new.size()-1];
                else if(zr<=this->zn_new[0])
                    nb = this->nbar_new[0];
                else
                   nb=  static_cast<real_prec>(gsl_spline_eval (spline, zr, acc));
#ifdef _USE_MASS_
                int indexm=get_bin(log10(this->prop[this->i_mass+i*this->NCOLS]),lmass_min,N_bin_lmass, deltam, true);
           //    real_prec  DeltaM= pow(10, lmass_min+(indexm+1)*deltam)-pow(10, lmass_min+(indexm)*deltam);


                dNdM[indexm]++;
                ULONG idmz=index_2d(indexm,iz_low,N_bin_z_low_res);
                real_prec  alpha=this->NXY[idmz]/(this->dNdz_low_res[iz_low]*deltaz_low_res);// this is N(z,M)/N(z) such that alpha * nb sim nbar(z,M)
                vdina<<static_cast<float>(prop[i_z+i*this->NCOLS])<<"\t"<<static_cast<float>(this->prop[this->i_ra+i*this->NCOLS])<<"\t"<<static_cast<float>(this->prop[this->i_dec+i*this->NCOLS])<<"\t"<<static_cast<float>(this->prop[this->i_mass+i*this->NCOLS])<<"\t"<<static_cast<float>(nb)<<"\t"<<static_cast<float>(nb*alpha)<<endl;

#else
                vdina<<prop[i_z+i*this->NCOLS]<<"\t"<<this->prop[this->i_ra+i*this->NCOLS]<<"\t"<<this->prop[this->i_dec+i*this->NCOLS]<<"\t"<<nb<<endl;

#endif


            count++;
            }
            }
  vdina.close();


  vector<real_prec>mass_aux(N_bin_lmass,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif

  for(auto i=0;i<dNdz_low_res.size();i++)
      this->dNdM[i]/=(deltam*static_cast<real_prec>(count));

  for(auto i=0;i<dNdz_low_res.size();i++)
        mass_aux[i]=lmass_min+(i+0.5)*deltam;
  //this->File.write_to_file(this->output_dir+"dndM.txt",mass_aux,dNdM);




  std::cout<<"with "<<count<<" objects"<<endl;
}


  
  
// ######################################################################
// ######################################################################

// void GALAXY::get_alpha_rsd(){


//     vector<vector<real_prec> >propa;
//     ofstream daa;daa.open("alpha_rsd_zbin1.txt");
//     Fo.read_file("/home/andres/data/Numerics/EUCLID/CODES/GitHub/Cl/2MPZ_dndz_zbin1_gaussian_sigmaz_var.dat", propa);
//     int nzz=propa.size();


//     vector<real_prec> alpha_rsd(nzz+1,0);
//     vector<real_prec> dnbar(nzz+1,0);
//     vector<real_prec> nbaR(nzz+1,0);
//     vector<real_prec> DNDZ(nzz+1,0);
//     vector<real_prec> ZZ(nzz+1,0);
//     for(int i=0;i<nzz;i++)ZZ[i]=propa[0+i*this->NCOLS];
//     for(int i=0;i<nzz;i++)DNDZ[i]=propa[i][2];


//     for(int i=0;i<nzz;i++)nbaR[i]=DNDZ[i]/gsl_inter_new(this->zv,this->Vshell,ZZ[i]);

//     for(int i=0;i<nzz;++i){
//       dnbar[i]=(nbaR[i+1]-nbaR[i])/(ZZ[1]-ZZ[0]);
//     }


//     for(int i=0;i<nzz;++i){
//       if(DNDZ[i]!=0){
//         real_prec  rr=gsl_inter_new(this->zv,this->rv,ZZ[i]);
//         alpha_rsd[i]=(dnbar[i]/nbaR[i])*gsl_inter_new(this->zv,this->Hv,ZZ[i])/Constants::speed_light   + 2./rr;
//         daa<<ZZ[i]<<"  "<<rr<<"  "<<alpha_rsd[i]<<endl;
//         cout<<ZZ[i]<<"  "<<rr<<"  "<<nbaR[i]<<"  "<<alpha_rsd[i]<<endl;

//       }
//     }
//     daa.close();

// }

// ######################################################################
// ######################################################################


// void GALAXY::get_new_cat_bcg(string file, string cat){

//   time_t start;
//   time(&start);
  
//   ofstream vdina;
//   vdina.precision(12);
//   vdina.setf(ios::showpoint);
//   vdina.setf(ios::scientific);
//   vdina.width(6);
//   vdina.open(file.c_str());
//   int count=0;
//   real_prec  aa=((real_prec)prop.size())/((real_prec)prop_r.size());
//   if(cat=="g"){
//     for(int i=0;i<this->prop.size();++i){
//       count++;
//       real_prec  nb=aa*gsl_inter_new(this->zn,this->nbar,this->prop[i_z+i*this->NCOLS]);
//       vdina<<this->prop[this->i_ra+i*this->NCOLS]<<"\t"<<this->prop[this->i_dec+i*this->NCOLS]<<"\t"<<prop[i_z+i*this->NCOLS]<<"\t"<<nb<<endl;
//     }
//   }
//   else{
//     for(int i=0;i<this->prop_r.size();++i){
//       count++;
//       real_prec  nb=aa*gsl_inter_new(this->zn,this->nbar,this->prop_r[i][i_z]);
//       vdina<<this->prop_r[i][this->i_ra]<<"\t"<<this->prop_r[i][this->i_dec]<<"\t"<<prop_r[i][i_z]<<"\t"<<nb<<endl;
//     }
//   }
//   vdina.close();
//   std::cout<<BLUE<<"New cat written in file "<<file<<endl;
//   std::cout<<"with "<<count<<" objects"<<RESET<<endl;
// }



// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
#define _USE_OMP_rangen_
  
void GALAXY::get_random_cat(string file, long n_random){
  




  real_prec  lka=-100.;
  real_prec  lkb=0.;
  for(int i=0;i<this->nbar.size();++i)
     {
        lkb=max(static_cast<real_prec>(this->nbar[i]), lka);
        lka=lkb;
     }

   this->prob.clear(); prob.shrink_to_fit();
   this->prob.resize(this->nbar.size(),0);
   for(int i=0;i<this->nbar.size();++i)
        this->prob[i]=this->nbar[i]/lka;

  std::cout<<"zmin = "<<z_min<<endl;
  std::cout<<"zmax = "<<z_max<<endl;

  deltaz=this->zn_new[1]-this->zn_new[0];// this is for the dNdz computed on the fly

  gsl_interp_accel *acc;

  acc = gsl_interp_accel_alloc ();
  gsl_spline *spline_zr  = gsl_spline_alloc (gsl_interp_linear, zv.size());
  gsl_spline_init (spline_zr, &(zv[0]), &(rv[0]), zv.size());

  real_prec  rmax=static_cast<real_prec>(gsl_spline_eval(spline_zr, z_max, acc)); // gsl_inter_new(zv,rv,this->z_max);
  real_prec  rmin=static_cast<real_prec>(gsl_spline_eval(spline_zr, z_min, acc));//  gsl_inter_new(zv,rv,this->z_min);
  gsl_spline_free(spline_zr);




  // get z(r) given r
  gsl_interp_accel *acc_rz;
  acc_rz = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, rv.size());
  gsl_spline_init (spline, &(this->rv[0]), &(this->zv[0]), this->rv.size());


// get P(z) given z
  gsl_interp_accel *acc_zp;
  acc_zp = gsl_interp_accel_alloc ();
  gsl_spline *spline_p    = gsl_spline_alloc (gsl_interp_linear, this->prob.size());
  gsl_spline_init (spline_p, &(this->zn_new[0]), &(this->prob[0]), this->prob.size());


  // get nbar(z) given z
  gsl_interp_accel *acc_zn;
  acc_zn = gsl_interp_accel_alloc ();
  gsl_spline *spline_n    = gsl_spline_alloc (gsl_interp_linear, this->nbar.size());
  gsl_spline_init (spline_n, &(this->zn_new[0]), &(this->nbar[0]), this->nbar.size());


  std::cout<<"r min = "<<rmin<<endl;
  std::cout<<"r max = "<<rmax<<endl;


//  this->dNdz_low_res.clear();
//  this->dNdz_low_res.shrink_to_fit();
//  this->dNdz_low_res.resize(this->N_bin_z_low_res,0);
//  this->dNdM.clear();
//  this->dNdM.shrink_to_fit();
//  this->dNdM.resize(this->N_bin_lmass,0);


#ifdef _USE_MASSa_
  real_prec  deltam=(lmass_max-lmass_min)/static_cast<real_prec>(N_bin_lmass);
#endif

#ifdef _USE_MASK_
  Healpix_Map<real_prec>map_aux(log2(nside), RING);
#endif

 real_prec  factor=M_PI/180.0;

 const gsl_rng_type * rng_t;
 gsl_rng * gBaseRand;
#ifndef _USE_OMP_rangen_
 gsl_rng_env_setup();
 gsl_rng_default_seed=224;
 rng_t = gsl_rng_ranlux;
 gBaseRand = gsl_rng_alloc (rng_t);
#endif


 int Nran_files=10;
 std::cout<<"Generating random catalog with "<<n_random*Nran_files<<" objects in chuncks of"<<Nran_files<<" files"<<endl;


#ifdef _USE_OMP_rangen_
 int jthread;

  omp_set_num_threads(Nran_files);
 vector<ULONG>vseeds(Nran_files,0);
 for(int i=0;i<vseeds.size();++i)
   vseeds[i]=3+static_cast<ULONG>(i+27*i*i)*56145;
#endif

#ifdef _USE_OMP_rangen_
#pragma omp parallel private (jthread, gBaseRand, rng_t)
   {


#pragma omp for
#endif
     for(int ir=0;ir<Nran_files;++ir)
     {
         string new_file=file+"_ranfile"+to_string(ir);
         ofstream vdina;
         vdina.open(new_file.c_str());
         vdina.precision(6);
         vdina.setf(ios::showpoint);
         vdina.setf(ios::scientific);
         vdina.width(3);

         vector<real_prec>dNdz_random(this->dNdz_new.size(),0);

         jthread=omp_get_thread_num();
         gsl_rng_default_seed=vseeds[jthread];
         rng_t = gsl_rng_mt19937;//_default;
         gBaseRand = gsl_rng_alloc (rng_t);

         ULONG count=0;

         while(count<n_random)
        {

        real_prec  rd  = rmax*pow(gsl_rng_uniform(gBaseRand), 1./3.);

        if(rd >=rmin)
         {
          real_prec  ra = (2.*M_PI)*gsl_rng_uniform(gBaseRand);       // ra
          real_prec  theta= acos(-1.0+2.0*gsl_rng_uniform(gBaseRand));     // theta
          real_prec  dec= 0.5*M_PI-theta;                          //delta
    #ifdef _USE_MASK_
          long ipix;
          point.phi=lg;
          point.theta=0.5*M_PI-bg;
          ipix=map_aux.ang2pix(point);
          if(mask[this->i_mask_flag+ipix*NCOLS_MASK]==1)
              {
    #else
           if((ra<=this->phi_max*factor && ra>= this->phi_min*factor) && (dec<=this->dec_max*factor && dec>=this->dec_min*factor) )
              {
#endif
               real_prec  zr   =  static_cast<real_prec>(gsl_spline_eval (spline, rd, acc_rz)); //
               real_prec  proba= 0 ;
               if(zr>=this->zn_new[zn_new.size()-1])
                   proba = this->prob[zn_new.size()-1];
               else if(zr<=this->zn_new[0])
                  proba = this->prob[0];
               else
                 proba =  static_cast<real_prec>(gsl_spline_eval (spline_p, zr, acc_zp)) ;
               real_prec  xr   = gsl_rng_uniform (gBaseRand);
                    if(proba>=xr)
                      {
                       real_prec  nb = 0;
                       if(zr>this->zn_new[zn_new.size()-1])
                            nb = this->nbar[zn_new.size()-1];
                       else if(zr<=this->zn_new[0])
                           nb = this->nbar[0];
                       else
                           nb= static_cast<real_prec>(gsl_spline_eval (spline_n, zr, acc_zn));

                       //     int iz_low=get_bin(zr,z_min,N_bin_z_low_res,deltaz_low_res, true);
                       int iz=    get_bin(zr,zn_new[0],zn_new.size(), deltaz, true);
                       dNdz_random[iz]+=1./deltaz;

            #ifdef _USE_MASSa_
                       bool accept=false;
                       while(false==accept)
                        {
                          real_prec  lmass=lmass_min+(lmass_max-lmass_min)*gsl_rng_uniform(ral);
                          int indexm=get_bin(lmass,lmass_min,N_bin_lmass, deltam, true);
                          ULONG id=index_2d(indexm,iz_low,N_bin_z_low_res);
                          real_prec  proba_mz=this->PXY[id];
                          real_prec  xm   = gsl_rng_uniform (gBaseRand);

                          if(proba_mz>=xm)
                          {
                            accept=true;
                            prop[i_zs+count*this->NCOLS]=static_cast<float>(zr);
                            prop[i_ra+count*this->NCOLS]=static_cast<float>(180*ra/M_PI);
                            prop[i_dec+count*this->NCOLS]=static_cast<float>(180*dec/M_PI);
                            prop[i_mass+count*this->NCOLS]=static_cast<float>(pow(10,lmass));
                            real_prec  alpha=this->NXY[id]/(this->dNdz_low_res[iz_low]*deltaz_low_res);// this is N(z,M)/N(z) such that alpha * nb \sim nbar(z,M)
                            vdina<<static_cast<float>(zr)<<"\t"<<static_cast<float>(180*ra/M_PI)<<"\t"<<static_cast<float>(180.*dec/M_PI)<<"\t"<<static_cast<float>(pow(10,lmass))<<"\t"<<"\t"<<static_cast<float>(nb)<<"\t"<<static_cast<float>(alpha*nb)<<endl;
                            dNdM[indexm]++;

                          }
                        }
    #else
               vdina<<180*ra/M_PI<<"\t"<<180.*dec/M_PI<<"\t"<<zr<<"\t"<<nb<<endl;
    #endif
               count++;
                }
              }
            }
   //   So.message_screen_flush("Number of random objects = ",static_cast<int>(count));
     } // closes while
         vdina.close();

//         this->File.write_to_file(this->output_dir+"dndz_random.txt_file"+to_string(ir), zn_new,dNdz_random);
         gsl_rng_free (gBaseRand);

     }// closes loop over files



#ifdef _USE_OMP_rangen_
   }
#endif


  this->So.DONE();


#ifdef _USE_MASSa_
  vector<real_prec>mass_aux(N_bin_lmass,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif

  for(auto i=0;i<dNdz_low_res.size();i++)
      this->dNdM[i]/=(deltam*static_cast<real_prec>(count));

  for(auto i=0;i<dNdz_low_res.size();i++)
        mass_aux[i]=lmass_min+(i+0.5)*deltam;
  //this->File.write_to_file(this->output_dir+"dndM_random.txt",mass_aux,dNdM);
#endif


  gsl_spline_free(spline);
  gsl_spline_free(spline_p);
  gsl_spline_free(spline_n);
  gsl_interp_accel_free(acc_rz);
  gsl_interp_accel_free(acc_zp);
  gsl_interp_accel_free(acc_zn);


  gsl_rng_free (gBaseRand);

  
}
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################

void GALAXY::bspline(int u,int n){
  std::cout<<CYAN<<"Smoothing dNdz and nbar  ... "<<RESET<<endl;

  int nn0=n;

  vector<gsl_real> new_zn0;
  vector<gsl_real> new_nb0;

  vector<gsl_real>nv(zn.size(),0);
  vector<gsl_real>za(zn.size(),0);
  za=zn;
  nv=(u == 1 ? this->nbar : this->dNdz);

  cout<<zn.size()<<endl;

  for(int i=0;i<1;i++)
    {
      new_zn0.resize(nn0);
      new_nb0.resize(nn0);
      gsl_bspline(za,nv, new_zn0, new_nb0);// smooth dndz
      za.clear();za.shrink_to_fit();
      za.resize(nn0,0);
      nv.clear();nv.shrink_to_fit();
      nv.resize(nn0,0);
      za=new_zn0;
      nv=new_nb0;
      nn0-=2;
    }
  
  for(int i=0;i<nv.size();i++)
      nv[i]=fabs(nv[i]);

  if(u==1)
    {
      nbar_new.clear();
      nbar_new.shrink_to_fit();
      zn_new.clear();
      zn_new.shrink_to_fit();
      nbar_new.resize(new_nb0.size(),0);
      zn_new.resize(za.size(),0);
      zn_new=za;
      for(int i=0;i<nv.size();i++)
          this->nbar_new[i]=fabs(new_nb0[i]);
      //this->File.write_to_file(this->output_dir+"nbar_smooth.txt",za,this->nbar_new);
    }
  
  if(u==2)
    {
      dNdz_new.resize(zn.size());
      zn_new.resize(zn.size());
      for(int i=0;i<zn.size();i++)this->dNdz_new[i]=gsl_inter_new(za,new_nb0,this->zn[i]);
      //this->File.write_to_file(this->output_dir+"dNdz_smooth.txt",zn,dNdz_new);
      zn_new=za;
    }
}



// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################

void GALAXY::get_nbar(){


  this->nbar.clear();
  this->nbar.shrink_to_fit();

  this->deltaz=this->zn_new[1]-this->zn_new[0];

  this->nbar.resize(dNdz_new.size(),0);
  So.message_screen("Computing Mean number density");
  for(int i=0;i<dNdz_new.size();i++)
    {
//      real_prec  dr_dz=Constants::speed_light/cosmology.Hubble_function(this->zn_new[i]);
      real_prec  distance_i=cosmology.comoving_distance(this->zn_new[i]-0.5*deltaz); // gsl_inter_new(this->zv,this->rv,zn[i]);
      real_prec  distance_j=cosmology.comoving_distance(this->zn_new[i]+0.5*deltaz); // gsl_inter_new(this->zv,this->rv,zn[i]);
//      real_prec Vshell=this->total_area*pow(distance,2)*(dr_dz)*deltaz;
      real_prec Vshell=this->total_area*(pow(distance_j,3)-pow(distance_i,3));
      this->nbar[i]=dNdz_new[i]/Vshell; // dz is not included here as dNdz does not included it neither
  }
 // string filezs=this->output_dir+this->name_cat+"_nbar_p_mk"+to_string(this->i_magnitude_limit)+"_zmax"+to_string(this->i_z_limit)+".txt";
  string filezs=this->output_dir+"nbar.txt";


 this->So.message_screen("Writting to file ", filezs);
  ofstream sal(filezs);
  for(int i=0;i<this->dNdz_new.size();++i)
    sal<<this->zn_new[i]<<"\t"<<this->nbar[i]<<endl;
 sal.close();
 So.DONE();


}


// ######################################################################
// ######################################################################

void GALAXY::bspline(string input_file,int nfactor){
  std::cout<<CYAN<<"Smoothing dNdz and nbar  ... "<<RESET<<endl;


  const gsl_rng_type * rng_t;
  gsl_rng * gBaseRand;
  gsl_rng_env_setup();


  int nlines=0;
  ifstream infile;
  infile.open(input_file.c_str());
  std::cout<<CYAN<<"Reading "<<input_file<<endl;

  // container for original input
  vector<gsl_real>nv;
  vector<gsl_real>za;

  real_prec  az, nz, aux;
  while(!infile.eof())
    {
      infile>>az>>nz>>aux;
      nv.push_back(nz);
      za.push_back(az);
  }
  za.pop_back();
  nv.pop_back();
  nlines=za.size();
  for(int i=0;i<nv.size();++i)nv[i]/=nv.size(); // Dividir histograma entre el nḿuero de bines para que al sumar de el

  ULONG Ntot=0;
  for(int i=0;i<nv.size();++i)Ntot+=nv[i];


  infile.close();
  std::cout<<CYAN<<"TOtal number of galaxies in hist "<<Ntot<<endl;
  std::cout<<CYAN<<"Done "<<nlines<<endl;

  // make a copy of the original zbins
  vector<gsl_real>za_original(za.size());
  for(int i=0;i<za.size();i++)za_original[i]=static_cast<real_prec>(za[i]); // ensure positiviness
  vector<gsl_real>nz_original(za.size());
  for(int i=0;i<nv.size();i++)nz_original[i]=static_cast<real_prec>(nv[i]); // ensure positiviness


  // containers whose size will be reduced
  vector<gsl_real> new_zn;
  vector<gsl_real> new_nb;

  int nn0=nlines; // this is key
  int fraction=static_cast<int>(floor(nlines/nfactor));

  cout<<fraction<<endl;
  nn0=fraction;

      new_zn.resize(nn0);
      new_nb.resize(nn0);

      cout<<"Smooth"<<endl;
      gsl_bspline(za,nv, new_zn, new_nb);// smooth the nv(zn) array to a new new_nb(new_zn)

      for(int i=0;i<nv.size();i++)new_nb[i]=abs(new_nb[i]); // ensure positiviness

      // now sample randomly the z interval using the real_prec of points
      cout<<"Interpol"<<endl;

      gsl_rng_default_seed=25225;
      rng_t = gsl_rng_ranlux;
      gBaseRand = gsl_rng_alloc (rng_t);

      vector<real_prec>zran(2*nlines,0);

      for(int ir=0;ir<zran.size();++ir)
          zran[ir]=this->z_min+(this->z_max-this->z_min)*gsl_rng_uniform(gBaseRand);

      vector<real_prec>zran_s(zran.size(),0);
      sort_1d_vectors<real_prec>(zran,zran_s);

      vector<real_prec>dNran(zran.size(),0);
      for(int ir=0;ir<zran_s.size();++ir)// INterpola el prodcto de lbspline en los nuevos z randoms sorted
        if(zran_s[ir]<new_zn[new_zn.size()-1] && zran_s[ir]>=new_zn[0])
            dNran[ir]=gsl_inter_pointers(&new_zn[0], &new_nb[0], new_zn.size(), zran_s[ir]);

    // bautiza de nuevo
      za.clear();za.shrink_to_fit(); za.resize(zran_s.size(),0);
      for(int i=0;i<za.size();i++)za[i]=zran_s[i]; // ensure positiviness
      nv.clear();nv.shrink_to_fit();nv.resize(zran_s.size(),0);
      for(int i=0;i<nv.size();i++)nv[i]=dNran[i]; // ensure positiviness

      this->dNdz.clear();this->dNdz.shrink_to_fit();
      this->dNdz.resize(dNdz_new.size(),0);
      this->dNdz_new.clear();dNdz_new.shrink_to_fit();
      this->zn_new.clear();this->zn_new.shrink_to_fit();
      gsl_interp_accel *acc = gsl_interp_accel_alloc ();
      gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear,  za.size());
      gsl_spline_init (spline, &za[0], &nv[0], za.size());
      for(int i=0;i<nlines;i++)
      {
        if((za_original[i]>=za[0] && za_original[i]>=this->z_min) && (za_original[i]<=za[za.size()-1] && za_original[i]<=this->z_max))
          {
             this->dNdz_new.push_back(gsl_spline_eval (spline, za_original[i], acc)) ; //gsl_inter_new(za,new_nb,za_original[i]);
             this->zn_new.push_back(za_original[i]);
           }
       }

//      this->File.write_to_file(this->output_dir+"dndz_smooth_bspline.txt",this->zn_new,this->dNdz_new);


      this->dNdz = this->dNdz_new;
  /*
//========================================
// Now we apply a second method: Fourier transform and a low pass filter:

  std::cout<<CYAN<<"Applying low pass filter "<<endl;

  ULONG Nft_HR=static_cast<ULONG>(nlines);
  ULONG Nft_LR=static_cast<ULONG>(floor(nlines/nfactor));

  complex_prec *in =  (complex_prec *)fftw_malloc(Nft_HR*sizeof(complex_prec));
  complex_prec *out = (complex_prec *)fftw_malloc(Nft_HR*sizeof(complex_prec));
  //vector<real_prec>inv(Nft_HR,0);

  for(int i = 0 ; i < Nft_HR; ++i )
  {
    in[i][REAL]=static_cast<real_prec>(nz_original[i]);
  //  in[i][REAL]=static_cast<real_prec>(i*i);
    in[i][IMAG]=0;
  //  inv[i]=static_cast<real_prec>(nz_original[i]);
  }


  do_fftw_1d(Nft_HR, in ,out, FFTW_FORWARD);
  //do_fftw_1d_r2c(Nft_HR, inv ,out);

  real_prec lenght=za_original[nlines-1]-za_original[0];
  real_prec delta_box_hr=za_original[1]-za_original[0];//   lenght/static_cast<real_prec>(Nft_HR);
  real_prec alpha=Nft_HR/Nft_LR; //static_cast<real_prec>(Nft_HR)/static_cast<real_prec>(Nft_LR);
  real_prec delta_shift= 0.5*(alpha-1.0)*delta_box_hr; //Yu

  real_prec delta_k= 2*M_PI/lenght;
  vector<real_prec> kmodes_lr(Nft_LR,0.0);

   for(int i = 0 ; i < Nft_LR; ++i )
    {
       int coords= i < Nft_LR/2? i: i-Nft_LR;
       kmodes_lr[i]=static_cast<real_prec>(coords)*delta_k;
    }

   complex_prec *out_lr = (complex_prec *)fftw_malloc(Nft_LR*sizeof(complex_prec));
   complex_prec *in_lr =  (complex_prec *)fftw_malloc(Nft_LR*sizeof(complex_prec));

  for(ULONG i=0;i<Nft_LR;++i)
    {
      ULONG i_hr = i < Nft_LR/2 +1 ? i : i+Nft_HR-Nft_LR;
      real_prec new_real=out[i_hr][REAL];
      real_prec new_imag=out[i_hr][IMAG];
      real_prec shift=  delta_shift*kmodes_lr[i];
      out_lr[i][REAL] = (cos(shift)*new_real - sin(shift)*new_imag); //This follows the convention delta_2-> exp(iks)*delta_2
      out_lr[i][IMAG] = (cos(shift)*new_imag + sin(shift)*new_real);
    }


  int ii=0; real_prec a,b;
   a  =  out_lr[ii][REAL];
   b  =  out_lr[ii][IMAG];
   out_lr[ii][REAL]=sqrt(a*a+b*b);
   out_lr[ii][IMAG]=0.0;
   ii=Nft_LR/2-1;
   a  =  out_lr[ii][REAL];
   b  =  out_lr[ii][IMAG];
   out_lr[ii][REAL]=sqrt(a*a+b*b);
   out_lr[ii][IMAG]=0.0;


   do_fftw_1d(Nft_LR, out_lr ,in_lr, FFTW_BACKWARD);
   cout<<alpha<<"  "<<delta_shift<<"   "<<za_original[0]<<endl;

   vector<gsl_real>new_z(Nft_LR,0);
   vector<gsl_real>new_dn(Nft_LR,0);
    for(int i = 0 ; i < Nft_LR; ++i )
      {
        new_z[i]=za_original[0]+(i+0.5)*delta_box_hr*alpha;
        new_dn[i]=in_lr[i][REAL]/alpha;
      }
    this->File.write_to_file(this->output_dir+"dndz_fourier.txt",new_z,new_dn);

    spline = gsl_spline_alloc (gsl_interp_linear, Nft_LR);
    gsl_spline_init (spline, &new_z[0], &new_dn[0], new_z.size());

    this->dNdz_new.clear();
    this->dNdz_new.shrink_to_fit();
    this->dNdz_new.resize(Nft_HR,0);
    for(int i=0;i<Nft_HR-1;i++)
    {
        if(za_original[i]>=new_z[0] && za_original[i]<=new_z[new_z.size()-1])
            this->dNdz_new[i]=gsl_spline_eval (spline, za_original[i], acc) ; //gsl_inter_new(za,new_nb,za_original[i]);
      }
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    cout<<this->dNdz_new.size()<<"  "<<Nft_HR<<endl;
    this->File.write_to_file(this->output_dir+"dndz_smooth_fourier"+to_string(nfactor)+".txt",za_original,this->dNdz_new);
    this->N_bin_z=this->dNdz_new.size();

*/

}
 
// ######################################################################
// ######################################################################
void GALAXY::get_jacknife_mask(){

  // REDEFINE THE WY YOU COUNT PIXELS IN THE CREATION OF MASKED PIXELS. NO ES NECESARIO CREAR MASCARA GRANDE
  // SIMPLEMENTE DECIR CUANDTOS PIXELSES ORIGINALES QUEREMOS QUITAR Y DARLE. NO NECESITAOS SER INDEPENDEIDNTES
  
  std::cout<<CYAN<<"Creating Jacknife masks"<<RESET<<std::endl;
  /// Create Nres masks, eacsh with resolution Nres and subtracting one pixel.
  /// The mask should follow 2MPZ mask.
  ofstream sal;
  
  
  //In the low resolution mask, each pixel contains Nside_new**2 / nside**2 pixels in the high res mask.
  
  
  /*
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    
    Healpix_Map<real_prec>mask_aux2(log2(this->nside), RING);
    int count=0;
    
    
    for(int count_regions=1;count_regions<=1000;++count_regions){ // loop over the low res mask
    
    vector<int>mask_pix(this->n_pixels, 0);
    for(int i=0;i<this->n_pixels;i++)mask_pix[i]=this->mask[i][this->i_mask_flag];
    
    gsl_rng_default_seed=78+5*count_regions;
    T = gsl_rng_ranlux;
    r = gsl_rng_alloc (T);
    
    pointing point_aux_high;
    count=0;
    
    do{
    init:
    real_prec  lg = (2.*M_PI)*gsl_rng_uniform(r);       // define galactic longitude
    real_prec  bg= acos(-1.0+2.0*gsl_rng_uniform(r));   // define galactic latitude (-PI, PI)
    bg =  0.5*M_PI-bg;                              // Convert to theta (0, PI)
    point_aux_high.phi=lg;
    point_aux_high.theta=0.5*M_PI-bg;
    long ipix=mask_aux2.ang2pix(point_aux_high); // It is the same for z1, z2, z3, just a map

    if(this->mask[ipix][this->i_mask_flag]==1){
    //        cout<<ipix<<endl;
    mask_pix[ipix]=0;
    count++;
    }
    else goto init;
    }while(count<256); //que hayan 64 pixeles, como si fuera un pixel grande de Nres ¡ 8

    string new_mask="../MASK/new_mask_jacknife_"+to_string(count_regions)+".txt";
    sal.open(new_mask);
    for(int ih=0;ih<this->n_pixels;++ih)sal<<ih<<"\t"<<this->mask[ih][this->i_lpix]<<"\t"<<this->mask[ih][this->i_bpix]<<"\t"<<mask_pix[ih]<<endl;
    sal.close();
    std::cout<<CYAN<<"Writing file "<<new_mask<<RESET<<endl;
    
    }
  */
  
  
  
  int Nside_new=4;
  Healpix_Map<real_prec>mask_aux_low(log2(Nside_new), RING);
  
  int n_pixels_new=12*Nside_new*Nside_new;
  real_prec  Npp=pow(this->nside,2)/pow(Nside_new,2);
  std::cout<<RED<<"In the low resolution mask, each pixel contains "<< Npp<<" pixels in the high res mask."<<RESET<<std::endl;
  
  int count_regions=0;
  
  
  string apix_file="../MASK/area_low_res_mask_jacknife_Nside_"+to_string(Nside_new)+".txt";
  ofstream apix;
  apix.open(apix_file.c_str());


  string low_res_file="../MASK/low_res_mask_Nside_"+to_string(Nside_new)+".txt";
  ofstream lres(low_res_file.c_str());

  
  vector<vector<int> > ihh;
  ihh.resize(this->n_pixels);


  vector<int> ihh_new;
  ihh_new.resize(n_pixels_new,1-observed_pixels_in_mask);

  vector<int> ihh_count;
  ihh_count.resize(n_pixels_new,-999);
  
 
  for(int il=0;il<n_pixels_new;++il) // loop over the low res mask
    { 
      
      vector<int>index(this->n_pixels,0);
      int count=0;//Counter on the number of pixels that are going to me masked
      int count_ori=0;  //Counter on the pixels in the original mask that are masked
      
      pointing point_aux_high;
      
      for(int ih=0;ih<this->n_pixels;++ih)
	{  // loop over the high resolution mask
	  point_aux_high.phi=fac*this->mask[this->i_lpix+ih*NCOLS_MASK];

	  point_aux_high.theta=0.5*M_PI-fac*this->mask[this->i_bpix+ih*NCOLS_MASK];
	  
	  int pixel_high= mask_aux_low.ang2pix(point_aux_high);
	  // get the pixel id of the high resolution map
	  
	  if(pixel_high==il)
	    {
	      index[ih]=static_cast<int>(fabs(1-this->observed_pixels_in_mask)); //If I'm in the low res pixel I want to subtract, mask all the hres pixels inside
	      
	      if(this->mask[this->i_mask_flag+ih*NCOLS_MASK]==static_cast<int>(fabs(1-this->observed_pixels_in_mask)))
		count_ori++;//count how many pixels in the original mask, contained in this low res pixel il, were already masked
	      
	      count++;  //Count how many pixels of this high_res pix are in the corresponding low_res
	    }
	  else
	    index[ih]=this->mask[this->i_mask_flag+ih*NCOLS_MASK]; // If this pixel is in other low res pixel, then assign the original value of the mask
	}
      
      
      
      int Npx_masked=count-count_ori;//number of pixels that are effectively masked from the low_res pixel il
      
      //205 256 for Nside_new = 4
      //814 1024 for Nside_new = 2
      if(Npx_masked>= 205 && Npx_masked<=256) // 10 percent variation in the area
	{
	  //if(Npx_masked>=205 && Npx_masked<=Npp){
	  //          if(Npx_masked>=0 && Npx_masked<=5555556){
	  
  
	  ihh_new[il]=this->observed_pixels_in_mask;
      
	  count_regions++;
	  ihh_count[il]=count_regions;
	  
	  
	  string new_mask="../MASK/new_mask_jacknife_Nside_"+to_string(Nside_new)+"_d1_"+to_string(count_regions)+".txt";
	  //     string new_mask="../MASK/new_mask_jacknife_all_regions_"+to_string(count_regions)+".txt";
	  
	  sal.open(new_mask);
	  
	  for(int ih=0;ih<this->n_pixels;++ih)
	    sal<<ih<<"\t"<<this->mask[this->i_lpix+ih*NCOLS_MASK]<<"\t"<<this->mask[this->i_bpix+ih*NCOLS_MASK]<<"\t"<<index[ih]<<endl;

	  for(int ih=0;ih<this->n_pixels;++ih)
	    ihh[ih].push_back(index[ih]);
	  
	  std::cout<<"Low res pixel "<<il<<CYAN<<":  "<<count_ori<<" originally masked pixels.    ";
	  std::cout<<count<<" High res pixels in the pixel "<<il <<RESET<<BLUE<<"  "<<Npx_masked<<" extra-masked pixels in region "<<count_regions<<RESET<<std::endl;
	  sal.close();
	  
	  
	  apix<<count_regions<<"\t"<<Npx_masked*this->area_pixel<<endl;
	  
	  
	  
	  // string new_fjn="../MASK/cat_low_res_mask_jacknife_Nside_"+to_string(Nside_new)+"_"+to_string(count_regions)+".txt"; // here I only write few columns, meant to be plot
	  // sal.open(new_fjn);
	  
	  // string new_fjn2="../MASK/cat_low_res_mask_jacknife_Nside_"+to_string(Nside_new)+"_"+to_string(count_regions)+"_fullinfo.txt"; // Here I write the same catalog in different low res pixels
	  // ofstream sal2; sal2.open(new_fjn2);



	  for(int i=0;i<this->NGAL;++i)
	    {
	      point.phi=prop[this->i_lgal+i*this->NCOLS]*fac;
	      point.theta=0.5*M_PI-this->prop[this->i_bgal+i*this->NCOLS]*fac;
	      if(this->gal_mask[i]==this->observed_pixels_in_mask)
		{
		  long pixel=mask_aux_low.ang2pix(point);
		  if(pixel==il){
            real_prec  alpha=this->prop[this->i_lgal+i*this->NCOLS]+180.0;real_prec  alpha2;

		    if(alpha>=360)alpha2=alpha-360.0; else alpha2=alpha;

		    // if(this->prop[this->i_zp+i*this->NCOLS]>=0)sal<<alpha2<<"  "<<prop[this->i_bgal+i*this->NCOLS]<<"  "<<prop[this->i_zp+i*this->NCOLS]<<endl;
		    // for(int ik=0;ik<this->NGAL;++ik)sal2<<prop[ik+i*this->NCOLS]<<"\t";
		    // sal2<<endl;
		    
		  }
		}
	    }
	  // sal.close();
	  // sal2.close();


	
	}


  


    }


  int cc=0;
  for(int i=0;i<n_pixels_new; ++i)
    {
      lres<<i<<" "<<ihh_new[i]<<" "<<ihh_count[i]<<endl;
      cc+=(1-ihh_new[i]);
    }
  cout<<RED<<cc<<endl;
  apix.close();





  /*

    for(int i=0;i<count_regions;++i){
    for(int j=i+1;j<count_regions;++j){
    if(i!=j){
    nmk++;
    string new_mask="../MASK/new_mask_jacknife_Nside_"+to_string(Nside_new)+"_d2_"+to_string(nmk)+".txt";
    sal.open(new_mask);
    for(int ih=0;ih<this->n_pixels;++ih)sal<<ih<<"\t"<<this->mask[ih][this->i_lpix]<<"\t"<<this->mask[ih][this->i_bpix]<<"\t"<<ihh[ih][i]*ihh[ih][j]<<endl;
    cout<<"Writting d2 mask "<<nmk<<" in file"<<new_mask<<endl;
    sal.close();
    }
    }
    }


    int nmk=0;

    for(int i=0;i<count_regions;++i){
    for(int j=i+1;j<count_regions;++j){
    if(i!j){
    for(int k=j+1;k<count_regions;++k){
    nmk++;
    string new_mask="../MASK/new_mask_jacknife_Nside_"+to_string(Nside_new)+"_d3_"+to_string(nmk)+".txt";
    sal.open(new_mask);
    for(int ih=0;ih<this->n_pixels;++ih)sal<<ih<<"\t"<<this->mask[ih][this->i_lpix]<<"\t"<<this->mask[ih][this->i_bpix]<<"\t"<<ihh[ih][i]*ihh[ih][j]*ihh[ih][k]<<endl;
    cout<<"Writting d3 mask "<<nmk<<" in file"<<new_mask<<endl;
    sal.close();
    }
    }
    }
    }

    cout<<nmk<<endl;
  */



  /*
    int nmk=0;

    for(int i=0;i<count_regions;++i){
    for(int j=i+1;j<count_regions;++j){
    for(int k=j+1;k<count_regions;++k){
    for(int l=k+1;l<count_regions;++l){
    for(int m=l+1;m<count_regions;++m){
    for(int n=m+1;n<count_regions;++n){
    nmk++;
    string new_mask="../MASK/new_mask_jacknife_Nside_"+to_string(Nside_new)+"_d6_"+to_string(nmk)+".txt";
    sal.open(new_mask);
    for(int ih=0;ih<this->n_pixels;++ih)sal<<ih<<"\t"<<this->mask[ih][this->i_lpix]<<"\t"<<this->mask[ih][this->i_bpix]<<"\t"<<ihh[ih][i]*ihh[ih][j]*ihh[ih][k]*ihh[ih][l]*ihh[ih][m]*ihh[ih][n]<<endl;
    cout<<"Writting d5 mask "<<nmk<<" in file"<<new_mask<<endl;
    sal.close();
    }
    }
    }
    }
    }
    }

  */



}
// ######################################################################
// ######################################################################

