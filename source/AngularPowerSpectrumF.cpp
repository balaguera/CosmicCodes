//##################################################################################
//##################################################################################
/** @file Bam.cpp
 *
 *  @brief Measurements of angular pwoer spectrum
 *  based on the BAM method.
 *  @author: Andrés Balaguera-Antolínez
 */
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

# include "../Headers/AngularPowerSpectrumF.h"
using namespace std;



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void message(){
  time_t rawtime;
  time (&rawtime);
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"ESTIMATES OF ANGULAR POWER SPECTRUM                   *"<<endl;
  std::cout<<"v1.2                                                  *"<<endl;
  std::cout<<"abalant@R3 2014-2016                                  *"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"*******************************************************"<<endl;
  std::cout<<"Starting time and date"<<std::endl;
  std::cout<<ctime (&rawtime)<<std::endl;
  std::cout<<"*******************************************************"<<endl;
}



void usage(string s){
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<RED<<"  Clgal a code to measure angular power spectrum of cosmological mass tracers"<<endl;
  cout<<RED<<"  Usage "<<s<<" [-option] [argument]"<<endl;
  cout<<RED<<"  Options:     -h for more information "<<endl;
  cout<<RED<<"               -a for information on the author"<<endl;
  cout<<RED<<"               -p parameter_file.ini to execture with input parameter file "<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;
  cout<<CYAN<<" ********************************************************************************"<<endl;}

void author(){
  cout<<CYAN<<"                                                                                 "<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" Copyright. 2017. Andres Balaguera Antolinez"<<endl;
  cout<<" Code developed for the power spectrum analysis of the 2MPZ sample"<<endl;
  cout<<" This code is public. If you want to use it for your research,     "<<endl;
  cout<<" contact the author at balaguera@iac.es"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<" ***************************************************************************************"<<endl;
  cout<<CYAN<<"                                                                                 "<<endl;
}



void message_time(time_t start_all){
    time_t end;
    time(&end);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
    healpix_real lapse=difftime(end,start_all);
    std::cout<<CYAN<<"Time elapsed: "<<lapse<<"  secs \r"<<RESET<<std::endl;
    time(&start_all);
    std::cout<<RED<<"=================================="<<std::endl;
    std::cout<<"=================================="<<RESET<<std::endl;
}




// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::get_Cl_gal()
{

  string file_Jlm_init_fixed=this->output_file_Jlm+".dat";
  this->output_file_Jlm=file_Jlm_init_fixed;
  
  if(this->file_type=="ascii")this->read_input_cats("g",this->input_file);


  if(!this->use_random_cat){
    // If using a mask instead of a random
    // catalog, use the nside (resolution) from this point on
    string mask_f=(this->sky=="full_sky"? this->input_file_mask_fs: this->input_file_mask);
    this->get_mask(mask_f);
    this->nside=floor(sqrt((healpix_real)this->n_pixels/12.));
    cout<<"Mask with nside = "<<this->nside<<endl;
  }

  // *********************************************************

  set_healpix_pars();
  get_pars_mask();
  set_vectors();     // Allocate memory for used vectors


  // *********************************************************


  if(this->n_z_bins!=-1){
    cout<<RED<<"SELECTION INFO "<<RESET<<endl;
    cout<<CYAN;
    if(selection=="fls"){
      cout<<"Flux limited sample selected. Bins in redshift"<<endl;
      if(this->define_z_bins=="delta")cout<<"zBins with constant width"<<endl;
      if(this->define_z_bins=="number")cout<<"zBins chosen with equal number of galaxies"<<endl;
    }
    
    // Get zmax
    this->zmax_cat=this->get_max('r',this->i_z);
    this->zmin_cat=this->get_min('d',this->i_z);
    cout<<RED<<"REDSHIFT INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim redshift found in catalogue = "<<this->zmax_cat<<endl;
    cout<<CYAN<<"Minimum redshift found in catalogue = "<<this->zmin_cat<<endl;
    if(!this->use_z_min_max_from_cat){
      cout<<CYAN<<"Maximim redshift selected in parameter file = "<<this->zmax<<endl;
      cout<<CYAN<<"Minimum redshift selected in parameter file = "<<this->zmin<<endl;
    }
    else{
      this->zmin=this->zmin_cat;
      this->zmax=this->zmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    // Get Magnitude limits
    this->MKmax_cat=this->get_max('d',this->i_M);
    this->MKmin_cat=this->get_min('d',this->i_M);
    cout<<RED<<"MAGNITUDE Mk INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim Mk found in catalogue = "<<this->MKmax_cat<<endl;
    cout<<CYAN<<"Minimum Mk found in catalogue = "<<this->MKmin_cat<<endl;
    if(!this->use_Mk_min_max_from_cat){
      cout<<CYAN<<"Maximim Mk selected in parameter file = "<<this->MKmax<<endl;
      cout<<CYAN<<"Minimum Mk selected in parameter file = "<<this->MKmin<<endl;
    }
    else{
      this->MKmin=this->MKmin_cat;
      this->MKmax=this->MKmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    
    // SOME WARNINGS:
    std::cout<<RED;
    if(this->MKmax > this->MKmax_cat)cout<<"Warning: Maximum Mk in parameter file ("<<this->MKmax<<") greater than maximum value found in the catalog ("<<this->MKmax_cat <<")"<<endl;
    if(this->MKmin < this->MKmin_cat)cout<<"Warning: Minimum Mk in parameter file ("<<this->MKmin<<") smaller than minimum value found in the catalog ("<<this->MKmin_cat <<")"<<endl;
    if(this->zmax > this->zmax_cat){
      //cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Setting value from catalog"<<endl;
      cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Doing nothing"<<endl;
      //   this->zmax= this->zmax_cat;
    }
    if(this->zmin < this->zmin_cat){
      cout<<"Warning: Minimum z in parameter file ("<<this->zmin<<") smaller than minimum value found in the catalog ("<<this->zmin_cat <<"). Doing nothing"<<endl;
    }
  }
  
  else{
    this->zmin=-1000;
    this->zmax=+1000;
  }
  
  
  // *******************************************************
  
  int NB;
  if(this->n_z_bins!=-1)NB = (this->n_z_bins);
  else NB=0;
  
  this->set_BLMZ(NB);
  
  // *******************************************************
  // Feed class with redshift bins
  
  if(this->n_z_bins!=-1)this->set_zbins();
  
  if(this->use_random_cat)this->read_input_cats("r", this->input_file_random);

  // ******************************************************************************************************************
  // *********************************************************
  // *********************************************************
  // *********************************************************
  // OPERATIONS IN HARMONIC SPACE
  // *********************************************************
  // *********************************************************
  // ******************************************************************************************************************

  
  // Estos los dejo definidos aca en directo...
  Alm<xcomplex <healpix_real> > Ilm(this->Lmax+1,this->Lmax+1);
  // Compute Jlm, Ilm. These objects are inizialized to the
  // full sky- case, Jlm=1, Ilm=0 (l>0)
  // such that if not computed in Map2ILM_Jlm, they remain full sky.
  // but if used, these have to be initialized to zero in that funcition
  // For full sky, Ilm os not set to zero by directly not used in the estimation

  this->set_ilm_jlm(); // So far here we only set Jlm
  if(this->use_random_cat)cout<<"Nothing to do here, read cat"<<endl;
  else if(!this->use_random_cat){  //IF WE HAVE A MASK instead of a random, get the Ilm, even IF DIRECT SUM ESTIMATION
    if(this->sky!="full_sky")this->Map2Ilm_Jlm(Ilm);
  }
  if(this->sky=="full_sky"){ // Intead if we have full sky, set Ilm to 1
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).real(1);//Just for full sky
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).imag(1);
  }

  cout<<"*******************************************************"<<endl;
  cout<<BLUE<<"Computing Alm in z bins"<<RESET<<endl;

  // *******************************************************
  int izf;
  int imf;
  int icut=0;

  if(this->n_z_bins!=-1){
    imf=0;
    izf=(this->n_z_bins==1? 0 : this->n_z_bins);
    this->set_mean_ngal_pix(this->n_z_bins);
  }
  else{
    this->set_mean_ngal_pix(1);
    izf=0;
  }
  // *******************************************************
  // *******************************************************

  this->set_Lbins();
  
  // START LOOP OVER THE REDSHIFT BINS ONLY IF FLS
  for(int IZ=0;IZ<=izf;++IZ){
      
    this->set_IZ(IZ);
      
    ///Get map from cat
    Healpix_Map<healpix_real>map(log2(this->nside), RING);
    this->get_map('d', this->file_type, map);
      
    // Define weights
    arr<healpix_real>weight_data(2*map.Nside());
    weight_data.fill(1.0);
      
    Healpix_Map<healpix_real>map_ran(log2(this->nside), RING);
    if(this->use_random_cat){
      this->get_map('r',this-> file_type, map_ran);
    }
      
    // ***************************************************************************************
    this->get_pars_cl_estimator();
    this->write_pars_cl_estimator();
    // ***************************************************************************************
    this->Mean_ngal_pix[IZ]=this->mean_number_galaxies_pix;
    this->Mean_ngal[IZ]=this->ngal_used/(this->area_survey);
    this->Shot_Noise[IZ]=this->shot_noise;
      
    // ***************************************************************************************
    Alm<xcomplex <healpix_real> > Alm(this->Lmax+1,this->Lmax+1);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).real(0);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).imag(0);
    weight_data.fill(1);
      
    cout<<RED<<"Computing Alm..."<<RESET<<endl;
    map2alm(map,Alm, weight_data,false); //HealPix
    cout<<"Done"<<endl;
      
    if(this->use_random_cat){
      arr<healpix_real>weight_random(4*map.Nside()); //Do not understand why 2*map.Nside()
      for(int i=0; i < 4*map.Nside(); ++i) weight_random[i]=1.0;
      map2alm(map_ran,Ilm, weight_random,false);
      this->Map2Jlm(map_ran);
    }
    cout<<"*******************************************************"<<endl;
      
    // ***********************************************************************************************
    // Fill the vectors Blm. These are the Alm in different redshift or Magnitude bins
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].real(Alm(il,im).real());
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].imag(Alm(il,im).imag());
      
  } //close loop over z, opened in sample=fls
    // ***********************************************************************************************
    
  cout<<BLUE;
  // Having allocated quantities in z bin (fls)
  // we now do a healpix_real loop (fls) over the zbins, or
  // we do only one loop (vls) (the second is currently running)
  // and compute things for the current Mbins
    
  cout<<"*******************************************************"<<endl;
  
  if(this->n_z_bins>1){
    cout<<"Computing Cross Cls"<<RESET<<endl;
    for(int IZ=1;IZ<=izf;IZ++){
      vector<healpix_real> ClzbinI(this->Lmax+1,0);
      this->Alm2Cl(IZ,IZ, Ilm, ClzbinI);//Auto Cl bin IZ
      for(int JZ=IZ;JZ<=izf;JZ++){
        vector<healpix_real> Cl(this->Lmax+1,0);
        vector<healpix_real> eCl(this->Lmax+1,0);
        vector<healpix_real> ClzbinJ(this->Lmax+1,0);
  
  // ***********************************************************************************************
  // We need to compute the autos within each loop in order to compute th
  // variance for the cross power spectrum (C_i+sn_j)*(C_j+sn_j).
  // NOTE THAT FIRST WE NEED TO COMPUTE THE AUTO. DO not change ordering here.
  // ***********************************************************************************************
  //Auto Cl in bin J:
  this->Alm2Cl(JZ,JZ,Ilm,ClzbinJ);
  //Cross:
  this->Alm2Cl(IZ,JZ,Ilm,Cl);
  //get variance of cross power:
  this->get_eCl(IZ,JZ,ClzbinI, ClzbinJ, eCl);
  
  // ***********************************************************************************************
  // Write raw estimates
  string rest;
  if(this->sampling=="H"){
    rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
  }
  else if(this->sampling=="DS"){
    rest="_"+this->sampling+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
    
  }
  string cfile=this->output_file_raw+rest;
        ////////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl,eCl);
  // ***********************************************************************************************
  // Cl in Bins
  this->Cl_bins(Cl, this->Clbin);
  // ***********************************************************************************************
  
  if(IZ==JZ){
    int ll=-1;
    do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[JZ]);
    if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<JZ<<". Select lmax = 100. "<<RESET<<std::endl;
    else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< JZ<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
  }
  // ***********************************************************************************************
  // Variance in Bins
  this->eCl_bins(eCl,this->eClbin);
  // Wl in Bins
  this->Cl_bins(this->Wl, this->Wlbin);
  // ***********************************************************************************************
  // Write Binned estimates
  string ofile=this->output_file+rest;
        ////////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
      }
    }
    
    {  // Get the Cl of the full sample
      vector<healpix_real> Cl(this->Lmax+1,0);
      vector<healpix_real> eCl(this->Lmax+1,0);
      this->Alm2Cl(0,0, Ilm,Cl);
      this->get_eCl(0,0,Cl,Cl,eCl);
      // ***********************************************************************************************
      string rest;
      rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(0)+"_"+to_string(0)+"_"+this->hemis+"_"+this->coord+".dat";
        
      string cfile=this->output_file_raw+rest;
      ////////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
      // ***********************************************************************************************
      this->Cl_bins(Cl,  this->Clbin);
      this->eCl_bins(eCl, this->eClbin);
      // ***********************************************************************************************
      this->Cl_bins(this->Wl, this->Wlbin);
  
  
      {
  int ll=-1;
  do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
  if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
  else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
      }
      // ***********************************************************************************************
      // ***********************************************************************************************
      string ofile=this->output_file+rest;
      ////////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
    }
      
  }

  else{   // Get the Cl when no redshift info is available
    vector<healpix_real> Cl(this->Lmax+1,0);
    vector<healpix_real> eCl(this->Lmax+1,0);
    this->Alm2Cl(0,0, Ilm,Cl);
    this->get_eCl(0,0,Cl,Cl,eCl);
    // ***********************************************************************************************
    string rest;
    rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_"+this->hemis+"_"+this->coord+".dat";
    
    string cfile=this->output_file_raw+rest;
    //////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
    // ***********************************************************************************************
    this->Cl_bins(Cl,  this->Clbin);
    this->eCl_bins(eCl, this->eClbin);
    // ***********************************************************************************************
    this->Cl_bins(this->Wl, this->Wlbin);
    {
      int ll=-1;
      do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
      if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
      else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
    }
    string ofile=this->output_file+rest;
    //////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
  }
  

  // Finally, get the mixig matrix, if desired
  if(this->compute_mixing_matrix)this->get_mixing_matrix();
    
  this->  prop.clear();
    
  return;



}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::get_my_power(){


  string file_Jlm_init_fixed=this->output_file_Jlm+".dat";
  this->output_file_Jlm=file_Jlm_init_fixed;
  
  if(this->file_type=="ascii")this->read_input_cats("g",this->input_file);


  if(!this->use_random_cat){
    // If using a mask instead of a random
    // catalog, use the nside (resolution) from this point on
    string mask_f=(this->sky=="full_sky"? this->input_file_mask_fs: this->input_file_mask);
    this->get_mask(mask_f);
    this->nside=floor(sqrt((healpix_real)this->n_pixels/12.));
    cout<<"Mask with nside = "<<this->nside<<endl;
  }

  // *********************************************************

  set_healpix_pars();
  get_pars_mask();
  set_vectors();     // Allocate memory for used vectors


  // *********************************************************


  if(this->n_z_bins!=-1){
    cout<<RED<<"SELECTION INFO "<<RESET<<endl;
    cout<<CYAN;
    if(selection=="fls"){
      cout<<"Flux limited sample selected. Bins in redshift"<<endl;
      if(this->define_z_bins=="delta")cout<<"zBins with constant width"<<endl;
      if(this->define_z_bins=="number")cout<<"zBins chosen with equal number of galaxies"<<endl;
    }
    
    // Get zmax
    this->zmax_cat=this->get_max('r',this->i_z);
    this->zmin_cat=this->get_min('d',this->i_z);
    cout<<RED<<"REDSHIFT INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim redshift found in catalogue = "<<this->zmax_cat<<endl;
    cout<<CYAN<<"Minimum redshift found in catalogue = "<<this->zmin_cat<<endl;
    if(!this->use_z_min_max_from_cat){
      cout<<CYAN<<"Maximim redshift selected in parameter file = "<<this->zmax<<endl;
      cout<<CYAN<<"Minimum redshift selected in parameter file = "<<this->zmin<<endl;
    }
    else{
      this->zmin=this->zmin_cat;
      this->zmax=this->zmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    // Get Magnitude limits
    this->MKmax_cat=this->get_max('d',this->i_M);
    this->MKmin_cat=this->get_min('d',this->i_M);
    cout<<RED<<"MAGNITUDE Mk INFO "<<RESET<<endl;
    cout<<CYAN<<"Maximim Mk found in catalogue = "<<this->MKmax_cat<<endl;
    cout<<CYAN<<"Minimum Mk found in catalogue = "<<this->MKmin_cat<<endl;
    if(!this->use_Mk_min_max_from_cat){
      cout<<CYAN<<"Maximim Mk selected in parameter file = "<<this->MKmax<<endl;
      cout<<CYAN<<"Minimum Mk selected in parameter file = "<<this->MKmin<<endl;
    }
    else{
      this->MKmin=this->MKmin_cat;
      this->MKmax=this->MKmax_cat;
    }
    std::cout<<"*******************************************************"<<endl;
    
    // SOME WARNINGS:
    std::cout<<RED;
    if(this->MKmax > this->MKmax_cat)cout<<"Warning: Maximum Mk in parameter file ("<<this->MKmax<<") greater than maximum value found in the catalog ("<<this->MKmax_cat <<")"<<endl;
    if(this->MKmin < this->MKmin_cat)cout<<"Warning: Minimum Mk in parameter file ("<<this->MKmin<<") smaller than minimum value found in the catalog ("<<this->MKmin_cat <<")"<<endl;
    if(this->zmax > this->zmax_cat){
      //cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Setting value from catalog"<<endl;
      cout<<"Warning: Maximum z in parameter file ("<<this->zmax<<") greater than maximum value found in the catalog ("<<this->zmax_cat <<"). Doing nothing"<<endl;
      //   this->zmax= this->zmax_cat;
    }
    if(this->zmin < this->zmin_cat){
      cout<<"Warning: Minimum z in parameter file ("<<this->zmin<<") smaller than minimum value found in the catalog ("<<this->zmin_cat <<"). Doing nothing"<<endl;
    }
  }
  
  else{
    this->zmin=-1000;
    this->zmax=+1000;
  }
  
  
  // *******************************************************
  
  int NB;
  if(this->n_z_bins!=-1)NB = (this->n_z_bins);
  else NB=0;
  
  this->set_BLMZ(NB);
  
  // *******************************************************
  // Feed class with redshift bins
  
  if(this->n_z_bins!=-1)this->set_zbins();
  
  if(this->use_random_cat)this->read_input_cats("r", this->input_file_random);

  // ******************************************************************************************************************
  // *********************************************************
  // *********************************************************
  // *********************************************************
  // OPERATIONS IN HARMONIC SPACE
  // *********************************************************
  // *********************************************************
  // ******************************************************************************************************************

  
  // Estos los dejo definidos aca en directo...
  Alm<xcomplex <healpix_real> > Ilm(this->Lmax+1,this->Lmax+1);
  // Compute Jlm, Ilm. These objects are inizialized to the
  // full sky- case, Jlm=1, Ilm=0 (l>0)
  // such that if not computed in Map2ILM_Jlm, they remain full sky.
  // but if used, these have to be initialized to zero in that funcition
  // For full sky, Ilm os not set to zero by directly not used in the estimation

  this->set_ilm_jlm(); // So far here we only set Jlm
  if(this->use_random_cat)cout<<"Nothing to do here, read cat"<<endl;
  else if(!this->use_random_cat){  //IF WE HAVE A MASK instead of a random, get the Ilm, even IF DIRECT SUM ESTIMATION
    if(this->sky!="full_sky")this->Map2Ilm_Jlm(Ilm);
  }
  if(this->sky=="full_sky"){ // Intead if we have full sky, set Ilm to 1
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).real(1);//Just for full sky
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Ilm(i,im).imag(1);
  }

  cout<<"*******************************************************"<<endl;
  cout<<BLUE<<"Computing Alm in z bins"<<RESET<<endl;

  // *******************************************************
  int izf;
  int imf;
  int icut=0;

  if(this->n_z_bins!=-1){
    imf=0;
    izf=(this->n_z_bins==1? 0 : this->n_z_bins);
    this->set_mean_ngal_pix(this->n_z_bins);
  }
  else{
    this->set_mean_ngal_pix(1);
    izf=0;
  }
  // *******************************************************
  // *******************************************************

  this->set_Lbins();
  
  // START LOOP OVER THE REDSHIFT BINS ONLY IF FLS
  for(int IZ=0;IZ<=izf;++IZ){
      
    this->set_IZ(IZ);
      
    ///Get map from cat
    Healpix_Map<healpix_real>map(log2(this->nside), RING);
    this->get_map('d', this->file_type, map);
      
    // Define weights
    arr<healpix_real>weight_data(2*map.Nside());
    weight_data.fill(1.0);
      
    Healpix_Map<healpix_real>map_ran(log2(this->nside), RING);
    if(this->use_random_cat){
      this->get_map('r',this-> file_type, map_ran);
    }
      
    // ***************************************************************************************
    this->get_pars_cl_estimator();
    this->write_pars_cl_estimator();
    // ***************************************************************************************
    this->Mean_ngal_pix[IZ]=this->mean_number_galaxies_pix;
    this->Mean_ngal[IZ]=this->ngal_used/(this->area_survey);
    this->Shot_Noise[IZ]=this->shot_noise;
      
    // ***************************************************************************************
    Alm<xcomplex <healpix_real> > Alm(this->Lmax+1,this->Lmax+1);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).real(0);
    for(int i=0;i<this->Lmax+1;i++)for(int im=0;im<this->Lmax+1;im++)Alm(i,im).imag(0);
    weight_data.fill(1);
      
    cout<<RED<<"Computing Alm..."<<RESET<<endl;
    map2alm(map,Alm, weight_data,false); //HealPix
    cout<<"Done"<<endl;
      
    if(this->use_random_cat){
      arr<healpix_real>weight_random(4*map.Nside()); //Do not understand why 2*map.Nside()
      for(int i=0; i < 4*map.Nside(); ++i) weight_random[i]=1.0;
      map2alm(map_ran,Ilm, weight_random,false);
      this->Map2Jlm(map_ran);
    }
    cout<<"*******************************************************"<<endl;
      
    // ***********************************************************************************************
    // Fill the vectors Blm. These are the Alm in different redshift or Magnitude bins
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].real(Alm(il,im).real());
    for(int il=0;il<this->Lmax+1;il++)for(int im=0;im<this->Lmax+1;im++)this->Blm[il][im][IZ].imag(Alm(il,im).imag());
      
  } //close loop over z, opened in sample=fls
    // ***********************************************************************************************
    
  cout<<BLUE;
  // Having allocated quantities in z bin (fls)
  // we now do a healpix_real loop (fls) over the zbins, or
  // we do only one loop (vls) (the second is currently running)
  // and compute things for the current Mbins
    
  cout<<"*******************************************************"<<endl;
  
  if(this->n_z_bins>1){
    cout<<"Computing Cross Cls"<<RESET<<endl;
    for(int IZ=1;IZ<=izf;IZ++){
      vector<healpix_real> ClzbinI(this->Lmax+1,0);
      this->Alm2Cl(IZ,IZ, Ilm, ClzbinI);//Auto Cl bin IZ
      for(int JZ=IZ;JZ<=izf;JZ++){
        vector<healpix_real> Cl(this->Lmax+1,0);
        vector<healpix_real> eCl(this->Lmax+1,0);
        vector<healpix_real> ClzbinJ(this->Lmax+1,0);
	
	// ***********************************************************************************************
	// We need to compute the autos within each loop in order to compute th
	// variance for the cross power spectrum (C_i+sn_j)*(C_j+sn_j).
	// NOTE THAT FIRST WE NEED TO COMPUTE THE AUTO. DO not change ordering here.
	// ***********************************************************************************************
	//Auto Cl in bin J:
	this->Alm2Cl(JZ,JZ,Ilm,ClzbinJ);
	//Cross:
	this->Alm2Cl(IZ,JZ,Ilm,Cl);
	//get variance of cross power:
	this->get_eCl(IZ,JZ,ClzbinI, ClzbinJ, eCl);
	
	// ***********************************************************************************************
	// Write raw estimates
	string rest;
	if(this->sampling=="H"){
	  rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
	}
	else if(this->sampling=="DS"){
	  rest="_"+this->sampling+"_zbins_"+to_string(IZ)+"_"+to_string(JZ)+"_"+this->hemis+"_"+this->coord+".dat";
	  
	}
	string cfile=this->output_file_raw+rest;
        ////////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl,eCl);
	// ***********************************************************************************************
	// Cl in Bins
	this->Cl_bins(Cl, this->Clbin);
	// ***********************************************************************************************
	
	if(IZ==JZ){
	  int ll=-1;
	  do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[JZ]);
	  if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<JZ<<". Select lmax = 100. "<<RESET<<std::endl;
	  else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< JZ<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
	}
	// ***********************************************************************************************
	// Variance in Bins
	this->eCl_bins(eCl,this->eClbin);
	// Wl in Bins
	this->Cl_bins(this->Wl, this->Wlbin);
	// ***********************************************************************************************
	// Write Binned estimates
	string ofile=this->output_file+rest;
        ////////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
      }
    }
    
    {  // Get the Cl of the full sample
      vector<healpix_real> Cl(this->Lmax+1,0);
      vector<healpix_real> eCl(this->Lmax+1,0);
      this->Alm2Cl(0,0, Ilm,Cl);
      this->get_eCl(0,0,Cl,Cl,eCl);
      // ***********************************************************************************************
      string rest;
      rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_zbins_"+to_string(0)+"_"+to_string(0)+"_"+this->hemis+"_"+this->coord+".dat";
      	
      string cfile=this->output_file_raw+rest;
      ////////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
      // ***********************************************************************************************
      this->Cl_bins(Cl,  this->Clbin);
      this->eCl_bins(eCl, this->eClbin);
      // ***********************************************************************************************
      this->Cl_bins(this->Wl, this->Wlbin);
	
	
      {
	int ll=-1;
	do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
	if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
	else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
      }
      // ***********************************************************************************************
      // ***********************************************************************************************
      string ofile=this->output_file+rest;
      ////////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
    }
      
  }

  else{   // Get the Cl when no redshift info is available
    vector<healpix_real> Cl(this->Lmax+1,0);
    vector<healpix_real> eCl(this->Lmax+1,0);
    this->Alm2Cl(0,0, Ilm,Cl);
    this->get_eCl(0,0,Cl,Cl,eCl);
    // ***********************************************************************************************
    string rest;
    rest="_"+this->sampling+"_Nside_"+to_string(this->nside)+"_"+this->hemis+"_"+this->coord+".dat";
    
    string cfile=this->output_file_raw+rest;
    //////////////////////////////////////////////this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
    // ***********************************************************************************************
    this->Cl_bins(Cl,  this->Clbin);
    this->eCl_bins(eCl, this->eClbin);
    // ***********************************************************************************************
    this->Cl_bins(this->Wl, this->Wlbin);
    {
      int ll=-1;
      do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
      if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
      else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
    }
    string ofile=this->output_file+rest;
    //////////////////////////////////////////////this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
  }
  

  // Finally, get the mixig matrix, if desired
  if(this->compute_mixing_matrix)this->get_mixing_matrix();
    
  this->  prop.clear();
    
  return;
    
}

// *****************************************************************************************************************
// *****************************************************************************************************************
// *****************************************************************************************************************
// *****************************************************************************************************************
// *****************************************************************************************************************
// *****************************************************************************************************************



void AngularPowerF::read_pars(string &file){

  // Intializing parameters:
  nzbins=0;
  code="cross";
  statistics = "Cl";
  name_catalog = "Galaxy_Catalog";
  name_output_dir ="../Output/";
  name_input_dir = "../Input/";
  input_file = "catalog.dat";
  file_type="ascii";
  file_type_mask="ascii";
  use_random_cat = false;
  sky="masked";
  coord="galactic";
  generate_fits_files=false;
  n_z_bins=-1;
  define_z_bins= "delta";
  bin_type="linear";
  sampling="H";
// ***********************************************************************************************

    // Reading from input file
  
  ifstream fin_parameters (file.c_str());
  if (!fin_parameters) { cerr <<"Error in opening the parameters file "<<file<<"!"<<endl; exit(1); }

  string line_in_file;
  string par_name;
  string equality;
  string par_value;

  
  while (getline(fin_parameters,line_in_file)) {
    // ignore lines starting with hashtag
    if(line_in_file[0] != '#' && line_in_file.empty()==0) {
      stringstream line_string (line_in_file);      // read parameter name
      line_string << line_in_file;
      line_string >> par_name;
      line_string >> equality;      // check that second word is "="
      if (equality != "=") {
	cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
	cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      line_string >> par_value;      // read parameter value
      if (par_value.empty()) {
	cout << "Value of " << par_name << " not specified in " <<file << endl;
	cout << "Assuming a default value for " << par_name << endl;
	continue;
      }
      if (par_name == "Lmax")Lmax = atoi(par_value.c_str());
      else if (par_name == "Lmin")Lmin = atoi(par_value.c_str());
      else if (par_name == "N_L_bins")N_L_bins = atoi(par_value.c_str());
      else if (par_name == "nside")nside = atoi(par_value.c_str());

      else if (par_name == "file_type")file_type = par_value;
      else if (par_name == "file_type_mask")file_type_mask = par_value;
      else if (par_name == "name_input_dir")name_input_dir = par_value;
      else if (par_name == "input_file")input_file = par_value;
      else if (par_name == "bin_type")bin_type = par_value;
      else if (par_name == "name_catalog")name_catalog = par_value;
      else if (par_name == "name_output_dir") name_output_dir = par_value;

      else if (par_name == "use_random_cat"){
        if(par_value=="yes")use_random_cat=true;
	else use_random_cat=false;
      }
      else if (par_name == "generate_fits_files"){
        if(par_value=="yes")generate_fits_files=true;
        else generate_fits_files=false;
      }

      else if (par_name == "compute_mixing_matrix"){
        if(par_value=="yes")compute_mixing_matrix=true;
	else compute_mixing_matrix=false;
      }
      else if (par_name == "mixing_matrix_exact"){
        if(par_value=="yes") mixing_matrix_exact=true;
	else mixing_matrix_exact=false;
      }
      else if (par_name == "use_z_min_max_from_cat"){
        if(par_value=="yes") use_z_min_max_from_cat=true;
	else use_z_min_max_from_cat=false;
      }
      else if (par_name == "use_weight"){
        if(par_value=="yes") use_weight=true;
        else use_weight=false;
      }
      else if (par_name == "compute_jlm"){
        if(par_value=="yes")compute_jlm=true;
        else compute_jlm=false;
      }
      else if (par_name == "compute_property_weighted_cl"){
        if(par_value=="yes")this->compute_property_weighted_cl=true;
        else this->compute_property_weighted_cl=false;
      }
      else if (par_name == "define_z_bins")this->define_z_bins= par_value;
      else if (par_name == "sky")this->sky= par_value;
      else if (par_name == "type_of_P_estimator")this->type_of_P_estimator= par_value;
      else if (par_name == "hemis")this->hemis = par_value;
      else if (par_name == "redshift")this->redshift = par_value;
      else if (par_name == "shot_noise_correction")this->shot_noise_correction = par_value;
      else if (par_name == "input_file_random")this->input_file_random = par_value;
      else if (par_name == "input_file_mask") this->input_file_mask = par_value;
      else if (par_name == "input_file_mask_fs") this->input_file_mask_fs = par_value;
      else if (par_name == "input_file_mask_north_equatorial") this->input_file_mask_north_equatorial = par_value;
      else if (par_name == "input_file_mask_south_equatorial") this->input_file_mask_south_equatorial = par_value;
      else if (par_name == "coord")this->coord = par_value.c_str();
      else if (par_name == "output_file_raw") this->output_file_raw = par_value;
      else if (par_name == "output_file")this->output_file = par_value;
      else if (par_name == "output_file_Jlm") this->output_file_Jlm = par_value;
      else if (par_name == "output_file_window") this->output_file_window = par_value;
      else if (par_name == "fits_map")this->fits_map = par_value;
      else if (par_name == "i_alpha") this->i_alpha= atoi(par_value.c_str());
      else if (par_name == "i_delta") this->i_delta= atoi(par_value.c_str());
      else if (par_name == "i_z") this->i_z= atoi(par_value.c_str());
      else if (par_name == "i_M") this->i_M= atoi(par_value.c_str());
      else if (par_name == "i_alpha") i_alpha_ran= atoi(par_value.c_str());
      else if (par_name == "i_w") i_w= atoi(par_value.c_str());
      else if (par_name == "i_w_ran") i_w_ran= atoi(par_value.c_str());
      else if (par_name == "i_delta_ran") this->i_delta_ran= atoi(par_value.c_str());
      else if (par_name == "i_z_ran") this->i_z_ran= atoi(par_value.c_str());
      else if (par_name == "i_M_ran") this->i_M_ran= atoi(par_value.c_str());
      else if (par_name == "i_mask_alpha") this->i_mask_alpha= atoi(par_value.c_str());
      else if (par_name == "i_mask_delta") this->i_mask_delta= atoi(par_value.c_str());
      else if (par_name == "i_mask_pixel") this->i_mask_pixel= atoi(par_value.c_str());
      else if (par_name == "i_mask_flag") this->i_mask_flag= atoi(par_value.c_str());
      else if (par_name == "zmin")this->zmin= (healpix_real)(atof(par_value.c_str()));
      else if (par_name == "zmax")this->zmax= (healpix_real)(atof(par_value.c_str()));
      else if (par_name == "zmin_bin")this->zmin_bin= (healpix_real)(atof(par_value.c_str()));
      else if (par_name == "zmax_bin")this->zmax_bin= (healpix_real)(atof(par_value.c_str()));



    }
  }

  // Ammend name of files:
  input_file=name_input_dir+input_file;
  output_file=name_output_dir+statistics+"_"+type_of_P_estimator+"_"+name_catalog+"_zbtype_"+this->define_z_bins+"_"+output_file;
  output_file_raw=name_output_dir+statistics+"_"+type_of_P_estimator+"_"+name_catalog+"_zbtype_"+this->define_z_bins+"_"+output_file_raw;
  output_file_window=name_output_dir+statistics+"_"+name_catalog+"_"+output_file_window+"_"+this->hemis+"_"+this->coord;
;
  fits_map=name_output_dir+fits_map;
  if(this->compute_property_weighted_cl)output_file=output_file+"_marked";
  if(this->compute_property_weighted_cl)output_file_raw=output_file_raw+"_marked";
  output_file_Jlm=name_output_dir+output_file_Jlm+"_2MPZ_mask_nside_"+to_string(this->nside)+"_Lmax_"+to_string(this->Lmax)+"_"+this->hemis+"_"+this->coord;


  cout<<RESET<<endl;
 }

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *****************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::set_BLMZ(int nb){
  // Resize Blm. I have to do it in an independend method
  this->Blm.resize(this->Lmax+1);
  for(int i=0; i<Blm.size();++i)this->Blm[i].resize(this->Lmax+1);
  for(int i=0; i<this->Blm.size();++i)for(int j=0;j<this->Blm[i].size();++j)this->Blm[i][j].resize(nb+1);

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::set_vectors(){


    this->z_min.resize(this->n_z_bins+1,0);
    this->z_max.resize(this->n_z_bins+1,0);
    this->zminmax.resize(this->n_z_bins+1);
    for(int i=0;i<=this->n_z_bins;i++)this->zminmax[i].resize(2,0);
    this->Wl.resize(this->Lmax+1);
    this->lvec.resize(this->Lmax+1);
    this->lbin.resize(this->N_L_bins,0);
    this->lbin_min.resize(this->N_L_bins,0);
    this->lbin_max.resize(this->N_L_bins,0);
    this->Clbin.resize(this->N_L_bins,0);
    this->eClbin.resize(this->N_L_bins,0);
    this->nmodes.resize(this->N_L_bins,0);
    this->Wlbin.resize(this->N_L_bins,0);


    this->R.resize(this->Lmax+1);
    for(int i=0;i<this->R.size();i++)this->R[i].resize(this->Lmax+1,0);

    // THIS IS ALSO USED WHTN WE READ THE MIXING MATRIX!!!

    if(this->compute_mixing_matrix){
      R.resize(this->Lmax+1);
      for(int i=0;i<this->R.size();i++)R[i].resize(this->Lmax+1,0);
      Rll_bin.resize(this->N_L_bins);
      for(int i=0;i<this->Rll_bin.size();i++)Rll_bin[i].resize(this->Lmax+1,0);
    }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::set_Lbins(){

   if(this->bin_type=="linear"){
    healpix_real deltal=((healpix_real)(this->Lmax-this->Lmin))/((healpix_real)lbin.size());
    //here and just here I define N_L_bins as a healpix_real instead of an integer
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
      healpix_real deltal=log10(Lmax/1.0)/((healpix_real)lbin.size());
      for(int i=0;i<lbin.size();i++)lbin[i]=pow(10,log10(1)+(i+0.5)*deltal);
      for(int i=0;i<lbin.size();i++)lbin_min[i]=pow(10,log10(1)+i*deltal);
      for(int i=0;i<lbin.size();i++)lbin_max[i]=pow(10,log10(1)+(i+1)*deltal);
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::set_IZ(int IZZ){IZ=IZZ;}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::set_zbins(){

  this->z_min[0]=this->zmin;
  this->z_max[0]=this->zmax;
  healpix_real Dz=(this->zmax-this->zmin)/(healpix_real(this->n_z_bins));


   if(this->define_z_bins=="number")get_zbins_same_ngal(this->n_z_bins,this->i_z,this->zmin,this->zmax,this->zminmax);

  for(int i=1;i<=this->n_z_bins;i++){
    z_min[i]= (this->define_z_bins=="number"? this->zminmax[i][0]:this->zmin+(i-1)*Dz);
    z_max[i]= (this->define_z_bins=="number"? this->zminmax[i][1]:this->zmin+(i)*Dz);
  }
  cout<<RED<<"REDSHIFT BINS: "<<endl;
  if(this->n_z_bins>1)for(int i=0;i<=this->n_z_bins;i++)cout<<i<<"  ("<<z_min[i]<<"-"<<z_max[i]<<")"<<endl;
  else cout<<"0  ("<<z_min[0]<<"-"<<z_max[0]<<")"<<endl;
  cout<<RESET<<endl;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::set_healpix_pars(){
  n_pixels=12*this->nside*this->nside;
  Nrings=4*this->nside-1;
}


// *******************************************************************************************************************************************************
int AngularPowerF::npix_ring(int nside_h, int ir){
  // Number of pixels contained in a ring labled with ir. with 0<ir< 4nside-1 rings
  int Npixels_ring;
  if(ir<nside_h-1)Npixels_ring=4*(ir+1);
  if(ir>=nside_h-1 && ir<3*nside_h)Npixels_ring=4*nside_h;
  if(ir>=3*nside_h)Npixels_ring=4*(4*nside_h-ir-1);
  return Npixels_ring;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::get_eCl(int iz, int jz, vector<healpix_real>cli,vector<healpix_real>clj, vector<healpix_real>&ecl){
  /// WARNING: THE CROSS POWER SPECTRUM NEEDS HAS UNDESTIMATED GAUSSIAN ERROR BARS, CHECN WHITE, SUNG, PERCIVAL

  for(int l=this->Lmin;l<=this->Lmax;++l){
    ecl[l]=sqrt(  (2./(this->sky_fraction*(2*l+1)))*(     (cli[l]+this->Shot_Noise[iz])*(clj[l]+this->Shot_Noise[jz])) );
  }

}


// *******************************************************************************************************************************************************
void AngularPowerF::get_mask(string input_file_mask){


  if(this->file_type_mask=="ascii"){

    vector<real_prec>mask2;
    this->n_columns_mask=Fmd.read_file(input_file_mask, mask2);

    // set the n_pixels of the mask
    this->n_pixels=mask2.size()/this->n_columns_mask;
    this->n_total_pixels=this->n_pixels;
    
  // set the nside of the mask
    this->nside=sqrt(this->n_pixels/12);

  // nr=nrings, have to compute it here for nside 
  // is going to be the size of the mask and has not been yet computed
    int nr=4*this->nside-1;
  
    this->theta_new.resize(nr);
    fill(this->theta_new.begin(), this->theta_new.end(), 0);

    this->phi_new.resize(mask2.size()/this->n_columns_mask);
    fill(phi_new.begin(), this->phi_new.end(), 0);


    this->theta_new_pix.resize(n_pixels);
    fill(theta_new_pix.begin(), this->theta_new_pix.end(), 0);


    this->pixmask.resize(mask2.size()/this->n_columns_mask);
    fill(this->pixmask.begin(), this->pixmask.end(), 0);
  
    n_observed_pixels=0;
  
  // Create an auxiliary Healpix map to
  // copy the mask and find the angles
    Healpix_Map<healpix_real>mask_aux(log2(this->nside), RING);

    
    for(int i=0;i<n_pixels;i++)mask_aux[i]=mask2[this->i_mask_flag+this->n_columns_mask*i];
    
    pointing point_rev;
    int ip=0;

    for(int ir=0;ir<nr;ir++){
      int Npixels_ring=npix_ring(nside, ir);
      for(int ipr=0;ipr<Npixels_ring;++ipr){
        point_rev=mask_aux.pix2ang(ip);
        this->theta_new[ir]=point_rev.theta;
        this->phi_new[ip]=point_rev.phi;
        this->theta_new_pix[ip]=point_rev.theta;
        int paux;

	int pmask=mask2[this->i_mask_flag+this->n_columns_mask*ip];
	if(this->hemis=="all"){
          paux=(pmask== 0 ? 0 : 1);  //what comes from a real mask
        }
        else{
          if(this->hemis=="north"){
            if(theta_new[ir]>=0.5*M_PI)paux=0; // mask south
	    //            else paux=mask[ip][this->i_mask_flag];
	    else paux=pmask;
          }
          else if(this->hemis=="south"){
            if(theta_new[ir]<0.5*M_PI)paux=0;  //mask north
	    else paux=pmask;
	    //            else paux = mask[ip][this->i_mask_flag];
          }
        }
        this->pixmask[ip]=paux;
        if(paux==1)++n_observed_pixels;
        ++ip;
      }
    }

    mask2.clear();
  }

  else if(this->file_type_mask=="fits"){
      /*
    cout<<CYAN<<"Reading mask from "<<endl;
    cout<<input_file_mask<<endl;
    cout<<RESET<<endl;
    Healpix_Map<healpix_real>mask_aux(log2(this->nside), RING);
    read_Healpix_map_from_fits(this->input_file_mask,mask_aux, 1);
    this->nside=mask_aux.Nside();
    gal.mask.resize(12*nside*nside);
    for(int i=0;i<mask.size();++i)mask[i].resize(3,0);
    for(int i=0;i<n_pixels;i++)gal.mask[i][this->i_mask_flag]=mask_aux[i];
    set_healpix_pars();
*/
  }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::get_map(char c, string file_type, Healpix_Map<healpix_real>&map){

  if(file_type=="ascii"){
    if(c=='d')cout<<CYAN<<"Building map for real catalogue..."<<endl;
    if(c=='r')cout<<"Building map for random catalogue..."<<RESET<<endl;

    if(this->n_z_bins!=-1)Cat2Map(c,map);
    else Cat2Map_noz(c,map);

    string map_fits;
    map_fits=fits_map+"_zbin_"+to_string(this->IZ)+".fits";
    if(this->generate_fits_files){
      cout<<CYAN<<"Writing map in fits format in file "<<map_fits<<RESET<<endl;
      write_Healpix_map_to_fits(map_fits, map, (PDT)1);
    }
    if(c=='d')ngal=this->  prop.size()/this->n_columns_gal;
    if(c=='r')nran=this->  prop_r.size()/this->n_columns_ran;
    cout<<CYAN<<"Done.   "<<ngal<<"  "<<this->n_columns_gal<<RESET<<endl;
  }
  else{
    if(file_type=="fits"){
      cout<<CYAN<<"Reading map from "<<endl;
      cout<<input_file<<endl;
      cout<<RESET<<endl;
  //    read_Healpix_map_from_fits(input_file,map, 1);
      this->nside=map.Nside();
      set_healpix_pars();
    }
  }

  if(c=='d'){
    rms_ngal=map.rms();
    this->mean_number_galaxies_pix=(healpix_real)ngal_used/((healpix_real)n_observed_pixels);
  }
  if(c=='r'){
    rms_nran=map.rms();
    this->mean_number_randoms_pix=(healpix_real)nran_used/((healpix_real)n_observed_pixels);
  }

  
  if(file_type=="ascii"){
    Healpix_Map<healpix_real>map_delta(log2(this->nside), RING);
    if(this->n_z_bins!=-1)Cat2Map_delta(c,map_delta);
    
    /*
      string map_fits_delta=fits_map+"_delta_zbin_"+to_string(this->IZ)+".fits";
      cout<<CYAN<<"Writing delta_map in fits format in file "<<map_fits_delta<<RESET<<endl;
      write_Healpix_map_to_fits(map_fits_delta, map_delta, (PDT)1);
      cout<<CYAN<<"Done."<<RESET<<endl;
    */
  
 }
 

}


// *******************************************************************************************************************************************************
void AngularPowerF::Cat2Map(char c, Healpix_Map<healpix_real>&map){
  
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int iz=(c=='r'? this->i_z_ran:this->i_z);
  int iw=(c=='r'? this->i_w_ran:this->i_w);
  int iM=(c=='r'? this->i_M_ran:this->i_M);
  int ncols = (c=='r'? this->n_columns_ran:this->n_columns_gal);
  if(c=='r')nran_used=0;
  if(c=='d')ngal_used=0;
  vector<real_prec> prop_a;
  prop_a = (c=='d'? this->  prop:this->  prop_r) ;
  vector<real_prec>w_property(prop_a.size()/ncols,0);
   int nau=0;  
  // Here we first allocate the normalized property, in case we want it
  if(this->compute_property_weighted_cl){
    healpix_real mean_p=0;
    for(int i=0;i<prop_a.size();++i){
      long pixel;
      point.phi=prop_a[ira+ncols*i]*fac;
      point.theta=0.5*M_PI-fac*prop_a[idec+ncols*i];
      healpix_real weight1=(this->use_weight? prop_a[iw+ncols*i]:1.0);
      healpix_real weight2=(this->compute_property_weighted_cl? w_property[i]:1.0);
      pixel=map.ang2pix(point);
      if(this->pixmask[pixel]==1){
	if(prop_a[iz+ncols*i]<this->z_max[this->IZ] && prop_a[iz+ncols*i]>=this->z_min[this->IZ]){     // Now do the map in the redshift bin:
	  ++nau;
	  if(this->property_i_M_type=="original"){
	    w_property[i]=prop_a[iM+ncols*i];
	    mean_p+=prop_a[iM+ncols*i];
	  }
	  else if (this->property_i_M_type=="log"){
	    w_property[i]=log10(prop_a[iM+ncols*i]);
	    mean_p+=w_property[i];
	  }
	  else if (this->property_i_M_type=="minus_power_ten"){
	    w_property[i]=pow(10,-2.*prop_a[iM+ncols*i]/5.);
	    mean_p+=w_property[i];
	  }
	}
      }
    }
    mean_p/=(healpix_real)nau;
    for(int i=0;i<w_property.size();++i)w_property[i]/=mean_p;
    cout<<CYAN<<"Mean value of property "<<mean_p<<RESET<<endl;  
  }
  else {
    for(int i=0;i<w_property.size();++i)w_property[i]=1.0;
  }
  vector<healpix_real>map_aux(map.Npix(),0);
  nau=0;
  map.fill(0.0);
  for(int i=0;i<prop_a.size()/ncols;i++){
    // These lines
    // allows us to fill the vector map_aux
    // which will tell us the number of empty pixels
    // for all redshifts, used to compute 
    // the total surveyed area.
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;

    healpix_real weight1=(this->use_weight? prop_a[iw+ncols*i]:1.0);
    healpix_real weight2=(this->compute_property_weighted_cl? w_property[i]:1.0);

    pixel=map.ang2pix(point);
    map_aux[pixel]++;

    if(pixmask[pixel]==1){ 
      //Just to be sure, we go throgh pixels allowed by the mask
      //If the catalogue is already masked with the same mask, fine
      //but there might be cases in which we use a non masked catalogue.
      if(prop_a[iz+ncols*i]< z_max[this->IZ] && prop_a[iz+ncols*i]>=z_min[this->IZ]){     // Now do the map in the redshift bin:
	++nau;
	map[pixel]+=weight1*weight2;
      }
    }
  }
  prop_a.clear();
  if(c=='d')this->ngal_used=nau;
  if(c=='r')this->nran_used=nau;
  

  // Build here theta_new and phi_new
  // from the map of the randoms
  if(this->use_random_cat){
    if(c=='r'){
      // Create an auxiliary Healpix map to
      // copy the mask and find the angles

      int nr=4*nside-1; 
      
      theta_new.resize(nr);
      fill(theta_new.begin(), theta_new.end(), 0);
      
      phi_new.resize(n_pixels);
      fill(phi_new.begin(), phi_new.end(), 0);
      
      Healpix_Map<healpix_real>map_aux(log2(nside), RING);

      for(int i=0;i<n_pixels;i++)map_aux[i]=map[i];
      pointing point_rev;
      int ip=0;
      int no_pix=0;
      for(int ir=0;ir<nr;ir++)
      {
      	int Npixels_ring=npix_ring(nside, ir);
	      for(int ipr=0;ipr<Npixels_ring;ipr++)
         {
          point_rev=map_aux.pix2ang(ip);
	        theta_new[ir]=point_rev.theta;
      	  phi_new[ip]=point_rev.phi;
           if(map_aux[ip]==0)
            ++no_pix;
      	  ++ip;
	       }
      }
      n_observed_pixels=map.Npix()-no_pix;
      n_total_pixels=map.Npix();
    }
  }
  cout<<"Cat2Map done"<<endl;

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::Cat2Map_noz(char c, Healpix_Map<healpix_real>&map){
  
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int ncols = (c=='r'? this->n_columns_ran:this->n_columns_gal);
  if(c=='r')nran_used=0;
  if(c=='d')ngal_used=0;
  vector<real_prec> prop_a;
  prop_a = (c=='d'? this->  prop:this->  prop_r) ;
  
  int nau=0;
  

  vector<healpix_real>map_aux(map.Npix(),0);
  nau=0;
  map.fill(0.0);
  cout<<prop_a.size()/ncols<<endl;
  for(int i=0;i<prop_a.size()/ncols;i++){
    // These lines
    // allows us to fill the vector map_aux
    // which will tell us the number of empty pixels
    // for all redshifts, used to compute 
    // the total surveyed area.
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;

    pixel=map.ang2pix(point);
    map_aux[pixel]++;

    if(pixmask[pixel]==1){ 
      //Just to be sure, we go throgh pixels allowed by the mask
      //If the catalogue is already masked with the same mask, fine
      //but there might be cases in which we use a non masked catalogue.
      ++nau;
      map[pixel]++;
    } 
  }
  prop_a.clear();
  if(c=='d')this->ngal_used=nau;
  if(c=='r')this->nran_used=nau;
  
  
  // Build here theta_new and phi_new
  // from the map of the randoms
  if(this->use_random_cat){
    if(c=='r'){
      // Create an auxiliary Healpix map to
      // copy the mask and find the angles

      int nr=4*nside-1; 
      
      theta_new.resize(nr);
      fill(theta_new.begin(), theta_new.end(), 0);
      
      phi_new.resize(n_pixels);
      fill(phi_new.begin(), phi_new.end(), 0);
      
      Healpix_Map<healpix_real>map_aux(log2(nside), RING);

      for(int i=0;i<n_pixels;i++)map_aux[i]=map[i];
      pointing point_rev;
      int ip=0;
      int no_pix=0;
      for(int ir=0;ir<nr;ir++){
	int Npixels_ring=npix_ring(nside, ir);
	for(int ipr=0;ipr<Npixels_ring;ipr++){
          point_rev=map_aux.pix2ang(ip);
	  theta_new[ir]=point_rev.theta;
	  phi_new[ip]=point_rev.phi;
          if(map_aux[ip]==0)++no_pix;
	  ++ip;
	}
      }
      n_observed_pixels=map.Npix()-no_pix;
      n_total_pixels=map.Npix();
    }
  }
  cout<<"Cat2Map_noz done"<<endl;
  
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// Compute the overdensity in pixels to write it in a map
void AngularPowerF::Cat2Map_delta(char c, Healpix_Map<healpix_real>&map){
  int ira=(c=='r'? this->i_alpha_ran:this->i_alpha);
  int idec=(c=='r'? this->i_delta_ran:this->i_delta);
  int iz=(c=='r'? this->i_z_ran:this->i_z);
  int ncols=(c=='r'? this->n_columns_ran:this->n_columns_gal);
  vector<real_prec>prop_a=(c=='r'? this->  prop_r:this->  prop);
  map.fill(0.0);
  for(int i=0;i<prop_a.size();i++){
    long pixel;
    point.phi=prop_a[ira+ncols*i]*fac;
    point.theta=0.5*M_PI-prop_a[idec+ncols*i]*fac;
    pixel=map.ang2pix(point);
    if(pixmask[pixel]==1){
      if(prop_a[iz+ncols*i]< z_max[this->IZ] && prop_a[iz+ncols*i]>=z_min[this->IZ]){     // Now do the map in the redshift bin:
	map[pixel]++;
      }
    }
  }
  prop_a.clear();
  
  string fileo= this->name_output_dir+"2MPZ_delta_pix_dist_"+this->hemis+"_"+this->coord+"_zbin_"+to_string(this->IZ)+".dat";
  ofstream sout; sout.open(fileo.c_str());
  healpix_real ngala=this->mean_number_galaxies_pix;
  for(int i=0;i<map.Npix();++i)map[i]=this->mean_number_galaxies_pix==0.0 ? 0.0: ((healpix_real)map[i]-ngala)/ngala ;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)sout<<log(1+map[i])<<endl;
  int nnpix=0;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)nnpix++;
  cout<<CYAN<<"File "<<fileo<<" written with log (1+delta)"<<RESET<<endl;

  healpix_real mean_lg=0;
  healpix_real var_lg=0;
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)mean_lg+=log(1+map[i])/((healpix_real)nnpix);
  for(int i=0;i<map.Npix();++i)if(map[i]!=-1)var_lg+=pow(log(1+map[i]) - mean_lg,2)/((healpix_real)nnpix);
  var_lg=sqrt(var_lg);

  cout<<CYAN<<"Mean log normal = "<<mean_lg<<"  Variance log normal = "<<var_lg<<RESET<<endl;

  cout<<"Cat2Map_delta done"<<endl;
  sout.close();
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::get_pars_cl_estimator(){
  cout<<RED;
  this->sky_fraction=(healpix_real)n_observed_pixels/(healpix_real)(n_total_pixels);
  this->area_pixel=4.*M_PI/(healpix_real)n_total_pixels;
  this->area_survey=((healpix_real)n_observed_pixels)*area_pixel;
  this->shot_noise=area_survey/((healpix_real)ngal_used);
  this->alpha=((healpix_real)ngal_used)/((healpix_real)nran_used);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::set_mean_ngal_pix(int nb){
  this->Mean_ngal_pix.resize(nb+1,0);
  this->Mean_ngal.resize(nb+1,0);
  this->Shot_Noise.resize(nb+1,0);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::get_pars_mask(){
  cout<<RED;
  this->area_pixel=4.*M_PI/(healpix_real)n_total_pixels;
  this->area_survey=(healpix_real)n_observed_pixels*area_pixel;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::write_pars_cl_estimator(){
  cout<<CYAN<<"*******************************************************"<<endl;  
  cout<<"Information computed from the galaxy catalog: "<<endl;
  if(this->n_z_bins!=-1){
    cout<<"Redshift bin "<<IZ<<endl;
    cout<<"z min = "<<this->z_min[IZ]<<endl;
    cout<<"z max = "<<this->z_max[IZ]<<endl;
  }
  cout<<"Number of galaxies in catalogue= "<<this->ngal<<endl;
  cout<<"Number of galaxies used = "<<ngal_used<<endl;
  if(use_random_cat)cout<<"Number of randoms in catalogue = "<<nran<<endl;
  if(use_random_cat)cout<<"Number of randoms used = "<<nran_used<<endl;
  if(use_random_cat)cout<<"alpha ="<<alpha<<endl;
  
  cout<<"Number of pixels in the mask = "<<n_total_pixels<<endl;
  cout<<"Number of pixels in the suveyed area = "<<n_observed_pixels<<endl;
  cout<<"Skyfraction = "<<sky_fraction<<endl;  
  cout<<"Area survey = "<<area_survey<<endl;  
  cout<<"Area pixel = "<<area_pixel<<endl;  
  cout<<"Shot-noise = "<<shot_noise<<endl;
  cout<<"Mean number of galaxies in pixels = "<<mean_number_galaxies_pix<<endl;
  cout<<"Mean number of galaxies= "<<ngal_used/area_survey<<endl;
  cout<<"RMS  = "<<rms_ngal<<endl;
  if(use_random_cat)cout<<"Mean number of randoms  in pixels = "<<mean_number_randoms_pix<<endl;
  if(use_random_cat) cout<<"RMS  = "<<rms_nran<<endl;
  std::cout<<"*******************************************************"<<endl;  
  cout<<RESET;

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::Map2Alm(Healpix_Map<healpix_real>map, Alm<xcomplex <healpix_real> >&alm,  Alm<xcomplex <healpix_real> >&Ilm, arr<arr<healpix_real> >&Jlm){

  time_t start;  
  time (&start);
  int iring=0;

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    healpix_real x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside, ir);
    for(int m=0;m<=Lmax;m++){
     healpix_real *Pl= new healpix_real[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
      	if(pixmask[ip]==1){	
          healpix_real cc=cos(healpix_real(m)*phi_new[ip]);
          healpix_real ss=sin(healpix_real(m)*phi_new[ip]);
   	  int il=0;
   	  for(int l=m;l<=Lmax;l++){	
            healpix_real Plm=Pl[il];
            healpix_real Yreal=Plm*cc;
            healpix_real Yimag=Plm*ss;
	    /*
	    alm(l,m).re+= area_pixel*Yreal*map[ip];
   	    alm(l,m).im+=-area_pixel*Yimag*map[ip];
	    Ilm(l,m).re+= area_pixel*Yreal;
	    Ilm(l,m).im+=-area_pixel*Yimag;*/
	    Jlm[l][m]  += area_pixel*pow(Plm,2);
	    ++il;
   	  }
   	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::Map2Alm(Healpix_Map<healpix_real>map, Alm<xcomplex <healpix_real> >&alm){

  // This function does the job of the map2alm pf Healpix in longer time

  time_t start;
  time (&start);
  int iring=0;

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    healpix_real x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside,ir);

    healpix_real a1=0,a2=0;
    for(int m=0;m<=Lmax;m++){
      healpix_real *Pl= new healpix_real[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
	if(pixmask[ip]==1){	
          healpix_real cc=cos(healpix_real(m)*phi_new[ip]);
          healpix_real ss=sin(healpix_real(m)*phi_new[ip]);
   	  int il=0;
   	  for(int l=m;l<=Lmax;l++){	
            healpix_real Plm=Pl[il];
            healpix_real Yreal=Plm*cc;
            healpix_real Yimag=Plm*ss;
	    a1+=area_pixel*Yreal*map[ip];
	    alm(l,m).real(a1); 
	    a2+=-area_pixel*Yimag*map[ip];
   	    alm(l,m).imag(a2);
   	    ++il;
   	  }
	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }



}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::set_ilm_jlm(){
  this->Jlm.resize(this->Lmax+1);
  for(int i=0;i<Jlm.size();++i)this->Jlm[i].resize(this->Lmax+1,1.0);
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::Map2Ilm_Jlm(Alm<xcomplex <healpix_real> >&Ilm){

  time_t start;  
  time (&start);

  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)Ilm(i,im).real(0);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)Ilm(i,im).imag(0);


  int iring=0;
  // ******************************************************************** 
  // In this part we get the Ilm using HEALPIX functions 
  
  cout<<CYAN<<"Computing Ilm...";
  Healpix_Map<healpix_real>pmask(log2(nside), RING);
  arr<healpix_real>weight_data(2*pmask.Nside(),1.0);
  for(int i=0;i<n_pixels;i++)pmask[i]=this->pixmask[i];
  map2alm(pmask,Ilm,weight_data,false);
//  map2alm_iter(pmask,Ilm,NUM_ITER_MAP2ALM, weight_data);


  cout<<"Done."<<RESET<<endl;
  
  // ********************************************************************

  if(this->type_of_P_estimator=="D"){
      cout<<CYAN<<"Computing Jlm";
      if(this->compute_jlm){
       for(int l=this->Lmin; l<this->Lmax;++l)for(int m=this->Lmin; m<this->Lmax;++m)this->Jlm[l][m]=0;
       cout<<this->compute_jlm<<"  "<<Nrings<<endl;

       for(int ir=0;ir<Nrings;ir++){
         healpix_real x=cos(theta_new[ir]);
         int Npixels_ring=npix_ring(nside, ir);
         for(int m=0;m<=Lmax;m++){
           healpix_real *Pl= new healpix_real[this->Lmax-m+1];
           gsl_sf_legendre_sphPlm_array(this->Lmax,m,x,Pl);
           for(int ipr=0;ipr<Npixels_ring;ipr++){
             int ip=ipr+iring;
             if(pixmask[ip]==1){
              int il=0;
              for(int l=m;l<=this->Lmax;l++){
                healpix_real Plm=Pl[il];
                this->Jlm[l][m]  += this->area_pixel*pow(Plm,2);
                ++il;
              }
            }
          }
          delete[] Pl;
        }
        iring+=Npixels_ring;
      }

      ofstream jlm_file;
      jlm_file.open(this->output_file_Jlm.c_str());
      for(int l=this->Lmin;l<=this->Lmax;++l){
        for(int m=0;m<=this->Lmax;++m){
          jlm_file<<l<<"\t"<<m<<"\t"<<this->Jlm[l][m]<<endl;
        }
      }
      jlm_file.close();
     }
  
    else if(!this->compute_jlm){
      vector<vector<real_prec> > propJ;
      Fmd.read_file(this->output_file_Jlm,propJ);
      int ik=-1;
      for(int l=this->Lmin;l<=this->Lmax;++l){
        for(int m=0;m<=this->Lmax;++m){
          ik++;this->Jlm[l][m]= propJ[ik][2];
        }
      }
      propJ.clear();
    }
    cout<<"Done"<<RESET<<endl;
    }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::Map2Jlm(Healpix_Map<healpix_real>&map_ran){

  time_t start;  
  time (&start);
  for(int i=0;i<Lmax+1;i++)for(int im=0;im<Lmax+1;im++)this->Jlm[i][im]=0;
  int iring=0;
  for(int ir=0;ir<Nrings;ir++){
    healpix_real x=cos(theta_new[ir]);
    int Npixels_ring=npix_ring(nside, ir);
    for(int m=0;m<=Lmax;m++){
      healpix_real *Pl= new healpix_real[Lmax-m+1];
      gsl_sf_legendre_sphPlm_array(Lmax,m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;      
        healpix_real cc=cos(healpix_real(m)*phi_new[ip]);
        healpix_real ss=sin(healpix_real(m)*phi_new[ip]);
	int il=0;
	for(int l=m;l<=Lmax;l++){	
          healpix_real Plm=Pl[il];
          healpix_real Yreal=Plm*cc;
          healpix_real Yimag=Plm*ss;
	  this->Jlm[l][m]+= area_pixel*pow(Plm,2)*map_ran[ip];
	  ++il;
	}
      }
      delete[] Pl;
    }
    iring+=Npixels_ring;
  }
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************



// This is the one

void  AngularPowerF::Alm2Cl(int iz, int jz,  Alm<xcomplex <healpix_real> >&Ilm, vector<healpix_real> &Cl){
  healpix_real sn=(iz==jz? Shot_Noise[iz]:0);
  healpix_real norm_fsky=1./sn; // The mean surface number density
  sn = (this->shot_noise_correction=="yes"? sn:0);
  healpix_real FSky = this->type_of_P_estimator=="D" ? 1.0 : this->sky_fraction;

  healpix_real Mean_ngal_pix_i=Mean_ngal_pix[iz];
  healpix_real Mean_ngal_pix_j=Mean_ngal_pix[jz];


  for(int l=this->Lmin;l<=this->Lmax;l++)this->Wl[l]=0;
  for(int l=this->Lmin;l<=this->Lmax;l++)Cl[l]=0;
  
  if(this->sky == "masked_sky"){
    if(this->sampling=="H"){
      
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){

          healpix_real jlm=this->type_of_P_estimator=="D" ? this->Jlm[l][m] : 1.0 ;
          healpix_real Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal_pix_i*Ilm(l,m).real())/Mean_ngal_pix_i;
          healpix_real Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal_pix_i*Ilm(l,m).imag())/Mean_ngal_pix_i;
          healpix_real Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal_pix_j*Ilm(l,m).real())/Mean_ngal_pix_j;
          healpix_real Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal_pix_j*Ilm(l,m).imag())/Mean_ngal_pix_j;
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0){Cl[l]/=2.; this->Wl[l]/=2.;}
	}
	
  //      healpix_real cor=pow(this-> pixel_window[l],-1);
        Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
        this->Wl[l]=Wl[l]/(l+0.5);
        this->lvec[l]=(healpix_real)l;
      }
    }    
    else{ //if direct sum. Here we are not sure about the implementation of the Ilm, which comes from a Masked. Should we use randoms?
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
          healpix_real jlm=this->type_of_P_estimator=="D" ? this->Jlm[l][m] : 1.0 ;
          healpix_real Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal[iz]*Ilm(l,m).real())/Mean_ngal[iz];
          healpix_real Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal[iz]*Ilm(l,m).imag())/Mean_ngal[iz];
          healpix_real Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal[jz]*Ilm(l,m).real())/Mean_ngal[jz];
          healpix_real Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal[jz]*Ilm(l,m).imag())/Mean_ngal[jz];
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
          this->Wl[l]+=norm(Ilm(l,m));
          if(m==0){Cl[l]/=2.;this->Wl[l]/=2.;}
	}
	Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
	this->Wl[l]/=(l+0.5);
        this->lvec[l]=(healpix_real)l;
      }
    }
  }

  
  else {  // If full sky. In this case, the estimator D and K are equivalent, Jlm=1, Ilm=0
    
    if(this->sampling=="H"){
      
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
          healpix_real Br_z1=this->Blm[l][m][iz].real()/Mean_ngal_pix[iz];
          healpix_real Bi_z1=this->Blm[l][m][iz].imag()/Mean_ngal_pix[iz];
          healpix_real Br_z2=this->Blm[l][m][jz].real()/Mean_ngal_pix[jz];
          healpix_real Bi_z2=this->Blm[l][m][jz].imag()/Mean_ngal_pix[jz];
	  
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)Cl[l]/=2.;
	  if(m==0) this->Wl[l]/=2.;
	}
        Cl[l]=(Cl[l]/(l+0.5))-sn;  // FOR D
        this->Wl[l]=Wl[l]/(l+0.5);
        this->lvec[l]=(healpix_real)l;
      }
    }
    
    else{ //if Direct summation && full sky
      for(int l=this->Lmin;l<=this->Lmax;l++){
	for(int m=0;m<=l;m++){
          healpix_real Br_z1=this->Blm[l][m][iz].real();
          healpix_real Bi_z1=this->Blm[l][m][iz].imag();
          healpix_real Br_z2=this->Blm[l][m][jz].real();
          healpix_real Bi_z2=this->Blm[l][m][jz].imag();
	  
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)Cl[l]/=2.;
	  if(m==0) this->Wl[l]/=2.;
	}
	Cl[l]=Cl[l]/(l+0.5)/pow(norm_fsky,2)-sn;  // FOR D
        this->Wl[l]/=(l+0.5);  // Applied only to the window function correctly?
        this->lvec[l]=(healpix_real)l;
      }
    }
  }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::Cl_bins(vector<healpix_real>Cl,  vector<healpix_real>&Clbin){

  cout<<BLUE<<endl;
  if(this->bin_type=="linear"){ 
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    for(int i=0;i<(signed)lbin.size();i++)Clbin[i]=0;
    
    healpix_real deltal=((healpix_real)(this->Lmax-this->Lmin))/((healpix_real)this->lbin.size());
    for(int l=this->Lmin;l<=this->Lmax;l++){
      int il=(int)floor((healpix_real)(l-this->Lmin)/((healpix_real)deltal));
      if(il==lbin.size())il--;
      Clbin[il]+=(l+0.5)*Cl[l];
      this->nmodes[il]+=(l+0.5);
    }
    for(int i=0;i<lbin.size();++i)Clbin[i]/=this->nmodes[i];
  }
  else{
    if(this->bin_type=="log"){
      healpix_real deltal=log10(Lmax/1.0)/((healpix_real)this->lbin.size());
      for(int i=0;i<lbin.size();++i)this->nmodes[i]=0;
      for(int l=0;l<=Lmax;l++){
        int il=(int)floor((log10(l)-log10(1.0))/deltal);
        if(il==lbin.size())il--;
        Clbin[il]+=(l+0.5)*Cl[l];
	this->nmodes[il]+=(l+0.5);
      }
      for(int i=0;i<lbin.size();++i)Clbin[i]/=this->nmodes[i];
    }
  }
  cout<<RESET;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::eCl_bins(vector<healpix_real>Cl,  vector<healpix_real>&Clbin){

  cout<<BLUE<<endl;
  if(this->bin_type=="linear"){
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    for(int i=0;i<(signed)lbin.size();i++)Clbin[i]=0;

    healpix_real deltal=((healpix_real)(this->Lmax-this->Lmin))/((healpix_real)this->lbin.size());
    for(int l=this->Lmin;l<=this->Lmax;l++){
      int il=(int)floor((healpix_real)(l-this->Lmin)/((healpix_real)deltal));
      if(il==lbin.size())il--;
      Clbin[il]+=pow((l+0.5)*Cl[l],2);
      this->nmodes[il]+=(l+0.5);
    }
    for(int i=0;i<lbin.size();i++)Clbin[i]=sqrt(Clbin[i])/this->nmodes[i];
  }
  else{
    if(this->bin_type=="log"){
      healpix_real deltal=log10(Lmax/1.0)/((healpix_real)this->lbin.size());
      for(int i=0;i<lbin.size();++i)this->nmodes[i]=0;
      for(int l=0;l<=Lmax;l++){
        int il=(int)floor((log10(l)-log10(1.0))/deltal);
        if(il==lbin.size())il--;
        Clbin[il]+=pow(Cl[l],2);
    this->nmodes[il]++;
      }
      for(int i=0;i<lbin.size();i++)Clbin[i]=sqrt(Clbin[i])/this->nmodes[i];
    }
  }
  cout<<RESET;
}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::W3J(){
  cout<<RED<<"Computing Wigner Symbols"<<endl;
  Wigner3J.resize(Lmax+1);
  for(int i=0;i<Lmax+1;i++){
    Wigner3J[i].resize(Lmax+1);
    for(int j=0;j<Lmax+1;j++){
      Wigner3J[i][j].resize(Lmax+1);
    }
  }

  time_t start;  
  time (&start);
  for(int i=0;i<Lmax+1;i++){
    for(int j=0;j<Lmax+1;j++){
      for(int k=0;k<Lmax+1;k++){
	int J=i+j+k;
	if(abs(i-j) <= k &&  k<=i+j){
	  if(J%2==0){
	    Wigner3J[i][j][k]= sqrt(factorial(J-2*i)*factorial(J-2*j)*factorial(J-2*k)/factorial(J+1))*(factorial(J/2))/(factorial(J/2-i)*factorial(J/2-j)*factorial(J/2-k));
	  }
	  else{
	    Wigner3J[i][j][k]=0;
	  }
	}
      }
    }
  }
  
  cout<<"DONE"<<RESET<<endl;								       
  
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerF::Mixing_matrix(){
  cout<<"Computing the mixing matrix ..."<<endl;
  for(int l1=Lmin;l1<Lmax;l1++){
    for(int l2=Lmin;l2<Lmax+1;l2++){
      healpix_real Rbis=0;
      for(int l3=Lmin;l3<Lmax+1;l3++){
       	if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
          healpix_real WC = Wigner3J[l1][l2][l3];
	  Rbis+=(2.*l3+1)*this->Wl[l3]*(WC*WC);
	}
      }
      this->R[l1][l2]=((2.*l2+1)/(4.*M_PI))*Rbis;
    }
  }
  cout<<"Done"<<RESET<<endl;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void AngularPowerF::Rll_bins(){


  // Only available for the Rll written in term of the Wigner symbols

  healpix_real deltal;
  
  if(this->bin_type=="linear"){ 
    
    for(int i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;
    
    deltal=(this->Lmax-this->Lmin)/((healpix_real)lbin.size()); //here and just here I define N_L_bins as a healpix_real instead of an integer
    cout<<RESET;

    for(int l1=this->Lmin;l1<this->Lmax;++l1){
      int lbina=floor((l1-this->Lmin)/deltal);

      if(lbina==lbin.size())lbina--;
      this->nmodes[lbina]+=(2.*l1+1);
      for(int l2=this->Lmin;l2<this->Lmax+1;++l2){
        healpix_real Rbis=0;
    	for(int l3=this->Lmin;l3<this->Lmax+1;l3++){
    	  if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
            healpix_real WC = this->Wigner3J[l1][l2][l3];
  	    Rbis+=(2*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	  }
  	}
        this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
      }
    }
    for(int i=0;i<lbin.size();i++)for(int l2=this->Lmin;l2<this->Lmax+1;l2++)Rll_bin[i][l2]/=this->nmodes[i];   }
  else{
    if(this->bin_type=="log"){
      
      deltal=log10(Lmax/1.0)/((healpix_real)lbin.size());
      cout<<RESET;
      for(int i=0;i<lbin.size();i++)this->nmodes[i]=0;
      
      for(int l1=this->Lmin;l1<this->Lmax;l1++){
  	int lbina=floor((log(l1)-log(1.0))/deltal);
	this->nmodes[lbina]+=(2.*l1+1);
  	for(int l2=this->Lmin;l2<this->Lmax+1;l2++){
          healpix_real Rbis=0;
  	  for(int l3=this->Lmin;l3<this->Lmax+1;l3++){
  	    if(abs(l1-l2) <= l3 &&  l3<=l1+l2){
              healpix_real WC = this->Wigner3J[l1][l2][l3];
              Rbis+=(2.*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	    }
  	  }
          this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
  	}
      }
      
      for(int i=0;i<lbin.size();i++)for(int l2=this->Lmin;l2<this->Lmax+1;l2++)Rll_bin[i][l2]/=this->nmodes[i]; 
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerF::get_mixing_matrix(){

  // compute the wigner symbols
  this->W3J();

  this->Mixing_matrix();
  healpix_real factor_fs = (this->mixing_matrix_exact? 1.0:this->sky_fraction);
  for(int l1=this->Lmin;l1<=this->Lmax;l1++)for(int l2=this->Lmin;l2<=this->Lmax;l2++)this->R[l1][l2]/=factor_fs;
  for(int l2=this->Lmin;l2<=this->Lmax;l2++){
    vector<healpix_real> slice (this->Lmax+1,0);
    string rfile=this->output_file_window+"_MixingMatrix_l_"+to_string(l2)+".dat";
    for(int l1=this->Lmin;l1<this->Lmax;l1++)slice[l1]=this->R[l1][l2];
    ////////////////////////////////////////////////this->Fmi.write_to_file(rfile, this->lvec, slice);
  }
  string rfile;
  if(this->mixing_matrix_exact)rfile=this->output_file_window+"_MixingMatrix_exact_nside_"+to_string(this->nside)+".dat";
  else rfile=this->output_file_window+"_MixingMatrix_nside_"+to_string(this->nside)+".dat";
  ////////////////////////////////////////////////this->Fmi.write_to_file(rfile, this->lvec, this->lvec, this->R);

  this->Rll_bins();
  for(int l1=0;l1<this->N_L_bins;l1++)for(int l2=this->Lmin;l2<=this->Lmax;l2++)this->Rll_bin[l1][l2]/=factor_fs;
  rfile=this->output_file_window+"_MixingMatrix_lbins_"+to_string(this->N_L_bins)+"_nside_"+to_string(this->nside)+".dat";
  ////////////////////////////////////////////////this->Fmi.write_to_file(rfile, this->lbin, this->lvec, this->Rll_bin);


  for(int l1=0;l1<this->N_L_bins;l1++){
    vector<healpix_real> slice3 (this->Lmax+1,0);
    for(int l2=this->Lmin;l2<=this->Lmax;l2++)slice3[l2]=this->Rll_bin[l1][l2];
    rfile=this->output_file_window+"_MixingMatrix_lbin_"+to_string(l1)+"_nside_"+to_string(this->nside)+".dat";
    ////////////////////////////////////////////////this->Fmi.write_to_file(rfile, this->lvec,slice3);
  }

}





// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// ######################################################################
// ######################################################################
// ######################################################################
// ######################################################################
healpix_real AngularPowerF::get_min(char c, int ii){

  vector<real_prec> pro=(c=='d'? this->  prop:this->  prop_r);
  int ncols=(c=='d'? this->n_columns_gal:this->n_columns_ran);
  real_prec ZX=10000;
  real_prec zmx;
  for(int i=0;i<pro.size()/ncols;i++){
    zmx=min(ZX,pro[ii+ncols*i]);
    ZX=zmx;
  }
  return zmx;
}
// ######################################################################
// ######################################################################
healpix_real AngularPowerF::get_max(char c, int ii ){
  vector<real_prec> pro=(c=='d'? this->  prop  :this->  prop_r);
  int ncols=(c=='d'? this->n_columns_gal:this->n_columns_ran);
  real_prec ZX=-100;
  real_prec zmx;
  for(int i=0;i<pro.size()/ncols;i++){
    zmx=max(ZX,pro[ii+ncols*i]);
    ZX=zmx;
  }
  return zmx;
}

// ######################################################################
// ######################################################################
void AngularPowerF::read_input_cats(string cat, string file){
   if(cat=="g")
       this->n_columns_gal=Fmd.read_file(file, this-> prop);
   else if(cat=="r")
       this->n_columns_ran=Fmd.read_file(file, this->prop_r);
}



// ######################################################################
// ######################################################################

void AngularPowerF::get_zbins_same_ngal(int nbins, int iiz, healpix_real zmn, healpix_real zmx, vector<vector<healpix_real> >&zzv){

  // ************************************
  // Construct the dNdz with many bins
  int na=20000; //This value is critical to avoid seg foults. The higher, the best.
  healpix_real delta=(zmx-zmn)/((healpix_real)na);
  vector<int>dn(na,0);
  for(int i=0;i<this->  prop.size()/this->n_columns_gal;++i){
    if(this->  prop[iiz+this->n_columns_gal*i]<zmx && this->  prop[iiz+this->n_columns_gal*i]>=zmn){
      int iza=floor((this->  prop[iiz+this->n_columns_gal*i]-zmn)/delta);
      dn[iza]++;
    }
  }

  int Nca=(int)(floor)(  prop.size()/this->n_columns_gal/((healpix_real)nbins)); //Desired number of galaxies per redshift bin:
  cout<<"Number of galaxies per redshift bin = "<<Nca<<endl;
  vector<healpix_real>zan(na,0);
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
      vector<healpix_real>zaux;
      for(int i=0;i<dn.size();++i){
        caa+=dn[i];  //Cumulative number of galaxies
        if((caa>=Nca*(ib-1)) &&  (caa<Nca*ib)) zaux.push_back(zan[i]);
      }
      zzv[ib][0]=zaux[0]-0.5*delta;  //Allocate the zmin of the ib zbin
      zzv[ib][1]=zaux[zaux.size()-1]+0.5*delta; //Allocate the zmax of the ib zbin
      zaux.clear();
    }
  }
 dn.clear();
 return ;
}

