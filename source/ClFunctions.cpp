# include "../headers/ClFunctions.h"

using namespace std;


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_Cl_bias()
{
  // This funciton computes the Catalog bias by taking ratios bvetween Cl from gaalxies and a
  // preduction of the Cl from a dark matter realization (model).
  // Catalog catalog
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
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::get_angular_power()
{
  
  string temp="lbins";
  string temp_raw="raw";
  
  this->params.set_output_file(this->params._name_output_dir()+this->params._statistics()+"_"+this->params._type_of_P_estimator()+"_"+this->params._name_catalog()+"_zbtype_"+this->params._define_z_bins()+"_"+temp);
  this->params.set_output_file_raw(this->params._name_output_dir()+this->params._statistics()+"_"+this->params._type_of_P_estimator()+"_"+this->params._name_catalog()+"_zbtype_"+this->params._define_z_bins()+"_"+temp_raw);
  
  if(this->params._input_catalog_type()!="ascii")
    cout<<"WARNING, catalog not in ascii, "<<endl;

  // Read the mask
  this->catalog.read_mask();
  this->catalog.read_catalog();
  this->set_healpix_pars();  // Compute Healpix numbers
  this->set_vectors();      // Allocate memory for used vectors
  // *********************************************************
  // *********************************************************
  
  if(this->params._n_z_bins()!=-1) // If we have more than one redshift bin, proceed
    {

#ifdef _FULL_VERBOSE_
      cout<<GREEN<<"SELECTION INFO: "<< this->params._selection()<<RESET<<endl;
      cout<<CYAN;
      if(this->params._selection()=="fls"){
        cout<<"Flux limited sample selected. Bins in redshift"<<endl;
      if(this->params._define_z_bins()=="delta")
         cout<<"zBins with constant width"<<endl;
      if(this->params._define_z_bins()=="number")
        cout<<"zBins chosen with equal number of galaxies"<<endl;
      }
#endif
      // Get this->params._zmax()
      real_prec zmax_cat=this->catalog._redshift_max();
      real_prec zmin_cat=this->catalog._redshift_min();
      this->catalog.set_redshift_min(zmin_cat);
      this->catalog.set_redshift_max(zmax_cat);
#ifdef _FULL_VERBOSE_
      cout<<GREEN<<"REDSHIFT INFO "<<RESET<<endl;
      cout<<CYAN<<"Maximum redshift found in catalogue = "<<this->catalog._redshift_max()<<endl;
      cout<<CYAN<<"Minimum redshift found in catalogue = "<<this->catalog._redshift_min()<<endl;
#endif
      
      if(false==this->params._use_z_min_max_from_cat())
        {
          cout<<CYAN<<"Maximim redshift selected in parameter file = "<<this->params._zmax()<<endl;
          cout<<CYAN<<"Minimum redshift selected in parameter file = "<<this->params._zmin()<<endl;
        }
      else
        {
          this->params.set_zmin(zmin_cat);
          this->params.set_zmax(zmax_cat);
        }
      
      // Get Magnitude limits if the MK is used
      if(this->params._i_M()>0)
        {
          this->MKmax_cat=this->catalog._MK_max();
          this->MKmax_cat=this->catalog._MK_min();
#ifdef _FULL_VERBOSE_
          cout<<RED<<"MAGNITUDE Mk INFO "<<RESET<<endl;
          cout<<CYAN<<"Maximim Mk found in catalogue = "<<this->catalog._MK_max()<<endl;
          cout<<CYAN<<"Minimum Mk found in catalogue = "<<this->catalog._MK_min()<<endl;
#endif
          if(false==this->params._use_Mk_min_max_from_cat()){
#ifdef _FULL_VERBOSE_
            cout<<CYAN<<"Maximim Mk selected in parameter file = "<<this->params._MKmax()<<endl;
            cout<<CYAN<<"Minimum Mk selected in parameter file = "<<this->params._MKmin()<<endl;
#endif
          }
	  else
	    {
          this->params.set_MKmin(this->catalog._MK_min());
          this->params.set_MKmax(this->catalog._MK_max());
	    }
#ifdef _FULL_VERBOSE_
	  std::cout<<"*******************************************************"<<endl;
#endif
	  // SOME WARNINGS:
	  std::cout<<RED;
#ifdef _FULL_VERBOSE_
      if(this->params._MKmax() > this->catalog._MK_max())
          cout<<"Warning: Maximum Mk in parameter file ("<<this->params._MKmax()<<") greater than maximum value found in the catalog ("<<this->catalog._MK_max()<<")"<<endl;
      if(this->params._MKmin() < this->catalog._MK_min())
          cout<<"Warning: Minimum Mk in parameter file ("<<this->params._MKmin()<<") smaller than minimum value found in the catalog ("<<this->catalog._MK_min() <<")"<<endl;
#endif
	}
      
      if(this->params._zmax() > this->catalog._redshift_max())
	{
#ifdef _FULL_VERBOSE_
      cout<<"Warning: Maximum z in parameter file ("<<this->params._zmax()<<") greater than maximum value found in the catalog ("<<this->catalog._redshift_max()<<"). Doing nothing"<<endl;
#endif
	  //   this->this->params._zmax()= this->this->params._zmax()_cat;
	}
      if(this->params._zmin() < this->catalog._redshift_min())
	{
#ifdef _FULL_VERBOSE_
      cout<<"Warning: Minimum z in parameter file ("<<this->params._zmin()<<") smaller than minimum value found in the catalog ("<<this->catalog._redshift_min()<<"). Doing nothing"<<endl;
#endif
	  
	}
    } // end if (n_zbins!=-1)
  else
    {
      this->params.set_zmin(-1000);
      this->params.set_zmax(1000);
    }
  
  // *******************************************************
  
  int NB=0;
  if(this->params._n_z_bins()>=1)
   {
      NB =this->params._n_z_bins();
      this->catalog.set_zbins();
  }
  this->set_BLMZ();
  
  
  // ******************************************************************************************************************
  // *********************************************************
  // OPERATIONS IN HARMONIC SPACE
  // *********************************************************
  // *********************************************************
  cout<<CYAN<<"Operations in Harmonic space"<<RESET<<endl;

  Alm<xcomplex <real_prec> > Ilm(this->params._Lmax(),this->params._Lmax());
  // Compute Jlm, Ilm. These objects are inizialized to the
  // full sky- case, Jlm=1, Ilm=0 (l>0)
  // such that if not computed in Map2ILM_Jlm, they remain full sky.
  // but if used, these have to be initialized to zero in that funcition
  // For full sky, Ilm os not set to zero by directly not used in the estimation
  this->Jlm.resize(this->params._Lmax()+1);
  for(ULONG i=0;i<Jlm.size();++i)this->Jlm[i].resize(this->params._Lmax()+1,1.0);
  this->Map2Ilm_Jlm(Ilm);
  // *********************************************************
  // *********************************************************
#ifdef _FULL_VERBOSE_
  cout<<"*******************************************************"<<endl;
  cout<<CYAN<<"Computing Alm in z bins"<<RESET<<endl;
#endif

  ULONG izf=0;
  izf=(this->params._n_z_bins()==1? 0 : this->params._n_z_bins());
  this->set_mean_ngal_pix();
  if(this->params._n_z_bins()==-1)
    izf=0;
  // *******************************************************
  this->set_Lbins();
  // *******************************************************
#ifdef _USE_GNUPLOT_
  this->gp_power << "set log\n";
  this->gp_power << "set xlabel 'l' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_power << "set ylabel 'C_l' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
  this->gp_power << "set size square\n";
#endif

  // *********************REDSHIFT BINS**********************************************************
  // *******************************************************
  // *******************************************************


  // START LOOP OVER THE REDSHIFT BINS ONLY IF FLS
  for(ULONG aIZ=0;aIZ<=izf;++aIZ)
    {
      this->set_IZ(aIZ);
      this->catalog.set_redshift_bin(aIZ);

      Healpix_Map<real_prec>healpix_map(ilog2(this->params._nside()), RING);
      this->catalog.Cat2Map(healpix_map,0, true);
      // ***************************************************************************************
      this->shot_noise=this->catalog._area_survey()/(static_cast<real_prec>(this->catalog._weighted_ngals_in_zbin()));
      this->Shot_Noise[aIZ]=this->shot_noise;// this already has the effect of the weights in the SN
      this->write_pars_cl_estimator();
      this->Mean_ngal_pix[aIZ]=this->catalog._mean_number_galaxies_pix();
      this->Mean_ngal[aIZ]=this->catalog._weighted_ngals_in_zbin()/(this->catalog._area_survey());
      // ***************************************************************************************
      cout<<GREEN<<"Computing Alm for current redshift bin:"<<RESET<<endl;
      Alm<xcomplex <real_prec> > Alm(this->params._Lmax(),this->params._Lmax());
      arr<real_prec>weight_data(2*healpix_map.Nside());
      weight_data.fill(1.0);
      map2alm(healpix_map,Alm, weight_data,false); //HealPix
      cout<<GREEN<<"DONE"<<RESET<<endl;
      // Fill the vectors Blm. These are the Alm in different redshift or Magnitude bins
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG il=0;il<this->params._Lmax()+1;il++)
        for(long im=0;im<this->params._Lmax()+1;im++)
          {
            this->Blm[il][im][aIZ].real(Alm(il,im).real());
            this->Blm[il][im][aIZ].imag(Alm(il,im).imag());
          }
     }
  // ***********************************************************************************************
  // Here we can estimate the shot noise as shuffled realziations of the map, and compute its power,
  // since the shot_noise for the cross is not needed
  // This is equivalent as to have a random catalog with the same angular and redhsft distribution
  
  struct zsn{vector<real_prec>Clvec_sn;};
  vector<zsn> shot_noise_ell(izf+1);
  if(true==this->params._shot_noise_correction())
    if(true==this->params._shot_noise_correction_randomized())
      {
    //Here we do a number of realizations randomizing the map, computing the Cl
    //allocating the results in the container Cl_sn.
        cout<<GREEN<<"Computing shot noise from random realization of gal coordinates using "<<this->params._N_random_shot_noise()<<" realiations"<<RESET<<endl;
        for(ULONG aIZ=0;aIZ<=izf;++aIZ)//Start again redshift bins to compute shot noise
          {
            this->set_IZ(aIZ);
            this->catalog.set_redshift_bin(aIZ);
            int Nran=this->params._N_random_shot_noise();
            shot_noise_ell[aIZ].Clvec_sn.resize(this->params._Lmax()+1,0);
            gsl_rng_env_setup();
            const gsl_rng_type *Tn;
            gsl_rng_default_seed=1014435;
            Tn = gsl_rng_default;
            gsl_rng *rn = gsl_rng_alloc (Tn);
            for(int ir=0;ir<Nran;ir++)
              {
                pointing point_r;
                Healpix_Map<real_prec>healpix_map(ilog2(this->params._nside()), RING);
                for(ULONG i=0;i<this->catalog._NOBJS();++i)
                {
                    ULONG npix=0;
                    bool inmask=false;
                    while(inmask==false)
                     {
                        point_r.phi=(360.0*gsl_rng_uniform(rn))*CFACTOR;
                        point_r.theta=acos(-1+2.0*gsl_rng_uniform(rn));
                        npix=static_cast<ULONG>(healpix_map.ang2pix(point_r));
                        if(1==this->catalog.pixmask[npix])
                            inmask=true;
                     }
                    this->catalog.tracer[i].galPIXEL=npix;
                }
                //Get map. We use again healpix_map as its original version has been saved in Blm
                this->catalog.Cat2Map(healpix_map,ir,false);
                //Get Alm
                arr<real_prec>weight_data(2*healpix_map.Nside());
                weight_data.fill(1.0);
                Alm<xcomplex <real_prec> > Alm(this->params._Lmax(),this->params._Lmax());
                map2alm(healpix_map, Alm, weight_data,false);
                // Compute power
                vector<real_prec>Cl_temp(this->params._Lmax()+1,0);
                this->Alm2Cl_sn(aIZ,Alm,Ilm,Cl_temp);
                for(int il=0;il<Cl_temp.size();++il)
                  shot_noise_ell[aIZ].Clvec_sn[il]+=(Cl_temp[il]/static_cast<real_prec>(Nran));
            }
          }
      }  //close loop over z, opened in sample=fls
  
  // ***********************************************************************************************
  
  cout<<BLUE;
  // Having allocated quantities in z bin (fls)
  // we now do a real_prec loop (fls) over the zbins, or
  // we do only one loop (vls) (the second is currently running)
  // and compute things for the current Mbins
  cout<<"*******************************************************"<<endl;
  
  if(this->params._n_z_bins()>1)
    {
#ifdef _FULL_VERBOSE_
      cout<<"Computing Cross Cls"<<RESET<<endl;
#endif

      for(int aIZ=1;aIZ<=izf; aIZ++) // Loop over the first bin
        {

          vector<real_prec> ClzbinI(this->params._Lmax()+1,0);
          this->Alm2Cl(aIZ,aIZ, Ilm, ClzbinI);//Get auto Cl at bin IZ
          if(true==this->params._shot_noise_correction_randomized() && true==this->params._shot_noise_correction())
              for(int il=0;il<ClzbinI.size();++il)
                ClzbinI[il]-=shot_noise_ell[aIZ].Clvec_sn[il];

          vector<real_prec> eClzbinI(this->params._Lmax()+1,0);
          this->get_eCl(aIZ,aIZ,ClzbinI,ClzbinI,eClzbinI);
          string rest="_"+this->params._sampling()+"_Nside_"+to_string(this->params._nside())+"_zbins_"+to_string(aIZ)+"_"+to_string(aIZ)+"_"+this->params._hemis()+"_"+this->params._coord()+".txt";
          string cfile=this->params._output_file_raw()+rest;
          this->Fmi.write_to_file(cfile, this->lvec, ClzbinI, this->Wl,eClzbinI);

          fill(this->Clbin.begin(),this->Clbin.end(),0);
          this->get_Cl_bins(ClzbinI, this->Clbin);
          fill(this->eClbin.begin(),this->eClbin.end(),0);
          this->get_eCl_bins(eClzbinI,this->eClbin);
          string ofile=this->params._output_file()+rest;
          this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
#ifdef _USE_GNUPLOT_
          vector<tuple<real_prec, real_prec, real_prec> > xy_pts_a;
          for(ULONG i=0; i<Clbin.size(); ++i)
            xy_pts_a.push_back(std::make_tuple(lbin[i], this->Clbin[i], this->eClbin[i]));
          if(aIZ==1)
            this->gp_power<<"plot[1:][] " << gp_power.file1d(xy_pts_a) << " u 1:2:3 w e ps 1 pt 7 lt "<<aIZ+1<<" title 'z-bin "<<aIZ<<"'"<<endl;
          else
            this->gp_power<<"replot " << gp_power.file1d(xy_pts_a) << " u 1:2:3 w e ps 1 pt 7 lt "<<aIZ+1<<" title 'z-bin "<<aIZ<<"'"<<endl;
#endif


#ifdef _GET_CROSS_POWER_
          for(ULONG aJZ=aIZ;aJZ<=izf; aJZ++) // Loop over the second bin
          {
              vector<real_prec> Cl(this->params._Lmax()+1,0); // COntainer for cross power
              vector<real_prec> eCl(this->params._Lmax()+1,0); // Container for errors of power
              vector<real_prec> ClzbinJ(this->params._Lmax()+1,0); // Container for auto power in bin J
              // ***********************************************************************************************
              // We need to compute the autos within each loop in order to compute th
              // variance for the cross power spectrum (C_i+sn_j)*(C_j+sn_j).
              // NOTE THAT FIRST WE NEED TO COMPUTE THE AUTO. DO not change ordering here.
              // ***********************************************************************************************

              if(aJZ==aIZ)//Auto Cl in bin J:
              {
              this->Alm2Cl(aJZ,aJZ,Ilm,ClzbinJ);
              if(true==this->params._shot_noise_correction_randomized())
                  for(int il=0;il<Clvec_sn.size();++il)
                    ClzbinJ[il]-=this->Clvec_sn[il];
               }
              //Cross:
              this->Alm2Cl(aIZ,aJZ,Ilm,Cl);

                  //get variance of cross power:
              this->get_eCl(aIZ,aJZ,ClzbinI, ClzbinJ, eCl);
              // ***********************************************************************************************
              // Write raw estimates
              string rest;
              if(this->params._sampling()=="H")
                rest="_"+this->params._sampling()+"_Nside_"+to_string(this-          vector<real_prec> eClzbinI(this->params._Lmax()+1,0);
>params._nside())+"_zbins_"+to_string(aIZ)+"_"+to_string(aJZ)+"_"+this->params._hemis()+"_"+this->params._coord()+".txt";
              else if(this->params._sampling()=="DS")
                rest="_"+this->params._sampling()+"_zbins_"+to_string(aIZ)+"_"+to_string(aJZ)+"_"+this->params._hemis()+"_"+this->params._coord()+".txt";
		  
              string cfile=this->params._output_file_raw()+rest;
              this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl,eCl);
              // ***********************************************************************************************
              // **************************************
              // Compute Cl in Bins of l
              this->Cl_bins(Cl, this->Clbin);
              // ***********************************************************************************************
              if(aIZ==aJZ) // For auto power get lmax before SN domination
                {
                  int ll=-1;
                  do{ll++;}
                  while(this->Clbin[ll]>this->Shot_Noise[aJZ]);
                  if(ll==this->lbin.size())
                    std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<aJZ<<". Select lmax = 100. "<<RESET<<std::endl;
                  else
                    std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< aJZ<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
                }
              // ***********************************************************************************************
              // Comute the Variance in Bins of L
              this->eCl_bins(eCl,this->eClbin);
              /*  Get power of mask in Bins */
              this->Cl_bins(this->Wl, this->Wlbin);
              // ***********************************************************************************************
              // Write Binned estimates
              string ofile=this->params._output_file()+rest;
              this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
            }
#endif // enf if cross is available
      }
	  

      {
	// Get the Cl of the full sample
	vector<real_prec> Cl(this->params._Lmax()+1,0);
	vector<real_prec> eCl(this->params._Lmax()+1,0);
	this->Alm2Cl(0,0, Ilm,Cl);
	this->get_eCl(0,0,Cl,Cl,eCl);
	// ***********************************************************************************************
    string rest="_"+this->params._sampling()+"_Nside_"+to_string(this->params._nside())+"_zbins_"+to_string(0)+"_"+to_string(0)+"_"+this->params._hemis()+"_"+this->params._coord()+".txt";
    string cfile=this->params._output_file_raw()+rest;
	this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
	// ***********************************************************************************************
    this->get_Cl_bins(Cl,  this->Clbin);
    this->get_eCl_bins(eCl, this->eClbin);
	// ***********************************************************************************************
    this->get_Cl_bins(this->Wl, this->Wlbin);

	{
	  int ll=-1;
	  do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
	  if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
	  else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
	}
	// ***********************************************************************************************
	// ***********************************************************************************************
	/* string ofile=this->params._output_file()+rest; */
	/* this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin); */
      }
    }      
    else if(this->params._n_z_bins()<=1)
    {   // Get the Cl when no redshift info is available
      vector<real_prec> Cl(this->params._Lmax()+1,0);
      vector<real_prec> eCl(this->params._Lmax()+1,0);
      this->Alm2Cl(0,0, Ilm,Cl);
      this->get_eCl(0,0,Cl,Cl,eCl);
      // ***********************************************************************************************
      string rest="_"+this->params._sampling()+"_Nside_"+to_string(this->params._nside())+"_"+this->params._hemis()+"_"+this->params._coord()+".txt";
      string cfile=this->params._output_file_raw()+rest;
      this->Fmi.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
      // ***********************************************************************************************
      this->get_Cl_bins(Cl,  this->Clbin);
      this->get_eCl_bins(eCl, this->eClbin);
      // ***********************************************************************************************
      this->get_Cl_bins(this->Wl, this->Wlbin);
      string ofile=this->params._output_file()+rest;
      this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
    }
      
  // Finally, get the mixig matrix
  if(true==this->params._compute_mixing_matrix())
    this->get_mixing_matrix();
      
  return;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *****************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_BLMZ(){
  // Resize Blm. I have to do it in an independend method
  this->Blm.resize(this->params._Lmax()+1);
  for(ULONG i=0; i<Blm.size();++i)this->Blm[i].resize(this->params._Lmax()+1);
  for(ULONG i=0; i<this->Blm.size();++i)
      for(ULONG j=0;j<this->Blm[i].size();++j)
          this->Blm[i][j].resize(this->params._n_z_bins()+1);

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_vectors()
{
  this->Wl.resize(this->params._Lmax()+1);
  this->lvec.resize(this->params._Lmax()+1);
  this->lbin.resize(this->params._N_L_bins(),0);
  this->lbin_min.resize(this->params._N_L_bins(),0);
  this->lbin_max.resize(this->params._N_L_bins(),0);
  this->Clbin.resize(this->params._N_L_bins(),0);
  this->eClbin.resize(this->params._N_L_bins(),0);
  this->nmodes.resize(this->params._N_L_bins(),0);
  this->Wlbin.resize(this->params._N_L_bins(),0);

  this->R.resize(this->params._Lmax()+1);
  for(ULONG i=0;i<this->R.size();i++)this->R[i].resize(this->params._Lmax()+1,0);

  // THIS IS ALSO USED WHTN WE READ THE MIXING MATRIX!!!

  if(this->params._compute_mixing_matrix()){
    this->R.resize(this->params._Lmax()+1);
    for(ULONG i=0;i<this->R.size();i++)this->R[i].resize(this->params._Lmax()+1,0);
    Rll_bin.resize(this->params._N_L_bins());
    for(ULONG i=0;i<this->Rll_bin.size();i++)Rll_bin[i].resize(this->params._Lmax()+1,0);
  }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::set_Lbins(){

  if(this->params._bin_type()=="linear")
    {
      real_prec deltal=0;
      if(this->params._use_deltaL()==true)
        deltal=this->params._DeltaL();
      else
        deltal=(static_cast<real_prec>(this->params._Lmax()-this->params._Lmin()))/(static_cast<real_prec>(lbin.size()));
      //here and just here I define N_L_bins as a real_prec instead of an integer
#ifdef _FULL_VERBOSE_
      cout<<CYAN<<endl;
      std::cout<<"*********************************************************"<<endl;
      cout<<BLUE<<"L-bins INFO "<<endl;
      cout<<"Bin-averaged power spectrum"<<endl;
      cout<<"Number of bins = "<<this->lbin.size()<<endl;
      cout<<"l-Bin width = "<<deltal<<RESET<<endl;
      std::cout<<"*********************************************************"<<endl;
      cout<<RESET<<endl;
#endif
          for(ULONG i=0;i<lbin.size();i++)this->lbin[i]=this->params._Lmin()+(i+0.5)*deltal;
          for(ULONG i=0;i<lbin.size();i++)this->lbin_min[i]=this->params._Lmin()+(i)*deltal;
         for(ULONG i=0;i<lbin.size();i++)this->lbin_max[i]=this->params._Lmin()+(i+1)*deltal;
  }
  else
    {
      if(this->params._bin_type()=="log")
	{
      real_prec deltal=log10(params._Lmax()/1.0)/(static_cast<real_prec>(lbin.size()));
          for(ULONG i=0;i<lbin.size();i++)lbin[i]=pow(10,log10(1)+(i+0.5)*deltal);
          for(ULONG i=0;i<lbin.size();i++)lbin_min[i]=pow(10,log10(1)+i*deltal);
          for(ULONG i=0;i<lbin.size();i++)lbin_max[i]=pow(10,log10(1)+(i+1)*deltal);
      }
    }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_healpix_pars(){
  this->n_pixels=12*this->params._nside()*this->params._nside();
  this->Nrings=4*this->params._nside()-1;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::get_eCl(int iz, int jz, vector<real_prec>cli,vector<real_prec>clj, vector<real_prec>&ecl){
  /// WARNING: THE CROSS POWER SPECTRUM HAS UNDESTIMATED GAUSSIAN ERROR BARS, CHECN WHITE, SUNG, PERCIVAL
  
#pragma omp parallel for
  for(long l=this->params._Lmin();l<=this->params._Lmax();++l)
    {
      real_prec cl1=cli[l] + this->Shot_Noise[iz];
      real_prec cl2=clj[l] + this->Shot_Noise[jz];
      real_prec factor=2./(this->catalog._sky_fraction()*(2*l+1));
      ecl[l]=sqrt(factor*cl1*cl2);
    }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::set_mean_ngal_pix(){

  int nb= this->params._n_z_bins()==-1 ? 1 : this->params._n_z_bins();
  this->Mean_ngal_pix.resize(nb+1,0);
  this->Mean_ngal.resize(nb+1,0);
  this->Shot_Noise.resize(nb+1,0);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::write_pars_cl_estimator(){

#ifdef _FULL_VERBOSE_
  cout<<endl;
  cout<<GREEN<<"Information computed from the Catalog catalog: "<<endl;
  cout<<endl;
  if(this->params._n_z_bins()!=-1)
  {
    cout<<"Redshift bin "<<this->IZ<<endl;
    cout<<"z min = "<<this->catalog._z_min(this->IZ)<<endl;
    cout<<"z max = "<<this->catalog._z_max(this->IZ)<<endl;
  }
  cout<<"Number of galaxies in catalog= "<<this->catalog._NOBJS()<<endl;
  cout<<"Number of galaxies in current redshift bin= "<<this->catalog._ngals_in_zbin()<<endl;
  cout<<"Check: Number of galaxies in current redshift bin (frmo map)= "<<this->catalog._ngals_in_zbin_from_map()<<endl;
  cout<<"Weighted (Gauss) number of galaxies in catalogue= "<<this->catalog._weighted_ngals_in_zbin()<<endl;
  cout<<"Number of pixels in the mask = "<<this->catalog._n_pixels()<<endl;
  cout<<"Number of pixels in the suveyed area = "<<this->catalog._n_observed_pixels()<<endl;
  cout<<"Skyfraction = "<<this->catalog._sky_fraction()<<endl;
  cout<<"Area survey = "<<this->catalog._area_survey()<<endl;
  cout<<"Area pixel = "<<this->catalog._area_pixel()<<endl;
  cout<<"Poisson Shot-noise = "<<this->shot_noise<<endl;
  cout<<"Mean number of galaxies in pixels = "<<this->catalog._mean_number_galaxies_pix()<<endl;
  cout<<"Mean number of galaxies= "<<this->catalog._ngals_in_zbin()/this->catalog._area_survey()<<endl;
  cout<<"RMS  = "<<this->catalog._rms_ngal()<<endl;
  cout<<endl;
  cout<<RESET;
#endif

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Map2Alm(Healpix_Map<real_prec>map, Alm<xcomplex <real_prec> >&alm,  Alm<xcomplex <real_prec> >&Ilm, arr<arr<real_prec> >&Jlm){

  time_t start;
  time (&start);
  int iring=0;

  for(long i=0;i<params._Lmax()+1;i++)
      for(long im=0;im<params._Lmax()+1;im++)
          alm(i,im).real(0);
  for(long i=0;i<params._Lmax()+1;i++)
      for(long im=0;im<params._Lmax()+1;im++)
          alm(i,im).imag(0);

  for(ULONG ir=0;ir<Nrings;ir++){
    comp_time(start,Nrings,ir);
    real_prec x=cos(this->catalog.theta_new[ir]);
    int Npixels_ring=npix_ring(this->params._nside(), ir);
    for(long m=0;m<=params._Lmax();m++){
      real_prec *Pl= new real_prec[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(ULONG ipr=0;ipr<Npixels_ring;ipr++)
      {
        ULONG ip=ipr+iring;
        if(this->catalog.pixmask[ip]==1){
//   	  real_prec cc=cos(real_prec(m)*this->catalog.phi_new[ip]);
//  	  real_prec ss=sin(real_prec(m)*this->catalog.phi_new[ip]);
          int il=0;
          for(long l=m;l<=params._Lmax();l++)
          {
            real_prec Plm=Pl[il];

            //real_prec Yreal=Plm*cc;
            //real_prec Yimag=Plm*ss;
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


void Cl_FUNCTIONS::Map2Alm(Healpix_Map<real_prec>map, Alm<xcomplex <real_prec> >&alm){

  // This function does the job of the map2alm pf Healpix in longer time

  time_t start;
  time (&start);
  int iring=0;

  for(long i=0;i<params._Lmax()+1;i++)for(long im=0;im<params._Lmax()+1;im++)alm(i,im).real(0);
  for(long i=0;i<params._Lmax()+1;i++)for(long im=0;im<params._Lmax()+1;im++)alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    real_prec x=cos(this->catalog.theta_new[ir]);
    int Npixels_ring=npix_ring(this->params._nside(),ir);

    real_prec a1=0,a2=0;
    for(long m=0;m<=params._Lmax();m++){
      real_prec *Pl= new real_prec[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++){
	int ip=ipr+iring;
    if(this->catalog.pixmask[ip]==1){
      real_prec cc=cos(real_prec(m)*this->catalog.phi_new[ip]);
      real_prec ss=sin(real_prec(m)*this->catalog.phi_new[ip]);
   	  int il=0;
      for(int l=m;l<=params._Lmax();l++){
  	    real_prec Plm=Pl[il];
   	    real_prec Yreal=Plm*cc;
   	    real_prec Yimag=Plm*ss;
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

void Cl_FUNCTIONS::Map2Ilm_Jlm(Alm<xcomplex <real_prec> >&Ilm){

  time_t start;
  time (&start);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(long i=0;i<params._Lmax()+1;i++)
      for(long im=0;im<params._Lmax()+1;im++)
        {
          Ilm(i,im).real(0);
          Ilm(i,im).imag(0);
        }


  int iring=0;
  // ********************************************************************
  // In this part we get the Ilm using HEALPIX functions
#ifdef _FULL_VERBOSE_
  cout<<CYAN<<"Computing Ilm...";
#endif
  Healpix_Map<real_prec>pmask(ilog2(this->params._nside()), RING);
  arr<real_prec>weight_data(2*pmask.Nside(),1.0);
  for(ULONG i=0;i<n_pixels;i++)
      pmask[i]=this->catalog.pixmask[i];
//  map2alm_iter(pmask,Ilm,NUM_ITER_MAP2ALM, weight_data);
  map2alm(pmask,Ilm,weight_data);

#ifdef _FULL_VERBOSE_
  cout<<"Done."<<RESET<<endl;
#endif
  // ********************************************************************

  if(this->params._type_of_P_estimator()=="D"){
      cout<<CYAN<<"Computing Jlm";
      if(this->params._compute_jlm()){
       for(long l=this->params._Lmin(); l<this->params._Lmax();++l)for(int m=this->params._Lmin(); m<this->params._Lmax();++m)this->Jlm[l][m]=0;
       cout<<this->params._compute_jlm()<<"  "<<Nrings<<endl;

       for(ULONG ir=0;ir<Nrings;ir++){
         comp_time(start,Nrings,ir);
         real_prec x=cos(this->catalog.theta_new[ir]);
         int Npixels_ring=npix_ring(this->params._nside(), ir);
         for(long m=0;m<=params._Lmax();m++){
           real_prec *Pl= new real_prec[this->params._Lmax()-m+1];
           gsl_sf_legendre_sphPlm_array(this->params._Lmax(),m,x,Pl);
           for(ULONG ipr=0;ipr<Npixels_ring;ipr++){
             ULONG ip=ipr+iring;
             if(this->catalog.pixmask[ip]==1){
              ULONG il=0;
              for(long l=m;l<=this->params._Lmax();l++){
                real_prec Plm=Pl[il];
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
      jlm_file.open(this->params._output_file_Jlm().c_str());
      for(long l=this->params._Lmin();l<=this->params._Lmax();++l)
      {
        for(long m=0;m<=this->params._Lmax();++m)
        {
          jlm_file<<l<<"\t"<<m<<"\t"<<this->Jlm[l][m]<<endl;
        }
      }
      jlm_file.close();
     }

    else if(!this->params._compute_jlm()){
      vector<vector<real_prec> > propJ;
      Fmd.read_file(this->params._output_file_Jlm(),propJ);
      int ik=-1;
      for(long l=this->params._Lmin();l<=this->params._Lmax();++l){
        for(long m=0;m<=this->params._Lmax();++m){
          ik++;this->Jlm[l][m]= propJ[ik][2];
        }
      }
      propJ.clear();
    }


  }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::Map2Jlm(Healpix_Map<real_prec>&map_ran){

  time_t start;
  time (&start);
  for(long i=0;i<params._Lmax()+1;i++)for(long im=0;im<params._Lmax()+1;im++)this->Jlm[i][im]=0;
  int iring=0;
  for(ULONG ir=0;ir<Nrings;ir++){
    comp_time(start,Nrings,ir);
    real_prec x=cos(this->catalog.theta_new[ir]);
    ULONG Npixels_ring=npix_ring(this->params._nside(), ir);
    for(long m=0;m<=params._Lmax();m++){
      real_prec *Pl= new real_prec[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(ULONG ipr=0;ipr<Npixels_ring;ipr++){
        ULONG ip=ipr+iring;
    real_prec cc=cos(real_prec(m)*this->catalog.phi_new[ip]);
    real_prec ss=sin(real_prec(m)*this->catalog.phi_new[ip]);
        ULONG il=0;
        for(long l=m;l<=params._Lmax();l++){
	  real_prec Plm=Pl[il];
	  real_prec Yreal=Plm*cc;
	  real_prec Yimag=Plm*ss;
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

void  Cl_FUNCTIONS::Alm2Cl(int iz, int jz,  Alm<xcomplex <real_prec> >&Ilm, vector<real_prec> &Cl){
  real_prec sn=(iz==jz? this->Shot_Noise[iz]:0);
  real_prec norm_fsky=1./sn; // The mean surface number density
  sn = (this->params._shot_noise_correction()==true ? sn:0);
  if(true==this->params._shot_noise_correction_randomized())
      sn=0; // set it to zero as SN will be removed after the call of this method

  real_prec FSky = this->params._type_of_P_estimator()=="D" ? 1.0 : this->catalog._sky_fraction();

  real_prec Mean_ngal_pix_i=0;
  real_prec Mean_ngal_pix_j=0;
  real_prec Mean_ngal_pix_i_aux=1;
  real_prec Mean_ngal_pix_j_aux=1;

  // OJo con esta parte que depende de si hemos noramlizado antes o no

  if(this->params._statistics()=="ARF")
    {
        // For ARF we need sto divide for the mean (weighted ) number of galaxies but we do not subtract the term Ilm (why?)
      Mean_ngal_pix_i=0;
      Mean_ngal_pix_j=0;
      Mean_ngal_pix_i_aux=this->Mean_ngal_pix[iz];
      Mean_ngal_pix_j_aux=this->Mean_ngal_pix[jz];
    }
    else
    {
      Mean_ngal_pix_i=this->Mean_ngal_pix[iz];
      Mean_ngal_pix_j=this->Mean_ngal_pix[jz];
      Mean_ngal_pix_i_aux=Mean_ngal_pix_i;
      Mean_ngal_pix_j_aux=Mean_ngal_pix_j;
    }
  


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int l=this->params._Lmin();l<=this->params._Lmax();l++)
    {
      this->Wl[l]=0;
      Cl[l]=0;
    }
  
  if(this->params._sky() == "masked_sky")
    {
      if(this->params._sampling()=="H")
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(int l=this->params._Lmin();l<=this->params._Lmax();l++)
	    {
	      for(long m=0;m<=l;m++)
          {
		  real_prec jlm=this->params._type_of_P_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
		  real_prec Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal_pix_i*Ilm(l,m).real())/Mean_ngal_pix_i_aux;
		  real_prec Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal_pix_i*Ilm(l,m).imag())/Mean_ngal_pix_i_aux;
		  real_prec Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal_pix_j*Ilm(l,m).real())/Mean_ngal_pix_j_aux;
		  real_prec Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal_pix_j*Ilm(l,m).imag())/Mean_ngal_pix_j_aux;
		  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          this->Wl[l]+=norm(Ilm(l,m));
		  if(m==0)
		    {
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              Cl[l]/=2.;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              this->Wl[l]/=2.;
		    }
		}
	      
          Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          this->Wl[l]/=(l+0.5);

          this->lvec[l]=static_cast<real_prec>(l);
	    }
	}
      else{ //if direct sum. Here we are not sure about the implementation of the Ilm, which comes from a Masked. Should we use randoms?
	for(int l=this->params._Lmin();l<=this->params._Lmax();l++){
	  for(long m=0;m<=l;m++){
	    real_prec jlm=this->params._type_of_P_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
	    real_prec Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal[iz]*Ilm(l,m).real())/Mean_ngal[iz];
	    real_prec Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal[iz]*Ilm(l,m).imag())/Mean_ngal[iz];
	    real_prec Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal[jz]*Ilm(l,m).real())/Mean_ngal[jz];
	    real_prec Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal[jz]*Ilm(l,m).imag())/Mean_ngal[jz];
	    Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
	    this->Wl[l]+=norm(Ilm(l,m));
	    if(m==0)
	      {
		Cl[l]/=2.;
		this->Wl[l]/=2.;
	      }
	  }
	  Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
	  this->Wl[l]/=(l+0.5);
	  this->lvec[l]=static_cast<real_prec>(l);
	}
      }
    }
  
  
  else {  // If full sky. In this case, the estimator D and K are equivalent, Jlm=1, Ilm=0
    
    if(this->params._sampling()=="H"){
      
      for(int l=this->params._Lmin();l<=this->params._Lmax();l++){
        for(long m=0;m<=l;m++){
	  real_prec Br_z1=this->Blm[l][m][iz].real()/Mean_ngal_pix[iz];
	  real_prec Bi_z1=this->Blm[l][m][iz].imag()/Mean_ngal_pix[iz];
	  real_prec Br_z2=this->Blm[l][m][jz].real()/Mean_ngal_pix[jz];
	  real_prec Bi_z2=this->Blm[l][m][jz].imag()/Mean_ngal_pix[jz];

	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)Cl[l]/=2.;
	  if(m==0) this->Wl[l]/=2.;
	}
        Cl[l]=(Cl[l]/(l+0.5))-sn;  // FOR D
        this->Wl[l]=Wl[l]/(l+0.5);
	this->lvec[l]=static_cast<real_prec>(l);
      }
    }

    else{ //if Direct summation && full sky
      for(int l=this->params._Lmin();l<=this->params._Lmax();l++){
        for(long m=0;m<=l;m++){
	  real_prec Br_z1=this->Blm[l][m][iz].real();
	  real_prec Bi_z1=this->Blm[l][m][iz].imag();
	  real_prec Br_z2=this->Blm[l][m][jz].real();
	  real_prec Bi_z2=this->Blm[l][m][jz].imag();

	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  this->Wl[l]+=norm(Ilm(l,m));
	  if(m==0)
	    {
	      Cl[l]/=2.;
	      this->Wl[l]/=2.;
	    }
	}
	Cl[l]=Cl[l]/(l+0.5)/pow(norm_fsky,2)-sn;  // FOR D
        this->Wl[l]/=(l+0.5);  // Applied only to the window function correctly?
	this->lvec[l]=static_cast<real_prec>(l);
      }
    }
  }
}

// ====================================================================================
// ====================================================================================
// ====================================================================================
void  Cl_FUNCTIONS::Alm2Cl_sn(int iz, Alm<xcomplex <real_prec> >&Alm_ran,Alm<xcomplex <real_prec> >&Ilm, vector<real_prec> &Cl){
  real_prec FSky = this->params._type_of_P_estimator()=="D" ? 1.0 : this->catalog._sky_fraction();

  real_prec Mean_ngal_pix_i=0;
  real_prec Mean_ngal_pix_i_aux=1;
  if(this->params._statistics()=="ARF")
    {
      Mean_ngal_pix_i=0;
      Mean_ngal_pix_i_aux=1;
    }
    else if(this->params._statistics()=="Cl")
   {
      Mean_ngal_pix_i=this->Mean_ngal_pix[iz];
      Mean_ngal_pix_i_aux=this->Mean_ngal_pix[iz];
    }

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int l=this->params._Lmin();l<=this->params._Lmax();l++)
      {
        for(long m=0;m<=l;m++)
          {
          real_prec jlm=this->params._type_of_P_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
          real_prec Br=(Alm_ran(l,m).real()-Mean_ngal_pix_i*Ilm(l,m).real())/Mean_ngal_pix_i_aux;
          real_prec Bi=(Alm_ran(l,m).imag()-Mean_ngal_pix_i*Ilm(l,m).imag())/Mean_ngal_pix_i_aux;
          Cl[l]+=(Br*Br+Bi*Bi)/jlm;
          if(m==0)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              Cl[l]/=2.;
        }
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Cl[l]=Cl[l]/(FSky*(l+0.5));
    }

}



// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void Cl_FUNCTIONS::get_Cl_bins(vector<real_prec>Cl,  vector<real_prec>&Clb){

  cout<<BLUE<<endl;
  if(this->params._bin_type()=="linear"){
    for(ULONG i=0;i<lbin.size();i++)
      this->nmodes[i]=0;
    for(ULONG i=0;i<lbin.size();i++)
      Clb[i]=0;

    real_prec deltal=static_cast<real_prec>(this->params._Lmax()-this->params._Lmin())/static_cast<real_prec>(this->lbin.size());
    if(this->params._use_deltaL()==true)
          deltal=this->params._DeltaL();
    
    for(long l=this->params._Lmin();l<=this->params._Lmax();l++)
      {
        ULONG il=static_cast<int>(floor((l-this->params._Lmin())/deltal));
        if(il==lbin.size())il--;
        Clb[il]+=(l+0.5)*Cl[l];
        this->nmodes[il]+=(l+0.5);
      }
    for(ULONG i=0;i<lbin.size();++i)
        Clb[i]/=this->nmodes[i];
  }
  else
    {
      if(this->params._bin_type()=="log")
      {
        real_prec deltal=log10(params._Lmax()/1.0)/((real_prec)this->lbin.size());
        for(ULONG i=0;i<lbin.size();++i)this->nmodes[i]=0;
        for(long l=0;l<=params._Lmax();l++)
        {
          ULONG il=static_cast<ULONG>(floor((log10(l)-log10(1.0))/deltal));
          if(il==lbin.size())il--;
          Clb[il]+=(l+0.5)*Cl[l];
          this->nmodes[il]+=(l+0.5);
         }
        for(ULONG i=0;i<lbin.size();++i)
            Clb[i]/=this->nmodes[i];
      }
    }
  cout<<RESET;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_eCl_bins(vector<real_prec>Cl,  vector<real_prec>&Clb)
{

  cout<<BLUE<<endl;
  if(this->params._bin_type()=="linear")
  {
    for(ULONG i=0;i<(signed)lbin.size();i++)
      this->nmodes[i]=0;
    for(ULONG i=0;i<(signed)lbin.size();i++)
      Clb[i]=0;
    
    real_prec deltal=(static_cast<real_prec>(this->params._Lmax()-this->params._Lmin()))/static_cast<real_prec>(this->lbin.size());

    if(this->params._use_deltaL()==true)
      deltal=this->params._DeltaL();

    for(long l=this->params._Lmin();l<=this->params._Lmax();l++)
      {
        ULONG il=static_cast<int>(floor((l-this->params._Lmin())/deltal));
	
        if(il==lbin.size())il--;
        Clb[il]+=pow((l+0.5)*Cl[l],2);
        this->nmodes[il]+=(l+0.5);
      }
    for(ULONG i=0;i<lbin.size();i++)
        Clb[i]=sqrt(Clb[i])/this->nmodes[i];
  }
  else{
    if(this->params._bin_type()=="log")
      {
	real_prec deltal=log10(params._Lmax()/1.0)/(static_cast<real_prec>(this->lbin.size()));
	for(ULONG i=0;i<lbin.size();++i)
          this->nmodes[i]=0;
	for(long l=1;l<=params._Lmax();l++)
	  {
	    ULONG il=static_cast<ULONG>(floor((log10(l)-log10(1.0))/deltal));
	    if(il==lbin.size())
	      il--;
        Clb[il]+=pow(Cl[l],2);
	    this->nmodes[il]++;
	  }
    for(ULONG i=0;i<lbin.size();i++)
        Clb[i]=sqrt(Clb[i])/this->nmodes[i];
      }
  }
  cout<<RESET;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::W3J(){
  cout<<RED<<"Computing Wigner Symbols"<<endl;
  Wigner3J.resize(params._Lmax()+1);
  for(long i=0;i<params._Lmax()+1;i++){
    Wigner3J[i].resize(params._Lmax()+1);
    for(int j=0;j<params._Lmax()+1;j++){
      Wigner3J[i][j].resize(params._Lmax()+1);
    }
  }

#pragma omp parallel for collapse(2)
  for(long i=0;i<params._Lmax()+1;i++){
    for(ULONG j=0;j<params._Lmax()+1;j++){
      for(ULONG k=0;k<params._Lmax()+1;k++){
        ULONG J=i+j+k;
        if(sqrt(pow(i-j,2)) <= k &&  k<=i+j){
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

void Cl_FUNCTIONS::Mixing_matrix(){
  cout<<"Computing mixing matrix:"<<endl;

#pragma omp parallel for collapse(2)
  for(ULONG l1=params._Lmin();l1<params._Lmax();l1++)
    {
      for(ULONG l2=params._Lmin();l2<params._Lmax()+1;l2++)
	{
	  real_prec Rbis=0;
	  for(int l3=params._Lmin();l3<params._Lmax()+1;l3++)
	    {
	      if(abs(sqrt(pow(l1-l2,2))) <= l3 &&  l3<=l1+l2)
		{
		  real_prec WC = Wigner3J[l1][l2][l3];
		  Rbis+=(2.*l3+1)*this->Wl[l3]*(WC*WC);
		}
	    }
	  this->R[l1][l2]=((2.*l2+1)/(4.*M_PI))*Rbis;
	}
    }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************


void Cl_FUNCTIONS::Rll_bins(){


  // Only available for the Rll written in term of the Wigner symbols

  real_prec deltal;

  if(this->params._bin_type()=="linear"){

    for(ULONG i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;

    deltal=(this->params._Lmax()-this->params._Lmin())/((real_prec)lbin.size()); //here and just here I define N_L_bins as a real_prec instead of an integer

    if(this->params._use_deltaL()==true)
      deltal=this->params._DeltaL();

    
    for(ULONG l1=this->params._Lmin();l1<this->params._Lmax();++l1)
      {
	int lbina=floor((l1-this->params._Lmin())/deltal);
	
	if(lbina==lbin.size())lbina--;
	this->nmodes[lbina]+=(2.*l1+1);
	for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;++l2){
	  real_prec Rbis=0;
	  for(int l3=this->params._Lmin();l3<this->params._Lmax()+1;l3++){
	    if(abs(static_cast<real_prec>(l1-l2)) <= l3 &&  l3<=l1+l2){
	      real_prec WC = this->Wigner3J[l1][l2][l3];
	      Rbis+=(2*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
	    }
	  }
	  this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
	}
      }
    for(ULONG i=0;i<lbin.size();i++)for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++)Rll_bin[i][l2]/=this->nmodes[i];   }
  else{
    if(this->params._bin_type()=="log"){
      
      deltal=log10(params._Lmax()/1.0)/((real_prec)lbin.size());
      
      for(ULONG i=0;i<lbin.size();i++)this->nmodes[i]=0;

      for(ULONG l1=this->params._Lmin();l1<this->params._Lmax();l1++){
  	int lbina=floor((log(l1)-log(1.0))/deltal);
	this->nmodes[lbina]+=(2.*l1+1);
        for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++){
  	  real_prec Rbis=0;
      for(int l3=this->params._Lmin();l3<this->params._Lmax()+1;l3++){
            if(abs(static_cast<real_prec>(l1-l2)) <= l3 &&  l3<=l1+l2){
  	      real_prec WC = this->Wigner3J[l1][l2][l3];
              Rbis+=(2.*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	    }
  	  }
          this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
  	}
      }

      for(ULONG i=0;i<lbin.size();i++)for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++)Rll_bin[i][l2]/=this->nmodes[i];
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_mixing_matrix(){

  // compute the wigner symbols
  this->W3J();
  
  this->Mixing_matrix();
  real_prec factor_fs = (this->params._mixing_matrix_exact()? 1.0:this->catalog._sky_fraction());
  for(ULONG l1=this->params._Lmin();l1<=this->params._Lmax();l1++)
    for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
      this->R[l1][l2]/=factor_fs;

  for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
    {
      vector<real_prec> slice (this->params._Lmax()+1,0);
      string rfile=this->params._output_file_window()+"_MixingMatrix_l_"+to_string(l2)+".txt";
      for(ULONG l1=this->params._Lmin();l1<this->params._Lmax();l1++)
        slice[l1]=this->R[l1][l2];
      this->Fmi.write_to_file(rfile, this->lvec, slice);
    }
  string rfile;
  if(this->params._mixing_matrix_exact())rfile=this->params._output_file_window()+"_MixingMatrix_exact_nside_"+to_string(this->params._nside())+".txt";
  else rfile=this->params._output_file_window()+"_MixingMatrix_nside_"+to_string(this->params._nside())+".txt";
  this->Fmi.write_to_file(rfile, this->lvec, this->lvec, this->R);

  this->Rll_bins();
  for(ULONG l1=0;l1<this->params._N_L_bins();l1++)
    for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
      this->Rll_bin[l1][l2]/=factor_fs;

  rfile=this->params._output_file_window()+"_MixingMatrix_lbins_"+to_string(this->params._N_L_bins())+"_nside_"+to_string(this->params._nside())+".txt";
  this->Fmi.write_to_file(rfile, this->lbin, this->lvec, this->Rll_bin);

  for(ULONG l1=0;l1<this->params._N_L_bins();l1++)
    {
      vector<real_prec> slicea (this->params._Lmax()+1,0);
      for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
          slicea[l2]=this->Rll_bin[l1][l2];
      rfile=this->params._output_file_window()+"_MixingMatrix_lbin_"+to_string(l1)+"_nside_"+to_string(this->params._nside())+".txt";
      this->Fmi.write_to_file(rfile, this->lvec,slicea);
  }

}





// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void Cl_FUNCTIONS::get_alm_gaussian(int seed, vector<real_prec>&cl_th, Alm<xcomplex <real_prec> >&alm){
  const gsl_rng_type * T;

  gsl_rng_env_setup();
  gsl_rng_default_seed=seed;
  T = gsl_rng_ranlux;
  r = gsl_rng_alloc (T);

  for(long l=0;l<=params._Lmax();++l){
    for(long m=0; m<=l;++m){
      if(m==0){
	alm(l,m).real(gsl_ran_gaussian(r, sqrt(cl_th[l])));
	alm(l,m).imag(0);
      }
      else{
	alm(l,m).real(gsl_ran_gaussian(r, sqrt(0.5*cl_th[l])));
        alm(l,m).imag(gsl_ran_gaussian(r, sqrt(0.5*cl_th[l])));
      }
    }
  }

}



// ######################################################################
// ######################################################################

void Cl_FUNCTIONS::get_pixel_window(){
  pixel_window.resize(this->params._Lmax()+1,0);
  for(ULONG i=0;i<this->params._Lmax()+1;++i)this->pixel_window[i]=1.0;
}

// ######################################################################
// ######################################################################


void Cl_FUNCTIONS::read_mixing_matrix_lbins(){
   vector< vector<real_prec> > pR;
   Fmd.read_file(this->params._input_mixing_matrix_lbins(),pR);
   ULONG ik=0;
   for(ULONG i=0;i<this->params._N_L_bins();i++)
     for(ULONG j=this->params._Lmin();j<this->params._Lmax()+1;++j){
       this->Rll_bin[i][j]=pR[ik][2];
       ++ik;
     }
   
  pR.clear();
  // Normalize mixing matrix
  for(long l=0;l<this->lbin.size();l++)
    {
      real_prec ch=0;
      for(long lp=this->params._Lmin();lp<this->params._Lmax()+1;lp++)
	ch+=this->Rll_bin[l][lp];
      for(long lp=this->params._Lmin();lp<this->params._Lmax()+1;lp++)
	this->Rll_bin[l][lp]/=ch;
    }
}

// ######################################################################
// ######################################################################

void Cl_FUNCTIONS::read_mixing_matrix_l(){
   vector< vector<real_prec> > pR;
   Fmd.read_file(this->params._input_mixing_matrix_l(),pR);

   int ik=0;
   for(ULONG i=0;i<=this->params._Lmax();i++)
     for(ULONG j=0;j<=this->params._Lmax();++j)
       {
	 this->R[i][j]=pR[ik][2];
	 ++ik;	 
       }
   pR.clear();
   // Normalize mixing matrix
}


// ######################################################################
// ######################################################################




void Cl_FUNCTIONS::get_lg_map(Healpix_Map<real_prec>map, Healpix_Map<real_prec>&lmap){

  // This returns a delta that is now log-normal distributed
  real_prec mean=0.0;
#pragma omp parallel for reduction(+:mean)
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    mean+=map[i];
  mean/=(real_prec(static_cast<ULONG>(map.Npix())));
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    map[i]-=mean;
  real_prec var=0.0;
#pragma omp parallel for reduction(+:var)
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    var+=pow(map[i]-mean,2);
  var/=(static_cast<ULONG>(map.Npix()));
  real_prec m=-0.5*log(1+var);
#pragma omp parallel for
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    lmap[i]=exp(map[i]+m)-1.;
}

// ######################################################################
// ######################################################################

void Cl_FUNCTIONS::get_cross_Cl()
{
  
  set_vectors();
  set_Lbins();
  
  FILE_MANAGER<real_prec>Fmd;
  vector<real_prec>mask;
  ULONG N_columns_mask=Fmd.read_file2(this->params._input_file_mask(), mask);
  ULONG N_pixels=static_cast<ULONG>(mask.size()/N_columns_mask);
  ULONG Nside=sqrt(N_pixels/12);
  this->params.set_nside(Nside);
  Healpix_Map<real_prec>mask_aux(ilog2(Nside), RING);
  for(ULONG i=0;i<N_pixels;i++)
    mask_aux[i]=mask[this->params._i_mask_flag()+N_columns_mask*i];
  
  arr<real_prec>weight_data(2*mask_aux.Nside());
  weight_data.fill(1.0);
  
  // ****** get Alm from first input***********************
  Alm<xcomplex <real_prec> > Alm_m1(this->params._Lmin(),this->params._Lmax());
  map2alm(mask_aux,Alm_m1, weight_data,false); //HealPix
  vector<real_prec>mask2;
  Healpix_Map<real_prec>mask_aux2(ilog2(Nside), RING);
  ULONG N_columns_mask2=Fmd.read_file2(this->params._input_file_mask2(), mask2);
  for(ULONG i=0;i<N_pixels;i++)mask_aux2[i]=mask2[this->params._i_mask_flag()+N_columns_mask2*i];

  // ****** get Alm from second input***********************
  Alm<xcomplex<real_prec> > Alm_m2(this->params._Lmin(),this->params._Lmax());
  map2alm(mask_aux2,Alm_m2, weight_data,false); //HealPix
  vector<real_prec> Cl(this->params._Lmax()+1,0);
  vector<ULONG> lvec(this->params._Lmax()+1,0);

  // ****** get Ilm ****************************************

  Healpix_Map<real_prec>binary_mask(ilog2(Nside), RING);
  vector<real_prec>mask_bin;
  ULONG N_columns_mask3=Fmd.read_file2(this->params._input_file_binary_mask(), mask_bin);
  for(ULONG i=0;i<N_pixels;i++)binary_mask[i]=mask_bin[this->params._i_mask_flag()+N_columns_mask3*i];

  Alm<xcomplex <real_prec> > Ilm(this->params._Lmin(),this->params._Lmax());
  Healpix_Map<real_prec>pmask(ilog2(Nside), RING);
  map2alm_iter(binary_mask,Ilm,NUM_ITER_MAP2ALM, weight_data);
  // *******************************************************
  
#pragma omp parallel for
  for(long l=this->params._Lmin(); l<=this->params._Lmax(); l++)
    {
      for(ULONG m=0;m<=l;m++)
	{
	  real_prec Br_z1=Alm_m1(l,m).real()-Ilm(l,m).real();
	  real_prec Bi_z1=Alm_m1(l,m).imag()-Ilm(l,m).imag();
	  real_prec Br_z2=Alm_m2(l,m).real()-Ilm(l,m).real();
	  real_prec Bi_z2=Alm_m2(l,m).imag()-Ilm(l,m).imag();
	  Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
	  if(m==0)
	    Cl[l]/=2.;
	}
#pragma omp atomic
      Cl[l]/=(l+0.5);
      lvec[l]=static_cast<real_prec>(l);
    }

  vector<real_prec> Clbin(this->params._N_L_bins(),0);
  this->get_Cl_bins(Cl,Clbin);
  string ofile=this->params._output_file()+"Crossed_power_masks.txt";
  this->Fmd.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,Clbin);

}













