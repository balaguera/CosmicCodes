
/** @file AngularPowerSpectrumF.cpp
 *
 * @brief This file contains the class AngularPowerSpectrum
 * @details Generates measurements of angular power spectrum
 * @author Andres Balaguera Antolinez
 */

# include "../include/AngularPowerSpectrumF.hpp"

using namespace std;


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

void AngularPowerSpectrum::measure_angular_power()
{
  
  string temp="lbins";
  string temp_raw="raw";

  Catalogue tracer(this->params, "TRACER");

  fs::path filetr = fs::path(params._Input_dir_cat()) / params._file_catalogue(); 
  tracer.read_catalog_new(filetr.string()); // Read the asciii file still

  HaloTools htools(params, tracer);

  htools.set_type_of_object("TRACER");

  this->params.set_output_file_power_bins(this->params._Output_directory()+this->params._statistics()+"_"+this->params._type_of_angular_power_estimator()+"_"+this->params._Name_survey()+"_zbtype_"+this->params._define_z_bins()+"_"+temp);
  this->params.set_output_file_power(this->params._Output_directory()+this->params._statistics()+"_"+this->params._type_of_angular_power_estimator()+"_"+this->params._Name_survey()+"_zbtype_"+this->params._define_z_bins()+"_"+temp_raw);
  
  htools.read_mask();    // Read the mask

  this->set_healpix_pars();  // Compute Healpix numbers

  this->set_vectors();      // Allocate memory for used vectors

  if(this->params._n_z_bins_tomography()!= NEGATIVE_INT ) // If we have more than one redshift bin, proceed
    {

#ifdef _FULL_VERBOSE_
      cout<<GREEN<<"SELECTION INFO: "<< this->params._type_of_sample()<<RESET<<endl;
      cout<<CYAN;

      if(this->params._type_of_sample()=="fls")
      {
        So.message_screen("Flux limited sample selected. Bins in redshift");

        if(this->params._define_z_bins()=="delta")
            So.message_screen("zBins with constant width");

        if(this->params._define_z_bins()=="number")
           So.message_screen("zBins chosen with equal number of galaxies");
      }
#endif

// Get this->params._zmax()
      healpix_real zmax_cat=get<static_cast<int>(PropStats::MAX)>(tracer.info_redshift)  ;
      healpix_real zmin_cat=get<static_cast<int>(PropStats::MIN)>(tracer.info_redshift)  ;


      htools.set_min_redshift(zmin_cat);
      htools.set_max_redshift(zmax_cat);

      #ifdef _FULL_VERBOSE_
       So.message_screen("REDSHIFT INFO ");
       So.message_screen("Maximum redshift found in catalogue = ", zmax_cat);
       So.message_screen("Minimum redshift found in catalogue = ", zmin_cat);
#endif
      
      if(!this->params._use_min_max_redshift_from_cat())
        {
           So.message_screen("Maximim redshift selected in parameter file = ", this->params._max_redshift_tomography());
           So.message_screen("Minimum redshift selected in parameter file = ", this->params._min_redshift_tomography());
        }
      else
        {
          this->params.set_min_redshift_tomography(zmin_cat);
          this->params.set_max_redshift_tomography(zmax_cat);
        }
      
      // Get Magnitude limits if the MK is used
      if(this->params._i_abs_mag_g()>0)
        {
#ifdef _FULL_VERBOSE_
           So.message_screen("MAGNITUDE Mk INFO ");
           So.message_screen("Maximim Mk found in catalogue = ",htools._max_MK());
           So.message_screen("Minimum Mk found in catalogue = ", htools._min_MK());
#endif
          if(!this->params._use_min_max_absMag_from_cat())
          {
#ifdef _FULL_VERBOSE_
             So.message_screen("Maximim Mk selected in parameter file = ", this->params._MKmax());
             So.message_screen("Minimum Mk selected in parameter file = ", this->params._MKmin());
#endif
          }
        else
          {
              this->params.set_MKmin(htools._min_MK());
              this->params.set_MKmax(htools._max_MK());
          }

     // SOME WARNINGS:
	  std::cout<<RED;
#ifdef _FULL_VERBOSE_
      if(this->params._MKmax() > htools._max_MK())
          So.message_screen("Warning: Maximum Mk in parameter file (", this->params._MKmax(), ") greater than maximum value found in the catalog", htools._max_MK(), ")");
      if(this->params._MKmin() < htools._min_MK())
           So.message_screen("Warning: Minimum Mk in parameter file (", this->params._MKmin(),") smaller than minimum value found in the catalog (",htools._min_MK() ,")");
#endif
 	}
      
     if(this->params._max_redshift_tomography() > htools._max_redshift())
      	{
#ifdef _FULL_VERBOSE_
        cout<<"Warning: Maximum z in parameter file ("<<this->params._min_redshift_tomography()<<") greater than maximum value found in the catalog ("<<htools._max_redshift()<<"). Doing nothing"<<endl;
#endif
	  //   this->this->params._zmax()= this->this->params._zmax()_cat;
    	}
      if(this->params._min_redshift_tomography() < htools._min_redshift())
	    {
#ifdef _FULL_VERBOSE_
        cout<<"Warning: Minimum z in parameter file ("<<this->params._min_redshift_tomography()<<") smaller than minimum value found in the catalog ("<<htools._min_redshift()<<"). Doing nothing"<<endl;
#endif
  	}
   } // end if (n_zbins!=-1)
  else
    {
      this->params.set_min_redshift_tomography(-1000);
      this->params.set_max_redshift_tomography(1000);
    }
  
  int NB=0;

  if(this->params._n_z_bins_tomography()>=1)
   {
      NB =this->params._n_z_bins_tomography();
      htools.set_zbins();
  }
  
  // ******************************************************************************************************************
  // *********************************************************
  // OPERATIONS IN HARMONIC SPACE
  // *********************************************************Ilm
  // *********************************************************
   So.message_screen("Operations in harmonic space");

  Alm<xcomplex <healpix_real> > Ilm(this->params._Lmax(),this->params._Lmax());
  // Compute Jlm, Ilm. These objects are inizialized to the
  // full sky- case, Jlm=1, Ilm=0 (l>0)
  // such that if not computed in Map2ILM_Jlm, they remain full sky.
  // but if used, these have to be initialized to zero in that funcition
  // For full sky, Ilm os not set to zero by directly not used in the estimation
  this->Jlm.resize(this->params._Lmax()+1);

  for (auto& vec : this->Jlm)
    vec.resize(this->params._Lmax() + 1, 1.0);

   this->Map2Ilm_Jlm(htools,Ilm);

  int izf=0;
  if(this->params._n_z_bins_tomography()==-1)
    izf=0;

  izf=(this->params._n_z_bins_tomography()==1? 0 : this->params._n_z_bins_tomography());

  this->set_vectors_mean_ngal_pix();

  this->set_Lbins();

  // *********************REDSHIFT BINS**********************************************************

  // START LOOP OVER THE REDSHIFT BINS ONLY IF FLS
   So.message_screen("Computing Alm in z bins");

  this->set_BLMZ();
  

   for(ULONG aIZ=0;aIZ<=izf;++aIZ)
    {
      this->set_index_tomographic_zbin(aIZ);

      htools.set_redshift_bin(aIZ);

      Healpix_Map<healpix_real>healpix_map(ilog2(this->params._nside()), RING);

      this->n_pixels=healpix_map.Npix(); 

      bool data = true;
      htools.Cat2Map(healpix_map,data);

      
      this->shot_noise=htools._area_survey()/(static_cast<healpix_real>(htools._weighted_ngals_in_zbin()));

      this->Shot_Noise[aIZ]=this->shot_noise;// this already has the effect of the weights in the SN

      this->write_pars_cl_estimator(htools);

      this->Mean_ngal_pix[aIZ]=htools._mean_number_galaxies_pix();

      this->Mean_ngal[aIZ]=htools._weighted_ngals_in_zbin()/(htools._area_survey());

      Alm<xcomplex <healpix_real> > Alm(this->params._Lmax(),this->params._Lmax());

      arr<healpix_real>weight_data(2*healpix_map.Nside());

      weight_data.fill(1.0);

      So.message_screen("Computing Alm for redshift bin ", aIZ);
      map2alm(healpix_map,Alm, weight_data,false); //HealPix
      So.DONE();

      // Fill the vectors Blm. These are the Alm in different redshift or Magnitude bins
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG il=0;il<this->params._Lmax()+1;il++)
      {
        for(long im=0;im< static_cast<long>(this->params._Lmax()+1);im++)
          {
            this->Blm[il][im][aIZ].real(Alm(il,im).real());
            this->Blm[il][im][aIZ].imag(Alm(il,im).imag());
          }
      }

    }

     
  // ***********************************************************************************************
  // Here we can estimate the shot noise as shuffled realziations of the map, and compute its power,
  // since the shot_noise for the cross is not needed
  // This is equivalent as to have a random catalog with the same angular and redhsft distribution
  
  struct zsn{
    vector<healpix_real>Clvec_sn;
  };

  vector<zsn> shot_noise_ell(izf+1);

  if(this->params._SN_correction())
  {
    if(this->params._shot_noise_correction_randomized())
      {
        //Here we do a number of realizations randomizing the map, computing the Cl
        //allocating the results in the container Cl_sn.
        So.message_screen("Computing shot noise using ", this->params._N_random_shot_noise()," realizations");
        
        for(ULONG aIZ=0;aIZ<=izf;++aIZ)//Start again redshift bins to compute shot noise
          {

            this->set_index_tomographic_zbin(aIZ);

            htools.set_redshift_bin(aIZ);

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
                Healpix_Map<healpix_real>healpix_map(ilog2(this->params._nside()), RING);

                for(ULONG i=0;i<tracer._NOBJS();++i)
                {
                    ULONG npix=0;
                    bool inmask=false;
                    while(inmask==false)
                     {
                        point_r.phi=(360.0*gsl_rng_uniform(rn))*CFACTOR;
                        point_r.theta=acos(-1+2.0*gsl_rng_uniform(rn));
                        npix=static_cast<ULONG>(healpix_map.ang2pix(point_r));
                        if(1==htools.pixmask[npix])
                            inmask=true;
                     }
                    tracer.set_galPIXEL(npix,i);
                }

                //Get map. We use again healpix_map as its original version has been saved in Blm

                htools.Cat2Map(healpix_map,false);

                //Get Alm
                arr<healpix_real>weight_data(2*healpix_map.Nside());
                weight_data.fill(1.0);

                Alm<xcomplex <healpix_real> > Alm(this->params._Lmax(),this->params._Lmax());

                map2alm(healpix_map, Alm, weight_data,false);

                // Compute power
                vector<healpix_real>Cl_temp(this->params._Lmax()+1,0);

                this->Alm2Cl_sn(htools,aIZ,Alm,Ilm,Cl_temp);

                const auto invNran = num_1 / static_cast<healpix_real>(Nran);
                for (size_t il = 0; il < Cl_temp.size(); ++il)
                  shot_noise_ell[aIZ].Clvec_sn[il] += Cl_temp[il] * invNran;
            }
         }
         So.DONE();
      }  
    }

   htools.catalogue.clear_mem();


  // Having allocated quantities in z bin (fls)
  // we now do a healpix_real loop (fls) over the zbins, or
  // we do only one loop (vls) (the second is currently running)
  // and compute things for the current Mbins

  if(this->params._n_z_bins_tomography()>1)
    {

      So.message_screen("Computing auto Cls");

      for(int aIZ=1;aIZ<=izf; aIZ++) // Loop over the first bin
        {
          vector<healpix_real> ClzbinI(this->params._Lmax()+1,0);

          this->Alm2Cl(htools, aIZ,aIZ, Ilm, ClzbinI);//Get auto Cl at bin IZ

          if(this->params._shot_noise_correction_randomized() && this->params._SN_correction())
            for(size_t il=0 ; il<ClzbinI.size(); ++il)
               ClzbinI[il]-=shot_noise_ell[aIZ].Clvec_sn[il];

          vector<healpix_real> eClzbinI(this->params._Lmax()+1,0);
  
          this->get_eCl(htools, aIZ,aIZ,ClzbinI,ClzbinI,eClzbinI);
  
          string rest="_"+this->params._harmonic_sampling()+"_Nside_"+to_string(this->params._nside())+"_zbins_"+to_string(aIZ)+"_"+to_string(aIZ)+"_"+this->params._define_hemisphere()+"_"+this->params._coordinate_system()+".txt";
  
          string cfile=this->params._file_power()+rest;
          vector<vector<healpix_real>> cols = {this->lvec, ClzbinI, this->Wl, eClzbinI}; 

          So.message_screen("Writing raw power");
          this->File.write_to_file(cfile, cols);

          fill(this->Clbin.begin(),this->Clbin.end(),0);
  
          this->get_Cl_bins(ClzbinI, this->Clbin);
  
          fill(this->eClbin.begin(),this->eClbin.end(),0);
  
          this->get_eCl_bins(eClzbinI,this->eClbin);
  
          cfile=this->params._file_power_lbins()+rest;
  
          So.message_screen("Writing binned power");
          vector<vector<healpix_real>> colsb = {this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin}; 
          this->File.write_to_file(cfile, colsb);

          this->params.output_files_hgaps.push_back(cfile);
          
          if(this->params._generate_fits_files())
            {//This is  atrick; fts maps are written by the methods of the class Catalog, but we save files names here in oprder to read them form hgaps.cpp as aprt of Cl.params method.
              string mfile = this->params._output_file_fits()+"_zbin_"+to_string(this->index_tomographic_zbin)+"_"+this->params._statistics()+".fits";
              this->set_index_tomographic_zbin(aIZ);
              this->params.output_files_hgaps_maps.push_back(mfile);
            }

          if(this->params._measure_cross())
          {
            So.message_screen("Computing Cross Cls");

            for(int aJZ=aIZ;aJZ<=izf; aJZ++) // Loop over the second bin
            {
                vector<healpix_real> Cl(this->params._Lmax()+1,0); // COntainer for cross power
                vector<healpix_real> eCl(this->params._Lmax()+1,0); // Container for errors of power
                vector<healpix_real> ClzbinJ(this->params._Lmax()+1,0); // Container for auto power in bin J
                // ***********************************************************************************************
                // We need to compute the autos within each loop in order to compute th
                // variance for the cross power spectrum (C_i+sn_j)*(C_j+sn_j).
                // NOTE THAT FIRST WE NEED TO COMPUTE THE AUTO. DO not change ordering here.
                // ***********************************************************************************************
                if(aJZ==aIZ)//Auto Cl in bin J:
                {
                  this->Alm2Cl(htools, aJZ,aJZ,Ilm,ClzbinJ);
                  if(true==this->params._shot_noise_correction_randomized())
                      for(int il=0;il<shot_noise_ell[aIZ].Clvec_sn.size();++il)
                        ClzbinJ[il]-=shot_noise_ell[aIZ].Clvec_sn[il];
                  }

                this->Alm2Cl(htools, aIZ,aJZ,Ilm,Cl); //Cross

                this->get_eCl(htools,aIZ,aJZ,ClzbinI, ClzbinJ, eCl); //get variance of cross power:

                // Write raw estimates
                string rest;

                if(this->params._harmonic_sampling()=="H")
                  rest="_"+this->params._harmonic_sampling()+"_Nside_"+to_string(this->params._nside())+"_zbins_"+to_string(aIZ)+"_"+to_string(aJZ)+"_"+this->params._define_hemisphere()+"_"+this->params._coordinate_system()+".txt";
                else if(this->params._harmonic_sampling()=="DS")
                  rest="_"+this->params._harmonic_sampling()+"_zbins_"+to_string(aIZ)+"_"+to_string(aJZ)+"_"+this->params._define_hemisphere()+"_"+this->params._coordinate_system()+".txt";
        
                string cfile=this->params._file_power()+rest;
                this->File.write_to_file(cfile, this->lvec, Cl, this->Wl,eCl);

                this->get_Cl_bins(Cl, this->Clbin);  // Compute Cl in Bins of l


                if(aIZ==aJZ) // For auto power get lmax before SN domination
                  {
                    int ll=-1;
                    do{ll++;}
                    while(this->Clbin[ll] > this->Shot_Noise[aJZ]);
                    if(ll==this->lbin.size())
                      So.message_screen("No shot noise domination reached for estimates in redshift bin ", aJZ, ". Select lmax = 100. ");
                    else
                      So.message_screen("Maxium scale before shot-noise regime for Auto Cl in zBin ", aJZ, " is lmax = ", this->lbin[ll]);
                  }

                this->get_eCl_bins(eCl,this->eClbin);   /* Comute the Variance in Bins of L */

                this->get_Cl_bins(this->Wl, this->Wlbin);  /*  Get power of mask in Bins */

               // Write Binned estimates
                string ofile=this->params._file_power()+rest;
                this->File.write_to_file(ofile, this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin);
              }
          } // enf if cross is available
      }
	  
      So.DONE();

      {
        So.message_screen("Computing power for full sample");
        // Get the Cl of the full sample
        
        vector<healpix_real> Cl(this->params._Lmax()+1,0);
        
        vector<healpix_real> eCl(this->params._Lmax()+1,0);
        
        this->Alm2Cl(htools,0,0, Ilm,Cl);
        
        this->get_eCl(htools,0,0,Cl,Cl,eCl);

        string rest="_"+this->params._harmonic_sampling()+"_Nside_"+to_string(this->params._nside())+"_zbins_"+to_string(0)+"_"+to_string(0)+"_"+this->params._define_hemisphere()+"_"+this->params._coordinate_system()+".txt";

        vector<vector<healpix_real>> cols = {lvec, Cl, Wl, eCl};
        this->File.write_to_file(this->params._file_power()+rest,cols);

        this->get_Cl_bins(Cl,  this->Clbin);

        this->get_eCl_bins(eCl, this->eClbin);

        this->get_Cl_bins(this->Wl, this->Wlbin);
        {
          int ll=-1;
          do{ll++;}while(this->Clbin[ll]>this->Shot_Noise[0]);
          if(ll==this->lbin.size())std::cout<<RED<<"No shot noise domination reached for estimates in redshift bin "<<0<<". Select lmax = 100. "<<RESET<<std::endl;
          else std::cout<<RED<<"Maxium scale before shot-noise regime for Auto Cl in zBin "<< 0<<" is lmax = "<<(int)this->lbin[ll]<<RESET<<std::endl;
        }

        vector<vector<healpix_real>> colsb = {this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin};
        this->File.write_to_file(this->params._file_power_lbins()+rest,colsb);

      }
    }      
    else if(this->params._n_z_bins_tomography()<=1)// Get the Cl when no redshift info is available
    {   
 
      vector<healpix_real> Cl(this->params._Lmax()+1,0);
 
      vector<healpix_real> eCl(this->params._Lmax()+1,0);
 
      this->Alm2Cl(htools,0,0, Ilm,Cl);
 
      this->get_eCl(htools,0,0,Cl,Cl,eCl);

      string rest="_"+this->params._harmonic_sampling()+"_Nside_"+to_string(this->params._nside())+"_"+this->params._define_hemisphere()+"_"+this->params._coordinate_system()+".txt";

      string cfile=this->params._file_power()+rest;
      
      this->File.write_to_file(cfile, this->lvec, Cl, this->Wl, eCl);
      
      this->get_Cl_bins(Cl,  this->Clbin);

      this->get_eCl_bins(eCl, this->eClbin);

      this->get_Cl_bins(this->Wl, this->Wlbin);

      string ofile=this->params._file_power()+rest;
      vector<vector<healpix_real>> colsb = {this->lbin, this->lbin_min, this->lbin_max,this->Clbin,this->eClbin};
      this->File.write_to_file(ofile, colsb);
   

    }

    So.DONE();

  // Finally, get the mixig matrix
  if(true==this->params._compute_mixing_matrix())
    this->get_mixing_matrix(htools);
      
  return;
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *****************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerSpectrum::set_BLMZ(){
  // Resize Blm. I have to do it in an independend method
  this->Blm.resize(this->params._Lmax()+1);

  for(ULONG i=0; i<Blm.size();++i)
   this->Blm[i].resize(this->params._Lmax()+1);
 
   for(ULONG i=0; i<this->Blm.size();++i)
      for(ULONG j=0;j<this->Blm[i].size();++j)
          this->Blm[i][j].resize(this->params._n_z_bins_tomography()+1,0);

}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::set_vectors()
{
  this->So.enter(__PRETTY_FUNCTION__);
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

  for(size_t i=0;i<this->R.size();i++)
    this->R[i].resize(this->params._Lmax()+1,0);

  // THIS IS ALSO USED WHTN WE READ THE MIXING MATRIX!!!

  if(this->params._compute_mixing_matrix())
    {
      this->R.resize(this->params._Lmax()+1);
      for(size_t i=0;i<this->R.size();i++)
        this->R[i].resize(this->params._Lmax()+1,0);
      Rll_bin.resize(this->params._N_L_bins());
      for(size_t i=0;i<this->Rll_bin.size();i++)
        Rll_bin[i].resize(this->params._Lmax()+1,0);
    }
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerSpectrum::set_Lbins(){

  healpix_real deltal=0;
  if(this->params._type_of_binning()==LINEAR)
    {
      if(this->params._use_deltaL())
        deltal=this->params._deltaL();
      else
        deltal=(static_cast<healpix_real>(this->params._Lmax()-this->params._Lmin()))/(static_cast<healpix_real>(lbin.size()));
        //here and just here I define N_L_bins as a healpix_real instead of an integer
      for(size_t i=0;i<lbin.size();i++)
        this->lbin[i]=this->params._Lmin()+(i+0.5)*deltal;

      for(size_t i=0;i<lbin.size();i++)
        this->lbin_min[i]=this->params._Lmin()+i*deltal;

      for(size_t i=0;i<lbin.size();i++)
        this->lbin_max[i]=this->params._Lmin()+(i+1)*deltal;
  }
  else
    {
      if(this->params._type_of_binning()==LOG)
    	{
          healpix_real deltal=log10(params._Lmax()/1.0)/(static_cast<healpix_real>(lbin.size()));
          for(size_t i=0;i<lbin.size();i++)
            lbin[i]=pow(10,log10(1)+(i+0.5)*deltal);
          for(size_t i=0;i<lbin.size();i++)
            lbin_min[i]=pow(10,log10(1)+i*deltal);
          for(size_t i=0;i<lbin.size();i++)
            lbin_max[i]=pow(10,log10(1)+(i+1)*deltal);
      }
    }

#ifdef _FULL_VERBOSE_
      So.message_screen("L-bins INFO");
      So.message_screen("Bin-averaged power spectrum");
      So.message_screen("Number of bins = ", this->lbin.size());
      So.message_screen("l-Bin width = ", deltal);
#endif

  }

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::set_healpix_pars(){
  this->n_pixels=12*this->params._nside()*this->params._nside();
  this->Nrings=4*this->params._nside()-1;
}


// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerSpectrum::get_eCl(HaloTools &htools, int iz, int jz, vector<healpix_real>&cli,vector<healpix_real>&clj, vector<healpix_real>&ecl){
  /// WARNING: THE CROSS POWER SPECTRUM HAS UNDESTIMATED GAUSSIAN ERROR BARS, CHECN WHITE, SUNG, PERCIVAL
  
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto l=this->params._Lmin();l<=this->params._Lmax();++l)
    {
      healpix_real cl1=cli[l] + this->Shot_Noise[iz];
      healpix_real cl2=clj[l] + this->Shot_Noise[jz];
      healpix_real factor=2./(htools._sky_fraction()*(2*l+1));
      ecl[l]=sqrt(factor*cl1*cl2);
    }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::set_vectors_mean_ngal_pix(){

  int nb= this->params._n_z_bins_tomography()==-1 ? 1 : this->params._n_z_bins_tomography();
  this->Mean_ngal_pix.resize(nb+1,0);
  this->Mean_ngal.resize(nb+1,0);
  this->Shot_Noise.resize(nb+1,0);
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::write_pars_cl_estimator(HaloTools &htools){
  So.message_screen("Information computed from the catalog: ");
  if(this->params._n_z_bins_tomography()!=-1)
  {
    So.message_screen("Redshift bin ", this->index_tomographic_zbin);
    So.message_screen("z min = ",htools._z_min(this->index_tomographic_zbin));
    So.message_screen("z max = ",htools._z_max(this->index_tomographic_zbin));
  }

  So.message_screen("Number of galaxies in current redshift bin = ",htools._ngals_in_zbin());
  So.message_screen("Check: Number of galaxies in current redshift bin (from map) = ",htools._ngals_in_zbin_from_map());
  So.message_screen("Weighted (Gauss) number of galaxies in catalogue = ",htools._weighted_ngals_in_zbin());
  So.message_screen("Number of pixels in the mask = ",htools._n_pixels());
  So.message_screen("Number of pixels in the surveyed area = ",htools._n_observed_pixels());
  So.message_screen("Sky fraction = ",htools._sky_fraction());
  this->sky_fraction=htools._sky_fraction();
  So.message_screen("Area survey = ",htools._area_survey());
  So.message_screen("Area pixel = ",htools._area_pixel());
  So.message_screen("Poisson Shot-noise = ",this->shot_noise);
  So.message_screen("Mean number of galaxies in pixels = ",htools._mean_number_galaxies_pix());
  So.message_screen("Mean number of galaxies = ", htools._ngals_in_zbin() / htools._area_survey());
  So.message_screen("RMS = ",htools._rms_ngal());
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::Map2Alm(HaloTools &htools, Healpix_Map<healpix_real>map, Alm<xcomplex <healpix_real> >&alm,  Alm<xcomplex <healpix_real> >&Ilm, arr<arr<healpix_real> >&Jlm){

  int iring=0;

  for(long i=0;i<params._Lmax()+1;i++)
      for(long im=0;im<params._Lmax()+1;im++)
          alm(i,im).real(0);
  for(long i=0;i<params._Lmax()+1;i++)
      for(long im=0;im<params._Lmax()+1;im++)
          alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++){
    healpix_real x=cos(htools.theta_new[ir]);
    int Npixels_ring=npix_ring(this->params._nside(), ir);
    for(long m=0;m<=params._Lmax();m++){
      gsl_real *Pl= new gsl_real[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++)
      {
        int ip=ipr+iring;
        if(htools.pixmask[ip]==1){
//   	  healpix_real cc=cos(healpix_real(m)*htools.phi_new[ip]);
//  	  healpix_real ss=sin(healpix_real(m)*htools.phi_new[ip]);
          int il=0;
          for(long l=m;l<=params._Lmax();l++)
          {
            healpix_real Plm=Pl[il];
            Jlm[l][m]  += htools._area_pixel()*pow(Plm,2);
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


void AngularPowerSpectrum::Map2Alm(HaloTools & htools, Healpix_Map<healpix_real>map, Alm<xcomplex <healpix_real> >&alm){

  // This function does the job of the map2alm pf Healpix in longer time

  time_t start;
  time (&start);
  int iring=0;

  for(auto i=0;i<params._Lmax()+1;i++)
    for(auto im=0;im<params._Lmax()+1;im++)
      alm(i,im).real(0);
  for(long i=0;i<params._Lmax()+1;i++)
    for(long im=0;im<params._Lmax()+1;im++)
      alm(i,im).imag(0);

  for(int ir=0;ir<Nrings;ir++)
  {
    healpix_real x=cos(htools.theta_new[ir]);
    int Npixels_ring=npix_ring(this->params._nside(),ir);

    healpix_real a1=0,a2=0;
    for(auto m=0;m<=params._Lmax();m++){
      gsl_real *Pl= new gsl_real[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(int ipr=0;ipr<Npixels_ring;ipr++)
      {
      	int ip=ipr+iring;
        if(htools.pixmask[ip]==1)
        {
          healpix_real cc=cos(healpix_real(m)*htools.phi_new[ip]);
          healpix_real ss=sin(healpix_real(m)*htools.phi_new[ip]);
          int il=0;
          for(int l=m;l<=params._Lmax();l++)
            {
              healpix_real Plm=Pl[il];
              healpix_real Yreal=Plm*cc;
              healpix_real Yimag=Plm*ss;
              a1+=htools._area_pixel()*Yreal*map[ip];
              alm(l,m).real(a1);
              a2+=-htools._area_pixel()*Yimag*map[ip];
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

void AngularPowerSpectrum::Map2Ilm_Jlm(HaloTools &htools, Alm<xcomplex <healpix_real> >&Ilm){

  time_t start;
  time (&start);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(auto i=0;i<this->params._Lmax()+1;i++)
      for(auto im=0;im<this->params._Lmax()+1;im++)
        {
          Ilm(i,im).real(0);
          Ilm(i,im).imag(0);
        }

  int iring=0;

  // In this part we get the Ilm using HEALPIX functions
  So.message_screen("Computing Ilm...");

  Healpix_Map<healpix_real>pmask(ilog2(this->params._nside()), RING);
  arr<healpix_real>weight_data(2*pmask.Nside(),1.0);

  for(ULONG i=0;i<n_pixels;i++)
      pmask[i]=static_cast<healpix_real>(htools.pixmask[i]);

  map2alm(pmask,Ilm,weight_data);
  So.DONE();

  // ********************************************************************

  if(this->params._type_of_angular_power_estimator()=="D")
    {
      So.message_screen("Computing Jlm");
      if(this->params._compute_jlm())
      {
       for(long l=this->params._Lmin(); l<this->params._Lmax();++l)for(int m=this->params._Lmin(); m<this->params._Lmax();++m)this->Jlm[l][m]=0;

       for(int ir=0;ir<Nrings;ir++){
         comp_time(start,Nrings,ir);
         healpix_real x=cos(htools.theta_new[ir]);
         int Npixels_ring=npix_ring(this->params._nside(), ir);
         for(auto m=0;m<=params._Lmax();m++){
           gsl_real *Pl= new gsl_real[this->params._Lmax()-m+1];
           gsl_sf_legendre_sphPlm_array(this->params._Lmax(),m,x,Pl);
           for(int ipr=0;ipr<Npixels_ring;ipr++){
             int ip=ipr+iring;
             if(htools.pixmask[ip]==1){
              ULONG il=0;
              for(auto l=m;l<=this->params._Lmax();l++){
                healpix_real Plm=Pl[il];
                this->Jlm[l][m]  += htools._area_pixel()*pow(Plm,2);
                ++il;
              }
            }
          }
          delete[] Pl;
        }
        iring+=Npixels_ring;
      }

      ofstream jlm_file;
      jlm_file.open(this->params._output_file_jlm().c_str());
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
      File.read_file(this->params._output_file_jlm(),propJ);
      int ik=-1;
      for(long l=this->params._Lmin();l<=this->params._Lmax();++l){
        for(long m=0;m<=this->params._Lmax();++m)
        {
          ik++;
          this->Jlm[l][m]=static_cast<healpix_real>(propJ[ik][2]);
        }
      }
      propJ.clear();
    }


  }

}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::Map2Jlm(HaloTools &htools, Healpix_Map<healpix_real>&map_ran){

  for(long i=0;i<params._Lmax()+1;i++)
    for(long im=0;im<params._Lmax()+1;im++)
      this->Jlm[i][im]=0;

  int iring=0;
  for(ULONG ir=0;ir<Nrings;ir++){
    healpix_real x=cos(htools.theta_new[ir]);
    ULONG Npixels_ring=npix_ring(this->params._nside(), ir);
    for(long m=0;m<=params._Lmax();m++){
      gsl_real *Pl= new gsl_real[params._Lmax()-m+1];
      gsl_sf_legendre_sphPlm_array(params._Lmax(),m,x,Pl);
      for(ULONG ipr=0;ipr<Npixels_ring;ipr++){
        ULONG ip=ipr+iring;
        healpix_real cc=cos(healpix_real(m)*htools.phi_new[ip]);
        healpix_real ss=sin(healpix_real(m)*htools.phi_new[ip]);
        ULONG il=0;
        for(long l=m;l<=params._Lmax();l++)
        {
          healpix_real Plm=Pl[il];
          healpix_real Yreal=Plm*cc;
          healpix_real Yimag=Plm*ss;
          this->Jlm[l][m]+= htools._area_pixel()*pow(Plm,2)*map_ran[ip];
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

void  AngularPowerSpectrum::Alm2Cl(HaloTools &htools, int iz, int jz,  Alm<xcomplex <healpix_real> >&Ilm, vector<healpix_real> &Cl)
{

  healpix_real sn=(iz==jz? this->Shot_Noise[iz]:0);

  healpix_real norm_fsky=num_1/sn; // The mean surface number density

  sn = (this->params._SN_correction()==true ? sn:0);

  if(true==this->params._shot_noise_correction_randomized())
      sn=0; // set it to zero as SN will be removed after the call of this method

  healpix_real FSky = this->params._type_of_angular_power_estimator()=="D" ? 1.0 : htools._sky_fraction();

  healpix_real Mean_ngal_pix_i=0;
  
  healpix_real Mean_ngal_pix_j=0;
 
  healpix_real Mean_ngal_pix_i_aux=1;
 
  healpix_real Mean_ngal_pix_j_aux=1;

  // OJo con esta parte que depende de si hemos noramlizado antes o no

  if(this->params._statistics()=="ARF")
    {
        // For ARF we need sto divide for the mean (weighted ) number of galaxies but we do not subtract the term Ilm (why?)
      Mean_ngal_pix_i=0;
      Mean_ngal_pix_j=0;
      Mean_ngal_pix_i_aux=this->Mean_ngal_pix[iz];
      Mean_ngal_pix_j_aux=this->Mean_ngal_pix[jz];
    }
    else if(this->params._statistics()=="Cl")
    {
      Mean_ngal_pix_i=this->Mean_ngal_pix[iz];
      Mean_ngal_pix_j=this->Mean_ngal_pix[jz];
      Mean_ngal_pix_i_aux=Mean_ngal_pix_i;
      Mean_ngal_pix_j_aux=Mean_ngal_pix_j;
    }
 
  healpix_real iMean_ngal_pix_i_aux=num_1/Mean_ngal_pix_i_aux;
  healpix_real iMean_ngal_pix_j_aux=num_1/Mean_ngal_pix_j_aux;

  #ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int l=0;l<=this->params._Lmax();l++)
    {
      this->Wl[l]=0;
      Cl[l]=0;
      lvec[l]=static_cast<healpix_real>(l);
    }
  
  if(this->params._type_of_sky_coverage() == "masked_sky")
    {
      if(this->params._harmonic_sampling()=="H")
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for(auto l=0;l<=this->params._Lmax();++l)
    	    {
	          for(auto m=0;m<=l;m++)
              {
                healpix_real jlm=this->params._type_of_angular_power_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
                healpix_real Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal_pix_i*Ilm(l,m).real())*iMean_ngal_pix_i_aux;
                healpix_real Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal_pix_i*Ilm(l,m).imag())*iMean_ngal_pix_i_aux;
                healpix_real Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal_pix_j*Ilm(l,m).real())*iMean_ngal_pix_j_aux;
                healpix_real Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal_pix_j*Ilm(l,m).imag())*iMean_ngal_pix_j_aux;
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
                    Cl[l]*=itwo;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
                    this->Wl[l]*=itwo;
		              }
            	}
           Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          this->Wl[l]/=(l+0.5);
	      }
	    }
    else{ //if direct sum. Here we are not sure about the implementation of the Ilm, which comes from a Masked. Should we use randoms?
	    for(int l=this->params._Lmin();l<=this->params._Lmax();l++)
        {
          for(long m=0;m<=l;m++)
            {
              healpix_real jlm=this->params._type_of_angular_power_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
              healpix_real Br_z1=(this->Blm[l][m][iz].real()-Mean_ngal_pix_i_aux*Ilm(l,m).real())*iMean_ngal_pix_i_aux;
              healpix_real Bi_z1=(this->Blm[l][m][iz].imag()-Mean_ngal_pix_i_aux*Ilm(l,m).imag())*iMean_ngal_pix_i_aux;
              healpix_real Br_z2=(this->Blm[l][m][jz].real()-Mean_ngal_pix_j_aux*Ilm(l,m).real())*iMean_ngal_pix_j_aux;
              healpix_real Bi_z2=(this->Blm[l][m][jz].imag()-Mean_ngal_pix_j_aux*Ilm(l,m).imag())*iMean_ngal_pix_j_aux;
              Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2)/jlm;
              this->Wl[l]+=norm(Ilm(l,m));
              if(m==0)
                {
                  Cl[l]*=itwo;
                  this->Wl[l]*=itwo;
                }
            }
          Cl[l]=Cl[l]/(FSky*(l+0.5))-sn;
          this->Wl[l]/=(l+0.5);
      	}
      }
    }
   
  else if(this->params._type_of_sky_coverage() == "full_sky")
   {  // If full sky. In this case, the estimator D and K are equivalent, Jlm=1, Ilm=0
    if(this->params._harmonic_sampling()=="H")
     {
       for(int l=0;l<=this->params._Lmax();l++)
          {
            for(long m=0;m<=l;m++)
              {
                healpix_real Br_z1=this->Blm[l][m][iz].real()*iMean_ngal_pix_i_aux;
                healpix_real Bi_z1=this->Blm[l][m][iz].imag()*iMean_ngal_pix_i_aux;
                healpix_real Br_z2=this->Blm[l][m][jz].real()*iMean_ngal_pix_j_aux;
                healpix_real Bi_z2=this->Blm[l][m][jz].imag()*iMean_ngal_pix_j_aux;
                Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
                this->Wl[l]+=norm(Ilm(l,m));
                if(m==0)Cl[l]*=itwo;
                if(m==0) this->Wl[l]*=itwo;
            	}
          Cl[l]=(Cl[l]/(l+0.5))-sn;  // FOR D
          this->Wl[l]=Wl[l]/(l+0.5);
         }
      }

    else
      { //if Direct summation && full sky
        for(int l=0;l<=this->params._Lmax();l++)
          {
           for(long m=0;m<=l;m++)
            {
              healpix_real Br_z1=this->Blm[l][m][iz].real();
              healpix_real Bi_z1=this->Blm[l][m][iz].imag();
              healpix_real Br_z2=this->Blm[l][m][jz].real();
              healpix_real Bi_z2=this->Blm[l][m][jz].imag();
              Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
              this->Wl[l]+=norm(Ilm(l,m));
              if(m==0)
                {
                  Cl[l]*=itwo;
                  this->Wl[l]*=itwo;
                }
            }
            Cl[l]=Cl[l]/(l+0.5)/pow(norm_fsky,2)-sn;  // FOR D
            this->Wl[l]/=(l+0.5);  // Applied only to the window function correctly?
          }
        }
    }
}

// ====================================================================================
// ====================================================================================
// ====================================================================================
void  AngularPowerSpectrum::Alm2Cl_sn(HaloTools &htools, int iz, Alm<xcomplex <healpix_real> >&Alm_ran,Alm<xcomplex <healpix_real> >&Ilm, vector<healpix_real> &Cl){

  healpix_real FSky = this->params._type_of_angular_power_estimator()=="D" ? 1.0 : htools._sky_fraction();

  healpix_real Mean_ngal_pix_i=0;
  healpix_real Mean_ngal_pix_i_aux=1;
 
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

  healpix_real iMean_ngal_pix_i_aux=num_1/Mean_ngal_pix_i_aux;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int l=this->params._Lmin();l<=this->params._Lmax();l++)
    {
      for(long m=0;m<=l;m++)
        {
          healpix_real jlm=this->params._type_of_angular_power_estimator()=="D" ? this->Jlm[l][m] : 1.0 ;
          healpix_real Br=(Alm_ran(l,m).real()-Mean_ngal_pix_i*Ilm(l,m).real())*iMean_ngal_pix_i_aux;
          healpix_real Bi=(Alm_ran(l,m).imag()-Mean_ngal_pix_i*Ilm(l,m).imag())*iMean_ngal_pix_i_aux;
          Cl[l]+=(Br*Br+Bi*Bi)/jlm;
          if(m==0)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              Cl[l]*=itwo;
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

void AngularPowerSpectrum::get_Cl_bins(vector<healpix_real>&Cl,  vector<healpix_real>&Clb){

  if(this->params._type_of_binning()==LINEAR)
    {
      for(ULONG i=0;i<lbin.size();i++)
        this->nmodes[i]=0;
      for(ULONG i=0;i<lbin.size();i++)
        Clb[i]=0;

      healpix_real deltal=static_cast<healpix_real>(this->params._Lmax()-this->params._Lmin())/static_cast<healpix_real>(this->lbin.size());
      healpix_real ideltal = num_1/deltal;

      if(this->params._use_deltaL()==true)
          deltal=this->params._deltaL();
    
      for(long l=this->params._Lmin();l<=this->params._Lmax();l++)
        {
          ULONG il=static_cast<int>(floor((l-this->params._Lmin())*ideltal));
          if(il==lbin.size())il--;
          Clb[il]+=(l+0.5)*Cl[l];
          this->nmodes[il]+=(l+0.5);
        }
      for(ULONG i=0;i<lbin.size();++i)
        Clb[i]/=this->nmodes[i];
    }
    else
      {
        if(this->params._type_of_binning()==LOG)
          {

            healpix_real deltal=log10(this->params._Lmax()/this->params._Lmin())/((healpix_real)this->lbin.size());
            healpix_real ideltal = num_1/deltal;

            for(ULONG i=0;i<lbin.size();++i)
              this->nmodes[i]=0;

            for(long l=0;l<=params._Lmax();l++)
              {
                ULONG il=static_cast<ULONG>(floor((log10(l)-log10(params._Lmin()))*ideltal));
                if(il==lbin.size())
                  il--;
                Clb[il]+=(l+0.5)*Cl[l];
                this->nmodes[il]+=(l+0.5);
              }

            for(ULONG i=0;i<lbin.size();++i)
              Clb[i]/=this->nmodes[i];
        }
      }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::get_eCl_bins(vector<healpix_real>&Cl,  vector<healpix_real>&Clb)
{

  if(this->params._type_of_binning()==LINEAR)
  {
    for(ULONG i=0;i<(signed)lbin.size();i++)
      this->nmodes[i]=0;

    for(ULONG i=0;i<(signed)lbin.size();i++)
      Clb[i]=0;
    
    healpix_real deltal=(static_cast<healpix_real>(this->params._Lmax()-this->params._Lmin()))/static_cast<healpix_real>(this->lbin.size());
    healpix_real ideltal = num_1/deltal;


    if(this->params._use_deltaL()==true)
      deltal=this->params._deltaL();

    for(long l=this->params._Lmin();l<=this->params._Lmax();l++)
      {
        ULONG il=static_cast<int>(floor((l-this->params._Lmin())*ideltal));

        if(il==lbin.size())
          il--;

        Clb[il]+=pow((l+0.5)*Cl[l],2);
        this->nmodes[il]+=(l+0.5);
      }

    for(ULONG i=0;i<lbin.size();i++)
      Clb[i]=sqrt(Clb[i])/this->nmodes[i];
  }
  else if(this->params._type_of_binning()==LOG)
    {
      healpix_real deltal=log10(params._Lmax()/params._Lmin())/(static_cast<healpix_real>(this->lbin.size()));
      healpix_real ideltal = num_1/deltal;
  
      for(ULONG i=0;i<lbin.size();++i)
        this->nmodes[i]=0;
    
      for(long l=params._Lmin();l<=params._Lmax();l++)
        {
          ULONG il=static_cast<ULONG>(floor((log10(l)-log10(params._Lmin()))*ideltal));

          if(il==lbin.size())
            il--;

          Clb[il]+=pow(Cl[l],2);

          this->nmodes[il]++;
        }

        for(ULONG i=0;i<lbin.size();i++)
            Clb[i]=sqrt(Clb[i])/this->nmodes[i];
    }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::W3J(){
  So.message_screen("Computing Wigner Symbols");
  Wigner3J.resize(params._Lmax()+1);
  for(ULONG i=0;i<=params._Lmax();i++)
    {
      Wigner3J[i].resize(params._Lmax()+1);
      for(int j=0;j<=params._Lmax();j++)
        {
          Wigner3J[i][j].resize(params._Lmax()+1);
        }
    }

#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(ULONG i=0;i<params._Lmax()+1;i++)
    {
     for(ULONG j=0;j<params._Lmax()+1;j++)
      {
       for(ULONG k=0;k<params._Lmax()+1;k++)
        {
         ULONG J=i+j+k;
          if(sqrt(pow(i-j,2)) <= k &&  k<=i+j)
            {
            if(J%2==0)
             {
              long double num = factorial(J-2*i) * factorial(J-2*j) * factorial(J-2*k);
              long double pref =factorial(J/2-i)*factorial(J/2-j)*factorial(J/2-k) ;
              Wigner3J[i][j][k]= sqrt(num/factorial(J+1))*(factorial(J/2))/pref;
             }
            else
            {
        	    Wigner3J[i][j][k]=0;
            }
          }
      }
    }
  }
  So.DONE();
}
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************

void AngularPowerSpectrum::Mixing_matrix()
{

  So.message_screen("Computing mixing matrix:");
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(auto l1=params._Lmin();l1<params._Lmax();l1++)
    {
      for(auto l2=params._Lmin();l2<params._Lmax()+1;l2++)
	      {
          healpix_real Rbis=0;
          for(auto l3=params._Lmin();l3<params._Lmax()+1;l3++)
            {
              if(abs(sqrt(pow(l1-l2,2))) <= l3 &&  l3<=l1+l2)
                {
                  healpix_real WC = Wigner3J[l1][l2][l3];
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


void AngularPowerSpectrum::Rll_bins(){


  // Only available for the Rll written in term of the Wigner symbols

  healpix_real deltal=1.;
  healpix_real ifpi = num_1/static_cast<healpix_real>(4.*M_PI);

  if(this->params._type_of_binning()==LINEAR){

    for(ULONG i=0;i<(signed)lbin.size();i++)this->nmodes[i]=0;

    healpix_real deltal=(this->params._Lmax()-this->params._Lmin())/((healpix_real)lbin.size()); //here and just here I define N_L_bins as a healpix_real instead of an integer

    if(this->params._use_deltaL()==true)
      deltal=this->params._deltaL();

    
    for(ULONG l1=this->params._Lmin();l1<this->params._Lmax();++l1)
      {
	int lbina=floor((l1-this->params._Lmin())/deltal);
	
	if(lbina==lbin.size())lbina--;
	this->nmodes[lbina]+=(2.*l1+1);
	for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;++l2){
	  healpix_real Rbis=0;
	  for(int l3=this->params._Lmin();l3<this->params._Lmax()+1;l3++){
	    if(abs(static_cast<healpix_real>(l1-l2)) <= l3 &&  l3<=l1+l2){
	      healpix_real WC = this->Wigner3J[l1][l2][l3];
	      Rbis+=(2*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
	    }
	  }
	  this->Rll_bin[lbina][l2]+=((2.*l2+1)/(4.*M_PI))*Rbis;
    	}
   }
    for(ULONG i=0;i<lbin.size();i++)for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++)Rll_bin[i][l2]/=this->nmodes[i];   }
  else{
    if(this->params._type_of_binning()==LOG){
      
      deltal=log10(params._Lmax()/1.0)/(static_cast<healpix_real>(this->lbin.size()));
      
      for(ULONG i=0;i<lbin.size();i++)this->nmodes[i]=0;

      for(ULONG l1=this->params._Lmin();l1<this->params._Lmax();l1++){
  	int lbina=floor((log(l1)-log(1.0))/deltal);
	this->nmodes[lbina]+=(2.*l1+1);
        for(ULONG l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++){
  	  healpix_real Rbis=0;
      for(int l3=this->params._Lmin();l3<this->params._Lmax()+1;l3++){
            if(abs(static_cast<healpix_real>(l1-l2)) <= l3 &&  l3<=l1+l2){
  	      healpix_real WC = this->Wigner3J[l1][l2][l3];
              Rbis+=(2.*l1+1)*(2.*l3+1)*this->Wl[l3]*(WC*WC);
  	    }
  	  }
          this->Rll_bin[lbina][l2]+=((2.*l2+1)*ifpi)*Rbis;
  	}
      }

      for(size_t i=0;i<lbin.size();i++)
       for(auto l2=this->params._Lmin();l2<this->params._Lmax()+1;l2++)
        Rll_bin[i][l2]/=this->nmodes[i];
    }
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::get_mixing_matrix(HaloTools &htools){

  // compute the wigner symbols
  this->W3J();
  
  this->Mixing_matrix();

  healpix_real factor_fs = (this->params._mixing_matrix_exact()? 1.0:htools._sky_fraction());
  healpix_real ifactor_fs = num_1/factor_fs;

  for(auto l1=this->params._Lmin();l1<=this->params._Lmax();++l1)
    for(auto l2=this->params._Lmin();l2<=this->params._Lmax();++l2)
      this->R[l1][l2]*=ifactor_fs;

  for(auto l2=this->params._Lmin();l2<=this->params._Lmax();++l2)
    {
      vector<healpix_real> slice (this->params._Lmax()+1,0);
      string sl2=to_string(l2);
      fs::path rfile=fs::path(this->params._Output_directory()) / (this->params._file_window() + "_MixingMatrix_l_" + sl2);
      rfile.replace_extension("txt");
      for(auto l1=this->params._Lmin();l1<this->params._Lmax();l1++)
        slice[l1]=this->R[l1][l2];
      this->File.write_to_file(rfile.string(), this->lvec, slice);
    }

    fs::path rfile;
    string snside = to_string(this->params._nside());
  if(this->params._mixing_matrix_exact())
    rfile=fs::path(this->params._Output_directory()) / (this->params._file_window() + "_MixingMatrix_exact_nside_" + snside);
  else 
    rfile=fs::path(this->params._Output_directory()) / (this->params._file_window()  + "_MixingMatrix_nside_" + snside);
  rfile.replace_extension("txt");
  this->File.write_to_file(rfile.string(), this->lvec, this->lvec, this->R);

  this->Rll_bins();
  for(ULONG l1=0;l1<this->params._N_L_bins();l1++)
    for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
      this->Rll_bin[l1][l2]*=ifactor_fs;

  string snlbins = to_string(this->params._N_L_bins());
  rfile=fs::path(this->params._Output_directory()) / (this->params._file_window() + "_MixingMatrix_lbins_" + snlbins + "_nside_" + snside );
  rfile.replace_extension("txt");
  this->File.write_to_file(rfile.string(), this->lbin, this->lvec, this->Rll_bin);

  for(ULONG l1=0;l1<this->params._N_L_bins();l1++)
    {
      vector<healpix_real> slicea (this->params._Lmax()+1,0);
      for(ULONG l2=this->params._Lmin();l2<=this->params._Lmax();l2++)
          slicea[l2]=this->Rll_bin[l1][l2];
      rfile= fs::path(this->params._Output_directory()) / (this->params._file_window() + "_MixingMatrix_lbin_" + to_string(l1) + "_nside_" + to_string(this->params._nside()));
      rfile.replace_extension("txt");
      this->params.output_files_mixing_matrix.push_back(rfile);
      this->File.write_to_file(rfile, this->lvec,slicea);
  }
}

// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
// *******************************************************************************************************************************************************
void AngularPowerSpectrum::get_alm_gaussian(int seed, vector<healpix_real>&cl_th, Alm<xcomplex <healpix_real> >&alm){

  const gsl_rng_type * T;
  gsl_rng_env_setup();
  gsl_rng_default_seed=seed;
  T = gsl_rng_ranlux;
  gsl_rng *r = gsl_rng_alloc (T);

  for(auto l=0;l<=params._Lmax();++l){
    for(auto m=0; m<=l;++m){
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

void AngularPowerSpectrum::read_mixing_matrix_lbins(){
   vector< vector<real_prec> > pR;
   File.read_file(this->params._input_mixing_matrix_lbins(),pR);
   ULONG ik=0;
   for(auto i=0;i<this->params._N_L_bins();i++)
     for(auto j=this->params._Lmin();j<this->params._Lmax()+1;++j){
       this->Rll_bin[i][j]=static_cast<healpix_real>(pR[ik][2]);
       ++ik;
     }
   
  pR.clear();
  // Normalize mixing matrix
  for(long l=0;l<this->lbin.size();l++)
    {
      healpix_real ch=0;
      for(auto lp=this->params._Lmin();lp<this->params._Lmax()+1;lp++)
      	ch+=this->Rll_bin[l][lp];
     healpix_real ich=num_1/ch;
     
      for(auto lp=this->params._Lmin();lp<this->params._Lmax()+1;lp++)
      	this->Rll_bin[l][lp]*=ich;
    }
}

// ######################################################################
// ######################################################################

void AngularPowerSpectrum::read_mixing_matrix_l(){
   vector< vector<real_prec> > pR;
   File.read_file(this->params._input_mixing_matrix_l(),pR);

   int ik=0;
   for(ULONG i=0;i<=this->params._Lmax();i++)
     for(ULONG j=0;j<=this->params._Lmax();++j)
       {
        this->R[i][j]=static_cast<healpix_real>(pR[ik][2]);
        ++ik;	 
       }
   pR.clear();
   // Normalize mixing matrix
}


// ######################################################################
// ######################################################################




void AngularPowerSpectrum::get_lg_map(Healpix_Map<healpix_real>&map, Healpix_Map<healpix_real>&lmap){

  // This returns a delta that is now log-normal distributed
  healpix_real mean=0.0;
#ifdef _USE_OMP_
  #pragma omp parallel for reduction(+:mean)
#endif
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    mean+=map[i];
  mean/=(healpix_real(static_cast<ULONG>(map.Npix())));
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    map[i]-=mean;
  healpix_real var=0.0;

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:var)
#endif
  for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    var+=pow(map[i]-mean,2);

  var/=(static_cast<ULONG>(map.Npix()));
  healpix_real m=-0.5*log(1+var);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for(ULONG i=0;i<static_cast<ULONG>(map.Npix());i++)
    lmap[i]=exp(map[i]+m)-1.;
}

// ######################################################################
// ######################################################################

void AngularPowerSpectrum::get_cross_Cl()
{
  
  set_vectors();
  set_Lbins();
  
  vector<real_prec>mask_aux_r;
  ULONG N_columns_mask=this->File.read_file(this->params._input_file_mask(), mask_aux_r,1);
  ULONG N_pixels=static_cast<ULONG>(mask_aux_r.size()/N_columns_mask);
  ULONG Nside=sqrt(N_pixels/12);
  this->params.set_nside(Nside);
  Healpix_Map<healpix_real>mask_aux(ilog2(Nside), RING);

  for(ULONG i=0;i<N_pixels;i++)
    mask_aux[i]= static_cast<healpix_real>(mask_aux_r[this->params._i_mask_flag()+N_columns_mask*i]);

    mask_aux_r.clear(); mask_aux_r.shrink_to_fit();  
  arr<healpix_real>weight_data(2*mask_aux.Nside());
  weight_data.fill(1.0);
  
  // ****** get Alm from first input***********************
  Alm<xcomplex <healpix_real> > Alm_m1(this->params._Lmin(),this->params._Lmax());
  map2alm(mask_aux,Alm_m1, weight_data,false); //HealPix
  vector<real_prec>mask2;
  Healpix_Map<healpix_real>mask_aux2(ilog2(Nside), RING);
  ULONG N_columns_mask2=File.read_file(this->params._input_file_mask2(), mask2);
  for(ULONG i=0;i<N_pixels;i++)
    mask_aux2[i]=static_cast<healpix_real>(mask2[this->params._i_mask_flag()+N_columns_mask2*i]);

  // ****** get Alm from second input***********************
  Alm<xcomplex<healpix_real> > Alm_m2(this->params._Lmin(),this->params._Lmax());
  map2alm(mask_aux2,Alm_m2, weight_data,false); //HealPix
  vector<healpix_real> Cl(this->params._Lmax()+1,0);
  vector<ULONG> lvec(this->params._Lmax()+1,0);

  // ****** get Ilm ****************************************

  Healpix_Map<healpix_real>mask(ilog2(Nside), RING);
  vector<real_prec>mask_bin;
  ULONG N_columns_mask3=File.read_file(this->params._input_file_mask(), mask_bin,1);
  for(ULONG i=0;i<N_pixels;i++)
    mask[i]=static_cast<healpix_real>(mask_bin[this->params._i_mask_flag()+N_columns_mask3*i]);
  mask_bin.clear(); mask_bin.shrink_to_fit(); 
  Alm<xcomplex <healpix_real> > Ilm(this->params._Lmin(),this->params._Lmax());
  Healpix_Map<healpix_real>pmask(ilog2(Nside), RING);
  map2alm(mask,Ilm,weight_data);
  // *******************************************************
#ifdef _USE_OMP_  
#pragma omp parallel for
#endif
  for(auto l=this->params._Lmin(); l<=this->params._Lmax(); l++)
    {
      for(auto m=0;m<=l;m++)
      {
        healpix_real Br_z1=Alm_m1(l,m).real()-Ilm(l,m).real();
        healpix_real Bi_z1=Alm_m1(l,m).imag()-Ilm(l,m).imag();
        healpix_real Br_z2=Alm_m2(l,m).real()-Ilm(l,m).real();
        healpix_real Bi_z2=Alm_m2(l,m).imag()-Ilm(l,m).imag();
        Cl[l]+=(Br_z1*Br_z2+Bi_z1*Bi_z2);
        if(m==0)
          Cl[l]*=itwo;
      }
#ifdef _USE_OMP_  
#pragma omp atomic
#endif
      Cl[l]/=(l+0.5);
      lvec[l]=static_cast<healpix_real>(l);
    }

  vector<healpix_real> Clbin(this->params._N_L_bins(),0);
  this->get_Cl_bins(Cl,Clbin);
  fs::path ofile=fs::path(this->params._file_power()) / "Crossed_power_masks.txt";
  this->File.write_to_file(ofile.string(), this->lbin, this->lbin_min, this->lbin_max,Clbin);

}











