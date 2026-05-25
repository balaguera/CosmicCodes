////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<HaloTools>
 * @file HaloTools.cpp
 * @brief Methods of the class HaloTools
 * @details The class HaloTools reads and analyses an input catalog od dark matter tracers
 * @author Andres Balaguera Antolinez 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../include/def.hpp"
#include "../include/HaloTools.hpp"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void HaloTools::select_mark(std::vector<real_prec>&mark,std::string property)
{    

       if (property=="_MASS_")
#pragma prallel for
       for(size_t i=0;i<mark.size();++i)
           mark[i]=log10(this->catalogue.mass_at(i));
       else if(property=="_VMAX_")
#pragma prallel for
         for(size_t i=0;i<mark.size();++i)
           mark[i]=log10(this->catalogue.vmax_at(i));
       else if  (property=="_RS_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.rs_at(i);
       else if  (property=="_RVIR_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.rvir_at(i);
       else if  (property=="_CONCENTRATION_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=log10(this->catalogue.concentration_at(i));
       else if  (property=="_SPIN_")
#pragma prallel for
       for(size_t i=0;i<mark.size();++i)
           mark[i]=log10(this->catalogue.spin_at(i));
       else if  (property=="_SPIN_BULLOCK_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=log10(this->catalogue.spin_bullock_at(i));
       else if  (property=="_VRMS_")
#pragma prallel for
         for(size_t i=0;i<mark.size();++i)
           mark[i]=log10(this->catalogue.vrms_at(i));
       else if (property=="_VIRIAL_")
#pragma prallel for
         for(size_t i=0;i<mark.size();++i)
           mark[i]=this->catalogue.virial_at(i);
       else if  (property=="_BTOA_")
#pragma prallel for
         for(size_t i=0;i<mark.size();++i)
           mark[i]=this->catalogue.b_to_a_at(i);
       else if  (property=="_CTOA_")
#pragma prallel for
         for(size_t i=0;i<mark.size();++i)
           mark[i]=this->catalogue.c_to_a_at(i);
       else if  (property=="_MACH_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.mach_number_at(i);
       else if  (property=="_BIAS_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.bias_at(i);
       else if  (property=="_TIDAL_ANISOTROPY_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.tidal_anisotropy_at(i);
       else if  (property=="_LOCAL_OVERDENSITY_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.local_overdensity_at(i);
       else if  (property=="_PEAK_HEIGHT_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.peak_height_at(i);
       else if  (property=="_DACH_")
#pragma prallel for
     for(size_t i=0;i<mark.size();++i)
       mark[i]=this->catalogue.dach_number_at(i);
     }



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_cosmo(){
  So.enter(__PRETTY_FUNCTION__);
  time_t start_all;
  time(&start_all);
  time_t end;
  So.message_screen("Interpolating cosmological functions");
  int nz=this->params._Nbins_redshift();
  this->rcosmo.resize(nz,0);
  this->zcosmo.resize(nz,0);
  real_prec  zmin_inter=0.0;
  real_prec  zmax_inter=this->params._redshift_max_sample();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<nz;++i)
    this->zcosmo[i]=zmin_inter+i*(zmax_inter-zmin_inter)/static_cast<double>(nz-1);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<nz;++i)
    this->rcosmo[i]=static_cast<gsl_real>(this->Cosmo.comoving_distance(this->zcosmo[i]));
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::analyze_cat(bool read){

  this->So.enter(__PRETTY_FUNCTION__);

  if(read)
    {
/*
      #if defined (_USE_ALL_PK_) || defined (_USE_MASS_CUTS_PK_)
      this->read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),pow(10,this->params._LOGMASSmin()));
#else
      this->read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),pow(10,this->params._LOGMASSmin()),pow(10,this->params._LOGMASSmax()));
#endif
*/

Catalogue catalogue(params, params._type_of_object());
      catalogue.read_catalog_new(this->params._Input_dir_cat()+this->params._file_catalogue());
  }
  //* COMPUTE LOCAL MACH NUMBER


  if(this->params._Get_cell_local_mach_number())
    this->get_local_mach_number(true);
  //* COMPUTE DISTRIBUTION OF SEPARATIONS
  // There is a bug related to memmory in jumping from get_neighbpuro to get_distribution. UNidentified.
  if(this->params._get_distribution_min_separations())
    {
      vector<s_nearest_cells> ncells_info;
      ncells_info.resize(this->params._NGRID());
      So.message_screen("Identifying neighbour cells");
      get_neighbour_cells_cat_analyze(this->params._Nft(), N_CELLS_BACK_FORTH, ncells_info, false);
      So.DONE();
      this->get_distribution_min_separations(ncells_info);
    }
  //* COMPUTE ABUNDANCE:
  // We must set sos ifs in order to aks what are we going to analyze from the catalog
  string fname_mass_function_Y = this->Output_directory+this->catalogue._type_of_object()+this->params._Name_survey()+"_abundance_R"+to_string(this->params._realization())+".txt";
  if(this->params._Get_prop_function_tracer() ||  this->params._Get_tracer_mean_number_density())
    this->get_property_function(fname_mass_function_Y);
  //* COMPUTE ABUNDANCE IN CWT:
  if(this->params._Get_prop_function_tracer_cwt())
    this->get_property_function_cosmic_web_types(fname_mass_function_Y);
 
  // ---------------------------------------------------------------------------- //
  // If requested, we assign an interpolation of the mean number density for each tracer given its mass.
  if(this->params._Get_tracer_mean_number_density())
  {
    vector<double> mfunc(this->mass_function.size(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0; i< mfunc.size(); ++i)
      mfunc[i]=this->mass_function[i]*pow(this->params._Lbox(),3);
  
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, mfunc.size());
    gsl_spline_init (spline, &(MBin[0]), &(mfunc[0]), this->mass_function.size());
    real_prec lm_min=this->params._LOGMASSmin();
    real_prec lm_max=this->params._LOGMASSmax();
    for(size_t i=0; i< this->catalogue._NOBJS(); ++i)
      this->catalogue.mean_number_density_at(i) = static_cast<real_prec>(gsl_spline_eval (spline, this->catalogue.mass_at(i), acc));
  }
 
 
  // ---------------------------------------------------------------------------- //
  // these are wainting for an input parameter
  /*
    string fname_mass_function_Y = this->Output_directory+this->catalogue._type_of_object()+this->params._Name_survey()+"_HOD_R"+to_string(this->params._realization())+".txt";
    this->get_HOD(fname_mass_function_Y);
    this->get_HOD_web_types(fname_mass_function_Y);
  */
  //this->get_pdf_vmax("VMAX");
  //this->get_sep_distribution(PROP_THRESHOLD_MULTI_SCALE_4);
  // ---------------------------------------------------------------------------- //
  //* COMPUTE NUMBER COUNTS IN AS MESH
  if(this->params._Get_tracer_number_counts())
    {
      string file=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft())+"_MAS"+to_string(this->box.masskernel);
      this->get_density_field_grid(_COUNTS_,file);
    }
  /*
    this->gp_abundance<<"set log x \n";
    this->gp_abundance<<"set border linewidth 2.2\n";
    this->gp_abundance<<"set xlabel 'Tracer Mass [Ms / h]' font 'Times-Roman,12'\n";
    this->gp_abundance<<"set ylabel 'Occupation Number in cells'  font 'Times-Roman,12'\n";
    vector<pair<real_prec, real_prec> > xy_pts_m;
    for(size_t i=0; i<this->catalogue._NOBJS(); ++i)
    xy_pts_m.push_back(std::make_pair(this->catalogue.mass_at(i), ncounts[this->catalogue.GridID]));
    this->gp_abundance<<"plot"<<this->gp_abundance.file1d(xy_pts_m) << " w p ps 0.2 pt 7 title 'Mass'"<<endl;
    xy_pts_m.clear(); xy_pts_m.shrink_to_fit();
  */
  // ---------------------------------------------------------------------------- //
  //* COMPUTE MASS WEIGHTED FIELD
  if(this->params._Get_tracer_mass_field())
    if(this->params._i_mass_g()>0 && this->params._i_mass_g()<this->NCOLS)
      {
        string file=this->params._Output_directory()+this->params._Name_survey()+"_MASS_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
        this->get_density_field_grid(_MASS_,file);
      }
  // ---------------------------------------------------------------------------- //
  //* COMPUTE VMAX WEIGHTED FIELD
  if(this->params._Get_tracer_vmax_field())
    if(this->params._i_vmax_g()>0 && this->params._i_vmax_g()<this->NCOLS)
      {
	string file=this->params._Output_directory()+this->params._Name_survey()+"_VMAX_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
	this->get_density_field_grid(_VMAX_,file);
      }
  // ---------------------------------------------------------------------------- //
  //* COMPUTE SPIN WEIGHTED FIELD
  if(this->params._Get_tracer_spin_field())
    if(this->params._i_spin_g()>0 && this->params._i_spin_g()<this->NCOLS)
      {
        string file=this->params._Output_directory()+this->params._Name_survey()+"_SPIN_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
        this->get_density_field_grid(_SPIN_,file);
      }
  // ---------------------------------------------------------------------------- //
  if(this->params._Get_tracer_bias_field())
    {
      vector<real_prec> DM_DEN_FIELD( this->params._NGRID(),0);
      this->File.read_array(this->params._Input_Directory_X()+this->params._Name_Catalog_X(),DM_DEN_FIELD);
      object_by_object_bias(this->params, this->catalogue, DM_DEN_FIELD);      // Esto es en Miscelanious.cpp
      string file=this->params._Output_directory()+this->params._Name_survey()+"_BIAS_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
      this->get_density_field_grid(_BIAS_,file);
    }

  // ---------------------------------------------------------------------------- //

  // these are waiting for an input parameter
  /*
    if(this->params._i_rs_g()>0 && this->params._i_rs_g()<this->NCOLS)
    {
    file=this->params._Output_directory()+"RS_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
    this->get_density_field_grid(_RS_,file);
    }

    if(this->params._i_virial_g()>0 && this->params._i_virial_g()<this->NCOLS)
    {
    file=this->params._Output_directory()+"VIR_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
    this->get_density_field_grid(_VIRIAL_,file);
    }
  */
  // these are wainting for an input parameter
  /*
  //These lines read a halo catalog and compute the cwt and local dm if the same DM -IC is available
  this->cwclass.set_params(this->params);
  string file_X=this->params._Input_Directory_X()+this->params._Name_Catalog_X();
  vector<real_prec>deltam(this->params._NGRID(),0);
  this->File.read_array(file_X,deltam);
  get_overdens(deltam, deltam);
  this->cwclass.get_CWC(deltam);
  #pragma omp parallel for
  for(size_t ig=0; ig<this->catalogue._NOBJS();++ig)
  {
  //        cout<<this->catalogue[ig].GridID<<endl;
  this->catalogue[ig].gal_cwt=this->cwclass.cwt_used[this->cwclass.get_Tclassification(this->catalogue[ig].GridID)]; // thisis written in the cats,
  this->catalogue.local_dm_at(ig)=deltam[this->catalogue[ig].GridID]; // thisis written in the cats,
  }
  string outputFileName=this->params._Output_directory()+"CAT_R"+to_string(this->params._realization())+"_"+this->params._Name_survey()+".txt_cwc";
  this->write_catalog(outputFileName);
  // ******************************************************************************************************
  //** This section is waiting for an input parameter:
  // These lines compute the Pearson coefficient between differnet halo properties in cosmic web types
  // and in in bis of DM
  this->cwclass.set_params(this->params);
  string file_X=this->params._Input_Directory_X()+this->params._Name_Catalog_X();
  vector<real_prec>deltam(this->params._NGRID(),0);
  this->File.read_array(file_X,deltam);
  get_overdens(deltam, deltam);
  real_prec dm_min=-1;// 0;
  real_prec dm_max=2;//0;
  this->cwclass.get_CWC(deltam);
  size_t NX=10;
  real_prec delta=(dm_max-dm_min)/static_cast<real_prec>(NX);
  size_t Ncwt=5;
  struct dm_cwt{
  vector<real_prec>mass;
  vector<real_prec>vel;
  vector<real_prec>den;
  vector<real_prec>idm;
  vector<int>ntr;
  };
  vector<dm_cwt> props(Ncwt);
  vector<dm_cwt> mean_props(Ncwt);
  vector<dm_cwt> std_props(Ncwt);
  vector<dm_cwt> corr_props(Ncwt);
  for(size_t ig=0; ig<this->catalogue._NOBJS();++ig)
  {
  size_t  IGRID=this->catalogue[ig].GridID;
  int ICWT=this->cwclass.cwt_used[this->cwclass.get_Tclassification(IGRID)];
  real_prec delta_dm=1+deltam[IGRID];
  int ddm_bin=get_bin(log10(delta_dm),dm_min,NX,delta, true);
  props[ICWT].mass.push_back(this->catalogue.mass_at(ig));
  props[ICWT].vel.push_back(this->catalogue.vmax_at(ig));
  props[ICWT].den.push_back(delta_dm);
  props[ICWT].idm.push_back(ddm_bin);
  props[0].mass.push_back(this->catalogue.mass_at(ig));
  props[0].vel.push_back(this->catalogue.vmax_at(ig));
  props[0].den.push_back(delta_dm);
  props[0].idm.push_back(ddm_bin);
  }
  for (int j=0;j<Ncwt;++j)
  {
  mean_props[j].mass.resize(NX,0);
  mean_props[j].vel.resize(NX,0);
  mean_props[j].den.resize(NX,0);
  mean_props[j].ntr.resize(NX,0);
  for(size_t i=0;i<props[j].mass.size();++i)
  {
  int ddm_bin=props[j].idm[i];
  mean_props[j].mass[ddm_bin]+=props[j].mass_at(i);
  mean_props[j].vel[ddm_bin]+=props[j].vel[i];
  mean_props[j].den[ddm_bin]+=props[j].den[i];
  mean_props[j].ntr[ddm_bin]++;
  }
  }
  for (int j=0;j<Ncwt;++j)
  {
  std_props[j].mass.resize(NX,0);
  std_props[j].vel.resize(NX,0);
  std_props[j].den.resize(NX,0);
  corr_props[j].mass.resize(NX,0);
  corr_props[j].vel.resize(NX,0);
  corr_props[j].den.resize(NX,0);
  for(size_t i=0;i<props[j].mass.size();++i)
  {
  int ddm_bin=props[j].idm[i];
  real_prec nnn=static_cast<real_prec>(mean_props[j].ntr[ddm_bin]);
  real_prec nnnv=static_cast<real_prec>(mean_props[j].ntr[ddm_bin]-1.);
  if(nnn>0)
  {
  real_prec mm=(props[j].mass_at(i)-mean_props[j].mass[ddm_bin]/nnn);
  real_prec vv=(props[j].vel[i]-mean_props[j].vel[ddm_bin]/nnn);
  real_prec dd=(props[j].den[i]-mean_props[j].den[ddm_bin]/nnn);
  std_props[j].mass[ddm_bin] +=mm*mm/nnnv;
  std_props[j].vel[ddm_bin]  +=vv*vv/nnnv;
  std_props[j].den[ddm_bin]  +=dd*dd/nnnv;
  corr_props[j].mass[ddm_bin]+=mm*vv/nnnv;
  corr_props[j].vel[ddm_bin] +=dd*vv/nnnv;
  corr_props[j].den[ddm_bin] +=dd*mm/nnnv;
  }
  }
  }
  for (int j=0;j<Ncwt;++j)
  for (int ddm_bin=0;ddm_bin<NX;++ddm_bin)
  {
  real_prec nnn=static_cast<real_prec>(mean_props[j].ntr[ddm_bin]-1);
  if(nnn>0){
  corr_props[j].mass[ddm_bin]/=sqrt(std_props[j].vel[ddm_bin]*std_props[j].mass[ddm_bin]);
  corr_props[j].vel[ddm_bin]/=sqrt(std_props[j].vel[ddm_bin]*std_props[j].den[ddm_bin]);
  corr_props[j].den[ddm_bin]/=sqrt(std_props[j].mass[ddm_bin]*std_props[j].den[ddm_bin]);
  }
  }
  ofstream sal;
  string outf_k=this->params._Output_directory()+"CorrelationMSIGMA_R"+to_string(this->params._realization())+".txt";
  sal.open(outf_k.c_str());
  So.message_screen("Writting correlation in bins of dm in file ", outf_k);
  for (size_t i=0;i<NX;++i){
  real_prec deltamm=dm_min+(i+0.5)*delta;
  sal<<deltamm<<"  "<<corr_props[0].mass_at(i)<<"   "<<corr_props[1].mass_at(i)<<"  "<<corr_props[2].mass_at(i)<<"   "<<corr_props[3].mass_at(i)<<"   "<<corr_props[4].mass_at(i)<<endl;
  cout<<deltamm<<"  "<<corr_props[0].mass_at(i)<<"   "<<corr_props[1].mass_at(i)<<"  "<<corr_props[2].mass_at(i)<<"   "<<corr_props[3].mass_at(i)<<"   "<<corr_props[4].mass_at(i)<<endl;
  }
  sal.close();
  outf_k=this->params._Output_directory()+"CorrelationDMSIGMA_R"+to_string(this->params._realization())+".txt";
  sal.open(outf_k.c_str());
  So.message_screen("Writting correlation in bins of dm in file ", outf_k);
  for (size_t i=0;i<NX;++i){
  real_prec deltamm=dm_min+(i+0.5)*delta;
  sal<<deltamm<<"  "<<corr_props[0].vel[i]<<"   "<<corr_props[1].vel[i]<<"  "<<corr_props[2].vel[i]<<"   "<<corr_props[3].vel[i]<<"   "<<corr_props[4].vel[i]<<endl;
  cout<<deltamm<<"  "<<corr_props[0].vel[i]<<"   "<<corr_props[1].vel[i]<<"  "<<corr_props[2].vel[i]<<"   "<<corr_props[3].vel[i]<<"   "<<corr_props[4].vel[i]<<endl;
  }
  sal.close();
  outf_k=this->params._Output_directory()+"CorrelationDMMAS_R"+to_string(this->params._realization())+".txt";
  sal.open(outf_k.c_str());
  So.message_screen("Writting correlation in bins of dm in file ", outf_k);
  for (size_t i=0;i<NX;++i){
  real_prec deltamm=dm_min+(i+0.5)*delta;
  sal<<deltamm<<"  "<<corr_props[0].den[i]<<"   "<<corr_props[1].den[i]<<"  "<<corr_props[2].den[i]<<"   "<<corr_props[3].den[i]<<"   "<<corr_props[4].den[i]<<endl;
  cout<<deltamm<<"  "<<corr_props[0].den[i]<<"   "<<corr_props[1].den[i]<<"  "<<corr_props[2].den[i]<<"   "<<corr_props[3].den[i]<<"   "<<corr_props[4].den[i]<<endl;
  }
  sal.close();
  */

  // The lines are meatn to determine tracer statisticla properties such as Lfunction, Mstar function
  // nethgods which demand the construction of the 1/Vmax estimtor etc
  // There will be also a request to build random catalogs as in Galxy.
  if(this->params._Get_Luminosity_function())
    {}
  if(this->params._Get_Mstellar_function())
    {}
  if(this->params._Get_Color_function())
    {}
  if(this->params._Get_Color_Mag_plane())
    {}
  if(this->params._Get_Random_Catalog())
    {
      this->get_random_catalog();
    }
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::write_catalog(string outputFileName)
{
  So.enter(__PRETTY_FUNCTION__);
  real_prec conversion_factor=1.0;
  int Nprop_file=MIN_N_PROP_CAT;
#if defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_ || defined _ASSIGN_VMAX_POST_
    Nprop_file++;
#ifdef _ASSIGN_MASS_POST_
    Nprop_file++;
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    Nprop_file++;
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
    Nprop_file++;
#endif
#endif
#endif
#ifdef _WRITE_BINARY_BAM_FORMAT_
    // The structure pof the binary file will be:
    //  size_t Nlines
    //  int Ncolumns = 10
    //   Nlines*Ncolumns floats x,y,z,vx,vy,vz,M,vmax,Rs,s
    ofstream outStream;
    this->So.message_screen("Writting to binary file", outputFileName);
    outStream.open(outputFileName.c_str(), ios::binary|ios::out);
    this->So.message_screen("Writting Number of tracers = ", this->catalogue._NOBJS());
    outStream.write(reinterpret_cast<char*>(&this->catalogue._NOBJS()), sizeof(size_t));
    this->So.DONE();
    this->So.message_screen("Writting number of columns = ", Nprop_file);
    outStream.write(reinterpret_cast<char*>(&Nprop_file), sizeof(size_t));
    this->So.DONE();
#ifdef _OUTPUT_WITH_HEADERS_
    vector<string> data_i(Nprop_file+1);
    data_i.push_back("#x(Mpc/h).");
    data_i.push_back("y(Mpc/h).");
    data_i.push_back("z(Mpc/h).");
    data_i.push_back("vx(km/s).");
    data_i.push_back("vy(km/s).");
    data_i.push_back("vz(km/s).");
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
#ifdef _ASSIGN_MASS_POST_
    data_i.push_back("Mass(Ms/h).");
    data_i.push_back("Vmax(km/s).");
#ifdef _USE_RS_AS_OBSERVABLE_POWER_
    data_i.push_back("Rs(kpc/h).");
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    data_i.push_back("Spin.");
#endif
#endif // end for assign_mass_post
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    data_i.push_back("Mass(Ms/h).");
#endif
    data_i.push_back("Ntracers");
    for(size_t i=0;i<Nprop_file+1;++i)
      {
    int sizes_i=data_i[i].size();
        outStream.write((char *)&sizes_i, sizeof(sizes_i));
        outStream.write(data_i[i].c_str(), sizes_i);
      }
#endif
    for(size_t i = 0; i< this->catalogue._NOBJS(); ++i)
      {
#ifdef _WRITE_COORDINATES_
    float x=static_cast<float>(this->catalogue.coord1_at(i));
    float y=static_cast<float>(this->catalogue.coord2_at(i));
    float z=static_cast<float>(this->catalogue.coord3_at(i));
    outStream.write(reinterpret_cast<char*>(&x), sizeof(x));
    outStream.write(reinterpret_cast<char*>(&y), sizeof(y));
    outStream.write(reinterpret_cast<char*>(&z), sizeof(z));
#endif
#ifdef _WRITE_VELOCITIES_
    float vx=static_cast<float>(this->catalogue.vel1_at(i)*conversion_factor);
    float vy=static_cast<float>(this->catalogue.vel2_at(i)*conversion_factor);
    float vz=static_cast<float>(this->catalogue.vel3_at(i)*conversion_factor);
    outStream.write(reinterpret_cast<char*>(&vx), sizeof(vx));
    outStream.write(reinterpret_cast<char*>(&vy), sizeof(vy));
    outStream.write(reinterpret_cast<char*>(&vz), sizeof(vz));
#endif
#if defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ || defined (_ASSIGN_MASS_POST_)
    float mass=static_cast<float>(this->catalogue.mass_at(i));
    outStream.write(reinterpret_cast<char*>(&mass), sizeof(mass));
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    float vmax=static_cast<float>(this->catalogue.vmax_at(i));
    outStream.write(reinterpret_cast<char*>(&vmax), sizeof(vmax));
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
    float rs=static_cast<float>(this->catalogue.rs_at(i));
        outStream.write(reinterpret_cast<char*>(&rs), sizeof(rs));
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    float spin=static_cast<float>(this->catalogue.spin_at(i));
    outStream.write(reinterpret_cast<char*>(&spin), sizeof(spin));
#endif
      }
#else  // else, write in ascii
    ofstream outStream;
    outStream.precision(6);
    outStream.setf(ios::showpoint);
    outStream.setf(ios::scientific);
    outStream.open(outputFileName.c_str(), ios::out);
     this->So.message_screen("Writting to ascii file", outputFileName);
    this->So.message_screen("with ", this->catalogue._NOBJS(), "objects");
#endif
    for(size_t i = 0; i< this->catalogue._NOBJS(); ++i)
    {
        if(this->catalogue.observed_at(i)==1)
#if defined _ASSIGN_PROPERTIES_
        outStream<<this->catalogue.coord1_at(i)<<"\t"<<this->catalogue.coord2<<"\t"<<this->catalogue.coord3<<"\t"<<this->catalogue.vel1*conversion_factor<<"\t"<<this->catalogue.vel2*conversion_factor<<"\t"<<this->catalogue.vel3*conversion_factor<<"\t"<<this->catalogue.mass<<"\t"<<this->catalogue.vmax<<"\t"<<this->catalogue.rs<<"\t"<<this->catalogue.spin<<"\t"<<this->catalogue.identity<<"\t"<<this->catalogue.gal_cwt<<"\t"<<this->catalogue.local_dm<<endl;
#elif defined _WRITE_COORDINATES_ || defined _WRITE_VELOCITIES_
        outStream<<this->catalogue.redshift_at(i)<<"\t"<<this->catalogue.mass_at(i)<<"\t"<<this->catalogue.color_at(i)<<"\t"<<this->catalogue.stellar_mass_at(i)<<"\t"<<this->catalogue.relative_bias_at(i)<<endl;
#endif
    }
    outStream.close();
#ifdef _VERBOSE_CAT_
    So.DONE();
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HaloTools::get_density_field_grid(string prop, string output_file)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _VERBOSE_CAT_
  So.message_screen("Interpolating on a grid using MAS= ", this->params._masskernel());
#endif
  vector<real_prec> deltaTR_counts(this->box.NGRID,0);
  MAS_NEW(&box,this->catalogue, _COUNTS_, deltaTR_counts);
  vector<real_prec> deltaTR;
  if(_COUNTS_ == prop)
    this->File.write_array(output_file, deltaTR_counts);
  else
    {
      deltaTR.resize(this->box.NGRID,0);
      MAS_NEW(&box,this->catalogue, prop, deltaTR);
      if(_MASS_ == prop || _RS_== prop || _VIRIAL_== prop || _SPIN_== prop || _VMAX_== prop || _BIAS_ == prop)
        {
      So.message_screen("Interpolating tracer on a grid weighting by", prop);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0 ; i< this->box.NGRID; ++i)
        if(deltaTR_counts[i]!=0)
              deltaTR[i]=static_cast<double>(deltaTR[i])/static_cast<double>(deltaTR_counts[i]); // Mean mass in cells
        else
          deltaTR[i]=0;
        }
      else if (_SAT_FRACTION_==prop)
    So.message_screen("Interpolating Number of satellites in a grid");
      this->File.write_array(output_file, deltaTR);
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HaloTools::get_density_field_grid(std::string prop, vector<real_prec>&deltaTR_counts,vector<real_prec>&out)
{
   vector<real_prec> deltaTR;
 deltaTR.resize(this->box.NGRID,0);
 So.message_screen("Interpolating tracer on a grid weighting");
 MAS_NEW(&box,this->catalogue, prop, deltaTR);
 for(size_t i=0 ; i< this->box.NGRID; ++i)
   if(deltaTR_counts[i]!=0)
    out[i]=static_cast<double>(deltaTR[i])/static_cast<double>(deltaTR_counts[i]); // Mean mass in cells
  else
   out[i]=0;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_density_field_grid(std::string prop, vector<real_prec>&out)
{
 So.message_screen("Interpolating tracer on a grid");
 MAS_NEW(&box,this->catalogue, prop, out);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::define_property_bins()
{
    So.enter(__PRETTY_FUNCTION__);
    cout<<__LINE__<<endl;
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  
   omp_set_num_threads(NTHREADS);
#endif
#ifdef _VERBOSE_CAT_
#ifdef MBINS
   So.message_screen("Defining bins for abundance in catalog type ", this->catalogue._type_of_object());
#elif defined MCUTS
  So.message_screen("Defining cuts for abundnace:");
#endif
#endif
  if(this->params._i_mass_g() >0)
    {
      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-mass");
#endif
      // to account for the zero bin which includes the full sample, except if one asks for only one bin, in which case the full is assumed
      real_prec MMin;
      real_prec MMax;
#ifdef _MASS_LOG_
      MMin = this->params._LOGMASSmin();
      MMax = this->params._LOGMASSmax();
#ifdef _VERBOSE_CAT_
      this->So.message_screen("Minimum Mass", pow(10,MMin), "Ms/h");
      this->So.message_screen("Maximum Mass", pow(10,MMax), "Ms/h");
#endif
#else
      MMin= pow(10,this->params._LOGMASSmin());
#ifdef _VERBOSE_CAT_
      this->So.message_screen("Minimum Mass", MMin, "Ms/h");
      this->So.message_screen("Maximum Mass", MMax, "Ms/h");
#endif
#endif
      this->logdeltaM=(MMax-MMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->MBmin.clear();this->MBmin.shrink_to_fit();
      this->MBmin.resize(this->NMBINS,0);
#ifdef _MASS_LOG_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->NMBINS;++i)
        this->MBmin[i]=pow(10,MMin+i*logdeltaM);
#else
      for(size_t i=0;i<this->NMBINS;++i)
            this->MBmin.push_back(MMin+i*logdeltaM);
#endif
#ifdef _MASS_LOG_
      this->MBmax.clear();this->MBmax.shrink_to_fit();
      this->MBin.clear();this->MBin.shrink_to_fit();
      this->MBin.resize(this->NMBINS,0);
      this->MBmax.resize(this->NMBINS,0);
#ifdef MBINS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->NMBINS;++i)
        {
          this->MBmax[i]=pow(10,MMin+(i+1)*logdeltaM);
          this->MBin[i]=this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf());
        }
#elif defined MCUTS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->NMBINS;++i)
        this->MBmax[i]=pow(10, MMax);
#endif
#else
#ifdef MBINS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NMBINS;++i)
        this->MBmax[i]=(MMin+(i)*logdeltaM);
#elif defined MCUTS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NMBINS;++i)
        this->MBmax[i]=(MMax);
#endif
#endif
#ifdef _VERBOSE_CAT_
      So.DONE();
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No halo mass info available for type ", this->catalogue._type_of_object());
#endif
      this->NMBINS=1;
      this->logdeltaM=1;
      this->MBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->MBmax.resize(NMBINS,100.0);  //
    }
#if defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_ || defined _ASSIGN_VMAX_POST_
  if(this->params._i_vmax_g() >0)
    {
      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-vmax");
#endif
      real_prec Min= log10(this->params._VMAXmin());
      real_prec Max= log10(this->params._VMAXmax());
#ifdef _VERBOSE_CAT_
      this->So.message_screen("Minimum VMAX", pow(10,Min), "km/s");
      this->So.message_screen("Maximum VMAX", pow(10,Max), "km/s");
#endif
      this->logdeltaVMAX=(Max-Min)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->VMAXBmin.clear();this->VMAXBmin.shrink_to_fit();this->VMAXBmin.resize(this->NMBINS,0);
      for(size_t i=0;i<this->NMBINS;++i)
        this->VMAXBmin[i]=pow(10,Min+i*logdeltaVMAX);
      this->VMAXBmax.clear();this->VMAXBmax.shrink_to_fit();this->VMAXBmax.resize(this->NMBINS,0);
      this->VMAXBin.clear();this->VMAXBin.shrink_to_fit();this->VMAXBin.resize(this->NMBINS,0);
#ifdef MBINS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->NMBINS;++i)
        {
          this->VMAXBmax[i]=pow(10,Min+(i+1)*logdeltaVMAX);
          this->VMAXBin[i]=(this->VMAXBmin[i]+(i+0.5)*(this->VMAXBmax[i]-this->VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));
        }
#elif defined MCUTS
      for(size_t i=0;i<this->NMBINS;++i)
        this->VMAXBmax.push_back(pow(10, MMax));
#endif
    }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No VMAX info available for type ", this->catalogue._type_of_object());
#endif
      this->NMBINS=1;
      this->logdeltaVMAX=1;
      this->VMAXBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->VMAXBmax.resize(NMBINS,100.0);  //
    }
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
  if(this->params._i_rs_g() >0)
    {
      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-vmax");
#endif
      real_prec RSMin;
      real_prec RSMax;
#ifdef _RS_LOG_
      RSMin= log10(this->params._RSmin());
      RSMax= log10(this->params._RSmax());
#ifdef _VERBOSE_CAT_
      this->So.message_screen("Minimum RS", pow(10,RSMin), "kpc/h");
      this->So.message_screen("Maximum RS", pow(10,RSMax), "kpc/h");
#endif
#else
      RSMin= pow(10,this->params._RSmin());
      RSMin= this->params._RSMASSmin();
#ifdef _VERBOSE_CAT_
      this->So.message_screen("Minimum RS", RSMin, "Ms/h");
      this->So.message_screen("Maximum RS", RSMax, "Ms/h");
#endif
#endif
      this->logdeltaRS=(RSMax-RSMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->RSBmin.clear();this->RSBmin.shrink_to_fit();
#ifdef _RS_LOG_
      for(size_t i=0;i<this->NMBINS;++i)
        this->RSBmin.push_back(pow(10,RSMin+i*this->logdeltaRS));
#else
      for(size_t i=0;i<this->NMBINS;++i)
        this->RSBmin.push_back(MMin+i*logdeltaVMAX);
#endif
#ifdef _RS_LOG_
      this->RSBmax.clear();this->RSBmax.shrink_to_fit();
      this->RSBin.clear();this->RSBin.shrink_to_fit();
#ifdef RSBINS
      for(size_t i=0;i<this->NMBINS;++i)
        this->RSBmax.push_back(pow(10,RSMin+(i+1)*logdeltaRS));
      for(size_t i=0;i<this->NMBINS;++i)
        this->RSBin.push_back(this->RSBmin[i]+(i+0.5)*(this->RSBmax[i]-this->RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));
#elif defined MCUTS
      for(size_t i=0;i<this->NMBINS;++i)
        this->RSBmax.push_back(pow(10, RSMax));
#endif
#else
#ifdef RSBINS
      for(size_t i=0;i<NMBINS;++i)
        this->RSBmax.push_back(RSMin+(i)*logdeltaRS);
#elif defined MCUTS
      for(size_t i=0;i<NMBINS;++i)
        this->RSmax.push_back(RSMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No RS info available for type ", this->catalogue._type_of_object());
#endif
      this->NMBINS=1;
      this->logdeltaRS=1;
      this->RSBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->RSBmax.resize(NMBINS,100.0);  //
    }
#endif
    //-----------------------Cvir
  bool use_cvir=false;
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
  use_cvir=true;
#endif
    if(true==use_cvir)
     {
      this->NMBINS = this->params._NMASSbins_mf() == 1 ? 1 : this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-cvir");
#endif
      real_prec SMin;
      real_prec SMax;
#ifdef _CONCENTRATION_LOG_
      SMin= log10(this->params._CONCENTRATIONmin());
      SMax= log10(this->params._CONCENTRATIONmax());
      this->So.message_screen("Minimum Cvir", pow(10,SMin));
      this->So.message_screen("Maximum Cvir", pow(10,SMax));
#endif
      this->logdeltaCONCENTRATION=(SMax-SMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->CONCENTRATIONBmin.clear();this->CONCENTRATIONBmin.shrink_to_fit();
#ifdef _CONCENTRATION_LOG_
      for(size_t i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmin.push_back(pow(10,SMin+i*this->logdeltaCONCENTRATION));
#else
      for(size_t i=0;i<this->NMBINS;++i)
        this->SPINBmin.push_back(SMin+i*logdeltaSPIN);
#endif
#ifdef _CONCENTRATION_LOG_
      this->CONCENTRATIONBmax.clear();this->CONCENTRATIONBmax.shrink_to_fit();
#ifdef CONCENTRATIONBINS
      for(size_t i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(pow(10,SMin+(i+1)*logdeltaCONCENTRATION));
      this->CONCENTRATIONBin.clear();this->CONCENTRATIONBin.shrink_to_fit();
      for(size_t i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBin.push_back(pow(10,SMin+(i+0.5)*this->logdeltaCONCENTRATION));
#elif defined CONCENTRATIONCUTS
      for(size_t i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(pow(10, SMax));
#endif
#else
#ifdef CONCENTRATIONBINS
      for(size_t i=0;i<NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(SMin+(i)*logdeltaSPIN);
#elif defined MCUTS
      for(size_t i=0;i<NMBINS;++i)
        this->CONCENTRATIONmax.push_back(SMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No Cvir info available for type ", this->catalogue._type_of_object());
#endif
      this->NMBINS=1;
      this->logdeltaCONCENTRATION=1;
      this->CONCENTRATIONBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->CONCENTRATIONBmax.resize(NMBINS,100.0);  //
    }
    //-----------------------spin
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
  if(this->params._i_spin_g() >0 ||this->params._i_spin_bullock_g() )
    {
      this->NMBINS = this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",this->NMBINS," bin(s) (zeroth -or  first bin- is the full sample) to tracer-spin");
#endif
      real_prec SMin;
      real_prec SMax;
#ifdef _SPIN_LOG_
      SMin= log10(this->params._SPINmin());
      SMax= log10(this->params._SPINmax());
      this->So.message_screen("Minimum Spin", pow(10,SMin));
      this->So.message_screen("Maximum Spin", pow(10,SMax));
#endif
      this->logdeltaSPIN=(SMax-SMin)/(static_cast<real_prec>(this->params._NMASSbins_mf()));
      this->SPINBmin.clear();this->SPINBmin.shrink_to_fit();
#ifdef _SPIN_LOG_
      for(size_t i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBmin.push_back(pow(10,SMin+i*this->logdeltaSPIN));
#else
      for(size_t i=0;i<this->NMBINS;++i)
        this->SPINBmin.push_back(SMin+i*logdeltaSPIN);
#endif
#ifdef _SPIN_LOG_
      this->SPINBmax.clear();this->SPINBmax.shrink_to_fit();
      this->SPINBin.clear();this->SPINBin.shrink_to_fit();
#ifdef SPINBINS
      for(size_t i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBmax.push_back(pow(10,SMin+(i+1)*this->logdeltaSPIN));
      for(size_t i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBin.push_back(pow(10,SMin+(i+0.5)*this->logdeltaSPIN));  
#elif defined SPINCUTS
      for(size_t i=0;i<this->NMBINS;++i)
        this->SPINBmax.push_back(pow(10, SMax));
#endif
#else
#ifdef SPINBINS
      for(size_t i=0;i<NMBINS;++i)
        this->SPINBmax.push_back(SMin+(i)*logdeltaSPIN);
#elif defined MCUTS
      for(size_t i=0;i<NMBINS;++i)
        this->SPINmax.push_back(SMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No SPIN info available for type ", this->catalogue._type_of_object());
#endif
      this->NMBINS=1;
      this->logdeltaSPIN=1;
      this->SPINBmin.resize(NMBINS,-100.0); //  if no mass is available, in the loop over the masses we set mass = 1,
      this->SPINBmax.resize(NMBINS,100.0);  //
    }
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::define_property_bins(size_t Nbins, string prop)
{
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif
    //this->NMBINS=this->params._NMASSbins_mf();
#ifdef _VERBOSE_CAT_
#ifdef MBINS
    So.message_screen("Defining bins for abundance:");

#elif defined MCUTS
  So.message_screen("Defining cuts for abundnace:");
#endif
#endif
  if(prop=="_MASS_")
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",Nbins," bin(s) (zeroth -or  first bin- is the full sample) to tracer-mass");
#endif
      // to account for the zero bin which includes the full sample, except if one asks for only one bin, in which case the full is assumed
      real_prec  MMin = (this->params._MASSbins_min(0));  // los this->params._MASSbins_min(0) ya vienen en log, ojo con eso; las demas propiedades no.
      real_prec  MMax = (this->params._MASSbins_max(this->params._NMASSbins_power()-1));
      this->logdeltaM=(MMax-MMin)/static_cast<real_prec>(Nbins);
      this->MBmin.clear();this->MBmin.shrink_to_fit();
      this->MBmin.resize(Nbins,0);
      this->MBmax.clear();this->MBmax.shrink_to_fit();
      this->MBmax.resize(Nbins,0);
      this->MBin.clear();this->MBin.shrink_to_fit();
      this->MBin.resize(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<Nbins;++i)
        {
          this->MBmin[i]=MMin+i*logdeltaM;
          this->MBmax[i]=MMin+(i+1)*logdeltaM;
          this->MBin[i]=this->MBmin[i]+(i+0.5)*logdeltaM;
        }
    }
  if(prop=="_VMAX_")
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating", Nbins," bin(s) (zeroth -or  first bin- is the full sample) to tracer-vmax");
#endif
      real_prec  MMin = log10(this->params._VMAXbins_min(0));
      real_prec  MMax = log10(this->params._VMAXbins_max(this->params._NVMAXbins_power()-1));
      this->logdeltaVMAX=(MMax-MMin)/(static_cast<real_prec>(Nbins));
      this->VMAXBmin.clear();this->VMAXBmin.shrink_to_fit();
      this->VMAXBmin.resize(Nbins,0);
      this->VMAXBmax.clear();this->VMAXBmax.shrink_to_fit();
      this->VMAXBmax.resize(Nbins,0);
      this->VMAXBin.clear();this->VMAXBin.shrink_to_fit();
      this->VMAXBin.resize(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<Nbins;++i)
        {
          this->VMAXBmin[i]=MMin+i*logdeltaVMAX;
          this->VMAXBmax[i]=MMin+(i+1)*logdeltaVMAX;
          this->VMAXBin[i]=this->VMAXBmin[i]+(i+0.5)*logdeltaVMAX;
        }
    }
  if(prop=="_CONCENTRATION_")
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",Nbins," bin(s) (zeroth -or  first bin- is the full sample) to tracer-vmax");
#endif
      real_prec  MMin = log10(this->params._CONCENTRATIONbins_min(0));
      real_prec  MMax = log10(this->params._CONCENTRATIONbins_max(this->params._NCONCENTRATIONbins_power()-1));
      this->logdeltaCONCENTRATION=(MMax-MMin)/(static_cast<real_prec>(Nbins));
      this->CONCENTRATIONBmin.clear();this->CONCENTRATIONBmin.shrink_to_fit();
      this->CONCENTRATIONBmin.resize(Nbins,0);
      this->CONCENTRATIONBmax.clear();this->CONCENTRATIONBmax.shrink_to_fit();
      this->CONCENTRATIONBmax.resize(Nbins,0);
      this->CONCENTRATIONBin.clear();this->CONCENTRATIONBin.shrink_to_fit();
      this->CONCENTRATIONBin.resize(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i< Nbins;++i)
        {
          this->CONCENTRATIONBmin[i]=MMin+i*logdeltaCONCENTRATION;
          this->CONCENTRATIONBmax[i]=MMin+(i+1)*logdeltaCONCENTRATION;
          this->CONCENTRATIONBin[i]=this->CONCENTRATIONBmin[i]+(i+0.5)*logdeltaCONCENTRATION;
        }
    }
  if(prop==_SPIN_ || prop==_SPIN_BULLOCK_)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Generating",Nbins," bin(s) (zeroth -or  first bin- is the full sample) to tracer-spin");
#endif
      real_prec  MMin = log10(this->params._SPINBULLOCKbins_min(0));
      real_prec  MMax = log10(this->params._SPINBULLOCKbins_max(this->params._NSPINBULLOCKbins_power()-1));
      this->logdeltaSPIN =(MMax-MMin)/(static_cast<real_prec>(Nbins));
      this->SPINBmin.clear();this->SPINBmin.shrink_to_fit();
      this->SPINBmin.resize(Nbins,0);
      this->SPINBmax.clear();this->SPINBmax.shrink_to_fit();
      this->SPINBmax.resize(Nbins,0);
      this->SPINBin.clear();this->SPINBin.shrink_to_fit();
      this->SPINBin.resize(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i< Nbins;++i)
        {
          this->SPINBmin[i]=MMin+i*logdeltaSPIN;
          this->SPINBmax[i]=MMin+(i+1)*logdeltaSPIN;
          this->SPINBin[i]=this->SPINBmin[i]+(i+0.5)*this->logdeltaSPIN;
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Estimates of Abundance as a function of halo properties
void HaloTools::get_property_function(string file)
{
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif
#ifdef accum
  bool accumulate=true;
#else
  bool accumulate=false;
#endif
  //Her we get by default the mass function and if vmax is used, the v,ax_function
  this->define_property_bins();
  // Mvir
  if(this->params._i_mass_g()>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring mass function");
#endif
      real_prec lm_min=this->params._LOGMASSmin();
      real_prec lm_max=this->params._LOGMASSmax();
      real_prec units_observable_m=this->params._MASS_units();
      this->mass_function.resize(this->params._NMASSbins_mf());
      vector<int >mf_counts(this->params._NMASSbins_mf(),0);
      int im=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(size_t i=0;i<this->catalogue._NOBJS();++i)
        {
         real_prec property=log10(this->catalogue.mass_at(i))+log10(units_observable_m);
#ifndef accum
          if(property >=lm_min && property<=lm_max)
           {
#endif
          im=get_bin(property,lm_min,this->params._NMASSbins_mf(),this->logdeltaM,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
         mf_counts[im]++;
#ifndef accum
          }
#endif
       }
      size_t naux=0;
#pragma omp parallel for reduction(+:naux)
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
       naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string file_m=file+"_mass";
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Mass function written in file ", file_m);
#endif
     for(size_t i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
           mout<<this->MBin[i]<<"\t"<<this->mass_function[i]<<"\t"<<mf_counts[i]<<endl;
         }
      mout.close();
#ifdef _VERBOSE_CAT_
      So.DONE();
#endif
      /*
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
   vector<pair<real_prec, real_prec> > xy_pts_ref;
   for(size_t i=0; i<this->mass_function.size(); ++i)
       xy_pts_ref.push_back(std::make_pair(this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf()), log10(this->mass_function[i])));
     this->gp_abundance<<"set log x \n";
     this->gp_abundance<<"set border linewidth 2.2\n";
     this->gp_abundance<<"set xlabel 'Mass [Ms / h]' font 'Times-Roman,12'\n";
     this->gp_abundance<<"set ylabel 'log n(M) h /Ms (h / Mpc)³]' font 'Times-Roman,12'\n";
     this->gp_abundance<<"plot"<<this->gp_abundance.file1d(xy_pts_ref) << "w l lw 2 lt 8 title 'Reference'"<<endl;
#endif
*/
    vector<real_prec>mcount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
     for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
      mcount[i]+=mf_counts[j]/static_cast<double>(this->catalogue._NOBJS());
    file_m+="_cumulative";
    this->File.write_to_file(file_m,this->MBin,mcount); 
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative Mass function written in file ", file_m);
      So.DONE();
#endif
    mf_counts.clear(); mf_counts.shrink_to_fit();
    mcount.clear();mcount.shrink_to_fit();
    }
  // Vmax
#if defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_ || defined _ASSIGN_VMAX_POST_
  if(this->params._i_vmax_g()>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring VMAX function");
#endif
      real_prec v_min=log10(this->params._VMAXmin());
      real_prec v_max=log10(this->params._VMAXmax());
      this->vmax_function.resize(this->params._NMASSbins_mf());
      vector<int >v_counts(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec property=log10(this->catalogue.vmax_at(i));
#ifndef accum
      if(property >=v_min && property<=v_max)
        {
#endif
          int im=get_bin(property,v_min,this->params._NMASSbins_mf(),this->logdeltaVMAX,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
           v_counts[im]++;
#ifndef accum
        }
#endif
    }
    int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
        naux+=v_counts[i];

#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Vmax) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Vmax) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string filev=file+"_vmax";
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting VMAX-function in file ", filev);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
        {
          real_prec mint=(this->VMAXBmax[i]-this->VMAXBmin[i]);
          this->vmax_function[i]=static_cast<real_prec>(v_counts[i])/(pow(this->params._Lbox(),3)*mint);
       }
      this->File.write_to_file(filev,this->VMAXBin,this->vmax_function,v_counts); 
      vector<real_prec>Vcount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i) // N(>Vmax)
        for(int j=i+1; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
            Vcount[i]+=v_counts[j]/static_cast<double>(this->catalogue._NOBJS());
      v_counts.clear();v_counts.shrink_to_fit();
      filev+="_cumulative";
      So.message_screen("Cumulative fraction of objects above Vmax N(>Vmax)/Ntot function written in file ", filev);
      this->File.write_to_file(filev,this->VMAXBin,Vcount); 
#ifdef _SET_THRESHOLD_FOR_PROP_RANDOM_TRACERS_
      if(this->aux_flag==true)
       {
      /** These lines are commneted for they were not giving a good solution to the finding of Vmax threshold.
          vector<gsl_real>Vaux;
          vector<gsl_real>Vbin;
          for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
          {
          gsl_real ppa=static_cast<gsl_real>(1.-Vcount[i]/static_cast<double>(this->catalogue._NOBJS()));
          //          gsl_real ppa=static_cast<gsl_real>(Vcount[i]/static_cast<double>(this->catalogue._NOBJS()));
          if(ppa > 0.01 && ppa < 0.93)
          {
              Vaux.push_back(ppa);
              Vbin.push_back(static_cast<gsl_real>(this->VMAXBin[i]));
          }
          }
          // Calculate the Vmax such that f(<=Vmax) = this->fraction_tracer_from_random;
          // from this, Vmax above the threshold will be assigned to DM
          // and values bvelow the threshold will be assigned to random placed tracers
          this->Prop_threshold_rand_dm=gsl_inter_new(Vaux,Vbin,this->fraction_tracer_from_random);
          //       this->Prop_threshold_rand_dm=gsl_inter_new(Vaux,Vbin,this->fraction_tracer_from_dm);
          cout<<endl;
          So.message_screen("Threshold for Vmax assigment to randoms = ", this->Prop_threshold_rand_dm);
          cout<<endl;
          Vbin.clear();Vbin.shrink_to_fit();
          Vaux.clear();Vaux.shrink_to_fit();
      **/
      // These lines order the ref catalog in ascending order with respect to vmax and associate to each tracer an index pin-poiting that order
          // such that in posterior loops, when the running index over the tracer reaches the number of tracers associated to randoms (in the assignment ala LPT),
      // the code will assign the los vmax to the randoms, and the high vmax to the dm particles
#ifdef _VERBOSE_CAT_
      cout<<endl;
#endif
      So.message_screen("Sorting values of Vmax in reference catalog:");
      gsl_vector *vmax_aux;
      vmax_aux=gsl_vector_alloc(this->catalogue._NOBJS());

      gsl_vector *id_gal;
      id_gal=gsl_vector_alloc(this->catalogue._NOBJS());

      for(size_t ig=0;ig<this->catalogue._NOBJS();++ig )
        {
          gsl_vector_set(vmax_aux,ig,this->catalogue.vmax_at(ig));
          gsl_vector_set(id_gal,ig,ig);
        }
      gsl_sort_vector2(vmax_aux,id_gal) ;   // sort the vmax and correspondingly their associated the gal id
      this->Prop_threshold_rand_dm=gsl_vector_get(vmax_aux, this->Ntracers_ran-1);// NOTE THAT this value won't be used in practice
#ifdef _VERBOSE_CAT_
      cout<<endl;
      So.message_screen("Threshold for Vmax assigment to randoms = ", this->Prop_threshold_rand_dm," km/s");
      cout<<endl;
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i ) //ORDERED loop in ascending order in vmax: i=0 is the lowest vmax value, i=NOBJS is the highest
        {
          // To each halo, associates the index i sorted in ascening order with respect to vmax
          // This will be used when assigning Vmax to the random tracers. Since we have Nrandom, we select the refernece tracers with this index <Nrandom to prodive their vmax to the random tracers
          // such athat these will properñly probe the Vmax end of the n(Vmax) function.
          size_t gal_id=gsl_vector_get(id_gal,i);
          this->catalogue[gal_id].vmax_index=i;
        }
      gsl_vector_free(vmax_aux);
      gsl_vector_free(id_gal);
      So.DONE();
    }
#endif
      Vcount.clear();Vcount.shrink_to_fit();
    }
#endif
  if(this->params._i_rs_g()>0) // Rs
    {
      int im=0;
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring RS function");
#endif
      real_prec rs_min=log10(this->params._RSmin());
      real_prec rs_max=log10(this->params._RSmax());
      real_prec units_observable_rs=1;
      this->rs_function.resize(this->params._NMASSbins_mf());
      vector<int >rs_counts(this->params._NMASSbins_mf(),0);

#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec property=log10(this->catalogue.rs_at(i)*units_observable_rs);
#ifndef accum
      if(property >=rs_min && property<=rs_max)
        {
#endif
          im=get_bin(property,rs_min,this->params._NMASSbins_mf(),this->logdeltaRS,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
              rs_counts[im]++;
#ifndef accum
        }
#endif
    }
    int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
      naux+=rs_counts[i];
#ifdef _VERBOSE_CAT_
    So.message_screen("Number of objects used in mass function n(Rs) = ", naux);
    if(this->catalogue._NOBJS()-naux>0)
      So.message_screen("Number of objects EXCLUDED in n(Rs) estimation = ", this->catalogue._NOBJS()-naux);
#endif
    string filers=file+"_rs";

#ifdef _VERBOSE_CAT_
    So.message_screen("Writting Rs-function in file ", filers);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->RSBmax[i]-this->RSBmin[i]);
      this->rs_function[i]=static_cast<real_prec>(rs_counts[i])/(pow(this->params._Lbox(),3)*mint);
    }
   this->File.write_to_file(filers,this->RSBin,this->rs_function,s_counts); 
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
      for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
    #pragma omp atomic
#endif
         RScount[i]+=this->rs_function[j]*(this->RSBmax[j]-this->RSBmin[j]);
    filers+="_cumulative";
    this->File.write_to_file(filers,this->RSBin,Rscount); 
    RScount.clear();RScount.shrink_to_fit();
#ifdef _VERBOSE_CAT_
    So.message_screen("Cumulative Rs-function written in file ", filers);
#endif
#endif
    }
//-----------
// SPIN
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
   if(this->params._i_spin_g()>0 || this->params._i_spin_bullock_g()>0)  // SPIN
    {
      int im=0;
      So.message_screen("Measuring SPIN-function");
      real_prec s_min=log10(this->params._SPINmin());
      real_prec s_max=log10(this->params._SPINmax());
      this->s_function.resize(this->params._NMASSbins_mf());
      vector<int >s_counts(this->params._NMASSbins_mf(),0);
      string filers=file+"_spin";
      ofstream mout;
      mout.open(filers.c_str());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<this->catalogue._NOBJS();++i)
      {
       real_prec property=log10(this->catalogue.spin_bullock_at(i));
#ifndef accum
      if(property >=s_min && property<=s_max)
        {
#endif
          im=get_bin(property,s_min,this->params._NMASSbins_mf(),this->logdeltaSPIN,accumulate);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          s_counts[im]++;
#ifndef accum
        }
#endif
      }
    int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
      naux+=s_counts[i];
#ifdef _USE_OMP_
#pragma omp parallel for 
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->SPINBmax[i]-this->SPINBmin[i]);
      this->s_function[i]=static_cast<real_prec>(s_counts[i])/(pow(this->params._Lbox(),3)*mint);
    }
    this->File.write_to_file(filers,this->SPINBin,this->s_function,s_counts); 
    s_counts.clear();s_counts.shrink_to_fit();
    vector<real_prec>Scount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
      for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
        Scount[i]+=this->s_function[j]*(this->SPINBmax[j]-this->SPINBmin[j]);
    filers+="_cumulative";
    this->File.write_to_file(filers,this->SPINBin,Scount); 
    Scount.clear();Scount.shrink_to_fit();
    So.message_screen("Cumulative Spin-function written in file ", filers);
    }
#endif
//-----------
// CVIR
  bool use_cvir=false;
#ifdef _USE_CONCENTRATION_AS_DERIVED_OBSERVABLE_
  use_cvir=true;
#endif
   if(true==use_cvir && this->params._i_rs_g()>0)  
    {
      int im=0;
      So.message_screen("Measuring CVIR-function");
      real_prec s_min=log10(this->params._CONCENTRATIONmin());
      real_prec s_max=log10(this->params._CONCENTRATIONmax());
      this->cvir_function.resize(this->params._NMASSbins_mf());
      vector<int >s_counts(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<this->catalogue._NOBJS();++i)
      {
       real_prec property=log10(this->catalogue.concentration_at(i));
#ifndef accum
      if(property >=s_min && property<=s_max)
        {
#endif
          im=get_bin(property,s_min,this->params._NMASSbins_mf(),this->logdeltaCONCENTRATION,accumulate);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          s_counts[im]++;
#ifndef accum
        }
#endif
      }
    int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
      naux+=s_counts[i];
#ifdef _VERBOSE_CAT_
    So.message_screen("Number of objects used in mass function n(cvir) = ", naux);
    if(this->catalogue._NOBJS()-naux>0)
      So.message_screen("Number of objects EXCLUDED in n(cvir) estimation = ", this->catalogue._NOBJS()-naux);
#endif
    string filers=file+"_cvir";
#ifdef _VERBOSE_CAT_
    So.message_screen("Writting Cvir-function in file ", filers);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->CONCENTRATIONBmax[i]-this->CONCENTRATIONBmin[i]);
      this->cvir_function[i]=static_cast<real_prec>(s_counts[i])/(pow(this->params._Lbox(),3)*mint);
    }
    this->File.write_to_file(filers,this->CONCENTRATIONBin,this->cvir_function,s_counts); 
    s_counts.clear();s_counts.shrink_to_fit();
    vector<real_prec>Scount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
      for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
        Scount[i]+=this->cvir_function[j]*(this->CONCENTRATIONBmax[j]-this->CONCENTRATIONBmin[j]);
    filers+="_cumulative";
    this->File.write_to_file(filers,this->CONCENTRATIONBin,Scount); 
    Scount.clear();Scount.shrink_to_fit();
    So.message_screen("Cumulative Spin-function written in file ", filers);
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Estimates of Abundance as a function of halo properties
void HaloTools::get_property_function_cosmic_web_types(string file)
{
#ifndef _USE_CWC_
  So.message_screen("Using ", __PRETTY_FUNCTION__);
  So.message_screen("Preprocessor directive _USE_CWC_ is not defined. Check in def.h file. Code ends here");
  exit(1);
#endif
  this->cwclass.set_params(this->params);
  string file_X=this->params._Input_Directory_X()+this->params._Name_Catalog_X();
  vector<real_prec>deltam(this->params._NGRID(),0);
  this->File.read_array(file_X,deltam);
  get_overdens(deltam, deltam);
  this->cwclass.get_CWC(deltam);
  So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
#endif
#ifdef accum
  bool accumulate=true;
#else
  bool accumulate=false;
#endif
  //Her we get by default the mass function and if vmax is used, the v,ax_function
  this->define_property_bins();

  // Mvir
  if(this->params._i_mass_g()>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring mass function");
#endif
      real_prec lm_min=this->params._LOGMASSmin();
      real_prec lm_max=this->params._LOGMASSmax();
      real_prec units_observable_m=this->params._MASS_units();
      for(int ict=0;ict<this->cwclass.cwt_used.size();ict++)
      {
        this->mass_function.resize(this->params._NMASSbins_mf());
        vector<int >mf_counts(this->params._NMASSbins_mf(),0);
        int im=0;
        size_t partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
        for(size_t i=0;i<this->catalogue._NOBJS();++i)
          {
           real_prec property=log10(this->catalogue.mass_at(i))+log10(units_observable_m);
#ifndef accum
            if(property >=lm_min && property<=lm_max && this->cwclass.get_Tclassification(this->catalogue.GridID_at(i))==static_cast<WebType>(ict))
             {
#endif
              partial_o++;
              im=get_bin(property,lm_min,this->params._NMASSbins_mf(),this->logdeltaM,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
             mf_counts[im]++;
#ifndef accum
            }
#endif  
          }
        size_t naux=0;
#pragma omp parallel for reduction(+:naux)
        for(size_t i=0;i< this->params._NMASSbins_mf();++i)
          naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
        So.message_screen("Number of objects used in mass function n(M) = ", naux);
        if(this->catalogue._NOBJS()-naux>0)
          So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->catalogue._NOBJS()-naux);
#endif
        string file_m=file+"_mass_cwt"+to_string(ict);
        ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(size_t i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
         mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();
#ifdef _VERBOSE_CAT_
      So.DONE();
#endif
      vector<real_prec>mcount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
    for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
          mcount[i]+=mf_counts[j];
     mf_counts.clear(); mf_counts.shrink_to_fit();

      file_m+="_cumulative_cwt"+to_string(ict);
      mout.open(file_m.c_str());
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
    mout<<this->MBin[i]<<"\t"<<mcount[i]/static_cast<double>(partial_o)<<endl;
      mout.close();
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative Mass function written in file ", file_m);
      So.DONE();
#endif
      mcount.clear();mcount.shrink_to_fit();
    }
}// closes loop pver cwtypes
  // Vmax
#if defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_ || defined _ASSIGN_VMAX_POST_
  if(this->params._i_vmax_g()>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring VMAX function");
#endif
      real_prec v_min=log10(this->params._VMAXmin());
      real_prec v_max=log10(this->params._VMAXmax());

      for(int ict=0;ict<this->cwclass.cwt_used.size();ict++)
      {
      real_prec units_observable_v=1;
      this->vmax_function.resize(this->params._NMASSbins_mf());
      vector<int >v_counts(this->params._NMASSbins_mf(),0);
     size_t partial_o=0;

#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec property=log10(this->catalogue.vmax_at(i)*units_observable_v);
#ifndef accum
      if(property >=v_min && property<=v_max && this->cwclass.get_Tclassification(this->catalogue.GridID_at(i))==static_cast<WebType>(ict))
        {
#endif
            partial_o++;
              int im=get_bin(property,v_min,this->params._NMASSbins_mf(),this->logdeltaVMAX,accumulate);
//              cout<<RED<<this->params._NMASSbins_mf()<<"   "<<im<<"   "<<v_min<<"   "<<v_max<<"  "<<property<<RESET<<endl;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
              v_counts[im]++;
#ifndef accum
        }
#endif
    }
    int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
        naux+=v_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Vmax) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Vmax) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string filev=file+"_vmax_cwt"+to_string(ict);
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting VMAX-function in file ", filev);
#endif
      ofstream mout;
      mout.open(filev.c_str());
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
        {
          real_prec mint=(this->VMAXBmax[i]-this->VMAXBmin[i]);
          this->vmax_function[i]=static_cast<real_prec>(v_counts[i])/(pow(this->params._Lbox(),3)*mint);
          mout<<this->VMAXBmin[i]+(i+0.5)*(this->VMAXBmax[i]-this->VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->vmax_function[i]<<"   "<<v_counts[i]<<endl;
       }
     mout.close();
      vector<real_prec>Vcount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i) // N(>Vmax)
        for(int j=i+1; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
            Vcount[i]+=v_counts[j];

      v_counts.clear();v_counts.shrink_to_fit();
      filev+="_cumulative_cwt"+to_string(ict);
      mout.open(filev.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative fraction of objects above Vmax N(>Vmax)/Ntot function written in file ", filev);
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
    mout<<this->VMAXBin[i]<<"\t"<<Vcount[i]/static_cast<double>(partial_o)<<endl;
      mout.close();
    }
}
#endif
// Rs
  if(this->params._i_rs_g()>0)
    {
      int im=0;
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring RS function");
#endif
      real_prec rs_min=log10(this->params._RSmin());
      real_prec rs_max=log10(this->params._RSmax());
      real_prec units_observable_rs=1;
      for(int ict=0;ict<this->cwclass.cwt_used.size();ict++)
      {
      this->rs_function.resize(this->params._NMASSbins_mf());
      vector<int >rs_counts(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec property=log10(this->catalogue.rs_at(i)*units_observable_rs);
#ifndef accum
      if(property >=rs_min && property<=rs_max)
        {
#endif
          im=get_bin(property,rs_min,this->params._NMASSbins_mf(),this->logdeltaRS,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
              rs_counts[im]++;
#ifndef accum
        }
#endif
    }
      int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)

    naux+=rs_counts[i];

#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Rs) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Rs) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string filers=file+"_rs_cwt"+to_string(ict);

#ifdef _VERBOSE_CAT_
      So.message_screen("Writting Rs-function in file ", filers);
#endif
      ofstream mout;
      mout.open(filers.c_str());
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->RSBmax[i]-this->RSBmin[i]);
      this->rs_function[i]=static_cast<real_prec>(rs_counts[i])/(pow(this->params._Lbox(),3)*mint);
      mout<<this->RSBmin[i]+(i+0.5)*(this->RSBmax[i]-this->RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->rs_function[i]<<"   "<<rs_counts[i]<<endl;
    }
      mout.close();
      rs_counts.clear();rs_counts.shrink_to_fit();
      vector<real_prec>RScount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
        for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#pragma omp atomic
          RScount[i]+=this->rs_function[j]*(this->RSBmax[j]-this->RSBmin[j]);

      filers+="_cumulative_cwt"+to_string(ict);
      mout.open(filers.c_str());
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
    mout<<this->RSBin[i]<<"\t"<<RScount[i]<<endl;
      mout.close();
      RScount.clear();RScount.shrink_to_fit();
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative Rs-function written in file ", filers);
#endif
}
#endif
}
  // SPIN
  if(this->params._i_spin_g()>0)
    {
      int im=0;
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring SPIN-function");
#endif
      real_prec s_min=log10(this->params._SPINmin());
      real_prec s_max=log10(this->params._SPINmax());
      real_prec units_observable_s=1;
      for(int ict=0;ict<this->cwclass.cwt_used.size();ict++)
      {
      this->s_function.resize(this->params._NMASSbins_mf());
      vector<int >s_counts(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec property=log10(this->catalogue.spin_at(i)*units_observable_s);
#ifndef accum
      if(property >=s_min && property<=s_max)
        {
#endif
          im=get_bin(property,s_min,this->params._NMASSbins_mf(),this->logdeltaSPIN,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
              s_counts[im]++;
#ifndef accum
        }
#endif
      }
      int naux=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:naux)
#endif
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    naux+=s_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(spin) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(spin) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string filers=file+"_spin_cwt"+to_string(ict);
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting Spin-function in file ", filers);
#endif
      ofstream mout;
      mout.open(filers.c_str());
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->SPINBmax[i]-this->SPINBmin[i]);
      this->s_function[i]=static_cast<real_prec>(s_counts[i])/(pow(this->params._Lbox(),3)*mint);
      mout<<this->SPINBmin[i]+(i+0.5)*(this->SPINBmax[i]-this->SPINBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->rs_function[i]<<"   "<<s_counts[i]<<endl;
    }
      mout.close();
      s_counts.clear();s_counts.shrink_to_fit();
      vector<real_prec>Scount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
        for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
            Scount[i]+=this->s_function[j]*(this->SPINBmax[j]-this->SPINBmin[j]);
      filers+="_cumulative_cwt"+to_string(ict);
      mout.open(filers.c_str());
      for(size_t i=0; i< this->params._NMASSbins_mf(); ++i)
    mout<<this->SPINBin[i]<<"\t"<<Scount[i]<<endl;
      mout.close();
      Scount.clear();Scount.shrink_to_fit();
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative Spin-function written in file ", filers);
#endif
}
#endif
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_HOD(string file)
{
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif
#ifdef accum
  bool accumulate=true;
#else
  bool accumulate=false;
#endif
  //Her we get by default the mass function and if vmax is used, the v,ax_function
  this->define_property_bins();

  if(this->params._i_mass_g()<0)
    So.message_warning("Mass nor properly specified in params file");
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring HOD");
#endif
      real_prec lm_min=this->params._LOGMASSmin();
      real_prec lm_max=this->params._LOGMASSmax();
      real_prec units_observable_m=this->params._MASS_units();

      this->mass_function.resize(this->params._NMASSbins_mf());
      vector<int >mf_counts(this->params._NMASSbins_mf(),0);

      int im=0;
      size_t partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
    for(size_t i=0;i<this->catalogue._NOBJS();++i)
        {
         real_prec property=log10(this->catalogue.mass_at(i))+log10(units_observable_m);
#ifndef accum
          if(property >=lm_min && property<=lm_max)
           {
#endif
            partial_o++;
           im=get_bin(property,lm_min,this->params._NMASSbins_mf(),this->logdeltaM,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
         mf_counts[im]++;
#ifndef accum
          }
#endif
       }
      size_t naux=0;
#pragma omp parallel for reduction(+:naux)
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
       naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string file_m=file+"_mass";
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(size_t i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
         mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();

  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_HOD_web_types(string file)
{
#ifndef _USE_CWC_
    So.message_screen("Using ", __PRETTY_FUNCTION__);
    So.message_screen("_USE_CWC_ not defined. Check. BAM ends here");
    exit(1);
#endif
    this->cwclass.set_params(this->params);
    string file_X=this->params._Input_Directory_X()+this->params._Name_Catalog_X();
    vector<real_prec>deltam(this->params._NGRID(),0);
    this->File.read_array(file_X,deltam);
    get_overdens(deltam, deltam);
    this->cwclass.get_CWC(deltam);
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  // omp_get_max_threads();
   omp_set_num_threads(NTHREADS);
#endif
#ifdef accum
  bool accumulate=true;
#else
  bool accumulate=false;
#endif
  //Her we get by default the mass function and if vmax is used, the v,ax_function
  this->define_property_bins();

  if(this->params._i_mass_g()<0)
    So.message_warning("Mass nor properly specified in params file");
#ifdef _VERBOSE_CAT_
      So.message_screen("Measuring HOD");
#endif
      real_prec lm_min=this->params._LOGMASSmin();
      real_prec lm_max=this->params._LOGMASSmax();
      real_prec units_observable_m=this->params._MASS_units();
      for(int ict=0;ict<this->cwclass.cwt_used.size();ict++)
      {

      this->mass_function.resize(this->params._NMASSbins_mf());
      vector<int >mf_counts(this->params._NMASSbins_mf(),0);

      int im=0;
      size_t partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
    for(size_t i=0;i<this->catalogue._NOBJS();++i)
        {
         real_prec property=log10(this->catalogue.mass_at(i))+log10(units_observable_m);
#ifndef accum
          if(property >=lm_min && property<=lm_max && this->cwclass.get_Tclassification(this->catalogue.GridID_at(i))==static_cast<WebType>(ict))
           {
#endif
            partial_o++;
           im=get_bin(property,lm_min,this->params._NMASSbins_mf(),this->logdeltaM,accumulate);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
         mf_counts[im]++;
#ifndef accum
          }
#endif
       }
      size_t naux=0;
#pragma omp parallel for reduction(+:naux)
      for(size_t i=0;i< this->params._NMASSbins_mf();++i)
    naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->catalogue._NOBJS()-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->catalogue._NOBJS()-naux);
#endif
      string file_m=file+"_mass_cwt"+to_string(ict);
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(size_t i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
         mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();
    }// closes loop pver cwtypes

  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools:: get_distribution_reduced_mass_in_cell(){
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
  So.message_screen("Measuring distribution of reduced mass in cells", MAXIMUM_DISTANCE_EXCLUSION, "Mpc/h");
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
#pragma omp parallel for
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      size_t ID=this->catalogue.GridID_at(i);
      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1_at(i));
      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2_at(i));
      cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3_at(i));
      cell_info_tr[ID].mass.push_back(this->catalogue.mass_at(i));
    }

  //  vector<int>dist_reduced_mass(N_BINS_REDUCED_MASS, 0);

  real_prec dmin=MINIMUM_DISTANCE_EXCLUSION;
  real_prec dmax=MAXIMUM_DISTANCE_EXCLUSION;
  vector<vector<int>> dist_reduced_mass(N_BINS_REDUCED_MASS,vector<int>(N_BINS_DIST_EX, 0));
  for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
    if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
      for(size_t i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
        for(int j=i+1;j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
          {
            real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
            int index_dist =get_bin(dist,dmin, N_BINS_DIST_EX,(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX),true);
            real_prec mu= cell_info_tr[id].mass[i]*cell_info_tr[id].mass[j]/(cell_info_tr[id].mass[i]+cell_info_tr[id].mass[j]);
            int index_mu =get_bin(mu,MIN_REDUCED_MASS,N_BINS_REDUCED_MASS,DELTA_REDUCED_MASS,true);
            dist_reduced_mass[index_mu][index_dist]++;
          }
  for(int id=0;id<N_BINS_DIST_EX;++id)
    {
      string file=this->Output_directory+"reduced_mass_dist_"+this->catalogue._type_of_object()+"_dbin"+to_string(id)+".txt";
      So.message_screen("Writting to file ", file);
      ofstream sal;
      sal.open(file.c_str());
      sal<<"#Bin info:     Dmin = "<<MINIMUM_DISTANCE_EXCLUSION+(id)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<"  Dmax = "<<MINIMUM_DISTANCE_EXCLUSION+(id+1.0)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<endl;
      for(size_t i=0;i<dist_reduced_mass.size();++i)
         sal<<MIN_REDUCED_MASS+(i+0.5)*DELTA_REDUCED_MASS<<" "<<dist_reduced_mass[i][id]<<endl;
      sal.close();
      So.DONE();
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function computes:
// minimum separation between pairs in cells
// mean separation between pairs in cells
// stdv separation between pairs in cells
// Allocates the result in a class member container min_separation_in_cell[ID]=
void HaloTools::get_stats_separation_in_cell()
{
#ifdef _FULL_VERBOSE_
    this->So.enter(__PRETTY_FUNCTION__);
#endif
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring minimum separation in cells (Mpc/h) ...");
  So.message_screen("Current type is ", this->catalogue._type_of_object());
  cout<<endl;
#endif
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
#ifdef _VERBOSE_CAT_
  So.message_screen("Getting positions inside cells:");
#endif
  if("TRACER_MOCK_ONLY_COORDS" != this->catalogue._type_of_object())
    {
      for(size_t i=0; i<this->catalogue._NOBJS(); ++i) //loop over the "observed objects", i.e, with cuts already set
    {
          size_t ID=this->catalogue.GridID_at(i);
      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1_at(i));
      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2_at(i));
      cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3_at(i));
#ifdef _test_ms_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      cell_info_tr[ID].property.push_back(this->catalogue.vmax_at(i));
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
      cell_info_tr[ID].property.push_back(this->catalogue.mass_at(i));
#endif // endif use_vmax_as _obs
#endif // endif _test_ms_
    }
    }
  else
    {
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed objects", i.e, with cuts already set
    {
          size_t ID=this->catalogue.GridID_at(i);
      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1_at(i));
      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2_at(i));
      cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3_at(i));
    }
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
  #ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
   this->min_separation_in_cell.clear();
   this->min_separation_in_cell.shrink_to_fit();
  this->min_separation_in_cell.resize(this->box.NGRID,0);
#endif
#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
   this->mean_separation_in_cell.clear();
   this->mean_separation_in_cell.shrink_to_fit();
   this->mean_separation_in_cell.resize(this->box.NGRID,0);
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
   this->stdv_separation_in_cell.clear();
   this->stdv_separation_in_cell.shrink_to_fit();
   this->stdv_separation_in_cell.resize(this->box.NGRID,0);
#endif
  // If the tracer has poperties, then ask for them when computing separations. If not, else below
//   if("TRACER_MOCK_ONLY_COORDS" != this->catalogue._type_of_object())
//     {
        for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
          {
            real_prec aux_min_distance_a;
            if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
            {

#if defined _USE_MEAN_SEPARATIONS_IN_CELLS_  || defined  _USE_STDV_SEPARATIONS_IN_CELLS_
                vector<real_prec>sep_cell;
#endif
                  for(size_t i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
            {
#if defined _USE_MIN_SEPARATIONS_IN_CELLS_
                      real_prec aux_min_distance_b=LARGE_NUMBER;
#endif
                          for(size_t j=i+1 ; j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
                {
#ifdef _USE_LOG_MSIC_
                              real_prec dist=  0.5*log10(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
#else
                              real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
#endif

#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
                              sep_cell.push_back(dist);
#endif
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
                              aux_min_distance_b=min(dist, aux_min_distance_a);
                              aux_min_distance_a=aux_min_distance_b;
#endif
                          }
            }

#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
                  this->mean_separation_in_cell[id]=get_mean(sep_cell);//allocate here the mean separations in cells
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
                  this->stdv_separation_in_cell[id]=get_var(sep_cell);//allocate here the mean separations in cells
#endif

#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
#ifdef _USE_LOG_MSIC_
                  this->min_separation_in_cell[id]=aux_min_distance_a >0 ? aux_min_distance_a:MIN_SEP_IN_CELLS;
#else
                  this->min_separation_in_cell[id]=aux_min_distance_a;
#endif
#endif
            }
            else{
#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
                this->mean_separation_in_cell[id]=MIN_MEAN_SEP_IN_CELLS;
#endif

#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
                this->min_separation_in_cell[id]=MIN_SEP_IN_CELLS;
#endif
            }
        }
/*	}
      else // This is menat for tracers with only coordinates
    {
      for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
        {
          real_prec aux_min_distance_a;
          if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
        {
          aux_min_distance_a=100.0;
          for(size_t i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
            {
              real_prec aux_min_distance_b;
              for(int j=i+1 ; j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
            {
              real_prec dist=  sqrt(pow(cell_info_tr[id].posx_p[i]-cell_info_tr[id].posx_p[j],2)+pow(cell_info_tr[id].posy_p[i]-cell_info_tr[id].posy_p[j],2)+pow(cell_info_tr[id].posz_p[i]-cell_info_tr[id].posz_p[j],2));
              aux_min_distance_b=min(dist, aux_min_distance_a);
              aux_min_distance_a=aux_min_distance_b;
            }
            }
#ifdef _USE_LOG_MSIC_
          this->min_separation_in_cell[id]=aux_min_distance_a>0? log10(aux_min_distance_a):MIN_SEP_IN_CELLS;
#else
          this->min_separation_in_cell[id]=aux_min_distance_a;
#endif
        }
          else
        this->min_separation_in_cell[id]=MIN_SEP_IN_CELLS;
        }
        }*/
      cell_info_tr.clear();cell_info_tr.shrink_to_fit();
      size_t nbins_new=50;
#ifdef _USE_LOG_MSIC_
#pragma omp parallel for
      for(size_t id=0; id<this->box.NGRID ;++id)
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
        min_separation_in_cell[id]=pow(10,min_separation_in_cell[id]);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
         mean_separation_in_cell[id]=pow(10,mean_separation_in_cell[id]);
#endif
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
      File.write_array(this->params._Output_directory()+"min_sep"+this->params._type_of_object(), min_separation_in_cell);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
      File.write_array(this->params._Output_directory()+"mean_sep"+this->params._type_of_object(), mean_separation_in_cell);
#endif
#pragma omp parallel for
    for(size_t id=0; id<this->box.NGRID ;++id)
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
        min_separation_in_cell[id]=log10(min_separation_in_cell[id]);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
        mean_separation_in_cell[id]=log10(mean_separation_in_cell[id]);
#endif

#endif // end for _USE_LOG_MISC
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
     vector<size_t>dist_min_sep(nbins_new,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
         if(this->min_separation_in_cell[id]>MIN_SEP_IN_CELLS)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              dist_min_sep[get_bin(this->min_separation_in_cell[id],MIN_SEP_IN_CELLS,nbins_new,DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/nbins_new,true)]++;
      So.DONE();

      this->min_halo_separation=get_min_nm(this->min_separation_in_cell, true);
/*
#ifdef _USE_LOG_MSIC_
      So.message_screen("Min value of log min_sep_in_cell= ", this->min_halo_separation);
      So.message_screen("Max value of log min_sep_in_cell= ", get_max<float>(this->min_separation_in_cell));
#else
      So.message_screen("Min value of min_sep_in_cell= ", this->min_halo_separation, "Mpc /h");
      So.message_screen("Max value of min_sep_in_cell= ", get_max<float>(this->min_separation_in_cell), "Mpc /h");
#endif // end for #ifdef _USE_LOG_MSIC_

*/
#endif // end for #ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
      vector<size_t>dist_mean_sep(nbins_new,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
          if(this->mean_separation_in_cell[id]>MIN_MEAN_SEP_IN_CELLS)
#ifdef _USE_OMP_
#pragma omp atomic
#endif
               dist_mean_sep[get_bin(this->mean_separation_in_cell[id],MIN_MEAN_SEP_IN_CELLS,nbins_new,DELTA_MEAN_SEP*N_BINS_MEAN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new),true)]++;
       So.DONE();
       this->min_mean_halo_separation=get_min_nm(this->mean_separation_in_cell);
       So.message_screen("Min value of mean_sep_in_cell= ", get_min_nm(this->mean_separation_in_cell));
       So.message_screen("Max value of mean_sep_in_cell= ", get_max_nm(this->mean_separation_in_cell));
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
       vector<size_t>dist_stdv_sep(nbins_new,0);
      for(size_t id=0; id<this->box.NGRID ;++id) // loop over cells
          if(this->stdv_separation_in_cell[id]>MIN_STDV_SEP_IN_CELLS)
               dist_stdv_sep[get_bin(this->stdv_separation_in_cell[id],MIN_STDV_SEP_IN_CELLS,nbins_new,DELTA_STDV_SEP*N_BINS_STDV_SEP_IN_CELLS/nbins_new,true)]++;
       So.DONE();
       So.message_screen("Min value of stdv_sep_in_cell= ", get_min_nm(this->stdv_separation_in_cell));
       So.message_screen("Max value of stdv_sep_in_cell= ", get_max_nm(this->stdv_separation_in_cell));
#endif
      ofstream sal;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
      string fileo=this->Output_directory+"min_sep_in_cells_dist_"+this->catalogue._type_of_object()+".txt";
      So.message_screen("Writting distribution of minimum separations within cells");
      So.message_screen("in file ", fileo);
      sal.open(fileo.c_str());
      for(size_t i=0;i<dist_min_sep.size();++i)
#ifdef _USE_LOG_MSIC_
          sal<<pow(10,MIN_SEP_IN_CELLS+(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new))<<"  "<<dist_min_sep[i]<<endl;
#else
          sal<<(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new)<<"  "<<dist_min_sep[i]<<endl;
#endif
      sal.close();
#
#endif
#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
      string filemean=this->Output_directory+"mean_sep_in_cells_dist_"+this->catalogue._type_of_object()+".txt";
      So.message_screen("Writting distribution of mean separations within cells");
      So.message_screen("in file ", filemean);
      sal.open(filemean.c_str());
      for(size_t i=0;i<dist_mean_sep.size();++i)
#ifdef _USE_LOG_MSIC_
          sal<<pow(10,MIN_SEP_IN_CELLS+(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new))<<"  "<<dist_min_sep[i]<<endl;
#else
          sal<<(i+0.5)*DELTA_MEAN_SEP*N_BINS_MEAN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new)<<"  "<<dist_mean_sep[i]<<endl;
#endif
      sal.close();
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
      string filestdv=this->Output_directory+"stdv_sep_in_cells_dist_"+this->catalogue._type_of_object()+".txt";
      So.message_screen("Writting distribution of stdv of separations within cells");
      So.message_screen("in file ", filestdv);
      sal.open(filestdv.c_str());
      for(size_t i=0;i<dist_stdv_sep.size();++i)
#ifdef _USE_LOG_MSIC_
          sal<<pow(10,MIN_STDV_SEP_IN_CELLS+(i+0.5)*DELTA_STDV_MIN_SEP*N_BINS_STDV_SEP_IN_CELLS/static_cast<real_prec>(nbins_new))<<"  "<<dist_stdv_sep[i]<<endl;
#else
          sal<<(i+0.5)*DELTA_STDV_SEP*N_BINS_STDV_SEP_IN_CELLS/static_cast<real_prec>(nbins_new)<<"  "<<dist_stdv_sep[i]<<endl;
#endif
      sal.close();
#endif

So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This ,ethod computes the distribution separation between pairs in one cell
// Allocates the result in a class member container min_separation_in_cell[ID]=
void HaloTools::get_sep_distribution(real_prec min_value_prop)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring separation distribution");
  So.message_screen("Current type is ", this->catalogue._type_of_object());
  cout<<endl;
#endif
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
int Nbins_sep=100;
vector<real_prec>dist_sep(Nbins_sep,0);
real_prec max_sep=log10(300.);
real_prec min_sep=log10(0.01);
real_prec delta_sep = (max_sep-min_sep)/static_cast<real_prec>(Nbins_sep);
#ifdef _VERBOSE_CAT_
  So.message_screen("Cutting sample");
  cout<<endl;
#endif
vector<real_prec>x;
vector<real_prec>y;
vector<real_prec>z;
 for(size_t ig=0; ig<this->catalogue._NOBJS() ;++ig) // loop over cells
  if(this->catalogue.vmax_at(ig)>min_value_prop)
   {
      x.push_back(this->catalogue.coord1_at(ig));
      y.push_back(this->catalogue.coord2_at(ig));
      z.push_back(this->catalogue.coord3_at(ig));
    }
 real_prec normal=1;//static_cast<double>(x.size())*(x.size()-1)*0.5;
#ifdef _VERBOSE_CAT_
  So.message_screen("Getting separations");
  cout<<endl;
#endif
size_t i=0,j=0;
#ifdef _USE_OMP_
#pragma omp parallel private (i,j)
  {
    vector<real_prec>h_par(Nbins_sep,0);
#pragma omp for nowait
#endif
   for(i=0; i<x.size() ;++i) // loop over cells
     for(j=i+1; j < x.size();++j) // loop overtracers in that cell
        {
          real_prec xa = x[i]-x[j];
          real_prec ya = y[i]-y[j];
          real_prec za = z[i]-z[j];
          real_prec dist= 0.5*log10(xa*xa + ya*ya + za*za);
          int bin_sep=get_bin(dist, min_sep, Nbins_sep,delta_sep,true);
#ifdef _USE_OMP_
           h_par[bin_sep]++;
#else
          dist_sep[bin_sep]++;
#endif
#ifdef bc
          xa = x[i]-(x[j]-this->params._Lbox());
          ya = y[i]-(y[j]-this->params._Lbox());
          za = z[i]-(z[j]-this->params._Lbox());
          dist= 0.5*log10(xa*xa + ya*ya + za*za);
          bin_sep=get_bin(dist, min_sep, Nbins_sep,delta_sep,true);
#ifdef _USE_OMP_
           h_par[bin_sep]++;
#else
          dist_sep[bin_sep]++;
#endif
          xa = x[i]-(x[j]+this->params._Lbox());
          ya = y[i]-(y[j]+this->params._Lbox());
          za = z[i]-(z[j]+this->params._Lbox());
          dist= 0.5*log10(xa*xa + ya*ya + za*za);
          bin_sep=get_bin(dist, min_sep, Nbins_sep,delta_sep,true);
#ifdef _USE_OMP_
           h_par[bin_sep]++;
#else
          dist_sep[bin_sep]++;
#endif
#endif
     }
#ifdef _USE_OMP_
#pragma omp critical
    for(size_t i=0;i<dist_sep.size();++i)
      dist_sep[i]+=h_par[i];
   }
#endif
    ofstream sal;
    string fileo;
    fileo=this->Output_directory+"sep_dist_"+this->catalogue._type_of_object()+this->params._Name_survey()+"_Real"+to_string(this->params._realization())+".txt";
    So.message_screen("Writting distribution of minimum separations within cells");
    So.message_screen("in file ", fileo);
    sal.open(fileo.c_str());
    for(size_t i=0;i<dist_sep.size()-1;++i) // avoid the last bin which has all info from bins beyond
       sal<< pow(10,min_sep+(i+0.5)*delta_sep)<<" "<<dist_sep[i]/normal<<endl;
    sal.close();
#ifdef _VERBOSE_CAT_
    So.DONE();
#endif
    dist_sep.clear();
    dist_sep.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This method identifies the set of tracers around each tracer within a radius of  SCALE_MAX_N_NEIGHBOURS Mpc/h
// The output is allocated in a class member container Number_of_tracers with dimension Nobjects
void HaloTools::get_neighbour_tracers(vector<s_nearest_cells>&nearest_cells_to_cell)
{
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
#ifdef _VERBOSE_CAT_
  So.message_screen("Identifying tracers in neighbouring cells");
#endif
  real_prec Rscale = SCALE_MAX_N_NEIGHBOURS;
  // we need to cunstrunc the tree.For each grid point, we can build a vector containing the 8 or 29 nearest cells
  int max_neigh_per_dim = 1+2*N_CELLS_BACK_FORTH; // Number of neighbouring cells pr dimension, including the same itself;
  int N_Neigh_cells=pow(max_neigh_per_dim,3); // total number of cells to explore around a cell
  So.message_screen("Maximum separation probed =",  max_neigh_per_dim*sqrt(3.0)*this->box.Lbox/static_cast<real_prec>(this->box.Nft),"Mpc /h");
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
#pragma omp parallel for
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      size_t ID=this->catalogue.GridID_at(i);
      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1_at(i));
      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2_at(i));
      cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3_at(i));
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
#endif
  // This arrawy will allcopate the number of neigh or the bin in distance in which the minn distance to pair falls
  this->Number_of_neighbours.clear();
  this->Number_of_neighbours.shrink_to_fit();
  this->Number_of_neighbours.resize(this->catalogue._NOBJS(),0);
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  So.message_screen("Computing number of neighbours for each tracer");
  So.message_screen("in spheres of radius ", Rscale, "Mpc/h");
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  So.message_screen("Computing minimum separation to other tracer for each tracer");
#endif
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->catalogue.coord1_at(i);
      real_prec y_coord=this->catalogue.coord2_at(i);
      real_prec z_coord=this->catalogue.coord3_at(i);
      size_t ID=this->catalogue.GridID_at(i);
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      vector<real_prec> min_distances_v(N_Neigh_cells,0);
#endif
      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
        {
          size_t ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
          aux_min_distance_b=BIG_NUMBER;
          real_prec min_distance_b;
#endif
          for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
            {
              size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
              if(index_gal!=i) //this simply avoids to get a 0 in the distance
                {
                  real_prec distx = x_coord - (cell_info_tr[ID_NEIGH].posx_p[k]+factor_bc_x);
                  real_prec disty = y_coord - (cell_info_tr[ID_NEIGH].posy_p[k]+factor_bc_y);
                  real_prec distz = z_coord - (cell_info_tr[ID_NEIGH].posz_p[k]+factor_bc_z);
                  real_prec dist  = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
#ifdef _USE_LOG_DIST_
                  dist=log10(dist);
#endif
#if defined (_USE_NUMBER_OF_NEIGHBOURS_) || defined (_USE_LOCAL_OVERDENSITY_ )
                  if(dist<Rscale)
            this->Number_of_neighbours_at(i)++;
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
                  min_distance_b=min(aux_min_distance_b,dist);
                  aux_min_distance_b=min_distance_b;
#endif
        }
            }
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      min_distances_v[j]=aux_min_distance_b;
#endif
        }

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      aux_min_distance_a=BIG_NUMBER;
      real_prec min_distance_a;
      for(int j=0;j<N_Neigh_cells; j++)
    if(min_distances_v[j]>0)
      {
        min_distance_a=min(min_distances_v[j],aux_min_distance_a);
        aux_min_distance_a=min_distance_a;
      }
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      int da_bin = get_bin(aux_min_distance_a,MIN_OF_MIN_SEPARATION,N_BINS_MIN_DIST_TO_NEI,MAX_OF_MIN_SEPARATION/static_cast<real_prec>(N_BINS_MIN_DIST_TO_NEI),true);
      this->Number_of_neighbours_at(i)=da_bin;
#endif
    }
#ifdef _VERBOSE_CAT_
   So.DONE();
#endif

#ifdef _USE_local_overdensity_
  So.message_screen("Computiong local clustering:");
  real_prec \tExpected_number_of_tracers = 4.0*M_PI*pow(SCALE_MAX_N_NEIGHBOURS,3)*this->mean_number_density/3.0;
  So.message_screen("\tExpected number of particles in the lcoal volume:", \tExpected_number_of_tracers);
  //  ofstream tea; tea.open("local_cl.txt");
#pragma omp parallel for
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i)
    {
      real_prec local_xi=log10(static_cast<real_prec>(this->Number_of_neighbours_at(i))/\tExpected_number_of_tracers+1);
      int index_xi=get_bin(local_xi,MIN_local_overdensity,N_BINS_MIN_DIST_TO_NEI, (MAX_local_overdensity-MIN_local_overdensity)/static_cast<real_prec>(N_BINS_MIN_DIST_TO_NEI),true);
      this->Number_of_neighbours_at(i)=index_xi;
    }
  //	tea.close();
  So.DONE();
#endif
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define  _get_m_m_prob_
void HaloTools::get_distribution_min_separations(vector<s_nearest_cells>&nearest_cells_to_cell)
{
    string data=this->catalogue._type_of_object()+"_"+this->params._Name_survey()+"_Real"+to_string(this->params._realization());

  if(this->catalogue._type_of_object()=="TRACER_MOCK")
  {
    #ifdef  _USE_CWC_
  data+="cwc_";
#endif
#ifdef _USE_BIAS_OBJECT_TO_OBJECT_
  data+="indiv_bias_";
#endif
#ifdef _USE_MACH_NUMBER_
  data+="mach5_";
#endif
#ifdef  _USE_LOCAL_OVERDENSITY_
  data+="delta5_";
#endif
#ifdef _MULTISCALE_
  data+="MS_";
#endif
  }

  // we need to cunstrunc the tree.For each grid point, we can build a vector containing the 8 or 29 nearest cells
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  int max_neigh_per_dim = 2*N_CELLS_BACK_FORTH+1 ; // Number of neighbouring cells pr dimension, including the same itself;
  size_t N_Neigh_cells=pow(max_neigh_per_dim,3)-2; // total number of cells to explore around a cell
  real_prec max_dist=(0.5*max_neigh_per_dim-0.5)*sqrt(3.0)*this->box.Lbox/static_cast<real_prec>(this->box.Nft);
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring distribution of pairs for ", this->catalogue._type_of_object());
  cout<<endl;
  So.message_screen("Maximum separation probed =", max_dist,"Mpc /h");
  cout<<endl;
#endif
  struct s_cell_info{
    vector<real_prec> posx_p;
    vector<real_prec> posy_p;
    vector<real_prec> posz_p;
    vector<size_t> gal_index;
  };
  //DEFINE VECTOR TO ALLOCATE THE NEIGHBOR CELLS OF EACH CELL
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
#ifdef _VERBOSE_CAT_
  So.message_screen("Identifying tracers in cells");
#endif
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      size_t ID=this->catalogue.GridID_at(i);
      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1_at(i));
      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2_at(i));
      cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3_at(i));
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
  int N_bins_dist=50;
  int N_mass_bins=10;
  vector<vector<real_prec>> min_distance_dist(N_bins_dist,vector<real_prec>(N_mass_bins,0));
  vector<real_prec> distance_dist_m(N_bins_dist*N_mass_bins*N_mass_bins,0);
  real_prec lm_min=this->params._LOGMASSmin();
  real_prec lm_max=this->params._LOGMASSmax();
  this->logdeltaM=(lm_max-lm_min)/static_cast<real_prec>(N_mass_bins);
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
  size_t pair_counting=0;
  So.message_screen("Computing minimum separations");
  // no paralelizar si hay un min o max() dentro del loop
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->catalogue.coord1_at(i);
      real_prec y_coord=this->catalogue.coord2_at(i);
      real_prec z_coord=this->catalogue.coord3_at(i);
      real_prec lmass=log10(this->catalogue.mass_at(i))+log10(this->params._MASS_units());
      int im=get_bin(lmass,lm_min,N_mass_bins,this->logdeltaM,true);
      size_t ID=this->catalogue.GridID_at(i);

      // Allocate the min_distances to i particle form neighbouring cells
      vector<real_prec> min_distances_v(N_Neigh_cells,0);

      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
       {
          size_t ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;
          aux_min_distance_b=1e5;
          real_prec min_distance_b=0;
          for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
        {
          pair_counting++;
          size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
          real_prec lmass_k=log10(this->catalogue.mass_at(index_gal))+log10(this->params._MASS_units());
          int im_k=get_bin(lmass_k,lm_min,N_mass_bins,this->logdeltaM,true);
          int index_mass=index_2d(im,im_k,N_mass_bins);
          if(index_gal!=i) //this simply avoids to get a 0 in the distance
            {
              real_prec distx = x_coord - (cell_info_tr[ID_NEIGH].posx_p[k]+factor_bc_x);
              real_prec disty = y_coord - (cell_info_tr[ID_NEIGH].posy_p[k]+factor_bc_y);
              real_prec distz = z_coord - (cell_info_tr[ID_NEIGH].posz_p[k]+factor_bc_z);
              real_prec dist  = sqrt(pow(distx,2)+pow(disty,2)+pow(distz,2));
              min_distance_b=min(aux_min_distance_b,dist);
              aux_min_distance_b=min_distance_b;
              int da_bin = get_bin(dist,0.,N_bins_dist,max_dist/static_cast<real_prec>(N_bins_dist),false);
              size_t sindex = index_2d(da_bin,index_mass,N_mass_bins);
              if(dist<max_dist)
               distance_dist_m[sindex]++;
            }
        }// at this point min_distance_b exits as the minimum within *one* neighbour cell
          min_distances_v[j]=aux_min_distance_b;
       }// at this point min_distance_a exits as the minimum among all neighbour cells
      aux_min_distance_a=1e5;
      real_prec min_distance_a;
      for(int j=0;j<N_Neigh_cells; j++)
        {
        if(min_distances_v[j]>0)
          {
            min_distance_a=min(min_distances_v[j],aux_min_distance_a);
            aux_min_distance_a=min_distance_a;
          }
        }
      int d_bin = get_bin(min_distance_a,0.,N_bins_dist,max_dist/static_cast<real_prec>(N_bins_dist),false);
      if(min_distance_a<max_dist)
        min_distance_dist[d_bin][im]++;
    }
#ifdef _VERBOSE_CAT_
  So.DONE();
  So.message_screen("Total number of pairs", pair_counting);
#endif
  ofstream sal;
  string ffa=this->params._Output_directory()+"minimum_separation_distribution_in_mass_bins_"+data+".txt";
  sal.open(ffa.c_str());
  So.message_screen("in file ", ffa);
  for(size_t i=0;i<min_distance_dist.size();++i)
    {
      sal<<(i+0.5)*max_dist/static_cast<real_prec>(N_bins_dist)<<" ";
      for(int j=0;j< N_mass_bins; ++j)
        sal<<min_distance_dist[i][j]/static_cast<real_prec>(pair_counting)<<" ";
      sal<<endl;
    }
  sal.close();
#ifdef _get_m_m_prob_
  So.message_screen("Writting P(r|M,M') ");
  for(size_t i=0;i<min_distance_dist.size();++i)
    {
      string ff=this->params._Output_directory()+"separation_distribution_mass_"+data+"_separation_bin"+to_string(i)+".txt"+data;
      So.message_screen("in file ", ff);
      sal.open(ff.c_str());
      for(int j=0;j< N_mass_bins; ++j)
        for(int k=0;k< N_mass_bins; ++k)
          {
            int index=index_2d(j,k,N_mass_bins);
            size_t sindex = index_2d(i,index,N_mass_bins);
            sal<<lm_min+(j+0.5)*logdeltaM<<" "<<lm_min+(k+0.5)*logdeltaM<<"  "<<static_cast<real_prec>(distance_dist_m[sindex])/static_cast<real_prec>(pair_counting)<<endl;
         }
      sal.close();
    }
#endif
  So.message_screen("Writting P(r) in mass bins ");
  string ff=this->params._Output_directory()+"separation_distribution_in_mass_bins_"+data+".txt";
  So.message_screen("in file ", ff);
  sal.open(ff.c_str());
  for(size_t i=0;i<min_distance_dist.size();++i)
    {
      sal<<(i+0.5)*max_dist/static_cast<real_prec>(N_bins_dist)<<" ";
      for(int j=0;j< N_mass_bins; ++j)
        sal<<static_cast<real_prec>(distance_dist_m[index_2d(i,j,N_mass_bins)])/static_cast<real_prec>(pair_counting)<<" ";
      sal<<endl;
    }
  sal.close();
#ifdef _VERBOSE_CAT_
  So.DONE();
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// for each bin in prop, we will allocate the list of masses (M1, M2) in pairs with separatIons d1-2 between d_min and 1.5Mpc
// WE do this from the reference as well as from the final mock with assigned masses
// the index im in this->masses_in_cells_min_sep_ref[ih].M1[im] denotes the mass1 of the im-th pair in the theta bin ih separated a distance [dmin, Dmax]
// frin the other partice having mass this->masses_in_cells_min_sep_ref[ih].M2[im]
//
#ifdef _CORRECT_EXCLUSION_
void HaloTools::get_masses_of_pairs_in_min_separation_bin_in_theta_bin(real_prec min_sep, vector<s_mass_members> & dm_properties_bins)
{
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);

  So.message_screen("**Getting list of pair of masses in bins of theta with separations");
  So.message_screen("**between the the minimum of the reference and", EXCLUSION_SCALE," Mpc/h:");
  //  So.message_screen("Using ", NTHREADS, "threads");

  //#pragma omp parallel for
  for(int ih=0;ih< dm_properties_bins.size();++ih)
    {
      int nobjs=  dm_properties_bins[ih].tracer_properties.size();
      if(nobjs>=2)
    {
      for(int j=0;j< nobjs;++j)
        {
          //                size_t Id1= dm_properties_bins[ih].GridID_bin_properties[j];
          for(int k=j+1;k< nobjs;++k)
        {
          //                  size_t Id2= dm_properties_bins[ih].GridID_bin_properties[k];
                  real_prec eta_dist=  pow(dm_properties_bins[ih].x_coord_bin_properties[j]-dm_properties_bins[ih].x_coord_bin_properties[k],2);
          eta_dist+=pow(dm_properties_bins[ih].y_coord_bin_properties[j]-dm_properties_bins[ih].y_coord_bin_properties[k],2);
          eta_dist+=pow(dm_properties_bins[ih].z_coord_bin_properties[j]-dm_properties_bins[ih].z_coord_bin_properties[k],2);
          eta_dist=sqrt(eta_dist);
                    if(eta_dist>= min_sep && eta_dist<=EXCLUSION_SCALE)
                     {
                       this->masses_in_cells_min_sep[ih].M1.push_back(dm_properties_bins[ih].tracer_properties[j]);
                       this->masses_in_cells_min_sep[ih].M2.push_back(dm_properties_bins[ih].tracer_properties[k]);
                     }
                  }
              }
          }
      }

    So.DONE();
}
#endif

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function tries to answer the question:
// at which value of Vmax do we start seeing cells of a given size occupied with *one* tracer
// of that vmax value?
struct pdf_prop{
  vector<real_prec>values;
  vector<size_t>counts;
};
void HaloTools::get_pdf_vmax(string property)
{

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  size_t NGRID=this->params._Nft()*this->params._Nft()*this->params._Nft();
  vector<pdf_prop> vmax_in_cells(NGRID);
  int Nbins=30; //this->params._NMASSbins_mf();
  int Nocc=30; //this->params._NMASSbins_mf();
  if("VMAX"==property)
    {
      So.message_warning("Allcoating properties");
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
       {
         size_t ID=this->catalogue.GridID_at(i);
         vmax_in_cells[ID].values.push_back(this->catalogue.vmax_at(i));
       }
       So.DONE();
      real_prec delta=log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<real_prec>(Nbins);
      So.message_warning("Getting histograms");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i< NGRID ;++i)
       {
         vmax_in_cells[i].counts.resize(Nbins,0);
         for(size_t j=0;j<vmax_in_cells[i].values.size();++j)
          {
            real_prec value=log10(vmax_in_cells[i].values[j]);
            int bin=get_bin(value, log10(this->params._VMAXmin()),Nbins, delta,false);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
            vmax_in_cells[i].counts[bin]++;
          }
        }
        So.DONE();
     vector<real_prec>pdf_c(Nocc*Nbins,0);
     for(int k=0;k<Nocc;++k)
       for(size_t i=0;i<Nbins;++i)
         for(size_t j=0;j<NGRID;++j)
           if(k==vmax_in_cells[j].counts[i])
             pdf_c[index_2d(k,i,Nbins)]++;

     ofstream sa;
     sa.open(this->params._Output_directory()+"pdf_vmax.txt");
     for(size_t i=0;i<Nbins;++i)
       {
       sa<<this->params._VMAXmin()*pow(10,(i+0.5)*delta)<<" ";
       for(int k=1;k<Nocc;++k)
          sa<<pdf_c[index_2d(k,i,Nbins)]<<"  ";
       sa<<endl;
      }
     sa.close();
      }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined (_USE_HYBRID_ASSIGNMENT_NEW_) || defined (_USE_HYBRID_ASSIGNMENT_NEW_NEW_)
void HaloTools::select_random_subsample(real_prec fraction){
#else
void HaloTools::select_random_subsample(string output_file, bool write_to_file, real_prec fraction){
#endif
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);

    size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));

   this->catalogue_random.clear_mem();
   size_t counter=0;
    do
     {
      size_t i= gsl_rng_uniform_int(gBaseRando,this->catalogue._NOBJS());
      this->catalogue_random.push_coord1(this->catalogue.coord1_at(i));
      this->catalogue_random.push_coord2(this->catalogue.coord2_at(i));
      this->catalogue_random.push_coord3(this->catalogue.coord3_at(i));
      this->catalogue_random.push_vmax(this->catalogue.vmax_at(i));
      this->catalogue_random.push_vmax_parent(this->catalogue.vmax_parent_at(i));
      this->catalogue_random.push_GridID(this->catalogue.GridID_at(i));
      this->catalogue_random.push_galThetaID(this->catalogue.galThetaID_at(i));
      counter++;
    }while(counter<Nobjs_fraction);
    this->set_NOBJS_subsample(counter);
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::select_random_subsample(real_prec fraction, string file){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);

   size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));
   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(3);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);
   size_t counter=0;
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
    do
     {
      int i= gsl_rng_uniform_int(gBaseRando,NOBJS);
      rcat<<log10(this->catalogue.mass_at(i))<<"\t"<<log10(this->catalogue.mean_number_density_at(i))<<"\t"
        <<log10(this->catalogue.vmax_at(i))<<"\t"<<log10(this->catalogue.rs_at(i))<<"\t"<<log10(this->catalogue.concentration_at(i))<<"\t"
        <<log10(this->catalogue.spin_bullock_at(i))<<"\t"<<log10(this->catalogue.vrms_at(i))<<"\t"<<this->catalogue.virial_at(i)<<"\t"
        <<this->catalogue.b_to_a_at(i)<<"\t"<<this->catalogue.c_to_a_at(i)<<"\t"<<this->catalogue.mach_number_at(i)<<"\t"
        <<this->catalogue.bias_at(i)<<"\t"<<this->catalogue.qbias_at(i)<<"\t"<<log10(1+this->catalogue.local_overdensity_at(i))<<"\t"
        <<this->catalogue.dach_number_at(i)<<"\t"<<this->catalogue.tidal_anisotropy_at(i)<<"\t"<<this->catalogue.peak_height_at(i)<<"\t"
        <<static_cast<int>(this->catalogue.gal_cwt_at(i))<<"\t"<<this->catalogue.relative_bias_at(i)<<"\t"<<this->catalogue.bias_rs_at(i)<<"\t"
        <<this->catalogue.rs_factor_at(i)<<"\t"<<this->catalogue.local_dm_at(i)<<"\t"<<this->catalogue.lambda1_at(i)<<"\t"
        <<this->catalogue.lambda2_at(i)<<"\t"<<this->catalogue.lambda3_at(i)<<"\t"<<log10(this->catalogue.mass_closest_neighbour_at(i))
        <<"\t"<<this->catalogue.distance_closest_neighbour_at(i)<<"\t"<<log10(this->catalogue.spin_closest_neighbour_at(i))<<"\t"<<log10(this->catalogue.concentration_closest_neighbour_at(i))<<"\t"<<log10(this->catalogue.most_massive_neighbour_at(i))<<"\t"<<this->catalogue.distance_to_most_massive_neighbour_at(i)<<endl;
      counter++;
    }while(counter<Nobjs_fraction);
    rcat.close();
   So.DONE();
   So.message_screen("\tWritting downsample catalog in file", file);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::select_random_subsample(real_prec fraction){
  this->So.enter(__PRETTY_FUNCTION__);
  const gsl_rng_type * rng_t;
  gsl_rng * gBaseRando;
  gsl_rng_env_setup();
  gsl_rng_default_seed=35;
  rng_t = gsl_rng_mt19937;//_default;
  gBaseRando = gsl_rng_alloc (rng_t);

  if(fraction <= 0 || fraction > 1)
     throw std::runtime_error("Fraction must be in (0,1]");

   size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));
   size_t counter=0;
   this->catalogue_random_subsample.clear_mem();
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%, ");

   vector<size_t> rind;
   So.message_screen("\tSelecting random ids");
   do
    {
      rind.push_back(gsl_rng_uniform_int(gBaseRando,this->catalogue._NOBJS()));
      counter++;
    }while(counter<Nobjs_fraction);

   gsl_rng_free(gBaseRando);


   So.message_screen("\tSelecting properties");
   if(params._i_coord1_g()>=0)
      for(size_t j=0; j<rind.size();++j)
      {
        real_prec val = this->catalogue.coord1_at(rind[j]);
        this->catalogue_random_subsample.push_coord1(val);
      }

   if(params._i_coord2_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_coord2(this->catalogue.coord2_at(rind[j]));
   if(params._i_coord3_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_coord3(this->catalogue.coord3_at(rind[j]));
   if(params._i_v1_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_vel1(this->catalogue.vel1_at(rind[j]));
   if(params._i_v2_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_vel2(this->catalogue.vel2_at(rind[j]));
   if(params._i_v3_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_vel3(this->catalogue.vel3_at(rind[j]));
   if(params._i_mass_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_mass(this->catalogue.mass_at(rind[j]));
   if(params._i_spin_bullock_g()>0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_spin_bullock(this->catalogue.spin_bullock_at(rind[j]));
   if(params._i_rs_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_rs(this->catalogue.rs_at(rind[j]));
   if(params._i_rs_g()>=0  && params._i_rvir_g()>=0)
      for(size_t j=0; j<rind.size();++j)
        this->catalogue_random_subsample.push_concentration(this->catalogue.concentration_at(rind[j]));
        
    if(this->params._get_cwc_properties()) 
      {
        for(size_t j=0; j<rind.size();++j)
        {
          this->catalogue_random_subsample.push_GridID(this->catalogue.GridID_at(rind[j]));
          this->catalogue_random_subsample.push_gal_cwt(this->catalogue.gal_cwt_at(rind[j]));
          this->catalogue_random_subsample.push_tidal_anisotropy(this->catalogue.tidal_anisotropy_at(rind[j]));
          this->catalogue_random_subsample.push_tidal_anisotropy_dm(this->catalogue.tidal_anisotropy_dm_at(rind[j]));
        }
     }

   this->catalogue_random_subsample.set_NOBJS(counter);
   So.DONE();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::select_random_subsample_bl(real_prec fraction, string file){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);

   size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));
   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(3);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);
   size_t counter=0;
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
    do
     {
      int i= gsl_rng_uniform_int(gBaseRando,NOBJS);
      rcat<<log10(this->catalogue.mass_at(i))<<"\t"<<"\t"<<log10(this->catalogue.vmax_at(i))<<"\t"<<log10(this->catalogue.rs_at(i))<<"\t"<<log10(this->catalogue.concentration_at(i))<<"\t"<<log10(this->catalogue.spin_bullock_at(i))<<"\t"<<this->catalogue.bias_at(i)<<"\t";
      for(int il=0;il<this->catalogue.bias_multipole_size();++il)rcat<<this->catalogue.bias_multipole_at(i,il)<<"\t";
      rcat<<endl;
      counter++;
    }while(counter<Nobjs_fraction);
  rcat.close();
  So.DONE();
  gsl_rng_free(gBaseRando);
  So.message_screen("\tWritting downsample catalog in file", file);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::select_random_subsample(real_prec fraction, int Nprop, vector<real_prec>&data, string file, string file_pca){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
   size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));
   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(3);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);
   ofstream rcatp;
   rcatp.open(file_pca.c_str());
   rcatp.precision(5);
   rcatp.setf(ios::showpoint);
   rcatp.setf(ios::scientific);
   size_t counter=0;
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
   while(counter<Nobjs_fraction)
     {
      int i= gsl_rng_uniform_int(gBaseRando,NOBJS);
      // write a fraction of original catalog
      rcat<<log10(this->catalogue.mass_at(i))<<"\t"<<log10(this->catalogue.vmax_at(i))<<"\t"<<log10(this->catalogue.rs_at(i))<<"\t"<<log10(this->catalogue.concentration_at(i))<<"\t"<<log10(this->catalogue.spin_bullock_at(i))<<"\t"<<log10(this->catalogue.vrms_at(i))<<"\t"<<this->catalogue.virial_at(i)<<"\t"<<this->catalogue.b_to_a_at(i)<<"\t"<<this->catalogue.c_to_a_at(i)<<"\t"<<this->catalogue.mach_number_at(i)<<"\t"<<this->catalogue.bias_at(i)<<"\t"<<this->catalogue.qbias_at(i)<<"\t"<<log10(1+this->catalogue.local_overdensity_at(i))<<"\t"<<this->catalogue.dach_number_at(i)<<"\t"<<this->catalogue.tidal_anisotropy_at(i)<<"\t"<<this->catalogue.peak_height_at(i)<<"\t"<<static_cast<int>(this->catalogue.gal_cwt_at(i))<<"\t"<<this->catalogue.relative_bias_at(i)<<"\t"<<this->catalogue.bias_rs_at(i)<<"\t"<<this->catalogue.rs_factor_at(i)<<"\t"<<this->catalogue.local_dm_at(i)<<"\t"<<this->catalogue.lambda1_at(i)<<"\t"<<this->catalogue.lambda2_at(i)<<"\t"<<this->catalogue.lambda3_at(i)<<"\t"<<log10(this->catalogue.mass_closest_neighbour_at(i))<<"\t"<<this->catalogue.distance_closest_neighbour_at(i)<<"\t"<<log10(this->catalogue.spin_closest_neighbour_at(i))<<"\t"<<log10(this->catalogue.concentration_closest_neighbour_at(i))<<"\t"<<log10(this->catalogue.most_massive_neighbour_at(i))<<"\t"<<this->catalogue.distance_to_most_massive_neighbour_at(i)<<endl;
      // write a fraction  of pca properties. Only 10 properties
      for(int j=0;j<Nprop;++j)
          rcatp<<data[index_2d(j,i,this->catalogue._NOBJS())]<<"\t";
      rcatp<<"\t"<<static_cast<int>(this->catalogue.gal_cwt_at(i))<<endl;
      counter++;
   }
   So.DONE();
   So.message_screen("\tDownsampled catalog written in file", file);
   So.message_screen("\tDownsampled PCA catalog with 10 PCA written in file", file_pca);
   rcat.close();
   rcatp.close();
   gsl_rng_free(gBaseRando);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::select_random_subsample_v(real_prec fraction){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
   size_t Nobjs_fraction=static_cast<size_t>(floor(fraction*this->catalogue._NOBJS()));
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
   So.message_screen("\tNumber of selected tracers: ",Nobjs_fraction);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<NOBJS;++i)
    this->catalogue.observed_at(i)=0;
   
    size_t counter=0;
   while(counter<Nobjs_fraction)
     {
       int i= gsl_rng_uniform_int(gBaseRando,NOBJS);
       if(0==this->catalogue.observed_at(i))// not to repeat
         {
          this->catalogue.set_observed(1,i);
          counter++;
        }
    }
   So.DONE();
   So.message_screen("\tDownsampled catalog allocated as observed");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::halos2galaxies_HOD()
{
 this->So.enter(__PRETTY_FUNCTION__);
    real_prec hod=0;
    // and this->gal the halo catalog
    // SET HOD parameters
    real_prec Mc=1e12;
    real_prec sigmaM=0.12;
    real_prec Mcut=1e12;
    real_prec alpha=1.10;
    // loop over the tracers to select central galaxies
    const gsl_rng_type *Tn;
    gsl_rng *rng ;
    Tn = gsl_rng_default;
    rng = gsl_rng_alloc (Tn);
    vector<real_prec>xnew;
    vector<real_prec>ynew;
    vector<real_prec>znew;
    vector<real_prec>vxnew;
    vector<real_prec>vynew;
    vector<real_prec>vznew;
    vector<real_prec>massg;
    vector<real_prec>vmaxg;
   for(size_t i=0;i<this->catalogue._NOBJS();++i)
    {
      real_prec mass=this->catalogue.mass_at(i);
      real_prec xp= log10(Mc/mass)/sigmaM;
      real_prec prob=gsl_sf_gamma(xp);
      real_prec xr=gsl_rng_uniform(rng);
      if(prob<xr) // select centrals
      {
          xnew.push_back(this->catalogue.coord1_at(i));
          ynew.push_back(this->catalogue.coord2_at(i));
          znew.push_back(this->catalogue.coord3_at(i));
          vxnew.push_back(this->catalogue.vel1_at(i));
          vynew.push_back(this->catalogue.vel2_at(i));
          vznew.push_back(this->catalogue.vel3_at(i));
          massg.push_back(this->catalogue.mass_at(i));
          vmaxg.push_back(this->catalogue.vmax_at(i));
      }
    }

    catalogue_central.resize_positions(vmaxg.size());
    catalogue_central.resize_velocities(vmaxg.size());
    catalogue_central.resize_mass(vmaxg.size());
    catalogue_central.resize_vmax(vmaxg.size());

    #ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<vmaxg.size();++i)
    {
      catalogue_central.set_coord1(xnew[i],i);
      catalogue_central.set_coord2(ynew[i],i);
      catalogue_central.set_coord3(znew[i],i);
      catalogue_central.set_vel1(vxnew[i],i);
      catalogue_central.set_vel2(vynew[i],i);
      catalogue_central.set_vel3(vznew[i],i);
      catalogue_central.set_mass(massg[i],i);
      catalogue_central.set_vmax(vmaxg[i],i);
    }
    size_t Ncentrals=vmaxg.size();
    xnew.clear();xnew.shrink_to_fit();
    ynew.clear();ynew.shrink_to_fit();
    znew.clear();znew.shrink_to_fit();
    vxnew.clear();vxnew.shrink_to_fit();
    vynew.clear();vxnew.shrink_to_fit();
    vznew.clear();vznew.shrink_to_fit();
    massg.clear();massg.shrink_to_fit();
    vmaxg.clear();vmaxg.shrink_to_fit();
    real_prec rmin=0.01; // minium radius used for normalization of NFR profile
    vector<size_t> Nsats_per_central(Ncentrals,0);
    DensityProfiles DenProf(this->s_cosmo_pars);

            // lopp over the centrals
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i<Ncentrals;++i)
          Nsats_per_central[i]=static_cast<size_t>(gsl_ran_poisson(rng,floor(pow(catalogue_central.mass_at(i) / Mcut , alpha))));  // draw number of satellites per central:

    size_t counter_sat=0;

    for(size_t i=0;i<Ncentrals;++i)
        for(size_t j=0;i<Nsats_per_central[i];++j)
        {
            real_prec phi=2.*M_PI*gsl_rng_uniform(rng);
            real_prec theta=acos(-1.0+2.*gsl_rng_uniform(rng));
            real_prec prob=-10.;
            real_prec radius=0.0;
            real_prec xra=0;
            DenProf.nfw_parameters(catalogue_central.mass_at(i));
            do {
                radius=DenProf._rvir()*pow(gsl_rng_uniform(rng),1./3.);
                prob=DenProf.DensityProfile_NFW_prob(radius,rmin);
                xra=gsl_rng_uniform(rng);
            }while(prob<xra);
            catalogue_satellite.set_coord1(radius*sin(theta)*cos(phi)+catalogue_central.coord1_at(i), counter_sat);
            catalogue_satellite.set_coord2(radius*cos(theta)*sin(phi)+catalogue_central.coord2_at(i), counter_sat);
            catalogue_satellite.set_coord3(radius*cos(phi)+catalogue_central.coord3_at(i), counter_sat);
            catalogue_satellite.set_vel1(catalogue_central.vel1_at(i)+gsl_ran_gaussian(rng,catalogue_central.vmax_at(i)), counter_sat) ;
            catalogue_satellite.set_vel2(catalogue_central.vel2_at(i)+gsl_ran_gaussian(rng,catalogue_central.vmax_at(i)), counter_sat) ;
            catalogue_satellite.set_vel3(catalogue_central.vel3_at(i)+gsl_ran_gaussian(rng,catalogue_central.vmax_at(i)),counter_sat) ;
            catalogue_satellite.set_mass(catalogue_central.mass_at(i), counter_sat) ;
            counter_sat++;
        }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 real_prec HaloTools::get_min(string prop){
  if(this->catalogue._NOBJS()==0)
     cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;

  real_prec lka=static_cast<real_prec>(LARGE_NUMBER);
  if(prop=="_MASS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.mass_at(i),lka);
    }
  else if(prop=="_VMAX_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
        lka=std::min(this->catalogue.vmax_at(i),lka);
    }
  else if(prop=="_RS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.rs_at(i),lka);
    }
  else if(prop=="_RVIR_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.rvir_at(i),lka);
    }
  else if(prop==_CONCENTRATION_)
    {
#pragma omp parallel for reduction(min:lka)
    for(size_t i=0;i<NOBJS;++i)
      lka=min(this->catalogue.concentration_at(i),lka);
    }
  else if(prop=="_SPIN_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.spin_at(i),lka);
    }
  else if(prop=="_SPIN_BULLOCK_")
    {
#pragma omp parallel for reduction(min:lka)
     for(size_t i=0;i<NOBJS;++i)
        lka=min(this->catalogue.spin_bullock_at(i),lka);
    }
  else if(prop== "_VIRIAL_")
    {
#pragma omp parallel for reduction(min:lka)
    for(size_t i=0;i<NOBJS;++i)
      lka=min(this->catalogue.virial_at(i),lka);
    }
  else if(prop== "_VRMS_")
    {
#pragma omp parallel for reduction(min:lka)
    for(size_t i=0;i<NOBJS;++i)
      lka=min(this->catalogue.vrms_at(i),lka);
    }
  else if(prop=="_BTOA_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.b_to_a_at(i),lka);
    }
  else if(prop=="_CTOA_")
    {
#pragma omp parallel for reduction(min:lka)
     for(size_t i=0;i<NOBJS;++i)
       lka=min(this->catalogue.c_to_a_at(i),lka);
    }
  else if(prop=="_MACH_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.mach_number_at(i),lka);
    }
  else if(prop=="_DACH_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.dach_number_at(i),lka);
    }
  else if(prop=="_BIAS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.bias_at(i),lka);
    }
  else if(prop=="_RELATIVE_BIAS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.relative_bias_at(i),lka);
    }
  else if(prop=="_LOCAL_OVERDENSITY_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.local_overdensity_at(i),lka);
    }
  else if(prop=="_TIDAL_ANISOTROPY_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.tidal_anisotropy_at(i),lka);
    }
  else if(prop=="_LOCALDM_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.local_dm_at(i),lka);
    }
  else if(prop=="_PEAK_HEIGHT_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.peak_height_at(i),lka);
    }
  else if(prop=="_XCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.coord1_at(i),lka);
    }
  else if(prop=="_YCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.coord2_at(i),lka);
    }
  else if(prop=="_ZCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.coord3_at(i),lka);
    }
  else if(prop=="_REDSHIFT_")
    {
#pragma omp parallel for reduction(min:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=min(this->catalogue.redshift_at(i),lka);
  }  
   
  else{
    this->So.message_screen("Property not found");
  }
  return lka;
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 real_prec HaloTools::get_max(string prop){
 if(this->catalogue._NOBJS()==0)
    cerr<<"Error. Array Halo is empty, in line "<<__PRETTY_FUNCTION__<<endl;
  real_prec lka=-static_cast<real_prec>(LARGE_NUMBER);
  real_prec lkb;
  if(prop=="_MASS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.mass_at(i),lka);
    }
  else if(prop=="_VMAX_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.vmax_at(i),lka);
    }
  else if(prop=="_RS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.rs_at(i),lka);
    }
  else if(prop=="_RVIR_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.rvir_at(i),lka);
    }
  else if(prop=="_CONCENTRATION_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.concentration_at(i),lka);
    }
  else if(prop=="_SPIN_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.spin_at(i),lka);
    }
  else if(prop=="_SPIN_BULLOCK_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.spin_bullock_at(i),lka);
    }
  else if(prop== "_VIRIAL_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.virial_at(i),lka);
    }
  else if(prop== "_VRMS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.vrms_at(i),lka);
    }
  else if(prop=="_BTOA_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.b_to_a_at(i),lka);
    }
   else if(prop=="_CTOA_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.c_to_a_at(i),lka);
     }
   else if(prop=="_MACH_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.mach_number_at(i),lka);
     }
  else if(prop=="_DACH_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
    lka=max(this->catalogue.dach_number_at(i),lka);
    }
   else if(prop=="_BIAS_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.bias_at(i),lka);
     }
   else if(prop=="_RELATIVE_BIAS_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.relative_bias_at(i),lka);
     }
   else if(prop=="_LOCAL_OVERDENSITY_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.local_overdensity_at(i),lka);
     }
   else if(prop=="_TIDAL_ANISOTROPY_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
       lka=max(this->catalogue.tidal_anisotropy_at(i),lka);
     }
  else if(prop=="_LOCALDM_")
    {
#pragma omp parallel for reduction(max:lka)
      for(size_t i=0;i<NOBJS;++i)
      lka=max(this->catalogue.local_dm_at(i),lka);
    }
   else if(prop=="_PEAK_HEIGHT_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.peak_height_at(i),lka);
     }
   else if(prop=="_XCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.coord1_at(i),lka);
     }
   else if(prop=="_YCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.coord2_at(i),lka);
     }
   else if(prop=="_ZCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.coord3_at(i),lka);
     }
   else if(prop=="_REDSHIFT_")
     {
#pragma omp parallel for reduction(max:lka)
       for(size_t i=0;i<NOBJS;++i)
     lka=max(this->catalogue.redshift_at(i),lka);
     }
   else
     this->So.message_screen("Property not found");

  return lka;
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  real_prec HaloTools::get_variance(real_prec meanp, string prop){
  if(NOBJS==0)
     cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;
   real_prec var=0;
  if(prop=="_VMAX_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(this->catalogue.vmax_at(i)-meanp,2);
   else if(prop=="_MASS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(log10(this->catalogue.mass_at(i))-meanp,2);
  else if(prop=="_RS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(this->catalogue.rs_at(i)-meanp,2);
  else if(prop=="_RVIR_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(this->catalogue.rvir_at(i)-meanp,2);
  else if(prop=="_CONCENTRATION_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(log10(this->catalogue.concentration_at(i))-meanp,2);
  else if(prop=="_SPIN_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(this->catalogue.spin_at(i)-meanp,2);
  else if(prop=="_SPIN_BULLOCK")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(size_t i=0;i<NOBJS;++i)
       var+=pow(this->catalogue.spin_bullock_at(i)-meanp,2);
   else if(prop== "_VIRIAL_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(size_t i=0;i<NOBJS;++i)
          var+=pow(this->catalogue.virial_at(i)-meanp,2);
   else if(prop== "_VRMS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(size_t i=0;i<NOBJS;++i)
          var+=pow(this->catalogue.vrms_at(i)-meanp,2);
   else if(prop=="_BTOA_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
          var+=pow(this->catalogue.b_to_a_at(i)-meanp,2);
   else if(prop=="_CTOA_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(size_t i=0;i<NOBJS;++i)
          var+=pow(this->catalogue.c_to_a_at(i)-meanp,2);
   else if(prop=="_MACH_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
        var+=pow(this->catalogue.mach_number_at(i)-meanp,2);
   else if(prop=="_BIAS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
        var+=pow(this->catalogue.bias_at(i)-meanp,2);
   else if(prop=="_RELATIVE_BIAS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
        var+=pow(this->catalogue.relative_bias_at(i)-meanp,2);
   else if(prop=="_LOCAL_OVERDENSITY_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
        var+=pow(this->catalogue.local_overdensity_at(i)-meanp,2);
  else if(prop=="_TIDAL_ANISOTROPY_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
      var+=pow(this->catalogue.tidal_anisotropy_at(i)-meanp,2);
  }
  else if(prop=="_LOCALDM_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
      var+=pow(this->catalogue.local_dm_at(i)-meanp,2);
  }
    else if(prop=="_PEAK_HEIGHT_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(size_t i=0;i<NOBJS;++i)
        var+=pow(this->catalogue.peak_height_at(i)-meanp,2);
    var/=static_cast<real_prec>(NOBJS);
  }
  else if(prop=="_DACH_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
  for(size_t i=0;i<NOBJS;++i)
      var+=pow(this->catalogue.dach_number_at(i)-meanp,2);
  var/=static_cast<real_prec>(NOBJS);
}
  else
      this->So.message_screen("Property not found");
  return var;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 pair<real_prec,real_prec> HaloTools::get_variance(string prop,bool factor){
     real_prec fac=0;
     if(factor==true)
         fac=1;
   if(NOBJS==0)
      cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;
    real_prec var=0;
    real_prec meanp=0;
   if(prop=="_VMAX_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
      for(size_t i=0;i<NOBJS;++i)
         meanp+=log10(this->catalogue.vmax_at(i));
      meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(size_t i=0;i<NOBJS;++i)
          var+=pow(log10(this->catalogue.vmax_at(i))-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
    }
    else if(prop=="_MASS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
      for(size_t i=0;i<NOBJS;++i)
         meanp+=log10(this->catalogue.mass_at(i));
      meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(size_t i=0;i<NOBJS;++i)
          var+=pow(log10(this->catalogue.mass_at(i))-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
    }
   else if(prop=="_RS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=log10(this->catalogue.rs_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.rs_at(i))-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_RVIR_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=log10(this->catalogue.rvir_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.rvir_at(i))-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_CONCENTRATION_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(log10(this->catalogue.concentration_at(i)));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.concentration_at(i))-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_SPIN_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=log10(this->catalogue.spin_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.spin_at(i))-fac*meanp,2);
     var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_SPIN_BULLOCK_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=log10(this->catalogue.spin_bullock_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.spin_bullock_at(i))-fac*meanp,2);
     var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_VRMS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=log10(this->catalogue.vrms_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(log10(this->catalogue.vrms_at(i))-fac*meanp,2);
     var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_VIRIAL_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.virial_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.virial_at(i)-fac*meanp,2);
      var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_BTOA_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.b_to_a_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.b_to_a_at(i)-fac*meanp,2);
     var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_CTOA_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.c_to_a_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.c_to_a_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_MACH_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.mach_number_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.mach_number_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_BIAS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.bias_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.bias_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_RELATIVE_BIAS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.relative_bias_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.relative_bias_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_LOCAL_OVERDENSITY_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.local_overdensity_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.local_overdensity_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_TIDAL_ANISOTROPY_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.tidal_anisotropy_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.tidal_anisotropy_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_LOCALDM_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.local_dm_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.local_dm_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_PEAK_HEIGHT_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)
        meanp+=(this->catalogue.peak_height_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.peak_height_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
   else if(prop=="_DACH_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(size_t i=0;i<NOBJS;++i)

        meanp+=(this->catalogue.dach_number_at(i));
     meanp/=static_cast<real_prec>(NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(size_t i=0;i<NOBJS;++i)
         var+=pow(this->catalogue.dach_number_at(i)-fac*meanp,2);
       var/=static_cast<real_prec>(NOBJS);
   }
    else
   this->So.message_screen("Property not found: ", prop);

   return std::make_pair(meanp,var);
   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_intervals_equal_number(string prop,vector<real_prec>&min_aux,vector<real_prec>&max_aux){

  this->So.enter(__PRETTY_FUNCTION__);

  size_t Ntracers=this->catalogue._NOBJS();
  int Nbins_prop=1;
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop!="_PEAK_HEIGHT_")
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers();
   else
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers_main_property();
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop!="_MASS_")
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers();
   else
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers_main_property();
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop!="_VMAX_")
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers();
   else
       Nbins_prop= this->params._Number_of_bins_equal_number_tracers_main_property();
#endif
   size_t Ntracers_bin=static_cast<size_t>(floor(Ntracers/Nbins_prop));
   So.message_screen("\tExpected (equal) number of tracers in bins of sec property =", Ntracers_bin);
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_PEAK_HEIGHT_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_ph_bins[0]=this->catalogue._NOBJS();
     this->Number_of_tracers_in_ph_bins[Nbins_prop]=this->catalogue._NOBJS()-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_ph_bins[i]=Ntracers_bin;
   }
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_MASS_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_mass_bins[0]=this->catalogue._NOBJS();
     this->Number_of_tracers_in_mass_bins[Nbins_prop]=this->catalogue._NOBJS()-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_mass_bins[i]=Ntracers_bin;
   }
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_VMAX_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_vmax_bins[0]=this->catalogue._NOBJS();
     this->Number_of_tracers_in_vmax_bins[Nbins_prop]=this->catalogue._NOBJS()-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_vmax_bins[i]=Ntracers_bin;
   }
#endif
   size_t NBins=100000;
   real_prec minp=get_min(prop);
   real_prec maxp=get_max(prop);
   if(prop=="_MASS_")
   {
      this->min_mass=minp;
      this->max_mass=maxp;
    }
   else if(prop=="_VMAX_")
    {
       this->min_vmax=minp;
       this->max_vmax=maxp;
   }
   else if(prop=="_RS_"){
      this->min_rs=minp;
      this->max_rs=maxp;
    }
   else if(prop=="_RVIR_"){
      this->min_rvir=minp;
      this->max_rvir=maxp;
    }
   else if(prop=="_CONCENTRATION_"){
      this->min_concentration=minp;
      this->max_concentration=maxp;
    }
   else if(prop=="_SPIN_")
   {
      this->min_spin=minp;
      this->max_spin=maxp;
    }
   else if(prop=="_SPIN_BULLOCK_")
   {
      this->min_spin_bullock=minp;
      this->max_spin_bullock=maxp;
    }
   else if(prop=="_VRMS_")
   {
      this->min_vrms=minp;
      this->max_vrms=maxp;
    }
   else if(prop=="_VIRIAL_"){
       this->min_virial=minp;
       this->max_virial=maxp;
   }
   else if(prop=="_BTOA_")
   {
      this->min_b_to_a=minp ;
      this->max_b_to_a=maxp;
   }
   else if(prop=="_CTOA_")
   {
      this->min_c_to_a =minp;
      this->max_c_to_a=maxp;
   }
   else if(prop=="_MACH_")
   {
      this->min_mach=minp ;
      this->max_mach=maxp;
   }
   else if(prop=="_BIAS_")
   {
      this->min_bias=minp ;
      this->max_bias=maxp;
   }
   else if(prop=="_RELATIVE_BIAS_")
   {
      this->min_bias=minp ;
      this->max_bias=maxp;
   }
   else if(prop=="_LOCAL_OVERDENSITY_")
   {
      this->min_local_overdensity=minp ;
      this->max_local_overdensity=maxp;
   }
   else if(prop=="_TIDAL_ANISOTROPY_")
   {
      this->min_tidal_anisotropy=minp ;
      this->max_tidal_anisotropy=maxp;
   }
   else if(prop=="_LOCALDM_")
   {
      this->min_local_dm=minp ;
      this->max_local_dm=maxp;
   }
   else if(prop=="_PEAK_HEIGHT_")
   {
      this->min_peak_height=minp;
      this->max_peak_height=maxp;
   }
   else if(prop=="_DACH_")
   {
      minp=this->min_dach=minp  ;
      maxp=this->max_dach=maxp;
   }
   else
       this->So.message_screen("Property not found");

   real_prec delta=(maxp-minp)/static_cast<real_prec>(NBins);
   vector<size_t>counts(NBins,0);
   if(prop=="_VMAX_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.vmax_at(i),minp,NBins,delta, true)]++;
}
   else if(prop=="_RS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.rs_at(i),minp,NBins,delta, true)]++;
  }
   else if(prop=="_RVIR_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.rvir_at(i),minp,NBins,delta, true)]++;
  }
   else if(prop=="_CONCENTRATION_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.concentration_at(i),minp,NBins,delta, true)]++;
  }
   else if(prop=="_SPIN_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.spin_at(i),minp,NBins,delta, true)]++;
}
   else if(prop=="_SPIN_BULLOCK_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.spin_bullock_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_VRMS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.vrms_at(i),minp,NBins,delta, true)]++;
 }

     else if(prop=="_VIRIAL_")
 {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[ get_bin(this->catalogue.virial_at(i),minp,NBins,delta, true)]++;
  }
     else if(prop=="_BTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.b_to_a_at(i),minp,NBins,delta,true)]++;
}
     else if(prop=="_CTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.c_to_a_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_MACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.mach_number_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.bias_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_RELATIVE_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.relative_bias_at(i), minp,NBins,delta, true)]++;
}

     else if(prop=="_LOCAL_OVERDENSITY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.local_overdensity_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_PEAK_HEIGHT_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.peak_height_at(i),minp,NBins,delta, true)]++;
}
   else if(prop=="_DACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.dach_number_at(i),minp,NBins,delta, true)]++;
}
   else if(prop=="_TIDAL_ANISOTROPY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.tidal_anisotropy_at(i),minp,NBins,delta, true)]++;
}
   else if(prop=="_LOCALDM_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.local_dm_at(i),minp,NBins,delta, true)]++;
}
     else if(prop=="_MASS_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.mass_at(i),minp,NBins,delta, true)]++;
    }
   min_aux[0]=minp;
   max_aux[0]=maxp;
   min_aux[1]=minp;
   max_aux[Nbins_prop]=maxp;
   size_t indi=0;
   for(int j=1;j<Nbins_prop;++j)// Go through the 4 bins for quartiles
     {
       real_prec new_counts=0;
       size_t i=indi;
       while(new_counts<=Ntracers_bin) // Go through the small bins untill we get the number of tracers expected for the quartile
         {
           new_counts+=counts[i];
           indi=i;
           i++;
         }
          max_aux[j]=minp+static_cast<real_prec>(indi-1)*delta;
          min_aux[j+1]=minp+static_cast<real_prec>(indi)*delta;
       }
 }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void HaloTools::get_intervals_equal_number_aux(string prop){
   this->So.enter(__PRETTY_FUNCTION__);
   int Nbins_prop= this->params._Number_of_bins_equal_number_tracers();
   size_t NBins=50000;
   size_t Ntracers=this->catalogue._NOBJS();
   real_prec minp=get_min(prop);
   real_prec maxp=get_max(prop);
   So.message_screen("Type of object:", this->catalogue._type_of_object());
   So.message_screen("Maximum of ", prop, "=", maxp);
   So.message_screen("Minimum of ", prop, "=", minp);
   if (std::isinf(maxp))
       maxp=LARGE_NUMBER;

   if(prop=="_MASS_")
   {
      this->min_mass=minp ;
      this->max_mass=maxp;
   }
   else if(prop=="_VMAX_")
     {
       this->min_vmax=minp;
       this->max_vmax=maxp;
     }
   else if(prop=="_RS_")
     {
       this->min_rs=minp;
       this->max_rs=maxp;
     }
   else if(prop=="_RVIR_")
   {
     this->min_rs=minp;
     this->max_rs=maxp;
   }
   else if(prop=="_CONCENTRATION_")
   {
     this->min_concentration=minp;
     this->max_concentration=maxp;
   }
   else if(prop=="_SPIN_")
     {
       this->min_spin=minp;
       this->max_spin=maxp;
     }
   else if(prop=="_SPIN_BULLOCK_")
     {
       this->min_spin_bullock=minp;
       this->max_spin_bullock=maxp;
     }
   else if(prop=="_VRMS_")
     {
       this->min_vrms=minp;
       this->max_vrms=maxp;
     }
   else if(prop=="_VIRIAL_")
     {
       this->min_virial=minp;
       this->max_virial=maxp;
     }
   else if(prop=="_BTOA_")
     {
       this->min_b_to_a=minp;
       this->max_b_to_a=maxp;
     }
   else if(prop=="_CTOA_")
     {
       this->min_c_to_a=minp;
       this->max_c_to_a=maxp;
     }
   else if(prop=="_MACH_")
     {
       this->min_mach=minp;
       this->max_mach=maxp;
     }
   else if(prop=="_BIAS_")
     {
       this->min_bias=minp;
       this->max_bias=maxp;
     }
   else if(prop=="_RELATIVE_BIAS_")
     {
       this->min_bias=minp;
       this->max_bias=maxp;
     }
   else if(prop=="_LOCAL_OVERDENSITY_")
     {
       this->min_local_overdensity=minp;
       this->max_local_overdensity=maxp;
     }
   else if(prop=="_TIDAL_ANISOTROPY_")
   {
      this->min_tidal_anisotropy=minp;
      this->max_tidal_anisotropy=maxp;
   }
   else if(prop=="_LOCALDM_")
   {
      this->min_local_dm=minp;
      this->max_local_dm=maxp;
   }
   else if(prop=="_PEAK_HEIGHT_")
   {
      this->min_peak_height=minp ;
      this->max_peak_height=maxp;
   }
   else if(prop=="_DACH_")
    {
       this->min_dach=minp;
       this->max_dach=maxp;
    }
   real_prec delta=(maxp-minp)/static_cast<real_prec>(NBins);
   vector<size_t>counts(NBins,0);
   if(prop=="_VMAX_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.vmax_at(i),minp,NBins,delta, true)]++;
    }
   else if(prop=="_RS_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.rs_at(i),minp,NBins,delta, true)]++;
    }
  else if(prop=="_RVIR_")
   {
   #ifdef _USE_OMP_
   #pragma omp parallel for
   #endif
      for(size_t i=0;i<Ntracers;++i)
   #ifdef _USE_OMP_
   #pragma atomic
   #endif
        counts[get_bin(this->catalogue.rvir_at(i),minp,NBins,delta, true)]++;
   }
   else if(prop==_CONCENTRATION_)
   {
   #ifdef _USE_OMP_
   #pragma omp parallel for
   #endif
      for(size_t i=0;i<Ntracers;++i)
   #ifdef _USE_OMP_
   #pragma atomic
   #endif
        counts[get_bin(this->catalogue.concentration_at(i),minp,NBins,delta, true)]++;
   }
  else if(prop=="_SPIN_") 
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.spin_at(i),minp,NBins,delta, true)]++;
  }
   else if(prop=="_SPIN_BULLOCK_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
   counts[get_bin(this->catalogue.spin_bullock_at(i),minp,NBins,delta, true)]++;
  }
  else if(prop=="_VRMS_")
  {
    for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
      counts[get_bin(this->catalogue.vrms_at(i),minp,NBins,delta, true)]++;
  }       
  else if(prop=="_VIRIAL_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[ get_bin(this->catalogue.virial_at(i),minp,NBins,delta, true)]++;
 }
     else if(prop=="_BTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.b_to_a_at(i),minp,NBins,delta,true)]++;
 }
     else if(prop=="_CTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.c_to_a_at(i),minp,NBins,delta, true)]++;
 }
     else if(prop=="_MACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.mach_number_at(i),minp,NBins,delta, true)]++;
 }
     else if(prop=="_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.bias_at(i),minp,NBins,delta, true)]++;
 }
      else if(prop=="_RELATIVE_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.relative_bias_at(i),minp,NBins,delta, true)]++;
 }
     else if(prop=="_LOCAL_OVERDENSITY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.local_overdensity_at(i),minp,NBins,delta, true)]++;
 }
     }
   else if(prop=="_TIDAL_ANISOTROPY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.tidal_anisotropy_at(i),minp,NBins,delta, true)]++;
     }
}

   else if(prop=="_LOCALDM_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.local_dm_at(i),minp,NBins,delta, true)]++;
     }
}
   else if(prop=="_PEAK_HEIGHT_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.peak_height_at(i),minp,NBins,delta, true)]++;
     }
 }
   else if(prop=="_DACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->catalogue.dach_number_at(i),minp,NBins,delta, true)]++;
}

     else if(prop=="_MASS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->catalogue.mass_at(i),minp,NBins,delta, true)]++;
 }
  // Get the intervals
   size_t Ntracers_bin=static_cast<size_t>(floor(Ntracers/Nbins_prop));

   vector<real_prec>min_aux(Nbins_prop+1,0);
   vector<real_prec>max_aux(Nbins_prop+1,0);
   min_aux[0]=minp;
   max_aux[0]=maxp;
   min_aux[1]=minp;
   max_aux[Nbins_prop]=maxp;
   size_t indi=0;
   for(int j=1;j<Nbins_prop;++j)
     {
       real_prec new_counts=0;
       size_t i=indi;
       while(new_counts<Ntracers_bin)
         {
           new_counts+=counts[i];
           indi=i;
           i++;
         }
       max_aux[j]=minp+static_cast<real_prec>(indi-1)*delta;
       min_aux[j+1]=minp+static_cast<real_prec>(indi)*delta;
     }

   if(prop=="_MASS_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_MASSbins_min(i,log10(min_aux[i]));
       this->params.set_MASSbins_max(i,log10(max_aux[i]));
     }
   else if(prop=="_VMAX_")
     for(size_t i=0; i<Nbins_prop+1;++i)
       {
     this->params.set_VMAXbins_min(i,log10(min_aux[i]));
     this->params.set_VMAXbins_max(i,log10(max_aux[i]));
       }
    else if(prop=="_RS_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_RSbins_min(i,min_aux[i]);
       this->params.set_RSbins_max(i,max_aux[i]);
     }
   else if(prop=="_RVIR_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_RVIRbins_min(i,min_aux[i]);
       this->params.set_RVIRbins_max(i,max_aux[i]);
     }
   else if(prop=="_CONCENTRATION_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_CONCENTRATIONbins_min(i,min_aux[i]);
       this->params.set_CONCENTRATIONbins_max(i,max_aux[i]);
     }
  else  if(prop=="_SPIN_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_SPINbins_min(i,min_aux[i]);
       this->params.set_SPINbins_max(i,max_aux[i]);
     }
    else if(prop=="_SPIN_BULLOCK_")
       for(size_t i=0; i<Nbins_prop+1;++i){
         this->params.set_SPINBULLOCKbins_min(i,min_aux[i]);
         this->params.set_SPINBULLOCKbins_max(i,max_aux[i]);
       }
    else if(prop=="_VIRIAL_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_VIRIALbins_min(i,min_aux[i]);
       this->params.set_VIRIALbins_max(i,max_aux[i]);
     }
   if(prop=="_VRMS_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_VRMSbins_min(i,min_aux[i]);
       this->params.set_VRMSbins_max(i,max_aux[i]);
     }
   if(prop=="_BTOA_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_BTOAbins_min(i,min_aux[i]);
       this->params.set_BTOAbins_max(i,max_aux[i]);
     }

   if(prop=="_CTOA_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_CTOAbins_min(i,min_aux[i]);
       this->params.set_CTOAbins_max(i,max_aux[i]);
     }

   if(prop=="_MACH_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_MACHbins_min(i,min_aux[i]);
       this->params.set_MACHbins_max(i,max_aux[i]);
     }
   if(prop=="_BIAS_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_BIASbins_min(i,min_aux[i]);
       this->params.set_BIASbins_max(i,max_aux[i]);
     }
   if(prop=="_RELATIVE_BIAS_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_RBIASbins_min(i,min_aux[i]);
       this->params.set_RBIASbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCAL_OVERDENSITY_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_LCbins_min(i,min_aux[i]);
       this->params.set_LCbins_max(i,max_aux[i]);
     }
   if(prop=="_TIDAL_ANISOTROPY_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_TAbins_min(i,min_aux[i]);
       this->params.set_TAbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCALDM_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_LOCALDMbins_min(i,min_aux[i]);
       this->params.set_LOCALDMbins_max(i,max_aux[i]);
     }
   if(prop=="_PEAK_HEIGHT_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_PHbins_min(i,min_aux[i]);
       this->params.set_PHbins_max(i,max_aux[i]);
     }
   if(prop=="_DACH_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_DACHbins_min(i,min_aux[i]);
       this->params.set_DACHbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCALDM_")
     for(size_t i=0; i<Nbins_prop+1;++i){
       this->params.set_LOCALDMbins_min(i,min_aux[i]);
       this->params.set_LOCALDMbins_max(i,max_aux[i]);
     }
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 struct proper{
    vector<real_prec> prop;
    string prop_name;

 };
void HaloTools::get_local_mach_number(bool write)
{
  So.enter(__PRETTY_FUNCTION__);
  int nprops=7;
  vector<proper>hprop(nprops);
  int IC=0;
  hprop[IC].prop_name="COUNTS";
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
  if(this->params._i_mass_g()>0)
    {
      IC++;
      hprop[IC].prop_name="MASS";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.mass_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  if(this->params._i_vmax_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VMAX";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.vmax_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  if(this->params._i_vrms_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VRMS";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.vrms_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  if(this->params._i_virial_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VIRIAL";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.virial_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  if(this->params._i_spin_g()>0)
    {
      IC++;
      hprop[IC].prop_name="SPIN";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.spin_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  if(this->params._i_rs_g()>0)
    {
      IC++;
      hprop[IC].prop_name="RS";
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->catalogue.rs_at(i)*sqrt(pow(this->catalogue.vel1_at(i),2)+pow(this->catalogue.vel2_at(i),2)+pow(this->catalogue.vel3_at(i),2)));
    }

  for(int iprop=0;iprop<nprops;++iprop)
    {
      vector<s_cell_info> cell_info_tr(this->box.NGRID);
      for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    cell_info_tr[this->catalogue.GridID_at(i)].property.push_back(hprop[iprop].prop[i]);

  So.DONE();

  this->mach_number.resize(this->box.NGRID,0);
  for(size_t i=0;i<this->box.NGRID ;++i) //loop over the "observed obejcts", i.e, with cuts already set
  {
    real_prec mean=0;
    real_prec var=0;
    for(size_t j=0;j< cell_info_tr[i].property.size();++j)
      mean+=cell_info_tr[i].property[j];
    mean/=static_cast<double>(cell_info_tr[i].property.size());
    for(size_t j=0;j<cell_info_tr[i].property.size();++j)
        var+=pow(cell_info_tr[i].property[j]-mean,2);
    var/=static_cast<double>(cell_info_tr[i].property.size());
    if(var>0)
        this->catalogue.set_mach_number(mean/sqrt(var),i);
  }
  if(write)
    this->File.write_array(this->params._Output_directory()+"mach_w_"+hprop[iprop].prop_name, this->mach_number);

 }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**************
// Enable this if relative velocities are to be used in the computation of the mach number per tracer.
// This applies also to the _chunck version of this method.
#define relative_mach
//**************
void HaloTools::get_local_mach_number(real_prec scale)
{
  omp_set_num_threads(_NTHREADS_);
  So.enter(__PRETTY_FUNCTION__);
  real_prec delta=this->params._d_delta_x(); // Fiducial size of cells
  int ncells_back_forth=static_cast<int>(floor(scale/delta))+1;
  // This localize for each tracer the teacer encoded in a distance Rscale from it
  vector<s_nearest_cells>nearest_cells_to_cell(this->params._NGRID());
#ifdef _VERBOSE_CAT_
  So.message_screen("\tNumber of cells back and forth used =", ncells_back_forth);
  So.message_screen("\tIdentifying neighbouring cells: ");
#endif
  get_neighbour_cells_cat_analyze(this->params._Nft(),ncells_back_forth,nearest_cells_to_cell, false);
  So.DONE();
   /*
  // These lines are meant ot write the info of the nearst cells to each cell, but it is heavy and slow
  // to write and read
  // here we write and read the nearst_grid informatio //
  string FNAME="nearest_grid_Nft"+to_string(this->params._Nft())+"_cellsf_bg"+to_string(ncells_back_forth)+".dat";
  ifstream infa(FNAME,ios::in |ios::binary); //we read from this
  if(infa.good()) // if the file exists, read it
  {
    So.message_screen("Reading nearest_grid info");
    for(size_t i=0;i<this->params._NGRID() ;++i)
    {
        for(size_t j=0;j<nearest_cells_to_cell[i].close_cell.size();++j)
        {
           int idnew;
           infa.read((char *)&idnew,sizeof(int));
           nearest_cells_to_cell[i].close_cell.push_back(idnew);
           int a1;
           infa.read((char *)&a1,sizeof(int));
           nearest_cells_to_cell[i].bc_x.push_back(a1);
           int a2;
           infa.read((char *)&a2,sizeof(int));
           nearest_cells_to_cell[i].bc_x.push_back(a2);
           int a3;
           infa.read((char *)&a3,sizeof(int));
           nearest_cells_to_cell[i].bc_x.push_back(a3);
        }
      }
    infa.close();
    So.DONE();
  }

  else // if file does not exist compute it and print it
  {
#ifdef _VERBOSE_CAT_
   So.message_screen("Identifying neighbouring cells");
   So.message_screen("Number of celLs back and forth used =", ncells_back_forth);
#endif
   get_neighbour_cells_cat_analyze(this->params._Nft(),ncells_back_forth,nearest_cells_to_cell, false);
   So.DONE();

   So.message_screen("Writing nearest_grid info");
   ofstream outf(FNAME,ios::out | ios::binary);

   for(size_t i=0;i<this->params._NGRID() ;++i)
   {
      for(size_t j=0;j<nearest_cells_to_cell[i].close_cell.size();++j)
      {
         int idnew=nearest_cells_to_cell[i].close_cell[j];
         outf.write((char *)&idnew,sizeof(int));
         int a1=nearest_cells_to_cell[i].bc_x[j];
         outf.write((char *)&a1,sizeof(int));
         int a2=nearest_cells_to_cell[i].bc_y[j];
         outf.write((char *)&a2,sizeof(int));
         int a3=nearest_cells_to_cell[i].bc_z[j];
         outf.write((char *)&a3,sizeof(int));
      }
    }
    outf.close();
    So.DONE();
  }
*/
  int max_neigh_per_dim = 1+2*ncells_back_forth; // Number of neighbouring cells pr dimension, including the same itself;
  int N_Neigh_cells=pow(max_neigh_per_dim,3); // total number of cells to explore around a cell, including the cell
#ifdef _VERBOSE_CAT_
  So.message_screen("\tMaximum separation probed =",  max_neigh_per_dim*sqrt(3.0)*this->params._d_delta_x(),"Mpc /h");
#endif
  //vector<s_cell_info> cell_info_tr(this->params._NGRID());
  vector<s_cell_info_reduced> cell_info_tr(this->params._NGRID());
#ifdef _VERBOSE_CAT_
  So.message_screen("\tIdentifying tracers in neighbouring cells:");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      size_t ID=this->catalogue.GridID_at(i);
//      cell_info_tr[ID].posx_p.push_back(this->catalogue.coord1);// esto se puede evitar escribiendo simplemente el indice y llamando al tracer.Halo en ese índice
//      cell_info_tr[ID].posy_p.push_back(this->catalogue.coord2);
 //     cell_info_tr[ID].posz_p.push_back(this->catalogue.coord3);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  size_t Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(size_t i=0;i<this->params._NGRID() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    Nt_check+=cell_info_tr[i].gal_index.size();
// expected number of tracer in a sphere of radius Scale
  real_prec Nbar=(4./3.)*M_PI*pow(scale/this->params._Lbox(), 3)*this->catalogue._NOBJS();
#ifdef _VERBOSE_CAT_
  So.message_screen("\tChecking. Number of tracers in structure: ", Nt_check);
  So.message_screen("\tNtracers Expected in structure: ", NOBJS);
  So.message_screen("\tMean number of tracers in a sphere of radius R :", Nbar);
  So.message_screen("\tcomputing...");
#endif
  scale=scale*scale; // use the square to avoid takin sqrt() inside loops
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->catalogue.coord1_at(i);
      real_prec y_coord=this->catalogue.coord2_at(i);
      real_prec z_coord=this->catalogue.coord3_at(i);
      size_t ID=this->catalogue.GridID_at(i);
      real_prec meanv=0; // variable used to compute mach number per each tracer
      size_t counter_gal=0;// number of tracers in the sphere of radius scale around the current tracer
      real_prec meand=0; // variable used to compute mach number per each tracer
      // Get the mean vel of tracers in a sphere of radii scale around the current tracer
      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
        {
          size_t ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->params._Lbox();
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->params._Lbox();
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->params._Lbox();
          for(int k=0; k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
            {
              size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
              real_prec distx = x_coord - (this->catalogue.coord1_at(index_gal)+factor_bc_x);
              real_prec disty = y_coord - (this->catalogue.coord2_at(index_gal)+factor_bc_y);
              real_prec distz = z_coord - (this->catalogue.coord3_at(index_gal)+factor_bc_z);
              real_prec dist  = distx*distx + disty*disty+ distz*distz;
              real_prec vx=this->catalogue.vel1_at(index_gal);
              real_prec vy=this->catalogue.vel2_at(index_gal);
              real_prec vz=this->catalogue.vel3_at(index_gal);
#ifdef relative_mach
              vx=this->catalogue.vel1_at(i)-vx;
              vy=this->catalogue.vel2_at(i)-vy;
              vz=this->catalogue.vel3_at(i)-vz;
#endif
              if(dist<=scale)
                {
                  meanv+=_get_modulo(vx,vy,vz);
                  meand+=sqrt(dist); // to compute mean saparation in the spheres
                  counter_gal++;
                }
            }
         }
        if(counter_gal>1)// IF we hae only one tracer (the ith) we will have a non-well defined variance
         {
            meanv/=static_cast<real_prec>(counter_gal);
            meand/=static_cast<real_prec>(counter_gal);
              // Get the variance vel of tracers in a sphere of radii Rscale around the current tracer
            real_prec varv=0;
            real_prec vard=0;
            for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
               {
                  size_t ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
                  real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->params._Lbox();
                  real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->params._Lbox();
                  real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->params._Lbox();
                  for(int k=0;k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
                    {
                      size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                      real_prec distx = x_coord - (this->catalogue.coord1_at(index_gal)+factor_bc_x);
                      real_prec disty = y_coord - (this->catalogue.coord2_at(index_gal)+factor_bc_y);
                      real_prec distz = z_coord - (this->catalogue.coord3_at(index_gal)+factor_bc_z);
                      real_prec dist  = distx*distx + disty*disty+ distz*distz;
                      real_prec vx=this->catalogue.vel1_at(index_gal);
                      real_prec vy=this->catalogue.vel2_at(index_gal);
                      real_prec vz=this->catalogue.vel3_at(index_gal);
#ifdef relative_mach
                      vx=this->catalogue.vel1_at(i)-vx;
                      vy=this->catalogue.vel2_at(i)-vy;
                      vz=this->catalogue.vel3_at(i)-vz;
#endif
                      if(dist<=scale)
                        {
                          varv+=pow(_get_modulo(vx,vy,vz)-meanv,2);
                          vard+=pow(sqrt(dist)-meand,2);
                         }
                     }
              }
            varv/=static_cast<real_prec>(counter_gal);
            vard/=static_cast<real_prec>(counter_gal);
            this->catalogue.set_mach_number(meanv/sqrt(varv), i);
            this->catalogue.set_dach_number(meand/sqrt(vard), i);
            this->catalogue.set_local_overdensity((static_cast<real_prec>(counter_gal)-Nbar)/Nbar , i);//meand/sqrt(vard);
            this->catalogue.set_number_of_neighbours(counter_gal,i);

         }
         else
         {
            this->catalogue.set_mach_number(0, i);
            this->catalogue.set_dach_number(0, i);
            this->catalogue.set_local_overdensity(0 , i);//meand/sqrt(vard);
            this->catalogue.set_number_of_neighbours(0,i);

         }
     }
  nearest_cells_to_cell.clear();nearest_cells_to_cell.shrink_to_fit();
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_local_mach_number_chuncks(real_prec scale)
{
  So.enter(__PRETTY_FUNCTION__);
  omp_set_num_threads(_NTHREADS_);
  real_prec delta=this->params._d_delta_x_low();
  int ncells_back_forth=static_cast<int>(floor(scale/delta))+1;
#ifdef _VERBOSE_CAT_
  So.message_screen("\tNumber of cells back and forth used =", ncells_back_forth);
  So.message_screen("\tIdentifying neighbouring cells: ");
#endif
  int max_neigh_per_dim = 1+2*ncells_back_forth; // Number of neighbouring cells pr dimension, including the current cell;
  int N_Neigh_cells=pow(max_neigh_per_dim,3); // total number of cells to explore around a cell, including the cell
#ifdef _VERBOSE_CAT_
  So.message_screen("\tNft used =",  this->params._Nft_low());
  So.message_screen("\tMaximum separation probed =",  max_neigh_per_dim*sqrt(3.0)*delta,"Mpc /h");
#endif
  real_prec Nbar=(4./3.)*M_PI*pow(scale/this->params._Lbox(), 3)*this->catalogue._NOBJS();
  So.message_screen("\tMean number of tracers in a sphere of radius R :", Nbar);
  int N_slices=static_cast<int>(this->params._Nft_low()/chunck_factor);
#ifdef _VERBOSE_CAT_
  So.message_screen("\tIdentifying tracers in cells:");
#endif
  // This vector structures contains the coords of the tracers in each cell identified with a 3D ID
  vector<s_cell_info_reduced> cell_info_tr(this->params._NGRID_low());
  // This vector structures contains the coords of the tracers in each of the N_slice slices.
  vector<s_cell_info_reduced> cell_info_tr_slice(N_slices);

  for(size_t i=0;i<NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      size_t ID=this->catalogue.GridID_n_at(i);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
      size_t xbin=get_bin(this->catalogue.coord1_at(i),0,N_slices, delta*chunck_factor,false);
      cell_info_tr_slice[xbin].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  size_t Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(size_t i=0;i<this->params._NGRID_low() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    Nt_check+=cell_info_tr[i].gal_index.size();
  So.message_screen("\tChecking: Number of tracers in structure: ", Nt_check);
  Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(size_t i=0;i< N_slices ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    Nt_check+=cell_info_tr_slice[i].gal_index.size();
  So.message_screen("\tChecking: Number of tracers in sliced structure: ", Nt_check);
  So.message_screen("\tcomputing in chuncks...");
  scale=scale*scale; // use the square to avoid takin sqrt() inside loops
  // Loop over chuncks, planes of y-z constant x. X can be one slide of the Nft, or a set composed of factor_chunck slides
  for(int icc=0;icc< N_slices ;++icc)
    {
      So.message_screen_flush("\t\tSlices completed: ",icc*100/static_cast<real_prec>(N_slices));
      //****************************************************************************************************
      // This structures contains the ID (3D) of the neighbouring cells beloinging the the y-z plane
      // abd hence has dimensions Nft² \times the number of planes to add (NC_factor)
      vector<s_nearest_cells>nearest_cells_to_cell(chunck_factor*this->params._Nft_low()*this->params._Nft_low());// 2D + CHUNCKSIZE
      get_neighbour_cells_cat_analyze_chuncks(icc,this->params._Nft_low(),ncells_back_forth,nearest_cells_to_cell, false);
      //****************************************************************************************************
      if(cell_info_tr_slice[icc].gal_index.size()>0)
        {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(size_t i=0;i<cell_info_tr_slice[icc].gal_index.size() ;++i) //loop over the "observed obejcts in the slice"
          {
          size_t index_gal_slice=cell_info_tr_slice[icc].gal_index[i];
          real_prec x_coord=this->catalogue.coord1_at(index_gal_slice);
          real_prec y_coord=this->catalogue.coord2_at(index_gal_slice);
          real_prec z_coord=this->catalogue.coord3_at(index_gal_slice);
          // This is the x-bin (e.g., 0 or 1 for chunck_factor=2) for each super/slice
          int xbinS=get_bin(x_coord, icc*delta*chunck_factor,  chunck_factor, delta,false);
          int ybin=get_bin(y_coord,0,this->params._Nft_low(),delta,false); // normal y-bin
          int zbin=get_bin(z_coord,0,this->params._Nft_low(),delta,false); // normal z-bin
          size_t ID_slice=index_3d(xbinS,ybin,zbin,this->params._Nft_low(),this->params._Nft_low());
          double meanv=0; // variable used to compute mach number per each tracer
          size_t counter_gal=0;// number of tracers in the sphere of radius scale around the current tracer
          double meand=0; // variable used to compute dach number per each tracer
          // Get the mean vel of tracers in a sphere of radii scale around the current tracer
          for(int j=0;j< N_Neigh_cells; ++j) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
            {
              size_t ID_NEIGH = nearest_cells_to_cell[ID_slice].close_cell[j];
              real_prec factor_bc_x=nearest_cells_to_cell[ID_slice].bc_x[j]*this->params._Lbox();
              real_prec factor_bc_y=nearest_cells_to_cell[ID_slice].bc_y[j]*this->params._Lbox();
              real_prec factor_bc_z=nearest_cells_to_cell[ID_slice].bc_z[j]*this->params._Lbox();
              vector<real_prec>mass_n;
              vector<real_prec>dist_n;
              vector<real_prec>spin_n;
              vector<real_prec>concentration_n;
              vector<size_t>id_n;
              for(int k=0; k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
                {
                  size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                  real_prec distx = x_coord - (this->catalogue.coord1_at(index_gal)+factor_bc_x);
                  real_prec disty = y_coord - (this->catalogue.coord2_at(index_gal)+factor_bc_y);
                  real_prec distz = z_coord - (this->catalogue.coord3_at(index_gal)+factor_bc_z);
                  
                  real_prec dist  = distx*distx+ disty*disty + distz*distz;
                  real_prec vx=this->catalogue.vel1_at(index_gal);
                  real_prec vy=this->catalogue.vel2_at(index_gal);
                  real_prec vz=this->catalogue.vel3_at(index_gal);
#ifdef relative_mach
                  vx-=this->catalogue.vel1_at(index_gal_slice);
                  vx-=this->catalogue.vel2_at(index_gal_slice);
                  vx-=this->catalogue.vel3_at(index_gal_slice);
#endif
                  if(dist<=scale)
                  {
                    meanv+=static_cast<double>(sqrt(vx*vx+vy*vy+vz*vz));
                    meand+=sqrt(dist);
                    counter_gal++;
                    mass_n.push_back(this->catalogue.mass_at(index_gal));
                    dist_n.push_back(sqrt(dist));
                    spin_n.push_back(this->catalogue.spin_at(index_gal));
                    id_n.push_back(index_gal);
                    concentration_n.push_back(this->catalogue.concentration_at(index_gal));
                  }
                }
                if(mass_n.size()>0)
                 {
                    size_t idmin_m, idmax_m;
                    size_t idmin_d, idmax_d;
                    // Get the id of the minimum and maximum mass
                    min_max_vector(mass_n, idmin_m,idmax_m);
                    // Get the id of the minimum and maximum separation
                    min_max_vector(dist_n, idmin_d,idmax_d);
                    this->catalogue.set_most_massive_neighbour(this->catalogue.mass_at(id_n[idmax_m]),index_gal_slice);
                    this->catalogue.set_distance_closest_neighbour(dist_n[idmin_d],index_gal_slice);
                    // Get the distance to the most massive neighbour
                    this->catalogue.set_distance_to_most_massive_neighbour(dist_n[idmax_m],index_gal_slice);
                    this->catalogue.set_mass_closest_neighbour(this->catalogue.mass_at(id_n[idmin_d]),index_gal_slice);
                    this->catalogue.set_spin(this->catalogue.spin_at(id_n[idmin_d]),index_gal_slice);
                    this->catalogue.set_concentration_closest_neighbour(this->catalogue.concentration_at(id_n[idmin_d]),index_gal_slice);
                 }
               else
                {
                  this->catalogue.set_most_massive_neighbour(0,index_gal_slice);
                  this->catalogue.set_distance_closest_neighbour(0,index_gal_slice);
                  this->catalogue.set_distance_to_most_massive_neighbour(0,index_gal_slice);
                  this->catalogue.set_mass_closest_neighbour(0,index_gal_slice);
                  this->catalogue.set_spin_closest_neighbour(0,index_gal_slice);
                  this->catalogue.set_concentration_closest_neighbour(0,index_gal_slice);
                }
            }
      // Either we set an if for only mach number or we use lambda/sigma lambda for separations. Open issue
          if(counter_gal>1)// IF we have only one tracer (the ith), we can have a non-well defined variance.
          {
            meanv/=static_cast<double>(counter_gal);
            meand/=static_cast<double>(counter_gal);
            // Get the variance vel of tracers in a sphere of radii Rscale around the current tracer
            double varv=0;
            double vard=0;
            for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
              {
                size_t ID_NEIGH = nearest_cells_to_cell[ID_slice].close_cell[j];
                real_prec factor_bc_x=nearest_cells_to_cell[ID_slice].bc_x[j]*this->params._Lbox();
                real_prec factor_bc_y=nearest_cells_to_cell[ID_slice].bc_y[j]*this->params._Lbox();
                real_prec factor_bc_z=nearest_cells_to_cell[ID_slice].bc_z[j]*this->params._Lbox();
                for(int k=0;k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
              {
                size_t index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                real_prec distx = x_coord - (this->catalogue.coord1_at(index_gal)+factor_bc_x);
                real_prec disty = y_coord - (this->catalogue.coord2_at(index_gal)+factor_bc_y);
                real_prec distz = z_coord - (this->catalogue.coord3_at(index_gal)+factor_bc_z);
                real_prec dist  = distx*distx+ disty*disty + distz*distz;
                real_prec vx=this->catalogue.vel1_at(index_gal);
                real_prec vy=this->catalogue.vel2_at(index_gal);
                real_prec vz=this->catalogue.vel3_at(index_gal);
    #ifdef relative_mach
                vx-=this->catalogue.vel1_at(index_gal_slice);
                vy-=this->catalogue.vel2_at(index_gal_slice);
                vz-=this->catalogue.vel3_at(index_gal_slice);
    #endif
                if(dist<=scale){
                  varv+=static_cast<double>(pow(sqrt(vx*vx+vy*vy+vz*vz)-meanv,2));
                  vard+=static_cast<double>(pow(sqrt(dist)-meand,2));
                }
              }
            }
            varv/=static_cast<real_prec>(counter_gal);
            vard/=static_cast<real_prec>(counter_gal);
            this->catalogue.mach_number_at(index_gal_slice)= varv ==0? 0.: meanv/sqrt(varv);
            this->catalogue.dach_number_at(index_gal_slice)= vard==0? 0.:  meand/sqrt(vard);
            this->catalogue.local_overdensity_at(index_gal_slice)=(static_cast<real_prec>(counter_gal)-Nbar)/Nbar;//meand/sqrt(vard);
            this->catalogue.set_number_of_neighbours(counter_gal, index_gal_slice);
          }
          else
          {
            this->catalogue.mach_number_at(index_gal_slice)=0.;
            this->catalogue.dach_number_at(index_gal_slice)=0.;
            this->catalogue.local_overdensity_at(index_gal_slice)=0.;
            this->catalogue.number_of_neighbours_at(index_gal_slice)=0;
        }
      }
    } // closes if (.gal.size()>0)
      nearest_cells_to_cell.clear();nearest_cells_to_cell.shrink_to_fit();
    }// closes loop over ic
  cout<<endl;
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_tracer_tidal_anisotropy(vector<real_prec>&tidal){
    So.message_screen("Computing tidal anisotropy at halo position");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(size_t i=0;i<this->catalogue._NOBJS() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
//     this->catalogue.tidal_anisotropy = linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->catalogue.coord1,this->catalogue.coord2,this->catalogue.coord3,tidal);
   this->catalogue.tidal_anisotropy_at(i) =tidal[this->catalogue.GridID_at(i)];
     So.DONE();
  this->min_tidal_anisotropy=this->get_min("_TIDAL_ANISOTROPY_");
  this->max_tidal_anisotropy=this->get_max("_TIDAL_ANISOTROPY_");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools:: get_peak_height_at_tracer(){
    So.message_screen("Computing peak height for tracers");
    this->ps.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    real_prec gr=this->Cosmo.growth_factor(this->params._redshift());
    real_prec gr_zero=this->Cosmo.growth_factor(0);
    this->s_cosmo_pars.growth_factor=gr/gr_zero;
    this->s_cosmo_pars.pk_normalization=ps.normalization();
    this->Cosmo.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    this->ps.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    this->stats.set_cosmo_pars(this->s_cosmo_pars); // update stats
    this->stats.compute_int_table_wavenumber(this->s_cosmo_pars.kmin_int,this->s_cosmo_pars.kmax_int,200);
    So.message_screen("        Normalization of power =",this->s_cosmo_pars.pk_normalization);
    So.message_screen("        Growth at current z    =",this->s_cosmo_pars.growth_factor);
     // ******************************************************************************
     int Nbinsm=200;
     vector<gsl_real>masses(Nbinsm,0);
     vector<gsl_real>nus(Nbinsm,0);
     for(size_t i=0;i< masses.size();++i)
     {
         masses[i]=log10(1e11)+(i+0.5)*log10(7e15/1e11)/static_cast<real_prec>(Nbinsm);
         nus[i]=this->stats.peak_height(masses[i], this->params._redshift(), &this->s_cosmo_pars);
    }
     gsl_interp_accel *acc = gsl_interp_accel_alloc();
     gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,masses.size());
     gsl_spline_init(spline,&masses[0],&nus[0],masses.size());
#ifdef _USE_OMP_
#pragma omp  parallel for
#endif
     for(size_t i=0;i< this->catalogue._NOBJS();++i)
       this->catalogue.peak_height_at(i)=log10(gsl_spline_eval(spline,log10(this->catalogue.mass_at(i)), acc));
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
      this->min_peak_height=get_min("_PEAK_HEIGHT_");
      this->max_peak_height=get_max("_PEAK_HEIGHT_");
      this->So.message_screen("\tMin NU  =", this->min_peak_height);
      this->So.message_screen("\tMax NU  =", this->max_peak_height);
      So.DONE();
      // This is in case the main property is biined in equal number bins, which is not often the case.
      if(this->params._set_bins_equal_number_tracers_main_property()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
       {
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       this->Number_of_tracers_in_ph_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
       this->get_intervals_equal_number("_PEAK_HEIGHT_",min_aux,max_aux);

       cout<<endl;
       So.message_screen("Mass bins with equal number of tracers identified");
       for(size_t i=0;i<=this->Number_of_tracers_in_ph_bins.size();++i)
       {
          this->params.set_PHbins_max(i,max_aux[i]);
          this->params.set_PHbins_min(i,min_aux[i]);
          So.message_screen("\tmin nu", min_aux[i]);
          So.message_screen("\tmax nu", max_aux[i]);
          So.message_screen("\tNtracers", this->Number_of_tracers_in_ph_bins[i]);
          cout<<endl;
       }
     }
      else // if bins in nu are to be define with constant width:
       {
           this->Number_of_tracers_in_ph_bins.resize(this->params._NPHbins_power(),0);
           So.message_screen("\tnu-bins defined fixed width (parameter file)");
       }


#endif
    this->s_cosmo_pars.cosmological_redshift = this->params._redshift() ;
     So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec HaloTools::pearson_correlation(string name_X, string name_Y){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  pair<real_prec,real_prec>stats_X=this->get_variance(name_X, true);
  real_prec Xmean=stats_X.first;
  real_prec Xvar=stats_X.second;
  pair<real_prec,real_prec>stats_Y=this->get_variance(name_Y, true);
  real_prec Ymean=stats_Y.first;
  real_prec Yvar=stats_Y.second;
  real_prec corr=0.;
// estoy hay que generalizarlo para otras propiedades primarias, ojo,
  if(Yvar>0 && Xvar>0)
    {
      if(name_X==_MASS_)
        {
          if(name_Y==_VMAX_)
          {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
       for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.vmax_at(i))-Ymean);
      }
      else if(name_Y=="_RS_")
        {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.rs_at(i))-Ymean);
        }
        else if(name_Y=="_RVIR_")
            {
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.rvir_at(i))-Ymean);
            }
          else if(name_Y=="_CONCENTRATION_")
              {
      #ifdef _USE_OMP_
      #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
      #endif
                for(size_t i=0;i<NOBJS;++i)
                 corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.concentration_at(i))-Ymean);
              }
      else if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.spin_at(i))-Ymean);
          }
          else if(name_Y=="_SPIN_BULLOCK_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
        }
      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
            {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
          else if(name_Y=="_LOCALDM_")
                {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
            }

      else if(name_Y=="_PEAK_HEIGHT_")
      {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);

      }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(log10(this->catalogue.mass_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }

    }
      // ----------------------------------------
      else if (name_X=="_VMAX_")
    {
      if(name_Y=="_RS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.rs_at(i))-Ymean);
        }
      else if(name_Y=="_RVIR_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.rvir_at(i))-Ymean);
        }
      else if(name_Y=="_CONCENTRATION_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.concentration_at(i))-Ymean);
        }
      else if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.spin_at(i))-Ymean);
        }
      else if(name_Y=="_SPIN_BULLOCK_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
        }
      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }

      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vmax_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
   else if (name_X=="_RS_")
    {
              if(name_Y=="_RVIR_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(log10(this->catalogue.rvir_at(i))-Ymean);
                }
             else if(name_Y=="_CONCENTRATION_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(log10(this->catalogue.concentration_at(i))-Ymean);
                }

      if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(log10(this->catalogue.spin_at(i))-Ymean);
        }
      if(name_Y=="_SPIN_BULLOCK_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
        }

      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=((this->catalogue.rs_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.rs_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
       // ----------------------------------------
           else if (name_X=="_RVIR_")
         {
                   if(name_Y=="_CONCENTRATION_")
                     {
             #ifdef  _USE_SIMD_OMP_
             #pragma omp simd reduction(+:corr)
             #else
             #ifdef _USE_OMP_

             #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
             #endif
             #endif
                       for(size_t i=0;i<NOBJS;++i)
                       corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(log10(this->catalogue.concentration_at(i))-Ymean);
                     }
           if(name_Y=="_SPIN_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(log10(this->catalogue.spin_at(i))-Ymean);
             }
           else if(name_Y=="_SPIN_BULLOCK_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
             }

           else if(name_Y=="_VRMS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
             }
           else if(name_Y=="_VIRIAL_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
             }
           else if(name_Y=="_BTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
             }
           else if(name_Y=="_CTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
             for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
             }
           else if(name_Y=="_MACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
             }
           else if(name_Y=="_BIAS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
             }
           else if(name_Y=="_LOCAL_OVERDENSITY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
             }
           else if(name_Y=="_TIDAL_ANISOTROPY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
             }
           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
             }
           else if(name_Y=="_PEAK_HEIGHT_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
                   corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);

             }
           else if(name_Y=="_DACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.rvir_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
             }
         }
     else if (name_X=="_CONCENTRATION_")
         {
           if(name_Y=="_SPIN_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(log10(this->catalogue.spin_at(i))-Ymean);
             }
           else if(name_Y=="_SPIN_BULLOCK_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_
     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
             }

           else if(name_Y=="_VRMS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
             }
           else if(name_Y=="_VIRIAL_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
             }
           else if(name_Y=="_BTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
             }
           else if(name_Y=="_CTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
             for(size_t i=0;i<NOBJS;++i)
               corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
             }
           else if(name_Y=="_MACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
             }
           else if(name_Y=="_BIAS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
             }
           else if(name_Y=="_LOCAL_OVERDENSITY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
             }
           else if(name_Y=="_TIDAL_ANISOTROPY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
             }

           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
             }


           else if(name_Y=="_PEAK_HEIGHT_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
                   corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);

             }
           else if(name_Y=="_DACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.concentration_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
             }

         }
          // ----------------------------------------
      else if (name_X=="_SPIN_")
    {

           if(name_Y=="_SPIN_BULLOCK_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(log10(this->catalogue.spin_bullock_at(i))-Ymean);
                }

      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
          }
           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
             }

      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
    }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      else if (name_X=="_SPIN_BULLOCK_")
    {
      if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(log10(this->catalogue.vrms_at(i))-Ymean);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
          }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
    }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.spin_bullock_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      else if (name_X=="_VRMS_")
    {
      if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.virial_at(i)-Ymean);
        }
        else if(name_Y=="_BTOA_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
          }
        else if(name_Y=="_CTOA_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
          }
        else if(name_Y=="_MACH_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
          }
        else if(name_Y=="_BIAS_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.bias_at(i)-Ymean);
          }
        else if(name_Y=="_LOCAL_OVERDENSITY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
          }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
          }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
          }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(log10(this->catalogue.vrms_at(i))-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
  }
      // ----------------------------------------
      else if (name_X=="_VIRIAL_")
    {
      if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.b_to_a_at(i)-Ymean);
        }

      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }

      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.virial_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      // ----------------------------------------
      else if (name_X=="_BTOA_")
    {
      if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.c_to_a_at(i)-Ymean);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }

      else if(name_Y=="_PEAK_HEIGHT_")

        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.b_to_a_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      // ----------------------------------------
      else if (name_X=="_CTOA_")
    {
      if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
      for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.mach_number_at(i)-Ymean);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
      for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
      for(size_t i=0;i<NOBJS;++i)
       corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
      }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.c_to_a_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
// ----------------------------------------
   else if (name_X=="_MACH_")
    {
      if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.bias_at(i)-Ymean);
        }
          else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
          else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.mach_number_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }

      }
      // ----------------------------------------
      else if (name_X=="_BIAS_")
    {
      if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.bias_at(i)-Xmean)*(this->catalogue.local_overdensity_at(i)-Ymean);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.bias_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.bias_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }


      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.bias_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.bias_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }

    }
      // ----------------------------------------
      else if (name_X=="_LOCAL_OVERDENSITY_")
    {
      if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif

    for(size_t i=0;i<NOBJS;++i)
          corr+=(this->catalogue.local_overdensity_at(i)-Xmean)*(this->catalogue.tidal_anisotropy_at(i)-Ymean);
      }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.local_overdensity_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
      {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.local_overdensity_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
      }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.local_overdensity_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }

    }
      // ----------------------------------------
      else if (name_X=="_TIDAL_ANISOTROPY_")
    {
          if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.tidal_anisotropy_at(i)-Xmean)*(this->catalogue.peak_height_at(i)-Ymean);
        }
       else if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.tidal_anisotropy_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
            }


      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=(this->catalogue.tidal_anisotropy_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      // ----------------------------------------
      else if (name_X=="_PEAK_HEIGHT_")
    {
       if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.peak_height_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
            }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=(this->catalogue.peak_height_at(i)-Xmean)*(this->catalogue.dach_number_at(i)-Ymean);
        }
    }
      else if (name_X=="_DACH_")
    {
       if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=(this->catalogue.dach_number_at(i)-Xmean)*(this->catalogue.local_dm_at(i)-Ymean);
            }
    }
     corr/=(sqrt(Xvar*Yvar)*static_cast<real_prec>(NOBJS));
    }
  return corr;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::PCA(vector<string>&name_props, vector<bool>&used_prop, string extra_info, bool write_cat){
  // steps:
  //    i) take var and mean, and standarize each variable for each halo.
  int dimen=name_props.size();
  vector<real_prec>mean(dimen, 0);
  vector<real_prec>sigma(dimen, 0);
  vector<string>names_used;
  real_prec sum_sigma=0;
  for(size_t i=0;i<dimen;++i)
    {
      if(true==used_prop[i])
        {
          pair<real_prec,real_prec>stats_X=this->get_variance(name_props[i], true);
          mean[i]=stats_X.first;
          sigma[i]=sqrt(stats_X.second);
          sum_sigma+=pow(sigma[i],2);// add the diagonal of the variances
      }
   }
  int Np_used=0;
  for(size_t i=0; i<dimen;++i)
    if(true==used_prop[i])
      Np_used++;
  vector<real_prec> standarized_data(this->catalogue._NOBJS()*Np_used, 0);
  int ip_label=0;// index over *ALL * properties
  int ip_label_used=0;// index over *USED * properties, again
  // Las matriz de datos estandarizados la llenamos de la forma C[i][j] donde i =0..Nobs y j=0, Np_used
  // Ojp que toca escribir esto a mano siguiendo el orden que viene de PowerSpectrum
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.mass_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("\\$\\log_{10} M_{vir}$\\");
  }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.vmax_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10} V_{max}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.rs_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$R_s$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.rvir_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10}Rvir$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.concentration_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10}c$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.spin_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10}\\lambda$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.spin_bullock_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10}\\lambda_b$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
        standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=(log10(this->catalogue.vrms_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\log_{10} \\sigma_{v}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.virial_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{V}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.b_to_a_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{T}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.c_to_a_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{E}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.mach_number_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{M}");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.bias_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$b_{h}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.local_overdensity_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\delta_{R}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.tidal_anisotropy_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{T}_{A}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.peak_height_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\nu$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.dach_number_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\\mathcal{D}$");
    }
  ip_label++;

  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0; i<this->catalogue._NOBJS();++i)
    standarized_data[index_2d(ip_label_used,i,this->catalogue._NOBJS())]=((this->catalogue.local_dm_at(i))-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$logDM$");
    }
  ip_label++;
  So.message_screen("\tNumber of properties:",ip_label);
  So.message_screen("\tNumber of used properties:",ip_label_used);
  this->So.message_screen("\tComputing covariance matrix of standarized properties (i.e, correlation coeff)");
  vector<real_prec> cova_matrix_b(ip_label_used*ip_label_used, 0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(int ip=0; ip<ip_label_used ;++ip)// loop over parameters
    for(int jp=0; jp<ip_label_used ;++jp)// loop over parameters
      for(size_t i=0; i<this->catalogue._NOBJS();++i) // loop over tracers
          cova_matrix_b[index_2d(ip,jp,ip_label_used)]+=(standarized_data[index_2d(ip,i,this->catalogue._NOBJS())]*standarized_data[index_2d(jp,i,this->catalogue._NOBJS())])/static_cast<real_prec>(this->catalogue._NOBJS());
  this->So.DONE();
  this->So.message_screen("\tComputing eigenvalues and eigenvectors");
  vector<s_eigen_vector>eigen_v (ip_label_used);
  for(size_t i=0; i<ip_label_used ;++i)
    eigen_v[i].eigen_vec.resize(ip_label_used,0);
  get_eigen(cova_matrix_b, eigen_v);
  this->So.DONE();
  // explained variance EV(lambda_j) = sum i=0,j (eigenalues)/ sum(eigenvalues)
  real_prec sum_eigen=0;
  for(size_t i=0; i<ip_label_used ;++i)
     sum_eigen+=abs(eigen_v[i].eigen_val);
  So.message_screen("\tTotal variation:",sum_eigen);
  vector<real_prec>iv(ip_label_used,0);//individual variance variance for each PC
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0; i<ip_label_used ;++i)
      iv[i]=eigen_v[i].eigen_val/sum_eigen;
  vector<real_prec>cv(ip_label_used,0); //cumulative variance for each PC
  for(size_t i=0; i<ip_label_used ;++i)
     for(int j=0; j<i+1 ;++j)
       cv[i]+=abs(eigen_v[j].eigen_val)/sum_eigen;
  string file_pca=this->params._Output_directory()+"PCA.txt"+extra_info;
  this->So.message_screen("\tWriting Eigenvalues for PCA in file", file_pca);
  ofstream pca;pca.open(file_pca.c_str());
  for(size_t j=0; j<ip_label_used ;++j)
      pca<<j+1<<"  "<<eigen_v[j].eigen_val<<"  "<<cv[j]<<"  "<<iv[j]<<"  "<<pow(sigma[j],2)/sum_sigma<<endl; // el signma estará bien acá solo cuando usemos todas las propiedades
  pca.close();
  this->So.DONE();
  vector<real_prec> new_data(this->catalogue._NOBJS()*ip_label_used, 0);
  //    ****************************************
  // Rotate the data to the new basis build by the eigen vectors: NEW_DATA= EIGEN_VEC * STANDARIZED_DATA
  // Ech object is now assigned a set of N PComponents, which are a linear combination of the original properties with weights given by the partial variance
  // (components of the eigenvectors)
  this->So.message_screen("\tComputing new data (PC)");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0; i<this->catalogue._NOBJS();++i)
    for(size_t ip=0; ip<ip_label_used ;++ip)// loop over the number of PC: each PC has its own eigenvalue, eigenvalue
        for(size_t jp=0; jp<ip_label_used ;++jp)//loop over the components of each eigenvector corresponding to each eigenvalue
              new_data[index_2d(ip,i,this->catalogue._NOBJS())]+=eigen_v[ip].eigen_vec[jp]*standarized_data[index_2d(jp,i,this->catalogue._NOBJS())];
  this->So.DONE();
   if(write_cat==true)
   {
    string rcatalog=this->params._Output_directory()+"catalog_reduced.txt";
    string rcatalog_pca=this->params._Output_directory()+"catalog_reduced_pca.txt";
    this->select_random_subsample(0.3, ip_label_used, new_data,rcatalog , rcatalog_pca);
  }
  this->So.message_screen("\tGetting mean and var of new data (PC)");
  vector<real_prec>mean_new_data(ip_label_used, 0);
  vector<real_prec>sigma_new_data(ip_label_used, 0);
  for(size_t ip=0; ip<ip_label_used ;++ip)
   {
     real_prec mean=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean)
#endif
     for(size_t i=0; i<this->catalogue._NOBJS();++i)
        mean+=new_data[index_2d(ip,i,this->catalogue._NOBJS())];
     mean_new_data[ip]=mean/static_cast<real_prec>(this->catalogue._NOBJS());
     mean=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean)
#endif
     for(size_t i=0; i<this->catalogue._NOBJS();++i)
        mean+=pow(new_data[index_2d(ip,i,this->catalogue._NOBJS())]-mean_new_data[ip],2);
     sigma_new_data[ip]=sqrt(mean/static_cast<real_prec>(this->catalogue._NOBJS()));
  }
  this->So.DONE();
  this->So.message_screen("\tComputing correlation of STD data and new data (PC)");
 // Get the correlation between the original parameter and the new ones: pick up a number of PC.
  int N_principal_comps=ip_label_used; // NUmber of chosen PC; so far we use all, but we can set a criteria based on the cumulative variance cv, e.g, chosse the first N containing 80% of the total variance
  vector<real_prec>correlation(ip_label_used*N_principal_comps, 0);
#ifdef _USE_OMP_
#pragma omp parallel for collapse(2)
#endif
  for(int ip=0; ip<ip_label_used ;++ip)// loop over the original parameters
    for(int jp=0; jp< N_principal_comps;++jp)// loop over the PC
      for(size_t i=0; i<this->catalogue._NOBJS();++i)// loop over the tracers
         correlation[index_2d(ip,jp,N_principal_comps)]+=standarized_data[index_2d(ip,i,this->catalogue._NOBJS())]*((new_data[index_2d(jp,i,this->catalogue._NOBJS())]-mean_new_data[jp])/sigma_new_data[jp])/static_cast<real_prec>(this->catalogue._NOBJS());
  this->So.DONE();

  // -Find which halo property correlats the most with PC1 (0)
   gsl_vector *max_correlation=gsl_vector_alloc(ip_label_used);
   for(int ip=0; ip<ip_label_used ;++ip)// loop over the original parameters
        gsl_vector_set(max_correlation, ip, correlation[index_2d(ip,0,N_principal_comps)]);
   gsl_vector *indices=gsl_vector_alloc(ip_label_used);
   for(int ip=0; ip<ip_label_used ;++ip)// loop over the original parameters
        gsl_vector_set(indices, ip,ip);
   gsl_sort_vector2(max_correlation,indices);// this sort in ascending order, so the maximum is the last element

   ofstream maxf;
   file_pca=this->params._Output_directory()+"Main_Halo_Properties.txt"+extra_info;
   maxf.open(file_pca.c_str());
   for(int ip=0; ip<ip_label_used ;++ip)// loop over the original parameters
    maxf<<names_used[gsl_vector_get(indices,ip_label_used-ip-1)]<<"   "<<gsl_vector_get(max_correlation, ip_label_used-ip-1)<<endl;
   maxf.close();
   So.message_screen("Writting ranked list of important halo properties correlatin with PC1 in file", file_pca);
   So.message_screen("Maximum correlation with PC1 obtained with property", names_used[gsl_vector_get(indices,ip_label_used-1)]);
   So.message_screen("Second Maximum correlation with PC1 obtained with property", names_used[gsl_vector_get(indices,ip_label_used-2)]);
   So.message_screen("Third Maximum correlation with PC1 obtained with property", names_used[gsl_vector_get(indices,ip_label_used-3)]);
  int N_eff_pc=0;
  real_prec threshold_var=0.70;
  while(cv[N_eff_pc]<=threshold_var)
    N_eff_pc++;
  So.message_screen("Number of PC components containing 70% of variance:", N_eff_pc);
  file_pca=this->params._Output_directory()+"PCA_correlations.txt"+extra_info;
  So.message_screen("Writting correlation between mass, vmax and rs with main PCs in file:",file_pca);
  ofstream pco;
  pca.open(file_pca.c_str());
  for(int ip=0; ip<ip_label_used ;++ip)
    {
       pca<<GREEN<<names_used[ip]<<RESET<<"\t";
       for(int jp=0; jp<N_eff_pc ;++jp)
           if(correlation[index_2d(ip,jp,N_principal_comps)]>0.8)
             pca<<BOLDBLUE<<correlation[index_2d(ip,jp,N_principal_comps)]<<RESET<<"\t"; // correlacion masa,vmax,rs,spin - PC1
            else
               pca<<correlation[index_2d(ip,jp,N_principal_comps)]<<"\t"; // correlacion masa,vmax,rs,spin - PC1
       pca<<endl;
   }
  pca.close();
  file_pca=this->params._Output_directory()+"PCA_correlations.tex"+extra_info;
  string hline="\\";
  So.message_screen("Writting correlation between mass, vmax and rs with main PCs in file (tex):",file_pca);
  pco.open(file_pca.c_str());
  for(int ip=0; ip<ip_label_used ;++ip)
    {
       pco<<names_used[ip]<<"&"<<"\t";
       for(int jp=0; jp<N_eff_pc ;++jp)
         pco<<"$"<<correlation[index_2d(ip,jp,N_principal_comps)]<<"$ & "; // correlacion masa,vmax,rs,spin - PC1
       pco<<hline<<endl;
   }
  pco.close();
  string file_importance=this->params._Output_directory()+"PCA_importance.txt"+extra_info;
  So.message_screen("Writting importance in file:",file_importance);
  pco.open(file_importance.c_str());
  for(int ip=0; ip<ip_label_used ;++ip)// loop over original parameters
    {
      real_prec importance=0;
      for(int jp=0; jp< N_principal_comps ;++jp) // loop over all PCs
        importance+=abs(correlation[index_2d(ip,jp,N_principal_comps)])*(eigen_v[jp].eigen_val/sum_eigen);
       pco<<importance<<endl;
   }
  pco.close();
  this->So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::Get_SO_tracer(){
    real_prec factor=0.001; //to convert Kpc to Mpc
    real_prec mean_overdensity=0;
    real_prec var_overdensity=0;
#pragma omp parallel for reduction(+:mean_overdensity)
    for(size_t i=0;i<this->catalogue._NOBJS();++i)
         mean_overdensity+=this->Cosmo.SO(this->params._redshift(), this->catalogue.mass_at(i), this->catalogue.rvir_at(i)*factor);
    mean_overdensity/=static_cast<real_prec>(this->catalogue._NOBJS());
#pragma omp parallel for reduction(+:mean_overdensity)
    for(size_t i=0;i<this->catalogue._NOBJS();++i)
         var_overdensity+=pow(this->Cosmo.SO(this->params._redshift(), this->catalogue.mass_at(i), this->catalogue.rvir_at(i)*factor)-mean_overdensity,2);
    var_overdensity/=static_cast<real_prec>(this->catalogue._NOBJS());
    So.message_screen("\tMean spherical overdensty of halos:", mean_overdensity);
    So.message_screen("\tSigma_SO of halos:", sqrt(var_overdensity));
    So.message_screen("\tBryan and Norman SO:", Cosmo.density_contrast_top_hat(this->params._redshift()));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::Get_Ranked_Props(string prop){

  gsl_vector *v=gsl_vector_alloc(this->catalogue._NOBJS());
  gsl_permutation * perm = gsl_permutation_alloc(this->catalogue._NOBJS());
  gsl_permutation * rank = gsl_permutation_alloc(this->catalogue._NOBJS());
 // The Halo_ranked[i].prop stores the rank in a sorted list of the proeprty prop.
  // The Halo_ranked[i].prop=0 corresponds to the lowest value of the quantity prop
  // The Halo_ranked[i].prop=NOBS-1 corresponds to the largest value of the quantity prop
  this->catalogue_ranked.set_params(params);

  if(prop=="_MASS_")
   {
      catalogue_ranked.resize_mass(this->catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i, log10(this->catalogue.mass_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.mass_at(i)=rank->data[i];
       }
    }

    if(prop=="_VMAX_")
    {
      catalogue_ranked.resize_vmax(this->catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.vmax_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.vmax_at(i)=rank->data[i];
       }
    }
    if(prop=="_RS_")
    {
      catalogue_ranked.resize_rs(this->catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.rs_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.rs_at(i)=rank->data[i];
       }
    }
    if(prop=="_RVIR_")
    {
      catalogue_ranked.resize_rvir(this->catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.rvir_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.rvir_at(i)=rank->data[i];
       }
    }
    if(prop=="_CONCENTRATION_")
    {
      
      catalogue_ranked.resize_concentration(this->catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.concentration_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.concentration_at(i)=rank->data[i];
       }
    }
    if(prop=="_SPIN_")
    {
      catalogue_ranked.resize_spin(NOBJS);
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.spin_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.spin_at(i)=rank->data[i];
       }
    }
    if(prop=="_SPIN_BULLOCK_")
    {
      catalogue_ranked.resize_spin_bullock(NOBJS);
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.spin_bullock_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.spin_bullock_at(i)=rank->data[i];
       }
    }
    if(prop=="_VRMS_")
    {
        catalogue_ranked.resize_vrms(NOBJS);
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,log10(this->catalogue.vrms_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_vrms(rank->data[i],i);
       }
    }

    if(prop=="_VIRIAL_")
    {
      catalogue_ranked.resize_virial(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.virial_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_virial(rank->data[i],i);
       }
    }
    if(prop=="_BTOA_")
    {
      catalogue_ranked.resize_b_to_a(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.b_to_a_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_b_to_a(rank->data[i],i);
       }
    }
    if(prop=="_CTOA_")
    {
      catalogue_ranked.resize_c_to_a(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.c_to_a_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_c_to_a(rank->data[i],i);
       }
    }

    if(prop=="_MACH_")
    {
      catalogue_ranked.resize_mach_number(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.mach_number_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_mach_number(rank->data[i],i);
       }
    }
    if(prop=="_BIAS_")
    {
      catalogue_ranked.resize_bias(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.bias_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         this->catalogue_ranked.set_bias(rank->data[i],i);
       }
    }
    if(prop=="_RELATIVE_BIAS_")
    {
       catalogue_ranked.resize_relative_bias(catalogue._NOBJS());
     for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.relative_bias_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_relative_bias(rank->data[i],i);
       }
    }
    if(prop=="_LOCAL_OVERDENSITY_")
    {
      catalogue_ranked.resize_local_overdensity(catalogue._NOBJS());
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.local_overdensity_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_local_overdensity(rank->data[i],i);
       }
    }
    if(prop=="_TIDAL_ANISOTROPY_")
    {
            catalogue_ranked.resize_tidal_anisotropy(NOBJS);

      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.tidal_anisotropy_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_tidal_anisotropy(rank->data[i], i);
       }
    }
    if(prop=="_PEAK_HEIGHT_")
    {
      catalogue_ranked.resize_peak_height(NOBJS);

      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v,i, this->catalogue.peak_height_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_peak_height(rank->data[i], i);
       }
    }
    if(prop=="_DACH_")
    {
            catalogue_ranked.resize_dach_number(NOBJS);

      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.dach_number_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_mach_number(rank->data[i], i);
       }
    }
    if(prop=="_LOCALDM_")
    {
            catalogue_ranked.resize_local_dm(NOBJS);

      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        gsl_vector_set(v, i,this->catalogue.local_dm_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < this->catalogue._NOBJS(); i++)
       {
         double vi = gsl_vector_get(v, i);
         this->catalogue_ranked.set_local_dm(rank->data[i],i);
       }
    }
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::Get_Ranked_Props(Catalogue & catalogue, string prop){

    size_t Nsize=this->catalogue._NOBJS();
    gsl_vector *v=gsl_vector_alloc(this->catalogue._NOBJS());
    gsl_permutation * perm = gsl_permutation_alloc(this->catalogue._NOBJS());
    gsl_permutation * rank = gsl_permutation_alloc(this->catalogue._NOBJS());

    if(prop=="_MASS_")
    {
      this->catalogue_ranked.resize_mass(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i, log10(this->catalogue.mass_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_mass(rank->data[i],i);
    }

    if(prop=="_VMAX_")
    {
      this->catalogue_ranked.resize_vmax(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.vmax_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_vmax(rank->data[i],i);
    }

    if(prop=="_RS_")
    {
      this->catalogue_ranked.resize_rs(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.rs_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_rs(rank->data[i],i);
    }
    if(prop=="_RVIR_")
    {
      this->catalogue_ranked.resize_rvir(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.rvir_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_rvir(rank->data[i],i);
    }
    if(prop=="_CONCENTRATION_")
    {
       this->catalogue_ranked.resize_concentration(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.concentration_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_concentration(rank->data[i],i);
    }

    if(prop=="_SPIN_")
    {
       this->catalogue_ranked.resize_spin(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.spin_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_spin(rank->data[i],i);
    }
    if(prop=="_SPIN_BULLOCK_")
    {
       this->catalogue_ranked.resize_spin_bullock(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.spin_bullock_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_spin_bullock(rank->data[i],i);
    }
    if(prop=="_VRMS_")
    {

      this->catalogue_ranked.resize_vrms(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(this->catalogue.vrms_at(i)));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_vrms(rank->data[i],i);
    }

    if(prop=="_VIRIAL_")
    {
      this->catalogue_ranked.resize_virial(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.virial_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_virial(rank->data[i],i);
    }
    if(prop=="_BTOA_")
    {

      this->catalogue_ranked.resize_b_to_a(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.b_to_a_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_b_to_a(rank->data[i],i);
    }

    if(prop=="_CTOA_")
    {
      this->catalogue_ranked.resize_c_to_a(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.c_to_a_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_c_to_a(rank->data[i],i);
    }

    if(prop=="_MACH_")
    {
      this->catalogue_ranked.resize_mach_number(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.mach_number_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_mach_number(rank->data[i],i);
    }
    if(prop=="_BIAS_")
    {
      this->catalogue_ranked.resize_bias(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.bias_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_bias(rank->data[i],i);
    }
    if(prop=="_RELATIVE_BIAS_")
    {
       this->catalogue_ranked.resize_relative_bias(Nsize);
      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.relative_bias_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_relative_bias(rank->data[i],i);
    }

    if(prop=="_LOCAL_OVERDENSITY_")
    {

      this->catalogue_ranked.resize_local_overdensity(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.local_overdensity_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_local_overdensity(rank->data[i],i);
    }

    if(prop=="_TIDAL_ANISOTROPY_")
    {
      this->catalogue_ranked.resize_tidal_anisotropy(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.tidal_anisotropy_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_tidal_anisotropy(rank->data[i],i);
    }
    if(prop=="_PEAK_HEIGHT_")
    {
             this->catalogue_ranked.resize_peak_height(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v,i, this->catalogue.peak_height_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_peak_height(rank->data[i],i);
    }

    if(prop=="_DACH_")
    {
             this->catalogue_ranked.resize_dach_number(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.dach_number_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_mach_number(rank->data[i],i);
    }
    if(prop=="_LOCALDM_")
    {
             this->catalogue_ranked.resize_local_dm(Nsize);

      for(size_t i=0;i<Nsize;++i)
        gsl_vector_set(v, i,this->catalogue.local_dm_at(i));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (size_t i = 0; i < Nsize; i++)
         this->catalogue_ranked.set_local_dm(rank->data[i],i);
    }
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec HaloTools::spearman_correlation(string name_X, string name_Y){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  real_prec corr=0.;

      if(name_X==_MASS_)
        {
          if(name_Y==_VMAX_)
          {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
       for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i) - this->catalogue_ranked.vmax_at(i),2);
      }
      else if(name_Y=="_RS_")
        {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.rs_at(i),2);
        }
        else if(name_Y=="_RVIR_")
            {
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.rvir_at(i),2);
            }
          else if(name_Y=="_CONCENTRATION_")
              {
      #ifdef _USE_OMP_
      #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
      #endif
                for(size_t i=0;i<NOBJS;++i)
                 corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.concentration_at(i),2);
              }
      else if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.spin_at(i),2);
          }
          else if(name_Y=="_SPIN_BULLOCK_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
        }
      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.vrms_at(i),2);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.virial_at(i),2);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i) -this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
            {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
           corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
          else if(name_Y=="_LOCALDM_")
                {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.local_dm_at(i),2);
            }

      else if(name_Y=="_PEAK_HEIGHT_")
      {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.peak_height_at(i),2);

      }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.mass_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      // ----------------------------------------
      else if (name_X=="_VMAX_")
    {
      if(name_Y=="_RS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.rs_at(i),2);
        }
      else if(name_Y=="_RVIR_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.rvir_at(i),2);
        }
      else if(name_Y=="_CONCENTRATION_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.concentration_at(i),2);
        }


      else if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.spin_at(i),2);
        }
      else if(name_Y=="_SPIN_BULLOCK_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
        }
      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.vrms_at(i),2);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.virial_at(i),2);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }

      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vmax_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }

    }

     // ----------------------------------------
      else if (name_X=="_RS_")
    {
              if(name_Y=="_RVIR_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.rvir_at(i),2);
                }
             else if(name_Y=="_CONCENTRATION_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.concentration_at(i),2);
                }

      if(name_Y=="_SPIN_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.spin_at(i),2);
        }
      if(name_Y=="_SPIN_BULLOCK_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
        }

      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.vrms_at(i),2);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.virial_at(i),2);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.peak_height_at(i),2);

        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.rs_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
     }
     else if (name_X=="_RVIR_")
         {
                   if(name_Y=="_CONCENTRATION_")
                     {
             #ifdef  _USE_SIMD_OMP_
             #pragma omp simd reduction(+:corr)
             #else
             #ifdef _USE_OMP_

             #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
             #endif
             #endif
                       for(size_t i=0;i<NOBJS;++i)
                       corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.concentration_at(i),2);
                     }

           if(name_Y=="_SPIN_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.spin_at(i),2);
             }
           else if(name_Y=="_SPIN_BULLOCK_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
             }

           else if(name_Y=="_VRMS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.vrms_at(i),2);
             }
           else if(name_Y=="_VIRIAL_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.virial_at(i),2);
             }
           else if(name_Y=="_BTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
             }
           else if(name_Y=="_CTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
             for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
             }
           else if(name_Y=="_MACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.mach_number_at(i),2);
             }
           else if(name_Y=="_BIAS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.bias_at(i),2);
             }
           else if(name_Y=="_LOCAL_OVERDENSITY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
             }
           else if(name_Y=="_TIDAL_ANISOTROPY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
             }
           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.local_dm_at(i),2);
             }



           else if(name_Y=="_PEAK_HEIGHT_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
                   corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.peak_height_at(i),2);

             }
           else if(name_Y=="_DACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.rvir_at(i)-this->catalogue_ranked.mach_number_at(i),2);
             }
         }
      else if (name_X=="_CONCENTRATION_")
         {

           if(name_Y=="_SPIN_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.spin_at(i),2);
             }
           else if(name_Y=="_SPIN_BULLOCK_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
             }

           else if(name_Y=="_VRMS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.vrms_at(i),2);
             }
           else if(name_Y=="_VIRIAL_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.virial_at(i),2);
             }
           else if(name_Y=="_BTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
             }
           else if(name_Y=="_CTOA_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
             for(size_t i=0;i<NOBJS;++i)
               corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
             }
           else if(name_Y=="_MACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.mach_number_at(i),2);
             }
           else if(name_Y=="_BIAS_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.bias_at(i),2);
             }
           else if(name_Y=="_LOCAL_OVERDENSITY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
             }
           else if(name_Y=="_TIDAL_ANISOTROPY_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
             }

           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.local_dm_at(i),2);
             }
           else if(name_Y=="_PEAK_HEIGHT_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
                   corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.peak_height_at(i),2);

             }
           else if(name_Y=="_DACH_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.concentration_at(i)-this->catalogue_ranked.mach_number_at(i),2);
             }

         }
      else if (name_X=="_SPIN_")
    {
           if(name_Y=="_SPIN_BULLOCK_")
                {
        #ifdef  _USE_SIMD_OMP_
        #pragma omp simd reduction(+:corr)
        #else
        #ifdef _USE_OMP_

        #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
        #endif
        #endif
                  for(size_t i=0;i<NOBJS;++i)
                  corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.spin_bullock_at(i),2);
                }

      else if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.vrms_at(i),2);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.virial_at(i),2);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
          }
           else if(name_Y=="_LOCALDM_")
             {
     #ifdef  _USE_SIMD_OMP_
     #pragma omp simd reduction(+:corr)
     #else
     #ifdef _USE_OMP_

     #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
     #endif
     #endif
               for(size_t i=0;i<NOBJS;++i)
             corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.local_dm_at(i),2);
             }

      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.peak_height_at(i),2);
    }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      else if (name_X=="_SPIN_BULLOCK_")
    {
      if(name_Y=="_VRMS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.vrms_at(i),2);
        }
      else if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.virial_at(i),2);
        }
      else if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }
      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
          }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }

      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.peak_height_at(i),2);
    }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.spin_bullock_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }

    }
      else if (name_X=="_VRMS_")
    {
      if(name_Y=="_VIRIAL_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.virial_at(i),2);
        }

        else if(name_Y=="_BTOA_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
          }
        else if(name_Y=="_CTOA_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
          }
        else if(name_Y=="_MACH_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.mach_number_at(i),2);
          }

        else if(name_Y=="_BIAS_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.bias_at(i),2);
          }
        else if(name_Y=="_LOCAL_OVERDENSITY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
          }
        else if(name_Y=="_TIDAL_ANISOTROPY_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
          }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
          {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.peak_height_at(i),2);
          }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.vrms_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
  }
      // ----------------------------------------
      else if (name_X=="_VIRIAL_")
    {
      if(name_Y=="_BTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.b_to_a_at(i),2);
        }

      else if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }

      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.virial_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      else if (name_X=="_BTOA_")
    {
      if(name_Y=="_CTOA_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.c_to_a_at(i),2);
        }
      else if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }

      else if(name_Y=="_PEAK_HEIGHT_")

        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.b_to_a_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }

    }
      else if (name_X=="_CTOA_")
    {
      if(name_Y=="_MACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
      else if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.bias_at(i),2);
        }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.c_to_a_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      else if (name_X=="_MACH_")
    {
      if(name_Y=="_BIAS_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
       for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.bias_at(i),2);
      }
      else if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
          else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
     }
      else if (name_X=="_BIAS_")
    {
      if(name_Y=="_LOCAL_OVERDENSITY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.bias_at(i)-this->catalogue_ranked.local_overdensity_at(i),2);
        }
      else if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.bias_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
        }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.bias_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }


      else if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.bias_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.bias_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      // ----------------------------------------
      else if (name_X=="_LOCAL_OVERDENSITY_")
    {
      if(name_Y=="_TIDAL_ANISOTROPY_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
    for(size_t i=0;i<NOBJS;++i)
          corr+=pow(this->catalogue_ranked.local_overdensity_at(i)-this->catalogue_ranked.tidal_anisotropy_at(i),2);
      }
      else if(name_Y=="_LOCALDM_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.local_overdensity_at(i)-this->catalogue_ranked.local_dm_at(i),2);
        }
      else if(name_Y=="_PEAK_HEIGHT_")
      {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
        for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.local_overdensity_at(i)-this->catalogue_ranked.peak_height_at(i),2);
      }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.local_overdensity_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      // ----------------------------------------
      else if (name_X=="_TIDAL_ANISOTROPY_")
    {
          if(name_Y=="_PEAK_HEIGHT_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.tidal_anisotropy_at(i)-this->catalogue_ranked.peak_height_at(i),2);
        }
       else if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.tidal_anisotropy_at(i)-this->catalogue_ranked.local_dm_at(i),2);
            }


      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
        corr+=pow(this->catalogue_ranked.tidal_anisotropy_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      // ----------------------------------------
      else if (name_X=="_PEAK_HEIGHT_")
    {
       if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_

    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.peak_height_at(i)-this->catalogue_ranked.local_dm_at(i),2);
            }
      else if(name_Y=="_DACH_")
        {
#ifdef  _USE_SIMD_OMP_
#pragma omp simd reduction(+:corr)
#else
#ifdef _USE_OMP_

#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
#endif
          for(size_t i=0;i<NOBJS;++i)
              corr+=pow(this->catalogue_ranked.peak_height_at(i)-this->catalogue_ranked.mach_number_at(i),2);
        }
    }
      else if (name_X=="_DACH_")
    {
       if(name_Y=="_LOCALDM_")
            {
    #ifdef  _USE_SIMD_OMP_
    #pragma omp simd reduction(+:corr)
    #else
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
    #endif
              for(size_t i=0;i<NOBJS;++i)
            corr+=pow(this->catalogue_ranked.mach_number_at(i)-this->catalogue_ranked.local_dm_at(i),2);
            }
    }
      real_prec nnn=static_cast<real_prec>(NOBJS);
      corr = 1-(6.0*corr)/(nnn*(nnn*nnn-1.));
      return corr;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function gets the scaling relation P(V|theta)
// This methods is replicated from Bam class, made here simpelrs as we only aim at measuring P(theta|M)
void HaloTools::get_scaling_relation_primary_property(string prop)
{
  this->So.enter(__PRETTY_FUNCTION__);
  size_t Nbins=this->params._NPROPbins_bam();
  size_t LENGHT_AB_ONE = Nbins*Nbins ;
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  real_prec lm_min=this->params._MASSbins_min(0);//este ya viene en log, no hay que sacarle log de nuevo
  So.message_screen("Getting scaling relation");
  real_prec lp_min = log10(this->get_min(prop));
  real_prec lp_max = log10(this->get_max(prop));
  real_prec delta_sec = (lp_max-lp_min)/static_cast<real_prec>(Nbins);
  vector<gsl_real> aux_min(Nbins,0);
  vector<gsl_real> aux_max(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop=="_VMAX_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
         size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
          size_t I_Y= get_bin(log10(this->catalogue.vmax_at(ig)),lp_min,Nbins,delta_sec, true);
         this->ABUNDANCE[index_2d(I_Y,I_X,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop=="_SPIN_" || prop=="_SPIN_BULLOCK_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.spin_bullock_at(ig)),lp_min,Nbins,delta_sec, true);
       this->ABUNDANCE[index_2d(I_Y,I_X,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop=="_CONCENTRATION_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.concentration_at(ig)),lp_min,Nbins,delta_sec  , true);
       this->ABUNDANCE[index_2d(I_Y,I_X,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->CONCENTRATIONBmin;
      aux_max=this->CONCENTRATIONBmax;
  }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
  long aux_a,aux_b;
  So.message_screen("Normalizing");
  for(size_t tr_x = 0; tr_x < Nbins; ++tr_x)
  {
      aux_a=-10000;
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins)]));
          aux_a=aux_b;
        }
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          size_t indexa=index_2d(tr_y, tr_x,Nbins);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
        }
   }
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->So.DONE();
  // Define bins in property. Primary bins are defined before calling this method.
  // Now assign
  So.message_screen("Assingning new property");
  for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
   {
      size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);//ojo con el numero de bies cá, arreglarlo
      real_prec prob=-10.0;
      real_prec ran=10.0;
      size_t i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        prob = this->ABUNDANCE_normalized[index_2d(i_halo_v_bin, I_X,Nbins)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop=="_VMAX_")
       this->catalogue.vmax_at(ig) = pow(10,lp_halo);
      else if(prop=="_CONCENTRATION_")
       this->catalogue.concentration_at(ig) = pow(10,lp_halo);
      else if(prop=="_SPIN_BULLOCK_")
       this->catalogue.spin_bullock_at(ig)= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_scaling_relation_primary_property(string prop_halo, string prop_env)
{
  this->So.enter(__PRETTY_FUNCTION__);
  size_t Nbins=this->params._NPROPbins_bam();
  size_t LENGHT_AB_ONE = Nbins*Nbins*Nbins ;
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  real_prec lm_min=this->params._MASSbins_min(0);//este ya viene en log, no hay que sacarle log de nuevo
  this->define_property_bins(Nbins, prop_halo);
  So.message_screen("Getting scaling relation");
  real_prec lenv_min = this->get_min(prop_env);
  real_prec lenv_max = this->get_max(prop_env);
  real_prec delta_env = (lenv_max-lenv_min)/static_cast<real_prec>(Nbins);
  real_prec lp_min = log10(this->get_min(prop_halo));
  real_prec lp_max = log10(this->get_max(prop_halo));
  real_prec delta_sec = (lp_max-lp_min)/static_cast<real_prec>(Nbins);
  vector<gsl_real> aux_min(Nbins,0);
  vector<gsl_real> aux_max(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop_halo=="_VMAX_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
         size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
         size_t I_Y= get_bin(log10(this->catalogue.vmax_at(ig)),lp_min,Nbins,delta_sec, true);
         size_t I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_VIRIAL_")
             I_Z= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_CTOA_")
             I_Z= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
          this->ABUNDANCE[index_3d(I_Y,I_X,I_Z,Nbins,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop_halo=="_SPIN_" || prop_halo=="_SPIN_BULLOCK_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.spin_bullock_at(ig)),lp_min,Nbins,delta_sec, true);
       size_t I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
       if(prop_env=="_DACH_")
          I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
       if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_VIRIAL_")
          I_Z= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_CTOA_")
          I_Z= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
       this->ABUNDANCE[index_3d(I_Y,I_X,I_Z,Nbins,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop_halo=="_CONCENTRATION_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.concentration_at(ig)),lp_min,Nbins,delta_sec  , true);
       size_t I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_DACH_")
          I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_VIRIAL_")
          I_Z= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_CTOA_")
          I_Z= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
       this->ABUNDANCE[index_3d(I_Y,I_X,I_Z,Nbins,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
     }
     aux_min=this->CONCENTRATIONBmin;
     aux_max=this->CONCENTRATIONBmax;
  }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
  long aux_a,aux_b;
  So.message_screen("Normalizing");
  for(size_t tr_x = 0; tr_x < Nbins*Nbins; ++tr_x)
  {
      aux_a=-10000;
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins)]));
          aux_a=aux_b;
        }
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          size_t indexa=index_2d(tr_y, tr_x,Nbins*Nbins);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
        }
   }
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->So.DONE();
  // Define bins in property. Primary bins are defined before calling this method.
  // Now assign
  So.message_screen("Assingning new property");
  for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
   {
      size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
      size_t I_Z= 0;
      if(prop_env=="_MACH_")
         I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_DACH_")
         I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCALDM_")
         I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_TIDAL_ANISOTROPY_")
         I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCAL_OVERDENSITY_")
         I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_VIRIAL_")
         I_Z= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_CTOA_")
         I_Z= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
      real_prec prob=-10.0;
      real_prec ran=10.0;
      size_t i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        prob = this->ABUNDANCE_normalized[index_3d(i_halo_v_bin, I_X,I_Z, Nbins,Nbins)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop_halo=="_VMAX_")
       this->catalogue.vmax_at(ig) = pow(10,lp_halo);
      else if(prop_halo=="_CONCENTRATION_")
       this->catalogue.concentration_at(ig) = pow(10,lp_halo);
      else if(prop_halo=="_SPIN_BULLOCK_")
       this->catalogue.spin_bullock_at(ig)= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_scaling_relation_primary_property(string prop_halo, string prop_env, string prop_extra)
{
  this->So.enter(__PRETTY_FUNCTION__);
  size_t Nbins=this->params._NPROPbins_bam();
  size_t Ncwt=4;
  size_t LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Ncwt ;
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  real_prec lm_min=this->params._MASSbins_min(0);//este ya viene en log, no hay que sacarle log de nuevo
  this->define_property_bins(Nbins, prop_halo);
  So.message_screen("Getting scaling relation");
  real_prec lenv_min = this->get_min(prop_env);
  real_prec lenv_max = this->get_max(prop_env);
  real_prec delta_env = (lenv_max-lenv_min)/static_cast<real_prec>(Nbins);
  real_prec lp_min = 0;
  real_prec lp_max = 0;
  if(prop_halo=="_BIAS_")
{
     lp_min = this->get_min(prop_halo);
     lp_max = this->get_max(prop_halo);
}
    else{
      lp_min = log10(this->get_min(prop_halo));
      lp_max = log10(this->get_max(prop_halo));
    }
  real_prec delta_sec = (lp_max-lp_min)/static_cast<real_prec>(Nbins);
  vector<gsl_real> aux_min(Nbins,0);
  vector<gsl_real> aux_max(Nbins,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop_halo=="_BIAS_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
         size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
         size_t I_Y= get_bin(this->catalogue.bias_at(ig),lp_min,Nbins,delta_sec, true);
         size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
         size_t I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
          size_t I_W=0;
           if(prop_extra=="_MACH_")
              I_W= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_DACH_")
              I_W= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCALDM_")
              I_W= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_TIDAL_ANISOTROPY_")
              I_W= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCAL_OVERDENSITY_")
              I_W= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_CTOA_")
              I_W= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_VIRIAL_")
              I_W= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
 
           std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
           this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
       }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  if(prop_halo=="_VMAX_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
         size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
         size_t I_Y= get_bin(log10(this->catalogue.vmax_at(ig)),lp_min,Nbins,delta_sec, true);
         size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
         size_t I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
          size_t I_W=0;
           if(prop_extra=="_MACH_")
              I_W= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_DACH_")
              I_W= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCALDM_")
              I_W= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_TIDAL_ANISOTROPY_")
              I_W= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCAL_OVERDENSITY_")
              I_W= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_CTOA_")
              I_W= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_VIRIAL_")
              I_W= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);

           std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
           this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop_halo=="_SPIN_" || prop_halo=="_SPIN_BULLOCK_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.spin_bullock_at(ig)),lp_min,Nbins,delta_sec, true);
       size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
       size_t I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
       if(prop_env=="_DACH_")
          I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
       if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
       size_t I_W=0;
        if(prop_extra=="_MACH_")
           I_W= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_DACH_")
           I_W= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCALDM_")
           I_W= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_env=="_TIDAL_ANISOTROPY_")
           I_W= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCAL_OVERDENSITY_")
           I_W= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_CTOA_")
           I_W= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_VIRIAL_")
           I_W= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);

           std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
           this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary

      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop_halo=="_CONCENTRATION_")
    {
      for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
      {
       size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
       size_t I_Y= get_bin(log10(this->catalogue.concentration_at(ig)),lp_min,Nbins,delta_sec  , true);
       size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
       size_t I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_DACH_")
          I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
       size_t I_W=0;
        if(prop_extra=="_MACH_")
           I_W= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_DACH_")
           I_W= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCALDM_")
           I_W= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_env=="_TIDAL_ANISOTROPY_")
           I_W= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCAL_OVERDENSITY_")
           I_W= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_CTOA_")
           I_W= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_VIRIAL_")
           I_W= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);

        std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
        std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
        this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
     }
     aux_min=this->CONCENTRATIONBmin;
     aux_max=this->CONCENTRATIONBmax;
  }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
  this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
  long aux_a,aux_b;
  So.message_screen("Normalizing");
  for(size_t tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
  {
      aux_a=-10000;
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
          aux_a=aux_b;
        }
      for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          size_t indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
          long AUX=this->ABUNDANCE[indexa];
          this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
        }
   }
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->So.DONE();
  // Define bins in property. Primary bins are defined before calling this method.
  // Now assign
  So.message_screen("Assingning new property");
  for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
   {
      size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lm_min,Nbins,this->logdeltaM, true);
      size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
      size_t I_Z= 0;
      if(prop_env=="_MACH_")
         I_Z= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_DACH_")
         I_Z= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCALDM_")
         I_Z= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_TIDAL_ANISOTROPY_")
         I_Z= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCAL_OVERDENSITY_")
         I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
      size_t I_W=0;
       if(prop_extra=="_MACH_")
          I_W= get_bin(this->catalogue.mach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_DACH_")
          I_W= get_bin(this->catalogue.dach_number_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_LOCALDM_")
          I_W= get_bin(this->catalogue.local_dm_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_TIDAL_ANISOTROPY_")
          I_W= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_LOCAL_OVERDENSITY_")
          I_W= get_bin(this->catalogue.local_overdensity_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_CTOA_")
          I_W= get_bin(this->catalogue.c_to_a_at(ig),lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_VIRIAL_")
          I_W= get_bin(this->catalogue.virial_at(ig),lenv_min,Nbins,delta_env, true);
      real_prec prob=-10.0;
      real_prec ran=10.0;
      size_t i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        std::vector<size_t>idxa={i_halo_v_bin, I_X,I_Z,I_W,I_CWT};
        std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
        prob = this->ABUNDANCE_normalized[index_nd(idxa,dimsa)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop_halo=="_VMAX_")
       this->catalogue.vmax_at(ig) = pow(10,lp_halo);
      if(prop_halo=="_BIAS_")
       this->catalogue.bias_at(ig) = lp_halo;
      else if(prop_halo=="_CONCENTRATION_")
       this->catalogue.concentration_at(ig)= pow(10,lp_halo);
      else if(prop_halo=="_SPIN_BULLOCK_")
       this->catalogue.spin_bullock_at(ig)= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_scaling_relation_bias()
{
  this->So.enter(__PRETTY_FUNCTION__);
  size_t Nbins=this->params._NPROPbins_bam();
  size_t Ncwt=4;
  size_t LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt ;
  this->ABUNDANCE.clear();
  this->ABUNDANCE.shrink_to_fit();
  this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
  So.message_screen("Getting scaling relation");
  real_prec lenv1_min = this->get_min("_LOCALDM_");
  real_prec lenv1_max = this->get_max("_LOCALDM_");
  real_prec delta_env1 = (lenv1_max-lenv1_min)/static_cast<real_prec>(Nbins);
  real_prec lenv2_min = this->get_min("_LOCAL_OVERDENSITY_");
  real_prec lenv2_max = this->get_max("_LOCAL_OVERDENSITY_");
  real_prec delta_env2 = (lenv2_max-lenv2_min)/static_cast<real_prec>(Nbins);
  real_prec lenv3_min = this->get_min("_MACH_");
  real_prec lenv3_max = this->get_max("_MACH_");
  real_prec delta_env3 = (lenv3_max-lenv3_min)/static_cast<real_prec>(Nbins);
  real_prec lenv4_min = this->get_min("_DACH_");
  real_prec lenv4_max = this->get_max("_DACH_");
  real_prec delta_env4 = (lenv4_max-lenv4_min)/static_cast<real_prec>(Nbins);
  real_prec lenv5_min = this->get_min("_TIDAL_ANISOTROPY_");
  real_prec lenv5_max = this->get_max("_TIDAL_ANISOTROPY_");
  real_prec delta_env5 = (lenv5_max-lenv5_min)/static_cast<real_prec>(Nbins);
  // Reconstruct bias from nonlocal properties P(b|dm, Delta5,M5,D5,CWT)
  So.message_screen("Learning bias");
  vector<gsl_real> aux_min(Nbins,0);
  vector<gsl_real> aux_max(Nbins,0);
  real_prec lp_min=get_min("_BIAS_");
  real_prec lp_max=get_max("_BIAS_");
  real_prec delta_sec=(lp_max-lp_min)/static_cast<real_prec>(Nbins);
  for(size_t i=0;i<Nbins;++i)
    {
      aux_min[i]=lp_min+i*delta_sec;
      aux_max[i]=lp_min+(i+1)*delta_sec;
    }
  for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
      size_t I_Y= get_bin(this->catalogue.bias_at(ig),lp_min,Nbins,delta_sec, true);// bin in bias
      size_t I_X= get_bin(this->catalogue.local_dm_at(ig),lenv1_min,Nbins,delta_env1, true);//bin of local dm
      size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property // cwt
      size_t I_Z= get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);//local overdenisty
      size_t I_W= get_bin(this->catalogue.mach_number_at(ig),lenv3_min,Nbins,delta_env3, true);//local overdenisty
      size_t I_D= get_bin(this->catalogue.dach_number_at(ig),lenv4_min,Nbins,delta_env4, true);//local overdenisty
      size_t I_T= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv5_min,Nbins,delta_env5, true);//local overdenisty
      std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_D,I_T,I_CWT};
      std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Nbins,Nbins,Ncwt};    
      this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
    }
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    long aux_a,aux_b;
    So.message_screen("Normalizing");
    for(size_t tr_x = 0; tr_x < Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         size_t indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
    }
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->So.DONE();
    // **********************************************************************************************
    // Assign bias
    So.message_screen("Assingning bias");// FOr the futire, the loop should be done over the mock
    So.message_warning("OJo con el loop en esta linea", __LINE__);
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
        this->catalogue.set_bias_aux(this->catalogue.bias_at(ig), ig); //keep track of the original bias
        size_t I_X= get_bin(this->catalogue.local_dm_at(ig),lenv1_min,Nbins,delta_env1, true);
        size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
        size_t I_Z = get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);
        size_t I_W= get_bin(this->catalogue.mach_number_at(ig),lenv3_min,Nbins,delta_env3, true);//local overdenisty
        size_t I_D= get_bin(this->catalogue.dach_number_at(ig),lenv4_min,Nbins,delta_env4, true);//local overdenisty
        size_t I_T= get_bin(this->catalogue.tidal_anisotropy_at(ig),lenv5_min,Nbins,delta_env5, true);//local overdenisty
        real_prec prob=-10.0;
        real_prec ran=10.0;
        size_t i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           std::vector<size_t>idxa={i_halo_v_bin, I_X,I_Z,I_W,I_D,I_T,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Nbins,Nbins,Ncwt};    
           prob = this->ABUNDANCE_normalized[index_nd(idxa,dimsa)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->catalogue.bias_at(ig) = lp_halo;
        }
    this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    So.message_screen("Learning mass");  // mass froom P(M|deltadm, Delta5, bias,cwt)
    //Now assign Mass based on BIAS and local dark matter and CWT
    LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Ncwt ;
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
    real_prec prop_min = log10(this->get_min(_MASS_));
    real_prec prop_max = log10(this->get_max(_MASS_));
    real_prec delta_prop = (prop_max-prop_min)/static_cast<real_prec>(Nbins);
    for(size_t i=0;i<Nbins;++i)
      {
        aux_min[i]=prop_min+i*delta_prop;
        aux_max[i]=prop_min+(i+1)*delta_prop;
      }
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
      size_t I_Y= get_bin(log10(this->catalogue.mass_at(ig)),prop_min,Nbins,delta_prop, true);
      size_t I_X= get_bin(this->catalogue.local_dm_at(ig),lenv1_min,Nbins,delta_env1, true);
      size_t I_Z = get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);
      size_t I_W = get_bin(this->catalogue.bias_aux_at(ig),lp_min,Nbins,delta_sec, true); // use the original bias to learn from, bias_aux
      size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
      std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
      std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
      this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
    }
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    So.message_screen("Normalizing");
    for(size_t tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         size_t indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
    }
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->So.DONE();
    // Assign new mass:
    So.message_screen("Assingning new mass");
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
        this->catalogue.set_mass_parent(this->catalogue.mass_at(ig), ig); // keep track of original mass
        size_t I_X= get_bin(this->catalogue.local_dm_at(ig),lenv1_min,Nbins,delta_env1, true);
        size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
        size_t I_Z = get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);
        size_t I_W = get_bin(this->catalogue.bias_at(ig),lp_min,Nbins,delta_sec, true);// use the new version of bias
        real_prec prob=-10.0;
        real_prec ran=10.0;
        size_t i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           std::vector<size_t>idxa={i_halo_v_bin, I_X,I_Z,I_W,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
           prob = this->ABUNDANCE_normalized[index_nd(idxa,dimsa)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->catalogue.mass_at(ig) = pow(10,lp_halo);
        }
        this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    ofstream tea; tea.open("test.txt");
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
        tea<<this->catalogue.mass_at(ig)<<"  "<<this->catalogue.mass_parent_at(ig)<<"  "<<this->catalogue.bias_at(ig)<<" "<<this->catalogue.bias_aux_at(ig)<<endl;
    tea.close();
    // **********************************************************************************************
    So.message_screen("Learning Vmax"); //P(Vmax|Masa,Delta5,bias,cwt)
    //Now assign Mass based on BIAS and local dark matter and CWT
    LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Ncwt;
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->ABUNDANCE.resize(LENGHT_AB_ONE, 0);
    prop_min = log10(this->get_min(_VMAX_));
    prop_max = log10(this->get_max(_VMAX_));
    delta_prop = (prop_max-prop_min)/static_cast<real_prec>(Nbins);
    for(size_t i=0;i<Nbins;++i)
      {
        aux_min[i]=prop_min+i*delta_prop;
        aux_max[i]=prop_min+(i+1)*delta_prop;
      }
    lenv1_min = log10(this->get_min(_MASS_));
    lenv1_max = log10(this->get_max(_MASS_));
    delta_env1 = (lenv1_max-lenv1_min)/static_cast<real_prec>(Nbins);
    cout<<prop_min<<"  "<<prop_max<<"  "<<delta_prop<<endl;
    cout<<lenv1_min<<"  "<<lenv1_max<<"  "<<delta_env1<<endl;
    cout<<lenv2_min<<"  "<<lenv2_max<<"  "<<delta_env2<<endl;
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
     size_t I_Y= get_bin(log10(this->catalogue.vmax_at(ig)),prop_min,Nbins,delta_prop, true);
     size_t I_X= get_bin(log10(this->catalogue.mass_parent_at(ig)),lenv1_min,Nbins,delta_env1, true);// use original mass to learn from
     size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
     size_t I_Z = get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);
     size_t I_W = get_bin(this->catalogue.bias_aux_at(ig),lp_min,Nbins,delta_sec, true); // use the original bias to learn from, bias2
     std::vector<size_t>idxa={I_Y,I_X,I_Z,I_W,I_CWT};
     std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
     this->ABUNDANCE[index_nd(idxa,dimsa)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
    }
    this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    So.message_screen("Normalizing");
    for(size_t tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(size_t tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         size_t indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
    }
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->So.DONE();
    // Assign new mass
    So.message_screen("Assingning new vmax");
    for(size_t ig = 0; ig< this->catalogue._NOBJS(); ++ig)
    {
        size_t I_X= get_bin(log10(this->catalogue.mass_at(ig)),lenv1_min,Nbins,delta_env1, true);/// use reconstructed mass
        size_t I_Z = get_bin(this->catalogue.local_overdensity_at(ig),lenv2_min,Nbins,delta_env2, true);
        size_t I_W = get_bin(this->catalogue.bias_at(ig),lp_min,Nbins,delta_sec, true);// use the new version of bias
        size_t I_CWT=static_cast<int>(this->catalogue.gal_cwt_at(ig))-1; // -1 to link it with bins in cwt property
        real_prec prob=-10.0;
        real_prec ran=10.0;
        size_t i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           std::vector<size_t>idxa={i_halo_v_bin, I_X,I_Z,I_W,I_CWT};
           std::vector<size_t>dimsa={0,Nbins,Nbins,Nbins,Ncwt};    
           prob = this->ABUNDANCE_normalized[index_nd(idxa,dimsa)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->catalogue.vmax_at(ig) = pow(10,lp_halo);
        }
    this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_superclusters(string realm, string tbias){
// real can be under, or over, referring to tracers with lowest bias (under) or hights bias
// tbias is the  BIAS or RELATIVE_BIAS
    real_prec mean_bias=0;
   if(tbias=="_BIAS_")
#pragma omp parallel for reduction(+:mean_bias)
    for (size_t i=0;i<this->catalogue._NOBJS(); ++i)
        mean_bias+=this->catalogue.bias_at(i);
   else if (tbias=="_RELATIVE_BIAS_")
#pragma omp parallel for reduction(+:mean_bias)
    for (size_t i=0;i<this->catalogue._NOBJS(); ++i)
        mean_bias+=this->catalogue.relative_bias_at(i);
    mean_bias/=static_cast<real_prec>(this->catalogue._NOBJS());
   So.message_screen("Mean bias = ", mean_bias);
//   this->get_intervals_equal_number_aux("_RELATIVE_BIAS_");
   this->get_intervals_equal_number_aux(tbias);
   real_prec quartile_bias=0;
   if(tbias=="_BIAS_")
       quartile_bias=this->params._BIASbins_min(4); //choose ethe last
   else if (tbias=="_RELATIVE_BIAS_")
       quartile_bias=this->params._RBIASbins_min(4); //choose ethe last
   So.message_screen("Upper quartile at bias = ", quartile_bias);
// ii) Select objects with bias above the mean
    this->catalogue_aux.clear_mem();
    // Rank tracers according to bias bottom-up. I.e, the tracer with lowest bias has this->catalogue_ranked.bias_at(i)=0
    // and the tracer with the highest bias gets this->catalogue_ranked.bias_at(i)=NOBJS-1.
   So.message_screen("Ranking");
   this->Get_Ranked_Props(tbias);
    // We nos select the  Ntracers_highb tracers with the highst or lowest bias
   size_t Ntracers_highb=200000;
   size_t bcounter=0;
   So.message_screen("Getting other numbers");
   if(tbias=="_BIAS_")
   {
     for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
       if(this->catalogue.bias_at(i)>quartile_bias)//selects the Ntracers_highb with the highest bias from the ranked list
       bcounter++;
     So.message_screen("Number of tracers above quartile = ",bcounter);
     bcounter=0;
     for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
     if(this->catalogue.bias_at(i)>mean_bias)//selects the Ntracers_highb with the highest bias from the ranked list
            bcounter++;
     So.message_screen("Number of tracers above mean = ",bcounter);
    }
   else if(tbias=="_RELATIVE_BIAS_")
    {
      for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
       if(this->catalogue.relative_bias_at(i)>quartile_bias)//selects the Ntracers_highb with the highest bias from the ranked list
         bcounter++;
      So.message_screen("Number of tracers above quartile = ",bcounter);
      bcounter=0;
      for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
        if(this->catalogue.relative_bias_at(i)>mean_bias)//selects the Ntracers_highb with the highest bias from the ranked list
         bcounter++;
       So.message_screen("Number of tracers above mean = ",bcounter);
    }
    vector<real_prec>bias_aux(this->catalogue._NOBJS(),0);
    vector<real_prec>bias_aux_ranked(this->catalogue._NOBJS(),0);
    if(tbias=="_BIAS_")
#pragma omp parallel for
       for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
       {
            bias_aux[i]=this->catalogue.bias_at(i);
            bias_aux_ranked[i]=this->catalogue_ranked.bias_at(i);
       }
    else if (tbias=="_RELATIVE_BIAS_")
#pragma omp parallel for
       for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
       {
            bias_aux[i]=this->catalogue.relative_bias_at(i);
            bias_aux_ranked[i]=this->catalogue_ranked.relative_bias_at(i);
       }

    size_t N_overbias=0;
    this->catalogue_aux.clear_mem();
    if(realm=="under")
     {
      for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
//       if(bias_aux_ranked[i]<Ntracers_highb &&  bias_aux[i]<mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
       if(bias_aux_ranked[i]<Ntracers_highb &&  bias_aux[i]<mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
        {
            this->catalogue_aux.push_coord1(this->catalogue.coord1_at(i));
            this->catalogue_aux.push_coord2(this->catalogue.coord2_at(i));
            this->catalogue_aux.push_coord3(this->catalogue.coord3_at(i));
            this->catalogue_aux.push_bias(this->catalogue.bias_at(i));
            this->catalogue_aux.push_spin(this->catalogue.spin_at(i));
            this->catalogue_aux.push_concentration(this->catalogue.concentration_at(i));
            this->catalogue_aux.push_mass(this->catalogue.mass_at(i));
            this->catalogue_aux.push_mach_number(this->catalogue.mach_number_at(i));
            this->catalogue_aux.push_b_to_a(this->catalogue.b_to_a_at(i));
            this->catalogue_aux.push_c_to_a(this->catalogue.c_to_a_at(i));
            this->catalogue_aux.push_relative_bias(this->catalogue.relative_bias_at(i));
            this->catalogue_aux.push_local_dm(this->catalogue.local_dm_at(i));
            this->catalogue_aux.push_tidal_anisotropy(this->catalogue.tidal_anisotropy_at(i));
            N_overbias++;
        }
    }
   else if(realm=="over")
     {
     for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
       if(bias_aux_ranked[i]>this->catalogue._NOBJS()-Ntracers_highb &&  bias_aux[i]>mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
         {
            this->catalogue_aux.push_coord1(this->catalogue.coord1_at(i));
            this->catalogue_aux.push_coord2(this->catalogue.coord2_at(i));
            this->catalogue_aux.push_coord3(this->catalogue.coord3_at(i));
            this->catalogue_aux.push_bias(this->catalogue.bias_at(i));
            this->catalogue_aux.push_spin(this->catalogue.spin_at(i));
            this->catalogue_aux.push_concentration(this->catalogue.concentration_at(i));
            this->catalogue_aux.push_mass(this->catalogue.mass_at(i));
            this->catalogue_aux.push_mach_number(this->catalogue.mach_number_at(i));
            this->catalogue_aux.push_b_to_a(this->catalogue.b_to_a_at(i));
            this->catalogue_aux.push_c_to_a(this->catalogue.c_to_a_at(i));
            this->catalogue_aux.push_relative_bias(this->catalogue.relative_bias_at(i));
            this->catalogue_aux.push_local_dm(this->catalogue.local_dm_at(i));
            this->catalogue_aux.push_tidal_anisotropy(this->catalogue.tidal_anisotropy_at(i));
            N_overbias++;
        }
   }
    // find the rank closer to the mean bias in sphere os bis around the mean
/*
    bool found=false;
    real_prec index_mean=0;
    real_prec counter_i=0;
    while(found==false)// nos asegiramos que hayan obvjetos en el radio de bias. Si no los hay, aumenta el radio.
    {
    real_prec brad=0.001;
    real_prec factor_increase=1;
    index_mean=0;
    counter_i=0;
    for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
        if(this->catalogue.relative_bias< mean_bias+brad*factor_increase && this->catalogue.relative_bias > mean_bias-brad*factor_increase)//selects the Ntracers_highb with the highest bias from the ranked list
        {
         index_mean+=this->catalogue_ranked.relative_bias_at(i);
         counter_i++;
        }

        if(counter_i>0)
            found=true;
        else
            factor_increase+=0.5;
    }
    size_t ind=static_cast<size_t>(index_mean/counter_i); // takje this index as the rank corresponding to the mean bias in the sample
    So.message_screen("Mean rank used:", ind);
   for (size_t i=0; i< this->catalogue._NOBJS(); ++i)
   {
       if(this->catalogue_ranked.relative_bias_at(i)> this->catalogue_ranked[ind].relative_bias-Ntracers_highb/2 &&   this->catalogue_ranked.relative_bias_at(i)<this->catalogue_ranked[ind].relative_bias+Ntracers_highb/2 )//selects the Ntracers_highb with the highest bias from the ranked list
        {
            this->catalogue_aux.push_back(s_Halo());
            this->catalogue_aux[N_overbias].coord1 = this->catalogue.coord1;
            this->catalogue_aux[N_overbias].coord2 = this->catalogue.coord2;
            this->catalogue_aux[N_overbias].coord3 = this->catalogue.coord3;
            this->catalogue_aux[N_overbias].bias = this->catalogue.bias;
            this->catalogue_aux[N_overbias].spin = this->catalogue.spin;
            this->catalogue_aux[N_overbias].concentration = this->catalogue.concentration;
            this->catalogue_aux[N_overbias].mass = this->catalogue.mass;
            this->catalogue_aux[N_overbias].mach_number = this->catalogue.mach_number;
            this->catalogue_aux[N_overbias].b_to_a = this->catalogue.b_to_a;
            this->catalogue_aux[N_overbias].c_to_a = this->catalogue.c_to_a;
            this->catalogue_aux[N_overbias].relative_bias = this->catalogue.relative_bias;
            this->catalogue_aux[N_overbias].local_dm = this->catalogue.local_dm;
            this->catalogue_aux[N_overbias].tidal_anisotropy = this->catalogue.tidal_anisotropy;
            N_overbias++;
        }
    }
    */
   this->So.message_screen("Number of tracers above Nlimit and above mean=", N_overbias);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If this is defined, the nbar is assgned by checking the bin in whcih the tracer z falls in the dndz vector. Only applied so far to CoordinateSystem::EQZ
//#define _use_simple_nbar_assignment_
#ifdef _USE_SEVERAL_RANDOM_FILES_
void HaloTools::ang_to_cart_coordinates(s_data_structure *s_data, int ir){
#else
void HaloTools::ang_to_cart_coordinates(s_data_structure *s_data){
#endif
    /**
   * Transform the coordinates of the input catalogues to cartessian coordinates
   * returning the same input vector in which the three first columns correspond to
   * to the X,Y,Z coordinates and the mean number density
   */
#ifdef _VERBOSE_POWER_
    So.enter(__PRETTY_FUNCTION__);
#endif
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
#else
  int NTHREADS=1;
#endif

  this->So.message_screen("Transform to cartesian coordinates in catalogue (and assigning nbar) for ", this->catalogue._type_of_object());

  real_prec aXMAX=-1e6,aYMAX=-1e6, aZMAX=-1e6;
  real_prec aXMIN=+1e6,aYMIN=+1e6, aZMIN=+1e6;

#ifdef _USE_SEVERAL_RANDOM_FILES_
 if(ir==0)
   {
     aXMAX=-1e6; aYMAX=-1e6; aZMAX=-1e6;
     aXMIN=+1e6; aYMIN=+1e6; aZMIN=+1e6;
   }
  else if (ir>0)
   {
     aXMAX=this->params._xmax(); aYMAX=this->params._ymax(); aZMAX=this->params._zmax();
     aXMIN=this->params._xmin(); aYMIN=this->params._ymin(); aZMIN=this->params._zmin();
    }
#endif

size_t nlines = this->catalogue._NOBJS();

#ifdef _use_simple_nbar_assignment_
  real_prec delta_z_nbar = 0;
  if(this->params._use_random_catalog())
    delta_z_nbar= s_data->zz_v[1]-s_data->zz_v[0];
#endif


  std::string angles_units = (this->catalogue._type_of_object()=="RANDOM" ? this->params._angles_units_r() : this->params._angles_units_g());
 
  CoordinateSystem sys_coord = (this->catalogue._type_of_object()=="RANDOM" ? this->params._sys_of_coord_r(): this->params._sys_of_coord_g());

  // Factor to transform deg to rad, as the function equ2cart expects angles in rads
  real_prec fac=1.0;
  if(angles_units=="D")
    fac=M_PI/static_cast<double>(180.);

    int n_rc=0;
  int n_zzv=0;
#ifndef _use_simple_nbar_assignment_
  // Preparing for interpolation of the relation nbar(z)
  gsl_spline *spline_nbar;
  gsl_interp_accel *spline_acc_nbar;
  if(this->params._use_random_catalog() && !this->params._nbar_tabulated() && this->params._use_file_nbar())
    {
      spline_nbar = gsl_spline_alloc (gsl_interp_linear,s_data->zz_v.size());
      gsl_spline_init (spline_nbar, &s_data->zz_v[0], &s_data->dndz_v[0], s_data->zz_v.size());// nbar(z)
    }
#endif

// Preparing for interpolation of the relation r(z):
  // if needed initialize gsl structure to interpolate Comoving distance as a function of redshift:
  gsl_spline *spline_zro;
  gsl_interp_accel *spline_acc_zro;
  if((this->params._use_random_catalog() && !this->params._nbar_tabulated()) || this->params._use_file_nbar() || this->params._nbar_tabulated())
    {
      n_rc = s_data->rr_c.size();
      spline_zro = gsl_spline_alloc (gsl_interp_linear, n_rc);
      switch(sys_coord)
      {
       case(CoordinateSystem::EQR): //in this coordinate system, we might need to get z from r and therefrom, nbar from z
         gsl_spline_init (spline_zro, &s_data->rr_c[0], &s_data->zz_c[0], n_rc); //z(r))
       break;
       case(CoordinateSystem::EQZ): //in this coordinates we have z, so we just need r(z) and nbar(z)
         gsl_spline_init (spline_zro, &s_data->zz_c[0], &s_data->rr_c[0], n_rc); //r(z)
       break;  
     }
    }


   switch(sys_coord)
    {
    case(CoordinateSystem::CART):   // Choose Cartesian coordinates:
#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
      {
        gsl_interp_accel *spline_acc_zro;
        spline_acc_zro = gsl_interp_accel_alloc ();
#ifndef _use_simple_nbar_assignment_
        gsl_interp_accel *spline_acc_nbar;
        if(this->params._use_random_catalog() && !this->params._nbar_tabulated())
          spline_acc_nbar = gsl_interp_accel_alloc ();// useless
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
#endif
      for(size_t i=0;i<nlines;++i)
        {
          real_prec nbar=1.0;
          real_prec x=this->catalogue.coord1_at(i);
          real_prec y=this->catalogue.coord2_at(i);
          real_prec z=this->catalogue.coord3_at(i);
          //compute nbar:
          if(this->params._use_random_catalog() && this->params._nbar_tabulated())
              nbar=this->catalogue.mean_density_at(i);
          else
            this->catalogue.set_mean_density(nbar,i);
          if(this->params._use_random_catalog())
            this->catalogue.set_mean_density(nbar,i);
          if(x > aXMAX)
            aXMAX = x;
          if(y > aYMAX)
            aYMAX = y;
          if(z > aZMAX)
            aZMAX = z;
          if(x < aXMIN)
            aXMIN = x;
          if(y < aYMIN)
            aYMIN = y;
          if(z < aZMIN)
            aZMIN = z;
        } // END parallel loop
      } // END parallel region
      break;
    case(CoordinateSystem::EQR):
#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
      {
    spline_acc_zro = gsl_interp_accel_alloc ();
    if(this->params._use_random_catalog() && !this->params._nbar_tabulated() && this->params._use_file_nbar() )
        gsl_interp_accel *spline_acc_nbar;
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
#endif
    for(size_t i=0;i<nlines;++i)
      {
        real_prec x,y,z, zro, nbar;
        real_prec ra_s = this->catalogue.coord1_at(i);
        real_prec dec_s= this->catalogue.coord2_at(i);
        real_prec rr=this->catalogue.coord3_at(i);
        equatorial_to_cartesian(ra_s,dec_s,rr,x, y, z); // Transform to cartesian coord:
        if(this->params._use_random_catalog()) 	  // Compute the mean number density if not tabulated :
          {
            if(this->params._nbar_tabulated())
              nbar=this->catalogue.mean_density_at(i);
           else
             {
                if(this->params._constant_depth() || this->params._use_file_nbar())
                  {
                    zro= gsl_spline_eval(spline_zro, this->catalogue.coord3_at(i), spline_acc_zro);
                    this->catalogue.set_redshift(zro,i);
#ifdef _use_simple_nbar_assignment_
                  size_t zbin=get_bin(zro ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                      nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar); // Compute the mean number density if not tabulated, either fropm file or from nbar measured from the randoms
#endif
                  }
                else
                {
                   zro= gsl_spline_eval(spline_zro, this->catalogue.coord3_at(i), spline_acc_zro); 		      // Interpolate to get redshift
               // Get ane stimate of the men number density:
#ifdef HEALPIX
              //	      my_get_mean_density_interpolated(map, this->params._new_n_dndz,this->params._redshift_min_sample, this->params._redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
#endif
                }
              }
           }
        // If a random catalog is not used, then use the
        // mean number density obtained from the information of L and Ngal
        // passed through the structure s_data
        else
          nbar=mean_density;
        // New coordinates:
        this->catalogue.set_coord1(x,i);
        this->catalogue.set_coord2(y,i);
        this->catalogue.set_coord3(z,i);

        if(!this->params._use_random_catalog()|| !this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->catalogue.set_mean_density(nbar,i);
        // Identify min and max coordinates in order to set size of cube.
        if(x > aXMAX)
          aXMAX = x;
        if(y > aYMAX)
          aYMAX = y;
        if(z > aZMAX)
          aZMAX = z;
        if(x < aXMIN)
          aXMIN = x;
        if(y < aYMIN)
          aYMIN = y;
        if(z < aZMIN)
          aZMIN = z;
       } // END parallel loop
     } // END parallel region
      break;
  case(CoordinateSystem::EQZ):
#pragma omp parallel num_threads(NTHREADS)
    {
      spline_acc_zro = gsl_interp_accel_alloc ();
#ifdef _use_simple_nbar_assignment_
      if(this->params._use_random_catalog() && !this->params._nbar_tabulated())
        spline_acc_nbar = gsl_interp_accel_alloc ();
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
      for(size_t i=0;i<nlines;++i)
       {
        real_prec x=0; real_prec y=0;
        real_prec z=0; real_prec zro=0; 
        real_prec ra_s=this->catalogue.coord1_at(i);
        real_prec dec_s=this->catalogue.coord2_at(i);
        this->catalogue.redshift_at(i)=this->catalogue.coord3_at(i); // just to keep track of the redshift of the tracer
        real_prec redshift_aux  = this->catalogue.redshift_at(i) < s_data->zz_c[0] ? s_data->zz_c[0]: this->catalogue.redshift_at(i);
        real_prec rr=gsl_spline_eval(spline_zro, redshift_aux, spline_acc_zro); // Transform to comoving distance *
        real_prec nbar=mean_density; 
        equatorial_to_cartesian(ra_s, dec_s, rr, x, y, z); 	    // Transform to cartesian coordinates

        if(this->params._use_random_catalog())
          if(!this->params._nbar_tabulated())
              if(this->params._constant_depth() || this->params._use_file_nbar() )
                {
#ifdef _use_simple_nbar_assignment_
                  size_t zbin=get_bin(redshift_aux ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                      nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar); // Compute the mean number density if not tabulated, either fropm file or from nbar measured from the randoms
#endif
                }
#ifdef HEALPIX
              else if(!this->params._constant_depth() && !this->params._use_file_nbar())
                 my_get_mean_density_interpolated(map, this->params._new_n_dndz,this->params._redshift_min_sample, this->params._redshift_max_sample,ra_s, dec_s, redshift_aux, dndz_m, &nbar);
#endif
#ifdef _USE_REDSHIFT_BINS_
        this->catalogue.observed=false;
        if(this->catalogue.coord3<=this->params._redshift_max_sample && this->catalogue.coord3>=this->params._redshift_min_sample)
          this->catalogue.observed=true;
#endif
        // New coordinates:
        this->catalogue.set_coord1(x,i);
        this->catalogue.set_coord2(y,i);
        this->catalogue.set_coord3(z,i);

        if(!this->params._use_random_catalog() || !this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->catalogue.set_mean_density(nbar,i);

          // Identify min and max coordinates in order to set size of cube.
        if(x > aXMAX)
          aXMAX = x;
        if(y > aYMAX)
          aYMAX = y;
        if(z > aZMAX)
          aZMAX = z;
        if(x < aXMIN)
          aXMIN = x;
        if(y < aYMIN)
          aYMIN = y;
        if(z < aZMIN)
          aZMIN = z;
        } // END parallel loop
      } // END parallel region
      break;
    case(CoordinateSystem::EQRZ):
#pragma omp parallel num_threads(NTHREADS)
      {
    spline_acc_zro = gsl_interp_accel_alloc ();
#ifndef _use_simple_nbar_assignment_
    if(this->params._use_random_catalog() && !this->params._nbar_tabulated())
      spline_acc_nbar = gsl_interp_accel_alloc ();
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
    for(size_t i=0;i<nlines;++i)
      {
        real_prec x,y,z,rr, zro, nbar, ra_s, dec_s;
        ra_s  = this->catalogue.coord1_at(i);
        dec_s = this->catalogue.coord2_at(i);
        rr    = this->catalogue.coord3_at(i);
        equatorial_to_cartesian(ra_s, dec_s, rr, x, y, z); 	  // Transform to cartesian coord. *
        if(this->params._use_random_catalog()) 	  // Compute the mean number density if not tabulated                         *
          {
            if(this->params._nbar_tabulated())
              nbar=this->catalogue.mean_density_at(i);
            else
             {
              if(this->params._constant_depth()  || this->params._use_file_nbar())
                {
                  zro=this->catalogue.coord3_at(i);
#ifdef _use_simple_nbar_assignment_
                  size_t zbin=get_bin(zro ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                  nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar);
#endif
                }
#ifdef HEALPIX
              else
                {
                  zro=this->catalogue.coord3_at(i);
            // my_get_mean_density_interpolated(map, this->params._new_n_dndz,this->params._redshift_min_sample, this->params._redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
                }
#endif
            }
         }
        else
          nbar=mean_density;
#ifdef _USE_REDSHIFT_BINS_
        this->catalogue.observed=false;
        if(this->catalogue.coord3<=this->params._redshift_max_sample && this->catalogue.coord3>=this->params._redshift_min_sample)
          this->catalogue.observed=true;
#endif
        // New coordinates:

        this->catalogue.set_coord1(x,i);
        this->catalogue.set_coord2(y,i);
        this->catalogue.set_coord3(z,i);

        if(!this->params._use_random_catalog()|| !this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->catalogue.set_mean_density(nbar,i);
  
          // Identify min and max coordinates in order to set size of cube.
        if(x > aXMAX)
          aXMAX = x;
        if(y > aYMAX)
          aYMAX = y;
        if(z > aZMAX)
          aZMAX = z;
        if(x < aXMIN)
          aXMIN = x;
        if(y < aYMIN)
          aYMIN = y;
        if(z < aZMIN)
          aZMIN = z;

      } // END parallel loop
      } // END parallel region
    } /// end switch
  // Determine dimension of box from the data
  // Get min and maxs:
  real_prec llx=0;
  real_prec lly=0;
  real_prec llz=0;
   if(true==params._new_Lbox())
    {
      llx=fabs(aXMAX-aXMIN);
      lly=fabs(aYMAX-aYMIN);
      llz=fabs(aZMAX-aZMIN);
      llx=max(llx,lly);
      lly=max(lly,llz);
      llz=max(llx,lly);
    }
  // Pass the new offsets to the params class only if not cartessian coordinates. If cartessian coords,
  // xmin remains the one read in parameter file and xmax is xmin+Lbox,m assuming that we are dealing witha cube.
  if(sys_coord>CoordinateSystem::CART)
  {
    this->params.set_Xoffset(0.5*(aXMAX+aXMIN));
    this->params.set_Yoffset(0.5*(aYMAX+aYMIN));
    this->params.set_Zoffset(0.5*(aZMAX+aZMIN));
    this->params.set_xmin(aXMIN);
    this->params.set_ymin(aYMIN);
    this->params.set_zmin(aZMIN);
    this->params.set_xmax(aXMAX);
    this->params.set_ymax(aYMAX);
    this->params.set_zmax(aZMAX);
 }
  real_prec shift_x= this->params._Xoffset() - 0.5*this->params._Lbox();
  real_prec shift_y= this->params._Yoffset() - 0.5*this->params._Lbox();
  real_prec shift_z= this->params._Zoffset() - 0.5*this->params._Lbox();
#ifdef _FULL_VERBOSE_
  if(this->catalogue._type_of_object()=="RANDOM" ){

  cout<<"\t"<<YELLOW<<"Range in x:[" << aXMIN << ":" << aXMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"Range in y:[" << aYMIN << ":" << aYMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"Range in z:[" << aZMIN << ":" << aZMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in x:[" << aXMIN-shift_x << ":" << aXMAX-shift_x << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in y:[" << aYMIN-shift_y << ":" << aYMAX-shift_y << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in z:[" << aZMIN-shift_z << ":" << aZMAX-shift_z << "]" << endl;
  So.message_screen("Xoffset =",this->params._Xoffset());
  So.message_screen("Yoffset =",this->params._Yoffset());
  So.message_screen("Zoffset =",this->params._Zoffset());
  }

  if(this->catalogue._type_of_object()=="RANDOM" )
    {
      if(this->params._new_Lbox())
        So.message_screen("Overwriting value ", this->params._Lbox(), " with ", llz );
      else
        So.message_screen("Keeping input box lengh", this->params._Lbox());
    }
#endif
#ifdef _USE_SEVERAL_RANDOM_FILES_
  if(true==params._new_Lbox() && ir ==0) // The box is defined with the randoms. If nreading seveal, do it from the first file
#else
  if(true==params._new_Lbox()) // The box is defined with the randoms
#endif
  {
  if(this->catalogue._type_of_object()=="RANDOM" ){ // this is not so necesasary, as in POwer we will feed Fft and power with random_this->catalogue.params
      So.message_screen("Computed box from random:", llz);
      this->params.set_Lbox(llz);
      this->params.derived_pars(); // This is only called if Lbox has changed
  }
   else
      So.message_screen("Computed box from data:", llz);
  }
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  _USE_SEVERAL_RANDOM_FILES_
  void HaloTools::get_interpolated_density_field(bool marked, string property, int ir) 
#else
  void HaloTools::get_interpolated_density_field(bool marked, string property) 
#endif
 {
   So.enter(__PRETTY_FUNCTION__);

#ifdef _USE_WEIGHTS_IN_POWER_
   // Identify columns in the corresponding catalog 
   int i_weight1     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
   int i_weight2     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
   int i_weight3     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
   int i_weight4     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
   bool use_weight1  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
   bool use_weight2  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
   bool use_weight3  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
   bool use_weight4  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
#endif
   size_t nlines= this->catalogue._NOBJS();
   vector<real_prec>field(this->params._NGRID(),0);
#ifdef _MASS_WEIGHT_POWER_
   vector<real_prec>field_mw(this->params._NGRID(),0);
#endif
#ifdef _USE_WEIGHTS_IN_POWER_
   vector<real_prec> ow(MAX_NUMBER_WEIGHTS,0);
#endif
   vector<real_prec>mark;



   if(true==marked)
     {
       mark.resize(this->catalogue._NOBJS(), 1);
       select_mark(mark, property);
     }

   // If rsd=true, then there is no need to shift coordinates even if we want redshift space pwoer spectrum, for the catalog might be already "distorted"

   real_prec rsd_x=1.0;
   real_prec rsd_y=1.0;
   real_prec rsd_z=1.0;
   real_prec conversion_factor=1;
   if(!this->params._redshift_space_coords_g() && "redshift_space" == this->params._clustering_space() )
     { 
      switch(this->params._direction_for_rsd_plane_parallel())
          {
            case(LineOfSight::X):rsd_x*=1.0;rsd_y*=0.0;rsd_z*=0.0; break;
            case(LineOfSight::Y):rsd_x*=0.0;rsd_y*=1.0;rsd_z*=0.0; break;
            case(LineOfSight::Z):rsd_x*=0.0;rsd_y*=0.0;rsd_z*=1.0; break;
          }

      if("kmps"==this->params._vel_units_g())
        conversion_factor=(1.+this->params._redshift())/(this->Cosmo.Hubble_function(this->params._redshift()));
      else if("alpt"==this->params._vel_units_g())
        conversion_factor= cgs_Mpc/(this->Cosmo.Hubble_function(this->params._redshift()));
      else if("Mpcph"==this->params._vel_units_g())
        conversion_factor=1;

      conversion_factor*=(VEL_BIAS_POWER);
      So.message_screen("Current redshift =",this->params._redshift());
      So.message_screen("Hubble function at current redshift =",this->Cosmo.Hubble_function(this->params._redshift()));
      So.message_screen("Conversion factor =",conversion_factor);
     
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i< nlines ;++i)
      {
        this->catalogue.coord1_at(i)+=rsd_x*this->catalogue.vel1_at(i)*conversion_factor;
        this->catalogue.coord2_at(i)+=rsd_y*this->catalogue.vel2_at(i)*conversion_factor;
        this->catalogue.coord3_at(i)+=rsd_z*this->catalogue.vel3_at(i)*conversion_factor;
      }
   }


   MassAssignment exp_mas= this->params._mass_assignment();
   So.message_screen("Interpolating", nlines,"tracers of type ", this->catalogue._type_of_object());
   // The selection of the MAS is done outside the loop over the tracers, so we change Ntracers ifs for only 4 ifs.
   double n_selected=0;
   double S_r_power=0;
   double normal_power=0;
   double W_r=0;
   double S_r1=0;
   double S_r2=0;
   double normal_bispectrum=0;
   real_prec vx=0;
   real_prec vy=0;
   real_prec vz=0;
   
#ifdef _USE_OMP_
#ifdef _GET_BISPECTRUM_NUMBERS_
#pragma omp parallel for reduction(+:n_selected,S_r_power,normal_power,W_r,S_r1,S_r2,normal_bispectrum)
#else
#pragma omp parallel for reduction(+:n_selected,S_r_power,normal_power,W_r)
#endif
#endif
   for(size_t i=0;i< nlines ;++i)
     {
       real_prec x=this->catalogue.coord1_at(i);
       real_prec y=this->catalogue.coord2_at(i);
       real_prec z=this->catalogue.coord3_at(i);
#ifdef _USE_MASS_AS_OBSERVABLE_
       real_prec property=this->catalogue.mass;
#elif defined _USE_VMAX_AS_OBSERVABLE_
       real_prec property=this->catalogue.vmax;
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
       if(property>=MINIMUM_PROP_CUT)
     {
#endif
#ifdef _MASS_WEIGHT_POWER_
       real_prec mass=this->catalogue.mass;
#endif
#ifdef _USE_REDSHIFT_BINS_
       if(this->catalogue.observed)
         {
#endif
           double nbar=static_cast<double>(this->mean_density);
           if(this->params._use_random_catalog())// nbar has been already tbulated or assigned, so it's enough asking if we use random cats or not
            nbar=static_cast<double>(this->catalogue.mean_density_at(i));
           double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
           ow[0] = (use_weight1 && (i_weight1<n_columns))? this->catalogue.weight1 : 1.0;
           ow[1] = (use_weight2 && (i_weight2<n_columns))? this->catalogue.weight2 : 1.0;
           ow[2] = (use_weight3 && (i_weight3<n_columns))? this->catalogue.weight3 : 1.0;
           ow[3] = (use_weight4 && (i_weight4<n_columns))? this->catalogue.weight4 : 1.0;
           ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
           double we_fkp=1.0;
           if(this->params._FKP_weight())
             we_fkp=1.0/(1.0+static_cast<double>(this->params._Pest())*nbar);
           ptotal_weight*=we_fkp;
           if(true==marked)
               ptotal_weight*=mark[i];
           // Parameters for the Pk
           n_selected++;                                //number of selected objects
           W_r += ptotal_weight;                          //weighted number of selected objects
           S_r_power += ptotal_weight*ptotal_weight;                   //sum of the squared of weight
           normal_power += nbar*ptotal_weight*ptotal_weight;    //normalization of power spectrum
           // Parameters for the bispectrum: shot noise
#ifdef _GET_BISPECTRUM_NUMBERS_
           real_prec ptotal_weight3 = ptotal_weight2 * ptotal_weight;
           S_r1 += ptotal_weight3*nbar;
           S_r2 += ptotal_weight3;
           normal_bispectrum+= nbar * nbar * ptotal_weight3;
#endif
           // Interpolation of the object density field:
#ifdef _MASS_WEIGHT_POWER_
#ifdef _USE_OMP_
#pragma atomic
#endif
           this->grid_assignment(x, y, z, ptotal_weight, mass, field, field_mw);
#else
           switch(exp_mas){
           case(MassAssignment::NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
             grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_TSC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::PCS):
#ifdef _USE_OMP_
#pragma atomic
#endif
             grid_assignment_PCS(&this->params, x, y, z, ptotal_weight, field);
         break;
           }
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
         }
#endif

#ifdef _USE_REDSHIFT_BINS_
     }
#endif
    }
   // When a random catalog is used, the
   // quantities below, computed from the random,
   // are used to get the normalization and the shot-noise.
   // If a simulation with no random is used,
   // the values below are computed in the function
   // raw sampling.


#ifdef  _USE_SEVERAL_RANDOM_FILES_
   this->n_gal+=n_selected;
   this->w_g+=W_r;
   this->normal_p+=normal_power;
   this->normal_b+=normal_bispectrum;
   this->s_g+=S_r_power;
   this->sg1+=S_r1;
   this->sg2+=S_r2;
   if(marked==true)
     {
       if(ir==0)
           this->field_external_marked.resize(field.size(),0);
       for(size_t i=0;i<this->field_external_marked.size();++i)
           this->field_external_marked[i]+=field[i];
     }
   else
    {
      if(ir==0)
        this->field_external.resize(field.size(),0);
      for(size_t i=0;i<this->field_external.size();++i)
        this->field_external[i]+=field[i];
     }
#else
   this->n_gal = static_cast<size_t>(n_selected);
   this->w_g = static_cast<real_prec>(W_r);
   this->normal_p = static_cast<real_prec>(normal_power);
   this->normal_b = static_cast<real_prec>(normal_bispectrum);
   this->s_g = static_cast<real_prec>(S_r_power);
   this->sg1 = static_cast<real_prec>(S_r1);
   this->sg2 = static_cast<real_prec>(S_r2);
   if(marked==true)
     {
       this->field_external_marked.resize(field.size(),0);
       this->field_external_marked=field;
     }
   else
    {
      this->field_external.resize(field.size(),0);
      this->field_external=field;
    }
#endif
  So.message_screen("\t\t Weighted number of tracers in mesh", get_nobjects(field));

  if(this->params._output_interpolated_field())
  {
    File.write_array(this->params._output_file_interpolated_field(), this->field_external);
      if(this->params._show_interpolated_field()){

        string json_file_plot ="plot_file_interpolated_field.json";
        std::ofstream jfile(json_file_plot);
        json j;
         j["Lbox"]=this->params._Lbox();
         j["show_interpolated_field"]=true;
         j["Nft"] = this->params._Nft() ;
         j["Nslices"] = 20;
         j["Initial_slice"] = 20;
         j["sample"] = params._Name_survey();
         j["name"] = "Density";
         j["redshift"] = params._redshift();
         j["output_file"] = params._output_file_interpolated_field()+".dat";
         jfile<<j.dump(4);
         cout<<"Generating json file "<<json_file_plot<<" for plotting"<<endl;  
         jfile.close();
         system("python3 ../python/cosmolib_plots.py plot_file_interpolated_field.json &");
      }


  }
   field.clear();
   field.shrink_to_fit(); // we liberate memme
#ifdef _MASS_WEIGHT_POWER_
   field_mw.clear();
   field_mw.shrink_to_fit();
#endif
    So.DONE();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void HaloTools::get_interpolated_density_field_real_space(bool marked, string property)
  {
  #ifdef _FULL_VERBOSE_
     So.enter(__PRETTY_FUNCTION__);
  #endif
    // Sampling of the random or real catalog and
    // interpolation of the object density field into a grid.
    // This method expects the coordinates of the input
    // catalogue already transformed into cartessian,
    // with the information of the mean number density written
    // in the corresponding slot as specified in the
    // parameter file or computed from the box.
    // Marks are expected to be normalized with respct to their mean
 #ifdef _USE_WEIGHTS_IN_POWER_
    int i_weight1     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
    int i_weight2     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
    int i_weight3     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
    int i_weight4     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
    bool use_weight1  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
    bool use_weight2  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
    bool use_weight3  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
    bool use_weight4  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
  #endif
    size_t nlines= NOBJS;
    real_prec nbar=static_cast<real_prec>(this->catalogue._NOBJS())/pow(this->params._Lbox(),3);
    size_t n_selected=0;
    double S_r_power=0;
    double normal_power=0;
    double W_r=0;
    double S_r1=0;
    double S_r2=0;
    double normal_bispectrum=0;
    vector<real_prec>field(this->params._NGRID(),0);
  #ifdef _MASS_WEIGHT_POWER_
    vector<real_prec>field_mw(this->params._NGRID(),0);
  #endif
  #ifdef _USE_WEIGHTS_IN_POWER_
    vector<real_prec> ow(MAX_NUMBER_WEIGHTS,0);
  #endif

   vector<real_prec>mark;
   if(true==marked)
    {
      mark.resize(this->catalogue._NOBJS(), 1);
      select_mark(mark, property);
    }


    MassAssignment exp_mas=this->params._mass_assignment();
   // The selection of the MAS is done outside the loop over the tracers, so we change Ntracers ifs for only 4 ifs.
        // Initialize variables used in the estimation of Pk / Bk
    n_selected=0;
    S_r_power=0;
    normal_power=0;
    W_r=0;
    S_r1=0;
    S_r2=0;
    normal_bispectrum=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:n_selected,S_r_power,normal_power,W_r)
#endif
     for (size_t i=0;i< nlines ;++i)
       {
            real_prec x=this->catalogue.coord1_at(i);
            real_prec y=this->catalogue.coord2_at(i);
            real_prec z=this->catalogue.coord3_at(i);
#ifdef _USE_MASS_AS_OBSERVABLE_
            real_prec property=this->catalogue.mass;
  #elif defined _USE_VMAX_AS_OBSERVABLE_
            real_prec property=this->catalogue.vmax;
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
            if(property>=MINIMUM_PROP_CUT)
              {
#endif
#ifdef _MASS_WEIGHT_POWER_
                real_prec mass=this->catalogue.mass;
#endif
#ifdef _USE_REDSHIFT_BINS_
                if(this->catalogue.observed)
                  {
#endif
                    double nbar=static_cast<double>(mean_density);
                    if(this->params._use_random_catalog())
                      nbar=this->catalogue.mean_density_at(i);
                    double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
                    ow[0] = (use_weight1 && (i_weight1<n_columns))? this->catalogue.weight1 : 1.0;
                    ow[1] = (use_weight2 && (i_weight2<n_columns))? this->catalogue.weight2 : 1.0;
                    ow[2] = (use_weight3 && (i_weight3<n_columns))? this->catalogue.weight3 : 1.0;
                    ow[3] = (use_weight4 && (i_weight4<n_columns))? this->catalogue.weight4 : 1.0;
                    ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
                    double we_fkp=1.0;
                    if(this->params._FKP_weight())
                      we_fkp=1.0/(1.0+this->params._Pest()*nbar);
                    ptotal_weight*=we_fkp;
                    if(true==marked)
                        ptotal_weight*=mark[i]/this->mean_property_value;
                    // Parameters for the Pk
                    n_selected++;                                //number of selected objects
                    W_r+=static_cast<double>(ptotal_weight);                          //weighted number of selected objects
                    double ptotal_weight2 = ptotal_weight * ptotal_weight;
                    S_r_power+= static_cast<double>(ptotal_weight2);                   //sum of the squared of weight
                    normal_power += static_cast<double>(nbar*ptotal_weight2);    //normalization of power spectrum
                    // Parameters for the bispectrum: shot noise
#ifdef _GET_BISPECTRUM_NUMBERS_
                    real_prec ptotal_weight3 = ptotal_weight2 * ptotal_weight;
                    S_r1 += ptotal_weight3*nbar;
                    S_r2 += ptotal_weight3;
                    normal_bispectrum+= nbar * nbar * ptotal_weight3;
#endif
#ifdef _MASS_WEIGHT_POWER_
#ifdef _USE_OMP_
#pragma atomic
#endif
           this->grid_assignment(x, y, z, ptotal_weight, mass, field, field_mw);
#else
           switch(exp_mas){
           case(MassAssignment::NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
            grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
         grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_TSC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(MassAssignment::PCS):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_PCS(&this->params, x, y, z, ptotal_weight, field);
         break;
           }
#endif
#ifdef _USE_REDSHIFT_BINS_
              }
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
                 }
#endif
    }
    this->n_gal=n_selected;
    this->w_g=W_r;
    this->normal_p=normal_power;
    this->normal_b=normal_bispectrum;
    this->s_g=S_r_power;
    this->sg1=S_r1;
    this->sg2=S_r2;
    if(marked==true)
    {
        this->field_external_marked.resize(field.size(),0);
        this->field_external_marked=field;
    }
    else{
      this->field_external.resize(field.size(),0);
      this->field_external=field;
    }
    field.clear();
    field.shrink_to_fit();
  #ifdef _MASS_WEIGHT_POWER_
    field_mw.clear();
    field_mw.shrink_to_fit();
  #endif
    So.DONE();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 void HaloTools::get_interpolated_density_field_real_and_redshift_space(bool marked, string property)
  {
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_WEIGHTS_IN_POWER_
    int i_weight1     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
    int i_weight2     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
    int i_weight3     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
    int i_weight4     = (this->catalogue._type_of_object()!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
    bool use_weight1  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
    bool use_weight2  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
    bool use_weight3  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
    bool use_weight4  = (this->catalogue._type_of_object()!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
#endif
    size_t nlines= this->catalogue._NOBJS();
    size_t n_selected=0;
    double S_r_power=0;
    double normal_power=0;
    double W_r=0;
    double S_r1=0;
    double S_r2=0;
    double normal_bispectrum=0;
    vector<real_prec>field(this->params._NGRID(),0); // FOR REAL SPACE
    vector<real_prec>field_s; // For redshift space
    field_s.resize(this->params._NGRID(),0); // For redshift space
#ifdef _MASS_WEIGHT_POWER_
    vector<real_prec>field_mw(this->params._NGRID(),0);
#endif
#ifdef _USE_WEIGHTS_IN_POWER_
    vector<real_prec> ow(MAX_NUMBER_WEIGHTS,0);
#endif

   vector<real_prec>mark;
   if(true==marked)
    {
      mark.resize(this->catalogue._NOBJS(), 1);
      select_mark(mark, property);
    }


    MassAssignment exp_mas=this->params._mass_assignment();// default, NGP

    real_prec rsd_x=1.0;
    real_prec rsd_y=1.0;
    real_prec rsd_z=1.0;

    if(!this->params._redshift_space_coords_g())

       switch(this->params._direction_for_rsd_plane_parallel())
        {
          case(LineOfSight::X):rsd_x*=1.0;rsd_y*=0.0;rsd_z*=0.0; break;
          case(LineOfSight::Y):rsd_x*=0.0;rsd_y*=1.0;rsd_z*=0.0; break;
          case(LineOfSight::Z):rsd_x*=0.0;rsd_y*=0.0;rsd_z*=1.0; break;
        }
    else // of coordinates have already the rsd included (as in a real catalog) and we still use _REDSHIFT_SPACE_ we have to set all rsd to 0
        {
         rsd_x=0.0;
         rsd_y=0.0;
         rsd_z=0.0;
       }
    // Conversion from km/s to Mpc/h for RSD
    real_prec conversion_factor=1;
    if(this->params._redshift_space_coords_g() == false)
    {
    if("kmps"==this->params._vel_units_g())
      conversion_factor=(1.+this->params.s_cosmo_pars.cosmological_redshift)/(this->Cosmo.Hubble_function(this->params.s_cosmo_pars.cosmological_redshift));
    else if("alpt"==this->params._vel_units_g())
      conversion_factor= cgs_Mpc/(this->Cosmo.Hubble_function(this->params.s_cosmo_pars.cosmological_redshift));
    else if("Mpcph"==this->params._vel_units_g())
      conversion_factor=1;
    }
        n_selected=0;
        S_r_power=0;
        normal_power=0;
        W_r=0;
        S_r1=0;
        S_r2=0;
        normal_bispectrum=0;
#ifdef _USE_OMP_
#ifdef _GET_BISPECTRUM_NUMBERS_
#pragma omp parallel for reduction(+:n_selected,S_r_power,normal_power,W_r,S_r1,S_r2,normal_bispectrum)
#else
#pragma omp parallel for reduction(+:n_selected,S_r_power,normal_power,W_r)
#endif
#endif
        for (size_t i=0;i< nlines ;++i)
          {
            real_prec x=this->catalogue.coord1_at(i);
            real_prec y=this->catalogue.coord2_at(i);
            real_prec z=this->catalogue.coord3_at(i);
            real_prec vx=rsd_x*this->catalogue.vel1_at(i)*conversion_factor;
            real_prec vy=rsd_y*this->catalogue.vel2_at(i)*conversion_factor;
            real_prec vz=rsd_z*this->catalogue.vel3_at(i)*conversion_factor;
            real_prec xs=x+vx;
            real_prec ys=y+vy;
            real_prec zs=z+vz;
#ifdef _USE_MASS_AS_OBSERVABLE_
            real_prec property=this->catalogue.mass_at(i);
#elif defined _USE_VMAX_AS_OBSERVABLE_
            real_prec property=this->catalogue.vmax_at(i);
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
            if(property>=MINIMUM_PROP_CUT)
              {
#endif
#ifdef _MASS_WEIGHT_POWER_
                real_prec mass=this->catalogue.mass_at(i);
#endif
#ifdef _USE_REDSHIFT_BINS_
                if(this->catalogue.observed_at(i))
                  {
#endif
                    double nbar=static_cast<double>(mean_density);
                    if(this->params._use_random_catalog())
                      nbar=this->catalogue.mean_density_at(i);
                    double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
                    ow[0] = (use_weight1 && (i_weight1<n_columns))? this->catalogue.weight1 : 1.0;
                    ow[1] = (use_weight2 && (i_weight2<n_columns))? this->catalogue.weight2 : 1.0;
                    ow[2] = (use_weight3 && (i_weight3<n_columns))? this->catalogue.weight3 : 1.0;
                    ow[3] = (use_weight4 && (i_weight4<n_columns))? this->catalogue.weight4 : 1.0;
                    ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
                    double we_fkp=1.0;
                    if(this->params._FKP_weight())
                      we_fkp=1.0/(1.0+this->params._Pest()*nbar);
                    ptotal_weight*=we_fkp;
                    if(true==marked)
                        ptotal_weight*=mark[i]/this->mean_property_value;
                    // Parameters for the Pk
                    n_selected++;                                //number of selected objects
                    W_r+=static_cast<double>(ptotal_weight);                          //weighted number of selected objects
                    double ptotal_weight2 = ptotal_weight * ptotal_weight;
                    S_r_power+= static_cast<double>(ptotal_weight2);                   //sum of the squared of weight
                    normal_power += static_cast<double>(nbar*ptotal_weight2);    //normalization of power spectrum
#ifdef _GET_BISPECTRUM_NUMBERS_
                    real_prec ptotal_weight3 = ptotal_weight2 * ptotal_weight;
                    S_r1 += ptotal_weight3*nbar;
                    S_r2 += ptotal_weight3;
                    normal_bispectrum+= nbar * nbar * ptotal_weight3;
#endif

#ifdef _MASS_WEIGHT_POWER_
                    this->grid_assignment(x, y, z, ptotal_weight, mass, field, field_mw);
#else
                    switch(exp_mas){
                    case(MassAssignment::NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
                     grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
                     grid_assignment_NGP(&this->params, xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(MassAssignment::CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
                        grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
                        grid_assignment_CIC(&this->params, xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(MassAssignment::TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
                        grid_assignment_TSC(&this->params,x, y, z, ptotal_weight, field);
                        grid_assignment_TSC(&this->params,xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(MassAssignment::PCS):
#ifdef _USE_OMP_
#pragma atomic
#endif
                  grid_assignment_PCS(&this->params, x, y, z, ptotal_weight, field);
                  grid_assignment_PCS(&this->params, xs, ys, zs, ptotal_weight, field_s);
                  break;
                    }
#endif
#ifdef _USE_REDSHIFT_BINS_
            }
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
                  }
  #endif
       }
    this->n_gal=n_selected;
    this->w_g=W_r;
    this->normal_p=normal_power;
    this->normal_b=normal_bispectrum;
    this->s_g=S_r_power;
    this->sg1=S_r1;
    this->sg2=S_r2;
    if(marked==true)
    {
        this->field_external_marked.resize(field.size(),0);
        this->field_external_marked=field;
        this->field_external_marked=field;
        this->field_external_marked_s=field_s;
    }
    else
    {
      this->field_external.resize(field.size(),0);
      this->field_external=field;
      this->field_external=field;
      this->field_external_s=field_s;
    }
    field_s.clear();
    field_s.shrink_to_fit();
#ifdef _MASS_WEIGHT_POWER_
    field_mw.clear();
    field_mw.shrink_to_fit();
#endif
    So.DONE();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_mean_number_density(real_prec alpha0, vector<gsl_real> &zz, vector<gsl_real>&rc, vector<gsl_real> &new_zz_v, vector<gsl_real>&new_dndz_v)
  {
   // Get the mean number density.
    So.enter(__PRETTY_FUNCTION__);
    this->So.message_screen("Preparing for interpolation of mean number density");
    size_t nlines= NOBJS;
    size_t nz = zz.size();
    size_t n_dndz = new_zz_v.size();
    real_prec redshift_min_sample=this->params._redshift_min_sample();
    real_prec redshift_max_sample=this->params._redshift_max_sample();
    real_prec area_survey=this->params._area_survey();
    CoordinateSystem sys_of_coord_r=this->params._sys_of_coord_r(); // Always refere to the random catalog, as these are smoothed
    string file_dndz=this->params._file_nbar()+"_calc";
    // WARNINGS
    cout<<RED<<"\t\tComputing dN/dz from the random catalog"<<RESET<<endl;
    cout<<RED<<"\t\t--->Warning 1: "<<CYAN<<" mean number density from redshift distribution computed"<<endl;
    cout<<"\t\tonly if positions of random objects are in spherical coordinates"<<endl;
    cout<<"\t\tIf positions of random objects are in cartesian coordinates,"<<endl;
    cout<<"\t\tmean number density set to unity. "<<RESET<<endl;
    cout<<RED<<"\t\t--->Warning 2: "<<CYAN<<" the area of the surveyed region must be well known"<<endl;
    cout<<"\t\tOtherwise the estimates of number density are biased. "<<endl;
    cout<<"\t\tConsider using a MASK"<<RESET<<endl;
    real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/static_cast<double>(n_dndz);
    size_t n_objects = this->catalogue._NOBJS();
    real_prec area=area_survey*pow(M_PI/180.,2); /*converting to strad*/

    for(size_t i=0;i<n_dndz;++i)
      new_zz_v[i]=redshift_min_sample+(i+0.5)*Delta_Z;

    if(CoordinateSystem::EQR==sys_of_coord_r)
      { // We must transform r to z
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
     for(size_t i=0;i<nlines;i++)
       {
         real_prec zro= gsl_inter_new(rc, zz, this->catalogue.coord3_at(i));
         int count= get_bin(zro, redshift_min_sample,n_dndz,Delta_Z,true); // this is not optimal
  #ifdef _USE_OMP_
  #pragma omp atomic
  #endif
          new_dndz_v[count]++;
        }
      }
     else  if(CoordinateSystem::EQZ==sys_of_coord_r)
      {
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
       for(size_t i=0;i<this->catalogue._NOBJS();i++)
        {
         int count= get_bin(this->catalogue.coord3_at(i), redshift_min_sample,n_dndz,Delta_Z,true);
  #ifdef _USE_OMP_
  #pragma omp atomic update
  #endif
         new_dndz_v[count]++;
      }
    }
    size_t n_rc = rc.size();
    vector<gsl_real>xa (n_rc,0);
    vector<gsl_real>ya (n_rc,0);
    for(int i = 0; i < n_rc; i++)
       {
         xa[i] = static_cast<gsl_real>(zz[i]);
         ya[i] = static_cast<gsl_real>(rc[i]);
     }
    gsl_interp_accel *spline_acc_zro=gsl_interp_accel_alloc ();
    gsl_spline *spline_zro = gsl_spline_alloc (gsl_interp_linear, n_rc);
    gsl_spline_init (spline_zro, &xa[0], &ya[0], n_rc); //r(z)
     // Compute Volume in redshift shell and divide by it to get nbar
    for(size_t i=0;i<n_dndz;i++)
      {
        real_prec rmax=gsl_spline_eval(spline_zro,new_zz_v[i]+0.5*Delta_Z, spline_acc_zro);
        real_prec rmin=gsl_spline_eval(spline_zro,new_zz_v[i]-0.5*Delta_Z, spline_acc_zro);
        real_prec Volume=(area_survey*M_PI/180.0)*(pow(rmax,3)-pow(rmin,3))/3.;
        new_dndz_v[i]*=(alpha0/Volume);  //do not divide by Delta Z.
      }
    /*Smooth the histogram with a cubic spline*/
  //  gsl_bspline(zz_v, dndz_v, new_zz_v, new_dndz_v);
    this->File.write_to_file(file_dndz,new_zz_v,new_dndz_v);
    So.DONE();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_random_catalog(){
    So.enter(__PRETTY_FUNCTION__);
    // FIrst get the distribuytion in z, color and stellar mass:
   try{
        if(this->params._sys_of_coord_g()>CoordinateSystem::CART)
            So.message_screen("");
        else
            throw (this->params._sys_of_coord_g());
    }
    catch (int soc){
        So.message_screen("System of coordinate not valid for this particular task");
        exit(1);
    }
    this->get_cosmo();
    size_t Nbins=this->params._Nbins_color()*this->params._Nbins_Mstellar()*this->params._Nbins_redshift();
    vector<real_prec>Pzcmdist(Nbins,0);// For z-c-M distributuion
    vector<real_prec>Pzdist(this->params._Nbins_redshift(),0);
    vector<real_prec>Pcdist(this->params._Nbins_color(),0);
    vector<real_prec>Pmdist(this->params._Nbins_Mstellar(),0);
    vector<real_prec>Pzdist_smoothed(this->params._Nbins_redshift(),0);
    real_prec delta_redshift=(this->params._redshift_max_sample()-this->params._redshift_min_sample())/static_cast<real_prec>(this->params._Nbins_redshift());
    real_prec delta_color=(this->params._Color_max()-this->params._Color_min())/static_cast<real_prec>(this->params._Nbins_color());
    real_prec delta_Mstellar=(this->params._Mstellar_max()-this->params._Mstellar_min())/static_cast<real_prec>(this->params._Nbins_Mstellar());
    vector<real_prec> zbins(Pzdist.size(), 0);
    vector<real_prec> cbins(Pcdist.size(), 0);
    vector<real_prec> mbins(Pmdist.size(), 0);

    for(size_t i=0;i<Pzdist.size();i++)
        zbins[i]=static_cast<gsl_real>(this->params._redshift_min_sample()+(i+0.5)*delta_redshift);
    
    for(size_t i=0;i<Pcdist.size();i++)
        cbins[i]=this->params._Color_min()+(i+0.5)*delta_color;
    
    for(size_t i=0;i<Pmdist.size();i++)
        mbins[i]=this->params._Mstellar_min()+(i+0.5)*delta_Mstellar;
    //So.message_screen("\tObtaining distributions in galaxy properties with a magnitud cut at m=", this->params._mK_max());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<this->catalogue._NOBJS();++i)
    {
  //      if(this->catalogue.app_mag <this->params._mK_max())
  //      {
        int ibin_z = get_bin(this->catalogue.redshift_at(i),this->params._redshift_min_sample(),this->params._Nbins_redshift(),delta_redshift,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pzdist[ibin_z]++;
        int ibin_color = get_bin(this->catalogue.color_at(i),this->params._Color_min(),this->params._Nbins_color(),delta_color,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pcdist[ibin_color]++;
        int ibin_Mstellar = get_bin(this->catalogue.stellar_mass_at(i),this->params._Mstellar_min(),this->params._Nbins_Mstellar(),delta_Mstellar,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pmdist[ibin_Mstellar]++;
        size_t index3d=index_3d(ibin_z, ibin_color, ibin_Mstellar, this->params._Nbins_color(), this->params._Nbins_Mstellar());
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pzcmdist[index3d]++;
    //}
        }
    So.DONE();
    So.message_screen("\tNormalizing distribution");
   //For each bin in redshift, normalize to maximum the remaining c-M distribution
    vector<real_prec>Pzcmdist_normal(Nbins,0);// For z-c-M distributuion
    for(int iz=0;iz<this->params._Nbins_redshift();iz++)
      {
        vector<real_prec>aux(this->params._Nbins_color()*this->params._Nbins_Mstellar(),0);// For z-c-M distributuion
        for(int ic=0;ic<this->params._Nbins_color();ic++)
          for(int im=0;im<this->params._Nbins_Mstellar();im++)
             aux[index_2d(ic,im, this->params._Nbins_Mstellar())]=Pzcmdist[index_3d(iz, ic, im, this->params._Nbins_color(), this->params._Nbins_Mstellar())];

        real_prec maxd=get_max_nm(aux);
        for(int ic=0;ic<this->params._Nbins_color();ic++)
          for(int im=0;im<this->params._Nbins_Mstellar();im++)
            {
              size_t indexa=index_3d(iz, ic, im, this->params._Nbins_color(), this->params._Nbins_Mstellar());
              Pzcmdist_normal[indexa]= maxd==0? 0: Pzcmdist[indexa]/maxd;
            }
      }
    So.DONE();
    string odn=this->params._Output_directory()+"dNdz.txt";
    this->File.write_to_file(odn, zbins, Pzdist);
    odn=this->params._Output_directory()+"dNdc.txt";
    this->File.write_to_file(odn, cbins, Pcdist);
    odn=this->params._Output_directory()+"dNdMs.txt";
    this->File.write_to_file(odn, mbins, Pmdist);
    // We now smooth the dNdz distribution and compute the mean number densiry, normalize by its maximum.
    this->So.message_screen("Smootring DnDz");
    smooth_distribution(zbins, Pzdist, Pzdist_smoothed);
    odn=this->params._Output_directory()+"dNdz_smooth.txt";
    this->File.write_to_file(odn, zbins, Pzdist_smoothed);
    this->So.DONE();
    vector<real_prec> nbar(Pzdist.size(), 0);
    vector<gsl_real> nbard(Pzdist.size(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i<Pzdist.size();i++)
      {
        real_prec  distance_i=this->Cosmo.comoving_distance(zbins[i]-0.5*delta_redshift); // gsl_inter_new(this->zv,this->rv,zn[i]);
        real_prec  distance_j=this->Cosmo.comoving_distance(zbins[i]+0.5*delta_redshift); // gsl_inter_new(this->zv,this->rv,zn[i]);
        real_prec Vshell=this->params._area_survey()*(pow(distance_j,3)-pow(distance_i,3));
        nbar[i]=Pzdist_smoothed[i]/Vshell; // dz is not included here as dNdz does not included it neither
        nbard[i]=static_cast<gsl_real>(nbar[i]);
    }
    real_prec max_nb=get_max_nm(nbar);
    vector<gsl_real>Pzdistd(this->params._Nbins_redshift(),0);
    vector<gsl_real>zbinsd(Pzdistd.size(), 0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(size_t i=0;i<nbar.size();++i)
     {
       Pzdist[i]=nbar[i]/max_nb;
       Pzdistd[i]=static_cast<gsl_real>(Pzdist[i]);
       zbinsd[i]=static_cast<gsl_real>(zbins[i]);
     }
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRand;
#ifndef _USE_OMP_
    gsl_rng_env_setup();
    gsl_rng_default_seed=224;
    rng_t = gsl_rng_ranlux;
    gBaseRand = gsl_rng_alloc (rng_t);
#endif
    // Now geenrate the random catalog in the following steps:
    //i) In the angular mask, throw a redshift according to the P(z)
    //ii)Throw two numbers, color and stellar mass, and accept then according to P(z,c,M)
    real_prec factor = (this->params._angles_units_g()=="D" ? M_PI/180.0: 1 );
    int Nran_files=this->params._Number_of_random_files();
    this->So.message_screen("Generating random catalog with ", this->catalogue._NOBJS()*Nran_files, " objects in chuncks of ", Nran_files, " files");
    string file=this->params._Output_directory()+"Random_cat";
    gsl_interp_accel *acc_rz;
    acc_rz = gsl_interp_accel_alloc ();
    gsl_spline *spline_rz    = gsl_spline_alloc (gsl_interp_linear, this->rcosmo.size());
    gsl_spline_init (spline_rz, &(this->rcosmo[0]), &(this->zcosmo[0]), this->rcosmo.size());
    gsl_interp_accel *acc_zr;
    acc_zr = gsl_interp_accel_alloc ();
    gsl_spline *spline_zr    = gsl_spline_alloc (gsl_interp_linear, this->rcosmo.size());
    gsl_spline_init (spline_zr, &(this->zcosmo[0]), &(this->rcosmo[0]),this->rcosmo.size());
  // get P(z) given z
    gsl_interp_accel *acc_zp;
    acc_zp = gsl_interp_accel_alloc ();
    gsl_spline *spline_p    = gsl_spline_alloc (gsl_interp_linear, Pzdistd.size());
    gsl_spline_init (spline_p, &(zbinsd[0]), &(Pzdistd[0]), Pzdistd.size());
    real_prec  rmax=static_cast<real_prec>(gsl_spline_eval(spline_zr, this->params._redshift_max_sample(), acc_zr)); // gsl_inter_new(zv,rv,this->z_max);
    real_prec  rmin=static_cast<real_prec>(gsl_spline_eval(spline_zr, this->params._redshift_min_sample(), acc_zr)); // gsl_inter_new(zv,rv,this->z_max);
    size_t n_random=this->catalogue._NOBJS();
#ifdef _USE_OMP_
    int jthread;
    omp_set_num_threads(Nran_files);
    vector<size_t>vseeds(Nran_files,0);
    for(size_t i=0;i<vseeds.size();++i)
      vseeds[i]=3+static_cast<size_t>(i+27*i*i)*56145;
#endif
#ifdef _USE_OMP_
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
        // To write them in binary
      //  ofstream outStream;
//        string out_bin=new_file+".dat";
 //       outStream.open(out_bin.c_str(), ios::binary|ios::out);
        vector<real_prec>aPzdist(this->params._Nbins_redshift(),0);
        vector<real_prec>aPcdist(this->params._Nbins_color(),0);
        vector<real_prec>aPmdist(this->params._Nbins_Mstellar(),0);
#ifdef _USE_OMP_
        jthread=omp_get_thread_num();
        gsl_rng_default_seed=vseeds[jthread];
#else
        gsl_rng_default_seed=12322;
#endif
        rng_t = gsl_rng_mt19937;//_default;
        gBaseRand = gsl_rng_alloc (rng_t);
        size_t count=0;
        while(count<=n_random)
           {
            real_prec rd = rmax*pow(gsl_rng_uniform(gBaseRand), 1./3.);
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
              if((ra<=this->params._RA_max()*factor && ra>= this->params._RA_min()*factor) && (dec<=this->params._DEC_max() *factor && dec>=this->params._DEC_min()*factor) )
                 {
#endif
                  real_prec zr = static_cast<double>(gsl_spline_eval (spline_rz, rd, acc_rz));
                  real_prec proba = 0 ;
                  if(zr>=zbins[zbins.size()-1])
                      proba = Pzdist[zbins.size()-1];
                  else if(zr<zbins[0])
                     proba = Pzdist[0];
                  else
                    proba =  static_cast<double>(gsl_spline_eval (spline_p, zr, acc_zp)) ;
                  real_prec xr = gsl_rng_uniform (gBaseRand);
                  if(proba>=xr)
                      {
                       int iz=    get_bin(zr,zbins[0],zbins.size(), delta_redshift, true);
                       aPzdist[iz]++;
                       bool accept=false;
                       while(false==accept)
                        {
                           real_prec new_color=this->params._Color_min()+(this->params._Color_max()-this->params._Color_min())*gsl_rng_uniform(gBaseRand);
                           real_prec new_ms=this->params._Mstellar_min()+(this->params._Mstellar_max()-this->params._Mstellar_min())*gsl_rng_uniform(gBaseRand);
                           int indexc=get_bin(new_color,this->params._Color_min(),this->params._Nbins_color(), delta_color, true);
                           int indexm=get_bin(new_ms,this->params._Mstellar_min(),this->params._Nbins_Mstellar(), delta_Mstellar, true);
                           size_t id=index_3d(iz, indexc, indexm,this->params._Nbins_color(), this->params._Nbins_Mstellar());
                           real_prec proba_mz=Pzcmdist_normal[id];
                           real_prec xm = gsl_rng_uniform (gBaseRand);
                           if(proba_mz>=xm)
                             {
                               real_prec ra_deg=180.*ra/M_PI;
                               real_prec dec_deg=180.*dec/M_PI;
                               accept=true;
                               aPcdist[indexc]++;
                               aPmdist[indexm]++;
                               //outStream.write(reinterpret_cast<char*>(&zr), sizeof(zr));
                               //outStream.write(reinterpret_cast<char*>(&ra_deg), sizeof(ra_deg));
                               //outStream.write(reinterpret_cast<char*>(&dec_deg), sizeof(dec_deg));
                               //outStream.write(reinterpret_cast<char*>(&new_color), sizeof(new_color));
                               //outStream.write(reinterpret_cast<char*>(&new_ms), sizeof(new_ms));
                               vdina<<zr<<"\t"<<ra_deg<<"\t"<<dec_deg<<"\t"<<new_color<<"\t"<<new_ms<<endl;
                           }
                       }
                       count++;
                    }
                 }
               }
        } // closes while
       vdina.close();
   //    outStream.close();
       odn=this->params._Output_directory()+"dNdz.txt_ranfile"+to_string(ir);
       this->File.write_to_file(odn, zbins, aPzdist);
       odn=this->params._Output_directory()+"dNdc.txt_ranfile"+to_string(ir);
       this->File.write_to_file(odn, cbins, aPcdist);
       odn=this->params._Output_directory()+"dNdMs.txt_ranfile"+to_string(ir);
       this->File.write_to_file(odn, mbins, aPmdist);
       gsl_rng_free (gBaseRand);
     }// closes loop over files
#ifdef _USE_OMP_
    }
#endif
     this->So.DONE();
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_intervals_multiscale(string prop)
{
    size_t Nft=this->params._Nft();
    gsl_vector *prop_aux=gsl_vector_alloc(this->catalogue._NOBJS());
    if(prop=="_VMAX_")
      for(size_t ig=0;ig<this->catalogue._NOBJS();++ig )
         gsl_vector_set(prop_aux,ig,log10(this->catalogue.vmax_at(ig)));
    if(prop=="_MASS_")
      for(size_t ig=0;ig<this->catalogue._NOBJS();++ig )
         gsl_vector_set(prop_aux,ig,log10(this->catalogue.mass_at(ig)));
    gsl_sort_vector(prop_aux) ;   // sort the vmax and correspondingly their associated the gal id
    size_t counter=0;
    for(int il=0;il<this->params._Number_of_MultiLevels();++il)
     {
        size_t N_level=pow(this->params.get_Nft_MultiLevels(il),3);
        this->params.set_PropThreshold_MultiLevels(il,pow(10,gsl_vector_get(prop_aux,this->catalogue._NOBJS()-(N_level+counter))));
        this->params.set_Ntracers_MultiLevels(il,N_level);
        So.message_screen("Number of reference properties to be assigned at individual level=",N_level);
        So.message_screen("Minimum property value in level",pow(10,gsl_vector_get(prop_aux,this->catalogue._NOBJS()-(N_level+counter))));
        counter=N_level;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_mean_bias_relation(s_info_in_bins &bias_info)
{
  So.enter(__PRETTY_FUNCTION__);
  int Nbins=bias_info.s_size();
  vector<int>ibin(NOBJS,0);
  vector<real_prec>prop(NOBJS,0);
  real_prec min_p=0;
  real_prec max_p=0;
  if(bias_info.name_info==_VMAX_)
   {
      min_p=log10(this->params._VMAXmin());
      max_p=log10(this->params._VMAXmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NOBJS;++i)
        prop[i]=log10(this->catalogue.vmax_at(i));
   }
  else if (bias_info.name_info==_MASS_)
   {
      min_p=this->params._LOGMASSmin();
      max_p=this->params._LOGMASSmax();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (size_t i=0;i<NOBJS;++i)
        prop[i]=log10(this->catalogue.mass_at(i));
   }
  else if(bias_info.name_info==_SPIN_ || bias_info.name_info==_SPIN_BULLOCK_)
   {
      min_p=log10(this->params._SPINmin());
      max_p=log10(this->params._SPINmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NOBJS;++i)
        prop[i]=log10(this->catalogue.spin_bullock_at(i));
   }
  else if(bias_info.name_info==_RS_ )
   {
      min_p=log10(this->params._RSmin());
      max_p=log10(this->params._RSmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NOBJS;++i)
        prop[i]=log10(this->catalogue.rs_at(i));
   }  
  else if(bias_info.name_info==_CONCENTRATION_ )
   {
      min_p=log10(this->params._CONCENTRATIONmin());
      max_p=log10(this->params._CONCENTRATIONmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NOBJS;++i)
        prop[i]=log10(this->catalogue.concentration_at(i));
   }  
  real_prec delta=(max_p-min_p)/static_cast<real_prec>(Nbins);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for (size_t i=0;i<Nbins;++i)
     bias_info.vbin[i]=min_p+(i+0.5)*delta;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<NOBJS;++i)
    if(prop[i]<max_p && prop[i]>= min_p)
      {
        int bin=get_bin(prop[i],min_p,Nbins,delta,true);
        ibin[i]=bin;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.vq1[bin]+=this->catalogue.bias_at(i);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.i_ncbin[bin]++;
     }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<Nbins;++i)  // mean
    if(bias_info.i_ncbin[i]>0)  
      bias_info.vq1[i]=bias_info.vq1[i]/static_cast<real_prec>(bias_info.i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<NOBJS;++i)
    if(prop[i]<max_p && prop[i]>= min_p)
      {
        int bin=ibin[i];
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.vq2[bin]+=pow(this->catalogue.bias_at(i)-bias_info.vq1[bin],2);
    }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<Nbins;++i)  // standard deviation
    if(bias_info.i_ncbin[i]>0)  
     bias_info.vq2[i]=sqrt(bias_info.vq2[i])/static_cast<real_prec>(bias_info.i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (size_t i=0;i<Nbins;++i) // error in the mean
    if(bias_info.i_ncbin[i]>0)  
      bias_info.vq3[i]=bias_info.vq2[i]/sqrt(static_cast<real_prec>(bias_info.i_ncbin[i]));
So.leaving(__PRETTY_FUNCTION__);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::get_mean_secondary_bias_relation(vector<s_info_in_bins> &sec_bias_info)
{
  So.enter(__PRETTY_FUNCTION__);
  this->params.set_Number_of_bins_equal_number_tracers(4);// divide in quartiles of thes econdary property
  string primary_prop=sec_bias_info[0].name_info;
  string secondary_prop=sec_bias_info[0].name_info_sec;
  if(secondary_prop==_SPIN_BULLOCK_)
    this->params._resizeSPINBULLOCKbins(4+1);
  else if(secondary_prop==_RS_)
    this->params._resizeRSbins(4+1);
  else if(secondary_prop==_CONCENTRATION_)
    this->params._resizeCONCENTRATIONbins(4+1);
  vector<int>quartiles={0,1,2,3,4};
  this->get_intervals_equal_number_aux(secondary_prop);
  int Nbins=sec_bias_info[0].vbin.size();
  vector<int>ibin(NOBJS,0);
  vector<bool>used(NOBJS,false);
  vector<real_prec>p_prop(NOBJS,0);
// ----------------------------------------------------------
// Get min and max of primeary proeprty, read fom parameter file:
  real_prec min_p, max_p, min_s, max_s;
// ----------------------------------------------------------
  if(primary_prop==_VMAX_)
   {
     min_p=log10(this->params._VMAXmin());
     max_p=log10(this->params._VMAXmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<NOBJS;++i)
     p_prop[i]=log10(this->catalogue.vmax_at(i));
   }
  else if (primary_prop==_MASS_)
   {
     min_p=this->params._LOGMASSmin();
     max_p=this->params._LOGMASSmax();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for (size_t i=0;i<NOBJS;++i)
       p_prop[i]=log10(this->catalogue.mass_at(i));
  }
     real_prec delta=(max_p-min_p)/static_cast<real_prec>(Nbins);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<Nbins;++i)
      sec_bias_info[0].vbin[i]=min_p+(i+0.5)*delta;
// ----------------------------------------------------------
  vector<real_prec>s_prop(NOBJS,0);
  if (secondary_prop ==_SPIN_ || secondary_prop==_SPIN_BULLOCK_)
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for (size_t i=0;i<NOBJS;++i)
       s_prop[i]=log10(this->catalogue.spin_bullock_at(i));
   }
  else if (secondary_prop==_RS_ )
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<NOBJS;++i)
      s_prop[i]=log10(this->catalogue.rs_at(i));
   } 
  else if (secondary_prop==_CONCENTRATION_ )
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<NOBJS;++i)
      s_prop[i]=log10(this->catalogue.concentration_at(i));
   } 
// ----------------------------------------------------------
  for(int iq=1;iq<quartiles.size();++iq)// start from 1; 0 is fo the full range
  {
    int quart=iq;//quartiles[iq];
    if (secondary_prop ==_SPIN_ || secondary_prop==_SPIN_BULLOCK_)
     {
       min_s=log10(this->params._SPINBULLOCKbins_min(quart));
       max_s=log10(this->params._SPINBULLOCKbins_max(quart));
     }
    else if (secondary_prop==_RS_ )
    {
      min_s=log10(this->params._RSbins_min(quart));
      max_s=log10(this->params._RSbins_max(quart));
    } 
    else if (secondary_prop==_CONCENTRATION_ )
    {
      min_s=log10(this->params._CONCENTRATIONbins_min(quart));
      max_s=log10(this->params._CONCENTRATIONbins_max(quart));
    } 

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (size_t i=0;i<NOBJS;++i)
      if(p_prop[i]<max_p && p_prop[i]>= min_p)
       if(s_prop[i]<max_s && s_prop[i]>= min_s)
        {
          int bin=get_bin(p_prop[i],min_p,Nbins,delta,true);
          ibin[i]=bin;      

#ifdef _USE_OMP_
#pragma omp atomic
#endif
          sec_bias_info[quart].vq1[bin]+=this->catalogue.bias_at(i);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          sec_bias_info[quart].i_ncbin[bin]++;
        }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (size_t i=0;i<Nbins;++i)  // mean
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq1[i]=sec_bias_info[quart].vq1[i]/static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<NOBJS;++i)
        if(p_prop[i]<max_p && p_prop[i]>= min_p)
          if(s_prop[i]<max_s && s_prop[i]>= min_s)
            {
              int bin=ibin[i];
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              sec_bias_info[quart].vq2[bin]+=pow(this->catalogue.bias_at(i)-sec_bias_info[quart].vq1[bin],2);
          }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (size_t i=0;i<Nbins;++i)  // standard deviation
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq2[i]=sqrt(sec_bias_info[quart].vq2[i])/static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (size_t i=0;i<Nbins;++i) // error in the mean
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq3[i]=sec_bias_info[quart].vq2[i]/sqrt(static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]));
 }//closing for(int iq=1;iq<quartiles.size();++iq)
So.leaving(__PRETTY_FUNCTION__);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This metods takes a box and converts it int a mock, although with no evolution
void HaloTools::snap_to_mock(bool set_dndz){

// convert equatorial to cartessian
// Assume a cosmological model and transform distance to redshift
// Use the infomration of the pculiar velocities along the line of sight and add pecular redshift
// Define a redshift distribution and select objects accordinly

int Nz=100;
vector<gsl_real>rrc(Nz,0);
vector<gsl_real>zzc(Nz,0);
this->Cosmo.set_cosmo_pars(this->s_cosmo_pars);
for(size_t i=0;i<rrc.size();++i)
{
  zzc[i]=0.5*(i+0.5)/static_cast<real_prec>(Nz);
  rrc[i]=this->Cosmo.comoving_distance(zzc[i]);
}

gsl_interp_accel *acc = gsl_interp_accel_alloc();
gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,rrc.size());
gsl_spline_init(spline,&rrc[0],&zzc[0],rrc.size());

// Asignar coordenadas esféricas y redshift cosmológico
real_prec ra, dec, r;
for(size_t i=0;i<this->catalogue._NOBJS();++i)
{
  // Place observer in the center of the box
  real_prec x=this->catalogue.coord1_at(i)-0.5*this->params._Lbox();
  real_prec y=this->catalogue.coord2_at(i)-0.5*this->params._Lbox();
  real_prec z=this->catalogue.coord3_at(i)-0.5*this->params._Lbox();
  // Read velocities
  real_prec vx=this->catalogue.vel1_at(i);
  real_prec vy=this->catalogue.vel2_at(i);
  real_prec vz=this->catalogue.vel3_at(i);
  // Convert cartesian to equatorial 
  cartesian_to_equatorial(x,y,z,ra, dec, r);
  this->catalogue.set_coord1(ra,i);
  this->catalogue.set_coord2(dec,i);
  this->catalogue.set_coord3(r,i);
  // Get cosmological redshift
  real_prec zc=gsl_spline_eval(spline,r, acc);
  // Add peculiar component to distance
  real_prec r_pec=r+(1+zc)*(x*vx+y*vy+z*vz)/this->Cosmo.Hubble_function(zc);
  // Obtain "observed" redshift
  real_prec zp= gsl_spline_eval(spline,r_pec, acc);
  this->catalogue.set_cosmological_redshift(zc,i);  // assign cosmological redhshift
  this->catalogue.set_redshift(zp,i);               // assign observed redshift
}


if(true==set_dndz){
  // aca tendré que elegir los redshifts con respecto a una dndz.
  // genera la nbar normalizada al máximo con alguna dndz
  // escoge el z de cada objeto son esa probabilidad
  //
}

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HaloTools::assign_idgrid_to_tracers()
{

  So.enter(__PRETTY_FUNCTION__);
  this->catalogue.resize_GridID(this->catalogue._NOBJS());

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0; i<this->catalogue._NOBJS();++i)
  {
    size_t idval = grid_ID(&this->box, this->catalogue.coord1_at(i),this->catalogue.coord2_at(i),this->catalogue.coord3_at(i));
    this->catalogue.set_GridID(idval, i);
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::assign_cwc_to_tracers(vector<WebType>&cwc)
{
 
  So.enter(__PRETTY_FUNCTION__);

  So.message_screen("Assigning CWT to tracers");
  this->catalogue.resize_gal_cwt(this->catalogue._NOBJS());

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<this->catalogue._NOBJS();++i)
    this->catalogue.set_gal_cwt(cwc[this->catalogue.GridID_at(i)], i);

  So.DONE();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::assign_tidal_anisotropy_to_tracers(vector<real_prec>&tidal, bool tr)
  {
    So.enter(__PRETTY_FUNCTION__);

    if(!tr) // false for tr
      {
      this->catalogue.resize_tidal_anisotropy(this->catalogue._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        this->catalogue.set_tidal_anisotropy(tidal[this->catalogue.GridID_at(i)], i);
      }
  else if(tr) // true for dm
    {
      this->catalogue.resize_tidal_anisotropy_dm(this->catalogue._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();++i)
        {
          ULONG index= this->catalogue.GridID_at(i);
          real_prec tidal_a=tidal[index];
          this->catalogue.set_tidal_anisotropy_dm(tidal_a,i);
        }
    }
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void HaloTools::Cat2Map(Healpix_Map<healpix_real>&healpix_map, bool data)
{
  if(data)
    So.enter(__PRETTY_FUNCTION__);
  
   real_prec nau=0;
   size_t ngal=0;
   healpix_map.fill(0);


  if(this->params._statistics()=="ARF")
    {
      real_prec mid_redshift_a=0;
      real_prec mid_redshift_b=0;
      real_prec bar_z=0;
      this->mid_redshift=this->z_min[this->index_zbin]+(this->z_max[this->index_zbin]-this->z_min[this->index_zbin])*0.5;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mid_redshift_a,mid_redshift_b)
#endif
      for(size_t i=0;i<this->catalogue._NOBJS();i++)
        {
          real_prec redshift=this->catalogue.redshift_at(i);
          real_prec xx=(redshift -this->mid_redshift)/this->params._sigma_arf();
          real_prec gauss=exp(-0.5*xx*xx);
          mid_redshift_a+=redshift*gauss;
          mid_redshift_b+=gauss;
        }
      this->bar_redshift=mid_redshift_a/mid_redshift_b;


      for(size_t i=0;i<this->catalogue._NOBJS();i++) // Loop over the HaloToolsCl
       {
          size_t pixel=this->catalogue.galPIXEL_at(i);
          if(1==this->pixmask[pixel])
          {
            real_prec gauss_weight=1.0;
            real_prec redshift=this->catalogue.redshift_at(i);
            if(redshift < this->z_max[this->index_zbin] && redshift >=this->z_min[this->index_zbin])
                {
                  real_prec xx=(this->catalogue.redshift_at(i)-this->mid_redshift)/this->params._sigma_arf();
                  gauss_weight=exp(-0.5*xx*xx);
                  nau+=gauss_weight;
                  ngal++;
                  healpix_map[pixel]+= (redshift - this->bar_redshift)*gauss_weight;
              }
        }
      }
  }

  else if("Cl"==this->params._statistics())
    {
      for(size_t i=0;i<this->catalogue._NOBJS();i++) // Loop over the HaloToolsCl
       {
         size_t pixel=this->catalogue.galPIXEL_at(i);
         if(1==this->pixmask[pixel])
          {
           real_prec redshift=this->catalogue.redshift_at(i);
           if(redshift < this->z_max[this->index_zbin] && redshift >=this->z_min[this->index_zbin])
             {
               healpix_real weight=1.0;
               ngal++;
               healpix_map[pixel]+=weight;
             }
	       }
	   }
  }

  this->ngals_in_zbin_from_map=0;
  real_prec nfm=0;

  #ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nfm)
#endif
  for(size_t i=0;i<healpix_map.Npix();++i)
      nfm+=healpix_map[i];

  if(data)
  {
  this->ngals_in_zbin=ngal; // Actual number o galaxies in the zbin

  this->ngals_in_zbin_from_map=nfm; // Check from the map
  
  if(this->params._statistics()=="ARF")
    this->weighted_ngals_in_zbin=nau; // Weighted number of galaxies. In bin-based tomographic analysis this equals ngals_in_zbins
  else
    this->weighted_ngals_in_zbin=ngal; // Weighted number of galaxies. In bin-based tomographic analysis this equals ngals_in_zbins


  this->mean_number_galaxies_pix=static_cast<real_prec>(this->weighted_ngals_in_zbin)/(static_cast<real_prec>(this->n_observed_pixels));

  this->rms_ngal=healpix_map.rms();

  if(this->params._generate_fits_files())
    {
      string file_map_fits;
      if(data)
          file_map_fits=this->params._output_file_fits()+"_zbin_"+to_string(this->index_zbin)+"_"+this->params._statistics()+".fits";
      else
          file_map_fits=this->params._output_file_fits()+"_zbin_"+to_string(this->index_zbin)+"_"+this->params._statistics()+"_rmaps.fits";

      So.message_screen("Writing map in fits format in file ", file_map_fits);
      if (std::filesystem::exists(file_map_fits)){
          std::filesystem::remove(file_map_fits);
        } 
      write_Healpix_map_to_fits(file_map_fits, healpix_map,(PDT)9);
    }
  }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void HaloTools::set_zbins()
{
  this->z_min.resize(this->params._n_z_bins_tomography()+1,0);
  this->z_max.resize(this->params._n_z_bins_tomography()+1,0);
  this->zminmax.resize(this->params._n_z_bins_tomography()+1);

  this->z_min[0]=this->params._min_redshift_tomography();
  this->z_max[0]=this->params._max_redshift_tomography();
  real_prec Dz=(this->params._max_redshift_tomography()-this->params._min_redshift_tomography())/(static_cast<real_prec>(this->params._n_z_bins_tomography()));

  if(this->params._define_z_bins()=="number")
    this->get_zbins_same_ngal(this->params._n_z_bins_tomography(),this->params._min_redshift_tomography(),this->params._max_redshift_tomography(),this->zminmax);

  for(size_t i=1;i<=this->params._n_z_bins_tomography();++i){
      z_min[i]= (this->params._define_z_bins()=="number"? this->zminmax[i][0]:this->params._min_redshift_tomography()+(i-1)*Dz);
      z_max[i]= (this->params._define_z_bins()=="number"? this->zminmax[i][1]:this->params._min_redshift_tomography()+(i)*Dz);
    }

#ifdef _FULL_VERBOSE_
  So.message_screen("REDSHIFT BINS SELECTED: ");
  if(this->params._n_z_bins_tomography()>1)
    for(size_t i=0;i<=this->params._n_z_bins_tomography();++i)
      cout<<i<<"  ("<<z_min[i]<<"-"<<z_max[i]<<")"<<endl;
  else cout<<"0  ("<<z_min[0]<<"-"<<z_max[0]<<")"<<endl;
  cout<<RESET<<endl;
#endif

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void HaloTools::get_zbins_same_ngal(int nbins, real_prec zmn, real_prec zmx, vector<vector<real_prec> >&zzv){
  // Construct the dNdz with many bins
  size_t na=20000; //This value is critical to avoid seg foults. The higher, the best.
  real_prec delta=(zmx-zmn)/(static_cast<real_prec>(na));
  vector<int>dn(na,0);
  for(size_t i=0;i<this->catalogue._NOBJS();++i){
    if(this->catalogue.redshift_at(i) <zmx && this->catalogue.redshift_at(i) >=zmn){
      int iza=floor((this->catalogue.redshift_at(i)-zmn)/delta);
      dn[iza]++;
    }
  }

  int Nca=static_cast<int>(floor(this->catalogue._NOBJS()/(static_cast<int>(nbins)))); //Desired number of galaxies per redshift bin:
#ifdef _FULL_VERBOSE_
  cout<<"Number of galaxies per redshift bin = "<<Nca<<endl;
#endif

  vector<real_prec>zan(na,0);
  for(size_t i=0;i<dn.size();++i)
    zan[i]=zmn+(i+0.5)*delta;

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
      for(size_t i=0;i<dn.size();++i){
        caa+=dn[i];  //Cumulative number of galaxies
        if((caa>=Nca*(ib-1)) &&  (caa<Nca*ib)) zaux.push_back(zan[i]);
      }
      zzv[ib][0]=zaux[0]-0.5*delta;  //Allocate the this->params._zmin() of the ib zbin
      zzv[ib][1]=zaux[zaux.size()-1]+0.5*delta; //Allocate the this->params._zmax() of the ib zbin
      zaux.clear();
    }
  }
  dn.clear();
  return ;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::read_nbar_photo(){  // TO BE MOVED OUT
  vector< vector<real_prec> > nbar;
  this->File.read_file(this->params._file_nbar(),nbar);

  for(size_t i=0;i<nbar.size();++i)this->dz.push_back(nbar[i][0]);
  for(size_t i=0;i<nbar.size();++i)this->nbar_photo.push_back(nbar[i][1]);
  if(this->params._use_min_max_redshift_from_nbar_file()){
    this->min_redshift_cl= nbar[0][0];
    this->max_redshift_cl= nbar[this->dz.size()-1][0];
  }
  else {
    if(this->max_redshift_cl > nbar[this->dz.size()-1][0]){
#ifdef _FULL_VERBOSE_
      std::cout<<RED<<"Warning: Maximum z in parameter file ("<<this->params._zmax()<<") greater than maximum value found in the dNdz ("<<this->dz[this->dz.size()-1] <<"). Setting value from HaloToolsCl"<<RESET<<endl;
#endif
      this->min_redshift_cl= nbar[this->  dz.size()-1][0];
    }
    else if(this->params._zmin() < nbar[0][0]){
#ifdef _FULL_VERBOSE_
      std::cout<<RED<<"Warning: Minimum z in parameter file ("<<this->params._zmin()<<") smaller than minimum value found in the dNdz ("<<nbar[0][2]<<"). Setting value from HaloToolsCl"<<RESET<<endl;
#endif
      this->min_redshift_cl= nbar[0][0];
    }
  }
  nbar.clear();
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void HaloTools::read_mask(){ // TO BE MOVED OUT

  std::string mask_f=(this->params._type_of_sky_coverage()=="full_sky"? this->params._input_file_mask_fullsky(): this->params._input_file_mask());
  std::vector<real_prec>mask2;
  this->n_pixels=this->File.read_file(mask_f, mask2,_NTHREADS_);
  this->n_columns_mask =static_cast<size_t>(mask2.size()/n_pixels);

  this->params.set_nside(static_cast<size_t>(sqrt(this->n_pixels/12)));
  So.message_screen("Using Nside =",this->params._nside() );

  // nr=nrings, have to compute it here for nside
  // is going to be the size of the mask and has not been yet computed
  size_t nr=4*this->params._nside()-1;

  // copy the mask and find the angles
  Healpix_Map<healpix_real>mask_aux(ilog2(this->params._nside()), RING);
  this->theta_new.resize(nr,0);
  this->phi_new.resize(mask_aux.Npix(),0);
  this->pixmask.resize(mask_aux.Npix(),0);

  n_observed_pixels=0;

  mask_aux.fill(0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(size_t i=0;i<n_pixels;i++)
        mask_aux[i]=static_cast<healpix_real>(mask2[this->params._i_mask_flag()+this->n_columns_mask*i]);

  if(this->params._generate_fits_files())
  {
    if (std::filesystem::exists(this->params._Output_directory()+"MASK.fits")){
       std::filesystem::remove(this->params._Output_directory()+"MASK.fits");
    } 
    write_Healpix_map_to_fits(this->params._Output_directory()+"MASK.fits", mask_aux, (PDT)9);
  }
  int ip=0;
  size_t nop=0;
  for(size_t ir=0;ir<nr;ir++)
    {
      size_t Npixels_ring=npix_ring(this->params._nside(), ir);
      for(size_t ipr=0;ipr<Npixels_ring;++ipr)
        {
          pointing point_rev;
          point_rev=mask_aux.pix2ang(ip);
          this->theta_new[ir]=point_rev.theta;
          this->phi_new[ip]=point_rev.phi;
          int pmask=static_cast<int>(mask2[this->params._i_mask_flag()+this->n_columns_mask*ip]);
          this->pixmask[ip]=pmask;
          ip++;
          if(1==pmask)
            nop++;
        }
    }
  So.message_screen("Number of observed pixels ", nop);
  mask2.clear(); mask2.shrink_to_fit();
  this->n_observed_pixels=nop;
  this->sky_fraction=static_cast<real_prec>(n_observed_pixels)/static_cast<real_prec>(n_pixels);
  this->area_pixel=4.*M_PI/static_cast<real_prec>(n_pixels);
  this->area_survey=static_cast<real_prec>(n_observed_pixels)*area_pixel;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void HaloTools::perform_extra_taks_on_tracers()
{
  
  this->So.enter(__PRETTY_FUNCTION__);

  CoordinateSystem systemc;
  if(this->catalogue._type_of_object()=="TRACER" || this->catalogue._type_of_object()=="TRACER_REF" || this->catalogue._type_of_object()=="TRACER_MOCK" || this->catalogue._type_of_object()=="TRACER_MOCK_ONLY_COORDS"|| this->catalogue._type_of_object()=="TRACER_REF_ONLY_COORDS")
    systemc = this->params._sys_of_coord_g();
  else if(this->catalogue._type_of_object()=="RANDOM")
    systemc = this->params._sys_of_coord_r();


  if(systemc ==CoordinateSystem::CART)
    {
      So.message_screen("Getting grid-ID from tracer coordinates");
      this->catalogue.resize_GridID(this->catalogue._NOBJS());
      this->catalogue.resize_GridID_n(this->catalogue._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for 
#endif
    for(size_t i =0; i< this->catalogue._NOBJS(); ++i)
      {
        real_prec x=this->catalogue.coord1_at(i);
        real_prec y=this->catalogue.coord2_at(i);
        real_prec z=this->catalogue.coord3_at(i);
        this->catalogue.set_GridID(grid_ID(&box, x,y,z),i);
        this->catalogue.set_GridID_n(grid_ID(&box_n, x,y,z),i);
      }
     So.DONE();
  }


  // ALLOCATE HEALPIUX INDEX
  if((this->params._statistics()=="Cl" || this->params._statistics()=="ARF") && systemc!=CoordinateSystem::CART) 
  {
     So.message_screen("Allocating healpix index");
     Healpix_Map<real_prec>map_aux(ilog2(this->params._nside()), RING);
     this->catalogue.resize_galPIXEL(this->catalogue._NOBJS());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(size_t i=0; i < this->catalogue._NOBJS() ;++i)
      {
        pointing point;
        point.phi=this->catalogue.coord1_at(i)*CFACTOR;
        point.theta=0.5*M_PI-this->catalogue.coord2_at(i)*CFACTOR;
        this->catalogue.set_galPIXEL(static_cast<size_t>(map_aux.ang2pix(point)), i);
      }
    So.DONE();
  }



  if(params.input_sections.HaloAnalysis)
   {


#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_

      if(this->params._set_bins_equal_number_tracers_main_property())// not often done, as it will be rather fixed to check diferent redshifts
      {
        int nbins_prop=1;
        vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
        vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
        nbins_prop=this->params._Number_of_bins_equal_number_tracers_main_property();
        this->Number_of_tracers_in_vmax_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
        this->get_intervals_equal_number("_VMAX_",min_aux,max_aux);

        for(size_t i=0;i<=nbins_prop;++i)
          this->params.set_VMAXbins_max(i,log10(max_aux[i]));
        for(size_t i=0;i<=nbins_prop;++i)
            this->params.set_VMAXbins_min(i,log10(min_aux[i]));

        this->So.message_screen("Number of bins in this property =", this->params._Number_of_bins_equal_number_tracers_main_property());
        this->So.message_screen("\tMin log Vmax  =", log10(get<static_cast<int>(PropStats::MIN)>(this->catalogue.info_vmax))));
        this->So.message_screen("\tMax log Vmax  =", log10(get<static_cast<int>(PropStats::MAX)>(this->catalogue.info_vmax))));
      }
      else
      {
        this->So.message_screen("\tMin log Vmax  =", log10(get<static_cast<int>(PropStats::MIN)>(this->catalogue.info_vmax))));
        this->So.message_screen("\tMax log Vmax  =", log10(get<static_cast<int>(PropStats::MAX)>(this->catalogue.info_vmax))));
        this->Number_of_tracers_in_vmax_bins.resize(this->params._NVMAXbins_power(),0);
        So.message_screen("\tVmax-bins defined fixed width (parameter file)");
      }
#else
      if(this->params._set_bins_equal_number_tracers())
      {
       int nbins_prop=1;
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       nbins_prop=this->params._Number_of_bins_equal_number_tracers_main_property();
       this->get_intervals_equal_number("_VMAX_",min_aux,max_aux);
       for(size_t i=0;i<=nbins_prop;++i)
         this->params.set_VMAXbins_max(i,max_aux[i]);
       for(size_t i=0;i<=nbins_prop;++i)
          this->params.set_VMAXbins_min(i,min_aux[i]);
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
       this->So.message_screen("\tMin Vmax  =", get<static_cast<int>(PropStats::MIN)>(this->catalogue.info_vmax));
       this->So.message_screen("\tMax Vmax  =", get<static_cast<int>(PropStats::MAX)>(this->catalogue.info_vmax));
     }
#endif

    

// ALLOCATE STUFF FOR MULTISCALING
#ifdef _MULTISCALE_
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
#ifndef _USE_FIXED_MULTISCALE_LEVELS_
if(type_of_object=="TRACER_REF")
    {
      vector<size_t> n_auxtr(this->params._Number_of_MultiLevels(),0);

      size_t part_level_counter=0;
      for(int il=0; il<this->params._Number_of_MultiLevels();++il )
        {
          size_t N_AUX=  n_auxtr[il];
          part_level_counter+=N_AUX;
          this->params.set_Ntracers_MultiLevels(il, N_AUX);
#ifdef _VERBOSE_CAT_
          So.message_screen("Number of reference tracers in multi-scale level", il+1,"=", N_AUX);
#endif
      }
      n_auxtr.clear();n_auxtr.shrink_to_fit();
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of reference properties to be assigned at individual level=",NLINES-part_level_counter);
      std::cout<<std::endl;
#endif
      }
      So.DONE();
#endif
#endif // end of _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
#endif // END OF _MULTISCALE_
   
  // ********************************************************************************************************************************

#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
      So.DONE();
      if(this->params._set_bins_equal_number_tracers_main_property()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
       {
         vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
         vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
         this->Number_of_tracers_in_mass_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
       this->get_intervals_equal_number("_MASS_",min_aux,max_aux);

       So.message_screen("Mass bins with equal number of tracers identified");
       for(size_t i=0;i<=this->Number_of_tracers_in_mass_bins.size();++i)
       {
          this->params.set_MASSbins_max(i,max_aux[i]);
          this->params.set_MASSbins_min(i,min_aux[i]);
          So.message_screen("\tmin mass", min_aux[i]);
          So.message_screen("\tmax mass", max_aux[i]);
          So.message_screen("\tNtracers", this->Number_of_tracers_in_mass_bins[i]);
          cout<<endl;
       }
     }
     else // if bins in mass are to bhe define with constant width:
      {
          this->Number_of_tracers_in_mass_bins.resize(this->params._NMASSbins_power(),0);
          So.message_screen("\tMass-bins defined fixed width (parameter file)");
      }

#else // If mass is not primary property then
      if(this->params._set_bins_equal_number_tracers()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
      {

       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       this->Number_of_tracers_in_mass_bins.resize(this->params._Number_of_bins_equal_number_tracers()+1,0);
       this->get_intervals_equal_number("_MASS_",min_aux,max_aux);
       for(size_t i=0;i<=this->Number_of_tracers_in_mass_bins.size();++i)
       {
          this->params.set_MASSbins_max(i,max_aux[i]);
          this->params.set_MASSbins_min(i,min_aux[i]);
       }
       this->So.message_screen("\tMin Mass  =", get<static_cast<int>(PropStats::MIN)>(this->catalogue.info_mass));
       this->So.message_screen("\tMax Mass  =", get<static_cast<int>(PropStats::MIN)>(this->catalogue.info_mass));
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
     }
    }
#endif

  // ********************************************************************************************************************************
  // THIS SECTION FILLS, IF REQUESTED, INFORMATION ON ENVIROMENTAL PROPERTIES SUCH AS MACH, DAH, TIDAL,
  if(this->params._Get_tracer_local_mach_number() || this->params._Get_local_overdensity() || this->params._Get_tracer_local_dach_number())
  {

 #ifndef _USE_VELOCITIES_TRACERS_
    So.message_warning("Option _USE_VELOCITIES_TRACERS_ is not defined. Please check def.h, clean and make");
    exit(1);
#endif
#ifdef _VERBOSE_CAT_
    So.message_screen("Computing statistiscs (mach number. local clustering) for tracers");
    So.message_screen("\tcomputing on sphers of radii ", this->params._Scale_mach_number(), "Mpc/h");
#endif
#ifdef _VERBOSE_CAT_
#endif
#ifdef _USE_CHUNCKS_NC_
    this->get_local_mach_number_chuncks(this->params._Scale_mach_number()); // local mach number at a radius of 8 Mpc/h
#else
    this->get_local_mach_number(this->params._Scale_mach_number()); // local mach number at a radius of 8 Mpc/h
#endif
    if(this->params._Get_tracer_local_mach_number())
     {
       this->min_mach=get_min("_MACH_");
       this->max_mach=get_max("_MACH_");
       this->So.message_screen("\tMin Mach  =", this->min_mach);
       this->So.message_screen("\tMax Mach =", this->max_mach);
     }
    if(this->params._Get_local_overdensity())
    {
       this->min_local_overdensity=get_min("_LOCAL_OVERDENSITY_");
       this->max_local_overdensity=get_max("_LOCAL_OVERDENSITY_");
       this->So.message_screen("\tMin Local Clustering =", this->min_local_overdensity);
       this->So.message_screen("\tMax Local Clustering =", this->max_local_overdensity);
    }
   if(this->params._Get_tracer_local_dach_number())
    {
       this->min_dach=get_min("_DACH_");
       this->max_dach=get_max("_DACH_");
       this->So.message_screen("\tMin Dach  =", this->min_dach);
       this->So.message_screen("\tMax Dach =", this->max_dach);
    }
  }
  
  this->mean_number_density=static_cast<real_prec>(this->catalogue._NOBJS())/pow(this->box.Lbox,3);
}

// *********************************************************************************************
#ifdef _MULTISCALE_
#ifdef _USE_FIXED_MULTISCALE_LEVELS_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
   this->get_intervals_multiscale("_MASS_");
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
   this->get_intervals_multiscale("_VMAX_");
#endif
#endif
#endif
  
}


