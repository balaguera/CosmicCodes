////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<Catalog>
 * @file Catalog.cpp
 * @brief Methods of the class Catalog
 * @details The class Catalog reads and analyses an input catalog od dark matter tracers
 * @author Andres Balaguera Antolinez 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../headers/def.h"
#include "Catalog.h"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_cosmo(){
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
  for(int i=0;i<nz;++i)
    this->zcosmo[i]=zmin_inter+i*(zmax_inter-zmin_inter)/static_cast<double>(nz-1);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<nz;++i)
    this->rcosmo[i]=static_cast<gsl_real>(this->Cosmo.comoving_distance(this->zcosmo[i]));
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::analyze_cat(bool read){
  this->So.enter(__PRETTY_FUNCTION__);
  // READ THE ASCII FILE WITH THE CATALOG
  if(true==read)
    {
#if defined (_USE_ALL_PK_) || defined (_USE_MASS_CUTS_PK_)
      this->read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),pow(10,this->params._LOGMASSmin()));
#else
      this->read_catalog(this->params._Input_dir_cat()+this->params._file_catalogue(),pow(10,this->params._LOGMASSmin()),pow(10,this->params._LOGMASSmax()));
#endif
    }
  //* COMPUTE LOCAL MACH NUMBER
  if(this->params._Get_cell_local_mach_number())
    this->get_local_mach_number(true);
  //* COMPUTE DISTRIBUTION OF SEPARATIONS
  // There is a bug related to memmory in jumping from get_neighbpuro to get_distribution. UNidentified.
  if(true==this->params._get_distribution_min_separations())
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
  string fname_mass_function_Y = this->Output_directory+this->type_of_object+this->params._Name_survey()+"_abundance_R"+to_string(this->params._realization())+".txt";
  if(true==this->params._Get_prop_function_tracer())
    this->get_property_function(fname_mass_function_Y);
  //* COMPUTE ABUNDANCE IN CWT:
  if(true==this->params._Get_prop_function_tracer_cwt())
    this->get_property_function_cosmic_web_types(fname_mass_function_Y);
  // ---------------------------------------------------------------------------- //
  // these are wainting for an input parameter
  /*
    string fname_mass_function_Y = this->Output_directory+this->type_of_object+this->params._Name_survey()+"_HOD_R"+to_string(this->params._realization())+".txt";
    this->get_HOD(fname_mass_function_Y);
    this->get_HOD_web_types(fname_mass_function_Y);
  */
  //this->get_pdf_vmax("VMAX");
  //this->get_sep_distribution(PROP_THRESHOLD_MULTI_SCALE_4);
  // ---------------------------------------------------------------------------- //
  //* COMPUTE NUMBER COUNTS IN AS MESH
  if(true==this->params._Get_tracer_number_counts())
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
    for(int i=0; i<this->NOBJS; ++i)
    xy_pts_m.push_back(std::make_pair(this->Halo[i].mass, ncounts[this->Halo[i].GridID]));
    this->gp_abundance<<"plot"<<this->gp_abundance.file1d(xy_pts_m) << " w p ps 0.2 pt 7 title 'Mass'"<<endl;
    xy_pts_m.clear(); xy_pts_m.shrink_to_fit();
  */
  // ---------------------------------------------------------------------------- //
  //* COMPUTE MASS WEIGHTED FIELD
  if(true==this->params._Get_tracer_mass_field())
    if(this->params._i_mass_g()>0 && this->params._i_mass_g()<this->NCOLS)
      {
        string file=this->params._Output_directory()+this->params._Name_survey()+"_MASS_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
        this->get_density_field_grid(_MASS_,file);
      }
  // ---------------------------------------------------------------------------- //
  //* COMPUTE VMAX WEIGHTED FIELD
  if(true==this->params._Get_tracer_vmax_field())
    if(this->params._i_vmax_g()>0 && this->params._i_vmax_g()<this->NCOLS)
      {
	string file=this->params._Output_directory()+this->params._Name_survey()+"_VMAX_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
	this->get_density_field_grid(_VMAX_,file);
      }
  // ---------------------------------------------------------------------------- //
  //* COMPUTE SPIN WEIGHTED FIELD
  if(true==this->params._Get_tracer_spin_field())
    if(this->params._i_spin_g()>0 && this->params._i_spin_g()<this->NCOLS)
      {
        string file=this->params._Output_directory()+this->params._Name_survey()+"_SPIN_Nft"+to_string(this->params._Nft())+"_"+this->params._Name_survey();
        this->get_density_field_grid(_SPIN_,file);
      }
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
  for(ULONG ig=0; ig<this->NOBJS;++ig)
  {
  //        cout<<this->Halo[ig].GridID<<endl;
  this->Halo[ig].gal_cwt=this->cwclass.cwt_used[this->cwclass.get_Tclassification(this->Halo[ig].GridID)]; // thisis written in the cats,
  this->Halo[ig].local_dm=deltam[this->Halo[ig].GridID]; // thisis written in the cats,
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
  ULONG NX=10;
  real_prec delta=(dm_max-dm_min)/static_cast<real_prec>(NX);
  ULONG Ncwt=5;
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
  for(ULONG ig=0; ig<this->NOBJS;++ig)
  {
  ULONG  IGRID=this->Halo[ig].GridID;
  int ICWT=this->cwclass.cwt_used[this->cwclass.get_Tclassification(IGRID)];
  real_prec delta_dm=1+deltam[IGRID];
  int ddm_bin=get_bin(log10(delta_dm),dm_min,NX,delta, true);
  props[ICWT].mass.push_back(this->Halo[ig].mass);
  props[ICWT].vel.push_back(this->Halo[ig].vmax);
  props[ICWT].den.push_back(delta_dm);
  props[ICWT].idm.push_back(ddm_bin);
  props[0].mass.push_back(this->Halo[ig].mass);
  props[0].vel.push_back(this->Halo[ig].vmax);
  props[0].den.push_back(delta_dm);
  props[0].idm.push_back(ddm_bin);
  }
  for (int j=0;j<Ncwt;++j)
  {
  mean_props[j].mass.resize(NX,0);
  mean_props[j].vel.resize(NX,0);
  mean_props[j].den.resize(NX,0);
  mean_props[j].ntr.resize(NX,0);
  for(ULONG i=0;i<props[j].mass.size();++i)
  {
  int ddm_bin=props[j].idm[i];
  mean_props[j].mass[ddm_bin]+=props[j].mass[i];
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
  for(ULONG i=0;i<props[j].mass.size();++i)
  {
  int ddm_bin=props[j].idm[i];
  real_prec nnn=static_cast<real_prec>(mean_props[j].ntr[ddm_bin]);
  real_prec nnnv=static_cast<real_prec>(mean_props[j].ntr[ddm_bin]-1.);
  if(nnn>0)
  {
  real_prec mm=(props[j].mass[i]-mean_props[j].mass[ddm_bin]/nnn);
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
  for (int i=0;i<NX;++i){
  real_prec deltamm=dm_min+(i+0.5)*delta;
  sal<<deltamm<<"  "<<corr_props[0].mass[i]<<"   "<<corr_props[1].mass[i]<<"  "<<corr_props[2].mass[i]<<"   "<<corr_props[3].mass[i]<<"   "<<corr_props[4].mass[i]<<endl;
  cout<<deltamm<<"  "<<corr_props[0].mass[i]<<"   "<<corr_props[1].mass[i]<<"  "<<corr_props[2].mass[i]<<"   "<<corr_props[3].mass[i]<<"   "<<corr_props[4].mass[i]<<endl;
  }
  sal.close();
  outf_k=this->params._Output_directory()+"CorrelationDMSIGMA_R"+to_string(this->params._realization())+".txt";
  sal.open(outf_k.c_str());
  So.message_screen("Writting correlation in bins of dm in file ", outf_k);
  for (int i=0;i<NX;++i){
  real_prec deltamm=dm_min+(i+0.5)*delta;
  sal<<deltamm<<"  "<<corr_props[0].vel[i]<<"   "<<corr_props[1].vel[i]<<"  "<<corr_props[2].vel[i]<<"   "<<corr_props[3].vel[i]<<"   "<<corr_props[4].vel[i]<<endl;
  cout<<deltamm<<"  "<<corr_props[0].vel[i]<<"   "<<corr_props[1].vel[i]<<"  "<<corr_props[2].vel[i]<<"   "<<corr_props[3].vel[i]<<"   "<<corr_props[4].vel[i]<<endl;
  }
  sal.close();
  outf_k=this->params._Output_directory()+"CorrelationDMMAS_R"+to_string(this->params._realization())+".txt";
  sal.open(outf_k.c_str());
  So.message_screen("Writting correlation in bins of dm in file ", outf_k);
  for (int i=0;i<NX;++i){
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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::write_catalog_bin(string outputFileName)
{
    So.enter(__PRETTY_FUNCTION__);
    real_prec conversion_factor=1.0;
    int Nprop_file=MIN_N_PROP_CAT;
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
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
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    Nprop_file++;
#endif
    // The structure pof the binary file will be:
    //  ULONG Nlines
    //  int Ncolumns = 10
   //   Nlines*Ncolumns floats x,y,z,vx,vy,vz,M,vmax,Rs,s
    ofstream outStream;
#ifdef _VERBOSE_CAT_
    this->So.message_screen("Writting to binary file", outputFileName);
#endif
    outStream.open(outputFileName.c_str(), ios::binary|ios::out);
#ifdef _VERBOSE_CAT_
    this->So.message_screen("Writting Number of tracers = ", this->NOBJS);
#endif
    outStream.write(reinterpret_cast<char*>(&this->NOBJS), sizeof(ULONG));
#ifdef _VERBOSE_CAT_
    this->So.DONE();
    this->So.message_screen("Writting number of columns = ", Nprop_file);
#endif
    outStream.write(reinterpret_cast<char*>(&Nprop_file), sizeof(ULONG));
#ifdef _VERBOSE_CAT_
    this->So.DONE();
#endif
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
    data_i.push_back("Mass(Ms/h).)";
    data_i.push_back("Vmax(km/s).)";
#ifdef _USE_RS_AS_OBSERVABLE_POWER_
    data_i.push_back("Rs(kpc/h).");
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    data_i.push_back("Spin.");
#endif
#endif // end for assign_mass_post

#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
    data_i.push_back("Mass(Ms/h).)";
#endif
    data_i.push_back("Ntracers");
    for(int i=0;i<Nprop_file+1;++i)
    {
        int sizes_i=data_i[i].size();
        outStream.write((char *)&sizes_i, sizeof(sizes_i));
        outStream.write(data_i[i].c_str(), sizes_i);
    }
#endif
    for(ULONG i = 0; i< this->NOBJS; ++i)
     {
#ifdef _WRITE_COORDINATES_
       float x=static_cast<float>(Halo[i].coord1);
       float y=static_cast<float>(Halo[i].coord2);
       float z=static_cast<float>(Halo[i].coord3);
       outStream.write(reinterpret_cast<char*>(&x), sizeof(x));
       outStream.write(reinterpret_cast<char*>(&y), sizeof(y));
       outStream.write(reinterpret_cast<char*>(&z), sizeof(z));
#endif
#ifdef _WRITE_VELOCITIES_
       float vx=static_cast<float>(Halo[i].vel1*conversion_factor);
       float vy=static_cast<float>(Halo[i].vel2*conversion_factor);
       float vz=static_cast<float>(Halo[i].vel3*conversion_factor);
       outStream.write(reinterpret_cast<char*>(&vx), sizeof(vx));
       outStream.write(reinterpret_cast<char*>(&vy), sizeof(vy));
       outStream.write(reinterpret_cast<char*>(&vz), sizeof(vz));
#endif

#if defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ || defined (_ASSIGN_MASS_POST_)
       float mass=static_cast<float>(Halo[i].mass);
       outStream.write(reinterpret_cast<char*>(&mass), sizeof(mass));
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
       float vmax=static_cast<float>(Halo[i].vmax);
       outStream.write(reinterpret_cast<char*>(&vmax), sizeof(vmax));
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
        float rs=static_cast<float>(Halo[i].rs);
        outStream.write(reinterpret_cast<char*>(&rs), sizeof(rs));
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
       float spin=static_cast<float>(Halo[i].spin);
       outStream.write(reinterpret_cast<char*>(&spin), sizeof(spin));
#endif
    }
    outStream.close();
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::write_catalog(string outputFileName)
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
    //  ULONG Nlines
    //  int Ncolumns = 10
    //   Nlines*Ncolumns floats x,y,z,vx,vy,vz,M,vmax,Rs,s
    ofstream outStream;
    this->So.message_screen("Writting to binary file", outputFileName);
    outStream.open(outputFileName.c_str(), ios::binary|ios::out);
    this->So.message_screen("Writting Number of tracers = ", this->NOBJS);
    outStream.write(reinterpret_cast<char*>(&this->NOBJS), sizeof(ULONG));
    this->So.DONE();
    this->So.message_screen("Writting number of columns = ", Nprop_file);
    outStream.write(reinterpret_cast<char*>(&Nprop_file), sizeof(ULONG));
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
    for(int i=0;i<Nprop_file+1;++i)
      {
    int sizes_i=data_i[i].size();
        outStream.write((char *)&sizes_i, sizeof(sizes_i));
        outStream.write(data_i[i].c_str(), sizes_i);
      }
#endif
    for(ULONG i = 0; i< this->NOBJS; ++i)
      {
#ifdef _WRITE_COORDINATES_
    float x=static_cast<float>(Halo[i].coord1);
    float y=static_cast<float>(Halo[i].coord2);
    float z=static_cast<float>(Halo[i].coord3);
    outStream.write(reinterpret_cast<char*>(&x), sizeof(x));
    outStream.write(reinterpret_cast<char*>(&y), sizeof(y));
    outStream.write(reinterpret_cast<char*>(&z), sizeof(z));
#endif
#ifdef _WRITE_VELOCITIES_
    float vx=static_cast<float>(Halo[i].vel1*conversion_factor);
    float vy=static_cast<float>(Halo[i].vel2*conversion_factor);
    float vz=static_cast<float>(Halo[i].vel3*conversion_factor);
    outStream.write(reinterpret_cast<char*>(&vx), sizeof(vx));
    outStream.write(reinterpret_cast<char*>(&vy), sizeof(vy));
    outStream.write(reinterpret_cast<char*>(&vz), sizeof(vz));
#endif
#if defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ || defined (_ASSIGN_MASS_POST_)
    float mass=static_cast<float>(Halo[i].mass);
    outStream.write(reinterpret_cast<char*>(&mass), sizeof(mass));
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    float vmax=static_cast<float>(Halo[i].vmax);
    outStream.write(reinterpret_cast<char*>(&vmax), sizeof(vmax));
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
    float rs=static_cast<float>(Halo[i].rs);
        outStream.write(reinterpret_cast<char*>(&rs), sizeof(rs));
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    float spin=static_cast<float>(Halo[i].spin);
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
    this->So.message_screen("with ", this->NOBJS, "objects");
#endif
    for(ULONG i = 0; i< this->NOBJS; ++i)
    {
        if(Halo[i].observed==1)
#if defined _ASSIGN_PROPERTIES_
        outStream<<Halo[i].coord1<<"\t"<<Halo[i].coord2<<"\t"<<Halo[i].coord3<<"\t"<<Halo[i].vel1*conversion_factor<<"\t"<<Halo[i].vel2*conversion_factor<<"\t"<<Halo[i].vel3*conversion_factor<<"\t"<<Halo[i].mass<<"\t"<<Halo[i].vmax<<"\t"<<Halo[i].rs<<"\t"<<Halo[i].spin<<"\t"<<Halo[i].identity<<"\t"<<Halo[i].gal_cwt<<"\t"<<Halo[i].local_dm<<endl;
#elif defined _WRITE_COORDINATES_ || defined _WRITE_VELOCITIES_
//        outStream<<Halo[i].coord1<<"\t"<<Halo[i].coord2<<"\t"<<Halo[i].coord3<<"\t"<<Halo[i].vel1*conversion_factor<<"\t"<<Halo[i].vel2*conversion_factor<<"\t"<<Halo[i].vel3*conversion_factor<<endl;
        outStream<<Halo[i].redshift<<"\t"<<Halo[i].mass<<"\t"<<Halo[i].color<<"\t"<<Halo[i].stellar_mass<<"\t"<<Halo[i].relative_bias<<endl;
#endif
    }
    outStream.close();
#ifdef _VERBOSE_CAT_
    So.DONE();
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::read_catalog_bin_tng()
{
#ifdef _VERBOSE_CATALOG_
  this->So.enter(__PRETTY_FUNCTION__);
#endif
    this->xgal.resize(this->params._N_lines_binary(),0);
    this->ygal.resize(this->params._N_lines_binary(),0);
    this->zgal.resize(this->params._N_lines_binary(),0);
    this->File.read_array(this->params._file_bin_x_coord(),&this->xgal[0],this->params._N_lines_binary());
    this->File.read_array(this->params._file_bin_y_coord(),&this->ygal[0],this->params._N_lines_binary());
    this->File.read_array(this->params._file_bin_z_coord(),&this->zgal[0],this->params._N_lines_binary());
  this->NOBJS=zgal.size();
  this->Halo.resize(this->NOBJS);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->NOBJS; ++i)
 {
     Halo[i].coord1=xgal[i];
     Halo[i].coord2=ygal[i];
     Halo[i].coord3=zgal[i];
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef _USE_VELOCITIES_
void Catalog::read_catalog_bin(ULONG n, string filex, string filey, string filez,string filevx, string filevy, string filevz)
#else
void Catalog::read_catalog_bin()
#endif
{
#ifdef _VERBOSE_CATALOG_
  this->So.enter(__PRETTY_FUNCTION__);
#endif
  this->xgal.clear();
  this->ygal.clear();
  this->zgal.clear();
  this->Pxgal.clear();
  this->Pygal.clear();
  this->Pzgal.clear();
  this->Mass.clear();
  this->property.clear();
//  string field=this->Output_directory;
  string Ofile=this->params._Output_directory()+this->params._Name_survey()+"_L"+to_string(static_cast<int>(this->params._Lbox()))+"_Nft"+to_string(this->box.Nft)+"_IC"+to_string(this->params._IC_index())+"_";
#ifndef _USE_VELOCITIES_
  ULONG n=this->params._N_lines_binary();
  this->So.message_warning("Why is this written here? we should add an if");
//  this->box.NGRID=n; //
  this->box.Nft=this->params._Nft();
  this->box.NGRID=this->params._NGRID();
#endif
  string filex=this->params._file_bin_x_coord();
  string filey=this->params._file_bin_y_coord();
  string filez=this->params._file_bin_z_coord();
#ifdef _VERBOSE_CAT_
#ifdef _USE_VELOCITIES_
  cout<<endl;
  So.message_screen("Reading binary files for pos and vel of DM, and interpolating into a grid.");
#else
  cout<<endl;
  So.message_screen("Reading binary files for pos of DM, and interpolating into a grid.");
#endif
#endif
#ifdef _GET_VEL_FIELD_
  vector<real_prec>vx_field(this->box.NGRID,0);
  vector<real_prec>vy_field(this->box.NGRID,0);
  vector<real_prec>vz_field(this->box.NGRID,0);
#endif
#ifdef _GET_NGP_DENS_FIELD_
  vector<real_prec>dens_field_ngp(this->box.NGRID,0);
#endif
#ifdef _GET_CIC_DENS_FIELD_
  vector<real_prec>dens_field_cic(this->box.NGRID,0);
#endif
#ifdef _GET_TSC_DENS_FIELD_
  vector<real_prec>dens_field_tsc(this->box.NGRID,0);
#endif
#ifdef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
  string bin_file;
  int N_bin_files=32; // this number comes from the number of hdf5 files converted to binary
  vector<ULONG>N_parts(N_bin_files,0);
  string path_hdfc="/net/deimos/scratch1/balaguera/data/Numerics/Minerva_DM/z1/001/";
  ifstream np; np.open(path_hdfc+"number_of_particles.txt");
  for(int ik=0; ik<N_bin_files;++ik)
    np>>N_parts[ik];
  np.close();
  for(int ik=0; ik<N_bin_files;++ik)
    {
      bin_file=path_hdfc+"minerva_dm_x_file"+to_string(ik)+".dat";
      ULONG Nn=N_parts[ik];
      So.message_screen("N particles =", Nn);
      this->xgal.clear();
      this->xgal.shrink_to_fit();
      this->xgal.resize(Nn,0);
      this->File.read_array(bin_file, this->xgal);			      /*   this might not be needed unless we use it for seome other calculations. If not, simply write it to output file
                                   this->power_in_bins[index_bins].modes=this->modes_g;
                                   this->power_in_bins[index_bins].kvector=this->kvector_data;
                                   this->power_in_bins[index_bins].pk0=this->pk0;
                                   this->power_in_bins[index_bins].pk2=this->pk2;
                                   this->power_in_bins[index_bins].pk4=this->pk4;
                              */
      bin_file=path_hdfc+"minerva_dm_y_file"+to_string(ik)+".dat";
      this->ygal.clear();
      this->ygal.shrink_to_fit();
      this->ygal.resize(Nn,0);
      this->File.read_array(bin_file, this->ygal);
      bin_file=path_hdfc+"minerva_dm_z_file"+to_string(ik)+".dat";
      this->zgal.clear();N_lines
      this->zgal.shrink_to_fit();
      this->zgal.resize(Nn,0);
      this->File.read_array(bin_file, this->zgal);
#ifdef _GET_VEL_FIELD_
      bin_file=path_hdfc+"minerva_dm_velx_file"+to_string(ik)+".dat";
      this->Pxgal.clear();
      this->Pxgal.shrink_to_fit();
      this->Pxgal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pxgal);
      bin_file=path_hdfc+"minerva_dm_vely_file"+to_string(ik)+".dat";
      this->Pygal.clear();
      this->Pygal.shrink_to_fit();
      this->Pygal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pygal);
      bin_file=path_hdfc+"minerva_dm_velz_file"+to_string(ik)+".dat";
      this->Pzgal.clear();
      this->Pzgal.shrink_to_fit();
      this->Pzgal.resize(Nn,0);
      this->File.read_array(bin_file, this->Pzgal);
#endif
#ifdef _GET_NGP_DENS_FIELD_
      So.message_screen("Interpolating density to grid ngp");
      getDensity_NGP(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_ngp,false);
      So.DONE();
#endif
#ifdef _GET_CIC_DENS_FIELD_
      So.message_screen("Interpolating density to grid CIC");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_cic,false);
      So.DONE();
#endif
#ifdef _GET_TSC_DENS_FIELD_
      So.message_screen("Interpolating density to grid TSC");
      getDensity_TSC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_tsc,false);
      So.DONE();
#endif
#ifdef _GET_VEL_FIELD_
      So.message_screen("Interpolating Vx to grid");N_lines
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pxgal,vx_field,true);
      So.DONE();
      So.message_screen("Interpolating Vy to grid");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pygal,vy_field,true);
      So.DONE();
      So.message_screen("Interpolating Vz to grid");
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->Pzgal,vz_field,true);
      So.DONE();
#endif
    }
#ifdef _GET_VEL_FIELD_
  ULONG ec=0;
#pragma omp parallel for reduction(+:ec)
  for(ULONG i=0;i<this->box.NGRID;++i)
    {
      if(dens_field[i]!=0)
        {
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          vx_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          vy_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          vz_field[i]/=static_cast<real_prec>(dens_field_cic[i]);
        }
    else{
      // this is not extrictly speaking right. We would have to assingn the velocity of the closest dm particle if the cell is empty.
      // However, for practical reasons, we will be using the velocity field in conjunction with the density (e.g., kinetic energy)
       ec++;
         vx_field[i]=0;
         vy_field[i]=0;
         vz_field[i]=0;
       }
  }
  So.message_screen("Number of empty cell=", ec);
#endif
#ifdef _GET_NGP_DENS_FIELD_
  File.write_array(Ofile+"density_NGP_Real1_z1_N"+to_string(this->box.Nft),dens_field_ngp);
  dens_field_ngp.clear();
  dens_field_ngp.shrink_to_fit();
#endif
#ifdef _GET_CIC_DENS_FIELD_
  File.write_array(Ofile+"density_CICp_Real1_z1_N"+to_string(this->box.Nft),dens_field_cic);
  dens_field_cic.clear();
  dens_field_cic.shrink_to_fit();
#endif
#ifdef _GET_TSC_DENS_FIELD_
  File.write_array(Ofile+"density_TSC_Real1_z1_N"+to_string(this->box.Nft),dens_field_tsc);
  dens_field_tsc.clear();
  dens_field_tsc.shrink_to_fit();
#endif
#ifdef _GET_VEL_FIELD_
  File.write_array(field+"Minerva_DM_vx_CIC_Real1_z1_N"+to_string(this->box.Nft),vx_field);
  File.write_array(field+"Minerva_DM_vy_CIC_Real1_z1_N"+to_string(this->box.Nft),vy_field);
  File.write_array(field+"Minerva_DM_vz_CIC_Real1_z1_N"+to_string(this->box.Nft),vz_field);
  vx_field.clear();N_lines
  vy_field.clear();
  vz_field.clear();
#endif
#else // if we only read from a single file
  vector<real_prec>dummy;
  if(true==this->params._dilute_dm_sample())
  {
      gsl_rng_env_setup();
      const gsl_rng_type *Tn;
      gsl_rng_default_seed=1015;
      Tn =  gsl_rng_mt19937 ;//gsl_rng_default;
      gsl_rng *rn = gsl_rng_alloc (Tn);
      dummy.resize(n,0);
      this->File.read_array(filex, dummy);
      ULONG Ndil=static_cast<ULONG>(n*this->params._fraction_dilute());
      vector<ULONG>aux_p;
      for(ULONG i=0;i<n;++i)
        aux_p.push_back(i);
      gsl_ran_shuffle(rn,&aux_p[0], aux_p.size() ,sizeof(ULONG));
      for(ULONG i=0;i<Ndil;++i)
        this->xgal.push_back(static_cast<double>(dummy[aux_p[i]]));
      dummy.resize(n,0);
      this->File.read_array(filey, dummy);
      for(ULONG i=0;i<Ndil;++i)
       this->ygal.push_back(static_cast<double>(dummy[aux_p[i]]));
      dummy.resize(n,0);
      this->File.read_array(filez, dummy);
      for(ULONG i=0;i<Ndil;++i)
       this->zgal.push_back(static_cast<double>(dummy[aux_p[i]]));
  }
else{
  this->xgal.resize(n,0);
  this->File.read_array(filex, this->xgal);
  this->ygal.resize(n,0);
  this->File.read_array(filey, this->ygal);
    this->zgal.resize(n,0);
    this->File.read_array(filez, this->zgal);
    this->params.set_N_lines_binary(this->zgal.size());
    }
 if(false==this->params._use_low_pass_filter())
 {
#if defined _GET_NGP_DENS_FIELD_
      getDensity_NGP(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_ngp,false);
      File.write_array(Ofile+"NGP",dens_field_ngp);
      dens_field_ngp.clear();
      dens_field_ngp.shrink_to_fit();
      So.DONE();
#endif
#if defined _GET_CIC_DENS_FIELD_
      getDensity_CIC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_cic,false);
      So.DONE();
      File.write_array(Ofile+"CIC",dens_field_cic);
      dens_field_cic.clear();
      dens_field_cic.shrink_to_fit();
#endif
#if defined _GET_TSC_DENS_FIELD_
      getDensity_TSC(this->box.Nft, this->box.Nft, this->box.Nft,this->box.Lbox, this->box.Lbox,this->box.Lbox, this->box.d1, this->box.d2, this->box.d3,0,0,0,this->xgal,this->ygal,this->zgal,this->zgal, dens_field_tsc,false);
      So.DONE();
      File.write_array(Ofile+"TSC",dens_field_tsc);
      dens_field_tsc.clear();
      dens_field_tsc.shrink_to_fit();
#endif
  }
#ifdef _USE_VELOCITIES_
  dummy.resize(n,0);
  this->File.read_array(filevx, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pxgal.push_back(static_cast<double>(dummy[i]));

  this->File.read_array(filevy, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pygal.push_back(static_cast<double>(dummy[i]));

  this->File.read_array(filevz, dummy);
  for(ULONG i=0;i<n;++i)
    this->Pzgal.push_back(static_cast<double>(dummy[i]));
#endif
  dummy.clear();
  dummy.shrink_to_fit();
#endif
 if(true==this->params._use_low_pass_filter())
  {
        string  file_hr=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft_HR())+"_MAS"+to_string(this->box.masskernel)+"_"+this->params._Name_survey();
        string  file=this->params._Output_directory()+this->params._Name_survey()+"_Ncounts_Nft"+to_string(this->params._Nft())+"_MAS"+to_string(this->box.masskernel)+"_"+this->params._Name_survey();
        vector<real_prec>field_hr(this->params._NGRID_HR() ,0);
        getDensity_CIC(this->params._Nft_HR(),this->params._Nft_HR(),this->params._Nft_HR(),this->params._Lbox(),this->params._Lbox(),this->params._Lbox(), this->params._d1_HR(), this->params._d2_HR(), this->params._d3_HR(),0,0,0,this->xgal,this->ygal,this->zgal,this->zgal,field_hr,false);
        this->File.write_array(file_hr, field_hr);
        vector<real_prec>field(this->params._NGRID() ,0);
        So.message_screen("Applying low pass filter");
        low_pass_filter(this->params._Nft_HR(), this->params._Nft(),this->params._masskernel(),false,field_hr, field,this->params._Lbox());
        field_hr.clear(); field_hr.shrink_to_fit();
        //get_overdens(field,field);
        this->File.write_array(file, field);
        field.clear(); field.shrink_to_fit();
    }
     this->xgal.clear(); this->xgal.shrink_to_fit();
     this->ygal.clear(); this->ygal.shrink_to_fit();
     this->zgal.clear(); this->zgal.shrink_to_fit();

  this->NOBJS=n;
  }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// In this function a catalog of DM or DM tracers, in ascii format, is read and their proeprtiess read, according to the request of the parameter file.
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
void Catalog::read_catalog(string input_file, real_prec prop_min)
#elif defined (_USE_MASS_BINS_PK_)
void Catalog::read_catalog(string input_file, real_prec prop_min, real_prec prop_max)
#endif
{
   this->So.enter(__PRETTY_FUNCTION__);
  // The input parameter min_mass cann be also VMAX min, according to preproc definitions
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
  real_prec prop_max=1e20;
#endif
  int i_x=-1;
  int i_y=-1;
  int i_z=-1;
  int i_vx=-1;
  int i_vy=-1;
  int i_vz=-1;
  int i_mass=-1; int i_weight=-1; int i_mean_density=-1; int i_sf=-1; int i_vmax=-1; int i_rs=-1;int i_rvir=-1; int i_virial=-1; int i_spin=-1;int i_spin_bullock=-1;int i_vrms=-1;int i_b_to_a=-1;int i_c_to_a=-1;
  int i_redshift=-1;
  int i_stellar_mass=-1;
  int i_color=-1;
  int i_abs_mag=-1;
  int i_app_mag=-1;
  if(this->type_of_object=="DM")
    {
      i_x= this->params._i_coord1_dm();
      i_y= this->params._i_coord2_dm();
      i_z= this->params._i_coord3_dm();
      i_vx= this->params._i_v1_dm();
      i_vy= this->params._i_v2_dm();
      i_vz= this->params._i_v3_dm();
      i_mass = this->params._i_mass_dm();
    }
  else if(this->type_of_object=="TRACER" || this->type_of_object=="TRACER_REF" || this->type_of_object=="TRACER_MOCK" || this->type_of_object=="TRACER_MOCK_ONLY_COORDS"|| this->type_of_object=="TRACER_REF_ONLY_COORDS")
    {
      i_x= this->params._i_coord1_g();
      i_y= this->params._i_coord2_g();
      i_z= this->params._i_coord3_g();
      i_vx= this->params._i_v1_g();
      i_vy= this->params._i_v2_g();
      i_vz= this->params._i_v3_g();
      i_mass = this->params._i_mass_g();
      i_vmax = this->params._i_vmax_g();
      i_vrms = this->params._i_vrms_g();
      i_weight = this->params._i_weight1_g();
      i_mean_density= this->params._i_mean_density_g();
      i_sf = this->params._i_sf_g();
      i_rs = this->params._i_rs_g();
      i_rvir = this->params._i_rvir_g();
      i_virial = this->params._i_virial_g();
      i_spin = this->params._i_spin_g();
      i_spin_bullock = this->params._i_spin_bullock_g();
      i_b_to_a = this->params._i_b_to_a_g();
      i_c_to_a = this->params._i_c_to_a_g();
      i_stellar_mass = this->params._i_stellar_mass_g();
      i_color = this->params._i_color_g();
      i_app_mag = this->params._i_app_mag_g();
      i_abs_mag = this->params._i_abs_mag_g();
      if(this->params._sys_of_coord_g()==2)
        i_redshift= this->params._i_coord3_g();// note that this applies only in the case in which coordsa are in pseudo*-equatorial
    }
  else if(this->type_of_object=="RANDOM")
    {
      i_x= this->params._i_coord1_r();
      i_y= this->params._i_coord2_r();
      i_z= this->params._i_coord3_r();
      i_weight = this->params._i_weight1_r();
      i_mean_density= this->params._i_mean_density_r();
      i_mass= this->params._i_mass_r();
      i_stellar_mass = this->params._i_stellar_mass_r();
      i_color = this->params._i_color_r();
      i_app_mag = this->params._i_app_mag_r();
      i_abs_mag = this->params._i_abs_mag_r();
      if(this->params._sys_of_coord_r()>=2)
        i_redshift= this->params._i_coord3_r();// note that this applies only in the case in which coords are in pseudo-equatorial
    }
  //We need to differentiate between observabloe and settings:
  //"Observable" will be using in case we want to do bins or cuts on a given property within a sample
  // -------------------------------------------------------------------------
  int i_observable=i_mass;
  int i_setting=i_mass;
  real_prec units_observable=1;
  real_prec units_settings=1;
  real_prec min_cut=0;
#ifdef _POWER_
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_
  i_observable=i_mass;
  i_setting=i_mass;
  units_observable=this->params._MASS_units();
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_
  i_observable=i_vmax;
  i_setting=i_mass;
  min_cut=pow(10, this->params._LOGMASSmin());
#elif defined _USE_RS_AS_OBSERVABLE_POWER_
  i_observable=i_rs;
  i_setting=i_mass;
#elif defined _USE_SPIN_AS_OBSERVABLE_POWER_
  i_observable=i_spin;
  i_setting=i_mass;
#endif
#else
#ifdef _USE_MASS_AS_PRIMARY_OBSERVABLE_
  i_observable=i_mass;
  units_observable=this->params._MASS_units();
#elif defined _USE_VMAX_AS_PRIMARY_OBSERVABLE_
  i_observable=i_vmax;
  units_observable=1;
#endif
  i_setting=i_mass;
  units_settings=this->params._MASS_units();
  min_cut=pow(10, this->params._LOGMASSmin());// Value of the settng-property defining the minimum setting property. Set mass by default.
#endif
#if defined (_USE_VMAX_AS_PRIMARY_OBSERVABLE_) || defined (_USE_VMAX_AS_PRIMARY_OBSERVABLE_POWER_)
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS" ||  this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )
    if(i_vmax<0)
        So.message_warning("No information for Vmax. Check .ini parameter");
#endif
#if defined (_USE_MASS_AS_PRIMARY_OBSERVABLE_) || defined (_USE_MASS_AS_PRIMARY_OBSERVABLE_POWER_)
  if(this->type_of_object!="TRACER_MOCK_ONLY_COORDS" ||this->type_of_object!="RANDOM" )
    if(i_mass<0)
        So.message_warning("No information for Mass. Check .ini parameter");
#endif
  vector<real_prec> prop;
  ULONG NLINES =0;
     // This is meant for ascie reading.
#ifdef _READ_BINARY_BAM_FORMAT_
   NLINES= this->File.read_binary_file(input_file, prop);
#else
  NLINES = this->File.read_file(input_file, prop,NTHREADS);
#endif
  this->NCOLS=(static_cast<ULONG>(prop.size()/NLINES));
#ifdef _ABACUS_
//* This section converts the fiducial input for abacus halo masses (given as number of particles) to Ms/h.  *//
  if(i_mass>=0)
  {
   So.message_screen("Converting number of particles to halo masses in ABACUS");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NLINES;++i)
        prop[i_mass+NCOLS*i] *= PART_TO_MASSES_ABACUS;
  }
#endif
  // We now proceed to allocate read and allcoate the different properties
  // taking into account the different limits on properties.
 // This is done individually, since putting Nob ifs inside a M loop costs more than a single "if" outside each M loop
  // The minimum mass (M or VMAX) defines the number of used tracers!!!!!!!!!!!!!
  // First stage: counting number of objects above limits. The corresponding loop is paralleized and tested
#ifdef _SET_CAT_WITH_CUT_
#ifdef _VERBOSE_CAT_
  if(i_mass>0)
    if(this->type_of_object!="RANDOM")
        So.message_screen("Catalog will be selected with mass cut at ", pow(10,this->params._LOGMASSmin()));
#endif
#endif

#ifdef _VERBOSE_CAT_
     So.message_screen("Counting number of", this->type_of_object);
#endif
 // **************************************************************************************************************
   ULONG count_new_nobj=0;
   if(this->type_of_object!="RANDOM")
    {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count_new_nobj)
#endif
          for(ULONG i=0;i<NLINES;++i)
            {
             real_prec obser=LARGE_NUMBER;
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
          if(this->type_of_object=="TRACER"  || this->type_of_object== "TRACER_REF" || this->type_of_object=="TRACER_MOCK" )
            {
#ifdef _POWER_
                  obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
                if(i_setting>0)//extra security ckeck
                   obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
           }

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
                if(prop[i_redshift+i*this->NCOLS] < this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
                    count_new_nobj++;

#ifdef _POWER_
#ifdef _SET_CAT_WITH_MASS_CUT_
        }
#endif
#endif
        }
#ifdef _POWER_
      So.message_screen("Minimum value of observable requested =", prop_min);
#ifdef _USE_REDSHIFT_BINS_
      So.message_screen("Minimum value of redshift  =", this->params._redshift_min_sample());
#endif
#ifdef _USE_MASS_BINS_PK_
      So.message_screen("Maximum value of observable requested =", prop_max);
#ifdef _USE_REDSHIFT_BINS_
      So.message_screen("Maximum value of redshift  =", this->params._redshift_max_sample());
#endif
      So.message_screen("Number of tracers in the prop-bin =", count_new_nobj);
#else
#ifdef _SET_CAT_WITH_MASS_CUT_
      So.message_screen("Minimum mass requested =", prop_min, "Ms/h");
#ifdef _USE_MASS_BINS_PK_
      So.message_screen("Maximum mass requested =", prop_max, "Ms/h");
      So.message_screen("Number of in the mass bin =", count_new_nobj);
#else
      So.message_screen("Number of tracers above mass limit =", count_new_nobj);
#endif
#elif defined (_SET_CAT_WITH_VMAX_CUT_)
      So.message_screen("Minimum VMAX requested =", prop_min, "Km/s");
#ifdef _USE_MASS_BINS_PK_
      So.message_screen("Maximum VMAX requested =", prop_max, "Km/s");
      So.message_screen("Number of in the mass bin =", count_new_nobj);
#else
      So.message_screen("Number of tracers above VMAX limit =", count_new_nobj);
#endif
#endif
#endif
#endif // end of def power
    }
    else if (this->type_of_object=="RANDOM")// ojo que aca asumimos que los randoms no tinen propiuedades para ser seleccionados
       count_new_nobj=NLINES;
    So.message_screen("Found ",count_new_nobj," objects");
    So.DONE();
   this->NOBJS=count_new_nobj;
   this->Halo.resize(count_new_nobj);
   if (this->type_of_object=="RANDOM")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NLINES;++i)  // This loop should not be paralelized, it generates problems.
       {
         this->Halo[i].coord1=prop[i_x+i*NCOLS];
         this->Halo[i].coord2=prop[i_y+i*NCOLS];
         this->Halo[i].coord3=prop[i_z+i*NCOLS];
         if(this->params._sys_of_coord_r()==I_EQZ)
             this->Halo[i].redshift=prop[i_z+i*NCOLS];
      }
    if(i_mean_density>0)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<NLINES;++i)  // This loop should not be paralelized, it generates problems.
            this->Halo[i].mean_density=prop[i_mean_density+i*NCOLS];
    }
  ULONG iN=0;
    // Now alocate the ID of the grid for each member This loop has problems with being parallelized.
#ifdef _VERBOSE_CAT_
#ifndef _POWER_   // This lines are only useful when analysing the catalog with BAM, not to measure the power spectrum
  if(this->params._sys_of_coord_g()==0)
      So.message_screen("Getting grid-ID from tracer coordinates");
#endif
#endif
 if (this->type_of_object!="RANDOM")
    {
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN)
#endif
  for(ULONG i=0;i<NLINES;++i)  // This loop should not be paralelized, it generates problems.
    {
      real_prec obser=LARGE_NUMBER;
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
      if(prop[i_setting+i*this->NCOLS]>min_cut)
    {
#endif
#endif
          if(this->type_of_object=="TRACER"  || this->type_of_object=="TRACER_REF" || this->type_of_object=="TRACER_MOCK"  )
            {
              //This if is inside, for we can lack of mass but still will to read coordinates
#ifdef _POWER_
             obser = (i_mass> 0 || i_vmax>0) ? prop[i_observable+i*this->NCOLS]*units_settings  : 0 ;  // the questio  must be generalized to other props
#else
              if(i_setting>0)//extra security ckeck
                 obser = prop[i_setting+i*this->NCOLS]*units_settings;
#endif
           }
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
           if(obser>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser>= prop_min && obser< prop_max)
#endif
          {
              real_prec x=prop[i_x+i*NCOLS];
              real_prec y=prop[i_y+i*NCOLS];
              real_prec z=prop[i_z+i*NCOLS];
              this->Halo[iN].coord1=x;
              this->Halo[iN].coord2=y;
              this->Halo[iN].coord3=z;
              if(this->params._sys_of_coord_r()==I_EQZ)
                  this->Halo[i].redshift=z;
#ifdef _USE_REDSHIFT_BINS_
                if(z< this->params._redshift_max_sample() && z>= this->params._redshift_min_sample())
                    {
                    this->Halo[iN].observed=true;
#endif
#ifndef _POWER_   // This lines are only useful when analysing the catalog with BAM, not to measure the power spectrum
                if(this->params._sys_of_coord_g()==0)
                {
                    this->Halo[iN].GridID = grid_ID(&box, x,y,z);
                    this->Halo[iN].GridID_n = grid_ID(&box_n, x,y,z);
                }
                this->Halo[iN].observed=true;
#endif
#ifdef _USE_NUMBER_OF_SATELLITES_
                this->Halo[iN].n_satellites=prop[i_sf+i*NCOLS];
#endif
                     iN++;
#ifdef _USE_REDSHIFT_BINS_
                }
#endif
                    }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
    }
#endif
#endif
    }
#ifdef _VERBOSE_CAT_a
   So.message_screen("Min x-coord", get_min("_XCOORD_"));
   So.message_screen("Min y-coord", get_min("_YCOORD_"));
   So.message_screen("Min z-coord", get_min("_ZCOORD_"));
   So.message_screen("Max x-coord", get_max("_XCOORD_"));
   So.message_screen("Max y-coord", get_max("_YCOORD_"));
   So.message_screen("Max z-coord", get_max("_ZCOORD_"));
#endif
  if(count_new_nobj!=iN)
    {
      So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
      So.message_warning("Line", __LINE__);
    }
  So.DONE();
}
 // ALLOCATE VMAX
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_vmax>0 && i_vmax<NCOLS) // This if is outside the loop, for we here want vmax and it applies only when this "if" is satisfied
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating Vmax from tracer");
#endif
      real_prec mean_m=0;
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
#ifdef _POWER_
             real_prec prop_sel =prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec prop_sel =prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop_sel>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(prop_sel>= prop_min && prop_sel < prop_max)
#endif
             #ifdef _USE_REDSHIFT_BINS_
                 if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
                 {
                   real_prec proper=prop[i_vmax+i*NCOLS]*units_observable;
                   this->Halo[iN].vmax=proper;
                   mean_m+=proper;
                   iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        }
#endif
#endif
      }
#ifdef _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
      this->min_vmax=get_min("_VMAX_");
      this->max_vmax=get_max("_VMAX_");
      if(true==this->params._set_bins_equal_number_tracers_main_property())// not often done, as it will be rather fixedto check diferent redshifts
      {
       int nbins_prop=1;
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       nbins_prop=this->params._Number_of_bins_equal_number_tracers_main_property();
       this->Number_of_tracers_in_vmax_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
       this->get_intervals_equal_number("_VMAX_",min_aux,max_aux);

       for(int i=0;i<=nbins_prop;++i)
         this->params.set_VMAXbins_max(i,log10(max_aux[i]));
       for(int i=0;i<=nbins_prop;++i)
          this->params.set_VMAXbins_min(i,log10(min_aux[i]));
       this->So.message_screen("Number of bins in this property =", this->params._Number_of_bins_equal_number_tracers_main_property());
       this->So.message_screen("\tMin log Vmax  =", log10(this->min_vmax));
       this->So.message_screen("\tMax log Vmax  =", log10(this->max_vmax));
      }
      else{
          this->So.message_screen("\tMin log Vmax  =", log10(this->min_vmax));
          this->So.message_screen("\tMax log Vmax  =", log10(this->max_vmax));
          this->Number_of_tracers_in_vmax_bins.resize(this->params._NVMAXbins_power(),0);
          So.message_screen("\tVmax-bins defined fixed width (parameter file)");
      }
#else
      if(true==this->params._set_bins_equal_number_tracers())
      {
       this->min_vmax=get_min("_VMAX_");
       this->max_vmax=get_max("_VMAX_");
       int nbins_prop=1;
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       nbins_prop=this->params._Number_of_bins_equal_number_tracers_main_property();
       this->get_intervals_equal_number("_VMAX_",min_aux,max_aux);
       for(int i=0;i<=nbins_prop;++i)
         this->params.set_VMAXbins_max(i,max_aux[i]);
       for(int i=0;i<=nbins_prop;++i)
          this->params.set_VMAXbins_min(i,min_aux[i]);
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
       this->So.message_screen("\tMin Vmax  =", this->min_vmax);
       this->So.message_screen("\tMax Vmax  =", this->max_vmax);
     }
#endif
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
        mean_m/=static_cast<real_prec>(this->NOBJS);
#ifdef _USE_LOG_MASS_
        this->So.message_screen("\tMean mass of tracer catalogue =", pow(10,mean_m)*this->params._MASS_units(), "km/s");
#else
        this->So.message_screen("\tMean VMAX of tracer catalogue =", mean_m*this->params._MASS_units(), "km/s");
#endif
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of vmax not available in catalog or column not properly set");
#endif
  So.DONE();
  // *********************************************************************
  // ALLOCATE VRMS:
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_vrms>0 && i_vrms<NCOLS) // This if is outside the loop, for we here want vmax and it applies only when this "if" is satisfied
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating Vrms from tracer");
#endif
      real_prec mean_m=0;
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
#ifdef _POWER_
             real_prec prop_sel =prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec prop_sel =prop[i_setting+i*this->NCOLS]*units_settings;
#endif

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop_sel>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(prop_sel>= prop_min && prop_sel < prop_max)
#endif
             #ifdef _USE_REDSHIFT_BINS_
                 if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
                 {
                   real_prec proper=prop[i_vrms+i*NCOLS];
                   this->Halo[iN].vrms=proper;
                   mean_m+=proper;
                   iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        }
#endif
#endif
      }
      if(true==this->params._set_bins_equal_number_tracers())
      {
      this->min_vrms=get_min("_VRMS_");
      this->max_vrms=get_max("_VRMS_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      this->get_intervals_equal_number("_VRMS_",min_aux,max_aux);
      for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_VRMSbins_max(i,max_aux[i]);
      for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_VRMSbins_min(i,min_aux[i]);
      this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
      this->So.message_screen("\tMin Vrms  =", this->min_vrms);
      this->So.message_screen("\tMax Vrms  =", this->max_vrms);
  }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
        mean_m/=static_cast<real_prec>(this->NOBJS);
        this->So.message_screen("\tMean Vrms of tracer catalogue =", mean_m*this->params._MASS_units(), "km/s");
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of vrms not available in catalog or column not properly set");
#endif
      So.DONE();
  // *********************************************************************
  // ALLOCATE SPIN:
if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )&&  i_spin>0 && i_spin<NCOLS)
{
#ifdef _VERBOSE_CAT_
  So.message_screen("Allocating Spin from tracers");
#endif
  iN=0;
  real_prec mean_m=0;

#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
  for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
      if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
#ifdef _POWER_
         real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
         real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
      if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
      {
            real_prec proper=prop[i_spin+i*NCOLS];
            this->Halo[iN].spin=proper;
    iN++;
            mean_m+=proper;
      }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        }
#endif
#endif
   }
  if(true==this->params._set_bins_equal_number_tracers())
  {
  this->min_spin=get_min("_SPIN_");
  this->max_spin=get_max("_SPIN_");
  vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
  vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
  this->get_intervals_equal_number("_SPIN_",min_aux,max_aux);
  for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
   {
     this->params.set_SPINbins_min(i,min_aux[i]);
     this->params.set_SPINbins_max(i,max_aux[i]);
   }
  this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
  this->So.message_screen("\tMin spin  =", this->min_spin);
  this->So.message_screen("\tMax spin  =", this->max_spin);
 }
 if(count_new_nobj!=iN)
   {
      So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
      So.message_warning("Line", __LINE__);
    }
#ifdef _VERBOSE_CAT_
     mean_m/=static_cast<real_prec>(this->NOBJS);
#ifdef _USE_LOG_MASS_
    this->So.message_screen("\tMean mass of tracer catalogue =", pow(10,mean_m)*this->params          this->So.message_screen("\tMin Vmax  =", this->min_vmax);
             this->So.message_screen("\tMax Vmax  =", this->max_vmax);
._MASS_units(), "km/s");
#else
    this->So.message_screen("\tMean spin of tracer catalogue =", mean_m*this->params._MASS_units());
#endif
#endif
}
#ifdef _VERBOSE_CAT_
else
So.message_warning("Information of spin not available in catalog or column not properly set");
So.DONE();
#endif
// ALLOCATE SPIN_BULLOCK
if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )&&  i_spin_bullock >0 && i_spin_bullock<NCOLS)
{
#ifdef _VERBOSE_CAT_
So.message_screen("Allocating Spin-Bullock from tracer catalog");
#endif
iN=0;
real_prec mean_m=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
for(ULONG i=0;i<NLINES;++i)
  {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
    if(prop[i_setting+i*this->NCOLS]>min_cut)
      {
#endif
#endif
#ifdef _POWER_
       real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
       real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
    if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
      if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
    {
    real_prec proper=prop[i_spin_bullock+i*NCOLS];
          this->Halo[iN].spin_bullock=proper;
  iN++;
          mean_m+=proper;
    }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
      }
#endif
#endif
 }
if(true==this->params._set_bins_equal_number_tracers())
{
this->min_spin=get_min("_SPIN_BULLOCK_");
this->max_spin=get_max("_SPIN_BULLOCK_");
vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
 vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
this->get_intervals_equal_number("_SPIN_BULLOCK_",min_aux,max_aux);
 for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
 {
   this->params.set_SPINBULLOCKbins_min(i,min_aux[i]);
   this->params.set_SPINBULLOCKbins_max(i,max_aux[i]);
 }
this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
this->So.message_screen("\tMin spin_b  =", this->min_spin_bullock);
this->So.message_screen("\tMax spin_b  =", this->max_spin_bullock);
}
if(count_new_nobj!=iN)
 {
    So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
    So.message_warning("Line", __LINE__);
  }
#ifdef _VERBOSE_CAT_
   mean_m/=static_cast<real_prec>(this->NOBJS);
#ifdef _USE_LOG_MASS_
  this->So.message_screen("\tMean mass of tracer catalogue =", pow(10,mean_m)*this->params          this->So.message_screen("\tMin Vmax  =", this->min_vmax);
           this->So.message_screen("\tMax Vmax  =", this->max_vmax);
._MASS_units(), "km/s");
#else
  this->So.message_screen("\tMean spin_b of tracer catalogue =", mean_m);
#endif
#endif
}
#ifdef _VERBOSE_CAT_
else
    So.message_warning("Information of spin_bullock not available in catalog or column not properly set");
So.DONE();
#endif
// ALLOCATE STUFF FOR MULTISCALING
#ifndef _POWER_
#ifdef _MULTISCALE_
#ifdef _USE_MULTISCALE_PROPERTY_ASSIGNMENT_NEW_
#ifndef _USE_FIXED_MULTISCALE_LEVELS_
if(type_of_object=="TRACER_REF")
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Computing numbers for multi-scale mass assgnment");
#endif
      vector<ULONG> n_auxtr(this->params._Number_of_MultiLevels(),0);
      for(ULONG i=0;i<NLINES;++i)
        {
         real_prec prop_sel = prop[i_setting+i*this->NCOLS]*units_settings;
         if(prop_sel>= prop_min)
          {
           real_prec obser= prop[i_observable+i*this->NCOLS]*units_observable;
           if(obser>= this->params.get_PropThreshold_MultiLevels(0))
             n_auxtr[0]++;
           for(int il=1; il<this->params._Number_of_MultiLevels();++il )
             if(obser>= this->params.get_PropThreshold_MultiLevels(il) &&  obser<this->params.get_PropThreshold_MultiLevels(il-1) )
               n_auxtr[il]++;
         }
      }
      ULONG part_level_counter=0;
      for(int il=0; il<this->params._Number_of_MultiLevels();++il )
        {
//          ULONG N_AUX=  static_cast<ULONG>(floor(this->params.get_Props_Tolerance_MultiLevels(il)*n_auxtr[il]));
          ULONG N_AUX=  n_auxtr[il];
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
#endif // END OF ifndef _power_
      // *********************************************************************
      // ALLOCATE WEIGHS
      // If information of weights is available, then
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") &&  i_weight>0 && i_weight<NCOLS)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating weighs from tracers");
#endif
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
          real_prec obser=LARGE_NUMBER;
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             obser= (i_mass> 0 || i_vmax>0) ? prop[i_observable+i*this->NCOLS]*units_settings : 0;
#else
             obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          real_prec propert= (i_mass> 0 || i_vmax>0) ?  prop[i_setting+i*this->NCOLS]*units_settings : 1.0;
           real_prec prop_min_new= (i_mass> 0 || i_vmax>0) ?  prop_min : 0.0;
          if(propert >= prop_min_new )
#elif defined (_USE_MASS_BINS_PK_)
            if( obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
          {
                    this->Halo[iN].weight1=prop[i_weight+i*NCOLS];
                    iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
          }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
    }
#ifdef _VERBOSE_CAT_
      So.DONE();
#endif
      // *********************************************************************
      // ALLOCATE MEAN DENSITY
  if(this->type_of_object!="RANDOM")
   if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") &&  i_mean_density>0 && i_mean_density<NCOLS )
     {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating number density from ", this->type_of_object);
#endif
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
          real_prec obser=LARGE_NUMBER;

#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             obser= (i_mass> 0 || i_vmax>0) ? prop[i_observable+i*this->NCOLS]*units_settings : 0;
#else
             obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
             real_prec propert= (i_mass> 0 || i_vmax>0) ?  prop[i_setting+i*this->NCOLS]*units_settings : 1.0;
              real_prec prop_min_new= (i_mass> 0 || i_vmax>0) ?  prop_min : 0.0;
             if(propert >= prop_min_new )
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
          {
        this->Halo[iN].mean_density=prop[i_mean_density+i*NCOLS];
        iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
              }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of mean_density not available in catalog or column not properly set");
#endif
      So.DONE();

// *********************************************************************
// ALLOCATE RS
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )&& i_rs>0 && i_rs<NCOLS )
    {
#ifdef _VERBOSE_CAT_
          So.message_screen("Allocating Rs from tracer catalog");
#endif
          real_prec mean_m=0;
          iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
             {
                real_prec proper= prop[i_rs+i*NCOLS]>0 ?   prop[i_rs+i*NCOLS]*units_observable :0.1 ;  //parche, porque luego tomo logaritmos y revienta
                this->Halo[iN].rs=proper;
                 mean_m+=proper;
                 iN++;
              }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
        }
      if(true==this->params._set_bins_equal_number_tracers())
      {
      this->min_rs=get_min("_RS_");
      this->max_rs=get_max("_RS_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       this->get_intervals_equal_number("_RS_",min_aux,max_aux);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_RSbins_max(i,max_aux[i]);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_RSbins_min(i,min_aux[i]);
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
       this->So.message_screen("\tMin Rs  =", this->min_rs);
       this->So.message_screen("\tMax Rs  =", this->max_rs);
        }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
#ifdef _USE_LOG_MASS_
        this->So.message_screen("\tMean mass of tracer catalogue =", pow(10,mean_m)*this->params._MASS_units(), "km/s");
#else
        this->So.message_screen("\tMean RS of tracer catalogue =", mean_m*this->params._MASS_units(), "kpc/h");
#endif
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of RS not available in catalog or column not properñy set");
      So.DONE();
#endif
// *********************************************************************
// ALLOCATE Rvir
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )&& (i_rvir>0 && i_rvir<NCOLS) )
    {
#ifdef _VERBOSE_CAT_
          So.message_screen("Allocating Rvir from tracers");
#endif
          real_prec mean_m=0;
          iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
             {
                real_prec proper=prop[i_rvir+i*NCOLS]*units_observable;
                this->Halo[iN].rvir=proper;
                 mean_m+=proper;
                 iN++;
              }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
        }
      if(true==this->params._set_bins_equal_number_tracers())
      {
      this->min_rs=get_min("_RVIR_");
      this->max_rs=get_max("_RVIR_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       this->get_intervals_equal_number("_RVIR_",min_aux,max_aux);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_RVIRbins_max(i,max_aux[i]);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_RVIRbins_min(i,min_aux[i]);
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
       this->So.message_screen("\tMin Rvir  =", this->min_rvir);
       this->So.message_screen("\tMax Rvir  =", this->max_rvir);
      }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
     this->So.message_screen("\tMean Rvir of tracer catalogue =", mean_m, "kpc/h");
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of Rvir not available in catalog or column not properñy set");
      So.DONE();
#endif
      // *********************************************************************
      // ALLOCATE CONCENTRATION
#ifdef _USE_CONCENTRATION_
      if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" )&& (i_rvir>0 && i_rvir<NCOLS) && (i_rs>0 && i_rs<NCOLS) )
       {
#ifdef _VERBOSE_CAT_
          So.message_screen("Computing concentration for tracers (c=Rvir/Rs)");
#endif
          real_prec mean_m=0;
          iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        if(prop[i_setting+i*this->NCOLS]>min_cut)
          {
#endif
#endif
#ifdef _POWER_
            real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
            real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
             {
               real_prec proper= prop[i_rs+i*NCOLS]>0 ? prop[i_rvir+i*NCOLS]/prop[i_rs+i*NCOLS]:1.0 ;
               this->Halo[iN].concentration=proper;
                mean_m+=proper;
                iN++;
              }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
        }
      if(true==this->params._set_bins_equal_number_tracers())
      {
        this->min_rs=get_min("_CONCENTRATION_");
        this->max_rs=get_max("_CONCENTRATION_");
        vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
        vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
        this->get_intervals_equal_number("_CONCENTRATION_",min_aux,max_aux);
        for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_CONCENTRATIONbins_max(i,max_aux[i]);
        for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
            this->params.set_CONCENTRATIONbins_min(i,min_aux[i]);
        this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
        this->So.message_screen("\tMin Concentration  =", this->min_concentration);
        this->So.message_screen("\tMax Concentration  =", this->max_concentration);
      }
     else 
     {
       this->min_concentration=get_min(_CONCENTRATION_);
       this->max_concentration=get_max(_CONCENTRATION_);
       this->So.message_screen("\tMin Concentration  =", this->min_concentration);
       this->So.message_screen("\tMax Concentration  =", this->max_concentration);
     }
     if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
     this->So.message_screen("\tMean Concentration of tracer catalog =", mean_m);
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of Concentration was not fully available");
      So.DONE();
#endif

#endif // ENFIF USE_CONCENTRATION
   // *********************************************************************
   // ALLOCATE VIRIAL
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_virial>0 && i_virial<NCOLS)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating T/|W| Viral from tracers");
#endif
     real_prec mean_m=0;
     iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser>= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
          {
        this->Halo[iN].virial=prop[i_virial+i*NCOLS];
        mean_m+=prop[i_virial+i*NCOLS];
        iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
    }
      if(true==this->params._set_bins_equal_number_tracers())
      {
       this->min_virial=get_min("_VIRIAL_");
       this->max_virial=get_max("_VIRIAL_");
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       this->get_intervals_equal_number("_VIRIAL_",min_aux,max_aux);
       for(int i=0;i<this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_VIRIALbins_max(i,max_aux[i]);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_VIRIALbins_min(i,min_aux[i]);
      this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
      this->So.message_screen("\tMin Virial  =", this->min_virial);
      this->So.message_screen("\tMax Virial  =", this->max_virial);
    }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
        this->So.message_screen("\tMean Viral of tracer catalogue =", mean_m, "");
#endif
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
    }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of virial not available in catalog or column not properly set.");
      So.DONE();
#endif
  // *********************************************************************
  // ALLOCATE B_TO_A
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_b_to_a>0 && i_b_to_a<NCOLS) // This if is outside the loop, for we here want vmax and it applies only when this "if" is satisfied
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating b_to_a from tracer");
#endif
      real_prec mean_m=0;
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
#ifdef _POWER_
             real_prec prop_sel =prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec prop_sel =prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop_sel>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(prop_sel>= prop_min && prop_sel < prop_max)
#endif
             #ifdef _USE_REDSHIFT_BINS_
                 if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
                 {
                   real_prec proper=prop[i_b_to_a+i*NCOLS];
                   this->Halo[iN].b_to_a=proper;
                   mean_m+=proper;
                   iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        }
#endif
#endif
      }
      if(true==this->params._set_bins_equal_number_tracers())
      {
      this->min_b_to_a=get_min("_BTOA_");
      this->max_b_to_a=get_max("_BTOA_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      this->get_intervals_equal_number("_BTOA_",min_aux,max_aux);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_BTOAbins_max(i,max_aux[i]);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_BTOAbins_min(i,min_aux[i]);
      this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
      this->So.message_screen("\tMin b_to_a  =", this->min_b_to_a);
      this->So.message_screen("\tMax b_to_a  =", this->max_b_to_a);

      }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
        mean_m/=static_cast<real_prec>(this->NOBJS);
        this->So.message_screen("\tMean b_to_a of tracer catalogue =", mean_m);
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of b_to_a not available in catalog or column not properly set");
#endif
      So.DONE();
  // *********************************************************************
  // ALLOCATE C_TO_A
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS"|| this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_c_to_a>0 && i_c_to_a<NCOLS) // This if is outside the loop, for we here want vmax and it applies only when this "if" is satisfied
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating c_to_a from tracer");
#endif
      real_prec mean_m=0;
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
        {
#endif
#endif
#ifdef _POWER_
             real_prec prop_sel =prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec prop_sel =prop[i_setting+i*this->NCOLS]*units_settings;
#endif

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(prop_sel>= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
        if(prop_sel>= prop_min && prop_sel < prop_max)
#endif
             #ifdef _USE_REDSHIFT_BINS_
                 if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
                 {
                   real_prec proper=prop[i_c_to_a+i*NCOLS];
                   this->Halo[iN].c_to_a=proper;
                   mean_m+=proper;
                   iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
        }
#endif
#endif
      }
      if(true==this->params._set_bins_equal_number_tracers())
      {
      this->min_c_to_a=get_min("_CTOA_");
      this->max_c_to_a=get_max("_CTOA_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      this->get_intervals_equal_number("_CTOA_",min_aux,max_aux);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
         this->params.set_CTOAbins_max(i,max_aux[i]);
       for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
          this->params.set_CTOAbins_min(i,min_aux[i]);
     this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
      this->So.message_screen("\tMin c_to_a  =", this->min_c_to_a);
      this->So.message_screen("\tMax c_to_a  =", this->max_c_to_a);
    }
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
#ifdef _VERBOSE_CAT_
        mean_m/=static_cast<real_prec>(this->NOBJS);
        this->So.message_screen("\tMean c_to_a of tracer catalogue =", mean_m);
#endif
      }
#ifdef _VERBOSE_CAT_
  else
    So.message_warning("Information of c_to_a not available in catalog or column not properly set");
#endif
      So.DONE();
  // ********************************************************************************************************************************
  // GET ELLIPTICITY AND PROLATNESS FROM SEMIAXIS
   if(this->params._i_b_to_a_g()>0 && this->params._i_c_to_a_g()>0)
    {
#ifdef _USE_HALO_PROL_ELL_
#ifdef _FULL_VERBOSE_
     So.message_screen("Converting from (q,s) to (T,e)");
      So.message_screen("   in the code/output, (T->b/a, and e->c/a)");
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
     {
        real_prec s=this->Halo[i].c_to_a;
        real_prec q=this->Halo[i].b_to_a;
        real_prec prol=(1.-q*q)/(1.-s*s);
        real_prec ell=(1.-s*s)/(1.+s*s+q*q);
        this->Halo[i].c_to_a=ell;
        this->Halo[i].b_to_a=prol;
     }
    {
      this->min_b_to_a=get_min("_BTOA_");
      this->max_b_to_a=get_max("_BTOA_");
      vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
      this->get_intervals_equal_number("_BTOA_",min_aux,max_aux);
      for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
        this->params.set_BTOAbins_max(i,max_aux[i]);
      for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
        this->params.set_BTOAbins_min(i,min_aux[i]);

      this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
      this->So.message_screen("\tMin T =", this->min_c_to_a);
      this->So.message_screen("\tMax T =", this->max_c_to_a);
    }
    {
    this->min_c_to_a=get_min("_CTOA_");
    this->max_c_to_a=get_max("_CTOA_");
    vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
    vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
    this->get_intervals_equal_number("_CTOA_",min_aux,max_aux);
    for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
      this->params.set_CTOAbins_max(i,max_aux[i]);
    for(int i=0;i<=this->params._Number_of_bins_equal_number_tracers();++i)
      this->params.set_CTOAbins_min(i,min_aux[i]);

    this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
    this->So.message_screen("\tMin e =", this->min_c_to_a);
    this->So.message_screen("\tMax e =", this->max_c_to_a);
    So.DONE();
    }
#endif
    }
   // *********************************************************************
    // ALLOCATE VELOCITIES
   // Here we have to implement an if statement in order to select objects within a "property" bin or cut
  // before assigning the cat to the structure this->Halo. The quantity this->NOBS must be then
  // the remaining number of tracers after the cuts.
if (this->type_of_object!="RANDOM")
{
#if defined (_USE_VELOCITIES_TRACERS_) || defined (_REDSHIFT_SPACE_)
#ifdef _VERBOSE_CAT_
  So.message_screen("Allocating velocities from tracers: ");
#endif
    // vx:
  if(i_vx>0) //Best to evaluate only one if than NOBS if's
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("\tVx");
#endif
      iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN)
#endif
     for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser  >= prop_min && obser < prop_max)
#endif
          {
        this->Halo[iN].vel1=prop[i_vx+i*NCOLS];
        iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
          }
     So.DONE();
  }
  // vy:
  if(i_vy>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("\tVy");
#endif
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
    {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
              {
        this->Halo[iN].vel2=prop[i_vy+i*NCOLS];
        iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
              }
      So.DONE();
  }
  // vz:
  if(i_vz>0)
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("\tVz");
#endif
      iN=0;
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif

#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
            if(obser >= prop_min && obser < prop_max)
#endif
              {
        this->Halo[iN].vel3=prop[i_vz+i*NCOLS];
        iN++;
          }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
        }
    }
#endif
    So.DONE();
}
// *********************************************************************
 // ALLOCATE MVIR
  if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM" ) && i_mass >0 && i_mass <NCOLS)
    {
#ifdef _USE_MASS_TRACERS_
      real_prec mean_m=0;
      iN=0;
#ifdef _VERBOSE_CAT_
      So.message_screen("Allocating Mvir of tracers");
#endif
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
      for(ULONG i=0;i<NLINES;++i)
        {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]*units_settings>min_cut)
            {
#endif
#endif
#ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
         if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
          if(obser >= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
            {
              real_prec proper=prop[i_mass+i*this->NCOLS];
#ifdef _USE_LOG_MASS_
#ifdef _VERBOSE_CAT_
              mean_m+=log10(proper);
#endif
              this->Halo[iN].mass=log10(proper);
#else
              this->Halo[iN].mass=proper;
              mean_m+=proper;
#endif
              iN++;
        }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
            }
#endif
#endif
      }
#ifdef _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
      this->min_mass=get_min("_MASS_");
      this->max_mass=get_max("_MASS_");
      this->So.message_screen("\tMin Mass  =", this->min_mass);
      this->So.message_screen("\tMax Mass  =", this->max_mass);
      So.DONE();
      if(true==this->params._set_bins_equal_number_tracers_main_property()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
       {
         vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
         vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
         this->Number_of_tracers_in_mass_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
       this->get_intervals_equal_number("_MASS_",min_aux,max_aux);

       cout<<endl;
       So.message_screen("Mass bins with equal number of tracers identified");
       for(int i=0;i<=this->Number_of_tracers_in_mass_bins.size();++i)
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
     if(count_new_nobj!=iN)
      {
         So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
         So.message_warning("Line", __LINE__);
      }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
      if(iN==0)
         So.message_warning("No mass found in the desidered limits");
#endif
#else // If mass is not primary property then
      if(true==this->params._set_bins_equal_number_tracers()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
      {
       this->min_mass=get_min("_MASS_");
       this->max_mass=get_max("_MASS_");

       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers(),0);
       this->Number_of_tracers_in_mass_bins.resize(this->params._Number_of_bins_equal_number_tracers()+1,0);
       this->get_intervals_equal_number("_MASS_",min_aux,max_aux);
       for(int i=0;i<=this->Number_of_tracers_in_mass_bins.size();++i)
       {
          this->params.set_MASSbins_max(i,max_aux[i]);
          this->params.set_MASSbins_min(i,min_aux[i]);
       }
       this->So.message_screen("\tMin Mass  =", this->min_mass);
       this->So.message_screen("\tMax Mass  =", this->max_mass);
       this->So.message_screen("\tNumber of bins in this property =", this->params._Number_of_bins_equal_number_tracers());
     }
     if(count_new_nobj!=iN)
      {
         So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
         So.message_warning("Line", __LINE__);
      }
#ifdef _VERBOSE_CAT_
      mean_m/=static_cast<real_prec>(this->NOBJS);
      if(iN==0)
         So.message_warning("No mass found in the desidered limits");
#endif
#endif
#ifdef _VERBOSE_CAT_
#ifdef _USE_LOG_MASS_
      this->So.message_screen("\tMean mass of tracer catalogue =", pow(10,mean_m)*this->params._MASS_units(), "Ms/h");
#else
      this->So.message_screen("\tMean mass of tracer catalogue =", mean_m*this->params._MASS_units(), "Ms/h");
#endif
#endif
    }
#ifdef _VERBOSE_CAT_
  else
      So.message_screen("No Information on the tracer mass available in", this->type_of_object);
  So.DONE();
#endif
#endif
    cout<<endl;
    // *********************************************************************
    // ALLOCATE STAR FORMATION
#ifdef _USE_SAT_FRACTION_
  if(i_sf>0)
    {
      iN=0;
      int iN_nosat=0;
      for(ULONG i=0;i<NLINES;++i)
    if(prop[i_mass+i*NCOLS]*params._MASS_units()>= prop_min)
#ifdef _USE_REDSHIFT_BINS_
          if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
           {
           this->Halo[iN].number_sub_structures=static_cast<int>(prop[i_sf+i*NCOLS]);
        iN++;
        if(prop[i_sf+i*NCOLS] < 1)
              iN_nosat++;
      }
      double mean_m=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean_m)
#endif
      for(ULONG i=0;i<this->NOBJS;++i)
    mean_m+=static_cast<double>(this->Halo[i].number_sub_structures);

      mean_m/=static_cast<double>(this->NOBJS);
      this->So.message_screen("Mean number of satellites per host =", mean_m);
      this->So.message_screen("Fraction of centrals without sub-structure =", 100.0*static_cast<double>(iN_nosat)/static_cast<double>(this->NOBJS),"%");
    }
#endif
  // *********************************************************************
  // ALLOCATE STELLAR_MASS
 if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_stellar_mass>0)
   {
#ifdef _VERBOSE_CAT_
     So.message_screen("Allocating stellar mass from tracers");
#endif
    real_prec mean_m=0;
    iN=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:iN,mean_m)
#endif
     for(ULONG i=0;i<NLINES;++i)
       {
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
         if(prop[i_setting+i*this->NCOLS]>min_cut)
           {
#endif
#endif
#ifdef _POWER_
            real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
#else
            real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
#endif
#if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
         if(obser >= prop_min)
#elif defined (_USE_MASS_BINS_PK_)
           if(obser>= prop_min && obser < prop_max)
#endif
#ifdef _USE_REDSHIFT_BINS_
   if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
#endif
         {
       this->Halo[iN].stellar_mass=prop[i_stellar_mass+i*NCOLS];
       mean_m+=prop[i_stellar_mass+i*NCOLS];
       iN++;
         }
#ifdef _POWER_
#ifdef _SET_CAT_WITH_CUT_
           }
#endif
#endif
   }
     mean_m/=static_cast<real_prec>(this->NOBJS);
       this->So.message_screen("\tMean Stellar mass in catalog =", mean_m, "");
     if(count_new_nobj!=iN)
       {
         So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
         So.message_warning("Line", __LINE__);
       }
   }
#ifdef _VERBOSE_CAT_
 else
   So.message_warning("Information of stellar mass not available in catalog or column not properly set.");
     So.DONE();
#endif
     // *********************************************************************
     // ALLOCATE COLOR
 if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_color>0)
    {
   #ifdef _VERBOSE_CAT_
        So.message_screen("Allocating color from tracers");
   #endif
       real_prec mean_m=0;
       iN=0;
   #ifdef _USE_OMP_TEST_CAT_
   #pragma omp parallel for reduction(+:iN,mean_m)
   #endif
        for(ULONG i=0;i<NLINES;++i)
          {
   #ifdef _POWER_
   #ifdef _SET_CAT_WITH_CUT_
            if(prop[i_setting+i*this->NCOLS]>min_cut)
              {
   #endif
   #endif
   #ifdef _POWER_
               real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
   #else
               real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
   #endif
   #if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
            if(obser >= prop_min)
   #elif defined (_USE_MASS_BINS_PK_)
              if(obser>= prop_min && obser < prop_max)
   #endif
   #ifdef _USE_REDSHIFT_BINS_
      if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
   #endif
         {
          this->Halo[iN].color=prop[i_color+i*NCOLS];
          mean_m+=prop[i_color+i*NCOLS];
          iN++;
         }
   #ifdef _POWER_
   #ifdef _SET_CAT_WITH_CUT_
           }
   #endif
   #endif
      }
      mean_m/=static_cast<real_prec>(this->NOBJS);
      this->So.message_screen("\tMean color in catalog =", mean_m, "");
      if(count_new_nobj!=iN)
        {
          So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
          So.message_warning("Line", __LINE__);
        }
      }
#ifdef _VERBOSE_CAT_
    else
     So.message_warning("Information of color not available in catalog or column not properly set.");
    So.DONE();
#endif

    // *********************************************************************
    // ALLOCATE APP_MAG
 if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_app_mag>0)
   {
  #ifdef _VERBOSE_CAT_
       So.message_screen("Allocating apparent magnitude from tracers");
  #endif
      real_prec mean_m=0;
      iN=0;
  #ifdef _USE_OMP_TEST_CAT_
  #pragma omp parallel for reduction(+:iN,mean_m)
  #endif
       for(ULONG i=0;i<NLINES;++i)
         {
  #ifdef _POWER_
  #ifdef _SET_CAT_WITH_CUT_
           if(prop[i_setting+i*this->NCOLS]>min_cut)
             {
  #endif
  #endif
  #ifdef _POWER_
              real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
  #else
              real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
  #endif
  #if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
           if(obser >= prop_min)
  #elif defined (_USE_MASS_BINS_PK_)
             if(obser>= prop_min && obser < prop_max)
  #endif
  #ifdef _USE_REDSHIFT_BINS_
     if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
  #endif
        {
         this->Halo[iN].app_mag=prop[i_app_mag+i*NCOLS];
         mean_m+=prop[i_app_mag+i*NCOLS];
         iN++;
        }
  #ifdef _POWER_
  #ifdef _SET_CAT_WITH_CUT_
          }
  #endif
  #endif
     }
     mean_m/=static_cast<real_prec>(this->NOBJS);
     this->So.message_screen("\tMean apparaemt magnitude in catalog =", mean_m, "");
     if(count_new_nobj!=iN)
       {
         So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
         So.message_warning("Line", __LINE__);
       }
     }
#ifdef _VERBOSE_CAT_
   else
    So.message_warning("Information of apparent magnitude not available in catalog or column not properly set.");
   So.DONE();
#endif
   // *********************************************************************
   // ALLOCATE ABS_MAG
if((this->type_of_object!="TRACER_MOCK_ONLY_COORDS" || this->type_of_object!="TRACER_REF_ONLY_COORDS" || this->type_of_object!="RANDOM") && i_abs_mag>0)
  {
 #ifdef _VERBOSE_CAT_
      So.message_screen("Allocating absoluite magnitude from tracers");
 #endif
     real_prec mean_m=0;
     iN=0;
 #ifdef _USE_OMP_TEST_CAT_
 #pragma omp parallel for reduction(+:iN,mean_m)
 #endif
      for(ULONG i=0;i<NLINES;++i)
        {
 #ifdef _POWER_
 #ifdef _SET_CAT_WITH_CUT_
          if(prop[i_setting+i*this->NCOLS]>min_cut)
            {
 #endif
 #endif
 #ifdef _POWER_
             real_prec obser=prop[i_observable+i*this->NCOLS]*units_settings;
 #else
             real_prec obser=prop[i_setting+i*this->NCOLS]*units_settings;
 #endif
 #if defined (_USE_MASS_CUTS_PK_) || defined (_USE_ALL_PK_)
          if(obser >= prop_min)
 #elif defined (_USE_MASS_BINS_PK_)
            if(obser>= prop_min && obser < prop_max)
 #endif
 #ifdef _USE_REDSHIFT_BINS_
    if(prop[i_redshift+i*this->NCOLS]< this->params._redshift_max_sample() && prop[i_redshift+i*this->NCOLS]>= this->params._redshift_min_sample())
 #endif
       {
        this->Halo[iN].app_mag=prop[i_abs_mag+i*NCOLS];
        mean_m+=prop[i_abs_mag+i*NCOLS];
        iN++;
       }
 #ifdef _POWER_
 #ifdef _SET_CAT_WITH_CUT_
         }
 #endif
 #endif
    }
    mean_m/=static_cast<real_prec>(this->NOBJS);
    this->So.message_screen("\tMean absolute magnitude in catalog =", mean_m, "");
    if(count_new_nobj!=iN)
      {
        So.message_warning("Counting not consistent in function", __PRETTY_FUNCTION__);
        So.message_warning("Line", __LINE__);
      }
    }
#ifdef _VERBOSE_CAT_
  else
   So.message_warning("Information of absolute magnitude not available in catalog or column not properly set.");
  So.DONE();
#endif
  // ********************************************************************************************************************************
  // THIS SECTION FILLS, IF REQUESTED, INFORMATION ON ENVIROMENTAL PROPERTIES SUCH AS MACH, DAH, TIDAL,
  if(true==this->params._Get_tracer_local_mach_number() || true==this->params._Get_local_overdensity() || true==this->params._Get_tracer_local_dach_number())
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
    if(true==this->params._Get_tracer_local_mach_number())
    {
      this->min_mach=get_min("_MACH_");
      this->max_mach=get_max("_MACH_");
      this->So.message_screen("\tMin Mach  =", this->min_mach);
      this->So.message_screen("\tMax Mach =", this->max_mach);
    }
    if(true==this->params._Get_local_overdensity())
    {
      this->min_local_overdensity=get_min("_LOCAL_OVERDENSITY_");
      this->max_local_overdensity=get_max("_LOCAL_OVERDENSITY_");
      this->So.message_screen("\tMin Local Clustering =", this->min_local_overdensity);
      this->So.message_screen("\tMax Local Clustering =", this->max_local_overdensity);
    }
   if(true==this->params._Get_tracer_local_dach_number())
    {
      this->min_dach=get_min("_DACH_");
      this->max_dach=get_max("_DACH_");
      this->So.message_screen("\tMin Dach  =", this->min_dach);
      this->So.message_screen("\tMax Dach =", this->max_dach);
    }
  }
  this->mean_number_density=static_cast<real_prec>(this->NOBJS)/pow(this->box.Lbox,3);
#ifdef _VERBOSE_FREEMEM_
 So.message_screen("Freeing memmory in ", __PRETTY_FUNCTION__);
#endif
 So.message_screen("Freeing memmory in ", __PRETTY_FUNCTION__);
 prop.clear(); prop.shrink_to_fit();
 So.DONE();

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
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This method reads and selects tracers in bins of props, with mins and max given as members of the param class
// Meant to be used in the PowerSpectrum Class, expexting randoms and galaxy catalogs with properties such
// as redshift, color, stellar mass.
// This is a reduced version of the class member read_catalog() above
void Catalog::read_catalog(string input_file, int i_bin1, int i_bin2, int i_bin3)
{
    //Hand made bins in redshift
    vector<real_prec> z_bins_min={0.0,0.0,0.3,0.6};
    vector<real_prec> z_bins_max={1.0,0.3,0.6,1.0};
    // in color
    vector<real_prec> c_bins_min={0.0,0.0,1.0,2.0};
    vector<real_prec> c_bins_max={3.0,1.0,2.0,3.0};
    // in log stellar mass
    vector<real_prec> ms_bins_min={ 8.5, 8.5, 9.5, 10.5};
    vector<real_prec> ms_bins_max={11.5, 9.5,10.5, 11.5};
    this->So.enter(__PRETTY_FUNCTION__);
    real_prec prop_min=0;
  int i_x=-1;
  int i_y=-1;
  int i_z=-1;
  int i_vx=-1;
  int i_vy=-1;
  int i_vz=-1;
  int i_mass=-1; int i_weight=-1; int i_mean_density=-1; int i_sf=-1; int i_vmax=-1; int i_rs=-1;int i_rvir=-1; int i_virial=-1; int i_spin=-1;int i_spin_bullock=-1;int i_vrms=-1;int i_b_to_a=-1;int i_c_to_a=-1;
  int i_redshift=-1;
  int i_stellar_mass=-1;
  int i_color=-1;
  int i_abs_mag=-1;
  int i_app_mag=-1;
  if(this->type_of_object=="TRACER" || this->type_of_object=="TRACER_REF" || this->type_of_object=="TRACER_MOCK" || this->type_of_object=="TRACER_MOCK_ONLY_COORDS"|| this->type_of_object=="TRACER_REF_ONLY_COORDS")
    {
      i_x= this->params._i_coord1_g();
      i_y= this->params._i_coord2_g();
      i_z= this->params._i_coord3_g();
      i_vx= this->params._i_v1_g();
      i_vy= this->params._i_v2_g();
      i_vz= this->params._i_v3_g();
      i_mass = this->params._i_mass_g();
      i_vmax = this->params._i_vmax_g();
      i_vrms = this->params._i_vrms_g();
      i_weight = this->params._i_weight1_g();
      i_mean_density= this->params._i_mean_density_g();
      i_sf = this->params._i_sf_g();
      i_rs = this->params._i_rs_g();
      i_rvir = this->params._i_rvir_g();
      i_virial = this->params._i_virial_g();
      i_spin = this->params._i_spin_g();
      i_spin_bullock = this->params._i_spin_bullock_g();
      i_b_to_a = this->params._i_b_to_a_g();
      i_c_to_a = this->params._i_c_to_a_g();
      i_stellar_mass = this->params._i_stellar_mass_g();
      i_color = this->params._i_color_g();
      i_app_mag = this->params._i_app_mag_g();
      i_abs_mag = this->params._i_abs_mag_g();
      if(this->params._sys_of_coord_g()==2)
        i_redshift= this->params._i_coord3_g();// note that this applies only in the case in which coordsa are in pseudo*-equatorial
    }
  else if(this->type_of_object=="RANDOM")
    {
      i_x= this->params._i_coord1_r();
      i_y= this->params._i_coord2_r();
      i_z= this->params._i_coord3_r();
      i_weight = this->params._i_weight1_r();
      i_mean_density= this->params._i_mean_density_r();
      i_mass= this->params._i_mass_r();
      i_stellar_mass = this->params._i_stellar_mass_r();
      i_color = this->params._i_color_r();
      i_app_mag = this->params._i_app_mag_r();
      i_abs_mag = this->params._i_abs_mag_r();
      if(this->params._sys_of_coord_r()>=2)
        i_redshift= this->params._i_coord3_r();// note that this applies only in the case in which coords are in pseudo-equatorial
    }
  real_prec units_observable=1;
  real_prec units_settings=1;
  int i_setting=i_mass;
  // **************************************************************************************************************
  vector<real_prec> prop;
  ULONG NLINES =0;
     // This is meant for ascie reading.
  NLINES = this->File.read_file(input_file, prop,_NTHREADS_);
  this->NCOLS=(static_cast<ULONG>(prop.size()/NLINES));
  So.message_screen("Counting number of", this->type_of_object);
  // **************************************************************************************************************
  // We now proceed to read and allcoate the different properties
  // taking into account the different limits on properties.
  // This is done individually, since putting Nob ifs inside a M loop costs more than a single "if" outside each M loop
   ULONG count_new_nobj=0;
   this->Halo.clear();this->Halo.shrink_to_fit();
   for(ULONG i=0;i<NLINES;++i)
      {
         real_prec zr=prop[i_redshift+i*this->NCOLS];
         real_prec cr=prop[i_color+i*this->NCOLS];
         real_prec mr=prop[i_stellar_mass+i*this->NCOLS];
         if( (zr>z_bins_min[i_bin1] && zr<=z_bins_max[i_bin1]) && ( cr>c_bins_min[i_bin2] && cr<=c_bins_max[i_bin2]) && (mr>ms_bins_min[i_bin3] && mr<=ms_bins_max[i_bin3]))
           {
             this->Halo.push_back(s_Halo());
             this->Halo[count_new_nobj].coord1=prop[i_x+i*NCOLS];
             this->Halo[count_new_nobj].coord2=prop[i_y+i*NCOLS];
             this->Halo[count_new_nobj].coord3=prop[i_z+i*NCOLS];
             count_new_nobj++;
            }
         }
   prop.clear();
   prop.shrink_to_fit();
   So.DONE();
   So.message_screen("Found ",count_new_nobj," objects");
   this->NOBJS=this->Halo.size();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_density_field_grid(string prop, string output_file)
{
  this->So.enter(__PRETTY_FUNCTION__);
#ifdef _VERBOSE_CAT_
  So.message_screen("Interpolating on a grid using MAS= ", this->params._masskernel());
#endif
  vector<real_prec> deltaTR_counts(this->box.NGRID,0);
  MAS_NEW(&box,this->Halo, _COUNTS_, deltaTR_counts);
  vector<real_prec> deltaTR;
  if(_COUNTS_ == prop)
    this->File.write_array(output_file, deltaTR_counts);
  else
    {
      deltaTR.resize(this->box.NGRID,0);
      MAS_NEW(&box,this->Halo, prop, deltaTR);
      if(_MASS_ == prop || _RS_== prop || _VIRIAL_== prop || _SPIN_== prop || _VMAX_== prop)
        {
      So.message_screen("Interpolating tracer on a grid weighting by", prop);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0 ; i< this->box.NGRID; ++i)
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
void Catalog::get_density_field_grid(string prop, vector<real_prec>&out)
{
  So.message_screen("Interpolating on a grid using MAS", this->params._masskernel());
  vector<real_prec> deltaTR_counts(this->box.NGRID,0);
  MAS_NEW(&box,this->Halo, _COUNTS_, deltaTR_counts);
  vector<real_prec> deltaTR;
  if(_COUNTS_ == prop)
    out=deltaTR_counts;
  else
    {
      deltaTR.resize(this->box.NGRID,0);
      MAS_NEW(&box,this->Halo, prop, deltaTR);
   if(_MASS_ == prop)
    {
      So.message_screen("Interpolating tracer on a grid weighting by mass");
      for(ULONG i=0 ; i< this->box.NGRID; ++i)
        if(deltaTR_counts[i]!=0)
          out[i]=static_cast<double>(deltaTR[i])/static_cast<double>(deltaTR_counts[i]); // Mean mass in cells
        else
          out[i]=0;
    }
    else if (_SAT_FRACTION_==prop)
    {
      So.message_screen("Interpolating Number of satellites in a grid");
      out=deltaTR;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_density_field_grid(s_params_box_mas sbox, string prop, vector<real_prec>&out)
{
  So.message_screen("Interpolating on a grid using MAS", this->params._masskernel());
  vector<real_prec> deltaTR_counts(this->box.NGRID,0);
  MAS_NEW(&sbox,this->Halo, _COUNTS_, deltaTR_counts);
  vector<real_prec> deltaTR;
  if(_COUNTS_ == prop)
    out=deltaTR_counts;
  else
    {
      deltaTR.resize(sbox.NGRID,0);
      MAS_NEW(&sbox,this->Halo, prop, deltaTR);
      if(_MASS_ == prop)
    {
      So.message_screen("Interpolating tracer on a grid weighting by mass");
      for(ULONG i=0 ; i< sbox.NGRID; ++i)
        if(deltaTR_counts[i]!=0)
          out[i]=static_cast<double>(deltaTR[i])/static_cast<double>(deltaTR_counts[i]); // Mean mass in cells
        else
          out[i]=0;
    }
    else if (_SAT_FRACTION_==prop)
    {
      So.message_screen("Interpolating Number of satellites in a grid");
      out=deltaTR;
    }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::define_property_bins()
{
    So.enter(__PRETTY_FUNCTION__);
    cout<<__LINE__<<endl;
#ifdef _USE_OMP_
   int NTHREADS = _NTHREADS_;  
   omp_set_num_threads(NTHREADS);
#endif
#ifdef _VERBOSE_CAT_
#ifdef MBINS
   So.message_screen("Defining bins for abundance in catalog type ", this->type_of_object);
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
      for(int i=0;i<this->NMBINS;++i)
        this->MBmin[i]=pow(10,MMin+i*logdeltaM);
#else
      for(int i=0;i<this->NMBINS;++i)
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
      for(int i=0;i<this->NMBINS;++i)
        {
          this->MBmax[i]=pow(10,MMin+(i+1)*logdeltaM);
          this->MBin[i]=this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf());
        }
#elif defined MCUTS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->NMBINS;++i)
        this->MBmax[i]=pow(10, MMax);
#endif
#else
#ifdef MBINS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<NMBINS;++i)
        this->MBmax[i]=(MMin+(i)*logdeltaM);
#elif defined MCUTS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<NMBINS;++i)
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
      So.message_screen("No halo mass info available for type ", this->type_of_object);
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
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmin[i]=pow(10,Min+i*logdeltaVMAX);
      this->VMAXBmax.clear();this->VMAXBmax.shrink_to_fit();this->VMAXBmax.resize(this->NMBINS,0);
      this->VMAXBin.clear();this->VMAXBin.shrink_to_fit();this->VMAXBin.resize(this->NMBINS,0);
#ifdef MBINS
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(int i=0;i<this->NMBINS;++i)
        {
          this->VMAXBmax[i]=pow(10,Min+(i+1)*logdeltaVMAX);
          this->VMAXBin[i]=(this->VMAXBmin[i]+(i+0.5)*(this->VMAXBmax[i]-this->VMAXBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));
        }
#elif defined MCUTS
      for(int i=0;i<this->NMBINS;++i)
        this->VMAXBmax.push_back(pow(10, MMax));
#endif
    }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No VMAX info available for type ", this->type_of_object);
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
      for(int i=0;i<this->NMBINS;++i)
        this->RSBmin.push_back(pow(10,RSMin+i*this->logdeltaRS));
#else
      for(int i=0;i<this->NMBINS;++i)
        this->RSBmin.push_back(MMin+i*logdeltaVMAX);
#endif
#ifdef _RS_LOG_
      this->RSBmax.clear();this->RSBmax.shrink_to_fit();
      this->RSBin.clear();this->RSBin.shrink_to_fit();
#ifdef RSBINS
      for(int i=0;i<this->NMBINS;++i)
        this->RSBmax.push_back(pow(10,RSMin+(i+1)*logdeltaRS));
      for(int i=0;i<this->NMBINS;++i)
        this->RSBin.push_back(this->RSBmin[i]+(i+0.5)*(this->RSBmax[i]-this->RSBmin[i])/static_cast<double>(this->params._NMASSbins_mf()));
#elif defined MCUTS
      for(int i=0;i<this->NMBINS;++i)
        this->RSBmax.push_back(pow(10, RSMax));
#endif
#else
#ifdef RSBINS
      for(int i=0;i<NMBINS;++i)
        this->RSBmax.push_back(RSMin+(i)*logdeltaRS);
#elif defined MCUTS
      for(int i=0;i<NMBINS;++i)
        this->RSmax.push_back(RSMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No RS info available for type ", this->type_of_object);
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
      for(int i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmin.push_back(pow(10,SMin+i*this->logdeltaCONCENTRATION));
#else
      for(int i=0;i<this->NMBINS;++i)
        this->SPINBmin.push_back(SMin+i*logdeltaSPIN);
#endif
#ifdef _CONCENTRATION_LOG_
      this->CONCENTRATIONBmax.clear();this->CONCENTRATIONBmax.shrink_to_fit();
#ifdef CONCENTRATIONBINS
      for(int i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(pow(10,SMin+(i+1)*logdeltaCONCENTRATION));
      this->CONCENTRATIONBin.clear();this->CONCENTRATIONBin.shrink_to_fit();
      for(int i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBin.push_back(pow(10,SMin+(i+0.5)*this->logdeltaCONCENTRATION));
#elif defined CONCENTRATIONCUTS
      for(int i=0;i<this->NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(pow(10, SMax));
#endif
#else
#ifdef CONCENTRATIONBINS
      for(int i=0;i<NMBINS;++i)
        this->CONCENTRATIONBmax.push_back(SMin+(i)*logdeltaSPIN);
#elif defined MCUTS
      for(int i=0;i<NMBINS;++i)
        this->CONCENTRATIONmax.push_back(SMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No Cvir info available for type ", this->type_of_object);
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
      for(int i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBmin.push_back(pow(10,SMin+i*this->logdeltaSPIN));
#else
      for(int i=0;i<this->NMBINS;++i)
        this->SPINBmin.push_back(SMin+i*logdeltaSPIN);
#endif
#ifdef _SPIN_LOG_
      this->SPINBmax.clear();this->SPINBmax.shrink_to_fit();
      this->SPINBin.clear();this->SPINBin.shrink_to_fit();
#ifdef SPINBINS
      for(int i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBmax.push_back(pow(10,SMin+(i+1)*this->logdeltaSPIN));
      for(int i=0;i<this->params._NMASSbins_mf();++i)
        this->SPINBin.push_back(pow(10,SMin+(i+0.5)*this->logdeltaSPIN));  
#elif defined SPINCUTS
      for(int i=0;i<this->NMBINS;++i)
        this->SPINBmax.push_back(pow(10, SMax));
#endif
#else
#ifdef SPINBINS
      for(int i=0;i<NMBINS;++i)
        this->SPINBmax.push_back(SMin+(i)*logdeltaSPIN);
#elif defined MCUTS
      for(int i=0;i<NMBINS;++i)
        this->SPINmax.push_back(SMax);
#endif
#endif
  }
  else
    {
#ifdef _VERBOSE_CAT_
      So.message_screen("No SPIN info available for type ", this->type_of_object);
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
void Catalog::define_property_bins(ULONG Nbins, string prop)
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
      for(int i=0;i<Nbins;++i)
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
      for(int i=0;i<Nbins;++i)
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
      for(int i=0;i< Nbins;++i)
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
      for(int i=0;i< Nbins;++i)
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
void Catalog::get_property_function(string file)
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
    for(ULONG i=0;i<this->NOBJS;++i)
        {
         real_prec property=log10(this->Halo[i].mass)+log10(units_observable_m);
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
      ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
       naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->NOBJS-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->NOBJS-naux);
#endif
      string file_m=file+"_mass";
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Mass function written in file ", file_m);
#endif
     for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
           mout<<this->MBin[i]<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();
#ifdef _VERBOSE_CAT_
      So.DONE();
#endif
      /*
#ifdef _USE_GNUPLOT_ABUNDANCE_PLOT_
   vector<pair<real_prec, real_prec> > xy_pts_ref;
   for(int i=0; i<this->mass_function.size(); ++i)
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
    for(int i=0; i< this->params._NMASSbins_mf(); ++i)
     for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
      mcount[i]+=mf_counts[j]/static_cast<double>(this->NOBJS);
    file_m+="_cumulative";
    this->File.write_to_file(file_m,this->MBin,mcount); 
#ifdef _VERBOSE_CAT_
      So.message_screen("Cumulative Mass function written in file ", file_m);
#endif
      So.DONE();
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
      for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec property=log10(this->Halo[i].vmax);
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
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
        naux+=v_counts[i];

#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Vmax) = ", naux);
      if(this->NOBJS-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Vmax) estimation = ", this->NOBJS-naux);
#endif
      string filev=file+"_vmax";
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting VMAX-function in file ", filev);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
        {
          real_prec mint=(this->VMAXBmax[i]-this->VMAXBmin[i]);
          this->vmax_function[i]=static_cast<real_prec>(v_counts[i])/(pow(this->params._Lbox(),3)*mint);
       }
      this->File.write_to_file(filev,this->VMAXBin,this->vmax_function,v_counts); 
      vector<real_prec>Vcount(this->params._NMASSbins_mf(),0);
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
      for(int i=0; i< this->params._NMASSbins_mf(); ++i) // N(>Vmax)
        for(int j=i+1; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
            Vcount[i]+=v_counts[j]/static_cast<double>(this->NOBJS);
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
          for(int i=0; i< this->params._NMASSbins_mf(); ++i)
          {
          gsl_real ppa=static_cast<gsl_real>(1.-Vcount[i]/static_cast<double>(this->NOBJS));
          //          gsl_real ppa=static_cast<gsl_real>(Vcount[i]/static_cast<double>(this->NOBJS));
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
      vmax_aux=gsl_vector_alloc(this->NOBJS);

      gsl_vector *id_gal;
      id_gal=gsl_vector_alloc(this->NOBJS);

      for(ULONG ig=0;ig<this->NOBJS;++ig )
        {
          gsl_vector_set(vmax_aux,ig,this->Halo[ig].vmax);
          gsl_vector_set(id_gal,ig,ig);
        }
      gsl_sort_vector2(vmax_aux,id_gal) ;   // sort the vmax and correspondingly their associated the gal id
      this->Prop_threshold_rand_dm=gsl_vector_get(vmax_aux, this->Ntracers_ran-1);// NOTE THAT this value won't be used in practice
#ifdef _VERBOSE_CAT_
      cout<<endl;
      So.message_screen("Threshold for Vmax assigment to randoms = ", this->Prop_threshold_rand_dm," km/s");
      cout<<endl;
#endif
      for(ULONG i=0;i<this->NOBJS;++i ) //ORDERED loop in ascending order in vmax: i=0 is the lowest vmax value, i=NOBJS is the highest
        {
          // To each halo, associates the index i sorted in ascening order with respect to vmax
          // This will be used when assigning Vmax to the random tracers. Since we have Nrandom, we select the refernece tracers with this index <Nrandom to prodive their vmax to the random tracers
          // such athat these will properñly probe the Vmax end of the n(Vmax) function.
          ULONG gal_id=gsl_vector_get(id_gal,i);
          this->Halo[gal_id].vmax_index=i;
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
      for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec property=log10(this->Halo[i].rs*units_observable_rs);
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
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
      naux+=rs_counts[i];
#ifdef _VERBOSE_CAT_
    So.message_screen("Number of objects used in mass function n(Rs) = ", naux);
    if(this->NOBJS-naux>0)
      So.message_screen("Number of objects EXCLUDED in n(Rs) estimation = ", this->NOBJS-naux);
#endif
    string filers=file+"_rs";

#ifdef _VERBOSE_CAT_
    So.message_screen("Writting Rs-function in file ", filers);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
    {
      real_prec mint=(this->RSBmax[i]-this->RSBmin[i]);
      this->rs_function[i]=static_cast<real_prec>(rs_counts[i])/(pow(this->params._Lbox(),3)*mint);
    }
   this->File.write_to_file(filers,this->RSBin,this->rs_function,s_counts); 
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for
#endif
    for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
     for(ULONG i=0;i<this->NOBJS;++i)
      {
       real_prec property=log10(this->Halo[i].spin_bullock);
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
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
      naux+=s_counts[i];
#ifdef _USE_OMP_
#pragma omp parallel for 
#endif
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
    for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
     for(ULONG i=0;i<this->NOBJS;++i)
      {
       real_prec property=log10(this->Halo[i].concentration);
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
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
      naux+=s_counts[i];
#ifdef _VERBOSE_CAT_
    So.message_screen("Number of objects used in mass function n(cvir) = ", naux);
    if(this->NOBJS-naux>0)
      So.message_screen("Number of objects EXCLUDED in n(cvir) estimation = ", this->NOBJS-naux);
#endif
    string filers=file+"_cvir";
#ifdef _VERBOSE_CAT_
    So.message_screen("Writting Cvir-function in file ", filers);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
    for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Estimates of Abundance as a function of halo properties
void Catalog::get_property_function_cosmic_web_types(string file)
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
        ULONG partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
        for(ULONG i=0;i<this->NOBJS;++i)
          {
           real_prec property=log10(this->Halo[i].mass)+log10(units_observable_m);
#ifndef accum
            if(property >=lm_min && property<=lm_max && this->cwclass.get_Tclassification(this->Halo[i].GridID)==ict)
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
        ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
        for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
          naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
        So.message_screen("Number of objects used in mass function n(M) = ", naux);
        if(this->NOBJS-naux>0)
          So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->NOBJS-naux);
#endif
        string file_m=file+"_mass_cwt"+to_string(ict);
        ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
    for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic update
#endif
          mcount[i]+=mf_counts[j];
     mf_counts.clear(); mf_counts.shrink_to_fit();

      file_m+="_cumulative_cwt"+to_string(ict);
      mout.open(file_m.c_str());
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
     ULONG partial_o=0;

#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
      for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec property=log10(this->Halo[i].vmax*units_observable_v);
#ifndef accum
      if(property >=v_min && property<=v_max && this->cwclass.get_Tclassification(this->Halo[i].GridID)==ict)
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
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
        naux+=v_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Vmax) = ", naux);
      if(this->NOBJS-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Vmax) estimation = ", this->NOBJS-naux);
#endif
      string filev=file+"_vmax_cwt"+to_string(ict);
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting VMAX-function in file ", filev);
#endif
      ofstream mout;
      mout.open(filev.c_str());
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
      for(int i=0; i< this->params._NMASSbins_mf(); ++i) // N(>Vmax)
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
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
      for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec property=log10(this->Halo[i].rs*units_observable_rs);
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
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
    naux+=rs_counts[i];

#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(Rs) = ", naux);
      if(this->NOBJS-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(Rs) estimation = ", this->NOBJS-naux);
#endif
      string filers=file+"_rs_cwt"+to_string(ict);

#ifdef _VERBOSE_CAT_
      So.message_screen("Writting Rs-function in file ", filers);
#endif
      ofstream mout;
      mout.open(filers.c_str());
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
        for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#pragma omp atomic
          RScount[i]+=this->rs_function[j]*(this->RSBmax[j]-this->RSBmin[j]);

      filers+="_cumulative_cwt"+to_string(ict);
      mout.open(filers.c_str());
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
      for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec property=log10(this->Halo[i].spin*units_observable_s);
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
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
    naux+=s_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(spin) = ", naux);
      if(this->NOBJS-naux>0)
         So.message_screen("Number of objects EXCLUDED in n(spin) estimation = ", this->NOBJS-naux);
#endif
      string filers=file+"_spin_cwt"+to_string(ict);
#ifdef _VERBOSE_CAT_
      So.message_screen("Writting Spin-function in file ", filers);
#endif
      ofstream mout;
      mout.open(filers.c_str());
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
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
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
        for(int j=i; j< this->params._NMASSbins_mf(); ++j)
#ifdef _USE_OMP_TEST_CAT_
#pragma omp atomic
#endif
            Scount[i]+=this->s_function[j]*(this->SPINBmax[j]-this->SPINBmin[j]);
      filers+="_cumulative_cwt"+to_string(ict);
      mout.open(filers.c_str());
      for(int i=0; i< this->params._NMASSbins_mf(); ++i)
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
void Catalog::get_HOD(string file)
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
      ULONG partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
    for(ULONG i=0;i<this->NOBJS;++i)
        {
         real_prec property=log10(this->Halo[i].mass)+log10(units_observable_m);
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
      ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
       naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->NOBJS-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->NOBJS-naux);
#endif
      string file_m=file+"_mass";
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
         mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();
      So.DONE();

  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_HOD_web_types(string file)
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
      ULONG partial_o=0;
#ifdef _USE_OMP_TEST_CAT_
#pragma omp parallel for reduction(+:partial_o)
#endif
    for(ULONG i=0;i<this->NOBJS;++i)
        {
         real_prec property=log10(this->Halo[i].mass)+log10(units_observable_m);
#ifndef accum
          if(property >=lm_min && property<=lm_max && this->cwclass.get_Tclassification(this->Halo[i].GridID)==ict)
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
      ULONG naux=0;
#pragma omp parallel for reduction(+:naux)
      for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
    naux+=mf_counts[i];
#ifdef _VERBOSE_CAT_
      So.message_screen("Number of objects used in mass function n(M) = ", naux);
      if(this->NOBJS-naux>0)
        So.message_screen("Number of objects EXCLUDED in n(M) estimation = ", this->NOBJS-naux);
#endif
      string file_m=file+"_mass_cwt"+to_string(ict);
      ofstream mout (file_m.c_str());
#ifdef _VERBOSE_CAT_
      So.message_screen("Writing mass function written in file ", file_m);
#endif
     for(ULONG i=0;i< this->params._NMASSbins_mf();++i)
         {
          real_prec mint=(this->MBmax[i]-this->MBmin[i]);
          this->mass_function[i]=static_cast<real_prec>(mf_counts[i])/(pow(this->params._Lbox(),3)*mint);
         mout<<this->MBmin[i]+(i+0.5)*(this->MBmax[i]-this->MBmin[i])/static_cast<double>(this->params._NMASSbins_mf())<<"\t"<<this->mass_function[i]<<"   "<<mf_counts[i]<<endl;
         }
      mout.close();
      So.DONE();
    }// closes loop pver cwtypes

  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog:: get_distribution_reduced_mass_in_cell(){
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
  So.message_screen("Measuring distribution of reduced mass in cells", MAXIMUM_DISTANCE_EXCLUSION, "Mpc/h");
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
#pragma omp parallel for
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].mass.push_back(this->Halo[i].mass);
    }

  //  vector<int>dist_reduced_mass(N_BINS_REDUCED_MASS, 0);

  real_prec dmin=MINIMUM_DISTANCE_EXCLUSION;
  real_prec dmax=MAXIMUM_DISTANCE_EXCLUSION;
  vector<vector<int>> dist_reduced_mass(N_BINS_REDUCED_MASS,vector<int>(N_BINS_DIST_EX, 0));
  for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
    if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
      for(int i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
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
      string file=this->Output_directory+"reduced_mass_dist_"+this->type_of_object+"_dbin"+to_string(id)+".txt";
      So.message_screen("Writting to file ", file);
      ofstream sal;
      sal.open(file.c_str());
      sal<<"#Bin info:     Dmin = "<<MINIMUM_DISTANCE_EXCLUSION+(id)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<"  Dmax = "<<MINIMUM_DISTANCE_EXCLUSION+(id+1.0)*(dmax-dmin)/static_cast<double>(N_BINS_DIST_EX)<<endl;
      for(int i=0;i<dist_reduced_mass.size();++i)
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
void Catalog::get_stats_separation_in_cell()
{
#ifdef _FULL_VERBOSE_
    this->So.enter(__PRETTY_FUNCTION__);
#endif
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring minimum separation in cells (Mpc/h) ...");
  So.message_screen("Current type is ", this->type_of_object);
  cout<<endl;
#endif
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
#ifdef _VERBOSE_CAT_
  So.message_screen("Getting positions inside cells:");
#endif
  if("TRACER_MOCK_ONLY_COORDS" != this->type_of_object)
    {
      for(ULONG i=0; i<this->NOBJS; ++i) //loop over the "observed objects", i.e, with cuts already set
    {
          ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
#ifdef _test_ms_
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
      cell_info_tr[ID].property.push_back(this->Halo[i].vmax);
#elif defined _USE_MASS_AS_PRIMARY_OBSERVABLE_
      cell_info_tr[ID].property.push_back(this->Halo[i].mass);
#endif // endif use_vmax_as _obs
#endif // endif _test_ms_
    }
    }
  else
    {
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed objects", i.e, with cuts already set
    {
          ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
    }
    }
  So.DONE();
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
//   if("TRACER_MOCK_ONLY_COORDS" != this->type_of_object)
//     {
        for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
          {
            real_prec aux_min_distance_a;
            if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
            {

#if defined _USE_MEAN_SEPARATIONS_IN_CELLS_  || defined  _USE_STDV_SEPARATIONS_IN_CELLS_
                vector<real_prec>sep_cell;
#endif
                  for(ULONG i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
            {
#if defined _USE_MIN_SEPARATIONS_IN_CELLS_
                      real_prec aux_min_distance_b=LARGE_NUMBER;
#endif
                          for(ULONG j=i+1 ; j < cell_info_tr[id].posx_p.size() ; ++j) // loop over ttracer in that very same cell
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
      for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
        {
          real_prec aux_min_distance_a;
          if(cell_info_tr[id].posx_p.size()>=2)//nos aseguramos de pasar por celdas que tengan al menos dos tracers para los cuales determinar la distancia. Si sólo hay uno, o ninguno, la separación es cero.
        {
          aux_min_distance_a=100.0;
          for(int i=0;i <cell_info_tr[id].posx_p.size();++i) // loop overtracers in that cell
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
      ULONG nbins_new=50;
#ifdef _USE_LOG_MSIC_
#pragma omp parallel for
      for(ULONG id=0; id<this->box.NGRID ;++id)
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
    for(ULONG id=0; id<this->box.NGRID ;++id)
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
        min_separation_in_cell[id]=log10(min_separation_in_cell[id]);
#elif defined _USE_MEAN_SEPARATIONS_IN_CELLS_
        mean_separation_in_cell[id]=log10(mean_separation_in_cell[id]);
#endif

#endif // end for _USE_LOG_MISC
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
     vector<ULONG>dist_min_sep(nbins_new,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
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
      vector<ULONG>dist_mean_sep(nbins_new,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
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
       vector<ULONG>dist_stdv_sep(nbins_new,0);
      for(ULONG id=0; id<this->box.NGRID ;++id) // loop over cells
          if(this->stdv_separation_in_cell[id]>MIN_STDV_SEP_IN_CELLS)
               dist_stdv_sep[get_bin(this->stdv_separation_in_cell[id],MIN_STDV_SEP_IN_CELLS,nbins_new,DELTA_STDV_SEP*N_BINS_STDV_SEP_IN_CELLS/nbins_new,true)]++;
       So.DONE();
       So.message_screen("Min value of stdv_sep_in_cell= ", get_min_nm(this->stdv_separation_in_cell));
       So.message_screen("Max value of stdv_sep_in_cell= ", get_max_nm(this->stdv_separation_in_cell));
#endif
      ofstream sal;
#ifdef _USE_MIN_SEPARATIONS_IN_CELLS_
      string fileo=this->Output_directory+"min_sep_in_cells_dist_"+this->type_of_object+".txt";
      So.message_screen("Writting distribution of minimum separations within cells");
      So.message_screen("in file ", fileo);
      sal.open(fileo.c_str());
      for(ULONG i=0;i<dist_min_sep.size();++i)
#ifdef _USE_LOG_MSIC_
          sal<<pow(10,MIN_SEP_IN_CELLS+(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new))<<"  "<<dist_min_sep[i]<<endl;
#else
          sal<<(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new)<<"  "<<dist_min_sep[i]<<endl;
#endif
      sal.close();
#
#endif
#ifdef _USE_MEAN_SEPARATIONS_IN_CELLS_
      string filemean=this->Output_directory+"mean_sep_in_cells_dist_"+this->type_of_object+".txt";
      So.message_screen("Writting distribution of mean separations within cells");
      So.message_screen("in file ", filemean);
      sal.open(filemean.c_str());
      for(ULONG i=0;i<dist_mean_sep.size();++i)
#ifdef _USE_LOG_MSIC_
          sal<<pow(10,MIN_SEP_IN_CELLS+(i+0.5)*DELTA_MIN_SEP*N_BINS_MIN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new))<<"  "<<dist_min_sep[i]<<endl;
#else
          sal<<(i+0.5)*DELTA_MEAN_SEP*N_BINS_MEAN_SEP_IN_CELLS/static_cast<real_prec>(nbins_new)<<"  "<<dist_mean_sep[i]<<endl;
#endif
      sal.close();
#endif
#ifdef _USE_STDV_SEPARATIONS_IN_CELLS_
      string filestdv=this->Output_directory+"stdv_sep_in_cells_dist_"+this->type_of_object+".txt";
      So.message_screen("Writting distribution of stdv of separations within cells");
      So.message_screen("in file ", filestdv);
      sal.open(filestdv.c_str());
      for(ULONG i=0;i<dist_stdv_sep.size();++i)
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
void Catalog::get_sep_distribution(real_prec min_value_prop)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring separation distribution");
  So.message_screen("Current type is ", this->type_of_object);
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
 for(ULONG ig=0; ig<this->NOBJS ;++ig) // loop over cells
  if(this->Halo[ig].vmax>min_value_prop)
   {
      x.push_back(this->Halo[ig].coord1);
      y.push_back(this->Halo[ig].coord2);
      z.push_back(this->Halo[ig].coord3);
    }
 real_prec normal=1;//static_cast<double>(x.size())*(x.size()-1)*0.5;
#ifdef _VERBOSE_CAT_
  So.message_screen("Getting separations");
  cout<<endl;
#endif
ULONG i=0,j=0;
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
    for(int i=0;i<dist_sep.size();++i)
      dist_sep[i]+=h_par[i];
   }
#endif
    ofstream sal;
    string fileo;
    fileo=this->Output_directory+"sep_dist_"+this->type_of_object+this->params._Name_survey()+"_Real"+to_string(this->params._realization())+".txt";
    So.message_screen("Writting distribution of minimum separations within cells");
    So.message_screen("in file ", fileo);
    sal.open(fileo.c_str());
    for(int i=0;i<dist_sep.size()-1;++i) // avoid the last bin which has all info from bins beyond
       sal<< pow(10,min_sep+(i+0.5)*delta_sep)<<" "<<dist_sep[i]/normal<<endl;
    sal.close();
    So.DONE();
    dist_sep.clear();
    dist_sep.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This method identifies the set of tracers around each tracer within a radius of  SCALE_MAX_N_NEIGHBOURS Mpc/h
// The output is allocated in a class member container Number_of_tracers with dimension Nobjects
void Catalog::get_neighbour_tracers(vector<s_nearest_cells>&nearest_cells_to_cell)
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
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
#endif
  // This arrawy will allcopate the number of neigh or the bin in distance in which the minn distance to pair falls
  this->Number_of_neighbours.clear();
  this->Number_of_neighbours.shrink_to_fit();
  this->Number_of_neighbours.resize(this->NOBJS,0);
#ifdef _USE_NUMBER_OF_NEIGHBOURS_
  So.message_screen("Computing number of neighbours for each tracer");
  So.message_screen("in spheres of radius ", Rscale, "Mpc/h");
#endif
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
  So.message_screen("Computing minimum separation to other tracer for each tracer");
#endif
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->Halo[i].coord1;
      real_prec y_coord=this->Halo[i].coord2;
      real_prec z_coord=this->Halo[i].coord3;
      ULONG ID=this->Halo[i].GridID;
#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
      vector<real_prec> min_distances_v(N_Neigh_cells,0);
#endif
      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
        {
          ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;

#ifdef _USE_MIN_DISTANCE_TO_NEIGHBOURS_
          aux_min_distance_b=BIG_NUMBER;
          real_prec min_distance_b;
#endif
          for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
            {
              ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
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
            this->Number_of_neighbours[i]++;
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
      this->Number_of_neighbours[i]=da_bin;
#endif
    }
  So.DONE();
#ifdef _USE_local_overdensity_
  So.message_screen("Computiong local clustering:");
  real_prec \tExpected_number_of_tracers = 4.0*M_PI*pow(SCALE_MAX_N_NEIGHBOURS,3)*this->mean_number_density/3.0;
  So.message_screen("\tExpected number of particles in the lcoal volume:", \tExpected_number_of_tracers);
  //  ofstream tea; tea.open("local_cl.txt");
#pragma omp parallel for
  for(ULONG i=0;i<this->NOBJS ;++i)
    {
      real_prec local_xi=log10(static_cast<real_prec>(this->Number_of_neighbours[i])/\tExpected_number_of_tracers+1);
      int index_xi=get_bin(local_xi,MIN_local_overdensity,N_BINS_MIN_DIST_TO_NEI, (MAX_local_overdensity-MIN_local_overdensity)/static_cast<real_prec>(N_BINS_MIN_DIST_TO_NEI),true);
      this->Number_of_neighbours[i]=index_xi;
    }
  //	tea.close();
  So.DONE();
#endif
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#define  _get_m_m_prob_
void Catalog::get_distribution_min_separations(vector<s_nearest_cells>&nearest_cells_to_cell)
{
    string data=this->type_of_object+"_"+this->params._Name_survey()+"_Real"+to_string(this->params._realization());

  if(this->type_of_object=="TRACER_MOCK")
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
  ULONG N_Neigh_cells=pow(max_neigh_per_dim,3)-2; // total number of cells to explore around a cell
  real_prec max_dist=(0.5*max_neigh_per_dim-0.5)*sqrt(3.0)*this->box.Lbox/static_cast<real_prec>(this->box.Nft);
#ifdef _VERBOSE_CAT_
  So.message_screen("Measuring distribution of pairs for ", this->type_of_object);
  cout<<endl;
  So.message_screen("Maximum separation probed =", max_dist,"Mpc /h");
  cout<<endl;
#endif
  struct s_cell_info{
    vector<real_prec> posx_p;
    vector<real_prec> posy_p;
    vector<real_prec> posz_p;
    vector<ULONG> gal_index;
  };
  //DEFINE VECTOR TO ALLOCATE THE NEIGHBOR CELLS OF EACH CELL
  // In order to gt ths miminum distance for each particle,
  // we nend to allocate for eachcell the coordintas of the particles belonging to it,
  // as we did when collapsing randoms:
  So.message_screen("Identifying tracers in cells");
  vector<s_cell_info> cell_info_tr(this->box.NGRID);
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);
      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
      cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  int N_bins_dist=50;
  int N_mass_bins=10;
  vector<vector<real_prec>> min_distance_dist(N_bins_dist,vector<real_prec>(N_mass_bins,0));
  vector<real_prec> distance_dist_m(N_bins_dist*N_mass_bins*N_mass_bins,0);
  real_prec lm_min=this->params._LOGMASSmin();
  real_prec lm_max=this->params._LOGMASSmax();
  this->logdeltaM=(lm_max-lm_min)/static_cast<real_prec>(N_mass_bins);
  real_prec aux_min_distance_a;
  real_prec aux_min_distance_b;
  ULONG pair_counting=0;
  So.message_screen("Computing minimum separations");
  // no paralelizar si hay un min o max() dentro del loop
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->Halo[i].coord1;
      real_prec y_coord=this->Halo[i].coord2;
      real_prec z_coord=this->Halo[i].coord3;
      real_prec lmass=log10(this->Halo[i].mass)+log10(this->params._MASS_units());
      int im=get_bin(lmass,lm_min,N_mass_bins,this->logdeltaM,true);
      ULONG ID=this->Halo[i].GridID;

      // Allocate the min_distances to i particle form neighbouring cells
      vector<real_prec> min_distances_v(N_Neigh_cells,0);

      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
       {
          ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->box.Lbox;
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->box.Lbox;
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->box.Lbox;
          aux_min_distance_b=1e5;
          real_prec min_distance_b=0;
          for(int k=0;k< cell_info_tr[ID_NEIGH].posx_p.size(); ++k) //loop sobre las particulas en las celdas cercanas
        {
          pair_counting++;
          ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
          real_prec lmass_k=log10(this->Halo[index_gal].mass)+log10(this->params._MASS_units());
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
              ULONG sindex = index_2d(da_bin,index_mass,N_mass_bins);
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
  for(int i=0;i<min_distance_dist.size();++i)
    {
      sal<<(i+0.5)*max_dist/static_cast<real_prec>(N_bins_dist)<<" ";
      for(int j=0;j< N_mass_bins; ++j)
        sal<<min_distance_dist[i][j]/static_cast<real_prec>(pair_counting)<<" ";
      sal<<endl;
    }
  sal.close();
#ifdef _get_m_m_prob_
  So.message_screen("Writting P(r|M,M') ");
  for(int i=0;i<min_distance_dist.size();++i)
    {
      string ff=this->params._Output_directory()+"separation_distribution_mass_"+data+"_separation_bin"+to_string(i)+".txt"+data;
      So.message_screen("in file ", ff);
      sal.open(ff.c_str());
      for(int j=0;j< N_mass_bins; ++j)
        for(int k=0;k< N_mass_bins; ++k)
          {
            int index=index_2d(j,k,N_mass_bins);
            ULONG sindex = index_2d(i,index,N_mass_bins);
            sal<<lm_min+(j+0.5)*logdeltaM<<" "<<lm_min+(k+0.5)*logdeltaM<<"  "<<static_cast<real_prec>(distance_dist_m[sindex])/static_cast<real_prec>(pair_counting)<<endl;
         }
      sal.close();
    }
#endif
  So.message_screen("Writting P(r) in mass bins ");
  string ff=this->params._Output_directory()+"separation_distribution_in_mass_bins_"+data+".txt";
  So.message_screen("in file ", ff);
  sal.open(ff.c_str());
  for(int i=0;i<min_distance_dist.size();++i)
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
void Catalog::get_masses_of_pairs_in_min_separation_bin_in_theta_bin(real_prec min_sep, vector<s_mass_members> & dm_properties_bins)
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
          //                ULONG Id1= dm_properties_bins[ih].GridID_bin_properties[j];
          for(int k=j+1;k< nobjs;++k)
        {
          //                  ULONG Id2= dm_properties_bins[ih].GridID_bin_properties[k];
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
  vector<ULONG>counts;
};
void Catalog::get_pdf_vmax(string property)
{

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG NGRID=this->params._Nft()*this->params._Nft()*this->params._Nft();
  vector<pdf_prop> vmax_in_cells(NGRID);
  int Nbins=30; //this->params._NMASSbins_mf();
  int Nocc=30; //this->params._NMASSbins_mf();
  if("VMAX"==property)
    {
      So.message_warning("Allcoating properties");
      for(ULONG i=0; i<this->NOBJS;++i)
       {
         ULONG ID=this->Halo[i].GridID;
         vmax_in_cells[ID].values.push_back(this->Halo[i].vmax);
       }
       So.DONE();
      real_prec delta=log10(this->params._VMAXmax()/this->params._VMAXmin())/static_cast<real_prec>(Nbins);
      So.message_warning("Getting histograms");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i< NGRID ;++i)
       {
         vmax_in_cells[i].counts.resize(Nbins,0);
         for(ULONG j=0;j<vmax_in_cells[i].values.size();++j)
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
       for(int i=0;i<Nbins;++i)
         for(ULONG j=0;j<NGRID;++j)
           if(k==vmax_in_cells[j].counts[i])
             pdf_c[index_2d(k,i,Nbins)]++;

     ofstream sa;
     sa.open(this->params._Output_directory()+"pdf_vmax.txt");
     for(int i=0;i<Nbins;++i)
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
void Catalog::select_random_subsample(real_prec fraction){
#else
void Catalog::select_random_subsample(string output_file, bool write_to_file, real_prec fraction){
#endif
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
   ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction*this->NOBJS));
   vector<real_prec>xnew;
   vector<real_prec>ynew;
   vector<real_prec>znew;
   vector<real_prec>vmax;
   vector<real_prec>vmax_parent;
   vector<ULONG>IDgrid;
   vector<ULONG>IDtheta;
   ULONG counter=0;
    do
     {
      ULONG i= gsl_rng_uniform_int(gBaseRando,this->NOBJS);
      xnew.push_back(this->Halo[i].coord1);
      ynew.push_back(this->Halo[i].coord2);
      znew.push_back(this->Halo[i].coord3);
      vmax.push_back(this->Halo[i].vmax);
      vmax_parent.push_back(this->Halo[i].vmax_parent);
      IDgrid.push_back(this->Halo[i].GridID);
      IDtheta.push_back(this->Halo[i].galThetaID);
      counter++;
    }while(counter<Nobjs_fraction);
    this->set_NOBJS_subsample(counter);
    this->Halo_random_subsample.resize(counter);
    for(ULONG i=0;i<counter;++i){
      this->Halo_random_subsample[i].coord1=xnew[i];
      this->Halo_random_subsample[i].coord2=ynew[i];
      this->Halo_random_subsample[i].coord3=znew[i];
      this->Halo_random_subsample[i].coord3=znew[i];
      this->Halo_random_subsample[i].vmax=vmax[i];
      this->Halo_random_subsample[i].vmax_parent=vmax_parent[i];
      this->Halo_random_subsample[i].GridID=IDgrid[i];
      this->Halo_random_subsample[i].galThetaID=IDtheta[i];
    }
    xnew.clear();xnew.shrink_to_fit();
    ynew.clear();ynew.shrink_to_fit();
    znew.clear();znew.shrink_to_fit();
    vmax.clear();vmax.shrink_to_fit();
    vmax_parent.clear();vmax_parent.shrink_to_fit();
    IDtheta.clear();IDtheta.shrink_to_fit();
    IDgrid.clear();IDgrid.shrink_to_fit();
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::select_random_subsample(real_prec fraction, string file){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);

   ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction*this->NOBJS));
   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(3);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);
   ULONG counter=0;
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
    do
     {
      int i= gsl_rng_uniform_int(gBaseRando,this->Halo.size());
      rcat<<log10(this->Halo[i].mass)<<"\t"<<log10(this->Halo[i].vmax)<<"\t"<<log10(this->Halo[i].concentration)<<"\t"<<log10(this->Halo[i].spin_bullock)<<"\t"<<this->Halo[i].mach_number<<"\t"<<this->Halo[i].bias<<"\t"<<log10(1+this->Halo[i].local_overdensity)<<endl;
//      rcat<<log10(this->Halo[i].mass)<<"\t"<<log10(this->Halo[i].vmax)<<"\t"<<log10(this->Halo[i].rs)<<"\t"<<log10(this->Halo[i].concentration)<<"\t"<<log10(this->Halo[i].spin_bullock)<<"\t"<<log10(this->Halo[i].vrms)<<"\t"<<this->Halo[i].virial<<"\t"<<this->Halo[i].b_to_a<<"\t"<<this->Halo[i].c_to_a<<"\t"<<this->Halo[i].mach_number<<"\t"<<this->Halo[i].bias<<"\t"<<this->Halo[i].qbias<<"\t"<<log10(1+this->Halo[i].local_overdensity)<<"\t"<<this->Halo[i].dach_number<<"\t"<<this->Halo[i].tidal_anisotropy<<"\t"<<this->Halo[i].peak_height<<"\t"<<this->Halo[i].gal_cwt<<"\t"<<this->Halo[i].relative_bias<<"\t"<<this->Halo[i].bias_rs<<"\t"<<this->Halo[i].rs_factor<<"\t"<<this->Halo[i].local_dm<<"\t"<<this->Halo[i].lambda1<<"\t"<<this->Halo[i].lambda2<<"\t"<<this->Halo[i].lambda3<<"\t"<<log10(this->Halo[i].mass_closest_neighbour)<<"\t"<<this->Halo[i].distance_closest_neighbour<<"\t"<<log10(this->Halo[i].spin_closest_neighbour)<<"\t"<<log10(this->Halo[i].concentration_closest_neighbour)<<"\t"<<log10(this->Halo[i].most_massive_neighbour)<<"\t"<<this->Halo[i].distance_to_most_massive_neighbour<<endl;
      counter++;
    }while(counter<Nobjs_fraction);
    rcat.close();
   So.DONE();
   So.message_screen("\tWritting downsample catalog in file", file);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::select_random_subsample(real_prec fraction, int Nprop, vector<real_prec>&data, string file, string file_pca){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
   ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction*this->NOBJS));
   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(8);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);
   ofstream rcatp;
   rcatp.open(file_pca.c_str());
   rcatp.precision(5);
   rcatp.setf(ios::showpoint);
   rcatp.setf(ios::scientific);
   ULONG counter=0;
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
   while(counter<Nobjs_fraction)
     {
      int i= gsl_rng_uniform_int(gBaseRando,this->Halo.size());
      // write a fraction of original catalog
      rcat<<log10(this->Halo[i].mass)<<"\t"<<log10(this->Halo[i].vmax)<<"\t"<<log10(this->Halo[i].rs)<<"\t"<<log10(this->Halo[i].concentration)<<"\t"<<log10(this->Halo[i].spin_bullock)<<"\t"<<log10(this->Halo[i].vrms)<<"\t"<<this->Halo[i].virial<<"\t"<<this->Halo[i].b_to_a<<"\t"<<this->Halo[i].c_to_a<<"\t"<<this->Halo[i].mach_number<<"\t"<<this->Halo[i].bias<<"\t"<<this->Halo[i].qbias<<"\t"<<"\t"<<log10(1+this->Halo[i].local_overdensity)<<"\t"<<this->Halo[i].dach_number<<"\t"<<this->Halo[i].tidal_anisotropy<<"\t"<<this->Halo[i].peak_height<<"\t"<<this->Halo[i].gal_cwt<<"\t"<<this->Halo[i].relative_bias<<"\t"<<this->Halo[i].bias_rs<<"\t"<<this->Halo[i].rs_factor<<"\t"<<this->Halo[i].local_dm<<"\t"<<this->Halo[i].lambda1<<"\t"<<this->Halo[i].lambda2<<"\t"<<this->Halo[i].lambda3<<"\t"<<log10(this->Halo[i].mass_closest_neighbour)<<"\t"<<this->Halo[i].distance_closest_neighbour<<"\t"<<log10(this->Halo[i].spin_closest_neighbour)<<"\t"<<log10(this->Halo[i].concentration_closest_neighbour)<<"\t"<<log10(this->Halo[i].most_massive_neighbour)<<"\t"<<this->Halo[i].distance_to_most_massive_neighbour<<endl;
      // write a fraction  of pca properties. Only 10 properties
      for(int j=0;j<Nprop;++j)
          rcatp<<data[index_2d(j,i,this->NOBJS)]<<"\t";
      rcatp<<"\t"<<this->Halo[i].gal_cwt<<endl;
      counter++;
   }
   So.DONE();
   So.message_screen("\tDownsampled catalog written in file", file);
   So.message_screen("\tDownsampled PCA catalog with 10 PCA written in file", file_pca);
   rcat.close();
   rcatp.close();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::select_random_subsample_v(real_prec fraction){
    this->So.enter(__PRETTY_FUNCTION__);
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRando;
    gsl_rng_env_setup();
    gsl_rng_default_seed=35;
    rng_t = gsl_rng_mt19937;//_default;
    gBaseRando = gsl_rng_alloc (rng_t);
   ULONG Nobjs_fraction=static_cast<ULONG>(floor(fraction*this->NOBJS));
   So.message_screen("\tSelecting a fraction: ",100*fraction, "%");
   So.message_screen("\tNumber of selected tracers: ",Nobjs_fraction);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->Halo.size();++i)this->Halo[i].observed=0;
   ULONG counter=0;
   while(counter<Nobjs_fraction)
     {
       int i= gsl_rng_uniform_int(gBaseRando,this->Halo.size());
       if(this->Halo[i].observed==0)// not to repeat
         {
          this->Halo[i].observed=1;
          counter++;
        }
    }
   So.DONE();
   So.message_screen("\tDownsampled catalog allocated as observed");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::halos2galaxies_HOD()
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
   for(ULONG i=0;i<this->NOBJS;++i)
    {
      real_prec mass=this->Halo[i].mass;
      real_prec xp= log10(Mc/mass)/sigmaM;
      real_prec prob=gsl_sf_gamma(xp);
      real_prec xr=gsl_rng_uniform(rng);
      if(prob<xr) // select centrals
      {
          xnew.push_back(this->Halo[i].coord1);
          ynew.push_back(this->Halo[i].coord2);
          znew.push_back(this->Halo[i].coord3);
          vxnew.push_back(this->Halo[i].vel1);
          vynew.push_back(this->Halo[i].vel2);
          vznew.push_back(this->Halo[i].vel3);
          massg.push_back(this->Halo[i].mass);
          vmaxg.push_back(this->Halo[i].vmax);
      }
    }
    this->galaxy_central.resize(vmaxg.size());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<vmaxg.size();++i)
    {
     galaxy_central[i].coord1=xnew[i];
     galaxy_central[i].coord2=xnew[i];
     galaxy_central[i].coord3=xnew[i];
     galaxy_central[i].vel1=vxnew[i];
     galaxy_central[i].vel2=ynew[i];
     galaxy_central[i].vel3=vznew[i];
     galaxy_central[i].mass=massg[i];
     galaxy_central[i].vmax=vmaxg[i];
    }
    ULONG Ncentrals=vmaxg.size();
    xnew.clear();xnew.shrink_to_fit();
    ynew.clear();ynew.shrink_to_fit();
    znew.clear();znew.shrink_to_fit();
    vxnew.clear();vxnew.shrink_to_fit();
    vynew.clear();vxnew.shrink_to_fit();
    vznew.clear();vznew.shrink_to_fit();
    massg.clear();massg.shrink_to_fit();
    vmaxg.clear();vmaxg.shrink_to_fit();
    real_prec rmin=0.01; // minium radius used for normalization of NFR profile
    vector<ULONG> Nsats_per_central(Ncentrals,0);
    DensityProfiles DenProf(this->s_cosmo_pars);
            // lopp over the centrals
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<Ncentrals;++i)
          Nsats_per_central[i]=static_cast<ULONG>(gsl_ran_poisson(rng,floor(pow(galaxy_central[i].mass / Mcut , alpha))));  // draw number of satellites per central:
    ULONG counter_sat=0;
    for(ULONG i=0;i<Ncentrals;++i)
        for(ULONG j=0;i<Nsats_per_central[i];++j)
        {
            real_prec phi=2.*M_PI*gsl_rng_uniform(rng);
            real_prec theta=acos(-1.0+2.*gsl_rng_uniform(rng));
            real_prec prob=-10.;
            real_prec radius=0.0;
            real_prec xra=0;
            DenProf.nfw_parameters(galaxy_central[i].mass);
            do {
                radius=DenProf._rvir()*pow(gsl_rng_uniform(rng),1./3.);
                prob=DenProf.DensityProfile_NFW_prob(radius,rmin);
                xra=gsl_rng_uniform(rng);
            }while(prob<xra);
            galaxy_satellite[counter_sat].coord1=radius*sin(theta)*cos(phi)+galaxy_central[i].coord1;
            galaxy_satellite[counter_sat].coord2=radius*cos(theta)*sin(phi)+galaxy_central[i].coord2;
            galaxy_satellite[counter_sat].coord3=radius*cos(phi)+galaxy_central[i].coord3;
            galaxy_satellite[counter_sat].vel1= galaxy_central[i].vel1+gsl_ran_gaussian(rng,galaxy_central[i].vmax) ;
            galaxy_satellite[counter_sat].vel2= galaxy_central[i].vel2+gsl_ran_gaussian(rng,galaxy_central[i].vmax) ;
            galaxy_satellite[counter_sat].vel3= galaxy_central[i].vel3+gsl_ran_gaussian(rng,galaxy_central[i].vmax) ;
            galaxy_satellite[counter_sat].mass= galaxy_central[i].mass ;
            counter_sat++;
        }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 real_prec Catalog::get_min(string prop){
  if(this->Halo.size()==0)
     cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;

  real_prec lka=static_cast<real_prec>(LARGE_NUMBER);
  if(prop=="_MASS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].mass,lka);
    }
  else if(prop=="_VMAX_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
        lka=std::min(this->Halo[i].vmax,lka);
    }
  else if(prop=="_RS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].rs,lka);
    }
  else if(prop=="_RVIR_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].rvir,lka);
    }
  else if(prop==_CONCENTRATION_)
    {
#pragma omp parallel for reduction(min:lka)
    for(ULONG i=0;i<this->Halo.size();++i)
      lka=min(this->Halo[i].concentration,lka);
    }
  else if(prop=="_SPIN_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].spin,lka);
    }
  else if(prop=="_SPIN_BULLOCK_")
    {
#pragma omp parallel for reduction(min:lka)
     for(ULONG i=0;i<this->Halo.size();++i)
        lka=min(this->Halo[i].spin_bullock,lka);
    }
  else if(prop== "_VIRIAL_")
    {
#pragma omp parallel for reduction(min:lka)
    for(ULONG i=0;i<this->Halo.size();++i)
      lka=min(this->Halo[i].virial,lka);
    }
  else if(prop== "_VRMS_")
    {
#pragma omp parallel for reduction(min:lka)
    for(ULONG i=0;i<this->Halo.size();++i)
      lka=min(this->Halo[i].vrms,lka);
    }
  else if(prop=="_BTOA_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].b_to_a,lka);
    }
  else if(prop=="_CTOA_")
    {
#pragma omp parallel for reduction(min:lka)
     for(ULONG i=0;i<this->Halo.size();++i)
       lka=min(this->Halo[i].c_to_a,lka);
    }
  else if(prop=="_MACH_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].mach_number,lka);
    }
  else if(prop=="_DACH_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].dach_number,lka);
    }
  else if(prop=="_BIAS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].bias,lka);
    }
  else if(prop=="_RELATIVE_BIAS_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].relative_bias,lka);
    }
  else if(prop=="_LOCAL_OVERDENSITY_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].local_overdensity,lka);
    }
  else if(prop=="_TIDAL_ANISOTROPY_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].tidal_anisotropy,lka);
    }
  else if(prop=="_LOCALDM_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].local_dm,lka);
    }
  else if(prop=="_PEAK_HEIGHT_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].peak_height,lka);
    }
  else if(prop=="_XCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].coord1,lka);
    }
  else if(prop=="_YCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].coord2,lka);
    }
  else if(prop=="_ZCOORD_")
    {
#pragma omp parallel for reduction(min:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=min(this->Halo[i].coord3,lka);
    }
  else{
  this->So.message_screen("Property not found");
  }
  return lka;
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 real_prec Catalog::get_max(string prop){
 if(this->Halo.size()==0)
    cerr<<"Error. Array Halo is empty, in line "<<__PRETTY_FUNCTION__<<endl;
  real_prec lka=-static_cast<real_prec>(LARGE_NUMBER);
  real_prec lkb;
  if(prop=="_MASS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].mass,lka);
    }
  else if(prop=="_VMAX_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].vmax,lka);
    }
  else if(prop=="_RS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].rs,lka);
    }
  else if(prop=="_RVIR_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].rvir,lka);
    }
  else if(prop=="_CONCENTRATION_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].concentration,lka);
    }
  else if(prop=="_SPIN_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].spin,lka);
    }
  else if(prop=="_SPIN_BULLOCK_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].spin_bullock,lka);
    }
  else if(prop== "_VIRIAL_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].virial,lka);
    }
  else if(prop== "_VRMS_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].vrms,lka);
    }
  else if(prop=="_BTOA_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].b_to_a,lka);
    }
   else if(prop=="_CTOA_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].c_to_a,lka);
     }
   else if(prop=="_MACH_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].mach_number,lka);
     }
  else if(prop=="_DACH_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
    lka=max(this->Halo[i].dach_number,lka);
    }
   else if(prop=="_BIAS_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].bias,lka);
     }
   else if(prop=="_RELATIVE_BIAS_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].relative_bias,lka);
     }
   else if(prop=="_LOCAL_OVERDENSITY_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].local_overdensity,lka);
     }
   else if(prop=="_TIDAL_ANISOTROPY_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
       lka=max(this->Halo[i].tidal_anisotropy,lka);
     }
  else if(prop=="_LOCALDM_")
    {
#pragma omp parallel for reduction(max:lka)
      for(ULONG i=0;i<this->Halo.size();++i)
      lka=max(this->Halo[i].local_dm,lka);
    }
   else if(prop=="_PEAK_HEIGHT_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].peak_height,lka);
     }
   else if(prop=="_XCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].coord1,lka);
     }
   else if(prop=="_YCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].coord2,lka);
     }
   else if(prop=="_ZCOORD_")
     {
#pragma omp parallel for reduction(max:lka)
       for(ULONG i=0;i<this->Halo.size();++i)
     lka=max(this->Halo[i].coord3,lka);
     }
   else
     this->So.message_screen("Property not found");

  return lka;
 }
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  real_prec Catalog::get_variance(real_prec meanp, string prop){
  if(this->Halo.size()==0)
     cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;
   real_prec var=0;
  if(prop=="_VMAX_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(this->Halo[i].vmax-meanp,2);
   else if(prop=="_MASS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(log10(this->Halo[i].mass)-meanp,2);
  else if(prop=="_RS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(this->Halo[i].rs-meanp,2);
  else if(prop=="_RVIR_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(this->Halo[i].rvir-meanp,2);
  else if(prop=="_CONCENTRATION_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(log10(this->Halo[i].concentration)-meanp,2);
  else if(prop=="_SPIN_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(this->Halo[i].spin-meanp,2);
  else if(prop=="_SPIN_BULLOCK")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
   for(ULONG i=0;i<this->Halo.size();++i)
       var+=pow(this->Halo[i].spin_bullock-meanp,2);
   else if(prop== "_VIRIAL_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(this->Halo[i].virial-meanp,2);
   else if(prop== "_VRMS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(this->Halo[i].vrms-meanp,2);
   else if(prop=="_BTOA_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(this->Halo[i].b_to_a-meanp,2);
   else if(prop=="_CTOA_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(this->Halo[i].c_to_a-meanp,2);
   else if(prop=="_MACH_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
        var+=pow(this->Halo[i].mach_number-meanp,2);
   else if(prop=="_BIAS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
        var+=pow(this->Halo[i].bias-meanp,2);
   else if(prop=="_RELATIVE_BIAS_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
        var+=pow(this->Halo[i].relative_bias-meanp,2);
   else if(prop=="_LOCAL_OVERDENSITY_")
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
        var+=pow(this->Halo[i].local_overdensity-meanp,2);
  else if(prop=="_TIDAL_ANISOTROPY_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
      var+=pow(this->Halo[i].tidal_anisotropy-meanp,2);
  }
  else if(prop=="_LOCALDM_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
      var+=pow(this->Halo[i].local_dm-meanp,2);
  }
    else if(prop=="_PEAK_HEIGHT_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
    for(ULONG i=0;i<this->Halo.size();++i)
        var+=pow(this->Halo[i].peak_height-meanp,2);
    var/=static_cast<real_prec>(this->Halo.size());
  }
  else if(prop=="_DACH_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
  for(ULONG i=0;i<this->Halo.size();++i)
      var+=pow(this->Halo[i].dach_number-meanp,2);
  var/=static_cast<real_prec>(this->Halo.size());
}
  else
      this->So.message_screen("Property not found");
  return var;
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 pair<real_prec,real_prec> Catalog::get_variance(string prop,bool factor){
     real_prec fac=0;
     if(factor==true)
         fac=1;
   if(this->Halo.size()==0)
      cerr<<"Error. Empty vector in function "<<__PRETTY_FUNCTION__<<endl;
    real_prec var=0;
    real_prec meanp=0;
   if(prop=="_VMAX_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
         meanp+=log10(this->Halo[i].vmax);
      meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(log10(this->Halo[i].vmax)-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
    }
    else if(prop=="_MASS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
         meanp+=log10(this->Halo[i].mass);
      meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
          var+=pow(log10(this->Halo[i].mass)-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
    }
   else if(prop=="_RS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=log10(this->Halo[i].rs);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].rs)-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_RVIR_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=log10(this->Halo[i].rvir);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].rvir)-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_CONCENTRATION_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(log10(this->Halo[i].concentration));
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].concentration)-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_SPIN_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=log10(this->Halo[i].spin);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].spin)-fac*meanp,2);
     var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_SPIN_BULLOCK_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=log10(this->Halo[i].spin_bullock);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].spin_bullock)-fac*meanp,2);
     var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_VRMS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=log10(this->Halo[i].vrms);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(log10(this->Halo[i].vrms)-fac*meanp,2);
     var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_VIRIAL_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].virial);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].virial-fac*meanp,2);
      var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_BTOA_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].b_to_a);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].b_to_a-fac*meanp,2);
     var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_CTOA_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].c_to_a);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].c_to_a-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_MACH_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].mach_number);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].mach_number-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_BIAS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].bias);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].bias-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_RELATIVE_BIAS_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].relative_bias);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].relative_bias-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_LOCAL_OVERDENSITY_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].local_overdensity);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].local_overdensity-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_TIDAL_ANISOTROPY_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].tidal_anisotropy);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].tidal_anisotropy-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_LOCALDM_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].local_dm);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].local_dm-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_PEAK_HEIGHT_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
        meanp+=(this->Halo[i].peak_height);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].peak_height-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
   else if(prop=="_DACH_")
   {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:meanp)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)

        meanp+=(this->Halo[i].dach_number);
     meanp/=static_cast<real_prec>(this->Halo.size());
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(_NTHREADS_) reduction(+:var)
#endif
     for(ULONG i=0;i<this->Halo.size();++i)
         var+=pow(this->Halo[i].dach_number-fac*meanp,2);
       var/=static_cast<real_prec>(this->Halo.size());
   }
    else
   this->So.message_screen("Property not found: ", prop);

   return std::make_pair(meanp,var);
   }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_intervals_equal_number(string prop,vector<real_prec>&min_aux,vector<real_prec>&max_aux){
  ULONG Ntracers=this->NOBJS;
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
   ULONG Ntracers_bin=static_cast<ULONG>(floor(Ntracers/Nbins_prop));
   So.message_screen("\tExpected (equal) number of tracers in bins of sec property =", Ntracers_bin);
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_PEAK_HEIGHT_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_ph_bins[0]=this->NOBJS;
     this->Number_of_tracers_in_ph_bins[Nbins_prop]=this->NOBJS-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_ph_bins[i]=Ntracers_bin;
   }
#elif defined  _USE_MASS_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_MASS_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_mass_bins[0]=this->NOBJS;
     this->Number_of_tracers_in_mass_bins[Nbins_prop]=this->NOBJS-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_mass_bins[i]=Ntracers_bin;
   }
#elif defined  _USE_VMAX_AS_PRIMARY_PROPERTY_HBIAS_
   if(prop=="_VMAX_")// ojo aca que quiza hace falta PH y VMAX tratado asi
   {
     this->Number_of_tracers_in_vmax_bins[0]=this->NOBJS;
     this->Number_of_tracers_in_vmax_bins[Nbins_prop]=this->NOBJS-(Nbins_prop-1)*Ntracers_bin;
     for(int i=1;i<Nbins_prop;++i)
       this->Number_of_tracers_in_vmax_bins[i]=Ntracers_bin;
   }
#endif
   ULONG NBins=100000;
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
   vector<ULONG>counts(NBins,0);
   if(prop=="_VMAX_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].vmax,minp,NBins,delta, true)]++;
}
   else if(prop=="_RS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].rs,minp,NBins,delta, true)]++;
  }
   else if(prop=="_RVIR_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].rvir,minp,NBins,delta, true)]++;
  }
   else if(prop=="_CONCENTRATION_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].concentration,minp,NBins,delta, true)]++;
  }
   else if(prop=="_SPIN_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].spin,minp,NBins,delta, true)]++;
}
   else if(prop=="_SPIN_BULLOCK_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].spin_bullock,minp,NBins,delta, true)]++;
}
     else if(prop=="_VRMS_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].vrms,minp,NBins,delta, true)]++;
 }

     else if(prop=="_VIRIAL_")
 {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[ get_bin(this->Halo[i].virial,minp,NBins,delta, true)]++;
  }
     else if(prop=="_BTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].b_to_a,minp,NBins,delta,true)]++;
}
     else if(prop=="_CTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].c_to_a,minp,NBins,delta, true)]++;
}
     else if(prop=="_MACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].mach_number,minp,NBins,delta, true)]++;
}
     else if(prop=="_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].bias,minp,NBins,delta, true)]++;
}
     else if(prop=="_RELATIVE_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].relative_bias, minp,NBins,delta, true)]++;
}

     else if(prop=="_LOCAL_OVERDENSITY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].local_overdensity,minp,NBins,delta, true)]++;
}
     else if(prop=="_PEAK_HEIGHT_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].peak_height,minp,NBins,delta, true)]++;
}
   else if(prop=="_DACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].dach_number,minp,NBins,delta, true)]++;
}
   else if(prop=="_TIDAL_ANISOTROPY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].tidal_anisotropy,minp,NBins,delta, true)]++;
}
   else if(prop=="_LOCALDM_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].local_dm,minp,NBins,delta, true)]++;
}
     else if(prop=="_MASS_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].mass,minp,NBins,delta, true)]++;
    }
   min_aux[0]=minp;
   max_aux[0]=maxp;
   min_aux[1]=minp;
   max_aux[Nbins_prop]=maxp;
   ULONG indi=0;
   for(int j=1;j<Nbins_prop;++j)// Go through the 4 bins for quartiles
     {
       real_prec new_counts=0;
       ULONG i=indi;
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
 void Catalog::get_intervals_equal_number_aux(string prop){
   this->So.enter(__PRETTY_FUNCTION__);
   int Nbins_prop= this->params._Number_of_bins_equal_number_tracers();
   ULONG NBins=50000;
   ULONG Ntracers=this->NOBJS;
   real_prec minp=get_min(prop);
   real_prec maxp=get_max(prop);
   So.message_screen("Type of object:", this->type_of_object);
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
   vector<ULONG>counts(NBins,0);
   if(prop=="_VMAX_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].vmax,minp,NBins,delta, true)]++;
    }
   else if(prop=="_RS_")
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].rs,minp,NBins,delta, true)]++;
    }
  else if(prop=="_RVIR_")
   {
   #ifdef _USE_OMP_
   #pragma omp parallel for
   #endif
      for(ULONG i=0;i<Ntracers;++i)
   #ifdef _USE_OMP_
   #pragma atomic
   #endif
        counts[get_bin(this->Halo[i].rvir,minp,NBins,delta, true)]++;
   }
   else if(prop==_CONCENTRATION_)
   {
   #ifdef _USE_OMP_
   #pragma omp parallel for
   #endif
      for(ULONG i=0;i<Ntracers;++i)
   #ifdef _USE_OMP_
   #pragma atomic
   #endif
        counts[get_bin(this->Halo[i].concentration,minp,NBins,delta, true)]++;
   }
  else if(prop=="_SPIN_") 
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].spin,minp,NBins,delta, true)]++;
  }
   else if(prop=="_SPIN_BULLOCK_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
 for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
   counts[get_bin(this->Halo[i].spin_bullock,minp,NBins,delta, true)]++;
  }
  else if(prop=="_VRMS_")
  {
    for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
      counts[get_bin(this->Halo[i].vrms,minp,NBins,delta, true)]++;
  }       
  else if(prop=="_VIRIAL_")
  {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[ get_bin(this->Halo[i].virial,minp,NBins,delta, true)]++;
 }
     else if(prop=="_BTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].b_to_a,minp,NBins,delta,true)]++;
 }
     else if(prop=="_CTOA_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].c_to_a,minp,NBins,delta, true)]++;
 }
     else if(prop=="_MACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].mach_number,minp,NBins,delta, true)]++;
 }
     else if(prop=="_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].bias,minp,NBins,delta, true)]++;
 }
      else if(prop=="_RELATIVE_BIAS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].relative_bias,minp,NBins,delta, true)]++;
 }
     else if(prop=="_LOCAL_OVERDENSITY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].local_overdensity,minp,NBins,delta, true)]++;
 }
     }
   else if(prop=="_TIDAL_ANISOTROPY_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].tidal_anisotropy,minp,NBins,delta, true)]++;
     }
}

   else if(prop=="_LOCALDM_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].local_dm,minp,NBins,delta, true)]++;
     }
}
   else if(prop=="_PEAK_HEIGHT_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i){
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].peak_height,minp,NBins,delta, true)]++;
     }
 }
   else if(prop=="_DACH_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
     counts[get_bin(this->Halo[i].dach_number,minp,NBins,delta, true)]++;
}

     else if(prop=="_MASS_")
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<Ntracers;++i)
#ifdef _USE_OMP_
#pragma atomic
#endif
       counts[get_bin(this->Halo[i].mass,minp,NBins,delta, true)]++;
 }
  // Get the intervals
   ULONG Ntracers_bin=static_cast<ULONG>(floor(Ntracers/Nbins_prop));

   vector<real_prec>min_aux(Nbins_prop+1,0);
   vector<real_prec>max_aux(Nbins_prop+1,0);
   min_aux[0]=minp;
   max_aux[0]=maxp;
   min_aux[1]=minp;
   max_aux[Nbins_prop]=maxp;
   ULONG indi=0;
   for(int j=1;j<Nbins_prop;++j)
     {
       real_prec new_counts=0;
       ULONG i=indi;
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
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_MASSbins_min(i,log10(min_aux[i]));
       this->params.set_MASSbins_max(i,log10(max_aux[i]));
     }
   else if(prop=="_VMAX_")
     for(int i=0; i<Nbins_prop+1;++i)
       {
     this->params.set_VMAXbins_min(i,log10(min_aux[i]));
     this->params.set_VMAXbins_max(i,log10(max_aux[i]));
       }
    else if(prop=="_RS_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_RSbins_min(i,min_aux[i]);
       this->params.set_RSbins_max(i,max_aux[i]);
     }
   else if(prop=="_RVIR_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_RVIRbins_min(i,min_aux[i]);
       this->params.set_RVIRbins_max(i,max_aux[i]);
     }
   else if(prop=="_CONCENTRATION_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_CONCENTRATIONbins_min(i,min_aux[i]);
       this->params.set_CONCENTRATIONbins_max(i,max_aux[i]);
     }
  else  if(prop=="_SPIN_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_SPINbins_min(i,min_aux[i]);
       this->params.set_SPINbins_max(i,max_aux[i]);
     }
    else if(prop=="_SPIN_BULLOCK_")
       for(int i=0; i<Nbins_prop+1;++i){
         this->params.set_SPINBULLOCKbins_min(i,min_aux[i]);
         this->params.set_SPINBULLOCKbins_max(i,max_aux[i]);
       }
    else if(prop=="_VIRIAL_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_VIRIALbins_min(i,min_aux[i]);
       this->params.set_VIRIALbins_max(i,max_aux[i]);
     }
   if(prop=="_VRMS_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_VRMSbins_min(i,min_aux[i]);
       this->params.set_VRMSbins_max(i,max_aux[i]);
     }
   if(prop=="_BTOA_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_BTOAbins_min(i,min_aux[i]);
       this->params.set_BTOAbins_max(i,max_aux[i]);
     }

   if(prop=="_CTOA_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_CTOAbins_min(i,min_aux[i]);
       this->params.set_CTOAbins_max(i,max_aux[i]);
     }

   if(prop=="_MACH_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_MACHbins_min(i,min_aux[i]);
       this->params.set_MACHbins_max(i,max_aux[i]);
     }
   if(prop=="_BIAS_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_BIASbins_min(i,min_aux[i]);
       this->params.set_BIASbins_max(i,max_aux[i]);
     }
   if(prop=="_RELATIVE_BIAS_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_RBIASbins_min(i,min_aux[i]);
       this->params.set_RBIASbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCAL_OVERDENSITY_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_LCbins_min(i,min_aux[i]);
       this->params.set_LCbins_max(i,max_aux[i]);
     }
   if(prop=="_TIDAL_ANISOTROPY_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_TAbins_min(i,min_aux[i]);
       this->params.set_TAbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCALDM_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_LOCALDMbins_min(i,min_aux[i]);
       this->params.set_LOCALDMbins_max(i,max_aux[i]);
     }
   if(prop=="_PEAK_HEIGHT_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_PHbins_min(i,min_aux[i]);
       this->params.set_PHbins_max(i,max_aux[i]);
     }
   if(prop=="_DACH_")
     for(int i=0; i<Nbins_prop+1;++i){
       this->params.set_DACHbins_min(i,min_aux[i]);
       this->params.set_DACHbins_max(i,max_aux[i]);
     }
   if(prop=="_LOCALDM_")
     for(int i=0; i<Nbins_prop+1;++i){
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
void Catalog::get_local_mach_number(bool write)
{
  So.enter(__PRETTY_FUNCTION__);
  int nprops=7;
  vector<proper>hprop(nprops);
  int IC=0;
  hprop[IC].prop_name="COUNTS";
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
  if(this->params._i_mass_g()>0)
    {
      IC++;
      hprop[IC].prop_name="MASS";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].mass*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  if(this->params._i_vmax_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VMAX";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].vmax*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  if(this->params._i_vrms_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VRMS";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].vrms*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  if(this->params._i_virial_g()>0)
    {
      IC++;
      hprop[IC].prop_name="VIRIAL";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].virial*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  if(this->params._i_spin_g()>0)
    {
      IC++;
      hprop[IC].prop_name="SPIN";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].spin*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  if(this->params._i_rs_g()>0)
    {
      IC++;
      hprop[IC].prop_name="RS";
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    hprop[IC].prop.push_back(this->Halo[i].rs*sqrt(pow(this->Halo[i].vel1,2)+pow(this->Halo[i].vel2,2)+pow(this->Halo[i].vel3,2)));
    }

  for(int iprop=0;iprop<nprops;++iprop)
    {
      vector<s_cell_info> cell_info_tr(this->box.NGRID);
      for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    cell_info_tr[this->Halo[i].GridID].property.push_back(hprop[iprop].prop[i]);

  So.DONE();

  this->mach_number.resize(this->box.NGRID,0);
  for(ULONG i=0;i<this->box.NGRID ;++i) //loop over the "observed obejcts", i.e, with cuts already set
  {
    real_prec mean=0;
    real_prec var=0;
    for(ULONG j=0;j< cell_info_tr[i].property.size();++j)
      mean+=cell_info_tr[i].property[j];
    mean/=static_cast<double>(cell_info_tr[i].property.size());
    for(ULONG j=0;j<cell_info_tr[i].property.size();++j)
        var+=pow(cell_info_tr[i].property[j]-mean,2);
    var/=static_cast<double>(cell_info_tr[i].property.size());
    if(var>0)
        this->mach_number[i]=mean/sqrt(var);
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
void Catalog::get_local_mach_number(real_prec scale)
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
    for(ULONG i=0;i<this->params._NGRID() ;++i)
    {
        for(ULONG j=0;j<nearest_cells_to_cell[i].close_cell.size();++j)
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

   for(ULONG i=0;i<this->params._NGRID() ;++i)
   {
      for(ULONG j=0;j<nearest_cells_to_cell[i].close_cell.size();++j)
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
   for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID;
//      cell_info_tr[ID].posx_p.push_back(this->Halo[i].coord1);// esto se puede evitar escribiendo simplemente el indice y llamando al tracer.Halo en ese índice
//      cell_info_tr[ID].posy_p.push_back(this->Halo[i].coord2);
 //     cell_info_tr[ID].posz_p.push_back(this->Halo[i].coord3);
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  ULONG Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(ULONG i=0;i<this->params._NGRID() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    Nt_check+=cell_info_tr[i].gal_index.size();
// expected number of tracer in a sphere of radius Scale
  real_prec Nbar=(4./3.)*M_PI*pow(scale/this->params._Lbox(), 3)*this->NOBJS;
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
  for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      real_prec x_coord=this->Halo[i].coord1;
      real_prec y_coord=this->Halo[i].coord2;
      real_prec z_coord=this->Halo[i].coord3;
      ULONG ID=this->Halo[i].GridID;
      real_prec meanv=0; // variable used to compute mach number per each tracer
      ULONG counter_gal=0;// number of tracers in the sphere of radius scale around the current tracer
      real_prec meand=0; // variable used to compute mach number per each tracer
      // Get the mean vel of tracers in a sphere of radii scale around the current tracer
      for(int j=0;j<N_Neigh_cells; j++) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
        {
          ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
          real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->params._Lbox();
          real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->params._Lbox();
          real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->params._Lbox();
          for(int k=0; k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
            {
              ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
              real_prec distx = x_coord - (this->Halo[index_gal].coord1+factor_bc_x);
              real_prec disty = y_coord - (this->Halo[index_gal].coord2+factor_bc_y);
              real_prec distz = z_coord - (this->Halo[index_gal].coord3+factor_bc_z);
              real_prec dist  = distx*distx + disty*disty+ distz*distz;
              real_prec vx=this->Halo[index_gal].vel1;
              real_prec vy=this->Halo[index_gal].vel2;
              real_prec vz=this->Halo[index_gal].vel3;
#ifdef relative_mach
              vx=this->Halo[i].vel1-vx;
              vy=this->Halo[i].vel2-vy;
              vz=this->Halo[i].vel3-vz;
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
                  ULONG ID_NEIGH = nearest_cells_to_cell[ID].close_cell[j];
                  real_prec factor_bc_x=nearest_cells_to_cell[ID].bc_x[j]*this->params._Lbox();
                  real_prec factor_bc_y=nearest_cells_to_cell[ID].bc_y[j]*this->params._Lbox();
                  real_prec factor_bc_z=nearest_cells_to_cell[ID].bc_z[j]*this->params._Lbox();
                  for(int k=0;k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
                    {
                      ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                      real_prec distx = x_coord - (this->Halo[index_gal].coord1+factor_bc_x);
                      real_prec disty = y_coord - (this->Halo[index_gal].coord2+factor_bc_y);
                      real_prec distz = z_coord - (this->Halo[index_gal].coord3+factor_bc_z);
                      real_prec dist  = distx*distx + disty*disty+ distz*distz;
                      real_prec vx=this->Halo[index_gal].vel1;
                      real_prec vy=this->Halo[index_gal].vel2;
                      real_prec vz=this->Halo[index_gal].vel3;
#ifdef relative_mach
                      vx=this->Halo[i].vel1-vx;
                      vy=this->Halo[i].vel2-vy;
                      vz=this->Halo[i].vel3-vz;
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
            this->Halo[i].mach_number=meanv/sqrt(varv);
            this->Halo[i].dach_number=meand/sqrt(vard);
            this->Halo[i].local_overdensity=(static_cast<real_prec>(counter_gal)-Nbar)/Nbar ;//meand/sqrt(vard);
            this->Halo[i].number_of_neighbours=counter_gal;

         }
         else
         {
           this->Halo[i].mach_number=0.;
           this->Halo[i].dach_number=0.;
           this->Halo[i].local_overdensity=0.;
           this->Halo[i].number_of_neighbours=0;
         }
     }
  nearest_cells_to_cell.clear();nearest_cells_to_cell.shrink_to_fit();
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_local_mach_number_chuncks(real_prec scale)
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
  real_prec Nbar=(4./3.)*M_PI*pow(scale/this->params._Lbox(), 3)*this->NOBJS;
  So.message_screen("\tMean number of tracers in a sphere of radius R :", Nbar);
  int N_slices=static_cast<int>(this->params._Nft_low()/chunck_factor);
#ifdef _VERBOSE_CAT_
  So.message_screen("\tIdentifying tracers in cells:");
#endif
  // This vector structures contains the coords of the tracers in each cell identified with a 3D ID
  vector<s_cell_info_reduced> cell_info_tr(this->params._NGRID_low());
  // This vector structures contains the coords of the tracers in each of the N_slice slices.
  vector<s_cell_info_reduced> cell_info_tr_slice(N_slices);

  for(ULONG i=0;i<this->Halo.size() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    {
      ULONG ID=this->Halo[i].GridID_n;
      cell_info_tr[ID].gal_index.push_back(i);//allocate the index of each tracer in that cell
      ULONG xbin=get_bin(this->Halo[i].coord1,0,N_slices, delta*chunck_factor,false);
      cell_info_tr_slice[xbin].gal_index.push_back(i);//allocate the index of each tracer in that cell
    }
  So.DONE();
  ULONG Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(ULONG i=0;i<this->params._NGRID_low() ;++i) //loop over the "observed obejcts", i.e, with cuts already set
    Nt_check+=cell_info_tr[i].gal_index.size();
  So.message_screen("\tChecking: Number of tracers in structure: ", Nt_check);
  Nt_check=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:Nt_check)
#endif
  for(ULONG i=0;i< N_slices ;++i) //loop over the "observed obejcts", i.e, with cuts already set
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
        for(ULONG i=0;i<cell_info_tr_slice[icc].gal_index.size() ;++i) //loop over the "observed obejcts in the slice"
          {
          ULONG index_gal_slice=cell_info_tr_slice[icc].gal_index[i];
          real_prec x_coord=this->Halo[index_gal_slice].coord1;
          real_prec y_coord=this->Halo[index_gal_slice].coord2;
          real_prec z_coord=this->Halo[index_gal_slice].coord3;
          // This is the x-bin (e.g., 0 or 1 for chunck_factor=2) for each super/slice
          int xbinS=get_bin(x_coord, icc*delta*chunck_factor,  chunck_factor, delta,false);
          int ybin=get_bin(y_coord,0,this->params._Nft_low(),delta,false); // normal y-bin
          int zbin=get_bin(z_coord,0,this->params._Nft_low(),delta,false); // normal z-bin
          ULONG ID_slice=index_3d(xbinS,ybin,zbin,this->params._Nft_low(),this->params._Nft_low());
          double meanv=0; // variable used to compute mach number per each tracer
          ULONG counter_gal=0;// number of tracers in the sphere of radius scale around the current tracer
          double meand=0; // variable used to compute dach number per each tracer
          // Get the mean vel of tracers in a sphere of radii scale around the current tracer
          for(int j=0;j< N_Neigh_cells; ++j) // loop sobre las celdas cercanas incluida la celda donde esta la particula i
            {
              ULONG ID_NEIGH = nearest_cells_to_cell[ID_slice].close_cell[j];
              real_prec factor_bc_x=nearest_cells_to_cell[ID_slice].bc_x[j]*this->params._Lbox();
              real_prec factor_bc_y=nearest_cells_to_cell[ID_slice].bc_y[j]*this->params._Lbox();
              real_prec factor_bc_z=nearest_cells_to_cell[ID_slice].bc_z[j]*this->params._Lbox();
              vector<real_prec>mass_n;
              vector<real_prec>dist_n;
              vector<real_prec>spin_n;
              vector<real_prec>concentration_n;
              vector<ULONG>id_n;
              for(int k=0; k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
                {
                  ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                  real_prec distx = x_coord - (this->Halo[index_gal].coord1+factor_bc_x);
                  real_prec disty = y_coord - (this->Halo[index_gal].coord2+factor_bc_y);
                  real_prec distz = z_coord - (this->Halo[index_gal].coord3+factor_bc_z);
                  
//                  cout<<index_gal_slice<<"  "<<distx<<"  "<<x_coord<<"  "<<this->Halo[index_gal].coord1<<"  "<<factor_bc_x<<"  "<<scale<<endl;
   //               cout<<index_gal<<"   "<<this->Halo[index_gal].GridID_n<<"  "<<factor_bc_x<<"  "<<scale<<endl;
                  
                  real_prec dist  = distx*distx+ disty*disty + distz*distz;
                  real_prec vx=this->Halo[index_gal].vel1;
                  real_prec vy=this->Halo[index_gal].vel2;
                  real_prec vz=this->Halo[index_gal].vel3;
#ifdef relative_mach
                  vx-=this->Halo[index_gal_slice].vel1;
                  vx-=this->Halo[index_gal_slice].vel2;
                  vx-=this->Halo[index_gal_slice].vel3;
#endif
                  if(dist<=scale)
                  {
                    meanv+=static_cast<double>(sqrt(vx*vx+vy*vy+vz*vz));
                    meand+=sqrt(dist);
                    counter_gal++;
                    mass_n.push_back(this->Halo[index_gal].mass);
                    dist_n.push_back(sqrt(dist));
                    spin_n.push_back(this->Halo[index_gal].spin);
                    id_n.push_back(index_gal);
                    concentration_n.push_back(this->Halo[index_gal].concentration);
                  }
                }
                if(mass_n.size()>0)
                 {
                    ULONG idmin_m, idmax_m;
                    ULONG idmin_d, idmax_d;
                    // Get the id of the minimum and maximum mass
                    min_max_vector(mass_n, idmin_m,idmax_m);
                    // Get the id of the minimum and maximum separation
                    min_max_vector(dist_n, idmin_d,idmax_d);
                    this->Halo[index_gal_slice].most_massive_neighbour=this->Halo[id_n[idmax_m]].mass;
                    this->Halo[index_gal_slice].distance_closest_neighbour=dist_n[idmin_d];
                    // Get the distance to the most massive neighbour
                    this->Halo[index_gal_slice].distance_to_most_massive_neighbour=dist_n[idmax_m];
                    this->Halo[index_gal_slice].mass_closest_neighbour=this->Halo[id_n[idmin_d]].mass;
                    this->Halo[index_gal_slice].spin_closest_neighbour=this->Halo[id_n[idmin_d]].spin;
                    this->Halo[index_gal_slice].concentration_closest_neighbour=this->Halo[id_n[idmin_d]].concentration;
                 }
               else
                {
                  this->Halo[index_gal_slice].most_massive_neighbour=0;
                  this->Halo[index_gal_slice].distance_closest_neighbour=0;
                  this->Halo[index_gal_slice].distance_to_most_massive_neighbour=0;
                  this->Halo[index_gal_slice].mass_closest_neighbour=0;
                  this->Halo[index_gal_slice].spin_closest_neighbour=0;
                  this->Halo[index_gal_slice].concentration_closest_neighbour=0;
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
                ULONG ID_NEIGH = nearest_cells_to_cell[ID_slice].close_cell[j];
                real_prec factor_bc_x=nearest_cells_to_cell[ID_slice].bc_x[j]*this->params._Lbox();
                real_prec factor_bc_y=nearest_cells_to_cell[ID_slice].bc_y[j]*this->params._Lbox();
                real_prec factor_bc_z=nearest_cells_to_cell[ID_slice].bc_z[j]*this->params._Lbox();
                for(int k=0;k< cell_info_tr[ID_NEIGH].gal_index.size(); ++k) //loop sobre las particulas en las celdas cercanas
              {
                ULONG index_gal=cell_info_tr[ID_NEIGH].gal_index[k];
                real_prec distx = x_coord - (this->Halo[index_gal].coord1+factor_bc_x);
                real_prec disty = y_coord - (this->Halo[index_gal].coord2+factor_bc_y);
                real_prec distz = z_coord - (this->Halo[index_gal].coord3+factor_bc_z);
                real_prec dist  = distx*distx+ disty*disty + distz*distz;
                real_prec vx=this->Halo[index_gal].vel1;
                real_prec vy=this->Halo[index_gal].vel2;
                real_prec vz=this->Halo[index_gal].vel3;
    #ifdef relative_mach
                vx-=this->Halo[index_gal_slice].vel1;
                vy-=this->Halo[index_gal_slice].vel2;
                vz-=this->Halo[index_gal_slice].vel3;
    #endif
                if(dist<=scale){
                  varv+=static_cast<double>(pow(sqrt(vx*vx+vy*vy+vz*vz)-meanv,2));
                  vard+=static_cast<double>(pow(sqrt(dist)-meand,2));
                }
              }
            }
            varv/=static_cast<real_prec>(counter_gal);
            vard/=static_cast<real_prec>(counter_gal);
            this->Halo[index_gal_slice].mach_number= varv ==0? 0.: meanv/sqrt(varv);
            this->Halo[index_gal_slice].dach_number= vard==0? 0.:  meand/sqrt(vard);
            this->Halo[index_gal_slice].local_overdensity=(static_cast<real_prec>(counter_gal)-Nbar)/Nbar;//meand/sqrt(vard);
            this->Halo[index_gal_slice].number_of_neighbours=counter_gal;
          }
          else
          {
            this->Halo[index_gal_slice].mach_number=0.;
            this->Halo[index_gal_slice].dach_number=0.;
            this->Halo[index_gal_slice].local_overdensity=0.;
            this->Halo[index_gal_slice].number_of_neighbours=0;
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
void Catalog::get_tracer_tidal_anisotropy(vector<real_prec>&tidal){
    So.message_screen("Computing tidal anisotropy at halo position");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<this->NOBJS ;++i) //loop over the "observed obejcts", i.e, with cuts already set
//     this->Halo[i].tidal_anisotropy = linInterpol(this->params._Nft(),this->params._Lbox(), this->params._d1(),this->Halo[i].coord1,this->Halo[i].coord2,this->Halo[i].coord3,tidal);
   this->Halo[i].tidal_anisotropy =tidal[this->Halo[i].GridID];
     So.DONE();
  this->min_tidal_anisotropy=this->get_min("_TIDAL_ANISOTROPY_");
  this->max_tidal_anisotropy=this->get_max("_TIDAL_ANISOTROPY_");
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog:: get_peak_height_at_tracer(){
    So.message_screen("Computing peak height for tracers");
    this->ps.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    real_prec gr=this->Cosmo.growth_factor(this->params._redshift());
    real_prec gr_zero=this->Cosmo.growth_factor(0);
    this->s_cosmo_pars.growth_factor=gr/gr_zero;
    this->s_cosmo_pars.pk_normalization=ps.normalization();
    this->Cosmo.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    this->ps.set_cosmo_pars(this->s_cosmo_pars); // update s_cosmo
    So.message_screen("        Normalization of power =",this->s_cosmo_pars.pk_normalization);
     // ******************************************************************************
     int Nbinsm=300;
     vector<gsl_real>masses(Nbinsm,0);
     vector<gsl_real>nus(Nbinsm,0);
     for(ULONG i=0;i< masses.size();++i)
     {
         masses[i]=log10(1e11)+(i+0.5)*log10(7e15/1e11)/static_cast<real_prec>(Nbinsm);
         nus[i]=this->stats.peak_height(masses[i], this->params._redshift(), &this->s_cosmo_pars);//
     }
     gsl_interp_accel *acc = gsl_interp_accel_alloc();
     gsl_spline *spline = gsl_spline_alloc(gsl_interp_linear,masses.size());
     gsl_spline_init(spline,&masses[0],&nus[0],masses.size());
#ifdef _USE_OMP_
#pragma omp  parallel for
#endif
     for(ULONG i=0;i< this->NOBJS;++i)
       this->Halo[i].peak_height=log10(gsl_spline_eval(spline,log10(this->Halo[i].mass), acc));
#ifdef _USE_PEAK_HEIGHT_AS_PRIMARY_PROPERTY_HBIAS_
      this->min_peak_height=get_min("_PEAK_HEIGHT_");
      this->max_peak_height=get_max("_PEAK_HEIGHT_");
      this->So.message_screen("\tMin NU  =", this->min_peak_height);
      this->So.message_screen("\tMax NU  =", this->max_peak_height);
      So.DONE();
      // This is in case the main property is biined in equal number bins, which is not often the case.
      if(true==this->params._set_bins_equal_number_tracers_main_property()) // bins in mass exit in log M to be consistent with the case in wich fixed bins are used and read from input file
       {
       vector<real_prec>min_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       vector<real_prec>max_aux(this->params._Number_of_bins_equal_number_tracers_main_property(),0);
       this->Number_of_tracers_in_ph_bins.resize(this->params._Number_of_bins_equal_number_tracers_main_property()+1,0);
       this->get_intervals_equal_number("_PEAK_HEIGHT_",min_aux,max_aux);

       cout<<endl;
       So.message_screen("Mass bins with equal number of tracers identified");
       for(int i=0;i<=this->Number_of_tracers_in_ph_bins.size();++i)
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
real_prec Catalog::pearson_correlation(string name_X, string name_Y){
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
       for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].vmax)-Ymean);
      }
      else if(name_Y=="_RS_")
        {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].rs)-Ymean);
        }
        else if(name_Y=="_RVIR_")
            {
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].rvir)-Ymean);
            }
          else if(name_Y=="_CONCENTRATION_")
              {
      #ifdef _USE_OMP_
      #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
      #endif
                for(ULONG i=0;i<this->Halo.size();++i)
                 corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].concentration)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].spin)-Ymean);
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
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].virial-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].peak_height-Ymean);

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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(log10(this->Halo[i].mass)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].rs)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].rvir)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].concentration)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].spin)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vmax)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].virial-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vmax)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=(log10(this->Halo[i].rs)-Xmean)*(log10(this->Halo[i].rvir)-Ymean);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=(log10(this->Halo[i].rs)-Xmean)*(log10(this->Halo[i].concentration)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].rs)-Xmean)*(log10(this->Halo[i].spin)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].rs)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=((this->Halo[i].rs)-Xmean)*(this->Halo[i].virial-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].rs)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
                       for(ULONG i=0;i<this->Halo.size();++i)
                       corr+=(log10(this->Halo[i].rvir)-Xmean)*(log10(this->Halo[i].concentration)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].rvir)-Xmean)*(log10(this->Halo[i].spin)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].rvir)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].virial-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
             for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].bias-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
                   corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].peak_height-Ymean);

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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].rvir)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].concentration)-Xmean)*(log10(this->Halo[i].spin)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].concentration)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].virial-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
             for(ULONG i=0;i<this->Halo.size();++i)
               corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].bias-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
                   corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].peak_height-Ymean);

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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].concentration)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=(log10(this->Halo[i].spin)-Xmean)*(log10(this->Halo[i].spin_bullock)-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].virial-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].bias-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(log10(this->Halo[i].vrms)-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].virial-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].bias-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].spin_bullock)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].virial-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].mach_number-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].bias-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].local_dm-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(log10(this->Halo[i].vrms)-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].b_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].virial-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].c_to_a-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].b_to_a-Xmean)*(this->Halo[i].dach_number-Ymean);
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
      for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].mach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
      for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
      for(ULONG i=0;i<this->Halo.size();++i)
       corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].c_to_a-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].bias-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].mach_number-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].bias-Xmean)*(this->Halo[i].local_overdensity-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].bias-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].bias-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].bias-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].bias-Xmean)*(this->Halo[i].dach_number-Ymean);
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

    for(ULONG i=0;i<this->Halo.size();++i)
          corr+=(this->Halo[i].local_overdensity-Xmean)*(this->Halo[i].tidal_anisotropy-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].local_overdensity-Xmean)*(this->Halo[i].local_dm-Ymean);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].local_overdensity-Xmean)*(this->Halo[i].peak_height-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].local_overdensity-Xmean)*(this->Halo[i].dach_number-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].tidal_anisotropy-Xmean)*(this->Halo[i].peak_height-Ymean);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].tidal_anisotropy-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=(this->Halo[i].tidal_anisotropy-Xmean)*(this->Halo[i].dach_number-Ymean);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].peak_height-Xmean)*(this->Halo[i].local_dm-Ymean);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=(this->Halo[i].peak_height-Xmean)*(this->Halo[i].dach_number-Ymean);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=(this->Halo[i].dach_number-Xmean)*(this->Halo[i].local_dm-Ymean);
            }
    }
     corr/=(sqrt(Xvar*Yvar)*static_cast<real_prec>(this->Halo.size()));
    }
  return corr;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::PCA(vector<string>&name_props, vector<bool>&used_prop, string extra_info, bool write_cat){
  // steps:
  //    i) take var and mean, and standarize each variable for each halo.
  int dimen=name_props.size();
  vector<real_prec>mean(dimen, 0);
  vector<real_prec>sigma(dimen, 0);
  vector<string>names_used;
  real_prec sum_sigma=0;
  for(int i=0;i<dimen;++i)
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
  for(int i=0; i<dimen;++i)
    if(true==used_prop[i])
      Np_used++;
  vector<real_prec> standarized_data(this->NOBJS*Np_used, 0);
  int ip_label=0;// index over *ALL * properties
  int ip_label_used=0;// index over *USED * properties, again
  // Las matriz de datos estandarizados la llenamos de la forma C[i][j] donde i =0..Nobs y j=0, Np_used
  // Ojp que toca escribir esto a mano siguiendo el orden que viene de PowerSpectrum
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].mass)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("\\$\log_{10} M_{vir}$\\");
  }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].vmax)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10} V_{max}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].rs)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$R_s$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].rvir)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10}Rvir$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].concentration)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10}c$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].spin)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10}\lambda$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].spin_bullock)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10}\lambda_b$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
        standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=(log10(this->Halo[i].vrms)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\log_{10} \sigma_{v}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].virial)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{V}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].b_to_a)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{T}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].c_to_a)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{E}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].mach_number)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{M}");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].bias)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$b_{h}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].local_overdensity)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\delta_{R}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].tidal_anisotropy)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{T}_{A}$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].peak_height)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\NU$");
    }
  ip_label++;
  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].dach_number)-mean[ip_label])/sigma[ip_label];
      ip_label_used++;
      names_used.push_back("$\mathcal{D}$");
    }
  ip_label++;

  if(true==used_prop[ip_label])
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<this->NOBJS;++i)
    standarized_data[index_2d(ip_label_used,i,this->NOBJS)]=((this->Halo[i].local_dm)-mean[ip_label])/sigma[ip_label];
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
      for(ULONG i=0; i<this->NOBJS;++i) // loop over tracers
          cova_matrix_b[index_2d(ip,jp,ip_label_used)]+=(standarized_data[index_2d(ip,i,this->NOBJS)]*standarized_data[index_2d(jp,i,this->NOBJS)])/static_cast<real_prec>(this->NOBJS);
  this->So.DONE();
  this->So.message_screen("\tComputing eigenvalues and eigenvectors");
  vector<s_eigen_vector>eigen_v (ip_label_used);
  for(int i=0; i<ip_label_used ;++i)
    eigen_v[i].eigen_vec.resize(ip_label_used,0);
  get_eigen(cova_matrix_b, eigen_v);
  this->So.DONE();
  // explained variance EV(lambda_j) = sum i=0,j (eigenalues)/ sum(eigenvalues)
  real_prec sum_eigen=0;
  for(int i=0; i<ip_label_used ;++i)
     sum_eigen+=abs(eigen_v[i].eigen_val);
  So.message_screen("\tTotal variation:",sum_eigen);
  vector<real_prec>iv(ip_label_used,0);//individual variance variance for each PC
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0; i<ip_label_used ;++i)
      iv[i]=eigen_v[i].eigen_val/sum_eigen;
  vector<real_prec>cv(ip_label_used,0); //cumulative variance for each PC
  for(int i=0; i<ip_label_used ;++i)
     for(int j=0; j<i+1 ;++j)
       cv[i]+=abs(eigen_v[j].eigen_val)/sum_eigen;
  string file_pca=this->params._Output_directory()+"PCA.txt"+extra_info;
  this->So.message_screen("\tWriting Eigenvalues for PCA in file", file_pca);
  ofstream pca;pca.open(file_pca.c_str());
  for(ULONG j=0; j<ip_label_used ;++j)
      pca<<j+1<<"  "<<eigen_v[j].eigen_val<<"  "<<cv[j]<<"  "<<iv[j]<<"  "<<pow(sigma[j],2)/sum_sigma<<endl; // el signma estará bien acá solo cuando usemos todas las propiedades
  pca.close();
  this->So.DONE();
  vector<real_prec> new_data(this->NOBJS*ip_label_used, 0);
  //    ****************************************
  // Rotate the data to the new basis build by the eigen vectors: NEW_DATA= EIGEN_VEC * STANDARIZED_DATA
  // Ech object is now assigned a set of N PComponents, which are a linear combination of the original properties with weights given by the partial variance
  // (components of the eigenvectors)
  this->So.message_screen("\tComputing new data (PC)");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0; i<this->NOBJS;++i)
    for(ULONG ip=0; ip<ip_label_used ;++ip)// loop over the number of PC: each PC has its own eigenvalue, eigenvalue
        for(ULONG jp=0; jp<ip_label_used ;++jp)//loop over the components of each eigenvector corresponding to each eigenvalue
              new_data[index_2d(ip,i,this->NOBJS)]+=eigen_v[ip].eigen_vec[jp]*standarized_data[index_2d(jp,i,this->NOBJS)];
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
  for(ULONG ip=0; ip<ip_label_used ;++ip)
   {
     real_prec mean=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean)
#endif
     for(ULONG i=0; i<this->NOBJS;++i)
        mean+=new_data[index_2d(ip,i,this->NOBJS)];
     mean_new_data[ip]=mean/static_cast<real_prec>(this->NOBJS);
     mean=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean)
#endif
     for(ULONG i=0; i<this->NOBJS;++i)
        mean+=pow(new_data[index_2d(ip,i,this->NOBJS)]-mean_new_data[ip],2);
     sigma_new_data[ip]=sqrt(mean/static_cast<real_prec>(this->NOBJS));
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
      for(ULONG i=0; i<this->NOBJS;++i)// loop over the tracers
         correlation[index_2d(ip,jp,N_principal_comps)]+=standarized_data[index_2d(ip,i,this->NOBJS)]*((new_data[index_2d(jp,i,this->NOBJS)]-mean_new_data[jp])/sigma_new_data[jp])/static_cast<real_prec>(this->NOBJS);
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
void Catalog::Get_SO_tracer(){
    real_prec factor=0.001; //to convert Kpc to Mpc
    real_prec mean_overdensity=0;
    real_prec var_overdensity=0;
#pragma omp parallel for reduction(+:mean_overdensity)
    for(ULONG i=0;i<this->NOBJS;++i)
         mean_overdensity+=this->Cosmo.SO(this->params._redshift(), this->Halo[i].mass, this->Halo[i].rvir*factor);
    mean_overdensity/=static_cast<real_prec>(this->NOBJS);
#pragma omp parallel for reduction(+:mean_overdensity)
    for(ULONG i=0;i<this->NOBJS;++i)
         var_overdensity+=pow(this->Cosmo.SO(this->params._redshift(), this->Halo[i].mass, this->Halo[i].rvir*factor)-mean_overdensity,2);
    var_overdensity/=static_cast<real_prec>(this->NOBJS);
    So.message_screen("\tMean spherical overdensty of halos:", mean_overdensity);
    So.message_screen("\tSigma_SO of halos:", sqrt(var_overdensity));
    So.message_screen("\tBryan and Norman SO:", Cosmo.density_contrast_top_hat(this->params._redshift()));
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::Get_Ranked_Props(string prop){
    this->Halo_ranked.resize(this->Halo.size());
    gsl_vector *v=gsl_vector_alloc(this->Halo.size());
    gsl_permutation * perm = gsl_permutation_alloc(this->Halo.size());
    gsl_permutation * rank = gsl_permutation_alloc(this->Halo.size());
    // The Halo_ranked[i].prop stores the rank in a sorted list of the proeprty prop.
    // The Halo_ranked[i].prop=0 corresponds to the lowest value of the quantity prop
    // The Halo_ranked[i].prop=NOBS-1 corresponds to the largest value of the quantity prop
    if(prop=="_MASS_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i, log10(this->Halo[i].mass));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].mass=rank->data[i];
       }
    }

    if(prop=="_VMAX_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].vmax));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].vmax=rank->data[i];
       }
    }
    if(prop=="_RS_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].rs));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].rs=rank->data[i];
       }
    }
    if(prop=="_RVIR_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].rvir));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].rvir=rank->data[i];
       }
    }
    if(prop=="_CONCENTRATION_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].concentration));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].concentration=rank->data[i];
       }
    }
    if(prop=="_SPIN_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].spin));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].spin=rank->data[i];
       }
    }
    if(prop=="_SPIN_BULLOCK_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].spin_bullock));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].spin_bullock=rank->data[i];
       }
    }
    if(prop=="_VRMS_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,log10(this->Halo[i].vrms));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].vrms=rank->data[i];
       }
    }

    if(prop=="_VIRIAL_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].virial);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].virial=rank->data[i];
       }
    }
    if(prop=="_BTOA_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].b_to_a);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].b_to_a=rank->data[i];
       }
    }
    if(prop=="_CTOA_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].c_to_a);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].c_to_a=rank->data[i];
       }
    }

    if(prop=="_MACH_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].mach_number);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].mach_number=rank->data[i];
       }
    }
    if(prop=="_BIAS_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].bias);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         this->Halo_ranked[i].bias=rank->data[i];
       }
    }
    if(prop=="_RELATIVE_BIAS_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].relative_bias);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].relative_bias=rank->data[i];
       }
    }
    if(prop=="_LOCAL_OVERDENSITY_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].local_overdensity);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].local_overdensity=rank->data[i];
       }
    }
    if(prop=="_TIDAL_ANISOTROPY_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].tidal_anisotropy);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].tidal_anisotropy=rank->data[i];
       }
    }
    if(prop=="_PEAK_HEIGHT_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v,i, this->Halo[i].peak_height);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].peak_height=rank->data[i];
       }
    }
    if(prop=="_DACH_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].dach_number);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].dach_number=rank->data[i];
       }
    }
    if(prop=="_LOCALDM_")
    {
      for(ULONG i=0;i<this->NOBJS;++i)
        gsl_vector_set(v, i,this->Halo[i].local_dm);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < this->NOBJS; i++)
       {
         double vi = gsl_vector_get(v, i);
         this->Halo_ranked[i].local_dm=rank->data[i];
       }
    }
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::Get_Ranked_Props(vector<s_Halo>&inHalo, string prop){
    this->Halo_ranked.resize(inHalo.size());
    gsl_vector *v=gsl_vector_alloc(inHalo.size());
    gsl_permutation * perm = gsl_permutation_alloc(inHalo.size());
    gsl_permutation * rank = gsl_permutation_alloc(inHalo.size());

    ULONG Nsize=inHalo.size();
    if(prop=="_MASS_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i, log10(inHalo[i].mass));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].mass=rank->data[i];
    }

    if(prop=="_VMAX_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].vmax));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].vmax=rank->data[i];
    }

    if(prop=="_RS_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].rs));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].rs=rank->data[i];
    }
    if(prop=="_RVIR_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].rvir));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].rvir=rank->data[i];
    }
    if(prop=="_CONCENTRATION_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].concentration));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].concentration=rank->data[i];
    }

    if(prop=="_SPIN_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].spin));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].spin=rank->data[i];
    }
    if(prop=="_SPIN_BULLOCK_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].spin_bullock));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].spin_bullock=rank->data[i];
    }
    if(prop=="_VRMS_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,log10(inHalo[i].vrms));
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].vrms=rank->data[i];
    }

    if(prop=="_VIRIAL_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].virial);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].virial=rank->data[i];
    }
    if(prop=="_BTOA_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].b_to_a);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].b_to_a=rank->data[i];
    }

    if(prop=="_CTOA_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].c_to_a);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].c_to_a=rank->data[i];
    }

    if(prop=="_MACH_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].mach_number);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].mach_number=rank->data[i];
    }
    if(prop=="_BIAS_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].bias);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].bias=rank->data[i];
    }
    if(prop=="_RELATIVE_BIAS_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].relative_bias);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].relative_bias=rank->data[i];
    }

    if(prop=="_LOCAL_OVERDENSITY_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].local_overdensity);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].local_overdensity=rank->data[i];
    }
    if(prop=="_TIDAL_ANISOTROPY_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].tidal_anisotropy);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].tidal_anisotropy=rank->data[i];
    }
    if(prop=="_PEAK_HEIGHT_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v,i, inHalo[i].peak_height);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].peak_height=rank->data[i];
    }

    if(prop=="_DACH_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].dach_number);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].dach_number=rank->data[i];
    }
    if(prop=="_LOCALDM_")
    {
      for(ULONG i=0;i<Nsize;++i)
        gsl_vector_set(v, i,inHalo[i].local_dm);
      gsl_sort_vector_index (perm, v);
      gsl_permutation_inverse (rank, perm);
      for (ULONG i = 0; i < Nsize; i++)
         this->Halo_ranked[i].local_dm=rank->data[i];
    }
    gsl_permutation_free (perm);
    gsl_permutation_free (rank);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
real_prec Catalog::spearman_correlation(string name_X, string name_Y){
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
       for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass - this->Halo_ranked[i].vmax,2);
      }
      else if(name_Y=="_RS_")
        {
#ifdef _USE_OMP_
#pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
#endif
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].rs,2);
        }
        else if(name_Y=="_RVIR_")
            {
    #ifdef _USE_OMP_
    #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
    #endif
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].rvir,2);
            }
          else if(name_Y=="_CONCENTRATION_")
              {
      #ifdef _USE_OMP_
      #pragma omp parallel for num_threads(NTHREADS) reduction(+:corr)
      #endif
                for(ULONG i=0;i<this->Halo.size();++i)
                 corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].concentration,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].spin,2);
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
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].spin_bullock,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].vrms,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].virial,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].b_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass -this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
           corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].tidal_anisotropy,2);
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
              for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].peak_height,2);

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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].mass-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].rs,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].rvir,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].concentration,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].spin,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].spin_bullock,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].vrms,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].virial,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].b_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vmax-this->Halo_ranked[i].dach_number,2);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].rvir,2);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].concentration,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].spin,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].spin_bullock,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].vrms,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].virial,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].b_to_a,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].peak_height,2);

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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].rs-this->Halo_ranked[i].dach_number,2);
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
                       for(ULONG i=0;i<this->Halo.size();++i)
                       corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].concentration,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].spin,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].spin_bullock,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].vrms,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].virial,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].b_to_a,2);
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
             for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].c_to_a,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].mach_number,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].bias,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].local_overdensity,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].tidal_anisotropy,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].local_dm,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
                   corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].peak_height,2);

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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].rvir-this->Halo_ranked[i].dach_number,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].spin,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].spin_bullock,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].vrms,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].virial,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].b_to_a,2);
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
             for(ULONG i=0;i<this->Halo.size();++i)
               corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].c_to_a,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].mach_number,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].bias,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].local_overdensity,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].tidal_anisotropy,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].local_dm,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
                   corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].peak_height,2);

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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].concentration-this->Halo_ranked[i].dach_number,2);
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
                  for(ULONG i=0;i<this->Halo.size();++i)
                  corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].spin_bullock,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].vrms,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].virial,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].b_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].bias,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].local_overdensity,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].tidal_anisotropy,2);
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
               for(ULONG i=0;i<this->Halo.size();++i)
             corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].local_dm,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin-this->Halo_ranked[i].dach_number,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].vrms,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].virial,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].b_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].bias,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].local_overdensity,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].local_dm,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].spin_bullock-this->Halo_ranked[i].dach_number,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].virial,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].b_to_a,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].c_to_a,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].mach_number,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].bias,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].local_overdensity,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].local_dm,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].vrms-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].b_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].virial-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].c_to_a,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].b_to_a-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].mach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].c_to_a-this->Halo_ranked[i].dach_number,2);
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
       for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].bias,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].mach_number-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].bias-this->Halo_ranked[i].local_overdensity,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].bias-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].bias-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].bias-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].bias-this->Halo_ranked[i].dach_number,2);
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
    for(ULONG i=0;i<this->Halo.size();++i)
          corr+=pow(this->Halo_ranked[i].local_overdensity-this->Halo_ranked[i].tidal_anisotropy,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].local_overdensity-this->Halo_ranked[i].local_dm,2);
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
        for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].local_overdensity-this->Halo_ranked[i].peak_height,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].local_overdensity-this->Halo_ranked[i].dach_number,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].tidal_anisotropy-this->Halo_ranked[i].peak_height,2);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].tidal_anisotropy-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
        corr+=pow(this->Halo_ranked[i].tidal_anisotropy-this->Halo_ranked[i].dach_number,2);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].peak_height-this->Halo_ranked[i].local_dm,2);
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
          for(ULONG i=0;i<this->Halo.size();++i)
              corr+=pow(this->Halo_ranked[i].peak_height-this->Halo_ranked[i].dach_number,2);
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
              for(ULONG i=0;i<this->Halo.size();++i)
            corr+=pow(this->Halo_ranked[i].dach_number-this->Halo_ranked[i].local_dm,2);
            }
    }
      real_prec nnn=static_cast<real_prec>(this->Halo.size());
      corr = 1-(6.0*corr)/(nnn*(nnn*nnn-1.));
      return corr;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This function gets the scaling relation P(V|theta)
// This methods is replicated from Bam class, made here simpelrs as we only aim at measuring P(theta|M)
void Catalog::get_scaling_relation_primary_property(string prop)
{
  this->So.enter(__PRETTY_FUNCTION__);
  ULONG Nbins=this->params._NPROPbins_bam();
  ULONG LENGHT_AB_ONE = Nbins*Nbins ;
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
      for(int i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop=="_VMAX_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
         ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
          ULONG I_Y= get_bin(log10(this->Halo[ig].vmax),lp_min,Nbins,delta_sec, true);
         this->ABUNDANCE[index_2d(I_Y,I_X,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop=="_SPIN_" || prop=="_SPIN_BULLOCK_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].spin_bullock),lp_min,Nbins,delta_sec, true);
       this->ABUNDANCE[index_2d(I_Y,I_X,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop=="_CONCENTRATION_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].concentration),lp_min,Nbins,delta_sec  , true);
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
  for(ULONG tr_x = 0; tr_x < Nbins; ++tr_x)
  {
      aux_a=-10000;
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins)]));
          aux_a=aux_b;
        }
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          ULONG indexa=index_2d(tr_y, tr_x,Nbins);
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
  for(ULONG ig = 0; ig< this->NOBJS; ++ig)
   {
      ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);//ojo con el numero de bies cá, arreglarlo
      real_prec prob=-10.0;
      real_prec ran=10.0;
      int i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        prob = this->ABUNDANCE_normalized[index_2d(i_halo_v_bin, I_X,Nbins)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop=="_VMAX_")
       this->Halo[ig].vmax = pow(10,lp_halo);
      else if(prop=="_CONCENTRATION_")
       this->Halo[ig].concentration = pow(10,lp_halo);
      else if(prop=="_SPIN_BULLOCK_")
       this->Halo[ig].spin_bullock= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_scaling_relation_primary_property(string prop_halo, string prop_env)
{
  this->So.enter(__PRETTY_FUNCTION__);
  ULONG Nbins=this->params._NPROPbins_bam();
  ULONG LENGHT_AB_ONE = Nbins*Nbins*Nbins ;
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
      for(int i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop_halo=="_VMAX_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
         ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
         ULONG I_Y= get_bin(log10(this->Halo[ig].vmax),lp_min,Nbins,delta_sec, true);
         ULONG I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_VIRIAL_")
             I_Z= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_CTOA_")
             I_Z= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
          this->ABUNDANCE[index_3d(I_Y,I_X,I_Z,Nbins,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop_halo=="_SPIN_" || prop_halo=="_SPIN_BULLOCK_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].spin_bullock),lp_min,Nbins,delta_sec, true);
       ULONG I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
       if(prop_env=="_DACH_")
          I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
       if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_VIRIAL_")
          I_Z= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_CTOA_")
          I_Z= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
       this->ABUNDANCE[index_3d(I_Y,I_X,I_Z,Nbins,Nbins)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop_halo=="_CONCENTRATION_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].concentration),lp_min,Nbins,delta_sec  , true);
       ULONG I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_DACH_")
          I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_VIRIAL_")
          I_Z= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_CTOA_")
          I_Z= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
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
  for(ULONG tr_x = 0; tr_x < Nbins*Nbins; ++tr_x)
  {
      aux_a=-10000;
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins)]));
          aux_a=aux_b;
        }
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          ULONG indexa=index_2d(tr_y, tr_x,Nbins*Nbins);
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
  for(ULONG ig = 0; ig< this->NOBJS; ++ig)
   {
      ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
      ULONG I_Z= 0;
      if(prop_env=="_MACH_")
         I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_DACH_")
         I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCALDM_")
         I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_TIDAL_ANISOTROPY_")
         I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCAL_OVERDENSITY_")
         I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_VIRIAL_")
         I_Z= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_CTOA_")
         I_Z= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
      real_prec prob=-10.0;
      real_prec ran=10.0;
      int i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        prob = this->ABUNDANCE_normalized[index_3d(i_halo_v_bin, I_X,I_Z, Nbins,Nbins)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop_halo=="_VMAX_")
       this->Halo[ig].vmax = pow(10,lp_halo);
      else if(prop_halo=="_CONCENTRATION_")
       this->Halo[ig].concentration = pow(10,lp_halo);
      else if(prop_halo=="_SPIN_BULLOCK_")
       this->Halo[ig].spin_bullock= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_scaling_relation_primary_property(string prop_halo, string prop_env, string prop_extra)
{
  this->So.enter(__PRETTY_FUNCTION__);
  ULONG Nbins=this->params._NPROPbins_bam();
  int Ncwt=4;
  ULONG LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Ncwt ;
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
      for(int i=0;i<Nbins;++i)
        {
          aux_min[i]=lp_min+i*delta_sec;
          aux_max[i]=lp_min+(i+1)*delta_sec;
        }
  if(prop_halo=="_BIAS_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
         ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
         ULONG I_Y= get_bin(this->Halo[ig].bias,lp_min,Nbins,delta_sec, true);
         int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
         ULONG I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
          ULONG I_W=0;
           if(prop_extra=="_MACH_")
              I_W= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_DACH_")
              I_W= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCALDM_")
              I_W= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_TIDAL_ANISOTROPY_")
              I_W= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCAL_OVERDENSITY_")
              I_W= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_CTOA_")
              I_W= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_VIRIAL_")
              I_W= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
          this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  if(prop_halo=="_VMAX_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
         ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
         ULONG I_Y= get_bin(log10(this->Halo[ig].vmax),lp_min,Nbins,delta_sec, true);
         int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
         ULONG I_Z=0;
          if(prop_env=="_MACH_")
             I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_DACH_")
             I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
          if(prop_env=="_LOCALDM_")
             I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_TIDAL_ANISOTROPY_")
             I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
          else if(prop_env=="_LOCAL_OVERDENSITY_")
             I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
          ULONG I_W=0;
           if(prop_extra=="_MACH_")
              I_W= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_DACH_")
              I_W= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCALDM_")
              I_W= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_TIDAL_ANISOTROPY_")
              I_W= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_LOCAL_OVERDENSITY_")
              I_W= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_CTOA_")
              I_W= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
           else if(prop_extra=="_VIRIAL_")
              I_W= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
          this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->VMAXBmin;
      aux_max=this->VMAXBmax;
  }
  else if(prop_halo=="_SPIN_" || prop_halo=="_SPIN_BULLOCK_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].spin_bullock),lp_min,Nbins,delta_sec, true);
       int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
       ULONG I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
       if(prop_env=="_DACH_")
          I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
       if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
       ULONG I_W=0;
        if(prop_extra=="_MACH_")
           I_W= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_DACH_")
           I_W= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCALDM_")
           I_W= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
        else if(prop_env=="_TIDAL_ANISOTROPY_")
           I_W= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCAL_OVERDENSITY_")
           I_W= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_CTOA_")
           I_W= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_VIRIAL_")
           I_W= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
        this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
      }
      aux_min=this->SPINBmin;
      aux_max=this->SPINBmax;
  }
  else if(prop_halo=="_CONCENTRATION_")
    {
      for(ULONG ig = 0; ig< this->NOBJS; ++ig)
      {
       ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
       ULONG I_Y= get_bin(log10(this->Halo[ig].concentration),lp_min,Nbins,delta_sec  , true);
       int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
       ULONG I_Z= 0;
       if(prop_env=="_MACH_")
          I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_DACH_")
          I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCALDM_")
          I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_TIDAL_ANISOTROPY_")
          I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
       else if(prop_env=="_LOCAL_OVERDENSITY_")
          I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
       ULONG I_W=0;
        if(prop_extra=="_MACH_")
           I_W= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_DACH_")
           I_W= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCALDM_")
           I_W= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
        else if(prop_env=="_TIDAL_ANISOTROPY_")
           I_W= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_LOCAL_OVERDENSITY_")
           I_W= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_CTOA_")
           I_W= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
        else if(prop_extra=="_VIRIAL_")
           I_W= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
        this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
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
  for(ULONG tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
  {
      aux_a=-10000;
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
          aux_a=aux_b;
        }
      for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
        {
          ULONG indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
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
  for(ULONG ig = 0; ig< this->NOBJS; ++ig)
   {
      ULONG I_X= get_bin(log10(this->Halo[ig].mass),lm_min,Nbins,this->logdeltaM, true);
      int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
      ULONG I_Z= 0;
      if(prop_env=="_MACH_")
         I_Z= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_DACH_")
         I_Z= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCALDM_")
         I_Z= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_TIDAL_ANISOTROPY_")
         I_Z= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
      else if(prop_env=="_LOCAL_OVERDENSITY_")
         I_Z= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
      ULONG I_W=0;
       if(prop_extra=="_MACH_")
          I_W= get_bin(this->Halo[ig].mach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_DACH_")
          I_W= get_bin(this->Halo[ig].dach_number,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_LOCALDM_")
          I_W= get_bin(this->Halo[ig].local_dm,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_TIDAL_ANISOTROPY_")
          I_W= get_bin(this->Halo[ig].tidal_anisotropy,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_LOCAL_OVERDENSITY_")
          I_W= get_bin(this->Halo[ig].local_overdensity,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_CTOA_")
          I_W= get_bin(this->Halo[ig].c_to_a,lenv_min,Nbins,delta_env, true);
       else if(prop_extra=="_VIRIAL_")
          I_W= get_bin(this->Halo[ig].virial,lenv_min,Nbins,delta_env, true);
      real_prec prob=-10.0;
      real_prec ran=10.0;
      int i_halo_v_bin=0;
      while(prob<ran)
      {
        i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
        prob = this->ABUNDANCE_normalized[index_5d(i_halo_v_bin, I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)];
        ran = gsl_rng_uniform(this->rn_cat);
      }
      real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
      real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
      if(prop_halo=="_VMAX_")
       this->Halo[ig].vmax = pow(10,lp_halo);
      if(prop_halo=="_BIAS_")
       this->Halo[ig].bias = lp_halo;
      else if(prop_halo=="_CONCENTRATION_")
       this->Halo[ig].concentration = pow(10,lp_halo);
      else if(prop_halo=="_SPIN_BULLOCK_")
       this->Halo[ig].spin_bullock= pow(10,lp_halo);
   }
  this->So.DONE();
  this->ABUNDANCE_normalized.clear();
  this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_scaling_relation_bias()
{
  this->So.enter(__PRETTY_FUNCTION__);
  ULONG Nbins=this->params._NPROPbins_bam();
  int Ncwt=4;
  ULONG LENGHT_AB_ONE = Nbins*Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt ;
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
  for(int i=0;i<Nbins;++i)
    {
      aux_min[i]=lp_min+i*delta_sec;
      aux_max[i]=lp_min+(i+1)*delta_sec;
    }
  for(ULONG ig = 0; ig< this->NOBJS; ++ig)
  {
     ULONG I_Y= get_bin(this->Halo[ig].bias,lp_min,Nbins,delta_sec, true);// bin in bias
     ULONG I_X= get_bin(this->Halo[ig].local_dm,lenv1_min,Nbins,delta_env1, true);//bin of local dm
     int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property // cwt
     ULONG I_Z= get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);//local overdenisty
     ULONG I_W= get_bin(this->Halo[ig].mach_number,lenv3_min,Nbins,delta_env3, true);//local overdenisty
     ULONG I_D= get_bin(this->Halo[ig].dach_number,lenv4_min,Nbins,delta_env4, true);//local overdenisty
     ULONG I_T= get_bin(this->Halo[ig].tidal_anisotropy,lenv5_min,Nbins,delta_env5, true);//local overdenisty
     this->ABUNDANCE[index_7d(I_Y,I_X,I_Z,I_W,I_D,I_T,I_CWT,Nbins,Nbins,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
  }
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    long aux_a,aux_b;
    So.message_screen("Normalizing");
    for(ULONG tr_x = 0; tr_x < Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         ULONG indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Nbins*Nbins*Ncwt);
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
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
    {
        this->Halo[ig].bias2 = this->Halo[ig].bias; //keep track of the original bias
        ULONG I_X= get_bin(this->Halo[ig].local_dm,lenv1_min,Nbins,delta_env1, true);
        int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
        ULONG I_Z = get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);
        ULONG I_W= get_bin(this->Halo[ig].mach_number,lenv3_min,Nbins,delta_env3, true);//local overdenisty
        ULONG I_D= get_bin(this->Halo[ig].dach_number,lenv4_min,Nbins,delta_env4, true);//local overdenisty
        ULONG I_T= get_bin(this->Halo[ig].tidal_anisotropy,lenv5_min,Nbins,delta_env5, true);//local overdenisty
        real_prec prob=-10.0;
        real_prec ran=10.0;
        int i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           prob = this->ABUNDANCE_normalized[index_7d(i_halo_v_bin, I_X,I_Z,I_W,I_D,I_T,I_CWT,Nbins,Nbins,Nbins,Nbins,Nbins,Ncwt)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->Halo[ig].bias = lp_halo;
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
    cout<<prop_min<<"  "<<prop_max<<"  "<<delta_prop<<endl;
    for(int i=0;i<Nbins;++i)
      {
        aux_min[i]=prop_min+i*delta_prop;
        aux_max[i]=prop_min+(i+1)*delta_prop;
      }
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
    {
     ULONG I_Y= get_bin(log10(this->Halo[ig].mass),prop_min,Nbins,delta_prop, true);
     ULONG I_X= get_bin(this->Halo[ig].local_dm,lenv1_min,Nbins,delta_env1, true);
     ULONG I_Z = get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);
     ULONG I_W = get_bin(this->Halo[ig].bias2,lp_min,Nbins,delta_sec, true); // use the original bias to learn from, bias2
     int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
     this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
    }
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    So.message_screen("Normalizing");
    for(ULONG tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         ULONG indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
    }
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->So.DONE();
    // Assign new mass:
    So.message_screen("Assingning new mass");
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
    {
        this->Halo[ig].mass_parent = this->Halo[ig].mass; // keep track of original mass
        ULONG I_X= get_bin(this->Halo[ig].local_dm,lenv1_min,Nbins,delta_env1, true);
        int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
        ULONG I_Z = get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);
        ULONG I_W = get_bin(this->Halo[ig].bias,lp_min,Nbins,delta_sec, true);// use the new version of bias
        real_prec prob=-10.0;
        real_prec ran=10.0;
        int i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           prob = this->ABUNDANCE_normalized[index_5d(i_halo_v_bin, I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->Halo[ig].mass = pow(10,lp_halo);
        }
        this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    ofstream tea; tea.open("test.txt");
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
        tea<<this->Halo[ig].mass<<"  "<<this->Halo[ig].mass_parent<<"  "<<this->Halo[ig].bias<<" "<<this->Halo[ig].bias2<<endl;
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
    for(int i=0;i<Nbins;++i)
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
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
    {
     ULONG I_Y= get_bin(log10(this->Halo[ig].vmax),prop_min,Nbins,delta_prop, true);
     ULONG I_X= get_bin(log10(this->Halo[ig].mass_parent),lenv1_min,Nbins,delta_env1, true);// use original mass to learn from
     int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
     ULONG I_Z = get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);
     ULONG I_W = get_bin(this->Halo[ig].bias2,lp_min,Nbins,delta_sec, true); // use the original bias to learn from, bias2
     this->ABUNDANCE[index_5d(I_Y,I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)]++; // P(M, theta), Scaling relation of interest is P(Y|X), Y = secondary prop,m X=primary
    }
    this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
    this->ABUNDANCE_normalized.resize(LENGHT_AB_ONE,0);
    So.message_screen("Normalizing");
    for(ULONG tr_x = 0; tr_x < Nbins*Nbins*Nbins*Ncwt; ++tr_x)
    {
     aux_a=-10000;
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         aux_b=max(aux_a, static_cast<long>(this->ABUNDANCE[index_2d(tr_y,tr_x,Nbins*Nbins*Nbins*Ncwt)]));
         aux_a=aux_b;
       }
     for(ULONG tr_y = 0; tr_y < Nbins; ++tr_y)
       {
         ULONG indexa=index_2d(tr_y, tr_x,Nbins*Nbins*Nbins*Ncwt);
         long AUX=this->ABUNDANCE[indexa];
         this->ABUNDANCE_normalized[indexa]=(aux_b==0 ? 0. : static_cast<real_prec>(AUX)/static_cast<real_prec>(aux_b));
       }
    }
    this->ABUNDANCE.clear();
    this->ABUNDANCE.shrink_to_fit();
    this->So.DONE();
    // Assign new mass
    So.message_screen("Assingning new vmax");
    for(ULONG ig = 0; ig< this->NOBJS; ++ig)
    {
        ULONG I_X= get_bin(log10(this->Halo[ig].mass),lenv1_min,Nbins,delta_env1, true);/// use reconstructed mass
        ULONG I_Z = get_bin(this->Halo[ig].local_overdensity,lenv2_min,Nbins,delta_env2, true);
        ULONG I_W = get_bin(this->Halo[ig].bias,lp_min,Nbins,delta_sec, true);// use the new version of bias
        int I_CWT=this->Halo[ig].gal_cwt-1; // -1 to link it with bins in cwt property
        real_prec prob=-10.0;
        real_prec ran=10.0;
        int i_halo_v_bin=0;
        while(prob<ran)
         {
           i_halo_v_bin= gsl_rng_uniform_int(this->rn_cat,Nbins);
           prob = this->ABUNDANCE_normalized[index_5d(i_halo_v_bin, I_X,I_Z,I_W,I_CWT,Nbins,Nbins,Nbins,Ncwt)];
           ran = gsl_rng_uniform(this->rn_cat);
         }
        real_prec xr = static_cast<real_prec>(gsl_rng_uniform(this->rn_cat));
        real_prec lp_halo = aux_min[i_halo_v_bin] + xr*(aux_max[i_halo_v_bin]-aux_min[i_halo_v_bin]);
        this->Halo[ig].vmax = pow(10,lp_halo);
        }
    this->So.DONE();
    this->ABUNDANCE_normalized.clear();
    this->ABUNDANCE_normalized.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_superclusters(string realm, string tbias){
// real can be under, or over, referring to tracers with lowest bias (under) or hights bias
// tbias is the  BIAS or RELATIVE_BIAS
    real_prec mean_bias=0;
   if(tbias=="_BIAS_")
#pragma omp parallel for reduction(+:mean_bias)
    for (ULONG i=0;i<this->NOBJS; ++i)
        mean_bias+=this->Halo[i].bias;
   else if (tbias=="_RELATIVE_BIAS_")
#pragma omp parallel for reduction(+:mean_bias)
    for (ULONG i=0;i<this->NOBJS; ++i)
        mean_bias+=this->Halo[i].relative_bias;
    mean_bias/=static_cast<real_prec>(this->NOBJS);
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
    this->tracer_aux.clear();
    this->tracer_aux.shrink_to_fit();
    // Rank tracers according to bias bottom-up. I.e, the tracer with lowest bias has this->Halo_ranked[i].bias=0
    // and the tracer with the highest bias gets this->Halo_ranked[i].bias=NOBJS-1.
   So.message_screen("Ranking");
   this->Get_Ranked_Props(tbias);
    // We nos select the  Ntracers_highb tracers with the highst or lowest bias
   ULONG Ntracers_highb=200000;
   ULONG bcounter=0;
   So.message_screen("Getting other numbers");
   if(tbias=="_BIAS_")
   {
     for (ULONG i=0; i< this->NOBJS; ++i)
       if(this->Halo[i].bias>quartile_bias)//selects the Ntracers_highb with the highest bias from the ranked list
       bcounter++;
     So.message_screen("Number of tracers above quartile = ",bcounter);
     bcounter=0;
     for (ULONG i=0; i< this->NOBJS; ++i)
     if(this->Halo[i].bias>mean_bias)//selects the Ntracers_highb with the highest bias from the ranked list
            bcounter++;
     So.message_screen("Number of tracers above mean = ",bcounter);
    }
   else if(tbias=="_RELATIVE_BIAS_")
    {
      for (ULONG i=0; i< this->NOBJS; ++i)
       if(this->Halo[i].relative_bias>quartile_bias)//selects the Ntracers_highb with the highest bias from the ranked list
         bcounter++;
      So.message_screen("Number of tracers above quartile = ",bcounter);
      bcounter=0;
      for (ULONG i=0; i< this->NOBJS; ++i)
        if(this->Halo[i].relative_bias>mean_bias)//selects the Ntracers_highb with the highest bias from the ranked list
         bcounter++;
       So.message_screen("Number of tracers above mean = ",bcounter);
    }
    vector<real_prec>bias_aux(this->NOBJS,0);
    vector<real_prec>bias_aux_ranked(this->NOBJS,0);
    if(tbias=="_BIAS_")
#pragma omp parallel for
       for (ULONG i=0; i< this->NOBJS; ++i)
       {
            bias_aux[i]=this->Halo[i].bias;
            bias_aux_ranked[i]=this->Halo_ranked[i].bias;
       }
    else if (tbias=="_RELATIVE_BIAS_")
#pragma omp parallel for
       for (ULONG i=0; i< this->NOBJS; ++i)
       {
            bias_aux[i]=this->Halo[i].relative_bias;
            bias_aux_ranked[i]=this->Halo_ranked[i].relative_bias;
       }

    ULONG N_overbias=0;
    if(realm=="under")
     {
      for (ULONG i=0; i< this->NOBJS; ++i)
//       if(bias_aux_ranked[i]<Ntracers_highb &&  bias_aux[i]<mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
       if(bias_aux_ranked[i]<Ntracers_highb &&  bias_aux[i]<mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
        {
            this->tracer_aux.push_back(s_Halo());
            this->tracer_aux[N_overbias].coord1 = this->Halo[i].coord1;
            this->tracer_aux[N_overbias].coord2 = this->Halo[i].coord2;
            this->tracer_aux[N_overbias].coord3 = this->Halo[i].coord3;
            this->tracer_aux[N_overbias].bias = this->Halo[i].bias;
            this->tracer_aux[N_overbias].spin = this->Halo[i].spin;
            this->tracer_aux[N_overbias].concentration = this->Halo[i].concentration;
            this->tracer_aux[N_overbias].mass = this->Halo[i].mass;
            this->tracer_aux[N_overbias].mach_number = this->Halo[i].mach_number;
            this->tracer_aux[N_overbias].b_to_a = this->Halo[i].b_to_a;
            this->tracer_aux[N_overbias].c_to_a = this->Halo[i].c_to_a;
            this->tracer_aux[N_overbias].relative_bias = this->Halo[i].relative_bias;
            this->tracer_aux[N_overbias].bias = this->Halo[i].bias;
            this->tracer_aux[N_overbias].local_dm = this->Halo[i].local_dm;
            this->tracer_aux[N_overbias].tidal_anisotropy = this->Halo[i].tidal_anisotropy;
            N_overbias++;
        }
    }
   else if(realm=="over")
     {
     for (ULONG i=0; i< this->NOBJS; ++i)
       if(bias_aux_ranked[i]>this->NOBJS-Ntracers_highb &&  bias_aux[i]>mean_bias)//selects the Ntracers_highb with the lowest  bias from the ranked list
         {
            this->tracer_aux.push_back(s_Halo());
            this->tracer_aux[N_overbias].coord1 = this->Halo[i].coord1;
            this->tracer_aux[N_overbias].coord2 = this->Halo[i].coord2;
            this->tracer_aux[N_overbias].coord3 = this->Halo[i].coord3;
            this->tracer_aux[N_overbias].bias = this->Halo[i].bias;
            this->tracer_aux[N_overbias].spin = this->Halo[i].spin;
            this->tracer_aux[N_overbias].concentration = this->Halo[i].concentration;
            this->tracer_aux[N_overbias].mass = this->Halo[i].mass;
            this->tracer_aux[N_overbias].mach_number = this->Halo[i].mach_number;
            this->tracer_aux[N_overbias].b_to_a = this->Halo[i].b_to_a;
            this->tracer_aux[N_overbias].c_to_a = this->Halo[i].c_to_a;
            this->tracer_aux[N_overbias].relative_bias = this->Halo[i].relative_bias;
            this->tracer_aux[N_overbias].bias = this->Halo[i].bias;
            this->tracer_aux[N_overbias].local_dm = this->Halo[i].local_dm;
            this->tracer_aux[N_overbias].tidal_anisotropy = this->Halo[i].tidal_anisotropy;
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
    for (ULONG i=0; i< this->NOBJS; ++i)
        if(this->Halo[i].relative_bias< mean_bias+brad*factor_increase && this->Halo[i].relative_bias > mean_bias-brad*factor_increase)//selects the Ntracers_highb with the highest bias from the ranked list
        {
         index_mean+=this->Halo_ranked[i].relative_bias;
         counter_i++;
        }

        if(counter_i>0)
            found=true;
        else
            factor_increase+=0.5;
    }
    ULONG ind=static_cast<ULONG>(index_mean/counter_i); // takje this index as the rank corresponding to the mean bias in the sample
    So.message_screen("Mean rank used:", ind);
   for (ULONG i=0; i< this->NOBJS; ++i)
   {
       if(this->Halo_ranked[i].relative_bias> this->Halo_ranked[ind].relative_bias-Ntracers_highb/2 &&   this->Halo_ranked[i].relative_bias<this->Halo_ranked[ind].relative_bias+Ntracers_highb/2 )//selects the Ntracers_highb with the highest bias from the ranked list
        {
            this->tracer_aux.push_back(s_Halo());
            this->tracer_aux[N_overbias].coord1 = this->Halo[i].coord1;
            this->tracer_aux[N_overbias].coord2 = this->Halo[i].coord2;
            this->tracer_aux[N_overbias].coord3 = this->Halo[i].coord3;
            this->tracer_aux[N_overbias].bias = this->Halo[i].bias;
            this->tracer_aux[N_overbias].spin = this->Halo[i].spin;
            this->tracer_aux[N_overbias].concentration = this->Halo[i].concentration;
            this->tracer_aux[N_overbias].mass = this->Halo[i].mass;
            this->tracer_aux[N_overbias].mach_number = this->Halo[i].mach_number;
            this->tracer_aux[N_overbias].b_to_a = this->Halo[i].b_to_a;
            this->tracer_aux[N_overbias].c_to_a = this->Halo[i].c_to_a;
            this->tracer_aux[N_overbias].relative_bias = this->Halo[i].relative_bias;
            this->tracer_aux[N_overbias].local_dm = this->Halo[i].local_dm;
            this->tracer_aux[N_overbias].tidal_anisotropy = this->Halo[i].tidal_anisotropy;
            N_overbias++;
        }
    }
    */
   this->So.message_screen("Number of tracers above Nlimit and above mean=", N_overbias);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// If this is defined, the nbar is assgned by checking the bin in whcih the tracer z falls in the dndz vector. Only applied so far to I_EQZ
#define _use_simple_nbar_assignment_
#ifdef _USE_SEVERAL_RANDOM_FILES_
void Catalog::ang_to_cart_coordinates(s_data_structure *s_data, int ir){
#else
void Catalog::ang_to_cart_coordinates(s_data_structure *s_data){
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
  this->So.message_screen("Transform to cartesian coordinates in catalogue (and assigning nbar) for ", this->type_of_object);
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
  ULONG nlines = this->NOBJS;
#ifdef _use_simple_nbar_assignment_
  real_prec delta_z_nbar = 0;
  if(true==this->params._use_random_catalog())
    delta_z_nbar= s_data->zz_v[1]-s_data->zz_v[0];
#endif
  string angles_units = (this->type_of_object=="RANDOM" ? this->params._angles_units_r() : this->params._angles_units_g());
  int sys_coord = (this->type_of_object=="RANDOM" ? this->params._sys_of_coord_r(): this->params._sys_of_coord_g());
 // Factor to transform deg to rad:
  real_prec fac=1.0;
  if(angles_units=="R")
    fac=180.0/M_PI;
  int n_rc=0;
  int n_zzv=0;
#ifndef _use_simple_nbar_assignment_
  // Preparing for interpolation of the relation nbar(z)
  gsl_spline *spline_nbar;
  gsl_interp_accel *spline_acc_nbar;
  if(true==this->params._use_random_catalog() && false==this->params._nbar_tabulated() && true==this->params._use_file_nbar())
    {
      spline_nbar = gsl_spline_alloc (gsl_interp_linear,s_data->zz_v.size());
      gsl_spline_init (spline_nbar, &s_data->zz_v[0], &s_data->dndz_v[0], s_data->zz_v.size());// nbar(z)
    }
#endif
  // Preparing for interpolation of the relation r(z):
  // if needed initialize gsl structure to interpolate Comoving distance as a function of redshift:
  gsl_spline *spline_zro;
  gsl_interp_accel *spline_acc_zro;
  if((true==this->params._use_random_catalog() && false==this->params._nbar_tabulated()) || true==this->params._use_file_nbar() || true==this->params._nbar_tabulated())
    {
      n_rc = s_data->rr_c.size();
      spline_zro = gsl_spline_alloc (gsl_interp_linear, n_rc);
      switch(sys_coord)
      {
       case(I_EQR): //in this coordinate system, we might need to get z from r and therefrom, nbar from z
         gsl_spline_init (spline_zro, &s_data->rr_c[0], &s_data->zz_c[0], n_rc); //z(r))
       break;
       case(I_EQZ): //in this coordinates we have z, so we just need r(z) and nbar(z)
         gsl_spline_init (spline_zro, &s_data->zz_c[0], &s_data->rr_c[0], n_rc); //r(z)
       break;  
     }
    }
  switch(sys_coord)
    {
    case(I_CART):   // Choose Cartesian coordinates:
#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
      {
        gsl_interp_accel *spline_acc_zro;
        spline_acc_zro = gsl_interp_accel_alloc ();
#ifndef _use_simple_nbar_assignment_
        gsl_interp_accel *spline_acc_nbar;
        if(true==this->params._use_random_catalog() && false==this->params._nbar_tabulated())
          spline_acc_nbar = gsl_interp_accel_alloc ();// useless
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
#endif
      for(ULONG i=0;i<nlines;++i)
        {
          real_prec nbar=1.0;
          real_prec x=this->Halo[i].coord1;
          real_prec y=this->Halo[i].coord2;
          real_prec z=this->Halo[i].coord3;
          //compute nbar:
          if(true==this->params._use_random_catalog() && true==this->params._nbar_tabulated())
              nbar=this->Halo[i].mean_density;
          else
            this->Halo[i].mean_density=nbar;
          if(true==this->params._use_random_catalog())
            this->Halo[i].mean_density=nbar;
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
    case(I_EQR):
#ifdef _USE_OMP_
#pragma omp parallel num_threads(NTHREADS)
      {
    spline_acc_zro = gsl_interp_accel_alloc ();
    if(true==this->params._use_random_catalog() && false==this->params._nbar_tabulated() && true == this->params._use_file_nbar() )
        gsl_interp_accel *spline_acc_nbar;
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
#endif
    for(ULONG i=0;i<nlines;++i)
      {
        real_prec x,y,z, zro, nbar;
        real_prec ra_s=this->Halo[i].coord1/fac;
        real_prec dec_s=this->Halo[i].coord2/fac;
        real_prec rr=this->Halo[i].coord3;
        equatorial_to_cartesian(ra_s,dec_s,rr,x, y, z); // Transform to cartesian coord:
        if(true==this->params._use_random_catalog()) 	  // Compute the mean number density if not tabulated :
          {
            if(true==this->params._nbar_tabulated())
              nbar=this->Halo[i].mean_density;
           else
             {
                if(true==this->params._constant_depth() || true == this->params._use_file_nbar())
                  {
                    zro= gsl_spline_eval(spline_zro, this->Halo[i].coord3, spline_acc_zro);
                    this->Halo[i].redshift=zro;
#ifdef _use_simple_nbar_assignment_
                  ULONG zbin=get_bin(zro ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                      nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar); // Compute the mean number density if not tabulated, either fropm file or from nbar measured from the randoms
#endif
                  }
            else
              {
                zro= gsl_spline_eval(spline_zro, this->Halo[i].coord3, spline_acc_zro); 		      // Interpolate to get redshift
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
        this->Halo[i].coord1=x;
        this->Halo[i].coord2=y;
        this->Halo[i].coord3=z;
        if(false==this->params._use_random_catalog()|| false == this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->Halo[i].mean_density=nbar;
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
  case(I_EQZ):
#pragma omp parallel num_threads(NTHREADS)
    {
      spline_acc_zro = gsl_interp_accel_alloc ();
#ifndef _use_simple_nbar_assignment_
      if(true==this->params._use_random_catalog() && false==this->params._nbar_tabulated())
        spline_acc_nbar = gsl_interp_accel_alloc ();
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
      for(ULONG i=0;i<nlines;++i)
       {
        real_prec x=0; real_prec y=0;
        real_prec z=0; real_prec zro=0; 
        real_prec ra_s=this->Halo[i].coord1/fac;
        real_prec dec_s=this->Halo[i].coord2/fac;
        this->Halo[i].redshift=this->Halo[i].coord3; // just to keep track of the redshift of the tracer
        real_prec redshift_aux  = this->Halo[i].redshift < s_data->zz_c[0] ? s_data->zz_c[0]: this->Halo[i].redshift;
        real_prec rr=gsl_spline_eval(spline_zro, redshift_aux, spline_acc_zro); // Transform to comoving distance *
        real_prec nbar=mean_density;
        equatorial_to_cartesian(ra_s, dec_s, rr, x, y, z); 	    // Transform to cartesian coordinates
        if(true==this->params._use_random_catalog())
          if(false==this->params._nbar_tabulated())
              if(true==this->params._constant_depth() || true == this->params._use_file_nbar() )
                {
#ifdef _use_simple_nbar_assignment_
                  ULONG zbin=get_bin(redshift_aux ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                      nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar); // Compute the mean number density if not tabulated, either fropm file or from nbar measured from the randoms
#endif
                }
#ifdef HEALPIX
              else if(false==this->params._constant_depth() && false == this->params._use_file_nbar())
                 my_get_mean_density_interpolated(map, this->params._new_n_dndz,this->params._redshift_min_sample, this->params._redshift_max_sample,ra_s, dec_s, redshift_aux, dndz_m, &nbar);
#endif
#ifdef _USE_REDSHIFT_BINS_
        this->Halo[i].observed=false;
        if(this->Halo[i].coord3<=this->params._redshift_max_sample && this->Halo[i].coord3>=this->params._redshift_min_sample)
          this->Halo[i].observed=true;
#endif
        // New coordinates:
        this->Halo[i].coord1=x;
        this->Halo[i].coord2=y;
        this->Halo[i].coord3=z;
        if(false==this->params._use_random_catalog() || false==this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->Halo[i].mean_density=nbar;
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
    case(I_EQRZ):
#pragma omp parallel num_threads(NTHREADS)
      {
    spline_acc_zro = gsl_interp_accel_alloc ();
#ifndef _use_simple_nbar_assignment_
    if(true==this->params._use_random_catalog() && false==this->params._nbar_tabulated())
      spline_acc_nbar = gsl_interp_accel_alloc ();
#endif
#pragma omp for schedule(static, NTHREADS) reduction(min:aXMIN, aYMIN, aZMIN) reduction(max:aXMAX, aYMAX, aZMAX) nowait
    for(ULONG i=0;i<nlines;++i)
      {
        real_prec x,y,z,rr, zro, nbar, ra_s, dec_s;
        ra_s=this->Halo[i].coord1/fac;
        dec_s=this->Halo[i].coord2/fac;
        rr=this->Halo[i].coord3;
        equatorial_to_cartesian(ra_s, dec_s, rr, x, y, z); 	  // Transform to cartesian coord. *
        if(true==this->params._use_random_catalog()) 	  // Compute the mean number density if not tabulated                         *
          {
            if(true==this->params._nbar_tabulated())
              nbar=this->Halo[i].mean_density;
            else
             {
              if(true==this->params._constant_depth()  || true == this->params._use_file_nbar())
                {
                  zro=this->Halo[i].coord3;
#ifdef _use_simple_nbar_assignment_
                  ULONG zbin=get_bin(zro ,s_data->zz_v[0],s_data->zz_v.size(),delta_z_nbar,false);
                  nbar=s_data->dndz_v[zbin];
#else
                  nbar = gsl_spline_eval(spline_nbar, zro, spline_acc_nbar);
#endif
                }
#ifdef HEALPIX
              else
                {
                  zro=this->Halo[i].coord3;
            // my_get_mean_density_interpolated(map, this->params._new_n_dndz,this->params._redshift_min_sample, this->params._redshift_max_sample,ra_s, dec_s, zro, dndz_m, &nbar);
                }
#endif
            }
         }
        else
          nbar=mean_density;
#ifdef _USE_REDSHIFT_BINS_
        this->Halo[i].observed=false;
        if(this->Halo[i].coord3<=this->params._redshift_max_sample && this->Halo[i].coord3>=this->params._redshift_min_sample)
          this->Halo[i].observed=true;
#endif
        // New coordinates:
        this->Halo[i].coord1=x;
        this->Halo[i].coord2=y;
        this->Halo[i].coord3=z;
        if(false==this->params._use_random_catalog()|| false == this->params._nbar_tabulated()) // If we had to compute nbar, assigne it now:
          this->Halo[i].mean_density=nbar;
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
  if(sys_coord>0)
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
  cout<<"\t"<<YELLOW<<"Range in x:[" << aXMIN << ":" << aXMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"Range in y:[" << aYMIN << ":" << aYMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"Range in z:[" << aZMIN << ":" << aZMAX << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in x:[" << aXMIN-shift_x << ":" << aXMAX-shift_x << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in y:[" << aYMIN-shift_y << ":" << aYMAX-shift_y << "]" << endl;
  cout<<"\t"<<YELLOW<<"New Range in z:[" << aZMIN-shift_z << ":" << aZMAX-shift_z << "]" << endl;
  So.message_screen("Xoffset =",this->params._Xoffset());
  So.message_screen("Yoffset =",this->params._Yoffset());
  So.message_screen("Zoffset =",this->params._Zoffset());
  if(this->type_of_object=="RANDOM" )
    {
      if(true==this->params._new_Lbox())
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
      So.message_screen("Computed box from random:");
      this->params.set_Lbox(llz);
      this->params.derived_pars(); // This is only called if Lbox has changed
  }
  So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  _USE_SEVERAL_RANDOM_FILES_
  void Catalog::get_interpolated_density_field(bool marked, string property, int ir) // meant for redshift space
#else
  void Catalog::get_interpolated_density_field(bool marked, string property) // meant for redshift space
#endif
 {
   So.enter(__PRETTY_FUNCTION__);
   // Sampling of the random or real catalog and
   // interpolation of the object density field into a grid.
   // This method expects the coordinates of the input
   // catalogue already transformed into cartessian,
   // with the information of the mean number density written
   // in the corresponding slot as specified in the
   // parameter file or computed from the box.
   // Identify columns in the corresponding catalog                               *
#ifdef _USE_WEIGHTS_IN_POWER_
   int i_weight1     = (this->type_of_object!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
   int i_weight2     = (this->type_of_object!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
   int i_weight3     = (this->type_of_object!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
   int i_weight4     = (this->type_of_object!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
   bool use_weight1  = (this->type_of_object!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
   bool use_weight2  = (this->type_of_object!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
   bool use_weight3  = (this->type_of_object!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
   bool use_weight4  = (this->type_of_object!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
#endif
   ULONG nlines= this->Halo.size();
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
       mark.resize(this->NOBJS, 1);
       if  (property=="_MASS_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].mass);
       else if(property=="_VMAX_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].vmax);
       else if  (property=="_RS_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].rs;
       else if  (property=="_RVIR_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].rvir;
       else if  (property=="_CONCENTRATION_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=log10(this->Halo[i].concentration);
       else if  (property=="_SPIN_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].spin);
       else if  (property=="_SPIN_BULLOCK_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=log10(this->Halo[i].spin_bullock);
       else if  (property=="_VRMS_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].vrms);
       else if (property=="_VIRIAL_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].virial;
       else if  (property=="_BTOA_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].b_to_a;
       else if  (property=="_CTOA_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].c_to_a;
       else if  (property=="_MACH_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].mach_number;
       else if  (property=="_BIAS_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].bias;
       else if  (property=="_TIDAL_ANISOTROPY_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].tidal_anisotropy;
       else if  (property=="_LOCAL_OVERDENSITY_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].local_overdensity;
       else if  (property=="_PEAK_HEIGHT_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].peak_height;
       else if  (property=="_DACH_")
     for(ULONG i=0;i<mark.size();++i)
       mark[i]=this->Halo[i].dach_number;
     }
   // If rsd=true, then there is no need to shift coordinates even if we want redshift space pwoer spectrum, for the catalog might be already "distorted"
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
   real_prec conversion_factor=1;
   if(this->params._redshift_space_coords_g() == false)
     {
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
     }
#endif
   int exp_mas=this->params._mass_assignment();
   So.message_screen("Interpolating", nlines,"tracers of type ", this->type_of_object);
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
   for(ULONG i=0;i< nlines ;++i)
     {
#ifdef _REDSHIFT_SPACE_
       vx=rsd_x*this->Halo[i].vel1*conversion_factor;
       vy=rsd_y*this->Halo[i].vel2*conversion_factor;
       vz=rsd_z*this->Halo[i].vel3*conversion_factor;
#endif
       real_prec x=this->Halo[i].coord1+vx;
       real_prec y=this->Halo[i].coord2+vy;
       real_prec z=this->Halo[i].coord3+vz;
#ifdef _USE_MASS_AS_OBSERVABLE_
       real_prec property=this->Halo[i].mass;
#elif defined _USE_VMAX_AS_OBSERVABLE_
       real_prec property=this->Halo[i].vmax;
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
       if(property>=MINIMUM_PROP_CUT)
     {
#endif
#ifdef _MASS_WEIGHT_POWER_
       real_prec mass=this->Halo[i].mass;
#endif
#ifdef _USE_REDSHIFT_BINS_
       if(true==this->Halo[i].observed)
         {
#endif
           double nbar=static_cast<double>(this->mean_density);
           if(true==this->params._use_random_catalog())// nbar has been already tbulated or assigned, so it's enough asking if we use random cats or not
            nbar=static_cast<double>(this->Halo[i].mean_density);
           double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
           ow[0] = (use_weight1 && (i_weight1<n_columns))? this->Halo[i].weight1 : 1.0;
           ow[1] = (use_weight2 && (i_weight2<n_columns))? this->Halo[i].weight2 : 1.0;
           ow[2] = (use_weight3 && (i_weight3<n_columns))? this->Halo[i].weight3 : 1.0;
           ow[3] = (use_weight4 && (i_weight4<n_columns))? this->Halo[i].weight4 : 1.0;
           ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
           double we_fkp=1.0;
           if(true==this->params._FKP_weight())
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
           case(I_NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
             grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_TSC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_PCS):
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
       for(ULONG i=0;i<this->field_external_marked.size();++i)
           this->field_external_marked[i]+=field[i];
     }
   else
    {
      if(ir==0)
        this->field_external.resize(field.size(),0);
      for(ULONG i=0;i<this->field_external.size();++i)
        this->field_external[i]+=field[i];
     }
#else
   this->n_gal = static_cast<ULONG>(n_selected);
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
 void Catalog::get_interpolated_density_field_real_space(bool marked, string property)
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
    int i_weight1     = (this->type_of_object!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
    int i_weight2     = (this->type_of_object!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
    int i_weight3     = (this->type_of_object!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
    int i_weight4     = (this->type_of_object!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
    bool use_weight1  = (this->type_of_object!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
    bool use_weight2  = (this->type_of_object!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
    bool use_weight3  = (this->type_of_object!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
    bool use_weight4  = (this->type_of_object!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
  #endif
    ULONG nlines= this->Halo.size();
    real_prec nbar=static_cast<real_prec>(this->NOBJS)/pow(this->params._Lbox(),3);
    ULONG n_selected=0;
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
      mark.resize(this->NOBJS,1.0);
      if  (property=="_MASS_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].mass);
      else if(property=="_VMAX_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].vmax);
      else if  (property=="_RS_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].rs;
      else if  (property=="_RVIR_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].rvir;
      else if  (property=="_CONCENTRATION_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=log10(this->Halo[i].concentration);
       else if  (property=="_SPIN_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].spin);
      else if  (property=="_SPIN_BULLOCK_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=log10(this->Halo[i].spin_bullock);
       else if  (property=="_VRMS_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].vrms);
       else if (property=="_VIRIAL_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].virial;
       else if  (property=="_BTOA_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].b_to_a;
       else if  (property=="_CTOA_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].c_to_a;
      else if  (property=="_MACH_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].mach_number;
      else if  (property=="_BIAS_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].bias;
      else if  (property=="_TIDAL_ANISOTROPY_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].tidal_anisotropy;
      else if  (property=="_LOCAL_OVERDENSITY_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].local_overdensity;
      else if  (property=="_PEAK_HEIGHT_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].peak_height;
      else if  (property=="_DACH_")
       for(ULONG i=0;i<mark.size();++i)
         mark[i]=this->Halo[i].dach_number;
   }
   int exp_mas=this->params._mass_assignment();
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
     for (ULONG i=0;i< nlines ;++i)
       {
            real_prec x=this->Halo[i].coord1;
            real_prec y=this->Halo[i].coord2;
            real_prec z=this->Halo[i].coord3;
#ifdef _USE_MASS_AS_OBSERVABLE_
            real_prec property=this->Halo[i].mass;
  #elif defined _USE_VMAX_AS_OBSERVABLE_
            real_prec property=this->Halo[i].vmax;
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
            if(property>=MINIMUM_PROP_CUT)
              {
#endif
#ifdef _MASS_WEIGHT_POWER_
                real_prec mass=this->Halo[i].mass;
#endif
#ifdef _USE_REDSHIFT_BINS_
                if(true==this->Halo[i].observed)
                  {
#endif
                    double nbar=static_cast<double>(mean_density);
                    if(true==this->params._use_random_catalog())
                      nbar=this->Halo[i].mean_density;
                    double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
                    ow[0] = (use_weight1 && (i_weight1<n_columns))? this->Halo[i].weight1 : 1.0;
                    ow[1] = (use_weight2 && (i_weight2<n_columns))? this->Halo[i].weight2 : 1.0;
                    ow[2] = (use_weight3 && (i_weight3<n_columns))? this->Halo[i].weight3 : 1.0;
                    ow[3] = (use_weight4 && (i_weight4<n_columns))? this->Halo[i].weight4 : 1.0;
                    ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
                    double we_fkp=1.0;
                    if(true==this->params._FKP_weight())
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
           case(I_NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
            grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
         grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
               grid_assignment_TSC(&this->params, x, y, z, ptotal_weight, field);
         break;
           case(I_PCS):
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
 void Catalog::get_interpolated_density_field_real_and_redshift_space(bool marked, string property)
  {
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_WEIGHTS_IN_POWER_
    int i_weight1     = (this->type_of_object!="RANDOM"? this->params._i_weight1_g() : this->params._i_weight1_r());
    int i_weight2     = (this->type_of_object!="RANDOM"? this->params._i_weight2_g() : this->params._i_weight2_r());
    int i_weight3     = (this->type_of_object!="RANDOM"? this->params._i_weight3_g() : this->params._i_weight3_r());
    int i_weight4     = (this->type_of_object!="RANDOM"? this->params._i_weight4_g() : this->params._i_weight4_r());
    bool use_weight1  = (this->type_of_object!="RANDOM"? this->params._use_weight1_g() : this->params._use_weight1_r());
    bool use_weight2  = (this->type_of_object!="RANDOM"? this->params._use_weight2_g() : this->params._use_weight2_r());
    bool use_weight3  = (this->type_of_object!="RANDOM"? this->params._use_weight3_g() : this->params._use_weight3_r());
    bool use_weight4  = (this->type_of_object!="RANDOM"? this->params._use_weight4_g() : this->params._use_weight4_r());
#endif
    ULONG nlines= this->Halo.size();
    ULONG n_selected=0;
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
       mark.resize(nlines,1.0);
       if  (property=="_MASS_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].mass);
       else if(property=="_VMAX_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].vmax);
       else if  (property=="_RS_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].rs;
       else if  (property=="_RVIR_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].rvir;
       else if  (property=="_CONCENTRATION_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].concentration);
        else if  (property=="_SPIN_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].spin);
       else if  (property=="_SPIN_BULLOCK_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=log10(this->Halo[i].spin_bullock);
        else if  (property=="_VRMS_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=log10(this->Halo[i].vrms);
        else if (property=="_VIRIAL_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].virial;
        else if  (property=="_BTOA_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].b_to_a;
        else if  (property=="_CTOA_")
         for(ULONG i=0;i<mark.size();++i)
           mark[i]=this->Halo[i].c_to_a;
       else if  (property=="_MACH_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].mach_number;
       else if  (property=="_BIAS_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].bias;
       else if  (property=="_TIDAL_ANISOTROPY_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].tidal_anisotropy;
       else if  (property=="_LOCAL_OVERDENSITY_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].local_overdensity;
       else if  (property=="_PEAK_HEIGHT_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].peak_height;
       else if  (property=="_DACH_")
        for(ULONG i=0;i<mark.size();++i)
          mark[i]=this->Halo[i].dach_number;
    }
   int exp_mas=this->params._mass_assignment();// default, NGP
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
        for (ULONG i=0;i< nlines ;++i)
          {
            real_prec x=this->Halo[i].coord1;
            real_prec y=this->Halo[i].coord2;
            real_prec z=this->Halo[i].coord3;
            real_prec vx=rsd_x*this->Halo[i].vel1*conversion_factor;
            real_prec vy=rsd_y*this->Halo[i].vel2*conversion_factor;
            real_prec vz=rsd_z*this->Halo[i].vel3*conversion_factor;
            real_prec xs=x+vx;
            real_prec ys=y+vy;
            real_prec zs=z+vz;
#ifdef _USE_MASS_AS_OBSERVABLE_
            real_prec property=this->Halo[i].mass;
#elif defined _USE_VMAX_AS_OBSERVABLE_
            real_prec property=this->Halo[i].vmax;
#endif
#ifdef _SET_GLOBAL_MASS_CUT_
            if(property>=MINIMUM_PROP_CUT)
              {
#endif
#ifdef _MASS_WEIGHT_POWER_
                real_prec mass=this->Halo[i].mass;
#endif
#ifdef _USE_REDSHIFT_BINS_
                if(true==this->Halo[i].observed)
                  {
#endif
                    double nbar=static_cast<double>(mean_density);
                    if(true==this->params._use_random_catalog())
                      nbar=this->Halo[i].mean_density;
                    double ptotal_weight = 1.0;
#ifdef _USE_WEIGHTS_IN_POWER_
                    ow[0] = (use_weight1 && (i_weight1<n_columns))? this->Halo[i].weight1 : 1.0;
                    ow[1] = (use_weight2 && (i_weight2<n_columns))? this->Halo[i].weight2 : 1.0;
                    ow[2] = (use_weight3 && (i_weight3<n_columns))? this->Halo[i].weight3 : 1.0;
                    ow[3] = (use_weight4 && (i_weight4<n_columns))? this->Halo[i].weight4 : 1.0;
                    ptotal_weight = ow[0] * ow[1] * ow[2] * ow[3];
#endif
                    double we_fkp=1.0;
                    if(true==this->params._FKP_weight())
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
                    case(I_NGP):
#ifdef _USE_OMP_
#pragma atomic
#endif
                     grid_assignment_NGP(&this->params, x, y, z, ptotal_weight, field);
                     grid_assignment_NGP(&this->params, xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(I_CIC):
#ifdef _USE_OMP_
#pragma atomic
#endif
                        grid_assignment_CIC(&this->params, x, y, z, ptotal_weight, field);
                        grid_assignment_CIC(&this->params, xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(I_TSC):
#ifdef _USE_OMP_
#pragma atomic
#endif
                        grid_assignment_TSC(&this->params,x, y, z, ptotal_weight, field);
                        grid_assignment_TSC(&this->params,xs, ys, zs, ptotal_weight, field_s);
                  break;
                    case(I_PCS):
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
void Catalog::get_mean_number_density(real_prec alpha0, vector<gsl_real> &zz, vector<gsl_real>&rc, vector<gsl_real> &new_zz_v, vector<gsl_real>&new_dndz_v)
  {
   // Get the mean number density.
    So.enter(__PRETTY_FUNCTION__);
    this->So.message_screen("Preparing for interpolation of mean number density");
    ULONG nlines= this->Halo.size();
    ULONG nz = zz.size();
    ULONG n_dndz = new_zz_v.size();
    real_prec redshift_min_sample=this->params._redshift_min_sample();
    real_prec redshift_max_sample=this->params._redshift_max_sample();
    real_prec area_survey=this->params._area_survey();
    int sys_of_coord_r=this->params._sys_of_coord_r(); // Always refere to the random catalog, as these are smoothed
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
    ULONG n_objects = this->NOBJS;
    real_prec area=area_survey*pow(M_PI/180.,2); /*converting to strad*/
    for(int i=0;i<n_dndz;++i)
      new_zz_v[i]=redshift_min_sample+(i+0.5)*Delta_Z;
    if(I_EQR==sys_of_coord_r)
      { // We must transform r to z
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
     for(int i=0;i<nlines;i++)
       {
         real_prec zro= gsl_inter_new(rc, zz, this->Halo[i].coord3);
         int count= get_bin(zro, redshift_min_sample,n_dndz,Delta_Z,true); // this is not optimal
  #ifdef _USE_OMP_
  #pragma omp atomic
  #endif
          new_dndz_v[count]++;
        }
      }
     else  if(I_EQZ==sys_of_coord_r)
      {
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
       for(int i=0;i<this->NOBJS;i++)
        {
         int count= get_bin(this->Halo[i].coord3, redshift_min_sample,n_dndz,Delta_Z,true);
  #ifdef _USE_OMP_
  #pragma omp atomic update
  #endif
         new_dndz_v[count]++;
      }
    }
    ULONG n_rc = rc.size();
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
    for(int i=0;i<n_dndz;i++)
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
void Catalog::get_random_catalog(){
    So.enter(__PRETTY_FUNCTION__);
    // FIrst get the distribuytion in z, color and stellar mass:
   try{
        if(this->params._sys_of_coord_g()>0)
            So.message_screen("");
        else
            throw (this->params._sys_of_coord_g());
    }
    catch (int soc){
        So.message_screen("System of coordinate not valid for this particular task");
        exit(1);
    }
    this->get_cosmo();
    ULONG Nbins=this->params._Nbins_color()*this->params._Nbins_Mstellar()*this->params._Nbins_redshift();
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
    for(int i=0;i<Pzdist.size();i++)
        zbins[i]=static_cast<gsl_real>(this->params._redshift_min_sample()+(i+0.5)*delta_redshift);
    for(int i=0;i<Pcdist.size();i++)
        cbins[i]=this->params._Color_min()+(i+0.5)*delta_color;
    for(int i=0;i<Pmdist.size();i++)
        mbins[i]=this->params._Mstellar_min()+(i+0.5)*delta_Mstellar;
    //So.message_screen("\tObtaining distributions in galaxy properties with a magnitud cut at m=", this->params._mK_max());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<this->NOBJS;++i)
    {
  //      if(this->Halo[i].app_mag <this->params._mK_max())
  //      {
        int ibin_z = get_bin(this->Halo[i].redshift,this->params._redshift_min_sample(),this->params._Nbins_redshift(),delta_redshift,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pzdist[ibin_z]++;
        int ibin_color = get_bin(this->Halo[i].color,this->params._Color_min(),this->params._Nbins_color(),delta_color,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pcdist[ibin_color]++;
        int ibin_Mstellar = get_bin(this->Halo[i].stellar_mass,this->params._Mstellar_min(),this->params._Nbins_Mstellar(),delta_Mstellar,true);
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        Pmdist[ibin_Mstellar]++;
        ULONG index3d=index_3d(ibin_z, ibin_color, ibin_Mstellar, this->params._Nbins_color(), this->params._Nbins_Mstellar());
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
              ULONG indexa=index_3d(iz, ic, im, this->params._Nbins_color(), this->params._Nbins_Mstellar());
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
    for(int i=0;i<Pzdist.size();i++)
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
    for(int i=0;i<nbar.size();++i)
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
    this->So.message_screen("Generating random catalog with ", this->NOBJS*Nran_files, " objects in chuncks of ", Nran_files, " files");
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
    ULONG n_random=this->NOBJS;
#ifdef _USE_OMP_
    int jthread;
    omp_set_num_threads(Nran_files);
    vector<ULONG>vseeds(Nran_files,0);
    for(int i=0;i<vseeds.size();++i)
      vseeds[i]=3+static_cast<ULONG>(i+27*i*i)*56145;
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
        ULONG count=0;
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
                           ULONG id=index_3d(iz, indexc, indexm,this->params._Nbins_color(), this->params._Nbins_Mstellar());
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
void Catalog::get_intervals_multiscale(string prop)
{
    ULONG Nft=this->params._Nft();
    gsl_vector *prop_aux=gsl_vector_alloc(this->NOBJS);
    if(prop=="_VMAX_")
      for(ULONG ig=0;ig<this->NOBJS;++ig )
         gsl_vector_set(prop_aux,ig,log10(this->Halo[ig].vmax));
    if(prop=="_MASS_")
      for(ULONG ig=0;ig<this->NOBJS;++ig )
         gsl_vector_set(prop_aux,ig,log10(this->Halo[ig].mass));
    gsl_sort_vector(prop_aux) ;   // sort the vmax and correspondingly their associated the gal id
    ULONG counter=0;
    for(int il=0;il<this->params._Number_of_MultiLevels();++il)
     {
        ULONG N_level=pow(this->params.get_Nft_MultiLevels(il),3);
        this->params.set_PropThreshold_MultiLevels(il,pow(10,gsl_vector_get(prop_aux,this->NOBJS-(N_level+counter))));
        this->params.set_Ntracers_MultiLevels(il,N_level);
        So.message_screen("Number of reference properties to be assigned at individual level=",N_level);
        So.message_screen("Minimum property value in level",pow(10,gsl_vector_get(prop_aux,this->NOBJS-(N_level+counter))));
        counter=N_level;
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_mean_bias_relation(s_info_in_bins &bias_info)
{
  So.enter(__PRETTY_FUNCTION__);
  int Nbins=bias_info.s_size();
  vector<int>ibin(this->Halo.size(),0);
  vector<real_prec>prop(this->Halo.size(),0);
  real_prec min_p=0;
  real_prec max_p=0;
  if(bias_info.name_info==_VMAX_)
   {
      min_p=log10(this->params._VMAXmin());
      max_p=log10(this->params._VMAXmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
        prop[i]=log10(this->Halo[i].vmax);
   }
  else if (bias_info.name_info==_MASS_)
   {
      min_p=this->params._LOGMASSmin();
      max_p=this->params._LOGMASSmax();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG i=0;i<this->Halo.size();++i)
        prop[i]=log10(this->Halo[i].mass);
   }
  else if(bias_info.name_info==_SPIN_ || bias_info.name_info==_SPIN_BULLOCK_)
   {
      min_p=log10(this->params._SPINmin());
      max_p=log10(this->params._SPINmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
        prop[i]=log10(this->Halo[i].spin_bullock);
   }
  else if(bias_info.name_info==_RS_ )
   {
      min_p=log10(this->params._RSmin());
      max_p=log10(this->params._RSmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
        prop[i]=log10(this->Halo[i].rs);
   }  
  else if(bias_info.name_info==_CONCENTRATION_ )
   {
      min_p=log10(this->params._CONCENTRATIONmin());
      max_p=log10(this->params._CONCENTRATIONmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
        prop[i]=log10(this->Halo[i].concentration);
   }  
  real_prec delta=(max_p-min_p)/static_cast<real_prec>(Nbins);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for (ULONG i=0;i<Nbins;++i)
     bias_info.vbin[i]=min_p+(i+0.5)*delta;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->Halo.size();++i)
    if(prop[i]<max_p && prop[i]>= min_p)
      {
        int bin=get_bin(prop[i],min_p,Nbins,delta,true);
        ibin[i]=bin;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.vq1[bin]+=this->Halo[i].bias;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.i_ncbin[bin]++;
     }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nbins;++i)  // mean
    if(bias_info.i_ncbin[i]>0)  
      bias_info.vq1[i]=bias_info.vq1[i]/static_cast<real_prec>(bias_info.i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->Halo.size();++i)
    if(prop[i]<max_p && prop[i]>= min_p)
      {
        int bin=ibin[i];
#ifdef _USE_OMP_
#pragma omp atomic
#endif
        bias_info.vq2[bin]+=pow(this->Halo[i].bias-bias_info.vq1[bin],2);
    }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<Nbins;++i)  // standard deviation
    if(bias_info.i_ncbin[i]>0)  
     bias_info.vq2[i]=sqrt(bias_info.vq2[i])/static_cast<real_prec>(bias_info.i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (int i=0;i<Nbins;++i) // error in the mean
    if(bias_info.i_ncbin[i]>0)  
      bias_info.vq3[i]=bias_info.vq2[i]/sqrt(static_cast<real_prec>(bias_info.i_ncbin[i]));
So.leaving(__PRETTY_FUNCTION__);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalog::get_mean_secondary_bias_relation(vector<s_info_in_bins> &sec_bias_info)
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
  vector<int>ibin(this->Halo.size(),0);
  vector<bool>used(this->Halo.size(),false);
  vector<real_prec>p_prop(this->Halo.size(),0);
// ----------------------------------------------------------
// Get min and max of primeary proeprty, read fom parameter file:
  real_prec min_p, max_p, min_s, max_s;
  if(primary_prop==_VMAX_)
   {
     min_p=log10(this->params._VMAXmin());
     max_p=log10(this->params._VMAXmax());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<this->Halo.size();++i)
     p_prop[i]=log10(this->Halo[i].vmax);
   }
  else if (primary_prop==_MASS_)
   {
     min_p=this->params._LOGMASSmin();
     max_p=this->params._LOGMASSmax();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for (ULONG i=0;i<this->Halo.size();++i)
       p_prop[i]=log10(this->Halo[i].mass);
  }
     real_prec delta=(max_p-min_p)/static_cast<real_prec>(Nbins);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<Nbins;++i)
      sec_bias_info[0].vbin[i]=min_p+(i+0.5)*delta;
// ----------------------------------------------------------
  vector<real_prec>s_prop(this->Halo.size(),0);
  if (secondary_prop ==_SPIN_ || secondary_prop==_SPIN_BULLOCK_)
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for (ULONG i=0;i<this->Halo.size();++i)
       s_prop[i]=log10(this->Halo[i].spin_bullock);
   }
  else if (secondary_prop==_RS_ )
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<this->Halo.size();++i)
      s_prop[i]=log10(this->Halo[i].rs);
   } 
  else if (secondary_prop==_CONCENTRATION_ )
   {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for (ULONG i=0;i<this->Halo.size();++i)
      s_prop[i]=log10(this->Halo[i].concentration);
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
    for (ULONG i=0;i<this->Halo.size();++i)
      if(p_prop[i]<max_p && p_prop[i]>= min_p)
       if(s_prop[i]<max_s && s_prop[i]>= min_s)
        {
          int bin=get_bin(p_prop[i],min_p,Nbins,delta,true);
          ibin[i]=bin;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          sec_bias_info[quart].vq1[bin]+=this->Halo[i].bias;
#ifdef _USE_OMP_
#pragma omp atomic
#endif
          sec_bias_info[quart].i_ncbin[bin]++;
        }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (int i=0;i<Nbins;++i)  // mean
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq1[i]=sec_bias_info[quart].vq1[i]/static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<this->Halo.size();++i)
        if(p_prop[i]<max_p && p_prop[i]>= min_p)
          if(s_prop[i]<max_s && s_prop[i]>= min_s)
            {
              int bin=ibin[i];
#ifdef _USE_OMP_
#pragma omp atomic
#endif
              sec_bias_info[quart].vq2[bin]+=pow(this->Halo[i].bias-sec_bias_info[quart].vq1[bin],2);
          }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (int i=0;i<Nbins;++i)  // standard deviation
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq2[i]=sqrt(sec_bias_info[quart].vq2[i])/static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (int i=0;i<Nbins;++i) // error in the mean
        if(sec_bias_info[quart].i_ncbin[i]>0)  
          sec_bias_info[quart].vq3[i]=sec_bias_info[quart].vq2[i]/sqrt(static_cast<real_prec>(sec_bias_info[quart].i_ncbin[i]));
 }
So.leaving(__PRETTY_FUNCTION__);
}