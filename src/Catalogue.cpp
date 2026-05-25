///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<Catalogue>
 * @file Catalogue.cpp
 * @brief Methods of the class Catalog
 * @details The class Catalogue reads and input catalog of dark matter tracers, and allocate mins, max and mean of each property.
 * The properties are listed as parameters in the json file, associated to Tracer section.
 * The typle with min, max and mean are public class members.
 * @author Andres Balaguera Antolinez 2007-2026
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "../include/def.hpp"
#include "../include/Catalogue.hpp"
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



template<typename T> void Catalogue::load_column(std::vector<T>& target, const T* prop, int column_index, std::tuple<T,T,T>& information)
{
  if(column_index < 0)
   return;
  else
  {
    
    target.resize(this->NOBJS);
    T minv=static_cast<T>(1e20);
    T maxv=static_cast<T>(-1e20);
    T meanv=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+: meanv) reduction(min:minv) reduction(max:maxv)
#endif
    for(size_t i = 0; i < this->NOBJS; ++i)
    {
        T val = static_cast<T>(prop[column_index + i * this->NCOLS]);
        target[i] = val;
        meanv+=val;
        minv=min(val,minv);
        maxv=max(val,maxv);
      }
      meanv/=static_cast<T>(this->NOBJS);
      information=std::make_tuple(minv, maxv, meanv);

  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalogue::write_catalog_bin(string outputFileName)
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
       float x=static_cast<float>(coord1[i]);
       float y=static_cast<float>(coord2[i]);
       float z=static_cast<float>(coord3[i]);
       outStream.write(reinterpret_cast<char*>(&x), sizeof(x));
       outStream.write(reinterpret_cast<char*>(&y), sizeof(y));
       outStream.write(reinterpret_cast<char*>(&z), sizeof(z));
#endif
#ifdef _WRITE_VELOCITIES_
       float vx=static_cast<float>(vel1[i]*conversion_factor);
       float vy=static_cast<float>(vel2[i]*conversion_factor);
       float vz=static_cast<float>(vel3[i]*conversion_factor);
       outStream.write(reinterpret_cast<char*>(&vx), sizeof(vx));
       outStream.write(reinterpret_cast<char*>(&vy), sizeof(vy));
       outStream.write(reinterpret_cast<char*>(&vz), sizeof(vz));
#endif

#if defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ || defined (_ASSIGN_MASS_POST_)
       float fmass=static_cast<float>(mass[i]);
       outStream.write(reinterpret_cast<char*>(&fmass), sizeof(mass));
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
       float fvmax=static_cast<float>(vmax[i]);
       outStream.write(reinterpret_cast<char*>(&fvmax), sizeof(vmax));
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
        float frs=static_cast<float>(rs[i]);
        outStream.write(reinterpret_cast<char*>(&frs), sizeof(rs));
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
       float fspin=static_cast<float>(spin[i]);
       outStream.write(reinterpret_cast<char*>(&fspin), sizeof(spin));
#endif
    }
    outStream.close();
    So.DONE();
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalogue::write_catalog(string outputFileName)
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
    float x=static_cast<float>(coord1[i]);
    float y=static_cast<float>(coord2[i]);
    float z=static_cast<float>(coord3[i]);
    outStream.write(reinterpret_cast<char*>(&x), sizeof(x));
    outStream.write(reinterpret_cast<char*>(&y), sizeof(y));
    outStream.write(reinterpret_cast<char*>(&z), sizeof(z));
#endif
#ifdef _WRITE_VELOCITIES_
    float vx=static_cast<float>(vel1[i]*conversion_factor);
    float vy=static_cast<float>(vel2[i]*conversion_factor);
    float vz=static_cast<float>(vel3[i]*conversion_factor);
    outStream.write(reinterpret_cast<char*>(&vx), sizeof(vx));
    outStream.write(reinterpret_cast<char*>(&vy), sizeof(vy));
    outStream.write(reinterpret_cast<char*>(&vz), sizeof(vz));
#endif
#if defined _USE_MASS_AS_PRIMARY_OBSERVABLE_ || defined (_ASSIGN_MASS_POST_)
    float fmass=static_cast<float>( mass[i]);
    outStream.write(reinterpret_cast<char*>(&mfass), sizeof(fmass));
#endif
#ifdef _USE_VMAX_AS_PRIMARY_OBSERVABLE_
    float fvmax=static_cast<float>(vmax[i]);
    outStream.write(reinterpret_cast<char*>(&fvmax), sizeof(fvmax));
#endif
#ifdef _USE_RS_AS_DERIVED_OBSERVABLE_
    float frs=static_cast<float>(rs[i]);
        outStream.write(reinterpret_cast<char*>(&frs), sizeof(frs));
#endif
#ifdef _USE_SPIN_AS_DERIVED_OBSERVABLE_
    float fspin=static_cast<float>(fspin[i]);
    outStream.write(reinterpret_cast<char*>(&fspin), sizeof(fspin));
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
        if( observed[i]==1)
#if defined _ASSIGN_PROPERTIES_
        outStream<< coord1<<"\t"<< coord2<<"\t"<< coord3<<"\t"<< vel1*conversion_factor<<"\t"<< vel2*conversion_factor<<"\t"<< vel3*conversion_factor<<"\t"<< mass<<"\t"<< vmax<<"\t"<< rs<<"\t"<< spin<<"\t"<< identity<<"\t"<< gal_cwt<<"\t"<< local_dm<<endl;
#elif defined _WRITE_COORDINATES_ || defined _WRITE_VELOCITIES_
        outStream<< redshift[i]<<"\t"<< mass[i]<<"\t"<< color[i]<<"\t"<< stellar_mass[i]<<"\t"<< relative_bias[i]<<endl;
#endif
    }
    outStream.close();
#ifdef _VERBOSE_CAT_
    So.DONE();
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Catalogue::read_catalog_bin_tng()
{
#ifdef _VERBOSE_CATALOG_
  this->So.enter(__PRETTY_FUNCTION__);
#endif
    coord1.resize(this->params._N_lines_binary(),0);
    coord2.resize(this->params._N_lines_binary(),0);
    coord3.resize(this->params._N_lines_binary(),0);
    this->File.read_array(this->params._file_bin_x_coord(),&coord1[0],this->params._N_lines_binary());
    this->File.read_array(this->params._file_bin_y_coord(),&coord2[0],this->params._N_lines_binary());
    this->File.read_array(this->params._file_bin_z_coord(),&coord3[0],this->params._N_lines_binary());
    this->NOBJS=coord1.size();

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef STAYHERE

#ifdef _USE_VELOCITIES_
void Catalogue::read_catalog_bin(ULONG n, string filex, string filey, string filez,string filevx, string filevy, string filevz)
#else
void Catalogue::read_catalog_bin()
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
  if(this->params._dilute_dm_sample())
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
 if(!this->params._use_low_pass_filter())
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
 if(this->params._use_low_pass_filter())
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


#endif  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// In this function a catalog of DM or DM tracers, in ascii format, is read and their proeprtiess read, according to the request of the parameter file.
void Catalogue::read_catalog_new(std::string input_file)
{

  this->So.enter(__PRETTY_FUNCTION__);
  
  using ColumnTuple = std::tuple<std::vector<real_prec>*, int, std::tuple<real_prec, real_prec, real_prec>* >;

  int i_coord1=-1;
  int i_coord2=-1;
  int i_coord3=-1;
  int i_vel1=-1;
  int i_vel2=-1;
  int i_vel3=-1;
  int i_mass=-1; 
  int i_weight=-1; 
  int i_mean_density=-1; 
  int i_sf=-1; 
  int i_vmax=-1; 
  int i_rs=-1;
  int i_rvir=-1; 
  int i_virial=-1; 
  int i_spin=-1;
  int i_spin_bullock=-1;
  int i_vrms=-1;
  int i_b_to_a=-1;
  int i_c_to_a=-1;
  int i_redshift=-1;
  int i_stellar_mass=-1;
  int i_color=-1;
  int i_abs_mag=-1;
  int i_app_mag=-1;
  CoordinateSystem systemc;

 if(this->type_of_object=="TRACER" || this->type_of_object=="TRACER_REF" || this->type_of_object=="TRACER_MOCK" || this->type_of_object=="TRACER_MOCK_ONLY_COORDS"|| this->type_of_object=="TRACER_REF_ONLY_COORDS")
    {
      i_coord1= this->params._i_coord1_g();
      i_coord2= this->params._i_coord2_g();
      i_coord3= this->params._i_coord3_g();
      i_vel1= this->params._i_v1_g();
      i_vel2= this->params._i_v2_g();
      i_vel3= this->params._i_v3_g();
      i_mass = this->params._i_mass_g();
      i_vmax = this->params._i_vmax_g();
      i_vrms = this->params._i_vrms_g();
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
      i_weight = this->params._i_weight1_g();
      systemc = this->params._sys_of_coord_g();
      if(systemc ==CoordinateSystem::EQZ)
        i_redshift= this->params._i_coord3_g();// note that this applies only in the case in which coordsa are in pseudo-equatorial, i.e, ra, dec, z.
    }
  else if(this->type_of_object=="RANDOM")
    {
      i_coord1= this->params._i_coord1_r();
      i_coord2= this->params._i_coord2_r();
      i_coord3= this->params._i_coord3_r();
      i_vel1= -1; 
      i_vel2= -1; 
      i_vel3= -1; 
      i_mass = -1; 
      i_vmax = -1; 
      i_vrms = -1; 
      i_weight = this->params._i_weight1_r();
      i_mean_density= this->params._i_mean_density_r();
      i_sf = -1;
      i_rs = 
      i_rvir = -1;
      i_virial = -1;
      i_spin = -1;
      i_spin_bullock = -1;
      i_b_to_a = -1;
      i_c_to_a = -1;
      i_stellar_mass = this->params._i_stellar_mass_r();
      i_color = this->params._i_color_r();
      i_app_mag = this->params._i_app_mag_r();
      i_abs_mag = this->params._i_abs_mag_r();
      i_redshift= -1;
      systemc = this->params._sys_of_coord_r();
      if(this->params._sys_of_coord_r()==CoordinateSystem::EQZ)
        i_redshift= this->params._i_coord3_r();// note that this applies only in the case in which coordsa are in pseudo*-equatorial, i.e, ra, dec, z.
    }
 
      // List of avaibale entries for proeprties of tracer catalogue.
      std::vector<ColumnTuple> columns_to_load = {
      {&this->coord1, i_coord1, &this->info_coord1},
      {&this->coord2, i_coord2, &this->info_coord2},
      {&this->coord3, i_coord3, &this->info_coord3},
      {&this->vel1,   i_vel1, &this->info_vel1},
      {&this->vel2,   i_vel2, &this->info_vel2},
      {&this->vel3,   i_vel3, &this->info_vel3},
      {&this->mass,   i_mass, &this->info_mass},
      {&this->vmax,   i_vmax, &this->info_vmax},
      {&this->vrms,   i_vrms, &this->info_vrms},
      {&this->mean_density,   i_mean_density, &this->info_mean_density},
      {&this->rs,   i_rs,     &this->info_rs},
      {&this->rvir,   i_rvir,  &this->info_rvir},
      {&this->spin,   i_spin, &this->info_spin},
      {&this->spin_bullock,   i_spin_bullock, &this->info_spin_bullock},
      {&this->b_to_a, i_b_to_a, &this->info_b_to_a},
      {&this->c_to_a, i_c_to_a, &this->info_c_to_a},
      {&this->stellar_mass, i_stellar_mass,&this->info_stellar_mass},
      {&this->color, i_color ,&this->info_color},
      {&this->app_mag, i_app_mag ,&this->info_app_mag},
      {&this->abs_mag, i_abs_mag,&this->info_abs_mag},
      {&this->redshift, i_redshift,&this->info_redshift},
      {&this->weight1, i_weight, &this->info_weight1}
    };

   
  vector<real_prec> prop;
  this->NOBJS = this->File.read_file(input_file, prop,_NTHREADS_);
  this->NCOLS=(static_cast<ULONG>(prop.size()/this->NOBJS));
  
  this->So.message_screen("Loading data");

  for(auto& [vec_ptr, col_index, information] : columns_to_load)
      load_column(*vec_ptr, prop.data(), col_index, *information);
  So.DONE();  




  prop.clear(); prop.shrink_to_fit();   

//      this->So.message_screen("\tMin redshift  =", std::get<static_cast<int>(PropStats::MIN)>(this->info_redshift));
//      this->So.message_screen("\tMax redshift  =", std::get<static_cast<int>(PropStats::MAX)>(this->info_redshift));
//      this->So.message_screen("\tMean redshift  =", std::get<static_cast<int>(PropStats::MEAN)>(this->info_redshift));


  if(i_rvir>0 && i_rs>0)
  {
    this->concentration.resize(NOBJS);
    So.message_screen("Assigning concentration");
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NOBJS;++i)
      {
        if(rs[i]>0)
          this->concentration[i]=this->rvir[i]/this->rs[i];
        else
          this->concentration[i]=0;
        }
        So.DONE();
      }


  if(i_b_to_a>0 && i_c_to_a >0)
  {
    if(params._convert_to_halo_ellipiticity_and_prolatness())
    {
  #ifdef _USE_OMP_
  #pragma omp parallel for
  #endif
      for(ULONG i=0;i<NOBJS;++i)
      {
          real_prec s=this->c_to_a[i];
          real_prec q=this->b_to_a[i];
          real_prec prol=(1.-q*q)/(1.-s*s);
          real_prec ell=(1.-s*s)/(1.+s*s+q*q);
          this->c_to_a[i]=ell;
          this->b_to_a[i]=prol;
      }
    }
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 void Catalogue::print_catalogue(string file){
 
   So.enter(__PRETTY_FUNCTION__);

   ofstream rcat;
   rcat.open(file.c_str());
   rcat.precision(3);
   rcat.setf(ios::showpoint);
   rcat.setf(ios::scientific);

   ULONG counter=0;
  
   this->So.message_screen("Writting downsampled version of input catalogue with bias in file ",file);

   if(true==this->params._get_cwc_properties())
   {
      if(true==this->params._Get_tracer_bias_multipoles())
        {
          rcat<<"#coord1, coord2, coord3, logM, log rs, log c, log spin, tidal_ani_dm, tidal_ani_tr, bias,bias_lm"<<endl;
          for (size_t i = 0; i < this->NOBJS; ++i)
           {
               rcat << coord1[i]<<"\t"
                <<coord2[i]<<"\t"
                <<coord3[i]<<"\t"
                <<log10(mass[i]) << "\t"
                << log10(rs[i]) << "\t"
                << log10(concentration[i]) << "\t"
                << log10(spin_bullock[i]) << "\t"
                << tidal_anisotropy_dm[i] << "\t"
                << tidal_anisotropy[i] << "\t"
                << bias[i] << "\t"
                << static_cast<int>(gal_cwt[i]) << "\t";
                for (size_t il = 0; il < bias_multipole.size(); ++il)
                  rcat << bias_multipole[i][il] << "\t";
                rcat << std::endl;
          }
        }
    else 
    {
      rcat<<"#coord1, coord2, coord3, logM, log rs, log c, log spin, tidal_anisotropy_dm, tidal_anisotropy_tr, bias, cwc"<<endl;
      for (size_t i = 0; i < this->NOBJS; ++i) {
        rcat << coord1[i]<<"\t"
              <<coord2[i]<<"\t"
              <<coord3[i]<<"\t"
              <<log10(mass[i]) << "\t"
              << log10(rs[i]) << "\t"
              << log10(concentration[i]) << "\t"
              << log10(spin_bullock[i]) << "\t"
              << tidal_anisotropy_dm[i] << "\t"
              << tidal_anisotropy[i] << "\t"
              << bias[i] << "\t"
              << static_cast<int>(gal_cwt[i]) << endl;
        }
  
    }
  }
  else
   {
    if(this->params._Get_tracer_bias_multipoles())
      {
        rcat<<"#coord1, coord2, coord3, logM, log rs, log c, log spin, bias,bias_lm"<<endl;
        for (size_t i = 0; i < this->NOBJS; ++i) {
          rcat << coord1[i]<<"\t"
              <<coord2[i]<<"\t"
              <<coord3[i]<<"\t"
              <<log10(mass[i]) << "\t"
              << log10(rs[i]) << "\t"
              << log10(concentration[i]) << "\t"
              << log10(spin_bullock[i]) << "\t"
              << bias[i] << "\t";
              for (size_t il = 0; il < bias_multipole.size(); ++il)
                rcat << bias_multipole[i][il] << "\t";
            rcat << std::endl;
        }
      }
    else
     {
        rcat<<"#coord1, coord2, coord3, logM, log rs, log c, log spin, bias"<<endl;
        for (size_t i = 0; i < this->NOBJS; ++i) {
          rcat << coord1[i]<<"\t"
              <<coord2[i]<<"\t"
              <<coord3[i]<<"\t"
              <<log10(mass[i]) << "\t"
              << log10(rs[i]) << "\t"
              << log10(concentration[i]) << "\t"
              << log10(spin_bullock[i]) << "\t"
              << bias[i]<<endl;
        }
   }
  }
  rcat.close();
  So.DONE();
}
