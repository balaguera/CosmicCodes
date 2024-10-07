////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
#include "../headers/def.h"
# include "CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
//  system("clear");
  time_t start_all;
  time(&start_all);
  char temp;
  string par_file_BiasMT;

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec everywhere. Plese desable and compile again.");
#endif
  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);
  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }
  while((temp =  getopt(argc, argv, "hadi:c:h:b:u:s:m:")) != -1)
    {
      if(temp=='h')
        So.usage(argv[0]);

      else if (temp=='a')  // Show authors
        So.author();

      else if(temp=='i') // displays input parameters
        {
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          params.show_params();
        }
      else if(temp=='d') // displays input parameters
        {
          So.show_preproc();
        }
      else if(temp=='c')   // to read input tracer catalg and analyze it. 
        {

          So.message_screen("cosmicatlass running under the option -c");
#ifdef _REDSHIFT_SPACE_
          So.message_warning("Option _REDSHIFT_SPACE_ in def.h is enabled. If velocities are not to be used, please desable it to save space. Code ends here.");
          So.message_screen("After editing the def.h file, please do make clean; make -j BiasMT;");
#endif
          So.message(start_all);
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog catalog(params); // this is redundant
          catalog.analyze_cat(true);
          So.message_time(start_all);
        }
      else if(temp=='h')   // Read halo catalog and generate galaxy catalogs using an HOD
        {
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog tracer(params);
          tracer.halos2galaxies_HOD();
        }
      else if('b'==temp)  // assign bias to tracers
        {
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog cat(params);
          cat.read_catalog_bin_tng();
          vector<double>dm_field(params._NGRID(),0);
          File.read_array_t<double>(params._Input_Directory_X()+params._Name_Catalog_X(), dm_field);
//          for (ULONG i = 0; i < dm_field.size(); i++)if(isinf(dm_field[i]))cout<<i<<"  "<<dm_field[i]<<endl;
          vector<real_prec>dm_field_f(params._NGRID(),0);
          for (ULONG i = 0; i < dm_field.size(); i++)
            dm_field_f[i]=static_cast<real_prec>(dm_field[i]);       
          for (size_t i = 0; i < 6; i++)cout<<dm_field_f[i]<<endl;
          dm_field.clear();dm_field.shrink_to_fit();
//          exit(0);
          PowerSpectrumF power(params);
          power.object_by_object_bias(cat.Halo,  dm_field_f);
      }
      else if('s'==temp)  // Conpute bias as a funciton of intrinsic properties
        {
          int Nqrts=4+1; // NUmber of qaurtiles plus the first
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog cat(params);
          cat.read_catalog(params._Input_dir_cat()+params._file_catalogue(), 0);
#ifndef _READ_BIAS_
          vector<double>dm_field(params._NGRID(),0);
          File.read_array_t<double>(params._Input_Directory_X()+params._Name_Catalog_X(), dm_field);
          vector<real_prec>dm_field_f(params._NGRID(),0);
          for (ULONG i = 0; i < dm_field.size(); i++)
            dm_field_f[i]=static_cast<real_prec>(dm_field[i]);       
          for (size_t i = 0; i < 6; i++)cout<<dm_field_f[i]<<endl;
          dm_field.clear();dm_field.shrink_to_fit();
          PowerSpectrumF power(params);
          power.object_by_object_bias(cat.Halo,  dm_field_f);
         So.message_screen("Reading bias from file file",out_bias);
#else
         string out_bias=params._Output_directory()+"individual_bias_invPhase.txt";
         ifstream bout; bout.open(out_bias.c_str());
         for(ULONG i=0;i<cat.Halo.size();++i)
          bout>>cat.Halo[i].bias>>cat.Halo[i].mach_number>>cat.Halo[i].local_overdensity;
        bout.close();
#endif 
        string primary_prop="_MASS_";
        vector<s_info_in_bins>secondary_bias_info_ref(Nqrts);  //container for seconday bias: 2  the number of quartiles to use 
        int Nbins=N_BINS_BIAS;
        s_info_in_bins bias_info;
        bias_info.allocate_all(Nbins);
        bias_info.name_info=primary_prop;
        cat.get_mean_bias_relation(bias_info);                // Get bias_prop relation from mock
        string file_o=params._Output_directory()+"halo_bias.txt";
        ofstream bal;
        So.message_screen("Writing bias in file", file_o);
        bal.open(file_o.c_str());
        for(int i=0;i<Nbins;++i)
          if(bias_info.vq1[i]>0)
            bal<<bias_info.vbin[i]<<"\t"<<bias_info.vq1[i]<<"\t"<<bias_info.vq3[i]<<endl;
        bal.close();



        string secondary_prop="_CONCENTRATION_";
        vector<s_info_in_bins>secondary_bias_info(Nqrts);  //container for seconday bias: 2  the number of quartiles to use 
        for(int i=0;i<secondary_bias_info.size();++i)
          secondary_bias_info[i].allocate_all(Nbins);  
        for(int i=0;i<secondary_bias_info.size();++i)
          secondary_bias_info[i].name_info=primary_prop;  
        for(int i=0;i<secondary_bias_info.size();++i)
          secondary_bias_info[i].name_info_sec=secondary_prop;  
        cat.get_mean_secondary_bias_relation(secondary_bias_info); // Get bias_prop relation from reference
  
        file_o=params._Output_directory()+"halo_bias_"+secondary_prop+".txt";
        for(int iq=1;iq<Nqrts;iq++)//loop over quartiles
        {
          string bbfile=file_o+"_quartile"+to_string(iq);
          ofstream sbal;
          So.message_screen("Writing bias in file", bbfile);
          sbal.open(bbfile.c_str());
          for(int i=0;i<Nbins;++i)
            if(secondary_bias_info[iq].vq1[i]>0)
              sbal<<secondary_bias_info[0].vbin[i]<<"\t"<<secondary_bias_info[iq].vq1[i]<<"\t"<<secondary_bias_info[iq].vq3[i]<<"\t"<<secondary_bias_info[iq].vq1[i]<<"\t"<<secondary_bias_info[iq].vq3[i]<<endl;
          sbal.close();
        }
 
      }
      else if('u'==temp)  // assign bias to tracers. 
        {
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog cat(params);
          cat.read_catalog(params._Input_dir_cat()+params._file_catalogue(), 0);
          vector<double>dm_field;
          vector<real_prec>dm_field_f;
          if(true==params._Get_tracer_bias())
          {
          dm_field.resize(params._NGRID(),0);
          File.read_array_t<double>(params._Input_Directory_X()+params._Name_Catalog_X(), dm_field);
          dm_field_f.resize(params._NGRID(),0);
          for (ULONG i = 0; i < dm_field.size(); i++)
            dm_field_f[i]=static_cast<real_prec>(dm_field[i]);       
          dm_field.clear();dm_field.shrink_to_fit();
          }
          PowerSpectrumF power(params);
          if(true==params._Get_tracer_bias())
          {
            power.object_by_object_bias(cat.Halo,dm_field_f);
            string out_bias=params._Output_directory()+"individual_bias"+"_"+params._Name_survey()+".txt";
            ofstream bout; bout.open(out_bias.c_str());
            cout<<"Writting bias in file "<<out_bias<<endl;
            for(ULONG i=0;i<cat.Halo.size();++i)
              bout<<cat.Halo[i].bias<<endl;
            bout.close();
          }
          if(true==params._Get_tracer_relative_bias())
          {
            // here we need to pass the number counts of the tracer field to get power . The class power has a member with the catalog power.object_by_object_rbias

            power.object_by_object_rbias();
            string out_bias=params._Output_directory()+"individual_bias"+"_"+params._Name_survey()+".txt";
            ofstream bout; bout.open(out_bias.c_str());
            cout<<"Writting relative bias in file "<<out_bias<<endl;
            for(ULONG i=0;i<cat.Halo.size();++i)
              bout<<cat.Halo[i].relative_bias<<endl;
            bout.close();
          }           
      }
      else if('m'==temp)  // to build a mock catalog based on a tracer sim (snapshot) using some dNdz
        {
          FileOutput File;
          So.welcome_message_c();
          par_file_BiasMT = argv[2];
          Params params(par_file_BiasMT);
          Catalog cat(params); // Read the N-body 
          cat.read_catalog(params._Input_dir_cat()+params._file_catalogue(), 0);
          cat.snap_to_mock();


      }

      else if('?'==temp)
       {
         cout<<endl;
         cout<<"Argument not recognized."<<endl;
         cout<<"Please run ./cosmicatlass.exe -h."<<endl;
         exit(1);
       }
    }
  exit(0);
  return 0;

}
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################
// ##################################################################################
// #################################################################################


