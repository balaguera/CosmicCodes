////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
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
  while((temp =  getopt(argc, argv, "hadi:c:h:b:")) != -1)
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


