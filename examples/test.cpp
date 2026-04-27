////////////////////////////////////////////////////////////////////////////
/**
 *  @file Main file for BMT.
 *  @brief Uses the ComiCalc library and the file BiasMappingTechnique
 *  @authors Andrés Balaguera-Antolinez
 */
////////////////////////////////////////////////////////////////////////////

#include "../include/CosmiCalcLIB.h"

////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  ScreenOutput So(start_all, "logfile_bmt.log");

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  FileOutput File; 
  string par_file = argv[1];
  Params params(par_file);
 
  if(false==params.input_sections.DarkMatterCatalogue)
    {
     throw std::runtime_error("Section DarkMatterCatalogue is not enabled in parameter file"); 
     exit(0);
  }

  if(false==params.input_sections.TracerCatalogue)
    {
     throw std::runtime_error("Section TracerCatalogue is not enabled in parameter file"); 
     exit(0);
    }
    
  Catalogue tracer(params, "TRACER");// This method reads catalog

  fs::path filetr = fs::path(params._Input_dir_cat()) / params._file_catalogue(); 
  tracer.read_catalog_new(filetr.string()); // Read the asciii file still

  HaloTools tools(params);  //This class still reads, but we doin't call that method
  tools.perform_extra_taks_on_tracers(tracer);

  exit(0);
  return 0;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
