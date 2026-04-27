////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
#include "../include/def.hpp"
#include "../include/CosmiCalcLIB.hpp"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{

  time_t start_all;
  time(&start_all);

#ifdef _USE_OMP_
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  ScreenOutput So(start_all,"logfile_ha.log");

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  string par_file=argv[2];
  Params params(par_file);

  if (getenv("RUNNING_IN_XTERM") == nullptr) {
      std::string cmd = "RUNNING_IN_XTERM=1 xterm -bg black -fg white -hold -e ";
    // Add program name
      cmd += argv[0];
      // Add all original arguments
          for (int i = 1; i < argc; ++i) {
              cmd += " ";
              cmd += argv[i];
        }
      system(cmd.c_str());
      return 0;
      }

  if(!params.input_sections.HaloAnalysis)
    {
      throw std::runtime_error("Section Halo Analysis is not enabled in parameter file"); 
      exit(0);
   }


  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"power", required_argument, 0, 'p'},
        {"window", required_argument, 0, 'w'},
        {"lpfilter",  required_argument, 0, 's'},
        {"grf",      no_argument,       0, 'g'},
        {"m",      required_argument,       0, 'm'},
        {0, 0, 0, 0}
    };

  while ((c = getopt_long(argc, argv, "c:h:b:u:s:m:",long_options, &option_index)) != -1)
  {
  
    switch (c) {
  
    case 'c': // to read input tracer catalg and analyze it.
    {

      Catalogue cat(params, "TRACER");
      HaloTools htools(params, cat); 
      htools.analyze_cat(true);
      So.message_time(start_all);
      break;
    }
    case 'h': // Read halo catalog and generate galaxy catalogs using an HOD
    {
      Catalogue cat(params, "TRACER");
      HaloTools htools(params, cat); 
      htools.halos2galaxies_HOD();
     break;
    }

    
    case 'm': // to build a mock catalog based on a tracer sim (snapshot) using some dNdz and angular selection function.
    {
      FileOutput File;
      So.welcome_message_c();
      Catalogue cat(params, "TRACER"); // Read the N-body
      HaloTools htools(params,cat);
      cat.read_catalog_new(params._Input_dir_cat() + params._file_catalogue());
      htools.snap_to_mock(false);
    break;
    }
    default:
    break;
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
