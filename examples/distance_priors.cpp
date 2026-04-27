// #################################################################################
/** @file 
 *
 *  @brief Example of mcmc using distance priors
 *  @authors Andrés Balaguera-Antolinez
 */

#include "../include/def.hpp"
#include "../include/CosmiCalcLIB.hpp"

using namespace std;

// ##################################################################################
// ##################################################################################

void run_mcmc(string par_file, int ch){
   auto start_time = chrono::high_resolution_clock::now();
   Params params(par_file);
   ScreenOutput So; 
   McmcFunctions mcmc(params);
   mcmc.set_seed_index(ch);
   mcmc.set_action("RUN");
   mcmc.set_mcmc_vectors();  

  mcmc.set_distance_priors();
  int j=0;
  size_t current_number_of_accepted_models=0;
  real_prec curr_H=2e15; // a High value to accept from here
  real_prec prop_H=0;

      //--------------------------------------------
     //NOTE: These numbers are used in the upgrade of the covariance matrix.
   size_t ns = params._number_of_steps_to_update_covariance();
   size_t nbi= params._number_of_burnin_phase_models();
   size_t nn = params._number_of_burnin_phase_models()+params._number_of_post_burnin_phase_models();

        //--------------------------------------------
        // NOTE: START OF MARKOFF CHAINS          
        do{
          if(j>0)
            mcmc.jump(current_number_of_accepted_models);

          j++;
          
          //tmodel.resize(params._number_of_time_bins());

          //--------------------------------------------
          // NOTE: Compute the chi squared 
          double sum_chis=0;
          sum_chis=mcmc.chi_squared_distance_priors();
          prop_H=sum_chis; //Update the chi squared

          //--------------------------------------------
          // NOTE: Metropolis-Hasting Algorithm: get acceptance probability MH with the prior information
          int weight_param=1; 
          mcmc.MHalgorithm(curr_H, prop_H, current_number_of_accepted_models, weight_param);
          //--------------------------------------------
          // NOTE: Save accepted parameters (or models)
          mcmc.write_accepted_models(j, current_number_of_accepted_models, weight_param, false); 
          //--------------------------------------------
          // NOTE: Compute the acceptance rate
#ifndef _OPTIMAL_PERFORMANCE_
          double acceptance_rate = 100.0*static_cast<double>(current_number_of_accepted_models)/static_cast<double>(j); 
           So.message_screen_flush("Current number of accepted models:", current_number_of_accepted_models, "Acceptance rate", acceptance_rate);    
#endif           
          //--------------------------------------------
          // NOTE: Update covariance after biurn in phase.
          if(params._update_covariance())
           {
             
             if ((current_number_of_accepted_models >= nbi && current_number_of_accepted_models< nn) && ((current_number_of_accepted_models - nbi) % ns ==0))
               mcmc.get_covariance_parameters();
           }          
        }while (current_number_of_accepted_models <= params._number_of_accepted_models());
        
        double total_acceptance_rate = 100.0*static_cast<double>(current_number_of_accepted_models)/static_cast<double>(j); 
        mcmc.lpout.close();
        So.message_screen("Chain", ch, "finished with acceptance rate [%]:", total_acceptance_rate);
#ifndef _RUN_PARALLEL_CHAINS_
        mcmc.get_stats();
        for(size_t ip=0;ip < params._number_of_fit_parameters();++ip)
          if(params.action_parameters[ip]!=1)
            mcmc.analyze_chains(ip);
        auto end_time = chrono::high_resolution_clock::now();
        double difft=chrono::duration<double>(end_time-start_time).count();
  //      mcmc.set_computing_time(difft);
//        mcmc.output_parameters();
#endif
#ifdef _TESTS_
      }        
#endif // closes if TEST




  }

// ##################################################################################
// ##################################################################################


int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  char temp;

  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);

  string logfile="logfile.log";
  ScreenOutput So(start_all, logfile);

  if(argc==1){
    So.usage(argv[0]);
    exit(1);
  }

  string par_file = argv[2];
  Params params(par_file);


  if(false==params.input_sections.MCMC)
    {
     throw std::runtime_error("Section MCMC is not enabled in parameter file"); 
     exit(0);
    }

  int option_index = 0;
  int c;
    static struct option long_options[] = {
        {"run", required_argument, 0, 'r'},
        {"analyze", required_argument, 0, 'a'},
        {0, 0, 0, 0}
    };

  Params params_aux(par_file);
  int n_chains=params_aux._number_of_chains();

 while ((c = getopt_long(argc, argv, "r:a:",long_options, &option_index)) != -1) 
   {
      switch (c) 
      {
        case 'r' :
          {
#pragma omp parallel num_threads(n_chains)
        {        
          int thread=omp_get_thread_num();
          run_mcmc(par_file, thread);
        }
          break;
       }
          
        case  'a':
      
        { 
          break;
        }
   }
  }
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


