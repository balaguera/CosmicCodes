// #################################################################################
/** @file 
 *
 *  @brief Example of mc
 *  @authors Andrés Balaguera-Antolinez
 */

#include "../include/def.hpp"
#include "../include/CosmiCalcLIB.hpp"

using namespace std;

// #################################################################################
// #################################################################################
// #################################################################################
// #################################################################################
// #################################################################################
// #################################################################################
/**
  @brief Example of CMB distance priors for cosmological parameter constraints
*/
void measured_distance_priors(Params & params, vector<real_prec> &priors_measured, vector<vector<real_prec>>&inv_cov){

  int I= params._dprior_ID();
  int J= params._dprior_JD();  
  size_t n_priors=0;

  if(I==1 && (J==1 || J==2 || J==3 || J==4))n_priors=1;
  if(I==1 && (J==5 || J==6 || J==7))n_priors=2;
  if(I==1 && J==8) n_priors=3;
  if(I==2) n_priors=4;
  if(I==3) n_priors=3;
  if(I==4) n_priors=5;
  priors_measured.resize(n_priors,0);

  inv_cov.resize(n_priors);
  for(size_t i=0;i<n_priors;++i) 
    inv_cov[i].resize(n_priors,0);

  vector<vector<real_prec> > cov(n_priors,vector<real_prec>(n_priors,0));
  
  switch(I){
  case 1:
    switch(J){
    case 1:      /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5,	1) l_A*/
      priors_measured[0]=302.10;
      inv_cov[0][0]=pow(1./0.86,2);
      break;
      
    case 2: /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5	1) R = shift parameter  */ 
      priors_measured[0]=1.710;
      inv_cov[0][0]=pow(1./0.019,2);
      break;
      
    case 3:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5    1) z_*/
      priors_measured[0]=1090.04;
      inv_cov[0][0]=pow(1./0.93,2);
      break;
      
    case 4: /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5*   
	      1) r_s = sond horizon*/
      priors_measured[0]=146.8;
      inv_cov[0][0]=pow(1./1.8,2);
      break;
      
    case 5:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_A = shift parameter
	     2) R  = shift parameter
	   */
      priors_measured[0]=302.10;
      priors_measured[1]=1.710;
      for(int i=0;i<n_priors;++i)
        for(int j=0;j<n_priors;++j)
          cov[i][j]=0;
      cov[0][0]=pow(0.86,2);
      cov[1][1]=pow(0.019,2);
      cov[0][1]=0.1109*(0.86*0.019);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,inv_cov);
      break;
      
    case 6:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_A = sound horizon angular scale
	     2) z  = redshift at decoupling
	   */
      priors_measured[0]=302.10;
      priors_measured[1]=1090.04;
      for(int i=0;i<n_priors;++i)
        for(int j=0;j<n_priors;++j)
          cov[i][j]=0;
      cov[0][0]=pow(0.86,2);
      cov[1][1]=pow(0.93,2);
      cov[0][1]=0.4215*(0.86*0.93);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,inv_cov);
      break; 

    case 7:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) R_ = shift parameter
	     2) z  = redshift at decoupling
	   */

      priors_measured[0]=1.710;
      priors_measured[1]=1090.04;
      
      for(int i=0;i<n_priors;++i)
        for(int j=0;j<n_priors;++j)
          cov[i][j]=0;

      cov[0][0]=pow(0.019,2);
      cov[1][1]=pow(0.93,2);
      cov[0][1]=0.6928*(0.019*0.93);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,inv_cov);
      break;

    case 8:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_ = sound horizon angular scale
	     2) R_ = shift parameter
	     3 z  = redshift at decoupling
	   */
      priors_measured[0]= 302.10;
      priors_measured[1]= 1.710;
      priors_measured[2]= 1090.04;
      inv_cov[0][0]= 1.800;
      inv_cov[0][1]= 27.968;
      inv_cov[0][2]=-1.103;
      inv_cov[1][0]= inv_cov[0][1];
      inv_cov[1][1]= 5667.577;
      inv_cov[1][2]=-92.263;
      inv_cov[2][0]= inv_cov[0][2];
      inv_cov[2][1]= inv_cov[1][2];
      inv_cov[2][2]= 2.923;
      // //***********************************************************
        //invirtiendo la covariancia: esto es solo para chequear
      //for(int i=0;i<=n_priors;++i)for(int j=0;j<=n_priors;++j)cov[i][j]=0;
        // cov[1][1]= pow(0.86,2);
        // cov[1][2]= 0.1109*0.86*0.019;
        // cov[1][3]= 0.4215*0.86*0.93;
        // cov[2][1]= cov_priors_cmb[1][2];
        // cov[2][2]= pow(0.019,2);
        // cov[2][3]= 0.6928*1.019*0.93;
        // cov[3][1]= cov_priors_cmb[1][3];
        // cov[3][2]= cov_priors_cmb[2][3];
        // cov[3][3]= pow(0.93,2);
	
      //   //***********************************************************      
      break;
    }

    break;

  case 2:
    /*EXTENDED DISTANCE PRIORS FROM CMB+SNia+Ho FROM KOMATSU et al 2010
      1) l_A(z_star),  
      2) R(z_star), 
      3) z_star = redshift  at decoupling 
      4) 100w_b  
    */
    priors_measured[0]=302.10;
    priors_measured[1]=1.710;
    priors_measured[2]=1090.04;
    priors_measured[3]=2.2765;
    
    inv_cov[0][0]= 31.001;
    inv_cov[0][1]=-5015.642;
    inv_cov[0][2]= 183.903;
    inv_cov[0][3]= 2337.977;
    
    inv_cov[1][0]=inv_cov[0][1];
    inv_cov[1][1]= 876807.166;
    inv_cov[1][2]=-32046.750;
    inv_cov[1][3]=-403818.837;
    
    
    inv_cov[2][0]=inv_cov[0][2];
    inv_cov[2][1]=inv_cov[1][2];
    inv_cov[2][2]=1175.054;
    inv_cov[2][3]=14812.579;
    
    inv_cov[3][0]=inv_cov[0][3];
    inv_cov[3][1]=inv_cov[1][3];
    inv_cov[3][2]=inv_cov[2][3];
    inv_cov[3][3]=187191.186;
    
    break;    
    
  case 3:
    /*DISTANCE PRIORS FROM KOMATSU et al 2010 WMAP7
      1) l_A(z_star),  
      2) R(z_star), 
      3) z_star = redshift  at decoupling */
    priors_measured[0] = 302.09;
    priors_measured[1] = 1.725;
    priors_measured[2] = 1091.3;
    inv_cov[0][0]= 2.305;
    inv_cov[0][1]= 29.698;
    inv_cov[0][2]= -1.333;
    inv_cov[1][0]= inv_cov[0][1];
    inv_cov[1][1]= 6825.270;
    inv_cov[1][2]=-113.180;
    inv_cov[2][0]= inv_cov[0][2];
    inv_cov[2][1]= inv_cov[1][2];
    inv_cov[2][2]= 3.414;
    
    break;
    
    
  case 4:
    /*DISTANCE PRIORS FROM SDSS(LSS)+CMB+SNia FROM SANCHEZ et al. 2009
      1) 100 \Omega_bar h^2, 
      2) z_star = redshift  of decoupling 
      3) l_A(z_star),  
      4) R(z_star), 
      5) G(z_m), z_m=0.35*/
    priors_measured[0]=2.28;
    priors_measured[1]=1090.12;
    priors_measured[2]=301.58;
    priors_measured[3]=1.701;
    priors_measured[4]=1175.;
    
    for(int i=1;i<n_priors;++i)
      for(int j=1;j<n_priors;++j)
        cov[i][j]=0;
    
    cov[0][0]= 2.99e-7;
    cov[0][1]=-4.2030e-4;
    cov[0][2]=-1.9988e-4;
    cov[0][3]=-3.4393e-6;
    cov[0][4]=-1.9978e-3;
    cov[1][0]=cov[0][1]; 
    cov[1][1]=8.6812e-1;
    cov[1][2]=0.2557;
    cov[1][3]=1.1999e-1;
    cov[1][4]=3.2355;
    cov[2][0]=cov[0][2];
    cov[2][1]=cov[1][2];
    cov[2][2]=0.4558;
    cov[2][3]=3.1265e-2;
    cov[2][4]=1.8247;
    cov[3][0]=cov[0][3];
    cov[3][1]=cov[1][3];
    cov[3][2]=cov[2][3];
    cov[3][3]=3.1460e-4;
    cov[3][4]=3.9649e-2;
    cov[4][0]=cov[0][4];
    cov[4][1]=cov[1][4];
    cov[4][2]=cov[2][4];
    cov[4][3]=cov[3][4];
    cov[4][4]=4.3784e2;
    matrix_inversion(cov,inv_cov);
    break;

   default:
    break;

  }
  
  return ;
}
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################
/**
  @brief Example of ussage of CMB distance priors for cosmological parameter constraints
*/

void distance_priors_cmb_model(Params &params, s_CosmologicalParameters &scp, vector<real_prec>&model_priors){

  int I= params._dprior_ID();
  int J= params._dprior_JD();  

  Cosmology CF(scp);
  /*Distance priors from Komatsu et al. 2010 (CMB)*/
  /*Calculamos z_dec segun Hu & Sujiyama, ver Komatsu et al. 2007 WMAP5 eq 65 */
  real_prec omega_b  = (scp.Om_baryons)*pow(scp.hubble,2);
  real_prec omega_m  = (scp.Om_matter)*pow(scp.hubble,2);

  // Baryon drag epoch:
  real_prec z_drag = CF.drag_redshift(); 
  real_prec da_zd  = CF.proper_angular_diameter_distance(z_drag);
  real_prec rsound = CF.comoving_sound_horizon(z_drag);
  real_prec R      = (1+z_drag)*da_zd*sqrt(scp.Om_matter)*scp.Hubble/Constants::speed_light;
  real_prec lA     = (1+z_drag)*da_zd*M_PI/rsound;

  if(I==1 && J==1)model_priors[0]=lA;
  if(I==1 && J==2)model_priors[0]=R;
  if(I==1 && J==3)model_priors[0]=z_drag;
  if(I==1 && J==4)model_priors[0]=rsound;
  if(I==1 && J==5){
    model_priors[0] =  lA;
    model_priors[1] =  R;
  }
  if(I==1 && J==6){
    model_priors[0] =  lA;
    model_priors[1] =  z_drag;
  }
  if(I==1 && J==7){
    model_priors[0] =  R;
    model_priors[1] =  z_drag;
  }
  if((I==1 && J==8) || I==3){
    model_priors[0] =  lA;
    model_priors[1] =  R;
    model_priors[2] =  z_drag;
  }
  if(I==2 && J==1){
    model_priors[0] =  lA;
    model_priors[1] =  R;
    model_priors[2] =  z_drag;
    model_priors[3] =  100.0*omega_b;
  }
  if(I==4){
    model_priors[0] =  100.0*omega_b;
    model_priors[1] =  z_drag;
    model_priors[2] =  lA;
    model_priors[3] =  R;
    model_priors[4] =  (1+Constants::z_mean_sdss)*CF.proper_angular_diameter_distance(Constants::z_mean_sdss)*pow(0.01*CF.Hubble_function(Constants::z_mean_sdss),0.8);
  }

  
}

// ##################################################################################
// ##################################################################################
// ##################################################################################
// ##################################################################################

void get_results(string par_file)
{
  ScreenOutput So;
  Params params(par_file);
  McmcFunctions mcmc(params);
  mcmc.set_action("ANALYZE");
  mcmc.set_mcmc_vectors();
  mcmc.read_nchains();
  cout<<endl;
  So.message_screen("Getting staitistics:");
  for(size_t ip=0; ip < params._number_of_fit_parameters(); ++ip) // Loop over the chains
    if(params.action_parameters[ip]!=1)
      mcmc.analyze_chains(ip);
   mcmc.output_parameters();
}


// ##################################################################################
// ##################################################################################

void run_mcmc(string par_file, int ch)
{
   auto start_time = chrono::high_resolution_clock::now();

   Params params(par_file);
   ScreenOutput So; 

   McmcFunctions mcmc(params);

   mcmc.set_seed_index(ch);
   mcmc.set_action("RUN");
   mcmc.set_mcmc_vectors();  
 
   vector<real_prec>measured_priors;
   vector<vector<real_prec>>inv_cova;

   measured_distance_priors(params, measured_priors, inv_cova);
   
   vector<real_prec>model_priors(measured_priors.size(),0);

   s_CosmologicalParameters scp=params.s_cosmo_pars;
 
   int j=0;

   size_t current_number_of_accepted_models=0;

   real_prec curr_H=2e15; // a High value to accept from here
   real_prec prop_H=0;

      //--------------------------------------------
     //NOTE: These numbers are used in the upgrade of the covariance matrix.
    size_t ns = params._number_of_steps_to_update_covariance();
    size_t nbi= params._number_of_burnin_phase_models();
    size_t nn = params._number_of_burnin_phase_models()+params._number_of_post_burnin_phase_models();

    int nchains = params._number_of_chains();
    bool write_sc=true;
    if(nchains>1)
      write_sc=false;
    //--------------------------------------------
        // NOTE: START OF MARKOFF CHAINS          
    do{
       if(j>0)
         mcmc.jump(current_number_of_accepted_models);

        j++;

        //--------------------------------------------
        scp.Om_matter  = mcmc.parameters[0];
        scp.Om_baryons = mcmc.parameters[1];
        scp.wde_eos    = mcmc.parameters[2];
        distance_priors_cmb_model(params, scp, model_priors);
//        for(int i=0;i<model_priors.size();++i)cout<<i<<"  "<<measured_priors[i]<<"  "<<model_priors[i]<<endl;
        //--------------------------------------------
        // NOTE: Compute the chi squared 
        real_prec sum_chis=0;
        mcmc.chi_squared(measured_priors, model_priors,inv_cova,sum_chis);
        prop_H=sum_chis; //Update the chi squared
        //--------------------------------------------
        // NOTE: Metropolis-Hasting Algorithm: get acceptance probability MH with the prior information
        int weight_param=1; 
        mcmc.MHalgorithm(curr_H, prop_H, current_number_of_accepted_models, weight_param);
          //--------------------------------------------
          // NOTE: Save accepted parameters (or models)
        mcmc.write_accepted_models(j, current_number_of_accepted_models, weight_param, write_sc); 
          //--------------------------------------------
          // NOTE: Compute the acceptance rate
        if(nchains<2 || !write_sc)
        {
          double acceptance_rate = 100.0*static_cast<double>(current_number_of_accepted_models)/static_cast<double>(j); 
//          So.message_screen_flush("Current number of accepted models:", current_number_of_accepted_models, "Acceptance rate", acceptance_rate);    
        }
        //--------------------------------------------
          // NOTE: Update covariance after biurn in phase.
        if(params._update_covariance())
          {
           if ((current_number_of_accepted_models >= nbi && current_number_of_accepted_models< nn) && ((current_number_of_accepted_models - nbi) % ns ==0))
             mcmc.get_covariance_parameters();
           }          
        //--------------------------------------------

        }while (current_number_of_accepted_models < params._number_of_accepted_models());
        

        double total_acceptance_rate = 100.0*static_cast<double>(current_number_of_accepted_models)/static_cast<double>(j); 
        mcmc.lpout.close();
        So.message_screen("Chain", ch, "finished with acceptance rate [%]:", total_acceptance_rate);

        if(nchains<2)
        {
          mcmc.get_stats();
          for(size_t ip=0;ip < params._number_of_fit_parameters();++ip)
            if(params.action_parameters[ip]!=1)
              mcmc.analyze_chains(ip);
          auto end_time = chrono::high_resolution_clock::now();
          double difft=chrono::duration<double>(end_time-start_time).count();
        }
        //        mcmc.output_parameters();


  }

// ##################################################################################
// ##################################################################################


int main(int argc, char *argv[]){

  time_t start_all;
  time(&start_all);
  char temp;

  int NTHREADS=omp_get_max_threads();
  omp_set_num_threads(NTHREADS);

  ScreenOutput So;
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
  size_t n_pars = params.action_parameters.size();
 while ((c = getopt_long(argc, argv, "r:a:",long_options, &option_index)) != -1) 
   {
      switch (c) 
      {
        case 'r' :
          {
            So.message_screen("MCMC for cosmological_distance_priors");
            So.message_screen("Number of parallel chains: ", n_chains);   
            So.message_screen("Number of fitting parameters: ", n_pars);   

#pragma omp parallel num_threads(n_chains)
            {        
              int thread=omp_get_thread_num();
              run_mcmc(par_file, thread);
            }
            So.message_screen("Analyzing chains:");
            get_results(par_file);

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


