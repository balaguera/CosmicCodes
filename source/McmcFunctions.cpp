////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<McmcFunctions>
 * @file McmcFunctions.cpp
 * @brief Methods of the class McmcFunctions
 * @details McmcFunctions Bayesian methods
 * @author Andres Balaguera Antolinez 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../headers/McmcFunctions.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec maf(string mas, real_prec x){
  /*Mas assigment*/
  real_prec ans;
  real_prec daux;
  if(mas=="NGP")
   {
    daux=fabs(x);
    if(daux<0.5)
      ans=1.0;

    if(daux==0.5)
      ans=0.5;

    if(daux>0.5)
     ans=0.0;
  }
  if(mas=="CIC")

    ans=(fabs(x)<1 ? 1.-fabs(x): 0);

  if(mas=="TSC")
  {
    if(fabs(x)<0.5)ans=(0.75-x*x);
    if(fabs(x)>=0.5 && fabs(x)<1.5) ans= (0.5*pow(1.5-fabs(x),2));
    if(fabs(x)>=1.5)ans=0.0;
  }
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void McmcFunctions::read_parameters(string parameters_file, string code){
  ifstream fin_parameters (parameters_file.c_str());
  if (!fin_parameters)
   cerr <<"Error in opening the parameters file "<<parameters_file<<"!"<<endl; exit(1);
  gsl_rng_env_setup(); 
  gsl_rng_default_seed=this->seed; 
  this->T = gsl_rng_ranlux; 
  this->r = gsl_rng_alloc (T); 
  string line_in_file;
  string par_name;
  string equality;
  string par_value;
  string par_value0;
  string par_value1;
  string par_value2;
  string par_value3;
  string par_value4;
  string par_value5;
  string par_value6;
  while (getline(fin_parameters,line_in_file)) {
    // The single values will start with a P in the parameter file
    if(line_in_file[0] != '#' && line_in_file.empty()==0 && line_in_file[0] =='P'){
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;
      
      // check that second word is "="
      line_string >> equality;
      if(equality != "=") {
        cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
        cerr << "Using a default value for " << par_name << endl; exit(1);
      }
      line_string >> par_value;
      if(par_value.empty()) {
        cout << "Value of " << par_name << " not specified in " << parameters_file << endl;
        continue;
      }
      else if (par_name == "P_use_distance_priors"){
        if(par_value=="yes")use_distance_priors=true;
        else use_distance_priors=false;
      }
      else if (par_name == "P_diagonal_covariance_matrix"){
        if(par_value=="yes")this->diagonal_covariance_matrix=true;
        else this->diagonal_covariance_matrix=false;
      }
      else if (par_name == "P_analytic_marginalization_wrt_amplitude"){
        if(par_value=="yes")this->analytic_marginalization_wrt_amplitude =true;
        else this->analytic_marginalization_wrt_amplitude =false;
      }
      if(par_name=="P_experiment")this->experiment=par_value;
      if(par_name=="P_output_dir")this->output_dir=par_value;

      if(par_name=="P_model")this->model=par_value;
      if(par_name=="P_observable")this->observable=par_value;
      if(par_name=="P_n_parameters")this->n_parameters=atoi(par_value.c_str());
      if(par_name=="P_sampling")this->sampling=par_value;
      if(par_name=="P_nbins_1D")this->nbin_1D=atoi(par_value.c_str());
      if(par_name=="P_nbins_2D")this->nbin_2D=atoi(par_value.c_str());
      if(par_name=="P_MAS")this->MAS=par_value;
      if(par_name=="P_nchains")this->nchains=atoi(par_value.c_str());
      if(par_name=="P_chain_ini")this->chain_ini=atoi(par_value.c_str());
      if(par_name=="P_chain_fin")this->chain_fin=atoi(par_value.c_str());
      if(par_name=="P_npp")this->n_steps=atoi(par_value.c_str());
      if(par_name=="P_Id")this->Id=atoi(par_value.c_str());
      if(par_name=="P_Jd")this->Jd=atoi(par_value.c_str());
      if(par_name=="P_msteps")this->msteps=atoi(par_value.c_str());
      if(par_name=="P_input_par_file_model")this->input_par_file_model=par_value;
      if(par_name=="P_fit_type")this->fit_type=par_value;
    }
        // ignore lines starting with hashtag and read the vectors
     // Here read parameter that do not start with the character P_
    else if (line_in_file[0] != '#' && line_in_file.empty()==0 && line_in_file[0] != 'P') {
      stringstream line_string (line_in_file);
      line_string << line_in_file;
      line_string >> par_name;
      // check that second word is "="
      line_string >> equality;
      if (equality != "=") {
      cerr << "Error in parameters file at " << par_name << " = " << "???" << endl;
      cerr << "Using a default value for " << par_name << endl; exit(1);
    }
     // read first parameter value
    line_string >> par_value0;
    if(par_value0.empty()) {
      cout << "First entry of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }

      // read second parameter value
    line_string >> par_value1;
    if (par_value1.empty()) {
      cout << "Second entry of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }

    line_string >> par_value2;
    if (par_value2.empty()) {
      cout << "Third value of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }


    line_string >> par_value3;
    if (par_value3.empty()) {
      cout << "Third value of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }

      // read 4 parameter value 
    line_string >> par_value4;
    if(par_value4.empty()) {
      cout << "Fourth value of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }

      
      // read 5 parameter value 
    line_string >> par_value5;
    if (par_value5.empty()) {
      cout << "Fifth value of " << par_name << " not specified in " << parameters_file << endl;
      continue;
    }
    // read 6 parameter value
  line_string >> par_value6;
  if (par_value6.empty()) {
    cout << "Fifth value of " << par_name << " not specified in " << parameters_file << endl;
    continue;
  }
    // Now get the numbers 
    for(int i=0;i<this->n_parameters;i++){
      string parname = "par"+to_string(i);
      if(par_name==parname){
        this->fixed_parameters.push_back(atoi(par_value0.c_str()));
        this->parameters_max.push_back(atof(par_value1.c_str()));
        this->parameters_min.push_back(atof(par_value2.c_str()));
        this->delta_parameters.push_back(atof(par_value3.c_str()));
        this->error_fixed_parameters.push_back(atof(par_value4.c_str()));
        this->parameters_max_2dplot.push_back(atof(par_value5.c_str()));
        this->parameters_min_2dplot.push_back(atof(par_value6.c_str()));
      }
    }
    // If no min and max specified for plots, use the initial min and max values
    if (par_value5.empty()){
        for(int i=0;i<this->n_parameters;++i)this->parameters_max_2dplot.push_back(this->parameters_max[i]);
        continue;
    }
    if (par_value6.empty()){
        for(int i=0;i<this->n_parameters;++i)this->parameters_min_2dplot.push_back(this->parameters_min[i]);
        continue;
    }
    }
    
  }
  parameters.resize(this->n_parameters,0);
  parameters_ini.resize(this->n_parameters,0);
  for(int i=0;i<this->n_parameters;i++)
    if(this->fixed_parameters[i]!=0)
        this->parameters_ini[i]=this->parameters_min[i]+(this->parameters_max[i]-this->parameters_min[i])*gsl_rng_uniform(this->r);
  // Append name of output chains
  string dir_out = this->output_dir;
 // if(code!="run" || code!="analyze")cout<<RED<<"PLEASE PROVIDE A VALID OPTION, run or analize"<<RESET<<endl;
  if(code=="run"){
  cout<<RED<<"WARNING: creation of output file commented in McmcFunctions.cpp"<<RESET<<endl;
  if(this->analytic_marginalization_wrt_amplitude)this->file_oo = "acc_par_"+this->experiment+"_"+this->observable+"_"+this->model+"_"+this->sampling+"_amp_marg_chain"+to_string(this->chain)+".txt";
  else this->file_oo = "acc_par_"+this->experiment+"_"+this->observable+"_"+this->model+"_"+this->sampling+"_chain"+to_string(this->chain)+".txt";
  this->output_file = dir_out+this->file_oo;
  this->lpout.open(output_file.c_str());
  this->lpout.precision(12);
  this->lpout.setf(ios::showpoint);
  this->lpout.setf(ios::scientific);  
  }
  this->nchains=this->chain_fin-this->chain_ini+1;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

void McmcFunctions::set_mcmc_read_vectors(){

  cout<<BLUE<<"Defining space for accepted models"<<RESET<<endl;
  this->acc_parameters_all.resize(n_parameters);

  for (int i = 0; i<n_parameters;i++) {
    this->acc_parameters_all[i].resize(n_steps);
    for (int j=0;j<n_steps;j++){
      this->acc_parameters_all[i][j].resize(nchains,0);
    }
  }

  this->weight_r.resize(n_steps);
  for(int i=0;i<weight_r.size();i++)weight_r[i].resize(nchains,0);

  this->chi_models.resize(n_steps);
  for(int i=0;i<chi_models.size();i++)chi_models[i].resize(nchains,0);

  this->mark.resize(nchains,0);
  this->acca.resize(nchains,0);
  this->R.resize(n_parameters);

}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
void McmcFunctions::get_cova_pars(int acc){
  // Compute the covariance of the parameters
  // at the step acc of the markov chain
  int nacc=(acc==1? acc:acc-1);
  for(int ip=0;ip<n_parameters;ip++)this->mean_parameters[ip]=0;
  for(int ip=0;ip<n_parameters;ip++)if(this->fixed_parameters[ip]!=0)for(int jp=0;jp<n_parameters;++jp)this->cov_masses[ip][jp]=0;
  for(int ip=0;ip<n_parameters;ip++)if(this->fixed_parameters[ip]!=0)for(int ai=0;ai<acc;++ai)this->mean_parameters[ip]+=this->acc_parameters[ip][ai]/acc;
  for(int ip=0;ip<n_parameters;ip++)for(int jp=0;jp<n_parameters;++jp)for(int ai=0;ai<acc;++ai)this->cov_masses[ip][jp]+=(this->acc_parameters[ip][ai]-this->mean_parameters[ip])*(this->acc_parameters[jp][ai]-this->mean_parameters[jp])/(acc-1);
  
}


// ************************************************************************************
// ************************************************************************************
//#define _USE_HMS_
// ************************************************************************************
void McmcFunctions::set_mcmc_vectors(){

  acc_parameters.resize(this->n_parameters);
  for(int i=0;i<acc_parameters.size();i++)
      acc_parameters[i].resize(this->n_steps,0);

  weight.resize(this->n_steps,0);
  chiss.resize(this->n_steps,0);

  mean_parameters.resize(this->n_parameters,0);
  stdev_parameters.resize(this->n_parameters,0);
  d_min.resize(this->n_parameters,0);
  d_max.resize(this->n_parameters,0);


#ifdef _USE_HMS_
  grad_U.resize(this->n_parameters,1);  
  momentum.resize(this->n_parameters,1);  
  masses.resize(this->n_parameters,0);  
  epsilon.resize(this->n_parameters,0);  

  cov_masses.resize(this->n_parameters);
  icov_masses.resize(this->n_parameters);
  for(int i=0;i<this->n_parameters;++i){
    cov_masses[i].resize(this->n_parameters,0);
    icov_masses[i].resize(this->n_parameters,0);
  }
#endif
}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
//Compute the kinetic energy in the Hamiltonian sampling
// Input 0 if the mass matrix is diagonal
// 1 if not
real_prec McmcFunctions::kinetic(int in){

  real_prec ans=0;
  if(in==0){ // For diagonal mass matrix
    for(int ip=0;ip<momentum.size();ip++)if(this->fixed_parameters[ip]!=0)ans+=0.5*pow(this->momentum[ip],2)/this->masses[ip];
  }
  else if(in==1){ //For full mass matrix
    for(int ip=0;ip<this->momentum.size();++ip)if(this->fixed_parameters[ip]!=0)for(int jp=0;jp<this->momentum.size();++jp)if(this->fixed_parameters[jp]!=0)ans+=0.5*this->momentum[ip]*this->cov_masses[ip][jp]*this->momentum[jp];
  }
  return ans;    
}
// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
// Trastearlo a FB
void McmcFunctions::get_model_FB(vector<vector<real_prec> >&Step,vector<vector<real_prec> >&R,vector<vector<real_prec> >&V,vector<real_prec>&Cmodel){

  int nmodes=R[0].size();
  int nmodes_n=Cmodel.size();

  for(int ii=0;ii<nmodes_n;ii++)Cmodel[ii]=0;

  for(int i=0;i<nmodes_n;i++){
    for(int k=1;k<Step.size();++k){
      for(int j=0;j<nmodes;j++){
	Cmodel[i]+=pow(R[i][j] + this->parameters[0]*V[i][j],2)*pow(10,this->parameters[k])*Step[k][j];
      }
    }
  }

}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************

void McmcFunctions::get_cova_FB_delta(vector<matrices>&VM, int l, vector<vector<real_prec> >&iCov){

  int nmodes=VM[l].R[0].size();
  int nmodes_n=VM[l].R.size();
  vector<vector<real_prec> > cova;
  cova.resize(nmodes_n);
  for(int i=0;i<cova.size();i++)cova[i].resize(nmodes_n,0);

  for(int i=0;i<nmodes_n;i++){
    for(int j=0;j<nmodes_n;j++){
      for(int k=1;k<VM[l].Step.size();++k){
        for(int p=0;p<nmodes;p++){
          cova[i][j]+=(VM[l].R[i][p] + this->parameters[0]*VM[l].V[i][p])*(VM[l].R[j][p] + this->parameters[0]*VM[l].V[j][p])*pow(10,this->parameters[k])*VM[l].Step[k][p] + VM[l].N[i][j];
        }
      }
    }
  }
  matrix_inversion(cova,iCov);
  
  
}
// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
// ************************************************************************************
void McmcFunctions::get_gradU(vector<real_prec>&Cmodel,vector<real_prec>&Cmeas,  vector<vector<real_prec> >&Step, vector<vector<real_prec> >&R,vector<vector<real_prec>>&V,vector<vector<real_prec> >&icov){
  // Get the gradient of the potential and also the model for this step.
  // The model is rerturned and then will be assigned as member of the FB class
  int nn=Cmodel.size();
  int nn_n=R[0].size();
  dCmod_dp.resize(nn);
  for(int i=0;i<nn;++i)dCmod_dp[i].resize(this->n_parameters,0);
  
  // // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  get_model_FB(Step,R,V,Cmodel);

  
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Compute derivative of the model in the current step
  for(int ii=0;ii<nn;ii++)for(int ip=0;ip<this->n_parameters;ip++)dCmod_dp[ii][ip]=0;

  // Derivative wrt parameter[0]
  if(this->fixed_parameters[0]!=0){
    for(int ii=0;ii<nn;ii++)for(int ij=0;ij<nn_n;ij++)for(int k=1;k<Step.size();k++)dCmod_dp[ii][0]+=2.*(R[ii][ij]+this->parameters[0]*V[ii][ij])*V[ii][ij]*pow(10,this->parameters[k])*Step[k][ij];
  }
  else for(int ii=0;ii<nn;ii++)dCmod_dp[ii][0]=0;
  
  // Derivative wrt parameter[k] k>=1
  for(int ii=0;ii<nn;ii++)for(int k=1;k<Step.size();k++)for(int ij=0;ij<nn_n;ij++)dCmod_dp[ii][k]+=pow(R[ii][ij]+this->parameters[0]*V[ii][ij],2.)*log(10.0)*pow(10,this->parameters[k])*Step[k][ij];
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Compute the grad U, i.e, the derivative of the chi squared:
  for(int ip=0;ip<parameters.size();ip++)this->grad_U[ip]=0;
  for(int ip=0;ip<parameters.size();ip++)if(this->fixed_parameters[ip]!=0)for(int ii=0;ii<nn;ii++)for(int ij=0;ij<nn;ij++)this->grad_U[ip]-=dCmod_dp[ii][ip]*icov[ii][ij]*(Cmeas[ij]-Cmodel[ij]);

}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************

void McmcFunctions::get_gradU_full(int lmin, int lmax, vector<matrices> VM,vector<real_prec>&Cmodel,vector<real_prec>&gradU_full){
  // Get the gradient of the potential and also the model for this step.
  // The model is rerturned and then will be assigned as member of the FB class

  for(int ip=0;ip<parameters.size();ip++)this->grad_U[ip]=0;

  real_prec full_like;
  
  for(int l=lmin;l<=lmax;l++){
    
    // // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    int nn=VM[l].R.size();
    int nn_n=VM[l].R[0].size();

    vector<real_prec> model(nn,0);
    get_model_FB(VM[l].Step,VM[l].R,VM[l].V,model);
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // Compute derivative of the model in the current step
    dCmod_dp.resize(nn);
    for(int i=0;i<nn;++i)dCmod_dp[i].resize(this->n_parameters,0);

    
    // Derivative wrt parameter[0]
    if(this->fixed_parameters[0]!=0){
      for(int ii=0;ii<nn;ii++)for(int ij=0;ij<nn_n;ij++)for(int k=1;k<VM[l].Step.size();k++)dCmod_dp[ii][0]+=2.*(VM[l].R[ii][ij]+this->parameters[0]*VM[l].V[ii][ij])*VM[l].V[ii][ij]*pow(10,this->parameters[k])*VM[l].Step[k][ij];
    }
    else for(int ii=0;ii<nn;ii++)dCmod_dp[ii][0]=0;
    
    // Derivative wrt parameter[k] k>=1
    for(int ii=0;ii<nn;ii++)for(int k=1;k<VM[l].Step.size();k++)for(int ij=0;ij<nn_n;ij++)dCmod_dp[ii][k]+=pow(VM[l].R[ii][ij]+this->parameters[0]*VM[l].V[ii][ij],2.)*log(10.0)*pow(10,this->parameters[k])*VM[l].Step[k][ij];


    // Compute now the chi square for each l
    // real_prec chiss;
    // chi_squared(VM[l].Cmed,model,VM[l].iCov,chiss);
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //Compute the grad U, i.e, the derivative of the chi squared for each l
    // and multiply it by the likelihood
    // Lineas comnetadas 
    // para hacer el analisis sumando las posterior
    // en vez de sumar las likelihoods

    //    for(int ip=0;ip<parameters.size();ip++)if(this->fixed_parameters[ip]!=0)for(int ii=0;ii<nn;ii++)for(int ij=0;ij<nn_n;ij++)this->grad_U[ip]-=(exp(-0.5*chiss)/sqrt(2.*M_PI*VM[l].det_matrix))*dCmod_dp[ii][ip]*VM[l].iCov[ii][ij]*(VM[l].Cmed[ij]-model[ij]);
    //full_like+=exp(-0.5*chiss)/sqrt(2.*M_PI*VM[l].det_matrix);

    for(int ip=0;ip<parameters.size();ip++)if(this->fixed_parameters[ip]!=0)for(int ii=0;ii<nn;ii++)for(int ij=0;ij<nn_n;ij++)this->grad_U[ip]-=dCmod_dp[ii][ip]*VM[l].iCov[ii][ij]*(VM[l].Cmed[ij]-model[ij]);
    
    

  }

  // Divide ny the total likelihood
  //  for(int ip=0;ip<parameters.size();ip++)grad_U[ip]/=full_like;
  
  

}



// ************************************************************************************
// ************************************************************************************
// ************************************************************************************

void McmcFunctions::jump(){
  
  So.message_screen("Jumping to a new position in parameter space");
  for(int ip=0;ip<this->parameters.size();ip++){
    switch(this->fixed_parameters[ip]){
    case(0):
      this->parameters[ip]=this->parameters_ini[ip];
      break;
    case(1):
      this->parameters[ip]+=gsl_ran_gaussian(r, delta_parameters[ip]);
      break;  
    }
  }
  So.DONE();
}

// ************************************************************************************
// ************************************************************************************
// ************************************************************************************



// ************************************************************************************
// ************************************************************************************
// ************************************************************************************

void McmcFunctions::get_mean(vector<real_prec>  &weight,int slab,vector<vector<real_prec> > &AA){
  
  // ******************************************************************************************************
  // Computes the mean parameter in a Markov-Chain using the weights computed
  // from the MH algorithm. Returns the mean and the standard deviation
  // for all the np parameters
  // weight = matrix, weights of the MCMC for each model evaluated
  // msteps = 1 ofe all steps are to be used, 2 if we want to take every two stepes and so on
  // acc   = number of accepted models
  // slab  = label of the first model taken after burn-in phase
  // np    = number of parameters to be constrained
  // AA    = matrix (parameters , model), 
  // meanA = vector (parameters), returns the mean of each parameters 
  // stdev = vector (parameters), returns the standard deviation within the chain
  // NOTE: COMPUTE ALSO THE MEDIAN
  // ******************************************************************************************************
  
  real_prec new_acc=0;
  // WARNING , AL MOVER ESTO A C++ EMPIEZO A LLENAR TODO DESDE 0
  // no necesito pasar aca el numero np, uso .size()

  for(int ip=0;ip<this->mean_parameters.size();ip++){
    new_acc=0;
    for(int j=slab;j<weight.size();j=j+this->msteps)if(AA[ip][j]!= -999)new_acc+=weight[j];
    real_prec mm=0;
    for(int j=slab;j<weight.size();j=j+this->msteps)if(AA[ip][j]!= -999)mm+=AA[ip][j]*weight[j]/new_acc;
    this->mean_parameters[ip]=mm;
  }
  /*compute standard deviation*/
  int new_acc_warning=0;
  for(int j=slab;j<AA[0].size();j=j+this->msteps)if(AA[0][j]!= -999)new_acc++;

  if(new_acc_warning==1){
    So.message_screen("Warning: only one accepted model in total!");
    exit(1);
  }

  for(int ip=0;ip<this->stdev_parameters.size();ip++){
    real_prec  stdeviation=0;
    for(int j=slab;j<AA[0].size();j=j+this->msteps){
     if(AA[ip][j]!= -999) stdeviation+=weight[j]*pow(AA[ip][j]-this->mean_parameters[ip],2);
    }
    this->stdev_parameters[ip]=sqrt(stdeviation/(new_acc-1.));
  }
  return;
}

/************************************************************************************************************************************************************************/
/************************************************************************************************************************************************************************/
/************************************************************************************************************************************************************************/
/************************************************************************************************************************************************************************/

void McmcFunctions::gelman_rubbin_diag(int n,vector<int> &mark, vector<int> &acca ,vector<vector<vector<real_prec> > > &acc_parameters){
  
  // ******************************************************************************************************
  // Gelman Rubbin diagnostics to check convergence of differnet Markov chains :                      
  // n           = number of free parameters
  // nchains_tot = number of chains to be analized
  // chain_ini   = label of the first chain 
  // chain_fin   = label of the last chain 
  // mark[nc]    = label from which the chain nc is going to be analized
  // acca[nc]    = number of points in the chain nc after the first passage 
  // through the point with chis=chir
  // acca_parameters [n][acca[nc]][nc]  = points in the chain 
  // R = vector(n), if R<1.03, chains have converged
  // ******************************************************************************************************
  
  So.message_screen("Checking for convergence (Gelman & Rubbin diagnostic): R<1.03");
  vector<real_prec> B(n,0);
  vector<real_prec> W(n,0);
  vector<real_prec> Wm(n,0);
  vector< vector<real_prec> > mean_chain;
  vector<real_prec> mean_all(n,0);
  mean_chain.resize(n);
  for(int i=0;i<n;++i)mean_chain[i].resize(this->nchains,0);
  cout<<this->nchains<<endl;
  /*Mean within each chain:*/


  for(int i=0;i<n;++i){
    for(int nc=0;nc<this->nchains;nc++){
     for(int j=mark[nc];j<acca[nc];++j)mean_chain[i][nc]+=acc_parameters[i][j-mark[nc]+1][nc]/(acca[nc]-mark[nc]);
    }
  }


  for(int i=0;i<R.size();++i)R[i]=0;
  for(int i=0;i<mean_all.size();++i)for(int nc=0;nc<this->nchains;nc++)mean_all[i]+=mean_chain[i][nc]/this->nchains;                               /*Mean of all chains:*/
  for(int i=0;i<B.size();++i)for(int nc=0;nc<this->nchains;nc++)B[i]+=pow(mean_all[i]-mean_chain[i][nc],2)/(this->nchains-1);               /*Variance between chains:*/

  /*Variance within each chain*/
  for(int i=0;i<W.size();++i){
    for(int nc=0;nc<this->nchains;nc++){
      for(int j=mark[nc];j<acca[nc];++j){
        W[i]+=pow(mean_chain[i][nc]-acc_parameters[i][j-mark[nc]+1][nc],2)/(this->nchains*(acca[nc]-mark[nc]-1));
        Wm[i]+=pow(mean_chain[i][nc]-acc_parameters[i][j-mark[nc]+1][nc],2)/(this->nchains*(acca[nc]-mark[nc]));
      }
    }
  }
  /*Potential factor*/
  for(int i=0;i<this->R.size();++i)this->R[i]=(Wm[i]+B[i])/W[i];

  Wm.clear();
  B.clear();
  W.clear();
  return ;
}

// *******************************************************************************************************
// *******************************************************************************************************

void McmcFunctions::chi_squared(vector<real_prec> &data, vector<real_prec> &model,vector< vector <real_prec> > &inv_cova, real_prec &chis_one){

  // Compute chi square given the data, the model and the inverse
  //  of the covariance matrix of dimension nXn
  chis_one=0;
  for(int i=0;i<model.size();++i){
    for(int k=0;k<model.size();++k){
      chis_one+=(data[i]-model[i])*inv_cova[i][k]*(data[k]-model[k]);
    }
  }
}

// *******************************************************************************************************
// *******************************************************************************************************

void McmcFunctions::chi_squared(const vector<real_prec> &data, const vector<real_prec> &model,vector<real_prec> &sigma, real_prec &chis_one){
  /*
    Compute chi square given the data, the model and the inverse
    of the covariance matrix of dimension nXn
  */
  chis_one=0;
  int ncount=0;
#pragma omp parallel for reduction(+:chis_one, ncount)
  for(int i=0;i<model.size();++i)
   if(sigma[i]>0)
      {
          chis_one+=pow(data[i]-model[i],2)/pow(sigma[i],2);
          ncount++;
      }
  chis_one/=(static_cast<double>(ncount));

}

// *******************************************************************************************************

void McmcFunctions::chi_squared_poisson(vector<real_prec> &data, vector<real_prec> &model, real_prec &chis_one){
  /*
    Compute chi square given the data, the model and the inverse
    of the covariance matrix of dimension nXn
  */
  chis_one=0;
  for(int i=0;i<model.size();++i)
      if(model[i]>0)
          chis_one+=data[i]*log(model[i])-model[i];

}


// *******************************************************************************************************
// *******************************************************************************************************
void McmcFunctions::get_loglikelihood_Poisson(vector<vector<real_prec>> &data, vector<vector<real_prec>> &model, real_prec &log_likelihood){
  /*
    Compute chi square given the data, the model and the inverse
    of the covariance matrix of dimension nXn
  */

  log_likelihood=0;
#pragma omp parallel for reduction(+:log_likelihood) collapse(2)
  for(int i=0;i<model.size();++i)
    for(int j=0;j<model[0].size();++j)
      if(model[i][j]>0)
        log_likelihood += (data[i][j]*log(model[i][j])-model[i][j]);

}


// *******************************************************************************************************
// *******************************************************************************************************
void McmcFunctions::chi_squared(vector<real_prec> &data, vector<real_prec> &model,vector<real_prec> &sigma, vector<real_prec>&xvar, real_prec xmin, real_prec xmax, real_prec &chis_one){

  //  Compute chi square given the data, the model and the inverse
 // of the covariance matrix of dimension nXn
  chis_one=0;
  for(int i=0;i<model.size();++i){
    if(xvar[i]<xmax && xvar[i]>=xmin){
      chis_one+=pow( (data[i]-model[i])/sigma[i],2);
    }
  }
}



// ***********************************************************************************************
// ***********************************************************************************************
void McmcFunctions::chi_squared_marginalized_amplitude(vector<real_prec> &data, vector<real_prec> &model,vector< vector<real_prec> > &inv_cova, real_prec &chis_one){
  /*
    Chi square eith analytic marginalization over an amplitud of the model
  */
  chis_one=0;
  int n=model.size();
  vector < vector<real_prec> > inv_cor_aux1(n, vector<real_prec >(n,0));
  vector < vector<real_prec> > inv_cor_aux2(n, vector<real_prec >(n,0));
  vector < vector<real_prec> > inv_cor_aux3(n, vector<real_prec >(n,0));

  for(int i=0;i<n;++i)for(int k=0;k<n;k++)inv_cor_aux1[i][k]=model[i]*model[k];
  for(int i=0;i<n;++i)for(int k=1;k<n;k++)for(int p=0;p<n;p++)inv_cor_aux2[i][k]+=inv_cor_aux1[i][p]*inv_cova[p][k];
  for(int i=0;i<n;++i)for(int k=1;k<n;k++)for(int p=0;p<n;p++)inv_cor_aux3[i][k]+=inv_cova[i][p]*inv_cor_aux2[p][k];
  for(int i=0;i<n;++i)for(int k=1;k<n;k++)chis_one+=model[i]*inv_cova[i][k]*model[k];
  for(int i=0;i<n;++i)for(int k=1;k<n;k++)inv_cor_aux2[i][k]=inv_cova[i][k]-inv_cor_aux3[i][k]/chis_one;
  chis_one=0;
  for(int i=0;i<n;++i)for(int k=0;k<n;k++)chis_one+=data[i]*inv_cor_aux2[i][k]*data[k];

}

// ***********************************************************************************************
// ***********************************************************************************************

void McmcFunctions::chi_squared_marginalized_amplitude(vector<real_prec> &data, vector<real_prec> &model, vector<real_prec> &sigma,  vector<real_prec>&xvar, real_prec xmin, real_prec xmax, real_prec &chis_one){
  /*
    Chi square eith analytic marginalization over an amplitud of the model
  */
    int n=model.size();
    real_prec alpha=0;
    real_prec beta=0;
    real_prec gamma=0;
    for(int i=0;i<n;++i){
      if(xvar[i]<xmax && xvar[i]>=xmin){
        alpha+= data[i]*data[i]/pow(sigma[i],2);
        beta += model[i]*data[i]/pow(sigma[i],2);
        gamma+= model[i]*model[i]/pow(sigma[i],2);
       }
    }
    chis_one=alpha-beta*beta/gamma;
 }

// ***********************************************************************************************
// ***********************************************************************************************

void McmcFunctions::likelihood_full(int lmin, int lmax, vector<matrices>&VM,real_prec &like){
  /*
    Compute chi square given the data, the model and the inverse
    of the covariance matrix of dimension nXn
     ONLY FOR FB analysis
  */
  like=0;
  for(int l=lmin;l<=lmax;++l){
    vector<real_prec>model(VM[l].R.size(),0);
    get_model_FB(VM[l].Step,VM[l].R,VM[l].V,model);
    real_prec chiss;
    chi_squared(VM[l].Cmed,model,VM[l].iCov,chiss);
    like+=0.5*chiss;
  }
}


// *******************************************************************************************************
// *******************************************************************************************************

// Metropolis-Hasting algorithm:
void McmcFunctions::MHalgorithm(real_prec &curr_loglike, real_prec &prop_loglike, int &acc, int &weight_here){
  real_prec MH;

  for(int kl=0;kl<parameters.size();kl++)
   if(fixed_parameters[kl]==1)
    if(parameters[kl]>parameters_max[kl])
      {
       So.message_screen("Warning: parameter ",kl," is taking values above the allowed range.");
       So.message_screen("Proposal set to zero in MHalgorithm");
       cout<<endl;
      }
    for(int kl=0;kl<parameters.size();kl++)
     if(fixed_parameters[kl]==1)
      if(parameters[kl]< parameters_min[kl])
       {
         So.message_screen("Warning: parameter ",kl," is taking values below the allowed range.");
         So.message_screen("Proposal set to zero in MHalgorithm");
         cout<<endl;
       }


    // Get the ra tio between new_posterior/old_posterior: this assumes the same priors
   real_prec posterior_ratio=exp(-0.5*prop_loglike+0.5*curr_loglike);
//   real_prec posterior_ratio=(prop_loglike/curr_loglike);

   So.message_screen("Ratio Proposed/Current posteriors =", posterior_ratio);


   int ip=0;
   while((this->parameters[ip]>=this->parameters_min[ip]) && (this->parameters[ip]<=this->parameters_max[ip]) && (ip<this->n_parameters))
    ip++;
   if(ip<this->n_parameters)
     MH=0;
   else
    MH=min(1.0,static_cast<double>(posterior_ratio));

  double xx=gsl_rng_uniform (r);

  if(xx<MH)
    {
     // Accept this step
      acc++;
      curr_loglike=prop_loglike;
      this->chiss[acc]=prop_loglike;
      this->weight[acc]=num_1;
      weight_here=num_1;
      for(int ip=0;ip<this->parameters.size();++ip)
        this->acc_parameters[ip][acc]=this->parameters[ip];
    }
  else
    {
    // // Do not accept this step, return to then last accepted step
      So.message_screen("**Step rejected**");
      for(int ip=0;ip<this->parameters.size();++ip)
      this->parameters[ip]=this->acc_parameters[ip][acc];
    weight_here=0;
    this->weight[acc]++;
  }

return;
}


// *******************************************************************************************************
// *******************************************************************************************************
void McmcFunctions::write_accepted_models(int j, int acc, int weight_here, string screen){

  int macc=acc-this->weight[acc];
  int wacc=acc-1;
  if(acc>1 && weight_here!=0)
    {
    //   // ***********************
    //   // Write in output file: *
    //   // ***********************
     this->lpout<<wacc<<"\t"<<static_cast<real_prec>(100.0*wacc)/static_cast<real_prec>(j)<<"\t"<<this->chiss[macc]<<"\t"<<this->weight[wacc]<<"\t";
     for(int jk=0;jk<this->parameters.size();jk++)
       lpout<<this->acc_parameters[jk][macc]<<"\t";
     this->lpout<<endl;
     if(screen=="screen")
      {
      // ***********************
    // Write in the screen:  *
    // ***********************
       cout<<endl;
        So.message_screen("Step accepted. Info:");
        So.message_screen("Step = ",j);
        So.message_screen("# of accepted = ",wacc);
        So.message_screen("Acceptance rate = ",100.0*wacc/static_cast<double>(j)," %");
        So.message_screen("log_likelihood = ",this->chiss[macc]);
        So.message_screen("Weight at this point = ",weight[wacc]);
        for(int jk=0;jk<this->n_parameters;jk++)
            So.message_screen("Accepted value",jk," = ",acc_parameters[jk][macc]);
        cout<<endl;
      }
   }
}


// ************************************************************************************************
// ************************************************************************************************

void  McmcFunctions::auto_correlation_mcmc(vector<int> &fixed_par,vector< vector<real_prec> > &acc_param, vector< vector<real_prec> > &auto_corr){
/*
  Computs the auto correlation of the different parameters along the chain
  As input, needs the number of parameters np, the number of accepted models
  AFTER BURN-IN PHASE acc, the msteps, the vector indicating which parameters are to keep fixed, 
  the matrix indicating the parameters accepted along the chain.
  Return the matrix auto_corr (np, acc)
*/
  cout<<"Calculating Auto-correlation"<<endl;
  int np=acc_param.size();
  int acc=acc_param[0].size();
  vector< vector<real_prec> > aux(np, vector<real_prec>(acc,0));
  vector< vector<real_prec> > daux(np, vector<real_prec>(acc,0));
  vector<real_prec> mean(np,0);
  for(int ip=0;ip<np;ip++)
    for(int k=0;k<acc;k=k+this->msteps)
      mean[ip]+=acc_param[ip][k]/((real_prec)acc/this->msteps);
  
  for(int ip=0;ip<np;ip++)
    {
      if(fixed_par[ip]==1)
	{
	  for(int k=0;k<acc;k=k+this->msteps)
	    {
	      for(int i=0;i<acc-k;i=i+this->msteps)
		{
		  aux[ip][k] +=(acc_param[ip][i]-mean[ip])*(acc_param[ip][i+k]-mean[ip]);
		  daux[ip][k]+=pow(acc_param[ip][i]-mean[ip],2);
		}
	      auto_corr[ip][k]=aux[ip][k]/daux[ip][k];
	    }
	}
    }
}

// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************


void McmcFunctions::entropy_mcmc(vector<vector<real_prec> > &pdf, vector<real_prec> &entro){
  // ********************************************************************
  //  Entropy of a MCMC chain for each parameter. Check the definition!
  // ********************************************************************

  for(int i=0;i<entro.size();++i)for(int n=0;n<pdf[0].size();n++)if(pdf[i][n]!=0)entro[i]+=pdf[i][n]*log(pdf[i][n]);
  return;
}

// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************
// ************************************************************************************************


void McmcFunctions::posterior1d(string fname_mean, string fname_pdf, vector<real_prec>  &weight, vector<vector<real_prec> > &acc_par,int perc){
  
  /*
    -------
    Wiki:
    Percentile: is the value of a variable below which a certain percent of observations fall.
    Thus, rounding to two decimals places, −3(sigma) is the 0.13th percentile, −2(sigma) the 2.28th percentile, −1(sigma) the 15.87th percentile, 
    0(sigma) the 50th percentile (both the mean and median of the distribution), +1(sigma) the 84.13th percentile, +2(sigma)
    the 97.72nd percentile, and +3 the 99.87th percentile.
    -------
    
    The corresponding percentil is specified by the integer perc: 
    perc=1 will give the hight of a 68.3 per cent 
    perc=2 gives the  hight of a 97.72 per cent 
    perc=3 will give the hight of a 99.87 per cent confidence intervals.  
    The vector percentil[4] below encodes the percentiles for each value of perc.
    
    **The matrix acc_parameters contain the parameters of the accepted models AFTER BURN-IN PHASE, hence slab=1;
    
    **If you want the posterior (1d, 2d etc) to be shown in a smaller range of parameters (maybe the prior was too broad)!, 
    redefine the minimum and maximum allowed values of the parameters in the par_min and par_max vector
    here or in the function mcmc_settings, 
    before calling this function
    
    
    **If you are trying to recover known parameteres, you can compute the bias with respect to the true parameters
    by passing the fid[] vector  d_max (d_min) is a vector containing the upper (lower) limit for each parameter for the user-given perc value
    
    to be done: 
    skip the passsage of np and define it as an const int np in the parameter file and use it as a global variable
    then initialize all arrays defined here with ={0}
    
    note: in previus version, this function was called errors()
  */

  cout.precision(5);
  //  cout.setf(ios::showpoint);
  // cout.setf(ios::scientific);
  
  /*--------------------------*/
  real_prec aal,pdfmax;
  int slab=0;
  int np=this->parameters_max.size();
  int acc=acc_par[0].size();
  /*--------------------------*/



  /*Define arrays*/
  real_prec percentil[4]={0,0.1587,0.0228, 0.0013};    /*corresponds to 15.87th percentil, 2.28 percentil and 0.13th percentil*/
  vector< vector<real_prec > > cov(np, vector<real_prec>(np,0));
  vector<int> JJ(np,0);



  vector<real_prec> new_acc(np,0);
  vector< vector<real_prec> > e_al_max(np, vector<real_prec> (this->nbin_1D, 0));                                            /*for ech aprameter, computes the cumultive number counts starting from the upper limit*/
  vector< vector<real_prec> > al_n(np, vector<real_prec> (this->nbin_1D, 0));                                        /*1d posterior in nbins for each parameter np */
  
  vector<real_prec> ddel(np,0);
  for(int i=0;i<np;++i)ddel[i]=(this->parameters_max_2dplot[i]-this->parameters_min_2dplot[i])/((real_prec)this->nbin_1D);




  /*
    Compute the marginalized posterior distribution for each parameter in nbins number of bins. Each parameter is weighted with the matrix 
  */
  
      
  for(int ip=0;ip<np;ip++){                                                /*loop over each paraameter*/
    if(this->fixed_parameters[ip]!=0){                                            /*only work with parameers to be fitted */
      for(int j=0;j<this->nbin_1D;++j)al_n[ip][j]=0;                                /*initialize posterior distribution*/
      for(int j=slab;j<acc;j=j+this->msteps){                                    /*loop over the accepted models, one each this->msteps models*/
      if(acc_par[ip][j]!= -999 && acc_par[ip][j]<=this->parameters_max_2dplot[ip] && (acc_par[ip][j]>=this->parameters_min_2dplot[ip]) ){  /*check that the parameter is within the prior or desired range (smaller than the prior)*/
        aal=(int)floor((acc_par[ip][j]-this->parameters_min_2dplot[ip])/((real_prec)ddel[ip]));            /*identify the label of the interval (parmin, par mx) where the model is located*/
        al_n[ip][(int)aal]+=weight[j];
        /*build the weighted posterior distribution*/

	}
      }
    }
  }


  
  /*---------------------------------------------------------------------------------------------------------------------------------------*/  
  /*
    We now normalize the posterior distribution in each parameter to its maximum and write to an output file
    The output is the middle of the bin in the corresponding units and the posterior in that bin.
  */

  
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]!=0){
      string nfname=fname_pdf+"_par_"+to_string(ip)+".txt";
      osali.open(nfname.c_str());
      osali.precision(12);
      osali.setf(ios::showpoint);
      osali.setf(ios::scientific);
      real_prec pd=-1000;

      for(int nb=0;nb<this->nbin_1D;++nb){
        pdfmax=max(pd,al_n[ip][nb]);
        pd=pdfmax;
      }

      for(int nb=0;nb<this->nbin_1D;nb++)osali<<this->parameters_min_2dplot[ip]+(nb+0.5)*ddel[ip]<<"   "<<al_n[ip][nb]/pd<<endl;
      osali.close();
      

      std::cout<<BLUE<<"Wrote 1d posterior in file "<<RESET<<GREEN<<nfname<<RESET<<std::endl;
    }
  }



  /*==========================================================================================================================================
    We compute now the percentiles
    ================================================*/
  /*
    The vector new_acc contains the total (weighted) number of models for each parameter. It will be used
    as an estimate of the *total area* with respect to which the percentiles are identified
  */
  


  fill(new_acc.begin(), new_acc.end(), 0);
  for(int ip=0;ip<np;ip++)for(int j=slab;j<acc;j=j+this->msteps)new_acc[ip]+=weight[j];

  
  /*
    Sup error bars:
    the vector e_al_max[nbin+1] counts the cumulative number of events starting from the upper limit
    and surveying the x-axis towards the mean of the posterior distribution, until it finds the bin at which the 
    cumultive number is a *DESIRED* fraction of the total area, (i.e, new_acc[]): 
    The fraction of the total area depends on the sigma we want to quote as errors.
    At this point, we identify the corresponding
    bin_eff and the mean value of the parameter in that bin, taken just as the value in the middle of the bin (e.g, bin-0.5).
    The value of bin_eff is the average between the bins in which the equality is satisfied
    That value of the parameter corresponds to the maximim estimated by the corresponding percentil
    The error in that direction is the difference between that value and the mean of the distribution.
    
  */

  for(int ip=0;ip<np;ip++){                                                 /*Loop over the parameters allowed to vary in the MC chain*/
    if(this->fixed_parameters[ip]!=0){                                                              
      for(int j=this->nbin_1D-1;j>=0;j--)for(int k=this->nbin_1D-1;k>=j;k--)e_al_max[ip][(this->nbin_1D-1)-j]+=al_n[ip][k]; /*Add events starting from the upper limit towards the minimum*/
    }
  }




  // **************************************************************************************************
  //    Now look for the bin at which the fraction (cumulative) is equal to the corresponding percentil
  // **************************************************************************************************

  for(int ip=0;ip<np;ip++){                                                 /*Loop over the parameters allowed to vary in the MC chain*/
    if(this->fixed_parameters[ip]!=0){
      for(int h=0;h<this->nbin_1D;h++)if(e_al_max[ip][h]<percentil[perc]*(new_acc[ip]) && e_al_max[ip][h+1]>=percentil[perc]*(new_acc[ip])){
	  JJ[ip]=h;                                                          /*  Identify Bin*/
	}  
      this->d_max[ip]=this->parameters_max_2dplot[ip]-JJ[ip]*ddel[ip];                           /*  Compute value of the parameter in that bin*/
    }  
  }


  /*
    inf errors:
    same procedure. Now start counting from the minimum 
    until we find the bin at which the enclosed area equals a fraction of the total area
  */
  for(int i=0;i<np;++i)for(int j=0;j<this->nbin_1D;++j)e_al_max[i][j]=0;
  for(int ip=0;ip<np;ip++){                                                      /*Loop over the parameters allowed to vary in the MC chain*/ 
    if(this->fixed_parameters[ip]!=0){
      for(int j=0;j<this->nbin_1D;++j)for(int k=0;k<=j;k++)e_al_max[ip][j]+=al_n[ip][k];  /*Add events starttin from the upper limit towards the minimum*/
    }
  }


  fill(JJ.begin(), JJ.end(), 0);
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]!=0){
      for(int h=0;h<this->nbin_1D;h++){
	if(e_al_max[ip][h+1]>=percentil[perc]*(new_acc[ip]) && e_al_max[ip][h]<percentil[perc]*(new_acc[ip])){
	  JJ[ip]=h;                                                            /*  Identify Bin*/    
	}  
      }
      this->d_min[ip]=this->parameters_min_2dplot[ip]+JJ[ip]*ddel[ip];                             /*  Compute value of the parameter in that bin*/
    }
  }

  // *********************************************
  // Output     
  // *********************************************
  get_mean(weight,slab,acc_par);                                          /*Compute mean and standard deviation of the distribution*/
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]==0){
      this->d_max[ip]=this->mean_parameters[ip];
      this->d_min[ip]=this->mean_parameters[ip];
    }  
  }
  fname_mean=fname_mean+".txt";
  osali.open(fname_mean.c_str());
  osali.precision(12);
  cout<<BLUE<<"Writing mean, 68.2 percent confidence interval and standard deviation in file "<<RESET<<GREEN<<fname_mean<<endl;
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)osali<<ip<<"\t"<<this->mean_parameters[ip]<<"\t"<<this->d_max[ip]-this->mean_parameters[ip]<<"\t"<<this->mean_parameters[ip]-this->d_min[ip]<<"\t"<<this->stdev_parameters[ip]<<endl;
  osali.close();
  std::cout<<CYAN<<"        Mean values                   +           68.2 CL      -                 Standard deviation" <<endl;
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)cout<<ip<<"\t"<<this->mean_parameters[ip]<<"\t"<<this->d_max[ip]-this->mean_parameters[ip]<<"\t"<<this->mean_parameters[ip]-this->d_min[ip]<<"\t"<<this->stdev_parameters[ip]<<endl;

  cout<<"Latex Style"<<endl;
  std::cout.precision(3);
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)cout<<ip<<"\t"<<this->mean_parameters[ip]<<"^{+"<<this->d_max[ip]-this->mean_parameters[ip]<<"}_{-"<<this->mean_parameters[ip]-this->d_min[ip]<<"}"<<endl;
  

  std::cout<<RESET<<endl;
  
  // ***********************************************************************************************
  // Compute the covariance matrix of the parameters
  // and write it in screen
  // ***********************************************************************************************
  
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]!=0){
      for(int jp=0;jp<np;jp++){
	if(this->fixed_parameters[jp]!=0){
	  for(int i=0;i<acc;++i){
	    cov[ip][jp]+=(acc_par[ip][i]-this->mean_parameters[ip])*(acc_par[jp][i]-this->mean_parameters[jp])/(acc-1);
	  }  
	}
      }
    }
  }
   cout<<"Covariance matrix, correlation, variance,  Delta Chis= error/variance"<<endl;
    for(int ip=0;ip<np;ip++){
     if(this->fixed_parameters[ip]==1){
       for(int jp=0;jp<np;jp++)if(this->fixed_parameters[jp]==1)cout<<cov[ip][jp]<<"\t";
        cout<<endl;
     }
   }

    //          cout<<ip<<"\t"<<jp<<"\t"<<cov[ip][jp]<<"\t"<<(real_prec)cov[ip][jp]/sqrt(cov[ip][ip]*cov[jp][jp])<<"\t"<<sqrt(cov[ip][ip])<<"\t"<<(this->d_max[ip]-this->mean_parameters[ip])/sqrt(cov[ip][ip])<<"\t"<<(this->mean_parameters[ip]-this->d_min[ip])/sqrt(cov[ip][ip])<<endl;

  
  /*==============
    Figure of merit for ellipses
    *================================================*/
  //el error de cada parametro se puede aproximar a error=sqrt(delta Chi)*variancia, 
  //donde delta Chi depende del numero de parametros en los que estemos interesados.
  //Para el CL de 68.2 por ciento y para un parametro, DeltaChi=1, ver las tablas del
  //numerical recipes
  // real_prec scale=2.30;
  // cout<<"  FoM =  "<<1./(3.14159*scale*st[1]*st[2]*sqrt(1.0-pow(cov[1][2]/(st[1]*st[2]),2.0)))<<endl;

  return;

}







// ************************************************************************************************
// ************************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************

//NOTE: MAYBE WE CAN PUT HERE THE 1D DISTRIBUTION AND LET THE FUNCTION ABOVE ONLY FOR THE percentiles
void McmcFunctions::posterior2d(string tdpost, string CR,vector<real_prec>  &weight, vector<vector<real_prec> > &acc_par){
  /*
    BEING
    D= DATA
    X= MODEL
    THIS ROUTINE  COMPUTES THE POSTERIOR DISTRIBUTION OF THE PARAMETERS X GIVEN A DATA D , P(X|D).
    FROM BAYES THEOREM, P(X|D)=P(D|X)P(X)/P(D)
    WHERE
    P(D|X)=THE LIKELIHOOD, THE PROBABILITY OF MEASURE A DATA SET {D} GIVEN A SET OF PARAMETERS {X}. (FORWARD MODELING)
    P(X|D)=IS THE POSTERIOR PROBABILITY OF THE PARAMETERS,  THE PROBABILITY THAT THE PARAMETERS OF A MODEL X TAKE CERTAIN VALUES AFTER DOING THE EXPERIMENT D. 
    THIS LAST QUANTITY IS WHAT WE WANT TO DETERMINE WITH THE MARKOV CHAINS
    P(X)  = PRIOR, THE DEGREE OF KNOWLEDGE OF THE PARAMETERS PRIOR TO THE EXPERIMENT, eG., TAKEN FROM OTHER EXPERIMENT
    P(D)  = EVIDENCE, PLAYING THE ROLE OF NORMALIZATION AND ONLY USEFUL IN MODDEL SELECTION ANALYSIS.
    BAYES' THEOREM -> THE PRIOR IS UPDATED BY THE LIKELIHOOD TO OBTAIN THE POSTERIOR, NORMALIZED BY THE EVIDENCE.
    
    The Metropolis-Hasting algorithm used assumes a Gaussian Likelihood
    which, while being evaluated within the boundaries defined by the prior, gives a sampling of the posterior
    distribution
    

    Given the accpeted models acc_alp, this return the 2d histogram in the array posterior[][], e.g, 
    the two dimensional posterior distribution marginalzed over the rest of the parameters.

  */
  
  //Re-define the range to compute the posterior
  // this->parameters_max[1]=0.27;
  // this->parameters_min[1]=0.235 ;
  // this->parameters_max[2]=0.80;
  // this->parameters_min[2]=0.745;
  

  this->So.message_screen("Building posterior and writing to files");

  real_prec poster_all,lkmin;

  int nll  = 100;        /*Divide the interval of posterior [0,1] in this number of bins*/
  int slab = 0;       /*This means that the models read here are already selected after burn-in phase*/

  int np=this->parameters_max.size();
  int acc=acc_par[0].size();



  // Vector containing the heights for the posterior
  vector<real_prec>sigmalevel;
  sigmalevel.push_back(0.683);
  sigmalevel.push_back(0.90);
  sigmalevel.push_back(0.954);
  sigmalevel.push_back(0.9900);
  sigmalevel.push_back(0.9973);
  sigmalevel.push_back(0.9999);

  /*Define Vectors: Do not forget to initialize them when used within loops*/
  vector<vector<real_prec> >poster(this->nbin_2D,vector<real_prec>(this->nbin_2D,0));
  vector<real_prec>com(nll);
  
  vector<vector<vector<real_prec> > > las;
  las.resize(nll);
  for (int i = 0; i<nll;++i) {
    las[i].resize(this->nbin_2D);
    for (int j=0;j<this->nbin_2D;++j)las[i][j].resize(this->nbin_2D);
  }
  for(int i=0;i<nll;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;k++)las[i][j][k]=0;


  
  vector<real_prec> nppi(nll,0); /*Number of bins that contribute to the nll-th bin of posterior */
  vector<real_prec> cs(nll,0);   /*Create a vector if dimension nll: each element denotes the cumulative enclosed volume */
  vector<int> lsl(7,0);       /*Used to allocate number if bin in posterior at which the enclosed volume is a fraction of that given by different sigmas*/
  vector<real_prec> sf(7,0);     /*Use to allocate heights of posterior for different sigmas*/
  

  vector<real_prec> ddel(np,0);
  for(int i=0;i<np;++i)ddel[i]=(this->parameters_max_2dplot[i]-this->parameters_min_2dplot[i])/((real_prec)this->nbin_2D);
  
  vector<real_prec> lla(nll);                                                             /*Use to define bins of posterior*/
  for(int li=0;li<nll;++li)lla[li]=(real_prec)(li+0.5)/(real_prec)(nll);                      /*Define center of the bins in the posterior */

  
  /* ==========================================================================================================================================
     Compute the 2d posterior probability  distribution by assigning the accepted models to a grid using a Triangular 
     shape cloud mass asssignement scheme and unsing the weights.
     ================================================*/
  
  for(int ip=0;ip<np;ip++){                                 /*Loops over the parameters*/
    for(int jp=0;jp<np;jp++){
//      if(jp>ip){
        if(this->fixed_parameters[ip]==1 && this->fixed_parameters[jp]==1){
          for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)poster[i][j]=0;
	  for(int i=slab;i<acc;i=i+this->msteps){
            if(acc_par[ip][i]!= -999 && (acc_par[ip][i]>=this->parameters_min_2dplot[ip] && acc_par[jp][i]>=this->parameters_min_2dplot[jp]) &&  (acc_par[ip][i]<this->parameters_max_2dplot[ip] && acc_par[jp][i]<this->parameters_max_2dplot[jp])){
              real_prec x=acc_par[ip][i]-this->parameters_min_2dplot[ip];
              real_prec y=acc_par[jp][i]-this->parameters_min_2dplot[jp];
              // The weight for ip is the same for ij, these are models
              this->mass_asg2d(this->MAS,x,y,ddel[ip],ddel[jp],weight[i], poster);
            }
	  }


        /*buscamos el maximo y normalizamos*/
	lkmin=-100.0;
        real_prec lkk=0.0;
        for(int i=0;i<this->nbin_2D;++i){
          for(int j=0;j<this->nbin_2D;++j){
	    lkk=max(lkmin,poster[i][j]);
	    lkmin=lkk;
	  }  
	}         
	
        for(int i=0;i<poster.size();++i)for(int j=0;j<poster[0].size();++j)poster[i][j]/=(real_prec)lkk;     /*Normalize the posterior to 1 */



	ofstream lout;
	string ntdpost=tdpost+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
	lout.open(ntdpost.c_str());
	lout.precision(12);
        for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)lout<<this->parameters_min_2dplot[ip]+(i+0.5)*(this->parameters_max_2dplot[ip]-this->parameters_min_2dplot[ip])/((real_prec)this->nbin_2D)<<"   "<<this->parameters_min_2dplot[jp]+(j+0.5)*(this->parameters_max_2dplot[jp]-this->parameters_min_2dplot[jp])/((real_prec)this->nbin_2D)<<"   "<<poster[i][j]<<endl;
    std::cout<<BLUE<<"Wrote 2d posterior in file "<<RESET<<GREEN<<ntdpost<<RESET<<std::endl;
	lout.close();
	

	//===============================================================================
	/*
	  Estimates of the levels defining confidence regions
	*/
		
	poster_all=0;                                                           
        for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)poster_all+=poster[i][j];      /*total volumen  enclosed by the 2d de posterior*/
	
	fill(com.begin(), com.end(), 0);
	fill(nppi.begin(), nppi.end(), 0);

        for(int li=0;li<nll;++li)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)las[li][i][j]=0;

	for(int li=0;li<nll;++li){                                                   /*Loop over bins of the posterior, normalized between 0 and 1*/
          for(int i=0;i<this->nbin_2D;++i){                                                    /*Loop over bins of parameter 1*/
            for(int j=0;j<this->nbin_2D;++j){                                                  /*Loop over bins of parameter 2*/
	      if(li<nll-1){
                if(poster[i][j]>= lla[li]-0.5/(real_prec)nll && poster[i][j]<lla[li]+0.5/(real_prec)nll){                     /*At this cell in parameter space, read the posterior  and */
		  las[li][i][j]=poster[i][j];                                            /*fill tensor las[][][]*/
		  nppi[li]++;                                                            /*Compute the number of cells that contribute to this posterior*/
		  com[li]     +=poster[i][j];                                            /*Compute the volume enclosed by the posterior within this specific bin*/
		}		                                                         /*The bins here were [)*/
	      }
	      if(li==nll-1){
                if(poster[i][j]>= lla[nll-1]-0.5/(real_prec)nll && poster[i][j]<=lla[nll-1]+0.5/(real_prec)nll){                   /*Take care of the last bin for it must be [] such that we count everything*/
	       	  las[li][i][j]=poster[i][j];                                            
	       	  nppi[li]++;                                                            
	       	  com[li]     +=poster[i][j];                                            
	       	}		
	      }
	    }
	  }
	}
	
	// ******************************************************
	// Checking for consistency*/
	// al=0;
	// for(int li=0;li<nll;li++)al+=com[li];
	// cout<<al<<" must be equal to "<<poster_all<<endl;
	// ******************************************************
	
	
	fill(cs.begin(), cs.end(), 0);
	for(int li=0;li<nll;li++)for(int lj=0;lj<nll;++lj)if(lj<=li)cs[li]+=com[lj];  /*In each bin of posterior, compute the cumulative volume enclosed until that bin*/
	for(int sl=0;sl<6;sl++)lsl[sl]=0;
	
	
	for(int sl=0;sl<6;sl++){                                                                          /*Loop over the values defining the confidence regions*/
	  for(int li=0;li<nll;li++){                                                                      /*Loop over the bins in posterior: find the bin (for each sigma)*/ 
	                                                                                                  /* at which the enclosed*/
	    if( (1.0-sigmalevel[sl])*poster_all>=cs[li]  && (1.0-sigmalevel[sl])*poster_all<cs[li+1]){     /*volume corresponds to a fraction (1-sigmalevel) of the total*/
	      lsl[sl]=li;                                                                                  /*enclosed volume poster_all. The term (1-sigmalevel) means that*/
	    }                                                                                              /*we are staring from below upwards and once the bin is found, a fraction*/
      	  }                                                                                                /*sigmalevel remains above that bin*/ 
	}
	
	
	// ***********************************************************************************************
	// Computing the height in posterior defining cofidence regions:
	// for each sigma, compute the average posterior at the 
	// corresponding bin lsl[sl] found in the previus step: that is, compute the sum of the posterior
	// read from the tensor las[lsls[sl]][i][j] divided by the number
	// of parameter_cells that contributed to that particular bin of posterior nppi[lsl[sl]]
	// ***********************************************************************************************


	fill(sf.begin(), sf.end(), 0);
	for(int sl=0;sl<6;sl++){
	  //	  cout<<nppi[lsl[sl]]<<endl;
          for(int i=0;i<this->nbin_2D;++i){
            for(int j=0;j<this->nbin_2D;++j){
	      sf[sl]+=las[lsl[sl]][i][j]/nppi[lsl[sl]]; 
	    }
	  }
	}
	/*  Output of height of normalized posterior distribution defining different confidence regions
	    One line, six columns. A file for every pair of parameters*/
	
        ofstream kout;
	string nCR=CR+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
	kout.open(nCR.c_str());
	kout.precision(12);
        kout.setf(ios::showpoint);
        kout.setf(ios::scientific);


        for(int sl=0;sl<6;sl++)kout<<sf[sl]<<"\t"; kout<<"\t"<<endl;
	this->So.message_screen("Wrote CR in file",nCR);
	kout.close();
        }
     //}

   }
  }
  std:: cout<<GREEN<<"Posterior distributions succesfully written!"<<RESET<<std::endl;
  return;
}

// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************



// This does the same as pòsterior2d,  taking the 2d historgram as input.
// The histogram is expected to be normalized to its maximum
void McmcFunctions::get_contour_levels(string CR,vector<vector<real_prec> > &poster){
  /*    Given the posterior or a 2D hisgoram normaized to its maximim, this gives the height for differnet confidence intervales
  */
  
  this->So.message_screen("Building probability contours and writing to files");
  this->nbin_2D=poster.size();
  real_prec poster_all;
  int nll  = 100;  /*Divide the interval of posterior [0,1] in this number of bins*/
  real_prec min_log_posterior = -10;  // El valor minimo de la posterior será siempre cero, Acá coloco un log 0 = -5 para andar
  real_prec max_log_posterior = 0.0; // La posterior esta normalizada,
  real_prec delta_log_posterior = (max_log_posterior-min_log_posterior)/static_cast<real_prec>(nll);
  int np=2; // Number of parameters, so far I have only two for this is mainly used for the delta delta scatter plot
  
  // Vector containing the heights for the posterior
  vector<real_prec>sigmalevel;
  sigmalevel.push_back(0.683);
  sigmalevel.push_back(0.90);
  sigmalevel.push_back(0.954);
  sigmalevel.push_back(0.9900);
  sigmalevel.push_back(0.9973);
  sigmalevel.push_back(0.9999);

  /*Define Vectors: Do not forget to initialize them when used within loops*/
  vector<real_prec>com(nll);
  
  vector<vector<vector<real_prec> > > las;
  las.resize(nll);
  
  for (int i = 0; i<nll;++i)
  {
    las[i].resize(this->nbin_2D);
    for (int j=0;j<this->nbin_2D;++j)
      las[i][j].resize(this->nbin_2D);
  }
  for(int i=0;i<nll;++i)
    for(int j=0;j<this->nbin_2D;++j)
      for(int k=0;k<this->nbin_2D;k++)
        las[i][j][k]=0;
  
  vector<real_prec> nppi(nll,0); /*Number of bins that contribute to the nll-th bin of posterior */
  vector<real_prec> cs(nll,0);   /*Create a vector if dimension nll: each element denotes the cumulative enclosed volume */
  vector<int> lsl(7,0);       /*Used to allocate number if bin in posterior at which the enclosed volume is a fraction of that given by different sigmas*/
  vector<real_prec> sf(7,0);     /*Use to allocate heights of posterior for different sigmas*/
  
    
  for(int ip=0;ip<np;ip++)
    {                                 /*Loops over the parameters*/
      for(int jp=0;jp<np;jp++)
	{
	  /*
	    Estimates of the levels defining confidence regions
	  */
	  
	  poster_all=0;
	  for(int i=0;i<this->nbin_2D;++i)
	    for(int j=0;j<this->nbin_2D;++j)
	      poster_all+=poster[i][j];      /*total volumen  enclosed by the 2d de posterior*/
	  
	  fill(com.begin(), com.end(), 0);
	  fill(nppi.begin(), nppi.end(), 1); // Initialize it with 1
	  
	  for(int li=0;li<nll;++li)
	    for(int i=0;i<this->nbin_2D;++i)
	      for(int j=0;j<this->nbin_2D;++j)las[li][i][j]=0;
	  
	  
	  for(int i=0;i<this->nbin_2D;++i)
	    {                                                    /*Loop over bins of parameter 1*/
	      for(int j=0;j<this->nbin_2D;++j)
		{
		  real_prec lpost = log10(poster[i][j]);
		  if(lpost >= min_log_posterior && lpost<=max_log_posterior)
		    {
		      /*Loop over bins of parameter 2*/
		      int li = floor((lpost-min_log_posterior)/delta_log_posterior);    // get bin in the posterior
		      if (li==nll) li--;
		      las[li][i][j]=static_cast<real_prec>(poster[i][j]);                                            /*fill tensor las[][][]*/
		      nppi[li]++;      /*Compute the number of cells that contribute to this posterior*/
		      com[li]+=poster[i][j];   /*Compute the volume enclosed by the posterior within this specific bin i-j*/
		    }		                                                         	    
		}
	    }
	  
	  
	  fill(cs.begin(), cs.end(), 0);
	  for(int li=0;li<nll;li++)
	    for(int lj=0;lj<=li;++lj)
	      cs[li]+=com[lj];  /*In each bin of posterior, compute the cumulative volume enclosed until that bin*/
	  
	  for(int sl=0;sl<sigmalevel.size();sl++)
	    lsl[sl]=0;
	  
	  for(int sl=0;sl<sigmalevel.size();sl++) /*Loop over the values defining the confidence regions*/
	    {                           
	      for(int li=0;li<nll;li++)
		{   /*Loop over the bins in posterior: find the bin (for each sigma)*/                                                                   
		  /* at which the enclosed*/
		  if( (1.0-sigmalevel[sl])*poster_all>=cs[li]  && (1.0-sigmalevel[sl])*poster_all<cs[li+1]){     /*volume corresponds to a fraction (1-sigmalevel) of the total*/
		    lsl[sl]=li;                                                                                  /*enclosed volume poster_all. The term (1-sigmalevel) means that*/
		  }                                                                                              /*we are staring from below upwards and once the bin is found, a fraction*/
		}                                                                                                /*sigmalevel remains above that bin*/
	    }
	  
	  // ***********************************************************************************************
	  // Computing the height in posterior defining cofidence regions:
	  // for each sigma, compute the average posterior at the
	  // corresponding bin lsl[sl] found in the previus step: that is, compute the sum of the posterior
	  // read from the tensor las[lsls[sl]][i][j] divided by the number
	  // of parameter_cells that contributed to that particular bin of posterior nppi[lsl[sl]]
	  // ***********************************************************************************************
	  
	  
	  fill(sf.begin(), sf.end(), 0);
	  for(int sl=0;sl<sigmalevel.size();sl++)// loop sobre los distintos valores de sigma que queremos escribir
	    {
	      int index=lsl[sl];
	      int number = nppi[index] == 0 ? 1: nppi[index];
	      for(int i=0;i<this->nbin_2D;++i)
		for(int j=0;j<this->nbin_2D;++j)
		  sf[sl]+=las[index][i][j]/static_cast<real_prec>(number);
	    }
	  /*  Output of height of normalized posterior distribution defining different confidence regions
	      One line, six columns. A file for every pair of parameters*/
	  
	  if(ip<jp)
	    {
	      ofstream kout;
	      string nCR=CR+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
	      kout.open(nCR.c_str());
	      kout.precision(12);
	      kout.setf(ios::showpoint);
	      kout.setf(ios::scientific);
	      
	      for(int sl=0;sl<sigmalevel.size();sl++)kout<<sf[sl]<<"\t"; kout<<"\t"<<endl;
	      So.message_screen("Wrote CR in file", nCR);
	      kout.close();
	    }
	}
    }
  
  //  std:: cout<<GREEN<<"Posterior distributions succesfully written!"<<RESET<<std::endl;
  return;
}
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************

void McmcFunctions::get_contour_levels(string CR,ULONG Nbins, vector<real_prec> &poster){
  /*    Given the posterior or a 2D hisgoram normaized to its maximim, this gives the height for differnet confidence intervales
  */
  
  this->So.message_screen("Building probability contours and writing to files");
  this->nbin_2D=Nbins;
  real_prec poster_all;
  int nll  = 100;  /*Divide the interval of posterior [0,1] in this number of bins*/
  real_prec min_log_posterior = -10;  // El valor minimo de la posterior será siempre cero, Acá coloco un log 0 = -5 para andar
  real_prec max_log_posterior = 0.0; // La posterior esta normalizada,
  real_prec delta_log_posterior = (max_log_posterior-min_log_posterior)/static_cast<real_prec>(nll);
  int np=2; // Number of parameters, so far I have only two for this is mainly used for the delta delta scatter plot
  
  // Vector containing the heights for the posterior
  vector<real_prec>sigmalevel;
  sigmalevel.push_back(0.683);
  sigmalevel.push_back(0.90);
  sigmalevel.push_back(0.954);
  sigmalevel.push_back(0.9900);
  sigmalevel.push_back(0.9973);
  sigmalevel.push_back(0.9999);

  /*Define Vectors: Do not forget to initialize them when used within loops*/
  vector<real_prec>com(nll);
  
  vector<vector<vector<real_prec> > > las;
  las.resize(nll);
  
  for (int i = 0; i<nll;++i)
  {
    las[i].resize(this->nbin_2D);
    for (int j=0;j<this->nbin_2D;++j)
      las[i][j].resize(this->nbin_2D);
  }
  for(int i=0;i<nll;++i)
    for(int j=0;j<this->nbin_2D;++j)
      for(int k=0;k<this->nbin_2D;k++)
        las[i][j][k]=0;
  
  vector<real_prec> nppi(nll,0); /*Number of bins that contribute to the nll-th bin of posterior */
  vector<real_prec> cs(nll,0);   /*Create a vector if dimension nll: each element denotes the cumulative enclosed volume */
  vector<int> lsl(7,0);       /*Used to allocate number if bin in posterior at which the enclosed volume is a fraction of that given by different sigmas*/
  vector<real_prec> sf(7,0);     /*Use to allocate heights of posterior for different sigmas*/
  
    
  for(int ip=0;ip<np;ip++)
    {                                 /*Loops over the parameters*/
      for(int jp=0;jp<np;jp++)
  {
    /*
      Estimates of the levels defining confidence regions
    */
    
    poster_all=0;
    for(int i=0;i<this->nbin_2D;++i)
      for(int j=0;j<this->nbin_2D;++j)
        poster_all+=poster[index_2d(i,j,Nbins)];      /*total volumen  enclosed by the 2d de posterior*/
    
    fill(com.begin(), com.end(), 0);
    fill(nppi.begin(), nppi.end(), 1); // Initialize it with 1
    
    for(int li=0;li<nll;++li)
      for(int i=0;i<this->nbin_2D;++i)
        for(int j=0;j<this->nbin_2D;++j)las[li][i][j]=0;
    
    
    for(int i=0;i<this->nbin_2D;++i)
      {                                                    /*Loop over bins of parameter 1*/
        for(int j=0;j<this->nbin_2D;++j)
    {
      real_prec lpost = log10(poster[index_2d(i,j,Nbins)]);
      if(lpost >= min_log_posterior && lpost<=max_log_posterior)
        {
          /*Loop over bins of parameter 2*/
          int li = floor((lpost-min_log_posterior)/delta_log_posterior);    // get bin in the posterior
          if (li==nll) li--;
          las[li][i][j]=static_cast<real_prec>(poster[index_2d(i,j,Nbins)]);                                            /*fill tensor las[][][]*/
          nppi[li]++;      /*Compute the number of cells that contribute to this posterior*/
          com[li]+=poster[index_2d(i,j,Nbins)];   /*Compute the volume enclosed by the posterior within this specific bin i-j*/
        }                                                                 
    }
      }
    
    
    fill(cs.begin(), cs.end(), 0);
    for(int li=0;li<nll;li++)
      for(int lj=0;lj<=li;++lj)
        cs[li]+=com[lj];  /*In each bin of posterior, compute the cumulative volume enclosed until that bin*/
    
    for(int sl=0;sl<sigmalevel.size();sl++)
      lsl[sl]=0;
    
    for(int sl=0;sl<sigmalevel.size();sl++) /*Loop over the values defining the confidence regions*/
      {                           
        for(int li=0;li<nll;li++)
    {   /*Loop over the bins in posterior: find the bin (for each sigma)*/                                                                   
      /* at which the enclosed*/
      if( (1.0-sigmalevel[sl])*poster_all>=cs[li]  && (1.0-sigmalevel[sl])*poster_all<cs[li+1]){     /*volume corresponds to a fraction (1-sigmalevel) of the total*/
        lsl[sl]=li;                                                                                  /*enclosed volume poster_all. The term (1-sigmalevel) means that*/
      }                                                                                              /*we are staring from below upwards and once the bin is found, a fraction*/
    }                                                                                                /*sigmalevel remains above that bin*/
      }
    
    // ***********************************************************************************************
    // Computing the height in posterior defining cofidence regions:
    // for each sigma, compute the average posterior at the
    // corresponding bin lsl[sl] found in the previus step: that is, compute the sum of the posterior
    // read from the tensor las[lsls[sl]][i][j] divided by the number
    // of parameter_cells that contributed to that particular bin of posterior nppi[lsl[sl]]
    // ***********************************************************************************************
    
    
    fill(sf.begin(), sf.end(), 0);
    for(int sl=0;sl<sigmalevel.size();sl++)// loop sobre los distintos valores de sigma que queremos escribir
      {
        int index=lsl[sl];
        int number = nppi[index] == 0 ? 1: nppi[index];
        for(int i=0;i<this->nbin_2D;++i)
    for(int j=0;j<this->nbin_2D;++j)
      sf[sl]+=las[index][i][j]/static_cast<real_prec>(number);
      }
    /*  Output of height of normalized posterior distribution defining different confidence regions
        One line, six columns. A file for every pair of parameters*/
    
    if(ip<jp)
      {
        ofstream kout;
        string nCR=CR+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
        kout.open(nCR.c_str());
        kout.precision(12);
        kout.setf(ios::showpoint);
        kout.setf(ios::scientific);
        
        for(int sl=0;sl<sigmalevel.size();sl++)kout<<sf[sl]<<"\t"; kout<<"\t"<<endl;
        So.message_screen("Wrote CR in file", nCR);
        kout.close();
      }
  }
    }
  
  //  std:: cout<<GREEN<<"Posterior distributions succesfully written!"<<RESET<<std::endl;
  return;
}
// ********************************************************************************************
// ********************************************************************************************
// ********************************************************************************************



//NOTE: MAYBE WE CAN PUT HERE THE 1D DISTRIBUTION AND LET THE FUNCTION ABOVE ONLY FOR THE percentiles
void McmcFunctions::posterior2d_combined_experiments(string fname_mean, string fname_pdf, string tdpost, string CR, experiments experiments){
  /*
    BEING
    D= DATA
    X= MODEL
    THIS ROUTINE  COMPUTES THE POSTERIOR DISTRIBUTION OF THE PARAMETERS X GIVEN A DATA D , P(X|D).
    FROM BAYES THEOREM, P(X|D)=P(D|X)P(X)/P(D)
    WHERE
    P(D|X)=THE LIKELIHOOD, THE PROBABILITY OF MEASURE A DATA SET {D} GIVEN A SET OF PARAMETERS {X}. (FORWARD MODELING)
    P(X|D)=IS THE POSTERIOR PROBABILITY OF THE PARAMETERS,  THE PROBABILITY THAT THE PARAMETERS OF A MODEL X TAKE CERTAIN VALUES AFTER DOING THE EXPERIMENT D.
    THIS LAST QUANTITY IS WHAT WE WANT TO DETERMINE WITH THE MARKOV CHAINS
    P(X)  = PRIOR, THE DEGREE OF KNOWLEDGE OF THE PARAMETERS PRIOR TO THE EXPERIMENT, eG., TAKEN FROM OTHER EXPERIMENT
    P(D)  = EVIDENCE, PLAYING THE ROLE OF NORMALIZATION AND ONLY USEFUL IN MODDEL SELECTION ANALYSIS.
    BAYES' THEOREM -> THE PRIOR IS UPDATED BY THE LIKELIHOOD TO OBTAIN THE POSTERIOR, NORMALIZED BY THE EVIDENCE.

    The Metropolis-Hasting algorithm used assumes a Gaussian Likelihood
    which, while being evaluated within the boundaries defined by the prior, gives a sampling of the posterior
    distribution


    If you want the posterior (1d, 2d etc) to be shown in a smaller range of parameters (maybe the prior was too broad)!,
    redefine the minimum and maximum allowed values of the parameters in the par_min and par_max vector, before calling this function


    Given the accpeted models acc_alp, this return the 2d histogram in the array posterior[][], e.g,
    the two dimensional posterior distribution marginalzed over the rest of the parameters.

    note: in previus versions, this was called histo_2d.

  */


  std::cout<<YELLOW<<"Building posterior and writing to files"<<RESET<<std::endl;

  real_prec poster_all,lkmin;

  int nll  = 25; // Divide the interval of posterior [0,1] in this number of bins
  int slab = 0;    // This means that the models read here are already selected after burn-in phase, so the cut hs been already done beforea calling this function



  // Vector containing the heights for the posterior
  vector<real_prec>sigmalevel;
  sigmalevel.push_back(0.683);
  sigmalevel.push_back(0.90);
  sigmalevel.push_back(0.954);
  sigmalevel.push_back(0.9900);
  sigmalevel.push_back(0.9973);
  sigmalevel.push_back(0.9999);



  int np=this->parameters_max.size();


  int acc1=experiments.acc_par1[0].size();
  int acc2=experiments.acc_par2[0].size();
  int acc3=experiments.acc_par3[0].size();
  int acc4=experiments.acc_par4[0].size();
  int acc5=experiments.acc_par5[0].size();
  int acc6=experiments.acc_par6[0].size();

  vector<vector<real_prec> >poster(this->nbin_2D,vector<real_prec>(this->nbin_2D,0));

  // Define vectors to get the confidence levels
  vector<real_prec>com(nll);
  vector<vector<vector<real_prec> > > las;
  las.resize(nll);
  for (int i = 0; i<nll;++i) {
    las[i].resize(this->nbin_2D);
    for (int j=0;j<this->nbin_2D;++j)las[i][j].resize(this->nbin_2D);
  }
  for(int i=0;i<nll;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;k++)las[i][j][k]=0;


  vector<real_prec> nppi(nll,0); /*Number of bins that contribute to the nll-th bin of posterior */
  vector<real_prec> cs(nll,0);   /*Create a vector if dimension nll: each element denotes the cumulative enclosed volume */
  vector<int> lsl(7,0);       /*Used to allocate number if bin in posterior at which the enclosed volume is a fraction of that given by different sigmas*/
  vector<real_prec> sf(7,0);     /*Use to allocate heights of posterior for different sigmas*/

  vector<real_prec> ddel(np,0);
  for(int i=0;i<np;++i)ddel[i]=(this->parameters_max_2dplot[i]-this->parameters_min_2dplot[i])/((real_prec)this->nbin_2D);

  vector<real_prec> lla(nll);                                                             /*Use to define bins of posterior*/
  for(int li=0;li<nll;++li)lla[li]=(real_prec)(li+0.5)/(real_prec)(nll);                      /*Define center of the bins in the posterior */


  // ==========================================================================================================================================
  /*  DEBEMOS ARMAR UN VECTOR DE CUATRO DIMENSIONES CON LONGITUD NBINS2D PARA METER LA JOINT POSTERIOR DE TODOS LOS PARAMEROS
   *  DE CADA UNO DE LOS EXPERIMENTOS, PARA LUEGO MULTIPLICARLOS Y COMENZAR A MARGINALIZAR DESDE LAS PROBABILIDIDADES, 
   * NO DESDE LA DISTRIBCUION DE OBJETOS. NO USAMOS LA INFO SOBRE EL BIAS, PUES PARA CADA BIEN SERA TRATADO COMO NUISSANCE
   * ESTO EVITA QUE AL MULTIPLICAR CON EL B INCLUIDO COMO PARAMETRO SALGA TODO CERO. PIENSA CHE PARA CADA Cl EL B SERA DISTINTO, DANDOI REGIONES DISYUNTAS.
   */
  vector<vector<vector<vector< real_prec > > > >  ALL(this->nbin_2D,vector<vector<vector<real_prec> > >(this->nbin_2D,vector<vector<real_prec> >(this->nbin_2D,vector<real_prec>(this->nbin_2D,0))));
  vector<vector<vector<vector< real_prec > > > > ALL1(this->nbin_2D,vector<vector<vector<real_prec> > >(this->nbin_2D,vector<vector<real_prec> >(this->nbin_2D,vector<real_prec>(this->nbin_2D,0))));


  int p1,p2,p3,p5;

  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]=0;
  cout<<"Mes 11"<<endl;
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;

  for(int i=slab;i<acc1;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par1[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par1[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par1[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par1[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par1[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par1[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par1[3][i]>=this->parameters_min_2dplot[2] && experiments.acc_par1[2][i]<this->parameters_max_2dplot[3])p3=floor((experiments.acc_par1[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par1[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par1[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par1[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par1[i];
  }
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]=ALL1[i][j][k][l];


  cout<<"Mes 22"<<endl;
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;
  for(int i=slab;i<acc2;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par2[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par2[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par2[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par2[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par2[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par2[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par2[2][i]>=this->parameters_min_2dplot[2] && experiments.acc_par2[2][i]<this->parameters_max_2dplot[2])p3=floor((experiments.acc_par2[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par2[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par2[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par2[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par2[i];
  }

  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]*=ALL1[i][j][k][l];


  cout<<"Mes 33"<<endl;
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;
  for(int i=slab;i<acc3;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par3[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par3[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par3[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par3[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par3[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par3[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par3[2][i]>=this->parameters_min_2dplot[2] && experiments.acc_par3[3][i]<this->parameters_max_2dplot[2])p3=floor((experiments.acc_par3[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par3[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par3[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par3[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par3[i];
  }
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]*=ALL1[i][j][k][l];


  cout<<"Mes 23"<<endl;
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;
  for(int i=slab;i<acc4;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par4[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par4[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par4[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par4[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par4[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par4[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par4[2][i]>=this->parameters_min_2dplot[2] && experiments.acc_par4[2][i]<this->parameters_max_2dplot[2])p3=floor((experiments.acc_par4[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par4[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par4[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par4[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par4[i];
  }
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]*=ALL1[i][j][k][l];

  cout<<"Mes 12"<<endl;

  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;
  for(int i=slab;i<acc5;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par5[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par5[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par5[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par5[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par5[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par5[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par5[2][i]>=this->parameters_min_2dplot[2] && experiments.acc_par5[2][i]<this->parameters_max_2dplot[2])p3=floor((experiments.acc_par5[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par5[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par5[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par5[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par5[i];
  }
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]*=ALL1[i][j][k][l];





  cout<<"Mes 13"<<endl;
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL1[i][j][k][l]=0;
  for(int i=slab;i<acc5;i=i+this->msteps){// Loop over the models in the merged chains
    if(experiments.acc_par6[0][i]>=this->parameters_min_2dplot[0] && experiments.acc_par6[0][i]<this->parameters_max_2dplot[0])p1=floor((experiments.acc_par6[0][i]-this->parameters_min_2dplot[0])/ddel[0]);//fb
    if(experiments.acc_par6[1][i]>=this->parameters_min_2dplot[1] && experiments.acc_par6[1][i]<this->parameters_max_2dplot[1])p2=floor((experiments.acc_par6[1][i]-this->parameters_min_2dplot[1])/ddel[1]);//Om
    if(experiments.acc_par6[2][i]>=this->parameters_min_2dplot[2] && experiments.acc_par6[2][i]<this->parameters_max_2dplot[2])p3=floor((experiments.acc_par6[2][i]-this->parameters_min_2dplot[2])/ddel[2]);//A
    if(experiments.acc_par6[4][i]>=this->parameters_min_2dplot[4] && experiments.acc_par6[4][i]<this->parameters_max_2dplot[4])p5=floor((experiments.acc_par6[4][i]-this->parameters_min_2dplot[4])/ddel[4]);//S8
    ALL1[p1][p2][p3][p5]+=experiments.weight_par6[i];
  }
  for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)ALL[i][j][k][l]*=ALL1[i][j][k][l];




  for(int ip=0;ip<np;++ip){
    for(int jp=0;jp<np;++jp){

      if(ip!=jp && ip!=2 && jp!=2){  //EXCLUYO LA INFO DE A


       for(int i=0;i<poster.size();++i)for(int j=0;j<poster[0].size();++j)poster[i][j]=0;     /*Normalize the posterior to 1 */

	if(ip==0 && jp==1)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[i][j][k][l];
	else if(ip==0 && jp==3)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[i][k][j][l];
	else if(ip==0 && jp==4)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[i][k][l][j];
	else if(ip==1 && jp==0)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[j][i][k][l];
	else if(ip==1 && jp==3)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][i][j][l];
	else if(ip==1 && jp==4)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][i][l][j];
	else if(ip==3 && jp==0)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[j][k][i][l];
	else if(ip==3 && jp==1)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][j][i][l];
	else if(ip==3 && jp==4)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][l][i][j];
	else if(ip==4 && jp==0)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[j][k][l][i];
	else if(ip==4 && jp==1)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][j][l][i];
	else if(ip==4 && jp==3)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)poster[i][j]+=ALL[k][l][j][i];




	// buscamos el maximo de la full y normalizamos
        lkmin=-1000.0;
        real_prec lkk=0.0;
        for(int i=0;i<this->nbin_2D;++i){
          for(int j=0;j<this->nbin_2D;++j){
	    lkk=max(lkmin,poster[i][j]);
            lkmin=lkk;
          }
        }

        cout<<"Maximum likelihood = "<<lkk<<endl;
        for(int i=0;i<poster.size();++i)for(int j=0;j<poster[0].size();++j)poster[i][j]/=(real_prec)lkk;     /*Normalize the posterior to 1 */



        ofstream lout;
        string ntdpost=tdpost+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
        lout.open(ntdpost.c_str());
        lout.precision(12);
        for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)lout<<this->parameters_min_2dplot[ip]+(i+0.5)*(this->parameters_max_2dplot[ip]-this->parameters_min_2dplot[ip])/((real_prec)this->nbin_2D)<<"   "<<this->parameters_min_2dplot[jp]+(j+0.5)*(this->parameters_max_2dplot[jp]-this->parameters_min_2dplot[jp])/((real_prec)this->nbin_2D)<<"   "<<poster[i][j]<<endl;
	So.message_screen("Wrote 2d posterior in file",ntdpost);
        lout.close();


        //===============================================================================
        /*
          Estimates of the levels defining confidence regions
        */

        poster_all=0;
        for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)poster_all+=poster[i][j];      /*total volumen  enclosed by the 2d de posterior*/

        fill(com.begin(), com.end(), 0);
        fill(nppi.begin(), nppi.end(), 0);

        for(int li=0;li<nll;++li)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)las[li][i][j]=0;

        for(int li=0;li<nll;++li){                                                   /*Loop over bins of the posterior, normalized between 0 and 1*/
          for(int i=0;i<this->nbin_2D;++i){                                                    /*Loop over bins of parameter 1*/
            for(int j=0;j<this->nbin_2D;++j){                                                  /*Loop over bins of parameter 2*/
              if(li<nll-1){
                if(poster[i][j]>= lla[li]-0.5/(real_prec)nll && poster[i][j]<lla[li]+0.5/(real_prec)nll){                     /*At this cell in parameter space, read the posterior  and */
                  las[li][i][j]=poster[i][j];                                            /*fill tensor las[][][]*/
                  nppi[li]++;                                                            /*Compute the number of cells that contribute to this posterior*/
                  com[li]     +=poster[i][j];                                            /*Compute the volume enclosed by the posterior within this specific bin*/
                }		                                                         /*The bins here were [)*/
              }
              if(li==nll-1){
                if(poster[i][j]>= lla[nll-1]-0.5/(real_prec)nll && poster[i][j]<=lla[nll-1]+0.5/(real_prec)nll){                   /*Take care of the last bin for it must be [] such that we count everything*/
                  las[li][i][j]=poster[i][j];
                  nppi[li]++;
                  com[li]     +=poster[i][j];
                }
              }
            }
          }
        }

        // ******************************************************
        // Checking for consistency*/
        // al=0;
        // for(int li=0;li<nll;li++)al+=com[li];
        // cout<<al<<" must be equal to "<<poster_all<<endl;
        // ******************************************************


        fill(cs.begin(), cs.end(), 0);
        for(int li=0;li<nll;li++)for(int lj=0;lj<nll;++lj)if(lj<=li)cs[li]+=com[lj];  /*In each bin of posterior, compute the cumulative volume enclosed until that bin*/
        for(int sl=0;sl<6;sl++)lsl[sl]=0;


        for(int sl=0;sl<6;sl++){                                                                          /*Loop over the values defining the confidence regions*/
          for(int li=0;li<nll;li++){                                                                      /*Loop over the bins in posterior: find the bin (for each sigma)*/
                                                                                                          /* at which the enclosed*/
            if( (1.0-sigmalevel[sl])*poster_all>=cs[li]  && (1.0-sigmalevel[sl])*poster_all<cs[li+1]){     /*volume corresponds to a fraction (1-sigmalevel) of the total*/
              lsl[sl]=li;                                                                                  /*enclosed volume poster_all. The term (1-sigmalevel) means that*/
            }                                                                                              /*we are staring from below upwards and once the bin is found, a fraction*/
          }                                                                                                /*sigmalevel remains above that bin*/
        }


        // ***********************************************************************************************
        // Computing the height in posterior defining cofidence regions:
        // for each sigma, compute the average posterior at the
        // corresponding bin lsl[sl] found in the previus step: that is, compute the sum of the posterior
        // read from the tensor las[lsls[sl]][i][j] divided by the number
        // of parameter_cells that contributed to that particular bin of posterior nppi[lsl[sl]]
        // ***********************************************************************************************


        fill(sf.begin(), sf.end(), 0);
        for(int sl=0;sl<6;sl++){
          for(int i=0;i<this->nbin_2D;++i){
            for(int j=0;j<this->nbin_2D;++j){
              sf[sl]+=las[lsl[sl]][i][j]/nppi[lsl[sl]];
            }
          }
        }
        /*  Output of height of normalized posterior distribution defining different confidence regions
            One line, six columns. A file for every pair of parameters*/

        ofstream kout;
        string nCR=CR+"_"+to_string(ip)+"_"+to_string(jp)+".txt";
        kout.open(nCR.c_str());
        kout.precision(12);
        kout.setf(ios::showpoint);
        kout.setf(ios::scientific);


        for(int sl=0;sl<6;sl++)kout<<sf[sl]<<"\t"; kout<<"\t"<<endl;
	this->So.message_screen("Wrote CR in file",nCR);
        kout.close();


      }
    }

  }


  std:: cout<<GREEN<<"Posterior distributions succesfully written!"<<RESET<<std::endl;



 // NOW WE GET THERE THE 1D POSTERIORS FOLLOWING THE POSTERIOR1D FUNCTION. HERE THE 1D ARRAYS WILL HAVE THE SAME DIMENTIONA STHE 2D ARRAYS
  int perc=1;


  vector<real_prec> new_acc(np,0);
  vector< vector<real_prec> > e_al_max(np, vector<real_prec> (this->nbin_2D, 0));                                            /*for ech aprameter, computes the cumultive number counts starting from the upper limit*/
  vector< vector<real_prec> > al_n(np, vector<real_prec> (this->nbin_2D, 0));

  for(int ip=0; ip<np;++ip){
    if(ip!=3){
      if(ip==0)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)al_n[ip][i]+=ALL[i][j][k][l];
      if(ip==1)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)al_n[ip][i]+=ALL[j][i][k][l];
      if(ip==2)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)al_n[ip][i]+=ALL[j][k][i][l];
      if(ip==4)for(int i=0;i<this->nbin_2D;++i)for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<this->nbin_2D;++k)for(int l=0;l<this->nbin_2D;++l)al_n[ip][i]+=ALL[j][k][l][i];
    }
  }
  real_prec pdfmax;
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]!=0){
      string nfname=fname_pdf+"_par_"+to_string(ip)+".txt";
      osali.open(nfname.c_str());
      osali.precision(12);
      osali.setf(ios::showpoint);
      osali.setf(ios::scientific);
      real_prec pd=-1000;

      for(int nb=0;nb<this->nbin_2D;++nb){
        pdfmax=max(pd,al_n[ip][nb]);
        pd=pdfmax;
      }

      for(int nb=0;nb<this->nbin_2D;nb++)osali<<this->parameters_min_2dplot[ip]+(nb+0.5)*ddel[ip]<<"   "<<al_n[ip][nb]/pd<<endl; //write normalized to the maximum
      osali.close();


      real_prec tot=0;
      for(int nb=0;nb<this->nbin_2D;nb++)tot+=al_n[ip][nb];
      for(int nb=0;nb<this->nbin_2D;nb++)al_n[ip][nb]/=tot;   // normalized to unit area

      this->So.message_screen("Wrote 1xd posterior of combined chi squared in file",nfname);
;
    }
  }

  /*Define arrays*/
  real_prec percentil[4]={0,0.1587,0.0228, 0.0013};    /*corresponds to 15.87th percentil, 2.28 percentil and 0.13th percentil*/
  vector<int> JJ(np,0);
  fill(new_acc.begin(), new_acc.end(), 0);
  for(int ip=0;ip<np;ip++)for(int nb=0;nb<this->nbin_2D;nb++)new_acc[ip]+=al_n[ip][nb];
  for(int ip=0;ip<np;ip++){                                                 /*Loop over the parameters allowed to vary in the MC chain*/
    if(this->fixed_parameters[ip]!=0){
      for(int j=this->nbin_2D-1;j>=0;j--)for(int k=this->nbin_2D-1;k>=j;k--)e_al_max[ip][(this->nbin_2D-1)-j]+=al_n[ip][k]; /*Add events starting from the upper limit towards the minimum*/
    }
  }
  for(int ip=0;ip<np;ip++){                                                 /*Loop over the parameters allowed to vary in the MC chain*/
    if(this->fixed_parameters[ip]!=0){
      for(int h=0;h<this->nbin_2D;h++)if(e_al_max[ip][h]<percentil[perc]*(new_acc[ip]) && e_al_max[ip][h+1]>=percentil[perc]*(new_acc[ip])){
          JJ[ip]=h;                                                          /*  Identify Bin*/
        }
      this->d_max[ip]=this->parameters_max_2dplot[ip]-JJ[ip]*ddel[ip];                           /*  Compute value of the parameter in that bin*/
    }
  }

  for(int i=0;i<np;++i)for(int j=0;j<this->nbin_2D;++j)e_al_max[i][j]=0;
  for(int ip=0;ip<np;ip++){                                                      /*Loop over the parameters allowed to vary in the MC chain*/
    if(this->fixed_parameters[ip]!=0){
      for(int j=0;j<this->nbin_2D;++j)for(int k=0;k<=j;k++)e_al_max[ip][j]+=al_n[ip][k];  /*Add events starttin from the upper limit towards the minimum*/
    }
  }

  fill(JJ.begin(), JJ.end(), 0);
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]!=0){
      for(int h=0;h<this->nbin_2D;h++){
        if(e_al_max[ip][h+1]>=percentil[perc]*(new_acc[ip]) && e_al_max[ip][h]<percentil[perc]*(new_acc[ip])){
          JJ[ip]=h;                                                            /*  Identify Bin*/
        }
      }
      this->d_min[ip]=this->parameters_min_2dplot[ip]+JJ[ip]*ddel[ip];                             /*  Compute value of the parameter in that bin*/
    }
  }
  // Output
  // *********************************************
  for(int ip=0;ip<np;ip++){
    if(this->fixed_parameters[ip]==0 || ip==3){  //excluyo la info sobre el b
      this->d_max[ip]=this->mean_parameters[ip];
      this->d_min[ip]=this->mean_parameters[ip];
    }
    else  for(int h=0;h<this->nbin_2D;++h)this->mean_parameters[ip]+= (this->parameters_min_2dplot[ip]+(h+0.5)*ddel[ip])*al_n[ip][h];
  }


  fname_mean=fname_mean+".txt";
  osali.open(fname_mean.c_str());
  osali.precision(12);
  cout<<BLUE<<"Writing mean, 68.2 percent confidence interval and standard deviation in file "<<RESET<<GREEN<<fname_mean<<endl;
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)osali<<ip<<"\t"<<this->mean_parameters[ip]<<"\t"<<this->d_max[ip]-this->mean_parameters[ip]<<"\t"<<this->mean_parameters[ip]-this->d_min[ip]<<endl;
  osali.close();
  std::cout<<CYAN<<"        Mean values                   +           68.2 CL      -                 Standard deviation" <<endl;
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)cout<<ip<<"\t"<<this->mean_parameters[ip]<<"\t"<<this->d_max[ip]-this->mean_parameters[ip]<<"\t"<<this->mean_parameters[ip]-this->d_min[ip]<<endl;

  cout<<"Latex Style"<<endl;
  std::cout.precision(3);
  for(int ip=0;ip<np;ip++)if(this->fixed_parameters[ip]==1)cout<<ip<<"\t"<<this->mean_parameters[ip]<<"^{+"<<this->d_max[ip]-this->mean_parameters[ip]<<"}_{-"<<this->mean_parameters[ip]-this->d_min[ip]<<"}"<<endl;






  return;
}



/************************************************************************************************************************************************************************/

void McmcFunctions::mass_asg2d(string MAS,real_prec x,real_prec y,real_prec deltax,real_prec deltay,real_prec wpart,vector< vector<real_prec> >&data){
  /*Determining the number of cell*/
  int N=this->nbin_2D;

  int xc=(int)floor(x/deltax);
  int yc=(int)floor(y/deltay);

  /*Location of the center of each cell*/
  real_prec xx  = deltax*(xc+0.5);
  real_prec yy  = deltay*(yc+0.5);
  /*Location of the center of each cell forward*/
  real_prec xxf = deltax*(xc+1.5);
  real_prec yyf = deltay*(yc+1.5);
  /*Location of the center of each cell backward*/
  real_prec xxb = deltax*(xc-0.5);
  real_prec yyb = deltay*(yc-0.5);

  real_prec Xb, Yb;   //siempre apareceran como Xb-1
  Xb=(xc==0 ? this->nbin_2D: xc);//es decir si xc=0, Xf-1=n1-1, i.e, salta hacia atras y termina el la ultima
  Yb=(yc==0 ? this->nbin_2D: yc);

  real_prec Xf, Yf;  //siempre apareceran como Xf+1
  Xf=(xc==this->nbin_2D-1 ? -1: xc);//es decir si xc=n1-1, Xf+1=0, i.e, salta hacia adelante y termina el la primera
  Yf=(yc==this->nbin_2D-1 ? -1: yc);


  /*Filling cells*/
  data[Xb-1][Yb-1]+=wpart*maf(MAS,(xxb-x)/deltax)*maf(MAS,(yyb-y)/deltay);
  data[Xb-1][yc]  +=wpart*maf(MAS,(xxb-x)/deltax)*maf(MAS,(yy-y)/deltay);
  data[Xb-1][Yf+1]+=wpart*maf(MAS,(xxb-x)/deltax)*maf(MAS,(yyf-y)/deltay);

  data[xc][Yb-1]+=wpart*maf(MAS,(xx-x)/deltax)*maf(MAS,(yyb-y)/deltay);
  data[xc][yc]  +=wpart*maf(MAS,(xx-x)/deltax)*maf(MAS,(yy-y)/deltay);
  data[xc][Yf+1]+=wpart*maf(MAS,(xx-x)/deltax)*maf(MAS,(yyf-y)/deltay);

  data[Xf+1][Yb-1]+=wpart*maf(MAS,(xxf-x)/deltax)*maf(MAS,(yyb-y)/deltay);
  data[Xf+1][yc]  +=wpart*maf(MAS,(xxf-x)/deltax)*maf(MAS,(yy-y)/deltay);
  data[Xf+1][Yf+1]+=wpart*maf(MAS,(xxf-x)/deltax)*maf(MAS,(yyf-y)/deltay);


}

// ************************************************************************** 
// ************************************************************************* *

/*
  LSS Distance priors obtained from  different LSS probe: CMB power spectum, BAOS
  returns the inverse of the covariance matrix assocaited to a set of distance priors.
  That inv_cov will be used to constrain cosmological parameters


  As input, I and J denotes the distance and WMAP prios to be used as follows

  ---
  I=1: WMAP5 5 (Komatsu et al.) with the following sub-cases:
    J=1: l_A
    J=2: R
    J=3: z_dec
    J=4: r_s;
    J=5: l_A+R
    J=6: l_A+z
    J=7: R+z
    J=8: l_A,R,z
  
  --
  I=2, J=1: l_A+R+z+ob ,extended distance priors, WMAP5, Komatsu et al. 2009
  
  --
  I=3, J=1: l_A,R,z, from WMAP7(Komatsu et al 2010)
  
  --
  I=4, J=1: CMB+LSS: 100om_b, z*, l_A(z_star),R,G  extended distance priors (Sanchez et al 2010): Ariel sent a warning: these are wrong.


  WARNING: INCLUDE ALWAS AFTER THE INCLUSON OF THE FUNCTION subroutine matrix_inversion
 */




void McmcFunctions::set_distance_priors(int I, int J){


  if(I==1 && (J==1 || J==2 || J==3 || J==4))this->n_priors=1;
  if(I==1 && (J==5 || J==6 || J==7))this->n_priors=2;
  if(I==1 && J==8) this->n_priors=3;
  if(I==2) this->n_priors=4;
  if(I==3) this->n_priors=3;
  if(I==4) this->n_priors=5;



  this->priors_measured.resize(this->n_priors,0);
  this->inv_cov.resize(this->n_priors);

  for(int i=0;i<this->n_priors;++i)this->inv_cov[i].resize(this->n_priors,0);

  vector<vector<real_prec> > cov(n_priors,vector<real_prec>(n_priors,0));
  
  switch(I){
  case 1:
    switch(J){
    case 1:      /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5,	1) l_A*/
      this->priors_measured[0]=302.10;
      this->inv_cov[0][0]=pow(1./0.86,2);
      break;
      
    case 2: /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5	1) R = shift parameter  */ 
      this->priors_measured[0]=1.710;
      this->inv_cov[0][0]=pow(1./0.019,2);
      break;
      
    case 3:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5    1) z_*/
      this->priors_measured[0]=1090.04;
      this->inv_cov[0][0]=pow(1./0.93,2);
      break;
      
    case 4: /*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5*   
	      1) r_s = sond horizon*/
      this->priors_measured[0]=146.8;
      this->inv_cov[0][0]=pow(1./1.8,2);
      break;
      
    case 5:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_A = shift parameter
	     2) R  = shift parameter
	   */
      this->priors_measured[0]=302.10;
      this->priors_measured[1]=1.710;
      for(int i=0;i<n_priors;++i)for(int j=0;j<n_priors;++j)cov[i][j]=0;
      cov[0][0]=pow(0.86,2);
      cov[1][1]=pow(0.019,2);
      cov[0][1]=0.1109*(0.86*0.019);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,this->inv_cov);
      break;
      
    case 6:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_A = sound horizon angular scale
	     2) z  = redshift at decoupling
	   */
      this->priors_measured[0]=302.10;
      this->priors_measured[1]=1090.04;
      for(int i=0;i<n_priors;++i)for(int j=0;j<n_priors;++j)cov[i][j]=0;
      cov[0][0]=pow(0.86,2);
      cov[1][1]=pow(0.93,2);
      cov[0][1]=0.4215*(0.86*0.93);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,this->inv_cov);
      break; 

    case 7:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) R_ = shift parameter
	     2) z  = redshift at decoupling
	   */

      this->priors_measured[0]=1.710;
      this->priors_measured[1]=1090.04;
      
      for(int i=0;i<n_priors;++i)for(int j=0;j<n_priors;++j)cov[i][j]=0;
      cov[0][0]=pow(0.019,2);
      cov[1][1]=pow(0.93,2);
      cov[0][1]=0.6928*(0.019*0.93);
      cov[1][0]=cov[0][1];
      matrix_inversion(cov,this->inv_cov);
      break;

    case 8:/*DISTANCE PRIORS FROM KOMATSU et al 2009 WMAP5
	     1) l_ = sound horizon angular scale
	     2) R_ = shift parameter
	     3 z  = redshift at decoupling
	   */
      this->priors_measured[0]= 302.10;
      this->priors_measured[1]= 1.710;
      this->priors_measured[2]= 1090.04;
      this->inv_cov[0][0]= 1.800;
      this->inv_cov[0][1]= 27.968;
      this->inv_cov[0][2]=-1.103;
      this->inv_cov[1][0]= inv_cov[0][1];
      this->inv_cov[1][1]= 5667.577;
      this->inv_cov[1][2]=-92.263;
      this->inv_cov[2][0]= this->inv_cov[0][2];
      this->inv_cov[2][1]= this->inv_cov[1][2];
      this->inv_cov[2][2]= 2.923;
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
    this->priors_measured[0]=302.10;
    this->priors_measured[1]=1.710;
    this->priors_measured[2]=1090.04;
    this->priors_measured[3]=2.2765;
    
    this->inv_cov[0][0]= 31.001;
    this->inv_cov[0][1]=-5015.642;
    this->inv_cov[0][2]= 183.903;
    this->inv_cov[0][3]= 2337.977;
    
    this->inv_cov[1][0]=this->inv_cov[0][1];
    this->inv_cov[1][1]= 876807.166;
    this->inv_cov[1][2]=-32046.750;
    this->inv_cov[1][3]=-403818.837;
    
    
    this->inv_cov[2][0]=this->inv_cov[0][2];
    this->inv_cov[2][1]=this->inv_cov[1][2];
    this->inv_cov[2][2]=1175.054;
    this->inv_cov[2][3]=14812.579;
    
    this->inv_cov[3][0]=this->inv_cov[0][3];
    this->inv_cov[3][1]=this->inv_cov[1][3];
    this->inv_cov[3][2]=this->inv_cov[3][4];
    this->inv_cov[3][3]=187191.186;
    
    break;    
    
  case 3:
    /*DISTANCE PRIORS FROM KOMATSU et al 2010 WMAP7
      1) l_A(z_star),  
      2) R(z_star), 
      3) z_star = redshift  at decoupling */
    this->priors_measured[0] = 302.09;
    this->priors_measured[1] = 1.725;
    this->priors_measured[2] = 1091.3;
    this->inv_cov[0][0]= 2.305;
    this->inv_cov[0][1]= 29.698;
    this->inv_cov[0][2]= -1.333;
    this->inv_cov[1][0]= this->inv_cov[0][1];
    this->inv_cov[1][1]= 6825.270;
    this->inv_cov[1][2]=-113.180;
    this->inv_cov[2][0]= this->inv_cov[0][2];
    this->inv_cov[2][1]= this->inv_cov[1][2];
    this->inv_cov[2][2]= 3.414;
    
    break;
    
    
  case 4:
    /*DISTANCE PRIORS FROM SDSS(LSS)+CMB+SNia FROM SANCHEZ et al. 2009
      1) 100 \Omega_bar h^2, 
      2) z_star = redshift  of decoupling 
      3) l_A(z_star),  
      4) R(z_star), 
      5) G(z_m), z_m=0.35*/
    this->priors_measured[0]=2.28;
    this->priors_measured[1]=1090.12;
    this->priors_measured[2]=301.58;
    this->priors_measured[3]=1.701;
    this->priors_measured[4]=1175.;
    for(int i=1;i<n_priors;++i)for(int j=1;j<n_priors;++j)cov[i][j]=0;
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
    matrix_inversion(cov,this->inv_cov);
  }
  
  return ;
}


// ************************************************************************** 
// ************************************************************************* *

real_prec McmcFunctions::chi_squared_distance_priors(){
  real_prec ans=0;


  // for(int i=0;i<this->priors_model.size();++i){
  //   cout<<i<<"  "<<this->priors_measured[i]<<"  "<<this->priors_model[i]<<endl;
  // }

  for(int i=0;i<this->priors_model.size();++i){
    for(int k=0;k<this->priors_model.size();k++){
      ans+=(this->priors_measured[i]-this->priors_model[i])*inv_cov[i][k]*(this->priors_measured[k]-this->priors_model[k]);
    }
  }
  return ans;
}

// ************************************************************************** 

// ************************************************************************* *

void McmcFunctions::distance_priors_cmb_model(int I, int J, s_CosmologicalParameters *scp){
  Cosmology CF(*scp);
  this->priors_model.resize(this->n_priors,0);
  /*Distance priors from Komatsu et al. 2010 (CMB)*/
  /*Calculamos z_dec segun Hu & Sujiyama, ver Komatsu et al. 2007 WMAP5 eq 65 */
  real_prec omega_b  = (scp->Om_baryons)*pow(scp->hubble,2);
  real_prec omega_m  = (scp->Om_matter)*pow(scp->hubble,2);

  // Baryon drag epoch:
  real_prec z_drag = CF.drag_redshift(); 
  real_prec da_zd  = CF.proper_angular_diameter_distance(z_drag);
  real_prec rsound = CF.comoving_sound_horizon(z_drag);
  real_prec R      = (1+z_drag)*da_zd*sqrt(scp->Om_matter)*scp->Hubble/Constants::speed_light;
  real_prec lA     = (1+z_drag)*da_zd*M_PI/rsound;

  if(I==1 && J==1)this->priors_model[0]=lA;
  if(I==1 && J==2)this->priors_model[0]=R;
  if(I==1 && J==3)this->priors_model[0]=z_drag;
  if(I==1 && J==4)this->priors_model[0]=rsound;
  if(I==1 && J==5){
    this->priors_model[0] =  lA;
    this->priors_model[1] =  R;
  }
  if(I==1 && J==6){
    this->priors_model[0] =  lA;
    this->priors_model[1] =  z_drag;
  }
  if(I==1 && J==7){
    this->priors_model[0] =  R;
    this->priors_model[1] =  z_drag;
  }
  if((I==1 && J==8) || I==3){
    this->priors_model[0] =  lA;
    this->priors_model[1] =  R;
    this->priors_model[2] =  z_drag;
  }
  if(I==2 && J==1){
    this->priors_model[0] =  lA;
    this->priors_model[1] =  R;
    this->priors_model[2] =  z_drag;
    this->priors_model[3] =  100*omega_b;
  }
  if(I==4){
    this->priors_model[0] =  100*omega_b;
    this->priors_model[1] =  z_drag;
    this->priors_model[2] =  lA;
    this->priors_model[3] =  R;
    this->priors_model[4] =  (1+Constants::z_mean_sdss)*this->Cf.proper_angular_diameter_distance(Constants::z_mean_sdss)*pow(0.01*this->Cf.Hubble_function(Constants::z_mean_sdss),0.8);
  }

  

  cout<<"Distance priors with fid. cosmology      Measured priors              difference"<<endl;
  for(int i=0;i<n_priors;i++)cout<<priors_measured[i]<<"                       "<<priors_model[i]<<"     "<<100*(priors_measured[i]-priors_model[i])/priors_measured[i]<<"%"<<endl;
  real_prec chis_cmb=chi_squared_distance_priors();
  
  for(int i=0;i<priors_model.size();++i)priors_model[i]=0;  
  cout<<"  "<<endl;
  cout <<"========================================================================================== "<<endl; 
  cout<<"chi^2 from fiducial parameters   "<<chis_cmb<<endl;
  cout <<"========================================================================================== "<<endl;
  

}

