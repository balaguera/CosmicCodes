////////////////////////////////////////////////////////////////////////////
/**
 * @class<McmcFunctions>
 * @brief Header file for the class McmcFunctions
 * @file McmcFunctions.cpp
 * @title Methods related to McmcFunctions
 * @author   ABA
 */
 ////////////////////////////////////////////////////////////////////////////

#ifndef _McmcFunctions_
#define _McmcFunctions_

# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <sstream>
# include <cassert>
# include <vector>
# include <limits>
# include <algorithm>
# include "CosmologicalFunctions.hpp"
# include "Miscelanious.hpp"
# include "FileOutput.hpp"
# include "ScreenOutput.hpp"

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

using namespace std;
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

class McmcFunctions{

private:

  //------------------------------------------
  /** 
   * @brief 
  */
  Params params;
  //------------------------------------------
  /** 
   * @brief 
  */
  Cosmology Cf;
  //------------------------------------------
  /** 
   * @brief 
  */
  vector<real_prec>priors_measured;
  //------------------------------------------
  /** 
   * @brief 
  */
  vector<real_prec>priors_model;
  //------------------------------------------
  /** 
   * @brief 
  */
  vector<vector<real_prec> >inv_cov;
  //------------------------------------------
  /** 
   * @brief 
  */
  ofstream osali;
  //------------------------------------------
  /** 
   * @brief 
  */
  ScreenOutput So;

//------------------------------------------
  /** 
   * @brief  Number of models accepted by the HM algorithm
   */
  string action;


 public:

  McmcFunctions() {}

  McmcFunctions(Params &_pars): params(_pars) {}

  ~McmcFunctions()
  {

    if(gBaseRand)
        gsl_rng_free(gBaseRand);
  }
  //////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////
  //------------------------------------------
  /** 
   * @brief 
  */
  const gsl_rng_type * T;
  //------------------------------------------
  /** 
   * @brief 
  */
  gsl_rng * gBaseRand;
  //------------------------------------------
  /** 
   * @brief 
  */
  int seed;
  //------------------------------------------
  /** 
   * @brief 
  */
  ofstream lpout;

  //------------------------------------------
  /** 
   * @brief 
  */
  vector<double> tail_probabilities;

 //------------------------------------------
  /** 
   * @brief 
  */
  double step_size_multiplier;

  /**
   * @brief   read it from from params
   */
  int seed_index;

  //////////////////////////////////////////////////////////
  /**
   * @brief   read it from from params
   */
  int n_parameters;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  int n_priors;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  int chain_ini;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  int chain_fin;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  int chain;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  bool use_distance_priors;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string sampling;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  bool  diagonal_covariance_matrix;
  //////////////////////////////////////////////////////////
  /**
   * @brief  read it from from params
   */
  bool  analytic_marginalization_wrt_amplitude;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  int Id;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  int Jd;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  int nbin_1D;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  int nbin_2D;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string MAS;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string experiment;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string model;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string observable;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string fit_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string input_par_file_model;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  string file_oo;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector < vector<real_prec> > acc_parameters;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> weight;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> chiss;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> prior_parameters_max_values;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> prior_parameters_max_values_2dplot;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> prior_parameters_min_values;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> prior_parameters_min_values_2dplot;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> proposal_parameters;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> parameters;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec> initial_parameters;
  //////////////////////////////////////////////////////////
  /** 
   * @brief 
  */
  std::vector<double>mean_parameter;
 //------------------------------------------
  /** 
   * @brief 
  */
  std::vector<double>stdev_parameter;
 //------------------------------------------
  /** 
   * @brief 
  */
  std::vector<double>upper_bound_sigma;
 //------------------------------------------
  /** 
   * @brief 
  */
  std::vector<double>lower_bound_sigma;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<real_prec>    R;
  //////////////////////////////////////////////////////////
  //Container with all steps of merged chains
  /**
   * @brief 
   */
  vector<vector <vector <real_prec> > > acc_parameters_all;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<int> mark;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<int> acca;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<vector<real_prec> > weight_r;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector<vector<real_prec> > chi_models;
  //////////////////////////////////////////////////////////
  //Containers for HMC
  /**
   * @brief 
   */
  vector<real_prec> grad_U;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector< vector<real_prec> >dCmod_dp;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <real_prec> momentum;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <real_prec> masses;
    //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <real_prec> epsilon;
    //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <real_prec> mean_parameters;
    //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <real_prec> stdev_parameters;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <vector<real_prec> > cov_masses;
  //////////////////////////////////////////////////////////
  /** 
   * @brief 
  */
  std::vector<std::vector<double>>cov_parameters_aux;
  //////////////////////////////////////////////////////////
  /** 
   * @brief 
  */
  std::vector<std::vector<double>>cov_parameters_decomposed;  
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  vector <vector<real_prec> > icov_masses;
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void set_mcmc_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void set_mcmc_read_vectors();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void gelman_rubbin_diag( vector<vector<vector<double>>> &acca_parameters);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void likelihood_full(int lmin, int lmax, vector<matrices>&VM,real_prec &chis_one);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_det_matrix(vector<vector<real_prec> >&mat, real_prec &determinant);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared(vector<real_prec>&,vector<real_prec>&, vector< vector <real_prec> > &, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared(const vector<real_prec>&,const vector<real_prec>&, vector<real_prec> &, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared_poisson(vector<real_prec>&,vector<real_prec>&, real_prec &); //change names to loglikelihood
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_loglikelihood_Poisson(vector<vector<real_prec> >&,vector<vector<real_prec>>&, real_prec &);

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared(vector<real_prec>&, vector<real_prec> &,vector<real_prec> &, vector<real_prec>&, real_prec, real_prec, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared_marginalized_amplitude(vector<real_prec> &, vector<real_prec> &, vector< vector<real_prec> > &, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void chi_squared_marginalized_amplitude(vector<real_prec> &, vector<real_prec> &, vector<real_prec> &, vector<real_prec>&, real_prec, real_prec, real_prec &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_cova_pars(int acc);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void jump(size_t nacc);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void MHalgorithm(real_prec&, real_prec&, size_t &, int&);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_accepted_models(int j, size_t acc, int weight_here, bool screen);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec kinetic(int);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_gradU(vector<real_prec>&model, vector<real_prec>&meas, vector<vector<real_prec> >&Step,vector<vector<real_prec> >&R,vector<vector<real_prec> >&V,vector<vector<real_prec> >&icov);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_gradU_full(int lmin, int lmax, vector<matrices> VM,vector<real_prec>&Cmodel,vector<real_prec>&gradU_full);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_model_FB(vector<vector<real_prec> >&Step,vector<vector<real_prec> >&R,vector<vector<real_prec> >&V,vector<real_prec>&Cm);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_cova_FB_delta(vector<matrices>&,int, vector<vector<real_prec> >&Cmodel);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void auto_correlation_mcmc(vector< vector<real_prec> > &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void entropy_mcmc(vector< vector<real_prec> > &,vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void posterior1d(string, string, vector<real_prec>  &, vector<vector<real_prec> > &,int);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void posterior2d(string,string, vector<real_prec> &, vector<vector<real_prec> > &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_contour_levels(string,vector<vector<real_prec> > &);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void get_contour_levels(string,ULONG, vector<real_prec> &);

  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @returns
   */
  void posterior2d_combined_experiments(string, string, string,string, experiments ex);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @returns
   */
  void posterior1d_combined_experiments(string,string, experiments ex);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @returns
   */
  void get_stats(int ip, double &median, double &std);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @returns
   */
  void get_stats();
  //////////////////////////////////////////////////////////
  /**
   * @brief
   * @returns
   */
  void get_covariance_parameters();

  //////////////////////////////////////////////////////////


  /**
   * @brief 
   * @details LSS Distance priors obtained from  different LSS probe: CMB power spectum, BAOS
  *  returns the inverse of the covariance matrix assocaited to a set of distance priors.
  *  That inv_cov will be used to constrain cosmological parameters
  *  As input, I and J denotes the distance and WMAP prios to be used as follows

  *   I=1: WMAP5 5 (Komatsu et al.) with the following sub-cases:
  *   J=1: l_A
  *   J=2: R
  *   J=3: z_dec
  *   J=4: r_s;
  *   J=5: l_A+R
  *   J=6: l_A+z
  *   J=7: R+z
  *   J=8: l_A,R,z
  *   I=2, J=1: l_A+R+z+ob ,extended distance priors, WMAP5, Komatsu et al. 2009
  *   I=3, J=1: l_A,R,z, from WMAP7(Komatsu et al 2010)
  *   I=4, J=1: CMB+LSS: 100om_b, z*, l_A(z_star),R,G  extended distance priors (Sanchez et al 2010).
  * @warning Include after subroutine matrix_inversion
  * @warning I=4, J=1 case (Sanchez et al 2010) are wrong.
*/
  void set_distance_priors();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void distance_priors_cmb_model(int I, int J, s_CosmologicalParameters *scp);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec chi_squared_distance_priors();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void run_distance_priors(string, int);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void analyze_chains(int ip);

  //////////////////////////////////////////////////////////
 /** 
   * @brief Retrieve private variable
  */
  string _action (){return this->action;}
  //////////////////////////////////////////////////////////
 /** 
   * @brief Retrieve private variable
  */
  void set_action(string newa){this->action=newa;}
  //////////////////////////////////////////////////////////
 /** 
   * @brief Set value of private variable
  */
  void set_seed_index(int newt){this->seed_index = newt;}

};






#endif
