////////////////////////////////////////////////////////////////////////////
/**
 * @class<HOD>
 * @brief Header file for the class HOD::
 * @file HOD.h
 * @title Functions related to HOD analysis
 * @author ABA
 */
////////////////////////////////////////////////////////////////////////////
#ifndef _HOD__
#define _HOD__
# include "NumericalMethods.h"
# include "Params.h"
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @brief The s_hods_par struct
 * @details Structure containing parameers of simple model of HOD
 */
struct s_hods_par{
  int hod_model;
  real_prec mmin;
  real_prec alpha;
  real_prec scatter;
  real_prec muno;
};
////////////////////////////////////////////////////////////////////////////

class HOD{
 private:
    s_CosmologicalParameters s_cosmo_pars;
    Params params;
 public:
    ////////////////////////////////////////////////////////////////////////////
/**
* @brief Defauilt constructor
*/
  HOD(){}
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Default constructor passing an intance of the class Params
  */
  HOD(Params par){
      this->params=par;
  };
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Default constructor passing an intance of the structure cosmolo_pars
  */
  HOD(s_CosmologicalParameters par){
      this->s_cosmo_pars=par;
  };
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Default destructor
  */
  ~HOD(){}
  ////////////////////////////////////////////////////////////////////////////
/**
  * @brief Mean of the halo occupation distribution for central galaxies as a function of the host halo mass
  * @arg M : Halo mass in Ms/h
*/
real_prec CENTRAL(real_prec M);
////////////////////////////////////////////////////////////////////////////
/**
* @brief Mean of the halo occupation distribution for satellite galaxies as a function of the host halo mass
  * @arg M : Halo mass in Ms/h
*/
  real_prec SATELLITE(real_prec M);
//##################################################################################
/**
* @brief Minimum (cut) mas son HOD as a function of halo mass, derived from the Halo-Mass stellar mass relation
* @arg M : Halo mass in Ms/h
*/
  real_prec Mmin_hod(real_prec M);
////////////////////////////////////////////////////////////////////////////
/**
* @brief Width of the normal distribution upon which the HOD for centrals is computed
* @arg M : Halo mass in Ms/h
*/
  real_prec sigma_hod(real_prec M);
////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Conversion from Halo mass to Stellar mass for central galaxies (Median stellar mass)
  * @arg M : Halo mass in Ms/h
  */
  real_prec Mhalo_to_Mstellar(real_prec Mh);
////////////////////////////////////////////////////////////////////////////
  /**
  * @brief Conversion from Halo mass to Stellar mass for central galaxies (Median stellar mass)
  * @arg M : Halo mass in Ms/h
  */
  real_prec alpha_par();
  ////////////////////////////////////////////////////////////////////////////
  /**
  * @brief
  */
  void set_cosmo_pars(s_CosmologicalParameters s_pars){this->s_cosmo_pars=s_pars;}
////////////////////////////////////////////////////////////////////////////
  /**
  * @brief
  */
  void set_pars(Params s_pars){this->params=s_pars;}
////////////////////////////////////////////////////////////////////////////
};
#endif
