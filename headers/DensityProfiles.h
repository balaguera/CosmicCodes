////////////////////////////////////////////////////////////////////////////
/**
* @class
* @file DensityProfiles.h
* @brief Header file for the class DensityProfiles::
* @title Functions related to the generation of density profiles for dark matter haloes
* @details Bias Assignment method for mock catalogs
* @author Andres Balaguera-Antolínez
*/
////////////////////////////////////////////////////////////////////////////

#ifndef _DensityProfiles_
#define _DensityProfiles_
# include "NumericalMethods.h"
# include "PowerSpectrumTH.h"
# include "CosmologicalFunctions.h"
# include "Statistics.h"


using namespace std;


class DensityProfiles{

private:
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief 
    */
  real_prec rs;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief 
    */
  real_prec alpha;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  real_prec rvir;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  real_prec rhos;
  //////////////////////////////////////////////////////////
     /**
  * @private
    * @brief
    */
  real_prec concentration;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  real_prec fc;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  s_CosmologicalParameters s_cosmo_pars;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  static gsl_real idensity_rc(gsl_real, void *);
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  static gsl_real idensity_fourier_c(gsl_real, void *);
  //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief
    */
  static gsl_real idensity_fourier(gsl_real, void *);


 public:
  //////////////////////////////////////////////////////////
   /**
    * @brief Default constructor
    */
  DensityProfiles(){}
  //////////////////////////////////////////////////////////
   /**
    * @brief Overlodaded constructor
    */
  DensityProfiles(s_CosmologicalParameters _spars):s_cosmo_pars(_spars){}
  //////////////////////////////////////////////////////////
   /**
    * @brief Default destructor
    */
  ~DensityProfiles(){}
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec density_r(real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief NFW Density profile
    */
  real_prec DensityProfile_NFW(real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec DensityProfile_NFW_prob(real_prec ,real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief Density profile in Fourier space
    */
  real_prec density_fourier(real_prec,real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec mass_concentration_dis(real_prec,real_prec,real_prec, void*);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec density_rc(real_prec,real_prec, real_prec,  void *);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec density_k(real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec density_kc(real_prec,real_prec, real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  void einasto_parameters(real_prec, real_prec, real_prec *, real_prec *, real_prec *,real_prec *,real_prec *);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  void einasto_parameters(real_prec);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  void nfw_parameters(real_prec , real_prec, real_prec *, real_prec *,real_prec *,real_prec *,real_prec *);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  void nfw_parameters(real_prec mass);
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec _rvir(){return this->rvir;}
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec _concentration(){return this->concentration;}
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec _rhos(){return this->rhos;}
  //////////////////////////////////////////////////////////
   /**
    * @brief 
    */
  real_prec _rs(){return this->rs;}

};





#endif
