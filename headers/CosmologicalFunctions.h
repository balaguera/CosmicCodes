//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
 * @class<Cosmology>
 * @brief Header file for the class Cosmology::
 * @file CosmologicalFunctions.h
 * @title Functions to compute cosmological dependent quantities
 * @author ABA
 * @callgraph
 * @date  2012-2023
 */
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#ifndef __COSMOLOGY__
#define __COSMOLOGY__
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include "Params.h"
# include "Constants.h"
using namespace std;
using namespace Constants;
//////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 * @struct<s_CosmoInfo>
 * @brief The s_CosmoInfo struct
 * @details Auxiliary structure to allocate cosmological quantities derived from cosmological parameter in class::Cosmology
 * @note This strucure will be merged with s_cosmo_pars.
 */
struct s_CosmoInfo
{
  /**
   *@brief Critical energy density of the Universe
   */
  real_prec critical_density;
  /**
   *@brief Density contrast for spherical collapse */
  real_prec density_contrast_top_hat;
  /**
   *@brief Hubble parameter at input cosmological redshift */
  real_prec Hubble_parameter;
  /**
   *@brief Comoving distance at input cosmological redshift */
  real_prec comoving_distance;
  /**
   *@brief Comoving angular diameter distance */
  real_prec comoving_angular_diameter_distance;
  /**
   *@brief Mean matter density of the Universe at input cosmological redshift */
  real_prec mean_matter_density;
  /**
   *@brief Age of the Universe at input cosmological redshift */
  real_prec age_universe;
  /**
   *@brief Comoving sound horizon at input cosmological redshift*/
  real_prec comoving_sound_horizon;
  /**
   *@brief Growth factor (valid for LCDM) at input cosmological redshift*/
  real_prec growth_factor;
  /**
   *@brief Auxiliary growth used by 2LPT */
  real_prec D2;
  /**
   *@brief Growth index at input cosmological redshift*/
  real_prec growth_index;
  /**
   *@brief Growth index at input cosmological redshift*/
  real_prec growth_index2;
  /**
   *@brief Halo dynamical time*/
  real_prec halo_dynamical_time;
  /**
   *@brief Matter energy density parameter at input cosmological redshift*/
  real_prec omega_matter;
  /**
   *@brief Distance modulos at input cosmological redshift*/
  real_prec Distance_Modulus;
  /**
   *@brief Cosmological scale factor at input cosmological redshift*/
  real_prec scale_factor;
  /**
   *@brief Density contrast Dv/Omega_m */
  real_prec Delta_Vir;
};


//////////////////////////////////////////////////////////
/**
 * @brief This demands 1% convergence precision for the Newton-Rhapson algorithm
 */
#define ERROR_LIMIT_NR static_cast<real_prec>(0.01)
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand for growth
 */
gsl_real gint(gsl_real, void *);
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand for comoving distance
 */
gsl_real i_rs(gsl_real, void *);
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand for comoving and age
 */
gsl_real Hinv(gsl_real, void *);
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand
 */
gsl_real Hinva(gsl_real, void *);
//////////////////////////////////////////////////////////
/**
 * @static
 * @brief Integrand
 */
gsl_real Froot(gsl_real , void*);
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand
 */
gsl_real dFroot(gsl_real , void *);
//////////////////////////////////////////////////////////
/**
 * @static functions
 * @brief Integrand
 */
void F_dF(gsl_real , void *, gsl_real *, gsl_real *);
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
class Cosmology{
  
 private:

   //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief Object of type CosmologicalFunctions
    */
   real_prec scale_factor(real_prec);
   //////////////////////////////////////////////////////////
   /**
    * @private
    * @brief Object of type CosmologicalFunctions
    */
   real_prec redshift;
    //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Number of points to gsl integration
   */
  ULONG NP;
  //////////////////////////////////////////////////////////
  /** 
   * @private
   * @brief Object of type CosmologicalFunctions
   */
  gsl_integration_glfixed_table *wf;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  s_CosmologicalParameters s_cosmo_pars;

  //////////////////////////////////////////////////////////
 public:
  
  /**
   * @brief Default constructor
   */
  Cosmology(){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief
   */
  Cosmology(s_CosmologicalParameters spar, real_prec red):s_cosmo_pars(spar), redshift(red){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Default constructor
   * @param Structure of type s_CosmologicalParameters
   */
  Cosmology(s_CosmologicalParameters spar):s_cosmo_pars(spar){}  // Default constructor
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Default destructor
   */
  ~Cosmology(){}
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Constructor
   * @param Np Number of points for the Gauss-Legendre integration of cosmological functions.
   * @return wf gsl_table for integration
   */
  Cosmology(int np){
    NP=np;
    wf=gsl_integration_glfixed_table_alloc(NP);
  } 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Hubble function at the input redshift, \f$ H(z)=H_{0}\sqrt{(1+z)^{3}\Omega_{m}+(1+z)^{2}\Omega_{rad}+(1+z)^{2}\Omega_{K}+(1+z)^{3(w+1)}\Omega_{de}} \f$
   */
  real_prec Hubble_function(real_prec z);
  real_prec Hubble_function();
    //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Comoving distance at the input redshift \f$ \eta(z)=\int_{0}^{z}\frac{d z}{H(z)}. \f$
   */
  real_prec comoving_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Proper_angular_diameter_distance at the input redshift
   */
  real_prec proper_angular_diameter_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Comoving_angular_diameter_distance at the input redshift
   */
  real_prec comoving_angular_diameter_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Proper_angular_diameter_distance to be interpolated
   */
  real_prec inter_proper_angular_diameter_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Luminosity Distance
   */
  real_prec luminosity_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Luminosity Distance  to be interpolated
   */
  real_prec inter_luminosity_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Baryon drag redshift
  */
  real_prec drag_redshift();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief R_drag
  */
  real_prec R_drag();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief R_at equalyu
  */
  real_prec R_equality();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief R_drag
  */
  real_prec sound_horizon();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief redshift at equality matter radiation
  */
  real_prec redshift_equality();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief redshift at equality matter radiation
  */
  real_prec k_equality();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Silk diffusion lenght
  */
  real_prec k_Silk();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @returns transverse_comoving_distance \f$ d_{T}(z)\equiv c
    \begin{cases}
    \mathrm{asin} \left( \int_{0}^{z}\frac{d z'}{H(z')} \right)& k=+1 \\
    \int_{0}^{z}\frac{d z'}{H(z')} &k=0\\
    \mathrm{asinh} \left( \int_{0}^{z}\frac{d z'}{H(z')} \right) & k=-1. \\
    \end{cases}\f$
   */
  real_prec transverse_comoving_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec inter_transverse_comoving_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec derivative_transverse_comoving_distance(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @details Mean background density as a function of redshift in units (Ms/h)/(Mpc/h)^(-3)
   * @param z Cosmological Redshift
   * @returns
   */
  real_prec mean_matter_density(real_prec z);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Spherical Overdensity of halos computed from Mvir and Rvir
   * @input r Halo radii in Mpc/h
   * @input M Halo mass in M/h
   * @input z Cosmological redshift
   * @param z Cosmological Redshift
   * @output SO Spherical oVerdensity
   */
  real_prec SO(real_prec, real_prec, real_prec);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Age of the Universe
   * @param z Cosmological Redshift
   * @returns Age of the Universe \f$T(a)=T_{H}\int_{0}^{a} \frac{1}{a'E(a')} d a' \f$
   */
  real_prec age_universe(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Comoving sound horizon
   * @param z Cosmological Redshift
   * @returns Comoving sound horizon in Mpc/h, \f$s(z)=\int_{z}^{\infty}\frac{c_{s}(z)}{H(z)}d z \f$
   */
  real_prec comoving_sound_horizon(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Growth factor
   * @param z Cosmological Redshift
   * @return Growth factor computed as \f$ D(a)=\frac{5}{2}\Omega_{\mathrm{mat}}E(a)\int_{0}^{a}[a' E(a')]^{-3} da'\f$
   */
  real_prec growth_factor(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Growth Index
   * @param z Cosmological Redshift
   */
  real_prec growth_index(real_prec );
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Secondary growth index
   * @param z Cosmological Redshift
   */
  real_prec growth_index2(real_prec);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Halo dynamical time
   * @param z Cosmological Redshift
   * @returns Halo dynamical time \f$\tau_{h}=\frac{0.1}{H(z)} \times \frac{Mp2km}{yrs2sec}\f$
   */
  real_prec halo_dynamical_time(real_prec);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec omega_matter(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec omega_radiation(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec omega_curvature(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec omega_dark_energy(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @returns Distance modulus \f$ DM(z)=25+5\log_{10} d_{L}(z)\f$
*/
  real_prec Distance_Modulus(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec inter_Distance_Modulus(real_prec z);
    //////////////////////////////////////////////////////////
  /** 
   * @public
   */
  real_prec zmax();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * details Given the limiting magnitude of the sample and the absolute magnitude of the object, find the zmax to which this object could have been observed.
   */
  real_prec zmax_old();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @details Given the limiting magnitude of the sample and the absolute
     magnitude of the object, find the zmax to which this object could have been
     observed.
   * @param z Cosmological Redshift
   */
  real_prec zmax_old(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @details  This is r/a, i.e, \f$rr(z)\f is a comoving distance, only used when the top-hat window function is needed. In the end is redshift independent because mean density sales as (1+z)³
   * @returns Comoving distance in Mpc/h
   */
  real_prec rr(real_prec , real_prec);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @param z Cosmological Redshift
   * @return Comoving distance at input redshift
   */
  real_prec rr_lag(real_prec , real_prec);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Top-hot Density contrast
   * @details  Density contrast from Bryan and Normanm 98, assuming that the redhsift observed is the same as the redshift of collapse. Notice that we divide over omega_matter: this is because the B&N equation is meant for critical density, not mean density. Here we compute it against the mean density.
   */
  real_prec density_contrast_top_hat(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return Critical overdensity linearly extrapolated
   */
  real_prec critical_overdensity(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return K-correction at input redshift
   */
  real_prec K_correction(real_prec z);
  //////////////////////////////////////////////////////////
  /**
  /**
   * @public
   * @param z Cosmological Redshift
   */
  real_prec K_correction(real_prec,real_prec,real_prec);

  //////////////////////////////////////////////////////////
  /** 
   * @brief Object of type CosmologicalFunctions 
   */
  real_prec dK_correction_dz(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   * @return e-correction at input redshift
   */
  real_prec e_correction(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @param z Cosmological Redshift
   */
  real_prec de_correction_dz(real_prec z);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  void free_gsl_table();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  void check_cosmo_pars();
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  real_prec gsl_cosmo_integration(gsl_real (*function)(gsl_real, void *) ,void *p, gsl_real LowLimit,gsl_real UpLimit);
  //////////////////////////////////////////////////////////
  /** 
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  void Comoving_distance_tabulated (real_prec, real_prec, vector<gsl_real>&, vector<gsl_real>&);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  void set_redshift(real_prec red){this->redshift=red;}
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Object of type CosmologicalFunctions
   */
  void set_cosmo_pars(s_CosmologicalParameters spars){this->s_cosmo_pars=spars;}

};
#endif 















