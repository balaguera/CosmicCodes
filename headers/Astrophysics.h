////////////////////////////////////////////////////////////////////////////
/**
 * @class<Astrophysics>
 * @brief Header file for the class Cwclass::
 * @title Functions related to astrophysical procesess
 * @author ABA
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef _Astrophysics_
#define _Astrophysics_
# include "ScalingRelations.h"
# include "CosmologicalFunctions.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 *@brief
 * @struct<s_astrophysical_parameters>
 * @brief The s_astrophysical_parameters struct
 * @details Structure containing parameters of gas-mass scaling relation in clusters
 */
struct s_astrophysical_parameters{
  real_prec A_gas;
  real_prec B_gas;
  real_prec mstar;
  real_prec sigma_red;
  real_prec sigma_ln;
  real_prec missing_flux;
};
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
class Astrophysics{
 private:

 public:
  Astrophysics(){};
  ~Astrophysics(){};

  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec mass2temp(real_prec, real_prec, void*);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec filtering_mass(real_prec, void*);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec cooling_function(real_prec, void*);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec vc(real_prec, real_prec, void*);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec cooling_radius(real_prec, real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec cooling_mass_rate(real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec infall_mass_rate(real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec baryon_fraction(real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec star_formation_rate(real_prec, real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec baryon_mass_virial_mass_distribution(real_prec, real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */

  real_prec total_mass_nb_mass_relation(real_prec, real_prec,void *, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  //  real_prec Mobs_Mnb_distribution(real_prec, real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec cluster_mass_function(real_prec, real_prec, real_prec, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec mass_lum_distribution(real_prec, real_prec, real_prec, void *, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec mass_lum_distribution_errors(real_prec, real_prec, real_prec, void *, void*);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec G(real_prec, real_prec, real_prec, void *, void *);
  //////////////////////////////////////////////////////////
  /**
    *@public
   * @brief
   */
  real_prec G_Lmin_Lmax(real_prec, real_prec, real_prec, void *, void *);
};
#endif
