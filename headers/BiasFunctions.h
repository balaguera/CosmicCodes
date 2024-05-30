////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 * @class
 * @file BiasFunctions.h
 * @brief Header file for the class MassBiasFunction::
 * @title Functions related to interpolations of halo bias as a function of halo mass
 * @author Andres Balaguera-Antolínez
 * @version  1.0
 * @date     2017-2019
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "NumericalMethods.h"
# include "CosmologicalFunctions.h"

#ifndef _MASS_BIAS_FUNCTIONS_
#define _MASS_BIAS_FUNCTIONS_

class MASS_BIAS_FUNCTIONS {
 public:
  MASS_BIAS_FUNCTIONS(){}
  /*HERE, nu=(delta_sc/sigma)^2. In some references, nu=delta_sc/sigma, so be careful
    Ojo, nu*f(nu)=0.5*f(sigma), y es nu*f(nu) lo que paso abajo como funcion de  masa*/
  /*HALO MASS FUNCTIONS*/
  real_prec mass_function(real_prec, real_prec, void*);
  real_prec dm_h_bias(real_prec, real_prec, void *);
  ~MASS_BIAS_FUNCTIONS(){}
};

#endif
