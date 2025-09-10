//////////////////////////////////////////////////////////
/**
 * @brief Header file
*  @file CoordinateSystem.h
 * @details  Functions used to transform coordiante system in the catalogues
 * @author Andres Balaguera-Antolínez (ABA)
 * @version 1.0
 * @date  2012-2023
 */
//////////////////////////////////////////////////////////
#ifndef __COORDINATE_SYSTEM__
#define __COORDINATE_SYSTEM__
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>
using namespace std;
# include "NumericalMethods.h"
# include "Miscelanious.h"
//////////////////////////////////////////////////////////
 /**
   * @brief Convert equatorial coordinates, i.e, right ascension, declination and distance  
   * to cartessian coordinates. 
   * @details Input values of RA and Dec are in degrees,          
   * the output is in the units of r       
   */
  void equatorial_to_cartesian(real_prec, real_prec, real_prec, real_prec &, real_prec &, real_prec &);
//////////////////////////////////////////////////////////
 /**
 /**
   * @brief Convert cartessian coordinates  to equatorial coordinates. 
   */
  void cartesian_to_equatorial(real_prec, real_prec, real_prec, real_prec &, real_prec &, real_prec &);
//////////////////////////////////////////////////////////
 /**
   * @brief Convert cartessian coordinates  to spherical coordinates. 
   */
  void cartesian_to_spherical(real_prec, real_prec, real_prec, real_prec &, real_prec &, real_prec &);
///////////////////////////////////////////////////////
   /**
   * @brief Convert cartesian coordinates into equatorial RA, DEC, r  
   * to a new spherical coordinate system 
  */
  void new_equatorial_to_cartesian(real_prec, real_prec, real_prec , real_prec, real_prec, real_prec &, real_prec &, real_prec &);
///////////////////////////////////////////////////////
  /**
   * @brief Convert equatorial coordinates, i.e, right ascension, declination and distance
   * to galactic coordinates . 
   * @details Input values of RA and Dec are in degrees,
   * the output are in radians      
   */
  void equatorial_to_galactic(real_prec, real_prec, real_prec *, real_prec *);
///////////////////////////////////////////////////////
  /**
   * @brief Convert equatorial coordinates, i.e, right ascension and declination           
   * to a new equatorial coordinate 
   * @details in which the axis north pole has angular coordinates  real_prec m_ra, real_prec m_dec.
   * @param RA in Degrees,
   * @param Dec in Degrees,
   */
 void equatorial_to_equatorial(real_prec, real_prec, real_prec, real_prec, real_prec *, real_prec *);
 ///////////////////////////////////////////////////////
   /**
    * @brief Convert galactic coordinates to a  equatorial coordinate
    * @details in which the axis north pole has angular coordinates  real_prec m_ra, real_prec m_dec.
    * @param RA in Degrees,
    * @param Dec in Degrees,
    */
  void galactic_to_equatorial(real_prec, real_prec,  real_prec &, real_prec &);
#endif
