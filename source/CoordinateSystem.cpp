//////////////////////////////////////////////////////////
/**
 * @brief Header file for the functions in CoordinateTransform.cpp
 *  @file CoordinateSystem.cpp
 *  @brief Methods of the class CoordinateSystem
 *  @author Andres Balaguera-Antolínez (ABA)
 *  @details This file contains functions used to transform coordiante system in the catalogues
 */
//////////////////////////////////////////////////////////
#include "../headers/CoordinateSystem.h"
//////////////////////////////////////////////////////////
void equatorial_to_cartesian(real_prec ra, real_prec dec, real_prec r, real_prec &x, real_prec &y, real_prec &z)
{
  real_prec fac=M_PI/180.;
  x  =r*cos(ra*fac)*sin(0.5*M_PI-dec*fac);
  y  =r*sin(ra*fac)*sin(0.5*M_PI-dec*fac);
  z  =r*cos(0.5*M_PI-dec*fac);
}
//////////////////////////////////////////////////////////
void new_equatorial_to_cartesian(real_prec ra, real_prec dec, real_prec r, real_prec m_ra, real_prec m_dec, real_prec &x, real_prec &y, real_prec &z){
  real_prec xx, yy, zz;
  real_prec new_dec,new_ra;
  equatorial_to_equatorial(ra, dec, m_ra,m_dec,&new_ra,&new_dec); 
  equatorial_to_cartesian(new_ra,new_dec,r, x, y, z);
}
//////////////////////////////////////////////////////////
void equatorial_to_galactic(real_prec alpha, real_prec delta, real_prec *b, real_prec *l)
{
  real_prec fac = M_PI/180.;                         /* Convierte grados a radianes*/
  real_prec lly,llx,bb;
  real_prec alpha_p=fac*192.859508;                        /*Right ascention of galactic north pole in radians*/
  real_prec delta_p=fac*27.128336;                         /*Declination of galactic north pole in radians*/
  real_prec  l_n   =fac*122.932;                           /*Galactic longitude of the celestial pole, in radians*/
  bb = sin(delta)*sin(delta_p)+cos(delta)*cos(delta_p)*cos(alpha-alpha_p);
  lly= cos(delta)*sin(alpha-alpha_p);
  llx=-cos(delta)*sin(delta_p)*cos(alpha-alpha_p)+sin(delta)*cos(delta_p);
  *l = l_n-atan2(lly,llx);
  *b = asin(bb);
}
//////////////////////////////////////////////////////////
void equatorial_to_equatorial(real_prec old_ra, real_prec old_dec, real_prec m_ra, real_prec m_dec, real_prec *new_ra, real_prec *new_dec)
{
  real_prec fac = M_PI/180.0;                         /* Convierte grados a radianes*/
  real_prec aux1=asin(sin(old_dec*fac)*sin(m_dec*fac)-cos(old_dec*fac)*cos(m_dec*fac)*sin(old_ra*fac-m_ra*fac));
  real_prec aux2=cos(old_dec*fac)*cos(old_ra*fac-m_ra*fac)/cos(aux1);
  real_prec aux3=0;
  aux2=aux3-acos(aux2);
  *new_dec=aux1/fac;
  *new_ra =aux2/fac;
}
//////////////////////////////////////////////////////////
void galactoc_to_equatorial(real_prec  b, real_prec  l, real_prec  &ra, real_prec  &dec){
  real_prec fac = M_PI/180.0;                         /* Convierte grados a radianes*/
  real_prec  alpha_p=fac*192.859508;                        /*Right ascention of galactic north pole in radians*/
  real_prec  delta_p=fac*27.128336;                         /*Declination of galactic north pole in radians*/
  real_prec  l_n   =fac*122.932;                           /*Galactic longitude of the celestial pole, in radians*/
  real_prec  delta=sin(b*fac)*sin(delta_p)+cos(b*fac)*cos(delta_p)*cos(fac*l-l_n);
  dec=asin(delta);
  real_prec  alpha1=sin(fac*l-l_n);
  real_prec  alpha2=cos(fac*l-l_n)*sin(delta_p)-tan(fac*b)*cos(delta_p);
  ra=atan2(alpha1, alpha2)+12.25;
}



