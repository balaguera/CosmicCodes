////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<ScalingRelations>
 * @file ScalingRelations.cpp
 * @brief Methods of the class ScalnigRELATIONS
 * @details The class Catalog reads and analyses an input catalog od dark matter tracers
 * @author Andres Balaguera Antolinez 2007-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# include "../headers/ScalingRelations.h"
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*MASS-TEMPERATURE RELATION CALIBRATED(?), MASSES IN UNITS OF 10^14 Ms/h: RETURNS TEMPERATURE IN keV*/
real_prec ScalingRelations::M2T(real_prec m, void *p){
  return 0.595*pow(m,(real_prec)0.34);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*Converting mass in 10^{14}/h to the NATURAL LOGARITHM of 
  luminosities (bolometric or band) in units of 10^{44}erg/s/h^2  */
real_prec ScalingRelations::RB_M2L(real_prec m, void *p){/*Reipricht & Bohringer 2002 */
  return log(0.1169*pow(m,(real_prec)1.462));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ScalingRelations::FED_M2L(real_prec m, void *p){/*Fedeli et al*/
  return log(0.0862*pow(m,(real_prec)1.554));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ScalingRelations::STA_M2L(real_prec m, void *p){/*Measured, for Omega_matter=0.24 and sigma_ln M=0.39 leading to the */
  return log(0.1139*pow(m,(real_prec)1.46));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ScalingRelations::MANTZ_BOL_M2L(real_prec m, void *p){
  /*Measured (bolometric) by Mantz et al 2009, also based on XLF of REFLEX data. The factors hubble^2 and 0.1
    come to transform from the units of Mantz ([M]=[10^15 Ms], [L]=10^44 erg/s) to our units ([M]=[10^14 Ms/h], [L]=10^44 erg/s/h^2) */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p;
  return log(pow(s_cp->hubble,2)*1.26*pow(0.1*m/s_cp->hubble,(real_prec)1.59));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ScalingRelations::MANTZ_BAND_M2L(real_prec m, void *p){/*Measured in the ROSAT energy band by Mantz et al 2009, also based on XLF REFLEX data */
  struct s_CosmologicalParameters * s_cp= (struct s_CosmologicalParameters *)p; 
  return log(pow(s_cp->hubble,2)*0.82*pow(0.1*m/s_cp->hubble,(real_prec)1.29));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ScalingRelations::MOCKS_M2L(real_prec m, void *p){ /*Calibrated for the reflex2 mocks to follow the reflex2 luminosity function. */
  // real_prec a_ml       = -1.10;
  // real_prec b_ml       =  1.70;
  // real_prec c_ml       = -0.23;
  real_prec a_ml       = -0.64;
  real_prec b_ml       =  1.27;
  real_prec c_ml       =  0.0;
  return  log(pow(10,a_ml+b_ml*log10(m)+c_ml*pow(log10(m),2)));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
