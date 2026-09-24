/**
 *  @class <ScalingRelations>
 *  @ingroup classes
 *  @brief This class contains methods related to the link between tracer properties
 *  @file DensityProfiles.h
 *  @author Andres Balaguera-Antolínez
*/
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef __SCALINGRELATIONS__
#define __SCALINGRELATIONS__
////////////////////////////////////////////////////////////////////////////
# include "CosmologicalFunctions.hpp"
////////////////////////////////////////////////////////////////////////////
class ScalingRelations{
private:
    real_prec m;
public:
  ScalingRelations(){};
  ~ScalingRelations(){};
  real_prec M2T(real_prec, void *);
  real_prec RB_M2L(real_prec, void *);
  real_prec FED_M2L(real_prec, void *);
  real_prec STA_M2L(real_prec, void *);
  real_prec MANTZ_BOL_M2L(real_prec, void *);
  real_prec MOCKS_M2L(real_prec, void *);
  real_prec MANTZ_BAND_M2L(real_prec , void *);
};
#endif
////////////////////////////////////////////////////////////////////////////
