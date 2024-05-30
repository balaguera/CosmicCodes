////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef __SCALINGRELATIONS__
#define __SCALINGRELATIONS__
////////////////////////////////////////////////////////////////////////////
# include "CosmologicalFunctions.h"
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
