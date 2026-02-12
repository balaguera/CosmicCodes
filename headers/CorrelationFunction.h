
#ifndef __CorrelationFunctionTH__
#define __CorrelationFunctionTH__

# include <iostream>
# include <math.h>
# include <cmath>
# include <cctype>
# include <stdio.h>
# include <fstream>
# include "NumericalMethods.h"
# include "CosmologicalFunctions.h"
# include "BiasFunctions.h"
# include "Hod.h"
# include "DensityProfiles.h"
# include "ScreenOutput.h"
# include "PowerSpectrumTH.h"

#define NPOINTS_EVAL 1000

using namespace std;

class CorrelationFunctionTH{
 private:
    //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  Params par;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  static gsl_real i_Linear_Matter_Correlation_Function(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  static gsl_real i_Non_Linear_Matter_Correlation_Function_Halo_Fit(gsl_real, void*);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  static gsl_real i_Galaxy_Correlation_Function_1h_ss(gsl_real,void *);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  static gsl_real i_Galaxy_Correlation_Function_1h_sc(gsl_real,void *);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  static gsl_real i_Galaxy_Correlation_Function_2h(gsl_real,void *);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  gsl_integration_workspace  *w;
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  gsl_integration_workspace  *wc;

 public:
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
   CorrelationFunctionTH(){}
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */

 CorrelationFunctionTH(Params _par) : par (_par)
    {
      this->w  =  gsl_integration_workspace_alloc  (NPOINTS_EVAL);
      this->wc =  gsl_integration_workspace_alloc  (NPOINTS_EVAL);
    }
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  ~CorrelationFunctionTH(){};
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Linear_Matter_Correlation_Function(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Non_Linear_Matter_Correlation_Function_Halo_Fit(s_CosmologicalParameters *, real_prec);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Galaxy_Correlation_Function_1h_ss(s_CosmologicalParameters *scp, real_prec r);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Galaxy_Correlation_Function_1h_sc(s_CosmologicalParameters *scp, real_prec r);
  //////////////////////////////////////////////////////////
  /**
   *  @name input/output
   */
  real_prec Galaxy_Correlation_Function_2h(s_CosmologicalParameters *scp, real_prec r);
  
};

#endif  
  



