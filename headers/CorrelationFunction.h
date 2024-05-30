
#ifndef __CorrelationFunctionTH__
#define __CorrelationFunctionTH__

# include <iostream>
# include <math.h>
# include <cmath>
# include <cctype>
# include <stdio.h>
# include <fstream>

# include "CosmologicalFunctions.h"
# include "PowerSpectrumTH.h"


using namespace std;

class CorrelationFunctionTH{
 private:
  
  static gsl_real i_Linear_Matter_Correlation_Function(gsl_real, void*);
  static gsl_real i_Non_Linear_Matter_Correlation_Function_Halo_Fit(gsl_real, void*);
  static gsl_real i_Galaxy_Correlation_Function_1h_ss(gsl_real,void *);
  static gsl_real i_Galaxy_Correlation_Function_1h_sc(gsl_real,void *);
  static gsl_real i_Galaxy_Correlation_Function_2h(gsl_real,void *);
  gsl_integration_workspace  *w;
  gsl_integration_workspace  *wc;

 public:
  CorrelationFunctionTH(){}
  CorrelationFunctionTH(int NIT,real_prec L)
    {
      this->w  =  gsl_integration_workspace_alloc  (NIT);
      this->wc =  gsl_integration_workspace_alloc  (NIT);
    }
  ~CorrelationFunctionTH(){}
  
  real_prec Linear_Matter_Correlation_Function(s_CosmologicalParameters *, real_prec);
  real_prec Non_Linear_Matter_Correlation_Function_Halo_Fit(s_CosmologicalParameters *, real_prec);
  real_prec Galaxy_Correlation_Function_1h_ss(s_CosmologicalParameters *scp, real_prec r);
  real_prec Galaxy_Correlation_Function_1h_sc(s_CosmologicalParameters *scp, real_prec r);
  real_prec Galaxy_Correlation_Function_2h(s_CosmologicalParameters *scp, real_prec r);
  
};

#endif  
  



