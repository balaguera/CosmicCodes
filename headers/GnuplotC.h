////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/**
 *  @file Gnuplot.h
 *  @brief The class Parameters. Reads the parameter
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef __GNUPLOTC__
#define __GNUPLOTC__
#include "Structures.h"
#include "../external_libs/gnuplot-iostream/gnuplot-iostream.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


class GnuplotC{
  
 private:
  
  Gnuplot gp;
  
 public:
  /**
   *  @brief default constructor
   *  @return object of class Parameters
   */
  GnuplotC(){}

  
  /**
   *  @brief default destructor
   *  @return none
   */
  ~GnuplotC(){}
  
  /**
   *  @brief default destructor
   *  @return none
   */
  void plot_scaling_relations(string);
  void plot_power_spectrum(vector<real_prec>&, vector<real_prec>&);
  void plot_power_spectrum_redshift_space(vector<real_prec>&, vector<real_prec>&,vector<real_prec>&, vector<real_prec>&);
  void plot_power_spectrum_redshift_space(vector<real_prec>&, vector<real_prec>&,vector<real_prec>&);
  /*
    void plot_abundance(vector<real_prec>&,vector<real_prec>&);
    void plot_abundance(vector<real_prec>&,vector<real_prec>&,vector<real_prec>&);
  
  */

};



#endif
