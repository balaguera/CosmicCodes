#ifndef _LSS_VECTORS_
#define _LSS_VECTORS_

# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <vector>

using namespace std;

// *******************************************************************
// *******************************************************************
/**
 *@brief This class is designed to manipultae density fields and its statistics.
 */

class lss_vector{
  
 protected:
  
 private:
  

  
  // *******************************************************************
  /**
   *@brief Size of the container
   */
  ULONG size;
  // *******************************************************************
  /**
   *@brief Type of the density field: 
   *@options: density (DENSITY), overdensity (OVERDENSITY)
   */
  string type;
  // *******************************************************************
  /**
   *@brief Size of the container
   */
  vector<real_prec> field;
  
  // *******************************************************************
  /**
   *@brief Size of the container
   */
  //  PowerSpectrumF PS;

  // *******************************************************************
  // *******************************************************************
  // *******************************************************************
  // *******************************************************************

 public:
  
  /**
   *@brief Constructor, passing the sise of the container
   */
 lss_vector(ULONG nsize, string ntype):size(nsize), type(ntype){
    this->field.resize(nsize,0);
  };
  // *******************************************************************
  /**
   *@brief Default constructor
   */
  ~lss_vector(){};
  // *******************************************************************
  /**
   *@brief 
   */
  void get_overdensity();
  // *******************************************************************
  /**
   *@brief 
   */
  int number_counts(int);
  // *******************************************************************
  /**
   *@brief 
   */
  void lss_power_spectrum();
  // *******************************************************************
  /**
   *@brief 
   */
  real_prec kvector(int);
};


#endif
