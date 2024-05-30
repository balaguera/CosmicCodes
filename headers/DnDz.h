


#ifndef __DNDZ__
#define __DNDZ__


# include "CoordinateSystem.h"
# include "NumericalMethods.h"
# include "FileOutput.h"
# include "Catalog.h"



/**
* @struct<s_dndz>
* @file  Type_Structures.h
* @brief The structure file
* @author Andres Balaguera-Antolínez
* @version   1.0
* @date      2020
*/

struct s_dndz{
  /// Container with the galaxy catalogue
  vector<s_Halo > properties;
  
  /// the number of columns in the galaxy catalogue
  int n_columns;
  
  /// First estimate of the parameter alpha, N_ran / N_gal
  real_prec alpha0;

  //////////////////////////////////////////////////////////

  /// Is the nbar tabulated?
  bool nbar_tabulated;
  //////////////////////////////////////////////////////////

  /// Measure the reshift distribution
  bool compute_dndz;

  //////////////////////////////////////////////////////////
  /// Has the survey a constnat depth
  bool constant_depth;

  //////////////////////////////////////////////////////////
  /// Vector of redshiftsa
  vector<gsl_real> zz;

  //////////////////////////////////////////////////////////
  /// Vector of comoving distance evaluated at zz
  vector<gsl_real> rc;

  //////////////////////////////////////////////////////////
  /// Number of redshift bins to measure dNdz
  int n_dndz;

  //////////////////////////////////////////////////////////
  /// Number of redshift bins to provide smooth version of dNdz
  int new_n_dndz;

  //////////////////////////////////////////////////////////
  /// Minimum redshift of the sample
  real_prec redshift_min_sample;

  //////////////////////////////////////////////////////////
  /// Maximum redshift of the sample
  real_prec redshift_max_sample;

  //////////////////////////////////////////////////////////
  /// Area of the survey
  real_prec area_survey;

  //////////////////////////////////////////////////////////
  /// System of coordinate of the random catalogue
  int sys_of_coord_r;

  //////////////////////////////////////////////////////////
  /// the column where the first random coordinate (according to the system of coordinates of the catalog) is written. Here first column must come as 0.
  int i_coord1_r;

  //////////////////////////////////////////////////////////
  /// the column where the second random coordinate (according to the system of coordinates of the catalog) is written. Here first column must come as 0.
  int i_coord2_r;

  //////////////////////////////////////////////////////////
  // the column where the second random coordinate (according to the system of coordinates of the catalog) is written. Here third column must come as 0.
  int i_coord3_r;

  //////////////////////////////////////////////////////////
  /// Area of one Healpix pixel
  real_prec area_pixel;

  //////////////////////////////////////////////////////////
  /// Number of pixels used in thee stimation of dNdz for depth varying sample
  long npixels;

  //////////////////////////////////////////////////////////
  /// Nside, healpix parameter
  long nside;

  //////////////////////////////////////////////////////////
  /// Output file for the dNdz
  string file_dndz;
};




/**
* @class<DnDz>
* @file DndDz.h
* @brief Header file for the class DnDz::
* @title Functions related to the generation and measuremetns of redshift distributions and mean number densities
* @author Andres Balaguera-Antolínez
* @version 1.0
* @date    2020
* @details: This is an example of a main function to call bam. A file called cosmicatlas.cpp
*/




class DnDz{

 private:

  /// Object of type FileOutput
  FileOutput Fm;
  
 
 public:
  
  /**
   * @brief Constructor
   */
  DnDz(){}

  /**
   * @brief Destructor
   */
  ~DnDz(){}
  
  
  /**
   * @brief DnDz
   * @details  This function computes the the mean number density
   * from the random catalogue, in order to be interpolated
   * at the position of the galaxies and randoms in the 
   * estimate of power spectrum. 
   * @param void pointer ps to be recasted in type struct s_dndz
   * @result zv vector with redshifts of smoothed mean number density
   * @result nzv vector with smooth mean number density
   * @result Mzv container with estimates of mean number density in pixels and redshift bins
   */
  void dndz(void *ps, vector<gsl_real>&zv, vector<gsl_real>&nzv, vector< vector<gsl_real> > &Mzv);
};



#endif 
