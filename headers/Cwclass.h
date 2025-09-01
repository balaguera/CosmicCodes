////////////////////////////////////////////////////////////////////////////
/**
 * @class<Cwclass>
 * @brief Header file for the class Cwclass::
 * @file Cwclass.h
 * @title Functions related to Cosmic Web and V-WWeb classification
 * @author   ABA
 */
////////////////////////////////////////////////////////////////////////////
#ifndef _Cwclass__
#define _Cwclass__

# include <iostream>
# include <fstream>
# include <iostream>
# include <vector>
# include <algorithm>
# include <math.h>
# include <sstream>
# include <iomanip>
# include <string>
# include "def.h"
# include "NumericalMethods.h"
# include "Miscelanious.h"
# include "ScreenOutput.h"
# include "FileOutput.h"
# include "Params.h"

using namespace std;

class Cwclass {
  
 private:

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief outoput object
   */
  ofstream sal;

  //////////////////////////////////////////////////////////
  /**
    *@private
   * @brief Params object
   */
  Params params;


  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of CWT used, computed as the size of the
   * Bam::cwt_used container
   */

  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Convolution with Kernel used in CWT.
   * @details The kernel used for convlution is a public containers
   */
  void Konvolve(vector<real_prec> &in, vector<real_prec>&out, string type);
  
  //////////////////////////////////////////////////////////
  /**
   * @public
   *  @brief Object of type ScreenOutput.
   */
  ScreenOutput So;

  //////////////////////////////////////////////////////////
  /**
   *  @brief Fraction of cells defined as knots
   */
  real_prec knots_fraction;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Fraction of cells defined as filaments
   */
  real_prec filaments_fraction;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Fraction of cells defined as sheets
   */
  real_prec sheets_fraction;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Fraction of cells defined as voids
   */
  real_prec voids_fraction;
  //////////////////////////////////////////////////////////
  /**
   *  @brief Kernel used in convolution
   */
  vector<real_prec>Kernel;
  //////////////////////////////////////////////////////////

  
 public:

  //////////////////////////////////////////////////////////
  /**
    * @public
   *  @brief default constructor
   */
  Cwclass(){}
  
  //////////////////////////////////////////////////////////
  /**
    * @public
   *  @brief Use this constructor to pass a Params objects
   *  @details This constructor initializes the value of NTT and NGRID
   *  @param Inizialization of private variables
   */
  
 Cwclass(Params _par)
  {
    this->params=_par;
    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;
    sal.setf(ios::showpoint); 
    sal.setf(ios::scientific); 
    this->cwv_used=this->params._cwv_used();
    this->cwt_used=this->params._cwt_used();
    this->s_cosmo_pars=this->params.s_cosmo_pars;
 }

  //////////////////////////////////////////////////////////
  /**
    * @public
   *  @brief Default constructor
   */
  ~Cwclass(){}
    //////////////////////////////////////////////////////////
  /**
   * @brief Auxiliary s_CosmologicalParameters
   */
  s_CosmologicalParameters s_cosmo_pars;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Auxiliary s_CosmologicalInfo structure
   */
  s_CosmoInfo s_cosmo_info;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container to allocate the cosmic web classification for each cell.
   * @brief Example: for knots
   * @code
    for (ULONG index=0;index<NGRID;++index)
     if (lambda1[index]> lambdath && lambda2[index]> lambdath && lambda3[index]>lambdath) // lambdath = Threshold
       CWClass[index]=I_KNOT;                                                             // I_KNOT   = 1 Preproc definition
    *@endcode
   */
  vector<ULONG>CWClass;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container allocating the different T-CW classifications
   * @brief required in the parameter file.
   */
  vector<ULONG> cwt_used;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container allocating the different V-CW classifications
   * @brief required in the parameter file.
   */
  vector<ULONG> cwv_used;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container to allocate the V-cosmic web classification for each cell.
   * @brief Example: for v-knots
   * @code
    for (ULONG index=0;index<NGRID;++index)
     if (lambda1_vs[index]> lambdath_v && lambda2_vs[index]> lambdath_v && lambda3_vs[index]>lambdath_v) // lambdath_v = V-Threshold
       CWClass_V[index]=I_KNOT;                                                                          // I_KNOT     = 1 Preproc definition
    *@endcode
   */
  vector<int>CWClass_V;
    //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_1 on the grid
   */
  vector<real_prec> lambda1;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_2 on the grid
   */
  vector<real_prec> lambda2;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_3 on the grid
   */
  vector<real_prec> lambda3;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the potential
   */
  vector<real_prec> potential;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the S2 term
   */
  vector<real_prec> S2;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the S3 term
   */
  vector<real_prec> S3;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the S2*delta term
   */
  vector<real_prec> S2DELTA;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the delta² term
   */
  vector<real_prec> DELTA2;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the delta³ term
   */
  vector<real_prec> DELTA3;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant of the tidal field I
    *@details The Invariant I of the tidal filed is the overdensity
   */
  vector<real_prec> Invariant_TF_I;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant of the tidal field II
    *@details The Invariant II of the T-field is the reduced sum of the product of two different eigenvalues
   */
  vector<real_prec> Invariant_TF_II;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant of the tidal field III
    *@details The Invariant II of the T-field is product of the tree  eigenvalues
   */
  vector<real_prec> Invariant_TF_III;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant of the tidal field III
    *@details The Invariant II of the T-field is product of the tree  eigenvalues
   */
  vector<real_prec> Invariant_TF_IV;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant of the shear of the velocity field I
    *@details In the linear theory approximation this is the overdensity   */

  vector<real_prec> Invariant_VS_I;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant II of the shear of the velocity field
   */
  vector<real_prec> Invariant_VS_II;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the Invariant III of the shear of the velocity field
   */
  vector<real_prec> Invariant_VS_III;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the tidal anisotropy
   */
  vector<real_prec> Tidal_Anisotropy;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the ELLIPTICITY
   */
  vector<real_prec> Ellipticity;
  //////////////////////////////////////////////////////////

  /**
   * @public
   * @brief Container for the ELLIPTICITY
   */
  vector<real_prec> Prolatness;


    //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for Nbala²delta on the grid
   */
  vector<real_prec> N2D;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_1 of the velocity shear on the grid
   */
  vector<real_prec> lambda1_vs;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_2 of the velocity shear on the grid
   */
  vector<real_prec> lambda2_vs;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the eigenvalue lambda_3 of the velocity shear on the grid
   */
  vector<real_prec> lambda3_vs;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for divergence of the velocity field on the mnesh
   */
  vector<real_prec>Divergence_VelField;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container for the dispersion of the velocity field on the mesh
   */
  vector<real_prec>Dispersion_VelField;


  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Minimum mass in percolated knots
   */
  real_prec knot_mass_min;
  //////////////////////////////////////////////////////////

  /**
   * @public
   * @brief Maximum mass in percolated knots
   */
  real_prec knot_mass_max;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container loaded by the function get_Mk_collapsing_regions()
   * containing the bin in Mk in which a given cell is found, according to the results of the FoF.
   */
  vector<ULONG> SKNOT_M_info;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Container loaded by the function get_V-Mk_collapsing_regions()
   */
  vector<ULONG> VDISP_KNOT_info;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Inform whether a cell is classified with a given CWT
   * @param int sua: index of CWC as requested from the param file, in the option v_CWT_used = ... END
   * @param ULONG i: index running through the mesh
   * @return true/false if the cell i is classified as Bam::cwt_used[sua]
   */
  bool get_cell_classified(int sua, ULONG i);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief This function provides the index identifying the T-CWT classification
   * as requiested in the parameter file, for each cell, and according
   * to the CWClass array. Example, if requested 1 2 34, the allowed indices are 0,1,2.
   * @param i: index (row-major order, i<Bam::NGRID) identifying the cell. Example:
   * @code
      for(ULONG ig = 0; ig< NOBJS ; ++ig)  // Loop over tracers
      {
        int I_CWT=0;
        ULONG ID= tracer.Halo[ig].GridID;    // Get the ID of tracers [0,NOBJS-1]. tracer is type Catalog
        I_CWT=wclass.get_Tclassification(ID);// wclass is type Cwclass
      }
    *@endcode
    *@return Container (size NGRID) with the index of classifiction (according to input par file) for each cell
   */
  int get_Tclassification(ULONG i);

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief This function provides the index identifying the V-CWT classification
   * as requiested in the parameter file, for each cell, and according
   * to the CWClass array
   * @param i: index (row-major order, i<Bam::NGRID) identifying the cell.
   */
  int get_Vclassification(int i);

  /////////////////////////////////////////////////////////
  /**
   * @public
   * @brief FoF for Knots
   * @param rho: the density (if No=0) or the overdensity (with mean No)
   * @param No: the number of objects (different of zero if parameter rho is the density field itself)
   * @return This function loads the container Bam::SKNOT_M_info of size Bam::NGRID, which contains the bin in MK in which
   * the cell i is identified according to the FoF
   */
  void get_Mk_collapsing_regions(vector<real_prec> &rho, real_prec No);
  void get_collapsing_regions(ULONG &, vector<real_prec>&sk);
  /////////////////////////////////////////////////////////

  /**
   * @public
   * @brief FoF for Knots
   * @param rho: the density (if No=0) or the overdensity (with mean No)
   * @param No: the number of objects (different of zero if parameter rho is the density field itself)
   * @return This function loads the container Bam::SKNOT_M_info of size Bam::NGRID, which contains the bin in MK in which
   * the cell i is identified according to the FoF
   */
  void get_SigmaVel_collapsing_regions(const vector<real_prec>&in, vector<real_prec>&Vx, vector<real_prec>&Vy,vector<real_prec>&Vz, real_prec Nmean);
  /////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Performs the Cosmic-Web classification
   * @param delta: Overdensity field of the dark matter distribution interpolated on a grid with size Bam::NGRID
   * @returns This function loads the container Bam::CWClass of size Bam::NGRID specifying in
   * @details Performs the Cosmic-Web classificatio by solving Poisson equation using FFTW, identifying CW using the eigenvalues of the tidal tensor and classifying the different CW given a threshold Bam::lambdath
   * @details each cell the type of CW with 1, 2, 3, 4, for knots, filaments, sheets and voids, respectively.
   */
  void get_CWC(vector<real_prec>&delta);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Performs the Cosmic-Web classification based on the shear of the velocity field
   * @details Performs the Cosmic-Web classificatio from the shear of the velocity field identifying CW using the eigenvalues of the tidal tensor and classifying the different CW given a threshold Bam::lambdath
   * @param delta: Overdensity field of the dark matter distribution interpolated on a grid with size Bam::NGRID
   * @returns This function loads the container Bam::CWClass_V of size Bam::NGRID specifying in
   * each cell the type of CW_V with 1, 2, 3, 4, for (dynamical) knots, filaments, sheets and voids, respectively.
   */
  void get_CWC_V(vector<real_prec>&Vx,vector<real_prec>&Vy,vector<real_prec>&Vz);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief In this function we compute terms associated to the bias expansion ala McDonald.
   * @param The matter overdensity (vector)
   * @return cvontainers for the terms S2 and Ð²ð
   * The other temrs such as ð² and ð³ are directly computed in Bam.cpp when called, i.e, in get_bias and get_mock
   */
  void get_bias_terms(vector<real_prec>&delta);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Get tidal anisotropy from input binary to write into output bin
   */
  void get_tidal_anisotropy(string inp, string out);
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Get tidal anisotropy from input density_field and write into output bin
   */
  void get_tidal_anisotropy(vector<real_prec>&den, string out);

  ////////////////////////////////////////////////////////////////////////////
  void get_tidal_anisotropy(vector<real_prec>&delta);

  ////////////////////////////////////////////////////////////////////////////
  void get_tidal_anisotropy(vector<real_prec>&dens, vector<real_prec>&tidal);




  /////////////////////////////////////////////////////////
  /**
   * @public
   * @brief This identifies the iteration within BAM
   */
  int step;

  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Volume filling factors for knots
   */
  real_prec volume_knots;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Volume filling factors for sheets
   */
  real_prec volume_sheets;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Volume filling factors for filaments
   */
  real_prec volume_filaments;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Volume filling factors for voids
   */
  real_prec volume_voids;
  //////////////////////////////////////////////////////////
  /**
   * @brief Set params. 
   * @detail I have to do this when defining this class as a class membr (e.g, a Bam class member) abnd I cannot pass params through the constructor
   */
  void set_params(Params new_params){
    this->params=new_params;
    this->cwv_used=this->params._cwv_used();
    this->cwt_used=this->params._cwt_used();
    this->s_cosmo_pars=this->params.s_cosmo_pars;
  }
//////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec _knots_fraction(){return this->knots_fraction;}  
  real_prec _filaments_fraction(){return this->filaments_fraction;}  
  real_prec _sheets_fraction(){return this->sheets_fraction;}  
  real_prec _voids_fraction(){return this->voids_fraction;}  
  //////////////////////////////////////////////////////////



};
  
#endif
  
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

