#ifdef _USE_LPT_
// **************************************************************
/**
 * @class<LPT>
 * @brief Header file for the class LPT::
 * @file LPT.h
 * @title LPT: Approximated gravity solvers based on Lagrangian perturbation theory
 * @author Francisco-Shu Kitaura
 * @author Andres Balaguera-Antolínez
 * @version   1.0
 * @date      2020`
*/
// **************************************************************
// **************************************************************
// **************************************************************
#ifndef _LPT_
#define _LPT_

#include "bstream.h"
#include <math.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string.h>
#include <cassert>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <netinet/in.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
#include "FileOutput.h"
#include "NumericalMethods.h"
#include "Miscelanious.h"

using namespace std;


//##################################################################################

// *****************************************************************************************************
/*!\class LPT

 * @details   LPT
 * @author    Francisco-Shu Kitaura
 * @author    adapted by Andres Balaguera-Antolinez

 * @brief CLASS LPT
 */


//##################################################################################

class LPT{

  private:

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /* /\** */
  /*  * @brief outoput object */
  /*  *\/ */
  /* ofstream sal; */

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Params params;
  //////////////////////////////////////////////////////////

  /**
   * @brief Object of class FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_power;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  Gnuplot gp_pdf;
    //////////////////////////////////////////////////////////
  /**
   * @brief Object of class Cosmology
   */
  Cosmology Cosmo;
  //////////////////////////////////////////////////////////
  /**
   * @brief Structure allocating cosmological information
   */
  s_CosmoInfo s_cosmo_info; 
    //////////////////////////////////////////////////////////
  /**
   * @brief Structure allocating cosmological parameters
   */
  s_CosmologicalParameters s_cosmo_pars;
  //////////////////////////////////////////////////////////
  /**
   * @brief Input parameter file adapted to LPT. Ideally, to be merged with that of BAM
   */
  string par_file;

  //////////////////////////////////////////////////////////
  /**
   *  @brief ScreenOutput obejct
   */
  string fname_MOCK_NCOUNTS;
  //////////////////////////////////////////////////////////
  /**
   *  @brief ScreenOutput obejct
   */
  string fname_MOCK_MASS;
  //////////////////////////////////////////////////////////
  /**
   *  @brief ScreenOutput obejct
   */
  ScreenOutput So;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec velbias;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec dkbin;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  int N_bin;
  
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string buffsf;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec L1;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec L2;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec L3;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG N1;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG N2;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG N3;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG NGRID;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  ULONG seed;


    ////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the DM produced by LPT
   */
  string fnamePOSX;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
 string fnamePOSX_RS;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string fnameVXpart;
  ////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the DM produced by LPT
   */
  string fnamePOSY;
   //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnamePOSY_RS;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string fnameVYpart;
  ////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the DM produced by LPT
   */
  string fnamePOSZ;
   //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnamePOSZ_RS;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string fnameVZpart;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the VX interpolated on the mesh produced by LPT
   */
  string fnameVX;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the DM produced by LPT
   */
  string fnameVY;
  //////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the DM produced by LPT
   */
  string fnameVZ;

  //////////////////////////////////////////////////////////
  /**
   * @brief Name of file containing the Fourier grid with the 3D power spectrum therein interpolated. No ".dat"
   */
   string fname3DPOWER;

   //////////////////////////////////////////////////////////
  /**
   * @brief Name of file the White noise in conf space
   */
   string fnameIC;

   //////////////////////////////////////////////////////////
   /**
    * @brief Name of file the initial density field displaying an initial power spectrum
    */
   string fnameICDELTA;

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */

   string fnameDM;
   //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
   string fnameDM_RS;

   //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
   string fname2LPTTERM;
   //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
   string fnameTHETA;
   //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
   string fnameDMNGP;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameDMNGP_RS;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameS2TERM;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameS2TERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameP2TERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameS3TERM;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameS3TERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameSTTERM;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameSTTERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnamePSITERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fname2LPTTERMEUL;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string stradd;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string stradd_bam;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string fnameTRACERCAT;

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   ULONG Number_of_Tracers;


  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool planepar;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  real_prec Initial_Redshift_DELTA;
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  bool Normalize_IC_to_initial_redshift;
  real_prec growth_ini=1;
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  vector<real_prec>tab;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////

 public:

  /**
   *  @brief default constructor
   *  @brief object of class LPT. I use this to defined a PACHY objects as a BAM class member
   *  @brief
   */
  LPT(){}
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *  @brief constructor
   *  @param parameters_file parameter file
   *  @return object of class LPT
   */


 LPT(Params _params, s_CosmoInfo _s_cosmo_info) :s_cosmo_info(_s_cosmo_info)
  {
   this->set_params(_params);
   this->NGRID= this->params.d_NGRID();
   time_t time_bam;
   time(&time_bam);
   this->So.initial_time=time_bam;
   this->gp_power<<"set border linewidth 1.5\n";
   this->gp_pdf<<"set border linewidth 1.5\n";
   this->Cosmo.set_cosmo_pars(this->params.s_cosmo_pars);

 }

 //////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////
  /**
     *  @brief Default destructor
     *  @return
   */
  ~LPT(){}

 //////////////////////////////////////////////////////////
 //////////////////////////////////////////////////////////
 /**
  * @brief Reading and initializing structures with params for LPT
  */

  void set_params(Params pars);
  void set_params();
 //////////////////////////////////////////////////////////
 /**
  * @brief Set names of files. Strings are private class members
  */
  void set_fnames();
  //////////////////////////////////////////////////////////

  /**
   *  @brief Generate DM density field according to some SF model and IC
   */

#ifdef OMPPARRAN
 void get_dm_field(gsl_rng ** gBaseRand);
#else
  void get_dm_field(gsl_rng * gBaseRand);
#endif


  //////////////////////////////////////////////////////////

  /**
   *  @brief Generate DM density field according to some SF model and IC
   */

#ifdef OMPPARRAN
 void get_displacement(gsl_rng ** gBaseRand, vector<real_prec>&, int it);
#else
 void get_displacement(gsl_rng * gBaseRand, vector<real_prec>&, int it);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  vector<real_prec>Displacement;


  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void  displace_part_comp(vector<real_prec> &posi, vector<real_prec>&psii, bool periodic,int comp);
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
#ifndef _DISPLACEMENTS_
  void Lag2Eul_comp(real_prec biasLAG2,real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp);
#else
  void Lag2Eul_comp(real_prec biasLAG2,real_prec kth,int ftype,bool periodic,int comp);
#endif
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void Lag2Eul_vel(real_prec biasLAG2,real_prec kth,int ftype,bool periodic);
  //////////////////////////////////////////////////////////
  /**
   *  @brief
   */
  void Lag2Eul_compB(real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp);
  //////////////////////////////////////////////////////////
  /**
   *  @brief LPT routine to generate catalog with positions and velocities
   */
#ifdef OMPPARRANRSD
  void makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir);
#else
  void makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand, int ir);
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
#ifdef _USE_OMP_
  void makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir);
  void makecat_withv_new(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand,int ir);
#else
  void makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand, int ir);
  void makecat_withv_new(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand,int ir);
#endif

  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void read_tabulated_power();
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  real_prec linInterp(real_prec xpos, real_prec ypos, real_prec zpos, const vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void theta2velcomp(vector<real_prec> & delta, vector<real_prec> &vei, bool zeropad, bool norm, int comp);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void comp_velbias(vector<real_prec> &delta, vector<real_prec>&out, bool zeropad, bool norm);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void normalize_df_z_ini(vector<real_prec>&, string type, real_prec target_z_ini, real_prec target_z_end);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void DM_to_RSS(int los);

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   void set_fname_MOCK_NCOUNTS(string name){this->fname_MOCK_NCOUNTS=name;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string _fname_MOCK_NCOUNTS(){return this->fname_MOCK_NCOUNTS;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   void set_fname_MOCK_MASS(string name){this->fname_MOCK_MASS=name;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  string _fname_MOCK_MASS(){return this->fname_MOCK_MASS;}




  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   void set_s_cosmo_info(s_CosmoInfo cinf){this->s_cosmo_info=cinf;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnamePOSX(){return this->fnamePOSX;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnamePOSY(){return this->fnamePOSY;}

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnamePOSZ(){return this->fnamePOSZ;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVXpart(){return this->fnameVXpart;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVYpart(){return this->fnameVYpart;}

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVZpart(){return this->fnameVZpart;}
   //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVX(){return this->fnameVXpart;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVY(){return this->fnameVYpart;}

  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameVZ(){return this->fnameVZpart;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameTRACERCAT(){return this->fnameTRACERCAT;}
   void set_fnameTRACERCAT(string new_f){this->fnameTRACERCAT=new_f;}
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _fnameDM(){return this->fnameDM;}
   void set_fnameDM(string new_f){this->fnameDM=new_f;}
//////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _stradd(){return this->stradd;}
   void set_stradd(string new_f){this->stradd=new_f;}
//////////////////////////////////////////////////////////
  /**
   * @brief  
   */
   string _stradd_bam(){return this->stradd_bam;}
   void set_stradd_bam(string new_f){this->stradd_bam=new_f;}


   int _seed(){return this->params._seed();}


};


#endif

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
#endif
