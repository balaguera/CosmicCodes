//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
* @file cosmo_parametrs.h
* @title Cosmological parameters from different LSS analysis
* @author Andres Balaguera-Antolínez
* @version 1.0
* @date    2020
*/
//////////////////////////////////////////////////////////
#include <math.h>
using namespace std;
//////////////////////////////////////////////////////////
#ifndef __COSMO_PARAMETERS__

#define __COSMO_PARAMETERS__


//////////////////////////////////////////////////////////
// These are the params of the Unitsim
/**
* @namespace<Cosmo_parameters_PLANCK>
* @brief Cosmological parameter from PLANCK mission 2016
*/
namespace Cosmo_parameters_PLANCK
{
/**
    * @brief Total matter energy density parameter UNITSIM
    */
  double Om_matter=0.3089;
  double N_eff=3.15;
  double Om_radiation=(1.3157e-5)*(1.+(7./8.)*std::pow(4./11.,4./3.)*N_eff);
  double Om_baryons=0.04859;
  double wde_eos=-1;
  double Om_vac=0.6911 ; //1.-Om_matter-Om_k-Om_radiation;
  double Om_k=0.0;// 1.-Om_matter-Om_vac-Om_radiation;
  double Om_cdm=1-Om_matter-Om_baryons-Om_radiation;
  // hubble constant H/100km/s/Mpc
  double hubble=0.6774;
  double Hubble=100.0;
  // fraction of baryonic energy density to matter
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.725;
  double spectral_index = 0.9667;
  double alpha_s = 0.;
  double sigma8= 0.8147;
  double RR=8.0;
  double Delta_SO = 200;
  bool use_wiggles = true; 
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
* @namespace<Cosmo_parameters_MG>
* @brief Cosmological parameters usedd in MOdified gravity simulations
*/
namespace Cosmo_parameters_MG
{
/**
    * @brief Total matter energy density parameter MOdified gravirt Sims
    */
  double Om_matter=0.31899;
  double N_eff=3.046;
  double Om_radiation=(1.3157e-5)*(1.+(7./8.)*std::pow(4./11.,4./3.)*N_eff);
  double Om_baryons=0.04899;
  double Om_vac=0.6795;
  double wde_eos=-0.9;
  double Om_k=1.-Om_matter-Om_vac-Om_radiation;
  double Om_cdm=Om_matter-Om_baryons;
  // hubble constant H/100km/s/Mpc
  double hubble=0.67;
  double Hubble=100.0;
  // fraction of baryonic energy density to matter
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.7255;
  double spectral_index = 0.96;
  double alpha_s = 0.;
  double sigma8= 0.8147;
  double RR=8.0;
  double Delta_SO = 200;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
* @namespace<Cosmo_parameters_SLICS>
* @brief Cosmological parameters usedd in SLICS simulations
*/
namespace Cosmo_parameters_SLICS
{
/**
    * @brief Total matter energy density parameter SLICS simulation
    */
  double Om_matter=0.2905;
  double N_eff=3.046;
  double Om_radiation=(1.3157e-5)*(1.+(7./8.)*std::pow(4./11.,4./3.)*N_eff);
  double Om_baryons=0.0473;
  double wde_eos=-1;
  double Om_vac=0.7095;
  double Om_k=1.-Om_matter-Om_vac-Om_radiation;
  // hubble constant H/100km/s/Mpc
  double hubble = 0.6898;
  double Hubble = 100.0;
  // fraction of baryonic energy density to matter
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.725;
  double spectral_index = 0.9667;
  double alpha_s = 0.;
  double sigma8= 0.826;
  double RR=8.0;
  double Delta_SO = 200;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
* @namespace<Cosmo_parameters_Flagship>
* @brief Cosmological parameters usedd in Flagship simulations
*/
namespace Cosmo_parameters_Flagship
{
/**
    * @brief Total matter energy density parameter Flagship Simulation
    */
  double Om_matter=0.319;
  double N_eff=3.046;
  double Om_radiation=0;//(1.3157e-5)*(1.+0.2271*N_eff);
  double Om_vac=0.681;
  double Om_baryons=0.049;
  double Om_k=0.00;//1.-Omega0-Omegavac-Omegarad;
  double Om_cdm=Om_matter-Om_baryons;
  double wde_eos=-1;
  double hubble = 0.67;
  double Hubble=100.0;
  double f_baryon = Om_baryons/Om_matter;
  double Tcmb = 2.725;
  double spectral_index = 0.96;
  double alpha_s = 0.;
  double sigma8= 0.83;
  double RR=8.0;
  double Delta_SO = 200;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
/**
* @namespace<Cosmo_parameters_abacus>
* @brief Cosmological parameters usedd in ABACUs simulations
*/
namespace Cosmo_parameters_abacus //AbacusCosmos Planck LCDM cosmology
{
  double Om_k=0.;
  double omega_bar = 0.02237;
  double omega_cdm = 0.1200;
  double Tcmb = 2.7255;
  double N_eff =  2.038;
  double omega_ndcm = 0.00064420;
  double omega_rad = 0.00000009137;
  double omega_nu = omega_rad*(7./8.)*std::pow(4./11.,4./3.)*N_eff ;
  double hubble = 0.6736;
  double Hubble =100.0;
  double wde_eos = -1.0;
  double spectral_index = 0.9649;
  double alpha_s = 0.;
  double sigma8 = 0.887952;
  double RR = 8.0;
  double Delta_SO = 200;
  /**
      * @brief Total matter energy density parameter Abacus
      */
  double Om_matter = (omega_bar+omega_cdm+omega_ndcm)/std::pow(hubble,2); // 0.315192;
  double Om_cdm = omega_cdm/std::pow(hubble,2);
  double Om_baryons = omega_bar/std::pow(hubble,2);
  double Om_radiation = (omega_rad+omega_nu)/std::pow(hubble,2);
  double Om_vac = 1.0-Om_matter-Om_radiation-Om_k;
  double f_baryon = Om_baryons/Om_matter;
}
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
#endif

