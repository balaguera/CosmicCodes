////////////////////////////////////////////////////////////////////////////
/**
* @namespace<Constants>
* @file Constants.h
* @title Constants
* @author Andres Balaguera-Antolínez
*/
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifndef __CONSTANTS__
#define __CONSTANTS__

namespace Constants
{
/**
    @brief Speed of light in km/s
*/
  const double speed_light            = 299792.458; // km/s
  /**
      @brief Gravitational constant in units units (km)³ /kg/s²
  */
  const double Gravitational_constant = 6.673e-20;
  /**
      @brief Solar mass in Kg
  */
  const double Solar_mass             = 1.988e30; /*in kg*/
  /**
      @brief Proton Mass in kg
  */
  const double proton_mass            = 1.6726e-27;
  /**
      @brief Bolzman constant in (km)² kg /s²/K
  */
  const double Boltzmann_constant     = 1.3806503e-28;
  // Conversion factors
  /**
      @brief One kev in Kg km²/s²
  */
  const double kev                    = 1.602e-24; /*One Kiloelectronvolt in */
  /**
      @brief From solar masses to geometrcal units
  */
  const double msol_to_mpc       = 4.78e-20;      /*factor que pasa de masas solares a mpc en unidades geometricas*/
  /**
      @brief 1 Mpc in km
  */
  const double Mpc_to_km         = 3.0856e19;     /*1 megaparsec in Km*/
  /**
      @brief Years to seconds
  */
  const double years_to_sec      = 31557600.0;    /*One year in seconds*/
  const double M_reference       = 1e12;
  const double z_mean_sdss       = 0.3;
  const double mean_mol_weight    = 0.59;
  const double sfr_efficiency    = 0.18;          /*star formation rate efficiency*/
  const double min_circular_vel  = 100.0;         /*Minimum circular velocity in km/s for SFR recipie*/
}

#endif
