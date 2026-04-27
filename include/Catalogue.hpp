////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef CATALOG__
#define CATALOG__
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * @class<Catalogue>
 * @brief Header file for the class Catalogue::
 * @file Catalogue.h
 * @title Class in charge of reading and analysing input catalogs
 * @details The methods in this class are desinged to generate different statistics from a catalog of dark matter tracers
 * @author Andres Balaguera-Antolínez
 * @version 2.0
 * @date    2008-2024
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <math.h>
#include <fstream>
#include "bstream.h"
#include <iomanip>
#include <iostream>
#include <string>
#include <string.h>
#include <cassert>
#include <cfloat>
#include <stdlib.h>
#include <stdio.h>
#include <fftw3.h>
#include <omp.h>
#include <sstream>
#include <netinet/in.h>
#include "Params.hpp"
#include "NumericalMethods.hpp"
#include "Miscelanious.hpp"
#include "CoordinateSystem.hpp"
#include "FileOutput.hpp"
#include "ScreenOutput.hpp"
//////////////////////////////////////////////////////////
using namespace std;

//////////////////////////////////////////////////////////

#define DEFINE_ACCESSOR(name, type) \
    type& name##_at(size_t i) { return name[i]; } \
    const type& name##_at(size_t i) const { return name[i]; }

//////////////////////////////////////////////////////////
class Catalogue{


  private:


  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Object of type Class::FileOutput
   */
  FileOutput File;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Class:ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Object of type Class:ScreenOutput
   */
  Params params;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Temnplate class member used to load tracer property ionformation, in the SoA data model.
   */
   template<typename T> void load_column(std::vector<T>& target, const T* prop, int column_index, std::tuple<T, T, T>&information );

  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  const gsl_rng_type *Tn_cat;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Random number generator used for selecting a random subsample.
  */
   gsl_rng *rn_cat;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of tracers in catalogue
  */
  size_t NOBJS;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of columns in input ascii catalogue
  */
  size_t NCOLS;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief Number of tracers in catalogue
  */
  string type_of_object;

 //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Class containing the input parameters
   */
  s_params_box_mas box;
  //////////////////////////////////////////////////////////
  /**
   * @private
   * @brief  Class containing the input parameters
   */
  s_params_box_mas box_HR;
  ///////////////////////////////////////////////////////
  /**
   * @private
   * @brief
   */
  s_params_box_mas box_n;  //new, auxiliary
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   *@brief  ID (in the file) of each tracer 
   */
  vector<ULONG> galID;
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  // PROPERTIES READ FROM AN INPUT CATALOGUE:
  // These are identified from the input parameter file.

  // Coordinates
  /**
   *@brief  Coordinate 1 (X,r or z) depending on params._sys_of_coords_g/r()  
   */
  vector<real_prec> coord1;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Coordinate 2 (Y,phi or phi) params._sys_of_coords_g/r() */
  vector<real_prec> coord2;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Coordinate 1 (Z or theta, ) params._sys_of_coords_g/r() */
  vector<real_prec> coord3;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  //   Velocities:
  /**
   *@brief  Velocity component in  coord1  
   */
  vector<real_prec> vel1;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Velocity component in  coord2 
   */
  vector<real_prec> vel2;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Velocity component in  coord3 
   */
  vector<real_prec> vel3;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // Intrinsic Properties
  /**
   *@brief  Tracer Mass, idally in Ms/h  
   */
  vector<real_prec> mass;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Mass Bin in which the current tracer mass is located  
   */
  vector<int> mass_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Maximum circular velocity   
   */
  vector<real_prec> vmax;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Bin of Maximum circular velocity  
   */
  vector<int> vmax_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Vrms   
   */
  vector<real_prec> vrms;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Vrms   
   */
  vector<int> vrms_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Scale radius Rs   */
  vector<real_prec> rs;
  //////////////////////////////////////////////////////////
  /**
   *@brief Bin of Rs  
   */
  vector<int> rs_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
  * @brief  Tracer Virial radius   
  */
  vector<real_prec> rvir;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Concentration   
   */
  vector<real_prec> concentration;
  //////////////////////////////////////////////////////////
  /**
     *@brief  Tracer virial ratio T/|W|   
  */
  vector<real_prec> virial;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Bin of virial 
   */
  vector<int> virial_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Spin (Peebles)  
   */
  vector<real_prec> spin;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Spin defined by Bullock  
   */
  vector<real_prec> spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   *@brief Bin in spin
  */
  vector<int> spin_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer b_to_a ratio   
   */
  vector<real_prec> b_to_a;
  //////////////////////////////////////////////////////////
  /**
   *@brief Bin in b_to_a
  */
  vector<int> btoa_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer b_to_a ratio   
   */
  vector<real_prec> c_to_a;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Bin in c_to_a ratio   
   */
  vector<int> ctoa_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
   *@brief  mean number density at current tracer mass  
   */
  vector<real_prec> mean_number_density;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Number of substructures (if the tracer is parent halo) 
   */
  vector<int> number_sub_structures;

  //////////////////////////////////////////////////////////

 //////////////////////////////////////////////////////////
   // Galaxy like properties

  /**
   *@brief  Color (for galaxies)   
   */
  vector<real_prec> color;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Stellar mass (for galaxies)  
   */
  vector<real_prec> stellar_mass;
  //////////////////////////////////////////////////////////
  /**
   *@brief Apparent manitude  
   */
  vector<real_prec> app_mag;
  //////////////////////////////////////////////////////////
  /**
   *@brief Absolute manitude   
   */
  vector<real_prec> abs_mag;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  //Statistical properties 
  /**
   *@brief  Weight (if any)
  */
  vector<real_prec> weight1;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Weight (if any) 
   */
  vector<real_prec> weight2;
  //////////////////////////////////////////////////////////
   /**
   *@brief   Weight (if any) 
   */
  vector<real_prec> weight3;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Weight (if any) 
   */
  vector<real_prec> weight4;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Mean density (if any) evaluated at the position of the tracer 
   */
  vector<real_prec> mean_density;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // Individial halo/galaxy 
  /**
   *@brief  Effective bias computed object by object
   */
  vector<real_prec>  bias;
  //////////////////////////////////////////////////////////
  // Individial halo/galaxy 
  /**
   *@brief  Effective bias computed object by object
   */
  vector<real_prec>  bias_aux;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Quadratic bias computed object by object
   */
  vector<real_prec>  qbias;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Relative bias computed object by object
   */
  vector<real_prec>  relative_bias;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Effective bias factor computed from a box in redshift space
   */
  vector<real_prec>  bias_rs;
  //////////////////////////////////////////////////////////
  /**
   *@brief  bias squared, computed using suto power spectrum for halos  and matter
   */
  vector<real_prec>  bias_squared;
  //////////////////////////////////////////////////////////
  /**
   *@brief  bias computed from a a box as a function of multipole l (Sheth) OJO con este
   */
  vector<vector<real_prec>>  bias_multipole;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Parameter linekd to the individual bias in redshift space
   */
  vector<real_prec>  rs_factor;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  // Environmental properties
  /**
   *@brief  Number of neighbours in a sphere of radius R
   */
  vector<real_prec>  number_neigh;
  //////////////////////////////////////////////////////////
  /**
   *@brief Mean separation to local neihbours in a sphere of radius R
   */
  vector<real_prec>  mean_sep_local_neigh;

  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  //  Used for BMT Applications
  /**
   *@brief  Tracer poistion in the set (0..Nobjects) ordererd by increasing order in Vmax. Used in BMT applications.
   */
  vector<ULONG> vmax_index;
  //////////////////////////////////////////////////////////
  /**
   *@brief Bin number in the Theta (set of BMT) properties binning in which the GridID of the galaxy has fallen 
   */
  vector<ULONG> galThetaID;
  //////////////////////////////////////////////////////////
  /**
   *@brief ID to identify the multiscale level in which the tracer has acqured the main property  
   */
  vector<int> multi_scale_level;
  //////////////////////////////////////////////////////////
  /**
   *@brief  ID (in the mesh) where this tracer is located given the resolution with l1  */
  vector<ULONG> GridID_l1;
  //////////////////////////////////////////////////////////
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l2  */
  vector<ULONG> GridID_l2;
  //////////////////////////////////////////////////////////
  /**
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l3  */
  vector<ULONG> GridID_l3;
  //////////////////////////////////////////////////////////
  /**
  /**
   *@brief   ID (in the mesh) where this tracer is located given the resolution with l4  */
  vector<ULONG> GridID_l4;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Identify whether this is a random object (-1) or a real tracer from DM particle sin LPT (1). Used in BMT applications  
   */
  vector<int> identity;
  //////////////////////////////////////////////////////////
  /**
  /**
   *@brief  Identify whether this has been observed or avoided. Used in BMT applications  */
  vector<int> observed;
  //////////////////////////////////////////////////////////
  /**
  /**
   *@brief  Identify whether this has been observed or avoided. Used in BMT applications */
  vector<int> observed_secondary;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  /**
     @brief  Cosmic-web clasification of tracer 
  */
  vector<WebType> gal_cwt;
  //////////////////////////////////////////////////////////
  /**
     @brief Interpolation of the DM fretrieved from the ID of the cell where the halo is found. 
     */
  vector<real_prec> local_dm;
  //////////////////////////////////////////////////////////
  /**
     @brief Eigenvalues of DM retrieved from the ID of the cell where the halo is found. 
  */
  vector<real_prec> lambda1;
  //////////////////////////////////////////////////////////
  /**
     @brief Eigenvalues of DM retrieved from the ID of the cell where the halo is found. 
   */
  vector<real_prec> lambda2;
  //////////////////////////////////////////////////////////
  /**
     @brief Eigenvalues of DM retrieved from the ID of the cell where the halo is found. 
  */
  vector<real_prec> lambda3;
  //////////////////////////////////////////////////////////
  /**
     @brief  Measuremet of the number of tracers around a tracer in asphere of R set in parameter file. 
  */     
  vector<real_prec> local_overdensity;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Mach number   
   */
  vector<real_prec> mach_number;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Bin of tracer Mach number
  */
  vector<int> mach_number_bin_tracer;
  //////////////////////////////////////////////////////////
  /**
     @brief  Analogous to mach number, for distances to tracer: mean/var of the separations around the tracer in a sphere on a scale given in parameter file 
  */
  vector<real_prec> dach_number;
  //////////////////////////////////////////////////////////
  /**
     @brief  Analogous to mach number, for distances to tracer: mean/var of the separations around the tracer in a sphere on a scale given in parameter file 
   */
  vector<real_prec> number_of_neighbours;
  //////////////////////////////////////////////////////////
  /**
     @brief  Analogous to mach number, for distances to tracer: mean/var of the separations around the tracer in a sphere on a scale given in parameter file 
   */
  vector<int> number_of_substructures;
  //////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////
  /**
     @brief  Tidal anisotropy evaluated at the tracer position 
  */
  vector<real_prec> tidal_anisotropy;
  //////////////////////////////////////////////////////////
  /**
     @brief  Tidal anisotropy evaluated at the tracer position based on the DM tidal field 
  */
  vector<real_prec> tidal_anisotropy_dm;
  //////////////////////////////////////////////////////////
  /**
     @brief  Peak height delta/sigma computed for each tracer based on their virial mass
  */
  vector<real_prec> peak_height;
  //////////////////////////////////////////////////////////
  /**
    @brief  Spherical overdensity associated to the mass and radii of the tracer
  */
  vector<real_prec> spherical_overdensity;
  //////////////////////////////////////////////////////////
  /**
    @brief  Mass of the closest neighbour
  */
  vector<real_prec> mass_closest_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Distance to the closest neighbour
  */
  vector<real_prec> distance_closest_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Mass of the most massive neighbour
  */
  vector<real_prec> most_massive_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Distnace ot the most massive neighbour
  */
  vector<real_prec> distance_to_most_massive_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Spin of the clisest neighbour
  */
  vector<real_prec> spin_closest_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Concentration of the closet neighbour
  */
  vector<real_prec> concentration_closest_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  ID of the closest neighbour
  */
  vector<ULONG> id_closest_neighbour;
  //////////////////////////////////////////////////////////
  /**
    @brief  Healpixel where the galaxy is located 
  */
  vector<ULONG> galPIXEL;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Redshift (Z or theta, ) params._sys_of_coords_g/r() 
   */
  vector<real_prec> cosmological_redshift;
  //////////////////////////////////////////////////////////
  /**
   *@brief   Redshift (Z or theta, ) params._sys_of_coords_g/r() 
   */
  vector<real_prec> redshift;
 //////////////////////////////////////////////////////////
  /**
   *@brief  Mass of tracer reassigned,  used for mocks and bmt 
   */
  vector<real_prec> mass_assigned;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Mass of reference tracer, used for mocks  
   */
  vector<real_prec> mass_parent;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Vmax assigend. Used in BMT applications.
   */
  vector<real_prec> vmax_assigned;
  //////////////////////////////////////////////////////////
  /**
   *@brief  Tracer Vmax from parent halo. Used in BMT applications   
  */
  vector<real_prec> vmax_parent;

  //////////////////////////////////////////////////////////
  /**
   *@brief ID (in the mesh) where this tracer is located given the nomial resolution () 
   */
  vector<ULONG> GridID;
  //////////////////////////////////////////////////////////
  /**
   *@brief ID (in the mesh) where this tracer is located given the nomial resolution () 
   */
  vector<ULONG> GridID_n;

  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////


public:
  //////////////////////////////////////////////////////////
  /**
   * @brief Constructor
   */
  Catalogue(){}

  //////////////////////////////////////////////////////////
  /**
   * @brief  Constructor (overloaded)
   * @details
   * @arg Object of type Params with the parameters read from the input parameter file
   */
  Catalogue(Params _params, string _type_of_object):params (_params), type_of_object(_type_of_object)
  {
    // This constructor can be used when the type is declared
    // and used directly. Otherwise, if a class has an object of
    // thi type, we define the object calling the default
    // constructor and calling set_params_catalog to assign values to this class, specifçically, to a structure

    time_t time_bam;
    time(&time_bam);
    this->So.initial_time=time_bam;

    gsl_rng_env_setup();
    gsl_rng_default_seed=1015;
    this->Tn_cat =  gsl_rng_mt19937 ;//gsl_rng_default;
    rn_cat = gsl_rng_alloc (this->Tn_cat);


  }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Default constructor
   */
  ~Catalogue(){}

  //////////////////////////////////////////////////////////
  /**
   * @brief Writes a catalogue produced in BMT  with N_PROP columns (see def.h) to a binary file
   * @author ABA 
   */
  void write_catalog_bin(string output_file);

  //////////////////////////////////////////////////////////
  /**
   * @brief Writes catalog with N_PROP columns (see def.h) to an ascii file
   */
  void write_catalog(string output_file);
  //////////////////////////////////////////////////////////
  /**
  * @brief Reads particle catalog in ascii format. 
  * @details This method reads the ascii file containing the catalog (with columns specified in the parameter file) and, contrary to the read_catalog approach, it associates the catalog to a single structure of vectors.
  */
   void read_catalog_new(string);
  //////////////////////////////////////////////////////////
  /**
  * @brief Print Catalogue
  * @details Used to prints out of catalogue with assigned
  */
  void print_catalogue(string file, bool bias);

  //////////////////////////////////////////////////////////
  /**
   * @brief Read input catalog in binary format, with x, y, z (and velocities if requested) in three separate files of lenght n
   */
#ifdef _USE_VELOCITIES_
  void read_catalog_bin(ULONG n, string fx, string fy, string fz,string fvx, string fvy, string fvz); // reads binary from fastpm
#else
  void read_catalog_bin(); // reads binary from fastpm
  void read_catalog_bin_tng(); // reads binary from fastpm
#endif
  //////////////////////////////////////////////////////////
  /**
   * @brief Function to set parameters
   */
  void set_params(Params &new_params){
    this->params=new_params;
  }

  //////////////////////////////////////////////////////////
  // Define Accesors to the elements 
  DEFINE_ACCESSOR(coord1, real_prec);
  DEFINE_ACCESSOR(coord2, real_prec);
  DEFINE_ACCESSOR(coord3, real_prec);

  DEFINE_ACCESSOR(vel1, real_prec);
  DEFINE_ACCESSOR(vel2, real_prec);
  DEFINE_ACCESSOR(vel3, real_prec);

  DEFINE_ACCESSOR(mass, real_prec);
  DEFINE_ACCESSOR(mass_bin_tracer, int);

  DEFINE_ACCESSOR(vmax, real_prec);
  DEFINE_ACCESSOR(vmax_bin_tracer, int);

  DEFINE_ACCESSOR(vrms, real_prec);
  DEFINE_ACCESSOR(vrms_bin_tracer, int);

  DEFINE_ACCESSOR(rs, real_prec);
  DEFINE_ACCESSOR(rs_bin_tracer, int);

  DEFINE_ACCESSOR(rvir, real_prec);
  DEFINE_ACCESSOR(concentration, real_prec);

  DEFINE_ACCESSOR(virial, real_prec);
  DEFINE_ACCESSOR(virial_bin_tracer, int);

  DEFINE_ACCESSOR(spin, real_prec);
  DEFINE_ACCESSOR(spin_bullock, real_prec);
  DEFINE_ACCESSOR(spin_bin_tracer, int);

  DEFINE_ACCESSOR(b_to_a, real_prec);
  DEFINE_ACCESSOR(btoa_bin_tracer, int);

  DEFINE_ACCESSOR(c_to_a, real_prec);
  DEFINE_ACCESSOR(ctoa_bin_tracer, int);

  DEFINE_ACCESSOR(mean_number_density, real_prec);
  DEFINE_ACCESSOR(mean_density, real_prec);
  DEFINE_ACCESSOR(number_of_substructures, int);

  DEFINE_ACCESSOR(color, real_prec);
  DEFINE_ACCESSOR(stellar_mass, real_prec);
  DEFINE_ACCESSOR(app_mag, real_prec);
  DEFINE_ACCESSOR(abs_mag, real_prec);

  DEFINE_ACCESSOR(weight1, real_prec);
  DEFINE_ACCESSOR(weight2, real_prec);
  DEFINE_ACCESSOR(weight3, real_prec);
  DEFINE_ACCESSOR(weight4, real_prec);

  DEFINE_ACCESSOR(bias, real_prec);
  DEFINE_ACCESSOR(bias_aux, real_prec);
  DEFINE_ACCESSOR(qbias, real_prec);
  DEFINE_ACCESSOR(relative_bias, real_prec);
  DEFINE_ACCESSOR(bias_rs, real_prec);
  DEFINE_ACCESSOR(bias_squared, real_prec);
  DEFINE_ACCESSOR(rs_factor, real_prec);

  DEFINE_ACCESSOR(number_neigh, real_prec);
  DEFINE_ACCESSOR(mean_sep_local_neigh, real_prec);

  DEFINE_ACCESSOR(vmax_index, ULONG);
  DEFINE_ACCESSOR(galThetaID, ULONG);
  DEFINE_ACCESSOR(multi_scale_level, int);
  DEFINE_ACCESSOR(gal_cwt, WebType);

  DEFINE_ACCESSOR(GridID, ULONG);
  DEFINE_ACCESSOR(GridID_n, ULONG);

  DEFINE_ACCESSOR(GridID_l1, ULONG);
  DEFINE_ACCESSOR(GridID_l2, ULONG);
  DEFINE_ACCESSOR(GridID_l3, ULONG);
  DEFINE_ACCESSOR(GridID_l4, ULONG);

  DEFINE_ACCESSOR(identity, int);
  DEFINE_ACCESSOR(observed, int);
  DEFINE_ACCESSOR(observed_secondary, int);

  DEFINE_ACCESSOR(local_dm, real_prec);
  DEFINE_ACCESSOR(lambda1, real_prec);
  DEFINE_ACCESSOR(lambda2, real_prec);
  DEFINE_ACCESSOR(lambda3, real_prec);
  DEFINE_ACCESSOR(local_overdensity, real_prec);
  DEFINE_ACCESSOR(mach_number, real_prec);
  DEFINE_ACCESSOR(dach_number, real_prec);
  DEFINE_ACCESSOR(mach_number_bin_tracer, int);
  DEFINE_ACCESSOR(number_of_neighbours, real_prec);

  DEFINE_ACCESSOR(tidal_anisotropy, real_prec);
  DEFINE_ACCESSOR(tidal_anisotropy_dm, real_prec);
  DEFINE_ACCESSOR(peak_height, real_prec);
  DEFINE_ACCESSOR(spherical_overdensity, real_prec);
  DEFINE_ACCESSOR(mass_closest_neighbour, real_prec);
  DEFINE_ACCESSOR(distance_closest_neighbour, real_prec);
  DEFINE_ACCESSOR(most_massive_neighbour, real_prec);
  DEFINE_ACCESSOR(distance_to_most_massive_neighbour, real_prec);
  DEFINE_ACCESSOR(spin_closest_neighbour, real_prec);
  DEFINE_ACCESSOR(concentration_closest_neighbour, real_prec);
  DEFINE_ACCESSOR(id_closest_neighbour, ULONG);
  DEFINE_ACCESSOR(galPIXEL, ULONG);

  DEFINE_ACCESSOR(cosmological_redshift, real_prec);
  DEFINE_ACCESSOR(redshift, real_prec);

  DEFINE_ACCESSOR(mass_assigned, real_prec);
  DEFINE_ACCESSOR(mass_parent, real_prec);
  DEFINE_ACCESSOR(vmax_assigned, real_prec);
  DEFINE_ACCESSOR(vmax_parent, real_prec);

  //////////////////////////////////////////////////////////
  /**
   * @brief Resize 
   */
  real_prec bias_multipole_at(size_t i, int mult){
    return this->bias_multipole[i][mult];
  }
  void set_bias_multipole(real_prec abc, size_t i, int mult){
    this->bias_multipole[i][mult]=abc;
  }
  void resize_bias_multipole(size_t N, int Nm){
    bias_multipole.resize(N);
    for(ULONG itr=0;itr<N;++itr)
      bias_multipole[itr].resize(Nm+1,0);
  }
  int bias_multipole_size(){
     return bias_multipole[0].size();
  }

  //////////////////////////////////////////////////////////
  /**
   * @brief Resize positions
   */
  void resize_positions(std::size_t N)
  {
    coord1.resize(N);
    coord2.resize(N);
    coord3.resize(N);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Resize velocities
   */
  void resize_velocities(std::size_t N)
  {
    vel1.resize(N);
    vel2.resize(N);
    vel3.resize(N);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief Resize intrinsic properties
   */
  void resize_halo_properties(std::size_t N)
  {
    mass.resize(N);
    vmax.resize(N);
    vrms.resize(N);
    rs.resize(N);
    concentration.resize(N);
    spin.resize(N);
    spin_bullock.resize(N);
    bias.resize(N);
    bias_aux.resize(N);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void resize_GridIDs(std::size_t N)
  {
    this->GridID.resize(N);
    this->GridID_n.resize(N);
  }
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */

  void resize_galaxy_properties(std::size_t N)
  {
  }
 // add more

  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  size_t _NOBJS(){return this->NOBJS;}
  void set_NOBJS(size_t N){this->NOBJS=N;}

  //////////////////////////////////////////////////////////
  /**
   * @brief  Clean all memmory
   */
 void clear_mem()
    {
     auto clean_vec = [](auto& v)
      {
        v.clear();
        v.shrink_to_fit();
      };
     clean_vec(coord1);
     clean_vec(coord2);
     clean_vec(coord3);
     clean_vec(cosmological_redshift);
     clean_vec(redshift);
     clean_vec(vel1);
     clean_vec(vel2);
     clean_vec(vel3);

     clean_vec(mass);
     clean_vec(mass_bin_tracer);
     clean_vec(mass_assigned);
     clean_vec(mass_parent);

     clean_vec(vmax);
     clean_vec(vmax_bin_tracer);
     clean_vec(vrms);
     clean_vec(vrms_bin_tracer);
     clean_vec(mach_number);
     clean_vec(mach_number_bin_tracer);

     clean_vec(vmax_assigned);
     clean_vec(vmax_parent);
     clean_vec(vmax_index);

     clean_vec(bias);
     clean_vec(bias_aux);
     clean_vec(qbias);
     clean_vec(relative_bias);
     clean_vec(bias_rs);
     clean_vec(bias_squared);
     clean_vec(bias_multipole);

     clean_vec(rs_factor);
     clean_vec(number_neigh);
     clean_vec(mean_sep_local_neigh);

     clean_vec(rs);
     clean_vec(rvir);
     clean_vec(concentration);
     clean_vec(rs_bin_tracer);

     clean_vec(virial);
     clean_vec(virial_bin_tracer);

     clean_vec(spin);
     clean_vec(spin_bullock);
     clean_vec(spin_bin_tracer);

     clean_vec(b_to_a);
     clean_vec(btoa_bin_tracer);

     clean_vec(c_to_a);
     clean_vec(ctoa_bin_tracer);

     clean_vec(mean_number_density);

     clean_vec(color);
     clean_vec(stellar_mass);
     clean_vec(app_mag);
     clean_vec(abs_mag);

     clean_vec(GridID);
     clean_vec(GridID_n);
     clean_vec(galThetaID);
     clean_vec(galID);
     clean_vec(multi_scale_level);

     clean_vec(GridID_l1);
     clean_vec(GridID_l2);
     clean_vec(GridID_l3);
     clean_vec(GridID_l4);

     clean_vec(number_of_substructures);
     clean_vec(identity);
     clean_vec(observed);
     clean_vec(observed_secondary);

     clean_vec(weight1);
     clean_vec(weight2);
     clean_vec(weight3);
     clean_vec(weight4);

     clean_vec(mean_density);
     clean_vec(gal_cwt);
     clean_vec(local_dm);
     clean_vec(lambda1);
     clean_vec(lambda2);
     clean_vec(lambda3);
     clean_vec(local_overdensity);
     clean_vec(dach_number);
     clean_vec(number_of_neighbours);
     clean_vec(tidal_anisotropy);
     clean_vec(tidal_anisotropy_dm);
     clean_vec(peak_height);
     clean_vec(spherical_overdensity);

     clean_vec(mass_closest_neighbour);
     clean_vec(distance_closest_neighbour);
     clean_vec(most_massive_neighbour);
     clean_vec(distance_to_most_massive_neighbour);

     clean_vec(spin_closest_neighbour);
     clean_vec(concentration_closest_neighbour);

     clean_vec(id_closest_neighbour);
     clean_vec(galPIXEL);
   };


  //////////////////////////////////////////////////////////
   /**
    * @brief Set the grid ID for each tracer
    */
   void set_GridID(size_t ID, size_t i){
    this->GridID[i]=ID;
   }
   void push_GridID(size_t ID){
    this->GridID.push_back(ID);
   }
   void resize_GridID(size_t N){
    this->GridID.resize(N);
   }
 
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
    */
   void set_GridID_n(size_t ID, size_t i){
    this->GridID_n[i]=ID;
   }
   void push_GridID_n(size_t ID){
    this->GridID_n.push_back(ID);
   }
   void resize_GridID_n(size_t N){
    this->GridID_n.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_galPIXEL(size_t ID, size_t i){
     this->galPIXEL[i]=ID;
    }
   void push_galPIXEL(size_t ID){
    this->galPIXEL.push_back(ID);
   }
   void resize_galPIXEL(size_t N){
    this->galPIXEL.resize(N);
   }
 
    //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_mean_density(real_prec abc, size_t i){
    this->mean_density[i]=abc;
   }
   void push_mean_density(real_prec abc){
    this->mean_density.push_back(abc);
   }
   void resize_mean_density(size_t N){
    this->mean_density.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_redshift(real_prec abc, size_t i){
    this->redshift[i]=abc;
   }
   void push_redshift(real_prec abc){
    this->redshift.push_back(abc);
   }
   void resize_redshift(size_t N){
    this->redshift.push_back(N);
   }
   //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_cosmological_redshift(real_prec abc, size_t i){
    this->cosmological_redshift[i]=abc;
   }
   void push_cosmological_redshift(real_prec abc){
    this->cosmological_redshift.push_back(abc);
   }
   void resize_cosmological_redshift(size_t N){
    this->cosmological_redshift.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_coord1(real_prec abc, size_t i){
    this->coord1[i]=abc;
   }
   void push_coord1(real_prec abc){
    this->coord1.push_back(abc);
   }
   void resize_coord1(size_t N){
    this->coord1.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vel1(real_prec abc, size_t i){
    this->vel1[i]=abc;
   }
   void push_vel1(real_prec abc){
    this->vel1.push_back(abc);
   }
   void resize_vel1(size_t N){
    this->vel1.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_coord2(real_prec abc, size_t i){
    this->coord2[i]=abc;
   }
   void push_coord2(real_prec abc){
    this->coord2.push_back(abc);
   }
   void resize_coord2(size_t N){
    this->coord1.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vel2(real_prec abc, size_t i){
    this->vel2[i]=abc;
   }
   void push_vel2(real_prec abc){
    this->vel2.push_back(abc);
   }
   void resize_vel2(size_t N){
    this->vel2.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_coord3(real_prec abc, size_t i){
    this->coord3[i]=abc;
   }
   void push_coord3(real_prec abc){
    this->coord3.push_back(abc);
   }
   void resize_coord3(size_t N){
    this->coord1.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vel3(real_prec abc, size_t i){
    this->vel3[i]=abc;
   }
   void push_vel3(real_prec abc){
    this->vel3.push_back(abc);
   }
   void resize_vel3(size_t N){
    this->vel3.push_back(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_observed(int abc, size_t i){
    this->observed[i]=abc;
   }
   void push_observed(int abc){
    this->observed.push_back(abc);
   }
   void resize_observed(size_t N){
    this->observed.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_number_of_substructures(int abc, size_t i){
    this->number_of_substructures[i]=abc;
   }
   void push_number_of_substructures(int abc){
    this->number_of_substructures.push_back(abc);
   }
   void resize_number_of_substructures(size_t N){
    this->number_of_substructures.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_gal_cwt(WebType abc, size_t i){
    this->gal_cwt[i]=abc;
   }
   void push_gal_cwt(WebType abc){
    this->gal_cwt.push_back(abc);
   }
   void resize_gal_cwt(size_t N){
    this->gal_cwt.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_tidal_anisotropy(real_prec abc, size_t i){
    this->tidal_anisotropy[i]=abc;
   }
   void push_tidal_anisotropy(real_prec abc){
    this->tidal_anisotropy.push_back(abc);
   }
   void resize_tidal_anisotropy(size_t N){
    this->tidal_anisotropy.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_tidal_anisotropy_dm(real_prec abc, size_t i){
    this->tidal_anisotropy_dm[i]=abc;
   }
   void push_tidal_anisotropy_dm(real_prec abc){
    this->tidal_anisotropy_dm.push_back(abc);
   }
   void resize_tidal_anisotropy_dm(size_t N){
    this->tidal_anisotropy_dm.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_lambda1(real_prec abc, size_t i){
    this->lambda1[i]=abc;
   }
   void push_lambda1(real_prec abc){
    this->lambda1.push_back(abc);
   }
   void resize_lambda1(size_t N){
    this->lambda1.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_lambda2(real_prec abc, size_t i){
    this->lambda2[i]=abc;
   }
   void push_lambda2(real_prec abc){
    this->lambda2.push_back(abc);
   }
   void resize_lambda2(size_t N){
    this->lambda2.resize(N);
   }
     //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_lambda3(real_prec abc, size_t i){
    this->lambda3[i]=abc;
   }
   void push_lambda3(real_prec abc){
    this->lambda3.push_back(abc);
   }
   void resize_lambda3(size_t N){
    this->lambda3.resize(N);
   }

  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_local_dm(real_prec abc, size_t i){
    this->local_dm[i]=abc;
   }
   void push_local_dm(real_prec abc){
    this->local_dm.push_back(abc);
   }
   void resize_local_dm(size_t N){
    this->local_dm.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_multi_scale_level(int abc, size_t i){
    this->multi_scale_level[i]=abc;
   }
   void resize_multi_scale_level(size_t N){
    this->multi_scale_level.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_concentration(real_prec abc, size_t i){
    this->concentration[i]=abc;
   }
   void push_concentration(real_prec abc){
    this->concentration.push_back(abc);
   }
   void resize_concentration(size_t N){
    this->concentration.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_concentration_closest_neighbour(real_prec abc, size_t i){
    this->concentration_closest_neighbour[i]=abc;
   }
   void push_concentration_closest_neighbour(real_prec abc){
    this->concentration_closest_neighbour.push_back(abc);
   }
   void resize_concentration_closest_neighbour(size_t N){
    this->concentration_closest_neighbour.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_spin_closest_neighbour(real_prec abc, size_t i){
    this->spin_closest_neighbour[i]=abc;
   }
   void push_spin_closest_neighbour(real_prec abc){
    this->spin_closest_neighbour.push_back(abc);
   }
   void resize_spin_closest_neighbour(size_t N){
    this->spin_closest_neighbour.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_galThetaID(ULONG abc, size_t i){
    this->galThetaID[i]=abc;
   }
   void push_galThetaID(ULONG abc){
    this->galThetaID.push_back(abc);
   }
   void resize_galThetaID(size_t N){
    this->galThetaID.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_mass_closest_neighbour(real_prec abc, size_t i){
    this->mass_closest_neighbour[i]=abc;
   }
   void push_mass_closest_neighbour(real_prec abc){
    this->mass_closest_neighbour.push_back(abc);
   }
   void resize_mass_closest_neighbour(size_t N){
    this->mass_closest_neighbour.resize(N);
   }
  //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_distance_closest_neighbour(real_prec abc, size_t i){
    this->distance_closest_neighbour[i]=abc;
   }
   void push_distance_closest_neighbour(real_prec abc){
    this->distance_closest_neighbour.push_back(abc);
   }
   void resize_distance_closest_neighbour(size_t N){
    this->distance_closest_neighbour.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_spin(real_prec abc, size_t i){
    this->spin[i]=abc;
   }
   void push_spin(real_prec abc){
    this->spin.push_back(abc);
   }
   void resize_spin(size_t N){
    this->spin.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_spin_bullock(real_prec abc, size_t i){
    this->spin_bullock[i]=abc;
   }
   void push_spin_bullock(real_prec abc){
    this->spin_bullock.push_back(abc);
   }
   void resize_spin_bullock(size_t N){
    this->spin_bullock.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_rs(real_prec abc, size_t i){
    this->rs[i]=abc;
   }
   void push_rs(real_prec abc){
    this->rs.push_back(abc);
   }
   void resize_rs(size_t N){
    this->rs.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_b_to_a(real_prec abc, size_t i){
    this->b_to_a[i]=abc;
   }
   void push_b_to_a(real_prec abc){
    this->b_to_a.push_back(abc);
   }
   void resize_b_to_a(size_t N){
    this->b_to_a.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_c_to_a(real_prec abc, size_t i){
    this->c_to_a[i]=abc;
   }
   void push_c_to_a(real_prec abc){
    this->c_to_a.push_back(abc);
   }
   void resize_c_to_a(size_t N){
    this->c_to_a.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_qbias(real_prec abc, size_t i){
    this->qbias[i]=abc;
   }
   void push_qbias(real_prec abc){
    this->qbias.push_back(abc);
   }
   void resize_qbias(size_t N){
    this->qbias.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_bias_squared(real_prec abc, size_t i){
    this->bias_squared[i]=abc;
   }
   void push_bias_squared(real_prec abc){
    this->bias_squared.push_back(abc);
   }
   void resize_bias_squared(size_t N){
    this->bias_squared.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_bias(real_prec abc, size_t i){
    this->bias[i]=abc;
   }
   void push_bias(real_prec abc){
    this->bias.push_back(abc);
   }
   void resize_bias(size_t N){
    this->bias.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_bias_rs(real_prec abc, size_t i){
    this->bias_rs[i]=abc;
   }
   void push_bias_rs(real_prec abc){
    this->bias_rs.push_back(abc);
   }
   void resize_bias_rs(size_t N){
    this->bias_rs.resize(N);
   }

 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_rs_factor(real_prec abc, size_t i){
    this->rs_factor[i]=abc;
   }
   void push_rs_factor(real_prec abc){
    this->rs_factor.push_back(abc);
   }
   void resize_rs_factor(size_t N){
    this->rs_factor.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_virial(real_prec abc, size_t i){
    this->virial[i]=abc;
   }
   void push_virial(real_prec abc){
    this->virial.push_back(abc);
   }
   void resize_virial(size_t N){
    this->virial.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_rvir(real_prec abc, size_t i){
    this->rvir[i]=abc;
   }
   void push_rvir(real_prec abc){
    this->rvir.push_back(abc);
   }
   void resize_rvir(size_t N){
    this->rvir.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vrms(real_prec abc, size_t i){
    this->vrms[i]=abc;
   }
   void push_vrms(real_prec abc){
    this->vrms.push_back(abc);
   }
   void resize_vrms(size_t N){
    this->vrms.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vmax(real_prec abc, size_t i){
    this->vmax[i]=abc;
   }
   void push_vmax(real_prec abc){
    this->vmax.push_back(abc);
   }
   void resize_vmax(size_t N){
    this->vmax.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_vmax_parent(real_prec abc, size_t i){
    this->vmax_parent[i]=abc;
   }
   void push_vmax_parent(real_prec abc){
    this->vmax_parent.push_back(abc);
   }
   void resize_vmax_parent(size_t N){
    this->vmax_parent.resize(N);
   }


 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_bias_aux(real_prec abc, size_t i){
    this->bias[i]=abc;
   }
   void push_bias_aux(real_prec abc){
    this->bias_aux.push_back(abc);
   }
   void resize_bias_aux(size_t N){
    this->bias_aux.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_relative_bias(real_prec abc, size_t i){
    this->relative_bias[i]=abc;
   }
   void push_relative_bias(real_prec abc){
    this->relative_bias.push_back(abc);
   }
   void resize_relative_bias(size_t N){
    this->relative_bias.resize(N);
   }

 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_number_of_neighbours(real_prec abc, size_t i){
    this->number_of_neighbours[i]=abc;
   }
   void push_number_of_neighbours(real_prec abc){
    this->number_of_neighbours.push_back(abc);
   }
   void resize_number_of_neighbours(size_t N){
    this->number_of_neighbours.push_back(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_mass(real_prec abc, size_t i){
    this->mass[i]=abc;
   }
   void push_mass(real_prec abc){
    this->mass.push_back(abc);
   }
   void resize_mass(size_t N){
    this->mass.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_weight1(real_prec abc, size_t i){
    this->weight1[i]=abc;
   }
   void push_weight1(real_prec abc){
    this->weight1.push_back(abc);
   }
   void resize_weight1(size_t N){
    this->weight1.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_mass_parent(real_prec abc, size_t i){
    this->mass_parent[i]=abc;
   }
   void push_mass_parent(real_prec abc){
    this->mass_parent.push_back(abc);
   }
   void resize_mass_parent(size_t N){
    this->mass_parent.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_mach_number(real_prec abc, size_t i){
    this->mach_number[i]=abc;
   }
   void push_mach_number(real_prec abc){
    this->mach_number.push_back(abc);
   }
   void resize_mach_number(size_t N){
    this->mach_number.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_dach_number(real_prec abc, size_t i){
    this->dach_number[i]=abc;
   }
   void push_dach_number(real_prec abc){
    this->dach_number.push_back(abc);
   }
   void resize_dach_number(size_t N){
    this->dach_number.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_local_overdensity(real_prec abc, size_t i){
    this->local_overdensity[i]=abc;
   }
   void push_local_overdensity(real_prec abc){
    this->local_overdensity.push_back(abc);
   }
   void resize_local_overdensity(size_t N){
    this->local_overdensity.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_peak_height(real_prec abc, size_t i){
    this->peak_height[i]=abc;
   }
   void push_peak_height(real_prec abc){
    this->local_overdensity.push_back(abc);
   }
   void resize_peak_height(size_t N){
    this->local_overdensity.resize(N);
   }

 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_distance_to_most_massive_neighbour(real_prec abc, size_t i){
    this->distance_to_most_massive_neighbour[i]=abc;
   }
   void push_distance_to_most_massive_neighbour(real_prec abc){
    this->distance_to_most_massive_neighbour.push_back(abc);
   }
   void resize_distance_to_most_massive_neighbour(size_t N){
    this->distance_to_most_massive_neighbour.resize(N);
   }
 //////////////////////////////////////////////////////////
   /**
   * @public
   * @brief Set the grid ID for each tracer
   */
   void set_most_massive_neighbour(real_prec abc, size_t i){
    this->most_massive_neighbour[i]=abc;
   }
   void push_most_massive_neighbour(real_prec abc){
    this->most_massive_neighbour.push_back(abc);
   }
   void resize_most_massive_neighbour(size_t N){
    this->most_massive_neighbour.resize(N);
   }


   //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_coord1;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_coord2;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_coord3;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_vel1;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_vel2;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_vel3;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_mass;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_vmax;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_vrms;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_mean_density;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_rs;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_spin;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_spin_bullock;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_b_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_c_to_a;
  //////////////////////////////////////////////////////////
  /**
   * @public
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_stellar_mass;
  //////////////////////////////////////////////////////////
  /**
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_color;
  //////////////////////////////////////////////////////////
  /**
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_app_mag;
  //////////////////////////////////////////////////////////
  /**
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_abs_mag;
  //////////////////////////////////////////////////////////
  /**
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_redshift;
  //////////////////////////////////////////////////////////
  /**
   * @brief Tuple for min, max and mean of tracer mass
   */
  std::tuple<real_prec, real_prec, real_prec> info_weight1;
  //////////////////////////////////////////////////////////
 /**
   * @brief Retrieve the number of columns read from the ascii input tracer catalogue
   */
  string _type_of_object(){return this->type_of_object;}
  void set_type_of_object(string nt){this->type_of_object=nt;}
  //////////////////////////////////////////////////////////
 /**
   * @brief Retrieve the number of columns read from the ascii input tracer catalogue
   */
  size_t _NCOLS(){return this->NCOLS;}
  //////////////////////////////////////////////////////////
  /**
   * @brief Access the params class
   */
   Params& _params(){
    return this->params;
   }


};


#endif
