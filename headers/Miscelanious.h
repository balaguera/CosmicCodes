#ifndef __MISC__
#define __MISC__
//////////////////////////////////////////////////////////
/**
 * @brief Header file for the functions in Miscelanious.cpp
 * @file Miscelanious.h
 * @details Several functions used in the analysis of dark matter tracers and fields.
 * @author Andres Balaguera-Antolínez (ABA)
 * @version 1.0
 * @date  2012-2023
 */
//////////////////////////////////////////////////////////
#include "NumericalMethods.h"
#include "Params.h"
#include "FftwFunctions.h"
#include "ScreenOutput.h"
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP mass assignment scheme.
 * @param x x-coodinate of tracer in Mpc/h
 * @param y y-coodinate of tracer in Mpc/h
 * @param z y-coodinate of tracer in Mpc/h
 * @param weight Property to use as weight
 * @returns field An update of the density field
 * @details  This method uses a transformation of the coordinates described as follows.
 * @note Transformation of coordinates to embed the sample in a cubical box of size \f$V=L^{3}\f$ is
 \f$ X_i - > X_i - X_{Offset} + L/2 \f$. The second term (\f$X_{Offset}\f$) is \f$ X_{Offset}= (X_{max} + X_{min})/2\f$
 places the mid point as the origin of coordinates. The third one places the center of the sample in the center of the cube
 For a simulation with the point origin (0,0,0) at the lower corner of the box, \f$X_{Offset}=Y_{Offset}=Z_{Offset}=L/2\f$
 no real translation is done. If the box has its center at the origin of coordinates, Xoffset=0
 and then transformation moves the center of the box to \f$(L/2,L/2,L/2)\f$
 If Lbox is chosen as xmax-xmin (since it is the largest side) then
 \f$  X_{Offset}- L/2 = xmin\f$ and we end with \f$x=x-x_{min}\f$
 but this is not the case for the other components, e,g.
  \f$  Y_{Offset}  -L/2 = y_{min} -(Lx-Ly)/2 \f$
 Therefore, if the sides of the box are not equal (but still we aim at embedding the survey inside a cube),
 we need to transform with \f$x - > x - (X_{Offset}-L/2)\f$ instead of subtracting the minimum.
 This embeds the sample in a box in which, if Lbox is from x, \f$x_{min} =0\f$ \f$y_{min}=0\f$, \f$z_{min}=0\f$ but the
 actual boundary of the survey, say, in the y direction, is not \fy_{min}\f, but at \f$y=(Lx-Ly)/2\f$.
 Same applies for CIC, TSC and PCS interpolation schemes.
* @author ABA
*/
void grid_assignment_NGP(Params *params, real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP mass assignment scheme.
 * @param x x-coodinate of tracer in Mpc/h
 * @param y y-coodinate of tracer in Mpc/h
 * @param z y-coodinate of tracer in Mpc/h
 * @param weight Property to use as weight
 * @returns field An update of the density field
*/
void grid_assignment_CIC(Params *params, real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP mass assignment scheme.
 * @param x x-coodinate of tracer in Mpc/h
 * @param y y-coodinate of tracer in Mpc/h
 * @param z y-coodinate of tracer in Mpc/h
 * @param weight Property to use as weight
 * @returns field An update of the density field
*/
void grid_assignment_TSC(Params *params, real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP mass assignment scheme.
 * @param x x-coodinate of tracer in Mpc/h
 * @param y y-coodinate of tracer in Mpc/h
 * @param z y-coodinate of tracer in Mpc/h
 * @param weight Property to use as weight
 * @returns field An update of the density field
*/
void grid_assignment_PCS(Params *params, real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP
 * @param params: structure of type s_params_box_mas
 * @param Halo: vector of structures of type s_Halo, with the coordinates
 * @param weightmass: string denoting the property to use as weight.
 * @details The size of the mesh is characterized by params.Nft The size of the volume is defined by params.Lbox.
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_NGP
 * @author ABA
 */
void getDensity_NGP(s_params_box_mas *params, vector<s_Halo>&Halo, vector<real_prec>&delta, string weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using NGP
 * @param params: N1: number of cells in x-direction
 * @param params: N2: number of cells in y-direction
 * @param params: N3: number of cells in z-direction
 * @param params: L1: Lenght of cuboid in x-direction
 * @param params: L2: Lenght of cuboid in y-direction
 * @param params: L3: Lenght of cuboid in z-direction
 * @param params: d1: Resolution in x-direction
 * @param params: d2: Resolution in y-direction
 * @param params: d3: Resolution in z-direction
 * @param params: min1: Minumum of x-coord
 * @param params: min2: Minumum of y-coord
 * @param params: min3: Minumum of z-coord
 * @param params: xp: vector with x-coordinates
 * @param params: yp: vector with y-coordinates
 * @param params: zp: vector with z-coordinates
 * @param params: Parmass: vector with weights
 * @param weightmass: bool asking to weight using the entries by Parmass
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_NGP
 * @author ABA & FSK.
 */
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec>&Par_mass, vector<real_prec>&delta, bool weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using CIC
 * @param params: N1: number of cells in x-direction
 * @param params: N2: number of cells in y-direction
 * @param params: N3: number of cells in z-direction
 * @param params: L1: Lenght of cuboid in x-direction
 * @param params: L2: Lenght of cuboid in y-direction
 * @param params: L3: Lenght of cuboid in z-direction
 * @param params: d1: Resolution in x-direction
 * @param params: d2: Resolution in y-direction
 * @param params: d3: Resolution in z-direction
 * @param params: min1: Minumum of x-coord
 * @param params: min2: Minumum of y-coord
 * @param params: min3: Minumum of z-coord
 * @param params: xp: vector with x-coordinates
 * @param params: yp: vector with y-coordinates
 * @param params: zp: vector with z-coordinates
 * @param params: Parmass: vector with weights
 * @param weightmass: bool asking to weight using the entries by Parmass
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_CIC
 * @author ABA & FSK.
 */
void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>& xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using CIC
 * @param params: structure of type s_params_box_mas
 * @param Halo: vector of structures of type s_Halo, with the coordinates
 * @param weightmass: string denoting the property to use as weight.
 * @details The size of the mesh is characterized by params.Nft The size of the volume is defined by params.Lbox.
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_CIC
 * @author ABA
 */
void getDensity_CIC(s_params_box_mas *params, const vector<s_Halo>& Halo, vector<real_prec>&delta, string weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief CIC + Phase-Space Mapping
 * @param params: N1: number of cells in x-direction
 * @param params: N2: number of cells in y-direction
 * @param params: N3: number of cells in z-direction
 * @param params: L1: Lenght of cuboid in x-direction
 * @param params: L2: Lenght of cuboid in y-direction
 * @param params: L3: Lenght of cuboid in z-direction
 * @param params: d1: Resolution in x-direction
 * @param params: d2: Resolution in y-direction
 * @param params: d3: Resolution in z-direction
 * @param params: min1: Minumum of x-coord
 * @param params: min2: Minumum of y-coord
 * @param params: min3: Minumum of z-coord
 * @param params: xp: vector with x-coordinates
 * @param params: yp: vector with y-coordinates
 * @param params: zp: vector with z-coordinates
 * @param NOBJS : Number of tracers interpolated  on the mnesh
 * @returns delta vector, \f$ \rho(\vec{x})\f$, modified by TET
 * @author KSK. Parallelized and optimzed by ABA.
 */
void getDensity_TETCIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec> &yp, const vector<real_prec>&zp, ULONG N_OBJ, vector<real_prec>&delta);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using TSC
 * @param params: N1: number of cells in x-direction
 * @param params: N2: number of cells in y-direction
 * @param params: N3: number of cells in z-direction
 * @param params: L1: Lenght of cuboid in x-direction
 * @param params: L2: Lenght of cuboid in y-direction
 * @param params: L3: Lenght of cuboid in z-direction
 * @param params: d1: Resolution in x-direction
 * @param params: d2: Resolution in y-direction
 * @param params: d3: Resolution in z-direction
 * @param params: min1: Minumum of x-coord
 * @param params: min2: Minumum of y-coord
 * @param params: min3: Minumum of z-coord
 * @param params: xp: vector with x-coordinates
 * @param params: yp: vector with y-coordinates
 * @param params: zp: vector with z-coordinates
 * @param params: Parmass: vector with weights
 * @param weightmass: bool asking to weight using the entries by Parmass
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_TSC
 * @author ABA & FSK.
 */
void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>& zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh using TSC
 * @param params: structure of type s_params_box_mas
 * @param Halo: vector of structures of type s_Halo, with the coordinates
 * @details The size of the mesh is characterized by params.Nft The size of the volume is defined by params.Lbox.
 * @param weightmass: string denoting the property to use as weight.
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh.
 * @warning This function will be either deprecated or modified to implement grid_assignment_TSC
 * @author ABA
 */
void getDensity_TSC(s_params_box_mas *params, const vector<s_Halo>& Halo, vector<real_prec>&delta, string weightmass);
//////////////////////////////////////////////////////////
/**
 * @brief Interpolation on a mesh
 * @param params: structure of type s_params_box_mas
 * @param Halo: vector of structures of type s_Halo, with the coordinates
 * @param weightmass: string denoting the property to use as weight.
 * @details THe size of the mesh is characterized by params.Nft The size of the volume is defined by params.Lbox. The MAS scheme is defined by the element params->masskernel whcih can be 0,1,2. Under this choice, the function calls the corresponding functions in this same file.
 * @returns delta vector, \f$ \rho(\vec{x})\f$, interpolation of number counts field on a mesh
 * @author ABA
 */
void MAS_NEW(s_params_box_mas *params, vector<s_Halo>&, string weight, vector<real_prec> &delta);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Solution of the Poisson equation
 * @details Using FFTW solves for the potential \f$\phi(\vec{x})\f$ from \f$\nabla^{2}\phi(\vec{x})=\delta(\vec{x})\f$, as \f$ \phi(\vec{x})=iFFTW (\delta(\vec{k})/k^{2})\f$.
 * @param Lbox : lenght-size of the box
 * @param Nft : Number of cellls per dimension
 * @param in: overdensity \f$\delta(\vec{x})\f$
 * @returns gravitational potential \f$\phi(\vec{x})\f$
 * @author ABA
*/
void PoissonSolver(real_prec Lbox, ULONG Nft, vector<real_prec>&in, vector<real_prec>&out);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Function used only in LPT.
 * @warning This application will be deprecated
 */
void cellbound(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, vector<real_prec>&v1, vector<real_prec>&v2, vector<real_prec>&v3);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Function used only in LPT.
 * @warning This application will be deprecated
 */
void cellboundcomp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>vi);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Function used only in LPT.
 * @warning This application will be deprecated
 */
void cellbound_comp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&,vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Nearest-gripd point interpolation function
 * @param Coordinate x
 * @returns weight in x
 * @author ABA
 */
real_prec MAS_NGP(real_prec x);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Cloud-in-cell interpolation function
 * @param Coordinate x
 * @returns weight in x
 * @author ABA
 */
real_prec MAS_CIC(real_prec x);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Triangular-shape-cloud interpolation function
 * @param Coordinate x
 * @returns weight in x
 * @author ABA
 */
real_prec MAS_TSC(real_prec x);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Piecewise Cubic Spline interpolation function
 * @param Coordinate x
 * @returns weight in x
 * @author ABA
 */
real_prec MAS_PCS(real_prec x);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Overdensity field
 * @param in: input vector of density field \f$ \rho(\vec{x}) \f$ of size N
 * @returns out: \f$ \delta (\vec{x}) = (\rho(\vec{x})-\bar{\rho})/\bar{\rho}  \f$
 * @details Internally computes the mean field
 * @author ABA
 */
void get_overdens(const vector<real_prec>&, vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Weighted overdensity field
 * @param in: input vector of density field \f$ \rho(\vec{x}) \f$ of size N
 * @param weight: input vector of weights \f$ w(\vec{x}) \f$ of size N
 * @details Computes \f$ \bar{\rho}_{w}=\sum_{i}\rho(x_{i})/ \sum_{i}w(x_{i})\f$
 * @returns out: \f$ \delta (\vec{x}) = (\rho(\vec{x})-w(\vec{x}))/\bar{\rho}_{w}\f$
 * @author ABA
 */
void get_overdens(const vector<real_prec>&, const vector<real_prec>&, vector<real_prec>&, bool);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Overdensity field
 * @param in: input vector of density field \f$ \rho(\vec{x}) \f$ of size N
 * @param mean: mean value of field \f$ \bar{\rho}=\sum_{i}\rho(x_i)/N\f$
 * @returns \f$ \delta (\vec{x}) = (\rho(\vec{x})-\bar{\rho})/\bar{\rho}  \f$
 * @author ABA
 */
void get_overdens(const vector<real_prec>&in, real_prec mean, vector<real_prec>&delta);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gets number of tracers in input field
 * @param in: input vector of density field \f$ \rho(\vec{x}) \f$ of size N
 * @returns \f$ \sum_{i}\rho(x_{i})\f$
 * @author ABA
 */
real_prec get_nobjects(const vector<real_prec>&);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Gets number of tracers in input field
 * @param in: input vector of density field \f$ \rho(\vec{x}) \f$ of size N
 * @returns \f$ \sum_{i}\rho(x_{i})\f$
 * @author ABA
 */
ULONG get_nobjects(const vector<ULONG>&);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Get bin from position
 * @param x coordinate
 * @param x_min minimum of x-coordinate
 * @param nbins number of bins in the x-coordinate
 * @param delta rsolution in the x-coordinate
 * @param w true of coordinates outside range are to be assigned  to the last bin
 * @author ABA
 */
ULONG get_bin(real_prec x, real_prec xmin,ULONG nbins, real_prec delta, bool w);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Tidal Anisotropy
 * @details Computes the tidal anisotropy
 * @param lambda 1: eigenvalue 1 (scalar)
 * @param lambda 2: eigenvalue 2 (scalar)
 * @param lambda 3: eigenvalue 3 (scalar)
  *@returns \f$ T_{A}= \frac{\sqrt{(\lambda_{3} -\lambda_{1})^{2}+(\lambda_{3}-\lambda_{2})^{2}+(\lambda_{2}-\lambda_{1})^{2}}}{2(1+\lambda_{1}+\lambda_{2}+\lambda_{3})} \f$
 * @author ABA
 */
real_prec tidal_anisotropy(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Invariant of tidal field
 * @details Get Invariant of the tidal field tensor
 * @param Eigenvalues of the tidal field lambda1, lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns \f$ I_{1}= \lambda_{1}+\lambda_{2}+\lambda_{3} \f$
 * @author ABA
*/
real_prec invariant_field_I(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief  Tidal field Invariant I2
 * @param Eigenvalues of the tidal field lambda1, lambda2, lambda3 = \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns  \f$ I_{2}= \lambda_{1}\lambda_{2}+\lambda_{2}\lambda_{3}+\lambda_{1}\lambda_{3} \f$
 * @author ABA
 */
real_prec invariant_field_II(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief  Tidal field invariant I3
 * @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns \f$ I_{3}= \lambda_{1}\lambda_{2}\lambda_{3} \f$
 */
real_prec invariant_field_III(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Tidal field invariant I4
 * @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
 * @returns \f$ I_{4}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
 */
real_prec invariant_field_IV(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

* @brief  Ellipticity
* @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
* @returns \f$ I_{3}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
 * @author ABA
*/
real_prec ellipticity(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

* @brief  Invariant IV
* @param Eigenvalues of the tidal field \f$ \lambda_{1},\lambda_{2},\lambda_{3}\f$
* @details \f$ I_{4}= \lambda_{1}^{2}+\lambda_{2}^{2}+\lambda_{3}^{2} \f$
 * @author ABA
*/
real_prec prolatness(real_prec lambda1, real_prec lambda2, real_prec lambda3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Computes Eigenvalues of the velocity shear
 * @returns
 * @author ABA
 */
void EigenValuesVweb(ULONG Nft, real_prec L1, vector<real_prec>&inx,vector<real_prec>&iny,vector<real_prec>&inz, vector<real_prec>&diver, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Computes Eigenvalues ofthe tidal tensor of density filed.
 * @param Lbox Lenght-size of the box
 * @param Nft  Number of cellls per dimension
 * @param delta Input overdensity
 * @param phi Gravitational potential
 * @returns Eigenvalues of the tidal field lambda1, lambda2, lambda3
 * @author ABA & FSK
 */
void EigenValuesTweb(ULONG Nft, real_prec L1,  vector<real_prec> &delta, vector<real_prec> &ohi, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Computes Eigenvalues ofthe tidal tensor of density filed.
 * @param Lbox Lenght-size of the box
 * @param Nft  Number of cellls per dimension
 * @param delta Input overdensity
 * @returns ou11, out2, out 3 Eigenvalues of the tidal field sorted out1>lout2>out3
 * @author ABA & FSK
 */
void EigenValuesTweb(ULONG Nft, real_prec L1,  vector<real_prec> &delta, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Computation of perturbative terms
 * @param Lbox Lenght-size of the box
 * @param Nft  Number of cellls per dimension
 * @param delta \f$\delta(vec{x})\f$ Input overdensity on the mesh
 * @details If requested via pre_proc directives, this function computes the terms \f$s^{2},s^{3}\f$ and \f$\nabla^{2}\delta \f$
 * @returns Perturbative terms
*/
void EigenValuesTweb_bias(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief 2D histogram
 * @param XMin minimum in the x-axis
 * @param XMax maximum in the x-axis
 * @param YMin minimum in the y-axis
 * @param YMax maximum in the y-axis
 * @param Nbin Number of bins per dimension
 * @param field_X Quantity X
 * @param field_Y Quantity Y
 * @param logs True if the histogram is deseired to be done with the variables \f$ x=\log(1+X)\f$ and \f$ x=\log(1+Y)\f$
 * @warning If logs is true, the values of Ymax, min and Xmax, min must be already computed in terms of the variables (x,y)
 * @returns HIST a 2D histogram \f$P(X,Y)\f$ encoded as a 1D container with column-major ordering. The histogram is normalized to its maximum.
 * @author ABA
 */
void get_2d_histogram(real_prec XMin,real_prec XMax,real_prec YMin,real_prec YMax, ULONG Nbins, vector<real_prec>&field_X,vector<real_prec>&field_Y,vector<real_prec>&HIST, bool logs);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief This function generates an integer in the range (0,Nmax-1)
 * @details It is in this set of functions becuase it uses some functions defiend here.
 * @author ABA
 */
int my_gsl_rng_uniform_(gsl_rng *r, int Nmax);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Rank ordering
 * @param seed inizialization for random generator
 * @param dens Vector defining the bins in overdensity
 * @param NK Number of bins
 * @param maxk Maximum value of overdensity
 * @param mink Minimum value of overdensity
 * @param in Original overdenisty field
 * @param pdfin Pdf of the origial field
 * @param pdfout Target pdf
 * @warning Input must be delta, and the limits must be max and mins in \f$log(1+\delta)\f$. The pdfs must be \f$P(log(1+delta))\f$
 * @returns In the same array "in", the new overdensity field is returned foollowing the target pdf.
 * @author FSK & ABA
 */
void rankorder(int seed, vector<real_prec>dens, ULONG Nk, real_prec maxk, real_prec mink, vector<real_prec>&in, vector<real_prec>&pdfin, vector<real_prec>&pdfout);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function computes the pdf of the input field in. builds the kernel from the ratio
 * @details If type is log, the pdf is measured in log(numinlog+in). The inputs min and max must be either in log(numinlog+"in") or in "in"
 * @param type lin or log, if  \f$log(1+\delta)\f$
 * @param maxk Maximum value of overdensity
 * @param mink Minimum value of overdensity
 * @param in Original overdenisty field
 * @warning Input must be delta, and the limits must be max and mins in \f$log(1+\delta)\f$. The pdfs must be \f$P(log(1+delta))\f$
 * @returns pdf Container with the density distribution
 * @author ABA
 */
void calc_pdf(string type, ULONG N, ULONG Nk, real_prec maxk, real_prec mink, const vector<real_prec>&in, vector<real_prec>&pdf);
///////////////////////////////////////////////////////////////////////////
/**

 * @brief Computes the ID of low resolution mesh of a cell defined in the same volume but on a high resolution mesh
 * @param Nft Number of cells/dimention in high resolution mesh
 * @param Nft_low Number of cells/dimention in low resolution mesh
 * @returns id_info containr of lengh  \f$ N_{ft}^{3} \f$
 * @author ABA
 * @author ABA
 */
void  get_high_res_id_from_low_res_id(int Nft, int Nft_low,vector<ULONG>&id_info);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief Computes the ID of a high resolution mesh of a cell defined in the same volume but on a low resolution mesh
 * @param Nft Number of cells/dimention in high resolution mesh
 * @param Nft_low Number of cells/dimention in low resolution mesh
 * @returns cell_id_info_low containr of lengh  \f$ N_{ft}^{3} \f$ and type s_cell_info_reduced
 * @author ABA
 */
void  get_low_res_id_from_high_res_id(int Nft, int Nft_low,vector<s_cell_info_reduced>&cell_info_low);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief  This function finds the index ID of the neighbour cells for all cells in a mesh
 * @param Nft Number of cells/dimension
 * @param N_cells_bf Number of cells (of the fiducial resolution) backward and forward to look for neighbour cells.
 * @returns ncells  Container of dimension \f$N_{ft}^{3}\f$ and type s_nearest_cells with information on neighbouring cells for each cell.
 * @author ABA
 */
void get_neighbour_cells(ULONG Nft, int N_cells_bf, vector<s_nearest_cells>& ncells);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief This function finds the index of the neighbour cells for a given set of emtpy cells located at "cell_index"
 * This is mainly implemented in the determination of the velocity field in the assignment campaigns of BAM
 * @param Nft Number of cells/dimension
 * @param cells_index ID of cells, container of size \f$N_{ft}^{3}\f$
 * @returns ncells  Container of dimension \f$N_{ft}^{3}\f$ and type s_nearest_cells with information on neighbouring cells for each cell.
 * @author ABA
 */
void get_neighbour_cells_of_cell(ULONG Nft, vector<ULONG>&cell_index, int ,vector<s_nearest_cells>&ncells);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief  This function finds the index of the neighbour cells for all cells in a mesh
 * @param Nft Number of cells/dimension
 * @param cells_index ID of cells, container of size \f$N_{ft}^{3}\f$
 * @param exclude if true, the current cell is ignored
 * @returns ncells  Container of dimension \f$N_{ft}^{3}\f$ and type s_nearest_cells with information on neighbouring cells for each cell.
 * @author ABA
 */
void get_neighbour_cells_cat_analyze(int Nft, int N_cells_bf, vector<s_nearest_cells>&, bool exclude);
////////////////////////////////////////////////////////////////////////////
/**

 * @brief  This function finds the index of the neighbour cells for all cells in a mesh
 * @param Nft Number of cells/dimension
 * @param cells_index ID of cells, container of size \f$N_{ft}^{3}\f$
 * @param exclude if true, the current cell is ignored
 * @details This function aims at being used inside a loop within the calculation of mach number and local overdensityies (see \refitem<Catalog> Class::Catalog.get_local_mach_number_chuncks())
 * in order to avoid high RAM consumption, which is the case of using laptops instead of a powerful machine or cluster.
 * The goal is to do the identification of neighbouring cells by traveling the full volume in chunks or planes (or slices) of constant x
 * @returns ncells  Container of dimension \f$N_{ft}^{3}\f$ and type s_nearest_cells with information on neighbouring cells for each cell.
 * @author ABA
 */
void get_neighbour_cells_cat_analyze_chuncks(int Nc, int Nft, int N_cells_bf, vector<s_nearest_cells>&, bool exclude);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Top-hat window function in Fourier space
 * @param k Wavenumber
 * @param R Scale
 * @author ABA
 */
real_prec window(real_prec k, real_prec R);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Gaussian window function in Fourier Space
 * @param k Wavenumber
 * @param R Scale
 * @author ABA
 */
real_prec windowg(real_prec k, real_prec R);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief This function assingns a velocity to an empty cell by averaging the velocities in the surrounding cells.
 * @param id_empty_cells  Container with the ID of empty cells
 * @param nearest_cells  Contaienr of structures with the information of the neighbouring cells to the current cells
 * @returns velf Container with the value of the velocity field at the emptu cell
 * @author ABA
 */
void get_vel_field_using_neigh_cells(vector<ULONG>&id_empty_cells,vector<s_nearest_cells>&nearest_cells,vector<real_prec>&velf);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Poisson Shot-noise Correcting mode-by-mode
 *
 * @param e Index denoting the MAS (0 for NGP, 1 for CIC, 2 for TSC, 3 for PCS)
 * @param xx kx-coordinate in units of the size of the cell in the kx-direction
 * @param yy ky-coordinate in units of the size of the cell in the ky-direction
 * @param zz kz-coordinate in units of the size of the cell in the kz-direction
 * @returns The product \f$S_{e}(x)S_{e}(y)S_{e}(z)\f$ where
  \f$ S(x)=\begin{cases}1 & e = 0 \\ 1-\frac{2}{3}x^{2} & e = 1\\ 1-x^{2}+\frac{2}{15}x^{4} & e = 2 \\
  1-\frac{4}{3}x^{2}+\frac{2}{5}x^{4}-\frac{4}{315}x^{6} & e = 3  \\
  \end{cases}\f$ and \f$ x_{i}=\sin(i_x\pi/N_{ft})/(i_x\pi/N_{ft}) \f$ used to correct the 3D power spectrum
 * @author ABA
 */
real_prec SN_correction_MAS(int e,real_prec xx,real_prec yy,real_prec zz);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Mass assignment correction mode-by-mode
 * @param e Index denoting the MAS (0 for NGP, 1 for CIC, 2 for TSC, 3 for PCS)
 * @param xx kx-coordinate in units of the size of the cell in the kx-direction
 * @param yy ky-coordinate in units of the size of the cell in the ky-direction
 * @param zz kz-coordinate in units of the size of the cell in the kz-direction
 * @returns The product \f$(S(xx)S(yy)S(zz))^{e+1}\f$ used to correct the 3D power spectrum
 * @author ABA
 */
real_prec correction_MAS(int exp,real_prec xx,real_prec yy,real_prec z);

////////////////////////////////////////////////////////////////////////////
/**
* @brief Assignmet of individual bias, 
* @details This funciton is meant to be used in Catlog, since there is a conflict defining the class PowerSpectrumF in the Catalog class.
*/
void object_by_object_bias(Params, vector<s_Halo>&, vector<real_prec>&);

////////////////////////////////////////////////////////////////////////////
/**

* @brief Mass assignment correction mode-by-mode
* @param e Index denoting the MAS (0 for NGP, 1 for CIC, 2 for TSC, 3 for PCS)
* @param i k-index  in the kx-direction
* @param j k-index  in the ky-direction
* @param i k-index  in the kz-direction
* @returns The product \f$(S(i\pi/N_{ft})S(jx\pi/N_{ft})S(k\pi/N_{ft}))^{e+1}\f$ used to correct the 3D power spectrum
* @author ABA
*/
real_prec correction_MAS(ULONG Nft, int expo, int i,int j,int k);
////////////////////////////////////////////////////////////////////////////
/**
 * @brief Generate an interpolated value of the mean number density
 * @details Given the position in the sky of the object, when the depth varies with the angular position. This interpolates the matrix computed on the method dndz
 * @param zmin Minimum redshift of the sample
 * @param zmax Maximum redshift of the sample
 * @param ra Right ascencion of the galaxy
 * @param ra Declination of the galaxy
 * @param zg Redshift of the galaxy
 * @param dndz_m container with the vaues of the mean number density in redshift bins and Healpix pixels tabulated.
 * @result nbar Mean number denisiy tabulated at the position of the galaxy
 */
void get_mean_density_interpolated(real_prec zmin, real_prec zmax, real_prec ra, real_prec dec, real_prec zg, vector< vector<real_prec> >&dndz_m, real_prec *nbar);
////////////////////////////////////////////////////////////////////////////
/**
* @brief Remap indices
* @details  This auxiliary function returns the indices that are
* arguments of the function ijk(), used to retrieve (a C++ note, ijk is  a method of same class, I do not need to create another FFTW_FUNCTION object)
* the amplitudes of the DTF as given by the FFTW. This takes into account  the fact that we are adding one more frequency
* the -Nyquist to the output and that we use Hermitian
   symmetry to retrieve the negative section
   of the third component, i.e,
   \f$ \delta(k_x, k_y, -k_z)=\delta{*}(-_kx, -k_y, k_z)\f$. sThe ofactor is +1 for the real part,
   -1 for the imaginary part when Hermitian symmetry is explicitely used.
¨* @author ABA
*/
void remap(int nn1, int nn2, int nn3, int i, int j, int k, int *ooi, int *ooj, int *ook, real_prec *ofactor);
////////////////////////////////////////////////////////////////////////////
#ifndef _MOCabc_
#define _MOCabc_
/**
* @brief Get mean Occupation number
 * @author ABA
*/
template<typename Type> Type get_mean_occupation_number(real_prec max, const vector<Type>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  ULONG Ntot=0;
  int Nd=static_cast<int>(max);
  vector<int>dist(Nd+1,0);
  real_prec delta=1;
  for(int i=0;i<in.size();++i)
      dist[get_bin(in[i],0,Nd,delta,true)]++;
 int mean=0;
 Ntot=0;
 for(int i=0;i<=Nd;++i)
    {
     mean+=i*dist[i];
     Ntot+=dist[i];
  }
 mean/=static_cast<real_prec>(Ntot);
 return static_cast<int>(mean);
}
#endif
////////////////////////////////////////////////////////////////////////////
#ifdef HEALPIX
inline void my_get_mean_density_interpolated(Healpix_Map<real_prec> &map, int new_n_dndz, real_prec redshift_min_sample, real_prec redshift_max_sample, real_prec ra, real_prec dec, real_prec redshift, vector< vector<real_prec> > &dndz_m, real_prec *nb)
{
// Function used when the mean number density is not tabulated
  pointing point;
  real_prec fac=M_PI/180.0;
  point.phi=ra*fac;
  point.theta=0.5*M_PI-dec*fac;
  long ipix=map.ang2pix(point);
  real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/dndz_m.size();
  int iz=(int)floor((float)((redshift-redshift_min_sample)/Delta_Z));
  *nb=dndz_m[iz][ipix];
}
#endif
////////////////////////////////////////////////////////////////////////////
#endif

////////////////////////////////////////////////////////////////////////////
void print_catalog(vector<s_Halo>&Halo, string file, bool);
////////////////////////////////////////////////////////////////////////////

void comp_time_MH(time_t start, long full, long step, real_prec H);
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
  /**
   * @brief  Add the elements of input container
   */
void comp_time(time_t start, long full, long step);
