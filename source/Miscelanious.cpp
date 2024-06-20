////////////////////////////////////////////////////////////////////////////
/** @file Miscelanious.cpp
 * @brief This file contains functions used in several applications
 * @author Andres Balaguera Antolinez 2017-2023
 */
////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <math.h>
#include <fstream>
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf.h>
using namespace std;
////////////////////////////////////////////////////////////////////////////
#include "Miscelanious.h"
////////////////////////////////////////////////////////////////////////////
real_prec window(real_prec k,real_prec R){
  /*Top hat window function, normalized suich that \int \dtx W(\xv) = 1*/
  real_prec y=k*R;
  return 3.*(sin(y)/pow(y,3)-cos(y)/pow(y,2));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec windowg(real_prec k,real_prec R){
  return exp(-0.5*pow(k*R,2));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void remap(int nn1, int nn2, int nn3, int i, int j, int k, int *ooi, int *ooj, int *ook, real_prec *ofactor)
{
  if(k<=nn3/2){
    if(i<=nn1/2)*ooi=i; else *ooi=i-1;
    if(j<=nn2/2)*ooj=j; else *ooj=j-1;
    *ook=k;
    *ofactor=1.00;
  }
  else{
    if(k>=nn3/2+1){  //negative components of z
      if(i>0 && i<=nn1/2)*ooi=nn1-i;else if(i>=nn1/2+1)*ooi=nn1-i+1;
      if(j>0 && j<=nn2/2)*ooj=nn2-j;else if(j>=nn2/2+1)*ooj=nn2-j+1;
      if(i==0)*ooi=0;
      if(j==0)*ooj=0;
      *ook=nn3-(k-1);
      *ofactor=-1.00;
    }
  }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//#define parallel_ngp
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass)
{
#ifdef parallel_ngp
#define DELTAnb(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAnb(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  ScreenOutput So;
  ULONG N_OBJ=xp.size();
#ifdef _FULL_VERBOSE_
  So.message_screen("Interpolating on a grid using NGP for", N_OBJ, " objects");
#endif
#ifndef _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
#endif
  ULONG NLOSS=0;
  ULONG NTOT=0;
#ifdef parallel_ngp
#pragma omp parallel
  {
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp for reduction(+:NTOT,NLOSS)
#endif
    for (ULONG n=0; n<N_OBJ; n++)
      {
	//check if particle is in selected Domain, else discard it
	real_prec xnew=xp[n];
	real_prec ynew=yp[n];
	real_prec znew=zp[n];
#ifdef _PERIODIC_BC_MAS_
	if (xnew<min1)
	  xnew+=L1;
	if (xnew>=L1)
	  xnew-=L1;
	if (ynew<min2)
	  ynew+=L1;
	if (ynew>=L1)
	  ynew-=L1;
	if (znew<min3)
	  znew+=L1;
	if (znew>=L1)
	  znew-=L1;
#endif
	if((xnew>=min1 && xnew<min1+L1) && (ynew>=min2 && ynew<min2+L2) && (znew>=min3 && znew<min3+L3))
	  {
	    ULONG i = get_bin(xnew,min1,N1,d1,false);
	    ULONG j = get_bin(ynew,min2,N2,d2,false);
	    ULONG k = get_bin(znew,min3,N3,d3,false);
	    real_prec mass=num_1;
	    if (true==weightmass)
	      mass=Par_mass[n];
	    DELTAnb(i,j,k) +=mass;
	    NTOT++;
	  }
	else
	  NLOSS++;
      }
#ifdef parallel_ngp
#pragma omp critical
    {
      for(ULONG i=0;i<delta.size();++i)
	delta[i]+=private_delta[i];
      private_delta.clear();private_delta.shrink_to_fit();
    }
  }
#endif
#ifndef  _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
  if(false==weightmass)
    {
#ifdef _FULL_VERBOSE_
      So.message_screen("Number of objects assigned to grid =",NTOT);
#endif
#endif
#ifdef _FULL_VERBOSE_
      if(NLOSS!=0)
	if(NLOSS!=0) So.message_screen("mass assignment found",NLOSS," particles outside mesh boundary");
#endif
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//#define parallel_ngp_b
void getDensity_NGP(s_params_box_mas *params, vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
#ifdef parallel_ngp_b
#define DELTAno(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAno(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  ScreenOutput So;
  ULONG N_OBJ=Halo.size(); 
#ifdef _FULL_VERBOSE_
  So.message_screen("Number of objects to be assigned =", N_OBJ);
#endif
  real_prec min1=params->min1;
  real_prec min2=params->min2;
  real_prec min3=params->min3;
  real_prec L1=params->Lbox;
  real_prec L2=params->Lbox;
  real_prec L3=params->Lbox;
  ULONG N1=params->Nft;
  ULONG N2=params->Nft;
  ULONG N3=params->Nft;
  real_prec d1=params->d1;
  real_prec d2=params->d2;
  real_prec d3=params->d3;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  vector<real_prec>pwe(N_OBJ,1);
  if(weight_prop!=_COUNTS_)
    pwe.resize(N_OBJ,1);
  if (weight_prop==_MASS_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].mass;
    }
  else if (weight_prop==_SAT_FRACTION_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].number_sub_structures;
    }
  else if (weight_prop==_RS_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].rs;
    }
  else if (weight_prop==_VMAX_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].vmax;
    }
  else if (weight_prop==_VIRIAL_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].virial;
    }
  else if (weight_prop==_SPIN_)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for (ULONG n=0; n<N_OBJ; n++)
	pwe[n]=Halo[n].spin;
    }
  ULONG NLOSS=0;
  ULONG NTOT=0;
#ifdef  parallel_ngp_b
#pragma omp parallel
  {
    vector<real_prec> private_delta(delta.size(),0);
#pragma omp for reduction(+:NTOT,NLOSS)
#endif
    for (ULONG n=0; n<N_OBJ; n++)
      { 
        real_prec xp=Halo[n].coord1-min1; // Since I add -min1, I implicitly set minimums to 0 below
        real_prec yp=Halo[n].coord2-min2;
        real_prec zp=Halo[n].coord3-min3;
	if((xp >=0 && xp <=L1) && (yp >=0 && yp <=L2) && (zp >=0 && zp <=L3))
	  {
            ULONG i = static_cast<ULONG>(floor(xp/d1 + BIN_SHIFT)); // indices of the cell of the particle
            ULONG j = static_cast<ULONG>(floor(yp/d2 + BIN_SHIFT));
            ULONG k = static_cast<ULONG>(floor(zp/d3 + BIN_SHIFT));
            i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
            j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
            k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
            if(weight_prop==_COUNTS_)
	      DELTAno(i,j,k)++;
            else
	      DELTAno(i,j,k)+= pwe[n];
            NTOT++;
	  }
        else
	  {
	    NLOSS++;
	  }
      }
#ifdef  parallel_ngp_b
#pragma omp critical
    {
      for(ULONG i=0;i<delta.size();++i)
        delta[i]+=private_delta[i];
    }
    private_delta.clear();private_delta.shrink_to_fit();
  }
#endif

  So.DONE();
  pwe.clear();pwe.shrink_to_fit();
  
  if(_COUNTS_ == weight_prop)
    {
      ULONG count=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count)
#endif
      for(ULONG i=0;i<delta.size();++i)
	count+=delta[i];
#ifdef _FULL_VERBOSE_
      So.message_screen("Number of objects assigned to grid =",count);
      So.message_screen("Number of objects assigned to grid =",NTOT);
      if(NLOSS>0)
	So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");
#endif
      for(ULONG i=0;i<delta.size();++i)
        if(delta[i]>50)cout<<i<<"  " <<delta[i]<<endl;
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_tetrahedron_centroid( real_prec d1, real_prec d2, real_prec d3, const vector<real_prec>&posx, const vector<real_prec>&posy, const vector<real_prec>&posz, const vector<int>&vertids, real_prec boxsize, real_prec boxhalf, real_prec &xp,real_prec &yp,real_prec &zp )
{
  real_prec x0_1=posx[vertids[0]]-num_0_5*d1;
  real_prec x0_2=posy[vertids[0]]-num_0_5*d2;
  real_prec x0_3=posz[vertids[0]]-num_0_5*d3;
  xp= 0.0;
  yp= 0.0;
  zp= 0.0;
  for(int i=1; i<4; ++i )
    {
      vector<real_prec> dx(3,0);
      dx[0]=posx[vertids[i]]-x0_1-num_0_5*d1;
      dx[1]=posy[vertids[i]]-x0_2-num_0_5*d2;
      dx[2]=posz[vertids[i]]-x0_3-num_0_5*d3;   
      
      if( dx[0] < -boxhalf )
	dx[0] += boxsize;
      else if( dx[0] > boxhalf )
	dx[0] -= boxsize;
      if( dx[1] < -boxhalf )
	dx[1] += boxsize;
      else if( dx[1] > boxhalf )
	dx[1] -= boxsize;
      if( dx[2] < -boxhalf )
	dx[2] += boxsize;
      else if( dx[2] > boxhalf )
	dx[2] -= boxsize;
      xp += dx[0];
      yp += dx[1];
      zp += dx[2];
    }
  xp = fmodf(0.25*xp+x0_1+boxsize, boxsize );
  yp = fmodf(0.25*yp+x0_2+boxsize, boxsize );
  zp = fmodf(0.25*zp+x0_3+boxsize, boxsize );
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define parallel_cictet
#ifdef parallel_cictet
//#define parallel_cictet_shared
#endif
void getDensity_TETCIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp,ULONG N_OBJ, vector<real_prec>&delta)
{
  
#if defined parallel_cictet && !defined (parallel_cictet_shared)
#define DELTAco(i,j,k) private_delta[k+N3*(j+N2*i)]
#else 
#define DELTAco(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  int NTHREADS=1;
#ifdef _USE_OMP_
  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  //----> connectivity for cubic grid, six tetrahedron decomposition
  real_prec boxhalf = num_0_5 * L1;
  const int vert[8][3] = { {0,0,0}, {0,0,1}, {0,1,0}, {0,1,1}, {1,0,0}, {1,0,1}, {1,1,0}, {1,1,1} };
  const int nbase = static_cast<int>(pow(static_cast<real_prec>(N_OBJ),1.0/3.0)+0.1);
  const int conn[6][4] = { {1,0,2,4}, {3,1,2,4}, {3,5,1,4}, {3,6,5,4}, {3,2,6,4}, {3,7,5,6} };
  ULONG Nbase_cube=nbase*nbase*nbase;
  real_prec pweight = static_cast<real_prec>(1.0/6.0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); ++i)
    delta[i]= 0.; 
  int ix,iy,iz,li;
#ifdef parallel_cictet
#ifndef parallel_cictet_shared
#pragma omp parallel private(ix,iy,iz,li)
  { 
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp for nowait
#else
#pragma omp parallel for collapse(4)
#endif
#endif
    for(ix=0; ix<nbase; ++ix)
      for(iy=0; iy<nbase; ++iy)
	for(iz=0; iz<nbase; ++iz)
	  for(li=0; li<6; ++li)
	    {
	      vector<int> vertids(4,0);
	      for(int m=0; m<4; ++m)
		vertids[m] = (((ix+vert[conn[li][m]][0])%nbase)*nbase + (iy+vert[conn[li][m]][1])%nbase)*nbase + (iz+vert[conn[li][m]][2])%nbase;
	      real_prec xpos, ypos, zpos;
	      get_tetrahedron_centroid(d1, d2, d3, xp, yp, zp, vertids, L1, boxhalf, xpos, ypos,zpos);
	      if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
		{
		  ULONG i = get_bin(xpos,min1,N1,d1,false);
		  ULONG j = get_bin(ypos,min2,N2,d2,false);
		  ULONG k = get_bin(zpos,min3,N3,d3,false);
		  ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
		  ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
		  ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
		  real_prec xc = static_cast<real_prec>(i);  
		  real_prec yc = static_cast<real_prec>(j);
		  real_prec zc = static_cast<real_prec>(k);
		  real_prec dx = (xpos-min1)/d1 - xc; 
		  real_prec dy = (ypos-min2)/d2 - yc;
		  real_prec dz = (zpos-min3)/d3 - zc;
		  real_prec tx = num_1 - dx;
		  real_prec ty = num_1 - dy;
		  real_prec tz = num_1 - dz;
#ifdef parallel_cictet_shared
#pragma omp task firstprivate(i,j,k) shared(delta)
                  {
#pragma omp atomic
#endif
		    DELTAco(i,j,k)   +=  pweight*tx*ty*tz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif

		    DELTAco(ii,j,k)  +=  pweight*dx*ty*tz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(i,jj,k)  +=  pweight*tx*dy*tz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(i,j,kk)  +=  pweight*tx*ty*dz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(ii,jj,k) +=  pweight*dx*dy*tz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(ii,j,kk) +=  pweight*dx*ty*dz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(i,jj,kk) +=  pweight*tx*dy*dz;
#ifdef parallel_cictet_shared
#pragma omp atomic
#endif
		    DELTAco(ii,jj,kk)+=  pweight*dx*dy*dz;
#ifdef parallel_cictet_shared
		  }
#endif
		}
	    }
#ifdef parallel_cictet
#ifndef parallel_cictet_shared
#pragma omp critical
    {
      for(ULONG i = 0; i < delta.size(); i++)
	delta[i]+=private_delta[i];
    }
  }
#endif 
#endif 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define parallel_cic
void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, const vector<real_prec> &Par_mass, vector<real_prec>&delta, bool weightmass)
{
#ifdef parallel_cic
#define DELTAc(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAc(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  ScreenOutput So;
  ULONG N_OBJ=xp.size();
#ifdef _FULL_VERBOSE_
  So.message_screen("Interpolating on a grid using CIC for", N_OBJ, " objects");
#endif
  int NTHREADS=1;
#ifdef _USE_OMP_
  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif  
#ifndef  _GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif  
  ULONG NTOT=0;
  ULONG NLOSS=0;
#ifdef parallel_cic
#pragma omp parallel 
  { // opens parallel section
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp for reduction(+:NTOT, NLOSS)
#endif
    for (ULONG n=0; n<N_OBJ; n++)
      {
	real_prec xpos=xp[n]-(num_0_5)*d1;
      	real_prec ypos=yp[n]-(num_0_5)*d2;
      	real_prec zpos=zp[n]-(num_0_5)*d3; 
	real_prec xnew=xpos;
	real_prec ynew=ypos;
	real_prec znew=zpos;
#ifdef _PERIODIC_BC_MAS_
        if (xnew<min1)
          xnew+=L1;
        if (xnew>=L1+min1)
          xnew-=L1;
        if (ynew<min2)
          ynew+=L1;
        if (ynew>=L1+min2)
          ynew-=L1;
        if (znew<min3)
          znew+=L1;
        if (znew>=L1+min3)
          znew-=L1;
        xpos=xnew;
	ypos=ynew;
	zpos=znew;
#endif
        if((xpos>=min1 && xpos<min1+L1) && (ypos>=min2 && ypos<min2+L2) && (zpos>=min3 && zpos<min3+L3))
	  {	
	    ULONG i = get_bin(xpos,min1,N1,d1,false);
	    ULONG j = get_bin(ypos,min2,N2,d2,false);
	    ULONG k = get_bin(zpos,min3,N3,d3,false);
	    ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	    ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	    ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	    real_prec xc = static_cast<real_prec>(i);
	    real_prec yc = static_cast<real_prec>(j);
	    real_prec zc = static_cast<real_prec>(k);
	    real_prec dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	    real_prec dy = (ypos-min2)/d2 - yc;
	    real_prec dz = (zpos-min3)/d3 - zc;
	    real_prec tx = num_1 - dx;
	    real_prec ty = num_1 - dy;
	    real_prec tz = num_1 - dz;
	    real_prec mass=num_1;
	    if (true==weightmass)
	      mass=Par_mass[n]; 
	    DELTAc(i,j,k)   += mass*tx*ty*tz;
	    DELTAc(ii,j,k)  += mass*dx*ty*tz;
	    DELTAc(i,jj,k)  += mass*tx*dy*tz;
	    DELTAc(i,j,kk)  += mass*tx*ty*dz;
	    DELTAc(ii,jj,k) += mass*dx*dy*tz;
	    DELTAc(ii,j,kk) += mass*dx*ty*dz;
	    DELTAc(i,jj,kk) += mass*tx*dy*dz;
	    DELTAc(ii,jj,kk)+= mass*dx*dy*dz;
	    NTOT++;
	  }
	else{
	  NLOSS++;
        }
      }
#ifdef parallel_cic
#pragma omp critical
    {
      for(ULONG i = 0; i < delta.size(); i++)
	delta[i] += private_delta[i];
    }
    private_delta.clear(); private_delta.shrink_to_fit();
  } 
#endif 
  
  So.message_screen("Number of objects assigned to grid =",NTOT);
  So.message_screen("Number of objects NOT  assigned to grid =",NLOSS);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void getDensity_CIC(s_params_box_mas *params,const vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
#ifdef parallel_cic
#define DELTAcc(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAcc(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  int NTHREADS=1;
#ifdef _USE_OMP_
  NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N_OBJ=Halo.size();
  real_prec min1=params->min1;
  real_prec min2=params->min2;
  real_prec min3=params->min3;
  real_prec L1=params->Lbox;
  real_prec L2=params->Lbox;
  real_prec L3=params->Lbox;
  real_prec N1=params->Nft;
  real_prec N2=params->Nft;
  real_prec N3=params->Nft;
  real_prec d1=params->d1;
  real_prec d2=params->d2;
  real_prec d3=params->d3;
#ifdef _USE_OMP_  
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  ULONG NLOSS=0;
  ULONG NTOT=0;
#ifdef parallel_cic
#pragma omp parallel
  {
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp parallel reduction(+:NLOSS,NTOT)    
#endif
    for (ULONG n=0; n<N_OBJ; n++)
      {
        real_prec xpos=Halo[n].coord1-num_0_5*d1-min1;
        real_prec ypos=Halo[n].coord2-num_0_5*d2-min2;
        real_prec zpos=Halo[n].coord3-num_0_5*d3-min3;
        real_prec xnew=xpos;
        real_prec ynew=ypos;
        real_prec znew=zpos;
#ifdef _PERIODIC_BC_MAS_
        if (xnew<0)
	  xnew+=L1;
        if (xnew>=L1)
	  xnew-=L1;
	if (ynew<0)
	  ynew+=L1;
        if (ynew>=L1)
	  ynew-=L1;
	if (znew<0)
	  znew+=L1;
        if (znew>=L1)
	  znew-=L1;
        xpos=xnew;
	ypos=ynew;
	zpos=znew;
#endif      
	//check if particle is in selected Domain, else discard it
        if((xpos>=0 && xpos<L1) && (ypos>=0 && ypos<L2) && (zpos>0 && zpos<L3))
	  {	
	    NTOT++;
      ULONG i = static_cast<ULONG>(floor((xpos)/d1+BIN_SHIFT)); // indices of the cell of the particle
      ULONG j = static_cast<ULONG>(floor((ypos)/d2+BIN_SHIFT));
      ULONG k = static_cast<ULONG>(floor((zpos)/d3+BIN_SHIFT));
	    i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	    j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	    k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	    ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	    ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	    ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	    real_prec xc = static_cast<real_prec>(i); 
	    real_prec yc = static_cast<real_prec>(j);
	    real_prec zc = static_cast<real_prec>(k);
      real_prec dx = (xpos)/d1 - xc; // distance of particle to center of the cell
      real_prec dy = (ypos)/d2 - yc;
      real_prec dz = (zpos)/d3 - zc;
	    real_prec tx = num_1 - dx;
	    real_prec ty = num_1 - dy;
	    real_prec tz = num_1 - dz;
	    real_prec tracer_weight=num_1;
	    if (weight_prop==_MASS_)
	      tracer_weight=Halo[n].mass;
	    else if (weight_prop==_SAT_FRACTION_)
	      tracer_weight=Halo[n].number_sub_structures;
	    else if (weight_prop==_RS_)
	      tracer_weight=Halo[n].rs;
	    else if (weight_prop==_VMAX_)
	      tracer_weight=Halo[n].vmax;
	    else if (weight_prop==_VIRIAL_)
	      tracer_weight=Halo[n].virial;
	    else if (weight_prop==_SPIN_)
	      tracer_weight=Halo[n].spin;
	    else if (weight_prop=="_BIAS_")
	      tracer_weight=Halo[n].bias;
	    DELTAcc(i,j,k)   += tracer_weight*tx*ty*tz;
	    DELTAcc(ii,j,k)  += tracer_weight*dx*ty*tz;
	    DELTAcc(i,jj,k)  += tracer_weight*tx*dy*tz;
	    DELTAcc(i,j,kk)  += tracer_weight*tx*ty*dz;
	    DELTAcc(ii,jj,k) += tracer_weight*dx*dy*tz;
	    DELTAcc(ii,j,kk) += tracer_weight*dx*ty*dz;
	    DELTAcc(i,jj,kk) += tracer_weight*tx*dy*dz;
	    DELTAcc(ii,jj,kk)+= tracer_weight*dx*dy*dz;
	  }
	else NLOSS++;	
      }
#ifdef parallel_cic
#pragma omp critical
    {
      for(ULONG i = 0; i < delta.size(); i++)
	delta[i] += private_delta[i];
    }
    private_delta.clear(); private_delta.shrink_to_fit();
  } 
#endif 
  ScreenOutput So;
  So.DONE();
  if(_COUNTS_ == weight_prop)
    {
      ULONG count=get_nobjects(delta);
      So.message_screen("Number of objects assigned to grid =",count);
      if(NLOSS!=0)
	So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");    
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#define parallel_tsc
void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3, const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, const vector<real_prec> & Par_mass, vector<real_prec>&delta, bool weightmass)
{
#ifdef parallel_tsc
#define DELTAt(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAt(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  ULONG N_OBJ=xp.size();
#if !defined(_GET_INTERPOLATED_FIELDS_FROM_BIN_FILES_)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i] = 0.;  //-1 if we want to calculate overdensity
#endif
  ULONG NLOSS=0;
  ULONG NTOT=0;
#ifdef parallel_tsc
#pragma omp parallel
  {
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp for reduction(+:NLOSS,NTOT)
#endif      
    for (ULONG n=0; n<N_OBJ; n++)
      {
	real_prec xpos=xp[n];
	real_prec ypos=yp[n];
	real_prec zpos=zp[n];
	real_prec xnew=xpos;
	real_prec ynew=ypos;
	real_prec znew=zpos;
#ifdef _PERIODIC_BC_MAS_      
	if (xnew<min1)
	  xnew+=L1;
	if (xnew>=L1+min1)
	  xnew-=L1;
	if (ynew<min2)
	  ynew+=L1;
	if (ynew>=L1+min2)
	  ynew-=L1;
	if (znew<min3)
	  znew+=L1;
	if (znew>=L1+min3)
	  znew-=L1;
	xpos=xnew;
	ypos=ynew;
	zpos=znew;
#endif      
	if((xpos>=min1 && xpos<=min1+L1) && (ypos>=min2 && ypos<=min2+L2) && (zpos>=min3 && zpos<=min3+L3))
	  {
	    ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
	    ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2));
	    ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3));
	    i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	    j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	    k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	    ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	    ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	    ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	    ULONG iii= static_cast<ULONG>(fmod(real_prec(i-1+N1),real_prec(N1)));
	    ULONG jjj= static_cast<ULONG>(fmod(real_prec(j-1+N2),real_prec(N2)));
	    ULONG kkk= static_cast<ULONG>(fmod(real_prec(k-1+N3),real_prec(N3)));
	    real_prec xc = static_cast<real_prec>(i)+static_cast<real_prec>(0.5); // centers of the cells
	    real_prec yc = static_cast<real_prec>(j)+static_cast<real_prec>(0.5);
	    real_prec zc = static_cast<real_prec>(k)+static_cast<real_prec>(0.5);
	    real_prec dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	    real_prec dy = (ypos-min2)/d2 - yc;
	    real_prec dz = (zpos-min3)/d3 - zc;
	    real_prec hx0=static_cast<real_prec>(0.75-dx*dx); // fraction of particle assigned
	    real_prec hy0=static_cast<real_prec>(0.75-dy*dy); // fraction of particle assigned
	    real_prec hz0=static_cast<real_prec>(0.75-dz*dz); // fraction of particle assigned
	    real_prec hxp1=static_cast<real_prec>(0.5*(0.5+dx)*(0.5+dx)); // fraction of particle assigned
	    real_prec hyp1=static_cast<real_prec>(0.5*(0.5+dy)*(0.5+dy)); // fraction of particle assigned
	    real_prec hzp1=static_cast<real_prec>(0.5*(0.5+dz)*(0.5+dz)); // fraction of particle assigned
	    real_prec hxm1=static_cast<real_prec>(0.5*(0.5-dx)*(0.5-dx)); // fraction of particle assigned
	    real_prec hym1=static_cast<real_prec>(0.5*(0.5-dy)*(0.5-dy)); // fraction of particle assigned
	    real_prec hzm1=static_cast<real_prec>(0.5*(0.5-dz)*(0.5-dz)); // fraction of particle assigned
	    real_prec mass=num_1;
	    if (weightmass==true)
	      mass=Par_mass[n];
	    NTOT++;
	    DELTAt(i,j,k)    	+= mass*hx0*hy0*hz0;
	    DELTAt(ii,jj,kk) 	+= mass*hxp1*hyp1*hzp1;
	    DELTAt(iii,jjj,kkk)+= mass*hxm1*hym1*hzm1;
	    DELTAt(ii,j,k)   	+= mass*hxp1*hy0*hz0;
	    DELTAt(i,jj,k)   	+= mass*hx0*hyp1*hz0;
	    DELTAt(i,j,kk)   	+= mass*hx0*hy0*hzp1;
	    DELTAt(ii,jj,k)  	+= mass*hxp1*hyp1*hz0;
	    DELTAt(ii,j,kk)  	+= mass*hxp1*hy0*hzp1;
	    DELTAt(i,jj,kk)  	+= mass*hx0*hyp1*hzp1;
	    DELTAt(iii,jj,kk)  += mass*hxm1*hyp1*hzp1;
	    DELTAt(ii,jjj,kk)  += mass*hxp1*hym1*hzp1;
	    DELTAt(ii,jj,kkk) += mass*hxp1*hyp1*hzm1;
	    DELTAt(iii,jjj,kk)+= mass*hxm1*hym1*hzp1;
	    DELTAt(iii,jj,kkk)+= mass*hxm1*hyp1*hzm1;
	    DELTAt(ii,jjj,kkk)+= mass*hxp1*hym1*hzm1;
	    DELTAt(i,jjj,kkk) += mass*hx0*hym1*hzm1;
	    DELTAt(iii,j,kkk) += mass*hxm1*hy0*hzm1;
	    DELTAt(iii,jjj,k) += mass*hxm1*hym1*hz0;
	    DELTAt(i,j,kkk)  	+= mass*hx0*hy0*hzm1;
	    DELTAt(i,jjj,k)  	+= mass*hx0*hym1*hz0;
	    DELTAt(iii,j,k)  	+= mass*hxm1*hy0*hz0;
	    DELTAt(i,jj,kkk)  += mass*hx0*hyp1*hzm1;
	    DELTAt(ii,j,kkk) += mass*hxp1*hy0*hzm1;
	    DELTAt(iii,jj,k)  += mass*hxm1*hyp1*hz0;
	    DELTAt(ii,jjj,k)  += mass*hxp1*hym1*hz0;
	    DELTAt(i,jjj,kk)  += mass*hx0*hym1*hzp1;
	    DELTAt(iii,j,kk)  += mass*hxm1*hy0*hzp1;
	  }	
      }
#ifdef parallel_tsc
#pragma omp critical
    {
      for(ULONG i=0;i<delta.size();++i)
	delta[i]+=private_delta[i];	
    }      
    private_delta.clear();private_delta.shrink_to_fit();
  }
#endif
  if(NLOSS!=0) cout << " >>> mass assignment found "<<NLOSS<<" particles outside mesh boundary...."<<endl<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void getDensity_TSC(s_params_box_mas *params, const vector<s_Halo>&Halo, vector<real_prec>&delta, string weight_prop)
{
#ifdef parallel_tsc
#define DELTAtt(i,j,k) private_delta[k+N3*(j+N2*i)]
#else
#define DELTAtt(i,j,k) delta[k+N3*(j+N2*i)]
#endif
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  ULONG N_OBJ=Halo.size();
  real_prec min1=params->min1;
  real_prec min2=params->min2;
  real_prec min3=params->min3;
  real_prec L1=params->Lbox;
  real_prec L2=params->Lbox;
  real_prec L3=params->Lbox;
  real_prec N1=params->Nft;
  real_prec N2=params->Nft;
  real_prec N3=params->Nft;
  real_prec d1=params->d1;
  real_prec d2=params->d2;
  real_prec d3=params->d3;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
  ULONG NLOSS=0;
  ULONG NTOT=0;
#ifdef parallel_tsc
#pragma omp parallel
  {
    vector<real_prec>private_delta(delta.size(),0);
#pragma omp for reduction(+:NTOT)
#endif
    for (ULONG n=0; n<N_OBJ; n++)
      {
	real_prec xpos=Halo[n].coord1;
	real_prec ypos=Halo[n].coord2;
	real_prec zpos=Halo[n].coord3;
	real_prec xnew=xpos;
	real_prec ynew=ypos;
	real_prec znew=zpos;
#ifdef _PERIODIC_BC_MAS_      
	if (xnew<min1)
	  xnew+=L1;
	if (xnew>=L1+min1)
	  xnew-=L1;
	if (ynew<min1)
	  ynew+=L1;
	if (ynew>=L2+min2)
	  ynew-=L2;
	if (znew<min3)
	  znew+=L3;
	if (znew>=L3+min3)
	  znew-=L1;
	xpos=xnew;
	ypos=ynew;
	zpos=znew;
#endif      
	if((xpos>=min1 && xpos<=min1+L1) && (ypos>=min2 && ypos<=min2+L2) && (zpos>=min3 && zpos<=min3+L3))
	  {
	    ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1+BIN_SHIFT)); // indices of the cell of the particle
	    ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2+BIN_SHIFT));
	    ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3+BIN_SHIFT));
	    i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
	    j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
	    k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
	    ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
	    ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
	    ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
	    ULONG iii= static_cast<ULONG>(fmod(real_prec(i-1+N1),real_prec(N1)));
	    ULONG jjj= static_cast<ULONG>(fmod(real_prec(j-1+N2),real_prec(N2)));
	    ULONG kkk= static_cast<ULONG>(fmod(real_prec(k-1+N3),real_prec(N3)));
	    real_prec xc = static_cast<real_prec>(i)+static_cast<real_prec>(0.5); // centers of the cells
	    real_prec yc = static_cast<real_prec>(j)+static_cast<real_prec>(0.5);
	    real_prec zc = static_cast<real_prec>(k)+static_cast<real_prec>(0.5);
	    real_prec dx = (xpos-min1)/d1 - xc; // distance of particle to center of the cell
	    real_prec dy = (ypos-min2)/d2 - yc;
	    real_prec dz = (zpos-min3)/d3 - zc;
	    real_prec hx0=static_cast<real_prec>(0.75-dx*dx); // fraction of particle assigned
	    real_prec hy0=static_cast<real_prec>(0.75-dy*dy); // fraction of particle assigned
	    real_prec hz0=static_cast<real_prec>(0.75-dz*dz); // fraction of particle assigned
	    real_prec hxp1=static_cast<real_prec>(0.5*(0.5+dx)*(0.5+dx)); // fraction of particle assigned
	    real_prec hyp1=static_cast<real_prec>(0.5*(0.5+dy)*(0.5+dy)); // fraction of particle assigned
	    real_prec hzp1=static_cast<real_prec>(0.5*(0.5+dz)*(0.5+dz)); // fraction of particle assigned
	    real_prec hxm1=static_cast<real_prec>(0.5*(0.5-dx)*(0.5-dx)); // fraction of particle assigned
	    real_prec hym1=static_cast<real_prec>(0.5*(0.5-dy)*(0.5-dy)); // fraction of particle assigned
	    real_prec hzm1=static_cast<real_prec>(0.5*(0.5-dz)*(0.5-dz)); // fraction of particle assigned
	    real_prec mass=num_1;
	    if (weight_prop==_MASS_)
	      mass=Halo[n].mass;
	    else if (weight_prop==_SAT_FRACTION_)
	      mass=Halo[n].number_sub_structures;
	    DELTAtt(i,j,k)    	+= mass*hx0*hy0*hz0;
	    DELTAtt(ii,jj,kk) 	+= mass*hxp1*hyp1*hzp1;
	    DELTAtt(iii,jjj,kkk) 	+= mass*hxm1*hym1*hzm1;
	    DELTAtt(ii,j,k)   	+= mass*hxp1*hy0*hz0;
	    DELTAtt(i,jj,k)   	+= mass*hx0*hyp1*hz0;
	    DELTAtt(i,j,kk)   	+= mass*hx0*hy0*hzp1;
	    DELTAtt(ii,jj,k)  	+= mass*hxp1*hyp1*hz0;
	    DELTAtt(ii,j,kk)  	+= mass*hxp1*hy0*hzp1;
	    DELTAtt(i,jj,kk)  	+= mass*hx0*hyp1*hzp1;
	    DELTAtt(iii,jj,kk)   	+= mass*hxm1*hyp1*hzp1;
	    DELTAtt(ii,jjj,kk)  	+= mass*hxp1*hym1*hzp1;
	    DELTAtt(ii,jj,kkk)  	+= mass*hxp1*hyp1*hzm1;
	    DELTAtt(iii,jjj,kk)  	+= mass*hxm1*hym1*hzp1;
	    DELTAtt(iii,jj,kkk)  	+= mass*hxm1*hyp1*hzm1;
	    DELTAtt(ii,jjj,kkk)  	+= mass*hxp1*hym1*hzm1;
	    DELTAtt(i,jjj,kkk)   	+= mass*hx0*hym1*hzm1;
	    DELTAtt(iii,j,kkk)  	+= mass*hxm1*hy0*hzm1;
	    DELTAtt(iii,jjj,k)  	+= mass*hxm1*hym1*hz0;
	    DELTAtt(i,j,kkk)  	+= mass*hx0*hy0*hzm1;
	    DELTAtt(i,jjj,k)  	+= mass*hx0*hym1*hz0;
	    DELTAtt(iii,j,k)  	+= mass*hxm1*hy0*hz0;
	    DELTAtt(i,jj,kkk)   	+= mass*hx0*hyp1*hzm1;
	    DELTAtt(ii,j,kkk)  	+= mass*hxp1*hy0*hzm1;
	    DELTAtt(iii,jj,k)   	+= mass*hxm1*hyp1*hz0;
	    DELTAtt(ii,jjj,k)  	+= mass*hxp1*hym1*hz0;
	    DELTAtt(i,jjj,kk)   	+= mass*hx0*hym1*hzp1;
	    DELTAtt(iii,j,kk)  	+= mass*hxm1*hy0*hzp1;
	  }	
      }
#ifdef parallel_tsc
#pragma omp critical
    {
      for(ULONG i=0;i<delta.size();i++)
	delta[i]+=private_delta[i];
    }
  }
#endif
  ScreenOutput So;
  So.DONE();
  if(_COUNTS_==weight_prop)
    {
      ULONG count=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:count)
#endif
      for(ULONG i=0;i<delta.size();++i)
	count+=delta[i];
      So.message_screen("Number of objects assigned to grid =",count);
    }
  if(NLOSS!=0) So.message_screen("mass assignment found",NLOSS,"particles outside mesh boundary");
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void cellbound(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, vector<real_prec>&v1, vector<real_prec>&v2, vector<real_prec>&v3)
{
  int NTHREADS=_NTHREADS_;
  omp_set_num_threads(NTHREADS);
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  /* start: interpolate from cell center to cell boundaries */
  vector<real_prec> vx(N,0),vy(N,0),vz(N,0);
  for (ULONG i=0;i<N1; i++)
    for (ULONG j=0; j<N2; j++)
      for (ULONG k=0; k<N3; k++)
	{     
	  ULONG l=k+N3*(j+N2*i);
	  ULONG m=k-1+N3*(j-1+N2*(i-1));
	  if (i>0 && j>0 && k>0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[m]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[m]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[m]+v3[l]);
	    }
	  /* periodic boundary conditions */
	  ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	  ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	  ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	  ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	  ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	  ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	  ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	  if (i==0 && j>0 && k>0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[ii]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[ii]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[ii]+v3[l]);
	    }
	  if (i==0 && j==0 && k>0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[ij]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[ij]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[ij]+v3[l]);
	    }
	  if (i==0 && j>0 && k==0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[ik]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[ik]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[ik]+v3[l]);
	    }
	  if (i==0 && j==0 && k==0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[ijk]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[ijk]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[ijk]+v3[l]);
	    }
	  if (i>0 && j==0 && k==0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[jk]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[jk]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[jk]+v3[l]);
	    }
	  if (i>0 && j==0 && k>0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[jj]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[jj]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[jj]+v3[l]);
	    }
	  if (i>0 && j>0 && k==0)
	    {
	      vx[l]=static_cast<real_prec>(0.5)*(v1[kk]+v1[l]);
	      vy[l]=static_cast<real_prec>(0.5)*(v2[kk]+v2[l]);
	      vz[l]=static_cast<real_prec>(0.5)*(v3[kk]+v3[l]);
	    }
	}
  for (ULONG l=0;l<N; l++)
    {
      v1[l]=vx[l];
      v2[l]=vy[l];
      v3[l]=vz[l];
    }
  /* end: interpolate from cell center to cell boundaries */   
}

void cellboundcomp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&vi)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  /* start: interpolate from cell center to cell boundaries */
  vector<real_prec> viout(N,0);
      
  for (ULONG i=0;i<N1; i++)
    for (ULONG j=0; j<N2; j++)
      for (ULONG k=0; k<N3; k++)
	{     
	  ULONG l=k+N3*(j+N2*i);
	  ULONG m=k-1+N3*(j-1+N2*(i-1));
	      
	  if (i>0 && j>0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[m]+vi[l]);
	    }
	      
	  /* periodic boundary conditions */
	      
	  ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	  ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	  ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	  ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	  ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	  ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	  ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	      
	  if (i==0 && j>0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[ii]+vi[l]);
	    }
	  if (i==0 && j==0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[ij]+vi[l]);
	    }
	  if (i==0 && j>0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[ik]+vi[l]);
	    }
	  if (i==0 && j==0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[ijk]+vi[l]);
	    }
	  if (i>0 && j==0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[jk]+vi[l]);
	    }
	  if (i>0 && j==0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[jj]+vi[l]);
	    }
	  if (i>0 && j>0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[kk]+vi[l]);
	    }
	}
	  
  for (ULONG l=0;l<N; l++)
    {
      vi[l]=viout[l];
    }
  /* end: interpolate from cell center to cell boundaries */   
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void cellbound_comp(ULONG N1,ULONG N2,ULONG N3,vector<real_prec>&vi, vector<real_prec>&viout)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  /* start: interpolate from cell center to cell boundaries */
  for (ULONG i=0;i<N1; i++)
    for (ULONG j=0; j<N2; j++)
      for (ULONG k=0; k<N3; k++)
	{     
	  ULONG l=k+N3*(j+N2*i);
	  ULONG m=k-1+N3*(j-1+N2*(i-1));
	      
	  if (i>0 && j>0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[m]+vi[l]);
	    }
	  /* periodic boundary conditions */
	  ULONG ii=k-1+N3*(j-1+N2*(N1-1));
	  ULONG jj=k-1+N3*(N2-1+N2*(i-1));
	  ULONG kk=N3-1+N3*(j-1+N2*(i-1));
	  ULONG ij=k-1+N3*(N2-1+N2*(N1-1));
	  ULONG ik=N3-1+N3*(j-1+N2*(N1-1));
	  ULONG jk=N3-1+N3*(N2-1+N2*(i-1));
	  ULONG ijk=N3-1+N3*(N2-1+N2*(N1-1));
	  if (i==0 && j>0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[ii]+vi[l]);
	    }
	  if (i==0 && j==0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[ij]+vi[l]);
	    }
	  if (i==0 && j>0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[ik]+vi[l]);
	    }
	  if (i==0 && j==0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[ijk]+vi[l]);
	    }
	  if (i>0 && j==0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[jk]+vi[l]);
	    }
	  if (i>0 && j==0 && k>0)
	    {
	      viout[l]=num_0_5*(vi[jj]+vi[l]);
	    }
	  if (i>0 && j>0 && k==0)
	    {
	      viout[l]=num_0_5*(vi[kk]+vi[l]);
	    }
	}
#pragma omp parallel for
  for (ULONG l=0;l<N; l++)
    {
      vi[l]=viout[l];
    }
  /* end: interpolate from cell center to cell boundaries */   
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void MAS_NEW(s_params_box_mas *params, vector<s_Halo>&Halo,  string weight_prop, vector<real_prec>&delta)
{
  switch (params->masskernel)
    {
    case 0:
      getDensity_NGP(params, Halo, delta,weight_prop);
      break;
    case 1:
      getDensity_CIC(params, Halo,delta,weight_prop);
      break;
    case 2:
      getDensity_TSC(params, Halo,delta,weight_prop);
      break;
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void PoissonSolver(real_prec Lbox, ULONG Nft,  vector<real_prec>&in, vector<real_prec>&out)
{
#ifdef _USE_OMP_
  omp_set_num_threads(_NTHREADS_);
#endif

#ifdef _USE_ZERO_PADDING_POT_
  int factor_ntt=_EXTRA_NFT_FACTOR_*_EXTRA_NFT_FACTOR_*_EXTRA_NFT_FACTOR_;
  ULONG NGRIDnew=factor_ntt*Nft*Nft*Nft;
  vector<real_prec>new_in(NGRIDnew,0);
#pragma omp parallel for
  for(ULONG i=0;i<_EXTRA_NFT_FACTOR_*Nft;++i)
    for(ULONG j=0;j<_EXTRA_NFT_FACTOR_* Nft;++j)
      for(ULONG k=0;k<_EXTRA_NFT_FACTOR_* Nft;++k)
	{
	  ULONG index_padd=index_3d(i,j,k,_EXTRA_NFT_FACTOR_*Nft,_EXTRA_NFT_FACTOR_*Nft);
	  ULONG index_or=index_3d(i,j,k,Nft,Nft);
	  if(index_or<Nft*Nft*Nft)
	    new_in[index_padd]=in[index_or];

	}
  Nft*=_EXTRA_NFT_FACTOR_;
#endif
  ULONG NTT = static_cast<ULONG>(Nft*Nft*(Nft/2+1));
  real_prec deltak=2.*M_PI/Lbox;
#ifdef DOUBLE_PREC
  complex_prec *data_out= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#else
  complex_prec *data_out= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NTT;++i)data_out[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NTT;++i)data_out[i][IMAG]=0;

#ifdef _USE_ZERO_PADDING_POT_
  do_fftw_r2c(Nft,new_in,data_out);
  new_in.clear();new_in.shrink_to_fit();
#else
  do_fftw_r2c(Nft,in,data_out);
#endif


  vector<real_prec> coords(Nft,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft ;++i)
    coords[i]=deltak*(i<=Nft/2? static_cast<real_prec>(i): -static_cast<real_prec>(Nft-i));
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< Nft ; i++)
    for(ULONG j=0;j< Nft ; j++)
      for(ULONG k=0;k< Nft/2+1; k++)
        {
          ULONG ind=index_3d(i,j,k, Nft, Nft/2+1);
          real_prec kmod2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);  // Get k**2
          real_prec term=0.;
          if(kmod2>0.)
            term=-static_cast<real_prec>(1./kmod2);
          data_out[ind][REAL]*=term;   // Divide by k**2
          data_out[ind][IMAG]*=term;
        }
  coords.clear();
  coords.shrink_to_fit();
#ifdef _USE_ZERO_PADDING_POT_
  vector<real_prec> newout(Nft*Nft*Nft,0);
  do_fftw_c2r(Nft ,data_out,newout);
#pragma omp parallel for
  for(ULONG i=0;i<Nft*Nft*Nft/factor_ntt;++i)out[i]=newout[i];
  newout.clear();newout.shrink_to_fit();
#else
  do_fftw_c2r(Nft ,data_out,out);
#endif
#ifdef _DOUBLE_PREC_
  fftw_free(data_out);
#else
  fftwf_free(data_out);
#endif

}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MAS_NGP(real_prec x)
{
  real_prec ans=0;
  if(fabs(x)<0.5)
    ans= 1.0;
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MAS_CIC(real_prec x)
{
  x=fabs(x);
  real_prec ans=0;
  if(x<1)
    ans= 1.-x;
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MAS_TSC(real_prec x)
{
  x=fabs(x);
  real_prec ans=0;
  if(x<0.5)
    ans=(0.75-x*x);
  else if(x>=0.5 && x<1.5)
    {
      real_prec r = 1.5 - x;
      r *= r;
      ans=0.5*r;
    }
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MAS_PCS(real_prec x)
{
  real_prec ans=0;
  real_prec y=fabs(x);
  if(y>=0 && y<1)
    ans=(1./6.)*(4.0- 6.0*x*x + 3.*pow(y,3));
  else if(y>=1 && y<2)
    ans= (1./6.)*pow(2.0-y, 3);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_overdens(const vector<real_prec>&in, real_prec mean, vector<real_prec>&out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\tConverting to Overdensity"<<RESET<<endl;
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/static_cast<double>(mean)-1.0;

#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;

#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_overdens(const vector<real_prec>&in,const vector<real_prec>&weight,vector<real_prec>&out, bool wind_binary)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\tConverting to Overdensity with weights"<<RESET<<endl;
#endif
  double mean=0;
  double veff=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:mean, veff)
#endif
  for(ULONG i=0;i< in.size();++i)
    {
      mean+=static_cast<real_prec>(in[i]);
      if(false==wind_binary)
	veff+=static_cast<real_prec>(weight[i]);
      else
	if (weight[i]>0.)
	  veff++;
    }
  mean/=static_cast<double>(veff);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/mean-weight[i];
#ifdef _FULL_VERBOSE_
  std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_overdens(const vector<real_prec>&in, vector<real_prec>&out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef _VERBOSE_OVERDENS_
  cout<<YELLOW<<"\tConverting to Overdensity"<<RESET<<endl;
#endif

  double mean=get_nobjects(in);
  mean/=static_cast<double>(in.size());

#ifdef _VERBOSE_OVERDENS_
  cout<<YELLOW<<"\tMean =\t"<<mean<<RESET<<endl;
#endif


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i< in.size();++i)
    out[i]=in[i]/mean-num_1;

#ifdef _VERBOSE_OVERDENS_
  std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_nobjects(const vector<real_prec> &density_field)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  real_prec ans=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ans)
#endif
  for(ULONG i=0;i< density_field.size();++i)
    ans+=density_field[i];
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG get_nobjects(const vector<ULONG> &density_field)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG ans=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ans)
#endif
  for(ULONG i=0;i< density_field.size();++i)
    ans+=density_field[i];
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG get_bin(real_prec x, real_prec xmin, ULONG nbins, real_prec delta, bool accumulate)
{
  ULONG ibin  =  static_cast<ULONG>(floor((static_cast<double>(x)-static_cast<double>(xmin))/static_cast<double>(delta)));
  if(ibin==nbins)ibin--;
  if(true==accumulate)
    {
      if (ibin>nbins)ibin=nbins-1; // All values above the maximum are placed in the last bin
      if (ibin<0)ibin=0; // All values below the minimum are placed in the first bin
    }
  return ibin;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec tidal_anisotropy(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  return sqrt(0.5*(pow(lambda3-lambda1,2)+pow(lambda3-lambda2,2)+pow(lambda2-lambda1,2)))/(2.+lambda1+lambda2+lambda3);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec invariant_field_I(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= lambda1+lambda2+lambda3;
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec invariant_field_II(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
  real_prec inv= lambda1;
#else
  real_prec inv= lambda1*lambda2 +lambda2*lambda3 + lambda1*lambda3;
#endif
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec invariant_field_III(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
  real_prec inv= lambda2;
#else
  real_prec inv= (lambda1*lambda2)*lambda3;
#endif
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec invariant_field_IV(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
#ifdef _USE_EIGENVALUES_
  real_prec inv= lambda3;
#else
  real_prec inv=  pow(lambda1,2) + pow(lambda2,2) + pow(lambda3,2);
#endif
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ellipticity(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= (lambda1-lambda3)/(2.*(lambda1+lambda2+lambda3));
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec prolat(real_prec lambda1, real_prec lambda2, real_prec lambda3)
{
  real_prec inv= (lambda1+lambda3-2.*lambda2)/(2.*(lambda1+lambda2+lambda3));
  return inv;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_2d_histogram(real_prec XMin,real_prec XMax, real_prec YMin,real_prec YMax, ULONG Nbins, vector<real_prec>&field_X,vector<real_prec>&field_Y,vector<real_prec>&HIST, bool log_scale){
  ULONG Ngrid=field_X.size();
  if(field_Y.size()!=Ngrid){
    cout<<"ERROR IN "<<__PRETTY_FUNCTION__<<": Meshes do not have the same dimentions"<<endl; exit(0);
  }
  // default is log, fields come normally in delta=rho/mea- 1
  if(true==log_scale)
    {
      real_prec deltaX=(XMax-XMin)/static_cast<double>(Nbins);
      real_prec deltaY=(YMax-YMin)/static_cast<double>(Nbins);
      for(ULONG iX = 0; iX<Ngrid ; ++iX)
	{
	  ULONG ixx=get_bin(log10(1.0+field_X[iX]),XMin,Nbins,deltaX,false);
	  ULONG iyy=get_bin(log10(1.0+field_Y[iX]),YMin,Nbins,deltaY,false);
	  if((ixx<Nbins && iyy <Nbins) && (ixx>=0 && iyy>=0))
	    HIST[index_2d(ixx,iyy,Nbins)]++;
	}
    }
  else
    {
      real_prec deltaX=(XMax-XMin)/static_cast<double>(Nbins);
      real_prec deltaY=(YMax-YMin)/static_cast<double>(Nbins);
      for(ULONG iX = 0; iX<Ngrid ; ++iX)
	{
	  ULONG ixx=get_bin(field_X[iX],XMin,Nbins,deltaX,false);
	  ULONG iyy=get_bin(field_Y[iX],YMin,Nbins,deltaY,false);
	  if((ixx<Nbins && iyy <Nbins) && (ixx>=0 && iyy>=0))
	    HIST[index_2d(ixx,iyy,Nbins)]++;
	}
    }
  real_prec max_hist=get_max(HIST);
  for(ULONG iX = 0; iX<HIST.size() ; ++iX)
    HIST[iX]/=max_hist;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int my_gsl_rng_uniform_(gsl_rng *r, int Nmax){
  real_prec xr=gsl_rng_uniform(r);
  real_prec delta= (1.0)/static_cast<real_prec>(Nmax);
  int val = get_bin(xr, 0., Nmax,delta, true);
  return val;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void rankorder(int seed, vector<real_prec>dens_bins, ULONG Nk, real_prec maxk, real_prec mink, vector<real_prec>&in, vector<real_prec>&pdfin, vector<real_prec>&pdfout)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N=in.size();
  double dk=(maxk-mink)/static_cast<double>(Nk);
  for (ULONG i=0;i<N;i++)
    in[i]=log10(NUM_IN_LOG+static_cast<real_prec>(in[i]));
  const gsl_rng_type *  Ta;
  gsl_rng * ra ;
  gsl_rng_env_setup();
  gsl_rng_default_seed=seed;
  Ta = gsl_rng_ranlux;
  ra = gsl_rng_alloc (Ta);
  for (ULONG i=0;i<N;i++)
    {
      ULONG ik=get_bin(in[i],mink,Nk,dk,true);
      real_prec pdfdi=0.;
      for (ULONG j=0;j<ik+1;j++)
	pdfdi+=pdfin[j];
      real_prec pdfdo=0.;
      ULONG kk=0;
      while (pdfdo<pdfdi)
	{
	  pdfdo+=pdfout[kk];
	  kk++;
	}
      ULONG l=kk;
      if (kk>0)
	l-=1;
      if (kk>Nk)
	l=Nk;
      //* / I defined it to lie a the center of the bin */
      // El factor 10 permite que la pdf final sea más suave
      real_prec xr= static_cast<real_prec>(gsl_rng_uniform(ra));
      real_prec denB=dens_bins[l]+ (xr-0.5)*dk;  //s ?? repartirlo aleatoriamente en este bin. OJO ACA VERIFICARLOS */
      in[i]=pow(10.0, denB)-NUM_IN_LOG; // Convert to the format of the input density field
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void calc_pdf(string type, ULONG N, ULONG Nk, real_prec maxk, real_prec mink, const vector<real_prec>&in, vector<real_prec>&pdfin)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  real_prec dk=(maxk-mink)/static_cast<double>(Nk);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<pdfin.size();++i)
    pdfin[i]=0;

  ULONG ntot=0;

  if(type=="log")
    {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ntot)
#endif
      for(ULONG i=0;i<N;i++)
        {
      real_prec den=log10(NUM_IN_LOG+static_cast<real_prec>(in[i]));
	  if (den>=mink && den<=maxk)
	    {
	      ULONG ik=static_cast<ULONG>(floor((den-mink)/dk));
	      if(ik==Nk)ik--;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
	      pdfin[ik]+=1.0;
	      ntot++;
	    }
        }
    }
  else
    {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:ntot)
#endif
      for(ULONG i=0;i<N;i++)
	{
	  real_prec den=in[i];
	  if(den>=mink && den<=maxk)
	    {
	      ULONG ik=static_cast<ULONG>(floor((den-mink)/dk));
	      if(ik==Nk)ik--;
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
	      pdfin[ik]+=1.0;
	      ntot++;
	    }
	}
    }
  for(ULONG i=0;i<Nk;i++)
    pdfin[i]/=static_cast<double>(ntot);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void EigenValuesTweb(ULONG Nft, real_prec L1,  vector<real_prec> &delta, vector<real_prec> &phi, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  real_prec L2=L1;
  real_prec L3=L1;
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
  ULONG Nhalf=static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft/2+1);
  vector<real_prec> LapPhivx(NGRID,0), LapPhivy(NGRID,0), LapPhivz(NGRID,0);
  vector<real_prec> LapPhivxy(NGRID,0), LapPhivxz(NGRID,0), LapPhivyz(NGRID,0);
  vector<real_prec> LapPhivzx(NGRID,0), LapPhivzy(NGRID,0), LapPhivyx(NGRID,0);
#ifdef _USE_GFINDIFF_EIGENV_
  vector<real_prec> dummy(NGRID,0);
  // Get gradient in x-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,1);
  // Get second derivatives of grad_x
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivx,1); // T_xx
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxy,2);// T_xy
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxz,3);// T_xz
  // Get gradient in y-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,2);
  // Get second derivatives of grad_y
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyx,1);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivy,2);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyz,3);
  // Get gradient in z-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,3);
  // Get second derivatives of grad_z
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzx,1);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzy,2);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivz,3);
  dummy.clear();
  dummy.shrink_to_fit();
#elif defined _USE_GFFT_EIGENV_
#ifdef DOUBLE_PREC
  complex_prec *Deltak= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *Deltak= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
  // Get Tij:
  do_fftw_r2c(Nft,delta,Deltak);  // Get delta(k)
  calc_LapPhiv(Nft,L1,Deltak,LapPhivx,1,1); //Get T11
  calc_LapPhiv(Nft,L1,Deltak,LapPhivxy,1,2);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivxz,1,3);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivyx,2,1);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivy,2,2); // Get T22
  calc_LapPhiv(Nft,L1,Deltak,LapPhivyz,2,3);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivzx,3,1);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivzy,3,2);
  calc_LapPhiv(Nft,L1,Deltak,LapPhivz,3,3);// Get T33
#endif
  ULONG ii=0;
#ifdef _USE_OMP_
#pragma omp parallel private(ii)
  {
#pragma omp parallel for
#endif
    for(ii=0; ii< NGRID;++ii)
      {
	// armamos la matriz
	vector<real_prec>in_matrix(9, 0);
	// THis is the right ordeing for this matrix: T11, T21, T31, T12, ....
	in_matrix[0]=LapPhivx[ii];
	in_matrix[1]=LapPhivyx[ii];
	in_matrix[2]=LapPhivzx[ii];
	in_matrix[3]=LapPhivxy[ii];
	in_matrix[4]=LapPhivy[ii];
	in_matrix[5]=LapPhivzy[ii];
	in_matrix[6]=LapPhivxz[ii];
	in_matrix[7]=LapPhivyz[ii];
	in_matrix[8]=LapPhivz[ii];
	vector<real_prec>eigenv(3,0);
	get_eigen(in_matrix, eigenv);
	out1[ii]=eigenv[0];
	out2[ii]=eigenv[1];
	out3[ii]=eigenv[2];
      }
#ifdef _USE_OMP_
  }//closes parallel region
#endif
  LapPhivx.clear();
  LapPhivx.shrink_to_fit();
  LapPhivy.clear();
  LapPhivy.shrink_to_fit();
  LapPhivz.clear();
  LapPhivz.shrink_to_fit();
  LapPhivxy.clear();
  LapPhivxy.shrink_to_fit();
  LapPhivyx.clear();
  LapPhivyx.shrink_to_fit();
  LapPhivzy.clear();
  LapPhivzy.shrink_to_fit();
  LapPhivyz.clear();
  LapPhivyz.shrink_to_fit();
#ifdef DOUBLE_PREC
  fftw_free(Deltak);
#else
  fftwf_free(Deltak);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void EigenValuesTweb(ULONG Nft, real_prec L1,  vector<real_prec> &delta, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
  ULONG Nhalf=static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft)*static_cast<ULONG>(Nft/2+1);
  vector<real_prec> LapPhivx(NGRID,0), LapPhivy(NGRID,0), LapPhivz(NGRID,0);
  vector<real_prec> LapPhivxy(NGRID,0), LapPhivxz(NGRID,0), LapPhivyz(NGRID,0);
  vector<real_prec> LapPhivzx(NGRID,0), LapPhivzy(NGRID,0), LapPhivyx(NGRID,0);
#ifdef DOUBLE_PREC
  complex_prec *Deltak= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *Deltak= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
  // Get Tij:
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\t\tComputing Tij"<<RESET<<endl;
#endif
  do_fftw_r2c(Nft,delta,Deltak);  // Get delta(k)
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T11"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivx,1,1); //Get T11
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T12"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivxy,1,2);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T13"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivxz,1,3);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T21"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivyx,2,1);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T22"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivy,2,2); // Get T22
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T23"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivyz,2,3);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T31"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivzx,3,1);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T32"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivzy,3,2);
#ifdef _FULL_VERBOSE_
//  cout<<YELLOW<<"\t\t\t\t T33"<<RESET<<endl;
#endif
  calc_LapPhiv(Nft,L1,Deltak,LapPhivz,3,3);// Get T33
#ifdef DOUBLE_PREC
  fftw_free(Deltak);
#else
  fftwf_free(Deltak);
#endif
  ULONG ii=0;
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\t\tComputing lambdas"<<RESET<<endl;
#endif
#ifdef _USE_OMP_
#pragma omp parallel for private(ii)
#endif
  for(ii=0; ii< NGRID;++ii)
    {
      // armamos la matriz
      vector<real_prec>in_matrix(9, 0);
      // This is the right ordering for this matrix: T11, T21, T31, T12, ....
      in_matrix[0]=LapPhivx[ii];
      in_matrix[1]=LapPhivyx[ii];
      in_matrix[2]=LapPhivzx[ii];
      in_matrix[3]=LapPhivxy[ii];
      in_matrix[4]=LapPhivy[ii];
      in_matrix[5]=LapPhivzy[ii];
      in_matrix[6]=LapPhivxz[ii];
      in_matrix[7]=LapPhivyz[ii];
      in_matrix[8]=LapPhivz[ii];
      vector<real_prec>eigenv(3,0);
      get_eigen(in_matrix, eigenv);//This gives only eigenvalues
      out1[ii]=eigenv[0];
      out2[ii]=eigenv[1];
      out3[ii]=eigenv[2];
      in_matrix.clear();in_matrix.shrink_to_fit();
      eigenv.clear();eigenv.shrink_to_fit();
    }
  LapPhivx.clear();
  LapPhivx.shrink_to_fit();
  LapPhivy.clear();
  LapPhivy.shrink_to_fit();
  LapPhivz.clear();
  LapPhivz.shrink_to_fit();
  LapPhivxy.clear();
  LapPhivxy.shrink_to_fit();
  LapPhivxz.clear();
  LapPhivxz.shrink_to_fit();
  LapPhivyx.clear();
  LapPhivyx.shrink_to_fit();
  LapPhivzx.clear();
  LapPhivzx.shrink_to_fit();
  LapPhivzy.clear();
  LapPhivzy.shrink_to_fit();
  LapPhivyz.clear();
  LapPhivyz.shrink_to_fit();
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void EigenValuesTweb_bias(ULONG Nft, real_prec L1, const vector<real_prec> &delta, const vector<real_prec> &phi, vector<real_prec> &S2, vector<real_prec> &S3, vector<real_prec> &N2D)
{
  real_prec L2=L1;
  real_prec L3=L2;
#ifdef _USE_OMP_
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
#endif
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
  vector<real_prec> LapPhivx(NGRID,0), LapPhivy(NGRID,0), LapPhivz(NGRID,0);
  vector<real_prec> LapPhivxy(NGRID,0), LapPhivxz(NGRID,0), LapPhivyz(NGRID,0);
  vector<real_prec> LapPhivzx(NGRID,0), LapPhivzy(NGRID,0), LapPhivyx(NGRID,0);
  vector<real_prec> dummy(NGRID,0);
  // Get gradient in x-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,1);
  // Get second derivatives of grad_x
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivx,1);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxy,2);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivxz,3);
  // Get gradient in y-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,2);
  // Get second derivatives of grad_y
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivy,2);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyx,1);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivyz,3);
  // Get gradient in z-direction
  gradfindif(Nft,Nft,Nft,L1,L2,L3,phi,dummy,3);
  // Get second derivatives of grad_z
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivz,3);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzx,1);
  gradfindif(Nft,Nft,Nft,L1,L2,L3,dummy,LapPhivzy,2);
#ifdef _USE_S2_
  // Get S2:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;++i)
    {
      real_prec S11=LapPhivx[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
      real_prec S12=LapPhivxy[i];
      real_prec S13=LapPhivxz[i];
      real_prec S21=LapPhivyx[i];
      real_prec S22=LapPhivy[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
      real_prec S23=LapPhivyz[i];
      real_prec S31=LapPhivzx[i];
      real_prec S32=LapPhivzy[i];
      real_prec S33=LapPhivz[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
      S2[i]= S11 *S11 +S12*S21 + S13*S31+ S21*S12 +S22*S22 + S23*S32 +S31*S13 +S23*S32 +S33*S33;
    }
#endif
#ifdef _USE_S3_
  // Get S3:
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;++i)
    {
      real_prec S11=LapPhivx[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
      real_prec S12=LapPhivxy[i];
      real_prec S13=LapPhivxz[i];
      real_prec S21=LapPhivyx[i];
      real_prec S22=LapPhivy[i]-(1./3.)*(La1[i]+La2[i]+La3[i]);
      real_prec S23=LapPhivyz[i];
      real_prec S31=LapPhivzx[i];
      real_prec S32=LapPhivzy[i];
      real_prec S33=LapPhivz[i]-(1./3.)*(La1[i]+La2[i]+La2[i]);
      real_prec A1=S11*(S11*S11+S12*S21+S12*S31)+S21*(S11*S12+S12*S21+S12*S31)+S31*(S11*S13+S12*S23+S13*S33);
      real_prec A2=S12*(S21*S11+S22*S21+S23*S31)+S22*(S21*S12+S22*S22+S23*S32)+S32*(S21*S13+S22*S23+S23*S33);
      real_prec A3=S13*(S31*S11+S32*S21+S33*S31)+S23*(S31*S12+S32*S22+S33*S32)+S33*(S31*S31+S32*S23+S33*S33);
      S3[i]=A1+A2+A3;
    }
#endif
#if defined (_USE_S3_) || defined (_USE_S3_)
  La1.clear(); La1.shrink_to_fit();
  La2.clear(); La2.shrink_to_fit();
  La3.clear(); La3.shrink_to_fit();
#endif
#ifdef _USE_NABLA2DELTA_
  vector<real_prec> NLapPhivx(NGRID,0), NLapPhivy(NGRID,0), NLapPhivz(NGRID,0);
  // Get partial delta/partial x
  gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivx,1);
  // Get partial delta/partial y
  gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivy,2);
  // Get partial delta/partial z
  gradfindif(Nft,Nft,Nft,L1,L2,L3,delta,NLapPhivz,3);
  // Get partial( delta/partial x) /partial x
  gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivx,NLapPhivx,1);
  // Get partial( delta/partial y) /partial y
  gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivy,NLapPhivy,2);
  // Get partial( delta/partial z) /partial z
  gradfindif(Nft,Nft,Nft,L1,L2,L3,NLapPhivz,NLapPhivz,3);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;++i)
    N2D[i]=NLapPhivx[i]+NLapPhivy[i]+NLapPhivz[i];
  NLapPhivx.clear();
  NLapPhivx.shrink_to_fit();
  NLapPhivy.clear();
  NLapPhivy.shrink_to_fit();
  NLapPhivz.clear();
  NLapPhivz.shrink_to_fit();

#endif
  dummy.clear();
  dummy.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void EigenValuesVweb(ULONG Nft, real_prec L1, vector<real_prec>&Vinx, vector<real_prec>&Viny,vector<real_prec>&Vinz, vector<real_prec>&diver, vector<real_prec> &out1, vector<real_prec> &out2, vector<real_prec> &out3)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
    real_prec L2=L1;
    real_prec L3=L2;
    ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
    vector<real_prec> vxy(NGRID,0),vyx(NGRID,0), vzy(NGRID,0);
    vector<real_prec> vxz(NGRID,0),vyz(NGRID,0), vzx(NGRID,0);
    // Sheer tensor, 0.5(dv_i/ dx_j + dv_j/dx_i), symmetric, only 6 independent components
    vector<real_prec> Svxx(NGRID,0), Svyy(NGRID,0),Svzz(NGRID,0);
    vector<real_prec> Svxy(NGRID,0), Svxz(NGRID,0),Svyz(NGRID,0);
    // Get gradient in x-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,Svxx,1); // dVx/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,vxy,2); //dVx/dy
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinx,vxz,3); // dVx/dz
    // Get gradient in y-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,Svyy,2); // dVy/dy
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,vyx,1);// dVy/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Viny,vyz,3);// dVy/dz
    // Get gradient in z-direction
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,Svzz,3); // dVz/dz
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,vzx,1); // dVz/dx
    gradfindif(Nft,Nft,Nft,L1,L2,L3,Vinz,vzy,2); // dVz/dy
    // Get divergence
    #ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0; i<NGRID;++i)
        diver[i]=Svxx[i]+Svyy[i]+Svzz[i];
    // Get shear
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG i=0;i<NGRID;i++)
      {
        Svxy[i]=-0.5*(vxy[i]+vyx[i]);
        Svxz[i]=-0.5*(vxz[i]+vzx[i]);
        Svyz[i]=-0.5*(vyz[i]+vzy[i]);
      }
    vxy.clear();
    vxy.shrink_to_fit();
    vxz.clear();
    vxz.shrink_to_fit();
    vyx.clear();
    vyx.shrink_to_fit();
    vzx.clear();
    vzx.shrink_to_fit();
    vzy.clear();
    vzy.shrink_to_fit();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(ULONG ii=0;ii<NGRID;ii++)
      {
        vector<real_prec>in_matrix(9, 0);
        // This is the right ordering for this matrix: T11, T21, T31, T12, ....
        in_matrix[0]=Svxx[ii];
        in_matrix[1]=Svxy[ii];
        in_matrix[2]=Svxz[ii];
        in_matrix[3]=Svxy[ii];
        in_matrix[4]=Svyy[ii];
        in_matrix[5]=Svyz[ii];
        in_matrix[6]=Svxz[ii];
        in_matrix[7]=Svyz[ii];
        in_matrix[8]=Svzz[ii];
        vector<real_prec>eigenv(3,0);
        get_eigen(in_matrix, eigenv);//This gives only eigenvalues
        out1[ii]=eigenv[0];
        out2[ii]=eigenv[1];
        out3[ii]=eigenv[2];
        in_matrix.clear();in_matrix.shrink_to_fit();
        eigenv.clear();eigenv.shrink_to_fit();
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void  get_high_res_id_from_low_res_id(int Nft_high, int Nft_low,vector<ULONG>&id_info){
  double delta_low=static_cast<double>(Nft_high)/static_cast<double>(Nft_low);
  for(int i=0;i<Nft_high ;++i)
    for(int j=0;j<Nft_high ;++j)
      for(int k=0;k<Nft_high ;++k)
	{
	  // now get the coords of the i, j , k in the lwo res:
	  ULONG il = static_cast<ULONG>(floor((i+0.5)/delta_low)); // indices of the cell of the particle
	  ULONG jl = static_cast<ULONG>(floor((j+0.5)/delta_low));
	  ULONG kl = static_cast<ULONG>(floor((k+0.5)/delta_low));
	  ULONG id_h   = index_3d(i, j, k, Nft_high, Nft_high);
	  id_info[id_h]= index_3d(il,jl,kl,Nft_low,Nft_low);
	}
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void  get_low_res_id_from_high_res_id(int Nft, int Nft_low, vector<s_cell_info_reduced>&cell_info_low){
  double delta_low=static_cast<double>(Nft)/static_cast<double>(Nft_low);
  for(int i=0;i<Nft ;++i)
    for(int j=0;j<Nft ;++j)
      for(int k=0;k<Nft ;++k)
	{
	  // now get the coords of the i, j , k in the lwo res:
	  ULONG il = static_cast<ULONG>(floor((i+0.5)/delta_low)); // indices of the cell of the particle
	  ULONG jl = static_cast<ULONG>(floor((j+0.5)/delta_low));
	  ULONG kl = static_cast<ULONG>(floor((k+0.5)/delta_low));
	  ULONG id_l = index_3d(il,jl,kl,Nft_low,Nft_low);
	  ULONG id   = index_3d(i, j, k, Nft, Nft);
	  cell_info_low[id_l].gal_index.push_back(id);  // allocate for each lowres ID all those high_res id living inside
	}
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_neighbour_cells(ULONG Nft, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#if defined _USE_NEIGHBOURS_ || defined _GET_DIST_MIN_SEP_REF_ || defined _GET_DIST_MIN_SEP_MOCK_
  int max_neigh_per_dim = 2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the centrak cell;
  int N_Neigh_cells=pow(max_neigh_per_dim+1,2)+ pow(max_neigh_per_dim+1,2)  + (pow(max_neigh_per_dim+1,2)); // total number of cells to explore around a cell, excluding the central
#else
  int max_neigh_per_dim =2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  int N_Neigh_cells=pow(max_neigh_per_dim,2)+ pow(max_neigh_per_dim,2)  + (pow(max_neigh_per_dim,2)-1); // total number of cells to explore around a cell, excluding the central
#endif
  vector<int> index_cells_dime(max_neigh_per_dim,0);
  for(int i=0;i< max_neigh_per_dim ;++i)
    index_cells_dime[i]= i-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
  for(ULONG i=0;i<Nft ;++i)
    for(ULONG j=0;j<Nft ;++j)
      for(ULONG k=0;k<Nft ;++k)
        {
          ULONG ID=index_3d(i,j,k,Nft,Nft);//get ID of cell
          nearest_cells_to_cell[ID].close_cell.resize(N_Neigh_cells,0);
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
          nearest_cells_to_cell[ID].bc_x.resize(N_Neigh_cells,0);
          nearest_cells_to_cell[ID].bc_y.resize(N_Neigh_cells,0);
          nearest_cells_to_cell[ID].bc_z.resize(N_Neigh_cells,0);
#endif
          ULONG count=0;
          for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
            {
              int new_ni = i - index_cells_dime[idx];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
	      int aux_bc_x=0;
#endif
	      if(new_ni<0)
		{
		  new_ni+= Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		  aux_bc_x=-1;
#endif
		}
	      if(new_ni>= Nft)
		{
		  new_ni-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		  aux_bc_x=1;
#endif
		}
	      for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
		{
		  int new_nj = j - index_cells_dime[idy];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		  int aux_bc_y=0;
#endif
		  if(new_nj<0)
		    {
		      new_nj+=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		      aux_bc_y=-1;
#endif
		    }//si la celda esta por detrñas, ponga -1
		  if(new_nj>= Nft)
		    {
		      new_nj-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		      aux_bc_y=1;
#endif
		    }
		  for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
		    {
		      int new_nk = k - index_cells_dime[idz];
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
		      int aux_bc_z=0;
#endif
		      if(new_nk<0)
			{
			  new_nk+=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
			  aux_bc_z=-1;
#endif
			}
		      if(new_nk>=Nft)
			{
			  new_nk-=Nft;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
			  aux_bc_z=1;
#endif
			}
		      ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
		      if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
			{
			  nearest_cells_to_cell[ID].close_cell[count]=new_id;
#if defined _USE_NEIGHBOURS_|| defined (_GET_DIST_MIN_SEP_REF_) || defined (_GET_DIST_MIN_SEP_MOCK_)
			  nearest_cells_to_cell[ID].bc_x[count]=aux_bc_x;
			  nearest_cells_to_cell[ID].bc_y[count]=aux_bc_y;
			  nearest_cells_to_cell[ID].bc_z[count]=aux_bc_z;
#endif
			  count++;
			}
		    }
		}
	    }
        }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_neighbour_cells_of_cell(ULONG Nft, vector<ULONG>&cell_index, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  int max_neigh_per_dim =2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  //  ULONG N_Neigh_cells=pow(max_neigh_per_dim,2)+ pow(max_neigh_per_dim,2)  + (pow(max_neigh_per_dim,2)-1); // total number of cells to explore around a cell, excluding the central
  ULONG N_Neigh_cells=pow(max_neigh_per_dim,3)-1; // total number of cells to explore around a cell, excluding the central
  vector<ULONG> index_cells_dime(max_neigh_per_dim,0);
  for(int i=0;i< max_neigh_per_dim ;++i)
    index_cells_dime[i]= i-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded
  for(ULONG id=0;id<cell_index.size();++id)
    {
      ULONG ID=cell_index[id];
      ULONG i=0, j=0, k=0;
      index2coords(ID,Nft,i,j,k);
      nearest_cells_to_cell[id].close_cell.resize(N_Neigh_cells,0);
      ULONG count=0;
      for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
        {
          int new_ni = i - index_cells_dime[idx];
          if(new_ni<0)
            new_ni+= Nft;
          if(new_ni>= Nft)
            new_ni-=Nft;
          for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
	    {
	      int new_nj = j - index_cells_dime[idy];
	      if(new_nj<0)
		new_nj+=Nft;//si la celda esta por detras, ponga -1
	      if(new_nj>= Nft)
		new_nj-=Nft;
	      for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
		{
		  int new_nk = k - index_cells_dime[idz];
		  if(new_nk<0)
		    new_nk+=Nft;
		  if(new_nk>=Nft)
		    new_nk-=Nft;
		  ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
		  if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
		    {
		      nearest_cells_to_cell[id].close_cell[count]=new_id;
		      count++;
                    }
		}
            }
	}
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_neighbour_cells_cat_analyze(int Nft, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell, bool exclude){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG max_neigh_per_dim = 2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  vector<int> index_cells_dime(max_neigh_per_dim,0);
  for(int i=0;i< max_neigh_per_dim ;++i)
    index_cells_dime[i]= i-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded
  // total number of cells to explore around a cell, including the central only once as it should be
  ULONG N_Neigh_cells=pow(max_neigh_per_dim,3);
  if(true==exclude)
    N_Neigh_cells-=1;  // subtract 1 to avoid counting 3 (one per dimention) times the central cell, which has been already counted once.
  ULONG NGRID = Nft*Nft*Nft;

#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
  for(int i=0;i<Nft ;++i)
    for(int j=0;j<Nft ;++j)
      for(int k=0;k<Nft ;++k)
        {
          ULONG ID=index_3d(i,j,k,Nft,Nft);//get ID of cell
          ULONG count=0;
          for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
            {
              int new_ni = i - index_cells_dime[idx];
	      int aux_bc_x=0;
	      if(new_ni<0)
		{
		  new_ni+= Nft;
		  aux_bc_x=-1;
		}
	      if(new_ni>= Nft)
		{
		  new_ni-=Nft;
		  aux_bc_x=1;
		}
	      for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
		{
		  int new_nj = j - index_cells_dime[idy];
		  int aux_bc_y=0;
		  if(new_nj<0)
		    {
		      new_nj+=Nft;
		      aux_bc_y=-1;
		    }//si la celda esta por detrñas, ponga -1
		  if(new_nj>= Nft)
		    {
		      new_nj-=Nft;
		      aux_bc_y=1;
		    }

		  for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
		    {
		      int new_nk = k - index_cells_dime[idz];
		      int aux_bc_z=0;
		      if(new_nk<0)
			{
			  new_nk+=Nft;
			  aux_bc_z=-1;
			}
		      if(new_nk>=Nft)
			{
			  new_nk-=Nft;
			  aux_bc_z=1;
			}

		      ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
		      if(exclude==true)
                        {
                          if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
			    {
                              nearest_cells_to_cell[ID].close_cell.push_back(new_id);
                              nearest_cells_to_cell[ID].bc_x.push_back(aux_bc_x);
                              nearest_cells_to_cell[ID].bc_y.push_back(aux_bc_y);
                              nearest_cells_to_cell[ID].bc_z.push_back(aux_bc_z);
                              count++;
			    }
                        }
		      else
                        {
			  nearest_cells_to_cell[ID].close_cell.push_back(new_id);
			  nearest_cells_to_cell[ID].bc_x.push_back(aux_bc_x);
			  nearest_cells_to_cell[ID].bc_y.push_back(aux_bc_y);
			  nearest_cells_to_cell[ID].bc_z.push_back(aux_bc_z);
			  count++;
                        }
		    }
		}
	    }
	}
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_neighbour_cells_cat_analyze_chuncks(int Nchunck, int Nft, int N_cells_back_forth,vector<s_nearest_cells>&nearest_cells_to_cell, bool exclude){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG max_neigh_per_dim = 2*N_cells_back_forth+1; // Number of neighbouring cells pr dimension, including the central cell;
  vector<int> index_cells_dime(max_neigh_per_dim,0);
  for(int ic=0;ic< max_neigh_per_dim ;++ic)
    index_cells_dime[ic]= ic-N_cells_back_forth;  // e.g., 1, 0, -1 if N_cells_back_forth=1 and the central excluded
  // total number of cells to explore around a cell, including the central only once as it should be
  ULONG N_Neigh_cells=pow(max_neigh_per_dim,3);
  if(true==exclude)
    N_Neigh_cells-=1;  // subtract 1 to avoid counting 3 (one per dimention) times the central cell, which has been already counted once.
  if(true==exclude)
    {
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
      // The loop over ic
 for(int ic=0; ic<chunck_factor ;++ic) // this loop is not necessary if chunck_factor=1
	for(int j=0;j<Nft ;++j)
	  for(int k=0;k<Nft ;++k)
	    {
	      // The variable "i" goes from 0 to Nft-1 when the full loop in Catalog is done;
	      // if chunck_factor is 2, this is 0,2,4,6,8...255.
	      int i=chunck_factor*Nchunck + ic;
	      ULONG ID=index_3d(ic,j,k,Nft,Nft);//get ID of cell
	      // ULONG count=0;
	      for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
		{
		  int new_ni = i - index_cells_dime[idx];
		  int aux_bc_x=0;
		  if(new_ni<0)
		    {
		      new_ni+= Nft;
		      aux_bc_x=-1;
		    }
		  if(new_ni>= Nft)
		    {
		      new_ni-=Nft;
		      aux_bc_x=1;
		    }
		  for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
		    {
		      int new_nj = j - index_cells_dime[idy];
		      int aux_bc_y=0;
		      if(new_nj<0)
			{
			  new_nj+=Nft;
			  aux_bc_y=-1;
			}//si la celda esta por detrñas, ponga -1
		      if(new_nj>= Nft)
			{
			  new_nj-=Nft;
			  aux_bc_y=1;
			}
		      for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
			{
			  int new_nk = k - index_cells_dime[idz];
			  int aux_bc_z=0;
			  if(new_nk<0)
			    {
			      new_nk+=Nft;
			      aux_bc_z=-1;
			    }
			  if(new_nk>=Nft)
			    {
			      new_nk -= Nft;
			      aux_bc_z = 1;
			    }
			  ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
			  if(new_id!=ID) // This explicitely excludes the ID of the same cell where we are now.
			    {
			      nearest_cells_to_cell[ID].close_cell.push_back(new_id);
			      nearest_cells_to_cell[ID].bc_x.push_back(aux_bc_x);
			      nearest_cells_to_cell[ID].bc_y.push_back(aux_bc_y);
			      nearest_cells_to_cell[ID].bc_z.push_back(aux_bc_z);
			      //count++;
			    }
			}
		    }
		}
	    }

    }
  else  // Explicitely separated to avoid if's in the deeper part of the loops
    {
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
  for(int ic=0; ic<chunck_factor ;++ic)
  	for(int j=0;j<Nft ;++j)
	    for(int k=0;k<Nft ;++k)
	      {
          ULONG ID=index_3d(ic,j,k,Nft,Nft);//get ID of cell
          int i=chunck_factor*Nchunck + ic;// if chunck_factor is 2, this is 0,2,4,6,8...254
          for(int idx =0; idx < max_neigh_per_dim; ++idx)  // loop in the x-direction
          {
            int new_ni = i - index_cells_dime[idx];
            int aux_bc_x=0;
            if(new_ni<0)
              {
                new_ni+= Nft;
                aux_bc_x=-1;
              }
            if(new_ni>= Nft)
              {
                new_ni-=Nft;
                aux_bc_x=1;
              }
              for(int idy =0; idy < max_neigh_per_dim; ++idy) // loop in the y-direction
                {
                  int new_nj = j - index_cells_dime[idy];
                  int aux_bc_y=0;
                  if(new_nj<0)
                    {
                      new_nj+=Nft;
                      aux_bc_y=-1;
                    }//si la celda esta por detrñas, ponga -1
                 if(new_nj>= Nft)
                  {
                    new_nj-=Nft;
                  aux_bc_y=1;
                  }
                  for(int idz =0; idz < max_neigh_per_dim; ++idz) // loop in the z-direction
                  {
                    int new_nk = k - index_cells_dime[idz];
                    int aux_bc_z=0;
                    if(new_nk<0)
                      {
                        new_nk+=Nft;
                        aux_bc_z=-1;
                      }
                    if(new_nk>=Nft)
                      {
                        new_nk-=Nft;
                        aux_bc_z=1;
                      }
                    ULONG new_id=index_3d(new_ni,new_nj,new_nk,Nft,Nft); // Get the ID of the neighbouring cells;
                    nearest_cells_to_cell[ID].close_cell.push_back(new_id);
                    nearest_cells_to_cell[ID].bc_x.push_back(aux_bc_x);
                    nearest_cells_to_cell[ID].bc_y.push_back(aux_bc_y);
                    nearest_cells_to_cell[ID].bc_z.push_back(aux_bc_z);
                  }
		           }
		        }
	      }
   }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_vel_field_using_neigh_cells(vector<ULONG>&id_empty_cells,vector<s_nearest_cells>&nearest_cells,vector<real_prec>&velf){
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG j=0;j<id_empty_cells.size();j++)
    {
      ULONG i=id_empty_cells[j]; // ID f emty cell
      real_prec vx_t=0;
      ULONG nc=nearest_cells[j].close_cell.size();
      for(ULONG k=0;k<nc;++k) // loop over the neighbour cells if the cell with ID empty_cells[i]
	vx_t+=velf[nearest_cells[j].close_cell[k]]; // Assign the average of the neighbour cells
      velf[i]= nc>0? vx_t/static_cast<real_prec>(nc) : 0;   //ojo con esto
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
inline real_prec Sq(int exp, real_prec x){
  real_prec ans;
  //////////////////////////////////////////////////////////
  // Correcting for the MAS mode-by-mode
  //////////////////////////////////////////////////////////
  switch(exp)
    {
    case(0):
      ans=1; break;
    case(1):
      ans=(1.0-(2./3.)*x*x); break;
    case(2):
      ans=(1.0-x*x+(2./15.)*x*x*x*x); break;
    case(3):
      ans=(1.0-(4./3.)*x*x+(2./5.)*x*x*x*x-(4./315.)*x*x*x*x*x*x); break;}
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec SN_correction_MAS(int exp,real_prec xx,real_prec yy,real_prec zz)
{
  return Sq(exp,xx)*Sq(exp,yy)*Sq(exp,zz);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec correction_MAS(int expo,real_prec xx,real_prec yy,real_prec zz)
{
  real_prec res = xx*yy*zz;
  real_prec ans=pow(res, expo+1);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec correction_MAS(ULONG Nft, int expo, int i,int j,int k)
{
  // Correcting for the MAS mode-by-mode
  real_prec xx=i*M_PI/Nft;
  real_prec yy=j*M_PI/Nft;
  real_prec zz=k*M_PI/Nft;
  real_prec sx=(i==0? 1.00: sin(xx)/xx);
  real_prec sy=(j==0? 1.00: sin(yy)/yy);
  real_prec sz=(k==0? 1.00: sin(zz)/zz);
  real_prec res = sx*sy*sz;
  real_prec ans = pow(res, expo+1);
  return ans;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_mean_density_interpolated(real_prec redshift_min_sample, real_prec redshift_max_sample, real_prec ra, real_prec dec, real_prec redshift, vector< vector<real_prec> > &dndz_m, real_prec *nb)
{
  long ipix=0;
#ifdef HEALPIX
  Healpix_Map<real_prec>map(log2(this->nside), RING);
  real_prec fac=M_PI/180.0;
  point.phi=ra*fac;
  point.theta=0.5*M_PI-dec*fac;
  ipix=map.ang2pix(point);
#endif
  real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/dndz_m.size();
  int iz=(int)floor((float)((redshift-redshift_min_sample)/Delta_Z));
  *nb=dndz_m[iz][ipix];
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void grid_assignment_NGP(Params *params, real_prec x, real_prec y, real_prec z,real_prec weight, vector<real_prec>& field)
{
  // Interpolation of the galaxy overdensity field
  // We apply periodic bounday conditions to remap objects
  // with coords. outside the range [0,Lside] within the box.
  // Same structure as TSC, only MAS(x) is different
  x -= (params->_Xoffset() - 0.5*params->_Lbox());
  y -= (params->_Yoffset() - 0.5*params->_Lbox());
  z -= (params->_Zoffset() - 0.5*params->_Lbox());
  // Apply boundary conditions to map tracers inside the box
  // Notice that now mins are = 0
  while(x< 0){x+=params->_Lbox();}
  while(y< 0){y+=params->_Lbox();}
  while(z< 0){z+=params->_Lbox();}
  while(x>=params->_Lbox()){x-=params->_Lbox();}
  while(y>=params->_Lbox()){y-=params->_Lbox();}
  while(z>=params->_Lbox()){z-=params->_Lbox();}
  // Determining the number of cell in each direction:   *
  // Output in the range [0,N-1]                         *
  // Notice that now mins are = 0
  ULONG xc=get_bin(x,0, params->_Nft(), params->_d_delta_x(), true);
  ULONG yc=get_bin(y,0, params->_Nft(), params->_d_delta_y(), true);
  ULONG zc=get_bin(z,0, params->_Nft(), params->_d_delta_z(), true);
  // This is just to be extremily sure that we do not go beyond a cell with index larger
  // than the number of cells per dimention. It should never happen by construction.
  if(xc==params->_Nft())
    xc-=1;
  if(yc==params->_Nft())
    yc-=1;
  if(zc==params->_Nft())
    zc-=1;
  ULONG index_aux = index_3d(xc, yc, zc,params->_Nft(),params->_Nft());
  field[index_aux]+= weight;
#ifdef _MASS_WEIGHT_POWER_
  field_mw[index_aux] +=weight_aux*weight_mass;
#endif
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void grid_assignment_CIC(Params *params,real_prec x, real_prec y, real_prec z,real_prec weight, vector<real_prec>& field)
{
  x -= (params->_Xoffset() - 0.5*params->_Lbox());
  y -= (params->_Yoffset() - 0.5*params->_Lbox());
  z -= (params->_Zoffset() - 0.5*params->_Lbox());
  while(x< 0){x+=params->_Lbox();}
  while(y< 0){y+=params->_Lbox();}
  while(z< 0){z+=params->_Lbox();}
  while(x>=params->_Lbox()){x-=params->_Lbox();}
  while(y>=params->_Lbox()){y-=params->_Lbox();}
  while(z>=params->_Lbox()){z-=params->_Lbox();}
  ULONG xc=get_bin(x,0, params->_Nft(), params->_d_delta_x(), true);
  ULONG yc=get_bin(y,0, params->_Nft(), params->_d_delta_y(), true);
  ULONG zc=get_bin(z,0, params->_Nft(), params->_d_delta_z(), true);
  if(xc==params->_Nft())
    xc-=1;
  if(yc==params->_Nft())
    yc-=1;
  if(zc==params->_Nft())
    zc-=1;
  real_prec xx  = params->_d_delta_x()*(static_cast<real_prec>(xc)+0.5);
  real_prec yy  = params->_d_delta_y()*(static_cast<real_prec>(yc)+0.5);
  real_prec zz  = params->_d_delta_z()*(static_cast<real_prec>(zc)+0.5);
  real_prec xxf = params->_d_delta_x()*(static_cast<real_prec>(xc)+1.5);
  real_prec yyf = params->_d_delta_y()*(static_cast<real_prec>(yc)+1.5);
  real_prec zzf = params->_d_delta_z()*(static_cast<real_prec>(zc)+1.5);
  real_prec xxb = params->_d_delta_x()*(static_cast<real_prec>(xc)-0.5);
  real_prec yyb = params->_d_delta_y()*(static_cast<real_prec>(yc)-0.5);
  real_prec zzb = params->_d_delta_z()*(static_cast<real_prec>(zc)-0.5);
  ULONG Xb=(xc==0 ? params->_Nft(): xc);
  ULONG Yb=(yc==0 ? params->_Nft(): yc);
  ULONG Zb=(zc==0 ? params->_Nft(): zc);
  ULONG Xf=(xc==params->_Nft()-1 ? -1: xc);
  ULONG Yf=(yc==params->_Nft()-1 ? -1: yc);
  ULONG Zf=(zc==params->_Nft()-1 ? -1: zc);
  ULONG i_idx[] = {Xb-1, xc, Xf+1};
  ULONG j_idx[] = {Yb-1, yc, Yf+1};
  ULONG k_idx[] = {Zb-1, zc, Zf+1};
  vector<real_prec> MAS_xx =
    {
      MAS_CIC((xxb  - x)/static_cast<double>(params->_d_delta_x())),
      MAS_CIC((xx   - x)/static_cast<double>(params->_d_delta_x())),
      MAS_CIC((xxf  - x)/static_cast<double>(params->_d_delta_x())),
    };
  vector<real_prec> MAS_yy =
    {
      MAS_CIC((yyb  - y)/static_cast<double>(params->_d_delta_y())),
      MAS_CIC((yy   - y)/static_cast<double>(params->_d_delta_y())),
      MAS_CIC((yyf  - y)/static_cast<double>(params->_d_delta_y())),
    };
  vector<real_prec> MAS_zz =
    {
      MAS_CIC((zzb  - z)/static_cast<double>(params->_d_delta_z())),
      MAS_CIC((zz   - z)/static_cast<double>(params->_d_delta_z())),
      MAS_CIC((zzf  - z)/static_cast<double>(params->_d_delta_z())),
    };
  for(int ih=0;ih<MAS_xx.size();++ih)
    for(int jh=0;jh<MAS_yy.size();++jh)
      for(int kh=0;kh<MAS_zz.size();++kh)
      	{
          ULONG index_aux = index_3d(i_idx[ih], j_idx[jh], k_idx[kh],params->_Nft(),params->_Nft());
          field[index_aux]+=  weight*MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
#ifdef _MASS_WEIGHT_POWER_
      	  field_mw[index_aux] +=weight_aux*weight_mass;
#endif
	      }
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void grid_assignment_TSC(Params *params,real_prec x, real_prec y, real_prec z,real_prec weight, vector<real_prec>& field)
{
  x -= params->_Xoffset() - 0.5*params->_Lbox();
  y -= params->_Yoffset() - 0.5*params->_Lbox();
  z -= params->_Zoffset() - 0.5*params->_Lbox();
  while(x< 0){x+=params->_Lbox();}
  while(y< 0){y+=params->_Lbox();}
  while(z< 0){z+=params->_Lbox();}
  while(x>=params->_Lbox()){x-=params->_Lbox();}
  while(y>=params->_Lbox()){y-=params->_Lbox();}
  while(z>=params->_Lbox()){z-=params->_Lbox();}
  ULONG xc=get_bin(x,0, params->_Nft(), params->_d_delta_x(), false);
  ULONG yc=get_bin(y,0, params->_Nft(), params->_d_delta_y(), false);
  ULONG zc=get_bin(z,0, params->_Nft(), params->_d_delta_z(), false);
  if(xc==params->_Nft())xc-=1;
  if(yc==params->_Nft())yc-=1;
  if(zc==params->_Nft())zc-=1;
  real_prec xx  = params->_d_delta_x()*(static_cast<real_prec>(xc)+0.5);
  real_prec yy  = params->_d_delta_y()*(static_cast<real_prec>(yc)+0.5);
  real_prec zz  = params->_d_delta_z()*(static_cast<real_prec>(zc)+0.5);
  real_prec xxf = params->_d_delta_x()*(static_cast<real_prec>(xc)+1.5);
  real_prec yyf = params->_d_delta_y()*(static_cast<real_prec>(yc)+1.5);
  real_prec zzf = params->_d_delta_z()*(static_cast<real_prec>(zc)+1.5);
  real_prec xxb = params->_d_delta_x()*(static_cast<real_prec>(xc)-0.5);
  real_prec yyb = params->_d_delta_y()*(static_cast<real_prec>(yc)-0.5);
  real_prec zzb = params->_d_delta_z()*(static_cast<real_prec>(zc)-0.5);
  ULONG Xb=(xc==0 ? params->_Nft(): xc);
  ULONG Yb=(yc==0 ? params->_Nft(): yc);
  ULONG Zb=(zc==0 ? params->_Nft(): zc);
  ULONG Xf=(xc==params->_Nft()-1 ? -1: xc);
  ULONG Yf=(yc==params->_Nft()-1 ? -1: yc);
  ULONG Zf=(zc==params->_Nft()-1 ? -1: zc);
  ULONG i_idx[] = {Xb-1, xc, Xf+1};
  ULONG j_idx[] = {Yb-1, yc, Yf+1};
  ULONG k_idx[] = {Zb-1, zc, Zf+1};
  vector<real_prec> MAS_xx =
    {
      MAS_TSC((xxb  - x)/static_cast<double>(params->_d_delta_x())),
      MAS_TSC((xx   - x)/static_cast<double>(params->_d_delta_x())),
      MAS_TSC((xxf  - x)/static_cast<double>(params->_d_delta_x())),
    };
  vector<real_prec> MAS_yy =
    {
      MAS_TSC((yyb  - y)/static_cast<double>(params->_d_delta_y())),
      MAS_TSC((yy   - y)/static_cast<double>(params->_d_delta_y())),
      MAS_TSC((yyf  - y)/static_cast<double>(params->_d_delta_y())),
    };
  vector<real_prec> MAS_zz =
    {
      MAS_TSC((zzb  - z)/static_cast<double>(params->_d_delta_z())),
      MAS_TSC((zz   - z)/static_cast<double>(params->_d_delta_z())),
      MAS_TSC((zzf  - z)/static_cast<double>(params->_d_delta_z())),
    };
  for(int ih=0;ih<MAS_xx.size();++ih)
    for(int jh=0;jh<MAS_yy.size();++jh)
      for(int kh=0;kh<MAS_zz.size();++kh)
	{
	  ULONG index_aux = index_3d(i_idx[ih], j_idx[jh], k_idx[kh],params->_Nft(),params->_Nft());
	  field[index_aux]+= weight*MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
#ifdef _MASS_WEIGHT_POWER_
	  field_mw[index_aux] +=weight_aux*weight_mass;
#endif
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void grid_assignment_PCS(Params *params,real_prec x, real_prec y, real_prec z, real_prec weight, vector<real_prec>&field)
{
  x -= params->_Xoffset() - 0.5*params->_Lbox();
  y -= params->_Yoffset() - 0.5*params->_Lbox();
  z -= params->_Zoffset() - 0.5*params->_Lbox();
  while(x< 0){x+=params->_Lbox();}
  while(y< 0){y+=params->_Lbox();}
  while(z< 0){z+=params->_Lbox();}
  while(x>=params->_Lbox()){x-=params->_Lbox();}
  while(y>=params->_Lbox()){y-=params->_Lbox();}
  while(z>=params->_Lbox()){z-=params->_Lbox();}
  long xc=static_cast<long>(get_bin(x,0, params->_Nft(), params->_d_delta_x(), true));
  long yc=static_cast<long>(get_bin(y,0, params->_Nft(), params->_d_delta_y(), true));
  long zc=static_cast<long>(get_bin(z,0, params->_Nft(), params->_d_delta_z(), true));
  if(xc==params->_Nft())xc=params->_Nft()-1;
  if(yc==params->_Nft())yc=params->_Nft()-1;
  if(zc==params->_Nft())zc=params->_Nft()-1;
  real_prec xx  = params->_d_delta_x()*(static_cast<real_prec>(xc)+0.5);//central cell
  real_prec yy  = params->_d_delta_y()*(static_cast<real_prec>(yc)+0.5);
  real_prec zz  = params->_d_delta_z()*(static_cast<real_prec>(zc)+0.5);
  real_prec xxf = params->_d_delta_x()*(static_cast<real_prec>(xc)+1.5);//one cell forward
  real_prec yyf = params->_d_delta_y()*(static_cast<real_prec>(yc)+1.5);
  real_prec zzf = params->_d_delta_z()*(static_cast<real_prec>(zc)+1.5);
  real_prec xxff = params->_d_delta_x()*(static_cast<real_prec>(xc)+2.5);//two cells forward
  real_prec yyff = params->_d_delta_y()*(static_cast<real_prec>(yc)+2.5);
  real_prec zzff = params->_d_delta_z()*(static_cast<real_prec>(zc)+2.5);
  real_prec xxb = params->_d_delta_x()*(static_cast<real_prec>(xc)-0.5);//once cell backward
  real_prec yyb = params->_d_delta_y()*(static_cast<real_prec>(yc)-0.5);
  real_prec zzb = params->_d_delta_z()*(static_cast<real_prec>(zc)-0.5);
  real_prec xxbb = params->_d_delta_x()*(static_cast<real_prec>(xc)-1.5);//two
  real_prec yybb = params->_d_delta_y()*(static_cast<real_prec>(yc)-1.5);
  real_prec zzbb = params->_d_delta_z()*(static_cast<real_prec>(zc)-1.5);
  long Xb=(xc==0 ? params->_Nft(): xc);// if at 0, when movng backward 1 cell, it ends at Nft-1
  long Yb=(yc==0 ? params->_Nft(): yc);
  long Zb=(zc==0 ? params->_Nft(): zc);
  long Xf=(xc==params->_Nft()-1 ? -1: xc); // if at Nft-1, when moving forward 1 cell, it ends at 0
  long Yf=(yc==params->_Nft()-1 ? -1: yc);
  long Zf=(zc==params->_Nft()-1 ? -1: zc);
  long Xbb=(xc==0 ? params->_Nft(): xc); // if at 0, when movnig backwards 2 cells, it goes to Nft-2
  long Ybb=(yc==0 ? params->_Nft(): yc);
  long Zbb=(zc==0 ? params->_Nft(): zc);
  if(xc!=0)
    Xbb=(xc==1 ? params->_Nft()+1: xc);  // if at 1, when subtracting 2, it goes to Nft-1
  if(yc!=0)
    Ybb=(yc==1 ? params->_Nft()+1: yc);
  if(zc!=0)
    Zbb=(zc==1 ? params->_Nft()+1: zc);
  long Xff=(xc==params->_Nft()-1 ? -1: xc);  //if at Nft-1, when adding 2, it goes to 1
  long Yff=(yc==params->_Nft()-1 ? -1: yc);
  long Zff=(zc==params->_Nft()-1 ? -1: zc);
  if(xc!=params->_Nft()-1)
    Xff=(xc==params->_Nft()-2 ? -2: xc); // if at Nft-2, if adding 2, it goes to zero
  if(yc!=params->_Nft()-1)
    Yff=(yc==params->_Nft()-2 ? -2: yc);
  if(zc!=params->_Nft()-1)
    Zff=(zc==params->_Nft()-2 ? -2: zc);
  long i_idx[MAX_MAS_DEG] = {Xbb-2, Xb-1, xc, Xf+1, Xff+2};
  long j_idx[MAX_MAS_DEG] = {Ybb-2, Yb-1, yc, Yf+1, Xff+2};
  long k_idx[MAX_MAS_DEG] = {Zbb-2, Zb-1, zc, Zf+1, Zff+2};
  vector<real_prec> MAS_xx=
    {
      MAS_PCS((xxbb - x)/static_cast<double>(params->_d_delta_x())),
      MAS_PCS((xxb  - x)/static_cast<double>(params->_d_delta_x())),
      MAS_PCS((xx   - x)/static_cast<double>(params->_d_delta_x())),
      MAS_PCS((xxf  - x)/static_cast<double>(params->_d_delta_x())),
      MAS_PCS((xxff - x)/static_cast<double>(params->_d_delta_x()))
    };
  vector<real_prec> MAS_yy=
    {
      MAS_PCS((yybb - y)/static_cast<double>(params->_d_delta_y())),
      MAS_PCS((yyb  - y)/static_cast<double>(params->_d_delta_y())),
      MAS_PCS((yy   - y)/static_cast<double>(params->_d_delta_y())),
      MAS_PCS((yyf  - y)/static_cast<double>(params->_d_delta_y())),
      MAS_PCS((yyff - y)/static_cast<double>(params->_d_delta_y()))
    };
  vector<real_prec> MAS_zz=
    {
      MAS_PCS((zzbb - z)/static_cast<double>(params->_d_delta_z())),
      MAS_PCS((zzb  - z)/static_cast<double>(params->_d_delta_z())),
      MAS_PCS((zz   - z)/static_cast<double>(params->_d_delta_z())),
      MAS_PCS((zzf  - z)/static_cast<double>(params->_d_delta_z())),
      MAS_PCS((zzff - z)/static_cast<double>(params->_d_delta_z()))
    };
  for(int ih=0;ih<MAS_xx.size();++ih)
    for(int jh=0;jh<MAS_yy.size();++jh)
      for(int kh=0;kh<MAS_zz.size();++kh)
        {
          real_prec weight_aux=weight*MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
          ULONG index_aux = index_3d(i_idx[ih], j_idx[jh], k_idx[kh],params->_Nft(),params->_Nft());
          field[index_aux] +=weight_aux;
#ifdef _MASS_WEIGHT_POWER_
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      	  field_mw[index_aux] +=weight_aux*weight_mass;
#endif
      	}
}
