/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
// This file contains the methods of the class FFTW_FUNCTION used in the estimates
// of clustering for EUCLID
// Developer:  Andres Balaguera Antolinez
// e-mail:     abalant@gmail.com
// Afiliation: INAF, OAR. Uni. Roma3
// 2013-2015
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


# include "../Headers/DnDz.h"
# include "../Headers/ScreenOutput.h"
# include <stdexcept>
# include <alm.h>
# include <alm_fitsio.h>
# include <healpix_map.h>
# include <healpix_map_fitsio.h>
# include <healpix_data_io.h>
# include <healpix_base.h>
//# include <healpix_base2.h>
# include <healpix_data_io.h>
# include <alm_powspec_tools.h>
# include <alm_healpix_tools.h>
//# include <cxxutils.h>


void DnDz::dndz(void *p, vector<gsl_real> &new_zz_v, vector<gsl_real>&new_dndz_v, vector< vector<gsl_real> > &new_dndz_matrix)
{

  ////////////////////////////////////////////////////////////////
  // This function computes the the mean number density
  // from the random catalogue, in order to be interpolated
  // at the position of the galaxies and randoms in the
  // estimate of power spectrum.
  ////////////////////////////////////////////////////////////////

  struct s_dndz * s_dndz_data = (struct s_dndz *)p;
  vector<s_Halo>prop = s_dndz_data->properties;
  vector<gsl_real>zz=s_dndz_data->zz;
  vector<gsl_real>rc=s_dndz_data->rc;
  int nz = zz.size();
  int n_dndz=s_dndz_data->n_dndz;
  int new_n_dndz=s_dndz_data->new_n_dndz;
  real_prec redshift_min_sample=s_dndz_data->redshift_min_sample;
  real_prec redshift_max_sample=s_dndz_data->redshift_max_sample;
  real_prec area_survey=s_dndz_data->area_survey;   //area survey in degrees
  int sys_of_coord_r=s_dndz_data->sys_of_coord_r;
  int i_coord1_r=s_dndz_data->i_coord1_r;
  int i_coord2_r=s_dndz_data->i_coord2_r;
  int i_coord3_r=s_dndz_data->i_coord3_r;
  int n_columns = s_dndz_data->n_columns;
  real_prec area_pixel=s_dndz_data->area_pixel;
  ULONG npixels=s_dndz_data->npixels;
  ULONG nside=s_dndz_data->nside;
  string file_dndz=s_dndz_data->file_dndz;

  real_prec area=area_survey*pow(M_PI/180.,2); /*converting to strad*/
  vector<gsl_real> zz_v(n_dndz,0);   /*define redshift vector for histogram*/
  vector<gsl_real> dndz_v(n_dndz,0); /*define dndz vector for histogram of full dNdz*/

  // *********************************************************************************
  // WARNINGS
  cout<<RED<<"Computing dN/dz from the random catalog"<<RESET<<endl;
  cout<<RED<<"--->Warning 1: "<<CYAN<<" mean number density from redshift distribution computed"<<endl;
  cout<<"only if positions of random objects are in spherical coordinates"<<endl;
  cout<<"If positions of random objects are in cartesian coordinates,"<<endl;
  cout<<"mean number density set to unity. "<<RESET<<endl;
  cout<<RED<<"--->Warning 2: "<<CYAN<<" the area of the surveyed region must be well known"<<endl;
  cout<<"Otherwise the estimates of number density are biased. "<<endl;
  cout<<"Better think of using a MASK"<<RESET<<endl;
  // *********************************************************************************

  real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/static_cast<double>(n_dndz);
  ULONG n_objects = prop.size();

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<n_dndz;i++)
    zz_v[i]=redshift_min_sample+(i+0.5)*Delta_Z;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<n_dndz;i++)dndz_v[i]=0;


  if(sys_of_coord_r==1){
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<n_objects;i++)
    {
      real_prec zro= gsl_inter_new(rc, zz, prop[i].coord3);
      int count= get_bin(zro, redshift_min_sample,n_dndz,Delta_Z,true);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
            dndz_v[count]++;
     }
   }
  else  if(sys_of_coord_r==2)
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<n_objects;i++)
    {
      int count= get_bin(prop[i].coord3, redshift_min_sample,n_dndz,Delta_Z,true);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
       dndz_v[count]++;
    }


   // Compute Volume in redshift shell and divide by it to get nbar
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<n_dndz;i++)
    {
      real_prec rmax=gsl_inter_new(zz, rc, zz_v[i]+0.5*Delta_Z);
      real_prec rmin=gsl_inter_new(zz, rc, zz_v[i]-0.5*Delta_Z);
      real_prec Volumen=(area_survey*M_PI/180.0)*(pow(rmax,3)-pow(rmin,3))/3.;
      dndz_v[i]=dndz_v[i]*(s_dndz_data->alpha0)/Volumen;  //do not divide by Delta Z.
    }

  /*Smooth the histogram with a cubic spline*/
  gsl_bspline(zz_v, dndz_v, new_zz_v, new_dndz_v);
    // or Pass the same historgram without smooth at all,
//  new_zz_v=zz_v;
//  new_dndz_v=dndz_v;
  Fm.write_to_file(file_dndz,new_zz_v,new_dndz_v);
    // **************************************************************************
    // If depth varies across the sky, go to pixelization and evaluate
    // the dNdz in each pixel, dividing the result by the respective volume
    // **************************************************************************
  real_prec fac=M_PI/180.0;

#ifdef USE_HEALPIX_DNDZ_
  if(false==s_dndz_data->constant_depth)
   {
    cout<<CYAN<<"Pixelizing to get nbar in a limiting flux-varying sample"<<endl;
    cout<<"using "<< npixels<< " pixels"<<RESET<<endl;

    vector< vector<real_prec> > dndz_v(n_dndz, vector<real_prec> (npixels, 0));

    // Create vector for the map
    pointing point;
    Healpix_Map<real_prec>map(log2(nside), RING);

    vector<int> pixel_counter(npixels,0);
    real_prec Delta_Z=(redshift_max_sample-redshift_min_sample)/n_dndz;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<n_dndz;i++)
      zz_v[i]=redshift_min_sample+(i+0.5)*Delta_Z;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<n_dndz;i++)
      for(int j=0;j<npixels;j++)
        dndz_v[i][j]=0;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
    for(int i=0;i<n_objects;i++)
      {
	      point.phi=prop[i*n_columns + i_coord1_r]*fac;   //Assuming that coordinates are coming in degrees
	      point.theta=0.5*M_PI-prop[i*n_columns + i_coord2_r]*fac;
	      long npix=0;
	      npix=map.ang2pix(point);
        pixel_counter[npix]++;
        real_prec zro=(sys_of_coord_r==1 ? gsl_inter_new(rc, zz, prop[i*n_columns + i_coord3_r]): prop[i*n_columns + i_coord3_r]);
	      if(zro<=redshift_max_sample && zro>=redshift_min_sample)
          {
	          int count =(int)floor(fabs(zro-redshift_min_sample)/Delta_Z);
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
            dndz_v[count][npix]++;
	       }
      }
      //////////////////////////////////////////////////////////////////
      // Here we are simply counting those empty pixels in order to give
      // an estimate of the syrveyed area,
      // but this is not fully correct. A mask should be used instead.
      // However, we only need the area of the pixel
      double nc=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:nc)
#endif
      for(int i=0;i<npixels;i++)
          if(pixel_counter[i]!=0)
              nc++;
      area=area_pixel*nc;
      cout<<"Area pixel [squared deg]    = "<<area_pixel*pow(M_PI/180.,-2)<<endl;
      cout<<"Computed area [squared deg] = "<<area*pow(acos(-1.0)/180.,-2)<<endl;
      //////////////////////////////////////////////////////////////////


      //ONCE THE dN/dz has been determined from the random catalog,
      //use the Volumen with the fiducial cosmology
      //to generate an estimate of nbar , which will be interpolated
      //below at the position of the objects

  for(int j=0;j<npixels;j++)
    {
      vector<gsl_real> dndz_vv(n_dndz,0); /*define dndz vector for histogram*/
	    for(int i=0;i<n_dndz;i++){
        real_prec rmax=      gsl_inter_new(zz, rc,new_zz_v[i]+0.5*Delta_Z);
        real_prec rmin=      gsl_inter_new(zz, rc,new_zz_v[i]-0.5*Delta_Z);
        real_prec Volumen=   area_pixel*(pow(rmax,3)-pow(rmin,3))/3.;
	      dndz_vv[i]=(s_dndz_data->alpha0)*dndz_v[i][j]/Volumen;
	     }
	// Spline in each pixel
	     gsl_bspline(zz_v,dndz_vv,new_zz_v,new_dndz_v);
	     for(int i=0;i<new_dndz_matrix.size();i++)
          new_dndz_matrix[i][j]=(sys_of_coord_r== 0 ? 1.0 : new_dndz_v[i]);
      }
    }
#endif


}
