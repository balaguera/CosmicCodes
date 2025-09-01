/** @file Catalogue.cpp
 *
 *  @brief Methods of the class Catalogue 
 *
 *  This file contains the implementation of the methods of the class
 *  Catalogue, used to handle catalogues of astronomical sources
 */

#include "../Include/Catalogue.h"
#include "../Include/CosmologicalFunctions.h"
#include <ctime>
#include <iostream>
#include <iterator>
#include <omp.h>


void Catalogue::read_catalogue (Params &params, bool type, int icatalogue)
{
  string file;
  int ic1, ic2, ic3, iw1, iw2, iw3, iw4;

  // set cosmolgical parameters

  s_CosmologicalParameters cosmo_par;
  cosmo_par.om_matter = params._om_matter();
  cosmo_par.om_radiation = params._om_radiation();
  cosmo_par.om_baryons = params._om_baryons();
  cosmo_par.om_vac = params._om_vac();
  cosmo_par.om_k = params._om_k();
  cosmo_par.Hubble = params._Hubble();
  cosmo_par.hubble = params._hubble();
  cosmo_par.spectral_index =  params._spectral_index();
  cosmo_par.w_eos  = params._w_eos();
  cosmo_par.N_eff = params._N_eff();
  cosmo_par.sigma8 = params._sigma8();
  cosmo_par.Tcmb = params._Tcmb();

  // data catalogue
  if (type == 0) {
    stringstream ss; ss << icatalogue;
    file = (params._n_catalogues()==1) ? params._file_catalogue() : params._file_catalogue()+"."+ss.str();
    ic1 = params._i_coord1_g()-1;
    ic2 = params._i_coord2_g()-1;
    ic3 = params._i_coord3_g()-1;
    iw1 = params._i_weight1_g()-1;
    iw2 = params._i_weight2_g()-1;
    iw3 = params._i_weight3_g()-1;
    iw4 = params._i_weight4_g()-1;
  }

  // random catalogue
  else {
    file = params._file_random();
    ic1 = params._i_coord1_r()-1;
    ic2 = params._i_coord2_r()-1;
    ic3 = params._i_coord3_r()-1;
    iw1 = params._i_weight1_r()-1;
    iw2 = params._i_weight2_r()-1;
    iw3 = params._i_weight3_r()-1;
    iw4 = params._i_weight4_r()-1;
  }

  time_t start, stop;
  char *BUFF;
  unsigned line_count = 0;

  ifstream inputf (file.c_str());
  if (!inputf) { cerr <<"Error in opening the input file "<<file<<"!"<<endl; exit(1); }
  BUFF = new char [20000];
  while (!inputf.eof())
    {
      inputf.read(BUFF, 20000);
      int count = inputf.gcount();
      for(int jj = 0; jj < count; jj++)
        if ( BUFF[jj] == '\n' )
          line_count++;
    }
  inputf.clear(); inputf.close();

  std::cout << line_count << " lines found" << std::endl;
  ifstream fin (file.c_str());
  resize_catalogue(line_count);

  unsigned int NTHREADS = omp_get_max_threads();
  string line[1000*NTHREADS]; 
  bool sstop = false;
  int line_idx = 0;
  unsigned int ii;


  time (&start);
#pragma omp parallel num_threads(NTHREADS)
  {
    //int myID = omp_get_thread_num();
    int buffer_size = 8;
    double aa, *cc;
    double ww = 1.0;

    cc = new double[buffer_size];

    while (!sstop) {      
#pragma omp single
      {
        for(ii = 0; ii < 1000 * NTHREADS && line_idx+ii < line_count; ii++)
          getline(fin,line[ii]);
        if(line_idx+ii == line_count)
          sstop = true;
      }
#pragma omp barrier      

#pragma omp for
      for(unsigned int ll = 0; ll < ii; ll++)
        {
          int jj = 0;
          stringstream ss(line[ll]);
          while (ss>>aa) cc[jj++] = aa;
          
          // const char *str;
          // char *nextstr;
          // str = line[ll].c_str();
          // while (nextstr != NULL) { cc[jj++] = strtod(str, &nextstr); str = nextstr;}
          ww = 1.0;
          if (0<=iw1 && iw1< jj) ww *= cc[iw1];
          if (0<=iw2 && iw2< jj) ww *= cc[iw2];
          if (0<=iw3 && iw3< jj) ww *= cc[iw3];
          if (0<=iw4 && iw4< jj) ww *= cc[iw4];

          int sc = type ? params._sys_of_coord_r() : params._sys_of_coord_g();
          Galaxy galaxy (cc[ic1], cc[ic2], cc[ic3], ww, sc, &cosmo_par);
          add_galaxy(galaxy, line_idx+ll);                  
        }
#pragma omp single
      line_idx += ii;
    }
  }
  time (&stop);
  std::cout << difftime(stop, start) << " seconds" << std::endl;

  fin.clear(); fin.close();
}



// void Catalogue::read_catalogue (Parameters &params, bool type)
// {
//   string file;
//   int ic1, ic2, ic3, iw1, iw2, iw3, iw4;

//   // set cosmolgical parameters

//   s_CosmologicalParameters cosmo_par;
//   cosmo_par.om_matter = params._om_matter();
//   cosmo_par.om_radiation = params._om_radiation();
//   cosmo_par.om_baryons = params._om_baryons();
//   cosmo_par.om_vac = params._om_vac();
//   cosmo_par.om_k = params._om_k();
//   cosmo_par.Hubble = params._Hubble();
//   cosmo_par.hubble = params._hubble();
//   cosmo_par.spectral_index = params._spectral_index();
//   cosmo_par.w_eos = params._w_eos();
//   cosmo_par.N_eff = params._N_eff();
//   cosmo_par.sigma8 = params._sigma8();
//   cosmo_par.Tcmb = params._Tcmb();

//   // data catalogue
//   if (type==0) {
//     file = params._file_catalogue();
//     ic1 = params._i_coord1_g()-1;
//     ic2 = params._i_coord2_g()-1;
//     ic3 = params._i_coord3_g()-1;
//     iw1 = params._i_weight1_g()-1;
//     iw2 = params._i_weight2_g()-1;
//     iw3 = params._i_weight3_g()-1;
//     iw4 = params._i_weight4_g()-1;
//   }

//   // random catalogue
//   else {
//     file = params._file_random();
//     ic1 = params._i_coord1_r()-1;
//     ic2 = params._i_coord2_r()-1;
//     ic3 = params._i_coord3_r()-1;
//     iw1 = params._i_weight1_r()-1;
//     iw2 = params._i_weight2_r()-1;
//     iw3 = params._i_weight3_r()-1;
//     iw4 = params._i_weight4_r()-1;
//   }


//   char *BUFF;
//   unsigned line_count = 0;
  
//   ifstream inputf (file.c_str());
//   if (!inputf) { cerr <<"Error in opening the input file "<<file<<"!"<<endl; exit(1); }
//   BUFF = new char [20000];
//   while (!inputf.eof())
//     {
//       inputf.read(BUFF, 20000);
//       int count = inputf.gcount();
//       for(int jj = 0; jj < count; jj++)
//         if ( BUFF[jj] == '\n' )
//           line_count++;
//     }
//   inputf.clear(); inputf.close();

//   std::cout << line_count << " lines found" << std::endl;
//   reserve_catalogue(line_count);
  
//   double aa;
//   string line; 
//   int gcount = 0;
  
//   ifstream fin (file.c_str());
//   if (!fin) { cerr <<"Error in opening the input file "<<file<<"!"<<endl; exit(1); }

//   while (getline(fin,line)) {
//     stringstream ss(line);
//     vector<double> cc;
//     while (ss>>aa) cc.push_back(aa);

//     double ww = 1.;
//     if (0<=iw1 && iw1<(int)cc.size()) ww *= cc[iw1];
//     if (0<=iw2 && iw2<(int)cc.size()) ww *= cc[iw2];
//     if (0<=iw3 && iw3<(int)cc.size()) ww *= cc[iw3];
//     if (0<=iw4 && iw4<(int)cc.size()) ww *= cc[iw4];

//     int sc = type ? params._sys_of_coord_r() : params._sys_of_coord_g();
//     Galaxy galaxy(cc[ic1], cc[ic2], cc[ic3], ww, sc, &cosmo_par);

//     add_galaxy(galaxy, gcount);
//     gcount++;
  
//   }
  
//   fin.clear(); fin.close();
// }


// void Catalogue::read_catalogue (Parameters &params, bool type)
// {
//   string file;
//   int ic1, ic2, ic3, iw1, iw2, iw3, iw4;

//   // set cosmolgical parameters

//   s_CosmologicalParameters cosmo_par = {
//     cosmo_par.om_matter = params._om_matter(),
//     cosmo_par.om_radiation = params._om_radiation(),
//     cosmo_par.om_baryons = params._om_baryons(),
//     cosmo_par.om_vac = params._om_vac(),
//     cosmo_par.om_k = params._om_k(),
//     cosmo_par.Hubble = params._Hubble(),
//     cosmo_par.hubble = params._hubble(),
//     cosmo_par.spectral_index =  params._spectral_index(),
//     cosmo_par.w_eos  = params._w_eos(),
//     cosmo_par.N_eff = params._N_eff(),
//     cosmo_par.sigma8 = params._sigma8(),
//     cosmo_par.Tcmb = params._Tcmb(),
//   };

//   // data catalogue
//   if (type) {
//     file = params._file_catalogue();
//     ic1 = params._i_coord1_g()-1;
//     ic2 = params._i_coord2_g()-1;
//     ic3 = params._i_coord3_g()-1;
//     iw1 = params._i_weight1_g()-1;
//     iw2 = params._i_weight2_g()-1;
//     iw3 = params._i_weight3_g()-1;
//     iw4 = params._i_weight4_g()-1;
//   }

//   // random catalogue
//   else {
//     file = params._file_random();
//     ic1 = params._i_coord1_r()-1;
//     ic2 = params._i_coord2_r()-1;
//     ic3 = params._i_coord3_r()-1;
//     iw1 = params._i_weight1_r()-1;
//     iw2 = params._i_weight2_r()-1;
//     iw3 = params._i_weight3_r()-1;
//     iw4 = params._i_weight4_r()-1;
//   }

//   time_t start, stop;
//   char *BUFF;
//   unsigned line_count = 0;
  
//   ifstream inputf (file.c_str());
//   if (!inputf) { cerr <<"Error in opening the input file "<<file<<"!"<<endl; exit(1); }
//   BUFF = new char [20000];
//   while (!inputf.eof())
//     {
//       inputf.read(BUFF, 20000);
//       int count = inputf.gcount();
//       for(int jj = 0; jj < count; jj++)
//         if ( BUFF[jj] == '\n' )
//           line_count++;
//     }
//   inputf.clear(); inputf.close();

//   std::cout << line_count << " lines found" << std::endl;
//   ifstream fin (file.c_str());
//   reserve_catalogue(line_count);

//   int NTHREADS = 1; //omp_get_max_threads();
//   string line[NTHREADS]; 
//   bool sstop = false;
//   int line_idx = 0;
//   int ii;

//   time (&start);
// #pragma omp parallel num_threads(NTHREADS)
//   {
//     int myID = omp_get_thread_num();
//     int buffer_size = 8;
//     double aa, *cc;
//     double ww = 1.;

//     cc = new double[buffer_size];

// #pragma omp single
//     ii = 0;
//     while (ii < NTHREADS && !fin.eof())
//       {
//          getline(fin,line[ii]);
//          ii++;
//       }
// #pragma omp barrier      

//     if(myID < ii)
//       {
//         stringstream ss(line[myID]);
//         int jj = 0;
//         while (ss>>aa) cc[jj++] = aa;
          
//         if (0<=iw1 && iw1< buffer_size) ww *= cc[iw1];
//         if (0<=iw2 && iw2< buffer_size) ww *= cc[iw2];
//         if (0<=iw3 && iw3< buffer_size) ww *= cc[iw3];
//         if (0<=iw4 && iw4< buffer_size) ww *= cc[iw4];
          
//         int sc = type ? params._sys_of_coord_r() : params._sys_of_coord_g();
//         Galaxy galaxy (cc[ic1], cc[ic2], cc[ic3], ww, sc, &cosmo_par);
//         add_galaxy(galaxy, line_idx+myID);                  
//       }
// #pragma omp single
//       line_idx += ii;
//     }

//   time (&stop);
//   std::cout << difftime(stop, start) << " seconds" << std::endl;

//   fin.clear(); fin.close();
// }



