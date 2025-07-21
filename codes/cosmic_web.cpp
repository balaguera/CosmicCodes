////////////////////////////////////////////////////////////////////////////
/** @file main_cosmicatlass.cpp
 *
 *  @brief
 *  @authors Andrés Balaguera-Antolinez
 */
#include "../headers/def.h"
#include "CosmiCalcLIB.h"
////////////////////////////////////////////////////////////////////////////
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  //  system("clear");
  time_t start_all;
  time(&start_all);
  char temp;
  string par_file;

#ifdef USE_GALAXY_TOOLS
  throw std::invalid_argument("USE_GALAXY_TOOLS is defined, enable double prec everywhere. Plese desable and compile again.");
#endif
  int NTHREADS = omp_get_max_threads();
  omp_set_num_threads(NTHREADS);
  string logfile = "logfile.log";
  ScreenOutput So(start_all, logfile);
  if (argc == 1)
  {
    So.usage(argv[0]);
    exit(1);
  }
  while ((temp = getopt(argc, argv, "hadi:w:")) != -1)
  {
    if (temp == 'h')
      So.usage(argv[0]);

    else if (temp == 'a') // Show authors
      So.author();

    else if (temp == 'i') // displays input parameters
    {
      par_file = argv[2];
      Params params(par_file);
      params.show_params();
    }
    else if (temp == 'd') // displays input parameters
    {
      So.show_preproc();
    }
    else if ( temp  == 'w') // assign bias to tracers.
    {
      FileOutput File;
      So.welcome_message_c();
      par_file = argv[2];
      Params params(par_file);

      //Read dm field
      vector<real_prec> dm_field(params._NGRID(), 0);
      File.read_array(params._Input_Directory_X() + params._Name_Catalog_X(), dm_field);


      // Read catalog
      Catalog cat(params);
      cat.read_catalog(params._Input_dir_cat() + params._file_catalogue(), 0);
      // Get density field from tracers


      vector<real_prec> tr_field_counts(params._NGRID(), 0);
      cat.get_density_field_grid(_COUNTS_,tr_field_counts);
      //Assign bias to tracers and get the bias field
      PowerSpectrumF power(params);
      power.object_by_object_bias(cat.Halo, dm_field);

      vector<real_prec> bias_field(params._NGRID(), 0);
      cat.get_density_field_grid(_COUNTS_,tr_field_counts,bias_field);
      string fileb=params._Output_directory()+"ls_bias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileb,bias_field);
      bias_field.clear();bias_field.shrink_to_fit();
  

     string fileo;
      vector<real_prec> tr_field(params._NGRID(), 0);
      cat.get_density_field_grid(_iBIAS_,tr_field);
      cat.get_density_field_grid(_COUNTS_,tr_field);
      fileo=params._Output_directory()+"field_hbias_p"+to_string(params._unitsim_plabel());

      File.write_array(fileo,tr_field);


      Cwclass cwclass(params); // for halos

      // Get tidal anisotropy from tracers
      vector<real_prec> tr_alpha(params._NGRID(), 0);
      cwclass.get_tidal_anisotropy(tr_field, tr_alpha);


   


      fileo=params._Output_directory()+"alpha_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,tr_alpha);
      fileo=params._Output_directory()+"l1_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda1);
      fileo=params._Output_directory()+"l2_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda2);
      fileo=params._Output_directory()+"l3_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.lambda3);

      fileo=params._Output_directory()+"cwc_hbias_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass.CWClass);

      vector<ULONG> cwc(params._NGRID(),0);

      // Get tidal anisotropy from dm

      vector<real_prec> dm_alpha(params._NGRID(), 0);
      Cwclass cwclass_dm(params); //para dm

      cwclass_dm.get_tidal_anisotropy(dm_field, dm_alpha);
      fileo=params._Output_directory()+"alpha_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,dm_alpha);

      vector<real_prec> alpha_bias(params._NGRID(), 0);
#pragma omp parallel for
    for(ULONG i=0;i<dm_field.size();++i)
      alpha_bias[i]=tr_alpha[i]/dm_alpha[i];



      fileo=params._Output_directory()+"bias_alpha_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,alpha_bias);
      fileo=params._Output_directory()+"l1_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda1);
      fileo=params._Output_directory()+"l2_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda2);
      fileo=params._Output_directory()+"l3_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.lambda3);
      fileo=params._Output_directory()+"cwc_dm_p"+to_string(params._unitsim_plabel());
      File.write_array(fileo,cwclass_dm.CWClass);

      // Asigna el cwc al ultimo calculado:


      // Here we rank orde the eigenvalues from the tracer tidal foeld to those form the dm
//      calc_pdf("log", params._NGRID(),NXn, this->Xmax, this->Xmin, this->delta_X_REF_PDF, this->pdf_ref);
//      rankorder(ZERO, xbaux, NXn,  this->Xmax, this->Xmin, delta_dm_aux_mem, pdf_in, this->pdf_ref);





      //cwc=cwclass.CWClass; // halos
      cwc=cwclass_dm.CWClass; // dm

      

      // Calculo de Histogramas


       int nb=200;
       real_prec bmin=-2; real_prec bmax=2; real_prec delta=(bmax-bmin)/static_cast<real_prec>(nb);

       for (int icwt=0; icwt<=4;++icwt)
        {
          vector<real_prec>hist(nb,0); 
          vector<real_prec>hist_alphah(nb,0); 
          vector<real_prec>hist_alphadm(nb,0); 

          if(icwt==0)
            {
              for(ULONG i=0;i<dm_field.size();++i)
                {
                  int ibin=get_bin(log10(alpha_bias[i]),bmin,nb,delta,0);
                  if(log10(alpha_bias[i])<bmax && log10(alpha_bias[i])>=bmin)
                    hist[ibin]++;

                  int abin=get_bin(log10(tr_alpha[i]),bmin,nb,delta,0);
                  if(log10(tr_alpha[i])<bmax && log10(tr_alpha[i])>=bmin)
                      hist_alphah[abin]++;

                      
                  int dbin=get_bin(log10(dm_alpha[i]),bmin,nb,delta,0);
                  if(log10(dm_alpha[i])<bmax && log10(dm_alpha[i])>=bmin)
                      hist_alphadm[dbin]++;


                  }
             }
          else
          {
            for(ULONG i=0;i<dm_field.size();++i)
            {
              if(cwc[i]==icwt)
                {
                  int ibin=get_bin(log10(alpha_bias[i]),bmin,nb,delta,0);
                  if(log10(alpha_bias[i])<bmax && log10(alpha_bias[i])>=bmin)
                    hist[ibin]++;

                  int abin=get_bin(log10(tr_alpha[i]),bmin,nb,delta,0);
                  if(log10(tr_alpha[i])<bmax && log10(tr_alpha[i])>=bmin)
                    hist_alphah[abin]++;

                  int dbin=get_bin(log10(dm_alpha[i]),bmin,nb,delta,0);
                  if(log10(dm_alpha[i])<bmax && log10(dm_alpha[i])>=bmin)
                    hist_alphadm[dbin]++;
  
                    
                }
            }
        }
        fileo=params._Output_directory()+"bias_alpha_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+".txt";
        cout<<"Writting alpha bias to file "<<fileo<<endl;
        ofstream off; off.open(fileo.c_str());
        for(int i=0; i<nb; ++i)
            off<<bmin+(i+0.5)*delta<<"\t"<<hist[i]/static_cast<real_prec>(dm_field.size())<<endl;
        off.close();

        fileo=params._Output_directory()+"alphah_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+".txt";
        cout<<"Writting alpha bias to file "<<fileo<<endl;
        off.open(fileo.c_str());
        for(int i=0; i<nb; ++i)
            off<<bmin+(i+0.5)*delta<<"\t"<<hist_alphah[i]/static_cast<real_prec>(dm_field.size())<<endl;
        off.close();

        fileo=params._Output_directory()+"alphadm_hbias_p"+to_string(params._unitsim_plabel())+"_cwt"+to_string(icwt)+".txt";
        cout<<"Writting alpha bias to file "<<fileo<<endl;
        off.open(fileo.c_str());
        for(int i=0; i<nb; ++i)
            off<<bmin+(i+0.5)*delta<<"\t"<<hist_alphadm[i]/static_cast<real_prec>(dm_field.size())<<endl;
        off.close();




        }

     dm_alpha.clear(); dm_alpha.shrink_to_fit();

     dm_field.clear();
     dm_field.shrink_to_fit();
     tr_field.clear();
     tr_field.shrink_to_fit();


    }
    else if ('?' == temp)
    {
      cout << endl;
      cout << "Argument not recognized." << endl;
      cout << "Please run ./cosmicatlass.exe -h." << endl;
      exit(1);
    }
  }
  exit(0);
  return 0;
}
// #################################################################################
// ##################################################################################
// #################################################################################
