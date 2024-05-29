#ifdef _USE_LPT_
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
/** @file LPT.cpp
   @brief Generation of  DM field using ALPT
   @author Andrés Balaguera, based on LPT code by F-S Kitaura
 */
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
# include "../Headers/LPT.h"
# include "../Headers/PowerSpectrumTH.h"
# include "../Headers/PowerSpectrumF.h"

inline void check_malloc(void *p, int line) {
  if (p == NULL) {
    fprintf(stderr, "malloc returned null at line %d", line);
    exit(EXIT_FAILURE);
  }

}
//##################################################################################
//##################################################################################




void LPT::set_params(Params _my_params)
{
#ifdef _FULL_VERBOSE_
  So.message_screen("Loading parameters for LPT");
#endif
  planepar=true;
  this->params=_my_params;  

  // For LPT we conserve the distinction among lengs and number of cells in the three directions.
  this->N1=this->params._N1();
  this->N2=N1;
  this->N3=N1;
  this->L1=this->params._Lbox();
  this->L2=L1;
  this->L3=L1;
  this->Initial_Redshift_DELTA=params._Initial_Redshift_DELTA();
  this->fnameIC=this->params._dir()+this->params._ic_WN_dir()+this->params._ic_WN_file();///string("deltaIC")+straddrefgen;
  this->Normalize_IC_to_initial_redshift=params._Normalize_IC_to_initial_redshift();
  // Assign values to the cosmological parameter stucture
  this->s_cosmo_pars.cosmological_redshift=params._redshift();
  this->s_cosmo_pars.Hubble=params._Hubble();
  this->s_cosmo_pars.hubble=params._hubble();
  this->s_cosmo_pars.Om_matter = params._om_matter();
  this->s_cosmo_pars.Om_cdm = params._om_cdm();
  this->s_cosmo_pars.Om_baryons = params._om_baryons();
  this->s_cosmo_pars.Om_radiation = params._om_radiation();
  this->s_cosmo_pars.Om_vac = params._om_vac();
  this->s_cosmo_pars.Om_k = params._om_k();
  this->s_cosmo_pars.n_s = params._spectral_index();
  this->s_cosmo_pars.w_eos = params._w_eos();
  this->s_cosmo_pars.sigma8 = params._sigma8();
  this->s_cosmo_pars.N_eff = params._N_eff();
  this->s_cosmo_pars.f_baryon = params._om_baryons()/params._om_matter();
  this->s_cosmo_pars.use_wiggles = params._use_wiggles();
  this->s_cosmo_pars.RR = params._RR();
  PowerSpectrum Pow(params.s_cosmo_pars);
  this->s_cosmo_pars.growth_factor=this->s_cosmo_info.growth_factor;
  this->s_cosmo_pars.pk_normalization=Pow.normalization();
  this->NGRID= this->params.d_NGRID();
  this->growth_ini=1.0;
  if(this->params._Initial_Redshift_TH_power_file() >0)
     this->growth_ini = this->Cosmo.growth_factor(this->params._Initial_Redshift_TH_power_file())/this->Cosmo.growth_factor(0.0);

  So.DONE();

}

//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

void LPT::set_fnames()
{
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Definig file names for LPT");
#endif
  real_prec kth=  this->params._slength();


  real_prec redshift=this->s_cosmo_pars.cosmological_redshift;


  int bmax=100;

  string buffz;
  string buffz_bam;

  char buffzc[bmax];

  sprintf(buffzc,"z%.3f",redshift);
  buffz=static_cast<string>(buffzc);
  buffz_bam=this->params._Name_survey();

  char buffvolc[bmax];

  sprintf(buffvolc,"V%.1f",L1);
  string buffvol=static_cast<string>(buffvolc);

  char buffgc[bmax];
  sprintf(buffgc,"G%d",static_cast<int>(N1));
  string buffg=static_cast<string>(buffgc);

  char buffs[bmax];
  //sprintf(buffs,"",static_cast<int>(seed));
  //if (inputmode==0)
  sprintf(buffs,"S%d",static_cast<int>(this->params._seed()));

  char buffsref[bmax];
   if (this->params._seed_ref()>0)
    sprintf(buffsref,"S%d",static_cast<int>(this->params._seed_ref()));
  else
    sprintf(buffsref,"S%d",static_cast<int>(this->params._seed()));


  string buffsl;
  //char buffslc[bmax];
  string buffslc=to_string(static_cast<int>(kth));
  //  sprintf(buffslc,"",static_cast<int>(kth));
  if (this->params._Structure_Formation_Model()==3) // for ALPT
    {
      char buffsl0[bmax];
      int slnum;
      buffslc=to_string(static_cast<int>(kth));
      //sprintf(buffslc,"r%d",static_cast<int>(kth));
      slnum=static_cast<int>(kth*static_cast<real_prec>(10.));
      slnum-=static_cast<int>(kth)*static_cast<int>(10);
      sprintf(buffsl0,"_%d",slnum);
      buffsl=static_cast<string>(buffslc)+static_cast<string>(buffsl0);
    }

  string bufftf=string("");

  string buffmk;
  switch (this->params._masskernel())
    {
    case 0:
      {
        buffmk=string("NGP");
      }
      break;
    case 1:
      {
        buffmk=string("CIC");
      }
      break;
    case 2:
      {
        buffmk=string("TSC");
      }
      break;
    }

  switch (this->params._Structure_Formation_Model())
    {
    case 1:
      {
        this->buffsf=string("ZELD");
      }
      break;
    case 2:
      {
        this->buffsf=string("2LPT");
      }
      break;
    case 3:
      {
        this->buffsf=string("ALPT");
        char buffalptc[bmax];
        sprintf(buffalptc,"rS%.1f",kth);
        string buffalpt=static_cast<string>(buffalptc);
        buffsf+=buffalpt;
      }
      break;
    }

#ifdef _USE_TET_
  this->buffsf+=string("TET");
#endif

  this->stradd=buffsf+buffmk+buffz+bufftf+buffg+buffvol+buffs;
  this->stradd_bam=buffs;

  string straddref=buffsf+buffmk+buffz+bufftf+buffg+buffvol+buffsref;
  string straddrefgen=buffg+buffvol+buffsref;

  this->fname3DPOWER=this->params._Output_directory()+string("3DpowerIC")+buffvol;
  this->fnameTHETA=this->params._Output_directory()+string("theta")+straddref;

#ifdef _MOVE_DM_TO_REDSHIFT_SPACE_
  this->fnameDM=this->params._Output_directory()+string("densDM_rss")+straddref;
  this->fnameDMNGP=this->params._Output_directory()+string("densDMNGP_rss")+straddref;
#else
  this->fnameDM=this->params._Output_directory()+string("densDM")+straddref;
  this->fnameDMNGP=this->params._Output_directory()+string("densDMNGP")+straddref;
#endif

#ifdef SAVECOMP
  this->fnameIC=this->params._dir()+this->params._ic_WN_file();///string("deltaIC")+straddrefgen;

  if(false==this->params._use_ic_file())
    this->fnameICDELTA=this->params._dir()+this->params._ic_WN_file()+"ICField";///string("deltaIC")+straddrefgen;
  else
    this->fnameICDELTA=this->params._Output_directory()+"ICDELTAField";

  this->fnameS2TERM=this->params._Output_directory()+string("stwoterm")+straddref;
  this->fnameS2TERMEUL=this->params._Output_directory()+string("stwotermeul")+straddref;
  this->fnameP2TERMEUL=this->params._Output_directory()+string("ptwotermeul")+straddref;
  this->fnameS3TERM=this->params._Output_directory()+string("sthreeterm")+straddref;
  this->fnameS3TERMEUL=this->params._Output_directory()+string("sthreetermeul")+straddref;
  this->fnameSTTERM=this->params._Output_directory()+string("stterm")+straddref;
  this->fnameSTTERMEUL=this->params._Output_directory()+string("sttermeul")+straddref;
  this->fnamePSITERMEUL=this->params._Output_directory()+string("psitermeul")+straddref;
  this->fname2LPTTERM=this->params._Output_directory()+string("twolptterm")+straddrefgen;
  this->fname2LPTTERMEUL=this->params._Output_directory()+string("twolpttermeul")+straddref;

  this->fnamePOSX=this->params._Output_directory()+string("posx")+straddref;
  this->fnamePOSY=this->params._Output_directory()+string("posy")+straddref;
  this->fnamePOSZ=this->params._Output_directory()+string("posz")+straddref;

  this->fnamePOSX_RS=this->params._Output_directory()+string("posx_rss")+straddref;
  this->fnamePOSY_RS=this->params._Output_directory()+string("posy_rss")+straddref;
  this->fnamePOSZ_RS=this->params._Output_directory()+string("posz_rss")+straddref;


  this->fnameVXpart=this->params._Output_directory()+string("vxpart")+straddref;
  this->fnameVYpart=this->params._Output_directory()+string("vypart")+straddref;
  this->fnameVZpart=this->params._Output_directory()+string("vzpart")+straddref;

  this->fnameVX=this->params._Output_directory()+string("vx")+straddref;
  this->fnameVY=this->params._Output_directory()+string("vy")+straddref;
  this->fnameVZ=this->params._Output_directory()+string("vz")+straddref;

#else
  string fnameIC=string("auxdeltaIC");

  this->fnamePOSX=string("auxposx");
  this->fnamePOSY=string("auxposy");
  this->fnamePOSZ=string("auxposz");

  this->fnameVX=string("auxvx");
  this->fnameVY=string("auxvy");
  this->fnameVZ=string("auxvz");
#endif
  this->So.DONE();

}


// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================

#ifdef _DISPLACEMENTS_
#ifdef _USE_OMP_
void LPT::get_displacement(gsl_rng ** gBaseRand, vector<real_prec>&DELTA_IC, int it)
#else
void LPT::get_displacement(gsl_rng * gBaseRand, vector<real_prec>&DELTA_IC, int it)
#endif
{
    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  real_prec kth=  this->params._slength();
  int ftype=1;// Smoothing kernel. 1 Is forGaussian
  kernelcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,kth,ftype,this->params._Output_directory());

  if(it==0)
   if(true==this->Normalize_IC_to_initial_redshift)
    normalize_df_z_ini(DELTA_IC,"delta",0,0);

  if(it==0)
    File.write_array(this->fnameICDELTA, DELTA_IC); // Lag2Eul reads thois file

  this->Lag2Eul_comp(1.0,kth,ftype,true,DISP_COORD);
}

#endif
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
// =================================================================================================
#ifndef _DISPLACEMENTS_
#ifdef _USE_OMP_
void LPT::get_dm_field(gsl_rng ** gBaseRand)
#else
void LPT::get_dm_field(gsl_rng * gBaseRand)
#endif
{

    So.enter(__PRETTY_FUNCTION__);
#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  cout<<endl;
#ifdef _FULL_VERBOSE_
  this->So.message_screen("OBTAINING DARK MATTER DENSITY FIELD with LPT");
#endif

  real_prec kth=  this->params._slength();
  real_prec kthv= this->params._slengthv();

  int ftype=1;// Smoothing kernel. 1 Is forGaussian
  ULONG N=this->NGRID;

  real_prec min1=this->params._xmin();
  real_prec min2=this->params._ymin();
  real_prec min3=this->params._zmin();

  int facL=1;//attention!!

  ULONG NL1p=this->N1*static_cast<ULONG>(facL);
  ULONG NL2=this->N2*static_cast<ULONG>(facL);
  ULONG NL3=this->N3*static_cast<ULONG>(facL);


  ULONG NLL=static_cast<ULONG>(NL1p)*static_cast<ULONG>(NL2)*static_cast<ULONG>(NL3);

  vector<real_prec> dummy(N,0);



  if (true==this->params._runsim())
    if(false==this->params._use_ic_file())
      if (true==this->params._readPS())
        {
          cout<<endl;
#ifdef _VERIFY_FILES_
          if(false==this->File.exists(this->fname3DPOWER))
            {
#endif
              cout<<endl;
              this->So.message_screen("Reading Initial P(k) from table and interpolating on Fourier grid");

              this->read_tabulated_power();
#ifdef _VERIFY_FILES_
            }
#endif

        }
      else
        {
          // if (inStream.is_open() == false )
          {
            cout<<endl;
            this->So.message_screen("Computing P(k) from fitting formulae");
            // PowerSpectrum Power;
            // real_prec Pa=Power.Linear_Matter_Power_Spectrum(&this->s_cosmo_pars, 0.1);

            // ATTENTION, to be placed in the right class, commented by ABA
            //initialize_pow_spec(s_LPT_pars,fnamePOWER);
          }
        }


  string fnamePOSXORD=this->params._Output_directory()+string("aux1");//+stradd;
  string fnamePOSYORD=this->params._Output_directory()+string("aux2");//+stradd;
  string fnamePOSZORD=this->params._Output_directory()+string("aux3");//+stradd;

  if (true==this->params._runv())
    kernelcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,kthv,ftype,this->params._Output_directory());

  // (1) start running simulation
  if (true==this->params._runsim())
    {


      kernelcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,kth,ftype,this->params._Output_directory());

      if(false==this->params._use_ic_file())
        {

          // if (vsmoo>0.)
          // 	kernelcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,vsmoo,ftype);

          switch (this->params._inputmode())
            {
            case 0:   // 0 for White Noise generated here
              {
                //	    vector<real_prec>dummy2(N,0);

#ifndef SAVECOMP
                this->So.message_screen("Generate initial gaussian field");

                // Read input power spectrum

                vector<real_prec>dummy2(N,0);

                this->File.read_array(this->fnamePOWER3D+".dat",dummy2);

#ifdef OMPPARRAN
#ifdef OMPPARRANGAR
                create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand);
#else
                create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand[0]);
#endif
#else
                create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand);
#endif

                this->So.Done();
#else
                bifstream inStream2(this->fnameIC.data(),file_is_natural);


                if (inStream2.is_open()==true)
                  {
                    this->File.read_array(this->fnameIC,dummy);  //read white noise if the fnameIC file does not exists (perhaps created in the case above with def SAVECOMP)
                  }
                  else
                  {
                    this->So.message_screen("Generate initial gaussian field");
                    vector<real_prec> dummy2(N);
                    this->File.read_array(this->fname3DPOWER+".dat",dummy2);

#ifdef OMPPARRAN
#ifdef OMPPARRANGAR
                    create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand);
#else
                    create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand[0]);
                    //create_GARFIELD_FIXED_AMP(N1,N2,N3,dummy,dummy2,gBaseRand[0]);

                    real_prec meanWN, sigma2D;
                    meanWN=get_mean(dummy);
#ifdef _FULL_VERBOSE_
                    this->So.message_screen("Mean WN =",meanWN);
#endif

                    sigma2D=get_var(meanWN, dummy);
#ifdef _FULL_VERBOSE_
                    this->So.message_screen("Sigma² WN =",sigma2D);
#endif


#endif
#else
                    create_GARFIELDR(N1,N2,N3,dummy,dummy2,gBaseRand);
#endif
                    this->File.write_array(this->fnameICDELTA,dummy);
                    this->So.DONE();

#ifdef _FULL_VERBOSE_
                    So.message_screen("Measuring power spectrum of IC");
#endif
                    this->params.set_SN_correction(false);
                    this->params.set_mass_assignment_scheme("CIC");
                    if(true==this->params._ic_alias_corrected())
                          this->params.set_MAS_correction(false);
                    else
                        this->params.set_MAS_correction(true);
                    this->params.set_input_type("delta");

                    this->params.set_Name_survey("IC_GRF");
                    PowerSpectrumF cPSF(this->params);
                    cPSF.compute_power_spectrum_grid(dummy,true);
                    this->So.DONE();
#ifdef _USE_GNUPLOT_

                     int Nk=600;
                     vector<real_prec>xbaux(Nk, 0);
                     vector<real_prec>pdf(Nk, 0);
                     real_prec max=get_max<real_prec>(dummy);
                     real_prec min=get_min<real_prec>(dummy);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                    for(int i=0;i<Nk; ++i)
                      xbaux[i]=min+static_cast<real_prec>(i+0.5)*(max-min)/static_cast<real_prec>(Nk);
                    string filex=this->params._Output_directory()+"pdf_IC"+to_string(this->params._seed())+".txt";
                    calc_pdf("lin", dummy.size(),Nk,max,min,dummy,pdf);
                    this->File.write_to_file(filex, xbaux,pdf);
                    vector<pair<real_prec, real_prec> > xy_pts_ref;
                    for(int i=1; i<Nk;++i) 
                      xy_pts_ref.push_back(std::make_pair(xbaux[i], pdf[i])); 
                    xbaux.clear(); xbaux.shrink_to_fit();
                    pdf.clear(); pdf.shrink_to_fit();
                    this->gp_pdf<<"set border linewidth 1.5\n";
                    this->gp_pdf<<"set border linecolor '"<<FG_COLOR<<"' \n";
                    this->gp_pdf<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";

                    this->gp_pdf<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
                    this->gp_pdf<<"set ytics textcolor '"<<FG_COLOR<<"' \n";

                    this->gp_pdf << "set log y\n";
                    this->gp_pdf << "set title 'R"<< this->params._seed()<<"'\n";
                    this->gp_pdf << "set xlabel '{/Symbol d}' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
                    this->gp_pdf << "set ylabel 'P({/Symbol d})' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
                    this->gp_pdf << "set size square\n";
                    this->gp_pdf << "plot" << gp_pdf.file1d(xy_pts_ref) << "w l lw 2 lt 1 title 'PDF IC'"<<endl;

                    xy_pts_ref.clear();
                    xy_pts_ref.shrink_to_fit();

                    for(int i=1; i<cPSF._kvector_data_size(); ++i) 
                      xy_pts_ref.push_back(std::make_pair(cPSF._kvector_data(i), cPSF._pk0(i))); 
                    this->gp_power<<"set size 1.0,0.5\n";
                    this->gp_power<<"set origin 0.,0.5\n";
                    this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
                    this->gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
                    this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n"; 
                    this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
                    this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";

                    this->gp_power << "set log\n";
                    this->gp_power << "set xlabel 'k [h / Mpc]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"'  \n";
                    this->gp_power<< "set ylabel 'P(k) [(Mpc / h)³]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
                    this->gp_power << "set size square\n";
                    this->gp_power << "plot" << gp_power.file1d(xy_pts_ref) << "w l lw 2 lt 1 title 'IC'"<<endl;
                    xy_pts_ref.clear();
                    xy_pts_ref.shrink_to_fit();
#endif


                    }
#endif
              }
              break;
            case 1:   // White Noise taken from some input file
              {

#define WHITE_NOISE dummy
#ifdef _FULL_VERBOSE_
                So.message_screen("White Noise:");
#endif
                this->File.read_array(this->fnameIC,WHITE_NOISE);               
                real_prec meanWN=get_mean(WHITE_NOISE);
                real_prec sigma2D=get_var(meanWN,WHITE_NOISE);
                  this->So.message_screen("Mean  WN =",meanWN);
                this->So.message_screen("Sigma² WN =",sigma2D);

                So.message_screen("Generating Initial DM field with input Power Spectrum");

#define DELTA_IC dummy  // dummy change names from WHITE_NOISE TO DELTA_IC

                create_GARFIELDR_from_WHITENOISE(this->fname3DPOWER+".dat",N1,N2,N3,DELTA_IC);
                this->So.DONE();


#ifdef _VERIFY_FILES_
                if(false==this->File.exists(this->fnameICDELTA))
#endif
                this->File.write_array(this->fnameICDELTA,DELTA_IC);

                sigma2D=get_var(DELTA_IC);
#ifdef _FULL_VERBOSE_
                this->So.message_screen("Sigma² IC with P(k) =",sigma2D);
                this->So.message_screen("Measuring power spectrum of IC: This must be in agreement with the input power spectrum");
#endif
                this->params.set_SN_correction(false);
                this->params.set_mass_assignment_scheme("CIC");
                if(true==this->params._ic_alias_corrected())
                      this->params.set_MAS_correction(false);
                else
                    this->params.set_MAS_correction(true);
                this->params.set_input_type("delta_grid");
                this->params.set_Name_survey("IC");
                PowerSpectrumF cPSF(this->params);
                cPSF.compute_power_spectrum_grid(DELTA_IC,true);
                this->So.DONE();

#ifdef _USE_GNUPLOT_
                int Nk=600;
                vector<real_prec>xbaux(Nk, 0);
                vector<real_prec>pdf(Nk, 0);
                real_prec max=get_max<real_prec>(dummy);
                real_prec min=get_min<real_prec>(dummy);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
                for(int i=0;i<Nk; ++i)
                  xbaux[i]=min+static_cast<real_prec>(i+0.5)*(max-min)/static_cast<real_prec>(Nk);
                string filex=this->params._Output_directory()+"pdf_IC"+to_string(this->params._seed())+".txt";
                calc_pdf("lin", dummy.size(),Nk,max,min,dummy,pdf);
                this->File.write_to_file(filex, xbaux,pdf);
                vector<pair<real_prec, real_prec> > xy_pts_ref;
                for(int i=0; i<Nk;++i) 
                  xy_pts_ref.push_back(std::make_pair(xbaux[i], pdf[i])); 
                xbaux.clear(); xbaux.shrink_to_fit();
                pdf.clear(); pdf.shrink_to_fit();
                    pdf.clear(); pdf.shrink_to_fit();
                this->gp_pdf<<"set border linewidth 1.5\n";
                this->gp_pdf<<"set border linecolor '"<<FG_COLOR<<"' \n";
                this->gp_pdf<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
                this->gp_pdf<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
                this->gp_pdf<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
                this->gp_pdf<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
                this->gp_pdf << "set title 'LOS"<< this->params._seed()<<"'\n";
                this->gp_pdf << "set log y\n";
                this->gp_pdf << "set xlabel 'delta '\n";
                this->gp_pdf << "set ylabel 'P(delta)'\n";
                this->gp_pdf << "set size square\n";
                this->gp_pdf << "plot" << gp_pdf.file1d(xy_pts_ref) << "w l lw 3 lt 1 title 'PDF IC'"<<endl;
                xy_pts_ref.clear();
                xy_pts_ref.shrink_to_fit();

                for(int i=1; i<cPSF._kvector_data_size(); ++i) 
                   xy_pts_ref.push_back(std::make_pair(cPSF._kvector_data(i), cPSF._pk0(i))); 
                this->gp_power << "set log\n";
                this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
                this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
                this->gp_power<<"set title textcolor rgb '"<<FG_COLOR<<"' \n";
                this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
                this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
                this->gp_power << "set xlabel 'k [h / Mpc]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
                this->gp_power << "set ylabel 'P(k) [(Mpc / h)³]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
                this->gp_power << "set size square\n";
                this->gp_power<<"set grid \n";
                this->gp_power << "plot" << gp_power.file1d(xy_pts_ref) << "w l lw 2 lt 1 title 'IC @ z=0'"<<endl;
                xy_pts_ref.clear();
                xy_pts_ref.shrink_to_fit();
#endif

              }

              break;
            }
        }
      else if(true==this->params._use_ic_file())
        {
#undef DELTA_IC
#define DELTA_IC dummy

#ifdef _FULL_VERBOSE_
          So.message_screen("Reading initial DM field with P(k)");
#endif
          this->File.read_array(this->params._ic_file(),DELTA_IC);
          //--------------------> This is meant to measure the power of the original input field
          this->params.set_SN_correction(false);
          if(this->params._ic_input_type()==DENSITY)
            this->params.set_input_type("density_grid");
          else
            this->params.set_input_type("delta_grid");
          this->params.set_mass_assignment_scheme("CIC");
          if(true==this->params._ic_alias_corrected())
                this->params.set_MAS_correction(false);
          else
              this->params.set_MAS_correction(true);
          this->params.set_Name_survey("IC_original");
          // This is meant to measure the power of the original input field
          PowerSpectrumF aPSF(this->params);
          aPSF.compute_power_spectrum_grid(DELTA_IC,true);
          //-------------------->


#ifdef _USE_GNUPLOT_
          vector<pair<real_prec, real_prec> > xy_pts_a;
          for(int i=1; i<aPSF._kvector_data_size(); ++i) 
           xy_pts_a.push_back(std::make_pair(aPSF._kvector_data(i), aPSF._pk0(i))); 
          this->gp_power<<"set border linewidth 1.5\n";
          this->gp_power<<"set border linecolor '"<<FG_COLOR<<"' \n";
          this->gp_power<<"set xtics textcolor '"<<FG_COLOR<<"' \n"; 
          this->gp_power<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
          this->gp_power<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
          this->gp_power << "set log\n";
          this->gp_power << "set xlabel 'k [h / Mpc]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
          this->gp_power << "set ylabel 'P(k) [(Mpc / h)³]' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"' \n";
          this->gp_power << "set size square\n";
          this->gp_power << "plot" << gp_power.file1d(xy_pts_a) << "w l lw 2 lt 8 title 'IC @ z = "<<this->params._Initial_Redshift_DELTA()<<"'"<<endl;

          xy_pts_a.clear();
          xy_pts_a.shrink_to_fit();
#endif


         //-------------------->We read a density field with the P(k) initial, we need to transform it to delta in order to pass it to the Grav Solver
         if(DENSITY==this->params._ic_input_type())
            get_overdens(DELTA_IC, DELTA_IC);
          //-------------------->

          //--------------------> This is meant to measure some statistis of the IC
          real_prec meanWN=0; 
          real_prec sigma2D=0;
          meanWN=get_mean(DELTA_IC);
#ifdef _FULL_VERBOSE_
          this->So.message_screen("Statistics of the original IC");
          this->So.message_screen("Mean delta_IC =",meanWN);
#endif
          sigma2D=get_var(meanWN, DELTA_IC);
#ifdef _FULL_VERBOSE_
          this->So.message_screen("Sigma² IC =",sigma2D);
#endif
          //-------------------->

          //--------------------> If the input overdensity IC has an amplitude at a z different from the IC of the simulations, we have to normalize to that IC
          // E.g, the IC can be normalized at z = 1, while the simulation has started at z=99
          if(true==this->params._Normalize_IC_to_initial_redshift())
            normalize_df_z_ini(DELTA_IC,"delta", this->params._Initial_Redshift_DELTA(), this->params._Initial_Redshift_SIM());

          //--------------------> Normalize now from the redshift of the initial conditions to redshift z=0. LPT starts from z=0
          normalize_df_z_ini(DELTA_IC,"delta",this->params._Initial_Redshift_SIM(), 0.0);

           //-------------------->  Here we write the alrady normalized to z=0 Initial DF and hereafter we will read it
          this->File.write_array(this->fnameICDELTA,DELTA_IC);
          //--------------------> This line will let us read the file written in the line above whenever we need the initial delta field below
          this->params.set_use_ic_file(false);

         //--------------------> This is meant to measure some statistis of the IC with amplitude at inital redshift of simulation
#ifdef _FULL_VERBOSE_
          So.message_screen("Measuring power spectrum of IC read in FILE",  this->fnameICDELTA);
#endif
          this->params.set_Name_survey("IC");
          this->params.set_input_type("delta_grid");
          PowerSpectrumF cPSF(this->params);
          cPSF.compute_power_spectrum_grid(DELTA_IC,true);
          this->So.DONE();
         //-------------------->

         //--------------------> GNUPLOT
#ifdef _USE_GNUPLOT_
          vector<pair<real_prec, real_prec> > xy_pts_ref;

#ifdef _USE_GNUPLOT_PDF_
          int Nk=300;
          vector<pair<real_prec, real_prec> > xy_pts_g;
          vector<real_prec>xbaux(Nk, 0);
          vector<real_prec>pdf(Nk, 0);
          real_prec max=get_max<real_prec>(DELTA_IC);
          real_prec min=get_min<real_prec>(DELTA_IC);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0;i<Nk; ++i)
            xbaux[i]=min+static_cast<real_prec>(i+0.5)*(max-min)/static_cast<real_prec>(Nk);
          string filex=this->params._Output_directory()+"pdf_IC"+to_string(this->params._seed())+".txt";
          calc_pdf("lin", dummy.size(),Nk,max,min,DELTA_IC,pdf);
          this->File.write_to_file(filex, xbaux,pdf);
          for(int i=1; i<Nk;++i) 
            xy_pts_ref.push_back(std::make_pair(xbaux[i], pdf[i])); 
          pdf.clear(); pdf.shrink_to_fit();

          for(int i=0; i<Nk;++i) 
            xy_pts_g.push_back(std::make_pair(xbaux[i], exp(-0.5*xbaux[i]*xbaux[i]/sigma2D)/sqrt(2.*M_PI*sigma2D))); 
          xbaux.clear(); xbaux.shrink_to_fit();

          this->gp_pdf<<"set border linewidth 1.5\n";
          this->gp_pdf<<"set border linecolor '"<<FG_COLOR<<"' \n";
          this->gp_pdf<<"set key textcolor rgb '"<<FG_COLOR<<"'\n";
          this->gp_pdf<<"set xtics textcolor '"<<FG_COLOR<<"' \n";
          this->gp_pdf<<"set ytics textcolor '"<<FG_COLOR<<"' \n";
          this->gp_pdf << "set log y\n";
          this->gp_pdf << "set xlabel 'log (1+{/Symbol d})' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"'\n";
          this->gp_pdf << "set ylabel 'P({/Symbol d})' font 'Times-Roman,10' textcolor rgb '"<<FG_COLOR<<"'\n";
          this->gp_pdf << "set size square\n";
          this->gp_pdf << "plot" << gp_pdf.file1d(xy_pts_ref) << "w l lw 2 lt 1 title 'PDF IC'"<<endl;
          this->gp_pdf << "replot" << gp_pdf.file1d(xy_pts_g) << "w l lw 2 lt 2 title 'N(0,1)'"<<endl;
          xy_pts_ref.clear();
          xy_pts_ref.shrink_to_fit();
          xy_pts_g.clear();
          xy_pts_g.shrink_to_fit();
#endif
         string power_th=this->params._dir()+this->params._ic_power_file();
         for(ULONG i=1; i<cPSF._kvector_data_size(); ++i) 
            xy_pts_ref.push_back(std::make_pair(cPSF._kvector_data(i), cPSF._pk0(i))); 
          this->gp_power << "replot" << gp_power.file1d(xy_pts_ref) << "w l lw 2 lt 1 title 'IC @ z = 0'"<<endl;
#ifdef _USE_IC_INPUT_POWER_DELTA_
          this->gp_power << "replot '"<<power_th<<"' u 1:($2*(2*acos(-1.)*acos(-1.0))/($1*$1*$1)) w l lw 2 lt 3 title 'Theoretical power'"<<endl;
#else
#ifdef _correct_shape_theoretica_power_
         this->gp_power << "replot '"<<power_th<<"' u 1:($2/(1.+0.05*$1**2)) w l lw 2 lt 3 title 'Theoretical power'"<<endl;
#else
         this->gp_power << "replot '"<<power_th<<"' w l lw 2 lt 3 title 'Theoretical power'"<<endl;

#endif
#endif          
          xy_pts_ref.clear();
          xy_pts_ref.shrink_to_fit();
#endif

        }

      real_prec biasLAG=static_cast<real_prec>(1.);
#ifdef _FULL_VERBOSE_
      So.message_screen("Computing Displacement Field");
#endif

#ifdef SAVEIC
      vector<real_prec>posx(NLL,0),posy(NLL,0),posz(NLL,0);
      vector<real_prec>posicx(NLL,0),posicy(NLL,0),posicz(NLL,0);
      vector<real_prec>dummy2(N,0);
#ifdef PARTRAN
      So.message_screen("Writting initial Positions");
      this->File.write_array(fnamePOSX,posx);
      this->File.write_array(fnamePOSY,posy);
      this->File.write_array(fnamePOSZ,posz);
      getDensity_NGP(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,posx,posy,posz,dummy,dummy2,false);
      this->File.write_array(fnameDMNGP,dummy2);
#endif //endif PARTRAN

#endif  //endif SAVEIC

#ifndef SAVEIC

#define POSX dummy
      this->Lag2Eul_comp(biasLAG,kth,ftype,true,POSX,1);
      // ABA: velocities are written inside, but we can do as with the positions
      this->File.write_array(this->fnamePOSX,POSX);

#define POSY dummy
      this->Lag2Eul_compB(kth,ftype,true,POSY,2);
      this->File.write_array(this->fnamePOSY, POSY);

#define POSZ dummy
      this->Lag2Eul_compB(kth,ftype,true,POSZ,3);
      this->File.write_array(this->fnamePOSZ, POSZ);

      //----------------> Compute the components of the dm tracer velocities. The differet components are writtn inside this function:
      this->Lag2Eul_vel(biasLAG,kth,ftype,true);


#ifdef _MOVE_DM_TO_REDSHIFT_SPACE_
      //----------------> This will modify POSX, POSY OR POSZ
      //----------------> ACCORDING TO LOS (line of sight)
     this->DM_to_RSS(LOS);
     //----------------> This has to be done snice below the Z coordinate is read from memory
     if(3==LOS)
       this->File.read_array(this->fnamePOSZ_RS+".dat",POSZ);
#endif

      vector<real_prec>dens_field_cic;// container for the interpoation of the dhnumber count used in the construction of the velocity field
      {
        vector<real_prec>dummy2(N,0);
        vector<real_prec>dummy3(N,0);
        vector<real_prec>dummy4(N,0);

#undef POSX
#define POSX dummy3

#undef POSY
#define POSY dummy4

#define DENS_FIELD dummy2

        // read to avoid keeping it in memmory. Z-compònent is the only one kept in memmory
#ifdef _MOVE_DM_TO_REDSHIFT_SPACE_
        if(1==LOS)
          {
          this->File.read_array(this->fnamePOSX_RS+".dat",POSX);
          this->File.read_array(this->fnamePOSY+".dat",POSY);
        }
        else if(2==LOS){
          this->File.read_array(this->fnamePOSX+".dat",POSX);
          this->File.read_array(this->fnamePOSY_RS+".dat",POSY);
        }
#else
        this->File.read_array(this->fnamePOSX+".dat",POSX);
        this->File.read_array(this->fnamePOSY+".dat",POSY);
#endif


        getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,POSZ,DENS_FIELD,false);
        this->File.write_array(this->fnameDMNGP, DENS_FIELD);
        {
          this->params.set_input_type("density_grid");
          this->params.set_MAS_correction(true);
          if(this->params._Structure_Formation_Model()==1)
            this->params.set_Name_survey("DM_ZELD");
          else if(this->params._Structure_Formation_Model()==2)
            this->params.set_Name_survey("DM_2LPT");
          else if(this->params._Structure_Formation_Model()==3)
           this->params.set_Name_survey("DM_ALPT");

          PowerSpectrumF cPSF(this->params);
          cPSF.compute_power_spectrum_grid(DENS_FIELD,true);
          this->So.DONE();

#ifdef _USE_GNUPLOT_
          vector<pair<real_prec, real_prec> > xy_pts_new;
          for(int i=1; i<cPSF._kvector_data_size(); ++i) 
            xy_pts_new.push_back(std::make_pair(cPSF._kvector_data(i), cPSF._pk0(i))); 
          this->gp_power << "set title 'LOS"<< this->params._seed()<<"'\n";
          if(this->params._Structure_Formation_Model()==1)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 4 title 'DM-ZELD z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          if(this->params._Structure_Formation_Model()==2)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 4 title 'DM-2LPT z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          if(this->params._Structure_Formation_Model()==3)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 4 title 'DM-ALPT z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          xy_pts_new.clear();
          xy_pts_new.shrink_to_fit();
#endif
        }




#ifdef _GET_DENSITY_FIELDS_
#ifdef _USE_TET_
#ifdef _FULL_VERBOSE_
        So.message_screen("Going TET density field");
#endif
        // Here I have explicitely exchanged x<->y when using TET.
        getDensity_TETCIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,NLL,DENS_FIELD);
#else

        switch (this->params._masskernel())
          {
          case 0:
            getDensity_NGP(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,POSZ,DENS_FIELD,false);
            break;
          case 1:
            getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSY,POSX,POSZ,POSZ,DENS_FIELD,false);
            break;
          case 2:
            getDensity_TSC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,POSZ,DENS_FIELD,false);
            break;
          }
#endif

        this->File.write_array(this->fnameDM, DENS_FIELD);  //write DMDF to binary file



        {
          this->params.set_input_type("density_grid");
          this->params.set_MAS_correction(true);
       if(this->params._Structure_Formation_Model()==1)
            this->params.set_Name_survey("DM_ZELD");
          else if(this->params._Structure_Formation_Model()==2)
            this->params.set_Name_survey("DM_2LPT");
          else if(this->params._Structure_Formation_Model()==3)
           this->params.set_Name_survey("DM_ALPT");
          PowerSpectrumF cPSF(this->params);
          cPSF.compute_power_spectrum_grid(DENS_FIELD,true);
          this->So.DONE();

#ifdef _USE_GNUPLOT_
          vector<pair<real_prec, real_prec> > xy_pts_new;
          string filex=this->params._Output_directory()+"pdf_DM"+this->buffsf+"_IC"+to_string(this->params._seed())+".txt";
          int Nk=300;
          vector<real_prec>xbaux(Nk, 0);
          vector<real_prec>pdf(Nk, 0);
          real_prec max=get_max<real_prec>(DENS_FIELD);
          real_prec min=get_min<real_prec>(DENS_FIELD);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(int i=0;i<Nk; ++i)
            xbaux[i]=min+static_cast<real_prec>(i+0.5)*(max-min)/static_cast<real_prec>(Nk);

          vector<pair<real_prec, real_prec> > xy_pts_g;
          calc_pdf("lin", dummy.size(),Nk,max,min,DENS_FIELD,pdf);
          this->File.write_to_file(filex, xbaux,pdf);
#ifdef _USE_GNUPLOT_PDF_
          for(int i=1; i<Nk;++i) 
            xy_pts_new.push_back(std::make_pair(xbaux[i], pdf[i])); 
          pdf.clear(); pdf.shrink_to_fit();
          this->gp_pdf << "replot" << gp_pdf.file1d(xy_pts_new) << "w l lw 2 lt 5 title 'DM-ALPT+TET'"<<endl;
          xy_pts_new.clear();
          xy_pts_new.shrink_to_fit();
#endif

          for(int i=1; i<cPSF._kvector_data_size(); ++i) 
            xy_pts_new.push_back(std::make_pair(cPSF._kvector_data(i), cPSF._pk0(i))); 
          if(this->params._Structure_Formation_Model()==1)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 5 title 'DM-ZELD-TET z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          if(this->params._Structure_Formation_Model()==2)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 5 title 'DM-2LPT-TET z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          if(this->params._Structure_Formation_Model()==3)
            this->gp_power << "replot" << gp_power.file1d(xy_pts_new) << "w l lw 2 lt 5 title 'DM-ALPT-TET z= "<<static_cast<float>(this->params._redshift())<<"'"<<endl;
          xy_pts_new.clear();
          xy_pts_new.shrink_to_fit();
#endif
        }



//        dens_field_cic.resize(this->NGRID,0);// container for the interpoation of the dhnumber count used in the construction of the velocity field
//        getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSY,POSX,POSZ,POSZ,dens_field_cic,false);

#endif




#ifdef _GET_VELOCITIES_
        // Now the velocities:
        vector<real_prec>dummy5(N,0);
        vector<real_prec>dummy6(N,0);

        // x:

#define VXp dummy5
#define VX  dummy6


        this->File.read_array(this->fnameVXpart+".dat",VXp);

        switch (this->params._masskernel_vel())
          {
          case 0:
            getDensity_NGP(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VXp,VX,true);
            break;
          case 1:
            getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSY,POSX,POSZ,VXp,VX,true);
            break;
          case 2:
            getDensity_TSC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VXp,VX,true);
            break;
          }

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          if(dens_field_cic[i]>0.)
            VX[i]/=static_cast<real_prec>(DENS_FIELD[i]);
          else
            VX[i]=0;  // Technically one would need to assign thevelocity of the nearest particle, for the V field can be different of zero even if there are no particles.
        this->File.write_array(this->fnameVX, VX);  //write DMDF to binary file

        //y:
#define VYp dummy5
#define VY  dummy6

        this->File.read_array(this->fnameVYpart+".dat",VYp);

        switch (this->params._masskernel_vel())
          {
          case 0:
            getDensity_NGP(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VYp,VY,true);
            break;
          case 1:
            getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSY,POSX,POSZ,VYp,VY,true);
            break;
          case 2:
            getDensity_TSC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VYp,VY,true);
            break;
          }

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          if(dens_field_cic[i]>0.)
            VY[i]/=static_cast<real_prec>(DENS_FIELD[i]);
          else
            VY[i]=0;
        this->File.write_array(this->fnameVY, VY);  //write DMDF to binary file



        // z:
#define VZp dummy5
#define VZ  dummy6

        this->File.read_array(this->fnameVZpart+".dat",VZp);

        switch (this->params._masskernel_vel())
          {
          case 0:
            getDensity_NGP(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VZp,VZ,true);
            break;
          case 1:
            getDensity_CIC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSY,POSX,POSZ,VZp,VZ,true);
            break;
          case 2:
            getDensity_TSC(N1,N2,N3,L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,POSX,POSY,POSZ,VZp,VZ,true);
            break;
          }

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          if(dens_field_cic[i]>0.)
            VZ[i]/=static_cast<real_prec>(DENS_FIELD[i]);
          else
            VZ[i]=0;
        this->File.write_array(this->fnameVZ, VZ);  //write DMDF to binary file

#endif


#undef DENS_FIELD

      }
#endif

      // #ifdef DMONLY
      //       exit(0);
      // #endif

    }
  // end running simulation (1)// Here it is produced a DM density field, that goes into BAM
  // The positions of the DM field have to be saved/printed in an output file that has to be read
  // by bam when assigning position withion cells.

}


#endif // end of ifndef displacements. I have to do it because some functions change arguments if displacements is used 


//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################
//##################################################################################

void LPT::displace_part_comp(vector<real_prec>&posi,vector<real_prec>&psii,bool periodic,int comp)
{

    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);

  this->So.enter(__PRETTY_FUNCTION__);

#ifdef _FULL_VERBOSE_
  this->So.message_screen("Moving particles from Lagrangian space to Eulerian space ");
#endif

  real_prec posii;

  ULONG nout=0;
  ULONG jobj=0;
  ULONG N=this->NGRID;

  ULONG i,j,k;
  real_prec posxi,posyi,poszi,psi,inew,rx,ry,rz;
  bool bout;

#ifdef _USE_OMP_
#pragma omp parallel for default(none) private(i,j,k,jobj,rx,ry,rz,posxi,posyi,poszi,posii,psi,inew,bout)  shared (periodic,comp,posi,psii)  reduction(+:nout)
#endif
  for (i=0;i<this->N1;i++)
    for (j=0;j<this->N2;j++)
      for (k=0;k<this->N3;k++)
        {
          jobj = index_3d(i,j,k,this->N2,this->N3);


          // This is where we fixe the initial positions (Lagrangian coordinates). Subtracting half bin here means that these particles are placed at the corner (lattice picture)
          posxi=this->params._d1()*(static_cast<real_prec>(i)+0.5-BIN_SHIFT);
          posyi=this->params._d2()*(static_cast<real_prec>(j)+0.5-BIN_SHIFT);
          poszi=this->params._d3()*(static_cast<real_prec>(k)+0.5-BIN_SHIFT);

          switch (comp)
            {
            case 1:
              posii=posxi;
              break;
            case 2:
              posii=posyi;
              break;
            case 3:
              posii=poszi;
              break;
            }
          //* These are the coordinates in Lagrangian space (Initial conditions) refered aither to  the center
          // (cell picture) or to the corners (lattice picture) * //
          //*  If ABACUS is enabled, 0.5-F=0 means that we use the lattice-like description.  * //
          //* The lininterp functio will be accordingly modified. *//

          psi=this->linInterp(posxi,posyi,poszi,psii);


          posi[jobj]=posii+psi;

          if (periodic==true)
            {
              inew=posi[jobj];

              if (inew<0.)
                inew+=this->L1;
              if (inew>=this->L1)
                inew-=this->L1;

              posi[jobj]=inew;
            }
          bout=false;
          if (posi[jobj]<0 || posi[jobj]>=this->L1)
            {
              nout++;
              bout=true;
            }
          //	  jobj++;
        }

  this->So.DONE();
#ifdef _FULL_VERBOSE_
  if(nout>0)
    this->So.message_screen("Number of objects outside boundaries =",nout);
#endif

}

// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ========================================================================== gf*=factor_growth; ============================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================





#ifndef _DISPLACEMENTS_ // When studying the displacement field, we do not need to move particles. We only want the displacement
void LPT::Lag2Eul_comp(real_prec biasLAG2,real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp)
#else
void LPT::Lag2Eul_comp(real_prec biasLAG2,real_prec kth,int ftype,bool periodic,int comp)
#endif

{

#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

 So.enter(__PRETTY_FUNCTION__);

#ifdef _FULL_VERBOSE_
  this->So.message_screen("Going Lagrangian to Eulerian, ", comp);
#endif
  real_prec D1=this->s_cosmo_info.growth_factor;

 #ifdef _UNITSIM_
    D1*=factor_growth;
#endif

  real_prec D2=this->s_cosmo_info.D2;
  real_prec f1=this->s_cosmo_info.growth_index;
  real_prec f2=this->s_cosmo_info.growth_index2;
  real_prec Dcurl=static_cast<real_prec>(-1./7.)*D1*D1*D1;
  real_prec D3a=static_cast<real_prec>(-1./3.)*D1*D1*D1;
  real_prec D3b=static_cast<real_prec>(10./21.)*D1*D1*D1;
  real_prec cvel=num_1/(cgs_km/cgs_Mpc);
  real_prec rsmolC=this->params._d1();//static_cast<real_prec>(1.0);//kth


  kernelcomp(this->L1, this->L2,this->L3,this->params._d1(),this->params._d2(),this->params._d3(),this->N1,this->N2,this->N3,rsmolC,ftype,this->params._Output_directory());

  string fname;

  bool zeld=false;
  bool twolpt=false;
  bool alpt=false;

  switch (this->params._Structure_Formation_Model())
    {
    case 1:
      zeld=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using Zel'dovich");
#endif
      break;
    case 2:
      twolpt=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using 2LPT");
#endif
      break;
    case 3:
      alpt=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using ALPT");
#endif
      break;
    }

  ULONG N=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3);
  ULONG Nhalf=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3/2+1);

  vector<real_prec>dummy2(N,0);
  vector<real_prec>dummy(N,0);

  // Here we read the initial overdenisty field. I add *dat, otherwise it will read the white noise (perhaps we can change names to avoid this confusion)

#define DELTA_IC dummy
#define POTENTIAL dummy


  if (true==zeld)
    {
      this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);// Get IC
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)               // get g(z) * delta 
        DELTA_IC[i]*=D1;

#define dPSI_TOT dummy   // when the sf ifs are close, the code expects a dPSI_TOT. BUt since this 

    }

  if (true==twolpt || true==alpt)
    {

      this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);

      // delta(1)-> Phi(1), needed to obtaine delta(2)
      PoissonSolver(this->L1,this->N1,DELTA_IC,POTENTIAL);// here dummy is overwritten, now is POTENTIAL

#ifdef PTCURL
      fname=string("auxs1");
      this->File.write_array(fname,POTENTIAL);
#endif

      // Phi(1)  -> delta(2)

#define DELTA_TWO dummy2

#ifndef _DISPLACEMENTS_
      ifstream inStream;
      string ffnn=this->fname2LPTTERM+".dat";
      inStream.open(ffnn.data());
      if (false==inStream.is_open())
        {
#endif
          calc_twolptterm(N1,N2,N3,L1,L2,L3,POTENTIAL,DELTA_TWO);
#ifndef _DISPLACEMENTS_
          this->File.write_array(this->fname2LPTTERM,DELTA_TWO);
        }
      else
        {
          fname=this->fname2LPTTERM+".dat";
          this->File.read_array(fname,DELTA_TWO);
        }
#endif

     this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);

#define PSI_LPT dummy
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        PSI_LPT[i]=D1*DELTA_IC[i]-D2*DELTA_TWO[i];	//attention!!!


#ifdef PTCURL
      fname=string("auxs2");

      this->File.read_array(fname,PSI_LPT);

#define POT_PSI_LPT dummy
      PoissonSolver(this->L1,this->N1,PSI_LPT,POT_PHI_LPT);


      vector<real_prec> dummy3(N,0);
#define CURL1 dummy3

      fname=string("auxs1");
      this->File.read_array(fname,POT_PHI_LPT);


      calc_curlcomp(N1,N2,N3,L1,L2,L3,POT_PHI_LPT,DELTA_TWO,CURL1,1);
      fname=string("curl1");
      this->File.write_array(fname,CURL1);

#define CURL2 dummy3
      calc_curlcomp(N1,N2,N3,L1,L2,L3,POT_PHI_LPT,DELTA_TWO,CURL2,2);
      fname=string("curl2");
      this->File.write_array(fname,CURL2);

#define CURL3 dummy3
      calc_curlcomp(N1,N2,N3,L1,L2,L3,POT_PHI_LPT,DELTA_TWO,CURL3,3);
      fname=string("curl3");
      this->File.write_array(fname,CURL3);


      fname=string("auxs2");
      this->File.read_array(fname,dummy);

      /**/
#ifdef LPT3B
      fname=string("auxs1");
      this->File.read_array(fname,dummy);

#define MU2 dummy3
      calc_mu2term(N1,N2,N3,L1,L2,L3,dummy,DELTA_TWO,MU2);

      fname=string("auxs2");
      this->File.read_array(fname,dummy);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        dummy[i]-=D3b*MU2[i];
      fname=string("auxs2");
      this->File.write_array(fname,dummy);
#endif // end LPT2B

#ifdef LPT3A
     fname=string("auxs1");
     this->File.read_array(fname,dummy);
     calc_Det(N1,N2,N3,L1,L2,L3,dummy,dummy3);
     fname=string("auxs2");
     this->File.read_array(fname,dummy);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for(ULONG i=0;i<N;i++)
       dummy[i]-=D3a*dummy3[i];
#endif // end LPT3A
    /**/
#endif // end LPT3B

//#define dPSI_TOT dummy   // this was alrady defined in if(true==zeld), ao that the compilor already saw this #def

    }

  if (true==alpt)
    {
#ifdef MUSCLE2LPT
      vector<real_prec> dummy3(N,0);

      int Nsmol=3;
      real_prec drsmol=static_cast<real_prec>(8.);
      real_prec rsmin=this->params._d1();

      for(int i=0;i<Nsmol;i++)
        {
          real_prec rsmol=static_cast<real_prec>(i)*drsmol+rsmin;
          kernelcomp(this->L1,this->L2,this->L3,this->params._d1(),this->params._d2(),this->params._d3(),this->N1,this->N2,this->N3,rsmol,1,this->params._Output_directory());
          convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dummy,dummy3,ftype,rsmol,this->params._Output_directory());
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
          for(ULONG i=0;i<N;i++)
            {
              real_prec psilin=dummy3[i];

              if (static_cast<real_prec>(1.+2./3.*psilin) < 0.)
                dummy[i]=3.;
            }
        }
#endif

      //Convolution of Kernel K with divergence of (Psi^2LPT)
      convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dummy,dummy,ftype,kth,this->params._Output_directory());



#undef DELTA_IC
#define DELTA_IC dummy2
      //      if(false==this->params._use_ic_file())
      this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);
      // else
      // 	{
      // 	  this->File.read_array(this->params._ic_file(),DELTA_IC);
      // 	  get_overdens(DELTA_IC, DELTA_IC);
      // 	  if(true==this->Normalize_initial_redshift)
      // 	    normalize_df_z_ini(DELTA_IC, DELTA_IC,"delta");
      // 	}

      /* Psi^tot=K o Psi^2LPT + Psi^SC - K o Psi^SC step by step*/

      // First: div Psi^SC

#define dPSI_SC dummy2

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        {
          real_prec psilin=-D1*DELTA_IC[i];
          real_prec psisc=0.;
          real_prec insidesqrt=static_cast<real_prec>(1.+2./3.*psilin);
          if (insidesqrt > 0.)
            psisc=static_cast<real_prec>(3.*(sqrt(insidesqrt)-1.));
          else
            psisc=-3.;

          psisc*=static_cast<real_prec>(-1.);
          dPSI_SC[i]=psisc;
        }


      //Second : K o div v^2LPT + div v^SC
#define K_CONV_dPSI_2LPT dummy

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        dummy[i]=K_CONV_dPSI_2LPT[i]+ dPSI_SC[i];

      //K o div v^SC
#define K_CONV_dPSI_SC dummy2
      convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dPSI_SC,K_CONV_dPSI_SC,ftype,kth,this->params._Output_directory());

//#define dPSI_TOT dummy // this was already defined in the if(true== zeld), so that the compilor already saw taht definition
      //K o div v^2LPT + div v^SC - K o div v^SC
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        dPSI_TOT[i]=dummy[i]- K_CONV_dPSI_SC[i];

      if (this->params._vslength()>0.)
        convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dPSI_TOT,dPSI_TOT,1,this->params._vslength(),this->params._Output_directory());

    } // close parenthesis over the type of sf model. From here dummy represents the Displacement field

/*
#ifdef _POWER_BOOST_ALPT_
#ifdef _FULL_VERBOSE_
  So.message_screen("Applying power-boost to div of Psi");
#endif
    vector<real_prec>KERNEL_POWER(Nhalf,1.0);
    kernelcomp_for_boost(L1,d1,N1,KERNEL_POWER, this->params._Output_directory());
    convolvek(N1,dPSI_TOT,KERNEL_POWER,dPSI_TOT);
    KERNEL_POWER.clear();KERNEL_POWER.shrink_to_fit();
#endif
*/


  // write divergence of the displacement field
#ifndef _DISPLACEMENTS_  // / When studying the displacement field, we do not need to write this file 
  this->File.write_array(this->fnameTHETA,dPSI_TOT);
#endif

  if (this->params._velbias()>0.)
    comp_velbias(dPSI_TOT,dPSI_TOT,false,true);


  // Now start component by component. This function deasl with the x-component, for
  // it also computes the theta field. The function *B deals with the other two components, reading the theta field

#ifndef _DISPLACEMENTS_
  if (comp==1)   // The other two components are computed in the function compB, for it reads files already computed, saving time.
    {
#endif

      /* Psi^tot x-component */
      //      this->So.message_screen("Getting Psi_tot, x-comp");
#define DISPLACEMENT dummy2
      theta2velcomp(dPSI_TOT,DISPLACEMENT,false,false,comp);
      this->So.DONE();



#ifdef CELLBOUND
      cellboundcomp(N1,N2,N3,VEL);
#endif

#ifdef PTCURL
      {
        vector<real_prec> dummy3(N,0);
        fname=string("curl1");
        this->File.read_array(fname,dummy);

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          dummy3[i]=Dcurl*PSI_TOT[i];

#ifdef SMOOTHCURL
        convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dummy3,dummy3,ftype,rsmolC,this->params._Output_directory());
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<N;i++)
          DISPLACEMENT[i]+=dummy3[i];
        dummy3.clear();dummy3.shrink_to_fit();
      }
#endif


#ifdef _DISPLACEMENTS_
      this->Displacement=DISPLACEMENT;
#endif


#ifndef _DISPLACEMENTS_
#ifdef _WRITE_DISPLACEMENTS_
      string disp_file=this->params._Output_directory()+"Displacement_comp_"+to_string(comp)+"_"+buffsf;
      this->File.write_array(disp_file, DISPLACEMENT);
#endif
#endif

#ifndef _DISPLACEMENTS_ // When studying the displacement field, we do not need to move particles. We only want the displacement
      this->displace_part_comp(out,DISPLACEMENT,periodic,comp);
#endif

#ifndef _DISPLACEMENTS_
   }
#endif

#undef DELTA_IC

  dummy.clear();dummy.shrink_to_fit();
  dummy2.clear();dummy2.shrink_to_fit();

}


// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================





void LPT::Lag2Eul_vel(real_prec biasLAG2,real_prec kth,int ftype,bool periodic)
{

#ifdef _USE_OMP_
    int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

#ifdef _FULL_VERBOSE_
  this->So.message_screen("Getting particle velocities");
#endif
  real_prec D1=this->s_cosmo_info.growth_factor;
 #ifdef _UNITSIM_
    D1*=factor_growth;
#endif

  real_prec D2=this->s_cosmo_info.D2;
  real_prec f1=this->s_cosmo_info.growth_index;
  real_prec f2=this->s_cosmo_info.growth_index2;
  real_prec Dcurl=static_cast<real_prec>(-1./7.)*D1*D1*D1;
  real_prec cvel=num_1/(cgs_km/cgs_Mpc);
  real_prec rsmolC=this->params._d1();//static_cast<real_prec>(1.0);//kth

  kernelcomp(this->L1, this->L2,this->L3,this->params._d1(),this->params._d2(),this->params._d3(),this->N1,this->N2,this->N3,rsmolC,ftype,this->params._Output_directory());

  string fname;

  bool zeld=false;
  bool twolpt=false;
  bool alpt=false;

  switch (this->params._Structure_Formation_Model())
    {
    case 1:
      zeld=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using Zel'dovich");
#endif
      break;
    case 2:
      twolpt=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using 2LPT");
#endif
      break;
    case 3:
      alpt=true;
#ifdef _FULL_VERBOSE_
      So.message_screen("Using ALPT");
#endif
      break;
    }

  ULONG N=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3);
  ULONG Nhalf=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3/2+1);

  vector<real_prec>dummy2(N,0);
  vector<real_prec>dummy(N,0);

  // Here we read the initial overdenisty field. I add *dat, otherwise it will read the white noise (perhaps we can change names to avoid this confusion)

#define DELTA_IC dummy
#define POTENTIAL dummy

  this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);

  if (true==zeld)
    {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        DELTA_IC[i]*=D1;
    }

  if (true==twolpt || true==alpt)
    {
      // delta(1)-> Phi(1), needed to obtaine delta(2)
      PoissonSolver(this->L1,this->N1,DELTA_IC,POTENTIAL);


      // Phi(1)  -> delta(2)
#define DELTA_TWO dummy2
      ifstream inStream;
      string ffnn=this->fname2LPTTERM+".dat";
      inStream.open(ffnn.data());

      if (false==inStream.is_open())
        {
          calc_twolptterm(N1,N2,N3,L1,L2,L3,POTENTIAL,DELTA_TWO);
          this->File.write_array(this->fname2LPTTERM,DELTA_TWO);
        }
      else
        {
          fname=this->fname2LPTTERM+".dat";
          this->File.read_array(fname,DELTA_TWO);
        }

        this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);


#define PSI_LPT dummy
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        PSI_LPT[i]=D1*DELTA_IC[i]-D2*(f2/f2)*DELTA_TWO[i];
        //attention!!! The original form for the velocity field has the form f1HaD1ð - f2HaD2ð2. Here we only compute (D1ð - f2D2ð2/f1). The multiplicative factor is added in the
        // function theta2velcomp(dPSI_TOT,VELX_part,false,true,1) with true, such that inside that function the term f1(z)H(z)a is added.
        // The term Vsc has also a multiplicative factor f1(z)H(z)a,
        // The term D1ð - D2ð2 has units of Mpc/h, such that when multiplying it by H /(f and a are dimensionless) we obtain units of km/s.

    }

  if (true==alpt)
    {
      //Convolution of Kernel K with divergence of (Psi^2LPT)
      convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,PSI_LPT,PSI_LPT,ftype,kth,this->params._Output_directory());

#undef DELTA_IC
#define DELTA_IC dummy2
      this->File.read_array(this->fnameICDELTA+".dat",DELTA_IC);

      /* Psi^tot=K o Psi^2LPT + Psi^SC - K o Psi^SC step by step*/
      // First: div Psi^SC

#define dPSI_SC dummy2
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        {
          if(D1*DELTA_IC[i] <3./2.)
            dPSI_SC[i]=D1*DELTA_IC[i]/static_cast<real_prec>(sqrt(1.-(2./3.)*D1*DELTA_IC[i]));
          else
            dPSI_SC[i]=0.0;
          }

      //Second : K o div v^2LPT + div v^SC
#define K_CONV_dPSI_2LPT dummy

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        dummy[i]=K_CONV_dPSI_2LPT[i] + dPSI_SC[i]; // dummy+=dummy2

      //K o div v^SC
#define K_CONV_dPSI_SC dummy2
      convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dPSI_SC,K_CONV_dPSI_SC,ftype,kth,this->params._Output_directory());


#define dPSI_TOT dummy
      //K o div v^2LPT + div v^SC - K o div v^SC
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        dPSI_TOT[i]=dummy[i]- K_CONV_dPSI_SC[i];  //dummy-=dummy2


      if (this->params._vslength()>0.)
        convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dPSI_TOT,dPSI_TOT,1,this->params._vslength(),this->params._Output_directory());

    } // close parenthesis over the type of sf model. From here dummy represents the Displacement field


  // write divergence of the displacement field
 // this->File.write_array(this->fnameTHETA,dPSI_TOT);

  if (this->params._velbias()>0.)
    comp_velbias(dPSI_TOT,dPSI_TOT,false,true);


  // Now start component by component. This function deasl with the x-component, for
  // it also computes the theta field. The function *B deals with the other two components, reading the theta field


  /* VEL x-component */
  //      this->So.message_screen("Getting Psi_tot, x-comp");
#define VELX_part dummy2
  theta2velcomp(dPSI_TOT,VELX_part,false,true,1);
  this->File.write_array(this->fnameVXpart, VELX_part);
#undef VELX_part

  /* VEL y-component */
  //      this->So.message_screen("Getting Psi_tot, x-comp");
#define VELY_part dummy2
  theta2velcomp(dPSI_TOT,VELY_part,false,true,2);
  this->File.write_array(this->fnameVYpart, VELY_part);
#undef VELY_part

  /* VEL z-component */
  //      this->So.message_screen("Getting Vel, x-comp");
#define VELZ_part dummy2
  theta2velcomp(dPSI_TOT,VELZ_part,false,true,3);
  this->File.write_array(this->fnameVZpart, VELZ_part);
#undef VELZ_part
  this->So.DONE();

  dummy.clear();dummy.shrink_to_fit();
  dummy2.clear();dummy2.shrink_to_fit();



}
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================






void LPT::Lag2Eul_compB(real_prec kth,int ftype,bool periodic,vector<real_prec> &out,int comp)
{

#ifdef _FULL_VERBOSE_
  this->So.message_screen("Going Lagrangian to Eulerian, ", comp);
#endif

#ifdef _USE_OMP_
       int NTHREADS=_NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  real_prec D1=this->s_cosmo_info.growth_factor;
 #ifdef _UNITSIM_
    D1*=factor_growth;
#endif

  real_prec D2=this->s_cosmo_info.D2;
  real_prec Dcurl=static_cast<real_prec>(-1./7.)*D1*D1*D1;
  real_prec D3a=static_cast<real_prec>(-1./3.)*D1*D1*D1;
  real_prec D3b=static_cast<real_prec>(10./21.)*D1*D1*D1;

  real_prec cvel=num_1/(cgs_km/cgs_Mpc);

  real_prec rsmolC=this->params._d1();//static_cast<real_prec>(1.0);//kth

  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  vector<real_prec>  dummy(N,0), dummy2(N,0);
  string fname;

  this->File.read_array(this->fnameTHETA+".dat",dummy);

  if (this->params._velbias()>0.)
    comp_velbias(dummy,dummy,false,true);

  if (comp==2)
    {
      /* Psi^tot y-component */
#define DISPLACEMENT2 dummy2

      this->theta2velcomp(dummy,DISPLACEMENT2,false,false,comp);
#ifdef CELLBOUND
      cellboundcomp(N1,N2,N3,dummy2);
#endif

#ifdef PTCURL
      {
        vector<real_prec> dummy3(N,0);

        fname=string("curl2");

        this->File.read_array(fname,dummy);

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          dummy3[i]=Dcurl*dummy[i];


#ifdef SMOOTHCURL
        convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dummy3,dummy3,ftype,rsmolC,this->params._Output_directory());
#endif

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          dummy2[i]+=dummy3[i];
      }
#endif

#ifndef _DISPLACEMENTS_
      this->displace_part_comp(out,DISPLACEMENT2,periodic,comp);
#endif


#ifdef _DISPLACEMENTS_
      this->Displacement=DISPLACEMENT2;
#endif

#ifdef _WRITE_DISPLACEMENTS_
      string disp_file=this->params._Output_directory()+"Displacement_comp_"+to_string(comp)+"_"+buffsf;
      this->File.write_array(disp_file, DISPLACEMENT2);
#endif



    }

  if (comp==3) // esto sobra, ya usando comp en todas partes se hace el trabajo si la función es llamada para cada componente por separado
    {
      /* Psi^tot z-component */
#define DISPLACEMENT3 dummy2
          theta2velcomp(dummy,DISPLACEMENT3,false,false,3);
#ifdef CELLBOUND
      cellboundcomp(N1,N2,N3,dummy2);
#endif

#ifdef PTCURL
      {
        vector<real_prec> dummy3(N,0);

        fname=string("curl3");
        this->File.read_array(fname,dummy);

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<N;i++)
          dummy3[i]=Dcurl*dummy[i];

#ifdef SMOOTHCURL
        convcomp(L1,L2,L3,this->params._d1(),this->params._d2(),this->params._d3(),N1,N2,N3,dummy3,dummy3,ftype,rsmolC,this->params._Output_directory());
#endif

#pragma omp parallel for
        for(ULONG i=0;i<N;i++)
          dummy2[i]+=dummy3[i];
      }
#endif

#ifndef _DISPLACEMENTS_
      this->displace_part_comp(out,DISPLACEMENT3,periodic,comp);


#endif



#ifdef _DISPLACEMENTS_
      this->Displacement=DISPLACEMENT3;
#endif
    }


#ifdef _WRITE_DISPLACEMENTS_
      string disp_file=this->params._Output_directory()+"Displacement_comp_"+to_string(comp)+"_"+buffsf;
      this->File.write_array(disp_file, DISPLACEMENT3);
#endif


  dummy.clear();dummy.shrink_to_fit();
  dummy2.clear();dummy2.shrink_to_fit();

}


// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================






// makecat cleaned!!!!!!!
#ifdef OMPPARRANRSD
void LPT::makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand, int ir)
{
#else
  void LPT::makecat(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand, int ir)
  {
#endif


    this->So.message_screen("Assigning position to tracers using DM particles:");
    cout<<endl;

    real_prec redshift=this->s_cosmo_pars.cosmological_redshift;

    ULONG N=this->NGRID;

  real_prec min1=this->params._xmin();
  real_prec min2=this->params._ymin();
  real_prec min3=this->params._zmin();

  real_prec ascale = 1./(1+redshift);
  real_prec Omega_L=this->s_cosmo_pars.Om_vac;
  real_prec Omega_M=this->s_cosmo_pars.Om_matter;
  real_prec hconst=this->s_cosmo_pars.hubble;
  real_prec Omega_c=this->s_cosmo_pars.Om_k;
  real_prec cvel=num_1/(cgs_km/cgs_Mpc);

  real_prec H0=static_cast<real_prec>(100.*hconst *cgs_km/cgs_Mpc/cgs_sec);
  real_prec Hub=H0*sqrt(Omega_M/ascale/ascale/ascale+Omega_L+Omega_c/ascale/ascale);
  real_prec Omega= this->s_cosmo_info.omega_matter;
  real_prec cpecvel=this->s_cosmo_info.growth_index*this->s_cosmo_pars.Hubble/(1.+redshift);

  int Nchunk = 1;
  ULONG NLLchunk=NLL/Nchunk;

  ULONG NNCH=1;

  string fname;

  string fnameNUMMRAN=string("auxnummran");//+stradd;

  string fnamePOSXORD=this->params._Output_directory()+string("aux1");//+stradd;
  string fnamePOSYORD=this->params._Output_directory()+string("aux2");//+stradd;
  string fnamePOSZORD=this->params._Output_directory()+string("aux3");//+stradd;

  vector<real_prec> dummy(N,0);
  vector<real_prec> dummy2(N,0);


  ULONG NOBJt=0;
  //4 perform RSD
  {

    {
      vector<real_prec> posx(NLL,0),posy(NLL,0),posz(NLL,0);
      {
        fname=this->fnamePOSX;
        this->File.read_array(fname+".dat",posx);
      }
      {
        fname=this->fnamePOSY;
        this->File.read_array(fname+".dat",posy);
      }
      {
        fname=this->fnamePOSZ;
        this->File.read_array(fname+".dat",posz);
      }

      //      getDensity_NGP(this->N1,this->N2,this->N3,this->L1,this->L2,this->L3,this->params._d1(),this->params._d2(),this->params._d3(),min1,min2,min3,posx,posy,posz,dummy2,NLL,dummy2,false);

      int bmax=100;
      char buffc[bmax];
      string buffchunk;

      vector<real_prec> dummychunk(NLLchunk*NNCH,0);

      //      fname=string("NG")+stradd;

#define MOCK_DEN_FIELD dummy2
      this->File.read_array(fnameMOCK+".dat",MOCK_DEN_FIELD);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NLLchunk*NNCH;i++)
        dummychunk[i]=0.0;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NLL;i++)
        dummy[i]=0.0;

      for(ULONG mm=0;mm<Nchunk;mm++)
        {
          sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
          buffchunk=static_cast<string>(buffc);

          for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
            {
              real_prec x=posx[ll];
              real_prec y=posy[ll];
              real_prec z=posz[ll];

              ULONG ix = static_cast<ULONG>(floor((x)/this->params._d1()));
              ULONG iy = static_cast<ULONG>(floor((y)/this->params._d2()));
              ULONG iz = static_cast<ULONG>(floor((z)/this->params._d3()));

              ix = static_cast<ULONG>(fmod(real_prec(ix),real_prec(N1)));
              iy = static_cast<ULONG>(fmod(real_prec(iy),real_prec(N2)));
              iz = static_cast<ULONG>(fmod(real_prec(iz),real_prec(N3)));

              if(ix >= this->N1)
                ix -= this->N1;
              if(iy >= this->N2)
                iy -= this->N2;
              if(iz >= this->N3)
                iz -= this->N3;

              if(ix < 0)
                ix += this->N1;
              if(iy < 0)
                iy += this->N2;
              if(iz < 0)
                iz += this->N3;

              ULONG lc=index_3d(ix,iy,iz,this->N2, this->N3);//iz+N3*(iy+N2*ix);

              ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
              ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));

              ULONG ichunklcc=static_cast<ULONG>(floor(lc/NLLchunk));
              ULONG llc=static_cast<ULONG>(floor(lc-ichunklcc*NLLchunk));

              if (static_cast<int>(floor(MOCK_DEN_FIELD[lc]))>0 && static_cast<int>(floor(dummy[lc])) < static_cast<int>(floor(MOCK_DEN_FIELD[lc])))
                {
                  dummychunk[lolc]=x;
                  dummy[lc]+=num_1;
                }
              else
                dummychunk[lolc]=-num_1;
            }
          {

            fname=fnamePOSXORD+static_cast<string>(buffchunk);
            this->File.write_array(fname,dummychunk);
          }
        }


  /* order simulated data in grid  */
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NLLchunk*NNCH;i++)
    dummychunk[i]=0.0;

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<NLL;i++)
    dummy[i]=0.0;

  for(ULONG mm=0;mm<Nchunk;mm++)
    {
      sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
      buffchunk=static_cast<string>(buffc);

      for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
        {
          real_prec x=posx[ll];
          real_prec y=posy[ll];
          real_prec z=posz[ll];

          ULONG ix = static_cast<ULONG>(floor((x)/this->params._d1()));
          ULONG iy = static_cast<ULONG>(floor((y)/this->params._d2()));
          ULONG iz = static_cast<ULONG>(floor((z)/this->params._d3()));

          ix = static_cast<ULONG>(fmod(real_prec(ix),real_prec(N1)));
          iy = static_cast<ULONG>(fmod(real_prec(iy),real_prec(N2)));
          iz = static_cast<ULONG>(fmod(real_prec(iz),real_prec(N3)));

          if(ix >= N1)
            ix -= N1;
          if(iy >= N2)
            iy -= N2;
          if(iz >= N3)
            iz -= N3;

          if(ix < 0)
            ix += N1;
          if(iy < 0)
            iy += N2;
          if(iz < 0)
            iz += N3;

          ULONG lc=iz+N3*(iy+N2*ix);

          ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
          ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));

          ULONG ichunklcc=static_cast<ULONG>(floor(lc/NLLchunk));
          ULONG llc=static_cast<ULONG>(floor(lc-ichunklcc*NLLchunk));

          if (static_cast<int>(floor(MOCK_DEN_FIELD[lc]))>0 && static_cast<int>(floor(dummy[lc])) < static_cast<int>(floor(MOCK_DEN_FIELD[lc])))
            {
              dummychunk[lolc]=y;
              dummy[lc]+=num_1;
            }
          else
            dummychunk[lolc]=-num_1;
        }
      {
        fname=fnamePOSYORD+static_cast<string>(buffchunk);
        this->File.write_array(fname,dummychunk);
      }
    }

  /* order simulated data in grid  */
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NLLchunk*NNCH;i++)
    dummychunk[i]=0.0;


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NLL;i++)
    dummy[i]=0.0;

  for(ULONG mm=0;mm<Nchunk;mm++)
    {
      sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
      buffchunk=static_cast<string>(buffc);

      for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
        {

          real_prec x=posx[ll];
          real_prec y=posy[ll];
          real_prec z=posz[ll];

          ULONG ix = static_cast<ULONG>(floor((x)/this->params._d1()));
          ULONG iy = static_cast<ULONG>(floor((y)/this->params._d2()));
          ULONG iz = static_cast<ULONG>(floor((z)/this->params._d3()));

          ix = static_cast<ULONG>(fmod(real_prec(ix),real_prec(N1)));
          iy = static_cast<ULONG>(fmod(real_prec(iy),real_prec(N2)));
          iz = static_cast<ULONG>(fmod(real_prec(iz),real_prec(N3)));

          if(ix >= N1)
            ix -= N1;
          if(iy >= N2)
            iy -= N2;
          if(iz >= N3)
            iz -= N3;

          if(ix < 0)
            ix += N1;
          if(iy < 0)
            iy += N2;
          if(iz < 0)
            iz += N3;

          ULONG lc=iz+N3*(iy+N2*ix);

          ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
          ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));

          ULONG ichunklcc=static_cast<ULONG>(floor(lc/NLLchunk));
          ULONG llc=static_cast<ULONG>(floor(lc-ichunklcc*NLLchunk));

#ifdef SAMPRAN
          if (static_cast<int>(floor(MOCK_DEN_FIELD[lc]))>0)
            {
              dummychunk[lolc]=z;
              dummy2[lc]-=num_1;
            }
          else
            dummychunk[lolc]=-num_1;
#else
          if (static_cast<int>(floor(MOCK_DEN_FIELD[lc]))>0 && static_cast<int>(floor(dummy[lc])) < static_cast<int>(floor(MOCK_DEN_FIELD[lc])))
            {
              dummychunk[lolc]=z;
              dummy[lc]+=num_1;
            }
          else
            dummychunk[lolc]=-num_1;

#endif

        }

      {
         fname=fnamePOSZORD+static_cast<string>(buffchunk);
         this->File.write_array(fname,dummychunk);
      }
    }

#ifdef SAVEMEM
  {
    string outputFileName0=this->params._Output_directory()+string("RANHALOS")+stradd+string(".txt");
    ofstream outStream0;
    outStream0.open(outputFileName0.data());
    assert(outStream0.is_open());
    this->So.message_screen("Writting in file ", outputFileName0);

    for(ULONG i=0;i<this->N1;i++)
      for(ULONG j=0;j<this->N2;j++)
        for(ULONG k=0;k<this->N3;k++)
          {
            ULONG ind = index_3d(i,j,k,N2,N3);
            if (MOCK_DEN_FIELD[ind]>0.)
              outStream0<<MOCK_DEN_FIELD[ind]<<" "<<ind<<endl;
          }
    outStream0.close();
  }
#endif

  //#ifdef SAMPRAN
  // ATTENTION: to be computed from the right class. Commented by ABA. fnameNUMRAN not defined
  //this->File.write_array(fnameNUMRAN,dummy2);
  //#endif

  }

    {
      fname=this->fnameDM+".dat";
      this->File.read_array(fname,dummy);
      get_overdens(dummy,dummy);
    }

    string outputFileName=this->params._Output_directory()+"CAT_realization"+to_string(ir)+"_"+this->stradd+string(".txt");
    this->fnameTRACERCAT=outputFileName;
    ofstream outStream;
    outStream.open(outputFileName.data());
    assert(outStream.is_open());
    this->So.message_screen("Ready to write in file ", outputFileName);


    int bmax_res=100;
    char buffer_res[bmax_res];

    {
      vector<real_prec> dummychunk1(NLLchunk*NNCH,0), dummychunk2(NLLchunk*NNCH,0), dummychunk3(NLLchunk*NNCH,0);

      ULONG ichunkold=1;

      int bmax=100;
      char buffc[bmax];
      string buffchunk;

      ULONG keyold=0;

#ifdef OMPPARRANRSD
      int nt=omp_get_max_threads();
  omp_set_num_threads(nt);

  {
    int jthread = omp_get_thread_num();

#pragma omp for
    for(ULONG mm=0;mm<Nchunk;mm++)
      {
#else
       for(ULONG mm=0;mm<Nchunk;mm++)
          {
#endif
            sprintf(buffer_res,"chunk%d",static_cast<int>(mm+1));

            {
              string fname_res1=fnamePOSXORD;
              string FileName_res1=string(fname_res1)+static_cast<string>(buffer_res)+string(".dat");
              this->File.read_array(FileName_res1,dummychunk1);
            }
            {
              string fname_res2=fnamePOSYORD;
              string FileName_res2=string(fname_res2)+static_cast<string>(buffer_res)+string(".dat");
              this->File.read_array(FileName_res2,dummychunk2);
            }
            {
              string fname_res3=fnamePOSZORD;
              string FileName_res3=string(fname_res3)+static_cast<string>(buffer_res)+string(".dat");
              this->File.read_array(FileName_res3,dummychunk3);
            }

        for(ULONG jj=NLLchunk*(mm);jj<NLLchunk*(1+mm);jj++)
          {
            ULONG jchunk=static_cast<ULONG>(floor(jj/NLLchunk));
            ULONG jjchunk=static_cast<ULONG>(floor(jj-jchunk*NLLchunk));

            if (dummychunk1[jjchunk]>=0. && dummychunk2[jjchunk]>=0. && dummychunk3[jjchunk]>=0. )
              {
                real_prec posxoj=dummychunk1[jjchunk];
                real_prec posyoj=dummychunk2[jjchunk];
                real_prec poszoj=dummychunk3[jjchunk];


                real_prec posxzj=posxoj;
                real_prec posyzj=posyoj;
                real_prec poszzj=poszoj;

                bool periodic=true;
                if (periodic==true)
                  {
                    real_prec xnew=posxzj;
                    real_prec ynew=posyzj;
                    real_prec znew=poszzj;

                    if (xnew<0.)
                      xnew+=L1;
                    if (xnew>=L1)
                      xnew-=L1;

                    if (ynew<0.)
                      ynew+=L1;
                    if (ynew>=L1)
                      ynew-=L1;

                    if (znew<0.)
                      znew+=L1;
                    if (znew>=L1)
                      znew-=L1;

                    posxzj=xnew;
                    posyzj=ynew;
                    poszzj=znew;
                  }

                real_prec dums=0.0;
                outStream << posxoj <<" "<< posyoj <<" "<< poszoj <<" "<< dums  <<" "<< dums <<" "<< dums <<" "<< dums <<endl;

                NOBJt++;
              }
          }
      }
    }
#ifdef OMPPARRANRSD
  }
#endif

#ifdef SAMPRAN
   {
#ifndef SAVEMEM

    fname=fnameNUMMRAN;
    this->File.read_array(fname,dummy2);

#else

    string outputFileName0=this->params._Output_directory()+string("RANHALOS")+stradd+string(".txt");

    ifstream inStream0;
    inStream0.open(outputFileName0.data());
    assert(inStream0.is_open());

    ULONG NCO=this->File.count_lines(outputFileName0);

    vector<real_prec> dumRH(NCO,0);

    vector<ULONG> dumRHind(NCO,0);

    for(ULONG i=0;i<NCO;i++)
      inStream0>>dumRH[i]>>dumRHind[i];

    //inStream0>>dumRH[i]>>dumi[i]>>dumj[i]>>dumk[i];commented by FSK

    inStream0.close();

    ULONG count=0;
#endif


#ifdef OMPPARRANRSD
  int nt=omp_get_max_threads();
  omp_set_num_threads(nt);

  {
   int jthread = omp_get_thread_num();

#pragma omp for
    for (ULONG i=0;i<N1;i++)
      for (ULONG j=0;j<N2;j++)
        for (ULONG k=0;k<N3;k++)
          {
#else
    for (ULONG i=0;i<N1;i++)
      for (ULONG j=0;j<N2;j++)
        for (ULONG k=0;k<N3;k++)
          {
#endif
            //            ULONG lc=k+N3*(j+N2*i);
            ULONG lc= index_3d(i,j,k,N2,N3);
#ifdef SAVEMEM
            if (lc==dumRHind[count])
            //if (i==static_cast<ULONG>(dumi[count]) && j==static_cast<ULONG>(dumj[count]) && k==static_cast<ULONG>(dumk[count]) )
            {
            for(ULONG kk=0;kk<static_cast<ULONG>(floor(dumRH[count]));kk++)
              {
#else
            for(ULONG kk=0;kk<static_cast<ULONG>(floor(MOCK_DEN_FIELD[lc]));kk++)
              {
#endif

#ifdef OMPPARRANRSD
                real_prec rx=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]))*this->params._d1();
                real_prec ry=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]))*this->params._d2();
                real_prec rz=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]))*this->params._d3();
#else
                real_prec rx=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*this->params._d1();
                real_prec ry=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*this->params._d2();
                real_prec rz=static_cast<real_prec>(gsl_rng_uniform(gBaseRand))*this->params._d3();
#endif


                real_prec posxoj=this->params._d1()*(static_cast<real_prec>(i))+rx;
                real_prec posyoj=this->params._d2()*(static_cast<real_prec>(j))+ry;
                real_prec poszoj=this->params._d3()*(static_cast<real_prec>(k))+rz;

                real_prec posxzj=posxoj;
                real_prec posyzj=posyoj;
                real_prec poszzj=poszoj;

                bool periodic=true;
                if (periodic==true)
                  {
                    real_prec xnew=posxzj;
                    real_prec ynew=posyzj;
                    real_prec znew=poszzj;

                    if (xnew<0.)
                      xnew+=L1;
                    if (xnew>=L1)
                      xnew-=L1;

                    if (ynew<0.)
                      ynew+=L1;
                    if (ynew>=L1)
                      ynew-=L1;

                    if (znew<0.)
                      znew+=L1;
                    if (znew>=L1)
                      znew-=L1;

                    posxzj=xnew;
                    posyzj=ynew;
                    poszzj=znew;
                  }

                real_prec dums=0.0;
                outStream << posxoj <<"\t"<< posyoj <<"\t"<< poszoj <<"\t"<< dums <<"\t"<< dums  <<"\t"<< dums <<"\t"<< dums <<endl;

                NOBJt++;
              }
#ifdef SAVEMEM
            count++;
              }
#endif
          }
    }
#ifdef OMPPARRANRSD
  }
#endif

#endif

   outStream.close();
  }
#ifdef _FULL_VERBOSE_
  this->So.message_screen("Number of tracers = ",NOBJt);
#endif
}



// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================




void LPT::theta2velcomp(vector<real_prec> &Theta, vector<real_prec> &vei, bool zeropad, bool norm, int comp)
 {
    // theta = nabla V. Getting V given theta is the main goal.
    // For norm=false, this function solves for the Displacement(vei) from the Divergencve of the Displacement (nabla delta) in order to obtain the component comp of the displacement.
    // For norm=true, this function solves for the Displacement from the Divergencve of the Displacement and multiply for
    // the factor fHa in order to get the component comp of te velocity filed.
    //

   ULONG N=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3);

   real_prec Lzp1=this->L1;
   real_prec Lzp2=this->L2;
   real_prec Lzp3=this->L3;

   ULONG Nzp1=this->N1;
   ULONG Nzp2=this->N2;
   ULONG Nzp3=this->N3;

   if (zeropad==true)
     {
       Lzp1=L1*num_2;
       Lzp2=L2*num_2;
       Lzp3=L3*num_2;

       Nzp1=this->N1*2;
       Nzp2=this->N2*2;
       Nzp3=this->N3*2;
     }


   real_prec cpecvel=1.;
   if(true==norm)
     cpecvel=this->s_cosmo_info.growth_index*this->s_cosmo_info.Hubble_parameter/(1.+this->s_cosmo_pars.cosmological_redshift);
   // f(z)*H(z)/(1+z). The scale factor is used such that we deal with proper velocities

   ULONG Nhalf=static_cast<ULONG>(this->N1)*static_cast<ULONG>(this->N2)*static_cast<ULONG>(this->N3/2+1);
#ifdef DOUBLE_PREC
   complex_prec *Theta_f= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
   complex_prec *veli= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
   complex_prec *Theta_f= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
   complex_prec *veli= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif

   // Fourier transform of the nabla X
   do_fftw_r2c(N1,Theta,Theta_f);

  vector<real_prec> coords(N1,0);
  real_prec deltak=2.*M_PI/this->params._Lbox();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N1 ;++i)
    coords[i]=deltak*(i<=N1/2? static_cast<real_prec>(i): -static_cast<real_prec>(N1-i));

  // Get X


/*real_prec kstar2=0;
  // This par tis commeted, for the convolution is to be made at Bam.cpp
if(true==this->params._use_vel_kernel())
  kstar2=pow(this->params._slengthv(),2.);  
*/

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for (ULONG i=0;i<N1;i++)
     for (ULONG j=0;j<N2;j++)
       for (ULONG k=0;k<N3/2+1;k++)
         {
           real_prec kx=calc_kx(i,Lzp1,Nzp1);
           real_prec ky=calc_ky(j,Lzp2,Nzp2);
           real_prec kz=calc_kz(k,Lzp3,Nzp3);
           real_prec kmod2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);  // Get k**2
           ULONG index=index_3d(i,j,k,this->N2,this->N3/2+1);
           real_prec vel_kernel = 1.0;//1./(1.+kmod2*kstar2);
           veli[index][REAL]=cpecvel*linearvel3d(comp, kx, ky, kz, -Theta_f[index][IMAG]*vel_kernel);
           veli[index][IMAG]=cpecvel*linearvel3d(comp, kx, ky, kz, Theta_f[index][REAL]*vel_kernel);
         }

   /*

     switch(comp)
     {
     case(1):
       {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for (ULONG i=0;i<N1;i++)
           for (ULONG j=0;j<N2;j++)
             for (ULONG k=0;k<N3/2+1;k++)
               {
                 real_prec kx=calc_kx(i,Lzp1,Nzp1);
                 real_prec ky=calc_ky(j,Lzp2,Nzp2);
                 real_prec kz=calc_kz(k,Lzp3,Nzp3);
                 ULONG index=index_3d(i,j,k,this->N2,this->N3/2+1);
                 veli[index][REAL]=cpecvel*linearvel3d(1, kx, ky, kz, -phi[index][IMAG]);
                 veli[index][IMAG]=cpecvel*linearvel3d(1, kx, ky, kz, phi[index][REAL]);
               }
       }
       break;
     case(2):
       {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for (ULONG i=0;i<N1;i++)
           for (ULONG j=0;j<N2;j++)
             for (ULONG k=0;k<N3/2+1;k++)
               {
                 real_prec kx=calc_kx(i,Lzp1,Nzp1);
                 real_prec ky=calc_ky(j,Lzp2,Nzp2);
                 real_prec kz=calc_kz(k,Lzp3,Nzp3);
                 ULONG index=index_3d(i,j,k,this->N2,this->N3/2+1);
                 veli[index][REAL]=cpecvel*linearvel3d(2, kx, ky, kz, -phi[index][IMAG]);
                 veli[index][IMAG]=cpecvel*linearvel3d(2, kx, ky, kz, phi[index][REAL]);
               }
       }
       break;
     case(3):
       {

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
         for (ULONG i=0;i<N1;i++)
           for (ULONG j=0;j<N2;j++)
             for (ULONG k=0;k<N3/2+1;k++)
               {
                 real_prec kx=calc_kx(i,Lzp1,Nzp1);
                 real_prec ky=calc_ky(j,Lzp2,Nzp2);
                 real_prec kz=calc_kz(k,Lzp3,Nzp3);
                 ULONG index=index_3d(i,j,k,this->N2,this->N3/2+1);
                 veli[index][REAL]=cpecvel*linearvel3d(3, kx, ky, kz, -phi[index][IMAG]);
                 veli[index][IMAG]=cpecvel*linearvel3d(3, kx, ky, kz, phi[index][REAL]);
               }
       }
     }

   */



   do_fftw_c2r(this->N1,veli,vei);
#ifdef DOUBLE_PREC
   fftw_free(veli);
   fftw_free(Theta_f);
#else
   fftwf_free(veli);
   fftwf_free(Theta_f);
#endif
 }





// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// THIS FUNCTION IS ORIGINALLY IN calc_power and is meant to interpolate an input
// power spectrum into a grid. We can put this into the Ffftwfunctions class
// where we use to do that as well.

 void LPT::read_tabulated_power()
{

  ULONG N=this->NGRID;
  int NTHREADS = omp_get_max_threads();
  vector<real_prec> Power3D(N,0);

  //read file input Pk

  real_prec Vol=L1*L2*L3;
  real_prec CONT_NORM=Vol;
  real_prec Norm=1.0;


#ifdef	  FOURIER_DEF_1
  Norm=static_cast<real_prec>(1./CONT_NORM);
#endif

#ifdef	  FOURIER_DEF_2
  Norm=static_cast<real_prec>(static_cast<real_prec>(N)*static_cast<real_prec>(N)/CONT_NORM);
#endif

  vector<real_prec>prop;
  ULONG nbinPS=this->File.read_file(this->params._dir()+this->params._ic_power_file(), prop, NTHREADS);
  int NCOLS=(static_cast<int>(prop.size()/nbinPS));
  vector<gsl_real> kPS(nbinPS,0);
  vector<gsl_real> powPS(nbinPS,0);



#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<nbinPS;i++)
    {
      kPS[i]=static_cast<real_prec>(prop[i*NCOLS]);
      powPS[i]=static_cast<real_prec>(prop[1+i*NCOLS]) / (this->growth_ini*this->growth_ini);
    }
  prop.clear();
  prop.shrink_to_fit();


  this->So.message_screen("Interpolating initial P(k) on the Fourier grid...");
  this->So.message_screen("Growth factor at initial redshift in input P(K)=", this->growth_ini);


  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, kPS.size());
  gsl_spline_init (spline, &kPS[0], &powPS[0], kPS.size());

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<this->N1;i++)
    for(ULONG j=0;j<this->N2;j++)
      for(ULONG k=0;k<this->N3;k++)
        {
          real_prec k2=k_squared(i,j,k,this->L1,this->L2,this->L3,this->N1,this->N2,this->N3);
          real_prec ktot=sqrt(k2);
          ULONG iid=index_3d(i,j,k,N2,N3);
          if(iid==0)
            Power3D[iid]=0.;
          else
            Power3D[iid]=static_cast<real_prec>(gsl_spline_eval (spline, ktot, acc)*Norm);
        }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  this->So.DONE();
  this->File.write_array(this->fname3DPOWER,Power3D);

}

 // ======================================================================================================

 real_prec LPT::linInterp(real_prec xpos, real_prec ypos, real_prec zpos, const vector<real_prec>&tab)
 {

#ifndef CELLBOUND
   xpos-=num_0_5*this->params._d1();
   ypos-=num_0_5*this->params._d2();
   zpos-=num_0_5*this->params._d3();

   {
     real_prec xnew=xpos;
     real_prec ynew=ypos;
     real_prec znew=zpos;

     if (xnew<0.)
       xnew+=L1;
     if (xnew>=L1)
       xnew-=L1;

     if (ynew<0.)
       ynew+=L1;
     if (ynew>=L1)
       ynew-=L1;

     if (znew<0.)
       znew+=L1;
     if (znew>=L1)
       znew-=L1;
     xpos=xnew;
     ypos=ynew;
     zpos=znew;
   }
#endif

  ULONG ix = static_cast<ULONG>(floor(xpos/this->params._d1()));
  ULONG iy = static_cast<ULONG>(floor(ypos/this->params._d2()));
  ULONG iz = static_cast<ULONG>(floor(zpos/this->params._d3()));


  real_prec tx = (xpos-static_cast<real_prec>(ix)*this->params._d1())/this->params._d1();
  real_prec ty = (ypos-static_cast<real_prec>(iy)*this->params._d2())/this->params._d2();
  real_prec tz = (zpos-static_cast<real_prec>(iz)*this->params._d3())/this->params._d3();


  ULONG shiftx=1;
  ULONG shifty=1;
  ULONG shiftz=1;

  long ixs=ix+shiftx;
  long iys=iy+shifty;
  long izs=iz+shiftz;

  if (ixs>=this->N1)
    ixs-=this->N1;
  if (iys>=this->N2)
    iys-=this->N2;
  if (izs>=this->N3)
    izs-=this->N3;

  real_prec y1 = tab[iz+this->N3*(iy+N2*ix)];
  real_prec y2 = tab[izs+this->N3*(iy+N2*ix)];
  real_prec y3 = tab[izs+this->N3*(iys+N2*ix)];
  real_prec y4 = tab[izs+this->N3*(iys+N2*ixs)];
  real_prec y5 = tab[iz+this->N3*(iys+N2*ix)];
  real_prec y6 = tab[iz+this->N3*((iys)+N2*ixs)];
  real_prec y7 = tab[iz+this->N3*(iy+N2*ixs)];
  real_prec y8 = tab[izs+this->N3*(iy+N2*ixs)];

  real_prec  out = static_cast<real_prec>((1.-tx)*(1.-ty)*(1.-tz)*y1 + tx*(1.-ty)*(1.-tz)*y2 + tx*ty*(1.-tz)*y3 + tx*ty*tz*y4 + (1.-tx)*ty*(1.-tz)*y5 + (1.-tx)*ty*tz*y6 + (1.-tx)*(1.-ty)*tz*y7 + tx*(1.-ty)*tz*y8);

  return out;
}




// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================


#ifdef _USE_OMP_
 void LPT::makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand,int ir)
#else
 void LPT::makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand,int ir)
#endif
 {
   // #else
   //  void LPT::makecat_withv(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng *  gBaseRand,int ir){
   // #endif


     //NLL is passed as NGRID from BAM


   this->So.message_screen("Assigning position and velocities to tracers using DM particles:");
   cout<<endl;

   real_prec redshift=this->s_cosmo_pars.cosmological_redshift;


   s_params_box_mas box;
   box.min1=this->params._xmin();
   box.min2=this->params._ymin();
   box.min3=this->params._zmin();
   box.Lbox=this->params._Lbox();
   box.Nft=this->N1;
   box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
   box.d2= box.d1;
   box.d3= box.d1;
   box.NGRID=this->NGRID;



   ULONG N=this->NGRID;

   real_prec min1=this->params._xmin();
   real_prec min2=this->params._ymin();
   real_prec min3=this->params._zmin();


   real_prec ascale = 1./(1+redshift);
   real_prec Omega_L=this->s_cosmo_pars.Om_vac;
   real_prec Omega_M=this->s_cosmo_pars.Om_matter;
   real_prec hconst=this->s_cosmo_pars.hubble;
   real_prec Omega_c=this->s_cosmo_pars.Om_k;
   real_prec cvel=num_1/(cgs_km/cgs_Mpc);

   real_prec H0=static_cast<real_prec>(100.*hconst *cgs_km/cgs_Mpc/cgs_sec);
   real_prec Hub=H0*sqrt(Omega_M/ascale/ascale/ascale+Omega_L+Omega_c/ascale/ascale);
   real_prec Omega= this->s_cosmo_info.omega_matter;
   real_prec cpecvel=this->s_cosmo_info.growth_index*this->s_cosmo_pars.Hubble/(1.+this->s_cosmo_pars.cosmological_redshift);


   ULONG Nchunk=1;   //original ok
   ULONG NLLchunk=this->NGRID/Nchunk;

   ULONG NNCH=1;


   string fname;

   string fnameNUMMRAN=string("auxnummran");//+stradd;

   string fnamePOSXORD=this->params._Output_directory()+string("aux1");//+stradd;
   string fnamePOSYORD=this->params._Output_directory()+string("aux2");//+stradd;
   string fnamePOSZORD=this->params._Output_directory()+string("aux3");//+stradd;

   vector<real_prec> dummy(this->NGRID,0);
   vector<real_prec> dummy2(this->NGRID,0);

   real_prec vsmoo=this->params._vslength();

   ULONG NOBJt=0;


   vector<real_prec> posx(this->NGRID,0),posy(this->NGRID,0),posz(this->NGRID,0);

   // Filenames of the bninary files containing the positions of the dark matter particles
   this->File.read_array(this->fnamePOSX+".dat",posx);
   this->File.read_array(this->fnamePOSY+".dat",posy);
   this->File.read_array(this->fnamePOSZ+".dat",posz);

   int bmax=100;
   char buffc[bmax];
   string buffchunk;

   vector<real_prec> dummychunk(NLLchunk*NNCH,0);


  //Read the Tracer number counts
#define MOCK_DEN_FIELD dummy2
   this->File.read_array(fnameMOCK+".dat",MOCK_DEN_FIELD);

   // abatest ****************
   vector<real_prec>NUMBER_COUNTS(MOCK_DEN_FIELD.size(),0);
   NUMBER_COUNTS=MOCK_DEN_FIELD;
   ULONG nob=get_nobjects(NUMBER_COUNTS);

   // **************************

#pragma omp parallel for
   for(ULONG i=0;i<NLLchunk*NNCH;i++)
     dummychunk[i]=0.0;

#pragma omp parallel for
   for(ULONG i=0;i<this->NGRID;i++)
     dummy[i]=0.0;

   for(ULONG mm=0;mm<Nchunk;mm++)
     {
       sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
       buffchunk=static_cast<string>(buffc);

       for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
         {
           real_prec x=posx[ll]; //positions of the dm particles
           real_prec y=posy[ll];
           real_prec z=posz[ll];
           ULONG lc =grid_ID(&box, x,y,z); //index_3d(ix,iy,iz,N2,N3);
           ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
           ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));

//	   ULONG ichunklcc=static_cast<ULONG>(floor(lc/NLLchunk));
//	   ULONG llc=static_cast<ULONG>(floor(lc-ichunklcc*NLLchunk));

           if (static_cast<int>(MOCK_DEN_FIELD[lc])>0 && ( static_cast<int>(floor(dummy[lc])) < static_cast<int>(MOCK_DEN_FIELD[lc]) ))
             {
               dummychunk[lolc]=x;     //save the x-component and write to binary
               dummy[lc]++;
             }
           else
             dummychunk[lolc]=-num_1;  //else mark them with -1
         }

       fname=fnamePOSXORD+static_cast<string>(buffchunk);
       this->File.write_array(fname,dummychunk);
     }

   /* order simulated data in grid  */

#pragma omp parallel for
   for(ULONG i=0;i<NLLchunk*NNCH;i++)
     dummychunk[i]=0.0;

#pragma omp parallel for
   for(ULONG i=0;i<this->NGRID;i++)
     dummy[i]=0.0;

   for(ULONG mm=0;mm<Nchunk;mm++)
     {
       sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
       buffchunk=static_cast<string>(buffc);

       for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
         {
           real_prec x=posx[ll];
           real_prec y=posy[ll];
           real_prec z=posz[ll];
           ULONG lc=grid_ID(&box, x,y,z);
           ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
           ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));
           if (static_cast<int>(MOCK_DEN_FIELD[lc])>0 && static_cast<int>(floor(dummy[lc])) < static_cast<int>(floor(MOCK_DEN_FIELD[lc])))
             {
               dummychunk[lolc]=y;   //save the y-component and write to binary
               dummy[lc]++;
             }
           else
             dummychunk[lolc]=-num_1;
         }

       fname=fnamePOSYORD+static_cast<string>(buffchunk);
       this->File.write_array(fname,dummychunk);
     }

   /* order simulated data in grid  */

#pragma omp parallel for
   for(ULONG i=0;i<NLLchunk*NNCH;i++)
     dummychunk[i]=0.0;


#pragma omp parallel for
   for(ULONG i=0;i<NLL;i++)
     dummy[i]=0.0;

   for(ULONG mm=0;mm<Nchunk;mm++)
     {
       sprintf(buffc,"chunk%d",static_cast<int>(mm+1));
       buffchunk=static_cast<string>(buffc);

       for(ULONG ll=NLLchunk*(mm);ll<NLLchunk*(1+mm);ll++)
         {
           real_prec x=posx[ll];
           real_prec y=posy[ll];
           real_prec z=posz[ll];
           ULONG lc=grid_ID(&box, x,y,z);
           ULONG ichunklc=static_cast<ULONG>(floor(ll/NLLchunk));
           ULONG lolc=static_cast<ULONG>(floor(ll-ichunklc*NLLchunk));
           if (static_cast<int>(MOCK_DEN_FIELD[lc])>0) // If we have tracers in that cell:
             {
               dummychunk[lolc]=z;	      // save the z-coordinate component to be written to binary
               MOCK_DEN_FIELD[lc]--;     // and subtract one particle from that cell. Note that this subtraction is only done here, for this is the last place calling the MOCK_NUM_DEN
             }
           else
             dummychunk[lolc]=-num_1;         // otherwise assign a -1
         }
       fname=fnamePOSZORD+static_cast<string>(buffchunk);
       this->File.write_array(fname,dummychunk);
     }



   string outputFileName0=this->params._Output_directory()+string("RANHALOS")+stradd+string(".txt");
   ofstream outStream0;
   outStream0.open(outputFileName0.data());
   assert(outStream0.is_open());

   // In the previous part the MOCK_DEN_FIELD container was reduced, leaving the number of random objects required per cell:
   // This is written now here:
   ULONG cc=0;
   for(ULONG ind=0;ind<this->NGRID;++ind)
     if(MOCK_DEN_FIELD[ind]>0.)
       {
         cc++;
         outStream0<<MOCK_DEN_FIELD[ind]<<" "<<ind<<endl;
         }
   outStream0.close();
   So.message_screen("fraction of random objects = ", 100.0*cc/static_cast<real_prec>(nob));

   // Open file to print catalog
   string outputFileName=this->params._Output_directory()+"CAT_realization"+to_string(ir)+"_"+this->stradd+string(".txt");
   this->fnameTRACERCAT=outputFileName;
   ofstream outStream;
   outStream.open(outputFileName.c_str());
   assert(outStream.is_open());
   outStream.precision(8);
   outStream.setf(ios::showpoint);
   outStream.setf(ios::scientific);


   vector<real_prec>vex (NLL,0);
   vector<real_prec>vey(NLL,0);



#ifdef NORSD
#pragma omp parallel for
   for(ULONG i=0;i<NLL;i++)
     {
       vex[i]=0.0;
       vey[i]=0.0;
       dummy2[i]=0.0;
     }
#else
   this->File.read_array(this->fnameVX+".dat",vex);
   this->File.read_array(this->fnameVY+".dat",vey);
   this->File.read_array(this->fnameVZ+".dat",dummy2);

#endif

   int bmax_res=100;
   char buffer_res[bmax_res];


   vector<real_prec> dummychunk1(NLLchunk*NNCH,0), dummychunk2(NLLchunk*NNCH,0), dummychunk3(NLLchunk*NNCH,0);

   ULONG ichunkold=1;

   //***** test by aba to verify number counts******************

   int nt=omp_get_max_threads();
   omp_set_num_threads(nt);


   {

     int jthread = omp_get_thread_num();   // original ok

    //Write to cat the positions of tracers identified with dark matter particles

#pragma omp for
     for(ULONG mm=0;mm<Nchunk;mm++)
       {
         sprintf(buffer_res,"chunk%d",static_cast<int>(mm+1));

         string fname_res1=fnamePOSXORD;
         string FileName_res1=string(fname_res1)+static_cast<string>(buffer_res)+string(".dat");
         this->File.read_array(FileName_res1,dummychunk1);


         string fname_res2=fnamePOSYORD;
         string FileName_res2=string(fname_res2)+static_cast<string>(buffer_res)+string(".dat");
         this->File.read_array(FileName_res2,dummychunk2);


         string fname_res3=fnamePOSZORD;
         string FileName_res3=string(fname_res3)+static_cast<string>(buffer_res)+string(".dat");
         this->File.read_array(FileName_res3,dummychunk3);


         for(ULONG jj=NLLchunk*(mm);jj<NLLchunk*(1+mm);jj++)
           {
             ULONG jchunk=static_cast<ULONG>(floor(jj/NLLchunk));
             ULONG jjchunk=static_cast<ULONG>(floor(jj-jchunk*NLLchunk));

             if (dummychunk1[jjchunk]>=0. && dummychunk2[jjchunk]>=0. && dummychunk3[jjchunk]>=0. )
               {
                 real_prec posxoj=dummychunk1[jjchunk];
                 real_prec posyoj=dummychunk2[jjchunk];
                 real_prec poszoj=dummychunk3[jjchunk];

                 real_prec vx=this->linInterp(posxoj,posyoj,poszoj,vex);
                 real_prec vy=this->linInterp(posxoj,posyoj,poszoj,vey);
                 real_prec vz=this->linInterp(posxoj,posyoj,poszoj,dummy2);

                 real_prec vexoj=vx;
                 real_prec veyoj=vy;
                 real_prec vezoj=vz;

                 real_prec dums=0.0;
                 int index=1;
                 outStream << posxoj <<"\t"<< posyoj <<"\t"<< poszoj <<"\t"<< vexoj <<"\t"<< veyoj <<"\t"<< vezoj <<"\t"<< dums <<"\t"<<index<<endl;

                 NOBJt++;
               }
           }
       }

   }


   outputFileName0=this->params._Output_directory()+string("RANHALOS")+stradd+string(".txt");
   ifstream inStream0;
   this->So.message_screen("Writting in file", outputFileName0);
   inStream0.open(outputFileName0.data());
   assert(inStream0.is_open());
   ULONG NCO=m_countLines(inStream0);
   vector<int> dumRH(NCO, 0);
   vector<ULONG> dumRHind(NCO, 0);

   // Read the number of random objects dumRH in a cell identified with the index dumRHind
   for(ULONG i=0;i<NCO;i++)
     inStream0>>dumRH[i]>>dumRHind[i];
   inStream0.close();

  // Write in this part the position of the tracers identified with random positions

   //#ifdef OMPPARRANRSD

   nt=omp_get_max_threads();
   omp_set_num_threads(nt);

   {
     int jthread = omp_get_thread_num();
     ULONG count=0;
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
     for (ULONG i=0;i<N1;i++)
       for (ULONG j=0;j<N2;j++)
         for (ULONG k=0;k<N3;k++)
           {
             ULONG lc=index_3d(i,j,k,N2,N3);
             if (lc==dumRHind[count])
               {
                 for(ULONG kk=0;kk < dumRH[count]; ++kk) // loop over the number of random objects in cells
                   {

#ifdef _USE_OMP_
                     real_prec rx=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                     real_prec ry=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                     real_prec rz=static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
#else
                         real_prec rx=static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec ry=static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                         real_prec rz=static_cast<real_prec>(gsl_rng_uniform(gBaseRand));

#endif

                     real_prec posxoj=this->params._d1()*(static_cast<real_prec>(i)+rx);
                     real_prec posyoj=this->params._d2()*(static_cast<real_prec>(j)+ry);
                     real_prec poszoj=this->params._d3()*(static_cast<real_prec>(k)+rz);

                     // esto las pone en el centro de la celda, y reduce de 876 a 680 el numero de celdas que no coinciden en number counts
                     // dejo la de arriba pero ultiplico el random por un factor 0.95 par asegurarnos que que quede dentro
                     // real_prec posxoj=d1*(static_cast<real_prec>(i)+0.5);
                     // real_prec posyoj=d2*(static_cast<real_prec>(j)+0.5);
                     // real_prec poszoj=d3*(static_cast<real_prec>(k)+0.5);



                     real_prec vx=this->linInterp(posxoj,posyoj,poszoj,vex);
                     real_prec vy=this->linInterp(posxoj,posyoj,poszoj,vey);
                     real_prec vz=this->linInterp(posxoj,posyoj,poszoj,dummy2);

                     real_prec vexoj=vx;
                     real_prec veyoj=vy;
                     real_prec vezoj=vz;
                     real_prec dums=0.0;
                     int index=-num_1;

                     outStream << posxoj <<"\t"<< posyoj <<"\t"<< poszoj <<"\t"<< vexoj <<"\t"<< veyoj <<"\t"<< vezoj <<"\t"<< dums <<"\t"<<index<<endl;
                     NOBJt++;
                   }
                 count++;
               }
           }
     outStream.close();
#ifdef _FULL_VERBOSE_
     this->So.message_screen("Number of tracers = ",NOBJt);
#endif
   }// close omp

   this->Number_of_Tracers = NOBJt;
 }

// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================


#ifdef _USE_OMP_
 void LPT::makecat_withv_new(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng ** gBaseRand,int ir)
#else
 void LPT::makecat_withv_new(string stradd,string fnameMOCK, ULONG NLL,ULONG NOBJ,gsl_rng * gBaseRand,int ir)
#endif
 {

   this->So.message_screen("Assigning position and velocities to tracers using DM particles:");
   cout<<endl;

   real_prec redshift=this->s_cosmo_pars.cosmological_redshift;


   s_params_box_mas box;
   box.min1=this->params._xmin();
   box.min2=this->params._ymin();
   box.min3=this->params._zmin();
   box.Lbox=this->params._Lbox();
   box.Nft=this->N1;
   box.d1= box.Lbox/static_cast<real_prec>(box.Nft);		/* grid spacing x-direction */
   box.d2= box.d1;
   box.d3= box.d1;
   box.NGRID=this->NGRID;


   ULONG N=this->NGRID;

   real_prec min1=this->params._xmin();
   real_prec min2=this->params._ymin();
   real_prec min3=this->params._zmin();


   real_prec ascale = 1./(1+redshift);
   real_prec Omega_L=this->s_cosmo_pars.Om_vac;
   real_prec Omega_M=this->s_cosmo_pars.Om_matter;
   real_prec hconst=this->s_cosmo_pars.hubble;
   real_prec Omega_c=this->s_cosmo_pars.Om_k;
   real_prec cvel=num_1/(cgs_km/cgs_Mpc);

   real_prec H0=static_cast<real_prec>(100.*hconst *cgs_km/cgs_Mpc/cgs_sec);
   real_prec Hub=H0*sqrt(Omega_M/ascale/ascale/ascale+Omega_L+Omega_c/ascale/ascale);
   real_prec Omega= this->s_cosmo_info.omega_matter;
   real_prec cpecvel=this->s_cosmo_info.growth_index*this->s_cosmo_pars.Hubble/(1.+this->s_cosmo_pars.cosmological_redshift);

   ULONG Nchunk=1;   //original ok
   ULONG NLLchunk=this->NGRID/Nchunk;

   ULONG NNCH=1;
   string fname;
   string fnameNUMMRAN=string("auxnummran");//+stradd;
   string fnamePOSXORD=this->params._Output_directory()+string("aux1");//+stradd;
   string fnamePOSYORD=this->params._Output_directory()+string("aux2");//+stradd;
   string fnamePOSZORD=this->params._Output_directory()+string("aux3");//+stradd;

   // ****************************DM POSITIONS *******************************************
   ULONG N_dm=this->NGRID;
   // NOte that hwre the numberof dm particles is taht odf the grid

   vector<real_prec> posx(N_dm,0),posy(N_dm,0),posz(N_dm,0);
   vector<ULONG>index(N_dm,0);
   // Filenames of the bninary files containing the positions of the dark matter particles
   this->File.read_array(this->fnamePOSX+".dat",posx);
   this->File.read_array(this->fnamePOSY+".dat",posy);
   this->File.read_array(this->fnamePOSZ+".dat",posz);

   vector<real_prec> velx(N_dm,0),vely(N_dm,0),velz(N_dm,0);
   this->File.read_array(this->fnameVX+".dat",velx);
   this->File.read_array(this->fnameVY+".dat",vely);
   this->File.read_array(this->fnameVZ+".dat",velz);

   struct s_cell_info{
     vector<real_prec> posx_p;
     vector<real_prec> posy_p;
     vector<real_prec> posz_p;
     vector<real_prec> velx_p;
     vector<real_prec> vely_p;
     vector<real_prec> velz_p;
     vector<real_prec> mass_p;
     vector<ULONG> id_p;
   };

   vector<s_cell_info> cell_inf_dm(this->NGRID);

   vector<real_prec> MOCK_DEN_FIELD(this->NGRID,0);
   this->File.read_array(fnameMOCK+".dat",MOCK_DEN_FIELD);

   vector<int> aux_cont(this->NGRID,0);
#pragma omp parallel for
   for(ULONG i=0;i<N_dm;++i)
     {
       real_prec x = static_cast<real_prec>(posx[i]);
       real_prec y = static_cast<real_prec>(posy[i]);
       real_prec z = static_cast<real_prec>(posz[i]);
       ULONG ind=grid_ID(&box, x,y,z);
       index[i]=ind;
#pragma omp atomic update
       aux_cont[ind]++;  //count the number of dm particle in a cell
     }

   vector<bool>dm_used(this->NGRID,false);
   vector<int>dm_cases(this->NGRID,0);

   //ACA DEBEMOS CONTEMPLAR TRES CASOS:

   // I) EN UNA CELDA HAY MAS TRACERS QUE DM
   // II) EN UNA CELDA HAY DM QUE TRACERS
   // III) En una celda hay Tracers pero no hay DM
   // IV) Celdas vacias, deben permanecer vacias


   ULONG Nobjects_mock=get_nobjects(MOCK_DEN_FIELD);

   ULONG empty_cells_original=0;
#pragma omp parallel for reduction(+:empty_cells_original)
   for(ULONG id=0;id<this->NGRID;++id)
     if(MOCK_DEN_FIELD[id]==0)
       {
         empty_cells_original++;
         dm_cases[id]=4;
         dm_used[id]=false;
       }

   //   So.message_screen("Original number of empty cells ", empty_cells_original);

   So.message_screen("Retrieving positions of dm particles");

   vector<int> aux_cont1(this->NGRID,0);
   vector<int> aux_cont2(this->NGRID,0);

   // caso I: mas trcers que dm
   for(ULONG i=0;i<N_dm;++i)
    {
      real_prec x = static_cast<real_prec>(posx[i]);
      real_prec y = static_cast<real_prec>(posy[i]);
      real_prec z = static_cast<real_prec>(posz[i]);
      ULONG id   =  index[i];    //identify the cell where this particle lives
      if(MOCK_DEN_FIELD[id]>0 &&  aux_cont[id]>0)  // si hay tracers en esta celda y dm tambien
        {
          if(MOCK_DEN_FIELD[id]> aux_cont[id]) //si  hay mas o igual número tracers que dm (y hay dm), tomar todas las dm que hay en cada celda. Al haber mas tracers que dm, vendrá la necesidad de tener randoms
            {
              real_prec vx=this->linInterp(x,y,z,velx);
              real_prec vy=this->linInterp(x,y,z,vely);
              real_prec vz=this->linInterp(x,y,z,velz);
              cell_inf_dm[id].posx_p.push_back(x);
              cell_inf_dm[id].posy_p.push_back(y);
              cell_inf_dm[id].posz_p.push_back(z);
              cell_inf_dm[id].velx_p.push_back(vx);
              cell_inf_dm[id].vely_p.push_back(vy);
              cell_inf_dm[id].velz_p.push_back(vz);
              cell_inf_dm[id].id_p.push_back(id);
              dm_used[id]=true;
              dm_cases[id]=1;
              aux_cont1[id]++;
            }
          else
            { // caso II: se necesitan menor numero de tracers que el numero de dm en la celda. Ntr<=Ndm
              if(aux_cont2[id]< MOCK_DEN_FIELD[id]) // el "if" es para tomar sólo los que necesitamos
                {
                  real_prec vx=this->linInterp(x,y,z,velx);
                  real_prec vy=this->linInterp(x,y,z,vely);
                  real_prec vz=this->linInterp(x,y,z,velz);
                  cell_inf_dm[id].posx_p.push_back(x);
                  cell_inf_dm[id].posy_p.push_back(y);
                  cell_inf_dm[id].posz_p.push_back(z);
                  cell_inf_dm[id].velx_p.push_back(vx);
                  cell_inf_dm[id].vely_p.push_back(vy);
                  cell_inf_dm[id].velz_p.push_back(vz);
                  cell_inf_dm[id].id_p.push_back(id);
                  dm_used[id]=true;
                  dm_cases[id]=2;
                  aux_cont2[id]++;
                }
            }
        }
    }
   index.clear();
   index.shrink_to_fit();

   // Get the number of randoms from case i and ii
   vector<int>Nrandom_tracers(this->NGRID,0);

#pragma omp parallel for
   for(ULONG id=0;id<this->NGRID;++id)
     if(MOCK_DEN_FIELD[id]>0)
       if(1==dm_cases[id] || 2==dm_cases[id])
         Nrandom_tracers[id]=MOCK_DEN_FIELD[id]-cell_inf_dm[id].posx_p.size();


   ULONG Ndm_used_I=0;
   ULONG Nrandoms1=0;
#pragma omp parallel for reduction(+:Ndm_used_I, Nrandoms1)
   for(ULONG id=0;id<this->NGRID;++id)
     if(dm_cases[id]==1)
       {
         Ndm_used_I+=aux_cont1[id];
         Nrandoms1+=Nrandom_tracers[id];
       }
   // So.message_screen("Number of dm used in case I",Ndm_used_I);
   // So.message_screen("Number of randoms demanded from case I ", Nrandoms1);


//    ULONG Ndm_used_II=0;
// #pragma omp parallel for reduction(+:Ndm_used_II)
//    for(ULONG id=0;id<this->NGRID;++id)
//      Ndm_used_II+=aux_cont2[id];
   // So.message_screen("Number of dm used in case II",Ndm_used_II);
   // So.message_screen("(case II demands no randoms)");

   posx.clear();posx.shrink_to_fit();
   posy.clear();posy.shrink_to_fit();
   posz.clear();posz.shrink_to_fit();

   // So.message_screen("Total Number of dm used",Ndm_used_I+Ndm_used_II);


   // caso 3, celdas con tracers pero sin dm particles -> Todas random
   //   So.message_screen("Getting randoms from case 3");
#pragma omp parallel for
   for(ULONG id=0;id<this->NGRID;++id)
     {
       if(MOCK_DEN_FIELD[id]>0 && 0==aux_cont[id])
         if(dm_cases[id]!=1 && dm_cases[id]!=2)
           {
             dm_used[id]=false;
             Nrandom_tracers[id]=MOCK_DEN_FIELD[id]; //todas randoms
             dm_cases[id]=3;
           }
     }

   ULONG Nrandoms3=0;
#pragma omp parallel for reduction(+:Nrandoms3)
   for(ULONG id=0;id<this->NGRID;++id)
     if(MOCK_DEN_FIELD[id]>0 && aux_cont[id]==0)
       if(3==dm_cases[id])
         Nrandoms3+=Nrandom_tracers[id];

   // So.message_screen("Total number of randoms requested from case 3", Nrandoms3);
   // So.message_screen("Total number of randoms requested", Nrandoms1+Nrandoms3);
   // So.message_screen("Randoms requested + Ndm", Nrandoms3+Nrandoms1+Ndm_used_I+Ndm_used_II);
   // So.message_screen("Number of original tracers",Nobjects_mock);
   Nrandoms1+=Nrandoms3;

   ULONG ncells_check=0;
#pragma omp parallel for reduction(+:ncells_check)
   for(ULONG id=0;id<this->NGRID;++id)
     {
       if(dm_cases[id]>=1 && dm_cases[id]<=4)
         ncells_check++;
     }
   if(ncells_check!=this->NGRID)
     {
       So.message_screen("Cells in cases", ncells_check);
       So.message_screen("Total number of cells", this->NGRID);
       exit(0);
     }


   int jthread=1;
   vector<s_cell_info> cell_inf_ran(this->NGRID);
   vector<bool>random_used(this->NGRID,false);

  if(Nrandoms1>0)
    {

      So.message_screen("Creating random positions");
      for(ULONG i=0; i<this->N1; ++i)
        for(ULONG j=0; j<this->N2; ++j)
          for(ULONG k=0; k<this->N3; ++k)
            {
              ULONG id=index_3d(i,j,k,N2,N3);
              if(Nrandom_tracers[id]>0)
                {
                  cell_inf_ran[id].posx_p.resize(Nrandom_tracers[id],0);
                  cell_inf_ran[id].posy_p.resize(Nrandom_tracers[id],0);
                  cell_inf_ran[id].posz_p.resize(Nrandom_tracers[id],0);
                  cell_inf_ran[id].velx_p.resize(Nrandom_tracers[id],0);
                  cell_inf_ran[id].vely_p.resize(Nrandom_tracers[id],0);
                  cell_inf_ran[id].velz_p.resize(Nrandom_tracers[id],0);
                  random_used[id]=true;
                  for(int ir=0;ir<Nrandom_tracers[id];++ir)
                    {
#ifdef _USE_OMP_
                      real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                      real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
                      real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand[jthread]));
#else

                      real_prec rx= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                      real_prec ry= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
                      real_prec rz= static_cast<real_prec>(gsl_rng_uniform(gBaseRand));
#endif
                      real_prec x = this->params._d1()*static_cast<real_prec>(i) + this->params._d1()*rx;
                      real_prec y = this->params._d2()*static_cast<real_prec>(j) + this->params._d2()*ry;
                      real_prec z = this->params._d3()*static_cast<real_prec>(k) + this->params._d3()*rz;
                      cell_inf_ran[id].posx_p[ir]=x;
                      cell_inf_ran[id].posy_p[ir]=y;
                      cell_inf_ran[id].posz_p[ir]=z;
                      cell_inf_ran[id].velx_p[ir]=this->linInterp(x,y,z,velx);
                      cell_inf_ran[id].vely_p[ir]=this->linInterp(x,y,z,vely);
                      cell_inf_ran[id].velz_p[ir]=this->linInterp(x,y,z,velz);
                      cell_inf_ran[id].id_p.push_back(id);
                    }
                }
            }
      So.DONE();
    }
  velx.clear();velx.shrink_to_fit();
  vely.clear();vely.shrink_to_fit();
  velz.clear();velz.shrink_to_fit();


  // Before write to file we assign masses. This replaces the call of the function in bamrunner



  string outputFileName=this->params._Output_directory()+"CAT_realization"+to_string(ir)+"_"+this->stradd+string(".txt");
  this->fnameTRACERCAT=outputFileName;
  ofstream outStream;
  outStream.open(outputFileName.c_str());
  assert(outStream.is_open());
  outStream.precision(8);
  outStream.setf(ios::showpoint);
  outStream.setf(ios::scientific);

  this->So.message_screen("Writting dm particles to file", outputFileName);

  //write dms pos to file:
  ULONG Nobjects_mock_rc=0;

  for(ULONG id=0;id<this->NGRID;++id)
    {
      if(MOCK_DEN_FIELD[id]>0)
        if(true==dm_used[id])
          for(int in=0;in < cell_inf_dm[id].posx_p.size(); ++in)
            {
              Nobjects_mock_rc++;
              outStream<<cell_inf_dm[id].posx_p[in]<<"\t"<<cell_inf_dm[id].posy_p[in]<<"\t"<<cell_inf_dm[id].posz_p[in]<<"\t"<<cell_inf_dm[id].velx_p[in]<<"\t"<<cell_inf_dm[id].vely_p[in]<<"\t"<<cell_inf_dm[id].velz_p[in]<<"\t"<<0.0<<"\t"<<1<<endl;
            }
    }
  So.DONE();

  // So.message_screen("Number dm objecs written = ", Nobjects_mock_rc);

  if(Nrandoms1>0)
    {
      this->So.message_screen("Adding random particles:");
      for(ULONG id=0;id<this->NGRID;++id)
        {
          if(MOCK_DEN_FIELD[id]>0)
            if(true==random_used[id])
              for(int in=0;in<cell_inf_ran[id].posx_p.size();++in)
                {
                  Nobjects_mock_rc++;
                  outStream<<cell_inf_ran[id].posx_p[in]<<"\t"<<cell_inf_ran[id].posy_p[in]<<"\t"<<cell_inf_ran[id].posz_p[in]<<"\t"<<cell_inf_ran[id].velx_p[in]<<"\t"<<cell_inf_ran[id].vely_p[in]<<"\t"<<cell_inf_ran[id].velz_p[in]<<"\t"<<0.0<<"\t"<<-1<<endl;
                }
        }
      So.DONE();
    }


  outStream.close();
  if(Nobjects_mock_rc!=Nobjects_mock)
    {
      So.message_screen("Number of objecs written = ", Nobjects_mock_rc);
      So.message_screen("Original Number of objects = ", Nobjects_mock);
      exit(0);
    }

 }

// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================


 void LPT::comp_velbias(vector<real_prec> &delta, vector<real_prec>&out, bool zeropad, bool norm)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);

  real_prec Lzp1=L1;
  real_prec Lzp2=L2;
  real_prec Lzp3=L3;

  ULONG Nzp1=N1;
  ULONG Nzp2=N2;
  ULONG Nzp3=N3;

  if (zeropad==true)
    {
      Lzp1=L1*num_2;
      Lzp2=L2*num_2;
      Lzp3=L3*num_2;

      Nzp1=N1*2;
      Nzp2=N2*2;
      Nzp3=N3*2;
    }

  ULONG Nzp=Nzp1*Nzp2*Nzp3;


  real_prec cpecvel=num_1;
  if (norm==true)
    cpecvel=this->s_cosmo_info.growth_index*this->s_cosmo_info.Hubble_parameter/(1.+this->s_cosmo_pars.cosmological_redshift);



  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
#ifdef DOUBLE_PREC
  complex_prec *phi= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *phi= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif


  do_fftw_r2c(N1, delta, phi);


#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<N1;i++)
    for (ULONG j=0;j<N2;j++)
      for (ULONG k=0;k<N3/2+1;k++)
        {
          real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3);
          ULONG ll=index_3d(i,j,k,this->N2,this->N3/2+1);
          phi[ll][REAL]=(num_1+this->params._velbias()*this->params._velbias()*k2)*phi[ll][REAL];
          phi[ll][IMAG]=(num_1+this->params._velbias()*this->params._velbias()*k2)*phi[ll][IMAG];

        }
  do_fftw_c2r(this->N1,phi,out);

#ifdef DOUBLE_PREC
       fftw_free(phi);
#else
       fftwf_free(phi);
#endif
}

// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
// ======================================================================================================
void LPT::normalize_df_z_ini(vector<real_prec>&field,  string type, real_prec target_red_ini,real_prec target_red_end)
 {
   real_prec Dz = this->Cosmo.growth_factor(target_red_ini);
   real_prec gf = Dz/this->Cosmo.growth_factor(target_red_end);

#ifdef _UNITSIM_
    gf*=factor_growth;
#endif


#ifdef _FULL_VERBOSE_
   So.message_screen("Normalizing Initial density field using growth D(ini) =", Dz);
   So.message_screen("@ z inital = ",target_red_ini);
   So.message_screen("@ z end = ",target_red_end);
   So.message_screen("D(z_ini)/D(z_end) =", gf);
#endif

/*
   // Get the PDF if the IC:
   int Nk=600;
   vector<real_prec>xbaux(Nk, 0);
   vector<real_prec>pdf(Nk, 0);
   real_prec max=get_max<real_prec>(in);
   real_prec min=get_min<real_prec>(in);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<Nk; ++i)
     xbaux[i]=min+static_cast<real_prec>(i+0.5)*(max-min)/static_cast<real_prec>(Nk);
   string filex="pdf_X_MASX_IC.txt";
   calc_pdf("lin", in.size(),Nk,max,min,in,pdf);
   this->File.write_to_file(filex, xbaux,pdf);
*/


   // -----------------------------------------------------------------------------------------


   // This test is meant to FIX THE ISSUE OF THE MISSING POWER IN THE Ic
//#define _DO_TEST_

#ifdef _DO_TEST_

   So.message_screen("THIS IS A TEST");
   ULONG N1=static_cast<ULONG>(pow(in.size(),1./3.))+1;
   ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);
   fftw_complex *AUX= (fftw_complex *)fftw_malloc(2*Nhalf*sizeof(real_prec));
   do_fftw_r2c(N1,in,AUX);

   int NTHREADS = omp_get_max_threads();
   vector<real_prec>prop;
   ULONG nbinPS=this->File.read_file(this->params._dir()+this->params._ic_power_file(), prop, NTHREADS);
   int NCOLS=(static_cast<int>(prop.size()/nbinPS));
   vector<real_prec> kPS(nbinPS,0);
   vector<real_prec> powPS(nbinPS,0);
   for(ULONG i=0;i<nbinPS;i++)
     {
       kPS[i]=static_cast<real_prec>(prop[i*NCOLS]);
       powPS[i]=static_cast<real_prec>(prop[1+i*NCOLS]);
     }
   prop.clear();
   prop.shrink_to_fit();

   gsl_interp_accel *acc = gsl_interp_accel_alloc ();
   gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, kPS.size());
   gsl_spline_init (spline, &kPS[0], &powPS[0], kPS.size());

   real_prec Vol=L1*L2*L3;
   real_prec CONT_NORM=Vol;
   real_prec Norm=1.0;

#ifdef	  FOURIER_DEF_1
 Norm=static_cast<real_prec>(1./CONT_NORM);
#endif

#ifdef	  FOURIER_DEF_2
 Norm=static_cast<real_prec>(real_prec(in.size())*real_prec(in.size())/CONT_NORM);
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N1;i++)
    for(ULONG j=0;j<N1;j++)
      for(ULONG k=0;k<=N1/2;k++)
        {
          ULONG ihalf= index_3d(i,j,k,N1,N1/2+1);

          real_prec k2=k_squared(i,j,k,this->L1,this->L2,this->L3,this->N1,this->N2,this->N3);
          real_prec ktot=sqrt(k2);

          real_prec power_shells;
          if(i==0 && j==0 && k==0)
            power_shells=0;
           else
              power_shells=static_cast<real_prec>(gsl_spline_eval (spline, ktot, acc));

          // espectro de las IC con fixed amplitude
          real_prec Power3D=pow(AUX[ihalf][REAL],2)+pow(AUX[ihalf][IMAG],2);

          AUX[ihalf][REAL]*=sqrt(power_shells/Power3D);
          AUX[ihalf][IMAG]*=sqrt(power_shells/Power3D);
          //AUX[ihalf][REAL]=sqrt(power_shells);
          //AUX[ihalf][IMAG]=sqrt(power_shells);


        }

  do_fftw_c2r(N1,AUX,in);

  real_prec meanWN=get_mean(in);
  this->So.message_screen("Mean WN =",meanWN);
  real_prec sigma2D=get_var(meanWN, in);
  this->So.message_screen("Sigma_corr WN =",sqrt(sigma2D));

  fftw_free(AUX);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  So.message_screen("END OF TEST");

#endif
  // -----------------------------------------------------------------------------------------




  double ginv=1./static_cast<double>(gf);
// NOrmalising IC from initial z to z=0

 if(type=="delta")
 {
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<field.size();++i)
        field[i]*=ginv;
    }
  else if(type=="density")
    {
      real_prec mean=static_cast<real_prec>(get_nobjects(field))/pow(this->N1,3);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<field.size();++i)
        field[i]=mean*(1.-ginv)+ginv*field[i];
    }

   So.DONE();
 }

 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
 void LPT::DM_to_RSS(int los){

  int NTHREADS=1;
#ifdef _USE_OMP_
  NTHREADS=_NTHREADS_;
#endif

   ULONG Ndm=this->NGRID;
   vector<real_prec>Position(Ndm,0);
   vector<real_prec>Velocity(Ndm,0);

   switch(los){
      case(1):
        this->File.read_array(this->fnamePOSX+".dat",Position);
        this->File.read_array(this->fnameVXpart+".dat",Velocity);
        break;
      case(2):
          this->File.read_array(this->fnamePOSY+".dat",Position);
          this->File.read_array(this->fnameVYpart+".dat",Velocity);
      break;
      case(3):
          this->File.read_array(this->fnamePOSZ+".dat",Position);
          this->File.read_array(this->fnameVZpart+".dat",Velocity);
      break;
      }

      // This converts km /s to Mpc/h
    real_prec conversion_factor=1./(this->Cosmo.Hubble_function(this->s_cosmo_pars.cosmological_redshift));



    // Add peculiar velocities 
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0; i<Ndm; ++i)
          Position[i]+=Velocity[i]*conversion_factor;

        Velocity.clear(); Velocity.shrink_to_fit();

//    Add fingers of God
#ifdef _USE_FOGS_

   vector<real_prec>Position1(Ndm,0);
   vector<real_prec>Position2(Ndm,0);

   switch(los){
      case(1):
        this->File.read_array(this->fnamePOSY+".dat",Position1);
        this->File.read_array(this->fnamePOSZ+".dat",Position2);
        break;
      case(2):
        this->File.read_array(this->fnamePOSX+".dat",Position1);
        this->File.read_array(this->fnamePOSZ+".dat",Position2);
      break;
      case(3):
        this->File.read_array(this->fnamePOSX+".dat",Position1);
        this->File.read_array(this->fnamePOSY+".dat",Position2);
      break;
      }



   s_params_box_mas box;
   box.min1=this->params._xmin();
   box.min2=this->params._ymin();
   box.min3=this->params._zmin();
   box.Lbox=this->params._Lbox();
   box.Nft=this->params._Nft();
   box.d1= box.Lbox/static_cast<real_prec>(box.Nft);   /* grid spacing x-direction */
   box.d2= box.d1;
   box.d3= box.d1;
   box.NGRID=(box.Nft*box.Nft*box.Nft);

 //Read the density field
  vector<real_prec>density(this->NGRID,0);
  this->File.read_array(this->fnameDM, density);  //write DMDF to binary file
  get_overdens(density, density);

  int jthread=0;
  const gsl_rng_type *Trn;
  gsl_rng *gBaseRand;
  vector<ULONG>vseeds(NTHREADS,0);
  for(int i=0;i<vseeds.size();++i)
    vseeds[i]=static_cast<ULONG>(i+14)*56045;
#ifdef _USE_OMP_
#pragma omp parallel private (jthread, Trn, gBaseRand)// this is causing problem
   {
#endif

    jthread=omp_get_thread_num();
    gsl_rng_default_seed=vseeds[jthread];
    Trn = gsl_rng_mt19937;
    gBaseRand= gsl_rng_alloc (Trn);

   switch(los){
    case(1):
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0; i<Ndm; ++i)
           {
              ULONG iindex=grid_ID(&box,Position[i],Position1[i],Position2[i]);
              real_prec deltaloc=density[iindex];
              real_prec sigma=pow(1+this->params._biasL()*deltaloc,this->params._ep())*this->params._sfac();
              Position[i]+= gsl_ran_gaussian(gBaseRand,sigma)*conversion_factor;
            }
      break;
     case(2):
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0; i<Ndm; ++i)
           {
              ULONG iindex=grid_ID(&box,Position1[i],Position[i],Position2[i]);
              real_prec deltaloc=density[iindex];
              real_prec sigma=pow(1+this->params._biasL()*deltaloc,this->params._ep())*this->params._sfac();
              Position[i]+= gsl_ran_gaussian(gBaseRand,sigma)*conversion_factor;
            }
     break;
     case(3):
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0; i<Ndm; ++i)
           {
              ULONG iindex=grid_ID(&box,Position1[i],Position2[i],Position[i]);
              real_prec deltaloc=density[iindex];
              real_prec sigma=pow(1+this->params._biasL()*deltaloc,this->params._ep())*this->params._sfac();
              Position[i]+= gsl_ran_gaussian(gBaseRand,sigma)*conversion_factor;
            }
     break;
      }

#ifdef _USE_OMP_
  }
#endif

 
#endif // end for use fogs



   string output_file;
   switch(los)
    {
      case(1):
         output_file= this->fnamePOSX_RS;
        break;
      case(2):
        output_file= this->fnamePOSY_RS;
        break;
      case(3):
        output_file= this->fnamePOSZ_RS;
        break;
    }

      this->File.write_array(output_file,Position);
 }
 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
 // ======================================================================================================
#endif
