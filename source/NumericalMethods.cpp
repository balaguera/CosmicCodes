////////////////////////////////////////////////////////////////////////////
/**
 * @brief This file contains functions for numerical analysis
 * @file NumericalMethods.cpp
 * @author Andres Balaguera Antolinez
 * @version 1.0
 * @date 2017-2024
 */
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
# include "../headers/NumericalMethods.h"
using namespace std;
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
inline ULONG factorial(int n){
  if(n<=0 || n==1 ) return 1;
  else return tgammal(n+1);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void randomize_vector(vector<ULONG>&data)
{
    const gsl_rng_type *Tn = gsl_rng_default;
    gsl_rng *rn = gsl_rng_alloc (Tn);
    gsl_rng_default_seed=155;
    int n_aux=data.size();
    vector<ULONG>cells_id_still_to_assign_aux(n_aux,0);  //container to allocate the ID of the cells with particles still to get masses
    vector<bool>cells_id_still_to_assign_chosen(n_aux,false);  //container to allocate the ID of the cells with particles still to get masses
    for(int ide=0; ide < n_aux ;++ide)
      {
        bool flag=false;
        while(false==flag)
          {
            int jk= gsl_rng_uniform_int(rn,n_aux);
            ULONG new_id=data[jk];
            if(cells_id_still_to_assign_chosen[jk]==false)
              {
                cells_id_still_to_assign_aux[ide]=new_id;  //container to allocate the ID of the cells with particles still to get masses
                cells_id_still_to_assign_chosen[jk]=true;  //container to allocate the ID of the cells with particles still to get masses
                flag=true;
              }
          }
      }
    for(ULONG ide=0; ide < n_aux ;++ide)
      data[ide]=cells_id_still_to_assign_aux[ide];  //container to allocate the ID of the cells with particles still to get masses
    cells_id_still_to_assign_aux.clear();
    cells_id_still_to_assign_aux.shrink_to_fit();
    cells_id_still_to_assign_chosen.clear();
    cells_id_still_to_assign_chosen.shrink_to_fit();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
// This is the same function write_array. Defined here to be used in the functions inherited from FSK codes
void dump_scalar(const vector<real_prec>&A_rm,ULONG N1,ULONG N2,ULONG N3,int sample_number,string fname)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
   fftw_array<real_prec> dummy(N);

#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
    dummy[i]=A_rm[i];
  
  string FNAME=fname+string(".dat");
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Writting binary file "<<CYAN<<FNAME<<RESET<<endl;
#endif
  bofstream outStream(FNAME.data(),file_is_natural);
  assert(outStream.is_open());
  outStream.put(dummy.data,N);
  outStream.close();
#ifdef _FULL_VERBOSE_
      std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
string dto_string (double Number)
{
  stringstream ss;
  ss << Number;
  return ss.str();
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *p,gsl_real LowLimit,gsl_real UpLimit)
{
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (1000);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  real_prec result=static_cast<real_prec>(gsl_integration_glfixed(&F,LowLimit,UpLimit,wf));
  gsl_integration_glfixed_table_free(wf);
  return result;
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration(int N, gsl_real(*function)(gsl_real, void *) ,void *p, gsl_real LowLimit, gsl_real UpLimit)
{
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (N);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  real_prec result=static_cast<real_prec>(gsl_integration_glfixed(&F,LowLimit,UpLimit,wf));
  gsl_integration_glfixed_table_free(wf);
  return result;
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *p,vector<gsl_real>&XX, vector<gsl_real>&WW, bool paral)
{
  real_prec result=0;
  int nn=WW.size();

    if(paral==true)
    {
#ifdef _USE_OMP_
  omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for reduction(+:result)
#endif
  for(ULONG i=0;i<nn;++i)
     result+=WW[i]*function(XX[i],p);
}
    else
    {
     for(ULONG i=0;i<nn;++i)
        result+=WW[i]*function(XX[i],p);
    }
  return static_cast<real_prec>(result);
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *p,vector<gsl_real>&XX, vector<gsl_real>&WW, real_prec MMin, bool paral)
{
  real_prec result=0;
  int nn=WW.size();

    if(paral==true)
    {
#ifdef _USE_OMP_
  omp_set_num_threads(omp_get_max_threads());
#pragma omp parallel for reduction(+:result)
#endif
  for(ULONG i=0;i<nn;++i)
     result+=WW[i]*function(XX[i],p);
}
    else
    {
     for(ULONG i=0;i<nn;++i)
         if(XX[i]>=MMin)
     {
        result+=WW[i]*function(XX[i],p);
      }
     }
  return static_cast<real_prec>(result);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration(gsl_real (*function)(gsl_real, void *) ,void *p,vector<gsl_real>&XX, vector<gsl_real>&WW)
{
  real_prec result=0;
  int nn=WW.size();
#ifdef _USE_OMP_
  omp_set_num_threads(omp_get_max_threads());
#endif

#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:result)
#endif
  for(ULONG i=0;i<nn;++i)
     result+=WW[i]*function(XX[i],p);

  return static_cast<real_prec>(result);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void gsl_get_GL_weights(gsl_real LowLimit,gsl_real UpLimit, vector<real_prec> &XX, vector<real_prec>&WW)
{
  ULONG nn=XX.size();
  gsl_integration_glfixed_table *wf = gsl_integration_glfixed_table_alloc (nn);
  for(int i=0;i<nn;++i)
  {
    gsl_real xi, wi;
    gsl_integration_glfixed_point(LowLimit,UpLimit, i, &xi, &wi, wf);
    XX[i]=static_cast<real_prec>(xi);
    WW[i]=static_cast<real_prec>(wi);
  }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void gsl_get_GL_weights(gsl_real LowLimit,gsl_real UpLimit,  gsl_integration_glfixed_table *wf, vector<gsl_real> &XX, vector<gsl_real>&WW)
{
  gsl_real xi, wi;
  ULONG nn=XX.size();
#ifdef _USE_OMP_
  omp_set_num_threads(omp_get_max_threads());
#endif

#ifdef _USE_OMP_
//#pragma omp parallel for
#endif
  for(int i=0;i<nn;++i){
    gsl_integration_glfixed_point(LowLimit,UpLimit, i, &xi, &wi, wf);
    XX[i]=xi;
    WW[i]=wi;
  }
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration_sin_kernel(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit, gsl_real UpLimit){
  int NIT=1e3;
  gsl_real result, error;
  gsl_real L=UpLimit-LowLimit;
  gsl_real relerr=0.1;
  gsl_real abserr=0.01;
  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);   
  // La subrutina qawo integra una funcion f peseada con un sin(omega x), 
  // en un intervalo definido
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)
  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;
  for(;;){
    if(gsl_integration_qawo(&F,LowLimit,abserr,relerr,NIT,w,wf,&result,&error)){
      abserr*=10;
    }
    else{
      break;
    }
  }
  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return static_cast<real_prec>(result);
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration_sin_kernel_lowlimit_inf2(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit, gsl_integration_workspace *w,gsl_integration_workspace *wc ){
  gsl_real result, error;
  int NIT=1e3;
  gsl_real L=1;
  gsl_real relerr=0.1;
  gsl_real abserr=0.01;
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);
  // La subrutina qawo integra una funcion f peseada con un sin(omega x),
  // en un intervalo definido
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)

  gsl_function F;
  F.params   = p;
  F.function = function;
  result=0;
  for(;;){
    if(gsl_integration_qawf(&F,LowLimit,relerr,NIT,w,wc,wf,&result,&error)){
      abserr*=10;
    }
    else{
      break;
    }
  }
  return static_cast<real_prec>(result);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_integration_sin_kernel_lowlimit_inf(gsl_real (*function)(gsl_real, void *), void *p, gsl_real omega, gsl_real LowLimit){
  int NIT=1e3;
  gsl_real result, error;
  gsl_real L=1;
  gsl_real relerr=1.;
  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT);
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT);
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,NIT);

  // La subrutina qawo integra una funcion f peseada con un sin( omega x), 
  // el cual es en nuestro caso el que viene de j0(x)
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)
  // Integramos en k para poder usar esta subrutina 
  // utiliza w = 1 con eta=kr, eta lo llamamos k arriba en las funcciones
  
  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;

   for(;;){
    if(gsl_integration_qawf(&F,LowLimit,relerr,NIT,w,wc,wf,&result,&error)){
      relerr*=10;
    }
    else{
      break;
    }
  }
  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return static_cast<real_prec>(result);
}    
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_inter_pointers(real_prec *x, real_prec *y, int n, real_prec xx){
  real_prec ans;
  gsl_real xa[n];
  gsl_real ya[n];
  for(int i=0;i<n;i++)xa[i]=*(x+i);
  for(int i=0;i<n;i++)ya[i]=*(y+i);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, xa, ya, n);
  ans=(xx<*x? *x : gsl_spline_eval (spline, xx, acc));
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return ans;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef SINGLE_PREC
real_prec gsl_inter_pointers(double *x, double *y, int n, double xx){
  double ans;
  gsl_real xa[n];
  gsl_real ya[n];
  for(int i=0;i<n;i++)xa[i]=*(x+i);
  for(int i=0;i<n;i++)ya[i]=*(y+i);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, xa, ya, n);
  ans=(xx<*x? *x : gsl_spline_eval (spline, xx, acc));
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return ans;
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void gsl_bspline(vector<gsl_real>&xx, vector<gsl_real>&yy, vector<gsl_real>&new_xx, vector<gsl_real> &new_yy)
{
  int n=xx.size();
  int new_n=new_yy.size();
  gsl_real x_ini=static_cast<gsl_real>(xx[0]);
  gsl_real x_fin=static_cast<gsl_real>(xx[n-1]);
  gsl_matrix *X, *cov;
  gsl_vector *c, *w;
  gsl_vector *x, *y;
  int k = 4;   /*cubic spline stands for k=4*/
  int nbreak=new_n+2-k;
  gsl_bspline_workspace *bw;
  gsl_multifit_linear_workspace *mw;
  gsl_real chisq;
  gsl_vector *B;
  X= gsl_matrix_alloc(n,new_n);
  cov= gsl_matrix_alloc(new_n,new_n);
  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  X = gsl_matrix_alloc(n, new_n);
  c = gsl_vector_alloc(new_n);
  w = gsl_vector_alloc(n);
  mw = gsl_multifit_linear_alloc(n, new_n);
  bw=gsl_bspline_alloc(k,nbreak);
  B=gsl_vector_alloc(new_n);
  gsl_bspline_knots_uniform(x_ini, x_fin,bw);
  for(int i=0;i<n;i++)gsl_vector_set(x,i,xx[i]);
  for(int i=0;i<n;i++)gsl_vector_set(y,i,yy[i]);
  for(int i=0;i<n;i++){
    gsl_vector_set(w,i,1e7); /**/
    gsl_real xi=gsl_vector_get(x,i);
    gsl_bspline_eval(xi,B,bw);
    for(int j=0;j<new_n;j++){
      gsl_real Bj=gsl_vector_get(B,j);
      gsl_matrix_set(X,i,j,Bj);
    }
  }
  gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
  {
    gsl_real xi, yi, yerr;
    for(int i=0;i<new_n; i++){
      xi=x_ini+i*(x_fin-x_ini)/(new_n-1);
      gsl_bspline_eval(xi, B, bw);
      gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
      new_xx[i]=xi;
      new_yy[i]=yi;
    }
  }
//  int dof = n - new_n;
//  gsl_real tss = gsl_stats_wtss(w->data, 1, y->data, 1, y->size);
//  gsl_real Rsq = 1.0 - chisq / tss;
//  cout<<"Rsq = "<<Rsq<<endl;
//  cout<<"chisq/dof = "<<chisq/dof<<endl;
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  return; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_inter(gsl_real *x, gsl_real *y, int n, gsl_real xn){
  real_prec ans;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_inter_new(const vector<gsl_real> &xx, const vector<gsl_real> &yy, gsl_real xn){
  int n=xx.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, &xx[0], &yy[0], n);
  real_prec ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef SINGLE_PREC
real_prec gsl_inter_new(const vector<float> &xx, const vector<float> &yy, gsl_real xn){
  int n=xx.size();
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear, n);
  vector<double>xxd(xx.size(),0); vector<double>yyd(yy.size(),0);
  for(int i =0;i<xx.size();++i)xxd[i]=static_cast<double>(xx[i]);
  for(int i =0;i<yy.size();++i)yyd[i]=static_cast<double>(yy[i]);
  gsl_spline_init(spline, &xxd[0], &yyd[0], n);
  real_prec ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans;
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec gsl_inter_new2(vector<real_prec> &xx, vector<vector<real_prec> > &yy, int li, real_prec xn){
  real_prec ans;
  int n=xx.size();
  gsl_real x[n],y[n];
  for(int i=0;i<n;++i)
  {
    x[i]=static_cast<gsl_real>(xx[i]);
    y[i]=static_cast<gsl_real>(yy[li][i]);
  }
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef SINGLE_PREC
real_prec gsl_inter_new2(vector<gsl_real> &xx, vector<vector<gsl_real> > &yy, int li, real_prec xn){
  gsl_real ans;
  int n=xx.size();
  gsl_real x[n],y[n];
  for(int i=0;i<n;++i){x[i]=static_cast<gsl_real>(xx[i]);y[i]=static_cast<gsl_real>(yy[li][i]);}
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  static_cast<real_prec>(ans); 
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void sort_index(int i, int j, int k, int *ii, int *jj, int *kk){
  vector<real_prec> maa;
  maa.push_back(i);
  maa.push_back(j);
  maa.push_back(k);
  sort(maa.begin(), maa.end());
  *ii=maa.at(0);
  *jj=maa.at(1);
  *kk=maa.at(2);
  return ;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void matrix_inversion(vector< vector<real_prec> > &G, vector< vector<real_prec> > &y){
  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  gsl_matrix *Ai = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)
      for(int j=0;j<n;++j)
	gsl_matrix_set(A, i,j,static_cast<gsl_real>(G[i][j]));
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;
  gsl_linalg_LU_decomp(A, p, &signum);
  gsl_linalg_LU_invert(A, p, Ai);
  gsl_permutation_free(p);
  for (int i=0;i<n;++i)for(int j=0;j<n;++j)y[i][j]=gsl_matrix_get(Ai, i,j);
  gsl_matrix_free (A);
  gsl_matrix_free (Ai);    
  return ;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void matrix_det(vector< vector<real_prec> > &G, real_prec &det){
  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)
      for(int j=0;j<n;++j)
          gsl_matrix_set(A, i,j,G[i][j]);
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;
  gsl_linalg_LU_decomp(A, p, &signum);
  det=gsl_linalg_LU_det(A, signum);
  gsl_permutation_free(p);
  gsl_matrix_free (A);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_eigen(vector<vector<real_prec>>&icov_masses, vector<real_prec>&masses){
  gsl_matrix *icova;
  gsl_vector *eig;
  gsl_eigen_symm_workspace *workspace;
  int n_par=masses.size();
  icova = gsl_matrix_alloc(n_par,n_par);
  eig = gsl_vector_alloc(n_par);
  for(int i=0;i<n_par;i++)
      for(int j=0;j<n_par;j++)
          gsl_matrix_set(icova, i,j, icov_masses[i][j]);
  workspace = gsl_eigen_symm_alloc(n_par);
  gsl_eigen_symm(icova, eig, workspace);
  gsl_sort_vector(eig);
  gsl_eigen_symm_free(workspace);
  for(int ip=0;ip<n_par;ip++)masses[ip]=sqrt(gsl_vector_get(eig,ip));
  gsl_matrix_free(icova);
  gsl_vector_free(eig);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_eigen(vector<real_prec>&in_matrix, vector<s_eigen_vector>&eigen){
  int n_par=eigen[0].eigen_vec.size();
  gsl_matrix *evec = gsl_matrix_alloc (n_par,n_par);
  gsl_matrix *aux_matrix = gsl_matrix_alloc(n_par,n_par);
  gsl_vector *eig = gsl_vector_alloc(n_par);
  gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(n_par);
  for(int i=0;i<n_par;i++)
      for(int j=0;j<n_par;j++)
          gsl_matrix_set(aux_matrix, i,j, in_matrix[i+j*n_par]);
  gsl_eigen_symmv(aux_matrix, eig, evec, workspace);
  gsl_eigen_symmv_sort(eig,evec,GSL_EIGEN_SORT_VAL_DESC);
  gsl_eigen_symmv_free(workspace);
  for(int ip=0;ip<n_par;ip++)
  {
    eigen[ip].eigen_val=gsl_vector_get(eig,ip); // sorterd in descending order
    for(int jp=0;jp<n_par;jp++)
      eigen[ip].eigen_vec[jp]=gsl_matrix_get(evec,ip,jp);
   }
  gsl_matrix_free(evec);
  gsl_matrix_free(aux_matrix);
  gsl_vector_free(eig);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_eigen(vector<real_prec>&in_matrix, vector<real_prec>&eigen){
  int n_par=eigen.size();
  gsl_matrix *evec = gsl_matrix_alloc (n_par,n_par);
  gsl_matrix *aux_matrix = gsl_matrix_alloc(n_par,n_par);
  gsl_vector *eig = gsl_vector_alloc(n_par);
  gsl_eigen_symmv_workspace *workspace = gsl_eigen_symmv_alloc(n_par);

  for(int i=0;i<n_par;i++)
      for(int j=0;j<n_par;j++)
          gsl_matrix_set(aux_matrix, i,j, in_matrix[i+j*n_par]);
  gsl_eigen_symmv(aux_matrix, eig, evec, workspace);
  gsl_eigen_symmv_sort(eig,evec,GSL_EIGEN_SORT_VAL_DESC);
  gsl_eigen_symmv_free(workspace);
  for(int ip=0;ip<n_par;ip++)
    eigen[ip]=gsl_vector_get(eig,ip); // sorterd in descending order

  gsl_matrix_free(evec);
  gsl_matrix_free(aux_matrix);
  gsl_vector_free(eig);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_det_matrix(vector<vector<real_prec>>&matriz, real_prec &determinant){
  vector<real_prec>det(matriz.size(),0);
  get_eigen(matriz, det);
  determinant=1;         ;//Para escalar
  for(int i=0;i<matriz.size();i++)
     determinant*=det[i];
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec real_sh(int l, int m, real_prec theta, real_prec phi){
  return  gsl_sf_legendre_sphPlm(l,m,cos(theta))*cos(m*phi);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec imag_sh(int l, int m, real_prec theta, real_prec phi){
  return gsl_sf_legendre_sphPlm(l,m,cos(theta))*sin(m*phi);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec bessel(real_prec x, int n){
  return gsl_sf_bessel_jl(n,x);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec dbessel(real_prec x, int n){
  return  (n*bessel(x,n-1)-(n+1)*bessel(x,n+1))/(2.*n+1);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec ddbessel(real_prec x, int n){
  return  -(2/x)*dbessel(x,n)-bessel(x,n)*(1.0-n*(n+1)/pow(x,2));
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
inline void SWAP_L(int &a, int &b) {
  int dum=a;
  a=b;
  b=dum;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
inline void SWAP_LU(ULONG &a, ULONG &b) {
  ULONG dum=a;
  a=b;
  b=dum;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void indexx(vector<int>&arr, vector<int>&indx)
{
  const int M=7,NSTACK=50;
  int n = indx.size();
  int i,indxt,ir,j,k,jstack=-1,l=0;
  int a;
  int istack[NSTACK];
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
	       a=arr[indxt];
	       for (i=j-1;i>=l;i--) 
          {
	          if (arr[indx[i]] <= a) break;
       	  indx[i+1]=indx[i];
          }
	      indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_L(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
        	SWAP_L(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
	      SWAP_L(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
       	SWAP_L(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
    	 do i++; while (arr[indx[i]] < a);
	     do j--; while (arr[indx[j]] > a);
	     if (j < i) break;
	       SWAP_L(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
	       cerr << "NSTACK too small in indexx." << endl;
      	exit(1);
      }
      if (ir-i+1 >= j-l) {
	      istack[jstack]=ir;
	      istack[jstack-1]=i;
      	ir=j-1;
      } 
      else {
	     istack[jstack]=j-1;
       istack[jstack-1]=l;
    	l=i;
     }
    }
  }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void indexx_ulong(vector<ULONG>&arr, vector<ULONG>&indx)
{
  const int M=7,NSTACK=50;
  ULONG n = indx.size();
  ULONG i,indxt,ir,j,k,jstack=-1,l=0;
  ULONG a;
  int istack[NSTACK];
  ir=n-1;
  for (j=0;j<n;j++) indx[j]=j;
  for (;;) {
    if (ir-l < M)
     {
       for (j=l+1;j<=ir;j++) 
       {
         indxt=indx[j];
         a=arr[indxt];
         for (i=j-1;i>=l;i--) 
          {
            if (arr[indx[i]] <= a) break;
          indx[i+1]=indx[i];
          }
        indx[i+1]=indxt;
       }
      if (jstack < 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
      } 
      else 
      {
       k=(l+ir) >> 1;
       SWAP_LU(indx[k],indx[l+1]);
        if (arr[indx[l]] > arr[indx[ir]]) {
          SWAP_LU(indx[l],indx[ir]);
        }
      if (arr[indx[l+1]] > arr[indx[ir]]) {
        SWAP_LU(indx[l+1],indx[ir]);
      }
      if (arr[indx[l]] > arr[indx[l+1]]) {
        SWAP_LU(indx[l],indx[l+1]);
      }
      i=l+1;
      j=ir;
      indxt=indx[l+1];
      a=arr[indxt];
      for (;;) {
       do i++; while (arr[indx[i]] < a);
       do j--; while (arr[indx[j]] > a);
       if (j < i) break;
         SWAP_LU(indx[i],indx[j]);
      }
      indx[l+1]=indx[j];
      indx[j]=indxt;
      jstack += 2;
      if (jstack >= NSTACK) {
         cerr << "NSTACK too small in indexx." << endl;
        exit(1);
      }
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } 
      else {
       istack[jstack]=j-1;
       istack[jstack-1]=l;
      l=i;
     }
    }
  }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void sort_vectors(vector<ULONG>& v1,vector<ULONG>&v2, vector<ULONG>&v3, vector<ULONG>& v4, vector<ULONG>&v5, vector<ULONG>& v6, vector<ULONG>& v7, vector<real_prec>& v8 ){
  // v1 is to be sorted. The elements of the other vectors
  // are shuffled accordingly.
  unsigned long j;
  long n=v1.size();
  vector<ULONG>iwksp(n,0);
  vector<float>wksp(n,0);
  indexx_ulong(v1,iwksp);
  for (j=0;j<n;++j) wksp[j]=v1[j];
  for (j=0;j<n;++j) v1[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v2[j];
  for (j=0;j<n;++j) v2[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v3[j];
  for (j=0;j<n;++j) v3[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v4[j];
  for (j=0;j<n;++j) v4[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v5[j];
  for (j=0;j<n;++j) v5[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v6[j];
  for (j=0;j<n;++j) v6[j]=wksp[iwksp[j]];
  for (j=0;j<n;++j) wksp[j]=v7[j];
  for (j=0;j<n;++j) v7[j]=wksp[iwksp[j]];
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
size_t m_getBuffer(istream& inStream) {
    // reset buffer
    inStream.clear();
    inStream.seekg(0, inStream.beg);
    string line = "";
    // get line from stream until is not a comment or empty
    while (false == inStream.eof() && ((true == line.empty()) || ('#' == line[0]))) {
      getline(inStream, line);
    }
    // even number of lines to read as 3MB buffer
    size_t numLines = 3*1024*1024/line.size();
    numLines += numLines&0x1;
    // reset buffer
    inStream.clear();
    inStream.seekg(0, inStream.beg);
    return numLines;
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
size_t m_countLines(istream& inStream) {
  // reset stream to beginning and clear status
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  // counter of lines
  size_t counter = -1;
  // line to read
  string line = "";
  while (false == inStream.eof()) {
    getline(inStream, line);
    counter++;
  }  
  // reset buffer
  inStream.clear();
  inStream.seekg(0, inStream.beg);
  return counter;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_mean(const vector<real_prec>& ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  double mean=0.;

  if(ini.size()>1)
  {
#ifdef _USE_OMP_
//#pragma omp parallel for reduction(+:mean)
#endif
  for(ULONG i=0;i<ini.size();++i)
  {
  if(isinf(ini[i]))cout<<i<<endl;
    mean+=static_cast<double>(ini[i]);
  }
  mean/=static_cast<double>(ini.size());
  }
 else
  mean=ini[0];
return static_cast<real_prec>(mean);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_var(const vector<real_prec> & ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  real_prec mean=get_mean(ini);
  double vari=0.;
#ifdef _USE_OMP_
#pragma opmp parallel for reduction(+:vari)
#endif
  for(ULONG i=0;i<ini.size();++i)
    vari+=pow(ini[i]-mean,2);
  vari/=(static_cast<double>(ini.size())-1.0);
  return static_cast<real_prec>(vari);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_var(real_prec mean, const vector<real_prec> & ini)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  real_prec vari=0.;

  if(ini.size()>0){
#ifdef _USE_OMP_
#pragma opmp parallel for reduction(+:vari)
#endif
  for(ULONG i=0;i<ini.size();++i)
    vari+=pow(ini[i]-mean,2);
  vari/=(static_cast<double>(ini.size())-1.0);
    }
  return vari;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec pearson_correlation(vector<real_prec>&X,vector<real_prec>&Y){
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif

  if(X.size()!=Y.size())
  {
     cout<<RED<<"Error in sze of files"<<endl;
  }

  real_prec corr=0.;
  if(X.size()>0)
  {
  real_prec Xmean=get_mean(X);
  real_prec Ymean=get_mean(Y);
  real_prec Xvar=get_var(Xmean,X);
  real_prec Yvar=get_var(Ymean,Y);

  if(Yvar>0 && Xvar>0)
  {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:corr)
#endif
  for(ULONG i=0;i<X.size();++i)
    corr+=(X[i]-Xmean)*(Y[i]-Ymean);
   corr/=(sqrt(Xvar*Yvar)*static_cast<real_prec>(X.size()));
    }
  }
   return corr;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec pearson_correlation(real_prec Xmean,vector<real_prec>&X,real_prec Ymean,vector<real_prec>&Y){
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  if(X.size()!=Y.size())
  {
     cout<<RED<<"Error in sze of files"<<endl;
  }
  real_prec corr=0.;
  if(X.size()>0)
  {
  real_prec Xvar=get_var(Xmean,X);
  real_prec Yvar=get_var(Ymean,Y);
  if(Yvar>0 && Xvar>0)
  {
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:corr)
#endif
  for(ULONG i=0;i<X.size();++i)
    corr+=(X[i]-Xmean)*(Y[i]-Ymean);
   corr/=(sqrt(Xvar*Yvar)*static_cast<real_prec>(X.size()));
    }
  }
   return corr;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_min_nm(const vector<real_prec>&in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  double lka=static_cast<double>(LARGE_NUMBER);
#pragma omp parallel for reduction(min:lka)
  for(ULONG i=0;i<in.size();++i)
      lka=min(static_cast<double>(in[i]),lka);
  return static_cast<real_prec>(lka);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_min_nm(const vector<real_prec>&in, bool zero)
{
 if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  double lka=static_cast<double>(LARGE_NUMBER);
  double lkb;
  if(zero){
  for(ULONG i=0;i<in.size();++i)
    {
      if(in[i]!=0)
        {
          lkb=min(static_cast<double>(in[i]),lka);
          lka=lkb;
            }
    }
  }
  else
  for(ULONG i=0;i<in.size();++i)
    {
          lkb=min(static_cast<double>(in[i]),lka);
          lka=lkb;
    }
  return static_cast<real_prec>(lka);  
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec get_max_nm(const vector<real_prec> &in)
{
  if(in.size()==0)
    cerr<<"Error. Empty vector"<<endl;
  double lka=-static_cast<double>(LARGE_NUMBER);
#pragma omp parallel for reduction(max:lka)
  for(ULONG i=0;i<in.size();++i)
      lka=max(static_cast<double>(in[i]), lka);
  return static_cast<real_prec>(lka);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void log_smooth(vector<real_prec>&xin, vector<real_prec>&vin)
{
  // steps:
  // i ) interpolate the raw ratio in log bins in k
  // ii) smooth the result
  // iii) interpolate the smoothed version and update power_ratio
  // Original modes, converted to log
  vector<gsl_real>lkvec(vin.size(),0);
  vector<gsl_real>vin_new(vin.size(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    lkvec[i]=static_cast<gsl_real>(log10(xin[i]));
  for(ULONG i = 0;i < vin.size(); ++i)
    vin_new[i]=static_cast<gsl_real>(vin[i]);
  int n_smooth=vin.size()*2;
  // Begin smooth ker_aux in log k
  // Log bins, in the same range as the original modes
  real_prec deltalk= static_cast<real_prec>(log10(xin[xin.size()-1.0]/xin[0])/static_cast<real_prec>(n_smooth-1.0));
  vector<gsl_real>aux_logk(n_smooth,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0; i < n_smooth-1; ++i)
    aux_logk[i]=log10(xin[0])+i*deltalk;
  aux_logk[n_smooth-1]=lkvec[xin.size()-1];
    // Interpolate original ratio in log k
  vector<gsl_real>aux_oratio(n_smooth,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < n_smooth; ++i)
    aux_oratio[i]=gsl_inter_new(lkvec,vin_new, aux_logk[i]);
  // Bspline smooth ratio, already interpolated in log
  int n_smooth_new=  static_cast<int>(floor(vin.size()/2));    //60;
  vector<gsl_real>aux_kvec(n_smooth_new, 0);
  vector<gsl_real>aux_ratio(n_smooth_new, 0);
  gsl_bspline(aux_logk, aux_oratio,  aux_kvec, aux_ratio);
  // Reassign the smoothed version to Kernel_bins via interpolation
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    vin[i]=gsl_inter_new(aux_kvec, aux_ratio, lkvec[i]);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void lin_smooth(vector<real_prec>&xin, vector<real_prec>&vin, int nr)
{
  // steps:
  // i ) interpolate the raw ratio in log bins in k
  // ii) smooth the result
  // iii) interpolate the smoothed version and update power_ratio
  // Original modes, converted to log
  vector<gsl_real>lkvec(vin.size(),0);
  vector<gsl_real>vin_new(vin.size(),0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < xin.size(); ++i)
    lkvec[i]=static_cast<gsl_real>(xin[i]);
  for(ULONG i = 0;i < vin.size(); ++i)
    vin_new[i]=static_cast<gsl_real>(vin[i]);
  /*
  int n_smooth=static_cast<int>(floor(vin.size()/3));
  // Begin smooth ker_aux in log k
  // Log bins, in the same range as the original modes
  real_prec deltalk= (xin[xin.size()-1]-xin[0])/static_cast<real_prec>(n_smooth-1);
  vector<gsl_real>aux_logk(n_smooth,xin[0]);
  aux_logk[n_smooth-1]=lkvec[xin.size()-1];
#pragma omp parallel for
  for(ULONG i = 1; i < n_smooth; ++i)
      aux_logk[i]=xin[0]+static_cast<real_prec>(i)*deltalk;
  // Interpolate original ratio in log k
  vector<gsl_real>aux_oratio(n_smooth,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < n_smooth; ++i)
    aux_oratio[i]=gsl_inter_new(lkvec,vin_new, aux_logk[i]);
  // Bspline smooth ratio, already interpolated in log
*/
  int n_smooth_new= static_cast<int>(floor(vin.size())/nr);
  vector<gsl_real>aux_kvec(n_smooth_new, 0);
  vector<gsl_real>aux_ratio(n_smooth_new, 0);
//  gsl_bspline(aux_logk, aux_oratio,  aux_kvec, aux_ratio);
 gsl_bspline(lkvec, vin_new,  aux_kvec, aux_ratio);
 // Reassign the smoothed version to Kernel_bins via interpolation
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i = 0;i < vin.size(); ++i)
    vin[i]=gsl_inter_new(aux_kvec, aux_ratio, lkvec[i]);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void gradfindif(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3,vector<real_prec>in,vector<real_prec>&out,int dim)
{
    // Gradient using finite differences O(∆x4) centered diﬀerence approximations
    real_prec deltaH=  12.0*(static_cast<real_prec>(L1)/static_cast<real_prec>(N1));
    real_prec ideltaH = 1.0/static_cast<real_prec>(deltaH);
    //    real_prec factor=static_cast<real_prec>(N1)/static_cast<real_prec>(num_2*L1);
  switch(dim){
  case(1):
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
    for(int x = 0; x < N1; x++)
    for(int y = 0; y < N2; y++)
      for(int z = 0; z < N3; z++)
	    {
          int xr = x + 1;
          int xl = x - 1;
          int xrr = x + 2;
          int xll = x - 2;
	      if(xr >= N1)
            xr -= N1;
	      if(xrr >= N1)
            xrr -= N1;
	      if(xl < 0)
            xl += N1;
	      if(xll < 0)
            xll += N1;
          ULONG iinl = index_3d(xl,y,z,N2,N3);
          ULONG iinr = index_3d(xr,y,z,N2,N3);
          ULONG iinll= index_3d(xll,y,z,N2,N3);
          ULONG iinrr= index_3d(xrr,y,z,N2,N3);
          ULONG iout = index_3d(x,y,z,N2,N3);
          out[iout] =(-in[iinrr]+8.0*in[iinr]-8.0*in[iinl]+in[iinll])*ideltaH;
//	      out[iout] = -static_cast<real_prec>((factor*((4.0/3.0)*(in[iinl]-in[iinr])-(1.0/6.0)*(in[iinll]-in[iinrr]))));
            }
      break;
  case(2):
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
      for(int x = 0; x < N1; x++)
    for(int y = 0; y < N2; y++)
      for(int z = 0; z < N3; z++)
	    {
          int yr = y + 1;
          int yl = y - 1;
          int yrr = y + 2;
          int yll = y - 2;
	      if(yr >= N2)
            yr -= N2;
	      if(yrr >= N2)
            yrr -= N2;
	      if(yl < 0)
            yl += N2;
	      if(yll < 0)
            yll += N2;
          ULONG iinl = index_3d(x,yl,z,N2,N3);
          ULONG iinr = index_3d(x,yr,z,N2,N3);
          ULONG iinll= index_3d(x,yll,z,N2,N3);
          ULONG iinrr= index_3d(x,yrr,z,N2,N3);
          ULONG iout = index_3d(x,y,z,N2,N3);
          out[iout] =(-in[iinrr]+8.0*in[iinr]-8.0*in[iinl]+in[iinll])*ideltaH;
//	      out[iout] = -static_cast<real_prec>((factor*((4.0/3.0)*(in[iinl]-in[iinr])-(1.0/6.0)*(in[iinll]-in[iinrr]))));
            }
        break;
     case(3):
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif
      for(int x = 0; x < N1; x++)
    for(int y = 0; y < N2; y++)
      for(int z = 0; z < N3; z++)
	    {
          int zr = z + 1;
          int zl = z - 1;
          int zrr = z + 2;
          int zll = z - 2;
	      if(zr >= N3)
		zr -= N3;
	      if(zrr >= N3)
		zrr -= N3;
	      if(zl < 0)
		zl += N3;
	      if(zll < 0)
		zll += N3;
          ULONG iinl = index_3d(x,y,zl,N2,N3);
          ULONG iinr = index_3d(x,y,zr,N2,N3);
          ULONG iinll= index_3d(x,y,zll,N2,N3);
          ULONG iinrr= index_3d(x,y,zrr,N2,N3);
          ULONG iout = index_3d(x,y,z,N2,N3);
          out[iout] =(-in[iinrr]+8.0*in[iinr]-8.0*in[iinl]+in[iinll])*ideltaH;
	    }
    break;
  }
 }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void low_pass_filter(ULONG Nft_HR, ULONG Nft_LR, int imas, bool correct, vector<real_prec>&HR_field, vector<real_prec>&LR_field, real_prec Lbox)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }

#ifdef _DOUBLE_PREC_
    complex_prec *HR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_HR*sizeof(real_prec));
#else
  complex_prec *HR_FOURIER= (complex_prec *)fftwf_malloc(2*NTT_HR*sizeof(real_prec));
#endif
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][REAL]=0;
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][IMAG]=0;
   // DO FFT of the high-res density field
   do_fftw_r2c(Nft_HR,HR_field,HR_FOURIER);
    // Get the MAS correction for the H-Res density field
   vector<real_prec> correction(Nft_HR,1.0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i = 1 ; i < Nft_HR; ++i )
    {
      int  coords= i<=Nft_HR/2? i: i-Nft_HR;
      real_prec xx=static_cast<real_prec>(coords)*M_PI/static_cast<real_prec>(Nft_HR);
      real_prec kernel= sin(xx)/static_cast<real_prec>(xx);
      correction[i]= true==correct ? pow(kernel,imas+1) : 1.0 ;
    }
#ifdef _DOUBLE_PREC_
  complex_prec *LR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_LR*sizeof(real_prec));
#else
   complex_prec *LR_FOURIER= (complex_prec *)fftwf_malloc(2*NTT_LR*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
  for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][IMAG]=0;
#endif
  vector<real_prec> kmodes_lr(Nft_LR,0.0);
  real_prec delta_k= 2*M_PI/Lbox;
#pragma omp parallel for
   for(int i = 0 ; i < Nft_LR; ++i )
   {
      int coords= i < Nft_LR/2+1? i: i-Nft_LR;
      kmodes_lr[i]=static_cast<real_prec>(coords)*delta_k;
   }
   // Get the FT of the Low-res density field
real_prec delta_box_hr= Lbox/static_cast<real_prec>(Nft_HR);
real_prec alpha=static_cast<double>(Nft_HR)/static_cast<double>(Nft_LR);
#ifdef _ABACUS_
real_prec delta_shift= 0;// Abacus uses lattice-like descriptions which do not require any shift at a low-pass filter
#else
real_prec delta_shift= 0.5*(alpha-1.0)*delta_box_hr; //Yu
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft_LR;++i)
    {
      ULONG i_hr = i < Nft_LR/2 +1 ? i : i+Nft_HR-Nft_LR;
      for(ULONG j=0;j< Nft_LR;++j)
        {
          ULONG j_hr = j < Nft_LR/2 +1 ? j : j+Nft_HR-Nft_LR;
          for(ULONG k=0;k<Nft_LR/2+1;++k)
            {
              ULONG k_hr=  k < Nft_LR/2 +1 ? k : k+Nft_HR-Nft_LR;
              ULONG index_lr=index_3d(i,   j,   k,   Nft_LR, Nft_LR/2+1);
              ULONG index_hr=index_3d(i_hr,j_hr,k_hr,Nft_HR ,Nft_HR/2+1);
              real_prec corr= 1.0;//correction[i_hr]*correction[j_hr]*correction[k_hr];
              real_prec new_real=HR_FOURIER[index_hr][REAL]/corr;
              real_prec new_imag=HR_FOURIER[index_hr][IMAG]/corr;
              real_prec shift=delta_shift*(kmodes_lr[i] + kmodes_lr[j] + kmodes_lr[k]);
              LR_FOURIER[index_lr][REAL] = (cos(shift)*new_real - sin(shift)*new_imag); //This follows the convention delta_2-> exp(iks)*delta_2
              LR_FOURIER[index_lr][IMAG] = (cos(shift)*new_imag + sin(shift)*new_real);
            }
        }
    }
  // Set Ny freq of the FT of the L-res density field to real:
  real_prec a=0, b=0;
  ULONG ii=0, jj=0, kk=0;
  ULONG ind=0;
  real_prec fac=1.0;
  ii=0; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=0; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=0; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=0; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=Nft_LR/2; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=Nft_LR/2; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=Nft_LR/2; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
  ii=Nft_LR/2; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
 // Transform to L-res density field.
  do_fftw_c2r(Nft_LR ,LR_FOURIER,LR_field);
#ifdef _DOUBLE_PREC_
  fftw_free(LR_FOURIER);
  fftw_free(HR_FOURIER);
#else
  fftwf_free(LR_FOURIER);
  fftwf_free(HR_FOURIER);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// Esta versión toma el campo, lo filtra, pero devuelve el campo en FOurier ya filtrado
void low_pass_filter(ULONG Nft_HR, ULONG Nft_LR, int imas, bool correct, vector<real_prec>&HR_field, complex_prec *LR_FOURIER, real_prec Lbox)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
  ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
  if(NTT_HR<NTT_LR)
    {
      cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
      exit(0);
    }
#ifdef _DOUBLE_PREC_
    complex_prec *HR_FOURIER= (complex_prec *)fftw_malloc(2*NTT_HR*sizeof(real_prec));
#else
  complex_prec *HR_FOURIER= (complex_prec *)fftwf_malloc(2*NTT_HR*sizeof(real_prec));
#endif
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][REAL]=0;
#pragma omp parallel for
   for(ULONG i=0;i<NTT_HR;++i)HR_FOURIER[i][IMAG]=0;
   // DO FFT of the high-res density field
   do_fftw_r2c(Nft_HR,HR_field,HR_FOURIER);
    // Get the MAS correction for the H-Res density field
   vector<real_prec> correction(Nft_HR,1.0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i = 1 ; i < Nft_HR; ++i )
    {
      int  coords= i<=Nft_HR/2? i: i-Nft_HR;
      real_prec xx=static_cast<real_prec>(coords)*M_PI/static_cast<real_prec>(Nft_HR);
      real_prec kernel= sin(xx)/static_cast<real_prec>(xx);
      correction[i]= true==correct ? pow(kernel,imas+1) : 1.0 ;
    }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][REAL]=0;
#ifdef _USE_OMP_
#pragma omp parallel for
  for(ULONG i=0;i<NTT_LR;++i)LR_FOURIER[i][IMAG]=0;
#endif
  vector<real_prec> kmodes_lr(Nft_LR,0.0);
  real_prec delta_k= 2*M_PI/Lbox;

#pragma omp parallel for
   for(int i = 0 ; i < Nft_LR; ++i )
   {
      int coords= i < Nft_LR/2+1? i: i-Nft_LR;
      kmodes_lr[i]=static_cast<real_prec>(coords)*delta_k;
   }
   // Get the FT of the Low-res density field
real_prec delta_box_hr= Lbox/static_cast<real_prec>(Nft_HR);
real_prec alpha=static_cast<double>(Nft_HR)/static_cast<double>(Nft_LR);

#ifdef _ABACUS_
real_prec delta_shift= 0;// Abacus uses lattice-like descriptions which do not require any shift at a low-pass filter
#else
real_prec delta_shift= 0.5*(alpha-1.0)*delta_box_hr; //Yu
#endif

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft_LR;++i)
    {
      ULONG i_hr = i < Nft_LR/2 +1 ? i : i+Nft_HR-Nft_LR;
      for(ULONG j=0;j< Nft_LR;++j)
        {
          ULONG j_hr = j < Nft_LR/2 +1 ? j : j+Nft_HR-Nft_LR;
          for(ULONG k=0;k<Nft_LR/2+1;++k)
            {
              ULONG k_hr=  k < Nft_LR/2 +1 ? k : k+Nft_HR-Nft_LR;
              ULONG index_lr=index_3d(i,   j,   k,   Nft_LR, Nft_LR/2+1);
              ULONG index_hr=index_3d(i_hr,j_hr,k_hr,Nft_HR ,Nft_HR/2+1);
              real_prec corr= 1.0;//correction[i_hr]*correction[j_hr]*correction[k_hr];
              real_prec new_real=HR_FOURIER[index_hr][REAL]/corr;
              real_prec new_imag=HR_FOURIER[index_hr][IMAG]/corr;
              real_prec shift=delta_shift*(kmodes_lr[i] + kmodes_lr[j] + kmodes_lr[k]);
              LR_FOURIER[index_lr][REAL] = (cos(shift)*new_real - sin(shift)*new_imag); //This follows the convention delta_2-> exp(iks)*delta_2
              LR_FOURIER[index_lr][IMAG] = (cos(shift)*new_imag + sin(shift)*new_real);
            }
        }
    }
  // Set Ny freq of the FT of the L-res density field to real:
  real_prec a=0, b=0;
  ULONG ii=0, jj=0, kk=0;
  ULONG ind=0;
  real_prec fac=1.0;

  ii=0; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=0; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=0; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=0; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=Nft_LR/2; kk=0;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;

  ii=Nft_LR/2; jj=Nft_LR/2; kk=Nft_LR/2;
  ind = index_3d(ii,jj,kk, Nft_LR, Nft_LR/2+1);
  a  =  LR_FOURIER[ind][REAL];
  b  =  fac*LR_FOURIER[ind][IMAG];
  LR_FOURIER[ind][REAL]=sqrt(a*a+b*b);
  LR_FOURIER[ind][IMAG]=0.0;
}
////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
/*
delta1 = Aexp(i phi1)
delta2 = Bexp(i phi2)
output delta2=A exp(i phi2), that is, conserves the phases but with a different amplitude
*/
void swap_amp_fourier(ULONG Nft, vector<real_prec>&in_ref, vector<real_prec>&out_field)
{
#ifdef _USE_OMP_
    int NTHREADS = _NTHREADS_;
    omp_set_num_threads(NTHREADS);
#endif
  ULONG NTT=static_cast<ULONG>(Nft*Nft*(Nft/2+1));
  ULONG NGRID=in_ref.size();
#ifdef _DOUBLE_PREC_
  complex_prec *in_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
  complex_prec *out_FOURIER= (complex_prec *)fftw_malloc(2*NTT*sizeof(real_prec));
#else
  complex_prec *in_FOURIER= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
  complex_prec *out_FOURIER= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
#endif

#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)in_FOURIER[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)in_FOURIER[i][IMAG]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)out_FOURIER[i][REAL]=0;
#pragma omp parallel for
  for(ULONG i=0;i<NTT;++i)out_FOURIER[i][IMAG]=0;
  do_fftw_r2c(Nft,in_ref,in_FOURIER);
  do_fftw_r2c(Nft,out_field,out_FOURIER);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NTT;++i)
   {
     real_prec Amp_in = sqrt(pow(in_FOURIER[i][REAL],2)+pow(in_FOURIER[i][IMAG],2));
     real_prec phi_o  = atan2(out_FOURIER[i][IMAG],out_FOURIER[i][REAL]);
     out_FOURIER[i][REAL] = Amp_in*cos(phi_o); 
     out_FOURIER[i][IMAG] = Amp_in*sin(phi_o);
   }
  do_fftw_c2r(Nft,out_FOURIER,out_field);
#ifdef _DOUBLE_PREC_
  fftw_free(in_FOURIER);
  fftw_free(out_FOURIER);
#else
  fftwf_free(in_FOURIER);
  fftwf_free(out_FOURIER);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void exchange_xy(int Nft,const vector<real_prec>&in, vector<real_prec>&out)
{
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
for(ULONG i=0;i<Nft;++i)
  for(ULONG j=0;j<Nft;++j)
    for(ULONG k=0;k<Nft;++k)
      out[index_3d(i,j,k,Nft,Nft)]=in[index_3d(j,i,k,Nft,Nft)];
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void exchange_xz(int Nft,const vector<real_prec>&in, vector<real_prec>&out)
{
 // this is alsme meant to be column (F) to major (c) order
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nft;++i)
      for(ULONG j=0;j<Nft;++j)
         for(ULONG k=0;k<Nft;++k)
             out[index_3d(i,j,k,Nft,Nft)]=in[index_3d(k,j,i,Nft,Nft)];
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_10d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr)
{
  return  static_cast<ULONG>(r) + static_cast<ULONG>(Nr*q) +  static_cast<ULONG>(Nr*Nq*p)+ static_cast<ULONG>(Nr*Nq*Np*o) + static_cast<ULONG>(Nr*Nq*Np*No*n)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*m)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*l)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*k)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_11d(ULONG i, ULONG j, ULONG k, ULONG l, ULONG m, ULONG n, ULONG o, ULONG p, ULONG q, ULONG r, ULONG s, ULONG Nj, ULONG Nk, ULONG Nl, ULONG Nm, ULONG Nn, ULONG No, ULONG Np, ULONG Nq, ULONG Nr, ULONG Ns)
{
 return  static_cast<ULONG>(s) + static_cast<ULONG>(Ns*r)+static_cast<ULONG>(Ns*Nq*q)+ static_cast<ULONG>(Ns*Nr*Nq*p)+ static_cast<ULONG>(Ns*Nr*Nq*Np*o) +static_cast<ULONG>(Ns*Nr*Nq*Np*No*n)
            +static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*m)+static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Ns*Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*k)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_11d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns)
{
 return  static_cast<ULONG>(s) + static_cast<ULONG>(Ns*r)+static_cast<ULONG>(Ns*Nq*q)+ static_cast<ULONG>(Ns*Nr*Nq*p)+ static_cast<ULONG>(Ns*Nr*Nq*Np*o) +static_cast<ULONG>(Ns*Nr*Nq*Np*No*n)
            +static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*m)+static_cast<ULONG>(Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Ns*Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*k)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_12d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int r, int s, int v, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq, int Nr, int Ns, int Nv)
{
 return static_cast<ULONG>(v) + static_cast<ULONG>(Nv*s) + static_cast<ULONG>(Nv*Ns*r)+static_cast<ULONG>(Nv*Ns*Nr*q)+static_cast<ULONG>(Nv*Ns*Nr*Nq*p)+ static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*o)+ static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No*n) +static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No*Nn*m)
            +static_cast<ULONG>(Nv*Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Nv*Ns*Nr*Nq)*static_cast<ULONG>(Np*No*Nn*Nm*Nl*k)+static_cast<ULONG>(Nv*Ns*Nr*Nq*Np*No)*static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nv*Ns*Nr*Nq*Np)*static_cast<ULONG>(No*Nn*Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_9d(int i, int j, int k, int l, int m, int n, int o, int p, int q, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np, int Nq)
{
  return  static_cast<ULONG>(q) + static_cast<ULONG>(Nq*p) +  static_cast<ULONG>(Nq*Np*o)+ static_cast<ULONG>(Nq*Np*No*n) + static_cast<ULONG>(Nq*Np*No*Nn*m)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*l)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*k)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*Nk*j)+static_cast<ULONG>(Nq*Np*No*Nn)*static_cast<ULONG>(Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_8d(int i, int j, int k, int l, int m, int n, int o, int p, int Nj, int Nk, int Nl, int Nm, int Nn, int No, int Np)
{
  return  static_cast<ULONG>(p) + static_cast<ULONG>(Np*o) +  static_cast<ULONG>(Np*No*n)+ static_cast<ULONG>(Np*No*Nn*m) + static_cast<ULONG>(Np*No*Nn*Nm*l)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*k)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*Nk*j)+static_cast<ULONG>(Np*No*Nn*Nm)*static_cast<ULONG>(Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_7d(int i, int j, int k, int l, int m, int n, int o, int Nj, int Nk, int Nl, int Nm, int Nn, int No)
{
  return  static_cast<ULONG>(o) + static_cast<ULONG>(No*n) +  static_cast<ULONG>(No*Nn*m)+ static_cast<ULONG>(No*Nn*Nm*l) + static_cast<ULONG>(No*Nn*Nm*Nl*k)+static_cast<ULONG>(No*Nn*Nm*Nl)*static_cast<ULONG>(Nk*j)+static_cast<ULONG>(No*Nn*Nm*Nl)*static_cast<ULONG>(Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_6d(int i, int j, int k, int l, int m, int n, int Nj, int Nk, int Nl, int Nm, int Nn)
{
  return  static_cast<ULONG>(n) + static_cast<ULONG>(Nn*m) +  static_cast<ULONG>(Nn*Nm*l)+ static_cast<ULONG>(Nn*Nm*Nl*k) + static_cast<ULONG>(Nn*Nm*Nl*Nk*j)+static_cast<ULONG>(Nn*Nm*Nl)*static_cast<ULONG>(Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_5d(int i, int j, int k, int l, int m, int Nj, int Nk, int Nl, int Nm)
{
  return static_cast<ULONG>(m) + static_cast<ULONG>(Nm*l)+ static_cast<ULONG>(Nm*Nl*k)+static_cast<ULONG>(Nm*Nl*Nk*j) +  static_cast<ULONG>(Nm*Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

ULONG index_4d(int i, int j, int k, int l, int Nj, int Nk, int Nl)
{
  return static_cast<ULONG>(l)+ static_cast<ULONG>(Nl*k)+static_cast<ULONG>(Nl*Nk*j)+ static_cast<ULONG>(Nl*Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

ULONG index_3d(ULONG i, ULONG j, ULONG k, int Nj, int Nk)
{ // if k denots the z coordinate, then this is row majort (c) order
  return k+ static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_3d(ULONG i, ULONG j, ULONG k, ULONG Nj, ULONG Nk)
{ // if k denots the z coordinate, then this is row majort (c) order
  return k+ static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG index_2d(ULONG i, ULONG j, ULONG Nj)
{
  return j+Nj*i;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// This funciton can be applied to both C like (row-major) or Fortran like (columnm major) order
void index2coords(ULONG index, ULONG N,  ULONG  &XG, ULONG &YG, ULONG &ZG )
{
//  see https://math.stackexchange.com/questions/3758576/how-to-convert-from-a-flattened-3d-index-to-a-set-of-coordinates

  ZG=index % N;  // Fast varying coordinate. Can be z(c++) or x (Fortran)
  index = static_cast<ULONG>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  YG=index % N;
  XG = static_cast<ULONG>(static_cast<real_prec>(index)/static_cast<real_prec>(N));  // Slow varying coordinate. Can be x(c++) or z(Fortran)
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
ULONG grid_ID(s_params_box_mas *params, const real_prec &x, const real_prec &y, const real_prec &z)
{
  ULONG i = static_cast<ULONG>(floor((x-params->min1)/params->d1)+ BIN_SHIFT); // indices of the cell of the particle
  ULONG j = static_cast<ULONG>(floor((y-params->min2)/params->d2)+ BIN_SHIFT);
  ULONG k = static_cast<ULONG>(floor((z-params->min3)/params->d3)+ BIN_SHIFT);

  if(i==params->Nft)
    i--;
  if(j==params->Nft)
    j--;
  if(k==params->Nft)
    k--;
  
  i = static_cast<ULONG>(fmod(real_prec(i),real_prec(params->Nft)));
  j = static_cast<ULONG>(fmod(real_prec(j),real_prec(params->Nft)));
  k = static_cast<ULONG>(fmod(real_prec(k),real_prec(params->Nft)));
  return index_3d(i,j,k,params->Nft,params->Nft);
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_r2c(ULONG Nft, vector<real_prec>&in, complex_prec *out)
{

#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

#ifdef DOUBLE_PREC
  int *n=(int *)fftw_malloc(ic_rank*sizeof(int));
#else
  int *n=(int *)fftwf_malloc(ic_rank*sizeof(int));
#endif
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef SINGLE_PREC
#ifdef _USE_OMP_
  fftwf_init_threads();
  fftwf_plan_with_nthreads(NTHREADS);
#endif
  fftwf_plan plan_r2c=fftwf_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftwf_execute(plan_r2c);
  fftwf_destroy_plan(plan_r2c);
  fftwf_free(n);
#endif
#ifdef DOUBLE_PREC
#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(NTHREADS);
#endif
  fftw_plan plan_r2c=fftw_plan_dft_r2c(ic_rank,n,&in[0],out,FFTW_ESTIMATE);
  fftw_execute(plan_r2c);
  fftw_destroy_plan(plan_r2c);
  fftw_free(n);
#endif
}
///////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_1d(ULONG Nft, complex_prec *in, complex_prec *out, int dir)
{
#ifdef SINGLE_PREC
  fftwf_init_threads();
  fftwf_plan_with_nthreads(_NTHREADS_);
  fftwf_plan plan_p=fftwf_plan_dft_1d(Nft,in,out,dir,FFTW_ESTIMATE);
  fftwf_execute(plan_p);
  fftwf_destroy_plan(plan_p);
#elif defined DOUBLE_PREC
  fftw_init_threads();
  fftw_plan_with_nthreads(_NTHREADS_);
  fftw_plan plan_p=fftw_plan_dft_1d(Nft, in ,out, dir,FFTW_ESTIMATE);
  fftw_execute(plan_p);
  fftw_destroy_plan(plan_p);
#endif
  if(FFTW_BACKWARD==dir)
  {
      for(ULONG i=0;i<Nft;i++)  // normalize the c2r transform
          out[i][REAL]/=static_cast<double>(Nft);
      for(ULONG i=0;i<Nft;i++)  // normalize the c2r transform
          out[i][IMAG]/=static_cast<double>(Nft);
    }
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_1d_r2c(ULONG Nft, vector<real_prec>&in, complex_prec *out)
{
#ifdef SINGLE_PREC
  fftwf_init_threads();
  fftwf_plan_with_nthreads(_NTHREADS_);
  fftwf_plan plan_p=fftwf_plan_dft_r2c_1d(Nft, &in[0] ,out,FFTW_ESTIMATE);
  fftwf_execute(plan_p);
//  fftwf_destroy_plan(plan_p);
#elif defined DOUBLE_PREC
  fftw_init_threads();
  fftw_plan_with_nthreads(_NTHREADS_);
  fftw_plan plan_p=fftw_plan_dft_r2c_1d(Nft, &in[0] ,out,FFTW_ESTIMATE);
  fftw_execute(plan_p);
  fftw_destroy_plan(plan_p);
#endif
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_1d_c2r(ULONG Nft, complex_prec *in, vector<real_prec>&out)
{
#ifdef SINGLE_PREC
  fftwf_init_threads();
  fftwf_plan_with_nthreads(_NTHREADS_);
  fftwf_plan plan_p=fftwf_plan_dft_c2r_1d(Nft, in ,&out[0],FFTW_ESTIMATE);
  fftwf_execute(plan_p);
//  fftwf_destroy_plan(plan_p);
#elif defined DOUBLE_PREC
  fftw_init_threads();
  fftw_plan_with_nthreads(_NTHREADS_);
  fftw_plan plan_p=fftw_plan_dft_c2r_1d(Nft, in ,&out[0],FFTW_ESTIMATE);
  fftw_execute(plan_p);
  fftw_destroy_plan(plan_p);
#endif
      for(ULONG i=0;i<Nft;i++)  // normalize the c2r transform
          out[i]/=static_cast<double>(Nft);
  }
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_c2r(int Nft, complex_prec *in, vector<real_prec>&out)
{
  ULONG NGRID = static_cast<ULONG>(Nft*Nft*Nft);
#ifdef DOUBLE_PREC
  int *n=(int *)fftw_malloc(ic_rank*sizeof(int));
#else
  int *n=(int *)fftwf_malloc(ic_rank*sizeof(int));
#endif
  for(int i=0;i<ic_rank;++i)n[i]=Nft;
#ifdef SINGLE_PREC
#ifdef _USE_OMP_
  fftwf_init_threads();
  fftwf_plan_with_nthreads(_NTHREADS_);
#endif
  fftwf_plan plan_c2r=fftwf_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftwf_execute(plan_c2r);
  fftwf_destroy_plan(plan_c2r);
  fftwf_free(n);
#endif
#ifdef DOUBLE_PREC
#ifdef _USE_OMP_
  fftw_init_threads();
  fftw_plan_with_nthreads(omp_get_max_threads());
#endif
  fftw_plan plan_c2r=fftw_plan_dft_c2r(ic_rank,n,in,&out[0], FFTW_ESTIMATE);
  fftw_execute(plan_c2r);
  fftw_destroy_plan(plan_c2r);
  fftw_free(n);
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<NGRID;i++)  // normalize the c2r transform
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      out[i]/=static_cast<double>(NGRID);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void do_fftw_3d(ULONG Nft, bool direction, complex_prec *in, complex_prec *out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG factor=Nft*Nft*Nft;
  int Ni1=static_cast<int>(Nft);
  int Ni2=static_cast<int>(Nft);
  int Ni3=static_cast<int>(Nft);
  if (direction == true)//from conf space to fourier
    {
#ifdef SINGLE_PREC
      fftwf_plan fftp;
      fftp = fftwf_plan_dft_3d(Ni1,Ni2,Ni3,in,out,FORWARD,FFTW_OPTION);	
      fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_plan fftp;
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,FORWARD,FFTW_OPTION);	
      fftw_execute(fftp);
#endif

    }
  else // from fourier to conf space
    {
#ifdef SINGLE_PREC
      fftwf_plan fftp;
      fftp = fftwf_plan_dft_3d(Ni1,Ni2,Ni3,in,out,BACKWARD,FFTW_OPTION);
      fftwf_execute(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_plan fftp;
      fftp = fftw_plan_dft_3d(Ni1,Ni2,Ni3,in,out,BACKWARD,FFTW_OPTION);	
      fftw_execute(fftp);
#endif
      for(ULONG i=0;i<factor;i++)
	{// normalize the c2r transform
          out[i][REAL]/=static_cast<real_prec>(factor);
          out[i][IMAG]/=static_cast<real_prec>(factor);
	}
#ifdef SINGLE_PREC
      fftwf_destroy_plan(fftp);
#endif
#ifdef DOUBLE_PREC
      fftw_destroy_plan(fftp);
#endif
    }
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void get_cumulative(const vector<real_prec> &dist, vector<real_prec> &cumu, unsigned long &NTOT_h)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  NTOT_h=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:NTOT_h)
#endif
  for(int i=0;i<dist.size();++i)
    NTOT_h+=dist[i];

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
   for(int i=0;i<dist.size();++i)
    for(int j=0;j<=i;++j)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
        cumu[i]+=dist[j];

#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i=0;i<dist.size();++i)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
      cumu[i]/=static_cast<double>(NTOT_h);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void convert_ngp_to_cic(vector<real_prec>&in, vector<real_prec>&out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N=in.size();
  ULONG N1=static_cast<ULONG>(pow(N,1./3.))+1;
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);
#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }

  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);
  vector<real_prec> coords(N1,0);
  vector<real_prec> correction(N1,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(int i = 0 ; i < N1; ++i )
    {
      coords[i]=static_cast<real_prec>( i<=N1/2? i: i-(N1));
      real_prec xx=coords[i]*M_PI/static_cast<real_prec>(N1);
      correction[i]= i==0 ? 1.0 : sin(xx)/static_cast<real_prec>(xx) ;
    }
  real_prec we=0;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:we)
#endif
  for(ULONG i=0;i< N1; i++)
    for(ULONG j=0;j< N1; j++)
      for(ULONG k=0;k< N1/2+1; k++)
    {
          real_prec cor=correction[i]*correction[j]*correction[k];
      ULONG ind=index_3d(i,j,k, N1, N1/2+1);
      AUX[ind][REAL]*=cor;
      AUX[ind][IMAG]*=cor;
      we+=cor;
    }

  do_fftw_c2r(N1,AUX,out);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N;i++)
    out[i]=(out[i]<0 ? 0: out[i]);//*(static_cast<real_prec>(N)/static_cast<real_prec>(2.0*we));

#ifdef DOUBLE_PREC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
real_prec MAS_CIC_public(real_prec x)
{
  //////////////////////////////////////////////////////////
  // CIC Mass assignment scheme.
  //////////////////////////////////////////////////////////
  x=fabs(x);
  if(x<1)
    return 1.-x;
  else
    return 0.0;
}
////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
#ifdef OMPPARRANGAR
real_prec GR_NUM(gsl_rng ** SEED, real_prec sigma ,int GR_METHOD, int jthread)
{
#else
real_prec GR_NUM(gsl_rng * SEED, real_prec sigma ,int GR_METHOD)
{
#endif

  real_prec val=0.;
  switch (GR_METHOD)
    {
    case 0:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian(SEED, sigma));
#endif
      }
      break;

    case 1:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian_ziggurat(SEED, sigma));
#endif
      }
      break;

    case 2:
      {
#ifdef OMPPARRANGAR
        val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED[jthread], sigma));
#else
        val=static_cast<real_prec>(gsl_ran_gaussian_ratio_method(SEED, sigma));
#endif
      }
      break;

    }
  return(val);
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
#ifdef OMPPARRANGAR
 void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng ** seed)
{
#else
  void create_GARFIELDR(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed)
{
#endif

  #ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif

  //  FileOutput File;
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  //  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
  ULONG Nhalf = static_cast<ULONG>(N1*N2*(N3/2+1));

  #ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\tComputing white noise "<<RESET<<endl;
#endif
#ifdef DOUBLE_PREC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif


  for (ULONG i=0;i<N;i++)
    delta[i] = static_cast<real_prec>(GR_NUM(seed,num_1,0));

  do_fftw_r2c(N1,delta,AUX);

#ifdef _FULL_VERBOSE_
    cout<<YELLOW<<"\tGoing to the Fourier grid "<<RESET<<endl;
#endif

#ifdef _USE_OMP_
#pragma omp parallel for// take care!!!
#endif
    for (ULONG i=0 ; i<N1;i++)
      for (ULONG j=0 ; j<N2;j++)
	for (ULONG k=0 ; k<=N3/2;k++)
	  {
	    ULONG ihalf= index_3d(i,j,k,N2,N3/2+1);
	    ULONG iind = index_3d(i,j,k,N2,N3);
	    real_prec mass=0.0;
#ifdef FOURIER_DEF_1
	    mass=Power[iind];//testing
#endif
#ifdef FOURIER_DEF_2
            mass=Power[iind]/static_cast<real_prec>(N);
#endif
            real_prec sigma = sqrt(mass);
            AUX[ihalf][REAL]*=sigma;
            AUX[ihalf][IMAG]*=sigma;
/*
            // usa esta si queremos hacer gausr af en F space,m en ves de crear el delta y tranformarlo
            real_prec phase =2.*M_PI*gsl_rng_uniform(seed);
            real_prec sigma=gsl_ran_rayleigh(seed, sqrt(0.5*Power[iind]));
            AUX[ihalf][REAL]=sigma*cos(phase);
            AUX[ihalf][IMAG]=-sigma*sin(phase);
*/
            }

    do_fftw_c2r(N1,AUX,delta); 
#ifdef DOUBLE_PREC
    fftw_free(AUX);
#else
    fftwf_free(AUX);
#endif
  }
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
void create_GARFIELD_FIXED_AMP(ULONG N1,ULONG N2,ULONG N3,  vector<real_prec> &delta, const vector<real_prec> &Power, gsl_rng * seed)
 {
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
   ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
   ULONG Nhalf = static_cast<ULONG>(N1*N2*(N3/2+1));
   cout<<YELLOW<<"Computing white noise "<<N2<<RESET<<endl;

#ifdef DOUBLE_PREC
   complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
   complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#ifdef _FULL_VERBOSE_
   cout<<YELLOW<<"Going to the Fourier grid "<<RESET<<endl;
#endif
 #ifdef _USE_OMP_
 #pragma omp parallel for// take care!!!
 #endif
     for (ULONG i=0 ; i<N1;i++)
       for (ULONG j=0 ; j<N2;j++)
         for (ULONG k=0 ; k<=N3/2;k++)
           {
             ULONG iind = index_3d(i,j,k,N2,N3);
             real_prec sigma = sqrt(Power[iind]); // In order to generate WN set sigma=1
             real_prec phase =2.*M_PI*gsl_rng_uniform(seed);

             ULONG ihalf= index_3d(i,j,k,N2,N3/2+1);
             AUX[ihalf][REAL]=sigma*cos(phase);
             AUX[ihalf][IMAG]=sigma*sin(phase);
           }


     do_fftw_c2r(N1,AUX,delta);
     real_prec meanWN=get_mean(delta);
     cout<<YELLOW<<"Mean WN = "<<CYAN<<meanWN<<RESET<<endl;
     real_prec sigma2D=get_var(meanWN, delta);
     cout<<YELLOW<<"Sigma_corr WN = "<<CYAN<<sqrt(sigma2D)<<RESET<<endl;
#ifdef DOUBLE_PREC
      fftw_free(AUX);
#else
     fftwf_free(AUX);
#endif
   }
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
 void kernelcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3, real_prec smol, int filtertype, string output_dir)
 {
#ifdef _USE_OMP_
   int NTHREADS=_NTHREADS_;
   omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
   cout<<YELLOW<<"\tComputing kernel at "<<__PRETTY_FUNCTION__<<RESET<<endl;
#endif
   string fnameR=output_dir+"kernel"+to_string(N1)+"V"+to_string(L1)+"r"+to_string(smol);
   ifstream inStream;
   string ffnn=fnameR+".dat";
   inStream.open(ffnn.data());
   if (inStream.is_open() == false )
    {
      bool gauss=false;
      bool errfunc=false;
      bool tophat=false;
       
      switch (filtertype)
     	 {
 	      case 1:
	      gauss=true;
	      break;
	      case 2:
        tophat=true;
	      break;
	      case 3:
	      errfunc=true;
	      break;
  	   }

      ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
      ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
      vector<real_prec> out(Nhalf,0);

#ifdef DOUBLE_PREC
      complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
      complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
     real_prec asmth=num_1;
     real_prec u; 
     real_prec rS=smol;
     real_prec rS2=rS*rS;
     real_prec kcut=smol;//2.*M_PI/rS;
     real_prec sigma=static_cast<real_prec>(.3);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for (ULONG i=0;i<N1;i++)
	 for (ULONG j=0;j<N2;j++)
	   for (ULONG k=0;k<N3/2+1;k++)
	     {

	       ULONG ii=index_3d(i,j,k,N2,N3/2+1);
	       real_prec k2=k_squared(i,j,k,L1,L2,L3,N1,N2,N3/2+1);
	       
	       if (tophat==true)
		 {
		   u = sqrt(k2);
		   
		   if (u>kcut)
		     AUX[ii][REAL]=0.0;
		   else
		     AUX[ii][REAL]=num_1;
		 }
	       
	       if (errfunc==true)
		 {
		   u = static_cast<real_prec>((sqrt(k2)-kcut)/(sqrt(2.)*sigma));
		   real_prec fac = static_cast<real_prec>(erfc(u));
		   AUX[ii][REAL]=fac;
		 }
	       
	       if (gauss==true)
		 AUX[ii][REAL]=static_cast<real_prec>(exp(-k2*rS2/2.));
	       
	       
	       AUX[ii][IMAG]=0.0;

	     }
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
       for(ULONG i=0;i<Nhalf;i++)
   	   out[i]=AUX[i][REAL];
       vector<real_prec> aux(N,0);
       do_fftw_c2r(N1, AUX, aux);
	
        real_prec wtotD=0.;
#ifdef _USE_OMP_
#pragma omp parallel for reduction(+:wtotD)
#endif
        for(ULONG i=0;i<N;i++)
          wtotD+=static_cast<real_prec>(aux[i]);

	aux.clear();
	aux.shrink_to_fit();

	// Normalize kernel in Fourier space
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
        for(ULONG i=0;i<Nhalf;i++)
#ifdef _USE_OMP_
#pragma omp atomic update
#endif
          out[i]/=static_cast<real_prec>(wtotD);

	// este out debe ser en Fourier espace
       dump_scalar(out,N1,N2,N3/2+1,0,fnameR);
#ifdef DOUBLE_PREC
       fftw_free(AUX);
#else
       fftwf_free(AUX);
#endif
       }
 }
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
void calc_twolptterm(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, vector<real_prec>&phiv, vector<real_prec> &m2v)
{  
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);
  string fnamef;
#ifndef COMPATCTWOLPT
  #ifdef DOUBLE_PREC
  complex_prec *philv= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *philv= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#pragma omp parallel for
  for(ULONG i=0;i<Nhalf;i++)
    {
      philv[i][REAL]=0.0;
      philv[i][IMAG]=0.0;
    }
   do_fftw_r2c(N1,phiv,philv);
  vector<real_prec> LapPhivx(N,0), LapPhivy(N,0), LapPhivz(N,0);
  vector<real_prec> LapPhivxy(N,0), LapPhivxz(N,0), LapPhivyz(N,0);
 #ifdef  GFFT
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivx,1,1);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivy,2,2);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivz,3,3);
 calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivxy,1,2);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivxz,1,3);
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,LapPhivyz,2,3);
#endif
#if defined(GFINDIFF) || defined (GFFT)
  vector<real_prec> dummy(N,0);
#endif
#ifdef GFINDIFF  
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivx,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivxy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivxz,3);
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivyz,3);
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,LapPhivz,3);
#endif
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)
     m2v[i]=LapPhivx[i]*LapPhivy[i]-LapPhivxy[i]*LapPhivxy[i]+LapPhivx[i]*LapPhivz[i]-LapPhivxz[i]*LapPhivxz[i]+LapPhivy[i]*LapPhivz[i]-LapPhivyz[i]*LapPhivyz[i];
#else 
#ifdef GFFT
  //  LapPhivx[i]*LapPhivy[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,m2v,1,1);    
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,2,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,0,fnamef);	
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
   // -LapPhivxy[i]*LapPhivxy[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,1,2);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
   // +LapPhivx[i]*LapPhivz[i]
  {
    string fname=string("aux1");
    get_scalar(fname,dummy,N1,N2,N3);
  }
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,m2v,3,3);
  fnamef=string("aux3");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
   // -LapPhivxz[i]*LapPhivxz[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,1,3);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
   // +LapPhivy[i]*LapPhivz[i]
  {
    string fname=string("aux2");
    get_scalar(fname,dummy,N1,N2,N3);
  }
   {
    string fname=string("aux3");
    get_scalar(fname,m2v,N1,N2,N3);
  }
 #pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
   // -LapPhivyz[i]*LapPhivyz[i]
  calc_LapPhiv(N1,N2,N3,L1,L2,L3,philv,dummy,2,3);
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
#endif  
 #ifdef GFINDIFF
  //  LapPhivx[i]*LapPhivy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,0,fnamef);	
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // -LapPhivxy[i]*LapPhivxy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  
  // +LapPhivx[i]*LapPhivz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux3");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  {
    string fname=string("aux1");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=dummy[i];
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // -LapPhivxz[i]*LapPhivxz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  {
    string fname=string("aux");
    get_scalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	

  // +LapPhivy[i]*LapPhivz[i]
  {
    string fname=string("aux2");
    getscalar(fname,dummy,N1,N2,N3);
  }
  {
    string fname=string("aux3");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,0,fnamef);	
  // -LapPhivyz[i]*LapPhivyz[i]  
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i]*dummy[i];
#endif 
#endif



#ifdef DOUBLE_PREC
       fftw_free(philv);
#else
       fftwf_free(philv);
#endif
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
void calc_LapPhiv(ULONG N1,ULONG N2,ULONG N3,real_prec L1,real_prec L2,real_prec L3, complex_prec *philv,vector<real_prec>&LapPhiv,int index1,int index2)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);	
#ifdef DOUBLE_PREC
  complex_prec *LapPhivl= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *LapPhivl= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
  real_prec k1=0.;
  real_prec k2=0.;
  real_prec deltak=2.*M_PI/L1;
 vector<real_prec> coords(N1,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N1 ;++i)
    coords[i]=deltak*(i<=N1/2? static_cast<real_prec>(i): -static_cast<real_prec>(N1-i));
#ifdef _USE_OMP_
#pragma omp parallel for collapse(3)
#endif  
  for (ULONG i=0;i<N1;i++)
    for (ULONG j=0;j<N2;j++)
      for (ULONG k=0;k<N3/2+1;k++)
	{	  
      ULONG ind=index_3d(i,j,k, N1, N1/2+1);
	  switch (index1)
	    {	      
	    case 1:
          k1=coords[i];
	      break;
	    case 2:
          k1=coords[j];
	      break;
	    case 3:
          k1=coords[k];
	      break;
	    }

	  switch (index2)
	    {	      
	    case 1:
          k2=coords[i];
	      break;
	    case 2:
          k2=coords[j];
	      break;
	    case 3:
          k2=coords[k];
	      break;
	    }
      LapPhivl[ind][REAL]=-k1*k2*(philv[ind][REAL]);
      LapPhivl[ind][IMAG]=-k1*k2*(philv[ind][IMAG]);
	}
  do_fftw_c2r(N1,LapPhivl,LapPhiv);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 void calc_LapPhiv(ULONG N1,real_prec L1, complex_prec *Delta,vector<real_prec>&LapPhiv,int index1,int index2)
{
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N1)*static_cast<ULONG>(N1/2+1);
#ifdef DOUBLE_PREC
  complex_prec *LapPhivl= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *LapPhivl= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
  real_prec k1=0.;
  real_prec k2=0.;
  real_prec deltak=2.*M_PI/L1;
 vector<real_prec> coords(N1,0);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N1 ;++i)
    coords[i]=deltak*(i<=N1/2? static_cast<real_prec>(i): -static_cast<real_prec>(N1-i));
  ULONG i, j, k;
#ifdef _USE_OMPa_
#pragma omp parallel private(i,j,k, ind) //esta esta dando problema, pero es la que menos se demora, dejala comentada
   {
#pragma omp parallel for
#endif
  for (i=0;i<N1;i++)
    for (j=0;j<N1;j++)
      for (k=0;k<N1/2+1;k++)
      {
        ULONG ind=index_3d(i,j,k, N1, N1/2+1);
        real_prec kv2=pow(coords[i],2)+pow(coords[j],2)+pow(coords[k],2);
        switch (index1)
         {
          case 1:
            k1=coords[i];
          break;
          case 2:
            k1=coords[j];
          break;
          case 3:
            k1=coords[k];
          break;
        }
       switch (index2)
        {
        case 1:
          k2=coords[i];
        break;
        case 2:
          k2=coords[j];
        break;
        case 3:
          k2=coords[k];
        break;
        }
      real_prec term=0;
      if(kv2>0)
           term=(k1*k2)/static_cast<real_prec>(kv2);
       LapPhivl[ind][REAL]=term*Delta[ind][REAL];
       LapPhivl[ind][IMAG]=term*Delta[ind][IMAG];
     }
#ifdef _USE_OMPa_
 }
#endif
  do_fftw_c2r(N1,LapPhivl,LapPhiv);
#ifdef DOUBLE_PREC
  fftw_free(LapPhivl);
#else
  fftwf_free(LapPhivl);
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
 void calc_curlcomp(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec> &m2v, int comp)
{  
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
 int pcomp, mcomp;
  switch(comp)
    {
    case 1:
      {
	pcomp=2;
	mcomp=3;
	break;
      }
    case 2:
      {
	pcomp=1;
	mcomp=3;
	break;
      }
    case 3:
      {
	pcomp=1;
	mcomp=2;
	break;
      }
    }
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  string fnamef;     
  vector<real_prec> dummy(N,0);
#define num_0_2_5 static_cast<real_prec>(0.25)
#ifdef GFINDIFF  
  //LapPhiv1l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,L1,L2,L3,0,fnamef);	
  //LapPhiv2pp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,pcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv1l[i]*LapPhiv2pp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhivl[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2mm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,mcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv1l[i]*LapPhiv2mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  {
    string fname=string("aux");
    getscalar(fname,dummy,N1,N2,N3);
  }
  //1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("auxs");
  //LapPhiv2l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1pp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,pcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv2l[i]*LapPhiv1pp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2l[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1mm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,mcomp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv2l[i]*LapPhiv1mm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }

#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("auxs");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
 
  /////

  //LapPhiv1p[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,pcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2lp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv1p[i]*LapPhiv2lp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv1m[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,mcomp);
  fnamef=string("aux1");
  // dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2lm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv1m[i]*LapPhiv2lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("auxs");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2p[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,pcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1lp[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,pcomp);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }

  //1/4*LapPhiv2p[i]*LapPhiv1lp[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	

  //LapPhiv2m[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,mcomp);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1lm[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,comp);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,mcomp);
 
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //1/4*LapPhiv2m[i]*LapPhiv1lm[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_2_5*dummy[i];

  {
    string fname=string("auxs");
    get_scalar(fname,dummy,N1,N2,N3);
  }
 
  // 1/4*LapPhiv1l[i]*LapPhiv2pp[i]+1/4*LapPhiv1l[i]*LapPhiv2mm[i]
  //-1/4*LapPhiv2l[i]*LapPhiv1pp[i]+1/4*LapPhiv2l[i]*LapPhiv1mm[i]
  //-1/4*LapPhiv1l[i]*LapPhiv2lp[i]+1/4*LapPhiv1l[i]*LapPhiv2lm[i]
  //+1/4*LapPhiv2l[i]*LapPhiv1lp[i]+1/4*LapPhiv2l[i]*LapPhiv1lm[i]
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<N;i++)
    m2v[i]+=dummy[i];

#endif
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
 void calc_mu2term(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, const vector<real_prec> &phiv, vector<real_prec> phiv2, vector<real_prec>&m2v)
{  
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  string fnamef;     
  vector<real_prec> dummy(N);
#ifdef GFINDIFF  
  //LapPhiv1xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv1xx[i]*LapPhiv2yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,2);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,1);
   {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
   //0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
   /////
  //LapPhiv1xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,1);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,3);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv1xx[i]*LapPhiv2zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("auxb");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2xx[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,1);
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
 {
    string fname=string("auxb");
    get_scalar(fname,dummy,N1,N2,N3);
  }
  //0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];  
  {
    string fname=string("aux");
    get_scalar(fname,dummy,N1,N2,N3);
  }
   // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  ////
  //LapPhiv1yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,2);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,3);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	

  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv1yy[i]*LapPhiv2zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
  fnamef=string("auxb");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv2yy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,m2v,3);
  fnamef=string("aux1");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
  //LapPhiv1zz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux2");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  {
    string fname=string("aux1");
    get_scalar(fname,m2v,N1,N2,N3);
  }
  //0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]*=num_0_5*dummy[i];
 {
    string fname=string("auxb");
    getscalar(fname,dummy,N1,N2,N3);
  }
  //0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];  
  {
    string fname=string("aux");//!!!
    getscalar(fname,dummy,N1,N2,N3);
  }
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]+=dummy[i];
  fnamef=string("aux");
  // LapPhiv1xy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2xy[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,2);

  {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }

 // LapPhiv1xy[i]*LapPhiv2xy[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
   // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
    ////
  // LapPhiv1xz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2xz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,1);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
   {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }
 // LapPhiv1xz[i]*LapPhiv2xz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
  // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
  //-LapPhiv1xz[i]*LapPhiv2xz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
  fnamef=string("aux");
  dump_scalar(m2v,N1,N2,N3,9,fnamef);	
    ////
  // LapPhiv1yz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
  fnamef=string("aux1");
  dump_scalar(dummy,N1,N2,N3,9,fnamef);	
  // LapPhiv2yz[i]
  gradfindif(N1,N2,N3,L1,L2,L3,phiv2,m2v,2);
  gradfindif(N1,N2,N3,L1,L2,L3,m2v,dummy,3);
   {
    string fname=string("aux1");
    getscalar(fname,m2v,N1,N2,N3);
  }
 // LapPhiv1yz[i]*LapPhiv2yz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy[i]*=m2v[i];
  {
    string fname=string("aux");
    getscalar(fname,m2v,N1,N2,N3);
  }
   // 0.5*LapPhiv1xx[i]*LapPhiv2yy[i]+0.5*LapPhiv2xx[i]*LapPhiv1yy[i]
  //+0.5*LapPhiv1xx[i]*LapPhiv2zz[i]+0.5*LapPhiv2xx[i]*LapPhiv1zz[i]
  //+0.5*LapPhiv1yy[i]*LapPhiv2zz[i]+0.5*LapPhiv2yy[i]*LapPhiv1zz[i]
  //-LapPhiv1xy[i]*LapPhiv2xy[i] 
  //-LapPhiv1xz[i]*LapPhiv2xz[i] 
  //-LapPhiv1yz[i]*LapPhiv2yz[i] 
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    m2v[i]-=dummy[i];
#endif 
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
 void calc_Det(ULONG N1, ULONG N2, ULONG N3, real_prec L1, real_prec L2, real_prec L3, const vector<real_prec>&in, vector<real_prec> &out)
 {  
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
  ULONG N=N1*N2*N3;
  string fname;
  vector<real_prec> dummy(N,0), dummy2(N,0);
  //phi,11
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,1);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
  //phi,22
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,2);
  //phi,33
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);
  //phi,22*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,22*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  fname=string("auxs");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
   //phi,23
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);
   //phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];
   fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
   fname=string("auxs");
  get_scalar(fname,out,N1,N2,N3);
  //phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    out[i]-=dummy2[i];
  //  dump_scalar(out,N1,N2,N3,9,fname);
   ////
   //phi,12
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,2);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
   //phi,12*phi,12
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];
   //phi,33
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,3);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);
   //phi,12*phi,12*phi,33
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  fname=string("aux2");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
 //phi,23
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);
  //phi,13
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,3);
  //phi,23*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i];
  fname=string("aux1");
  get_scalar(fname,out,N1,N2,N3);
  //phi,12*phi,23*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=num_2*out[i];//this term appears twice in the determinant
  fname=string("aux2");
  get_scalar(fname,out,N1,N2,N3);
  //-phi,12*phi,12*phi,33+phi,12*phi,23*phi,13//sign must be changed below
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]-=out[i];
  fname=string("auxs");
  get_scalar(fname,dummy2,N1,N2,N3);
  //    phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
  //-(2*phi,12*phi,12*phi,33-phi,12*phi,23*phi,13)
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]+=out[i];//here comes a plus sign see above
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
  ////
  //phi,13
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,1);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,dummy2,3);
  fname=string("aux1");
  //  dump_scalar(dummy2,N1,N2,N3,9,fname);
  //phi,13*phi,13
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=dummy2[i];
  //phi,22
  gradfindif(N1,N2,N3,L1,L2,L3,in,dummy,2);
  gradfindif(N1,N2,N3,L1,L2,L3,dummy,out,2);
  //phi,13*phi,13*phi,22
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    dummy2[i]*=out[i]; 
  fname=string("auxs");
  get_scalar(fname,out,N1,N2,N3);
  //  phi,11*phi,22*phi,33-phi,11*phi,23*phi,23
  //-(phi,12*phi,12*phi,33-phi,12*phi,23*phi,13)
#pragma omp parallel for
  for(ULONG i=0;i<N;i++)    
    out[i]-=dummy2[i];
}
#endif
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
#ifdef _USE_LPT_
 void convcomp(real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3,ULONG N1, ULONG N2, ULONG N3,  vector<real_prec>&in, vector<real_prec> &out, int filtertype,real_prec smol, string file_kernel)
{
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"\tConvolution"<<RESET<<endl;
#endif
  ULONG N=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3);
  ULONG Nhalf=static_cast<ULONG>(N1)*static_cast<ULONG>(N2)*static_cast<ULONG>(N3/2+1);  	
#ifdef DOUBLE_PRC
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }
  // Convert input field to Fourier space
  //read the kernel in 3D
  string fname=file_kernel+"kernel"+to_string(N1)+"V"+to_string(L1)+"r"+to_string(smol);
  // read kernel in Fourier space
  vector<real_prec> kern(Nhalf,0);
  get_scalar(fname,kern,N1,N2,N3/2+1);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<Nhalf;i++)
    {
      real_prec cor=kern[i];
      AUX[i][REAL]*=cor;
      AUX[i][IMAG]*=cor;
    }

  do_fftw_r2c(N1,in,AUX);

  kern.clear();
  kern.shrink_to_fit();

  do_fftw_c2r(N1,AUX,out);

#ifdef DOUBLE_PRC
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif
 }
#endif
 ////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 void convolvek(ULONG N1, vector<real_prec>&in, vector<real_prec> &kernel, vector<real_prec> &out)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifdef _FULL_VERBOSE_
cout<<YELLOW<<"\tConvolution"<<RESET<<endl;
#endif
  ULONG N=in.size();
  ULONG Nhalf=kernel.size();
#ifdef DOUBLE
  complex_prec *AUX= (complex_prec *)fftw_malloc(2*Nhalf*sizeof(real_prec));
#else
  complex_prec *AUX= (complex_prec *)fftwf_malloc(2*Nhalf*sizeof(real_prec));
#endif
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for(ULONG i=0;i<Nhalf;++i)
    {
      AUX[i][REAL]=0;
      AUX[i][IMAG]=0;
    }
  // Convert input field to Fouroer space
  do_fftw_r2c(N1,in,AUX);
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
  for (ULONG i=0;i<Nhalf;i++)
    {
      real_prec cor=kernel[i];
      AUX[i][REAL]*=cor;
      AUX[i][IMAG]*=cor;
      }
  do_fftw_c2r(N1,AUX,out);
#ifdef DOUBLE
  fftw_free(AUX);
#else
  fftwf_free(AUX);
#endif
#ifdef _FULL_VERBOSE_
       std::cout<<BOLDGREEN<<"\t\t\t\t\t\t["<<BOLDBLUE<<"DONE"<<BOLDGREEN<<"]"<<RESET<<endl;
#endif
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 real_prec linearvel3d(int index, real_prec kx, real_prec ky, real_prec kz, real_prec phi)
{
  real_prec out, kl;
  switch(index)
    {
    case 1:
      {
	kl=kx;
	break;
      }
    case 2:
      {
	kl=ky;
	break;
      }
    case 3:
      {
	kl=kz;
	break;
      }
    }
  real_prec kmod2=kx*kx+ky*ky+kz*kz;

  real_prec fackern=0.0;
  if (kmod2>eps)
    fackern = kl/kmod2;
  out=fackern*phi;
  return out;
}
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 void sort_2vectors(vector<vector<ULONG> >& v1,vector< vector<ULONG> >&v2){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
   ULONG i, j;
   ULONG NN=v1.size();
   ULONG n=v1[0].size();
   for (i=0;i<NN;++i)
     {
       vector<ULONG>iwksp(n,0);
       vector<float>wksp(n,0);
       vector<ULONG>av1(n,0);
       for (j=0;j<n;++j) av1[j]=v1[i][j];
       indexx_ulong(av1,iwksp); //av1 must be int
       
       for (j=0;j<n;++j) wksp[j]=av1[j];
       for (j=0;j<n;++j) av1[j]=wksp[iwksp[j]];
       for (j=0;j<n;++j) v1[i][j]=av1[j];
       
       for (j=0;j<n;++j) wksp[j]=v2[i][j];
       for (j=0;j<n;++j) v2[i][j]=wksp[iwksp[j]];
     }
 }
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 void sort_1dvectors(vector<ULONG> & v1, vector<ULONG> &v2){
   // v1 is to be sorted. The elements of the vector v2 are shuffled accordingly. 
   ULONG i, j;
   ULONG n=v1.size();
   vector<ULONG>iwksp(n,0);
   vector<ULONG>wksp(n,0);
   vector<ULONG>av1(n,0);
   // Prepare working space
   for (j=0;j<n;++j) av1[j]=v1[j];
  // sort v1. iwksp contains the order of
   //v1= [ 10 , 8, 25, 54, 1] -> v1=  returns the smae and
  //        0   1   2   3  4    iwksp= [4,1,0,2,3] is the rank
   indexx_ulong(v1,iwksp); 
    // order av1
   for (j=0;j<n;++j) wksp[j]=av1[j];
   for (j=0;j<n;++j) v1[j]=wksp[iwksp[j]];
    // order v2
   for (j=0;j<n;++j) wksp[j]=v2[j];
    // allocate values of v2 sorted
   for (j=0;j<n;++j) v2[j]=wksp[iwksp[j]];
 }
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
  void min_max_vector(vector<real_prec>vec, size_t &imin, size_t &imax)
  {
    gsl_vector * v = gsl_vector_alloc (vec.size());
    for (ULONG i= 0; i < vec.size(); i++)
      gsl_vector_set (v, i, vec[i]);
    gsl_vector_minmax_index(v, &imin, &imax);
  }

 //##################################################################################
 void sort_1dvectors_v2(vector<ULONG> &v1, vector<ULONG> &v2, ULONG &v1cero, ULONG &v2cero){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
   ULONG n=v1.size();
   vector<ULONG>iwksp(n,0);
   indexx_ulong(v1,iwksp); //av1 must be int
   v1cero=v1[iwksp[0]];
   v2cero=v2[iwksp[0]];
  }
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
 void sort_1dvectors_iv2(vector<int> &v1, vector<int> &v2, int &v1cero, int &v2cero){
   // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
   ULONG n=v1.size();
   vector<int>iwksp(n,0);
   indexx(v1,iwksp); //av1 must be int
   v1cero=v1[iwksp[0]];
   v2cero=v2[iwksp[0]];
  }
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
 void sort_1dvectors_v3(vector<ULONG> & v1, vector<ULONG> &v2,  ULONG &v2cero){
  // v1 is to be sorted. The elements of the other vectors
   // are shuffled accordingly.
  // No ordered vector is returned, only the zero element of v2 sorted according to v1
   vector<ULONG>iwksp(v1.size(),0);
   indexx_ulong(v1,iwksp); 
   v2cero=v2[iwksp[0]];
 }
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
 void  reduce_resolution(real_prec L, ULONG Nft_h, ULONG Nft_l, vector<real_prec>&field_hr, vector<real_prec>&field_lr){
#ifdef _FULL_VERBOSE_
  cout<<YELLOW<<"Reducing resolution in dm"<<RESET<<endl;
#endif

     ULONG NLOSS=0;
     ULONG count=0;

     real_prec dl=L/static_cast<real_prec>(Nft_l);
     for (ULONG ip=0; ip<Nft_h; ip++)
     {
       real_prec xp=ip*dl;
       for (ULONG jp=0; jp<Nft_h; jp++)
       {
           real_prec yp=jp*dl;
          for (ULONG kp=0; kp<Nft_h; kp++)
          {
              real_prec zp=kp*dl;
              //check if particle is in selected Domain, else discard it
           ULONG i = static_cast<ULONG>(floor(xp/dl)); // indices of the cell of the particle
           ULONG j = static_cast<ULONG>(floor(yp/dl));
           ULONG k = static_cast<ULONG>(floor(zp/dl));

           i = static_cast<ULONG>(fmod(real_prec(i),real_prec(Nft_l)));
           j = static_cast<ULONG>(fmod(real_prec(j),real_prec(Nft_l)));
           k = static_cast<ULONG>(fmod(real_prec(k),real_prec(Nft_l)));
           ULONG lp_r=index_3d(i,j,k,Nft_l,Nft_l/2+1);;
           ULONG lp_h=index_3d(ip,jp,kp,Nft_h,Nft_h/2+1);;
           field_lr[lp_r] +=field_hr[lp_h];
        }
       }
     }
 }
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 void average(real_prec Lside, ULONG Nft_HR, ULONG Nft_LR, vector<real_prec>&HR_field, vector<real_prec>&LR_field)
 {
   cout<<"Averaging from high to low res"<<endl;
   ULONG NTT_HR=static_cast<ULONG>(Nft_HR*Nft_HR*(Nft_HR/2+1));
   ULONG NTT_LR=static_cast<ULONG>(Nft_LR*Nft_LR*(Nft_LR/2+1));
   ULONG NGRID_HR=static_cast<ULONG>(Nft_HR*Nft_HR*Nft_HR);
   ULONG NGRID_LR=static_cast<ULONG>(Nft_LR*Nft_LR*Nft_LR);
   if(NTT_HR<NTT_LR)
     {
       cout<<RED<<"Error. High resolution mesh must have larger number of Nft cells than low resolution."<<endl;
       exit(0);
     }
   real_prec delta=Lside/static_cast<real_prec>(Nft_LR);
   real_prec delta_hr=Lside/static_cast<real_prec>(Nft_HR);
   getDensity_NGP(Nft_LR,Nft_HR,Lside,delta,delta_hr,HR_field,LR_field);

   for(ULONG i=0;i<NGRID_LR; ++i)
     LR_field[i]/=static_cast<double>(pow(Nft_HR/Nft_LR,3));   //average
 }
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
 real_prec linInterpol(ULONG N, real_prec L, real_prec delta, real_prec xpos, real_prec ypos, real_prec zpos, const vector<real_prec>&tab)
 {
#ifndef CELLBOUND
   //* If ABACUS is enabled, num_0_5 -BIN_SHIFT =0 * //
   // * and there is no need to redefine particle posistions, as these already come referred to the corners*/
   xpos-=num_0_5*delta;
   ypos-=num_0_5*delta;
   zpos-=num_0_5*delta;
   {
     real_prec xnew=xpos;
     real_prec ynew=ypos;
     real_prec znew=zpos;

     if (xnew<0.)
       xnew+=L;
     if (xnew>=L)
       xnew-=L;

     if (ynew<0.)
       ynew+=L;
     if (ynew>=L)
       ynew-=L;

     if (znew<0.)
       znew+=L;
     if (znew>=L)
       znew-=L;
     xpos=xnew;
     ypos=ynew;
     zpos=znew;
   }
#endif
  ULONG ix = static_cast<ULONG>(floor(xpos/delta+BIN_SHIFT));
  ULONG iy = static_cast<ULONG>(floor(ypos/delta+BIN_SHIFT));
  ULONG iz = static_cast<ULONG>(floor(zpos/delta+BIN_SHIFT));
  real_prec tx = (xpos-static_cast<real_prec>(ix)*delta)/delta;
  real_prec ty = (ypos-static_cast<real_prec>(iy)*delta)/delta;
  real_prec tz = (zpos-static_cast<real_prec>(iz)*delta)/delta;
  ULONG shiftx=1;
  ULONG shifty=1;
  ULONG shiftz=1;
  long ixs=ix+shiftx;
  long iys=iy+shifty;
  long izs=iz+shiftz;
  if (ixs>=N)
    ixs-=N;
  if (iys>=N)
    iys-=N;
  if (izs>=N)
    izs-=N;
  real_prec y1 = tab[iz+N*(iy+N*ix)];
  real_prec y2 = tab[izs+N*(iy+N*ix)];
  real_prec y3 = tab[izs+N*(iys+N*ix)];
  real_prec y4 = tab[izs+N*(iys+N*ixs)];
  real_prec y5 = tab[iz+N*(iys+N*ix)];
  real_prec y6 = tab[iz+N*((iys)+N*ixs)];
  real_prec y7 = tab[iz+N*(iy+N*ixs)];
  real_prec y8 = tab[izs+N*(iy+N*ixs)];
  real_prec  out = static_cast<real_prec>((1.-tx)*(1.-ty)*(1.-tz)*y1 + tx*(1.-ty)*(1.-tz)*y2 + tx*ty*(1.-tz)*y3 + tx*ty*tz*y4 + (1.-tx)*ty*(1.-tz)*y5 + (1.-tx)*ty*tz*y6 + (1.-tx)*(1.-ty)*tz*y7 + tx*(1.-ty)*tz*y8);
  return out;
}
 ////////////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////////////
void getDensity_NGP(ULONG Nft, ULONG Nft_HR, real_prec L1, real_prec d1, real_prec delta_hr, const vector<real_prec>&we, vector<real_prec>&delta)
{
#ifdef _USE_OMP_
  int NTHREADS = _NTHREADS_;
  omp_set_num_threads(NTHREADS);
#endif
#ifndef _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
#endif
  ULONG count_previous=0;
#pragma omp parallel for reduction(+:count_previous)
  for(ULONG ip=0;ip<delta.size();++ip)
    count_previous+=delta[ip];
  ULONG counter_a=0;
#pragma omp parallel for reduction(+:counter_a)
  for (ULONG ig=0; ig<we.size(); ++ig)
    {
      ULONG xc,yc,zc;
      index2coords(Nft_HR,ig,xc,yc,zc);
      real_prec xpos=(static_cast<real_prec>(zc)+0.5)*delta_hr;
      real_prec ypos=(static_cast<real_prec>(yc)+0.5)*delta_hr;
      real_prec zpos=(static_cast<real_prec>(xc)+0.5)*delta_hr;
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L1;
      if(zpos< 0)zpos+=L1;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L1)ypos-=L1;
      if(zpos>=L1)zpos-=L1;
      ULONG i = static_cast<ULONG>(floor((xpos)/d1)); // indices of the cell of the particle
      ULONG j = static_cast<ULONG>(floor((ypos)/d1));
      ULONG k = static_cast<ULONG>(floor((zpos)/d1));
#pragma omp atomic update
        delta[k+Nft*j+Nft*Nft*i]+=we[ig];
        counter_a++;
    }
  ULONG count=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count)
#endif
  for(ULONG ip=0;ip<delta.size();++ip)
    count+=delta[ip];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid  = "<<count<<endl;
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// THis is inherited from Galaxy::bspline
void smooth_distribution(vector<real_prec>&x, vector<real_prec>&F,vector<real_prec>&Fs)
{
    int nfactor=30;
    const gsl_rng_type * rng_t;
    gsl_rng * gBaseRand;
    gsl_rng_env_setup();
    real_prec xmax=x[x.size()-1];
    real_prec xmin=x[0];
    // make a copy of the original zbins
    vector<gsl_real>za_original(x.size());
    for(int i=0;i<x.size();i++)za_original[i]=static_cast<gsl_real>(x[i]);
    vector<gsl_real>nz_original(F.size());
    for(int i=0;i<F.size();i++)nz_original[i]=static_cast<gsl_real>(F[i]);
    // Transform input type to gsl_real
    vector<gsl_real>xx(x.size());
    for(int i=0;i<x.size();i++)xx[i]=static_cast<gsl_real>(x[i]);
    vector<gsl_real>FF(F.size());
    for(int i=0;i<F.size();i++)FF[i]=static_cast<gsl_real>(F[i]);
    // containers whose size will be reduced
    vector<gsl_real> new_zn;
    vector<gsl_real> new_nb;
    int nn0=static_cast<int>(floor(x.size()/nfactor));
    new_zn.resize(nn0);
    new_nb.resize(nn0);
    gsl_bspline(xx,FF,new_zn, new_nb);// smooth the nv(zn) array to a new new_nb(new_zn)
    for(int i=0;i<x.size();i++)
        new_nb[i]=abs(new_nb[i]); // ensure positiviness
    // Redefine mins and max from the smoothed array
    xmax=new_zn[new_zn.size()-1];
    xmin=new_zn[0];
    gsl_rng_default_seed=25225;
    rng_t = gsl_rng_ranlux;
    gBaseRand = gsl_rng_alloc (rng_t);
    // I need to change types
    vector<real_prec>new_nbf(new_nb.size());
    for(int i=0;i<new_nb.size();i++)new_nbf[i]=static_cast<real_prec>(new_nb[i]);
    vector<real_prec>new_znf(new_zn.size());
    for(int i=0;i<new_nb.size();i++)new_znf[i]=static_cast<real_prec>(new_zn[i]);
    vector<gsl_real>zran(2*x.size(),0);
    for(int ir=0;ir<zran.size();++ir)
       zran[ir]=xmin+(xmax-xmin)*gsl_rng_uniform(gBaseRand);
    vector<gsl_real>zran_s(zran.size(),0);
    sort_1d_vectors<gsl_real>(zran,zran_s);

    // FOrce the random generated vector has the sae limits as the original one
    zran_s[0]=za_original[0];
    zran_s[zran_s.size()-1]=za_original[za_original.size()-1];

    vector<gsl_real>dNran(zran.size(),0);
    for(int ir=0;ir<zran_s.size();++ir)// INterpola el prodcto de lbspline en los nuevos z randoms sorted
      dNran[ir]=gsl_inter_pointers(&new_znf[0], &new_nbf[0], new_znf.size(), zran_s[ir]);

    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_linear,  zran_s.size());
    gsl_spline_init (spline, &zran_s[0], &dNran[0], zran_s.size());
    for(int i=0;i<za_original.size();i++)
        Fs[i]=gsl_spline_eval (spline, za_original[i], acc);
}
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
