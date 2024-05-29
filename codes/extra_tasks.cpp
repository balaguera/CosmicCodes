
# include "CosmiCalcLIB.h"
/////////////////////////////////////////////////////////////////////////////////////////
void get_tidal_anisotropy(string parfile){
    Params params(parfile);
    Cwclass cwc(params);
    string infile=params._Input_Directory_X()+params._Name_Catalog_X();
    string outfile=params._Output_directory()+"tidal";
    string outfile_pot=params._Output_directory()+"potential";
    FileOutput File;
    vector<real_prec>delta(params._NGRID(),0);
    File.read_array(infile, delta);
    real_prec mean = get_mean(delta);
    get_overdens(delta,mean, delta);
    vector<real_prec>potential(params._NGRID(),0);
    ULONG NTT = static_cast<ULONG>(params._Nft()*params._Nft()*(params._Nft()/2+1));
    complex_prec *data_out= (complex_prec *)fftwf_malloc(2*NTT*sizeof(real_prec));
    do_fftw_r2c(params._Nft(),delta,data_out);
    cout<<data_out[0][REAL]<<endl;
    cout<<data_out[0][IMAG]<<endl;
    cout<<data_out[1][REAL]<<endl;
    cout<<data_out[1][IMAG]<<endl;
    cout<<data_out[2][REAL]<<endl;
    cout<<data_out[2][IMAG]<<endl;

/*
    PoissonSolver(params._Lbox(),params._Nft(),delta,potential);
    File.write_array(outfile_pot, potential);
    cwc.get_CWC(delta);
#pragma omp parallel for
    for(ULONG i=0;i<delta.size();++i)
        delta[i]=tidal_anisotropy(cwc.lambda1[i],cwc.lambda2[i],cwc.lambda3[i]);
    File.write_array(outfile, delta);
*/

}
/////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]){
    time_t start_all;
    time(&start_all);
    string pf=argv[1];
    get_tidal_anisotropy(pf);
}
















