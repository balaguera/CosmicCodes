/**
 * @class<FileOutput>
 * @brief Header file for the class FileOutput::
 * @file FileOutput.h
 * @title Methods for input/output in BAM
 * @author Andres Balaguera-Antolínez
 * @version   1.0
 * @date      2020
*/

#ifndef __FILE_OUTPUT__
#define __FILE_OUTPUT__
//# include <CCfits>
# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <iomanip>
# include <sstream>
# include <cassert>
# include <vector>
# include "fftw_array.h"
# include "ScreenOutput.h"
# include <bstream.h>
// ***************************************************************************************************************
// ***************************************************************************************************************
using namespace std;
// ***************************************************************************************************************
// ***************************************************************************************************************

class FileOutput{
 private:
    /**
     * @brief Object of type bifstream
     */
  bifstream inStream;

 public:
  
  //////////////////////////////////////////////////////////
  /** 
   * @brief Constructor
   */
  FileOutput(){
#ifdef SINGLE_PREC
    this->out.precision(10);
#else
    this->out.precision(10);
#endif
    this->out.setf(ios::showpoint);
    this->out.width(8);
   }
  //////////////////////////////////////////////////////////
  /**
   * @brief  Destructor
   */
  ~FileOutput(){}
  //////////////////////////////////////////////////////////
  /**
   * @brief   String
   */
  string input_type;
  //////////////////////////////////////////////////////////
  /**
   * @brief  Object of type ScreenOutput
   */
  ScreenOutput So;
  //////////////////////////////////////////////////////////
  /**
   * @brief Reads input file and allocate it in a vector
   * @param ifile input file (ascii)
   * @authors Daniele Tavagnacco
   * @authors Andres Balaguera
   */
  int read_file(string ifile,vector<real_prec > &prop);
 //////////////////////////////////////////////////////////
  /**
   * @brief Object of type std::ifstream
   */
  ifstream in;
  //////////////////////////////////////////////////////////
  /**
   * @brief Object of type std::ofstream
   */
  ofstream out;
  //////////////////////////////////////////////////////////
  /**
   * @brief Reads input file and allocate it in a matrix
   * @authors Daniele Tavagnacco
   * @authors Andres Balaguera
   * @param ifile input file (ascii)
   * @result prop 2d container with the number of lines equal to the number of objects in the
   * catalogue and the number of rows equal to the number of columns in the catalogue
   */
  void read_file(string ifile, vector< vector< real_prec > > &prop);
  //////////////////////////////////////////////////////////
  /**
   * @brief Verify if file exists
   * @params file: string name of file
   */
  bool exists(string file);
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   * @authors Daniele Tavagnacco
   * @authors Andres Balaguera
   */
  int read_file_one(string, vector< real_prec > &); // For files with one column
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
//  int read_file_one(string, vector< float > &); // For files with one column
  //////////////////////////////////////////////////////////
  /**
   * @brief  
   */
  void read_file_one_N(long N, string, vector< real_prec > &); // For files with one column
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   * @authors Daniele Tavagnacco
   * @authors Andres Balaguera
   */
  ULONG read_file(string fname, vector<real_prec> &prop_ob, int NTHREADS);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  ULONG read_binary_file(string fname, vector<real_prec> &prop_ob);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_to_file(string,  vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void write_to_file_i(string, vector<ULONG>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief
   */
  void write_to_file_i(string,  vector<int>&);
  //////////////////////////////////////////////////////////
  /**
   * @brief 
   */
  void write_to_file(string, vector<real_prec>& , vector<real_prec>&);
  //////////////////////////////////////////////////////////
#ifdef SINGLE_PREC
  /**
   * @brief 
   */
  void write_to_file(string fname, vector<gsl_real>& kve, vector<real_prec> &pk);
#endif
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&,vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec>&, vector<real_prec> &,vector<real_prec> &,vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<int>& nmod);
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<int>& nmod);
  //////////////////////////////////////////////////////////
  /** 
   *@brief  Write output
   */
  void write_to_file(string fname, vector<real_prec>& kve, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<real_prec>& nmod);
  //////////////////////////////////////////////////////////
  /**
   *@brief  Write output
   */
  void write_to_file(string fname, vector<real_prec>&kve, vector<real_prec> &pk0, vector<real_prec> &pk, vector<real_prec> &pkk,vector<real_prec> &pkkk, vector<real_prec> &si, vector<int>& nmod);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector< real_prec >&, int);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec> &, vector<real_prec>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec> &, vector<int>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<int>&, vector< vector<real_prec> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<real_prec>&, vector<real_prec>&, vector<real_prec> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<real_prec>&, vector<real_prec> &, vector<int> &);
#ifdef SINGLE_PREC
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<real_prec> &, vector<int> &);
#endif
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file2(string, vector<real_prec>&, vector<real_prec> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file2(string, vector<real_prec>&, vector<real_prec> &, vector<int> &, bool);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string fname, vector<real_prec> &kve, vector< real_prec> &bis, vector< real_prec> &sn_bis, vector<int> &mod);

  //////////////////////////////////////////////////////////  
  /** 
   *@brief   Write output
   */
  void read_array(string, float *, ULONG);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void read_array(string, vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  void read_array(string, vector<ULONG>&);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  void read_array_fast(string, vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_array(string, float *, ULONG);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  template<class Type> void read_array(string, vector<real_prec>&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_array(string fname, vector<real_prec>&out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  void write_array(string fname, vector<ULONG>&out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write complex output
   */
  void write_array(string fname, int Nft, complex_prec *out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  void write_array(string fname, vector<int>&out);
  //////////////////////////////////////////////////////////
  /**
   *@brief   Write output
   */
  ULONG count_lines(string);
  //////////////////////////////////////////////////////////
  /**
   * template<read_array_t>
   *@brief Read binary file
   */
template<class Type> void read_array_t(string fname, vector<real_prec>&OUT){
      ULONG N=OUT.size();
      fftw_array<Type> auxil(N);
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Reading binary file", fname);
#endif
      this->inStream.open(fname.data(),file_is_natural);
      if(!this->inStream.is_open())
        {
          cout<<RED<<"File not found!."<<RESET<<endl;
          cout<<RED<<"Check input parameter file."<<RESET<<endl;
          cout<<RED<<"Code exits here."<<RESET<<endl;
          cout<<endl;
          exit(0);
        }
      assert(this->inStream.is_open());
      this->inStream.get(auxil.data,N);
      this->inStream.close();
#pragma omp parallel for
      for(ULONG i=0;i<N;i++)
        OUT[i]=static_cast<real_prec>(auxil[i]);
      So.DONE();
     return ;
  }

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
  /**
   * template<read_array_t>
   *@brief Read binary file
   */
template<class Type> void read_array_t(string fname, vector<Type>&OUT){
#ifdef _USE_OMP_
     int NTHREADS=_NTHREADS_;
 omp_set_num_threads(NTHREADS);
#endif

      ULONG N=OUT.size();
      fftw_array<Type> aux(N);
#ifdef _FULL_VERBOSE_
      this->So.message_screen("Reading binary file", fname);
#endif
      this->inStream.open(fname.data(),file_is_natural);
      if(!this->inStream.is_open())
        {
          cout<<RED<<"File not found!."<<RESET<<endl;
          cout<<RED<<"Check input parameter file."<<RESET<<endl;
          cout<<RED<<"Code exits here."<<RESET<<endl;
          cout<<endl;
          exit(0);
        }
      assert(this->inStream.is_open());
      this->inStream.get(aux.data,N);
      this->inStream.close();
#ifdef _USE_OMP_
#pragma omp parallel for
#endif
      for(ULONG i=0;i<N;i++)
        OUT[i]=static_cast<real_prec>(aux[i]);
      So.DONE();
     return ;
  }
  //////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////
  /**
   * template<read_array_t>
   *@brief Read binary file
   */
 template<class Type> void read_array_tt(string fname, vector<real_prec>&OUT){
  #ifdef _USE_OMP_
       int NTHREADS=_NTHREADS_;
       omp_set_num_threads(NTHREADS);
  #endif

  #ifdef _FULL_VERBOSE_
        this->So.message_screen("Reading binary file", fname);
  #endif
        ifstream fileb(fname.c_str(), ios::binary| ios::in);
        while(!fileb.eof())
        {
           Type dat;
#ifdef _TNG_
           fileb.read((char*)&dat, 2*sizeof(Type));
#elif defined _TNG_GAL_
           fileb.read((char*)&dat, sizeof(Type));
#endif
           OUT.push_back(static_cast<real_prec>(dat));
        }
        this->So.message_screen("Read", OUT.size(), " lines");
        cout<<OUT[0]<<"  "<<OUT[OUT.size()-1]<<endl;
        fileb.close();
        So.DONE();
       return ;
    }

//////////////////////////////////////////////////////////
  /**
   * template<read_array_t>
   *@brief Read binary file
   */
template<typename T> void write_to_file(const std::string& fname,const std::vector<std::vector<T>>& columns)
{
    std::ofstream out(fname, std::ios::out);
    if (!out)
        throw std::runtime_error("Cannot open file " + fname);
    if (columns.empty())
        return;
    size_t nrows = columns[0].size();
    size_t ncols = columns.size();
    // Safety check
    for (const auto& col : columns)
        if (col.size() != nrows)
            throw std::runtime_error("Column size mismatch");

    So.message_output_file(fname, nrows, ncols);

    for (size_t i = 0; i < nrows; ++i) {
        for (size_t j = 0; j < ncols; ++j)
            out << columns[j][i] << "\t";
        out << "\n";
    }

    out.close();
    So.DONE();
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
 
#if !defined(USE_GALAXY_TOOLS)
  /**
   * @brief Write to output file
   * @details Outputin ASCII format
   * @param file: <string> name of output file
   * @param data: <vector> vector to print out
   */
 void write_to_file(string file,  vector<gsl_real>&data);
//////////////////////////////////////////////////////////
  /**
   * @brief Write output
   */
  void write_to_file(string, vector<gsl_real>& , vector<gsl_real>&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&,vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&,vector<gsl_real> &,vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real> &,vector<gsl_real> &,vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /**
   *@brief Write output
   */
  void write_to_file(string fname, vector<gsl_real>& kve, vector<gsl_real> &pk, vector<gsl_real> &pkk,vector<gsl_real> &pkkk, vector<gsl_real> &si, vector<int>& nmod);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector< gsl_real >&, int);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real> &, vector<gsl_real>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real> &, vector<int>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<int>&, vector< vector<gsl_real> >&);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<int> &, vector<gsl_real>&, vector<gsl_real>&, vector<gsl_real> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string, vector<gsl_real>&, vector<gsl_real> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file2(string, vector<gsl_real>&, vector<gsl_real> &, vector<int> &);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  void write_to_file(string fname, vector<gsl_real> &kve, vector< gsl_real> &bis, vector< gsl_real> &sn_bis, vector<int> &mod);
  //////////////////////////////////////////////////////////
  /** 
   *@brief   Write output
   */
  bifstream inStreamp;
#endif
};



#endif

