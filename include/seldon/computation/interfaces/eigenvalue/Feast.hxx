#ifndef SELDON_FILE_FEAST_HXX
#define SELDON_FILE_FEAST_HXX

extern "C"
{
  void feastinit_(int*);
  
  void dfeast_srci_(int*, const int*, void*, double*, void*, double*, double*,
                   int*, double*, int*, double*, double*, int*,
                   double*, double*, int*, double*, int*);

  void zfeast_hrci_(int*, int*, void*, void*, void*, void*, void*,
                    int*, double*, int*, double*, double*, int*,
                    double*, void*, int*, double*, int*);
  
}

namespace Seldon
{

  void CallFeast(int& ijob, int n, complex<double>& ze,
		 Matrix<double, General, ColMajor>& work,
		 Matrix<complex<double>, General, ColMajor>& workc,
		 Vector<double>& aq, Vector<double>& sq,
		 IVect& fpm, double& epsout, int& loop,
		 double emin, double emax, int m0, Vector<double>& lambda,
		 Matrix<double, General, ColMajor>& eigen_vectors,
		 int& m, Vector<double>& res, int& info);
  
  void CallFeast(int& ijob, int n, complex<double>& ze,
		 Matrix<complex<double>, General, ColMajor>& work,
		 Matrix<complex<double>, General, ColMajor>& workc,
		 Vector<complex<double> >& aq, Vector<complex<double> >& sq,
		 IVect& fpm, double& epsout, int& loop,
		 double emin, double emax, int m0, Vector<double>& lambda,
		 Matrix<complex<double>, General, ColMajor>& eigen_vectors,
		 int& m, Vector<double>& res, int& info);
 
  // main function to find eigenvalues and eigenvectors with Feast (MKL implementation)
#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  void FindEigenvaluesFeast(EigenProblem_Base<T>& var,
                            Vector<T>& eigen_values,
                            Vector<T>& lambda_imag,
                            Matrix<T, General, ColMajor>& eigen_vectors);
#else
  template<class EigenProblem, class T, class Allocator1,
	   class Allocator2, class Al3>
  void FindEigenvaluesFeast(EigenProblem& var,
                            Vector<T, VectFull, Allocator1>& eigen_values,
                            Vector<T, VectFull, Allocator2>& lambda_imag,
                            Matrix<T, General, ColMajor, Al3>& eigen_vectors);
#endif
  
}

#endif
