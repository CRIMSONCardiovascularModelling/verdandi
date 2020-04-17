#ifndef SELDON_FILE_ANASAZI_HXX

#include "AnasaziConfigDefs.hpp"
#include "AnasaziBasicEigenproblem.hpp"
#include "AnasaziLOBPCG.hpp"
#include "AnasaziLOBPCGSolMgr.hpp"
#include "AnasaziBlockDavidson.hpp"
#include "AnasaziBlockDavidsonSolMgr.hpp"
#include "AnasaziBlockKrylovSchur.hpp"
#include "AnasaziBlockKrylovSchurSolMgr.hpp"
#include "AnasaziBasicSort.hpp"
#include "AnasaziBasicOutputManager.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
// #include "AnasaziEpetraAdapter.hpp"
#include "AnasaziMVOPTester.hpp"

#ifdef SELDON_WITH_MPI
//#include "Epetra_MpiComm.h"
#else
//#include "Epetra_SerialComm.h"
#endif

// #include "Epetra_Map.h"

// templated multivector
#include "AnasaziMultiVec.hpp"
#include "MyMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include "Teuchos_BLAS.hpp"

namespace Anasazi
{
  using namespace std;
  
  template<class EigenPb, class T>
  class OperatorAnasaziEigen : public Operator<T>
  {
  public :
    //! reference to the object containing eigenproblem data    
    EigenPb& var_eig;
    //! operator stored (M, A, Op)
    int type_operator;
    enum {OPERATOR_A, OPERATOR_M, OPERATOR_OP};
    
    OperatorAnasaziEigen(EigenPb& var, int type);

    void Apply(const MultiVec<T>& X0, MultiVec<T>& Y0) const;
    
  };

}

namespace Seldon
{
  
  template<class T>
  void SetComplexEigenvalue(const T& x, const T& y, T& z, T& zimag);
  
  template<class T>
  void SetComplexEigenvalue(const T& x, const T& y, complex<T>& z, complex<T>& zimag);

#ifdef SELDON_WITH_VIRTUAL
  template<class T>
  void FindEigenvaluesAnasazi(EigenProblem_Base<T>& var,
			      Vector<T>& eigen_values,
			      Vector<T>& lambda_imag,
			      Matrix<T, General, ColMajor>& eigen_vectors);
#else
  template<class EigenProblem, class T, class Allocator1,
           class Allocator2, class Allocator3>
  void FindEigenvaluesAnasazi(EigenProblem& var,
			      Vector<T, VectFull, Allocator1>& eigen_values,
			      Vector<T, VectFull, Allocator2>& lambda_imag,
			      Matrix<T, General, ColMajor, Allocator3>& eigen_vectors);
#endif
  
}

#define SELDON_FILE_ANASAZI_HXX
#endif
