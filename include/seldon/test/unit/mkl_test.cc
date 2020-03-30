// NewAlloc as default allocator in order to avoid problems
// with vectors of complex types
#define SELDON_DEFAULT_ALLOCATOR NewAlloc
// seldon will call abort() when encountering an exception
#define SELDON_WITH_ABORT
// no call of srand by Seldon
#define SELDON_WITHOUT_REINIT_RANDOM

// C library for time function and for randomization
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <stdint.h>

#ifdef SELDON_WITH_MPI
#include "mpi.h"
#endif

#define SELDON_WITH_PRECONDITIONING

#include "Seldon.hxx"
#include "SeldonComplexMatrix.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

double threshold;

namespace std
{
  inline bool isnan(const complex<double>& x)
  {
    if (isnan(real(x)))
      return true;
    
    if (isnan(imag(x)))
      return true;
    
    return false;
  }
}

template<class T>
void GetRandNumber(T& x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRandNumber(complex<T>& x)
{
  int type = rand()%3;
  if (type == 0)
    x = complex<T>(0, rand())/T(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/T(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
}

template<class T>
void GenerateRandomVector(Vector<T, VectSparse>& x, int n, int nnz)
{
  x.Clear();
  for (int i = 0; i < nnz; i++)
    {
      int j = rand()%n;
      GetRandNumber(x.Get(j));
    }
}

template<class T>
bool EqualVector(const Vector<T>& x, const Vector<T>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > threshold) || isnan(abs(x(i)-y(i))))
      {
        DISP(i); DISP(x(i)); DISP(y(i));
        return false;
      }
  
  return true;
}

template<class T>
bool EqualVector(const Vector<T, VectSparse>& x, const Vector<T, VectSparse>& y)
{
  if (x.GetM() != y.GetM())
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ( (x.Index(i) != y.Index(i)) || (abs(x.Value(i) - y.Value(i)) > threshold) 
         || isnan(abs(x.Value(i) - y.Value(i))) )
      {
        DISP(i); DISP(x.Value(i)); DISP(y.Value(i));
        return false;
      }
  
  return true;
}

template<class T, class Prop, class Storage,
         class T2, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop, Storage>& A,
                 const Matrix<T2, Prop2, Storage2>& B)
{
  if (A.GetM() != B.GetM())
    return false;

  if (A.GetN() != B.GetN())
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j)-B(i, j))))
        {
          DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
          return false;
        }
  
  return true;
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n, int nnz)
{
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  A.Clear();
  A.Reallocate(m, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%m;
      int j = rand()%n;
      GetRandNumber(x);
      A.Set(i, j, x);
    }
}


template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomMatrix(Matrix<T, Prop, Storage, Allocator>& A,
                          int m, int n)
{
  A.Reallocate(m, n);
  T x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        GetRandNumber(x);
        A(i, j) = x;
      }
}

template<class T>
void CheckSparseVector(Vector<T>& y)
{
  int n = 40, nnz = 8;
  Vector<T, VectSparse> x;
  
  GenerateRandomVector(x, n, nnz);
  GenerateRandomVector(y, n);
  
  // testing Add
  Vector<T> y0(y);
  Vector<T, VectSparse> x0(x);

  T alpha;
  GetRandNumber(alpha);

  Add(alpha, x, y);
  
  Vector<T> xd(n); xd.Fill(0);
  for (int i = 0; i < x.GetM(); i++)
    xd(x.Index(i)) = x.Value(i);
  
  Vector<T> ybis(y0);
  Add(alpha, xd, ybis);

  if (!EqualVector(y, ybis))
    {
      cout << "Add incorrect" << endl;
      abort();
    }
  
  
  // testing DotProd
  ybis = y0; y = y0;
  T scal = DotProd(x, y);
  T scal_ref; SetComplexZero(scal_ref);
  for (int i = 0; i < n; i++)
    scal_ref += xd(i)*y(i);
  
  if ((abs(scal - scal_ref) > threshold) || isnan(scal))
    {
      cout << "DotProd incorrect" << endl;
      abort();
    }

  // testing DotProdConj
  scal = DotProdConj(x, y);
  SetComplexZero(scal_ref);
  for (int i = 0; i < n; i++)
    scal_ref += conj(xd(i))*y(i);
  
  if ((abs(scal - scal_ref) > threshold) || isnan(scal))
    {
      cout << "DotProdConj incorrect" << endl;
      abort();
    }

  // testing GatherSparseEntry
  Vector<T, VectSparse> xbis(x);
  xbis.Fill(0);
  for (int i = 0; i < x.GetM(); i++)
    xbis.Value(i) = y(x.Index(i));
  
  GatherSparseEntry(y, x);
  
  if (!EqualVector(x, xbis))
    {
      cout << "GatherSparseEntry incorrect" << endl;
      abort();
    }

  // testing GatherSparseEntryZero
  xbis.Fill(0); x.Fill(0);
  for (int i = 0; i < x.GetM(); i++)
    {
      xbis.Value(i) = y(x.Index(i));
      SetComplexZero(y(x.Index(i)));
    }
  
  ybis = y0;
  GatherSparseEntryZero(ybis, x);
  
  if ( (!EqualVector(y, ybis)) || (!EqualVector(x, xbis)))
    {
      cout << "GatherSparseEntryZero incorrect" << endl;
      abort();
    }  
  
  // testing ScatterSparseEntry
  x = x0; y.Fill(0);
  xbis = x0; ybis.Fill(0);
  for (int i = 0; i < x.GetM(); i++)
    ybis(x.Index(i)) = x.Value(i);
  
  ScatterSparseEntry(x, y);
  
  if (!EqualVector(y, ybis))
    {
      cout << "ScatterSparseEntry incorrect" << endl;
      abort();
    }

  // testing Add with dense vectors
  /*Vector<T> z, zbis;
  GenerateRandomVector(y, n);
  GenerateRandomVector(z, n);

  T beta;
  GetRandNumber(beta);
  
  DISP(alpha); DISP(beta);
  ybis = y; zbis = z;
  Add(alpha, y, beta, z);
  
  for (int i = 0; i < y.GetM(); i++)
    zbis(i) = beta*zbis(i) + alpha*ybis(i);
  
  DISP(z); DISP(zbis);
  if (!EqualVector(z, zbis))
    {
      cout << "Add incorrect " << endl;
      abort();
      }*/
}


template<class T>
void CheckSparseVectorReal(Vector<T>& y)
{
  int n = 40, nnz = 8;
  Vector<T, VectSparse> x;
  
  GenerateRandomVector(x, n, nnz);
  GenerateRandomVector(y, n);
  
  Vector<T> y0(y), ybis, xd(y);
  Vector<T, VectSparse> x0(x), xbis;

  // testing ApplyRot
  x = x0; y = y0;
  xd.Fill(0);
  for (int i = 0; i < x.GetM(); i++)
    xd(x.Index(i)) = x.Value(i);
  
  T teta, c, s;
  GetRandNumber(teta);
  c = cos(teta); s = sin(teta);
  ApplyRot(x, y, c, s);
  
  xbis = x0; ybis = y0;
  ApplyRot(xbis, ybis, c, s);
  
  if ( (!EqualVector(y, ybis)) || (!EqualVector(x, xbis)))
    {
      cout << "ApplyRot incorrect" << endl;
      abort();
    }  
}  

template<class T, class Prop, class Storage>
void CheckSparseMatrix(Matrix<T, Prop, Storage>& A)
{
  int m = 40, n = 45, nnz = 200;
  if (IsSymmetricMatrix(A))
    n = m;
  
  GenerateRandomMatrix(A, m, n, nnz);
  
  Vector<T> x, y, ybis, x0, y0;
  GenerateRandomVector(x, n);
  
  x0 = x;
  y.Reallocate(m);
  y.Fill(0);
  ybis.Reallocate(m);
  ybis.Fill(0);
  
  Matrix<T, General, RowMajor> B(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      B(i, j) = A(i, j);
  
  // testing Mlt
  Mlt(A, x, y);
  Mlt(B, x, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  T alpha; GetRandNumber(alpha);
  Mlt(alpha, A, x, y);
  Mlt(alpha, B, x, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Mlt(SeldonNoTrans, A, x, y);
  Mlt(SeldonNoTrans, B, x, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  Vector<T> xbis(x);
  GenerateRandomVector(y, m);
  ybis = y;
  Mlt(SeldonTrans, A, y, x);
  Mlt(SeldonTrans, B, ybis, xbis);
  
  if (!EqualVector(x, xbis))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }

  Mlt(SeldonConjTrans, A, y, x);
  Mlt(SeldonConjTrans, B, ybis, xbis);
  
  if (!EqualVector(x, xbis))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
  // testing MltAdd
  T beta; GetRandNumber(beta); 
  x = x0; y0 = y;
  MltAdd(alpha, A, x, beta, y);
  MltAdd(alpha, B, x, beta, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

  y = y0; ybis = y0;
  MltAdd(alpha, SeldonNoTrans, A, x, beta, y);
  MltAdd(alpha, SeldonNoTrans, B, x, beta, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

  GenerateRandomVector(y0, n);
  GenerateRandomVector(x, m);
  y = y0; ybis = y0;
  MltAdd(alpha, SeldonTrans, A, x, beta, y);
  MltAdd(alpha, SeldonTrans, B, x, beta, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

  y = y0; ybis = y0;
  MltAdd(alpha, SeldonConjTrans, A, x, beta, y);
  MltAdd(alpha, SeldonConjTrans, B, x, beta, ybis);
  
  if (!EqualVector(y, ybis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

} 

template<class T, class Prop, class Storage>
void CheckGeneralMatrix(Matrix<T, Prop, Storage>& A)
{
  int m = 40, n = m, nnz = 200;
  T val;
  
  GenerateRandomMatrix(A, m, n, nnz);
  for (int i = 0; i < m; i++)
    {
      GetRandNumber(val);
      A.AddInteraction(i, i, T(1)+val);
    }
  
  Vector<T> x, y, ybis;
  GenerateRandomVector(x, n);
  
  y.Reallocate(m);
  y.Fill(0);
  ybis.Reallocate(m);
  ybis.Fill(0);
  
  Matrix<T, General, RowMajor> B(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      B(i, j) = A(i, j);
  
  // testing Solve with lower part
  Solve(SeldonLower, SeldonNoTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonLower, SeldonNoTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonLower, SeldonTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonLower, SeldonTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonLower, SeldonConjTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val/conj(A(i, i));
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonLower, SeldonConjTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  // testing Solve with upper part
  Solve(SeldonUpper, SeldonNoTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonUpper, SeldonNoTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonUpper, SeldonTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonUpper, SeldonTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonUpper, SeldonConjTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val/conj(A(i, i));
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonUpper, SeldonConjTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

} 


template<class T, class Prop, class Storage>
void CheckSymmetricMatrix(Matrix<T, Prop, Storage>& A)
{
  int m = 40, n = m, nnz = 200;
  T val;
  
  GenerateRandomMatrix(A, m, n, nnz);
  for (int i = 0; i < m; i++)
    {
      GetRandNumber(val);
      A.AddInteraction(i, i, T(1)+val);
    }
  
  Vector<T> x, y, ybis;
  GenerateRandomVector(x, n);
  
  y.Reallocate(m);
  y.Fill(0);
  ybis.Reallocate(m);
  ybis.Fill(0);
  
  Matrix<T, General, RowMajor> B(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      B(i, j) = A(i, j);
  
  // testing Solve with upper part
  Solve(SeldonNoTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonNoTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = m-1; i >= 0; i--)
    {
      val = ybis(i);
      for (int j = i+1; j < m; j++)
        val -= A(i, j)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val/A(i, i);
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= A(j, i)*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonConjTrans, SeldonNonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val/conj(A(i, i));
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  Solve(SeldonConjTrans, SeldonUnit, A, x, y);
  
  ybis = x;
  for (int i = 0; i < m; i++)
    {
      val = ybis(i);
      for (int j = 0; j < i; j++)
        val -= conj(A(j, i))*ybis(j);
      
      ybis(i) = val;
    }
  
  if (!EqualVector(y, ybis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

} 


template<class T, class Prop, class Storage, class PropDense, class StorageDense>
void CheckSparseMatrixProd(Matrix<T, Prop, Storage>& A,
                           Matrix<T, PropDense, StorageDense>& B)
{
  int m = 40, n = 30, k = 34, nnz = 200;
  if (IsSymmetricMatrix(A))
    n = m;
  
  Matrix<T, PropDense, StorageDense> C, Cbis, Ad;
  
  // testing MltAdd
  GenerateRandomMatrix(A, m, n, nnz);
  GenerateRandomMatrix(B, n, k);
  GenerateRandomMatrix(C, m, k); 
  
  Ad.Reallocate(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Ad(i, j) = A(i, j);
  
  Cbis = C;
  
  T alpha, beta;
  GetRandNumber(alpha); GetRandNumber(beta);
  
  MltAdd(alpha, SeldonNoTrans, A, SeldonNoTrans, B, beta, C);
  MltAdd(alpha, SeldonNoTrans, Ad, SeldonNoTrans, B, beta, Cbis);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(A, n, m, nnz);
  Ad.Reallocate(n, m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      Ad(i, j) = A(i, j);

  GenerateRandomMatrix(C, m, k); 
  Cbis = C;
  
  MltAdd(alpha, SeldonTrans, A, SeldonNoTrans, B, beta, C);
  MltAdd(alpha, SeldonTrans, Ad, SeldonNoTrans, B, beta, Cbis);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }

  GenerateRandomMatrix(C, m, k); 
  Cbis = C;
  
  MltAdd(alpha, SeldonConjTrans, A, SeldonNoTrans, B, beta, C);
  MltAdd(alpha, SeldonConjTrans, Ad, SeldonNoTrans, B, beta, Cbis);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "MltAdd incorrect" << endl;
      abort();
    }
}

template<class T, class Prop, class Storage, class PropDense, class StorageDense>
void CheckGeneralMatrixProd(Matrix<T, Prop, Storage>& A,
                            Matrix<T, PropDense, StorageDense>& B)
{
  int m = 40, n = m, k = 34, nnz = 200;
  
  Matrix<T, PropDense, StorageDense> C, Cbis, Ad;
  Vector<T> y(n); T val;
  
  // testing Solve
  GenerateRandomMatrix(A, n, n, nnz);
  GenerateRandomMatrix(B, n, k);
  C.Reallocate(n, k);   Cbis.Reallocate(n, k);
  for (int i = 0; i < m; i++)
    {
      GetRandNumber(val);
      A.AddInteraction(i, i, T(1)+val);
    }
  
  Ad.Reallocate(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Ad(i, j) = A(i, j);

  // testing Solve with lower part  
  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < n; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(i, j)*y(j);
          
          y(i) = val/Ad(i, i);
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonLower, SeldonNoTrans, SeldonNonUnit, A, B, C);
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < n; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(i, j)*y(j);
          
          y(i) = val;
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
  
  Solve(SeldonLower, SeldonNoTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val/Ad(i, i);
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
  
  Solve(SeldonLower, SeldonTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val;
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
  
  Solve(SeldonLower, SeldonTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val/conj(Ad(i, i));
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
  
  Solve(SeldonLower, SeldonConjTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val;
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
  
  Solve(SeldonLower, SeldonConjTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }


  // testing Solve with upper part
  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);

      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(i, j)*y(j);
      
          y(i) = val/Ad(i, i);
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
      
  Solve(SeldonUpper, SeldonNoTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);

      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(i, j)*y(j);
      
          y(i) = val;
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
      
  Solve(SeldonUpper, SeldonNoTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val/Ad(i, i);
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonUpper, SeldonTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val;
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonUpper, SeldonTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val/conj(Ad(i, i));
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonUpper, SeldonConjTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val;
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonUpper, SeldonConjTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  // testing Add
  Matrix<T, Prop, Storage> Bs, Cs;
  Matrix<T, PropDense, StorageDense> Bd;
  
  GenerateRandomMatrix(A, m, k, nnz-5);
  GenerateRandomMatrix(Bs, m, k, nnz+10);
  
  Bd.Reallocate(m, k); Ad.Reallocate(m, k);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < k; j++)
      {
        Ad(i, j) = A(i, j);
        Bd(i, j) = Bs(i, j);
      }

  T alpha; GetRandNumber(alpha);
  Add(alpha, Bs, A);
  Add(alpha, Bd, Ad);
  
  if (!EqualMatrix(A, Ad))
    {
      cout << "Add incorrect" << endl;
      abort();
    }
  
  // testing Mlt
  GenerateRandomMatrix(A, m, n, nnz-5);
  GenerateRandomMatrix(Bs, n, k, nnz-5);
  Ad.Reallocate(m, n); Bd.Reallocate(n, k);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      Ad(i, j) = A(i, j);

  for (int i = 0; i < n; i++)
    for (int j = 0; j < k; j++)
      Bd(i, j) = Bs(i, j);
  
  C.Reallocate(m, k);
  Mlt(A, Bs, Cs);
  Mlt(Ad, Bd, C);
  
  if (!EqualMatrix(C, Cs))
    {
      cout << "Mlt incorrect" << endl;
      abort();
    }
  
}


template<class T, class Prop, class Storage, class PropDense, class StorageDense>
void CheckSymmetricMatrixProd(Matrix<T, Prop, Storage>& A,
                              Matrix<T, PropDense, StorageDense>& B)
{
  int m = 40, n = m, k = 34, nnz = 200;
  
  Matrix<T, PropDense, StorageDense> C, Cbis, Ad;
  Vector<T> y(n); T val;
  
  // testing Solve
  GenerateRandomMatrix(A, n, n, nnz);
  GenerateRandomMatrix(B, n, k);
  C.Reallocate(n, k);   Cbis.Reallocate(n, k);
  for (int i = 0; i < m; i++)
    {
      GetRandNumber(val);
      A.AddInteraction(i, i, T(2)+val);
    }
  
  Ad.Reallocate(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      Ad(i, j) = A(i, j);

  // testing Solve with upper part
  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);

      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(i, j)*y(j);
      
          y(i) = val/Ad(i, i);
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
      
  Solve(SeldonNoTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);

      for (int i = m-1; i >= 0; i--)
        {
          val = y(i);
          for (int j = i+1; j < m; j++)
            val -= Ad(i, j)*y(j);
      
          y(i) = val;
        }
      
      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }
      
  Solve(SeldonNoTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val/Ad(i, i);
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= Ad(j, i)*y(j);
          
          y(i) = val;
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val/conj(Ad(i, i));
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonConjTrans, SeldonNonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  for (int p = 0; p < k; p++)
    {
      for (int i = 0; i < n; i++)
        y(i) = B(i, p);
      
      for (int i = 0; i < m; i++)
        {
          val = y(i);
          for (int j = 0; j < i; j++)
            val -= conj(Ad(j, i))*y(j);
          
          y(i) = val;
        }

      for (int i = 0; i < n; i++)
        Cbis(i, p) = y(i);
    }

  Solve(SeldonConjTrans, SeldonUnit, A, B, C);
  
  if (!EqualMatrix(C, Cbis))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
}

template<class T>
void CheckDotProd(Vector<T>& x)
{
  int n = 20;
  Vector<T> y;

  // checking DotProd/DotProdConj
  GenerateRandomVector(x, n);
  GenerateRandomVector(y, n);
  
  T scal = DotProd(x, y);
  T scal_ref; SetComplexZero(scal_ref);
  for (int i = 0; i < n; i++)
    scal_ref += x(i)*y(i);
  
  DISP(scal); DISP(scal_ref);

  scal = DotProdConj(x, y);
  SetComplexZero(scal_ref);
  for (int i = 0; i < n; i++)
    scal_ref += conj(x(i))*y(i);
  
  DISP(scal); DISP(scal_ref);
  
}


int main(int argc, char** argv)
{
  cout.precision(15);
  
  {
    threshold = 1e-5;
    Vector<float> x;
    CheckDotProd(x);
  }

  {
    threshold = 1e-12;
    Vector<double> x;
    CheckDotProd(x);
  }

  {
    threshold = 1e-5;
    Vector<complex<float> > x;
    CheckDotProd(x);
  }

  {
    threshold = 1e-12;
    Vector<complex<double> > x;
    CheckDotProd(x);
  }
  
  {
    threshold = 1e-5;
    Vector<float> x;
    CheckSparseVector(x);
    CheckSparseVectorReal(x);
  }

  {
    threshold = 1e-12;
    Vector<double> x;
    CheckSparseVector(x);
    CheckSparseVectorReal(x);
  }

  {
    threshold = 1e-5;
    Vector<complex<float> > x;
    CheckSparseVector(x);
  }

  {
    threshold = 1e-12;
    Vector<complex<double> > x;
    CheckSparseVector(x);
  }

  {
    threshold = 1e-5;
    Matrix<float, General, RowSparse> A;
    CheckSparseMatrix(A);
    CheckGeneralMatrix(A);
    
    Matrix<float, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad);
    CheckGeneralMatrixProd(A, Ad);
  }

  {
    threshold = 1e-12;
    Matrix<double, General, RowSparse> A;
    CheckSparseMatrix(A);
    CheckGeneralMatrix(A);

    Matrix<double, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad);
    CheckGeneralMatrixProd(A, Ad);
  }

  {
    threshold = 1e-5;
    Matrix<complex<float>, General, RowSparse> A;
    CheckSparseMatrix(A);
    CheckGeneralMatrix(A);

    Matrix<complex<float>, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad); 
    CheckGeneralMatrixProd(A, Ad);   
  }

  {
    threshold = 1e-12;
    Matrix<complex<double>, General, RowSparse> A;
    CheckSparseMatrix(A);
    CheckGeneralMatrix(A);

    Matrix<complex<double>, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad);
    CheckGeneralMatrixProd(A, Ad);
  }

  {
    threshold = 1e-5;
    Matrix<float, Symmetric, RowSymSparse> A;
    CheckSparseMatrix(A);
    CheckSymmetricMatrix(A);

    Matrix<float, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad); threshold = 1e-5;
    CheckSymmetricMatrixProd(A, Ad);
  }

  {
    threshold = 1e-12;
    Matrix<double, Symmetric, RowSymSparse> A;
    CheckSparseMatrix(A);
    CheckSymmetricMatrix(A);

    Matrix<double, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad);
    CheckSymmetricMatrixProd(A, Ad);
  }

  {
    threshold = 1e-5;
    Matrix<complex<float>, Symmetric, RowSymSparse> A;
    CheckSparseMatrix(A);
    CheckSymmetricMatrix(A);

    Matrix<complex<float>, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad); threshold = 1e-5;
    CheckSymmetricMatrixProd(A, Ad);
  }

  {
    threshold = 1e-12;
    Matrix<complex<double>, Symmetric, RowSymSparse> A;
    CheckSparseMatrix(A);
    CheckSymmetricMatrix(A);
    
    Matrix<complex<double>, General, RowMajor> Ad;
    CheckSparseMatrixProd(A, Ad);
    CheckSymmetricMatrixProd(A, Ad);
  }

  cout << "All tests passed successfully" << endl;
  return 0;
}
