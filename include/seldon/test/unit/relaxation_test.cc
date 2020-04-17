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

#include "Seldon.hxx"
#include "SeldonComplexMatrix.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

typedef double Real_wp;
typedef complex<double> Complex_wp;
typedef Vector<double> VectReal_wp;

Real_wp threshold = 1e-12;

template<class T1, class Prop, class Storage, class Allocator,
	 class T2, class T3>
void SorTest(const Matrix<T1, Prop, Storage, Allocator>& A,
	     Vector<T2>& X, const Vector<T2>& B,
	     const T3& omega, int iter, int type)
{
  // forward sweep
  T2 val, zero, one;
  SetComplexZero(zero);
  SetComplexOne(one);
  int n = A.GetM();
  if (type%2 == 0)
    for (int num = 0; num < iter; num++)
      {
        for (int i = 0; i < n; i++)
          {
            val = B(i);
            for (int j = 0; j < n; j++)
              if (i != j)
                val -= A(i, j)*X(j);
            
            X(i) = (one - omega) * X(i) + omega * val / A(i, i);
          }
      }
  
  if (type%3 == 0)
    for (int num = 0; num < iter; num++)
      for (int i = n-1; i >= 0; i--)
	{
	  val = B(i);
	  for (int j = 0; j < n; j++)
	    if (i != j)
	      val -= A(i, j)*X(j);
	  
	  X(i) = (one - omega) * X(i) + omega * val / A(i, i);
	}
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y)
{
  if ( (x.GetM() == 0) || (x.GetM() != y.GetM()))
    return false;
  
  Real_wp cte = Norm2(x);
  if (cte < threshold)
    return false;
  
  for (int i = 0; i < x.GetM(); i++)
    if ((abs(x(i) - y(i)) > threshold*cte) || isnan(abs(x(i) - y(i))))
      {
        DISP(i); DISP(x(i)); DISP(y(i)); DISP(x(i)-y(i));
        return false;
      }
  
  return true;
}

template<class T>
void GetRand(T& x)
{
  x = rand();
}

template<class T>
void GetRand(complex<T>& x)
{
  int type = rand()%3;
  switch (type)
    {
    case 0 : x = complex<T>(rand(), 0); return;
    case 1 : x = complex<T>(0, rand()); return;
    case 2 : x = complex<T>(rand(), rand()); return;
    }
  
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSorFunction(Matrix<T, Prop, Storage, Allocator>& A)
{
  typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type Tcplx;
  typedef typename ClassComplexType<T>::Treal Treal;
  Vector<Tcplx> x, b, y;

  int nb_iter = 2;
  Treal omega(0.5);
  int n = 37, nnz = 200;
  Matrix<Tcplx, General, ColSparse> Acsr(n, n);
  Acsr.FillRand(nnz);
  for (int i = 0; i < Acsr.GetDataSize(); i++)
    GetRand(Acsr.GetData()[i]);
  
  for (int i = 0; i < n; i++)
    {
      T val; SetComplexReal(Treal(2)*Treal(RAND_MAX) + Treal(rand()), val);
      Acsr.Get(i, i) += val;
    }
  
  Mlt(Treal(1)/Treal(RAND_MAX), Acsr);

  x.Reallocate(n); 
  y.Reallocate(n); 
  b.Reallocate(n);
  b.FillRand(); Mlt(Treal(1)/Treal(RAND_MAX), b);
  x.Fill(0); y.Fill(0);
  
  Copy(Acsr, A);
  SOR(A, x, b, omega, nb_iter, 2);
  SorTest(A, y, b, omega, nb_iter, 2);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  x.Fill(0); y.Fill(0);
  SOR(A, x, b, omega, nb_iter, 3);
  SorTest(A, y, b, omega, nb_iter, 3);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  x.Fill(0); y.Fill(0);
  SOR(A, x, b, omega, nb_iter, 6);
  SorTest(A, y, b, omega, nb_iter, 6);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }
}


template<class T, class Prop, class Storage, class Allocator>
void CheckSorSymFunction(Matrix<T, Prop, Storage, Allocator>& A)
{
  typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type Tcplx;
  typedef typename ClassComplexType<T>::Treal Treal;
  Vector<Tcplx> x, b, y;

  int nb_iter = 2;
  Treal omega(0.5);
  int n = 29, nnz = 100;
  Matrix<Tcplx, Symmetric, RowSymSparse> Acsr(n, n);
  for (int k = 0; k < nnz; k++)
    {
      int i = rand()%n;
      int j = rand()%n;
      GetRand(Acsr.Get(i, j));
    }

  for (int i = 0; i < Acsr.GetDataSize(); i++)
    GetRand(Acsr.GetData()[i]);

  for (int i = 0; i < n; i++)
    {
      T val; SetComplexReal(Treal(2)*Treal(RAND_MAX) + Treal(rand()), val);
      Acsr.Get(i, i) += val;
    }
  
  Mlt(Treal(1)/Treal(RAND_MAX), Acsr);
  
  x.Reallocate(n); 
  y.Reallocate(n); 
  b.Reallocate(n);
  b.FillRand(); Mlt(Treal(1)/Treal(RAND_MAX), b);
  x.Fill(0); y.Fill(0);
  
  Copy(Acsr, A);
  SOR(A, x, b, omega, nb_iter, 2);
  SorTest(A, y, b, omega, nb_iter, 2);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }
  return;
  
  x.Fill(0); y.Fill(0);
  SOR(A, x, b, omega, nb_iter, 3);
  SorTest(A, y, b, omega, nb_iter, 3);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }

  x.Fill(0); y.Fill(0);
  SOR(A, x, b, omega, nb_iter, 6);
  SorTest(A, y, b, omega, nb_iter, 6);
  
  if (!EqualVector(x, y))
    {
      cout << "SOR incorrect" << endl;
      abort();
    }  

}


int main(int argc, char** argv)
{
  threshold = 1e-12;
  
  //srand(time(NULL));
  //Complex_wp x(2, 5), y(3, -7);
  //DISP(x*y);
  //DISP(x*y - Complex_wp(41, 1));

  //Complex_wp a(3, 0), b(0, -5);
  //DISP(a*x); DISP(a*y); DISP(a*b); DISP(b*x); DISP(b*y);
  
  //DISP(x/y); DISP(a/x); DISP(a/y); DISP(a/b); DISP(b/x); DISP(b/y);
  //DISP(x/a); DISP(y/a); DISP(b/a); DISP(x/b); DISP(y/b);
  
  //Real_wp ar(3);
  //DISP(ar/x); DISP(ar/y); DISP(ar/a); DISP(ar/b);
  //DISP(x/ar); DISP(y/ar); DISP(b/ar); DISP(a/ar);

  //return 0;
  
  {
    Matrix<Real_wp, General, RowSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Real_wp, General, ColSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Real_wp, General, ArrayRowSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Real_wp, General, ArrayColSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, RowSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ColSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ArrayRowSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ArrayColSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, RowComplexSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ColComplexSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ArrayRowComplexSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Complex_wp, General, ArrayColComplexSparse> A;
    CheckSorFunction(A);
  }

  {
    Matrix<Real_wp, Symmetric, RowSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Real_wp, Symmetric, ColSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Real_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Real_wp, Symmetric, ArrayColSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayColSymSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymComplexSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayRowSymComplexSparse> A;
    CheckSorSymFunction(A);
  }
  
  {
    Matrix<Complex_wp, Symmetric, ColSymComplexSparse> A;
    CheckSorSymFunction(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ArrayColSymComplexSparse> A;
    CheckSorSymFunction(A);
  }
  
  cout << "All tests passed successfully" << endl;  

  return 0;
}
