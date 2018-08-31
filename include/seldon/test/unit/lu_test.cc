// Copyright (C) 2010 Marc Durufle
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.

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

Real_wp threshold = 1e-10;

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
    x = complex<T>(0, rand())/Real_wp(RAND_MAX);
  else if (type == 1)
    x = complex<T>(rand(), 0)/Real_wp(RAND_MAX);
  else
    x = complex<T>(rand(), rand())/Real_wp(RAND_MAX);
}

template<class T>
void GenerateRandomVector(Vector<T>& x, int n)
{
  x.Reallocate(n);
  for (int i = 0; i < n; i++)
    GetRandNumber(x(i));
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
        A.Set(i, j, x);
      }
}

template<class T, class Prop, class Storage, class Alloc>
bool CheckIdentity(const Matrix<T, Prop, Storage, Alloc>& A)
{
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      {
        if (i == j)
          {
            if (abs(A(i, j) - 1.0) > threshold)
              return false;
          }
        else
          {
            if (abs(A(i, j)) > threshold)
              return false;
          }
      }
  
  return true;
}

template<class T, class T2>
bool EqualVector(const Vector<T>& x, const Vector<T2>& y, Real_wp eps = threshold)
{
  if (x.GetM() != y.GetM())
    return false;
  
  if (Norm2(x) <= eps)
    return false;

  for (int i = 0; i < x.GetM(); i++)
    if (abs(x(i) - y(i)) > eps)
      return false;
  
  return true;
}

template<class T>
bool EqualConditionNumber(const T& c, const T& c2)
{
  if ((c > 3.0*c2) || (c < c2/3.0))
    return false;
  
  return true;
}

template<class T1, class Prop1, class Storage1, class Allocator1,
         class T2, class Prop2, class Storage2, class Allocator2,
         class T3, class Prop3, class Storage3, class Allocator3>
void MltTest(const Matrix<T1, Prop1, Storage1, Allocator1>& A,
             const Matrix<T2, Prop2, Storage2, Allocator2>& B,
             Matrix<T3, Prop3, Storage3, Allocator3>& C)
{
  C.Fill(0);
  T3 val;
  for (int i = 0; i < C.GetM(); i++)
    for (int j = 0; j < C.GetN(); j++)
      {
        SetComplexZero(val);
        for (int k = 0; k < A.GetN(); k++)
          val += A(i, k)*B(k, j);
        
        C.Set(i, j, val);
      }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckGeneralMatrix(Matrix<T, Prop, Storage, Allocator>& M)
{
  int n = 23;
  // testing LU resolution of general real matrices
  Matrix<T, Prop, Storage, Allocator> A(n, n), invA(n, n), Identity(n, n);
  Vector<T> x(n), b(n), y;
  IVect pivot(n);
  T zero; SetComplexZero(zero);
  
  GenerateRandomMatrix(A, n, n);
  GenerateRandomVector(x, n);
  y = x;
  b.Fill(zero);
  Mlt(A, x, b);
  Real_wp anorm = Norm1(A);
  Real_wp anorm_inf = NormInf(A);
    
  /*Matrix<Real_wp> Ar(n, n), Ai(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	Ar(i, j) = real(A(i, j));
	Ai(i, j) = imag(A(i, j));
      }
  Ar.WriteText("Ar.dat");
  Ai.WriteText("Ai.dat");
  */
  Matrix<T, Prop, Storage, Allocator> B(A);
  invA = A;
  x = b;
  GetLU(A, pivot);
  SolveLU(A, pivot, x);
  
  if (! EqualVector(x, y))
    {
      cout << "GetLU/SolveLU incorrect" << endl;
      abort();
    }
  
  x = b;
  SolveLU(SeldonNoTrans, A, pivot, x);
  
  if (! EqualVector(x, y))
    {
      cout << "GetLU/SolveLU incorrect" << endl;
      abort();
    }
  
  Real_wp ferr, berr;
  RefineSolutionLU(invA, A, pivot, x, b, ferr, berr);
  if (! EqualVector(x, y))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
    
  Real_wp eval_cond = ReciprocalConditionNumber(A, pivot, SeldonNorm1, anorm);
  Real_wp eval_cond_inf = ReciprocalConditionNumber(A, pivot, SeldonNormInf, anorm_inf);
  
  A = invA;
  GetInverse(invA);
  Mlt(invA, A, Identity);
  
  //DISP(Norm1(invA));
  //DISP(NormInf(invA));
  
  bool test_inverse = CheckIdentity(Identity);
  if (!test_inverse)
    {
      cout << "GetInverse incorrect" << endl;
      abort();
    }
  
  Real_wp cond_ref = anorm * Norm1(invA);
  if (!EqualConditionNumber(1.0/eval_cond, cond_ref))
    {
      DISP(1.0/eval_cond); DISP(1.0/eval_cond_inf); DISP(cond_ref);
      cout << "ReciprocalConditionNumber-Norm1 incorrect" << endl;
      abort();
    }
  
  Real_wp cond_ref_inf = anorm_inf * NormInf(invA);
  if (!EqualConditionNumber(1.0/eval_cond_inf, cond_ref_inf))
    {
      DISP(1.0/eval_cond_inf); DISP(cond_ref_inf);
      cout << "ReciprocalConditionNumber-NormInf incorrect" << endl;
      abort();
    }
  
  x = y;
  A = B;
  Mlt(SeldonTrans, A, x, b);
    
  x = b;
  GetLU(A, pivot);
  SolveLU(SeldonTrans, A, pivot, x);
    
  if (!EqualVector(x, y))
    {
      cout << "GetLU/SolveLU incorrect" << endl;
      abort();
    }
  
  RefineSolutionLU(SeldonTrans, B, A, pivot, x, b, ferr, berr);
  if (! EqualVector(x, y))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  


  if (IsComplexMatrix(A))
    {
      A = B;
      x = y;
      Mlt(SeldonConjTrans, A, x, b);

      x = b;
      GetLU(A, pivot);
      SolveLU(SeldonConjTrans, A, pivot, x);
      
      if (!EqualVector(x, y))
	{
	  cout << "GetLU/SolveLU incorrect" << endl;
	  abort();
	}
            
      RefineSolutionLU(SeldonConjTrans, B, A, pivot, x, b, ferr, berr);
      if (! EqualVector(x, y))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}  
      
      if ((ferr > 1e3*threshold) || (berr > threshold))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}  

    }
  
  Real_wp row_cond, col_cond, amax;
  VectReal_wp row_scale(n), col_scale(n);
  GetScalingFactors(A, row_scale, col_scale, row_cond, col_cond, amax);
  //DISP(row_scale); DISP(col_scale);
  //DISP(row_cond); DISP(col_cond); DISP(amax);
}

template<class T, class Prop, class Storage, class Allocator>
void CheckSymmetricMatrix(Matrix<T, Prop, Storage, Allocator>& M, bool herm = false)
{
  int n = 24;
  // testing LU resolution of general real matrices
  Matrix<T, Prop, Storage, Allocator> A(n, n), invA(n, n), Identity(n, n);
  Vector<T> x(n), b(n), y;
  IVect pivot(n);
  T zero; SetComplexZero(zero);
  
  GenerateRandomMatrix(A, n, n);
  GenerateRandomVector(x, n);
  if (herm)
    {
      for (int i = 0; i < n; i++)
	A.Set(i, i, real(A(i, i)));
    }
  
  y = x;
  b.Fill(zero);
  Mlt(A, x, b);
  Real_wp anorm = Norm1(A);
  Real_wp anorm_inf = NormInf(A);
  
  invA = A;
  x = b;
  GetLU(A, pivot);
  SolveLU(A, pivot, x);
  
  if (! EqualVector(x, y))
    {
      cout << "GetLU/SolveLU incorrect" << endl;
      abort();
    }
  
  Real_wp ferr, berr;
  RefineSolutionLU(invA, A, pivot, x, b, ferr, berr);
  if (! EqualVector(x, y))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  

  Real_wp eval_cond = ReciprocalConditionNumber(A, pivot, SeldonNorm1, anorm);
  Real_wp eval_cond_inf = ReciprocalConditionNumber(A, pivot, SeldonNormInf, anorm_inf);
  
  A = invA;
  GetInverse(invA);
  MltTest(invA, A, Identity);

  bool test_inverse = CheckIdentity(Identity);
  if (!test_inverse)
    {
      cout << "GetInverse incorrect" << endl;
      abort();
    }

  Real_wp cond_ref = anorm * Norm1(invA);
  if (!EqualConditionNumber(1.0/eval_cond, cond_ref))
    {
      DISP(1.0/eval_cond); DISP(1.0/eval_cond_inf); DISP(cond_ref);
      cout << "ReciprocalConditionNumber-Norm1 incorrect" << endl;
      abort();
    }
  
  Real_wp cond_ref_inf = anorm_inf * NormInf(invA);
  if (!EqualConditionNumber(1.0/eval_cond_inf, cond_ref_inf))
    {
      DISP(1.0/eval_cond_inf); DISP(cond_ref_inf);
      cout << "ReciprocalConditionNumber-NormInf incorrect" << endl;
      abort();
    }  
}

template<class T, class Prop, class Storage, class Allocator>
void GenerateRandomTriangular(Matrix<T, Prop, Storage, Allocator>& A,
			      int m, int n, bool low)
{
  A.Reallocate(m, n);
  typename Matrix<T, Prop, Storage, Allocator>::entry_type x;
  if (low)
    for (int i = 0; i < m; i++)
      for (int j = 0; j <= i; j++)
        {
          GetRandNumber(x);
          A.Set(i, j, x);
        }
  else
    for (int i = 0; i < m; i++)
      for (int j = i; j < n; j++)
        {
          GetRandNumber(x);
          A.Set(i, j, x);
        }

  for (int i = 0; i < m; i++)
    {
      T sum;
      SetComplexZero(sum);
      for (int j = 0; j < n; j++)
	if (j != i)
	  sum += abs(A(i, j));
      
      A.Get(i, i) += sum;
    }
}

template<class T, class Prop, class Storage, class Allocator>
void CheckTriangularMatrix(Matrix<T, Prop, Storage, Allocator>& A, bool low)
{
  int n = 22;
  GenerateRandomTriangular(A, n, n, low);
  
  Real_wp anorm = Norm1(A);
  Real_wp anorm_inf = NormInf(A);

  Real_wp eval_cond = ReciprocalConditionNumber(A, SeldonNorm1);
  Real_wp eval_cond_inf = ReciprocalConditionNumber(A, SeldonNormInf);
  
  // testing SolveLU of triangular systems (Lapack interface)
  Vector<T> X, Y, X0, B;
  GenerateRandomVector(X, n);
  Real_wp ferr, berr;
  
  X0 = X;
  Y = X;
  Mlt(A, X);  B = X;
  SolveLU(A, X);

  if (!EqualVector(X, X0, threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  RefineSolutionLU(A, X, B, ferr, berr);
  if (! EqualVector(X, X0))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  X = X0;
  Y = X;
  Mlt(SeldonNoTrans, SeldonNonUnit, A, X);
  SolveLU(SeldonNoTrans, SeldonNonUnit, A, X);

  if (!EqualVector(X, X0, threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  RefineSolutionLU(SeldonNoTrans, SeldonNonUnit, A, X, B, ferr, berr);
  if (! EqualVector(X, X0))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  

  X = X0;
  Y = X;
  Mlt(SeldonTrans, SeldonNonUnit, A, X); B = X;
  SolveLU(SeldonTrans, SeldonNonUnit, A, X);

  if (!EqualVector(X, X0, threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  RefineSolutionLU(SeldonTrans, SeldonNonUnit, A, X, B, ferr, berr);
  if (! EqualVector(X, X0))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  


  if (IsComplexMatrix(A))
    {
      X = X0;
      Y = X;
      Mlt(SeldonConjTrans, SeldonNonUnit, A, X); B = X;
      SolveLU(SeldonConjTrans, SeldonNonUnit, A, X);
      
      if (!EqualVector(X, X0, threshold))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}

      RefineSolutionLU(SeldonConjTrans, SeldonNonUnit, A, X, B, ferr, berr);
      if (! EqualVector(X, X0))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}  
      
      if ((ferr > 1e3*threshold) || (berr > threshold))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}        
    }

  T one;
  SetComplexOne(one);
  
  X = X0;
  Y = X;
  Mlt(SeldonNoTrans, SeldonUnit, A, X); B = X;
  SolveLU(SeldonNoTrans, SeldonUnit, A, X);
  
  if (!EqualVector(X, X0, threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }

  RefineSolutionLU(SeldonNoTrans, SeldonUnit, A, X, B, ferr, berr);
  if (! EqualVector(X, X0))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }        
  
  X = X0;
  Y = X;
  Mlt(SeldonTrans, SeldonUnit, A, X); B = X;
  SolveLU(SeldonTrans, SeldonUnit, A, X);

  if (!EqualVector(X, X0, 1e3*threshold))
    {
      cout << "Solve incorrect" << endl;
      abort();
    }
  
  RefineSolutionLU(SeldonTrans, SeldonUnit, A, X, B, ferr, berr);
  if (! EqualVector(X, X0))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }  
  
  if ((ferr > 1e3*threshold) || (berr > threshold))
    {
      cout << "RefineSolutionLU incorrect" << endl;
      abort();
    }        
  
  if (IsComplexMatrix(A))
    {
      X = X0;
      Y = X;
      Mlt(SeldonConjTrans, SeldonUnit, A, X); B = X;
      SolveLU(SeldonConjTrans, SeldonUnit, A, X);
      
      if (!EqualVector(X, X0, 1e3*threshold))
	{
	  cout << "Solve incorrect" << endl;
	  abort();
	}
      
      RefineSolutionLU(SeldonConjTrans, SeldonUnit, A, X, B, ferr, berr);
      if (! EqualVector(X, X0))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}  
      
      if ((ferr > 1e3*threshold) || (berr > threshold))
	{
	  cout << "RefineSolutionLU incorrect" << endl;
	  abort();
	}        

    }

  Matrix<T, Prop, Storage, Allocator> invA(A);
  Matrix<T, General, RowMajor, Allocator> Identity(n, n);
  GetInverse(invA);
  Mlt(invA, A, Identity);
  
  bool test_inverse = CheckIdentity(Identity);
  if (!test_inverse)
    {
      cout << "GetInverse incorrect" << endl;
      abort();
    }

  Real_wp cond_ref = anorm * Norm1(invA);
  if (!EqualConditionNumber(1.0/eval_cond, cond_ref))
    {
      DISP(1.0/eval_cond); DISP(1.0/eval_cond_inf); DISP(cond_ref);
      cout << "ReciprocalConditionNumber-Norm1 incorrect" << endl;
      abort();
    }
  
  Real_wp cond_ref_inf = anorm_inf * NormInf(invA);
  if (!EqualConditionNumber(1.0/eval_cond_inf, cond_ref_inf))
    {
      DISP(1.0/eval_cond_inf); DISP(cond_ref_inf);
      cout << "ReciprocalConditionNumber-NormInf incorrect" << endl;
      abort();
    }  
}

int main(int argc, char** argv)
{
  threshold = 1e-10;
  
  //srand(time(NULL));
    
  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckGeneralMatrix(A);
  }
  
  {
    Matrix<Complex_wp, General, RowMajor> A;
    CheckGeneralMatrix(A);
  }

  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckGeneralMatrix(A);
  }
  
  {
    Matrix<Complex_wp, General, ColMajor> A;
    CheckGeneralMatrix(A);
  }
  
  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymPacked> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymPacked> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSym> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    CheckSymmetricMatrix(A);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSym> A;
    CheckSymmetricMatrix(A);
  }
  
  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    CheckSymmetricMatrix(A, true);
  }

  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    CheckSymmetricMatrix(A, true);
  }

  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    CheckSymmetricMatrix(A, true);
  }

  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    CheckSymmetricMatrix(A, true);
  }

  {
    Matrix<Real_wp, General, ColLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Real_wp, General, ColLoTriang> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Real_wp, General, RowLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Real_wp, General, RowLoTriang> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriang> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriangPacked> A;
    CheckTriangularMatrix(A, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriang> A;
    CheckTriangularMatrix(A, true);
  }
  
  {
    Matrix<Real_wp, General, ColUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Real_wp, General, ColUpTriang> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Real_wp, General, RowUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Real_wp, General, RowUpTriang> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriang> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriangPacked> A;
    CheckTriangularMatrix(A, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriang> A;
    CheckTriangularMatrix(A, false);
  }
  
  cout << "All tests passed successfully" << endl;
  
  return 0;
}
