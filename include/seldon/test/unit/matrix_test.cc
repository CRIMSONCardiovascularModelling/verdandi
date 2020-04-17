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
void FillRand(Vector<T> & x)
{
  x.FillRand();
  Mlt(1e-9, x);
}

template<class T>
void FillRand(Vector<complex<T> > & x)
{
  for (int i = 0; i < x.GetM(); i++)
    x(i) = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T, class Prop, class Storage>
void FillRand(Matrix<T, Prop, Storage> & A, bool sym)
{
  A.FillRand();
  Mlt(1e-9, A);
}

template<class T, class Prop, class Storage>
void FillRand(Matrix<complex<T>, Prop, Storage> & A, bool sym)
{
  if (sym)
    {
      for (int i = 0; i < A.GetM(); i++)
	for (int j = i; j < A.GetN(); j++)
	  A.Val(i, j) = complex<T>(rand(), rand())/T(RAND_MAX);
    }
  else
    {
      for (int i = 0; i < A.GetM(); i++)
	for (int j = 0; j < A.GetN(); j++)
	  A.Val(i, j) = complex<T>(rand(), rand())/T(RAND_MAX);
    }  
}

template<class T, class Prop, class Storage>
void FillRandTri(Matrix<T, Prop, Storage> & A)
{
  A.FillRand();
  Mlt(1e-9, A);
}

template<class T, class Prop, class Storage>
void FillRandTri(Matrix<complex<T>, Prop, Storage> & A)
{
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j <= i; j++)
      A.Val(i, j) = complex<T>(rand(), rand())/T(RAND_MAX);
  
}

template<class T>
void GetRand(T & x)
{
  x = T(rand())/RAND_MAX;
}

template<class T>
void GetRand(complex<T> & x)
{
  x = complex<T>(rand(), rand())/T(RAND_MAX);
}

template<class T, class Prop1, class Storage1, class Prop2, class Storage2>
bool EqualMatrix(const Matrix<T, Prop1, Storage1>& A,
		 const Matrix<T, Prop2, Storage2>& B)
{
  if ( (A.GetM() != B.GetM()) || (A.GetN() != B.GetN()) )
    return false;
  
  for (int i = 0; i < A.GetM(); i++)
    for (int j = 0; j < A.GetN(); j++)
      if ( abs(A(i, j) - B(i, j) ) > threshold)
	return false;
  
  return true;
}

template<class T, class Prop, class Storage>
void CheckConstOperatorGeneral(const Matrix<T, Prop, Storage>& A, const Vector<T>& x)
{
  int m = A.GetM();
  int n = A.GetN();
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      {
        if (( (A.Val(i, j) != x(i*n + j)) || (A(i, j) != x(i*n + j))
             || (A.Get(i, j) != x(i*n + j)) )
            || isnan(A.Val(i, j)) || isnan(x(i*n+j)) || isnan(A(i, j)) || isnan(A.Get(i, j)))
          {
            cout << "Const operator incorrect" << endl;
            abort();
          }
      }
}


template<class T, class Prop, class Storage>
void CheckOperatorValModif(Matrix<T, Prop, Storage>& A)
{
  int m = 30, n = 20;
  A.Reallocate(m, n);
  Vector<T> x(m*n);
  T zero; SetComplexZero(zero);
  FillRand(x);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A(i, j) = x(i*n + j);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != x(i*n + j)) || isnan(A(i, j)) || isnan(x(i*n+j)))
	{
	  cout << "operator () incorrect" << endl;
	  abort();
	}
  
  A.Fill(zero);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A.Val(i, j) = x(i*n + j);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A.Val(i, j) != x(i*n + j)) || isnan(A.Val(i, j)) || isnan(x(i*n+j)))
	{
	  cout << "Val incorrect" << endl;
	  abort();
	}

  A.Fill(zero);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A.Get(i, j) = x(i*n + j);
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A.Get(i, j) != x(i*n + j)) || isnan(A.Get(i, j)) || isnan(x(i*n+j)))
	{
	  cout << "Get incorrect" << endl;
	  abort();
	}

  A.Fill(zero);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      A.Set(i, j, x(i*n + j));
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != x(i*n + j)) || isnan(A(i, j)) || isnan(x(i*n+j)))
	{
	  cout << "Set incorrect" << endl;
	  abort();
	}
  
  CheckConstOperatorGeneral(A, x);  
}

template<class T, class Prop, class Storage>
void CheckConstOperatorSym(const Matrix<T, Prop, Storage>& A,
                           const Matrix<T, General, RowMajor>& B, bool herm)
{
  int n = A.GetN();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
          {
            DISP(i); DISP(j); DISP(A(i, j)); DISP(B(i, j));
            cout << "Const operator incorrect" << endl;
            abort();
          }
        
        if (j >= i)
          if ((A.Val(i, j) != B(i, j)) || isnan(A.Val(i, j)) || isnan(B(i, j)))
            {
              cout << "Const operator incorrect" << endl;
              abort();
            }        
        
        if (!herm || (j >= i))
          if ((A.Get(i, j) != B(i, j)) || isnan(A.Get(i, j)) || isnan(B(i, j)))
            {
              cout << "Const operator incorrect" << endl;
              abort();
            }
      }
}

template<class T, class Prop, class Storage>
void CheckOperatorVal(Matrix<T, Prop, Storage>& A, bool herm = false)
{
  int n = 4;
  A.Reallocate(n, n);
  Vector<T> x(n*n);
  FillRand(x);
  T zero; SetComplexZero(zero);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      A.Val(i, j) = x(i*n + j);
  
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((A.Val(i, j) != x(i*n + j)) || isnan(A.Val(i, j)) || isnan(x(i*n + j)))
	{
	  cout << "operator Val(i, j) incorrect" << endl;
	  abort();
	}
  
  A.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      A.Val(i, j) = x(i*n + j);
  
  T y;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if (j < i)
	  {
	    if (herm)
	      y = conj(x(j*n + i));
	    else
	      y = x(j*n + i);
	  }
	else
	  y = x(i*n + j);
	
	if ((A(i, j) != y) || isnan(A(i, j)) || isnan(y))
	  {
	    cout << "Operator() incorrect" << endl;
	    abort();
	  }
      }

  Matrix<T, Prop, Storage> B(n, n);
  A.Fill(zero); FillRand(B, true);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      A.Set(i, j, B(i, j));
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
	  {
	    cout << "Set incorrect" << endl;
	    abort();
	  }
      }
  
  if (!herm)
    {
      for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	  A.Get(i, j) = B(i, j);
      
      for (int i = 0; i < n; i++)
	for (int j = 0; j < n; j++)
	  {
	    if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
	      {
		cout << "Get incorrect" << endl;
		abort();
	      }
	  }  
    }
  
  Matrix<T, General, RowMajor> C(n, n);
  FillRand(C, false);
  if (herm)
    {
      for (int i = 0; i < n; i++)
        {
          for (int j = i+1; j < n; j++)
            {
              C(i, j) += conj(C(j, i));
              C(j, i) = conj(C(i, j));
            }
          
          SetComplexReal(real(C(i, i)), C(i, i));
        }
    }
  else
    {
      for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
          {
            C(i, j) += C(j, i);
            C(j, i) = C(i, j);
          }
    }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      A.Set(i, j, C(i, j));
  
  CheckConstOperatorSym(A, C, herm);
}

template<class T, class Prop, class Storage>
void CheckMatrix(Matrix<T, Prop, Storage>& A, bool sym, bool packed, bool herm = false)
{
  if (Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  int m = 32, n = 25;
  if (sym)
    m = n;
  
  A.Reallocate(m, n);
  int size = m*n;
  if ((sym) && (packed))
    size = n*(n+1)/2;
  
  if ( (A.GetM() != m) || (A.GetN() != n) || (A.GetDataSize() != size) || (A.GetSize() != m*n) )
    {
      cout << "GetM, GetN, GetDataSize or GetSize incorrect" << endl;
      abort();
    }
  
  T x;
  GetRand(x);
  A.Fill(x);

  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
	  if ( herm && (i > j) )
	    {
	    }
	  else
	    {
	      cout << "Fill incorrect" <<endl;
	      abort();
	    }
	}
  
  A = x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
	  if ( herm && (i > j) )
	    {
	    }
	  else
	    {
	      cout << "operator = incorrect" <<endl;
	      abort();
	    }
	}
  
  FillRand(A, sym);
  Matrix<T, Prop, Storage> B(A);
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  B.Clear();
  if ( (B.GetM() != 0) || (B.GetN() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  Matrix<T, Prop, Storage> C(A);
  B.SetData(m, n, C.GetData()); C.Nullify();
  if (!EqualMatrix(A, B))
    {
      cout << "SetData incorrect" << endl;
      abort();
    }
  
  B.Clear();
  B = A;
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }

  int m2 = 44, n2 = 33;
  if (sym)
    m2 = n2;
  
  A.Resize(m2, n2);
  if ( (A.GetM() != m2) || (A.GetN() != n2))
    {
      cout << "Resize incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  A.Resize(m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  if (herm)
    SetComplexReal(real(x), x);
  
  A *= x;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      if ((abs(A(i, j) - x*B(i, j)) > threshold) || isnan(abs(A(i, j) - x*B(i, j))))
	{
	  cout << "Operator *= incorrect" << endl;
	  abort();
	}
  
  B.WriteText("toto.dat");
  A.Clear();
  A.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "ReadText/WriteText incorrect" << endl;
      abort();
    }

  B.Write("totob.dat");
  A.Clear();
  A.Read("totob.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
  
  A.Reallocate(n, n);
  A.SetIdentity();
  T zero, one;
  SetComplexOne(one);
  SetComplexZero(zero); 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if (i == j)
	  x = one;
	else
	  x = zero;
	
	if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	  {
	    cout << "SetIdentity incorrect" << endl;
	    abort();
	  }
      }	

  Matrix<T, General, RowMajor> D(m, n);
  FillRand(D, false);
  if (herm)
    {
      for (int i = 0; i < n; i++)
        {
          for (int j = i+1; j < n; j++)
            {
              D(i, j) += conj(D(j, i));
              D(j, i) = conj(D(i, j));
            }
          
          SetComplexReal(real(D(i, i)), D(i, i));
        }
    }
  else if (sym)
    {
      for (int i = 0; i < n; i++)
        for (int j = i+1; j < n; j++)
          {
            D(i, j) += D(j, i);
            D(j, i) = D(i, j);
          }
    }
  
  A.Reallocate(m, n);
  Vector<T> vec_row(n), vec_col(m);
  for (int i = 0; i < m; i++)
    {
      GetRow(D, i, vec_row);
      SetRow(vec_row, i, A);
    }
  
  if (!EqualMatrix(A, D))
    {
      cout << "GetRow/SetRow incorrect" << endl;
      abort();
    }

  A.Fill(zero);
  for (int i = 0; i < n; i++)
    {
      GetCol(D, i, vec_col);
      SetCol(vec_col, i, A);
    }
  
  if (!EqualMatrix(A, D))
    {
      cout << "GetCol/SetCol incorrect" << endl;
      abort();
    }

  for (int i = 0; i < n; i++)
    {
      GetCol(A, i, vec_col);
      for (int j = 0; j < m; j++)
        if ((D(j, i) != vec_col(j)) || isnan(D(j, i)) || isnan(vec_col(j)))
          {
            cout << "GetCol incorrect" << endl;
            abort();
          }
    }

  for (int i = 0; i < m; i++)
    {
      GetRow(A, i, vec_row);
      for (int j = 0; j < n; j++)
        if ((D(i, j) != vec_row(j)) || isnan(D(i, j)) || isnan(vec_row(j)))
          {
            cout << "GetRow incorrect" << endl;
            abort();
          }
    }
  
}

template<class T, class Prop, class Storage>
void CheckConstOperatorTriangUp(const Matrix<T, Prop, Storage>& A,
                                const Matrix<T, General, RowMajor>& B)
{
  int n = A.GetN();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
          {
            cout << "Const operator incorrect" << endl;
            abort();
          }
        
        if (j >= i)
          if (( (A.Val(i, j) != B(i, j)) || (A.Get(i, j) != B(i, j)) )
              || isnan(A.Val(i, j))  || isnan(B(i, j)) || isnan(A.Get(i, j)) || isnan(B(i, j)) )
            {
              cout << "Const operator incorrect" << endl;
              abort();
            }        
      }
}

template<class T, class Prop, class Storage>
void CheckTriangularUp(Matrix<T, Prop, Storage>& A, bool packed)
{
  if (Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }

  int n = 28;
  A.Reallocate(n, n);
  int size = n*(n+1)/2;
  if (! packed)
    size = n*n;
  
  if ( (A.GetM() != n) || (A.GetN() != n) || (A.GetDataSize() != size) || (A.GetSize() != n*n) )
    {
      cout << "GetM, GetN, GetDataSize or GetSize incorrect" << endl;
      abort();
    }

  Vector<T> x2(n*n);
  FillRand(x2);
  T zero; SetComplexZero(zero);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      A.Val(i, j) = x2(i*n + j);
  
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if (((A.Val(i, j) != x2(i*n + j)) || (A.Get(i, j) != x2(i*n + j)))
          || isnan(A.Val(i, j)) || isnan(x2(i*n + j)) || isnan(A.Get(i, j)))
	{
	  cout << "operator Val(i, j) incorrect" << endl;
	  abort();
	}
  
  A.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      A.Val(i, j) = x2(i*n + j);
  
  T y;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        if (i > j)
          y = zero;
        else
          y = x2(i*n + j);
	
	if ((A(i, j) != y) || isnan(A(i, j)) || isnan(y))
	  {
	    cout << "Operator() incorrect" << endl;
	    abort();
	  }
      }

  Matrix<T, General, RowMajor> Cd(n, n);
  FillRand(Cd, false);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < i; j++)
      Cd(i, j) = zero;
  
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      A.Set(i, j, Cd(i, j));
  
  CheckConstOperatorTriangUp(A, Cd);
  
  T x;
  GetRand(x);
  A.Fill(x);

  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
          cout << "Fill incorrect" <<endl;
          abort();
        }
  
  A = x;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
          cout << "operator = incorrect" <<endl;
          abort();
      	}
  
  FillRand(A, true);
  Matrix<T, Prop, Storage> B(A);
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  B.Clear();
  if ( (B.GetM() != 0) || (B.GetN() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  Matrix<T, Prop, Storage> C(A);
  B.SetData(n, n, C.GetData()); C.Nullify();
  if (!EqualMatrix(A, B))
    {
      cout << "SetData incorrect" << endl;
      abort();
    }
  
  B.Clear();
  B = A;
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }

  int n2 = 33;

  A.Resize(n2, n2);
  if ( (A.GetM() != n2) || (A.GetN() != n2))
    {
      cout << "Resize incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  A.Resize(n, n);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  A *= x;
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if ((abs(A(i, j) - x*B(i, j)) > threshold) || isnan(abs(A(i, j) - x*B(i, j))))
	{
	  cout << "Operator *= incorrect" << endl;
	  abort();
	}
  
  B.WriteText("toto.dat");
  A.Clear();
  A.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "ReadText/WriteText incorrect" << endl;
      abort();
    }

  B.Write("totob.dat");
  A.Clear();
  A.Read("totob.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
  
  A.Reallocate(n, n);
  A.SetIdentity();
  T one;
  SetComplexOne(one);
  SetComplexZero(zero); 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if (i == j)
	  x = one;
	else
	  x = zero;
	
	if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	  {
	    cout << "SetIdentity incorrect" << endl;
	    abort();
	  }
      }	
}

template<class T, class Prop, class Storage>
void CheckConstOperatorTriangLo(const Matrix<T, Prop, Storage>& A,
                                const Matrix<T, General, RowMajor>& B)
{
  int n = A.GetN();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        if ((A(i, j) != B(i, j)) || isnan(A(i, j)) || isnan(B(i, j)))
          {
            cout << "Const operator incorrect" << endl;
            abort();
          }
        
        if (i >= j)
          if ( (A.Val(i, j) != B(i, j)) || (A.Get(i, j) != B(i, j)) 
               || isnan(A.Val(i, j)) || isnan(B(i, j)) || isnan(A.Get(i, j)))
            {
              cout << "Const operator incorrect" << endl;
              abort();
            }        
      }
}

template<class T, class Prop, class Storage>
void CheckTriangularLo(Matrix<T, Prop, Storage>& A, bool packed)
{
  if (Storage::Sparse)
    {
      cout << "not a sparse matrix" << endl;
      abort();
    }
  
  int n = 28;
  A.Reallocate(n, n);
  int size = n*(n+1)/2;
  if (! packed)
    size = n*n;
  
  if ( (A.GetM() != n) || (A.GetN() != n) || (A.GetDataSize() != size) || (A.GetSize() != n*n) )
    {
      cout << "GetM, GetN, GetDataSize or GetSize incorrect" << endl;
      abort();
    }

  Vector<T> x2(n*n);
  FillRand(x2);
  T zero; SetComplexZero(zero);
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      A.Val(i, j) = x2(i*n + j);
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if (((A.Val(i, j) != x2(i*n + j)) || (A.Get(i, j) != x2(i*n + j)))
          || isnan(A.Val(i, j))  || isnan(x2(i*n + j)) || isnan(A.Get(i, j)))
	{
	  cout << "operator Val(i, j) incorrect" << endl;
	  abort();
	}
  
  A.Fill(zero);
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      A.Val(i, j) = x2(i*n + j);
  
  T y;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
        if (i < j)
          y = zero;
        else
          y = x2(i*n + j);
	
	if ((A(i, j) != y) || isnan(A(i, j)) || isnan(y))
	  {
	    cout << "Operator() incorrect" << endl;
	    abort();
	  }
      }

  Matrix<T, General, RowMajor> Cd(n, n);
  FillRand(Cd, false);
  for (int i = 0; i < n; i++)
    for (int j = i+1; j < n; j++)
      Cd(i, j) = zero;
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      A.Set(i, j, Cd(i, j));
  
  CheckConstOperatorTriangLo(A, Cd);
  
  T x;
  GetRand(x);
  A.Fill(x);

  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
          cout << "Fill incorrect" <<endl;
          abort();
        }
  
  A = x;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	{
          cout << "operator = incorrect" <<endl;
          abort();
      	}
  
  FillRandTri(A);
  Matrix<T, Prop, Storage> B(A);
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }
  
  B.Clear();
  if ( (B.GetM() != 0) || (B.GetN() != 0) )
    {
      cout << "Clear incorrect" << endl;
      abort();
    }

  Matrix<T, Prop, Storage> C(A);
  B.SetData(n, n, C.GetData()); C.Nullify();
  if (!EqualMatrix(A, B))
    {
      cout << "SetData incorrect" << endl;
      abort();
    }
  
  B.Clear();
  B = A;
  if (!EqualMatrix(A, B))
    {
      cout << "Copy constructor incorrect" << endl;
      abort();
    }

  int n2 = 33;

  A.Resize(n2, n2);
  if ( (A.GetM() != n2) || (A.GetN() != n2))
    {
      cout << "Resize incorrect" << endl;
      abort();
    }
  
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold)
          || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  A.Resize(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if ((abs(A(i, j) - B(i, j)) > threshold) || isnan(abs(A(i, j) - B(i, j))))
	{
	  cout << "Resize incorrect" << endl;
	  abort();
	}
  
  A *= x;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++)
      if ((abs(A(i, j) - x*B(i, j)) > threshold) || isnan(abs(A(i, j) - x*B(i, j))))
	{
	  cout << "Operator *= incorrect" << endl;
	  abort();
	}
  
  B.WriteText("toto.dat");
  A.Clear();
  A.ReadText("toto.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "ReadText/WriteText incorrect" << endl;
      abort();
    }

  B.Write("totob.dat");
  A.Clear();
  A.Read("totob.dat");
  if (!EqualMatrix(A, B))
    {
      cout << "Read/Write incorrect" << endl;
      abort();
    }
  
  A.Reallocate(n, n);
  A.SetIdentity();
  T one;
  SetComplexOne(one);
  SetComplexZero(zero); 
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      {
	if (i == j)
	  x = one;
	else
	  x = zero;
	
	if ((A(i, j) != x) || isnan(A(i, j)) || isnan(x))
	  {
	    cout << "SetIdentity incorrect" << endl;
	    abort();
	  }
      }	
}

int main(int argc, char** argv)
{
  threshold = 1e-12;
  srand(0);

  cout.precision(15);
  
  {
    Matrix<Real_wp, General, RowMajor> A;
    CheckOperatorValModif(A);
    CheckMatrix(A, false, false);
  }

  {
    Matrix<Real_wp, General, ColMajor> A;
    CheckOperatorValModif(A);
    CheckMatrix(A, false, false);
  }
  
  {
    Matrix<Complex_wp, General, RowMajor> A;
    CheckOperatorValModif(A);
    CheckMatrix(A, false, false);
  }

  {
    Matrix<Complex_wp, General, ColMajor> A;
    CheckOperatorValModif(A);
    CheckMatrix(A, false, false);
  }

  {
    Matrix<Real_wp, Symmetric, RowSym> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, false);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSym> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, false);
  }

  {
    Matrix<Real_wp, Symmetric, ColSym> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, false);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSym> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, false);
  }

  {
    Matrix<Real_wp, Symmetric, RowSymPacked> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, true);
  }

  {
    Matrix<Complex_wp, Symmetric, RowSymPacked> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, true);
  }

  {
    Matrix<Real_wp, Symmetric, ColSymPacked> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, true);
  }

  {
    Matrix<Complex_wp, Symmetric, ColSymPacked> A;
    CheckOperatorVal(A);
    CheckMatrix(A, true, true);
  }

  {
    Matrix<Complex_wp, Hermitian, RowHermPacked> A;
    CheckOperatorVal(A, true);
    CheckMatrix(A, true, true, true);
  }

  {
    Matrix<Complex_wp, Hermitian, ColHermPacked> A;
    CheckOperatorVal(A, true);
    CheckMatrix(A, true, true, true);
  }

  {
    Matrix<Complex_wp, Hermitian, RowHerm> A;
    CheckOperatorVal(A, true);
    CheckMatrix(A, true, false, true);
  }

  {
    Matrix<Complex_wp, Hermitian, ColHerm> A;
    CheckOperatorVal(A, true);
    CheckMatrix(A, true, false, true);
  }
  
  {
    Matrix<Real_wp, General, RowUpTriang> A;
    CheckTriangularUp(A, false);
  }

  {
    Matrix<Complex_wp, General, RowUpTriang> A;
    CheckTriangularUp(A, false);
  }

  {
    Matrix<Real_wp, General, ColUpTriang> A;
    CheckTriangularUp(A, false);
  }

  {
    Matrix<Complex_wp, General, ColUpTriang> A;
    CheckTriangularUp(A, false);
  }

  {
    Matrix<Real_wp, General, RowLoTriang> A;
    CheckTriangularLo(A, false);
  }

  {
    Matrix<Complex_wp, General, RowLoTriang> A;
    CheckTriangularLo(A, false);
  }

  {
    Matrix<Real_wp, General, ColLoTriang> A;
    CheckTriangularLo(A, false);
  }

  {
    Matrix<Complex_wp, General, ColLoTriang> A;
    CheckTriangularLo(A, false);
  }
  
  {
    Matrix<Real_wp, General, RowUpTriangPacked> A;
    CheckTriangularUp(A, true);
  }

  {
    Matrix<Complex_wp, General, RowUpTriangPacked> A;
    CheckTriangularUp(A, true);
  }

  {
    Matrix<Real_wp, General, ColUpTriangPacked> A;
    CheckTriangularUp(A, true);
  }

  {
    Matrix<Complex_wp, General, ColUpTriangPacked> A;
    CheckTriangularUp(A, true);
  }

  {
    Matrix<Real_wp, General, RowLoTriangPacked> A;
    CheckTriangularLo(A, true);
  }

  {
    Matrix<Complex_wp, General, RowLoTriangPacked> A;
    CheckTriangularLo(A, true);
  }

  {
    Matrix<Real_wp, General, ColLoTriangPacked> A;
    CheckTriangularLo(A, true);
  }

  {
    Matrix<Complex_wp, General, ColLoTriangPacked> A;
    CheckTriangularLo(A, true);
  }
  
  std::remove("toto.dat");
  std::remove("totob.dat");
  
  cout << "All tests passed successfully" << endl;  
  
  return 0;
}
