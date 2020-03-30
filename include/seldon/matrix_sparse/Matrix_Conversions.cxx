// Copyright (C) 2003-2011 Marc Duruflé
// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX


#include "Matrix_Conversions.hxx"


namespace Seldon
{

  /*
    From CSR formats to "Matlab" coordinate format.
    index => starting index (usually 0 or 1)
    sym = true => the upper part and lower part are both generated
  */
  
  
  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndRow(j) = i + index;
	  IndCol(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < n; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndCol(j) = i + index;
	  IndRow(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }


  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  if (ptr[i] < ptr[i+1])
            if (ind[ptr[i]] == i)
              nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ... and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }


  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    if (ind[j] == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;

	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }

	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
 	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }

      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndCol(j) = i + index;
	      IndRow(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }

  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  IndRow(nb) = i + index;
	  IndCol(nb) = A.Index(i, j) + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndRow(nb) = A.Index(i, j) + index;
	  IndCol(nb) = i + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }


  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;

	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }

        // Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }


  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;
	
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndCol(nb) = i + index;
	      IndRow(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;

	      if (A.Index(i, j) != i)
		{
		  IndCol(nb) = A.Index(i, j) + index;
		  IndRow(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }
	
	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);

	// ...and the column numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(nb) = A.Index(i, j) + index;
	      IndCol(nb) = i + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }

  
  /*
    From "Matlab" coordinate format to CSR formats.
  */


  //! Conversion from coordinate format to RowSparse.
  /*! Contrary to the other conversion functions
    ConvertMatrix_from_Coordinates, this one accepts duplicates.
    \param[in] IndRow_ row indexes of the non-zero elements.
    \param[in] IndCol_ column indexes of the non-zero elements.
    \param[in] Val values of the non-zero elements.
    \param[out] A matrix defined by \a IndRow, \a IndCol and \a Val.
    \param[in] index index of the first column and of the first row.
  */
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSparse, Allocator3>& A,
				 int index)
  {
    int Nelement = IndRow_.GetLength();

    Vector<int> IndRow(Nelement), IndCol(Nelement);
    
    for (int i = 0; i < Nelement; i++)
      {
	IndRow(i) = IndRow_(i) - index;
	IndCol(i) = IndCol_(i) - index;
      }
    
    IndRow_.Clear();
    IndCol_.Clear();
    
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    
    int m = row_max + 1;
    int n = col_max + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    Sort(IndRow, IndCol, Val);
    
    // Construction of array 'Ptr'.
    Vector<int> Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < Nelement; i++)
      Ptr(IndRow(i)+1)++;
    
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);
    
    A.SetData(m, n, Val, Ptr, IndCol);
  }


#ifndef SWIG

  //! Conversion from coordinate format to ColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSparse, Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(n + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }

    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndRow'
    for (int i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);

    A.SetData(m, n, Val, Ptr, IndRow);
  }


  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;

    int nnz = IndRow.GetM();
    
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);

    A.SetData(m, n, Val, Ptr, IndCol);
  }


  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }

    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);

    // Sorts 'IndRow'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);

    A.SetData(m, n, Val, Ptr, IndRow);
  }

  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */


  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator3>& A,
				 int index)
  {
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }


  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(n);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < n; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Clear(); A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
          // sorting column numbers
          Sort(offset, offset+Ptr(i)-1, IndCol, Val);

          // putting values in A
	  A.ReallocateRow(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
            }

	  offset += Ptr(i);
	}
  }
  
  
  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow_,
				 Vector<int, VectFull, Allocator2>& IndCol_,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymSparse,
				 Allocator3>& A,
				 int index)
  {
    // Assuming that there is no duplicate value.
    if (IndRow_.GetM() <= 0)
      return;

    int nnz = IndRow_.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) = IndRow_(i);
	IndCol(i) = IndCol_(i);
      }
    
    IndRow_.Clear();
    IndCol_.Clear();

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (int i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;

    if (nb_low > 0)
      {
	int nb = 0;
	for (int i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }

	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    Vector<int> Ptr(m);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (int j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }

#endif


  /*
    From Sparse formats to CSC format
  */

  
  //! Conversion from RowSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    
    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;
    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              k++;

            for (int j = ptr_[i]; j < ptr_[i+1]; j++)
              {
                int irow = ind_[j];
                while (k < IndCol.GetM() && IndCol(k) == i
                       && IndRow(k) < irow)
                  k++;

                if (k < IndCol.GetM() && IndCol(k) == i && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
          }
      }
    
    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
	IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        T zero; SetComplexZero(zero);
        for (int i = 0; i <= n; i++)
          {
	    while (k < IndCol.GetM() && IndCol(k) < i)
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
		nb++;
                k++;
              }
	    
	    if (i < n)
	      for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		{
		  int irow = ind_[j];
		  while (k < IndCol.GetM() && IndCol(k) == i
			 && OldInd(k) < irow)
		    {
		      IndRow(nb) = OldInd(k);
		      Val(nb) = OldVal(k);
		      nb++;
		      k++;
		    }
		  
		  if (k < IndCol.GetM() && IndCol(k) == i && OldInd(k) == irow)
		    {
		      // Already existing entry.
		      IndRow(nb) = OldInd(k);
		      Val(nb) = OldVal(k);
		      nb++;
		      k++;
		    }
		  else
		    {
		      // New entry (null).
		      IndRow(nb) = irow;
		      Val(nb) = zero;
		      nb++;
		    }
		}
          }
      }
  }


  //! Conversion from ArrayRowSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);

    // Counting non-zero entries per column.
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    int nb_new_val = 0;
    if (sym_pat)
      {
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        int k = 0;
        for (int i = 0; i < n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              k++;

	    for (int j = 0; j < A.GetRowSize(i); j++)
              {
                int irow = A.Index(i, j);
                while (k < IndCol.GetM() && IndCol(k) == i
                       && IndRow(k) < irow)
                  k++;

                if (k < IndCol.GetM() && IndCol(k) == i && IndRow(k) == irow)
                  // Already existing entry.
                  k++;
                else
                  {
                    // New entry.
                    Ptr(i + 1)++;
                    nb_new_val++;
                  }
              }
	  }
      }
    
    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
    if (sym_pat && (nb_new_val > 0))
      {
        // Changing 'IndRow' and 'Val', and assembling the pattern.
        Vector<Tint, VectFull, Alloc3> OldInd(IndRow);
        Vector<T, VectFull, Alloc4> OldVal(Val);
        IndRow.Reallocate(nnz + nb_new_val);
        Val.Reallocate(nnz + nb_new_val);
        int k = 0, nb = 0;
        T zero; SetComplexZero(zero);
        for (int i = 0; i <= n; i++)
          {
            while (k < IndCol.GetM() && IndCol(k) < i)
              {
                IndRow(nb) = OldInd(k);
                Val(nb) = OldVal(k);
                nb++;
                k++;
              }
	    
	    if (i < n)
	      for (int j = 0; j < A.GetRowSize(i); j++)
		{
		  int irow = A.Index(i, j);
		  while (k < IndCol.GetM() && IndCol(k) == i
			 && OldInd(k) < irow)
		    {
		      IndRow(nb) = OldInd(k);
		      Val(nb) = OldVal(k);
		      nb++;
		      k++;
		    }
		  
		  if (k < IndCol.GetM() && IndCol(k) == i && OldInd(k) == irow)
		    {
		      // Already existing entry.
		      IndRow(nb) = OldInd(k);
		      Val(nb) = OldVal(k);
		      nb++;
		      k++;
		    }
		  else
		    {
		      // New entry (null).
		      IndRow(nb) = irow;
		      Val(nb) = zero;
		      nb++;
		    }
		}
          }
      }
  }

  
  //! Conversion from ColSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    if (!sym_pat)
      {
	// direct conversion
	Ptr.Reallocate(n+1);
	IndRow.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i <= n; i++)
	  Ptr(i) = ptr_[i];
	
	for (int i = 0; i < nnz; i++)
	  {
	    IndRow(i) = ind_[i];
	    Val(i) = data_[i];
	  }
      }
    else
      {
	// Conversion in coordinate format.
	Vector<Tint> IndCol;
	ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
	
	// Sorting with respect to row numbers.
	Sort(IndRow, IndCol, Val);
	
	// Constructing pointer array 'Ptr'.
	Ptr.Reallocate(n + 1);
	Ptr.Fill(0);
	
	// Counting non-zero entries per column.
	for (int i = 0; i < nnz; i++)
	  Ptr(IndRow(i) + 1)++;
	
	int nb_new_val = 0;
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        {
	  int k = 0;
	  for (int i = 0; i < n; i++)
	    {
	      while (k < IndRow.GetM() && IndRow(k) < i)
		k++;
	      
	      for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		{
		  int icol = ind_[j];
		  while (k < IndRow.GetM() && IndRow(k) == i
			 && IndCol(k) < icol)
		    k++;
		  
		  if (k < IndRow.GetM() && IndRow(k) == i && IndCol(k) == icol)
		    // Already existing entry.
		    k++;
		  else
		    {
		      // New entry.
		      Ptr(i + 1)++;
		      nb_new_val++;
		    }
		}
	    }
	}
	
	// Accumulation to get pointer array.
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  Ptr(i + 1) += Ptr(i);
	
	if (nb_new_val > 0)
	  {
	    // Changing 'IndRow' and 'Val', and assembling the pattern.
	    Vector<Tint> OldInd(IndCol);
	    Vector<T, VectFull, Alloc4> OldVal(Val);
	    IndCol.Reallocate(nnz + nb_new_val);
	    Val.Reallocate(nnz + nb_new_val);
	    int k = 0, nb = 0;
            T zero; SetComplexZero(zero);
	    for (int i = 0; i <= n; i++)
	      {
		while (k < IndRow.GetM() && IndRow(k) < i)
		  {
		    IndCol(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		
		if (i < n)
		  for (int j = ptr_[i]; j < ptr_[i+1]; j++)
		    {
		      int icol = ind_[j];
		      while (k < IndRow.GetM() && IndRow(k) == i
			     && OldInd(k) < icol)
			{
			  IndCol(nb) = OldInd(k);
			  Val(nb) = OldVal(k);
			  nb++;
			  k++;
			}
		      
		      if (k < IndRow.GetM() && IndRow(k) == i && OldInd(k) == icol)
			{
			  // Already existing entry.
			  IndCol(nb) = OldInd(k);
			  Val(nb) = OldVal(k);
			  nb++;
			  k++;
			}
		      else
			{
			  // New entry
			  IndCol(nb) = icol;
			  Val(nb) = zero;
			  nb++;
			}
		    }
	      }
	    
	    IndRow.Reallocate(nnz + nb_new_val);
	    for (int i = 0; i < n; i++)
	      for (int j = Ptr(i); j < Ptr(i+1); j++)
		IndRow(j) = i;
	  }
	
	// sorting by columns
	Sort(IndCol, IndRow, Val);
      }    
  }
  
  
  //! Conversion from ArrayColSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with 'nnz' entries.
    int n = A.GetN();
        
    if (!sym_pat)
      {
	// direct conversion
	Ptr.Reallocate(n+1);
	int nnz = 0;
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  {
	    nnz += A.GetColumnSize(i);
	    Ptr(i+1) = nnz;
	  }
	
	IndRow.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(Ptr(i) + j) = A.Index(i, j);
	      Val(Ptr(i) + j) = A.Value(i, j);
	    }
      }
    else
      {
	// Conversion in coordinate format.
	Vector<Tint> IndCol;
	ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
	
	int nnz = IndRow.GetM();
	// Sorting with respect to row numbers.
	Sort(IndRow, IndCol, Val);
	
	// Constructing pointer array 'Ptr'.
	Ptr.Reallocate(n + 1);
	Ptr.Fill(0);
	
	// Counting non-zero entries per column.
	for (int i = 0; i < nnz; i++)
	  Ptr(IndRow(i) + 1)++;
	
	int nb_new_val = 0;
        // Counting entries that are on the symmetrized pattern without being
        // in the original pattern.
        {
	  int k = 0;
	  for (int i = 0; i < n; i++)
	    {
	      while (k < IndRow.GetM() && IndRow(k) < i)
		k++;
	      
	      for (int j = 0; j < A.GetColumnSize(i); j++)
		{
		  int icol = A.Index(i, j);
		  while (k < IndRow.GetM() && IndRow(k) == i
			 && IndCol(k) < icol)
		    k++;
		  
		  if (k < IndRow.GetM() && IndRow(k) == i && IndCol(k) == icol)
		    // Already existing entry.
		    k++;
		  else
		    {
		      // New entry.
		      Ptr(i + 1)++;
                    nb_new_val++;
		    }
		}
	    }
	}
	
	// Accumulation to get pointer array.
	Ptr(0) = 0;
	for (int i = 0; i < n; i++)
	  Ptr(i + 1) += Ptr(i);
	
	if (nb_new_val > 0)
	  {
	    // Changing 'IndRow' and 'Val', and assembling the pattern.
	    Vector<Tint> OldInd(IndCol);
	    Vector<T, VectFull, Alloc4> OldVal(Val);
	    IndCol.Reallocate(nnz + nb_new_val);
	    Val.Reallocate(nnz + nb_new_val);
	    int k = 0, nb = 0;
            T zero; SetComplexZero(zero);
	    for (int i = 0; i <= n; i++)
	      {
		while (k < IndRow.GetM() && IndRow(k) < i)
		  {
		    IndCol(nb) = OldInd(k);
		    Val(nb) = OldVal(k);
		    nb++;
		    k++;
		  }
		
		if (i < n)
		  for (int j = 0; j < A.GetColumnSize(i); j++)
		    {
		      int icol = A.Index(i, j);
		      while (k < IndRow.GetM() && IndRow(k) == i
			     && OldInd(k) < icol)
			{
			  IndCol(nb) = OldInd(k);
			  Val(nb) = OldVal(k);
			  nb++;
			  k++;
			}
		      
		      if (k < IndRow.GetM() && IndRow(k) == i && OldInd(k) == icol)
			{
			  // Already existing entry.
			  IndCol(nb) = OldInd(k);
			  Val(nb) = OldVal(k);
			  nb++;
			  k++;
			}
		      else
			{
			  // New entry
			  IndCol(nb) = icol;
			  Val(nb) = zero;
			  nb++;
			}
		    }
	      }
	    
	    IndRow.Reallocate(nnz + nb_new_val);
	    for (int i = 0; i < n; i++)
	      for (int j = Ptr(i); j < Ptr(i+1); j++)
		IndRow(j) = i;
	  }
	
	// sorting by columns
	Sort(IndCol, IndRow, Val);
      }    
  }
  
  
  //! Conversion from ColSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();    
    int nnz = A.GetDataSize();
    
    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    for (int i = 0; i <= n; i++)
      Ptr(i) = ptr_[i];
    
    for (int i = 0; i < nnz; i++)
      {
	Ind(i) = ind_[i];
	Value(i) = data_[i];
      }
  }

  
  //! Conversion from ColSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint, VectFull, Alloc3> IndCol;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by columns
    Sort(IndCol, IndRow, Value);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayColSymSparse to symmetric CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    int nnz = A.GetDataSize();
    
    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz);
    Value.Reallocate(nnz);
    
    Ptr(0) = 0;
    for (int i = 1; i <= n; i++)
      Ptr(i) = Ptr(i-1) + A.GetColumnSize(i-1);
    
    int nb = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	{
	  Ind(nb) = A.Index(i, j);
	  Value(nb) = A.Value(i, j);
	  nb++;
	}
  }

  
  //! Conversion from ArrayColSymSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint, VectFull, Alloc3> IndCol;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by columns
    Sort(IndCol, IndRow, Value);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from RowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    Vector<Tint, VectFull, Alloc3> IndCol;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);
    
    // sorting by columns
    Sort(IndCol, IndRow, Value);
    
    int n = A.GetN();
    int nnz = A.GetDataSize();
    
    // creating pointer array
    Ptr.Reallocate(n+1);
    Ptr.Fill(0);
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;
    
    for (int i = 0; i < n; i++)
      Ptr(i+1) += Ptr(i);
  }

  
  //! Conversion from RowSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    int n = A.GetN();
    
    Vector<Tint, VectFull, Alloc2> IndCol;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by columns
    Sort(IndCol, IndRow, Value);
    
    Ptr.Reallocate(n+1);
    Ptr.Zero();
    // counting number of non-zero entries
    int nnz = 0;
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	Ptr(IndCol(i) + 1)++;
	nnz++;
      }
    
    // incrementing Ptr
    for (int i = 2; i <= n; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayRowSymSparse to symmetric CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Value, bool sym_pat)
  {
    Vector<Tint, VectFull, Alloc3> IndCol;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);
    
    // sorting by columns
    Sort(IndCol, IndRow, Value);
    
    int n = A.GetN();
    int nnz = A.GetDataSize();
    
    // creating pointer array
    Ptr.Reallocate(n+1);
    Ptr.Fill(0);
    for (int i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;
    
    for (int i = 0; i < n; i++)
      Ptr(i+1) += Ptr(i);
  }

  
  //! Conversion from ArrayRowSymSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& Ind,
                    Vector<T, VectFull, Alloc4>& AllVal, bool sym_pat)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }

    Sort(ind, IndCol, IndRow, Val);

    Ptr.Reallocate(n+1);
    Ind.Reallocate(nnz + ind);
    AllVal.Reallocate(nnz+ind);
    nnz = ind;
    ind = 0;

    int offset = 0; Ptr(0) = 0;
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;

        int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;

	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    Ind(offset+j) = IndRow(ind);
	    AllVal(offset+j) = Val(ind);
            ind++;
	  }

	for (j = 0; j < size_upper; j++)
	  {
	    Ind(offset + size_lower + j) = A.Index(i, j);
	    AllVal(offset + size_lower + j) = A.Value(i, j);
          }

        offset += size_row; Ptr(i+1) = offset;
      }
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowMajor, Allocator>& A,
		  Matrix<T, Prop, RowMajor, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymPacked, Allocator>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColMajor, Allocator>& A,
		  Matrix<T, Prop, ColMajor, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymPacked, Allocator>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator>& A,
		  Matrix<T, Prop, RowSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator>& A,
		  Matrix<T, Prop, RowSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSparse, Allocator>& A,
		  Matrix<T, Prop, ColSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ColSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSymSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayRowSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Allocator>& A,
		  Matrix<T, Prop, ArrayColSymSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;

    General sym;
    ConvertToCSC(mat_array, sym, IndCol, IndRow, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();

    mat_csc.SetData(m, n, Val, IndCol, IndRow);
  }

  
  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }


  //! Conversion from RowSparse to ArrayColSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Alloc1>& A,
		  Matrix<T, Prop, ArrayColSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	B.ReallocateColumn(i, size_col);
	for (int j = Ptr(i); j < Ptr(i+1); j++)
	  {
	    B.Index(i, j-Ptr(i)) = Ind(j);
	    B.Value(i, j-Ptr(i)) = Val(j);
	  }
      }
  }

  
  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Val;

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    General sym;
    ConvertToCSC(mat_array, sym, Ptr, IndRow, Val);

    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }
  
  
  //! Conversion from RowSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayRowSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSymSparse to ColSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, ColSymSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from RowSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayRowSymSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T1, VectFull, Allocator1> AllVal;

    int n = A.GetM();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, AllVal);

    B.SetData(n, n, AllVal, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSymSparse to ColSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    Vector<int> Ptr;
    Vector<int> Ind;
    Vector<T, VectFull, Alloc2> Val;

    int m = A.GetM(), n = A.GetN();
    General sym;
    ConvertToCSC(A, sym, Ptr, Ind, Val);

    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  /*
    From Sparse formats to CSR format
  */
  
  
  //! Conversion from RowSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    for (int i = 0; i <= m; i++)
      Ptr(i) = ptr_[i];
    
    for (int i = 0; i < nnz; i++)
      {
        IndCol(i) = ind_[i];
        Value(i) = data_[i];
      }
  }
  
  
  //! Conversion from ColSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (int i = 0; i < nnz; i++)
      Ptr(ind_[i])++;

    // Incrementing in order to get the row indexes.
    int increment = 0, size, num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    
    // Last index.
    Ptr(m) = increment;
    
    // 'Offset' will be used to get current positions of new entries.
    Vector<Tint, VectFull, Alloc2> Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (int i = ptr_[j]; i < ptr_[j + 1]; i++)
	{
	  num_row = ind_[i];
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = data_[i];
	  Offset(num_row)++;
	}
  }
  
  
  //! Conversion from ArrayColSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int n = A.GetN();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    // Computation of the indexes of the beginning of rows.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    // Counting the number of entries per row.
    for (int i = 0; i < n; i++)
      for (int j = 0; j < A.GetColumnSize(i); j++)
	Ptr(A.Index(i, j))++;
    
    // Incrementing in order to get the row indexes.
    int increment = 0, size, num_row;
    for (int i = 0; i < m; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    
    // Last index.
    Ptr(m) = increment;
    
    // 'Offset' will be used to get current positions of new entries.
    Vector<Tint, VectFull, Alloc2> Offset(Ptr);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);

    // Loop over the columns.
    for (int j = 0; j < n; j++)
      for (int i = 0; i < A.GetColumnSize(j); i++)
	{
	  num_row = A.Index(j, i);
	  IndCol(Offset(num_row)) = j;
	  Value(Offset(num_row)) = A.Value(j, i);
	  Offset(num_row)++;
	}
  }
  
  
  //! Conversion from ArrayRowSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    int  nnz = A.GetDataSize();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0; Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  {
	    IndCol(nnz + j) = A.Index(i, j);
	    Value(nnz + j) = A.Value(i, j);
	  }
	
	nnz += A.GetRowSize(i);
	Ptr(i+1) = nnz;
      }
  }


  //! Conversion from ArrayRowSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    int m = A.GetM();
    if (m <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    int  nnz = 0;
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) >= i)
	  nnz++;
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(nnz);
    Value.Reallocate(nnz);
    nnz = 0; Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int j = 0; j < A.GetRowSize(i); j++)
	  if (A.Index(i, j) >= i)
	    {
	      IndCol(nnz) = A.Index(i, j);
	      Value(nnz) = A.Value(i, j);
	      nnz++;
	    }
	
	Ptr(i+1) = nnz;
      }
  }

  
  //! Conversion from ColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }
  
  
  //! Conversion from ColSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayColSymSparse to symmetric CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;

    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, false);

    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }
  
  
  //! Conversion from ArrayColSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
  }
  
  
  //! Conversion from RowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = ptr_[i]; k < ptr_[i+1]; k++)
	  {
	    IndCol(ind) = ind_[k];
	    Val(ind) = data_[k];
	    ind++;
	  }
	
	IndRow(i + 1) = ind;
      }
  }

  
  //! Conversion from RowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }
  
  
  //! Conversion from ArrayRowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& IndRow,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Number of rows and non-zero entries.
    int nnz = A.GetDataSize();
    int m = A.GetM();

    // Allocation of arrays for CSR format.
    Val.Reallocate(nnz);
    IndRow.Reallocate(m + 1);
    IndCol.Reallocate(nnz);

    int ind = 0;
    IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
	for (int k = 0; k < A.GetRowSize(i); k++)
	  {
	    IndCol(ind) = A.Index(i, k);
	    Val(ind) = A.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
  }

  
  //! Conversion from ArrayRowSymSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value)
  {
    Vector<Tint, VectFull, Alloc3> IndRow;
    
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);
    
    int m = A.GetM();
    Ptr.Reallocate(m+1);
    Ptr.Zero();

    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndRow(i) + 1)++;
    
    // incrementing Ptr
    for (int i = 2; i <= m; i++)
      Ptr(i) += Ptr(i-1);
    
  }

  
  //! Conversion from column-oriented sparse to row-oriented sparse.
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T1, VectFull, Alloc2> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }

  
  //! Conversion from ArrayColSparse to RowSparse
  /*!
    \param[in] A matrix to be converted.
    \param[out] B converted matrix.
  */
  template<class T1, class T2, class Prop1, class Prop2,
           class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayColSparse, Alloc1>& A,
		  Matrix<T2, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T1, VectFull, Alloc2> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }
  
  
  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;

    General sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }


  //! Conversion from ArrayRowSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;
    
    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }

  
  //! Conversion from ColSymSparse to RowSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ColSymSparse to RowSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    General unsym;
    ConvertToCSR(A, unsym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayColSymSparse to RowSymSparse
  template<class T, class Prop, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop, RowSymSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from ArrayColSymSparse to RowSparse
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, ArrayColSymSparse, Alloc1>& A,
		  Matrix<T, Prop2, RowSparse, Alloc2>& B)
  {
    Vector<int> Ptr, Ind;
    Vector<T, VectFull, Alloc2> Value;
    
    General unsym;
    ConvertToCSR(A, unsym, Ptr, Ind, Value);

    // creating the matrix
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, Ind);
  }


  //! Conversion from RowSymSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;

    General unsym;
    ConvertToCSR(mat_array, unsym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }
  
  
  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;

    Symmetric sym;
    ConvertToCSR(mat_array, sym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  //! Conversion from ArrayRowSymSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& mat_array,
		  Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<T1, VectFull, Allocator1> Val;
    Vector<int> IndRow;
    Vector<int> IndCol;

    General unsym;
    ConvertToCSR(mat_array, unsym, IndRow, IndCol, Val);

    int m = mat_array.GetM();
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }


  /******************************
   * From Sparse to ArraySparse *
   ******************************/
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<int> Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<int> Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	B.ReallocateColumn(i, size_col);
	for (int j = 0; j < size_col; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T1, VectFull, Allocator1> Val(nnz);
    int* indA = A.GetInd();
    int* ptrA = A.GetPtr();
    T0* dataA = A.GetData();
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = ptrA[i]; j < ptrA[i+1]; j++)
	if (indA[j] != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = indA[j];
	    Val(ind) = dataA[j];
	    ind++;
	  }
    
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
        
	int size_lower = ind - first_index;
	int size_upper = ptrA[i+1] - ptrA[i];
	int size_row = size_lower + size_upper;
	B.ResizeRow(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = indA[ptrA[i]+j];
	    B.Value(i, size_lower + j) = dataA[ptrA[i]+j];
	  }
        
	B.AssembleRow(i);
      }
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, RowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
	int size_col = ptr_[i+1] - ptr_[i];
	B.ReallocateColumn(i, size_col);
	for (int j = 0; j < size_col; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }    
  }

  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {    
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    CopyMatrix(Acsc, A);
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }
  
  
  /***********************************
   * From ArraySparse to ArraySparse *
   ***********************************/
  
  
  //! From ArrayRowSymSparse to ArrayRowSymSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }
  
  
  //! From ArrayRowSparse to ArrayRowSparse (T0 and T1 different)
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRowSize(i);
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = A.Index(i, j);
            B.Value(i, j) = A.Value(i, j);
          }
      }
  }
  
  
  //! upper part of A is used to obtain a symmetric matrix
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    B.Reallocate(A.GetM(), A.GetN());
    for (int i = 0; i < A.GetM(); i++)
      {
	int k = 0;
	while ( (k < A.GetRowSize(i)) && (A.Index(i, k) < i))
	  k++;
	
	if (k < A.GetRowSize(i))
	  {
	    int size_row = A.GetRowSize(i) - k;
	    B.ReallocateRow(i, size_row);
	    for (int j = k; j < A.GetRowSize(i); j++)
	      {
		B.Index(i, j-k) = A.Index(i, j);
		B.Value(i, j-k) = A.Value(i, j);
	      }
	  }
      }
  }

  
  //! conversion from ArrayColSymSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int n = A.GetM();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    Vector<int> Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    B.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      {
	int size_row = Ptr(i+1) - Ptr(i);
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = Ind(Ptr(i) + j);
	    B.Value(i, j) = Value(Ptr(i) + j);
	  }
      }    
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& Acsc,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {    
    int m = Acsc.GetM();
    int n = Acsc.GetN();
    if (n <= 0)
      {
	B.Clear();
	return;
      }
    
    // conversion to RowSparse
    Matrix<T0, Prop0, RowSparse, Allocator0> A;
    CopyMatrix(Acsc, A);
    
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
	int size_row = ptr_[i+1] - ptr_[i];
	B.ReallocateRow(i, size_row);
	for (int j = 0; j < size_row; j++)
	  {
	    B.Index(i, j) = ind_[ptr_[i] + j];
	    B.Value(i, j) = data_[ptr_[i] + j];
	  }
      }
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int> IndRow(nnz),IndCol(nnz);
    Vector<T1, VectFull, Allocator1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeRow(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleRow(i);
      }
  }


  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    Vector<int> IndRow(nnz), IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;

    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeColumn(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleColumn(i);
      }
  }

  
  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void CopyMatrix(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
		  Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i;

    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();

    // Conversion in coordinate format.
    Vector<T1> Val;
    Vector<int> IndRow, IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);

    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);

    // Constructing pointer array 'Ptr'.
    Vector<int> Ptr(n + 1);
    Ptr.Zero();
    
    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);

    // we fill matrix B
    B.Reallocate(A.GetM(), n);
    for (i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	if (size_col > 0)
	  {
	    B.ReallocateColumn(i, size_col);
	    for (int j = Ptr(i); j < Ptr(i+1); j++)
	      {
		B.Index(i, j-Ptr(i)) = IndRow(j);
		B.Value(i, j-Ptr(i)) = Val(j);
	      }
	  }
      }
  }

  
  /***********************
   * GetSymmetricPattern *
   ***********************/


  //! Returns pattern of A + A' in CSR format
  /*!
    From a sparse matrix, we compute the pattern of A + A'
    so that this pattern is symmetric even if A is non-symmetric    
   */
  template<class T, class Prop, class Storage, class Allocator,
           class Tint, class Allocator2, class Allocator3>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Vector<Tint, VectFull, Allocator2>& Ptr,
                           Vector<Tint, VectFull, Allocator3>& Ind)
  {
    typedef typename Matrix<T, Prop, Storage, Allocator>::entry_type T0;
    int n = A.GetM();

    // Converting to coordinates.
    Vector<int> IndRow, IndCol;
    Vector<T0> Value;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Value);

    // clearing values
    Value.Clear();

    // Sorting columns too.
    Vector<int> IndRow2, IndCol2, Index(2*n);
    IndRow2 = IndRow;
    IndCol2 = IndCol;
    Sort(IndCol2.GetM(), IndCol2, IndRow2);

    Tint max_nnz = 0;
    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndRow(i) <= IndCol(i))
        max_nnz++;

    for (int i = 0; i < IndRow.GetM(); i++)
      if (IndCol2(i) <= IndRow2(i))
        max_nnz++;

    // then symmetrization of pattern and conversion to csr.
    Ptr.Reallocate(n+1);
    Ind.Reallocate(max_nnz);
    Tint j_end = 0;
    int size_row = 0;
    Tint j2_end = 0;
    Ptr(0) = 0;
    for (int i = 0; i < A.GetM(); i++)
      {
        size_row = 0;
        // We retrieve column numbers.
        while ( (j_end < IndRow.GetM()) && (IndRow(j_end) == i))
          {
            if (IndRow(j_end) <= IndCol(j_end))
              {
                Index(size_row) = IndCol(j_end);
                size_row++;
              }

            j_end++;
          }

        while ( (j2_end < IndCol2.GetM()) && (IndCol2(j2_end) == i))
          {
            if (IndCol2(j2_end) <= IndRow2(j2_end))
              {
                Index(size_row) = IndRow2(j2_end);
                size_row++;
              }

            j2_end++;
          }

        // Sorting indexes.
        Assemble(size_row, Index);

        // Updating Ptr, Ind.
        for (int j = 0; j < size_row; j++)
	  Ind(Ptr(i) + j) = Index(j);

        Ptr(i+1) = Ptr(i) + size_row;
      }

    IndRow2.Clear(); IndCol2.Clear();
    IndRow.Clear(); IndCol.Clear();
    Ind.Resize(Ptr(n));
  }


  template<class T, class Prop, class Storage, class Allocator, class AllocI>
  void GetSymmetricPattern(const Matrix<T, Prop, Storage, Allocator>& A,
                           Matrix<int, Symmetric, RowSymSparse, AllocI>& B)
  {
    Vector<int> Ptr, Ind;

    GetSymmetricPattern(A, Ptr, Ind);

    int n = A.GetM();
    Vector<int, VectFull, AllocI> Val(Ptr(n));
    // We put Ptr and Ind into the matrix B.
    B.SetData(n, n, Val, Ptr, Ind);
  }

  
  /*****************************************************
   * Conversion from sparse matrices to dense matrices *
   *****************************************************/
  
  
  //! conversion from RowSparse to RowMajor
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSparse, Allocator1>& A,
		  Matrix<T, Prop, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = ptr[i]; j < ptr[i+1]; j++)
        B(i, ind[j]) = data[j];
    
  }

  
  //! conversion from ArrayRowSparse to RowMajor
  template<class T1, class T2, class Prop1, 
	   class Prop2, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& A,
		  Matrix<T2, Prop2, RowMajor, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    
    B.Reallocate(m, n);
    T2 zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
        B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from RowSymSparse to RowSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, RowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* data = A.GetData();
    
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = ptr[i]; j < ptr[i+1]; j++)
        B(i, ind[j]) = data[j];
    
  }

  
  //! conversion from ArrayRowSymSparse to RowSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to ColSymPacked
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSymPacked, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to RowSym
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, RowSym, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B.Val(i, A.Index(i, j)) = A.Value(i, j);
    
  }


  //! conversion from ArrayRowSymSparse to ColSym
  template<class T, class Prop, class Allocator1, class Allocator2>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymSparse, Allocator1>& A,
		  Matrix<T, Prop, ColSym, Allocator2>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    T zero;
    SetComplexZero(zero);
    B.Fill(zero);
    for (int i = 0; i < m; i++)
      for (int j = 0; j < A.GetRowSize(i); j++)
	B.Val(i, A.Index(i, j)) = A.Value(i, j);
    
  }
  

  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, PETScMPIDense, Alloc1>& A,
		  Matrix<T, Prop2, RowMajor, Alloc2>& B)
  {
    int m, n;
    MatGetLocalSize(A.GetPetscMatrix(), &m, &n);
    n = A.GetN();
    B.Reallocate(m, n);
    T *local_a;
    MatGetArray(A.GetPetscMatrix(), &local_a);
    for (int i = 0; i < m; i++)
      for(int j = 0; j < n; j++)
        B(i, j) = local_a[i + j * m];
    MatRestoreArray(A.GetPetscMatrix(), &local_a);
  }


  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void CopyMatrix(const Matrix<T, Prop1, RowMajor, Alloc1>& A,
		  Matrix<T, Prop2, PETScMPIDense, Alloc2>& B)
  {
    T *local_data;
    MatGetArray(B.GetPetscMatrix(), &local_data);
    int mlocal, nlocal;
    MatGetLocalSize(B.GetPetscMatrix(), &mlocal, &nlocal);
    Matrix<T, Prop1, ColMajor, Alloc1> local_D;
    local_D.SetData(mlocal, B.GetN(), local_data);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
        local_D(i, j) = A(i, j);
    local_D.Nullify();
    MatRestoreArray(B.GetPetscMatrix(), &local_data);
  }

  
  /*****************************************************
   * Conversion from dense matrices to sparse matrices *
   *****************************************************/
  
  
  //! conversion from RowSymPacked to RowSymSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, Symmetric, RowSymPacked>& A,
                       Matrix<T, Symmetric, RowSymSparse>& B,
		       const T& threshold)
  {
    int nnz = 0;
    int n = A.GetM();
    for (int i = 0; i < n; i++)
      for (int j = i; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;
    
    Vector<int> IndCol(nnz), IndRow(n+1); 
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < n; i++)
      {
        for (int j = i; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }
        
        IndRow(i+1) = nnz;
      }
    
    B.SetData(n, n, Value, IndRow, IndCol);
    
  }
  
  
  //! conversion from RowMajor to ArrayRowSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, ArrayRowSparse>& B,
		       const T& threshold)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            size_row++;
        
        B.ReallocateRow(i, size_row);
        
        size_row = 0;
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              B.Index(i, size_row) = j;
              B.Value(i, size_row) = A(i, j);
              size_row++;
            }
      }
  }
  
  
  //! conversion from RowMajor to RowSparse
  template<class T>
  void ConvertToSparse(const Matrix<T, General, RowMajor>& A,
                       Matrix<T, General, RowSparse>& B,
		       const T& threshold)
  {
    int nnz = 0;
    int m = A.GetM();
    int n = A.GetN();
    for (int i = 0; i < m; i++)
      for (int j = 0; j < n; j++)
        if (abs(A(i, j)) > threshold)
          nnz++;
    
    Vector<int> IndCol(nnz), IndRow(m+1); 
    Vector<T> Value(nnz);
    nnz = 0; IndRow(0) = 0;
    for (int i = 0; i < m; i++)
      {
        for (int j = 0; j < n; j++)
          if (abs(A(i, j)) > threshold)
            {
              IndCol(nnz) = j;
              Value(nnz) = A(i, j);
              nnz++;
            }
        
        IndRow(i+1) = nnz;
      }
    
    B.SetData(m, n, Value, IndRow, IndCol);
    
  }
    
} // namespace Seldon.

#define SELDON_FILE_MATRIX_CONVERSIONS_CXX
#endif
