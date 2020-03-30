// Copyright (C) 2003-2011 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_CXX


#include "Matrix_ComplexConversions.hxx"

/*
  Same functions as in Matrix_Conversions.cxx 
  for complex matrices (RowComplexSparse, etc)
 */

namespace Seldon
{
  
  //! Conversion from RowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    IVect col; Vector<T> value;
    for (int i = 0; i < m; i++)
      {
        int nb_r = real_ptr[i+1] - real_ptr[i];
        int nb_i = imag_ptr[i+1] - imag_ptr[i];
        int size_row = nb_r + nb_i;
        if (size_row > col.GetM())
          {
            col.Reallocate(size_row);
            value.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
            col(nb) = real_ind[j] + index;
            value(nb) = T(real_data[j], 0);
            nb++;
          }

        for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
          {
            col(nb) = imag_ind[j] + index;
            value(nb) = T(0, imag_data[j]);
            nb++;
          }
        
        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = index + i;
            IndCol(nnz + j) = col(j);
            Val(nnz + j) = value(j);
          }
        
        nnz += nb;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int n = A.GetN();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    IVect col; Vector<T> value;
    for (int i = 0; i < n; i++)
      {
        int nb_r = real_ptr[i+1] - real_ptr[i];
        int nb_i = imag_ptr[i+1] - imag_ptr[i];
        int size_col = nb_r + nb_i;
        if (size_col > col.GetM())
          {
            col.Reallocate(size_col);
            value.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
          {
            col(nb) = real_ind[j] + index;
            value(nb) = T(real_data[j], 0);
            nb++;
          }

        for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
          {
            col(nb) = imag_ind[j] + index;
            value(nb) = T(0, imag_data[j]);
            nb++;
          }
        
        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = col(j);
            IndCol(nnz + j) = index + i;
            Val(nnz + j) = value(j);
          }
        
        nnz += nb;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    
    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect col; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_row = nb_r + nb_i;
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }
            
            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                col(nb) = real_ind[j];
                value(nb) = T(real_data[j], 0);
                nb++;
              }
            
            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                col(nb) = imag_ind[j];
                value(nb) = T(0, imag_data[j]);
                nb++;
              }
            
            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = col(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;
                
                if (col(j) != i)
                  {
                    IndRow(nnz) = col(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(col(j))++;
                    nnz++;
                  }
              }
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
        
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
        nnz = 0;
        IVect col; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_row = nb_r + nb_i;
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }
            
            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                col(nb) = real_ind[j] + index;
                value(nb) = T(real_data[j], 0);
                nb++;
              }
            
            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                col(nb) = imag_ind[j] + index;
                value(nb) = T(0, imag_data[j]);
                nb++;
              }
            
            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz + j) = index + i;
                IndCol(nnz + j) = col(j);
                Val(nnz + j) = value(j);
              }
            
            nnz += nb;
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }

  
  //! Conversion from ColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSymComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* real_ptr = A.GetRealPtr();
    int* imag_ptr = A.GetImagPtr();
    int* real_ind = A.GetRealInd();
    int* imag_ind = A.GetImagInd();
    Treal* real_data = A.GetRealData();
    Treal* imag_data = A.GetImagData();
    
    if (sym)
      {
	nnz *= 2;
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect row; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_col = nb_r + nb_i;
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }
            
            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                row(nb) = real_ind[j];
                value(nb) = T(real_data[j], 0);
                nb++;
              }
            
            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                row(nb) = imag_ind[j];
                value(nb) = T(0, imag_data[j]);
                nb++;
              }
            
            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = row(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;
                
                if (row(j) != i)
                  {
                    IndRow(nnz) = row(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(row(j))++;
                    nnz++;
                  }
              }
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
        
        // Sorting the column numbers...
	Sort(IndCol, IndRow, Val);

	// ...and the row numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndRow, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect row; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int nb_r = real_ptr[i+1] - real_ptr[i];
            int nb_i = imag_ptr[i+1] - imag_ptr[i];
            int size_col = nb_r + nb_i;
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }
            
            int nb = 0;
            for (int j = real_ptr[i]; j < real_ptr[i+1]; j++)
              {
                row(nb) = real_ind[j] + index;
                value(nb) = T(real_data[j], 0);
                nb++;
              }
            
            for (int j = imag_ptr[i]; j < imag_ptr[i+1]; j++)
              {
                row(nb) = imag_ind[j] + index;
                value(nb) = T(0, imag_data[j]);
                nb++;
              }
            
            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndCol(nnz + j) = index + i;
                IndRow(nnz + j) = row(j);
                Val(nnz + j) = value(j);
              }
            
            nnz += nb;
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }

  
  //! Conversion from ArrayRowComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int m = A.GetM();
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    IVect col; Vector<T> value;
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
        if (size_row > col.GetM())
          {
            col.Reallocate(size_row);
            value.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            col(nb) = A.IndexReal(i, j) + index;
            value(nb) = T(A.ValueReal(i, j), 0);
            nb++;
          }

        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            col(nb) = A.IndexImag(i, j) + index;
            value(nb) = T(0, A.ValueImag(i, j));
            nb++;
          }
        
        Assemble(nb, col, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = index + i;
            IndCol(nnz + j) = col(j);
            Val(nnz + j) = value(j);
          }
        
        nnz += nb;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayColComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColComplexSparse,
			       Allocator1>& A,
			       Vector<Tint, VectFull, Allocator2>& IndRow,
			       Vector<Tint, VectFull, Allocator3>& IndCol,
			       Vector<T, VectFull, Allocator4>& Val,
			       int index, bool sym)
  {
    int nnz = A.GetRealDataSize() + A.GetImagDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    nnz = 0;
    IVect row; Vector<T> value;
    for (int i = 0; i < A.GetN(); i++)
      {
        int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
        if (size_col > row.GetM())
          {
            row.Reallocate(size_col);
            value.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = 0; j < A.GetRealColumnSize(i); j++)
          {
            row(nb) = A.IndexReal(i, j) + index;
            value(nb) = T(A.ValueReal(i, j), 0);
            nb++;
          }

        for (int j = 0; j < A.GetImagColumnSize(i); j++)
          {
            row(nb) = A.IndexImag(i, j) + index;
            value(nb) = T(0, A.ValueImag(i, j));
            nb++;
          }
        
        Assemble(nb, row, value);
        for (int j = 0; j < nb; j++)
          {
            IndRow(nnz + j) = row(j);
            IndCol(nnz + j) = index + i;
            Val(nnz + j) = value(j);
          }
        
        nnz += nb;
      }
    
    IndRow.Resize(nnz);
    IndCol.Resize(nnz);
    Val.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayRowSymComplexSparse,
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
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect col; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }
            
            int nb = 0;
            for (int j = 0; j < A.GetRealRowSize(i); j++)
              {
                col(nb) = A.IndexReal(i, j);
                value(nb) = T(A.ValueReal(i, j), 0);
                nb++;
              }
            
            for (int j = 0; j < A.GetImagRowSize(i); j++)
              {
                col(nb) = A.IndexImag(i, j);
                value(nb) = T(0, A.ValueImag(i, j));
                nb++;
              }
            
            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = col(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;
                
                if (col(j) != i)
                  {
                    IndRow(nnz) = col(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(col(j))++;
                    nnz++;
                  }
              }
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
        
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
        nnz = 0;
        IVect col; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
            if (size_row > col.GetM())
              {
                col.Reallocate(size_row);
                value.Reallocate(size_row);
              }
            
            int nb = 0;
            for (int j = 0; j < A.GetRealRowSize(i); j++)
              {
                col(nb) = A.IndexReal(i, j) + index;
                value(nb) = T(A.ValueReal(i, j), 0);
                nb++;
              }
            
            for (int j = 0; j < A.GetImagRowSize(i); j++)
              {
                col(nb) = A.IndexImag(i, j) + index;
                value(nb) = T(0, A.ValueImag(i, j));
                nb++;
              }
            
            Assemble(nb, col, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz + j) = index + i;
                IndCol(nnz + j) = col(j);
                Val(nnz + j) = value(j);
              }
            
            nnz += nb;
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }
  
  
  //! Conversion from ArrayColSymComplexSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2,
	   class Tint, class Allocator3, class Allocator4>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ArrayColSymComplexSparse,
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
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	Vector<int> Ptr(m);
	Ptr.Zero();
        nnz = 0;
        IVect row; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }
            
            int nb = 0;
            for (int j = 0; j < A.GetRealColumnSize(i); j++)
              {
                row(nb) = A.IndexReal(i, j);
                value(nb) = T(A.ValueReal(i, j), 0);
                nb++;
              }
            
            for (int j = 0; j < A.GetImagColumnSize(i); j++)
              {
                row(nb) = A.IndexImag(i, j);
                value(nb) = T(0, A.ValueImag(i, j));
                nb++;
              }
            
            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndRow(nnz) = i + index;
                IndCol(nnz) = row(j) + index;
                Val(nnz) = value(j);
                Ptr(i)++;
                nnz++;
                
                if (row(j) != i)
                  {
                    IndRow(nnz) = row(j) + index;
                    IndCol(nnz) = i + index;
                    Val(nnz) = value(j);
                    Ptr(row(j))++;
                    nnz++;
                  }
              }
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
        
        // Sorting the column numbers...
	Sort(IndCol, IndRow, Val);

	// ...and the row numbers.
	int offset = 0;
	for (int i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndRow, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
        nnz = 0;
        IVect row; Vector<T> value;
	for (int i = 0; i < m; i++)
	  {
            int size_col = A.GetRealColumnSize(i) + A.GetImagColumnSize(i);
            if (size_col > row.GetM())
              {
                row.Reallocate(size_col);
                value.Reallocate(size_col);
              }
            
            int nb = 0;
            for (int j = 0; j < A.GetRealColumnSize(i); j++)
              {
                row(nb) = A.IndexReal(i, j) + index;
                value(nb) = T(A.ValueReal(i, j), 0);
                nb++;
              }
            
            for (int j = 0; j < A.GetImagColumnSize(i); j++)
              {
                row(nb) = A.IndexImag(i, j) + index;
                value(nb) = T(0, A.ValueImag(i, j));
                nb++;
              }
            
            Assemble(nb, row, value);
            for (int j = 0; j < nb; j++)
              {
                IndCol(nnz + j) = index + i;
                IndRow(nnz + j) = row(j);
                Val(nnz + j) = value(j);
              }
            
            nnz += nb;
          }
        
        IndRow.Resize(nnz);
        IndCol.Resize(nnz);
        Val.Resize(nnz);
      }
  }

  
  //! Conversion from coordinate format to RowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(m+1), PtrImag(m+1), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndRow(i)+1)++;
        
	if (imag(Val(i)) != zero)
          PtrImag(IndRow(i)+1)++;
      }

    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    int real_nz = PtrReal(m), imag_nz = PtrImag(m);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (real(Val(offset + j)) != zero)
            {
              IndReal(nb) = IndCol(offset + j);
              ValReal(nb) = real(Val(offset + j));
              nb++;
            }
        
        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (imag(Val(offset + j)) != zero)
            {
              IndImag(nb) = IndCol(offset + j);
              ValImag(nb) = imag(Val(offset + j));
              nb++;
            }
        
        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
        
        offset += Ptr(i);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }
  
  
  //! Conversion from coordinate format to ColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);
    
    // Number of elements per column.
    Vector<int> PtrReal(n+1), PtrImag(n+1), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndCol(i)+1)++;
        
	if (imag(Val(i)) != zero)
          PtrImag(IndCol(i)+1)++;
      }

    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    int real_nz = PtrReal(n), imag_nz = PtrImag(n);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (real(Val(offset + j)) != zero)
            {
              IndReal(nb) = IndRow(offset + j);
              ValReal(nb) = real(Val(offset + j));
              nb++;
            }
        
        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (imag(Val(offset + j)) != zero)
            {
              IndImag(nb) = IndRow(offset + j);
              ValImag(nb) = imag(Val(offset + j));
              nb++;
            }
        
        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
        
        offset += Ptr(i);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to RowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, RowSymComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(m+1), PtrImag(m+1), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndRow(i)+1)++;
            
            if (imag(Val(i)) != zero)
              PtrImag(IndRow(i)+1)++;
          }
      }

    for (int i = 0; i < m; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    int real_nz = PtrReal(m), imag_nz = PtrImag(m);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (i <= IndCol(offset+j))
            if (real(Val(offset + j)) != zero)
              {
                IndReal(nb) = IndCol(offset + j);
                ValReal(nb) = real(Val(offset + j));
                nb++;
              }
        
        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (i <= IndCol(offset+j))
            if (imag(Val(offset + j)) != zero)
              {
                IndImag(nb) = IndCol(offset + j);
                ValImag(nb) = imag(Val(offset + j));
                nb++;
            }
        
        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
        
        offset += Ptr(i);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to ColSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ColSymComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());
    
    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);
    
    // Number of elements per column.
    Vector<int> PtrReal(n+1), PtrImag(n+1), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndCol.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndCol(i)+1)++;
            
            if (imag(Val(i)) != zero)
              PtrImag(IndCol(i)+1)++;
          }
      }

    for (int i = 0; i < n; i++)
      {
        PtrReal(i+1) += PtrReal(i);
        PtrImag(i+1) += PtrImag(i);
      }
    
    int real_nz = PtrReal(n), imag_nz = PtrImag(n);
    
    // Fills matrix 'A'.
    Vector<int> IndReal(real_nz), IndImag(imag_nz);
    Vector<Treal, VectFull, Allocator4> ValReal(real_nz), ValImag(imag_nz);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        int nb = PtrReal(i);
        for (int j = 0; j < Ptr(i); j++)
          if (IndRow(offset+j) <= i)
            if (real(Val(offset + j)) != zero)
              {
                IndReal(nb) = IndRow(offset + j);
                ValReal(nb) = real(Val(offset + j));
                nb++;
              }
        
        nb = PtrImag(i);
        for (int j = 0; j < Ptr(i); j++)
          if (IndRow(offset+j) <= i)
            if (imag(Val(offset + j)) != zero)
              {
                IndImag(nb) = IndRow(offset + j);
                ValImag(nb) = imag(Val(offset + j));
                nb++;
              }
        
        // sorting column numbers
        Sort(PtrReal(i), PtrReal(i+1)-1, IndReal, ValReal);
        Sort(PtrImag(i), PtrImag(i+1)-1, IndImag, ValImag);
        
        offset += Ptr(i);
      }
    
    // providing pointers to A
    A.SetData(m, n, ValReal, PtrReal, IndReal, ValImag, PtrImag, IndImag);
  }

  
  //! Conversion from coordinate format to ArrayRowComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndRow(i))++;
        
	if (imag(Val(i)) != zero)
          PtrImag(IndRow(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealRow(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  A.IndexReal(i, nb) = IndCol(offset + j);
                  A.ValueReal(i, nb) = real(Val(offset + j));
                  nb++;
                }
          }
        
        if (PtrImag(i) > 0)
          {
            A.ReallocateImagRow(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  A.IndexImag(i, nb) = IndCol(offset + j);
                  A.ValueImag(i, nb) = imag(Val(offset + j));
                  nb++;
                }
          }
        
        offset += Ptr(i);
      }

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayColComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColComplexSparse,
				 Allocator4>& A,
				 int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(n), PtrImag(n), Ptr(n);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
	if (real(Val(i)) != zero)
          PtrReal(IndCol(i))++;
        
	if (imag(Val(i)) != zero)
          PtrImag(IndCol(i))++;
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < n; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealColumn(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  A.IndexReal(i, nb) = IndRow(offset + j);
                  A.ValueReal(i, nb) = real(Val(offset + j));
                  nb++;
                }
          }
        
        if (PtrImag(i) > 0)
          {
            A.ReallocateImagColumn(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  A.IndexImag(i, nb) = IndRow(offset + j);
                  A.ValueImag(i, nb) = imag(Val(offset + j));
                  nb++;
                }
          }
        
        offset += Ptr(i);
      }

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayRowSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndRow(i))++;
        if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndRow(i))++;
            
            if (imag(Val(i)) != zero)
              PtrImag(IndRow(i))++;
          }
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealRow(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  if (i <= IndCol(offset+j))
                    {
                      A.IndexReal(i, nb) = IndCol(offset + j);
                      A.ValueReal(i, nb) = real(Val(offset + j));
                      nb++;
                    }
                }
          }
        
        if (PtrImag(i) > 0)
          {
            A.ReallocateImagRow(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  if (i <= IndCol(offset+j))
                    {
                      A.IndexImag(i, nb) = IndCol(offset + j);
                      A.ValueImag(i, nb) = imag(Val(offset + j));
                      nb++;
                    }
                }
          }
        
        offset += Ptr(i);
      }

    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }

  
  //! Conversion from coordinate format to ArrayRowSymComplexSparse.
  template<class T, class Prop, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4>
  void
  ConvertMatrix_from_Coordinates(Vector<int, VectFull, Allocator1>& IndRow,
				 Vector<int, VectFull, Allocator2>& IndCol,
				 Vector<T, VectFull, Allocator3>& Val,
				 Matrix<T, Prop, ArrayColSymComplexSparse,
				 Allocator4>& A, int index)
  {
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    m = max(m, A.GetM());
    n = max(n, A.GetN());

    // Sorts the array 'IndRow'.
    Sort(IndCol, IndRow, Val);
    
    // Number of elements per row.
    Vector<int> PtrReal(m), PtrImag(m), Ptr(m);
    PtrReal.Zero(); PtrImag.Zero(); Ptr.Zero();
    for (int i = 0; i < IndRow.GetM(); i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
        Ptr(IndCol(i))++;
        if (IndRow(i) <= IndCol(i))
          {
            if (real(Val(i)) != zero)
              PtrReal(IndCol(i))++;
            
            if (imag(Val(i)) != zero)
              PtrImag(IndCol(i))++;
          }
      }

    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (int i = 0; i < m; i++)
      {
        if (PtrReal(i) > 0)
          {
            A.ReallocateRealColumn(i, PtrReal(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (real(Val(offset + j)) != zero)
                {
                  if (IndRow(offset+j) <= i)
                    {
                      A.IndexReal(i, nb) = IndRow(offset + j);
                      A.ValueReal(i, nb) = real(Val(offset + j));
                      nb++;
                    }
                }
          }
        
        if (PtrImag(i) > 0)
          {
            A.ReallocateImagColumn(i, PtrImag(i));
            int nb = 0;
            for (int j = 0; j < Ptr(i); j++)
              if (imag(Val(offset + j)) != zero)
                {
                  if (IndRow(offset+j) <= i)
                    {
                      A.IndexImag(i, nb) = IndRow(offset + j);
                      A.ValueImag(i, nb) = imag(Val(offset + j));
                      nb++;
                    }
                }
          }
        
        offset += Ptr(i);
      }

    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }

  
  //! Conversion from RowComplexSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
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

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    int* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < m; i++)
      {
        int size_real = real_ptr_[i+1] - real_ptr_[i];
        int size_imag = imag_ptr_[i+1] - imag_ptr_[i];
        int size_row = size_real + size_imag;
        if (size_row > val.GetM())
          {
            col.Reallocate(size_row);
            val.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
            val(nb) = T(real_data_[j], zero);
            col(nb) = real_ind_[j];
            nb++;
          }
        
        for (int j = imag_ptr_[i]; j < imag_ptr_[i+1]; j++)
          {
            val(nb) = T(zero, imag_data_[j]);
            col(nb) = imag_ind_[j];
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ColComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from RowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
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

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    int* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < m; i++)
      {
        int size_real = real_ptr_[i+1] - real_ptr_[i];
        int size_imag = imag_ptr_[i+1] - imag_ptr_[i];
        int size_row = size_real + size_imag;
        if (size_row > val.GetM())
          {
            col.Reallocate(size_row);
            val.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
            val(nb) = T(real_data_[j], zero);
            col(nb) = real_ind_[j];
            nb++;
          }
        
        for (int j = imag_ptr_[i]; j < imag_ptr_[i+1]; j++)
          {
            val(nb) = T(zero, imag_data_[j]);
            col(nb) = imag_ind_[j];
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from RowSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
  }

  
  //! Conversion from ColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
  }


  //! Conversion from ColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);    
  }

  
  //! Conversion from ArrayRowComplexSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
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
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
        if (size_row > val.GetM())
          {
            col.Reallocate(size_row);
            val.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            val(nb) = T(A.ValueReal(i, j), zero);
            col(nb) = A.IndexReal(i, j);
            nb++;
          }
        
        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            val(nb) = T(zero, A.ValueImag(i, j));
            col(nb) = A.IndexImag(i, j);
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ArrayColComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint, VectFull, Alloc3> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from ArrayRowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
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
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(m+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < m; i++)
      {
        int size_row = A.GetRealRowSize(i) + A.GetImagRowSize(i);
        if (size_row > val.GetM())
          {
            col.Reallocate(size_row);
            val.Reallocate(size_row);
          }
        
        int nb = 0;
        for (int j = 0; j < A.GetRealRowSize(i); j++)
          {
            val(nb) = T(A.ValueReal(i, j), zero);
            col(nb) = A.IndexReal(i, j);
            nb++;
          }
        
        for (int j = 0; j < A.GetImagRowSize(i); j++)
          {
            val(nb) = T(zero, A.ValueImag(i, j));
            col(nb) = A.IndexImag(i, j);
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from ArrayRowSymComplexSparse to CSR
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
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
  
    
  //! Conversion from ArrayColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from ArrayColSymComplexSparse to CSR format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSR(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Val)
  {
    // Matrix (m,n) with nnz entries.
    int m = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndRow;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to row numbers.
    Sort(IndRow, IndCol, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(m + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per row.
    for (int i = 0; i < IndRow.GetM(); i++)
      Ptr(IndRow(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from RowComplexSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      " implemented for RowComplexSparse storage");
    
    // Matrix (m,n) with nnz entries.
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
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }
  
  
  //! Conversion from ColComplexSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      " implemented for ColComplexSparse storage");
    
    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }

    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    int* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < n; i++)
      {
        int size_real = real_ptr_[i+1] - real_ptr_[i];
        int size_imag = imag_ptr_[i+1] - imag_ptr_[i];
        int size_col = size_real + size_imag;
        if (size_col > val.GetM())
          {
            col.Reallocate(size_col);
            val.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
            val(nb) = T(real_data_[j], zero);
            col(nb) = real_ind_[j];
            nb++;
          }
        
        for (int j = imag_ptr_[i]; j < imag_ptr_[i+1]; j++)
          {
            val(nb) = T(zero, imag_data_[j]);
            col(nb) = imag_ind_[j];
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }
    
    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from RowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
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
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }
  

  //! Conversion from RowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, RowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int n = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per column.
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from ColSymComplexSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    int* real_ptr_ = A.GetRealPtr();
    int* real_ind_ = A.GetRealInd();
    Treal* real_data_ = A.GetRealData();

    int* imag_ptr_ = A.GetImagPtr();
    int* imag_ind_ = A.GetImagInd();
    Treal* imag_data_ = A.GetImagData();
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < n; i++)
      {
        int size_real = real_ptr_[i+1] - real_ptr_[i];
        int size_imag = imag_ptr_[i+1] - imag_ptr_[i];
        int size_col = size_real + size_imag;
        if (size_col > val.GetM())
          {
            col.Reallocate(size_col);
            val.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = real_ptr_[i]; j < real_ptr_[i+1]; j++)
          {
            val(nb) = T(real_data_[j], zero);
            col(nb) = real_ind_[j];
            nb++;
          }
        
        for (int j = imag_ptr_[i]; j < imag_ptr_[i+1]; j++)
          {
            val(nb) = T(zero, imag_data_[j]);
            col(nb) = imag_ind_[j];
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }

    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ColSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int n = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per column.
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from ArrayRowComplexSparse to CSC format
  /*!
    if sym_pat is equal to true, the pattern is symmetrized
    by adding artificial null entries
   */
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      "implemented for RowComplexSparse storage");
    
    // Matrix (m,n) with nnz entries.
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
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }
  
  
  //! Conversion from ArrayColComplexSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    if (sym_pat)
      throw Undefined("ConvertToCSC", "Symmetrization of pattern not"
                      "implemented for ColComplexSparse storage");
    
    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < n; i++)
      {
        int size_real = A.GetRealColumnSize(i);
        int size_imag = A.GetImagColumnSize(i);
        int size_col = size_real + size_imag;
        if (size_col > val.GetM())
          {
            col.Reallocate(size_col);
            val.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = 0; j < size_real; j++)
          {
            val(nb) = T(A.ValueReal(i, j), zero);
            col(nb) = A.IndexReal(i, j);
            nb++;
          }
        
        for (int j = 0; j < size_imag; j++)
          {
            val(nb) = T(zero, A.ValueImag(i, j));
            col(nb) = A.IndexImag(i, j);
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }

    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }

  
  //! Conversion from ArrayRowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
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
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }
  

  //! Conversion from ArrayRowSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayRowSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int n = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per column.
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }

  
  //! Conversion from ArrayColSymComplexSparse to CSC
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    Symmetric& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndCol,
                    Vector<T, VectFull, Alloc4>& Value,
                    bool sym_pat)
  {
    int n = A.GetN();
    if (n <= 0)
      {
	Ptr.Clear();
	IndCol.Clear();
	Value.Clear();
	return;
      }
    
    typedef typename ClassComplexType<T>::Treal Treal;
    Treal zero(0);
    
    int real_nnz = A.GetRealDataSize();
    int imag_nnz = A.GetImagDataSize();
    
    Ptr.Reallocate(n+1);
    IndCol.Reallocate(real_nnz + imag_nnz);
    Value.Reallocate(real_nnz + imag_nnz);
    int nnz = 0;
    Ptr(0) = 0;
    IVect col;
    Vector<T> val;
    for (int i = 0; i < n; i++)
      {
        int size_real = A.GetRealColumnSize(i);
        int size_imag = A.GetImagColumnSize(i);
        int size_col = size_real + size_imag;
        if (size_col > val.GetM())
          {
            col.Reallocate(size_col);
            val.Reallocate(size_col);
          }
        
        int nb = 0;
        for (int j = 0; j < size_real; j++)
          {
            val(nb) = T(A.ValueReal(i, j), zero);
            col(nb) = A.IndexReal(i, j);
            nb++;
          }
        
        for (int j = 0; j < size_imag; j++)
          {
            val(nb) = T(zero, A.ValueImag(i, j));
            col(nb) = A.IndexImag(i, j);
            nb++;
          }
        
        Assemble(nb, col, val);
        
        for (int j = 0; j < nb; j++)
          {
            IndCol(Ptr(i) + j) = col(j);
            Value(Ptr(i) + j) = val(j);
          }
        
        nnz += nb;
        Ptr(i+1) = nnz;
      }

    IndCol.Resize(nnz);
    Value.Resize(nnz);
  }


  //! Conversion from ArrayColSymComplexSparse to CSC format
  template<class T, class Prop, class Alloc1,
           class Tint, class Alloc2, class Alloc3, class Alloc4>
  void ConvertToCSC(const Matrix<T, Prop, ArrayColSymComplexSparse, Alloc1>& A,
                    General& sym, Vector<Tint, VectFull, Alloc2>& Ptr,
                    Vector<Tint, VectFull, Alloc3>& IndRow,
                    Vector<T, VectFull, Alloc4>& Val, bool sym_pat)
  {
    // Matrix (m,n) with nnz entries.
    int n = A.GetN();
    
    // Conversion in coordinate format.
    Vector<Tint> IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val, 0, true);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    Ptr.Reallocate(n + 1);
    Ptr.Fill(0);
    
    // Counting non-zero entries per column.
    for (int i = 0; i < IndCol.GetM(); i++)
      Ptr(IndCol(i) + 1)++;

    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (int i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayColComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from RowComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColSymComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColSymComplexSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ColComplexSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    Vector<int> Ptr, IndCol;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSR(mat_array, prop, Ptr, IndCol, Value);
    
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    mat_csr.SetData(m, n, Value, Ptr, IndCol);
  }


  //! Conversion from ArrayRowComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayRowSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayRowSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from ArrayColSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayColSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from RowComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }


  //! Conversion from RowSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  
  

  //! Conversion from RowSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  


  //! Conversion from ColComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }
  
  
  //! Conversion from ColSymComplexSparse to ColSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    General prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  
  

  //! Conversion from ColSymComplexSparse to ColSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSymSparse, Allocator1>& B)
  {
    Vector<int> Ptr, IndRow;
    Vector<T1, VectFull, Allocator1> Value;

    Symmetric prop;
    ConvertToCSC(A, prop, Ptr, IndRow, Value);
    
    int m = A.GetM();
    int n = A.GetN();
    B.SetData(m, n, Value, Ptr, IndRow);
  }  

  
  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowComplexSparse, Allocator>& A,
	    Matrix<T, Prop, RowComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, RowSymComplexSparse, Allocator>& A,
	    Matrix<T, Prop, RowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ColComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ColSymComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ColSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayRowSymComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayRowComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayColComplexSparse, Allocator>& B)
  {
    B = A;
  }


  //! B = A.
  template<class T, class Prop, class Allocator>
  void CopyMatrix(const Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& A,
	    Matrix<T, Prop, ArrayColSymComplexSparse, Allocator>& B)
  {
    B = A;
  }

  
  //! Conversion from ColSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Vector<int> Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    General sym;
    ConvertToCSR(A, sym, Ptr, Ind, Value);
    
    // counting number of non-zero elements
    typedef typename ClassComplexType<T1>::Treal Treal;
    int imag_nnz = 0, real_nnz = 0;
    Treal zero(0);
    for (int i = 0; i < Value.GetM(); i++)
      {
        if (real(Value(i)) != zero)
          real_nnz++;
        
        if (imag(Value(i)) != zero)
          imag_nnz++;
      }
    
    Vector<int> Ptr_real(m+1), Ptr_imag(m+1),
      Ind_real(real_nnz), Ind_imag(imag_nnz);
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < m; i++)
      {
        for (int j = Ptr(i); j < Ptr(i+1); j++)
          {
            if (real(Value(j)) != zero)
              {
                Ind_real(real_nnz) = Ind(j);
                Val_real(real_nnz) = real(Value(j));
                real_nnz++;
              }
            
            if (imag(Value(j)) != zero)
              {
                Ind_imag(imag_nnz) = Ind(j);
                Val_imag(imag_nnz) = imag(Value(j));
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! conversion from ColSparse to ColComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    int nnz = A.GetDataSize();
    
    // counting number of non-zero elements
    int imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (int i = 0; i < nnz; i++)
      {
        if (real(data_[i]) != zero)
          real_nnz++;
        
        if (imag(data_[i]) != zero)
          imag_nnz++;
      }
    
    Vector<int> Ptr_real(n+1), Ptr_imag(n+1),
      Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < n; i++)
      {
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                Ind_real(real_nnz) = ind_[j];
                Val_real(real_nnz) = real(data_[j]);
                real_nnz++;
              }
            
            if (imag(data_[j]) != zero)
              {
                Ind_imag(imag_nnz) = ind_[j];
                Val_imag(imag_nnz) = imag(data_[j]);
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! Conversion from RowSymSparse to ColSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
       Matrix<T1, Prop1, ColSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Vector<int> Ptr, Ind;
    Vector<T0, VectFull, Allocator0> Value;
    
    Symmetric sym;
    ConvertToCSC(A, sym, Ptr, Ind, Value);
    
    // counting number of non-zero elements
    int imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (int i = 0; i < Value.GetM(); i++)
      {
        if (real(Value(i)) != zero)
          real_nnz++;
        
        if (imag(Value(i)) != zero)
          imag_nnz++;
      }
    
    Vector<int> Ptr_real(m+1), Ptr_imag(m+1),
      Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < m; i++)
      {
        for (int j = Ptr(i); j < Ptr(i+1); j++)
          {
            if (real(Value(j)) != zero)
              {
                Ind_real(real_nnz) = Ind(j);
                Val_real(real_nnz) = real(Value(j));
                real_nnz++;
              }
            
            if (imag(Value(j)) != zero)
              {
                Ind_imag(imag_nnz) = Ind(j);
                Val_imag(imag_nnz) = imag(Value(j));
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }
  
  
  //! conversion from RowSymSparse to RowSymComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    int nnz = A.GetDataSize();
    
    // counting number of non-zero elements
    int imag_nnz = 0, real_nnz = 0;
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    for (int i = 0; i < nnz; i++)
      {
        if (real(data_[i]) != zero)
          real_nnz++;
        
        if (imag(data_[i]) != zero)
          imag_nnz++;
      }
    
    Vector<int> Ptr_real(n+1), Ptr_imag(n+1),
      Ind_real(real_nnz), Ind_imag(imag_nnz);
    
    Vector<Treal, VectFull, Allocator1> Val_real(real_nnz), Val_imag(imag_nnz);
    
    // filling arrays Ind_real, Ind_imag, Val_real, Val_imag
    Ptr_real(0) = 0; Ptr_imag(0) = 0;
    real_nnz = 0; imag_nnz = 0;
    for (int i = 0; i < n; i++)
      {
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                Ind_real(real_nnz) = ind_[j];
                Val_real(real_nnz) = real(data_[j]);
                real_nnz++;
              }
            
            if (imag(data_[j]) != zero)
              {
                Ind_imag(imag_nnz) = ind_[j];
                Val_imag(imag_nnz) = imag(data_[j]);
                imag_nnz++;
              }
          }
        
        Ptr_real(i+1) = real_nnz;
        Ptr_imag(i+1) = imag_nnz;
      }

    
    // creating the matrix
    B.SetData(m, n, Val_real, Ptr_real, Ind_real,
              Val_imag, Ptr_imag, Ind_imag);
  }

  
  //! Conversion from ColSparse to ArrayRowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Matrix<T0, Prop0, ArrayRowSparse, Allocator0> Ar;
    CopyMatrix(A, Ar);

    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = 0; j < Ar.GetRowSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              size_real++;
            
            if (imag(Ar.Value(i, j)) != zero)
              size_imag++;
          }
        
        B.ReallocateRealRow(i, size_real);
        B.ReallocateImagRow(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = 0; j < Ar.GetRowSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              {
                B.IndexReal(i, size_real) = Ar.Index(i, j);
                B.ValueReal(i, size_real) = real(Ar.Value(i, j));
                size_real++;
              }
            
            if (imag(Ar.Value(i, j)) != zero)
              {
                B.IndexImag(i, size_imag) = Ar.Index(i, j);
                B.ValueImag(i, size_imag) = imag(Ar.Value(i, j));
                size_imag++;
              }
          }
        
        Ar.ClearRow(i);
      }
  }
  
  
  //! conversion from ColSparse to ArrayColComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ColSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayColComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              size_real++;
            
            if (imag(data_[j]) != zero)
              size_imag++;
          }
        
        B.ReallocateRealColumn(i, size_real);
        B.ReallocateImagColumn(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                B.IndexReal(i, size_real) = ind_[j];
                B.ValueReal(i, size_real) = real(data_[j]);
                size_real++;
              }
            
            if (imag(data_[j]) != zero)
              {
                B.IndexImag(i, size_imag) = ind_[j];
                B.ValueImag(i, size_imag) = imag(data_[j]);
                size_imag++;
              }
          }
      }
  }
  
  
  //! Conversion from RowSymSparse to ArrayColSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayColSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM(), n = A.GetN();
    Matrix<T0, Prop0, ArrayColSymSparse, Allocator0> Ar;
    CopyMatrix(A, Ar);

    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);    
    B.Reallocate(m, n);
    for (int i = 0; i < m; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = 0; j < Ar.GetColumnSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              size_real++;
            
            if (imag(Ar.Value(i, j)) != zero)
              size_imag++;
          }
        
        B.ReallocateRealColumn(i, size_real);
        B.ReallocateImagColumn(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = 0; j < Ar.GetColumnSize(i); j++)
          {
            if (real(Ar.Value(i, j)) != zero)
              {
                B.IndexReal(i, size_real) = Ar.Index(i, j);
                B.ValueReal(i, size_real) = real(Ar.Value(i, j));
                size_real++;
              }
            
            if (imag(Ar.Value(i, j)) != zero)
              {
                B.IndexImag(i, size_imag) = Ar.Index(i, j);
                B.ValueImag(i, size_imag) = imag(Ar.Value(i, j));
                size_imag++;
              }
          }
        
        Ar.ClearColumn(i);
      }
  }
  
  
  //! conversion from RowSymSparse to ArrayRowSymComplexSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowSymComplexSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T0* data_ = A.GetData();
    
    typedef typename ClassComplexType<T1>::Treal Treal;
    Treal zero(0);
    B.Reallocate(m, n);
    for (int i = 0; i < n; i++)
      {
        int size_real = 0, size_imag = 0;
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              size_real++;
            
            if (imag(data_[j]) != zero)
              size_imag++;
          }
        
        B.ReallocateRealRow(i, size_real);
        B.ReallocateImagRow(i, size_imag);
        size_real = 0;
        size_imag = 0;
        for (int j = ptr_[i]; j < ptr_[i+1]; j++)
          {
            if (real(data_[j]) != zero)
              {
                B.IndexReal(i, size_real) = ind_[j];
                B.ValueReal(i, size_real) = real(data_[j]);
                size_real++;
              }
            
            if (imag(data_[j]) != zero)
              {
                B.IndexImag(i, size_imag) = ind_[j];
                B.ValueImag(i, size_imag) = imag(data_[j]);
                size_imag++;
              }
          }       
      }
  }

  
  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    typedef typename ClassComplexType<T1>::Treal Treal;
    Vector<Treal, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int> IndRow_real(m + 1);
    Vector<int> IndRow_imag(m + 1);
    Vector<int> IndCol_real(nnz_real);
    Vector<int> IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }


  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();

    // Allocation of arrays for CSR.
    typedef typename ClassComplexType<T1>::Treal Treal;
    Vector<Treal, VectFull, Allocator1> Val_real(nnz_real), Val_imag(nnz_imag);
    Vector<int> IndRow_real(m + 1);
    Vector<int> IndRow_imag(m + 1);
    Vector<int> IndCol_real(nnz_real);
    Vector<int> IndCol_imag(nnz_imag);

    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }

	IndRow_real(i + 1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }

	IndRow_imag(i + 1) = ind_imag;
      }

    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }
  
  
  //! Conversion from ArrayRowSymComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr)
  {
    Vector<int> IndRow, IndCol;
    Vector<T1, VectFull, Allocator1> Value;
    
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = mat_array.GetM();
    int n = mat_array.GetN();    
    mat_csr.Reallocate(m, n);
    int k = 0;
    for (int i = 0; i < m; i++)
      {
        int k0 = k;
        while ((k < IndRow.GetM()) && (IndRow(k) <= i))
          k++;

        int size_row = k - k0;
        if (size_row > 0)
          {
            mat_csr.ReallocateRow(i, size_row);
            for (int j = 0; j < size_row; j++)
              {
                mat_csr.Index(i, j) = IndCol(k0 + j);
                mat_csr.Value(i, j) = Value(k0 + j);
              }
          }
        else
          mat_csr.ClearRow(i);
      }
  }



  //! Conversion from ArrayRowSymComplexSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int size_real = A.GetRealRowSize(i);
        int size_imag = A.GetImagRowSize(i);
        int jr = 0, ji = 0;
        int size_row = 0;
        while (jr < size_real)
          {
            int col = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < col))
              {
                col_num(size_row) = A.IndexImag(i, ji);
                value(size_row) = T0(0, A.ValueImag(i, ji));
                ji++; size_row++;
              }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr),
				     A.ValueImag(i, ji));
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr), 0);
                size_row++;
              }
            
            jr++;
          }

        while (ji < size_imag)
          {
            col_num(size_row) = A.IndexImag(i, ji);
            value(size_row) = T0(0, A.ValueImag(i, ji));
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from ArrayRowComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int m = A.GetM();
    int n = A.GetN();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int size_real = A.GetRealRowSize(i);
        int size_imag = A.GetImagRowSize(i);
        int jr = 0, ji = 0;
        int size_row = 0;
        while (jr < size_real)
          {
            int col = A.IndexReal(i, jr);
            while ((ji < size_imag) && (A.IndexImag(i, ji) < col))
              {
                col_num(size_row) = A.IndexImag(i, ji);
                value(size_row) = T0(0, A.ValueImag(i, ji));
                ji++; size_row++;
              }
            
            if ((ji < size_imag) && (A.IndexImag(i, ji) == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr),
				     A.ValueImag(i, ji));
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(A.ValueReal(i, jr), 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < size_imag)
          {
            col_num(size_row) = A.IndexImag(i, ji);
            value(size_row) = T0(0, A.ValueImag(i, ji));
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from RowSymComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& mat_csr)
  {
    Vector<int> IndRow, IndCol;
    Vector<T1, VectFull, Allocator1> Value;
    
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Value, 0, true);
    
    // sorting by rows
    Sort(IndRow, IndCol, Value);

    int m = mat_array.GetM();
    int n = mat_array.GetN();    
    mat_csr.Reallocate(m, n);
    int k = 0;
    for (int i = 0; i < m; i++)
      {
        int k0 = k;
        while ((k < IndRow.GetM()) && (IndRow(k) <= i))
          k++;

        int size_row = k - k0;
        if (size_row > 0)
          {
            mat_csr.ReallocateRow(i, size_row);
            for (int j = 0; j < size_row; j++)
              {
                mat_csr.Index(i, j) = IndCol(k0 + j);
                mat_csr.Value(i, j) = Value(k0 + j);
              }
          }
        else
          mat_csr.ClearRow(i);
      }
  }


  //! Conversion from RowComplexSparse to ArrayRowSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {  
    typedef typename ClassComplexType<T0>::Treal Treal;  
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int jr = ptr_real[i], ji = ptr_imag[i];
        int size_row = 0;
        while (jr < ptr_real[i+1])
          {
            int col = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < col))
              {
                col_num(size_row) = ind_imag[ji];
                value(size_row) = T0(0, data_imag[ji]);
                ji++; size_row++;
              }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr],
				     data_imag[ji]);
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr], 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            col_num(size_row) = ind_imag[ji];
            value(size_row) = T0(0, data_imag[ji]);
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }


  //! Conversion from RowSymComplexSparse to ArrayRowSymSparse
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  CopyMatrix(const Matrix<T0, Prop0, RowSymComplexSparse, Allocator0>& A,
       Matrix<T1, Prop1, ArrayRowSymSparse, Allocator1>& B)
  { 
    typedef typename ClassComplexType<T0>::Treal Treal;   
    int m = A.GetM();
    int n = A.GetN();
    int* ptr_real = A.GetRealPtr();
    int* ind_real = A.GetRealInd();
    Treal* data_real = A.GetRealData();
    int* ptr_imag = A.GetImagPtr();
    int* ind_imag = A.GetImagInd();
    Treal* data_imag = A.GetImagData();
    B.Reallocate(m, n);
    Vector<T0> value(n);
    Vector<int> col_num(n);
    for (int i = 0; i < m; i++)
      {
        int jr = ptr_real[i], ji = ptr_imag[i];
        int size_row = 0;
        while (jr < ptr_real[i+1])
          {
            int col = ind_real[jr];
            while ((ji < ptr_imag[i+1]) && (ind_imag[ji] < col))
              {
                col_num(size_row) = ind_imag[ji];
                value(size_row) = T0(0, data_imag[ji]);
                ji++; size_row++;
              }
            
            if ((ji < ptr_imag[i+1]) && (ind_imag[ji] == col))
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr],
				     data_imag[ji]);
                ji++; size_row++;
              }
            else
              {
                col_num(size_row) = col;
                value(size_row) = T0(data_real[jr], 0);
                size_row++;
              }
            
            jr++;
          }
        
        while (ji < ptr_imag[i+1])
          {
            col_num(size_row) = ind_imag[ji];
            value(size_row) = T0(0, data_imag[ji]);
            ji++; size_row++;
          }
        
        B.ReallocateRow(i, size_row);
        for (int j = 0; j < size_row; j++)
          {
            B.Index(i, j) = col_num(j);
            B.Value(i, j) = value(j);
          }
      }
  }

}


#define SELDON_FILE_MATRIX_COMPLEX_CONVERSIONS_CXX
#endif
