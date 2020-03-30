// Copyright (C) 2001-2011 Vivien Mallet
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


#ifndef SELDON_FILE_MATRIX_HERMPACKED_INLINE_CXX

#include "Matrix_HermPacked.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>::Matrix_HermPacked():
    Matrix_Base<T, Allocator>()
  {
  }


  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>::~Matrix_HermPacked()
  {
    Clear();
  }


  /*******************
   * BASIC FUNCTIONS *
   *******************/


  //! Returns the number of elements stored in memory.
  /*!
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline int Matrix_HermPacked<T, Prop, Storage, Allocator>::GetDataSize() const
  {
    return (this->m_ * (this->m_ + 1)) / 2;
  }


  //! Returns size of A in bytes used to store the matrix.
  template <class T, class Prop, class Storage, class Allocator>
  inline int64_t Matrix_HermPacked<T, Prop, Storage, Allocator>::GetMemorySize() const
  {
    int64_t taille = sizeof(*this) + int64_t(GetDataSize())*sizeof(T);
    return taille;
  }
  

  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/
  
  
  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline const typename 
  Matrix_HermPacked<T, Prop, Storage, Allocator>::value_type
  Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, j, this->m_, this->n_, "Matrix_HermPacked");
#endif

    if (i > j)
      return conj(this->data_[Storage::GetFirst(j * this->m_
						- (j*(j+1)) / 2 + i,
						(i*(i+1)) / 2 + j)]);
    else
      return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					   (j*(j+1)) / 2 + i)];
  }


  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Val(int i, int j)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_HermPacked");
#endif

    return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					 (j*(j+1)) / 2 + i)];
  }


  //! Direct access method.
  /*!
    This method allows access to elements stored in memory, i.e. elements
    from the upper part. i <= j must be satisfied.
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Val(int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBoundsSym(i, j, this->m_, this->n_, "Matrix_HermPacked");
#endif

    return this->data_[Storage::GetFirst(i * this->n_ - (i*(i+1)) / 2 + j,
					 (j*(j+1)) / 2 + i)];
  }


  //! Returns access to an element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Get(int i, int j)
  {
    return this->Val(i, j);
  }


  //! Returns access to an element (i, j)
  /*!
    \param i row index.
    \param j column index.
    \return The value of the matrix at (i, j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::Get(int i, int j) const
  {
    return this->Val(i, j);
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>::reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator[] (int i)
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_HermPacked");
#endif

    return this->data_[i];
  }


  //! Access to elements of the data array.
  /*!
    Provides a direct access to the data array.
    \param i index.
    \return i-th element of the data array.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline typename Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::const_reference
  Matrix_HermPacked<T, Prop, Storage, Allocator>::operator[] (int i) const
  {

#ifdef SELDON_CHECK_BOUNDS
    CheckBounds(i, this->GetDataSize(), "Matrix_HermPacked");
#endif

    return this->data_[i];
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_HermPacked<T, Prop, Storage, Allocator>&
  Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_HermPacked<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }


  //! Sets an element of the matrix.
  /*!
    \param i row index.
    \param j column index.
    \param x new value for the matrix element (\a i, \a j).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Set(int i, int j, const T& x)
  {
    if (i > j)
      this->Val(j, i) = conj(x);
    else
      this->Val(i, j) = x;
  }


  //! Duplicates a matrix.
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::Copy(const Matrix_HermPacked<T, Prop, Storage, Allocator>& A)
  {
    this->Reallocate(A.GetM(), A.GetN());

    Allocator::memorycpy(this->data_, A.GetData(), this->GetDataSize());
  }
  

#ifdef SELDON_WITH_VIRTUAL
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::ApplySor(Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::ApplySor(const class_SeldonTrans& trans, Vector<T>& x, const Vector<T>& r,
	     const typename ClassComplexType<T>::Treal& omega,
	     int nb_iter, int stage_ssor) const
  {
    SOR(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	x, r, omega, nb_iter, stage_ssor);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::MltAddVector(const T& alpha, const Vector<T>& x,
		 const T& beta, Vector<T>& y) const
  {
    MltAdd(alpha,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::MltAddVector(const T& alpha, const SeldonTranspose& trans,
		 const Vector<T>& x,
		 const T& beta, Vector<T>& y) const
  {
    MltAdd(alpha, trans,
	   static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this),
	   x, beta, y);
  }
  
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>
  ::MltVector(const Vector<T>& x, Vector<T>& y) const
  {
    Mlt(static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_HermPacked<T, Prop, Storage, Allocator>  
  ::MltVector(const SeldonTranspose& trans,
	      const Vector<T>& x, Vector<T>& y) const
  {
    Mlt(trans,
	static_cast<const Matrix<T, Prop, Storage, Allocator>& >(*this), x, y);
  }

  template <class T, class Prop, class Storage, class Allocator>
  inline bool Matrix_HermPacked<T, Prop, Storage, Allocator>  
  ::IsSymmetric() const
  {
    return false;
  }
#endif

  
  ///////////////////////////
  // MATRIX<COLHERMPACKED> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHermPacked, Allocator>::Matrix():
    Matrix_HermPacked<T, Prop, ColHermPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j column-major hermitian matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHermPacked, Allocator>::Matrix(int i, int j):
    Matrix_HermPacked<T, Prop, ColHermPacked, Allocator>(i, j)
  {
  }


  /*******************
   * OTHER FUNCTIONS *
   *******************/


  //! Fills the matrix with a given value.
  /*!
    \param x the value with which the matrix is to be filled.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, ColHermPacked, Allocator>&
  Matrix<T, Prop, ColHermPacked, Allocator>
  ::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, ColHermPacked, Allocator>&
  Matrix<T, Prop, ColHermPacked, Allocator>::operator= (const Matrix<T, Prop,
                                                       ColHermPacked,
                                                       Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }

  
  
  ///////////////////////////
  // MATRIX<ROWHERMPACKED> //
  ///////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    On exit, the matrix is an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHermPacked, Allocator>::Matrix():
    Matrix_HermPacked<T, Prop, RowHermPacked, Allocator>()
  {
  }


  //! Main constructor.
  /*! Builds a i x j row-major hermitian matrix in packed form.
    \param i number of rows.
    \param j number of columns.
    \note 'j' is assumed to be equal to 'i' and is therefore discarded.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHermPacked, Allocator>::Matrix(int i, int j):
    Matrix_HermPacked<T, Prop, RowHermPacked, Allocator>(i, j)
  {
  }


  /*******************
   * OTHER FUNCTIONS *
   *******************/


  //! Fills the matrix with a given value.
  /*!
    \param x the value with which the matrix is to be filled.
  */
  template <class T, class Prop, class Allocator>
  template <class T0>
  inline Matrix<T, Prop, RowHermPacked, Allocator>&
  Matrix<T, Prop, RowHermPacked, Allocator>
  ::operator= (const T0& x)
  {
    this->Fill(x);

    return *this;
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: \a A is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Allocator>
  inline Matrix<T, Prop, RowHermPacked, Allocator>&
  Matrix<T, Prop, RowHermPacked, Allocator>::operator= (const Matrix<T, Prop,
                                                       RowHermPacked,
                                                       Allocator>& A)
  {
    this->Copy(A);
    return *this;
  }

  
} // namespace Seldon.

#define SELDON_FILE_MATRIX_HERMPACKED_INLINE_CXX
#endif
