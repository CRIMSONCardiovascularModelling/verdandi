// Copyright (C) 2010 Vivien Mallet
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


#include "SeldonHeader.hxx"


#ifndef SELDON_WITH_COMPILED_LIBRARY
#include "share/StorageInline.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#endif


#include <complex>
typedef std::complex<float> complexfloat;
typedef std::complex<double> complexdouble;


namespace Seldon
{

#ifdef SELDON_WITH_SLOW_COMPILED_LIBRARY

  SELDON_EXTERN template class Matrix_Base<@scalar, MallocAlloc<@scalar> >;
  SELDON_EXTERN template class Matrix_Pointers<@scalar, General, @storage_blasGE, MallocAlloc<@scalar> >;
  SELDON_EXTERN template class Matrix<@scalar, General, @storage_blasGE, MallocAlloc<@scalar> >;
  SELDON_EXTERN template void Matrix_Pointers<@real_complex, General, @storage_blasGE, MallocAlloc<@real_complex> >::Fill(const @real_complex&);
  SELDON_EXTERN template void Matrix_Pointers<@scalar, General, @storage_blasGE, MallocAlloc<@scalar> >::Fill(const int&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<@scalar, General, @storage_blasGE, MallocAlloc<@scalar> >&);
#endif


  SELDON_EXTERN template class Matrix_Base<@scalar, NewAlloc<@scalar> >;
  SELDON_EXTERN template class Matrix_Pointers<@scalar, General, @storage_blasGE, NewAlloc<@scalar> >;
  SELDON_EXTERN template class Matrix<@scalar, General, @storage_blasGE, NewAlloc<@scalar> >;
  SELDON_EXTERN template void Matrix_Pointers<@real_complex, General, @storage_blasGE, NewAlloc<@real_complex> >::Fill(const @real_complex&);
  SELDON_EXTERN template void Matrix_Pointers<@scalar, General, @storage_blasGE, NewAlloc<@scalar> >::Fill(const int&);

#ifndef SWIG
  SELDON_EXTERN template ostream& operator << (ostream&, const Matrix<@scalar, General, @storage_blasGE, NewAlloc<@scalar> >&);
#endif

#endif

}
