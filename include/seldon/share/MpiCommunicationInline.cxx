// Copyright (C) 2014 INRIA
// Author(s): Marc Duruflé
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

#ifndef SELDON_FILE_MPI_COMMUNICATION_INLINE_CXX

namespace Seldon
{

  inline const MPI::Datatype& GetMpiDataType(const Vector<bool>&)
  {
    return MPI::BOOL;
  }

  inline const MPI::Datatype& GetMpiDataType(const Vector<int>&)
  {
    return MPI::INTEGER;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<float>&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<complex<float> >&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<double>&)
  {
    return MPI::DOUBLE;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const Vector<complex<double> >&)
  {
    return MPI::DOUBLE;
  }

  inline const MPI::Datatype& GetMpiDataType(const int&)
  {
    return MPI::INTEGER;
  }

  inline const MPI::Datatype& GetMpiDataType(const bool&)
  {
    return MPI::BOOL;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const float&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const complex<float>&)
  {
    return MPI::FLOAT;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const double&)
  {
    return MPI::DOUBLE;
  }
  
  inline const MPI::Datatype& GetMpiDataType(const complex<double>&)
  {
    return MPI::DOUBLE;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(const Vector<T>&)
  {
    return 1;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(const Vector<complex<T> >&)
  {
    return 2;
  }

  template<class T>
  inline int GetRatioMpiDataType(const T&)
  {
    return 1;
  }
  
  template<class T>
  inline int GetRatioMpiDataType(complex<T>&)
  {
    return 2;
  }
  
}

#define SELDON_FILE_MPI_COMMUNICATION_INLINE_CXX
#endif
