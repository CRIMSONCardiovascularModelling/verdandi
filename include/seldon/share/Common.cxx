// Copyright (C) 2001-2012 Vivien Mallet
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


#ifndef SELDON_FILE_COMMON_CXX

#include "Common.hxx"

namespace Seldon
{
  //! returns extension of a string 
  string GetExtension(const string& nom)
  {
    size_t index = nom.find_last_of('.');
    if (index == string::npos)
      return string("");
    
    string extension = nom.substr(index+1, nom.size()-index);  
    return extension;
  }

  
  //! returns base of a string
  string GetBaseString(const string& nom)
  {
    size_t index = nom.find_last_of('.');
    if (index == string::npos)
      return nom;
    
    string base = nom.substr(0, index);
    return base;
  }
  
    
#ifdef SELDON_WITH_HDF5
  //! Gives for most C types the corresponding HDF5 memory datatype.
  /*!
    \param[in] input variable to analyze.
    \return HDF5 memory type of \a input.
  */
  template <class T>
  hid_t GetH5Type(T& input)
  {
    double d;
    float f;
    int i;
    long l;
    char c;
    unsigned char uc;
    long long ll;
    unsigned int ui;
    unsigned short us;
    unsigned long ul;
    unsigned long long ull;

    if (typeid(input) == typeid(d))
      return H5T_NATIVE_DOUBLE;
    if (typeid(input) == typeid(f))
      return H5T_NATIVE_FLOAT;
    if (typeid(input) == typeid(i))
      return H5T_NATIVE_INT;
    if (typeid(input) == typeid(l))
      return H5T_NATIVE_LONG;
    if (typeid(input) == typeid(c))
      return H5T_NATIVE_CHAR;
    if (typeid(input) == typeid(uc))
      return H5T_NATIVE_UCHAR;
    if (typeid(input) == typeid(ll))
      return H5T_NATIVE_LLONG;
    if (typeid(input) == typeid(ui))
      return H5T_NATIVE_UINT;
    if (typeid(input) == typeid(us))
      return H5T_NATIVE_USHORT;
    if (typeid(input) == typeid(ul))
      return H5T_NATIVE_ULONG;
    if (typeid(input) == typeid(ull))
      return H5T_NATIVE_ULLONG;
    else
      throw Error("hid_t GetH5Type(T& input)",
                  "Type has no corresponding native HDF5 datatype.");
  }
#endif


}  // namespace Seldon.

#define SELDON_FILE_COMMON_CXX
#endif
