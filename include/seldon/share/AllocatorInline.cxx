// Copyright (C) 2001-2010 Vivien Mallet, INRIA
// Author(s): Vivien Mallet, Marc Fragu
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


#ifndef SELDON_FILE_ALLOCATOR_INLINE_CXX

#include "malloc.h"
#include "Allocator.hxx"
#include <stdexcept>

namespace Seldon
{


  /////////////////
  // MALLOCALLOC //
  /////////////////


  template <class T>
  inline typename MallocAlloc<T>::pointer
  MallocAlloc<T>::allocate(int num, void* h)
  {
    return static_cast<pointer>( malloc(num * sizeof(T)) );
  }

  template <class T>
  inline void MallocAlloc<T>::deallocate(pointer data, int num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void* MallocAlloc<T>::reallocate(pointer data, int num, void* h)
  {
    return realloc(reinterpret_cast<void*>(data), num * sizeof(T));
  }

  template <class T>
  inline void MallocAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  MallocAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


  /////////////////
  // CALLOCALLOC //
  /////////////////


  template <class T>
  inline typename CallocAlloc<T>::pointer
  CallocAlloc<T>::allocate(int num, void* h)
  {
    return static_cast<pointer>( calloc(num, sizeof(T)) );
  }

  template <class T>
  inline void CallocAlloc<T>::deallocate(pointer data, int num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void* CallocAlloc<T>::reallocate(pointer data, int num, void* h)
  {
    void* data_pointer_cast_void = reinterpret_cast<void*>(data);

    // For CallocAlloc, the code requires that the allocated memory is initialised.
    //
    // This requires some work, as realloc doesn't know how to initialise
    // any new memory it allocates (which is required here if the new array length
    // is longer than the old).
    //
    // The current fix for this won't be portable (possibly GNU Linux only), but 
    // hopefully we'll always be compiling on GNU Linux...
    //
    // Anyway, it's preferable to the previous status quo of using
    // uninitialised memory. 
    //
    // There's hundreds of calls to reallocate() throught this library;
    // the correct fix would be to go through them and initialise the
    // memory properly in each case, but that will take a long time.
    // So unfortunately, this sub-optimal, probably-unstable fix
    // is all I can do for now.
    //
    // To quote the man page:
    //
    // "The value returned by malloc_usable_size() may be greater than the
    // requested size of the allocation because of alignment and minimum size
    // constraints. Although the excess bytes can be overwritten by the
    // application without ill effects, this is not good programming practice:
    // the number of excess bytes in an allocation depends on the underlying
    // implementation."
    //
    // So this usage is OK (on GNU), but not a great solution.
    // Unfortunately, I haven't found a better one beyond redesigning all
    // calls to reallocate() in the library.

    #ifndef __GNUC__
    #error malloc_usable_size() may not be available outside GNU C. Find an alternative for your compiler.
    #endif
    const size_t old_array_size = malloc_usable_size(data_pointer_cast_void);

    const size_t new_array_size = num * sizeof(T);
    void* reallocated_array = realloc(reinterpret_cast<void*>(data_pointer_cast_void), new_array_size);

    if (!reallocated_array) {
      throw std::runtime_error("Failed to reallocate Seldon array in AllocatiorInline.cxx.");
    }

    // Pointer arithmetic on void* is not defined by the standard,
    // but it is a GNU extension, see 
    // 
    // https://gcc.gnu.org/onlinedocs/gcc-4.4.2/gcc/Pointer-Arith.html#Pointer-Arith
    //
    // So in GNU C, sizeof(void) = 1 (byte).
    // So can increment the pointer a byte at a time...
    #ifndef __GNUC__
    #error void pointer arithmetic may not be available outside of GNU C. Investigate alternatives for your compiler.
    #endif
    if (new_array_size > old_array_size) {
      void* start_of_newly_allocated_array_region = reallocated_array + old_array_size;

      const size_t length_of_newly_allocated_part_of_array = new_array_size - old_array_size;
      memset(start_of_newly_allocated_array_region, 0, length_of_newly_allocated_part_of_array);
    }

    return reallocated_array;
  }

  template <class T>
  inline void CallocAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  CallocAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


  //////////////
  // NEWALLOC //
  //////////////


  template <class T>
  inline typename NewAlloc<T>::pointer
  NewAlloc<T>::allocate(int num, void* h)
  {
    return static_cast<pointer>(new T[num]);
  }

  template <class T>
  inline void NewAlloc<T>::deallocate(pointer data, int num, void* h)
  {
    delete [] data;
  }

  template <class T>
  inline void* NewAlloc<T>::reallocate(pointer data, int num, void* h)
  {
    if (data != NULL)
      delete [] data;
    return (new T[num]);
  }

  template <class T>
  inline void NewAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  NewAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    for (size_t i = 0; i < num; i++)
      datat[i] = datas[i];
  }


  //////////////
  // NANALLOC //
  //////////////


  template <class T>
  inline void NaNAlloc<T>::deallocate(pointer data, int num, void* h)
  {
    free(data);
  }

  template <class T>
  inline void NaNAlloc<T>::memoryset(pointer data, char c, size_t num)
  {
    memset(reinterpret_cast<void*>(data), c, num);
  }

  template <class T>
  inline void
  NaNAlloc<T>::memorycpy(pointer datat, pointer datas, size_t num)
  {
    memcpy(reinterpret_cast<void*>(datat), reinterpret_cast<void*>(datas),
	   num * sizeof(T));
  }


} // namespace Seldon.

#define SELDON_FILE_ALLOCATOR_INLINE_CXX
#endif
