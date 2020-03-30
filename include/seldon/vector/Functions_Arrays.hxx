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


#ifndef SELDON_FILE_FUNCTIONS_ARRAYS_HXX

namespace Seldon
{

  template<class T, class Storage, class Allocator>
  int PartitionQuickSort(int m, int n,
			 Vector<T, Storage, Allocator>& t);

  template<class T, class Storage, class Allocator>
  void QuickSort(int m, int n,
		 Vector<T, Storage, Allocator>& t);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  int PartitionQuickSort(int m, int n,
			 Vector<T1, Storage1, Allocator1>& t1,
			 Vector<T2, Storage2, Allocator2>& t2,
			 Vector<T3, Storage3, Allocator3>& t3);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void QuickSort(int m, int n,
		 Vector<T1, Storage1, Allocator1>& t1,
		 Vector<T2, Storage2, Allocator2>& t2,
		 Vector<T3, Storage3, Allocator3>& t3);
  
  template<class T, class Storage, class Allocator>
  void MergeSort(int m, int n, Vector<T, Storage, Allocator>& tab1);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2);
  
  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void MergeSort(int m, int n, Vector<T1, Storage1, Allocator1>& tab1,
		 Vector<T2, Storage2, Allocator2>& tab2,
		 Vector<T3, Storage3, Allocator3>& tab3);
  
  template<class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2 >
  void Assemble(int& n, Vector<int, Storage1, Allocator1>& Node,
		Vector<T2, Storage2, Allocator2>& Vect);
  
  template<class T, class Storage1, class Allocator1>
  void Assemble(int& n, Vector<T, Storage1, Allocator1>& Node);
  
  template<class T, class Storage1, class Allocator1>
  void Assemble(Vector<T, Storage1, Allocator1>& Node);

  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2);

  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(int& n, Vector<T, Storage1, Allocator1>& Node);

  template<class T, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node,
		       Vector<T2, Storage2, Allocator2>& Node2);

  template<class T, class Storage1, class Allocator1>
  void RemoveDuplicate(Vector<T, Storage1, Allocator1>& Node);

  template<class T, class Storage, class Allocator>
  void Sort(int m, int n, Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int m, int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);

  template<class T, class Storage, class Allocator>
  void Sort(int n, Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(int n, Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);
  
  template<class T, class Storage, class Allocator>
  void Sort(Vector<T, Storage, Allocator>& V);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2);

  template<class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2,
	   class T3, class Storage3, class Allocator3>
  void Sort(Vector<T1, Storage1, Allocator1>& V,
	    Vector<T2, Storage2, Allocator2>& V2,
	    Vector<T3, Storage3, Allocator3>& V3);

  template<class T, class Storage, class Allocator>
  bool HasElement(Vector<T, Storage, Allocator>& X, T& a);
  
} // namespace Seldon

#define SELDON_FILE_FUNCTIONS_ARRAYS_HXX
#endif
