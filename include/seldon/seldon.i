// Copyright (C) 2001-2009 Vivien Mallet
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


%module seldon
%{
#include "SeldonHeader.hxx"
namespace Seldon
{
  void skip_vector_double(istream& input_stream);
  void skip_matrix_double(istream& input_stream);
}
  %}

%include "std_iostream.i"
%include "std_string.i"

using namespace std;

namespace std
{
  class ifstream: public istream
  {
  public:
    ifstream(const char *fname);
    ~ifstream();
    bool is_open();
    void close();
  };

  class ofstream: public ostream
  {
  public:
    ofstream(const char *fname);
    ~ofstream();
    bool is_open();
    void close();
  };

}

%include "share/Errors.hxx"
%exception
{
  try
    {
      $action
	}
  catch(Seldon::WrongIndex& e)
    {
      PyErr_SetString(PyExc_IndexError, e.What().c_str());
      return NULL;
    }
  catch(Seldon::WrongRow& e)
    {
      PyErr_SetString(PyExc_IndexError, e.What().c_str());
      return NULL;
    }
  catch(Seldon::WrongCol& e)
    {
      PyErr_SetString(PyExc_IndexError, e.What().c_str());
      return NULL;
    }
  catch(Seldon::WrongDim& e)
    {
      PyErr_SetString(PyExc_IndexError, e.What().c_str());
      return NULL;
    }
  catch(Seldon::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "SeldonHeader.hxx"
%include "share/Common.hxx"
%include "share/Storage.hxx"
%include "share/Properties.hxx"
%include "array/Array3D.hxx"
%include "array/Array4D.hxx"
%include "vector/Vector.hxx"
%include "vector/VectorCollection.hxx"
%include "vector/SparseVector.hxx"
%include "matrix/Matrix_Base.hxx"
%include "matrix/Matrix_Pointers.hxx"
%include "matrix_sparse/Matrix_Sparse.hxx"
%include "share/Allocator.hxx"

namespace Seldon
{
  %extend Vector<int, VectFull, MallocAlloc<int> >
  {
    int __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, int value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }
  %extend Vector<float, VectFull, MallocAlloc<float> >
  {
    float __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, float value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }
  %extend Vector<double, VectFull, MallocAlloc<double> >
  {
    double __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, double value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }

  %extend Vector<Vector<int, VectFull, MallocAlloc<int> >, Collection, MallocAlloc<Vector<int, VectFull, MallocAlloc<int> > > >
  {
    int __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, int value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }
  %extend Vector<Vector<double, VectFull, MallocAlloc<double> >, Collection, MallocAlloc<Vector<double, VectFull, MallocAlloc<double> > > >
  {
    double __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, double value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }

  %extend Matrix<int, General, RowMajor, MallocAlloc<int> >
  {
    int __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > __getitem__(int i)
    {
      Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, int value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }
  %extend Matrix<float, General, RowMajor, MallocAlloc<float> >
  {
    float __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > __getitem__(int i)
    {
      Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, float value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }
  %extend Matrix<double, General, RowMajor, MallocAlloc<double> >
  {
    double __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > __getitem__(int i)
    {
      Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, double value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }

  %extend Array3D<int, MallocAlloc<int> >
  {
    int __getitem__(PyObject* args)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k);
    }
    Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > __getitem__(int i, int j)
    {
      Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > v(self->GetLength3());
      for (int k = 0; k < self->GetLength3(); k++)
	v(k) = (*self)(i, j, k);
      return v;
    }
    Seldon::Matrix<int, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<int> >
      __getitem__(int i)
    {
      Seldon::Matrix<int, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<int> >
        m(self->GetLength2(), self->GetLength3());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          m(j, k) = (*self)(i, j, k);
      return m;
    }
    void __setitem__(PyObject* args, int value)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }
  %extend Array3D<float, MallocAlloc<float> >
  {
    float __getitem__(PyObject* args)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k);
    }
    Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > __getitem__(int i, int j)
    {
      Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > v(self->GetLength3());
      for (int k = 0; k < self->GetLength3(); k++)
	v(k) = (*self)(i, j, k);
      return v;
    }
    Seldon::Matrix<float, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<float> >
      __getitem__(int i)
    {
      Seldon::Matrix<float, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<float> >
        m(self->GetLength2(), self->GetLength3());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          m(j, k) = (*self)(i, j, k);
      return m;
    }
    void __setitem__(PyObject* args, float value)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }
  %extend Array3D<double, MallocAlloc<double> >
  {
    double __getitem__(PyObject* args)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k);
    }
    Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > __getitem__(int i, int j)
    {
      Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > v(self->GetLength3());
      for (int k = 0; k < self->GetLength3(); k++)
	v(k) = (*self)(i, j, k);
      return v;
    }
    Seldon::Matrix<double, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<double> >
      __getitem__(int i)
    {
      Seldon::Matrix<double, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<double> >
        m(self->GetLength2(), self->GetLength3());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          m(j, k) = (*self)(i, j, k);
      return m;
    }
    void __setitem__(PyObject* args, double value)
    {
      int i, j, k;
      int success = PyArg_ParseTuple(args, "iii", &i, &j, &k);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }

  %extend Array4D<int, MallocAlloc<int> >
  {
    int __getitem__(PyObject* args)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k, l);
    }
    Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > __getitem__(int i, int j, int k)
    {
      Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > v(self->GetLength4());
      for (int l = 0; l < self->GetLength4(); l++)
	v(l) = (*self)(i, j, k, l);
      return v;
    }
    Seldon::Matrix<int, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<int> >
      __getitem__(int i, int j)
    {
      Seldon::Matrix<int, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<int> >
        m(self->GetLength3(), self->GetLength4());
      for (int k = 0; k < self->GetLength3(); k++)
        for (int l = 0; l < self->GetLength4(); l++)
          m(k, l) = (*self)(i, j, k, l);
      return m;
    }
    Seldon::Array3D<int, Seldon::MallocAlloc<int> > __getitem__(int i)
    {
      Seldon::Array3D<int, Seldon::MallocAlloc<int> >
        a(self->GetLength2(), self->GetLength3(), self->GetLength4());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          for (int l = 0; l < self->GetLength4(); l++)
            a(j, k, l) = (*self)(i, j, k, l);
      return a;
    }
    void __setitem__(PyObject* args, int value)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k, l) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }
  %extend Array4D<float, MallocAlloc<float> >
  {
    float __getitem__(PyObject* args)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k, l);
    }
    Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > __getitem__(int i, int j, int k)
    {
      Seldon::Vector<float, Seldon::VectFull, Seldon::MallocAlloc<float> > v(self->GetLength4());
      for (int l = 0; l < self->GetLength4(); l++)
	v(l) = (*self)(i, j, k, l);
      return v;
    }
    Seldon::Matrix<float, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<float> >
      __getitem__(int i, int j)
    {
      Seldon::Matrix<float, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<float> >
        m(self->GetLength3(), self->GetLength4());
      for (int k = 0; k < self->GetLength3(); k++)
        for (int l = 0; l < self->GetLength4(); l++)
          m(k, l) = (*self)(i, j, k, l);
      return m;
    }
    Seldon::Array3D<float, Seldon::MallocAlloc<float> > __getitem__(int i)
    {
      Seldon::Array3D<float, Seldon::MallocAlloc<float> >
        a(self->GetLength2(), self->GetLength3(), self->GetLength4());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          for (int l = 0; l < self->GetLength4(); l++)
            a(j, k, l) = (*self)(i, j, k, l);
      return a;
    }
    void __setitem__(PyObject* args, float value)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k, l) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }
  %extend Array4D<double, MallocAlloc<double> >
  {
    double __getitem__(PyObject* args)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j, k, l);
    }
    Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > __getitem__(int i, int j, int k)
    {
      Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > v(self->GetLength4());
      for (int l = 0; l < self->GetLength4(); l++)
	v(l) = (*self)(i, j, k, l);
      return v;
    }
    Seldon::Matrix<double, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<double> >
      __getitem__(int i, int j)
    {
      Seldon::Matrix<double, Seldon::General, Seldon::RowMajor, Seldon::MallocAlloc<double> >
        m(self->GetLength3(), self->GetLength4());
      for (int k = 0; k < self->GetLength3(); k++)
        for (int l = 0; l < self->GetLength4(); l++)
          m(k, l) = (*self)(i, j, k, l);
      return m;
    }
    Seldon::Array3D<double, Seldon::MallocAlloc<double> > __getitem__(int i)
    {
      Seldon::Array3D<double, Seldon::MallocAlloc<double> >
        a(self->GetLength2(), self->GetLength3(), self->GetLength4());
      for (int j = 0; j < self->GetLength2(); j++)
        for (int k = 0; k < self->GetLength3(); k++)
          for (int l = 0; l < self->GetLength4(); l++)
            a(j, k, l) = (*self)(i, j, k, l);
      return a;
    }
    void __setitem__(PyObject* args, double value)
    {
      int i, j, k, l;
      int success = PyArg_ParseTuple(args, "iiii", &i, &j, &k, &l);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j, k, l) = value;
    }
    unsigned long __len__() {
      return self->GetLength1();
    }
  }

  %template(IntMalloc) MallocAlloc<int>;
  
  %template(BaseSeldonVectorInt) Vector_Base<int, MallocAlloc<int> >;
  %template(VectorInt) Vector<int, VectFull, MallocAlloc<int> >;
  %template(FloatMalloc) MallocAlloc<float>;
  %template(BaseSeldonVectorFloat) Vector_Base<float, MallocAlloc<float> >;
  %template(VectorFloat) Vector<float, VectFull, MallocAlloc<float> >;
  %template(DoubleMalloc) MallocAlloc<double>;
  %template(BaseSeldonVectorDouble) Vector_Base<double, MallocAlloc<double> >;
  %template(VectorDouble) Vector<double, VectFull, MallocAlloc<double> >;

  %template(MatrixBaseInt) Matrix_Base<int, MallocAlloc<int> >;
  %template(MatrixPointersInt) Matrix_Pointers<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixInt) Matrix<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixBaseFloat) Matrix_Base<float, MallocAlloc<float> >;
  %template(MatrixPointersFloat) Matrix_Pointers<float, General, RowMajor, MallocAlloc<float> >;
  %template(MatrixFloat) Matrix<float, General, RowMajor, MallocAlloc<float> >;
  %template(MatrixBaseDouble) Matrix_Base<double, MallocAlloc<double> >;
  %template(MatrixPointersDouble) Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
  %template(MatrixDouble) Matrix<double, General, RowMajor, MallocAlloc<double> >;

  %template(Array3DInt) Array3D<int, MallocAlloc<int> >;
  %template(Array3DFloat) Array3D<float, MallocAlloc<float> >;
  %template(Array3DDouble) Array3D<double, MallocAlloc<double> >;
  %template(Array4DInt) Array4D<int, MallocAlloc<int> >;
  %template(Array4DFloat) Array4D<float, MallocAlloc<float> >;
  %template(Array4DDouble) Array4D<double, MallocAlloc<double> >;

  %template(VectorSparseDouble) Vector<double, VectSparse, MallocAlloc<double> >;

  %template(BaseMatrixSparseDouble) Matrix_Sparse<double, General, RowSparse, MallocAlloc<double> >;
  %template(MatrixSparseDouble) Matrix<double, General, RowSparse, MallocAlloc<double> >;

  void skip_vector_double(istream& input_stream);
  void skip_matrix_double(istream& input_stream);
}


// For conversions from Seldon to Numpy, and from Numpy to Seldon.
%pythoncode %{
def load_vector(filename, array = False):
    """
    Loads a Seldon vector (in double precision) from a file. If 'array' is set
    to True, the vector is converted to a numpy array.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    vector = seldon.VectorDouble()
    vector.Read(stream, True)

    if array:
        import numpy
        vector = numpy.array(vector)

    if isinstance(filename, str):
        stream.close()

    return vector


def load_vector_list(filename, array = False, begin = 0, N = 0):
    """
    Loads a list of Seldon vectors (in double precision) from a file, skipping
    the first 'begin' vectors. If 'array' is set to True, the vectors are
    converted to numpy arrays. 'N' is the number of vectors to read; all
    vectors are read if 'N' is 0 or negative.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    begin = max(begin, 0)
    count = 0
    while stream.peek() != -1 and count != begin:
        skip_vector_double(stream)
        count += 1

    vector_list = []
    count = 0
    while stream.peek() != -1 and (N <= 0 or count < N):
        vector = seldon.VectorDouble()
        vector.Read(stream, True)
        vector_list.append(vector)
        count += 1

    if array:
        import numpy
        vector_list = [numpy.array(x) for x in vector_list]

    if isinstance(filename, str):
        stream.close()

    return vector_list


def to_vector(v):
    """
    Converts a list or a numpy array to a Seldon vector (in double precision).
    """
    import seldon
    out = seldon.VectorDouble(len(v))
    for i in range(len(v)):
        out[i] = v[i]
    return out


def load_matrix(filename, array = False):
    """
    Loads a Seldon matrix (in double precision) from a file. If 'array' is set
    to True, the matrix is converted to a numpy array.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    matrix = seldon.MatrixDouble()
    matrix.Read(stream, True)

    if array:
        import numpy
        matrix = numpy.array(matrix)

    if isinstance(filename, str):
        stream.close()

    return matrix


def load_matrix_list(filename, array = False, begin = 0, N = 0):
    """
    Loads a list of Seldon matrices (in double precision) from a file,
    skipping the first 'begin' matrices.. If 'array' is set to True, the
    matrices are converted to numpy arrays. 'N' is the number of matrices to
    read; all matrices are read if 'N' is 0 or negative.
    """
    import seldon
    if isinstance(filename, str):
        stream = seldon.ifstream(filename)
    else: # assuming 'filename' is already a stream.
        stream = filename

    begin = max(begin, 0)
    count = 0
    while stream.peek() != -1 and count != begin:
        skip_matrix_double(stream)
        count += 1

    matrix_list = []
    count = 0
    while stream.peek() != -1 and (N <= 0 or count < N):
        matrix = seldon.MatrixDouble()
        matrix.Read(stream, True)
        matrix_list.append(matrix)
        count += 1

    if array:
        import numpy
        matrix_list = [numpy.array(x) for x in matrix_list]

    if isinstance(filename, str):
        stream.close()

    return matrix_list


def to_matrix(m):
    """
    Converts a numpy array to a Seldon matrix (in double precision).
    """
    import seldon
    out = seldon.MatrixDouble(m.shape[0], m.shape[1])
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            out[i, j] = m[i, j]
    return out
%}
