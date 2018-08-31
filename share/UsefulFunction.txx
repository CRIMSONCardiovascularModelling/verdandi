// Copyright (C) 2008, INRIA
// Author(s): Vivien Mallet, Anne Tilloy, Marc Fragu
//
// This file is part of the data assimilation library Verdandi.
//
// Verdandi is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Verdandi. If not, see http://www.gnu.org/licenses/.
//
// For more information, visit the Verdandi web site:
//      http://verdandi.gforge.inria.fr/


#ifndef VERDANDI_FILE_SHARE_USEFULFUNCTION_TXX
#define VERDANDI_FILE_SHARE_USEFULFUNCTION_TXX


namespace Verdandi
{


    //! Returns a value interpolated from a 2D field.
    /*! The output value is produced by bilinear interpolation of the field \a
      input at (\a x, \a y).
      \param x_min abscissa of the first point in \a input.
      \param Delta_x step along x.
      \param y_min ordinate of the first point in \a input.
      \param Delta_y step along y.
      \param input values of the field: \a input(i, j) is the value at
      \f$(x_{min} + i \times \Delta x, y_{min} + j \times \Delta y)\f$.
      \param x abscissa of the point where the field is interpolated.
      \param y ordinate of the point where the field is interpolated.
      \return The value of the field interpolated at (\a x, \a y).
    */
    template <class T, class TM>
    T interpolate(T x_min, T Delta_x, T y_min, T Delta_y,
                  const Matrix<TM>& input, T x, T y)
    {
        int Nx = input.GetN();
        int Ny = input.GetM();

        T distance_x = (x - x_min) / Delta_x;
        int pos_x = int(distance_x);
        int one_pos_x = 1 + pos_x;
        if (pos_x < 0 || pos_x >= Nx - 1)
            throw ErrorConfiguration("interpolate");
        T weight_x = distance_x - T(pos_x);
        T one_weight_x = 1. - weight_x;

        T distance_y = (y - y_min) / Delta_y;
        int pos_y = int(distance_y);
        int one_pos_y = 1 + pos_y;
        if (pos_y < 0 || pos_y >= Ny - 1)
            throw ErrorConfiguration("interpolate");
        T weight_y = distance_y - T(pos_y);
        T one_weight_y = 1. - weight_y;

        return one_weight_y * one_weight_x * input(pos_y, pos_x)
            + one_weight_y * weight_x * input(pos_y, one_pos_x)
            + weight_y * weight_x * input(one_pos_y, one_pos_x)
            + weight_y * one_weight_x * input(one_pos_y, pos_x);
    }


    /*! \brief Returns the coordinate of a point identified with a global
      index in a multidimensional grid. */
    /*! A global index gives the position in a grid with a single integer. For
      example, in 2D, the global index of the grid point \f$(i, j)\f$ is \f$i
      \times N + j\f$ if there are N points along the second dimension. This
      function would return the coordinates \f$(x, y)\f$ of the point: \f$x =
      x_0 + i \times \Delta_x\f$ and \f$y = y_0 + j \times \Delta_y\f$. It is
      assumed that the grid has a constant space step along each direction.
      \param[in] index the global index.
      \param[in] minimum coordinates of the lower left corner of the grid.
      \param[in] step space step along each dimension.
      \param[in] shape dimensions of the grid.
      \param[out] coordinate coordinates.
    */
    template <class T>
    void get_coordinate(int index, const Vector<T>& minimum,
                        const Vector<T>& step, const Vector<int>& shape,
                        Vector<T>& coordinate)
    {
        int Ndimension = shape.GetLength();
        Vector<int> position(Ndimension);
        get_position(index, shape, position);

        coordinate.Reallocate(Ndimension);
        for (int i = 0; i < Ndimension; i++)
            coordinate(i) = minimum(i) + T(position(i)) * step(i);
    }


    //! Converts strings to most types.
    /*!
      \param[in] s string to be converted.
      \param[out] out \a s converted to 'T'.
    */
    template <class T>
    void convert(const string& s, T& out)
    {
        istringstream str(s);
        str >> out;
    }


    //! Splits a string.
    /*!
      The string is split according to delimiters and elements are stored
      in the vector 'vect'.
      \param[in] str string to be split.
      \param[out] vect (output) vector containing elements of the string.
      \param[in] delimiters (optional) delimiters. Default: " \n\t".
    */
    template <class T>
    void split(string str, vector<T>& vect, string delimiters)
    {
        vect.clear();

        T tmp;
        string::size_type index_beg, index_end;

        index_beg = str.find_first_not_of(delimiters);

        while (index_beg != string::npos)
        {
            index_end = str.find_first_of(delimiters, index_beg);
            convert(str.substr(index_beg, index_end == string::npos ?
                               string::npos : (index_end - index_beg)), tmp);
            vect.push_back(tmp);
            index_beg = str.find_first_not_of(delimiters, index_end);
        }
    }


    //! Checks equality with a tolerance interval of epsilon.
    /*!
      \param x first number.
      \param y second number.
      \param epsilon relative tolerance.
      \return True if \a x equals \a y with a relative tolerance of \a
      epsilon.
    */
    template<class T>
    bool is_equal(T x, T y, T epsilon)
    {
        return abs(x - y) <=  0.5 * epsilon * (abs(x) + abs(y));
    }


    //! Checks whether a number is multiple of another, with given tolerance.
    /*!
      \param x possible multiple.
      \param d base number.
      \param epsilon relative tolerance.
      \return True if \a x is a multiple of \a d with a relative tolerance of
      epsilon.
    */
    template<class T>
    bool is_multiple(T x, T d, T epsilon)
    {
        int i = int(x / d + .5);
        return is_equal(x, T(i) * d, epsilon);
    }


    //! Builds a diagonal sparse matrix.
    /*!
      \param[in] size number of lines (or columns).
      \param[in] diagonal_value value on the diagonal of the matrix.
      \param[out] matrix the initialized diagonal sparse matrix.
    */
    template <class T>
    void build_diagonal_sparse_matrix(int size, T diagonal_value,
                                      Matrix<T, General, RowSparse>& matrix)
    {
        Vector<int> column(size);
        Vector<int> pointer(size + 1);
        Vector<T> value(size);

        value.Fill(diagonal_value);

        for (int i = 0; i < size; i++)
        {
            column(i) = i;
            pointer(i) = i;
        }
        pointer(size) = size;

        matrix.SetData(size, size, value, pointer, column);
    }


    ////////////////////
    // LINEAR ALGEBRA //
    ////////////////////


    //! Computes the Cholesky decomposition.
    /*!
      \param[in] A on entry, a symmetric definite positive matrix. On exit,
      the square root \f$ S \f$ of A as given by the Cholesky decomposition:
      \f$ A = S S^T \f$, with \a S is a lower triangular matrix.
    */
    template <class T, class Allocator>
    void GetCholesky(Matrix<T, General, RowMajor, Allocator>& A)
    {
        Matrix<T, General, RowSymPacked> A_sympacked(A.GetM(), A.GetN());
        for (int i = 0; i < A.GetM(); i++)
            for (int j = i; j < A.GetN(); j++)
                A_sympacked(i, j) = A(i, j);

        GetCholesky(A_sympacked);

        A.Zero();
        for (int i = 0; i < A.GetM(); i++)
            for (int j = 0; j <= i; j++)
                A(i, j) = A_sympacked(i, j);
    }


    //! Computes the Cholesky decomposition.
    /*!
      \param[in] A on entry, a symmetric definite positive matrix. On exit,
      the square root \f$ S \f$ of A as given by the Cholesky decomposition:
      \f$ A = S S^T \f$, with \a S is a lower triangular matrix.
    */
    template <class T, class Allocator>
    void GetCholesky(Matrix<T, General, RowSparse, Allocator>& A)
    {
        int m = A.GetM();
        int n = A.GetN();

        Matrix<T, General, ArrayRowSymSparse, Allocator> A_array_sym(m, n);

        T* data = A.GetData();
        int* ptr = A.GetPtr();
        int* column = A.GetInd();

        for (int i = 0; i < m; i++)
            for (int j = ptr[i]; j < ptr[i + 1]; j++)
                if (column[j] >= i)
                    A_array_sym.AddInteraction(i, column[j], data[j]);

        GetCholesky(A_array_sym);

        Matrix<T, General, ArrayRowSparse> A_array(m, n);
        for (int i = 0; i < m; i++)
            for (int j = 0; j <= i; j++)
                A_array.AddInteraction(i, j, A_array_sym(i, j));

        for (int i = 0; i < n; i++)
            A_array(i, i) = T(1) / A_array(i, i);

        A.Clear();
        Copy(A_array, A);
    }



    //! This function overwrites a sparse matrix with its inverse.
    /*!
      \param[in,out] A the matrix to be inverted.
    */
    template <class T, class Allocator>
    void GetInverse(Matrix<T, General, RowSparse, Allocator>& A)
    {
        Matrix<T, General, RowMajor, Allocator> A_dense;
        Copy(A, A_dense);

        GetInverse(A_dense);

        Matrix<T, General, ArrayRowSparse, Allocator> A_array;
        ConvertDenseToArrayRowSparse(A_dense, A_array);

        Copy(A_array, A);
    }


    //! Conversion from 'RowMajor' to 'RowSymPacked' format.
    /*!
      \param[in] A the 'RowMajor' matrix to be converted.
      \param[out] B the matrix \a A  stored in the 'RowSymPacked' format.
    */
    template <class T0, class Allocator0,
              class T1, class Allocator1>
    void Copy(const Matrix<T0, General, RowMajor, Allocator0>& A,
              Matrix<T1, General, RowSymPacked, Allocator1>& B)
    {
        int m = A.GetM();
        int n = A.GetN();

        B.Reallocate(m, n);

        for (int i = 0; i < A.GetM(); i++)
            for (int j = i; j < A.GetN(); j++)
                B(i, j) = A(i, j);
    }


    //! Conversion from 'RowSparse' to 'RowMajor' format.
    /*!
      \param[in] A the 'RowSparse' matrix to be converted.
      \param[out] A_dense the matrix \a A  stored in the 'RowMajor' format.
    */
    template <class T, class Allocator>
    void Copy(const Matrix<T, General, RowSparse, Allocator>& A,
              Matrix<T, General, RowMajor, Allocator>& A_dense)
    {
        int m = A.GetM();
        int n = A.GetN();

        A_dense.Reallocate(m, n);

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A_dense(i, j) = A(i, j);
    }


    //! Adds a 'vector collection' to a 'full vector'.
    /*!
      \param[in] alpha a given scalar.
      \param[in] X a given collection vector.
      \param[out] Y the result \f$ Y = Y + \alpha X\f$.
    */
    template <class T0,
              class T1, class Allocator1,
              class T2, class Allocator2>
    void Add(const T0 alpha, const Vector<T1, Collection, Allocator1>& X,
             Vector<T2, VectFull, Allocator2>& Y)
    {
        if (alpha != T0(0))
        {
            T0 alpha_ = alpha;
            int ma = X.GetM();
#ifdef VERDANDI_CHECK_DIMENSIONS
            CheckDim(X, Y, "Add(alpha, X, Y)");
#endif
            for (int i = 0; i < ma; i++)
                Y(i) += alpha_ * X(i);
        }
    }


#ifdef VERDANDI_WITH_DEPRECATED_SELDON
    //! Solves a sparse linear system using LU factorization.
    /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
      \f$ X \f$ and \f$ Y \f$ are vectors.
      \param[in] M the sparse matrix of the linear system, to be factorized in
      LU form by UMFPACK, SuperLU or Mumps. On exit, \a M is cleared.
      \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
      solution \f$ X \f$ of the system.
    */
    template <class T, class Prop0, class Allocator0, class Allocator1>
    void GetAndSolveLU(Matrix<T, Prop0, ColSparse, Allocator0>& M,
                       Vector<T, VectFull, Allocator1>& Y)
    {
        Solve(M, Y);
    }


    //! Solves a sparse linear system using LU factorization.
    /*! This function solves \f$ M X = Y \f$ where \f$ M \f$ is a matrix, and
      \f$ X \f$ and \f$ Y \f$ are vectors.
      \param[in] M the sparse matrix of the linear system, to be factorized in
      LU form by UMFPACK, SuperLU or Mumps. On exit, \a M is cleared.
      \param[in,out] Y on entry, the right-hand side \f$ Y \f$; on exit, the
      solution \f$ X \f$ of the system.
    */
    template <class T, class Prop0, class Allocator0, class Allocator1>
    void GetAndSolveLU(Matrix<T, Prop0, RowSparse, Allocator0>& M,
                       Vector<T, VectFull, Allocator1>& Y)
    {
        Solve(M, Y);
    }
#endif


    //! Sets the column of a given matrix.
    /*!
      \param[in] X a given column vector.
      \param[in] i a given index.
      \param[out] M a given matrix.
    */
    template <class T0, class Allocator0, class T1, class Allocator1>
    void SetCol(const Vector<T1, VectFull, Allocator1>& X,
                int i, Matrix<T0, General, RowSparse, Allocator0>& M)
    {
        Vector<T1, VectSparse, Allocator1> X_sparse;
        ConvertDenseToSparse(X, X_sparse);
        SetCol(X_sparse, i, M);
    }


    //! Conversion from 'RowMajor' to 'ArrayRowSparse' format.
    /*!
      \param[in] A_dense the 'RowMajor' matrix to be converted.
      \param[out] A_array the matrix \a A_dense  stored in the
      'ArrayRowSparse' format.
    */
    template <class T, class Allocator>
    void ConvertDenseToArrayRowSparse(
        const Matrix<T, General, RowMajor, Allocator>& A_dense,
        Matrix<T, General, ArrayRowSparse, Allocator>& A_array)
    {
        int m = A_dense.GetM();
        int n = A_dense.GetN();

        A_array.Reallocate(m, n);

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
            {
                T value = A_dense(i, j);
                if (value != T(0))
                    A_array.AddInteraction(i, j, value);
            }
    }


    //! Conversion from 'ArrayRowSparse' to 'RowMajor' format.
    /*!
      \param[in] A_array the matrix to be converted.
      \param[out] A_dense the matrix \a A_array stored in the
      'RowMajor' format.
    */
    template <class T, class Allocator>
    void ConvertArrayRowSparseToDense(
        const Matrix<T, General, ArrayRowSparse, Allocator>& A_array,
        Matrix<T, General, RowMajor, Allocator>& A_dense)
    {
        int m = A_array.GetM();
        int n = A_array.GetN();

        A_dense.Reallocate(m, n);

        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                A_dense(i, j) = A_array(i, j);
    }


    //! Conversion from 'RowSparse' to 'ArrayRowSparse' format.
    /*!
      \param[in] A the 'RowSparse' matrix to be converted.
      \param[out] A_array the matrix \a A stored in the
      'ArrayRowSparse' format.
    */
    template <class T, class Allocator>
    void ConvertRowSparseToArrayRowSparse(
        const Matrix<T, General, RowSparse, Allocator>& A,
        Matrix<T, General, ArrayRowSparse, Allocator>& A_array)
    {
        int m = A.GetM();
        int n = A.GetN();

        A_array.Reallocate(m, n);

        T* data = A.GetData();
        int* ptr = A.GetPtr();
        int* column = A.GetInd();

        for (int i = 0; i < m; i++)
            for (int j = ptr[i]; j < ptr[i + 1]; j++)
                A_array.AddInteraction(i, column[j], data[j]);
    }


    //! Conversion from 'VectSparse' to 'VectFull' format.
    /*!
      \param[in] V_sparse the 'VectSparse' vector to be converted.
      \param[out] V_dense the vector \a V_sparse  stored in the 'VectFull'
      format.
    */
    template <class T, class Allocator>
    void ConvertSparseToDense(const Vector<T, VectSparse, Allocator>&
                              V_sparse, Vector<T, VectFull, Allocator>&
                              V_dense)
    {
        V_dense.Fill(T(0.));
        for (int k = 0; k < V_sparse.GetM(); k++)
            V_dense(V_sparse.Index(k)) = V_sparse.Value(k);
    }


    //! Conversion from 'VectFull' to 'VectSparse' format.
    /*!
      \param[in] V_dense the 'VectFull' vector to be converted.
      \param[out] V_sparse the vector \a V_dense  stored in the 'VectSparse'
      format.
    */
    template <class T, class Allocator>
    void ConvertDenseToSparse(const Vector<T, VectFull, Allocator> V_dense,
                              Vector<T, VectSparse, Allocator>& V_sparse)
    {
        V_sparse.Clear();
        for (int k = 0; k < V_dense.GetLength(); k++)
        {
            T value = V_dense(k);
            if (value != T(0.))
                V_sparse.AddInteraction(k, value);
        }
    }


    //! Fills a given matrix with a given a value.
    /*!
      \param[in] value the value to fill the matrix with.
      \param[out] M the matrix to be filled.
    */
    template <class T>
    void Fill(T value, Matrix<T, Symmetric, RowSymSparse>& M)
    {
        Vector<T> working_vector;
        working_vector.SetData(M.GetDataSize(), M.GetData());
        working_vector.Fill(T(value));
        working_vector.Nullify();
    }



    //! Adds a small matrix to a bigger one at a given position.
    /*!
      \param[in] c the value to multiply the matrix to insert with.
      \param[in] pi the row to place the matrix in.
      \param[in] pj the column to place the matrix in.
      \param[in] B the matrix to be inserted.
      \param[out] A the matrix to be filled.
    */
    template <class T0,
              class T1, class Prop1, class Storage1, class Allocator1,
              class T2, class Storage2, class Allocator2>
    void AddMatrixPosition(T0 c, int pi, int pj,
                           const Matrix<T1, Prop1, Storage1, Allocator1>& B,
                           Matrix<T2, Symmetric, Storage2, Allocator2>& A)
    {
        if (A.GetM() - pi  < B.GetM())
            throw ErrorArgument("void AddMatrixPosition(T c, int i, int j, "
                                "Matrix B, Matrix A)",
                                string("Number of rows ") + to_str(A.GetM()) +
                                " in matrix A is not enough to copy another "
                                + to_str(B.GetM()) + " rows.");
        if (A.GetN() - pj  < B.GetN())
            throw ErrorArgument("void AddMatrixPosition(T c, int i, int j, "
                                "Matrix B, Matrix A)",
                                string("Number of columns ") +
                                to_str(A.GetN()) + " in matrix A is not"
                                " enough to copy another "
                                + to_str(B.GetN()) + " columns.");
        for (int i = 0; i < B.GetM(); i++)
            for (int j = i; j < B.GetN(); j++)
                A.Val(pi + i, pj + j) += c * B.Val(i, j);
    }


    //! Adds a small matrix to a bigger one at a given position.
    /*!
      \param[in] c the value to multiply the matrix to insert with.
      \param[in] pi the row to place the matrix in.
      \param[in] pj the column to place the matrix in.
      \param[in] B the matrix to be inserted.
      \param[out] A the matrix to be filled.
    */
    template <class T0,
              class T1, class Prop1, class Storage1, class Allocator1,
              class T2, class Prop2, class Storage2, class Allocator2>
    void AddMatrixPosition(T0 c, int pi, int pj,
                           const Matrix<T1, Prop1, Storage1, Allocator1>& B,
                           Matrix<T2, Prop2, Storage2, Allocator2>& A)
    {
        if (A.GetM() - pi  < B.GetM())
            throw ErrorArgument("void AddMatrixPosition(T c, int i, int j, "
                                "Matrix B, Matrix A)",
                                string("Number of rows ") + to_str(A.GetM()) +
                                " in matrix A is not enough to copy another "
                                + to_str(B.GetM()) + " rows.");
        if (A.GetN() - pj  < B.GetN())
            throw ErrorArgument("void AddMatrixPosition(T c, int i, int j, "
                                "Matrix B, Matrix A)",
                                string("Number of columns ") +
                                to_str(A.GetN()) + " in matrix A is not"
                                " enough to copy another "
                                + to_str(B.GetN()) + " columns.");

        for (int i = 0; i < B.GetM(); i++)
            for (int j = 0; j < B.GetN(); j++)
                A.Val(pi + i, pj + j) += c * B.Val(i, j);
    }


    //! Gets a vector referencing a line of a given matrix.
    /*!
      \param[in] M a given matrix.
      \param[in] i a given line index.
      \return a vector pointing on the line \a i of \a M.
    */
    template <class T, template <class U> class Allocator>
    void GetRowPointer(const Matrix<T, General, RowMajor, Allocator<T> >& M,
                       int i, Vector<T, VectFull, Allocator<T> >& V)
    {
        if (i < 0 || i >= M.GetM())
            throw ErrorArgument("void GetLineVector(Matrix<T, General, "
                                "RowMajor, Allocator<T> >& M, int i, "
                                "Vector<T, VectFull, Allocator<T> >& V)",
                                string("Index should be in [0, ")
                                + to_str(M.GetM() - 1)
                                + "], but is equal to " + to_str(i) + ".");
        V.SetData(M.GetN(), &M.GetData()[i * M.GetN()]);
    }


#ifdef VERDANDI_WITH_PETSC
    //! Reallocates memory of a given matrix.
    /*!
      \param[in, out] A a given matrix. On exit, the matrix is a \a i x \a j
      matrix.
      \param[in] i new number of rows.
      \param[in] j new number of columns.
      \param[in] model the model that gives, if \a A is distributed, the
      number of local rows.
    */
    template <class T0, class Allocator0,
              class Model>
    void Reallocate(Matrix<T0, General, PETScMPIDense, Allocator0>& A, int i,
                    int j, const Model& model)
    {
        A.Reallocate(i, j, model.GetLocalNstate());
    }
#endif


    //! Reallocates memory of a given matrix.
    /*!
      \param[in, out] A a given matrix. On exit, the matrix is a \a i x \a j
      matrix.
      \param[in] i new number of rows.
      \param[in] j new number of columns.
      \param[in] model this argument is not used by this function.
    */
    template <class T0, class Prop0, class Storage0, class Allocator0,
              class Model>
    void Reallocate(Matrix<T0, Prop0, Storage0, Allocator0>& A, int i, int j,
                    const Model& model)
    {
        A.Reallocate(i, j);
    }


#ifdef VERDANDI_WITH_PETSC
    //! Reallocates memory of a given vector.
    /*!
      \param[in, out] V a given vector. On exit, the vector is of size i.
      \param[in] i new vector size.
      \param[in] model the model that gives, if \a V is distributed, the
      number of local elements.
    */
    template <class T0, class Allocator0,
              class Model>
    void Reallocate(Vector<T0, PETScPar, Allocator0>& V, int i,
                    const Model& model)
    {
        V.Reallocate(i, model.GetLocalNstate());
    }
#endif


    //! Reallocates memory of a given vector.
    /*!
      \param[in, out] V a given vector. On exit, the vector is of size i.
      \param[in] i new vector size.
      \param[in] model this argument is not used by this function.
    */
    template <class T0, class Storage0, class Allocator0,
              class Model>
    void Reallocate(Vector<T0, Storage0, Allocator0>& V, int i,
                    const Model& model)
    {
        V.Reallocate(i);
    }


} // namespace Verdandi.


#endif
