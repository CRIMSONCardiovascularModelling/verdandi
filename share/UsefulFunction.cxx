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


#ifndef VERDANDI_FILE_SHARE_USEFULFUNCTION_CXX
#define VERDANDI_FILE_SHARE_USEFULFUNCTION_CXX


#include "VerdandiHeader.hxx"


namespace Verdandi
{



    /*! \brief Returns the position in a multidimensional grid that is
      associated with a global index. */
    /*! A global index gives the position in a grid with a single integer. For
      example, in 2D, the global index of the grid point \f$(i, j)\f$ is \f$i
      \times N + j\f$ if there are N points along the second dimension. This
      function would return \f$(i, j)\f$ in \a position from (\a index =) \f$i
      \times N + j\f$, with \a shape set to \f$(M, N)\f$.
      \param index the global index.
      \param shape dimensions of the grid.
      \param position position in the grid.
    */
    void get_position(int index, const Vector<int>& shape,
                      Vector<int>& position)
    {
        int d;

        int length = shape.GetLength();

        if (index < 0)
            throw ErrorArgument("get_position(int, Vector<int>, Vector<int>)",
                                "Wrong index: " + to_str(index) + ".");
        if (length == 0)
            throw ErrorArgument("get_position(int, Vector<int>, Vector<int>)",
                                "The shape vector is empty.");

        position.Reallocate(length);

        if (length == 1)
            if (index >= shape(0))
                throw ErrorArgument
                    ("get_position(int, Vector<int>, Vector<int>&)",
                     "The shape vector is " + Seldon::to_str(shape)
                     + ", but the index is " + to_str(index) + ".");
            else
            {
                position(0) = index;
                return;
            }

        Vector<int> size(length - 1);
        size(length - 2) = shape(length - 1);
        for (d = length - 3; d >= 0; d--)
            size(d) = size(d + 1) * shape(d + 1);

        for (d = 0; d < length - 1; d++)
        {
            position(d) = index / size(d);
            index = index - position(d) * size(d);
        }
        position(length - 1) = index;
    }


    /*! \brief Returns the global index in a multidimensional grid that is
      associated with a local position. */
    /*! A global index gives the position in a grid with a single integer. For
      example, in 2D, the global index of the grid point \f$(i, j)\f$ is \f$i
      \times N + j\f$ if there are N points along the second dimension. This
      function returns \f$i \times N + j\f$ from \a position set to \f$(i,
      j)\f$, if \a shape is \f$(M, N)\f$.
      \param[in] shape dimensions of the grid.
      \param[in] position position in the grid.
      \return index The global index.
    */
    int get_position(const Vector<int>& shape, const Vector<int>& position)
    {
        if (position.GetLength() == 0)
            return 0;

        int index = position(0);
        for (int d = 1; d < position.GetLength(); d++)
            index = shape(d) * index + position(d);

        return index;
    }


    //! Sets a string.
    /*!
      \param[in] s input string.
      \param[out] out output string, equal to \a s on exit.
    */
    void convert(const string& s, string& out)
    {
        out = s;
    }


    //! Checks whether a string is a number.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is a number, false otherwise.
    */
    bool is_num(const string& str)
    {
        if (str == "")
            return false;

        bool mant, mant_a, mant_b, exp;
        string::size_type pos;
        string m, e, m_a, m_b;

        pos = str.find_first_of("eE");
        // Mantissa.
        m = str.substr(0, pos);
        // Exponent.
        e = pos == string::npos ? "" : str.substr(pos + 1);

        exp = pos != string::npos;

        pos = m.find_first_of(".");
        // Mantissa in the form: [m_a].[m_b].
        m_a = m.substr(0, pos);
        // Exponent.
        m_b = pos == string::npos ? "" : m.substr(pos + 1);

        mant = m != "" && m != "-" && m != "+";
        mant_a = m_a != "" && m_a != "-" && m_a != "+";
        mant_b = m_b != "";

        return (mant
                && ((mant_a || mant_b)
                    && (!mant_a || is_integer(m_a))
                    && (!mant_b || is_unsigned_integer(m_b)))
                && (!exp || is_integer(e)));
    }


    //! Checks whether a string is an integer.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is an integer, false otherwise.
    */
    bool is_integer(const string& str)
    {
        bool ans;

        ans = (str.size() > 0 && isdigit(str[0]))
            || (str.size() > 1 && (str[0] == '+' || str[0] == '-'));

        unsigned int i(1);
        while (i < str.size() && ans)
        {
            ans = ans && isdigit(str[i]);
            i++;
        }

        return ans;
    }


    //! Checks whether a string is an unsigned integer.
    /*!
      \param[in] str string to be checked.
      \return True if \a str is an unsigned integer, false otherwise.
    */
    bool is_unsigned_integer(const string& str)
    {
        bool ans(str.size() > 0);

        unsigned int i(0);
        while (i < str.size() && ans)
        {
            ans = ans && isdigit(str[i]);
            i++;
        }

        return ans;
    }


    //! Trims off a string.
    /*!
      Removes delimiters at each edge of the string.
      \param[in] str string to be trimmed off.
      \param[in] delimiters characters to be removed.
      \return \a str trimmed off.
    */
    string trim(string str, string delimiters)
    {
        string::size_type index_end = str.find_last_not_of(delimiters);
        string::size_type index_beg = str.find_first_not_of(delimiters);

        if (index_beg == string::npos)
            return "";

        return str.substr(index_beg, index_end - index_beg + 1);
    }


    //! Splits a string.
    /*!
      The string is split according to delimiters.
      \param[in] str string to be split.
      \param[in] delimiters (optional) delimiters. Default: " \n\t".
      \return A vector containing elements of the string.
    */
    vector<string> split(string str, string delimiters)
    {
        vector<string> vect;
        split(str, vect, delimiters);
        return vect;
    }


    //! Finds and replace a substring.
    /*!
      \param str base string.
      \param old_str substring to be replaced.
      \param new_str substring to be put in place of 'old_str'.
      \return 'str' where 'old_str' was replaced by 'new'str'.
    */
    string find_replace(string str, string old_str, string new_str)
    {
        string::size_type index = str.find(old_str);

        while (index != string::npos)
        {
            str.replace(index, old_str.size(), new_str);
            index = str.find(old_str, index + new_str.size());
        }

        return str;
    }


    //! Converts a string to upper-case string.
    /*!
      \param str string to be converted.
      \return \a str in upper case.
    */
    string upper_case(string str)
    {
        string upper(str);
        std::transform(upper.begin(), upper.end(), upper.begin(),
                       (int(*)(int))toupper);
        return upper;
    }


#ifdef VERDANDI_WITH_MPI
    //! Builds the row and the column communicator of the current process.
    /*! Builds the row and the column communicator of the current process in a
      grid of size \a Nrow * \a Ncol.
      \param[in] Nrow the number of row of the grid.
      \param[in] Ncol the number of column of the grid.
      \param[out] row_communicator current process MPI row communicator.
      \param[out] col_communicator current process MPI column communicator.
    */
    void SetGridCommunicator(int Nrow, int Ncol, MPI_Comm *row_communicator,
                             MPI_Comm *col_communicator)
   {
       int irow, icol, color, key, world_rank;
       // Establishes the row and column to which this processor belongs.
       MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
       irow = world_rank % Nrow;
       icol = world_rank / Nrow;
       color = irow;
       key = world_rank;
       MPI_Comm_split(MPI_COMM_WORLD, color, key, row_communicator);
       color = icol;
       MPI_Comm_split(MPI_COMM_WORLD, color, key, col_communicator);
   }
#endif


} // namespace Verdandi.


#endif
