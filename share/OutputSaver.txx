// Copyright (C) 2008-2009 INRIA
// Author(s): Marc Fragu, Vivien Mallet
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


#ifndef VERDANDI_FILE_OUTPUTSAVER_OUTPUT_SAVER_TXX
#define VERDANDI_FILE_OUTPUTSAVER_OUTPUT_SAVER_TXX

#include "Variable.hxx"


namespace Verdandi
{

    /////////////
    // METHODS //
    /////////////


    //! Writes the variable in the file with the chosen format.
    /*! The variable \a variable_name is saved in its dedicated file with its
      saving mode, providing \a time satisfies the proper conditions.
      \param[in] x value of the variable.
      \param[in] time time.
      \param[in] variable_name name of the variable to be saved.
    */
    template <class S>
    void OutputSaver::Save(const S& x, double time, string variable_name)
    {
        if (is_active_ && (save_period_ == 0.
            || is_multiple(time, save_period_, time_tolerance_)))
        {
            Save(x, variable_name);
#ifdef VERDANDI_WITH_HDF5
            if (mode_ == "HDF" && IsSaved(variable_name))
            {
                map<string, Variable>::iterator im;

                im = variable_list_.find(variable_name);
                Variable& variable = im->second;
                string dataset_time_name =
                    variable_name + "_time";
                Vector<double> vector_time(1);
                vector_time(0) = time;
                WriteHDF5(vector_time, variable.GetFile(), group_,
                          dataset_time_name);
            }
#endif
        }
    }


    //! Writes the variable in the file with the chosen format.
    /*! The variable \a variable_name is saved in its dedicated file with its
      saving mode.
      \param[in] x value of the variable.
      \param[in] variable_name name of the variable to be saved.
    */
    template <class S>
    void OutputSaver::Save(const S& x, string variable_name)
    {
        if (!is_active_ || !IsSaved(variable_name))
            return;

#ifdef VERDANDI_WITH_PETSC
        /* In this case, the simulation is parallel, but x is sequential
         (duplicated on each processes): only process 0 perform the saving. */
        int rank;
        int ierr;
        ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        CHKERRABORT(MPI_COMM_WORLD, ierr);
        if (rank != 0)
            return;
#endif

        map<string, Variable>::iterator im;
        im = variable_list_.find(variable_name);
        Variable& variable = im->second;

        // In case the mode has not been set yet.
        SetVariable<S>(variable);
        if (variable.HasToEmptyFile())
            Empty(variable_name);
        if (variable.GetMode() == "text")
            WriteText(x, variable.GetFile());
        else if (variable.GetMode() == "binary")
            WriteBinary(x, variable.GetFile());
#ifdef VERDANDI_WITH_HDF5
        else if (variable.GetMode() == "HDF")
        {
            WriteHDF5(x, variable.GetFile(), group_, variable_name);
        }
#endif
    }


#ifdef VERDANDI_WITH_PETSC
    //! Writes the variable in the file with the chosen format.
    /*! The variable \a variable_name is saved in its dedicated file with its
      saving mode.
      \param[in] x value of the variable.
      \param[in] variable_name name of the variable to be saved.
    */
    template <class T, class Allocator>
    void OutputSaver::Save(const Vector<T, PETScPar, Allocator>& x,
                           string variable_name)
    {
        if (!is_active_ || !IsSaved(variable_name))
            return;

        map<string, Variable>::iterator im;
        im = variable_list_.find(variable_name);
        Variable& variable = im->second;

        // In case the mode has not been set yet.
        SetVariable<Vector<T, PETScPar, Allocator> >(variable);

        if (variable.HasToEmptyFile())
        {
            int rank;
            int ierr;
            ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            CHKERRABORT(MPI_COMM_WORLD, ierr);
            variable.SetFile(variable.GetFile()
                             + "-processor_" + to_str(rank));
            Empty(variable_name);
        }

        if (variable.GetMode() == "text")
            WriteText(x, variable.GetFile());
        else if (variable.GetMode() == "binary")
            WriteBinary(x, variable.GetFile());
    }
#endif



    //! Writes \a x in a text file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <class S>
    void OutputSaver::WriteText(const S& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteText(const S& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        x.WriteText(file);
        file << endl;
        file.close();
    }

#ifdef VERDANDI_WITH_HDF5
    //! Writes \a x in a HDF5 file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
      \param[in] group_name name of the group \a x must be stored in.
      \param[in] dataset_name name of the dataset \a x must be stored in.
      \param[in] time corresponding time of the variable.
    */
    template <class S>
    void OutputSaver::WriteHDF5(const S& x, string file_name,
                                string group_name, string dataset_name) const
    {
        DISP(file_name);
        x.WriteHDF5(file_name, group_name, dataset_name);
    }
#endif

    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
    */
    template <class S>
    void OutputSaver::WriteBinary(const S& x, string file_name) const
    {
        ofstream file(file_name.c_str(), ofstream::app | ofstream::binary);
#ifdef VERDANDI_CHECK_IO
        if (!file)
            throw ErrorIO("WriteBinary(const S& x , string file_name)",
                          "Cannot open file \"" + file_name + "\"." );
#endif
        x.Write(file);
        file.close();
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
      \warning This method is undefined: it throws an exception in any case.
    */
    template <class T, class Prop, class Allocator>
    void OutputSaver
    ::WriteBinary(const Matrix<T, Prop, RowSparse, Allocator>& x,
                  string file_name) const
    {
        throw ErrorUndefined("WriteBinary(const Matrix<T, Prop, RowSparse, "
                             "Allocator>& x, string file_name)");
    }


    //! Writes \a x in a binary file.
    /*!
      \param[in] x variable to be written.
      \param[in] file_name output filename.
      \warning This method is undefined: it throws an exception in any case.
    */
    template <class T, class Prop, class Allocator>
    void OutputSaver
    ::WriteBinary(const Matrix<T, Prop, ColSparse, Allocator>& x,
                  string file_name) const
    {
        throw ErrorUndefined("WriteBinary(const Matrix<T, Prop, ColSparse, "
                             "Allocator>& x, string file_name)");
    }


    //! Empties the output file associated with a variable.
    /*!
      \param[in] variable_name the name of the variable whose file should be
      emptied.
    */
    template <class S>
    void OutputSaver::Empty(string variable_name)
    {
        if (!is_active_ || !IsSaved(variable_name))
            return;

        map<string, Variable>::iterator im;
        im = variable_list_.find(variable_name);

        if (im == variable_list_.end())
            return;

        if (im->second.GetMode().empty())
            SetVariable<S>(im->second);

        ofstream output_stream(im->second.GetFile().c_str());
        output_stream.close();
        im->second.HasToEmptyFile(false);
    }


    ////////////////////
    // PRIVATE METHOD //
    ////////////////////


    //! Sets the parameters of a variable.
    /*!
      \param[in,out] variable the variable whose parameters are to be
      adjusted.
    */
    template <class S>
    void OutputSaver::SetVariable(Variable& variable)
    {
        if (variable.GetMode().empty())
            // In this case, the mode is not set yet, so the default mode
            // 'mode_' is selected.
        {
            variable.SetMode(mode_);
            SetVariableFile(variable);
        }
    }


} // namespace Verdandi.



#endif
