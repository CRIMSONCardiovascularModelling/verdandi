// Copyright (C) 2011-2012 INRIA
// Author(s): Vivien Mallet, Kévin Charpentier
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


#ifndef VERDANDI_FILE_OBSERVATION_MANAGER_PYTHONOBSERVATIONMANAGER_HXX

#define PYTHON_EVAL_HELPER <python
#define EVAL(v) v
#include EVAL(PYTHON_EVAL_HELPER)EVAL(VERDANDI_PYTHON_VERSION)/Python.h>

#include<numpy/arrayobject.h>


namespace Verdandi
{


    //! This class is an interface for Python observation managers.
    class PythonObservationManager: public VerdandiBase
    {
    public:
        //! Type of the tangent linear operator.
        typedef Matrix<double> tangent_linear_operator;
        //! Type of a row of the tangent linear operator.
        typedef Vector<double> tangent_linear_operator_row;
        //! Type of the observation error covariance matrix.
        typedef Matrix<double> error_variance;
        //! Type of the observation vector.
        typedef Vector<double> observation;

    protected:

        //! Instance of the Python model class.
        PyObject *pyObservationManagerInstance_;
        //! Name of the python model file.
        PyObject *pyObservationManagerFile_;
        //! Imported model module from Python.
        PyObject *pyObservationManagerModule_;
        //! Dictionary object that describes the module's namespace.
        PyObject *pyObservationManagerDict_;
        //! Model class from the module.
        PyObject *pyObservationManagerClass_;
        //! Has the Python module been initialized?
        bool is_module_initialized_;

        //! Name of the Python module that contains the observation manager.
        string module_;
        /*! \brief Directory to include in "sys.path" for the module to be
          found. If no directory is to be added to "sys.path", this
          attribute remains empty. */
        string directory_;
        //! Name of the Python observation manager class.
        string class_name_;

        //! Observation currently stored.
        observation observation_;
        //! Innovation currently stored.
        observation innovation_;

        /*** Observation operator ***/

        //! Tangent operator matrix (H).
        tangent_linear_operator tangent_operator_matrix_;

        //! Observation error covariance matrix (R).
        error_variance error_variance_;
        //! Inverse of the observation error covariance matrix (R).
        error_variance error_variance_inverse_;

        //! Index of the row of H  currently stored.
        int current_row_;
        //! Value of the row of H currently stored.
        tangent_linear_operator_row tangent_operator_row_;

    public:
        // Constructor and destructor.
        PythonObservationManager();
        ~PythonObservationManager();

        // Initialization.
        template <class Model>
        void Initialize(const Model& model, string configuration_file);

        void DiscardObservation(bool discard_observation);
        template <class Model>
        void SetTime(const Model& model, double time);


        /////////////////
        // OBSERVATION //
        /////////////////


        observation& GetObservation();


        ////////////////
        // INNOVATION //
        ////////////////


        template <class state>
        observation& GetInnovation(const state& x);


        ////////////
        // ACCESS //
        ////////////


        bool HasObservation() const;
        int GetNobservation() const;


        ///////////////
        // OPERATORS //
        ///////////////


        template <class state>
        void ApplyOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyTangentLinearOperator(const state& x, observation& y) const;
        double GetTangentLinearOperator(int i, int j) const;
        tangent_linear_operator_row& GetTangentLinearOperatorRow(int row);
        const tangent_linear_operator& GetTangentLinearOperator();

        template <class state>
        void ApplyAdjointOperator(const state& x, observation& y) const;

        double GetErrorVariance(int i, int j) const;
        const error_variance& GetErrorVariance();
        const error_variance& GetErrorVarianceInverse();

        string GetName() const;
        void Message(string message);

    private:
        string ErrorMessageNotContiguous(string function_name) const;
    };


} // namespace Verdandi


#define VERDANDI_FILE_OBSERVATION_MANAGER_PYTHONOBSERVATIONMANAGER_HXX
#endif
