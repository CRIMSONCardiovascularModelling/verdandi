// Copyright (C) 2010 INRIA
// Author(s): Marc Fragu, Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_OBSERVATIONMANAGER_OBSERVATIONMANAGERTEMPLATE_HXX


namespace Verdandi
{


    /////////////////////////////////
    // OBSERVATIONMANAGER TEMPLATE //
    /////////////////////////////////


    //! This class is a template of observation manager.
    class ObservationManagerTemplate: public VerdandiBase
    {
    public:
#ifdef VERDANDI_TANGENT_LINEAR_OPERATOR_SPARSE
        //! Type of the tangent linear operator.
        typedef Matrix<double, General, RowSparse> tangent_linear_operator;
#else
        //! Type of the tangent linear operator.
        typedef Matrix<double> tangent_linear_operator;
#endif
        //! Type of a row of the tangent linear operator.
        typedef Vector<double> tangent_linear_operator_row;

#ifdef VERDANDI_OBSERVATION_ERROR_SPARSE
        //! Type of the observation error covariance matrix.
        typedef Matrix<double, General, RowSparse> error_variance;
#else
        //! Type of the observation error covariance matrix.
        typedef Matrix<double> error_variance;
#endif

        //! Type of the observation vector.
        typedef Vector<double> observation;


    public:
        // Constructors and destructor.
        ObservationManagerTemplate();
        ~ObservationManagerTemplate();

        // Initialization.
        template <class Model>
        void Initialize(Model& model, string configuration_file);

        void DiscardObservation(bool discard_observation);
        template <class Model>
        void SetTime(Model& model, double time);


        /////////////////
        // PARALLEL OM //
        /////////////////

#ifdef VERDANDI_WITH_MPI
        void SetMPICommunicator(MPI_Comm& mpi_communicator);
#endif

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
        template <class state, class mat>
        void GetNudgingMatrix(const state& x, mat& M) const;

        ///////////////
        // OPERATORS //
        ///////////////


        template <class state>
        void ApplyOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyTangentLinearOperator(const state& x, observation& y) const;
        double GetTangentLinearOperator(int i, int j) const;
        tangent_linear_operator_row& GetTangentLinearOperatorRow(int row);
        const tangent_linear_operator& GetTangentLinearOperator() const;

        template <class state>
        void ApplyAdjointOperator(const state& x, observation& y) const;

        double GetErrorVariance(int i, int j) const;
        const error_variance& GetErrorVariance() const;
        const error_variance& GetErrorVarianceInverse() const;

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_OBSERVATIONMANAGER_OBSERVATIONMANAGERTEMPLATE_HXX
#endif
