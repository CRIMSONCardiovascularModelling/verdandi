// Copyright (C) 2015 INRIA
// Author(s): Gautier Bureau
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


#ifndef VERDANDI_FILE_OBSERVATIONMANAGER_LEVELSETOBSERVATIONMANAGER_HXX
#define VERDANDI_FILE_OBSERVATIONMANAGER_LEVELSETOBSERVATIONMANAGER_HXX


namespace Verdandi
{


    /////////////////////////////////
    // LEVELSETOBSERVATIONMANAGER  //
    /////////////////////////////////


    //! This class is a template of observation manager.
    template <class Model>
    class LevelSetObservationManager: public VerdandiBase
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
        //! Type of the observation vector.
        typedef Vector<double> observation_vector;
        //! Type of the observation vector 2.
        typedef Vector2<double> observation_vector2;

        //! Type of the time vector.
        typedef Vector<double> time_vector;

    private:

        Model* model_ = nullptr;

    protected:

        /*** Observation file structure ***/

        //! File that stores the observations.
        string observation_file_;
        //! Type of the file.
        string observation_file_type_;
        //! Path to the dataset where observations are stored (HDF5 filetype).
        string observation_dataset_path_;
        //! How are stored the observations.
        string observation_type_;
        //! Total number of observations at current time.
        int Nobservation_;
        //! Size in bytes of an observation vector.
        size_t Nbyte_observation_;
        //! Is the time interval between two observations constant?
        bool is_delta_t_constant_;
        //! Period with which observations are available (if constant).
        double Delta_t_;
        /*! Times at which observations are available (if the time interval
         between two observations is not constant). */
        Vector<double> observation_time_;

        /*! Path to the observation times (needed if the time interval between
         two observations is not constant). */
        string observation_time_file_;
        //! Period with which available observations are actually loaded.
        int Nskip_;
        //! First time at which observations are available.
        double initial_time_;
        //! Final time at which observations are available.
        double final_time_;

        /*** Observation times ***/

        //! Requested time.
        double time_;
        //! Available observation time of the time interval.
        time_vector available_time_;

        /*** Model domain ***/

        //! The size of a model state.
        int Nstate_model_;

        //! Threshold activation of the model.
        double threshold_activation_;
        //! Epsilon ofr Dirac approximation
        double epsilon_;
        //! With topological gradient ?
        bool with_topological_gradient_;
        //! Topological gradient coefficient.
        double topological_gradient_coeff_;

        //! Stimulation radius in 1D.
        double stimulation_radius_1D_;
        //! Delta_x in 1D.
        double delta_x_1D_;
        //! Is the problem in 1D?
        bool case_in_1D_;

        //! If observations not availabe at each time step inteporlate
        //observations.
        bool interpolate_observations_;

        //! Observations currently stored.
        observation_vector2 observation2_;
        //! Innovation currently stored.
        observation innovation_;
        //! Innovation contributions.
        observation_vector2 innovation2_;
        //! Vector that stores the current level set.
        observation level_set_;
        //! Vector that stores the current approximation of delta.
        observation delta_;
        //! Current sign mult for the topological gradient.
        observation sign_mult_;
        //! Vector filled with the threshold value.
        observation threshold_activation_vector_;
        //! Topological gradient currently stored.
        observation topological_gradient_;

        int iteration_;


    public:
        // Constructors and destructor.
        LevelSetObservationManager();
        ~LevelSetObservationManager();

        // Initialization.
        void Initialize(Model& model, string configuration_file);

        void SetTime(Model& model, double time);
        void SetTime(double time);
        void SetAvailableTime(double time, time_vector& available_time);
        double GetTime() const;
        time_vector& GetAvailableTime();


        /////////////////
        // PARALLEL OM //
        /////////////////


#ifdef VERDANDI_WITH_MPI
        void SetMPICommunicator(MPI_Comm& mpi_communicator);
#endif


        /////////////////
        // OBSERVATION //
        /////////////////


        observation_vector2& GetObservation();
        observation_vector2& GetVectorObservation2();
        observation_vector2& GetVectorInnovation2();
        observation& GetLevelSet();
        const observation& GetThresholdActivationVector() const;
        observation& GetSignMult();
        observation& GetDelta();
        observation& GetTopologicalGradient();


        ////////////////
        // INNOVATION //
        ////////////////


        template <class state>
        observation& GetInnovation(const state& x);
        observation& GetVectorInnovation();


        ////////////
        // ACCESS //
        ////////////

        bool HasObservation() const;
        bool InterpolateObservations() const;
        int GetNobservation() const;
        double GetTresholdActivation() const;
        bool CaseIn1D() const;
        double GetStimulationRadius1D() const;
        double GetDeltaX1D() const;
        double GetEpsilon() const;
        bool WithTopologicalGradient() const;
        double GetTopologicalGradientCoefficient() const;
        template <class state, class mat>
        void GetNudgingMatrix(const state& x, mat& M) const;
        const Model& GetModel() const;
        Model& GetNonCstModel();

        ///////////////
        // OPERATORS //
        ///////////////

        template <class state>
        void ApplyOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyTangentLinearOperator(const state& x, observation& y) const;

        template <class state>
        void ApplyAdjointOperator(const state& x, observation& y) const;

        void DiracApproximation(const observation& level_set,
                                const double epsilon);

        template <class state>
        double GetMeanValueOnDomain(const state& x,
                                    const observation& observation);

        template <class state>
        double
        GetMeanValueOnComplementaryDomain(const state& x,
                                          const observation& observation);

        void ReadObservation(ifstream& file_stream, double time, int variable,
                             observation_vector& observation);

        void ReadObservation(const time_vector& available_time,
                             observation_vector2& observation2);

        string GetName() const;

        void Message(string message);

    };


} // namespace Verdandi.


#endif // VERDANDI_FILE_OBSERVATIONMANAGER_LEVELSETOBSERVATIONMANAGER_HXX
