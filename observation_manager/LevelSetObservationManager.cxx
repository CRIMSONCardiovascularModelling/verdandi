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


#ifndef VERDANDI_FILE_OBSERVATIONMANAGER_LEVELSETOBSERVATIONMANAGER_CXX
#define VERDANDI_FILE_OBSERVATIONMANAGER_LEVELSETOBSERVATIONMANAGER_CXX


#include "LevelSetObservationManager.hxx"
#include <algorithm>
#include <cmath>


namespace Verdandi
{
    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Default constructor.
    template <class Model>
    LevelSetObservationManager<Model>::LevelSetObservationManager()
    {
    }


    //! Destructor.
    template <class Model>
    LevelSetObservationManager<Model>::~LevelSetObservationManager()
    {
        // Operations to be performed when the object is destroyed.
    }


    ////////////////////
    // INITIALIZATION //
    ////////////////////


    //! Initializes the observation manager.
    /*!
      \param[in] model model.
      \param[in] configuration_file configuration file.
      \tparam Model the model type; e.g. ShallowWater<double>
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::Initialize(Model& model, string configuration_file)
    {
        model_ = &model;

        VerdandiOps configuration(configuration_file);
        configuration.SetPrefix("observation.");

        Nstate_model_ = model.GetNstate();
        configuration.Set("level_set_observation_manager."
                          "threshold_activation",
                          threshold_activation_);

        Nobservation_ = Nstate_model_;

        iteration_ = 0;

        //! Computes once and for all a vector filled with the threshold that
        //will be needed.
        double threshold_activation = GetThresholdActivation();
        threshold_activation_vector_.Reallocate(Nobservation_);
        threshold_activation_vector_.Fill(threshold_activation);

        bool with_observation;
        configuration.Set("option.with_observation", with_observation);

        configuration.Set("Nskip", "v > 0", Nskip_);
        configuration.Set("initial_time", "", 0., initial_time_);
        configuration.Set("final_time", "", numeric_limits<double>::max(),
                          final_time_);

        configuration.Set("Delta_t_constant", is_delta_t_constant_);

        if (!is_delta_t_constant_)
        {
            configuration.Set("observation_time_file",
                              observation_time_file_);
            observation_time_.Read(observation_time_file_);
        }
        else
            configuration.Set("Delta_t", "v > 0", Delta_t_);

        configuration.SetPrefix("observation.level_set_observation_manager.");

        configuration.Set("case_in_1D", case_in_1D_);
        configuration.Set("stimulation_radius_1D", stimulation_radius_1D_);
        configuration.Set("delta_x_1D", delta_x_1D_);
        configuration.Set("epsilon", epsilon_);
        configuration.Set("with_topological_gradient",
                          with_topological_gradient_);
        configuration.Set("topological_gradient_coeff",
                          topological_gradient_coeff_);
        configuration.Set("interpolate_observations",
                          interpolate_observations_);

        if (!with_observation)
            return;

        if (WithTopologicalGradient())
            sign_mult_.Reallocate(Nobservation_);
        topological_gradient_.Reallocate(Nobservation_);

        configuration.SetPrefix("observation.");

        configuration.Set("file", observation_file_);
        configuration.Set("file_type", "", "binary", observation_file_type_);
        configuration.Set("observation_dataset_path", "", "",
                          observation_dataset_path_);
        configuration.Set("type", "", "state", observation_type_);

        time_ = numeric_limits<double>::min();

        if (observation_type_ == "state")
            Nbyte_observation_ = Nstate_model_ * sizeof(double) + sizeof(int);

        if (observation_type_ == "observation")
            Nbyte_observation_ = Nobservation_ * sizeof(double) + sizeof(int);

        if (is_delta_t_constant_)
        {
            int expected_file_size;
            expected_file_size = Nbyte_observation_
                * int((final_time_ - initial_time_)
                      / (Delta_t_ * double(Nskip_)) + 1.);

            int file_size;
            ifstream file_stream;
            file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("LevelSetObservationManager<Model>"
                              "::Initialize(model, configuration_file)",
                              "Unable to open file \""
                              + observation_file_ + "\".");
#endif

            file_stream.seekg(0, ios_base::end);
            file_size = file_stream.tellg() ;
            file_stream.close();

            if (expected_file_size > file_size)
                throw ErrorIO("LevelSetObservationManager<Model>"
                              "::Initialize(model, configuration_file)",
                              "Too few available observations, the size of \""
                              + observation_file_ + "\" must be greater than "
                              + to_str(expected_file_size) + " B.");
        }
        else
        {
            int file_size;
            ifstream file_stream;
            file_stream.open(observation_file_.c_str());

#ifdef VERDANDI_CHECK_IO
            // Checks if the file was opened.
            if (!file_stream.is_open())
                throw ErrorIO("LinearObservationManager"
                              "::Initialize(model, configuration_file)",
                              "Unable to open file \""
                              + observation_file_ + "\".");
#endif
            file_stream.close();
        }
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] model the model.
      \param[in] time a given time.
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::SetTime(Model& model, double time)
    {
        SetTime(time);
    }


    //! Sets the time of observations to be loaded.
    /*!
      \param[in] time a given time.
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::SetTime(double time)
    {
        time_ = time;
        SetAvailableTime(time_, available_time_);
    }


    //! Sets the available observation times at a given time.
    /*!
      \param[in] time the given time.
      \param[out] available_time the available observation times.
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::SetAvailableTime(double time,
                       time_vector& available_time)
    {
        available_time.Clear();

        if (is_delta_t_constant_)
        {
            double period = Delta_t_ * Nskip_;

            if (InterpolateObservations())
            {
                double t1, t2;

                if (is_multiple(time - initial_time_, period))
                {
                    t1 = time;

                    if (t1 <= final_time_)
                        available_time.PushBack(t1);

                    return;
                }
                else
                    t1 = initial_time_
                        + floor((time - initial_time_) / period) * period;

                t2 = t1 + period;

                if (t1 <= final_time_)
                    available_time.PushBack(t1);
                if (t2 <= final_time_)
                    available_time.PushBack(t2);

                return;

            }
            else
            {
                if (is_multiple(time - initial_time_, period))
                    available_time.PushBack(time);

                return;
            }
        }
        else
        {
            if (observation_time_(0) > time)
                return;

            if (InterpolateObservations())
            {
                for (int i = 1; i < observation_time_.GetM(); i++)
                    if (observation_time_(i) >= time)
                    {
                        available_time.PushBack(observation_time_(i - 1));
                        available_time.PushBack(observation_time_(i));
                        return;
                    }
                return;
            }
            else
            {
                for (int i = 1; i < observation_time_.GetM(); i++)
                    if ((observation_time_(i) - time) <= 1e-10)
                    {
                        available_time.PushBack(observation_time_(i));
                        return;
                    }
                return;
            }
        }

        Logger::Log<3>(*this, to_str(time) +
                       ", {" + to_str(available_time) + "}\n");
    }


    /////////////////
    // PARALLEL OM //
    /////////////////


#ifdef VERDANDI_WITH_MPI
    //! Sets the MPI communicator.
    /*!
      \param[in] mpi_communicator the MPI communicator to be set.
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::SetMPICommunicator(MPI_Comm& mpi_communicator)
    {
        throw ErrorUndefined("void LevelSetObservationManager<Model>::"
                             "SetMPICommunicator(MPI_Comm&)");
    }
#endif


    /////////////////
    // OBSERVATION //
    /////////////////


    //! Returns the observations.
    /*! This method is called after 'SetTime' set the time at which the
      observations are requested.
      \return The observation vector.
    */
    template <class Model>
    typename LevelSetObservationManager<Model>::observation_vector2&
    LevelSetObservationManager<Model>::GetObservation()
    {
        ReadObservation(available_time_, observation2_);

        return observation2_;
    }


    //! Builds observations associated with given times.
    /*!
      \param[in] available_time the given observation time vector.
      \param[out] observation2 the observations.
    */
    template <class Model>
    void LevelSetObservationManager<Model>::ReadObservation(
        const time_vector& available_time,
        observation_vector2& observation2)
    {
        ifstream file_stream;
        file_stream.open(observation_file_.c_str());
#ifdef VERDANDI_CHECK_IO
        // Checks if the file was opened.
        if (!file_stream.is_open())
            throw ErrorIO("LinearObservationManager"
                          "::LoadObservation(model)",
                          "Unable to open file \""
                          + observation_file_ + "\".");
#endif
        int Nt = available_time.GetSize();
        observation2.Reallocate(Nt);
        if (observation_file_type_ == "binary")
            for (int h = 0; h < Nt; h++)
                ReadObservation(file_stream, available_time(h), 0,
                                observation2(h));
#ifdef VERDANDI_WITH_HDF5
        else if (observation_file_type_ == "HDF")
            for (int h = 0; h < Nt; h++)
                ReadObservation(file_stream, available_time(h), 0,
                                observation_dataset_path_,
                                observation2(h));
#endif

        file_stream.close();
    }


    //! Reads observation from observation file given a time and a variable.
    /*!
      \param[in] file_stream the observation file stream.
      \param[in] time the time.
      \param[in] variable the variable.
      \param[out] observation the observations.
    */
    template <class Model>
    void LevelSetObservationManager<Model>::ReadObservation(
        ifstream& file_stream, double time, int variable,
        observation_vector&
        observation)
    {
        observation.Reallocate(Nobservation_);
        observation_vector input_data;

        streampos position;

        if (is_delta_t_constant_)
            position = (floor((time - initial_time_)
                              / Delta_t_ + 0.5) + variable)
                * Nbyte_observation_;
        else
        {
            int index;
            bool observation_available = false;
            for (index = 0; index < observation_time_.GetM(); index++)
                if (observation_time_(index) == time)
                {
                    observation_available = true;
                    break;
                }
            if (!observation_available)
                throw ErrorIO("LevelSetObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "No observation available at time "
                              + to_str(time) + ".");
            position = index * Nbyte_observation_;
        }
        file_stream.seekg(position);
        input_data.Read(file_stream);

        if (observation_type_ == "state")
        {
            if (input_data.GetSize() != Nstate_model_)
                throw ErrorIO("LevelSetObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "The observation type is 'state', so the whole"
                              " model state is supposed to be stored, but "
                              "the size of the observation read at time "
                              + to_str(time_) + " (Nread = "
                              + to_str(input_data.GetSize()) +
                              ") mismatches with the expected size (Nstate = "
                              + to_str(Nstate_model_) + ").");
            ApplyOperator(input_data, observation);
        }
        else
        {
            if (input_data.GetSize() != Nobservation_)
                throw ErrorIO("LevelSetObservationManager::ReadObservation"
                              "(ifstream& file_stream, double time, "
                              "int variable, LinearObservationManager"
                              "::observation_vector& observation) const",
                              "The observation type is 'observation', so "
                              "only observations are stored in the file, but"
                              " the size of the observation read at time "
                              + to_str(time_) + " (Nread = "
                              + to_str(input_data.GetSize())
                              + ") mismatches with the "
                              "expected size (Nobservation = "
                              + to_str(Nobservation_ ) + ").");
            Copy(input_data, observation);
        }
    }


    ////////////////
    // INNOVATION //
    ////////////////


    //! Returns an innovation.
    /*! This method is called after 'SetTime' set the time at which the
      innovation is requested.
      \param[in] state state vector.
      \return The innovation vector.
    */
    template <class Model>
    template <class state>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetInnovation(const state& x)
    {
        observation& innovation = GetVectorInnovation();
        innovation.Reallocate(Nobservation_);

        observation_vector2& innovation2 = GetVectorInnovation2();

        double alpha = 0.;

        if (InterpolateObservations())
        {
            time_vector available_time = GetAvailableTime();
            double time = GetTime();

            if (available_time.GetSize() > 1)
                alpha = (available_time(1) - time)
                    / (available_time(1) - available_time(0));
            else if ((time - available_time(0)) < 1.e-12)
                alpha = 1;
            else
                alpha = time / available_time(0);

            innovation2.Reallocate(2);
        }
        else
            innovation2.Reallocate(1);

        observation& level_set = GetLevelSet();
        level_set.Reallocate(Nobservation_);
        level_set.Copy(x);

        const observation& threshold_activation_vector =
            GetThresholdActivationVector();

        Add(-1., threshold_activation_vector, level_set);

        double epsilon = GetEpsilon();

        observation_vector2& observation_vector2 = GetObservation();

        int obs_num = observation_vector2.GetSize();

        for (int i_obs = 0; i_obs < obs_num ; ++i_obs)
        {
            observation& observation_vector = observation_vector2(i_obs);

            observation& innovation_vector = innovation2(i_obs);
            innovation_vector.Reallocate(Nobservation_);
            innovation_vector.Zero();

#if defined(CARTESIAN_GRID)
            GetDelta().Reallocate(Nobservation_);

            observation& delta = GetDelta();

            DiracApproximation(level_set, epsilon);

            double mean_value_on_domain =
                GetMeanValueOnDomain(x, observation_vector);

            double mean_value_on_complementary_domain =
                GetMeanValueOnComplementaryDomain(x, observation_vector);

            if (WithTopologicalGradient())
            {
                observation& sign_mult = GetSignMult();
                observation& topological_gradient = GetTopologicalGradient();

                double topological_gradient_coeff =
                    GetTopologicalGradientCoefficient();

                for (int i = 0 ; i < Nobservation_ ; ++i)
                {
                    topological_gradient(i) =
                        (observation_vector(i) - mean_value_on_domain) *
                        (observation_vector(i) - mean_value_on_domain) -
                        (observation_vector(i) -
                         mean_value_on_complementary_domain) *
                        (observation_vector(i) -
                         mean_value_on_complementary_domain);

                    if (level_set(i) * topological_gradient(i) > 0)
                        sign_mult(i) = 2.;
                    else if (level_set(i) * topological_gradient(i) < 0)
                        sign_mult(i) = 0.;
                    else
                        sign_mult(i) = 0.;

                    sign_mult(i) = -1. * topological_gradient_coeff
                        * sign_mult(i) * topological_gradient(i);
                }

                innovation_vector.Copy(sign_mult);
            }

            for (int i = 0 ; i < Nobservation_ ; ++i)
            {
                innovation_vector(i) -= delta(i) *
                    ((observation_vector(i) - mean_value_on_domain) *
                     (observation_vector(i) - mean_value_on_domain) -
                     (observation_vector(i) -
                      mean_value_on_complementary_domain) *
                     (observation_vector(i) -
                      mean_value_on_complementary_domain));
            }
#else
            Model& model = GetNonCstModel();

            const double threshold_activation = GetThresholdActivation();

            std::array<double,2> c1_c2 =
                model.GetMeanValues(observation_vector, threshold_activation);

            bool case_in_1D = CaseIn1D();

            if (WithTopologicalGradient())
            {
                double topological_gradient_coeff =
                    GetTopologicalGradientCoefficient();

                model.ComputeDirectCorrection(innovation_vector,
                                              observation_vector,
                                              level_set,
                                              epsilon,
                                              c1_c2[0],
                                              c1_c2[1],
                                              topological_gradient_coeff,
                                              true,
                                              case_in_1D);
            }
            else
                model.ComputeDirectCorrection(innovation_vector,
                                              observation_vector,
                                              level_set,
                                              epsilon,
                                              c1_c2[0],
                                              c1_c2[1],
                                              0.,
                                              false,
                                              case_in_1D);
#endif
        }

        innovation.Copy(innovation2(0));

        if (InterpolateObservations() && obs_num != 1)
        {
            Mlt(alpha, innovation);
            Add(1. - alpha, innovation2(1), innovation);
        }

        ++iteration_;

        return innovation;
    }


    //! Returns the approximation of a Dirac.
    template <class Model>
    void LevelSetObservationManager<Model>::DiracApproximation(
        const observation& level_set,
        const double epsilon)
    {
        delta_.Zero();

        if (CaseIn1D())
        {
            double epsilon_tilde = epsilon * 3.;

            for (int i = 0 ; i < Nobservation_ ; i++)
            {
                const double PI = 3.1415926535897;
                delta_(i) = 0.5 *
                    (1. + std::cos(PI * level_set(i) / epsilon_tilde));

                if (fabs(level_set(i)) > epsilon_tilde)
                    delta_(i) = 0.;
            }

            Mlt(1. / epsilon_tilde, delta_);
        }
        else
        {
            for (int i = 0 ; i < Nobservation_ ; i++)
            {
                const double PI = 3.1415926535897;
                delta_(i) = epsilon / (std::sqrt(PI) *
                                       (level_set(i) *
                                        level_set(i) + epsilon * epsilon));
            }
        }
    }


    template <class Model>
    template <class state>
    double LevelSetObservationManager<Model>::GetMeanValueOnDomain(
        const state& x,
        const observation& observation)
    {
        double domain, domain_size = 0;
        double threshold_activation = GetThresholdActivation();
        double mean_value = 0;

        for (int i = 0 ; i < Nobservation_ ; ++i)
        {
            if (x(i) > threshold_activation)
                domain = 1;
            else
                domain = 0;

            mean_value += domain * observation(i);
            domain_size += domain;
        }

        domain_size = max(domain_size, 1e-16); // to avoid getting 0.

        mean_value /= domain_size;

        return mean_value;
    }


    template <class Model>
    template <class state>
    double LevelSetObservationManager<Model>
    ::GetMeanValueOnComplementaryDomain(
        const state& x,
        const observation& observation)
    {
        double domain, domain_size = 0;
        double threshold_activation = GetThresholdActivation();
        double mean_value = 0;

        for (int i = 0 ; i < Nobservation_ ; ++i)
        {
            if (x(i) > threshold_activation)
                domain = 0;
            else
                domain = 1;

            mean_value += domain * observation(i);
            domain_size += domain;
        }

        domain_size = max(domain_size, 1e-16); // to avoid getting 0.

        mean_value /= domain_size;

        return mean_value;
    }


    ////////////
    // ACCESS //
    ////////////


    //! Returns the nudging matrix.
    /*!
      \param[in] x a vector representing the state.
      \param[out] output the nudging matrix at \a x.
    */
    template <class Model>
    template <class state, class mat>
    void LevelSetObservationManager<Model>
    ::GetNudgingMatrix(const state& x, mat& M) const
    {
        int state_size = x.GetSize();
        M.Reallocate(Nobservation_, state_size);
        M.SetIdentity();
    }


    //! Indicates if some observations are available at current time.
    template <class Model>
    bool LevelSetObservationManager<Model>::HasObservation() const
    {
        return available_time_.GetSize() != 0;
    }


    //! Returns the number of available observations.
    /*!
      \return The total number of observation at current time.
    */
    template <class Model>
    int LevelSetObservationManager<Model>::GetNobservation() const
    {
        return Nobservation_;
    }


    ///////////////
    // OPERATORS //
    ///////////////


    //! Applies the observation operator to a given vector.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] x a vector.
      \param[out] y the value of the operator applied to \a x. It is resized
      if needed.
    */
    template <class Model>
    template <class state>
    void LevelSetObservationManager<Model>
    ::ApplyOperator(const state& x, observation& y) const
    {
        double max_potential = x.GetNormInf();

        double threshold_activation = GetThresholdActivation();

        if (CaseIn1D())
        {
            if (max_potential > threshold_activation)
            {
                const observation& threshold_activation_vector =
                    GetThresholdActivationVector();

                state sorted_state(x);

                Add(-1.0, threshold_activation_vector, sorted_state);

                for (int i = 0 ; i < Nstate_model_ ; ++i)
                {
                    sorted_state(i) = fabs(sorted_state(i));
                }

                state indices(Nobservation_);
                indices.Zero();
                indices.Fill(); // Needed for the sort to work.

                Sort(sorted_state, indices);

                double x0 = indices(0);
                int i_count = 1;
                double x1 = indices(i_count);

                double model_radius = GetStimulationRadius1D();
                double epsilon = GetEpsilon();
                double delta_x_1D = GetDeltaX1D();

                while (fabs(indices(i_count) - x0) < model_radius)
                {
                    x1 = indices(i_count);
                    ++i_count;
                }

                for (int i = 0 ; i < Nobservation_ ; ++i)
                {
                    vector<double> v;
                    v.push_back(fabs(i - x0));
                    v.push_back(fabs(i - x1));
                    vector<int> swap_indices;
                    swap_indices.push_back(0);
                    swap_indices.push_back(1);
                    vector<double> distances = v;

                    if (v[0] > v[1])
                    {
                        swap_indices[0] = 1;
                        swap_indices[1] = 0;

                        distances[0] =
                            v[static_cast<unsigned int>(swap_indices[0])];
                        distances[1] =
                            v[static_cast<unsigned int>(swap_indices[1])];
                    }

                    double distance = distances[0];

                    double xg;

                    if (swap_indices[0] == 0)
                        xg = x0;
                    else
                        xg = x1;

                    double distance_compare = epsilon * fabs(x0 - x1);

                    if (x(i) > threshold_activation)
                        if (distance >= distance_compare)
                            y(i) = 1;
                        else
                            y(i) = fabs(i - xg) / distance_compare;
                    else if (x(i) < threshold_activation)
                        if (distance >= epsilon *
                            (time_ / delta_x_1D - fabs(x0 - x1)))
                            y(i) = -1;
                        else
                            y(i) = -fabs(i - xg) /
                                (epsilon *
                                 (time_ / delta_x_1D - fabs(x0 - x1)));
                    else
                        y(i) = 0;
                }
            }
        }
        else
            if (max_potential > threshold_activation)
                for (int i = 0 ; i < Nobservation_ ; ++i)
                {
                    if (x(i) > threshold_activation)
                        y(i) = 1;
                    else if (x(i) < threshold_activation)
                        y(i) = -1;
                    else
                        y(i) = 0;
                }
    }


    //! Applies the tangent linear operator to a given vector.
    /*! This method is called after 'SetTime' set the time at which the
      operator is defined.
      \param[in] x a vector.
      \param[out] y the value of the tangent linear operator applied to \a
      x. It is resized if needed.
    */
    template <class Model>
    template <class state>
    void LevelSetObservationManager<Model>
    ::ApplyTangentLinearOperator(const state& x, observation& y) const
    {
        throw ErrorUndefined(
            "void LevelSetObservationManager<Model>::"
            "ApplyTangentLinearOperator"
            "(const state& x, observation& y) const");
    }


    //! Applies the adjoint operator to a given vector.
    /*!
      \param[in] x a vector.
      \param[out] y the value of the operator at \a x. It is resized if
      needed.
    */
    template <class Model>
    template <class state>
    void LevelSetObservationManager<Model>
    ::ApplyAdjointOperator(const state& x, observation& y) const
    {
        throw ErrorUndefined(
            "void LevelSetObservationManager<Model>"
            "::ApplyAdjointOperator(const state& x,"
            "observation& y) const");
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model>
    string LevelSetObservationManager<Model>::GetName() const
    {
        return "LevelSetObservationManager";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model>
    void LevelSetObservationManager<Model>
    ::Message(string message)
    {
        // Put here any processing you need.
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation_vector2&
    LevelSetObservationManager<Model>::GetVectorObservation2()
    {
        return observation2_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation_vector2&
    LevelSetObservationManager<Model>::GetVectorInnovation2()
    {
        return innovation2_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetVectorInnovation()
    {
        return innovation_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::GetThresholdActivation() const
    {
        return threshold_activation_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::GetStimulationRadius1D() const
    {
        return stimulation_radius_1D_;
    }


    template <class Model>
    bool LevelSetObservationManager<Model>::CaseIn1D() const
    {
        return case_in_1D_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::GetDeltaX1D() const
    {
        return delta_x_1D_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::GetEpsilon() const
    {
        return epsilon_;
    }


    template <class Model>
    bool LevelSetObservationManager<Model>::WithTopologicalGradient() const
    {
        return with_topological_gradient_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetLevelSet()
    {
        return level_set_;
    }


    template <class Model>
    const typename  LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetThresholdActivationVector() const
    {
        return threshold_activation_vector_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetSignMult()
    {
        return sign_mult_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetDelta()
    {
        return delta_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::observation&
    LevelSetObservationManager<Model>::GetTopologicalGradient()
    {
        return topological_gradient_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::
    GetTopologicalGradientCoefficient() const
    {
        return topological_gradient_coeff_;
    }


    template <class Model>
    inline const Model& LevelSetObservationManager<Model>::GetModel() const
    {
        assert(!(!model_));
        return *model_;
    }

    template <class Model>
    inline Model& LevelSetObservationManager<Model>::GetNonCstModel()
    {
        return const_cast<Model&>(LevelSetObservationManager<Model>
                                  ::GetModel());
    }


    template <class Model>
    bool LevelSetObservationManager<Model>::InterpolateObservations() const
    {
        return interpolate_observations_;
    }


    template <class Model>
    typename LevelSetObservationManager<Model>::time_vector&
    LevelSetObservationManager<Model>::GetAvailableTime()
    {
        return available_time_;
    }


    template <class Model>
    double LevelSetObservationManager<Model>::GetTime() const
    {
        return time_;
    }

}

#endif
