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


#ifndef VERDANDI_FILE_EXTENDEDKALMANFILTER_CXX


#include "ExtendedKalmanFilter.hxx"

#include "BLUE.cxx"

namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Main constructor.
    /*! Builds the driver.
     */
    template <class Model, class ObservationManager>
    ExtendedKalmanFilter<Model, ObservationManager>
    ::ExtendedKalmanFilter(): iteration_(-1)
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ExtendedKalmanFilter::StaticMessage);
    }


    //! Destructor.
    template <class Model, class ObservationManager>
    ExtendedKalmanFilter<Model, ObservationManager>
    ::~ExtendedKalmanFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////



    //! Initializes the extended Kalman filter driver.
    /*! Initializes the model and the observation manager. */
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);
    }


    //! Initializes the extended Kalman filter driver.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("extended_kalman_filter.");


        /*********************************
         * Model and observation manager *
         *********************************/


        iteration_ = 0;

        configuration.Set("model.configuration_file", "",
                          configuration_file_,
                          model_configuration_file_);
        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        Nstate_ = model_.GetNstate();

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file_,
                          observation_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);

        Nobservation_  = observation_manager_.GetNobservation();

        Copy(model_.GetStateErrorVariance(), state_error_variance_);


        /***************************
         * Reads the configuration *
         ***************************/


        /*** Display options ***/

        // Should iterations be displayed on screen?
        configuration.Set("display.iteration", option_display_["iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("display.time", option_display_["time"]);
        // Should the analysis times be displayed on screen?
        configuration.Set("display.analysis_time",
                          option_display_["analysis_time"]);

        /*** Assimilation options ***/

        configuration.Set("data_assimilation.analyze_first_step",
                          analyze_first_step_);

        configuration.Set("BLUE_computation",
                          "ops_in(v, {'vector', 'matrix'})",
                          blue_computation_);
        configuration.Set("covariance_computation",
                          "ops_in(v, {'vector', 'matrix'})",
                          covariance_computation_);

        /*** Ouput saver ***/

        configuration.SetPrefix("extended_kalman_filter.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("forecast_time");
        output_saver_.Empty("forecast_state");
        output_saver_.Empty("analysis_time");
        output_saver_.Empty("analysis_state");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("extended_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        if (option_display_["iteration"])
            Logger::StdOut(*this, "Initialization");
        else
            Logger::Log<-3>(*this, "Initialization");
        if (option_display_["time"])
            Logger::StdOut(*this, "Initial time: "
                           + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Initial time: " + to_str(model_.GetTime()));

        /*** Assimilation ***/

        // To save the initial condition before assimilation.
        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "driver", "initial condition");

        if (analyze_first_step_)
            Analyze();

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the model.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

        if (option_display_["iteration"])
            Logger::StdOut(*this, "Starting iteration "
                           + to_str(iteration_)
                           + " -> " + to_str(iteration_ + 1));
        else
            Logger::Log<-3>(*this, "Starting iteration "
                            + to_str(iteration_)
                            + " -> " + to_str(iteration_ + 1));
        if (option_display_["time"])
            Logger::StdOut(*this, "Starting iteration at time "
                           + to_str(model_.GetTime()));
        else
            Logger::Log<-3>(*this,
                            "Starting iteration at time "
                            + to_str(model_.GetTime()));

        model_.InitializeStep();

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward, with extended Kalman filter.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::Forward()
    {
        Prediction();
        Analyze();
    }


    //! Performs a forecast step.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::Prediction()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        time_ = model_.GetTime();

        model_.Forward();

        PropagateCovarianceMatrix();

        ++iteration_;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    /*! Whenever observations are available, it computes BLUE.
     */
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::Analyze()
    {

        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            if (option_display_["analysis_time"])
                Logger::StdOut(*this, "Computing an analysis at time "
                               + to_str(model_.GetTime()));
            else
                Logger::Log<-3>(*this,"Computing an analysis at time "
                                + to_str(model_.GetTime()));

            Nstate_ = model_.GetNstate();
            model_state& x = model_.GetState();
            observation& innovation = observation_manager_.GetInnovation(x);
            Nobservation_ = innovation.GetSize();

            ComputeBLUE(innovation, x);

            model_.StateUpdated();

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Finalizes a step for the model.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Computes covariance.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::PropagateCovarianceMatrix_vector()
    {
        double saved_time = model_.GetTime();
        model_.SetTime(time_);

        // One column of covariance matrix P.
        model_state_error_variance_row error_covariance_column(Nstate_);

        for (int j = 0; j < Nstate_; j++)
        {
            GetCol(state_error_variance_, j,
                   error_covariance_column);
            model_.ApplyTangentLinearOperator(error_covariance_column);
            SetCol(error_covariance_column, j, state_error_variance_);
        }

        Transpose(state_error_variance_);

        for (int j = 0; j < Nstate_; j++)
        {
            GetCol(state_error_variance_, j, error_covariance_column);
            model_.ApplyTangentLinearOperator(error_covariance_column);
            SetCol(error_covariance_column, j, state_error_variance_);
        }

        if (model_.GetErrorVariance().GetM() != 0
            && model_.GetErrorVariance().GetN() != 0)
            Add(Ts(1.), model_.GetErrorVariance(), state_error_variance_);

        model_.SetTime(saved_time);
    }


    //! Computes covariance.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::PropagateCovarianceMatrix_matrix()
    {
        double saved_time = model_.GetTime();
        model_.SetTime(time_);

        model_tangent_linear_operator& A = model_.GetTangentLinearOperator();

        MltAdd(Ts(1.), A, state_error_variance_, Ts(0.),
               state_error_variance_);

        MltAdd(Ts(1.), SeldonNoTrans, state_error_variance_,
               SeldonTrans, A, Ts(0.), state_error_variance_);

        if (model_.GetErrorVariance().GetM() != 0
            && model_.GetErrorVariance().GetN() != 0)
            Add(Ts(1.), model_.GetErrorVariance(), state_error_variance_);

        model_.SetTime(saved_time);
    }


    //! Computes BLUE for Extended Kalman Filter.
    /*! The state is updated by the combination of background state and
      innovation. It computes the BLUE (best linear unbiased estimator).
      \param[in] state_vector the state vector to analyze.
    */
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::ComputeBLUE(const observation& innovation, model_state& state)
    {
        if (blue_computation_ == "vector")
            throw ErrorUndefined("ExtendedKalmanFilter"
                                 "::ComputeBLUE()");
        else
            ComputeBLUE_matrix(state_error_variance_,
                               observation_manager_.
                               GetTangentLinearOperator(),
                               innovation,
                               observation_manager_.GetErrorVariance(),
                               state, true, true);
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class Model, class ObservationManager>
    bool ExtendedKalmanFilter<Model, ObservationManager>::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Computes Covariance.
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::PropagateCovarianceMatrix()
    {
        if (covariance_computation_ == "vector")
            PropagateCovarianceMatrix_vector();
        else
            PropagateCovarianceMatrix_matrix();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class Model, class ObservationManager>
    Model&
    ExtendedKalmanFilter<Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class Model, class ObservationManager>
    ObservationManager&
    ExtendedKalmanFilter<Model, ObservationManager>
    ::GetObservationManager()
    {
        return observation_manager_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class Model, class ObservationManager>
    OutputSaver&
    ExtendedKalmanFilter<Model, ObservationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model, class ObservationManager>
    string ExtendedKalmanFilter<Model, ObservationManager>::GetName() const
    {
        return "ExtendedKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model, class ObservationManager>
    void ExtendedKalmanFilter<Model, ObservationManager>
    ::Message(string message)
    {
        if (message.find("initial condition") != string::npos)
        {
            output_saver_.Save(model_.GetTime(), model_.GetTime(),
                               "forecast_time");
            output_saver_.Save(model_.GetState(), model_.GetTime(),
                               "forecast_state");
        }
        if (message.find("forecast") != string::npos)
        {
            output_saver_.Save(model_.GetTime(), model_.GetTime(),
                               "forecast_time");
            output_saver_.Save(model_.GetState(), model_.GetTime(),
                               "forecast_state");
        }
        if (message.find("analysis") != string::npos)
        {
            output_saver_.Save(model_.GetTime(), model_.GetTime(),
                               "analysis_time");
            output_saver_.Save(model_.GetState(), model_.GetTime(),
                               "analysis_state");
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_EXTENDEDKALMANFILTER_CXX
#endif
