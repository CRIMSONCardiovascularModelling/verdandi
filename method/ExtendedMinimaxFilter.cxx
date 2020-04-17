// Copyright (C) 2013 INRIA
// Author(s): Vivien Mallet
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


#ifndef VERDANDI_FILE_EXTENDEDMINIMAXFILTER_CXX


#include "ExtendedMinimaxFilter.hxx"


namespace Verdandi
{


    /////////////////////////////////
    // CONSTRUCTORS AND DESTRUCTOR //
    /////////////////////////////////


    //! Main constructor.
    /*! Builds the driver.
     */
    template <class Model, class ObservationManager>
    ExtendedMinimaxFilter<Model, ObservationManager>
    ::ExtendedMinimaxFilter(): iteration_(-1)
    {

        /*** Initializations ***/

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ExtendedMinimaxFilter::StaticMessage);
    }


    //! Destructor.
    template <class Model, class ObservationManager>
    ExtendedMinimaxFilter<Model, ObservationManager>
    ::~ExtendedMinimaxFilter()
    {
    }


    /////////////
    // METHODS //
    /////////////



    //! Initializes the extended minimax filter.
    /*! Initializes the model and the observation manager. */
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);
    }


    //! Initializes the extended minimax filter.
    /*! Initializes the model and the observation manager. Optionally computes
      the analysis of the first step. */
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("extended_minimax_filter.");


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

        Copy(model_.GetStateErrorVariance(), G_);
        // TODO: Add conversion from variance to gain (tolerance region).
        GetInverse(G_);


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

        /*** Ouput saver ***/

        configuration.SetPrefix("extended_minimax_filter.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("forecast_time");
        output_saver_.Empty("forecast_state");
        output_saver_.Empty("forecast_gain");
        output_saver_.Empty("analysis_time");
        output_saver_.Empty("analysis_state");
        output_saver_.Empty("analysis_gain");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("extended_minimax_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        /*** Assimilation at the first step ***/

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

        // To save the initial condition before assimilation.
        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "driver", "initial condition");

        if (analyze_first_step_)
        {
            // Data assimilation per se.
            observation_manager_.SetTime(model_, model_.GetTime());
            if (observation_manager_.HasObservation())
            {
                Nobservation_ = observation_manager_.GetNobservation();

                R_inv_ = observation_manager_.GetErrorVariance();
                // TODO: Add conversion from variance to gain (tolerance
                // region).
                GetInverse(R_inv_);
            }
            UpdateGain();
            ComputeEstimator();

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes a step for the model.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::InitializeStep()
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


    //! Performs a step forward, with analysis and gain update.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::Forward()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        time_ = model_.GetTime();

        PropagateGain();
        model_.Forward();

        ++iteration_;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Propagates the gain matrix over one iteration.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::PropagateGain()
    {
        MessageHandler::Send(*this, "all", "::PropagateGain begin");

        // Operations are carried out on the inverse of the gain.
        GetInverse(G_);

        // Tangent linear model.
        model_tangent_linear_operator& M = model_.GetTangentLinearOperator();

        // Propagates the gain, without model error: computes $M G^{-1} M^T$.
        model_state_error_variance M_Ginv(Nstate_, Nstate_);
        MltAdd(Ts(1), M, G_, Ts(0), M_Ginv);
        MltAdd(Ts(1), SeldonNoTrans, M_Ginv, SeldonTrans, M, Ts(0), G_);

        // Adds the model error.
        // TODO: Add conversion from variance to gain (tolerance region).
        if (model_.GetErrorVariance().GetM() != 0
            && model_.GetErrorVariance().GetN() != 0)
            Add(Ts(1), model_.GetErrorVariance(), G_);

        // Back to the gain.
        GetInverse(G_);

        MessageHandler::Send(*this, "all", "::PropagateGain end");
    }


    //! Performs a step forward, with analysis and gain update.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::Analyze()
    {
        MessageHandler::Send(*this, "all", "::Analyze begin");

        // Data assimilation per se.
        observation_manager_.SetTime(model_, model_.GetTime());
        if (observation_manager_.HasObservation())
        {
            if (option_display_["analysis_time"])
                Logger::StdOut(*this, "Computing an analysis at time "
                               + to_str(model_.GetTime()));
            else
                Logger::Log<-3>(*this,"Computing an analysis at time "
                                + to_str(model_.GetTime()));

            Nobservation_ = observation_manager_.GetNobservation();

            R_inv_ = observation_manager_.GetErrorVariance();
            // TODO: Add conversion from variance to gain (tolerance region).
            GetInverse(R_inv_);

            UpdateGain();
            ComputeEstimator();

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Updates the gain matrix when observations are available.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::UpdateGain()
    {
        MessageHandler::Send(*this, "all", "::UpdateGain begin");

        if (!observation_manager_.HasObservation())
            return;

        const observation_tangent_linear_operator& H
            = observation_manager_.GetTangentLinearOperator();

        // Adds $H^T R^{-1} H$ to the gain.
        matrix_state_observation Ht_Rinv(Nstate_, Nobservation_);
        MltAdd(Ts(1), SeldonTrans, H, SeldonNoTrans, R_inv_, Ts(0), Ht_Rinv);
        MltAdd(Ts(1), Ht_Rinv, H, Ts(1), G_);

        MessageHandler::Send(*this, "all", "::UpdateGain end");
    }


    //! Computes the estimator for the current time step.
    /*!
      \warning This method must be called after 'UpdateGain', so that
      'Ht_Rinv_' has been already computed.
    */
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::ComputeEstimator()
    {
        MessageHandler::Send(*this, "all", "::ComputeEstimator begin");

        if (!observation_manager_.HasObservation())
            return;

        const observation_tangent_linear_operator& H
            = observation_manager_.GetTangentLinearOperator();

        // Innovation $(y - H x)$.
        observation innovation = observation_manager_.GetObservation();
        MltAdd(Ts(-1), H, model_.GetState(), Ts(1), innovation);

        // Adds $H^T R^{-1} H$ to the gain.
        matrix_state_observation Ht_Rinv(Nstate_, Nobservation_);
        MltAdd(Ts(1), SeldonTrans, H, SeldonNoTrans, R_inv_, Ts(0), Ht_Rinv);

        // $H^T R^{-1} (y - H x)$.
        model_state Ht_Rinv_innovation(Nstate_), correction(Nstate_);
        Mlt(Ht_Rinv, innovation, Ht_Rinv_innovation);

        // Gets $G^{-1}$.
        model_state_error_variance Ginv = G_;
        GetInverse(Ginv);

        // Adds the correction $G^{-1} H^T R^{-1} (y - H x)$ to the state
        // vector.
        MltAdd(Ts(1), Ginv, Ht_Rinv_innovation, Ts(1), model_.GetState());

        MessageHandler::Send(*this, "all", "::ComputeEstimator end");
    }


    //! Finalizes a step for the model.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class Model, class ObservationManager>
    bool ExtendedMinimaxFilter<Model, ObservationManager>::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class Model, class ObservationManager>
    Model&
    ExtendedMinimaxFilter<Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class Model, class ObservationManager>
    ObservationManager&
    ExtendedMinimaxFilter<Model, ObservationManager>
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
    ExtendedMinimaxFilter<Model, ObservationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model, class ObservationManager>
    string ExtendedMinimaxFilter<Model, ObservationManager>::GetName() const
    {
        return "ExtendedMinimaxFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model, class ObservationManager>
    void ExtendedMinimaxFilter<Model, ObservationManager>
    ::Message(string message)
    {
        if (message.find("initial condition") != string::npos
            || message.find("forecast") != string::npos)
        {
            output_saver_.Save(model_.GetTime(), model_.GetTime(),
                               "forecast_time");
            output_saver_.Save(model_.GetState(), model_.GetTime(),
                               "forecast_state");
            output_saver_.Save(G_, model_.GetTime(), "forecast_gain");
        }
        if (message.find("analysis") != string::npos)
        {
            output_saver_.Save(model_.GetTime(), model_.GetTime(),
                               "analysis_time");
            output_saver_.Save(model_.GetState(), model_.GetTime(),
                               "analysis_state");
            output_saver_.Save(G_, model_.GetTime(), "analysis_gain");
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_EXTENDEDMINIMAXFILTER_CXX
#endif
