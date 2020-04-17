// Copyright (C) 2014 INRIA
// Author(s): Claude Nicolas
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


#ifndef VERDANDI_FILE_METHOD_NUDGING_CXX


#include "Nudging.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    template <class Model, class ObservationManager>
    Nudging<Model, ObservationManager>::Nudging(): iteration_(-1)
    {
        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     Nudging::StaticMessage);
    }


    //! Destructor.
    template <class Model, class ObservationManager>
    Nudging<Model, ObservationManager>::~Nudging()
    {
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the nudging.
    /*! Initializes the model and the observation manager. */
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager);
    }


    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");
        configuration_file_ = configuration.GetFilePath();

        iteration_ = 0;
        previous_time_ = 0;

        /***************************
         * Reads the configuration *
         ***************************/

        configuration.SetPrefix("nudging.");

        /*** Model ***/

        configuration.Set("model.configuration_file", "",
                          configuration_file_,
                          model_configuration_file_);

        /*** Observation manager ***/

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file_,
                          observation_configuration_file_);

        /*** Display options ***/

        configuration.SetPrefix("nudging.display.");
        // Should iterations be displayed on screen?
        configuration.Set("iteration", option_display_["iteration"]);
        // Should current time be displayed on screen?
        configuration.Set("time", option_display_["time"]);
        // Should the analysis times be displayed on screen?
        configuration.Set("analysis_time", option_display_["analysis_time"]);

        /*** Assimilation options ***/

        configuration.
            SetPrefix("nudging.data_assimilation.");
        configuration.Set("analyze_first_step", analyze_first_step_);
        configuration.Set("nudging_type",
                          "ops_in(v, {'standard', 'dt', 'source', 'BNF'})",
                          nudging_type_);
        configuration.Set("matrix_fixed", matrix_fixed_);
        if (matrix_fixed_)
        {
            int m, n;
            configuration.Set("nudging_matrix.m", m);
            configuration.Set("nudging_matrix.n", n);
            nudging_matrix_.Reallocate(m, n);
            configuration.Set("nudging_matrix.matrix", nudging_matrix_);
        }

        if (analyze_first_step_ && nudging_type_ == "dt")
            throw  ErrorProcessing("Nudging::Initialize(VerdandiOps&, bool,"
                                   " bool)", "analyze_first_step can't be"
                                   "\"true\" when nudging_type = \"dt\".");

        configuration.Set("nudging_gain", nudging_gain_);

        /*** Ouput saver ***/

        configuration.SetPrefix("nudging.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("forecast_time");
        output_saver_.Empty("forecast_state");
        output_saver_.Empty("analysis_time");
        output_saver_.Empty("analysis_state");
        output_saver_.Empty("analysis_variance_diagonal");

        /*** Logger and read configuration ***/

        configuration.SetPrefix("nudging.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        /*** Initializations ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);
        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);

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


    //! Initializes a step for the nudging.
    /*! Initializes a step for the model.
     */
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::InitializeStep()
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


    //! Performs a step forward with nudging.
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::Prediction()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        model_.Forward();

        // If the nudging type is set to "dt", we compute the delta_t between
        // this time and the previous.
        if (nudging_type_ == "dt")
        {
            Ts new_time = model_.GetTime();
            delta_t_ = new_time - previous_time_;
            previous_time_ = new_time;
        }

        ++iteration_;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");
        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Performs a step forward then an analysis.
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::Forward()
    {
        Prediction();
        Analyze();
    }


    //! Computes an analysis.
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::Analyze()
    {
        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());

        if (observation_manager_.HasObservation())
        {
            model_state& state = model_.GetState();
            observation& innovation =
                observation_manager_.GetInnovation(state);
            if (!matrix_fixed_)
                observation_manager_.GetNudgingMatrix(state, nudging_matrix_);
            Ts gain = nudging_gain_;

            if (nudging_type_ == "source")
            {
                // In the "source" nudging, the discretization of the source
                // is done inside the model.
                model_state source;
                source.Resize(state.GetSize());
                MltAdd(gain, nudging_matrix_, innovation, Ts(0), source);
                model_.SetSource(source);
                return;
            }

            if (option_display_["analysis_time"])
                Logger::StdOut(*this, "Computing an analysis at time "
                               + to_str(model_.GetTime()));
            else
                Logger::Log<-3>(*this,"Computing an analysis at time "
                                + to_str(model_.GetTime()));


            // If the 'nudging_type_' is "dt", we consider that the time
            // discretization for the source is done here, hence the nudging
            // gain is multiplied by 'dt'.
            if (nudging_type_ == "dt")
                gain *= delta_t_;

            // The new state is a weighted average between the state and the
            // innovation. The weight is given by the 'nudging_matrix_' and
            // the 'nudging_gain_': "state = state + gain * matrix *
            // innovation".
            MltAdd(gain, nudging_matrix_, innovation, Ts(1), state);
            model_.StateUpdated();

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Finalizes a step for the model.
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>::Finalize()
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
    bool Nudging<Model, ObservationManager>::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class Model, class ObservationManager>
    Model& Nudging<Model, ObservationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class Model, class ObservationManager>
    OutputSaver& Nudging<Model, ObservationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model, class ObservationManager>
    string Nudging<Model, ObservationManager>::GetName() const
    {
        return "Nudging";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model, class ObservationManager>
    void Nudging<Model, ObservationManager>
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


#define VERDANDI_FILE_METHOD_NUDGING_CXX
#endif
