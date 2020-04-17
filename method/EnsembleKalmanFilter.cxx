// Copyright (C) 2011 INRIA
// Author(s): Kévin Charpentier, Vivien Mallet
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


#ifndef VERDANDI_FILE_METHOD_ENSEMBLEKALMANFILTER_CXX

#include "EnsembleKalmanFilter.hxx"

namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver.
     */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    EnsembleKalmanFilter<Model, ObservationManager, PerturbationManager>
    ::EnsembleKalmanFilter(): iteration_(-1)
    {

        /*** Initializations ***/

#if defined(VERDANDI_WITH_MPI)
        int initialized;
        MPI_Initialized(&initialized);
        if (!initialized)
            MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &Nprocess_);
#endif

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("observation_manager",
                                     observation_manager_,
                                     ObservationManager::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     EnsembleKalmanFilter::StaticMessage);
    }


    //! Destructor.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    EnsembleKalmanFilter<Model, ObservationManager, PerturbationManager>
    ::~EnsembleKalmanFilter()
    {
#if defined(VERDANDI_WITH_MPI)
        int finalized;
        MPI_Finalized(&finalized);
        if (!finalized)
            MPI_Finalize();
#endif
    }


    /////////////
    // METHODS //
    /////////////


    //! Initializes the ensemble Kalman filter.
    /*! It reads the configuration and initializes the model and the
      observation manager. It can also compute an analysis with the model's
      initial condition.
      \param[in] configuration_file configuration file for the method.
      \param[in] initialize_model should the model be initialized with a call
      to Model::Initialize(string)?
      \param[in] initialize_observation_manager should the observation manager
      be initialized with a call to ObservationManager::Initialize(Model&,
      string)?
      \param[in] initialize_perturbation_manager should the perturbation
      manager be initialized with a call to
      PerturbationManager::Initialize(string)?
      \warning If \a initialize_model is set to false, the model should be
      initialized before calling this function.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>
    ::Initialize(string configuration_file,
                 bool initialize_model, bool initialize_observation_manager,
                 bool initialize_perturbation_manager)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration,
                   initialize_model, initialize_observation_manager,
                   initialize_perturbation_manager);

    }


    //! Initializes the ensemble Kalman filter.
    /*! It reads the configuration and initializes the model and the
      observation manager. It can also compute an analysis with the model's
      initial condition.
      \param[in] configuration configuration for the method.
      \param[in] initialize_model should the model be initialized with a call
      to Model::Initialize(string)?
      \param[in] initialize_observation_manager should the observation manager
      be initialized with a call to ObservationManager::Initialize(Model&,
      string)?
      \param[in] initialize_perturbation_manager should the perturbation
      manager be initialized with a call to
      PerturbationManager::Initialize(string)?
      \warning If \a initialize_model is set to false, the model should be
      initialized before calling this function.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager, PerturbationManager>
    ::Initialize(VerdandiOps& configuration,
                 bool initialize_model, bool initialize_observation_manager,
                 bool initialize_perturbation_manager)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("ensemble_kalman_filter.");

        iteration_ = 0;

        /*******************************************************
         * Model, perturbation manager and observation manager *
         *******************************************************/


        configuration.Set("model.configuration_file", "",
                          configuration_file_, model_configuration_file_);

        configuration.Set("perturbation_manager.configuration_file", "",
                          configuration_file_,
                          perturbation_manager_configuration_file_);

        configuration.Set("observation_manager.configuration_file", "",
                          configuration_file_,
                          observation_configuration_file_);


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

        /*** Ensemble data ***/

        configuration.Set("Nmember", Nmember_);
        Nlocal_member_ = Nmember_;
        int global_member_number = 0;

#if defined(VERDANDI_WITH_MPI)
        Nlocal_member_ = Nmember_ / Nprocess_;
        if (rank_ < Nmember_ % Nprocess_)
            Nlocal_member_++;
        int div = Nmember_ % Nprocess_;
        if (rank_ < div)
            global_member_number = rank_ * Nlocal_member_;
        else
            global_member_number = div * (Nlocal_member_ + 1)
                + (rank_ - div) * Nlocal_member_;
#endif

        /*** Ouput saver ***/

        configuration.SetPrefix("ensemble_kalman_filter.output_saver.");
        output_saver_.Initialize(configuration);

        output_saver_.Empty("forecast_time");
        output_saver_.Empty("forecast_state");
        output_saver_.Empty("analysis_time");
        output_saver_.Empty("analysis_state");

        for (int k = 0; k < Nlocal_member_; k++)
        {
            output_saver_.Empty("forecast_state-"
                                + to_str(k + global_member_number));
            output_saver_.Empty("analysis_state-"
                                + to_str(k + global_member_number));
        }

        /*** Logger and read configuration ***/

        configuration.SetPrefix("ensemble_kalman_filter.");

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration",
                              output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

        /*** Initialization of model, perturbation manager and observation
             manager ***/

        if (initialize_model)
            model_.Initialize(model_configuration_file_);

        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "driver", "initial condition");

        if (initialize_perturbation_manager)
            perturbation_manager_
                .Initialize(perturbation_manager_configuration_file_);

        if (initialize_observation_manager)
            observation_manager_.Initialize(model_,
                                            observation_configuration_file_);

        configuration.SetPrefix("ensemble_kalman_filter.");

        configuration.Set("observation_tangent_linear_operator_access",
                          "ops_in(v, {'element', 'matrix'})",
                          observation_tangent_linear_operator_access_);

#ifdef VERDANDI_WITH_MPI
        if (world_rank_ == 0)
        {
#endif
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
#ifdef VERDANDI_WITH_MPI
        }
#endif

        /*** Ensemble initialization ***/

        model_state state;
        Nstate_ = model_.GetNstate();
        Nfull_state_ = model_.GetNfull_state();
        Nparameter_ = model_.GetNparameter();
        ensemble_.resize(Nlocal_member_);
        parameter_.resize(Nparameter_);

        for (int p = 0; p < Nparameter_; p++)
            parameter_[p].resize(Nlocal_member_);

        for (int m = 0; m < Nlocal_member_; m++)
            ensemble_[m] = model_.GetFullState();

        InitializeEnsemble();

        /*** Assimilation ***/

        if (analyze_first_step_)
            Analyze();

        perturbation_manager_.Finalize();

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initialization of the perturbations.
    /*! \brief The perturbations of the parameters are generated independently
      for each member.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::InitializeEnsemble()
    {
        uncertain_parameter reference_parameter;

        for (int i = 0; i < Nparameter_; i++)
        {
            if (model_.GetParameterOption(i) == "init_step")
            {
                // Parameter to be perturbed.
                reference_parameter = model_.GetParameter(i);

                for (int m = 0; m < Nlocal_member_ ; m++)
                {
                    uncertain_parameter sample;
                    SetDimension(model_.GetParameter(i), sample);
                    Fill(sample, model_.GetParameterPDF(i));

                    if (model_.GetParameterPDF(i) == "Normal"
                        || model_.GetParameterPDF(i) == "LogNormal"
                        || model_.GetParameterPDF(i) == "BlockNormal"
                        || model_.GetParameterPDF(i) == "BlockLogNormal")
                        perturbation_manager_
                            .Sample(model_.GetParameterPDF(i),
                                    model_.GetParameterVariance(i),
                                    model_.GetParameterPDFData(i),
                                    model_.GetParameterCorrelation(i),
                                    sample);
                    else if (model_.GetParameterPDF(i) == "NormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "LogNormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "BlockNormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "BlockLogNormalHomogeneous")
                        perturbation_manager_
                            .Sample(model_.GetParameterPDF(i),
                                    model_.GetParameterVariance(i)(0, 0),
                                    model_.GetParameterPDFData(i),
                                    model_.GetParameterCorrelation(i),
                                    sample);

                    if (model_.GetParameterPDF(i) == "Normal"
                        || model_.GetParameterPDF(i) == "BlockNormal"
                        || model_.GetParameterPDF(i) == "NormalHomogeneous"
                        || model_.GetParameterPDF(i) ==
                        "BlockNormalHomogeneous")
                        Add(1., sample, model_.GetParameter(i));
                    else if (model_.GetParameterPDF(0) == "LogNormal"
                             || model_.GetParameterPDF(0) ==
                             "BlockLogNormal"
                             || model_.GetParameterPDF(0) ==
                             "LogNormalHomogeneous"
                             || model_.GetParameterPDF(0) ==
                             "BlockLogNormalHomogeneous")
                        for (int j = 0; j < sample.GetM(); j++)
                            model_.GetParameter(i)(j) *= sample(j);

                    parameter_[i][m] = model_.GetParameter(i);
                    // Puts back the reference parameter into the model.
                    model_.GetParameter(i) = reference_parameter;
                    model_.ParameterUpdated(i);
                }

                MessageHandler::Send(*this, "model", "perturbation");
            }
        }
    }


    //! Initializes a step for the method.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::InitializeStep()
    {
        MessageHandler::Send(*this, "all", "::InitializeStep begin");

#ifdef VERDANDI_WITH_MPI
        if (world_rank_ == 0)
        {
#endif
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
#ifdef VERDANDI_WITH_MPI
        }
#endif

        model_.InitializeStep();

        uncertain_parameter reference_parameter;

        for (int i = 0; i < Nparameter_; i++)
        {
            if (model_.GetParameterOption(i) == "every_step")
            {
                // Parameter to be perturbed.
                reference_parameter = model_.GetParameter(i);

                for (int m = 0; m < Nlocal_member_ ; m++)
                {
                    uncertain_parameter sample;
                    SetDimension(model_.GetParameter(i), sample);
                    Fill(sample, model_.GetParameterPDF(i));

                    if (model_.GetParameterPDF(i) == "Normal"
                        || model_.GetParameterPDF(i) == "LogNormal"
                        || model_.GetParameterPDF(i) == "BlockNormal"
                        || model_.GetParameterPDF(i) == "BlockLogNormal")
                        perturbation_manager_
                            .Sample(model_.GetParameterPDF(i),
                                    model_.GetParameterVariance(i),
                                    model_.GetParameterPDFData(i),
                                    model_.GetParameterCorrelation(i),
                                    sample);
                    else if (model_.GetParameterPDF(i) == "NormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "LogNormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "BlockNormalHomogeneous"
                             || model_.GetParameterPDF(i) ==
                             "BlockLogNormalHomogeneous")
                        perturbation_manager_
                            .Sample(model_.GetParameterPDF(i),
                                    model_.GetParameterVariance(i)(0, 0),
                                    model_.GetParameterPDFData(i),
                                    model_.GetParameterCorrelation(i),
                                    sample);
                    if (model_.GetParameterPDF(i) == "Normal"
                        || model_.GetParameterPDF(i) == "BlockNormal"
                        || model_.GetParameterPDF(i) == "NormalHomogeneous"
                        || model_.GetParameterPDF(i) ==
                        "BlockNormalHomogeneous")
                        Add(1., sample, model_.GetParameter(i));
                    else if (model_.GetParameterPDF(0) == "LogNormal"
                             || model_.GetParameterPDF(0) ==
                             "BlockLogNormal"
                             || model_.GetParameterPDF(0) ==
                             "LogNormalHomogeneous"
                             || model_.GetParameterPDF(0) ==
                             "BlockLogNormalHomogeneous")
                        for (int j = 0; j < sample.GetM(); j++)
                            model_.GetParameter(i)(j) *= sample(j);

                    parameter_[i][m] = model_.GetParameter(i);
                    // Puts back the reference parameter into the model.
                    model_.GetParameter(i) = reference_parameter;
                    model_.ParameterUpdated(i);
                }

                MessageHandler::Send(*this, "model", "perturbation");
            }
        }

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward for the ensemble.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::Forward()
    {
        Prediction();
        Analyze();
    }


    //! Performs a forecast step for the ensemble.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::Prediction()
    {
        MessageHandler::Send(*this, "all", "::Forward begin");

        time_ = model_.GetTime();
        vector<uncertain_parameter> reference_parameter;
        for (int i = 0; i < Nparameter_; i++)
            reference_parameter.push_back(model_.GetParameter(i));

        for (int m = 0; m < Nlocal_member_; m++)
        {
            for (int i = 0; i < Nparameter_; i++)
            {
                model_.GetParameter(i) = parameter_[i][m];
                model_.ParameterUpdated(i);
            }
            model_.GetFullState() = ensemble_[m];
            model_.FullStateUpdated();

            model_.Forward();

            ensemble_[m] = model_.GetFullState();

            if (m < Nlocal_member_ - 1)
                model_.SetTime(time_);
        }

        model_state mean_state_vector(model_.GetNstate());
        // Sets state to ensemble mean.
        mean_state_vector = 0.;
        for (int m = 0; m < Nlocal_member_; m++)
        {
            model_.GetFullState() = ensemble_[m];
            model_.FullStateUpdated();
            Add(Ts(1), model_.GetState(), mean_state_vector);
        }

#if defined(VERDANDI_WITH_MPI)
        model_state mean_state_global_vector(model_.GetNstate());
        MPI_Allreduce(mean_state_vector.GetData(), mean_state_global_vector.
                      GetData(), model_.GetNstate(), MPI_DOUBLE, MPI_SUM,
                      MPI_COMM_WORLD);
        mean_state_vector = mean_state_global_vector;
#endif

        Mlt(Ts(1) / Ts(Nmember_), mean_state_vector);
        model_.GetState() = mean_state_vector;
        model_.StateUpdated();

        // Puts back the reference parameters.
        for (int i = 0; i < Nparameter_; i++)
        {
            model_.GetParameter(i) = reference_parameter[i];
            model_.ParameterUpdated(i);
        }

        ++iteration_;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "observation_manager", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");
        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Computes an analysis.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::Analyze()
    {
        MessageHandler::Send(*this, "all", "::Analyze begin");

        observation_manager_.SetTime(model_, model_.GetTime());
        if (observation_manager_.HasObservation())
        {
#ifdef VERDANDI_WITH_MPI
            if (world_rank_ == 0)
            {
#endif
                if (option_display_["analysis_time"])
                    Logger::StdOut(*this, "Computing an analysis at time "
                                   + to_str(model_.GetTime()));
                else
                    Logger::Log<-3>(*this,"Computing an analysis at time "
                                    + to_str(model_.GetTime()));
#ifdef VERDANDI_WITH_MPI
            }
#endif

            observation& obs = observation_manager_.GetObservation();
            Nobservation_ = obs.GetLength();

            model_state mean_state_vector(model_.GetNstate());
            // Sets state to ensemble mean.
            mean_state_vector = model_.GetState();

            // Computes the innovation vectors d = y - Hx where H is the
            // observation operator.
            Matrix<To> innovation_matrix(Nobservation_, Nlocal_member_);
            observation Hx(Nobservation_);
            model_state state(Nstate_);
            for (int m = 0; m < Nlocal_member_; m++)
            {
                model_.GetFullState() = ensemble_[m];
                model_.FullStateUpdated();
                observation_manager_.ApplyOperator(model_.GetState(), Hx);
                Mlt(To(-1), Hx);
                Add(To(1), obs, Hx);
                SetCol(Hx, m, innovation_matrix);
            }

            // Constructs the square root of the empirical variance.
            Matrix<Ts> L(Nstate_, Nlocal_member_);
            for (int l = 0; l < Nlocal_member_; l++)
                for (int k = 0; k < Nstate_; k++)
                    L(k, l) = Ts(1) / sqrt(Ts(Nmember_ - 1))
                        * (ensemble_[l](k) - mean_state_vector(k));

            // Computes H times L.
            Matrix<To> HL(Nobservation_, Nlocal_member_);

            if (observation_tangent_linear_operator_access_ == "matrix")
                MltAdd(To(1), observation_manager_.GetTangentLinearOperator(),
                       L, To(0), HL);
            else // "element".
            {
                HL.Fill(To(0));
                for (int i = 0; i < Nobservation_; i++)
                    for (int j = 0; j < Nlocal_member_; j++)
                        for (int k = 0; k < Nstate_; k++)
                            HL(i, j) +=
                                observation_manager_.
                                GetTangentLinearOperator(i, k) *
                                L(k, j);
            }

            // Reads R.
            Matrix<To> working_matrix(Nobservation_,Nobservation_);
            working_matrix = observation_manager_.GetErrorVariance();

            // 'working_matrix' stores HLL'H' + R.
            MltAdd(To(1), SeldonNoTrans, HL, SeldonTrans, HL,
                   To(1), working_matrix);

            // Computes (HLL'H' + R)^{-1} d.
            Matrix<To> correction(Nobservation_, Nlocal_member_);
            GetInverse(working_matrix);
            MltAdd(To(1), working_matrix, innovation_matrix,
                   To(0), correction);

            // Computes LL'H'
            Matrix<Ts> LLH(Nstate_, Nobservation_);
            MltAdd(Ts(1), SeldonNoTrans, L, SeldonTrans, HL, Ts(0), LLH);

            // Computes LL'H' (HLL'H' + R)^{-1} d.
            Matrix<Ts> Kd(Nstate_, Nlocal_member_);
            MltAdd(Ts(1), LLH, correction, Ts(0), Kd);

            // Updates the ensemble A += K * d.
            for (int m = 0; m < Nlocal_member_; m++)
                for (int l = 0; l < Nstate_; l++)
                    ensemble_[m](l) += Kd(l, m);

            // Sets state to ensemble mean.
            mean_state_vector = 0.;
            for (int m = 0; m < Nlocal_member_; m++)
            {
                model_.GetFullState() = ensemble_[m];
                model_.FullStateUpdated();
                Add(Ts(1), model_.GetState(), mean_state_vector);
            }

#if defined(VERDANDI_WITH_MPI)
            model_state mean_state_global_vector(model_.GetNstate());
            MPI_Allreduce(mean_state_vector.GetData(),
                          mean_state_global_vector.GetData(),
                          model_.GetNstate(), MPI_DOUBLE, MPI_SUM,
                          MPI_COMM_WORLD);
            mean_state_vector = mean_state_global_vector;
#endif

            Mlt(Ts(1) / Ts(Nmember_), mean_state_vector);
            model_.GetState() = mean_state_vector;
            model_.StateUpdated();

            MessageHandler::Send(*this, "model", "analysis");
            MessageHandler::Send(*this, "observation_manager", "analysis");
            MessageHandler::Send(*this, "driver", "analysis");
        }

        MessageHandler::Send(*this, "all", "::Analyze end");
    }


    //! Finalizes a step for the model.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    bool EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::HasFinished()
        const
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    Model&
    EnsembleKalmanFilter<Model, ObservationManager,
                         PerturbationManager>::GetModel()
    {
        return model_;
    }


    //! Returns the observation manager.
    /*!
      \return The observation manager.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    ObservationManager&
    EnsembleKalmanFilter<Model, ObservationManager,
                         PerturbationManager>::GetObservationManager()
    {
        return observation_manager_;
    }


    //! Returns the perturbation manager.
    /*!
      \return The perturbation manager.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    PerturbationManager&
    EnsembleKalmanFilter<Model, ObservationManager,
                         PerturbationManager>::GetPerturbationManager()
    {
        return perturbation_manager_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    OutputSaver& EnsembleKalmanFilter<Model, ObservationManager,
                                      PerturbationManager>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    string EnsembleKalmanFilter<Model, ObservationManager,
                                PerturbationManager>::GetName() const
    {
        return "EnsembleKalmanFilter";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>::Message(string message)
    {

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            if (message.find("initial condition") != string::npos)
            {
                output_saver_.Save(model_.GetTime(), model_.GetTime(),
                                   "forecast_time");
                output_saver_.Save(model_.GetState(), model_.GetTime(),
                                   "forecast_state");
            }
#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            if (message.find("forecast") != string::npos)
            {
                output_saver_.Save(model_.GetTime(), model_.GetTime(),
                                   "forecast_time");
                output_saver_.Save(model_.GetState(), model_.GetTime(),
                                   "forecast_state");
            }
        int global_member_number = 0;

#if defined(VERDANDI_WITH_MPI)
        int div = Nmember_ % Nprocess_;
        if (rank_ < div)
            global_member_number = rank_ * Nlocal_member_;
        else
            global_member_number = div * (Nlocal_member_ + 1)
                + (rank_ - div) * Nlocal_member_;
#endif

        if (message.find("forecast") != string::npos)
        {
            model_state mean_state_vector(model_.GetNfull_state());
            mean_state_vector = model_.GetFullState();

            for (int m = 0; m < Nlocal_member_; m++)
                if (output_saver_
                    .IsVariable("forecast_state-"
                                + to_str(m + global_member_number)))
                {
                    model_.GetFullState() = ensemble_[m];
                    model_.FullStateUpdated();
                    output_saver_.Save(model_.GetState(),
                                       model_.GetTime(),
                                       "forecast_state-"
                                       + to_str(m + global_member_number));
                }

            model_.GetFullState() = mean_state_vector;
            model_.FullStateUpdated();
        }

#if defined(VERDANDI_WITH_MPI)
        if (rank_ == 0)
#endif
            if (message.find("analysis") != string::npos)
            {
                output_saver_.Save(model_.GetTime(), model_.GetTime(),
                                   "analysis_time");
                output_saver_.Save(model_.GetState(), model_.GetTime(),
                                   "analysis_state");
            }

        if (message.find("analysis") != string::npos)
        {
            model_state mean_state_vector(model_.GetNfull_state());
            mean_state_vector = model_.GetFullState();

            for (int m = 0; m < Nlocal_member_; m++)
                if (output_saver_
                    .IsVariable("analysis_state-"
                                + to_str(m + global_member_number)))
                {
                    model_.GetFullState() = ensemble_[m];
                    model_.FullStateUpdated();
                    output_saver_.Save(model_.GetState(),
                                       model_.GetTime(),
                                       "analysis_state-" +
                                       to_str(m + global_member_number));
                }
            model_.GetFullState() = mean_state_vector;
            model_.FullStateUpdated();
        }
    }


    ///////////////////////
    // PROTECTED METHODS //
    ///////////////////////


    /*! \brief Fills an input vector collection according to its probability
      distribution. */
    /*!
      \param[in,out] in input collection vector.
      \param[in] pdf probability density function: Normal, NormalHomogeneous,
      LogNormal or LogNormalHomogeneous.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    template <class T0, class Allocator0>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>
    ::Fill(Vector<T0, Collection, Allocator0>& in,
           string pdf)
    {
        if (pdf == "Normal" || pdf == "NormalHomogeneous")
            for (int i = 0; i < in.GetNvector(); i++)
                in.GetVector(i).Fill(typename T0::value_type(0));
        else if (pdf == "LogNormal" || pdf == "LogNormalHomogeneous")
            for (int i = 0; i < in.GetNvector(); i++)
                in.GetVector(i).Fill(typename T0::value_type(1));
    }


    /*! \brief Fills an input vector according to its probability
      distribution. */
    /*!
      \param[in,out] in input vector.
      \param[in] pdf probability density function: Normal, NormalHomogeneous,
      LogNormal or LogNormalHomogeneous.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    template <class T0, class Storage0, class Allocator0>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>
    ::Fill(Vector<T0, Storage0, Allocator0>& in, string pdf)
    {
        if (pdf == "Normal" || pdf == "NormalHomogeneous")
            in.Fill(T0(0));
        else if (pdf == "LogNormal" || pdf == "LogNormalHomogeneous")
            in.Fill(T0(1));
    }


    //! Allocates an output vector to the dimension of the input vector.
    /*!
      \param[in] in input vector.
      \param[out] out output vector.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    template <class T0, class Storage0, class Allocator0>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>
    ::SetDimension(Vector<T0, Storage0, Allocator0>& in,
                   Vector<T0, Storage0, Allocator0>& out)
    {
        out.Reallocate(in.GetLength());
    }


    /*! \brief Allocates an output vector collection to the dimension of the
      input vector collection. */
    /*!
      \param[in] in input collection vector.
      \param[out] out output collection vector.
    */
    template <class Model, class ObservationManager,
              class PerturbationManager>
    template <class T0, class Allocator0>
    void EnsembleKalmanFilter<Model, ObservationManager,
                              PerturbationManager>
    ::SetDimension(Vector<T0, Collection, Allocator0>& in,
                   Vector<T0, Collection, Allocator0>& out)
    {
        T0 suboutput;
        for (int i = 0; i < in.GetNvector(); i++)
        {
            suboutput.Reallocate(in.GetVector(i).GetLength());
            out.AddVector(suboutput);
            suboutput.Nullify();
        }
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_ENSEMBLEKALMANFILTER_CXX
#endif
