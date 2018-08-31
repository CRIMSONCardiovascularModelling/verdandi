// Copyright (C) 2009 INRIA
// Author(s): Vivien Mallet, Claire Mouton
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


#ifndef VERDANDI_FILE_METHOD_FORWARDDRIVER_CXX


#include "ForwardDriver.hxx"


namespace Verdandi
{


    ////////////////////////////////
    // CONSTRUCTOR AND DESTRUCTOR //
    ////////////////////////////////


    //! Main constructor.
    /*! Builds the driver and reads option keys in the configuration file.
      \param[in] configuration configuration file.
    */
    template <class Model>
    ForwardDriver<Model>::ForwardDriver():
        iteration_(-1)
    {

        /*** Initializations ***/

#if defined(VERDANDI_WITH_MPI)
        int initialized;
        MPI_Initialized(&initialized);
        if (!initialized)
            MPI_Init(NULL, NULL);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &Nprocess_);
#endif

        MessageHandler::AddRecipient("model", model_, Model::StaticMessage);
        MessageHandler::AddRecipient("driver", *this,
                                     ForwardDriver::StaticMessage);
    }


    //! Destructor.
    template <class Model>
    ForwardDriver<Model>::~ForwardDriver()
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


    //! Initializes the simulation.
    /*! Initializes the model.
      \param[in] configuration_file configuration file to be given to the
      model initialization method.
    */
    template <class Model>
    void ForwardDriver<Model>::Initialize(string configuration_file,
                                          bool initialize_model)
    {
        VerdandiOps configuration(configuration_file);
        Initialize(configuration, initialize_model);
    }


    //! Initializes the simulation.
    /*! Initializes the model.
      \param[in] configuration configuration file to be given to the
      model initialization method.
    */
    template <class Model>
    void ForwardDriver<Model>::Initialize(VerdandiOps& configuration,
                                          bool initialize_model)
    {
        MessageHandler::Send(*this, "all", "::Initialize begin");

        configuration_file_ = configuration.GetFilePath();
        configuration.SetPrefix("forward.");

#ifdef VERDANDI_WITH_MPI
        configuration.Set("mpi_grid.Nrow", Nrow_);
        configuration.Set("mpi_grid.Ncol", Ncol_);

        if (Nprocess_ != Nrow_ * Ncol_)
            throw ErrorConfiguration("ForwardDriver<Model>::Initialize",
                                     "Wrong number of processes: "
                                     + to_str(Nprocess_)
                                     + ". The dimension of the MPI grid ("
                                     + to_str(Nrow_) + ", " + to_str(Ncol_) +
                                     ") requires " + to_str(Nrow_ * Ncol_)
                                     + " processes.");

        SetGridCommunicator(Nrow_, Ncol_, &row_communicator_,
                            &col_communicator_);
        MPI_Comm_rank(row_communicator_, &model_task_);
#endif

        /*********
         * Model *
         *********/


        configuration.Set("model.configuration_file", "",
                          configuration_file_,
                          model_configuration_file_);
        if (initialize_model)
        {
#ifdef VERDANDI_WITH_MPI
            model_.SetMPICommunicator(col_communicator_);
#endif
            model_.Initialize(model_configuration_file_);
        }

        /***************************
         * Reads the configuration *
         ***************************/


        /*** Display options ***/

        // Should the iteration be displayed on screen?
        configuration.Set("display.iteration", option_display_["iteration"]);
        // Should the time be displayed on screen?
        configuration.Set("display.time", option_display_["time"]);
#ifdef VERDANDI_WITH_MPI
        // Should the MPI grid be displayed on screen?
        configuration.Set("display.mpi_grid", option_display_["mpi_grid"]);
#endif

        iteration_ = 0;

        /*** Ouput saver ***/

        configuration.SetPrefix("forward.output_saver.");
        output_saver_.Initialize(configuration);
        output_saver_.Empty("forecast_time");
        output_saver_.Empty("forecast_state");
        configuration.SetPrefix("forward.");

        /*** Logger and read configuration ***/

        if (configuration.Exists("output.log"))
            Logger::SetFileName(configuration.Get<string>("output.log"));

        if (configuration.Exists("output.configuration"))
        {
            string output_configuration;
            configuration.Set("output.configuration", output_configuration);
            configuration.WriteLuaDefinition(output_configuration);
        }

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

#ifdef VERDANDI_WITH_MPI
        if (world_rank_ == 0)
        {
            if (option_display_["mpi_grid"])
                Logger::StdOut(*this, "world rank\tmodel task\tmodel rank");
            else
                Logger::Log<-3>(*this,
                                "world rank\tmodel task\tmodel rank");
        }

        MPI_Barrier(MPI_COMM_WORLD);
        int model_rank;
        MPI_Comm_rank(col_communicator_, &model_rank);
        if (option_display_["mpi_grid"])
            Logger::StdOut(*this, to_str(world_rank_) + "\t\t"
                           + to_str(model_task_) + "\t\t" +
                           to_str(model_rank));
        else
            Logger::Log<-3>(*this, to_str(world_rank_) + "\t\t"
                            + to_str(model_task_) + "\t\t" +
                            to_str(model_rank));
#endif

        MessageHandler::Send(*this, "model", "initial condition");
        MessageHandler::Send(*this, "driver", "initial condition");

        MessageHandler::Send(*this, "all", "::Initialize end");
    }


    //! Initializes the model before a time step.
    template <class Model>
    void ForwardDriver<Model>::InitializeStep()
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

        MessageHandler::Send(*this, "all", "::InitializeStep end");
    }


    //! Performs a step forward without optimal interpolation.
    template <class Model>
    void ForwardDriver<Model>::Forward()
    {
        time_.PushBack(model_.GetTime());

        MessageHandler::Send(*this, "all", "::Forward begin");

        model_.Forward();
        iteration_++;

        MessageHandler::Send(*this, "model", "forecast");
        MessageHandler::Send(*this, "driver", "forecast");

        MessageHandler::Send(*this, "all", "::Forward end");
    }


    //! Finalizes a step for the model.
    template <class Model>
    void ForwardDriver<Model>::FinalizeStep()
    {
        MessageHandler::Send(*this, "all", "::FinalizeStep begin");

        model_.FinalizeStep();

        MessageHandler::Send(*this, "all", "::FinalizeStep end");
    }


    //! Finalizes the model.
    template <class Model>
    void ForwardDriver<Model>::Finalize()
    {
        MessageHandler::Send(*this, "all", "::Finalize begin");

        model_.Finalize();

        MessageHandler::Send(*this, "all", "::Finalize end");
    }


    //! Checks whether the model has finished.
    /*!
      \return True if no more data assimilation is required, false otherwise.
    */
    template <class Model>
    bool ForwardDriver<Model>::HasFinished()
    {
        return model_.HasFinished();
    }


    //! Returns the model.
    /*!
      \return The model.
    */
    template <class Model>
    Model& ForwardDriver<Model>::GetModel()
    {
        return model_;
    }


    //! Returns the output saver.
    /*!
      \return The output saver.
    */
    template <class Model>
    OutputSaver&
    ForwardDriver<Model>::GetOutputSaver()
    {
        return output_saver_;
    }


    //! Returns the name of the class.
    /*!
      \return The name of the class.
    */
    template <class Model>
    string ForwardDriver<Model>::GetName() const
    {
        return "ForwardDriver";
    }


    //! Receives and handles a message.
    /*
      \param[in] message the received message.
    */
    template <class Model>
    void  ForwardDriver<Model>::Message(string message)
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
    }


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_FORWARDDRIVER_CXX
#endif
