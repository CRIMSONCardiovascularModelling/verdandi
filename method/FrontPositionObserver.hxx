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


#ifndef VERDANDI_FILE_METHOD_FRONTPOSITIONOBSERVER_HXX
#define VERDANDI_FILE_METHOD_FRONTPOSITIONOBSERVER_HXX


namespace Verdandi
{


    ///////////////////////////
    // FRONTPOSITIONOBSERVER //
    ///////////////////////////


    //! This class simply performs a forward simulation.
    template <class Model, class ObservationManager>
    class FrontPositionObserver: public VerdandiBase
    {

    public:
        //! Type of the model state vector.
        typedef typename Model::state model_state;
        //! Type of the observation vector.
        typedef typename ObservationManager::observation
        observation;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        //! Iteration.
        int iteration_;

        /*** Configuration ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the observation manager.
        string observation_configuration_file_;

        //! Display options.
        map<string, bool> option_display_;

        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;
        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;
        //! Current time.
        double time_;

#ifdef VERDANDI_WITH_MPI
        //! The rank in MPI_COMM_WORLD.
        int world_rank_;
        //! The rank in the model communicator.
        int model_task_;
        //! The number of rows of the MPI grid.
        int Nrow_;
        //! The number of columns of the MPI grid.
        int Ncol_;
        //! The number of processes in MPI_COMM_WORLD.
        int Nprocess_;
        //! The MPI grid row communicator of the current process.
        MPI_Comm row_communicator_;
        //! The MPI grid column communicator of the current process.
        MPI_Comm col_communicator_;
#endif

        /*** Output saver ***/

        //! Output saver.
        OutputSaver output_saver_;

        /*** Gain Filter ***/
        double filter_gain_;

    public:

        /*** Constructor and destructor ***/

        FrontPositionObserver();
        ~FrontPositionObserver();

        /*** Methods ***/

        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);

        void InitializeStep();

        void Prediction();
        void Forward();
        void Analyze();

        void FinalizeStep();
        void Finalize();

        void ComputeCorrection(const observation& innovation,
                               model_state& state);

        bool HasFinished();

        // Access methods.
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();

        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#endif // define VERDANDI_FILE_METHOD_FRONTPOSITIONOBSERVER_HXX
