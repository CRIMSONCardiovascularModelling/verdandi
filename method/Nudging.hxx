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


#ifndef VERDANDI_FILE_METHOD_NUDGING_HXX


namespace Verdandi
{


    /////////////
    // NUDGING //
    /////////////


    template <class Model, class ObservationManager>
    class Nudging: public VerdandiBase
    {
    public:
        typedef typename Model::state::value_type Ts;
        typedef typename Model::state model_state;
        typedef typename ObservationManager::observation observation;
        typedef Matrix<Ts, General, RowMajor, MallocAlloc<Ts> > gain_matrix;

    protected:

        /*** Main components ***/

        //! Underlying model.
        Model model_;

        //! Observation manager.
        ObservationManager observation_manager_;

        //! Output saver.
        OutputSaver output_saver_;

        /*** Configuration paths ***/

        //! Path to the configuration file.
        string configuration_file_;
        //! Path to the model configuration file.
        string model_configuration_file_;
        //! Path to the configuration file for the observation manager.
        string observation_configuration_file_;

        /*** Options ***/

        //! Display options.
        map<string, bool> option_display_;

        //! Should an analysis be computed at the first step?
        bool analyze_first_step_;

        /*** Assimilation variables ***/

        //! Iteration.
        int iteration_;

        //! Dimension of the state.
        int Nstate_;
        //! Number of observations.
        int Nobservation_;

        //! Type of nudging: "standard", "dt" or "source".
        string nudging_type_;

        //! Gain matrix for the nudging term.
        gain_matrix nudging_matrix_;
        //! Is the gain matrix stationary?
        bool matrix_fixed_;
        //! Multiplicative factor applied to the nudging term.
        Ts nudging_gain_;
        //! Time of the previous time step.
        Ts previous_time_;
        //! Duration of the time step.
        Ts delta_t_;

    public:

        Nudging();
        ~Nudging();
        void Initialize(string configuration_file,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void Initialize(VerdandiOps& configuration,
                        bool initialize_model = true,
                        bool initialize_observation_manager = true);
        void InitializeStep();
        void Forward();
        void Analyze();
        void Prediction();

        void FinalizeStep();
        void Finalize();
        bool HasFinished();
        Model& GetModel();
        ObservationManager& GetObservationManager();
        OutputSaver& GetOutputSaver();
        string GetName() const;
        void Message(string message);
    };


} // namespace Verdandi.


#define VERDANDI_FILE_METHOD_NUDGING_HXX
#endif
