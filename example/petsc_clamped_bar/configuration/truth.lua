----------------------------------- GLOBAL -----------------------------------


Delta_t_petsc_clamped_bar = 0.01
final_time_petsc_clamped_bar = 10
-- Saving period.
Nskip_save = 1

output_directory = "result/"
output_mode = "binary"
output_mode_scalar = "text"

observation_file = output_directory .. "truth-observation.bin"


----------------------------------- MODEL ------------------------------------


dofile("configuration/model.lua")


-------------------------------- OBSERVATION ---------------------------------


dofile("configuration/observation.lua")
-- In order to deactivate all observations.
observation.option.with_observation = false


----------------------------------- METHOD -----------------------------------


forward = {

   output_saver = {

      variable_list = {"forecast_time", "forecast_state"},
      file = output_directory .. "truth-%{name}.%{extension}",
      time = "step " .. Delta_t_petsc_clamped_bar * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   display = {

      iteration = false,
      time = true,
      mpi_grid = true

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   },

  mpi_grid = {

      -- The number of processes for each model task.
      Nrow = 3,
      -- The number of model tasks.
      Ncol = 1
   }

}

observation_generator = forward
observation_generator.output_saver.variable_list
    = {"observation_time", "observation"}
