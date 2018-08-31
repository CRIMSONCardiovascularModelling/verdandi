----------------------------------- GLOBAL -----------------------------------


output_directory = "result/"


----------------------------------- MODEL ------------------------------------


dofile("configuration/model.lua")


----------------------------------- METHOD -----------------------------------


-- Forward simulation.
forward = {

   display = {

      iteration = false,
      time = true

   },

   output_saver = {

      variable_list = {"forecast_time", "forecast_state"},
      file = output_directory .. "truth-%{name}.%{extension}",
      mode = "binary",
      mode_scalar = "text"

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }

}
