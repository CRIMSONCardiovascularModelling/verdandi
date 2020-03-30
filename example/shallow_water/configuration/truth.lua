----------------------------------- GLOBAL -----------------------------------


Delta_t_shallow_water = 0.03
final_time_shallow_water = 1500 * Delta_t_shallow_water
-- Saving period.
Nskip_save = 10

output_directory = "result/"
output_mode = "binary"
output_mode_scalar = "text"


----------------------------------- MODEL ------------------------------------


dofile("configuration/model.lua")


----------------------------------- METHOD -----------------------------------


forward = {

   display = {

      iteration = false,
      time = true

   },

   output_saver = {

      variable_list = {"forecast_time", "forecast_state"},
      file = output_directory .. "truth-%{name}.%{extension}",
      time = "step " .. Delta_t_shallow_water * Nskip_save .. " 1.e-6",
      mode = output_mode,
      mode_scalar = output_mode_scalar

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }

}
