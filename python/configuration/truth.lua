----------------------------------- GLOBAL -----------------------------------


output_directory = "result/"


----------------------------------- MODEL ------------------------------------


dofile("configuration/quadratic_model.lua")


python_model = {

   module = "QuadraticModel",
   directory = "../../model/",
   class_name = "QuadraticModel"

}


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

   },

   output = {

      configuration = output_directory .. "truth.lua",
      log = output_directory .. "truth.log"

   }
}

forward.display.show_iteration = false
forward.display.show_time = false
python_model.directory = "../model/"