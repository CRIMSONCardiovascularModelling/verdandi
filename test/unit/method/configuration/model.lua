----------------------------------- MODEL ------------------------------------


quadratic_model = {

   error = {

      -- Is the model error variance a scaled identity matrix?
      scaled_identity = true,
      -- If so, put the diagonal value:
      diagonal_value = 0.0,
      -- Otherwise, the operator value (file name or table):
      value = {}

   },

   error_sqrt = {

      -- Columns of the square root of the model error variance.
      value = {}

   },

   state_error = {

      -- Is the state error variance a scaled identity matrix?
      scaled_identity = true,
      -- If so, put the diagonal value:
      diagonal_value = 1.,
      -- Otherwise, the operator value (file name or table):
      value = {}

   },

   state_error_sqrt = {

      -- Columns of the square root of the state error variance.
      value = {}

   },

   output_saver = {

      variable_list = {"state", "S", "L", "b"},
      file = output_directory .. "model-%{name}.bin",

   },




--- 2D version.

definition = {

   initial_state = {0.5,
                    0.5,
                 },

   with_quadratic_term = false,
   with_linear_term = true,
   with_constant_term = true,

   quadratic_term = {1., 0.,
                     0., 1.,
                     1., 0.,
                     0., 0.},
   linear_term = {-1., 1.,
                  0., -1.},
   constant = {0.,
               0.},

   Delta_t = 0.0015,
   initial_time = 0.,
   final_time = 2.5

},

uncertainty = {

   uncertain_parameter_list = {"quadratic_term", "linear_term",
                               "constant"},

   quadratic_term = {

      distribution = "NormalHomogeneous",
      mean = {0., 0.},
      variance = {0.1^2, 0.,
                  0., 0.1^2},
      parameter = {-1., 1.}
   },

   linear_term = {

      distribution = "NormalHomogeneous",
      mean = {0., 0.},
      variance = {0.1^2, 0.,
                  0., 0.1^2},
      parameter = {-1., 1.}
   },

   constant = {

      distribution = "NormalHomogeneous",
      mean = {0., 0.},
      variance = {0.1^2, 0.,
                  0., 0.1^2},
      parameter = {-1., 1.}
   }

}

}