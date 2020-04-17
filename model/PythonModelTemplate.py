# Copyright (C) 2011-2012, INRIA
# Author(s): Kévin Charpentier, Vivien Mallet
#
# This file is an example part of the data assimilation library Verdandi.
#
# Verdandi is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# Verdandi is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
# more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Verdandi. If not, see http://www.gnu.org/licenses/.
#
# For more information, visit the Verdandi web site:
#      http://verdandi.gforge.inria.fr/


## This class is a template model in Python.
class PythonModelTemplate:
    ## Constructs the model.
    def __init__(self):
        return


    ## Initializes the model.
    # @param[in] path The path to the configuration file.
    def Initialize(self, path):
        return


    ## Initializes the current time step for the model.
    def InitializeStep(self):
        return


    # ## Processing ## #


    ## Advances one step forward in time.
    # \f[x^f_{h+1} = \mathcal{M}_h(x^a_h, p_h)\,.\f]
    def Forward(self):
        return


    ## Performs one step backward in adjoint model.
    # @param[in] observation_term \f$ H^T R^{-1}(y - Hx) \f$.
    def BackwardAdjoint(self, observation_term):
        return


    ## Applies the model to a given vector.
    # The current state of the model is modified.
    # @param[in, out] state a vector.
    # @param[in] preserve_state Boolean to indicate if the model state has to
    #  be preserved.
    def ApplyOperator(self, state, preserve_state = True):
        return


    ## Applies the tangent linear model to a given vector.
    # @param[in,out] x on entry, a vector to which the tangent linear model
    #  should be applied; on exit, the result.
    def ApplyTangentLinearOperator(self, x):
        return


    ## Returns the tangent linear model.
    # @param[out] M the matrix of the tangent linear model.
    def GetTangentLinearOperator(self):
        return self.M_


    ## Checks whether the model has finished.
    # @return True if the simulation is done, false otherwise.
    def HasFinished(self):
        return self.time_ >= self.final_time_


    ## Finalizes the current time step for the model.
    def FinalizeStep(self):
        return

    ## Finalizes the model.
    def Finalize(self):
        return


    # ## Access methods ## #


    ## Returns the current time.
    # @return The current time.
    def GetTime(self):
        return self.time_


    ## Sets the current time.
    # @param[in] time the current time.
    def SetTime(self, time):
        self.time_ = time


    ## Returns the dimension of the state.
    # @return The dimension of the state.
    def GetNstate(self):
        return self.Nstate_


    ## Returns the dimension of the full state.
    # @return The dimension of the full state.
    def GetNfull_state(self):
        return self.NFull_state_


    ## Provides the controlled state vector.
    # @return The controlled state vector.
    def GetState(self):
        return self.state_

    ## Performs some calculations when the update of the model state is done.
    def StateUpdated(self):
        return


    ## Provides the state lower bound.
    # @return lower_bound the state lower bound (componentwise).
    def GetStateLowerBound(self):
        return self.state_lower_bound_


    ## Provides the state upper bound.
    # @return lower_bound the state upper bound (componentwise).
    def GetStateUpperBound(self, state_upper_bound):
        return self.state_upper_bound_


    ## Provides the full state vector.
    # @return The full state vector.
    def GetFullState(self):
        return self.full_state_


    ## Performs some calculations when the update of the model state is done.
    def FullStateUpdated(self, full_state):
        return


    ## Provides the adjoint state vector.
    # @return The adjoint state vector.
    def GetAdjointState(self):
        return self.adjoint_state_


    ## Performs some calculations when the update of the adjoint state
    #  is done.
    def AdjointStateUpdated(self):
        return


    ## Returns the number of parameters to be perturbed.
    # @return The number of parameters to be perturbed.
    def GetNparameter(self):
        return self.Nparameter_


    ## Returns the PDF of the i-th parameter.
    # @param[in] i parameter index.
    # @return The PDF of the i-th parameter.
    def GetParameterPDF(self, i):
        return self.parameter_pdf_[i]


    ## Gets the i-th uncertain parameter.
    # @param[in] i parameter index.
    # @return The vector associated with the i-th parameter.
    def GetParameter(self, i):
        return self.parameter_[i]


    ## Triggers the actions after an uncertain parameter has been updated.
    # @param[in] i parameter index.
    def ParameterUpdated(self, i):
        return


    ## Returns the name of a parameter to be perturbed.
    # @param[in] i parameter index.
    # @return The name of the parameter.
    def GetParameterName(self, i):
        return self.parameter_name_[i]


    ## Returns the correlation between the uncertain parameters.
    # @param[in] parameter index.
    # @return The correlation between the uncertain parameters.
    def GetParameterCorrelation(self, i):
        return self.parameter_correlation_[i]


    ## Returns parameters associated with the PDF of some model parameter.
    # @param[in] i parameter index.
    # @return The parameters associated with the i-th parameter.
    def GetParameterPDFData(self, i):
        return self.parameter_pdf_data_[i]


    ## Returns the covariance matrix associated with the i-th parameter.
    # @param[in] i parameter index.
    # @return The covariance matrix associated with the i-th parameter.
    def GetParameterVariance(self, i):
        return self.parameter_variance_[i]


    ## Returns the perturbation option of the i-th parameter.
    # @param[in] i parameter index.
    # @return The perturbation option of the i-th parameter.
    def GetParameterOption(self, i):
        return self.parameter_option_[i]


    # ## Errors ## #


    ## Returns the state error variance.
    # @return The state error variance.
    def GetStateErrorVariance(self):
        return self.P_


    ## Returns a row of the state error variance.
    # @param[in] row row index.
    # @return The row with index \a row in the state error variance.
    def GetStateErrorVarianceRow(self, row):
        return self.P_[row]


    ## Returns the matrix L in the decomposition of the state error covariance
    # matrix (\f$B\f$) as a product \f$LUL^T\f$.
    # @return The matrix \f$L\f$.
    def GetStateErrorVarianceProjector(self):
        return self.L_


    ## Returns the matrix U in the decomposition of the state error covariance
    # matrix (\f$B\f$) as a product \f$LUL^T\f$.
    # @return The matrix \f$U\f$.
    def GetStateErrorVarianceReduced(self):
        return self.U_


    ## Returns the inverse of the background error variance (\f$B^{-1}\f$).
    # @return The inverse of the background error variance (\f$B^{-1}\f$
    def GetStateErrorVarianceInverse(self):
        return self.state_error_variance_inverse_


    ## Returns the name of the class.
    # @return The name of the class.
    def GetName(self):
        return "PythonModelTemplate"

    ## Receives and handles a message.
    # @param[in] message the received message
    def Message(self, message):
        return
