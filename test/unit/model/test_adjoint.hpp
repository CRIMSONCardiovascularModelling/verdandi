// Copyright (C) 2014 INRIA
// Author(s): Nicolas Claude
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

#include <cmath>

using namespace Verdandi;


//! This class allows to test the adjoint model.
/*!
  This class uses the gtest framework to perform two tests on the adjoint:
  1. Verification that <M v; M v> = <v; M^T M v>.
  2. Verification of the tangent linear operator with finite differences.

  M: tangent linear model. M^T: adjoint model. v: state.
  < ; >: scalar product.
*/
class AdjointTest: public testing::Test
{
public:
    typedef double real;
    typedef Verdandi::Vector<double> state;

protected:
    //! The accuracy of the test.
    double accuracy_;
    //! The size of the state.
    unsigned int Nstate_;
    //! Restricts the values of the states below a bound.
    double  state_upper_bound_;
    //! Restricts the values of the states above a bound.
    double  state_lower_bound_;
    //! Generates random numbers.
    PerturbationManager rng_;
    //! The working state, randomly generated.
    state crafted_state_;
    //! Tested model.
    VERDANDI_GTEST_MODEL model_;

public:
    // This method initializes the model and reads the configuration for the
    // tests, then randomizes a state which will be used in the test.
    void SetUp()
    {
        // Reads the configuration.
        Verdandi::VerdandiOps configuration(VERDANDI_GTEST_CONFIG_PATH);
        configuration.Set("accuracy", accuracy_);
        configuration.Set("state_lower_bound", state_lower_bound_);
        configuration.Set("state_upper_bound", state_upper_bound_);

        // Initializes the model and the model adjoint.
        model_.Initialize(VERDANDI_GTEST_CONFIG_PATH);
        model_.InitializeAdjoint();

        // Gets the size of the state.
        Nstate_ = model_.GetNstate();

        // Creates a random generator.
        crafted_state_.Reallocate(Nstate_);
        // Creates a randomized state: 'crafted_state_'. This state's values
        // are between 'state_lower_bound_' and 'state_upper_bound_'.
        model_.GetState().Reallocate(Nstate_);
        for (unsigned int i = 0; i < Nstate_; i++)
            crafted_state_(i) =
                rng_.Uniform(state_lower_bound_, state_upper_bound_);

        model_.GetState() = crafted_state_;
    }
};


// This test verifies that <M v; M v> = <v; M^T M v>.
TEST_F(AdjointTest, TestModelAdjoint)
{
    // Creates an null observation term for the adjoint. It permits to go
    // backward without observations.
    state observation_term;
    observation_term.Reallocate(Nstate_);
    observation_term.Fill(0);

    // Declares working variables.
    double output_left = 0;
    double output_right = 0;

    // Declares a copy of the randomized state.
    state crafted_state_copy;
    crafted_state_copy.Reallocate(Nstate_);
    crafted_state_copy = crafted_state_;

    // Sets 'crafted_state_' to  M v.
    model_.ApplyTangentLinearOperator(crafted_state_);

    // Computes M^T M v.
    model_.GetAdjointState() = crafted_state_;
    model_.BackwardAdjoint(observation_term);

    // Computes <v; M^T M v>.
    for (int i = 0; i < crafted_state_.GetM(); i++)
        output_right += crafted_state_copy(i) *
            (model_.GetAdjointState()(i) - crafted_state_(i));

    // Computes <M v; M v>.
    for (int i = 0; i < crafted_state_.GetM(); i++)
        output_left += crafted_state_(i) * crafted_state_(i);

    // Test itself: is <M v; M v> = <v; M^T M v> as it should be?
    ASSERT_NEAR(output_left, output_right, accuracy_);
}


// This test checks the tangent linear operator using finites differences.
TEST_F(AdjointTest, TestLinearTangent)
{
    // 'c' is a randomized point where we do the test.
    state c = crafted_state_;

    // Randomizes i and j. j determines the perturbation direction. i is the
    // index of the model output under consideration.
    int j = rng_.UniformInt(0, Nstate_);
    int i = rng_.UniformInt(0, Nstate_);

    // Denominator: (Model(c + Epsilon * v)_i - Model(c)_i) / Epsilon, where
    // Model(c)_i is the i-th component of the model output Model(c).
    double denominator = 1.0;

    state working_state;
    working_state.Reallocate(Nstate_);

    // Creates the perturbation direction 'v' (0 everywhere except at j-th
    // component where is set to 1).
    state v;
    v.Reallocate(Nstate_);
    v.Fill(0);
    v(j) = 1;

    // Epsilon is a scaling factor for the perturbation along 'v'. It
    // approaches 0.
    double epsilon = 0.1;

    // 'mv' is M times v, with M the tangent linear model at point c.
    state mv = v;
    model_.ApplyTangentLinearOperator(mv);

    // 'mvi' is the i-th component of M v.
    double mvi = mv(i);

    // 'model_mt_ci' is the i-th component of the model's output at point c.
    working_state = c;
    model_.ApplyOperator(working_state);
    double model_at_ci = working_state(i);
    unsigned count = 0;
    do
    {
        // Determines when the test should be considered as a failure.
        ASSERT_TRUE(count < 50);
        working_state = c;

        // Computes working_state += epsilon * v.
        Add(epsilon, v, working_state);

        // Computes Model(c + epsilon * v).
        model_.ApplyOperator(working_state);

        // Computes the finite difference.
        denominator = (working_state(i) - model_at_ci) / epsilon;

        // Forces epsilon to approach 0.
        epsilon /= 2.;

        count++;
    }
    // This test succeeds if 'mvi' approaches 'denominator' (down 'accuracy_').
    while (std::abs((mvi / denominator) - 1) > accuracy_);
}


// This test checks the consistency between 'GetTangentLinearOperator()' and
// 'ApplyTangentLinearOperator'.
TEST_F(AdjointTest, TestTangentLinearOperator)
{
    state FromMatrixOperator = crafted_state_;
    state FromApplyOperator = crafted_state_;
    Matrix<double> TL_model = model_.GetTangentLinearOperator();

    for(int i = 0 ; i < TL_model.GetM(); i++)
        for (int j = 0 ; j < TL_model.GetN(); j++)
            FromMatrixOperator(i) = TL_model(i,j) * crafted_state_(i);

    model_.ApplyTangentLinearOperator(FromApplyOperator);

    for (int i = 0; i < crafted_state_.GetM(); i++)
        ASSERT_NEAR(FromApplyOperator(i), FromMatrixOperator(i), accuracy_);
}
