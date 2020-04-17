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


//! This class performs a set of basic tests on the model.
/*!
  This classes uses the gtest framework to perform some tests on the model:
  1. Basic definitions (name, state defined, state dimension...);
  2. 'GetState' followed by 'StateUpdated' works correctly;
  3. Checks that 'Forward' seems correct (consistency between methods).
*/
class BasicModelTest: public testing::Test
{
public:
    typedef double real;
    typedef Verdandi::Vector<double> state;

protected:
    //! Tested model.
    VERDANDI_GTEST_MODEL model_;
    //! Generates random numbers, using Mersenne twister.
    PerturbationManager rng_;
    //! Restricts the values of the states below a bound.
    double  state_upper_bound_;
    //! Restricts the values of the states above a bound.
    double  state_lower_bound_;
    //! Maximum number of steps that the model is allowed to finish.
    int Nstep_;

public:
    // Reads the configuration and set up the model.
    void SetUp()
    {
        Verdandi::VerdandiOps configuration(VERDANDI_GTEST_CONFIG_PATH);
        configuration.Set("state_upper_bound", state_upper_bound_);
        configuration.Set("state_lower_bound", state_lower_bound_);
        configuration.Set("Nstep_init", Nstep_);
        model_.Initialize(VERDANDI_GTEST_CONFIG_PATH);
    }
};


// Checks that the state dimension is set to a strictly positive value.
TEST_F(BasicModelTest, CheckIfNstateIsDefined)
{
    ASSERT_GT(model_.GetNstate(), 0);
}


// Checks that the size of the state is consistent across the interface.
TEST_F(BasicModelTest, CheckIfStateSizeIsCorrect)
{
    ASSERT_EQ(model_.GetNstate(), model_.GetState().GetM());
}


// Checks that, after a call to 'StateUpdated', 'GetState' returns the
// previously modified state.
TEST_F(BasicModelTest, CheckSetStateCorrectness)
{


    // Generates an integer lower than Nstate at which we insert a random
    // value.
    int position = rng_.UniformInt(0, model_.GetNstate());
    // Generates a random value.
    double value = rng_.Uniform(state_lower_bound_, state_upper_bound_);
    // Sets the generated value into the model.
    model_.GetState()(position) = value;
    model_.StateUpdated();
    // Checks that modified state component has not changed.
    ASSERT_EQ(value, model_.GetState()(position));
}


// Checks that the model got a defined name.
TEST_F(BasicModelTest, CheckGetName)
{
    string name = model_.GetName();
    ASSERT_NE(name, "");
}


// Checks that the model can go forward, and checks that the size of the state
// remains consistent across the interface.
TEST_F(BasicModelTest,CheckForward)
{
    model_.Forward();
    ASSERT_EQ(model_.GetNstate(), model_.GetState().GetM());
}


// Checks that 'Forward' and 'ApplyOperator' give the same result.
TEST_F(BasicModelTest, CheckForwardAndApplyOperator)
{
    state state_copy = model_.GetState();
    model_.ApplyOperator(state_copy);
    model_.Forward();
    for (int i = 0; i < state_copy.GetM(); i++)
    {
        ASSERT_EQ(state_copy(i), model_.GetState()(i));
    }
}


// Checks that the model ever stops.
TEST_F(BasicModelTest, CheckFinish)
{
    int count = 0;
    do
    {
        ASSERT_TRUE(count < Nstep_);
        model_.Forward();
        count++;
    }
    while (!model_.HasFinished());
}
