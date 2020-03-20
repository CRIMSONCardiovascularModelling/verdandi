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

#include <iostream>

using namespace Verdandi;


//! This class allows to perform two time-based tests.
/*!
  This class uses the gtest framework to perform 2 tests:
  1. tries to apply the model from time h-1 to h with h=0..Nstep_time
  2. randomly samples two points in time, and applies the model from the first
     time to the second.

  Each test will succeed if the stored state equals the computed one.
*/
class TestTime: public testing::Test
{
protected:
    typedef double real;
    typedef Verdandi::Vector<double> state;

    //! The accuracy of the test.
    double accuracy_;
    //! The size of the state.
    unsigned int Nstate_;
    //! Generates random numbers, using Mersenne twister.
    PerturbationManager rng_;
    //! Tested model.
    VERDANDI_GTEST_MODEL model_;
    //! Stores the states to check the consistency of the model.
    vector<state> state_storage_;
    //! Times corresponding to the stored states.
    vector<double> time_storage_;
    //! Span of the storage.
    int Nstep_;

public:
    // This method initializes the model and reads the configuration for the
    // tests, then stores all the states and times needed for the tests.
    virtual  void SetUp()
    {
        // Initializes the model, reads the configuration, initializes
        // variables.
        model_.Initialize(VERDANDI_GTEST_CONFIG_PATH);
        VerdandiOps configuration(VERDANDI_GTEST_CONFIG_PATH);
        configuration.Set("Nstep_time", Nstep_);
        configuration.Set("accuracy", accuracy_);
        Nstate_ = model_.GetNstate();

        // Computes and stores the states and times.
        for (int i = 0; i < Nstep_; i++)
        {
            state_storage_.push_back(model_.GetState());
            time_storage_.push_back(model_.GetTime());
            model_.Forward();
        }
    }
};


// This test checks that the model will give the correct state,
// given the state N-1.
TEST_F(TestTime, TestTimeBackward)
{
    // For each time step, asks the model to go forward in time once and
    // checks if the given state is correct.
    for (int i = 0; i < Nstep_ - 1; i++)
    {
        model_.GetState() = state_storage_[Nstep_ - i - 2];
        model_.StateUpdated();
        model_.SetTime(time_storage_[Nstep_ - i - 2]);
        model_.Forward();
        for (unsigned int j = 0; j < Nstate_; j++)
            ASSERT_NEAR(model_.GetState()(j),
                        state_storage_[Nstep_ - i - 1](j), accuracy_);
    }
}


// This test takes two random points in time and tries to apply the model from
// the first to the second.
TEST_F(TestTime, TestTimeRandom)
{
    // Randomly samples two steps.
    unsigned int step1 = 0;
    unsigned int step2 = 0;
    while (step1 == step2 || step2 < step1)
    {
        step1 = rng_.UniformInt(0, Nstep_);
        step2 = rng_.UniformInt(0, Nstep_ - 1);
    }

    // Inserts the state at 'step1' into the model.
    model_.GetState() = state_storage_[step1];
    model_.StateUpdated();

    // Inserts the corresponding time into the model.
    model_.SetTime(time_storage_[step1]);

    // Applies the model forward in time.
    for (unsigned int i = step1; i < step2; i++)
        model_.Forward();

    // Checks if the computed state is correct.
    for (unsigned int j = 0; j < Nstate_; j++)
        ASSERT_NEAR(model_.GetState()(j), state_storage_[step2](j),
                    accuracy_);
}
