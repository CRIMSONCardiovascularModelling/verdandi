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


#include "method/MonteCarlo.cxx"

#include "method/UnscentedKalmanFilter.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "model/QuadraticModel.cxx"
#include "method/ExtendedKalmanFilter.cxx"
#include "method/OptimalInterpolation.cxx"
#include "method/ForwardDriver.cxx"
#include "method/ReducedOrderExtendedKalmanFilter.cxx"
#include "method/Nudging.cxx"


using namespace Verdandi;
typedef Verdandi::Vector<double> state;

//! This test allows to compare methods.
/*! This class allows to compare methods in order to verify that under the
  right hypothesis they produce the same output.
*/
class MethodCompare: public testing::Test
{
protected:
    typedef double real;
    typedef Verdandi::Vector<double> state;

public:
    //! This is the shared state used to compare states across tests.
    static state shared_state_;

    //! This method allows to compute an output from EKF and to store it.
    static void SetUpTestCase()
    {
        // First we produce an output with a forward driver, which will be
        // used as observations.
        ForwardDriver<VERDANDI_GTEST_MODEL> driver;
        driver.Initialize("configuration/truth.lua");
        while (!driver.HasFinished())
        {
            driver.InitializeStep();
            driver.Forward();
            driver.FinalizeStep();
        }
        driver.Finalize();

        // Those observations are used to run the EKF method. The final state
        // is stored in 'shared_state_'.
        Verdandi::ExtendedKalmanFilter<Verdandi::QuadraticModel<real>,
                                       Verdandi::LinearObservationManager
                                       <real> > driver2;
        driver2.Initialize(VERDANDI_GTEST_CONFIG_PATH);
         while (!driver2.HasFinished())
        {
            driver2.InitializeStep();
            driver2.Forward();
            driver2.FinalizeStep();
        }
        shared_state_ = driver2.GetModel().GetState();
    }
};


//! Uses Google Test assertion to compare two states.
void AssertStateEqual(state A, state B)
{
    double accuracy = 0.00001;
    int sizeA = A.GetM(), sizeB = B.GetM();
    ASSERT_FLOAT_EQ(sizeA, sizeB);
    for (int i = 0; i < sizeA; i++)
        ASSERT_NEAR(A(i), B(i), accuracy);
}


//! This state is the final state after being corrected by the EFK method.
state MethodCompare::shared_state_;


//! This test checks that the UKF produces the same output as EKF.
TEST_F(MethodCompare, test_UKF)
{
    // The UFK method is run using the same configuration and observations as
    // the EKF in the 'SetUp()' method.
    Verdandi::UnscentedKalmanFilter<Verdandi::QuadraticModel<real>,
                                    Verdandi::LinearObservationManager
                                    <real> > driver;

    driver.Initialize(VERDANDI_GTEST_CONFIG_PATH);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.FinalizeStep();
     }

    state compared_state = driver.GetModel().GetState();
    // With the hypothesis of model linearity, we should obtain the same exact
    // state.
    AssertStateEqual(compared_state, shared_state_);
    driver.Finalize();
}


//! This test checks that the ROEKF produces the same output as EKF.
TEST_F(MethodCompare, test_ROEKF)
{
    // The ROEKF method is run using the same configuration and observations
    // as the EKF in the 'SetUp()' method.
    Verdandi::ReducedOrderExtendedKalmanFilter<Verdandi::QuadraticModel<real>,
                                    Verdandi::LinearObservationManager
                                    <real> > driver;

    driver.Initialize(VERDANDI_GTEST_CONFIG_PATH);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.FinalizeStep();
     }

    state compared_state = driver.GetModel().GetState();
    // With the hypothesis of model linearity and a non-reduced state, we
    // should obtain the same exact state.
    AssertStateEqual(compared_state, shared_state_);
    driver.Finalize();
}


TEST_F(MethodCompare, test_nudging)
{
    // The nudging method is run using the same configuration and observations
    // as the EKF in the 'SetUp()' method.
    Verdandi::Nudging<Verdandi::QuadraticModel<real>,
                                    Verdandi::LinearObservationManager
                                    <real> > driver;

    driver.Initialize(VERDANDI_GTEST_CONFIG_PATH);
    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
        driver.FinalizeStep();
     }

    state compared_state = driver.GetModel().GetState();
    // With the hypothesis of model linearity we should obtain the same exact
    // state.
    AssertStateEqual(compared_state, shared_state_);
    driver.Finalize();
}
