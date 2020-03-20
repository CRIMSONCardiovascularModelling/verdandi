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


#include "module_performance/perf_timer.cxx"
#include "method/UnscentedKalmanFilter.cxx"
#include "method/ExtendedKalmanFilter.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/ForwardDriver.cxx"
using namespace Verdandi;


//! This class is an example of a performance test.
/*! This class provides an example of how to use the macros found in
  "perf_timer.cxx" in order to provide some useful information about
  computational time into a gtest fixture.
*/
class PerformanceTest: public testing::Test
{
protected:
    typedef double real;
    typedef Verdandi::Vector<double> state;

public:
    void SetUp()
    {
        ForwardDriver<VERDANDI_GTEST_MODEL> driver;
        driver.Initialize("configuration/truth.lua");
        while (!driver.HasFinished())
        {
            driver.InitializeStep();
            driver.Forward();
            driver.FinalizeStep();
        }
        driver.Finalize();
        // For each test a timer called "main" is started at the beginning.
        TIMER("main");
    }

    void TearDown()
    {
        // The main timer is stopped.
        TIMER("main");

        // Both timers are printed in a gtest formulation.
        PRINT_TIMER("main");
        PRINT_TIMER("forward");

        // In order to use the same names in the next test, the timers are
        // reset.
        RESET_TIMER("main");
        RESET_TIMER("forward");
    }
};


/*! \brief This first test uses the extended Kalman filter and looks further
    for the time spent in the Forward function of 'ExtendedKalmanFilter'. */
TEST_F(PerformanceTest, PerformanceEKF)
{
    // Declares the driver.
    Verdandi::ExtendedKalmanFilter<VERDANDI_GTEST_MODEL,
                                   Verdandi::LinearObservationManager<real> >
        driver;

    // Initializes the model, the driver and the observation manager.
    driver.Initialize(VERDANDI_GTEST_CONFIG_PATH);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        // A timer named "forward" is started before the call to
        // 'driver.Forward'.
        TIMER("forward");
        driver.Forward();
        // The timer is stopped right after the call, in order to get the time
        // spent in it.
        TIMER("forward");
        driver.Analyze();
        driver.FinalizeStep();
    }
    driver.Finalize();
}


/*! \brief This test is similar to 'PerformanceEKF' but looks at the
  'UnscentedKalmanFilter' method. */
TEST_F(PerformanceTest, PerformanceUKF)
{
    Verdandi::UnscentedKalmanFilter<VERDANDI_GTEST_MODEL,
                                    Verdandi::LinearObservationManager<real> >
        driver;
    driver.Initialize(VERDANDI_GTEST_CONFIG_PATH);
    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        TIMER("forward");
        driver.Forward();
        TIMER("forward");
        driver.Analyze();
        driver.FinalizeStep();
    }
}
