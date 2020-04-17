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


/*! This file includes tests for basics functions of a perturbation
  manager. Those functions are used for example by "test_TR1.hpp". As any test
  on a random number generator, those tests can fail once in a while, but they
  should pass most of the time.
*/


using namespace Verdandi;


// This test checks that the output of 'Uniform(min, max)' is correctly
// between 'min' and 'max' excluded.
template <class PerturbationManager>
void TestRange(PerturbationManager rng)
{
    double min = 13.2;
    double max = 17.5;
    double value = rng.Uniform(min, max);
    ASSERT_GE(value, min);
    ASSERT_LT(value, max);
}


// This test checks that the output of 'Uniform(min, max)' is correctly
// between 'min' and 'max' excluded.
template <class PerturbationManager>
void TestRangeInt(PerturbationManager rng)
{
    int min = 10;
    int max = 20;
    int value = rng.UniformInt(min, max);

    ASSERT_GE(value, min);
    ASSERT_LT(value, max);
}


// This test checks that there is a real variety of outputs from the
// 'UniformInt()' function.
template <class PerturbationManager>
void TestDiversity(PerturbationManager rng)
{
    int min = 0;
    unsigned int max = 170;
    unsigned int count = 100;
    int* stack = new int[max];
    for (unsigned int i = 0; i < max; i++)
        stack[i] = 0;
    int value;

    for (unsigned int i = 0; i < count; i++)
    {
        value = rng.UniformInt(min, max);
        stack[value]++;
    }

    for (unsigned int i = 0; i < max; i++)
        ASSERT_LT(stack[i], 6);
    delete[] stack;
}


// This test checks that the 'Normal(mean, variance, parameter)' function will
// produce an output equal to mean on average after a lot of randoms samples.
template <class PerturbationManager>
void TestMean(PerturbationManager rng)
{
    Vector<double, VectFull> parameter;
    unsigned int count = 1000000;
    double  mean = 17;
    double variance = 5;
    double value = 0;
    double accuracy = 0.01;
    for(unsigned int i = 0; i < count; i++)
        value+= rng.Normal(mean, variance, parameter);
    value /= count;
    ASSERT_NEAR(value, mean, accuracy);
}


// This tests perform a chi^2 test in order to check the independence of the
// samples.
template <class PerturbationManager>
void TestChiSquared(PerturbationManager rng)
{
    int min = 0;
    unsigned int max = 6;
    unsigned int count = 6000;
    int* stack = new int[max];
    int value;
    double critical_value = 11.07;
    for (unsigned int i = 0; i < max; i++)
        stack[i] = 0;
    for (unsigned int i = 0; i < count; i++)
    {
        value = rng.Uniform(min, max);
         stack[value]++;
    }
    double expected = count / max;
    double criteria = 0;
    for (unsigned int i = 0; i < max; i++)
         criteria += (stack[i] - expected) * (stack[i] - expected) / expected;
    ASSERT_LT(criteria, critical_value);
    delete[] stack;
}
