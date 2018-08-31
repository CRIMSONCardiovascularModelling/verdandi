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
#include "RandomPerturbationManager.cxx"

using namespace Verdandi;


class TestRandom: public testing::Test
{
protected:
    RandomPerturbationManager rng_;
public:
    void SetUp()
    {
    }
    void TearDown()
    {
    }
};


TEST_F(TestRandom, TestUniformRange)
{
    TestRange(rng_);
}


TEST_F(TestRandom, TestUniformRangeInt)
{
    TestRangeInt(rng_);
}


TEST_F(TestRandom, TestDiversityOfNumbers)
{
    TestDiversity(rng_);
}


TEST_F(TestRandom, TestNormalDistribution)
{
    TestMean(rng_);
}


TEST_F(TestRandom, TestChi2)
{
    TestChiSquared(rng_);
}
