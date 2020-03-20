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

#define VERDANDI_DEBUG_LEVEL_4
#define VERDANDI_WITH_ABORT
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK
#define VERDANDI_DENSE
#define VERDANDI_WITH_TRAJECTORY_MANAGER

#include "Verdandi.hxx"
#include VERDANDI_GTEST_MODEL_PATH
#include "seldon/computation/optimization/NLoptSolver.cxx"
#include "seldon/SeldonSolver.hxx"
#include "gtest/gtest.h"

#include "test_perf.hpp"

// Main function used to launch the gtest framework.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
