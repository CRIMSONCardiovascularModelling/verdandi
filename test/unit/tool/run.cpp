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
#include "seldon/SeldonSolver.hxx"

#include "gtest/gtest.h"
#include "test_helper.hpp"

#ifndef __clang__
#include "test_TR1.hpp"
#endif

#include "test_function_inverse.hpp"

#ifdef VERDANDI_HAS_CXX11
#include "test_random.hpp"
#endif

// Main function used to launch the Google Test framework.
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
