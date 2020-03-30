#ifndef SELDON_FILE_SELDON_DISTRIBUTED_HEADER_HXX

#include "share/MpiCommunication.hxx"

// including distributed vectors
#include "vector/DistributedVector.hxx"

// including distributed sparse matrices
#include "matrix_sparse/DistributedMatrix.hxx"

// including distributed solver if SeldonSolverHeader.hxx has been included
#ifdef SELDON_FILE_SELDON_SOLVER_HEADER_HXX
#include "computation/solver/DistributedSolver.hxx"
#endif

#define SELDON_FILE_SELDON_DISTRIBUTED_HEADER_HXX
#endif
