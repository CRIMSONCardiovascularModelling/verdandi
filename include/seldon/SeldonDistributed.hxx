#ifndef SELDON_FILE_SELDON_DISTRIBUTED_HXX

#include "SeldonDistributedHeader.hxx"
#include "SeldonDistributedInline.hxx"

#include "share/MpiCommunication.cxx"
#include "vector/DistributedVector.cxx"
#include "matrix_sparse/DistributedMatrix.cxx"

// including distributed solver if SeldonSolver.hxx has been included
#ifdef SELDON_FILE_SELDON_SOLVER_HXX
#include "computation/solver/DistributedSolver.cxx"
#endif

#define SELDON_FILE_SELDON_DISTRIBUTED_HXX
#endif
