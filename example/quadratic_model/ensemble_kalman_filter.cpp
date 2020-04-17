#define VERDANDI_DEBUG_LEVEL_4
#define SELDON_WITH_BLAS
#define SELDON_WITH_LAPACK

#define VERDANDI_DENSE
#define VERDANDI_WITH_ABORT

//#define VERDANDI_WITH_MPI

#if defined(VERDANDI_WITH_MPI)
#include <mpi.h>
#endif

#include "Verdandi.hxx"

#include "model/QuadraticModel.cxx"
#include "observation_manager/LinearObservationManager.cxx"
#include "method/EnsembleKalmanFilter.cxx"

#ifdef VERDANDI_HAS_CXX11
#include "method/RandomPerturbationManager.cxx"
#define RNG RandomPerturbationManager
#else
#include "method/TR1PerturbationManager.cxx"
#define RNG TR1PerturbationManager
#endif

int main(int argc, char** argv)
{

    TRY;

    if (argc != 2)
    {
        string mesg  = "Usage:\n";
        mesg += string("  ") + argv[0] + " [configuration file]";
        std::cout << mesg << std::endl;
        return 1;
    }

    typedef double real;

    Verdandi::EnsembleKalmanFilter<Verdandi::QuadraticModel<real>,
        Verdandi::LinearObservationManager<real>, Verdandi::RNG> driver;

    driver.Initialize(argv[1]);

    while (!driver.HasFinished())
    {
        driver.InitializeStep();
        driver.Forward();
    }

    END;

    return 0;

}
