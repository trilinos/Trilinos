#include <gtest/gtest.h>

#include <mpi.h>

#include <stk_util/environment/WallTime.hpp>

#if defined(_OPENMP) && !defined(__INTEL_COMPILER)
#include <omp.h>

namespace
{

TEST(MpiPlusOpenMp, HelloWorld)
{
#pragma omp parallel
    {
        std::cout << "Hello World\n";
    }
}

int getNumProcs()
{
    int numProcs = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    return numProcs;
}

int getProcId()
{
    int procId = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &procId);
    return procId;
}

TEST(MpiPlusOpenMp, Reduction)
{
    int numThreads = 0;
    double threadSum = 0;

    #pragma omp parallel reduction(+:threadSum)
    {
        numThreads = omp_get_num_threads();
        threadSum = 1.0;
    }

    double mpiSum = -1;
    const int numItemsPerProc = 1;
    MPI_Allreduce(&threadSum, &mpiSum, numItemsPerProc, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    int numProcs = getNumProcs();
    EXPECT_EQ(numThreads*numProcs, mpiSum);
}

TEST(MpiPlusOpenMp, VectorSumReduction)
{
    int numProcs = getNumProcs();
    int numThreads = 0;
//    size_t sizeOfGlobalVector = 4032000000;
    size_t sizeOfGlobalVector = 40320000;
    ASSERT_EQ(0u, sizeOfGlobalVector%numProcs);
    size_t sizeOfLocalVector = sizeOfGlobalVector / numProcs;
    double initVal = 1.0;

    double *vec = new double[sizeOfLocalVector];
    #pragma omp parallel for
    for (size_t i = 0; i < sizeOfLocalVector; i++)
    {
        vec[i] = initVal;
    }

    double start_time_wall = stk::wall_time();

    double threadSum = 0;
    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
        #pragma omp for reduction(+:threadSum)
        for (size_t i = 0; i < sizeOfLocalVector; i++)
        {
            threadSum += vec[i];
        }
    }

    double mpiSum = -1;
    const int numItemsPerProc = 1;
    MPI_Allreduce(&threadSum, &mpiSum, numItemsPerProc,
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    double end_time_wall = stk::wall_time();

    int procId = getProcId();
    if(procId == 0)
    {
        std::stringstream ss;
        ss << "Num MPI processes: " << numProcs << " ";
        ss << "Num threads per process: " << numThreads <<  " ";
        ss << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;
        std::cerr << ss.str();
    }

    EXPECT_EQ(sizeOfGlobalVector*initVal, mpiSum);
    delete [] vec;
}

}
#endif
