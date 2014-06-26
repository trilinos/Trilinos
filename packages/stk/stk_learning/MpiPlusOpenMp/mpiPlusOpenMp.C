#include <gtest/gtest.h>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
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
  return stk::parallel_machine_size(MPI_COMM_WORLD);
}

int getProcId()
{
  return stk::parallel_machine_rank(MPI_COMM_WORLD);
}

int getNumThreads()
{
    int numThreads = -1;
    #pragma omp parallel
    {
        #pragma omp single
        numThreads = omp_get_num_threads();
    }
    return numThreads;
}

TEST(MpiPlusOpenMp, Reduction)
{
    omp_set_num_threads(8);

    double threadSum = 0;
    #pragma omp parallel reduction(+:threadSum)
    {
        threadSum = 1.0;
    }

    double mpiSum = -1;
    const int numItemsPerProc = 1;
    stk::all_reduce_sum(MPI_COMM_WORLD, &threadSum, &mpiSum, numItemsPerProc);

    int numProcs = getNumProcs();
    int numThreads = getNumThreads();
    EXPECT_EQ(numThreads*numProcs, mpiSum);
}

TEST(MpiPlusOpenMp, VectorSumReduction)
{
    int numProcs = getNumProcs();
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
        #pragma omp for reduction(+:threadSum)
        for (size_t i = 0; i < sizeOfLocalVector; i++)
        {
            threadSum += vec[i];
        }
    }

    double mpiSum = -1;
    const int numItemsPerProc = 1;
    stk::all_reduce_sum(MPI_COMM_WORLD, &threadSum, &mpiSum, numItemsPerProc);

    double end_time_wall = stk::wall_time();

    int procId = getProcId();
    if(procId == 0)
    {
        int numThreads = getNumThreads();
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
