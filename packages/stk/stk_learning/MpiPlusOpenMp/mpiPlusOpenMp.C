// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

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
