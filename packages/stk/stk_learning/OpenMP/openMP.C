// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#if defined(_OPENMP) && !defined(__INTEL_COMPILER)
// there seems to be an issue with OpenMP combined with GoogleTest macros with Intel compilers
// example of error:
//    openMP.C(206): internal error: assertion failed at: "shared/cfe/edgcpfe/checkdir.c", line 5531
#include <omp.h>

namespace
{

//DocTest1
TEST(OPENMP, HelloWorldDontSetNumThreadsInsideCode)
{
    #pragma omp parallel
    {
        std::cout << "Hello, World!\n";
    }
}

//DocTest2
TEST(OpenMp, HelloWorldSetNumThreadsstart_timeInsideCode)
{
    int numThreads = 8;
    omp_set_num_threads(numThreads);
    #pragma omp parallel
    {
        std::cout << "Hello, World!\n";
    }
}

//DocTest3
TEST(OPENMP, HelloWorldUsingPrivate)
{
    int numThreads = 8;
    omp_set_num_threads(numThreads);
    int threadId = -1;
    #pragma omp parallel private(threadId)
    {
        threadId = omp_get_thread_num();
        std::stringstream avoidScrambledOutputByHavingSingleCoutCommand;
        avoidScrambledOutputByHavingSingleCoutCommand << "Hello, World! (from thread " << threadId << "/" << numThreads << ")" << std::endl;
        std::cout << avoidScrambledOutputByHavingSingleCoutCommand.str();
    }
}

//DocTest4
// #define DO_OUTPUT
TEST(OpenMp, MatrixVectorMultiplyUsingThreads)
{
    int numThreads = 8;
    omp_set_num_threads(numThreads);
//    size_t numRows = 200000;
//    size_t numCols = 20000;
    size_t numRows = 2000;
    size_t numCols = 2000;
    std::vector<std::vector<double> > matrix(numRows);
    std::vector<double> vecIn(numCols, 1);
    std::vector<double> vecOut(numRows, 0);
    for(size_t i = 0; i < numRows; i++)
    {
        matrix[i].resize(numCols, i + 1);
    }

    double startTime = stk::cpu_time();
    double startWallTime = stk::wall_time();

    #pragma omp parallel
    {
        #if defined(DO_OUTPUT)
        int tid = omp_get_thread_num();
        std::ostringstream oss;
        oss << "Thread_" << tid;
        std::string filename = oss.str();
        std::ofstream output(filename.c_str());
        #endif
        numThreads = omp_get_num_threads();
        #pragma omp for
        for(size_t i = 0; i < numRows; i++)
        {
            #if defined(DO_OUTPUT)
            output << "Thread " << tid << " working on row " << i << std::endl;
            #endif
            for(size_t j = 0; j < numCols; j++)
            {
                vecOut[i] += matrix[i][j] * vecIn[j];
            }
        }
        #if defined(DO_OUTPUT)
        output.close();
        #endif
    }

    double endTime = stk::cpu_time();
    double endWallTime = stk::wall_time();

    std::cerr << "Num threads: " << numThreads << std::endl;
    std::cerr << "Time: " << (endTime - startTime)/static_cast<double>(numThreads) << std::endl;
    std::cerr << "Wall Time: " << (endWallTime - startWallTime) << std::endl;

    for(size_t i = 0; i < numRows; i++)
    {
        double goldAnswer = (i + 1) * numCols;
        EXPECT_EQ(goldAnswer, vecOut[i]);
    }
}

//DocTest5
TEST(OpenMp, SumOverVector)
{
    int numThreads = 8;
    omp_set_num_threads(numThreads);
    //size_t sizeOfVector = 4000000000;
    size_t sizeOfVector = 1000000;
    double initVal = 1.0;
    std::vector<double> vec(sizeOfVector, initVal);

    double sum = 0;

    double startTime = stk::cpu_time();
    double startWallTime = stk::wall_time();

    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
        #pragma omp for reduction(+:sum)
        for (size_t i = 0; i < sizeOfVector; i++)
        {
            sum += vec[i];
        }
    }

    double endTime = stk::cpu_time();
    double endWallTime = stk::wall_time();

    std::cerr << "Num threads: " << numThreads << std::endl;
    std::cerr << "Time: " << (endTime - startTime)/static_cast<double>(numThreads) << std::endl;
    std::cerr << "Wall Time: " << (endWallTime - startWallTime) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}

//DocTest6
TEST(OpenMp, SumUsingSections)
{
    int numThreads = 2;
    omp_set_num_threads(numThreads);

    size_t sizeOfVector = 100000;
    double initVal = 1.0;
    std::vector<double> vec(sizeOfVector,initVal);

    double sum1 = 0, sum2 = 0;

    double startTime = stk::cpu_time();
    double startWallTime = stk::wall_time();

    size_t halfLoopSize = sizeOfVector/2;

    #pragma omp parallel
    {
        numThreads = omp_get_num_threads();
        #pragma omp sections
        {
            #pragma omp section
            {
                std::cerr << "Loop 1: hello from thread: " << omp_get_thread_num() << std::endl;
                for (size_t i = 0; i < halfLoopSize; i++)
                {
                    sum1 += vec[i];
                }
            }
            #pragma omp section
            {
                std::cerr << "Loop 2: hello from thread: " << omp_get_thread_num() << std::endl;
                for (size_t i = halfLoopSize; i < sizeOfVector; i++)
                {
                    sum2 += vec[i];
                }
            }
        }
    }

    double sum = sum1 + sum2;

    double endTime = stk::cpu_time();
    double endWallTime = stk::wall_time();

    std::cerr << "Num threads: " << numThreads << std::endl;
    std::cerr << "Time: " << (endTime - startTime)/static_cast<double>(numThreads) << std::endl;
    std::cerr << "Wall Time: " << (endWallTime - startWallTime) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}
//EndDocTest

//DocTestForLearningAboutPrivates
struct SimpleDefaultedInt
{
    SimpleDefaultedInt() : value(-1) {}
    SimpleDefaultedInt(int defaultValue) : value(defaultValue) {}
    int value;
};

TEST(OpenMp, learningAboutPrivates)
{
    SimpleDefaultedInt a(13);
    SimpleDefaultedInt b(14);
    SimpleDefaultedInt c(15);

    const int numberOfIterations = 10;
    #pragma omp parallel private(a) firstprivate(b) shared(c)
    {
        EXPECT_NE(13, a.value);
        EXPECT_EQ(14, b.value);
        EXPECT_EQ(15, c.value);
        #pragma omp barrier
        #pragma omp for lastprivate(c)
        for(int i =0; i <= numberOfIterations; i++)
        {
            EXPECT_NE(15, c.value);  // with excruciatingly high probability

            a.value = i;
            b.value = i;
            c.value = i;
        }
    }
    EXPECT_EQ(13, a.value);
    EXPECT_EQ(14, b.value);
    EXPECT_EQ(numberOfIterations, c.value);
}
//EndDocTest

} // end namespace
#endif
