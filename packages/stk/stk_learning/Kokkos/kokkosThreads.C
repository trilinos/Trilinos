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

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <Kokkos_Core.hpp>

namespace
{

struct MatVecFunctor
{
    typedef typename Kokkos::Threads device_type;

    MatVecFunctor(const std::vector<std::vector<double> > &matrix,
                  const std::vector<double> &vecIn,
                  std::vector<double> &vecOut) :
            mMatrix(matrix),
            mVecIn(vecIn),
            mVecOut(vecOut)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t rowIndex) const
    {
        const size_t numCols = mMatrix[rowIndex].size();
        for(size_t j = 0; j < numCols; j++)
        {
            mVecOut[rowIndex] += mMatrix[rowIndex][j] * mVecIn[j];
        }
    }

private:
    const std::vector<std::vector<double> > &mMatrix;
    const std::vector<double> &mVecIn;
    std::vector<double> &mVecOut;
};

void runMatVecThreadedTest(const int numThreads)
{
    size_t numRows = 2000;
    size_t numCols = 2000;
    std::vector<std::vector<double> > matrix(numRows);
    std::vector<double> vec_in(numCols, 1);
    std::vector<double> vec_out(numRows, 0);
    for(size_t i = 0; i < numRows; i++)
    {
        matrix[i].resize(numCols, i + 1);
    }

    double start_time = stk::cpu_time();
    double start_time_wall = stk::wall_time();

    Kokkos::Threads::initialize( numThreads );
    MatVecFunctor matVecFunctor(matrix, vec_in, vec_out);
    Kokkos::parallel_for(numRows, matVecFunctor);
    Kokkos::Threads::finalize();

    double end_time = stk::cpu_time();
    double end_time_wall = stk::wall_time();

    std::cerr << "Num threads: " << numThreads << std::endl;
    std::cerr << "Time: " << (end_time - start_time)/double(numThreads) << std::endl;
    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    for(size_t i = 0; i < numRows; i++)
    {
        double goldAnswer = (i + 1) * numCols;
        EXPECT_EQ(goldAnswer, vec_out[i]);
    }
}

TEST(KokkosThreads, MatrixVectorMultiplyUsingOneThread)
{
    const int numThreads = 1;
    runMatVecThreadedTest(numThreads);
}

TEST(KokkosThreads, MatrixVectorMultiplyUsingEightThreads)
{
    const int numThreads = 8;
    runMatVecThreadedTest(numThreads);
}




struct SumOverVectorFunctor
{
    typedef typename Kokkos::Threads device_type;
    typedef double value_type;

    SumOverVectorFunctor(const std::vector<value_type> &vector) :
        mVector(vector)
    {}

    KOKKOS_INLINE_FUNCTION
    void init(value_type &threadLocalSum) const
    {
        threadLocalSum = 0.0;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t vectorIndex, value_type &threadLocalSum) const
    {
        threadLocalSum += mVector[vectorIndex];
    }

    KOKKOS_INLINE_FUNCTION
    void join(volatile value_type & joinedThreadSum, volatile value_type const& threadLocalSum) const
    {
        joinedThreadSum += threadLocalSum;
    }

private:
    const std::vector<value_type> mVector;
};

void runSumOverVectorTest(const int numThreads)
{
    //size_t sizeOfVector = 4000000000;
    size_t sizeOfVector = 1000000;
    double initVal = 1.0;
    std::vector<double> vec(sizeOfVector,initVal);

    double start_time = stk::cpu_time();
    double start_time_wall = stk::wall_time();

    double sum = 0;

    Kokkos::Threads::initialize( numThreads );
    SumOverVectorFunctor sumOverVectorFunctor(vec);
    Kokkos::parallel_reduce(sizeOfVector, sumOverVectorFunctor, sum);
    Kokkos::Threads::finalize();

    double end_time = stk::cpu_time();
    double end_time_wall = stk::wall_time();

    std::cerr << "Num threads: " << numThreads << std::endl;
    std::cerr << "Time: " << (end_time - start_time)/double(numThreads) << std::endl;
    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}

TEST(KokkosThreads, SumOverVectorUsingOneThread)
{
    const int numThreads = 1;
    runSumOverVectorTest(numThreads);
}

TEST(KokkosThreads, SumOverVectorUsingEightThreads)
{
    const int numThreads = 8;
    runSumOverVectorTest(numThreads);
}

} // end namespace
