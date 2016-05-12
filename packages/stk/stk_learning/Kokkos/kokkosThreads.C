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

#include "mtk_kokkos.h"

#include <Kokkos_Core.hpp>

struct MatVecFunctor
{
    MatVecFunctor(Kokkos::View<double**> matrix, Kokkos::View<double*> vecIn, Kokkos::View<double*> vecOut) :
            matrix(matrix),
            vecIn(vecIn),
            vecOut(vecOut)
    {
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t row) const
    {
        const size_t numCols = matrix.extent(1);
        for(size_t col = 0; col < numCols; col++)
            vecOut(row) += matrix(row, col) * vecIn(col);
    }

private:
    Kokkos::View<double**> matrix;
    Kokkos::View<double*> vecIn;
    Kokkos::View<double*> vecOut;
};

void runMatVecThreadedTest()
{
    const size_t numRows = 3;
    const size_t numCols = 2;
    Kokkos::View<double**> matrix("matrix", numRows, numCols);
    Kokkos::View<double*> vecIn("vecIn", numCols);
    Kokkos::View<double*> vecOut("vecOut", numRows);

    Kokkos::deep_copy(vecIn, 1);
    Kokkos::deep_copy(vecOut, 0);

    Kokkos::parallel_for(numRows,
        KOKKOS_LAMBDA(size_t row)
        {
            for(size_t col=0; col<numCols; col++)
                matrix(row, col) = row + 1;
        }
    );

    double start_time_wall = stk::wall_time();

    MatVecFunctor matVecFunctor(matrix, vecIn, vecOut);
    Kokkos::parallel_for(numRows, matVecFunctor);

    double end_time_wall = stk::wall_time();

    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    Kokkos::View<double*>::HostMirror vecOutHost =  Kokkos::create_mirror_view(vecOut);
    Kokkos::deep_copy(vecOutHost, vecOut);
    for(size_t i = 0; i < numRows; i++)
    {
        double goldAnswer = (i + 1) * numCols;
        EXPECT_EQ(goldAnswer, vecOutHost(i));
    }
}

TEST_F(MTK_Kokkos, MatrixVectorMultiply)
{
    runMatVecThreadedTest();
}




struct SumOverVectorFunctor
{
    SumOverVectorFunctor(Kokkos::View<double*> vector) :
        mVector(vector)
    {}

    KOKKOS_INLINE_FUNCTION
    void operator()(size_t vectorIndex, double &threadLocalSum) const
    {
        threadLocalSum += mVector(vectorIndex);
    }

private:
    Kokkos::View<double*> mVector;
};

void runSumOverVectorTest()
{
    //size_t sizeOfVector = 4000000000;
    size_t sizeOfVector = 100000;
    Kokkos::View<double*> vec("vec", sizeOfVector);

    double initVal = 1.0;
    Kokkos::deep_copy(vec, initVal);

    double start_time_wall = stk::wall_time();

    double sum = 0;
    SumOverVectorFunctor sumOverVectorFunctor(vec);
    Kokkos::parallel_reduce(sizeOfVector, sumOverVectorFunctor, sum);

    double end_time_wall = stk::wall_time();

    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}

TEST_F(MTK_Kokkos, SumOverVector)
{
    runSumOverVectorTest();
}

