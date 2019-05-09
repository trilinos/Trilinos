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

namespace stk
{

#define STK_LAMBDA KOKKOS_LAMBDA
#define STK_INLINE KOKKOS_INLINE_FUNCTION

#ifdef KOKKOS_ENABLE_CUDA
using DeviceSpace = Kokkos::CudaSpace;
#else
using DeviceSpace = Kokkos::HostSpace;
#endif

class NgpVector
{
public:
    NgpVector(size_t size) :
        vector("NgpVector", size),
        hostVector(Kokkos::create_mirror_view(vector))
    {
    }
    NgpVector(size_t size, double initialValue) :
        NgpVector(size)
    {
        Kokkos::deep_copy(vector, initialValue);
    }

    STK_INLINE double & device_get(size_t i) const
    {
        return vector(i);
    }
    double & host_get(size_t i) const
    {
        return hostVector(i);
    }
    void copy_device_to_host()
    {
        Kokkos::deep_copy(hostVector, vector);
    }
private:
  typedef Kokkos::View<double *, DeviceSpace> KokkosVector;
    KokkosVector vector;
    KokkosVector::HostMirror hostVector;
};

class NgpMatrix
{
public:
    NgpMatrix(size_t numRows, size_t numCols) :
        matrix("NgpMatrix", numRows, numCols),
        hostMatrix(Kokkos::create_mirror_view(matrix))
    {
    }
    NgpMatrix(size_t numRows, size_t numCols, double initialValue) :
        NgpMatrix(numRows, numCols)
    {
        Kokkos::deep_copy(matrix, initialValue);
    }

    STK_INLINE double & device_get(size_t row, size_t col) const
    {
        return matrix(row, col);
    }
    double & host_get(size_t row, size_t col) const
    {
        return hostMatrix(row, col);
    }

    STK_INLINE size_t num_rows() const
    {
        return matrix.extent(0);
    }
    STK_INLINE size_t num_cols() const
    {
        return matrix.extent(1);
    }

    void copy_device_to_host()
    {
        Kokkos::deep_copy(hostMatrix, matrix);
    }
private:
  typedef Kokkos::View<double **, DeviceSpace> KokkosMatrix;
    KokkosMatrix matrix;
    KokkosMatrix::HostMirror hostMatrix;
};

template <typename FUNCTOR>
void parallel_for(size_t n, const FUNCTOR &functor)
{
    Kokkos::parallel_for(n, functor);
}
template <typename FUNCTOR, typename VAL_TYPE>
void parallel_reduce(size_t n, const FUNCTOR &functor, VAL_TYPE &sum)
{
    Kokkos::parallel_reduce(n, functor, sum);
}

} //namespace stk


namespace NgpMatrixComputations
{
    void initialize_matrix_to_row_values(stk::NgpMatrix matrix)
    {
        stk::parallel_for(matrix.num_rows(),
            STK_LAMBDA(size_t row)
            {
                const size_t numCols = matrix.num_cols();
                for(size_t col=0; col<numCols; col++)
                    matrix.device_get(row, col) = row + 1;
            }
        );
    }
    void compute_mat_vec(stk::NgpMatrix matrix, stk::NgpVector vecIn, stk::NgpVector vecOut)
    {
        stk::parallel_for(matrix.num_rows(),
            STK_LAMBDA(size_t row)
            {
                const size_t numCols = matrix.num_cols();
                for(size_t col = 0; col < numCols; col++)
                    vecOut.device_get(row) += matrix.device_get(row, col) * vecIn.device_get(col);
            }
        );
    }
}


void report_bandwidth(double numDoubles, double time)
{
    double matVecSizeGb = numDoubles * sizeof(double) / (1<<30);
    std::cerr << "\t\t\t\t\tBandwidth: " << matVecSizeGb / time << " GB/s" << std::endl;
}

class EncapsulateKokkos : public MTK_Kokkos {};
TEST_F(EncapsulateKokkos, MatrixVectorMultiply)
{
    const size_t numRows = 3000;
    const size_t numCols = 200;
    stk::NgpMatrix matrix(numRows, numCols);
    stk::NgpVector vecIn(numCols, 1);
    stk::NgpVector vecOut(numRows);

    NgpMatrixComputations::initialize_matrix_to_row_values(matrix);

    double startTime = stk::wall_time();
    NgpMatrixComputations::compute_mat_vec(matrix, vecIn, vecOut);
    double endTime = stk::wall_time();

    report_bandwidth(numRows*numCols + numRows + numCols, endTime - startTime);

    vecOut.copy_device_to_host();
    for(size_t i = 0; i < numRows; i++)
    {
        double goldAnswer = (i + 1) * numCols;
        EXPECT_EQ(goldAnswer, vecOut.host_get(i));
    }
}




struct SumOverVectorFunctor
{
    SumOverVectorFunctor(stk::NgpVector vector) :
        mVector(vector)
    {}

    STK_INLINE
    void operator()(size_t vectorIndex, double &threadLocalSum) const
    {
        threadLocalSum += mVector.device_get(vectorIndex);
    }

private:
    stk::NgpVector mVector;
};
TEST_F(EncapsulateKokkos, SumOverVector)
{
    size_t sizeOfVector = 1000000;
    double initVal = 1.0;
    stk::NgpVector vec(sizeOfVector, initVal);

    double startTime = stk::wall_time();
    double sum = 0;
    stk::parallel_reduce(sizeOfVector, SumOverVectorFunctor(vec), sum);
    double endTime = stk::wall_time();

    report_bandwidth(sizeOfVector, endTime - startTime);

    EXPECT_EQ(sizeOfVector*initVal, sum);
}


class CopyHostToDevice : public MTK_Kokkos {};
TEST_F(CopyHostToDevice, CopyDataToDeviceAndMirrorBack)
{
    size_t sizeOfVector = 1000000;
    double initVal = 1.0;
    Kokkos::View<double*, Kokkos::HostSpace> hostVec("hostVec", sizeOfVector);
    for(size_t i=0; i<sizeOfVector; i++)
        hostVec(i) = initVal;

    Kokkos::View<double*> vec("vec", sizeOfVector);
    Kokkos::deep_copy(vec, hostVec);

    Kokkos::View<double*>::HostMirror hostMirrorVector = Kokkos::create_mirror_view(vec);
    Kokkos::deep_copy(hostMirrorVector, vec);
    for(size_t i=0; i<sizeOfVector; i++)
        EXPECT_EQ(initVal, hostMirrorVector(i));
}

