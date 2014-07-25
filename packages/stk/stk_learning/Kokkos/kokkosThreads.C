#include <gtest/gtest.h>
#include <iostream>
#include <vector>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <Kokkos_Threads.hpp>
#include <Kokkos_hwloc.hpp>
#include <Kokkos_View.hpp>

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
