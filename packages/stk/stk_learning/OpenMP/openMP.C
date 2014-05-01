#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <stk_util/environment/WallTime.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <omp.h>

namespace
{

//DocTest1
TEST(OPENMP, HelloWorldDontSetNumThreadsInsideCode)
{
#pragma omp parallel
    {
        std::cout << "Hello World\n";
    }
}

//DocTest2
TEST(OPENMP, HelloWorldSetNumThreadsstart_timeInsideCode)
{
    int numThreads = 8;
    omp_set_num_threads(numThreads);
#pragma omp parallel
    {
        std::cout << "Hello World\n";
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
        std::cout << "Hello World from thread " << threadId << "/" << numThreads << std::endl;
    }
}

//DocTest4
// #define DO_OUTPUT

TEST(OPENMP, MatrixVectorMultiplyUsingThreads)
{
//    size_t numRows = 200000;
//    size_t numCols = 20000;
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

    size_t i = 0, j = 0;
    int num_threads = 0;
    int tid = -1;

#pragma omp parallel private (tid)
    {
        tid = omp_get_thread_num();
#if defined(DO_OUTPUT)
        std::ostringstream oss;
        oss << "Thread_" << tid;
        std::string filename = oss.str();
        std::ofstream output(filename.c_str());
#endif
        num_threads = omp_get_num_threads();
        #pragma omp for
        for(i = 0; i < numRows; i++)
        {
#if defined(DO_OUTPUT)
            output << "Thread " << tid << " working on row " << i << std::endl;
#endif
            for(j = 0; j < numCols; j++)
            {
                vec_out[i] += matrix[i][j] * vec_in[j];
            }
        }
#if defined(DO_OUTPUT)
        output.close();
#endif
    }

    double end_time = stk::cpu_time();
    double end_time_wall = stk::wall_time();

    std::cerr << "Num threads: " << num_threads << std::endl;
    std::cerr << "Time: " << (end_time - start_time)/double(num_threads) << std::endl;
    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    for(size_t i = 0; i < numRows; i++)
    {
        double goldAnswer = (i + 1) * numCols;
        EXPECT_EQ(goldAnswer, vec_out[i]);
    }
}

//DocTest5
TEST(OPENMP, SumOverVector)
{
    //size_t sizeOfVector = 4000000000;
    size_t sizeOfVector = 1000000;
    double initVal = 1.0;
    std::vector<double> vec(sizeOfVector,initVal);
    int num_threads = 0;

    double sum = 0;

    double start_time = stk::cpu_time();
    double start_time_wall = stk::wall_time();

    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp for reduction(+:sum)
        for (size_t i=0;i<sizeOfVector;i++)
        {
            sum += vec[i];
        }
    }

    double end_time = stk::cpu_time();
    double end_time_wall = stk::wall_time();

    std::cerr << "Num threads: " << num_threads << std::endl;
    std::cerr << "Time: " << (end_time - start_time)/double(num_threads) << std::endl;
    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}

//DocTest6
TEST(OPENMP, SumUsingSections)
{
    int numThreads = 2;
    omp_set_num_threads(numThreads);

    size_t sizeOfVector = 100000;
    double initVal = 1.0;
    std::vector<double> vec(sizeOfVector,initVal);
    int num_threads = 0;

    double sum1 = 0, sum2 = 0;

    double start_time = stk::cpu_time();
    double start_time_wall = stk::wall_time();

    size_t halfLoopSize = sizeOfVector/2;

    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
        #pragma omp sections
        {
            #pragma omp section
            {
                std::cerr << "Loop 1: hello from thread: " << omp_get_thread_num() << std::endl;
                for (size_t i=0;i<halfLoopSize;i++)
                {
                    sum1 += vec[i];
                }
            }
            #pragma omp section
            {
                std::cerr << "Loop 2: hello from thread: " << omp_get_thread_num() << std::endl;
                for (size_t i=halfLoopSize;i<sizeOfVector;i++)
                {
                    sum2 += vec[i];
                }
            }
        }
    }

    double sum = sum1+sum2;

    double end_time = stk::cpu_time();
    double end_time_wall = stk::wall_time();

    std::cerr << "Num threads: " << num_threads << std::endl;
    std::cerr << "Time: " << (end_time - start_time)/double(num_threads) << std::endl;
    std::cerr << "Wall Time: " << (end_time_wall - start_time_wall) << std::endl;

    double goldAnswer = sizeOfVector*initVal;
    EXPECT_EQ(goldAnswer, sum);
}
//EndDocTest

} // end namespace
