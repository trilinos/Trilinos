#include "gtest/gtest.h"
#include "stk_util/parallel/IallreduceReplacement.hpp"
#include "stk_util/parallel/Parallel.hpp"

namespace 
{

void sum_op(void* invec, void* inoutvec, int* len, MPI_Datatype* /*datatype*/)
{
  int* invec_int = reinterpret_cast<int*>(invec);
  int* inoutvec_int = reinterpret_cast<int*>(inoutvec);
  for (int i=0; i < *len; ++i)
    inoutvec_int[i] = invec_int[i] + inoutvec_int[i];
}

double test_allreduce_mpi(int numValues)
{
  int num_iters = 100;
  double t_start = 0;

  std::vector<int> input_data(numValues, 1), output_data(numValues, 42);
  MPI_Op mpi_op;
  MPI_Op_create(&sum_op, true, &mpi_op);
  for (int i=0; i < num_iters; ++i)
  {
    if (i == 1)
      t_start = MPI_Wtime();

    MPI_Request req;
    //MPI_Ibarrier(MPI_COMM_WORLD, &req);
    MPI_Iallreduce(input_data.data(), output_data.data(), numValues, MPI_INT, mpi_op, MPI_COMM_WORLD, &req);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
  }

  double elapsed_time = (MPI_Wtime() - t_start) / (num_iters - 1);
  double max_time = 0;
  MPI_Allreduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  return max_time;
}


void test_allreduce(stk::impl::IAllreduceReplacement<int>& allreduce, int numValues)
{
  auto op = [](const int& lhs, const int& rhs) { return lhs + rhs; };
  std::vector<int> input_data(numValues, stk::parallel_machine_rank(MPI_COMM_WORLD)), output_data(numValues, 42),
                   expected_data(numValues);

  int commRank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int commSize = stk::parallel_machine_size(MPI_COMM_WORLD);

  if (commRank == 0)
    std::cout << "\nTesting size " << numValues << std::endl;


  for (int i=0; i < numValues; ++i)
  {
    input_data[i] = commRank + i;
    expected_data[i] = i*commSize + commSize * (commSize - 1) / 2;
  }

  int num_iters = 100;
  double t_start = 0;

  for (int i=0; i < num_iters; ++i)
  {
    if (i == 1)
      t_start = MPI_Wtime();

    allreduce.startReduction(op, input_data, output_data);
    allreduce.finishReduction();
  }

  double elapsed_time = (MPI_Wtime() - t_start) / (num_iters - 1);
  double max_time = 0;
  MPI_Allreduce(&elapsed_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (commRank == 0)
    std::cout << "New algorithm time = " << max_time << " per iteration" << std::endl;

  for (int i=0; i < numValues; ++i)
  {
    EXPECT_EQ(output_data[i], expected_data[i]);
  }

  double max_time2 = test_allreduce_mpi(numValues);

  if (commRank == 0)
    std::cout << "MPI algorithm time = " << max_time2 << " per iteration" << ", ratio " << max_time / max_time2 << std::endl;
}
}




TEST(IallreduceReplacement, ReductionFactor2)
{
  stk::impl::IAllreduceReplacement<int> allreduce(MPI_COMM_WORLD, MPI_INT, 666, 2);
  MPI_Barrier(MPI_COMM_WORLD);
  test_allreduce(allreduce, 1);
  test_allreduce(allreduce, 1);
  test_allreduce(allreduce, 2);
  test_allreduce(allreduce, 10);
  test_allreduce(allreduce, 100);
  test_allreduce(allreduce, 1000);
  test_allreduce(allreduce, 10000);
}

TEST(IallreduceReplacement, ReductionFactor5)
{
  stk::impl::IAllreduceReplacement<int> allreduce(MPI_COMM_WORLD, MPI_INT, 666, 5);
  MPI_Barrier(MPI_COMM_WORLD);
  test_allreduce(allreduce, 1);
  test_allreduce(allreduce, 1);
  test_allreduce(allreduce, 2);
  test_allreduce(allreduce, 10);
  test_allreduce(allreduce, 100);
  test_allreduce(allreduce, 1000);
  test_allreduce(allreduce, 10000);
}

TEST(IallreduceReplacement, ReductionFactorError)
{
  EXPECT_ANY_THROW(stk::impl::IAllreduceReplacement<int> allreduce(MPI_COMM_WORLD, MPI_INT, 666, 1));
}

TEST(IallreduceReplacement, InputSizeError)
{
  auto op = [](const int& lhs, const int& rhs) { return lhs + rhs; };
  stk::impl::IAllreduceReplacement<int> allreduce(MPI_COMM_WORLD, MPI_INT, 666, 2);
  std::vector<int> input_data(2), output_data(3);
  EXPECT_ANY_THROW(allreduce.startReduction(op, input_data, output_data));
}

TEST(IallreduceReplacement, ReuseBeforeCompletedError)
{
  if (stk::parallel_machine_size(MPI_COMM_WORLD) < 4)
  {
    GTEST_SKIP();
  }

  auto op = [](const int& lhs, const int& rhs) { return lhs + rhs; };
  stk::impl::IAllreduceReplacement<int> allreduce(MPI_COMM_WORLD, MPI_INT, 666, 2); 
  std::vector<int> input_data(2, 0), output_data(2, 0);
  allreduce.startReduction(op, input_data, output_data);

  EXPECT_ANY_THROW(allreduce.startReduction(op, input_data, output_data));

  allreduce.finishReduction();
}