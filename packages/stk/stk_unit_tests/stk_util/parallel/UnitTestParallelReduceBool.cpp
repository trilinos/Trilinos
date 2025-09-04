#include "gtest/gtest.h"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"

#ifdef STK_HAS_MPI

//-----------------------------------------------------------------------------
// is_true_on_any_proc test

TEST(ParallelReduceBool, is_true_on_any_proc_all_true)
{
  EXPECT_TRUE(stk::is_true_on_any_proc(MPI_COMM_WORLD, true));
}

TEST(ParallelReduceBool, is_true_on_any_proc_all_false)
{
  EXPECT_FALSE(stk::is_true_on_any_proc(MPI_COMM_WORLD, false));
}

TEST(ParallelReduceBool, is_true_on_any_proc_one_true)
{
  int comm_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int comm_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  for (int i=0; i < comm_size; ++i)
  {
    bool val = comm_rank == i ? true : false;
    EXPECT_TRUE(stk::is_true_on_any_proc(MPI_COMM_WORLD, val));
  }
}

//-----------------------------------------------------------------------------
// is_true_on_all_procs tests

TEST(ParallelReduceBool, is_true_on_all_procs_all_true)
{
  EXPECT_TRUE(stk::is_true_on_all_procs(MPI_COMM_WORLD, true));
}

TEST(ParallelReduceBool, is_true_on_all_procs_all_false)
{
  EXPECT_FALSE(stk::is_true_on_all_procs(MPI_COMM_WORLD, false));
}

TEST(ParallelReduceBool, is_true_on_all_procs_one_true)
{
  int comm_rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
  int comm_size = stk::parallel_machine_size(MPI_COMM_WORLD);

  for (int i=0; i < comm_size; ++i)
  {
    bool val = comm_rank == i ? true : false;
    if (comm_size > 1) {
      EXPECT_FALSE(stk::is_true_on_all_procs(MPI_COMM_WORLD, val));
    } else {
      EXPECT_TRUE(stk::is_true_on_all_procs(MPI_COMM_WORLD, val));
    }
  }
}

#endif
