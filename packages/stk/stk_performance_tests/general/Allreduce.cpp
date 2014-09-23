#include <gtest/gtest.h>

#include <vector>

#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/WallTime.hpp>

TEST( Parallel, allreduce)
{
  stk::ParallelMachine comm = MPI_COMM_WORLD;

  int num_procs = stk::parallel_machine_size(comm);
  int local_proc = stk::parallel_machine_rank(comm);

  stk::parallel_machine_barrier(comm);

  double start_time = stk::wall_time();

  std::vector<double> local_data(num_procs, 0);
  std::vector<double> global_data(num_procs, 0);

  int num_iters = 100;

  for(int i=0; i<num_iters; ++i) {
    stk::all_reduce_max(comm, &local_data[0], &global_data[0], local_data.size());
  }

  double elapsed_time = stk::wall_time() - start_time;

  if (local_proc == 0) {
    std::cout << "Allreduce time (proc 0 of "<<num_procs<<"): " << elapsed_time << " seconds." << std::endl;
  }
}
