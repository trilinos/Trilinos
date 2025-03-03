#include "gtest/gtest.h"
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>

struct ReducedDependencyGeometricTransferCommTest : public ::testing::Test
{
  using data_t = size_t;

  void SetUp() override{

    auto rank = stk::parallel_machine_rank(MPI_COMM_WORLD);
    auto num_ranks = stk::parallel_machine_size(MPI_COMM_WORLD);

    // e.g. comm pattern 0 -> 1 -> 2 -> 0
    auto dest_rank = (rank + 1) % num_ranks;
    auto source_rank = (rank - 1) < 0 ? num_ranks - 1 : (rank - 1) % num_ranks;

    comm.numToMeshCommunications = 1;
    comm.uniqueToProcVec.push_back(source_rank);
    comm.offset_and_num_keys_to_mesh.emplace_back(0, buffer_size);

    comm.numFromMeshCommunications = 1;
    comm.uniqueFromProcVec.push_back(dest_rank);
    comm.offset_and_num_keys_from_mesh.emplace_back(0, buffer_size);
  }

  stk::transfer::ReducedDependencyCommData comm;
  const size_t max_int = std::numeric_limits<int>::max();
  const int buffer_size = 1 + max_int / sizeof(data_t);
};

TEST_F(ReducedDependencyGeometricTransferCommTest, do_comm_greater_than_2GB)
{
  EXPECT_TRUE(buffer_size * sizeof(data_t) > max_int);

  std::vector<data_t> send(buffer_size);
  for (size_t i = 0; i < send.size(); ++i) {
    send[i] = i;
  }

  std::vector<data_t> recv(buffer_size);

  do_communication(comm, send, recv);

  for (size_t i = 0; i < recv.size(); ++i) {
    EXPECT_TRUE(recv[i] == i);
  }
}



TEST_F(ReducedDependencyGeometricTransferCommTest, do_reverse_comm_greater_than_2GB)
{
  EXPECT_TRUE(buffer_size * sizeof(data_t) > max_int);

  std::vector<data_t> send(buffer_size);
  for (size_t i = 0; i < send.size(); ++i) {
    send[i] = i;
  }

  std::vector<data_t> recv(buffer_size);

  do_reverse_communication(comm, recv, send);

  for (size_t i = 0; i < recv.size(); ++i) {
    EXPECT_TRUE(recv[i] == i);
  }
}


