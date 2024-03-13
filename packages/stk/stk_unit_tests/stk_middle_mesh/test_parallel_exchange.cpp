#include "stk_middle_mesh/parallel_exchange.hpp"
#include "gtest/gtest.h"

namespace stk {
namespace middle_mesh {
namespace impl {

namespace {

using utils::impl::comm_rank;
using utils::impl::comm_size;

int get_value(int sendRank, int recvRank, int idx)
{
  return sendRank + 2 * recvRank + idx;
}
} // namespace

TEST(ParallelExchange, AllToAll)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  utils::impl::ParallelExchange<int> exchanger(comm, 42);

  int myrank = comm_rank(comm);
  int nvals  = 100;
  for (int destRank = 0; destRank < comm_size(comm); ++destRank)
  {
    exchanger.set_recv_buffer_size(destRank, nvals);

    auto& buf = exchanger.get_send_buffer(destRank);
    buf.resize(nvals);
    for (int i = 0; i < nvals; ++i)
      buf[i] = get_value(myrank, destRank, i);

    exchanger.start_send(destRank);
    exchanger.start_recv(destRank);
  }

  exchanger.complete_recvs();

  for (int sendRank = 0; sendRank < comm_size(comm); ++sendRank)
  {
    const auto& recvBuf = exchanger.get_recv_buffer(sendRank);
    EXPECT_EQ(recvBuf.size(), size_t(nvals));
    for (int i = 0; i < nvals; ++i)
      EXPECT_EQ(recvBuf[i], get_value(sendRank, myrank, i));
  }

  exchanger.complete_sends();
  MPI_Barrier(comm);
}

TEST(ParallelExchange, UpwardNeighbor)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  utils::impl::ParallelExchange<int> exchanger(comm, 42);

  int myrank   = comm_rank(comm);
  int destRank = (myrank + 1) % comm_size(comm);
  int srcRank  = (myrank - 1 + comm_size(comm)) % comm_size(comm);

  int nvals = 100;
  exchanger.set_recv_buffer_size(srcRank, nvals);

  auto& buf = exchanger.get_send_buffer(destRank);
  buf.resize(nvals);
  for (int i = 0; i < nvals; ++i)
    buf[i] = get_value(myrank, destRank, i);

  exchanger.start_send(destRank);
  exchanger.start_recv(srcRank);

  exchanger.complete_recvs();

  const auto& recvBuf = exchanger.get_recv_buffer(srcRank);
  EXPECT_EQ(recvBuf.size(), size_t(nvals));
  for (int i = 0; i < nvals; ++i)
    EXPECT_EQ(recvBuf[i], get_value(srcRank, myrank, i));

  exchanger.complete_sends();
  MPI_Barrier(comm);
}

} // namespace impl
} // namespace middle_mesh
} // namespace stk
