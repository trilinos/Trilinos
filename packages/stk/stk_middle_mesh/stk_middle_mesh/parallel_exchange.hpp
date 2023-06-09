#ifndef PARALLEL_EXCHANGE_H
#define PARALLEL_EXCHANGE_H

#include "mpi.h"
#include "utils.hpp"
#include <stdexcept>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

template <typename T>
class ParallelExchange
{
  public:
    ParallelExchange(MPI_Comm comm, int tag)
      : m_comm(comm)
      , m_tag(tag)
      , m_sendBufs(comm_size(comm))
      , m_recvBufs(comm_size(comm))
      , m_sendReqs(comm_size(comm), MPI_REQUEST_NULL)
      , m_recvReqs(comm_size(comm), MPI_REQUEST_NULL)
    {}

    ~ParallelExchange() { complete_sends(); }

    std::vector<T>& get_send_buffer(int rank) { return m_sendBufs[rank]; }

    void set_recv_buffer_size(int rank, int numel) { m_recvBufs[rank].resize(numel); }

    const std::vector<T>& get_recv_buffer(int rank) { return m_recvBufs[rank]; }

    void start_send(int rank)
    {
      if (m_sendReqs[rank] != MPI_REQUEST_NULL)
        throw std::runtime_error("previous send not completed, cannot start a new one");

      MPI_Isend(m_sendBufs[rank].data(), m_sendBufs[rank].size() * sizeof(T), MPI_BYTE, rank, m_tag, m_comm,
                &(m_sendReqs[rank]));
    }

    void start_sends()
    {
      for (std::size_t i = 0; i < m_sendBufs.size(); ++i)
        if (m_sendBufs[i].size() > 0)
          start_send(i);
    }

    void start_recv(int rank)
    {
      if (m_recvReqs[rank] != MPI_REQUEST_NULL)
        throw std::runtime_error("previous recv not completed, cannot start a new one");

      MPI_Irecv(m_recvBufs[rank].data(), m_recvBufs[rank].size() * sizeof(T), MPI_BYTE, rank, m_tag, m_comm,
                &(m_recvReqs[rank]));
    }

    void start_recvs()
    {
      for (std::size_t i = 0; i < m_recvBufs.size(); ++i)
        if (m_recvBufs[i].size() > 0)
          start_recv(i);
    }

    void complete_sends() { MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE); }

    void complete_recvs() { MPI_Waitall(m_recvReqs.size(), m_recvReqs.data(), MPI_STATUSES_IGNORE); }

  private:
    MPI_Comm m_comm;
    int m_tag;
    std::vector<std::vector<T>> m_sendBufs;
    std::vector<std::vector<T>> m_recvBufs;
    std::vector<MPI_Request> m_sendReqs;
    std::vector<MPI_Request> m_recvReqs;
};

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif