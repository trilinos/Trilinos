#include "ReceiveSizeCounter.hpp"

#ifdef STK_HAS_MPI

namespace stk {

bool ReceiveSizeCounter::is_complete()
{
  if (m_numReceived == m_nrecvs)
  {
    return true;
  } else
  {
    int flag = false;
    do {
      MPI_Status status;          
      int idx = 0;
      MPI_Testany(m_nrecvs, m_recvReqs.data(), &idx, &flag, &status);
      if (flag)          
      {
        m_recvCounts[status.MPI_SOURCE] = m_recvCountsCompressed[idx];
        m_numReceived++;
      }          
    } while (flag && m_numReceived < m_nrecvs);   
  }
  
  return false;
}

const std::vector<ReceiveSizeCounter::ULL>& ReceiveSizeCounter::get_receive_sizes() const
{
  assert(m_numReceived == m_nrecvs); 
  return m_recvCounts;
}

void ReceiveSizeCounter::finish_receives()
{
  for (int i=0; i < m_nrecvs; ++i)
  {
    MPI_Status status;
    int idx = 0;
    MPI_Waitany(m_recvReqs.size(), m_recvReqs.data(), &idx, &status);
    m_recvCounts[status.MPI_SOURCE] = m_recvCountsCompressed[idx];
    m_numReceived++;    
  }
}

}

#endif