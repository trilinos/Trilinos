#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include <thread>

namespace stk {


void DataExchangeUnknownPatternNonBlocking::reset()
{
    // check that communication finished from previous iteration
    ThrowRequireMsg(!m_areRecvsInProgress, "Previous receive must have completed before starting a new one");
    ThrowRequireMsg(!m_areSendsInProgress, "Previous send must have completed before starting a new one");

    // setup for new iteration
    m_sendRankMap.resize(0);
    m_recvRankMap.resize(0);
    m_sendReqs.resize(0);
    m_recvReqs.resize(0);
    m_recvcount = 0;

    MPI_Comm_rank(m_comm, &m_myrank);
    m_tag = m_tag == m_tag1 ? m_tag2 : m_tag1;
}

void DataExchangeUnknownPatternNonBlocking::yield()
{
  // Note: sleep_for would be better for this, but its minimum sleep time is
  // too long
  std::this_thread::yield();
}

}  // namespace
