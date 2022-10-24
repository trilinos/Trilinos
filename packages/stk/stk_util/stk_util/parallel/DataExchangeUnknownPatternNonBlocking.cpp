#include "stk_util/parallel/DataExchangeUnknownPatternNonBlocking.hpp"
#include <thread>
#include "stk_util/parallel/MPITagManager.hpp"

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

    m_tag = get_mpi_tag_manager().get_tag(m_comm, m_tagHint);
}

void DataExchangeUnknownPatternNonBlocking::yield()
{
  // Note: sleep_for would be better for this, but its minimum sleep time is
  // too long
  //std::this_thread::yield();
}

}  // namespace
