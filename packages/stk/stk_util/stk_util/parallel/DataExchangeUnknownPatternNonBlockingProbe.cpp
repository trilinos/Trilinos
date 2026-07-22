#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingProbe.hpp"
#include <thread>
#include "stk_util/parallel/MPITagManager.hpp"

namespace stk {


void DataExchangeUnknownPatternNonBlockingProbe::reset()
{
    // check that communication finished from previous iteration
    STK_ThrowRequireMsg(!m_areRecvsInProgress, "Previous receive must have completed before starting a new one");
    STK_ThrowRequireMsg(!m_areSendsInProgress, "Previous send must have completed before starting a new one");

    // setup for new iteration
    m_recvReqRanks.resize(0);
    m_extraRecvBufs.resize(0);
    for (std::vector<int> recvBufsIdxs : m_rankExtraRecvBufs)
    {
      recvBufsIdxs.clear();
    }
    m_sendReqs.resize(0);
    m_recvReqs.resize(0);
    m_recvcount = 0;

    m_tag = get_mpi_tag_manager().get_tag(m_comm, m_tagHint);
}

void DataExchangeUnknownPatternNonBlockingProbe::yield()
{
  // Note: sleep_for would be better for this, but its minimum sleep time is
  // too long
  //std::this_thread::yield();
}

}  // namespace
