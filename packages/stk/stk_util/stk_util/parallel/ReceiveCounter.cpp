#include "stk_util/parallel/ReceiveCounter.hpp"
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequireMsg
#include <cassert>                          // for assert

#include <iostream>  //TODO: DEBUGGING

namespace stk {

void ReceiveCounter::start_receive_count(const std::vector<int> &sendCounts)
{
    int commSize;
    MPI_Comm_size(m_comm, &commSize);
    STK_ThrowRequireMsg(sendCounts.size() == static_cast<size_t>(commSize), "send counts must have same length as MPI Communicator size");
    STK_ThrowRequireMsg(m_recvFinished, "Previous receive count must have completed before starting a new one");

    int commRank;
    MPI_Comm_rank(m_comm, &commRank);

    m_recvFinished = false;
    m_sendCount.resize(sendCounts.size());
    for (unsigned int i=0; i < sendCounts.size(); ++i)
    {
      m_sendCount[i] = sendCounts[i] > 0 ? 1 : 0;
    }

    MPI_Ireduce_scatter_block(m_sendCount.data(), &m_nrecv, 1, MPI_INT, MPI_SUM, m_comm, &m_recvReq);
}

bool ReceiveCounter::is_complete()
{

    if (m_recvFinished) {
        return true;
    } else {
        int isComplete = false;
        MPI_Test(&m_recvReq, &isComplete, MPI_STATUS_IGNORE);
        if (isComplete)
            m_recvFinished = true;

        return isComplete;
    }
}

int ReceiveCounter::get_receive_count()
{
    assert(m_recvFinished);
    return m_nrecv;
}

}
