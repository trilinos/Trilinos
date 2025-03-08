#include "DataExchangeKnownPatternUserDataNonBlocking.hpp"

namespace stk {

DataExchangeKnownPatternUserDataNonBlocking::~DataExchangeKnownPatternUserDataNonBlocking()
{
  if (m_areSendsInProgress) {
    complete_sends();
  }
}

void DataExchangeKnownPatternUserDataNonBlocking::start_nonblocking(const std::vector<PointerAndSize>& sendData, const std::vector<int>& sendRanks,
                        std::vector<PointerAndSize>& recvData, const std::vector<int>& recvRanks)
{
  start_recvs(recvData, recvRanks);
  start_sends(sendData, sendRanks);
}

void DataExchangeKnownPatternUserDataNonBlocking::start_recvs(std::vector<PointerAndSize> &recvData, const std::vector<int>& recvRanks)
{
  STK_ThrowRequireMsg(recvData.size() == recvRanks.size(), "recvData and recvRanks must have same length");
  STK_ThrowRequireMsg(!m_areRecvsInProgress,
                  "cannot start new round of communication until the recvs from the previous round are complete");

  m_recvIndexMap.resize(0);
  m_recvReqs.resize(0);
  m_recvCounts.resize(0);
  for (size_t i=0; i < recvRanks.size(); ++i)
  {
    int rank = recvRanks[i];
    PointerAndSize recvBuf = recvData[i];
    int numSegments = computeNumSegments(recvBuf.size, MAX_MESSAGE_SIZE);
    size_t bytesReceived = 0;
    m_recvCounts.emplace_back();
    
    for (int j=0; j < numSegments; ++j)
    {
      m_recvIndexMap.push_back(i);        
      m_recvReqs.emplace_back();
      m_recvCounts.back().numPosted++;
      size_t thisMsgSize = std::min(MAX_MESSAGE_SIZE, recvBuf.size - bytesReceived);
      MPI_Irecv(recvBuf.ptr + bytesReceived, thisMsgSize, MPI_BYTE, rank, m_tag, m_comm, &(m_recvReqs.back()));      
      bytesReceived += thisMsgSize;
    }
  }

  m_areRecvsInProgress = true;
}

void DataExchangeKnownPatternUserDataNonBlocking::start_sends(const std::vector<PointerAndSize>& sendData, const std::vector<int>& sendRanks)
{
  STK_ThrowRequireMsg(sendData.size() == sendRanks.size(), "sendData and sendRanks must have same length");
  STK_ThrowRequireMsg(!m_areSendsInProgress, 
                  "cannot start new round of communication until the sends from the previous round are complete");
  
  m_sendReqs.resize(0);
  for (size_t i=0; i < sendRanks.size(); ++i)
  {
    int rank = sendRanks[i];
    PointerAndSize sendBuf = sendData[i];
    int numSegments = computeNumSegments(sendBuf.size, MAX_MESSAGE_SIZE);
    size_t bytesSent = 0;
    
    for (int j=0; j < numSegments; ++j)
    {
      size_t thisMsgSize = std::min(MAX_MESSAGE_SIZE, sendBuf.size - bytesSent);
      m_sendReqs.emplace_back();       
      MPI_Isend(sendBuf.ptr + bytesSent, thisMsgSize, MPI_BYTE, rank, m_tag, m_comm, &(m_sendReqs.back()));
      bytesSent += thisMsgSize;        
    }
  }

  m_areSendsInProgress = true;
}


void DataExchangeKnownPatternUserDataNonBlocking::complete_sends()
{
  STK_ThrowRequireMsg(m_areSendsInProgress, "sends must be started before they can be completed");
  MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
  m_areSendsInProgress = false;
}

}