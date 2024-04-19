#include "DataExchangeKnownPatternNonBlocking.hpp"

namespace stk {

DataExchangeKnownPatternNonBlocking::DataExchangeKnownPatternNonBlocking(MPI_Comm comm, int tag_hint) :
  m_comm(comm),
  m_tag(stk::get_mpi_tag_manager().get_tag(comm, tag_hint))
{
  m_sendReqs.reserve(stk::parallel_machine_size(comm));
  m_recvReqs.reserve(stk::parallel_machine_size(comm));
}


DataExchangeKnownPatternNonBlocking::~DataExchangeKnownPatternNonBlocking()
{
  if (m_areSendsInProgress) {
    complete_sends();
  }
}

void DataExchangeKnownPatternNonBlocking::complete_sends()
{
  STK_ThrowRequireMsg(m_areSendsInProgress, "sends must be started before they can be completed");
  MPI_Waitall(m_sendReqs.size(), m_sendReqs.data(), MPI_STATUSES_IGNORE);
  m_areSendsInProgress = false;
}


}
