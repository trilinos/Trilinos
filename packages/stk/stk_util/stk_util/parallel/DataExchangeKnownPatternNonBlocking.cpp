#include "DataExchangeKnownPatternNonBlocking.hpp"

namespace stk {

DataExchangeKnownPatternNonBlocking::DataExchangeKnownPatternNonBlocking(MPI_Comm comm, int tag_hint) :
  m_exchanger(comm, tag_hint)
{}


void DataExchangeKnownPatternNonBlocking::complete_sends()
{
  m_exchanger.complete_sends();
}

}
