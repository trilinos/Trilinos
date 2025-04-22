#include "DataExchangeUnknownPatternNonBlocking.hpp"

#ifdef STK_HAS_MPI

namespace stk {

UnknownPatternExchanger stringToUnknownPatternExchangerEnum(const std::string& str)
{
  if (str == "Default")
  {
    return UnknownPatternExchanger::Default;
  } else if (str == "Probe")
  {
    return UnknownPatternExchanger::Probe;
  } else if (str == "Prepost")
  {
    return UnknownPatternExchanger::Prepost;
  } else
  {
    throw std::runtime_error(std::string("unable to convert string ") + str + " to UnknownPatternExchanger enum");
  }
}

DataExchangeUnknownPatternNonBlocking::DataExchangeUnknownPatternNonBlocking(MPI_Comm comm, int tag_hint, UnknownPatternExchanger exchanger_type) :
  m_exchanger_type(exchanger_type)
{
  if (exchanger_type == UnknownPatternExchanger::Default)
  {
    char* env_contents = std::getenv("STK_UNKNOWN_PATTERN_EXCHANGER");
    if (env_contents)
    {
      m_exchanger_type = stringToUnknownPatternExchangerEnum(env_contents);          
    } else
    {
      m_exchanger_type = UnknownPatternExchanger::Probe;
    }        
  }

  if (m_exchanger_type == UnknownPatternExchanger::Probe)
  {
    m_exchanger_probe = std::make_shared<stk::DataExchangeUnknownPatternNonBlockingProbe>(comm, tag_hint);
  } else if (m_exchanger_type == UnknownPatternExchanger::Prepost)
  {
    m_exchanger_prepost = std::make_shared<stk::DataExchangeUnknownPatternNonBlockingPrepost>(comm, tag_hint);
  } else
  {
    throw std::runtime_error("unhandled enum value in DataExchangeUnknownPatternNonBlocking");        
  }
}

MPI_Comm DataExchangeUnknownPatternNonBlocking::get_comm() const
{
  if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
    return m_exchanger_probe->get_comm();
  } else
  {
    return m_exchanger_prepost->get_comm();
  }
}


void DataExchangeUnknownPatternNonBlocking::complete_sends()
{
  if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
    return m_exchanger_probe->complete_sends();
  } else
  {
    return m_exchanger_prepost->complete_sends();
  }      
}

bool DataExchangeUnknownPatternNonBlocking::are_sends_in_progress() const
{
  if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
    return m_exchanger_probe->are_sends_in_progress();
  } else
  {
    return m_exchanger_prepost->are_sends_in_progress();
  }       
}

bool DataExchangeUnknownPatternNonBlocking::are_recvs_in_progress() const
{ 
  if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
    return m_exchanger_probe->are_recvs_in_progress();
  } else
  {
    return m_exchanger_prepost->are_recvs_in_progress();
  }
}

}

#endif