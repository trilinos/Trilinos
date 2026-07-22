#ifndef stk_util_parallel_DataExchangeUnknownPatternNonBlocking_hpp
#define stk_util_parallel_DataExchangeUnknownPatternNonBlocking_hpp

#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingPrepost.hpp"
#include "stk_util/parallel/DataExchangeUnknownPatternNonBlockingProbe.hpp"
#include "stk_util/parallel/Parallel.hpp"   // for MPI
#include "DataExchangeUnknownPatternNonBlockingProbe.hpp"
#include "DataExchangeUnknownPatternNonBlockingPrepost.hpp"
#include <string>

#ifdef STK_HAS_MPI

namespace stk {
  
enum class UnknownPatternExchanger
{
  Default,
  Probe,
  Prepost,
};

UnknownPatternExchanger stringToUnknownPatternExchangerEnum(const std::string& str);


class DataExchangeUnknownPatternNonBlocking
{
  public:
    static constexpr int Unknown = DataExchangeUnknownPatternNonBlockingProbe::Unknown;
    static constexpr size_t MAX_MESSAGE_SIZE = std::numeric_limits<int>::max();

    DataExchangeUnknownPatternNonBlocking(MPI_Comm comm, int tag_hint=11173, UnknownPatternExchanger exchanger_type=UnknownPatternExchanger::Default);

    ~DataExchangeUnknownPatternNonBlocking() {}

    MPI_Comm get_comm() const;

    template <typename T>
    void start_nonblocking(std::vector< std::vector<T> > &sendLists,
                           std::vector< std::vector<T> > &recvLists,
                           int numRecvsExpected=Unknown)
    {
      if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
        return m_exchanger_probe->start_nonblocking(sendLists, recvLists, numRecvsExpected);
      } else
      {
        return m_exchanger_prepost->start_nonblocking(sendLists, recvLists, numRecvsExpected);
      } 
    }

    template <typename T>
    void post_nonblocking_receives(std::vector< std::vector<T> > &recvLists)
    {
      if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
        return m_exchanger_probe->post_nonblocking_receives(recvLists);
      } else
      {
        return m_exchanger_prepost->post_nonblocking_receives(recvLists);
      }      
    }

    // Wait for the receives to finish, calling func on each receive
    // buffer as it becomes ready.
    // func must be callable as func(int rank, std::vector<T> recv_buf)
    // (although taking the recv_buf by reference is usually better)
    template <typename T, typename Tfunc>
    void complete_receives(std::vector< std::vector<T> >& recvLists, Tfunc func)
    {
      if (m_exchanger_type == UnknownPatternExchanger::Probe) { 
        return m_exchanger_probe->complete_receives(recvLists, func);
      } else
      {
        return m_exchanger_prepost->complete_receives(recvLists, func);
      }
    }

    // wait for sends to complete, after which time the caller can overwrite the send
    // buffers
    void complete_sends();

    bool are_sends_in_progress() const;

    bool are_recvs_in_progress() const;

  private:
    UnknownPatternExchanger m_exchanger_type;
    std::shared_ptr<DataExchangeUnknownPatternNonBlockingProbe> m_exchanger_probe;
    std::shared_ptr<DataExchangeUnknownPatternNonBlockingPrepost> m_exchanger_prepost;
};

  
}

#endif
#endif