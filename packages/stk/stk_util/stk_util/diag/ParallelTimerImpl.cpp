#include "stk_util/stk_config.h"
#include "ParallelTimerImpl.hpp"
#include "stk_util/util/Marshal.hpp"

namespace stk::diag::impl {

ParallelTimer::ParallelTimer()
  : m_name(),
    m_timerMask(0),
    m_subtimerLapCount(0),
    m_lapCount(),
    m_cpuTime(),
    m_wallTime(),
    m_MPICount(),
    m_MPIByteCount(),
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
    m_heapAlloc(),
#endif
    m_subtimerList()
{}

ParallelTimer::ParallelTimer(const ParallelTimer &parallel_timer)
  : m_name(parallel_timer.m_name),
    m_timerMask(parallel_timer.m_timerMask),
    m_subtimerLapCount(parallel_timer.m_subtimerLapCount),
    m_lapCount(parallel_timer.m_lapCount),
    m_cpuTime(parallel_timer.m_cpuTime),
    m_wallTime(parallel_timer.m_wallTime),
    m_MPICount(parallel_timer.m_MPICount),
    m_MPIByteCount(parallel_timer.m_MPIByteCount),
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
    m_heapAlloc(parallel_timer.m_heapAlloc),
#endif
    m_subtimerList(parallel_timer.m_subtimerList)
{}

ParallelTimer &ParallelTimer::operator=(const ParallelTimer &parallel_timer) {
  m_name = parallel_timer.m_name;
  m_timerMask = parallel_timer.m_timerMask;
  m_subtimerLapCount = parallel_timer.m_subtimerLapCount;
  m_lapCount = parallel_timer.m_lapCount;
  m_cpuTime = parallel_timer.m_cpuTime;
  m_wallTime = parallel_timer.m_wallTime;
  m_MPICount = parallel_timer.m_MPICount;
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  m_heapAlloc = parallel_timer.m_heapAlloc;
#endif
  m_subtimerList = parallel_timer.m_subtimerList;

  return *this;
}


Writer &
ParallelTimer::dump(Writer &dout) const {
  if (dout.shouldPrint()) {
    dout << "ParallelTimer " << m_name << push << dendl;
    dout << "m_name " << m_name << dendl;
    dout << "m_timerMask " << hex << m_timerMask << dendl;
    dout << "m_subtimerLapCount " << m_subtimerLapCount << dendl;
    dout << "m_lapCount " << m_lapCount << dendl;
    dout << "m_cpuTime " << m_cpuTime << dendl;
    dout << "m_wallTime " << m_wallTime << dendl;
    dout << "m_MPICount " << m_MPICount << dendl;
    dout << "m_MPIByteCount " << m_MPIByteCount << dendl;
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
    dout << "m_heapAlloc " << m_heapAlloc << dendl;
#endif
    dout << "m_subtimerList " << m_subtimerList << dendl;
    dout << pop;
  }
  return dout;
}

void
merge_parallel_timer(
  ParallelTimer &       p0,
  const ParallelTimer & p1,
  bool                  checkpoint)
{
  p0.m_timerMask = p1.m_timerMask;
  p0.m_subtimerLapCount += p1.m_subtimerLapCount;
  p0.m_lapCount.accumulate(p1.m_lapCount, checkpoint);
  p0.m_cpuTime.accumulate(p1.m_cpuTime, checkpoint);
  p0.m_wallTime.accumulate(p1.m_wallTime, checkpoint);
  p0.m_MPICount.accumulate(p1.m_MPICount, checkpoint);
  p0.m_MPIByteCount.accumulate(p1.m_MPIByteCount, checkpoint);
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  p0.m_heapAlloc.accumulate(p1.m_heapAlloc, checkpoint);
#endif

  for (std::list<ParallelTimer>::const_iterator p1_it = p1.m_subtimerList.begin(); p1_it != p1.m_subtimerList.end(); ++p1_it) {
    std::list<ParallelTimer>::iterator p0_it = std::find_if(p0.m_subtimerList.begin(), p0.m_subtimerList.end(), finder((*p1_it).m_name));
    if (p0_it == p0.m_subtimerList.end()) {
      p0.m_subtimerList.push_back((*p1_it));
    }
    else
      merge_parallel_timer(*p0_it, *p1_it, checkpoint);
  }
}

stk::Marshal &operator>>(stk::Marshal &min, ParallelTimer &t) {
  min >> t.m_name >> t.m_timerMask >> t.m_subtimerLapCount
      >> t.m_lapCount.m_value
      >> t.m_lapCount.m_checkpoint
      >> t.m_cpuTime.m_value
      >> t.m_cpuTime.m_checkpoint
      >> t.m_wallTime.m_value
      >> t.m_wallTime.m_checkpoint
      >> t.m_MPICount.m_value
      >> t.m_MPICount.m_checkpoint
      >> t.m_MPIByteCount.m_value
      >> t.m_MPIByteCount.m_checkpoint
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
      >> t.m_heapAlloc.m_value
      >> t.m_heapAlloc.m_checkpoint
#endif
      ;

  min >> t.m_subtimerList;

  return min;
}


#ifdef STK_HAS_MPI
size_t round_up_to_next_word(size_t value)
{
  const size_t SIZE_OF_WORD = 4;
  size_t remainder = value % SIZE_OF_WORD;
  if (remainder == 0) {
      return value;
  }
  return value + SIZE_OF_WORD - remainder;
}
#endif

ParallelTimer
collect_timers(
  const Timer &         root_timer,
  bool                  checkpoint,
  ParallelMachine       comm,
  const int max_procs_per_gather)
{
  Marshal mout;
  mout << root_timer;
  impl::ParallelTimer root_parallel_timer;

#ifdef STK_HAS_MPI
  const int parallel_root = 0 ;
  const int parallel_size = parallel_machine_size(comm);
  const int parallel_rank = parallel_machine_rank(comm);

  // Gather the send counts on root processor
  std::string send_string(mout.str());
  int send_count = send_string.size();
  send_string.resize(round_up_to_next_word(send_count));
  int padded_send_count = send_string.size();


  //We need to gather the timer data in a number of 'cycles' where we
  //only receive from a portion of the other processors each cycle.
  //This is because buffer allocation-failures have been observed for
  //runs on very large numbers of processors if the 'root' processor tries
  //to allocate a buffer large enough to hold timing data from all other
  //procesors.
  //We will set an arbitrary limit for now, making sure that no more than
  //a given number of processors' worth of timer data is gathered at a time.
  int num_cycles = parallel_size/max_procs_per_gather;
  if (parallel_size < max_procs_per_gather || num_cycles < 1) {
    num_cycles = 1;
  }

  std::vector<char> recv_buffer;

  for(int ii=0; ii<num_cycles; ++ii) {
    bool send_this_cycle = (parallel_rank+ii)%num_cycles == 0;

    //send_count is the amount of data this processor needs to send.
    int send_count_this_cycle        = send_this_cycle ? send_count : 0;
    int padded_send_count_this_cycle = send_this_cycle ? padded_send_count : 0;

    std::vector<int> recv_count(parallel_size, 0);
    std::vector<int> padded_recv_count(parallel_size, 0);

    {
      int result = MPI_Gather(&send_count_this_cycle, 1, MPI_INT,
                              recv_count.data(), 1, MPI_INT,
                              parallel_root, comm);
      if (MPI_SUCCESS != result) {
        std::ostringstream message ;
        message << "stk::diag::collect_timers FAILED: send_count MPI_Gather = " << result ;
        throw std::runtime_error(message.str());
      }
    }

    {
      int result = MPI_Gather(&padded_send_count_this_cycle, 1, MPI_INT,
                              padded_recv_count.data(), 1, MPI_INT,
                              parallel_root, comm);
      if (MPI_SUCCESS != result) {
        std::ostringstream message ;
        message << "stk::diag::collect_timers FAILED: padded_send_count MPI_Gather = " << result ;
        throw std::runtime_error(message.str());
      }
    }

    // Receive counts are only non-zero on the root processor:
    std::vector<int> recv_displ(parallel_size + 1, 0);
    std::vector<int> recv_end(parallel_size + 1, 0);

    for (int i = 0 ; i < parallel_size ; ++i) {
      recv_displ[i + 1] = recv_displ[i] + padded_recv_count[i] ;
      recv_end[i] = recv_displ[i] + recv_count[i] ;
    }

    const int recv_size = recv_displ[parallel_size] ;

    recv_buffer.assign(recv_size, 0);

    {
      int result = MPI_Gatherv(send_string.data(), padded_send_count_this_cycle, MPI_CHAR,
                               recv_buffer.data(), padded_recv_count.data(), recv_displ.data(), MPI_CHAR,
                               parallel_root, comm);
      if (MPI_SUCCESS != result) {
        std::ostringstream message ;
        message << "stk::diag::collect_timers FAILED: MPI_Gatherv = " << result ;
        throw std::runtime_error(message.str());
      }

      std::vector<impl::ParallelTimer> parallel_timer_vector;
      parallel_timer_vector.reserve(parallel_size);

      if (parallel_rank == parallel_root) {
        for (int j = 0; j < parallel_size; ++j) {
          int received_count = recv_displ[j+1] - recv_displ[j];
          if (received_count > 0) {
            //grow parallel_timer_vector by 1:
            parallel_timer_vector.resize(parallel_timer_vector.size()+1);
            Marshal min(std::string(recv_buffer.data() + recv_displ[j], recv_buffer.data() + recv_end[j]));
            //put this data into the last entry of parallel_timer_vector:
            min >> parallel_timer_vector[parallel_timer_vector.size()-1];
          }
        }

        if (parallel_rank==parallel_root && send_count_this_cycle>0)
        {
          root_parallel_timer = parallel_timer_vector[0];
        }

        for (size_t j = 0; j < parallel_timer_vector.size(); ++j)
        {
          merge_parallel_timer(root_parallel_timer, parallel_timer_vector[j], checkpoint);
        }
      }
    }
  }
#else
  Marshal min(mout.str());
  min >> root_parallel_timer;
  merge_parallel_timer(root_parallel_timer, root_parallel_timer, checkpoint);
#endif

  return root_parallel_timer;
}

}
