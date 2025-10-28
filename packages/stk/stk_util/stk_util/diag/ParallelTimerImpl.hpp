// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#ifndef STK_UTIL_DIAG_ParallelTimerImpl_hpp
#define STK_UTIL_DIAG_ParallelTimerImpl_hpp

#include "stk_util/stk_config.h"
#include "stk_util/diag/Timer.hpp"
#include "stk_util/util/Writer.hpp"
#include "WriterExt.hpp"
#include "stk_util/util/string_case_compare.hpp"  // for equal_case
#include "TimerMetricTraits.hpp"
#include <limits>
#include <typeinfo>
#include <list>

namespace stk { struct Marshal; }

namespace stk::diag {
namespace impl {

struct ParallelTimer
{
  template <typename T>
  struct Metric
  {
    Metric()
      : m_value(0),
        m_sum(0.0),
        m_min(std::numeric_limits<double>::max()),
        m_max(0.0)
    {}

    typename MetricTraits<T>::Type  m_value;  ///< Metric value
    typename MetricTraits<T>::Type  m_checkpoint;  ///< Metric checkpointed value
    double                          m_sum;    ///< Reduction sum
    double                              m_min;    ///< Reduction min
    double                    m_max;          ///< Reduction max

    void accumulate(const Metric<T> &metric, bool checkpoint) {
      double value = static_cast<double>(metric.m_value);
      if (checkpoint)
        value -= static_cast<double>(metric.m_checkpoint);

      m_sum += value;
      m_min = std::min(m_min, value);
      m_max = std::max(m_max, value);
    }

    Writer &dump(Writer &dout) const {
      if (dout.shouldPrint()) {
        dout << "Metric<" << typeid(typename MetricTraits<T>::Type) << ">" << push << dendl;
        dout << "m_value " << m_value << dendl;
        dout << "m_checkpoint " << m_value << dendl;
        dout << "m_sum " << m_sum << dendl;
        dout << "m_min " << m_min << dendl;
        dout << "m_max " << m_max << dendl;
        dout << pop;
      }
      return dout;
    }
  };

  ParallelTimer();

  ParallelTimer(const ParallelTimer &parallel_timer);

  ParallelTimer &operator=(const ParallelTimer &parallel_timer);

  template <class T>
  const Metric<T> &getMetric() const;

  std::string                   m_name;                 ///< Name of the timer
  TimerMask                     m_timerMask;
  double                        m_subtimerLapCount;     ///< Sum of subtimer lap counts and m_lapCount

  Metric<LapCount>              m_lapCount;             ///< Number of laps accumulated
  Metric<CPUTime>               m_cpuTime;              ///< CPU time
  Metric<WallTime>              m_wallTime;             ///< Wall time
  Metric<MPICount>              m_MPICount;             ///< MPI call count
  Metric<MPIByteCount>          m_MPIByteCount;	        ///< MPI byte count
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  Metric<HeapAlloc>             m_heapAlloc;            ///< MPI byte count
#endif

  std::list<ParallelTimer>      m_subtimerList;         ///< Sub timers

  Writer &dump(Writer &dout) const;
};

template<>
inline const ParallelTimer::Metric<LapCount> &
ParallelTimer::getMetric<LapCount>() const {
  return m_lapCount;
}


template<>
inline const ParallelTimer::Metric<CPUTime> &
ParallelTimer::getMetric<CPUTime>() const {
  return m_cpuTime;
}


template<>
inline const ParallelTimer::Metric<WallTime> &
ParallelTimer::getMetric<WallTime>() const {
  return m_wallTime;
}


template<>
inline const ParallelTimer::Metric<MPICount> &
ParallelTimer::getMetric<MPICount>() const {
  return m_MPICount;
}


template<>
inline const ParallelTimer::Metric<MPIByteCount> &
ParallelTimer::getMetric<MPIByteCount>() const {
  return m_MPIByteCount;
}

#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
template<>
STK_DEPRECATED inline const ParallelTimer::Metric<HeapAlloc> &
ParallelTimer::getMetric<HeapAlloc>() const {
  return m_heapAlloc;
}
#endif


template <typename T>
Writer &operator<<(Writer &dout, const ParallelTimer::Metric<T> &t) {
  return t.dump(dout);
}

inline Writer &operator<<(Writer &dout, const ParallelTimer &parallel_timer) {
  return parallel_timer.dump(dout);
}

stk::Marshal &operator>>(stk::Marshal &min, ParallelTimer &t);

#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 444)
#endif
class finder
{
public:
  finder(const std::string &name)
    : m_name(name)
  {}

  bool operator()(const ParallelTimer &parallel_timer) const {
    return equal_case(parallel_timer.m_name, m_name);
  }

private:
  std::string           m_name;
};
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

void
merge_parallel_timer(
  ParallelTimer &       p0,
  const ParallelTimer & p1,
  bool                  checkpoint);

ParallelTimer
collect_timers(
  const Timer &               root_timer,
  bool                  checkpoint,
  ParallelMachine       comm,
  const int max_procs_per_gather = 64);

}
}

#endif
