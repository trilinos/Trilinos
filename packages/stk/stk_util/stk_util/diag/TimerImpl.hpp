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

#ifndef STK_UTIL_DIAG_TimerImpl_hpp
#define STK_UTIL_DIAG_TimerImpl_hpp

#include "stk_util/stk_config.h"
#include "stk_util/diag/TimerMetricTraits.hpp"
#include "stk_util/util/string_case_compare.hpp"  // for equal_case
#include "stk_util/diag/Timer.hpp"
#include "stk_util/util/Writer.hpp"               // for operator<<, Writer, dendl, pop, push
#include "stk_util/diag/WriterExt.hpp"            // for operator<<


namespace stk::diag {


/**
 * Class <b>TimerImpl</b> is the core timer class.  The Timer class is a
 * wrapper around TimerImpl so that the buried references can be constructed more easily.
 *
 * Each timer has a lap counter, cpu timer, wall timer and other metrics.  Each time a timer is
 * started, the cpu start time, wall start time and other metrics, set to the process' current
 * values.  When the timer is stopped, the lap counter is incremented, and the cpu, wall, and other
 * values are accumulated with the difference between now and the start time.
 *
 * Each timer may have a list of subordinate timers.  The relationship is purely
 * hierarchical in that a there is no timing relationship assumed between the timers other
 * than the grouping.  There is no relation between the starting and stopping of parent
 * and subordinate timers.
 *
 * The subordinate timers are stored as pointers to a new timer on the heap, since the
 * calling function will be receiving a reference to this memory which can never change
 * location.  The subordinate timers are not sorted in the list as they should very
 * rarely be created or looked up by name, rather the calling function stores the
 * reference via the Timer class.
 *
 */
class TimerImpl
{
  friend class Timer;
  friend class TimerTester;

public:
  static void updateRootTimer(TimerImpl *root_timer);

  static Timer createRootTimer(const std::string &name, const TimerSet &timer_set);

  static void deleteRootTimer(TimerImpl *root_timer);

private:
  /**
   * Static function <b>reg</b> returns a reference to an existing timer or newly
   * created timer of the specified <b>name</b> which is subordinate to the
   * <b>parent</b> timer.
   *
   * @return      a <b>TimerImpl</b> reference to the timer with the
   *        specified name that is subordinate to the
   *        <b>parent</b> timer.
   */
  static TimerImpl *reg(const std::string &name, TimerMask timer_mask, TimerImpl *parent_timer, const TimerSet &timer_set) {
    return parent_timer->addSubtimer(name, timer_mask, timer_set);
  }

  /**
   * Creates a new <b>Timer</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the timer.
   *
   */
  TimerImpl(const std::string &name, TimerMask timer_mask, TimerImpl *parent_timer, const TimerSet &timer_set);

  /**
   * Destroys a <b>TimerImpl</b> instance.
   *
   */
  ~TimerImpl();

  TimerImpl(const TimerImpl &TimerImpl);
  TimerImpl &operator=(const TimerImpl &TimerImpl);

  /**
   * Class <b>finder</b> is a binary predicate for finding a subordinate timer.
   *
   * Note that the subordinate timer is an unsorted list as there are very few timers
   * created and should rarely be looked up by name.
   */
#ifdef __INTEL_COMPILER
#pragma warning(push)
#pragma warning(disable: 444)
#endif
  class finder
  {
  public:
    explicit finder(const std::string &name)
      : m_name(name)
    {}

    bool operator()(Timer timer) const {
      return equal_case(timer.getName(), m_name);
    }

  private:
    std::string    m_name;
  };
#ifdef __INTEL_COMPILER
#pragma warning(pop)
#endif

public:
  /**
   * Member function <b>getName</b> returns the name of the timer.
   *
   * @return      a <b>std::string</b> const reference to the timer's
   *        name.
   */
  const std::string &getName() const {
    return m_name;
  }

  /**
   * Member function <b>getTimerMask</b> returns the timer mask of the timer.
   *
   * @return      a <b>TimerMask</b> value to the timer mask.
   */
  TimerMask getTimerMask() const {
    return m_timerMask;
  }

  /**
   * Member function <b>getTimerSet</b> returns the timer set of the timer.
   *
   * @return      a <b>TimerSet</b> const reference to the timer set.
   */
  const TimerSet &getTimerSet() const {
    return m_timerSet;
  }

  /**
   * Member function <b>shouldRecord</b> returns true if any of the specified timer
   * bit masks are set in the enable timer bit mask.
   */
  bool shouldRecord() const;

  /**
   * Member function <b>getSubtimerLapCount</b> returns the subtimer lap counter.
   *
   * @return      a <b>Counter</b> value of the subtimer lap counter.
   */
  double getSubtimerLapCount() const {
    return m_subtimerLapCount;
  }

  void setSubtimerLapCount(double value) {
    m_subtimerLapCount = value;
  }

  /**
   * Member function <b>getLapCount</b> returns the lap counter metric.  The lap
   * count metric is the number of times the stop function has been executed.
   *
   * @return      a <b>CounterMetric</b> const reference of the lap counter
   *        metric.
   */
  template <class T>
  const Timer::Metric<T> &getMetric() const;

  /**
   * Member function <b>getTimerList</b> returns the subtimers associated with
   * this timer.
   *
   * @return      a <b>TimerList</b> const reference to the sub
   *        time list.
   */
  const TimerList &getTimerList() const {
    return m_subtimerList;
  }

  TimerList::iterator begin() {
    return m_subtimerList.begin();
  }

  TimerList::const_iterator begin() const {
    return m_subtimerList.begin();
  }

  TimerList::iterator end() {
    return m_subtimerList.end();
  }

  TimerList::const_iterator end() const {
    return m_subtimerList.end();
  }

  /**
   * Member function <b>reset</b> resets the accumulated time and lap times.
   *
   */
  void reset();

  /**
   * Member function <b>checkpoint</b> checkpoints the timer and all subtimers.
   *
   */
  void checkpoint() const;

  /**
   * Member function <b>start</b> sets the start timer.
   *
   * @return      a <b>TimerImpl</b> reference to the timer.
   */
  TimerImpl &start();

  /**
   * Member function <b>lap</b> sets the stop timer.
   *
   * @return      a <b>TimerImpl</b> reference to the timer.
   */
  TimerImpl &lap();

  /**
   * Member function <b>stop</b> sets the stop timer and sums the just completed lap
   * time to the timer.
   *
   * @return      a <b>TimerImpl</b> reference to the timer.
   */
  TimerImpl &stop();

  /**
   * Member function <b>accumulateSubtimerLapCounts</b> sums the lap counter of all
   * subordinate timers.  This is used to determin which timers have been activated at all.
   *
   * @return      an <b>int</b> value of the number of subordinate
   *        timer laps.
   */
  double accumulateSubtimerLapCounts() const;

  Timer getSubtimer(const std::string &name);

private:
  /**
   * Member function <b>addSubtimer</b> returns a reference to an existing or new
   * subtimer with the specified name.
   *
   * @param name    a <b>std::string</b> value of the timer's name.
   *
   * @param timer_mask    a <b>TimerMask</b> value of the class of the timer.
   *
   * @return      a <b>TimerImpl</b> reference to the timer with
   *        specified name.
   */
  TimerImpl *addSubtimer(const std::string &name, TimerMask timer_mask, const TimerSet &timer_set);
  TimerImpl & child_notifies_of_start();
  TimerImpl & child_notifies_of_stop();

private:
  std::string           m_name;                 ///< Name of the timer
  TimerMask             m_timerMask;            ///< Bit mask to enable timer
  TimerImpl *           m_parentTimer;          ///< Parent timer
  mutable double        m_subtimerLapCount;     ///< Sum of subtimer lap counts and m_lapCount
  unsigned              m_lapStartCount;        ///< Number of pending lap stops
  unsigned              m_activeChildCount;     ///< How many children timers have been started
  bool                  m_childCausedStart;     ///< Was this timer started because a child was started?

  TimerList             m_subtimerList;         ///< List of subordinate timers

  const TimerSet &              m_timerSet;     ///< Timer enabled mask
  Timer::Metric<LapCount>       m_lapCount;     ///< Number of laps accumulated
  Timer::Metric<CPUTime>        m_cpuTime;      ///< CPU time
  Timer::Metric<WallTime>       m_wallTime;     ///< Wall time
  Timer::Metric<MPICount>       m_MPICount;     ///< MPI call count
  Timer::Metric<MPIByteCount>   m_MPIByteCount; ///< MPI byte count
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  Timer::Metric<HeapAlloc>      m_heapAlloc;    ///< Heap allocated
#endif
};

template<>
inline const Timer::Metric<LapCount> &
TimerImpl::getMetric<LapCount>() const {
  return m_lapCount;
}


template<>
inline const Timer::Metric<CPUTime> &
TimerImpl::getMetric<CPUTime>() const {
  return m_cpuTime;
}


template<>
inline const Timer::Metric<WallTime> &
TimerImpl::getMetric<WallTime>() const {
  return m_wallTime;
}


template<>
inline const Timer::Metric<MPICount> &
TimerImpl::getMetric<MPICount>() const {
  return m_MPICount;
}


template<>
inline const Timer::Metric<MPIByteCount> &
TimerImpl::getMetric<MPIByteCount>() const {
  return m_MPIByteCount;
}


#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
template<>
STK_DEPRECATED
inline const Timer::Metric<HeapAlloc> &
TimerImpl::getMetric<HeapAlloc>() const {
  return m_heapAlloc;
}
#endif


}

#endif
