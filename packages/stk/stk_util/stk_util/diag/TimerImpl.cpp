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

#include "stk_util/stk_config.h"
#include "stk_util/diag/TimerImpl.hpp"
#include "stk_util/diag/Timer.hpp"

namespace stk::diag {

namespace {

template <class T>
typename MetricTraits<T>::Type
value_now() {
  if (MetricTraits<T>::METRIC & getEnabledTimerMetricsMask())
    return MetricTraits<T>::value_now();
  else
    return 0;
}

} // namespace <empty>


TimerImpl::TimerImpl(
  const std::string &  name,
  TimerMask    timer_mask,
  TimerImpl *    parent_timer,
  const TimerSet &      timer_set)
  : m_name(name),
    m_timerMask(timer_mask),
    m_parentTimer(parent_timer),
    m_subtimerLapCount(0.0),
    m_lapStartCount(0),
    m_activeChildCount(0),
    m_childCausedStart(false),
    m_subtimerList(),
    m_timerSet(timer_set)
{}


TimerImpl::~TimerImpl()
{
  try {
    for (TimerList::iterator it = m_subtimerList.begin(); it != m_subtimerList.end(); ++it)
      delete (*it).m_timerImpl;
  }
  catch (std::exception &) {
  }
}

bool TimerImpl::shouldRecord() const {
  return m_timerSet.shouldRecord(m_timerMask) && getEnabledTimerMetricsMask();
}

void
TimerImpl::reset()
{
  m_lapStartCount = 0;
  m_childCausedStart = false;
  m_activeChildCount = 0;

  m_lapCount.reset();
  m_cpuTime.reset();
  m_wallTime.reset();
  m_MPICount.reset();
  m_MPIByteCount.reset();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  m_heapAlloc.reset();
#endif
}


Timer
TimerImpl::getSubtimer(
  const std::string &  name)
{
  TimerList::iterator it = std::find_if(m_subtimerList.begin(), m_subtimerList.end(), finder(name));

  if (it == m_subtimerList.end())
    throw std::runtime_error("Timer not found");
  else
    return *it;
}


TimerImpl *
TimerImpl::addSubtimer(
  const std::string &  name,
  TimerMask          timer_mask,
  const TimerSet &      timer_set)
{
  TimerList::iterator it = std::find_if(m_subtimerList.begin(), m_subtimerList.end(), finder(name));

  if (it == m_subtimerList.end()) {
    TimerImpl *timer_impl = new TimerImpl(name, timer_mask, this, timer_set);
    m_subtimerList.push_back(Timer(timer_impl));
    return timer_impl;
  }
  else
    return (*it).m_timerImpl;
}


TimerImpl &
TimerImpl::start()
{
  if (shouldRecord()) {
    if (m_lapStartCount == 0) {
      ++m_lapStartCount;
      m_lapCount.m_lapStart = m_lapCount.m_lapStop;

      m_cpuTime.m_lapStop = m_cpuTime.m_lapStart = value_now<CPUTime>();
      m_wallTime.m_lapStop = m_wallTime.m_lapStart = value_now<WallTime>();
      m_MPICount.m_lapStop = m_MPICount.m_lapStart = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = m_MPIByteCount.m_lapStart = value_now<MPIByteCount>();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
      m_heapAlloc.m_lapStop = m_heapAlloc.m_lapStart = value_now<HeapAlloc>();
#endif
      if(m_parentTimer)
        m_parentTimer->child_notifies_of_start();
    }
  }

  return *this;
}


TimerImpl &
TimerImpl::lap()
{
  if (shouldRecord()) {
    if (m_lapStartCount > 0) {
      m_cpuTime.m_lapStop = value_now<CPUTime>();
      m_wallTime.m_lapStop = value_now<WallTime>();
      m_MPICount.m_lapStop = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();
#endif
    }
  }

  return *this;
}

TimerImpl & TimerImpl::child_notifies_of_start()
{
  //Start only if not already started and this isn't a root timer
  if(m_lapStartCount == 0 && m_parentTimer)
  {
    start();
    m_childCausedStart = true;
  }
  m_activeChildCount++;

  return *this;
}

TimerImpl & TimerImpl::child_notifies_of_stop()
{
  m_activeChildCount--;
  if(m_activeChildCount == 0 && m_childCausedStart)
  {
    stop();
  }
  return *this;
}

TimerImpl &
TimerImpl::stop()
{
  if (shouldRecord()) {
    if (m_lapStartCount > 0) {
      m_lapStartCount = 0;
      m_lapCount.m_lapStop++;
      m_childCausedStart = false;
      m_activeChildCount = 0;

      m_cpuTime.m_lapStop = value_now<CPUTime>();
      m_wallTime.m_lapStop = value_now<WallTime>();
      m_MPICount.m_lapStop = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();
#endif

      m_lapCount.addLap();
      m_cpuTime.addLap();
      m_wallTime.addLap();
      m_MPICount.addLap();
      m_MPIByteCount.addLap();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
      m_heapAlloc.addLap();
#endif
      if(m_parentTimer)
        m_parentTimer->child_notifies_of_stop();
    }
  }

  return *this;
}


double
TimerImpl::accumulateSubtimerLapCounts() const
{
  m_subtimerLapCount = m_lapCount.getAccumulatedLap(false);

  for (TimerList::const_iterator it = m_subtimerList.begin(); it != m_subtimerList.end(); ++it)
    (*it).m_timerImpl->accumulateSubtimerLapCounts();

  for (TimerList::const_iterator it = m_subtimerList.begin(); it != m_subtimerList.end(); ++it)
    m_subtimerLapCount += (*it).m_timerImpl->m_subtimerLapCount;

  return m_subtimerLapCount;
}


void
TimerImpl::checkpoint() const
{
  m_lapCount.checkpoint();
  m_cpuTime.checkpoint();
  m_wallTime.checkpoint();
  m_MPICount.checkpoint();
  m_MPIByteCount.checkpoint();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  m_heapAlloc.checkpoint();
#endif

  for (TimerList::const_iterator it = m_subtimerList.begin(); it != m_subtimerList.end(); ++it)
    (*it).m_timerImpl->checkpoint();
}


void
TimerImpl::updateRootTimer(TimerImpl *root_timer)
{
  root_timer->m_lapCount.m_lapStop = value_now<LapCount>();
  root_timer->m_cpuTime.m_lapStop = value_now<CPUTime>();
  root_timer->m_wallTime.m_lapStop = value_now<WallTime>();
  root_timer->m_MPICount.m_lapStop = value_now<MPICount>();
  root_timer->m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  root_timer->m_heapAlloc.m_lapStop = value_now<HeapAlloc>();
#endif

  root_timer->m_lapCount.m_accumulatedLap = root_timer->m_lapCount.m_lapStop - root_timer->m_lapCount.m_lapStart;
  root_timer->m_cpuTime.m_accumulatedLap = root_timer->m_cpuTime.m_lapStop - root_timer->m_cpuTime.m_lapStart;
  root_timer->m_wallTime.m_accumulatedLap = root_timer->m_wallTime.m_lapStop - root_timer->m_wallTime.m_lapStart;
  root_timer->m_MPICount.m_accumulatedLap = root_timer->m_MPICount.m_lapStop - root_timer->m_MPICount.m_lapStart;
  root_timer->m_MPIByteCount.m_accumulatedLap = root_timer->m_MPIByteCount.m_lapStop - root_timer->m_MPIByteCount.m_lapStart;
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
  root_timer->m_heapAlloc.m_accumulatedLap = root_timer->m_heapAlloc.m_lapStop - root_timer->m_heapAlloc.m_lapStart;
#endif
}



Timer
TimerImpl::createRootTimer(
  const std::string &   name,
  const TimerSet &      timer_set)
{
  TimerImpl *timer_impl = new TimerImpl(name, 0, 0, timer_set);
  return Timer(timer_impl);
}


void
TimerImpl::deleteRootTimer(
  TimerImpl *           root_timer)
{
  delete root_timer;
}


}
