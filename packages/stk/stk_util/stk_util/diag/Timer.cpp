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

#include "stk_util/diag/Timer.hpp"
#include "stk_util/diag/TimerImpl.hpp"
#include "stk_util/diag/WriterExt.hpp"            // for operator<<
#include "stk_util/stk_config.h"                  // for STK_HAS_MPI
#include "stk_util/util/Writer.hpp"               // for operator<<, Writer, dendl, pop, push
#include <algorithm>                              // for find_if
#include <exception>                              // for exception
#include <functional>                             // for unary_function
#include <memory>                                 // for shared_ptr, __shared_ptr_access
#include <stdexcept>                              // for runtime_error
#include <vector>                                 // for vector

namespace stk {
namespace diag {

MetricsMask s_enabledMetricsMask = METRICS_LAP_COUNT | METRICS_CPU_TIME | METRICS_WALL_TIME;        ///< Bit mask of enabled metrics

MetricsMask getEnabledTimerMetricsMask()
{
  return s_enabledMetricsMask;
}

void
setEnabledTimerMetricsMask(
  MetricsMask   timer_mask)
{
  s_enabledMetricsMask = timer_mask | METRICS_LAP_COUNT;
}

void
updateRootTimer(
  Timer                 root_timer)
{
  TimerImpl::updateRootTimer(root_timer.m_timerImpl);
}


Timer
createRootTimer(
  const std::string &   name,
  const TimerSet &      timer_set)
{
  return TimerImpl::createRootTimer(name, timer_set);
}


void
deleteRootTimer(
  Timer                 timer)
{
  TimerImpl::deleteRootTimer(timer.m_timerImpl);
  timer.m_timerImpl = nullptr;
}




Timer::~Timer() {}

Timer::Timer(const std::string &name, const Timer parent)
  : m_timerImpl(TimerImpl::reg(name, parent.getTimerMask(), parent.m_timerImpl, parent.getTimerSet()))
{}

Timer::Timer(const std::string &name, TimerMask timer_mask, const Timer parent)
  : m_timerImpl(TimerImpl::reg(name, timer_mask, parent.m_timerImpl, parent.getTimerSet()))
{}

Timer::Timer(const std::string &name, const Timer parent, const TimerSet &timer_set)
  : m_timerImpl(TimerImpl::reg(name, parent.getTimerMask(), parent.m_timerImpl, timer_set))
{}

Timer::Timer(const std::string &name, TimerMask timer_mask, const Timer parent, const TimerSet &timer_set)
  : m_timerImpl(TimerImpl::reg(name, timer_mask, parent.m_timerImpl, timer_set))
{}


const std::string &
Timer::getName() const {
  return m_timerImpl->m_name;
}

TimerMask
Timer::getTimerMask() const {
  return m_timerImpl->getTimerMask();
}

const TimerSet &
Timer::getTimerSet() const {
  return m_timerImpl->getTimerSet();
}

double
Timer::getSubtimerLapCount() const {
  return m_timerImpl->getSubtimerLapCount();
}

const TimerList &
Timer::getTimerList() const {
  return m_timerImpl->getTimerList();
}

template<class T>
const Timer::Metric<T> &
Timer::getMetric() const {
  return m_timerImpl->getMetric<T>();
}

template const Timer::Metric<LapCount> &Timer::getMetric<LapCount>() const;
template const Timer::Metric<CPUTime> &Timer::getMetric<CPUTime>() const;
template const Timer::Metric<WallTime> &Timer::getMetric<WallTime>() const;
template const Timer::Metric<MPICount> &Timer::getMetric<MPICount>() const;
template const Timer::Metric<MPIByteCount> &Timer::getMetric<MPIByteCount>() const;
#ifndef STK_HIDE_DEPRECATED_CODE // Delete after Nov 2025
template const Timer::Metric<HeapAlloc> &Timer::getMetric<HeapAlloc>() const;
#endif


bool
Timer::shouldRecord() const
{
  return m_timerImpl->shouldRecord();
}

TimerList::iterator
Timer::begin()
{
  return m_timerImpl->begin();
}

TimerList::const_iterator
Timer::begin() const
{
  return m_timerImpl->begin();
}

TimerList::iterator
Timer::end()
{
  return m_timerImpl->end();
}

TimerList::const_iterator
Timer::end() const
{
  return m_timerImpl->end();
}

double
Timer::accumulateSubtimerLapCounts() const {
  return m_timerImpl->accumulateSubtimerLapCounts();
}

Timer &
Timer::start() {
  m_timerImpl->start();
  return *this;
}

Timer &
Timer::lap() {
  m_timerImpl->lap();
  return *this;
}

Timer &
Timer::stop() {
  m_timerImpl->stop();
  return *this;
}

void
Timer::checkpoint() const {
  m_timerImpl->checkpoint();
}

TimeBlockSynchronized::TimeBlockSynchronized(
  Timer &    timer,
  MPI_Comm    mpi_comm,
  bool      start_timer)
  : m_timer(timer),
    m_mpiComm(mpi_comm),
    m_started(start_timer)
{
  if (m_timer.m_timerImpl->shouldRecord()) {
#ifdef STK_HAS_MPI
    if (mpi_comm != MPI_COMM_NULL)
      MPI_Barrier(mpi_comm);
#endif
    if (start_timer) {
      m_timer.start();
    }
  }
}

TimeBlockSynchronized::~TimeBlockSynchronized()
{
  if (m_started) {
    m_timer.stop();
  }
}

void
TimeBlockSynchronized::start()
{
  // Place barrier here
  MPI_Barrier(m_mpiComm);
  m_started = true;
  m_timer.start();
}

void
TimeBlockSynchronized::stop()
{
  m_started = false;
  m_timer.stop();
  // Does a barrier need to be here?
  // MPI_Barrier(Env::parallel_comm());
}

} // namespace diag
} // namespace stk

namespace sierra {
namespace Diag {

//
// SierraRootTimer member functions:
//
SierraRootTimer::SierraRootTimer()
  : m_sierraTimer(stk::diag::createRootTimer("Sierra", sierraTimerSet()))
{ }


SierraRootTimer::~SierraRootTimer()
{
  stk::diag::deleteRootTimer(m_sierraTimer);
}


stk::diag::Timer & SierraRootTimer::sierraTimer()
{
  return m_sierraTimer;
}


TimerSet &
sierraTimerSet()
{
  static TimerSet s_sierraTimerSet(TIMER_PROCEDURE | TIMER_REGION);

  return s_sierraTimerSet;
}


std::shared_ptr<SierraRootTimer> sierraRootTimer()
{
  static std::shared_ptr<SierraRootTimer> s_sierraRootTimer(new SierraRootTimer());
  if ( ! s_sierraRootTimer ) {
    s_sierraRootTimer.reset(new SierraRootTimer());
  }
  return s_sierraRootTimer;
}


Timer &
sierraTimer()
{
  return sierraRootTimer()->sierraTimer();
}

void
sierraTimerDestroy()
{
  sierraRootTimer().reset();
}


void
setEnabledTimerMask(
  TimerMask     timer_mask)
{
  sierraTimerSet().setEnabledTimerMask(timer_mask);
}


TimerMask
getEnabledTimerMask()
{
  return sierraTimerSet().getEnabledTimerMask();
}


void
setTimeFormat(int time_format) {
  stk::diag::setTimerTimeFormat(time_format);
}


void
setTimeFormatMillis()
{
  if ((getTimeFormat() & stk::TIMEFORMAT_STYLE_MASK ) == stk::TIMEFORMAT_HMS) {
    if (getSierraWallTime() > 3600.0)
      setTimeFormat(getTimeFormat() & ~stk::TIMEFORMAT_MILLIS);
    else
      setTimeFormat(getTimeFormat() | stk::TIMEFORMAT_MILLIS);
  }
  else if ((getTimeFormat() & stk::TIMEFORMAT_STYLE_MASK ) == stk::TIMEFORMAT_SECONDS) {
    if (getSierraWallTime() > 1000.0)
      setTimeFormat(getTimeFormat() & ~stk::TIMEFORMAT_MILLIS);
    else
      setTimeFormat(getTimeFormat() | stk::TIMEFORMAT_MILLIS);
  }
}

int
getTimeFormat()
{
  return stk::diag::getTimerTimeFormat();
}

void
setTimerNameMaxWidth(
  size_t        /*width*/)
{
}

stk::diag::MetricTraits<stk::diag::WallTime>::Type
getSierraWallTime()
{
  return sierraTimer().getMetric<stk::diag::WallTime>().getAccumulatedLap(false);
}

stk::diag::MetricTraits<stk::diag::CPUTime>::Type
getCPULapTime(Timer timer) {
  return timer.getMetric<stk::diag::CPUTime>().getLap();
}

TimerParser &
theTimerParser()
{
  static TimerParser parser;

  return parser;
}


TimerParser::TimerParser()
  : sierra::OptionMaskParser(),
    m_metricsSetMask(0),
    m_metricsMask(0)
{
  mask("cpu", 0, "Display CPU times");
  mask("wall", 0, "Display wall times");

  mask("hms", 0, "Display times in HH:MM:SS format");
  mask("seconds", 0, "Display times in seconds");

  mask("all", TIMER_ALL, "Enable all metrics");
  mask("none", TIMER_NONE, "Disable all timers");

  mask("domain", TIMER_DOMAIN, "Enable metrics on the domain");
  mask("region", TIMER_REGION, "Enable metrics on regions");
  mask("procedure", TIMER_PROCEDURE, "Enable metrics on procedures");
  mask("mechanics", TIMER_MECHANICS, "Enable metrics on mechanics");
  mask("algorithm", TIMER_ALGORITHM, "Enable metrics on algorithms");
  mask("solver", TIMER_SOLVER, "Enable metrics on solvers");
  mask("contact", TIMER_CONTACT, "Enable metrics on contact");
  mask("material", TIMER_MATERIAL, "Enable metrics on materials");
  mask("search", TIMER_SEARCH, "Enable metrics on searches");
  mask("transfer", TIMER_TRANSFER, "Enable metrics on user functions");
  mask("adaptivity", TIMER_ADAPTIVITY, "Enable metrics on adaptivity");
  mask("recovery", TIMER_RECOVERY, "Enable metrics on encore recovery");
  mask("profile1", TIMER_PROFILE_1, "Enable app defined profiling metrics");
  mask("profile2", TIMER_PROFILE_2, "Enable app defined profiling metrics");
  mask("profile3", TIMER_PROFILE_3, "Enable app defined profiling metrics");
  mask("profile4", TIMER_PROFILE_4, "Enable app defined profiling metrics");
  mask("app1", TIMER_APP_1, "Enable app defined metrics");
  mask("app2", TIMER_APP_2, "Enable app defined metrics");
  mask("app3", TIMER_APP_3, "Enable app defined metrics");
  mask("app4", TIMER_APP_4, "Enable app defined metrics");
}


OptionMaskParser::Mask
TimerParser::parse(
  const char *          option_mask) const
{
  m_metricsSetMask = 0;
  m_metricsMask = 0;
  m_optionMask = getEnabledTimerMask();

  m_optionMask = OptionMaskParser::parse(option_mask);

  setEnabledTimerMask(m_optionMask);

  if (m_metricsSetMask != 0)
    stk::diag::setEnabledTimerMetricsMask(m_metricsMask);

  return m_optionMask;
}


void
TimerParser::parseArg(
  const std::string &   name,
  const std::string &   arg) const
{
  if (name == "cpu") {
    m_metricsMask |= stk::diag::METRICS_CPU_TIME;
    m_metricsSetMask |= stk::diag::METRICS_CPU_TIME;
  }
  else if (name == "wall") {
    m_metricsMask |= stk::diag::METRICS_WALL_TIME;
    m_metricsSetMask |= stk::diag::METRICS_WALL_TIME;
  }
  else if (name == "heap") {
    m_metricsMask |= stk::diag::METRICS_HEAP_ALLOC;
    m_metricsSetMask |= stk::diag::METRICS_HEAP_ALLOC;
  }
  else if (name == "none") {
    m_optionMask = 0;
    m_metricsSetMask = stk::diag::METRICS_WALL_TIME | stk::diag::METRICS_CPU_TIME;
  }

  else if (name == "hms") {
    Diag::setTimeFormat(stk::TIMEFORMAT_HMS);
  }
  else if (name == "seconds") {
    Diag::setTimeFormat(stk::TIMEFORMAT_SECONDS);
  }

  else
    OptionMaskParser::parseArg(name, arg);
}

} // namespace Diag
} // namespace sierra
