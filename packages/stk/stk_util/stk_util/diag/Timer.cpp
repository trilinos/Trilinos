/*------------------------------------------------------------------------*/
/*                 Copyright 2010 - 2011 Sandia Corporation.              */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <math.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <limits>

#include <stk_util/util/string_case_compare.hpp>

#include <stk_util/diag/Timer.hpp>
#include <stk_util/diag/Writer.hpp>
#include <stk_util/diag/WriterExt.hpp>

namespace stk {
namespace diag {

namespace {

MetricsMask s_enabledMetricsMask = METRICS_LAP_COUNT | METRICS_CPU_TIME | METRICS_WALL_TIME;        ///< Bit mask of enabled metrics

template <class T>
typename MetricTraits<T>::Type
value_now() {
  if (MetricTraits<T>::METRIC & getEnabledTimerMetricsMask())
    return MetricTraits<T>::value_now();
  else
    return 0;
}

std::vector<std::string> &
split(
  const std::string &           path,
  char                          separator,
  std::vector<std::string> &    path_vector)
{
  for (std::string::const_iterator it = path.begin(); ; ) {
    std::string::const_iterator it2 = std::find(it, path.end(), separator);
    path_vector.push_back(std::string(it, it2));
    if (it2 == path.end())
      break;
    it = it2 + 1;
  }
  
  return path_vector;
}

} // namespace <empty>


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
  
public:
  static void updateRootTimer(TimerImpl *root_timer);
  
  static Timer createRootTimer(const std::string &name, const TimerSet &timer_set);
    
  static void deleteRootTimer(TimerImpl *root_timer);

  static std::vector<Timer> &findTimers(TimerImpl *root_timer, const std::string &path_tail, std::vector<Timer> &found_timers);

  static void findTimer(TimerImpl *timer, std::vector<std::string> &path_tail_vector, std::vector<Timer> &found_timers);
  
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
  class finder : private std::unary_function<Timer, bool>
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
   *
   * @param timer_mask    a <b>TimerMask</b> value to test the enable timer
   *        bit mask against.
   *
   */
  bool shouldRecord() const {
    return m_timerSet.shouldRecord(m_timerMask) && s_enabledMetricsMask;
  }

  /**
   * Member function <b>getSubtimerLapCount</b> returns the subtimer lap counter.
   *
   * @return      a <b>Counter</b> value of the subtimer lap
   *        counter.
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

public:
  /**
   * Member function <b>dump</b> writes the timer to the specified
   * diagnostic writer.
   *
   * @param dout    a <b>Writer</b> variable reference to write the timer to.
   *
   * @return      a <b>Writer</b> reference to <it>dout</it>.
   */
  Writer &dump(Writer &dout) const;

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

private:
  std::string           m_name;                 ///< Name of the timer
  TimerMask             m_timerMask;            ///< Bit mask to enable timer
  TimerImpl *           m_parentTimer;          ///< Parent timer
  mutable double        m_subtimerLapCount;     ///< Sum of subtimer lap counts and m_lapCount
  unsigned              m_lapStartCount;        ///< Number of pending lap stops

  TimerList             m_subtimerList;         ///< List of subordinate timers

  const TimerSet &              m_timerSet;     ///< Timer enabled mask
  Timer::Metric<LapCount>       m_lapCount;     ///< Number of laps accumulated
  Timer::Metric<CPUTime>        m_cpuTime;      ///< CPU time
  Timer::Metric<WallTime>       m_wallTime;     ///< Wall time
  Timer::Metric<MPICount>       m_MPICount;     ///< MPI call count
  Timer::Metric<MPIByteCount>   m_MPIByteCount; ///< MPI byte count
  Timer::Metric<HeapAlloc>      m_heapAlloc;    ///< Heap allocated
};


/**
 * Member function <b>operator&lt;&lt;</b> ...
 *
 * @param dout      a <b>Writer</b> variable ...
 *
 * @param timer      a <b>TimerImpl</b> variable ...
 *
 * @return      a <b>Writer</b> ...
 */
inline Writer &operator<<(Writer &dout, const TimerImpl &timer) {
  return timer.dump(dout);
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
  timer.m_timerImpl = 0;
}


std::vector<Timer> &
findTimers(Timer root_timer, const std::string &path_tail, std::vector<Timer> &found_timers) {
  TimerImpl::findTimers(root_timer.m_timerImpl, path_tail, found_timers);
  return found_timers;
}


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


template<>
const Timer::Metric<LapCount> &
TimerImpl::getMetric<LapCount>() const {
  return m_lapCount;
}


template<>
const Timer::Metric<CPUTime> &
TimerImpl::getMetric<CPUTime>() const {
  return m_cpuTime;
}


template<>
const Timer::Metric<WallTime> &
TimerImpl::getMetric<WallTime>() const {
  return m_wallTime;
}


template<>
const Timer::Metric<MPICount> &
TimerImpl::getMetric<MPICount>() const {
  return m_MPICount;
}


template<>
const Timer::Metric<MPIByteCount> &
TimerImpl::getMetric<MPIByteCount>() const {
  return m_MPIByteCount;
}


template<>
const Timer::Metric<HeapAlloc> &
TimerImpl::getMetric<HeapAlloc>() const {
  return m_heapAlloc;
}


void
TimerImpl::reset()
{
  m_lapStartCount = 0;

  m_lapCount.reset();
  m_cpuTime.reset();
  m_wallTime.reset();
  m_MPICount.reset();
  m_MPIByteCount.reset();
  m_heapAlloc.reset();
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
    if (m_lapStartCount++ == 0) {
      m_lapCount.m_lapStart = m_lapCount.m_lapStop;

      m_cpuTime.m_lapStop = m_cpuTime.m_lapStart = value_now<CPUTime>();
      m_wallTime.m_lapStop = m_wallTime.m_lapStart = value_now<WallTime>();
      m_MPICount.m_lapStop = m_MPICount.m_lapStart = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = m_MPIByteCount.m_lapStart = value_now<MPIByteCount>();
      m_heapAlloc.m_lapStop = m_heapAlloc.m_lapStart = value_now<HeapAlloc>();
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
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();
    }
  }

  return *this;
}


TimerImpl &
TimerImpl::stop()
{
  if (shouldRecord()) {
    if (--m_lapStartCount <= 0) {
      m_lapStartCount = 0;
      m_lapCount.m_lapStop++;

      m_cpuTime.m_lapStop = value_now<CPUTime>();
      m_wallTime.m_lapStop = value_now<WallTime>();
      m_MPICount.m_lapStop = value_now<MPICount>();
      m_MPIByteCount.m_lapStop = value_now<MPIByteCount>();
      m_heapAlloc.m_lapStop = value_now<HeapAlloc>();

      m_lapCount.addLap();
      m_cpuTime.addLap();
      m_wallTime.addLap();
      m_MPICount.addLap();
      m_MPIByteCount.addLap();
      m_heapAlloc.addLap();
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
  m_heapAlloc.checkpoint();

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
  root_timer->m_heapAlloc.m_lapStop = value_now<HeapAlloc>();

  root_timer->m_lapCount.m_accumulatedLap = root_timer->m_lapCount.m_lapStop - root_timer->m_lapCount.m_lapStart;
  root_timer->m_cpuTime.m_accumulatedLap = root_timer->m_cpuTime.m_lapStop - root_timer->m_cpuTime.m_lapStart;
  root_timer->m_wallTime.m_accumulatedLap = root_timer->m_wallTime.m_lapStop - root_timer->m_wallTime.m_lapStart;
  root_timer->m_MPICount.m_accumulatedLap = root_timer->m_MPICount.m_lapStop - root_timer->m_MPICount.m_lapStart;
  root_timer->m_MPIByteCount.m_accumulatedLap = root_timer->m_MPIByteCount.m_lapStop - root_timer->m_MPIByteCount.m_lapStart;
  root_timer->m_heapAlloc.m_accumulatedLap = root_timer->m_heapAlloc.m_lapStop - root_timer->m_heapAlloc.m_lapStart;
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


void
TimerImpl::findTimer(
  TimerImpl *                   timer,
  std::vector<std::string> &    path_tail_vector,
  std::vector<Timer> &          found_timers)
{
  if (timer->begin() == timer->end()) { // at leaf
    
  }
  else
    for (TimerList::const_iterator it = timer->begin(); it != timer->end(); ++it)
      findTimer((*it).m_timerImpl, path_tail_vector, found_timers);
}


std::vector<Timer> &
TimerImpl::findTimers(
  TimerImpl *                   root_timer,
  const std::string &           path_tail,
  std::vector<Timer> &          found_timers)
{
  std::vector<std::string> path_tail_vector;
  
  findTimer(root_timer, split(path_tail, '.', path_tail_vector), found_timers);

  return found_timers;
}



Writer &
TimerImpl::dump(
  Writer &    dout) const
{
  if (dout.shouldPrint()) {
    dout << "TimerImpl" << push << dendl;
    dout << "m_name, " << m_name << dendl;
    dout << "m_timerMask, " << m_timerMask << dendl;
//    dout << "m_parentTimer, " << c_ptr_name(m_parentTimer) << dendl;
    dout << "m_subtimerLapCount, " << m_subtimerLapCount << dendl;
    dout << "m_lapStartCount, " << m_lapStartCount << dendl;

    dout << "m_lapCount, " << m_lapCount << dendl;
    dout << "m_cpuTime, " << m_cpuTime << dendl;
    dout << "m_wallTime, " << m_wallTime << dendl;
    dout << "m_MPICount, " << m_MPICount << dendl;
    dout << "m_MPIByteCount, " << m_MPIByteCount << dendl;
    dout << "m_heapAlloc, " << m_heapAlloc << dendl;

    dout << "m_subtimerList, " << m_subtimerList << dendl;
    dout << pop;
  }

  return dout;
}


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
template const Timer::Metric<HeapAlloc> &Timer::getMetric<HeapAlloc>() const;


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

Writer &
Timer::dump(Writer& dout) const {
  return m_timerImpl->dump(dout);
}

template <class T>
Writer &
Timer::Metric<T>::dump(
  Writer &    dout) const
{
  if (dout.shouldPrint()) {
    dout << "Timer::Metric<T>" << push << dendl;
    dout << "m_lapStart, " << m_lapStart << dendl;
    dout << "m_lapStop, " << m_lapStop << dendl;
    dout << "m_accumulatedLap, " << m_accumulatedLap << dendl;
    dout << "m_checkpoint, " << m_checkpoint << dendl;
    dout << pop;
  }
  return dout;
}

template Writer &Timer::Metric<LapCount>::dump(Writer &) const;
template Writer &Timer::Metric<CPUTime>::dump(Writer &) const;
template Writer &Timer::Metric<WallTime>::dump(Writer &) const;
template Writer &Timer::Metric<MPICount>::dump(Writer &) const;
template Writer &Timer::Metric<MPIByteCount>::dump(Writer &) const;
template Writer &Timer::Metric<HeapAlloc>::dump(Writer &) const;


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

    if (start_timer)
      m_timer.start();
  }
}


TimeBlockSynchronized::~TimeBlockSynchronized()
{
  if (m_started) {
    try {
      m_timer.stop();
    }
    catch (...) {
    }
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
