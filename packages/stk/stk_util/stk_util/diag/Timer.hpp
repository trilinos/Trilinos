/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_UTIL_DIAG_Timer_hpp
#define STK_UTIL_DIAG_Timer_hpp

#include <iosfwd>
#include <vector>
#include <list>
#include <string>

#include <mpi.h>

#include <stk_util/diag/TimerMetricTraits.hpp>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_util/environment/FormatTime.hpp>
#include <stk_util/diag/Writer_fwd.hpp>

#include <stk_util/diag/String.hpp>
#include <stk_util/diag/WriterParser.hpp>
#include <stk_util/diag/Option.hpp>


///
/// @addtogroup DiagTimerDetail
/// @{
///

namespace stk {
namespace diag {

class Timer;
class TimerSet;
class TimerImpl;

typedef unsigned TimerMask;        ///< Timer classification mask

/**
 * Function <b>getEnabledMetricsMask</b> retruns the timer enable bit mask.
 *
 * @return      a <b>MetricsMask</b> value of the timer enable bit
 *        mask.
 */
MetricsMask getEnabledTimerMetricsMask();

/**
 * Function <b>setEnabledMetricsMask</b> set the timer enable bit mask to
 * <b>timer_mask</b>.
 *
 * @param timer_mask    a <b>MetricsMask</b> value to set the timer enable bit
 *        mask to.
 *
 */
void setEnabledTimerMetricsMask(MetricsMask timer_mask);

/**
 * Function <b>updateRootTimer</b> updates the root timers stop and total
 * metric values with the current time.
 *
 * @param root_timer      a <b>Timer</b> reference to the root timer.
 *
 */
void updateRootTimer(Timer root_timer);

/**
 * Function <b>createRootTimer</b> creates a root timer.  Root timers are the root of a timer
 * hierarchy.  The timer_set specifies the timer groupings for this root timer.  The percentage of a
 * child timer is the ratio of that timer the its root.
 *
 * @param name                  a <b>std::string</b> const reference to the name of the new root
 *                              timer.
 *
 * @param timer_set           a <b>TimerSet</b> const reference of the timer set of the new root
 *                              timer.
 *
 * @return      a <b>Timer</b> value of the new root timer.
 */
Timer createRootTimer(const std::string &name, const TimerSet &timer_set);

/**
 * Function <b>deleteRootTimer</b> deletes a root timer and all of it's children timers.  All
 * children Timers are invalidated and can no longer be used.
 *
 * @param                       a <b>Timer</b> value of the root timer to delete.
 */
void deleteRootTimer(Timer timer);

/**
 * @brief Member function <code>findTimer</code> return a vector of timers whose tail of the dot
 * separated name from root_time to leaf matches the specified path_tail.
 *
 * @param root_timer    a <code>Timer</code> value of the root to begin search.
 *
 * @param path_tail    a <code>std::string</code> const reference to the dot separated tail
 *                              to match.
 *
 * @param found_timer    a <code>std::vector<Timer></code> reference to the vector to store
 *                              matching timers.
 *
 * @return      a <code>std::vector<Timer></code> reference to found_timer.
 */
std::vector<Timer> &findTimers(Timer root_timer, const std::string &path_tail, std::vector<Timer> &found_timers);

/**
 * @brief Class <b>TimerSet</b> implements a set of timer classifications.  A time classification
 * consists of a bit mask set TimerMask
 *
 */
class TimerSet
{
public:
  explicit TimerSet(TimerMask enabled_timer_mask)
    : m_enabledTimerMask(enabled_timer_mask)
  {}

private:
  TimerSet(const TimerSet &timer_set)
    : m_enabledTimerMask(timer_set.m_enabledTimerMask)
  {}

  TimerSet &operator=(TimerSet &timer_set) {
    m_enabledTimerMask = timer_set.m_enabledTimerMask;

    return *this;
  }

public:
  ~TimerSet()
  {}

  /**
   * Member function <b>getEnabledTimerMask</b> returns the timer enable bit mask.
   *
   * @return      a <b>TimerMask</b> value of the timer enable bit
   *        mask.
   */
  TimerMask getEnabledTimerMask() const {
    return m_enabledTimerMask;
  }

  /**
   * Member function <b>setEnabledTimerMask</b> set the timer enable bit mask to
   * <b>timer_mask</b>.
   *
   * @param timer_mask    a <b>TimerMask</b> value to set the timer enable bit
   *        mask to.
   *
   */
  void setEnabledTimerMask(TimerMask timer_mask) {
    m_enabledTimerMask = timer_mask;
  }

  /**
   * Member function <b>shouldRecord</b> returns true if any of the specified timer
   * bit masks are set in the enable timer bit mask.
   *
   * @param timer_mask    a <b>TimerMask</b> value to test the enable timer
   *        bit mask against.
   *
   */
  bool shouldRecord(TimerMask timer_mask) const {
    return (timer_mask == 0 || (m_enabledTimerMask & timer_mask));
  }

private:
  TimerMask    m_enabledTimerMask;  ///< Bit mask of enabled timer
};


typedef std::list<Timer> TimerList;    ///< A vector of subordinate timers.

/**
 * @brief Class <b>Timer</b> implements a diagnostic timer and timer container for the
 * collection and display of execution times.
 *
 */
class Timer
{
  friend class TimerImpl;
  friend class TimeBlock;
  friend class TimeBlockSynchronized;
  friend void updateRootTimer(Timer);
  friend Timer createRootTimer(const std::string &, const TimerSet &);
  friend void deleteRootTimer(Timer);
  friend std::vector<Timer> &findTimers(Timer, const std::string &, std::vector<Timer> &);

public:
  /**
   * Class <b>Metric</b> maintains the metric data for the timer or counter.  The
   * start and stop times maintain the current lap time.  When a lap completes, its
   * time/count is accumlated to the total.  The total time/count can be stored in the
   * checkpoint member variable.  The total can be retrieved as either absolute time/count
   * the diffence from the checkpoint value.
   *
   */
  template <typename T>
  struct Metric
  {
    Metric()
      : m_lapStart(0),
        m_lapStop(0),
        m_accumulatedLap(0),
        m_checkpoint(0)
    {}

    /**
     * Member function <b>reset</b> resets the metric values to zero.
     *
     */
    void reset() {
      m_lapStart = m_lapStop = m_accumulatedLap = m_checkpoint = 0;
    }

    /**
     * Member function <b>addLap</b> adds the most recently completed lap to the total.
     *
     * @return      a <b>T</b> value of the total.
     */
    typename MetricTraits<T>::Type addLap() {
      return m_accumulatedLap += m_lapStop - m_lapStart;
    }

    /**
     * Member function <b>checkpoint</b> checkpoints the metrics by storing the
     * total time in the checkpoint value.
     *
     */
    void checkpoint() const {
      m_checkpoint = m_accumulatedLap;
    }

    /**
     * Member function <b>getLap</b> returns the value of the most recently
     * lap.
     *
     * @return      a <b>T</b> value of the most recent lap.
     */
    typename MetricTraits<T>::Type getLap() const {
      return m_lapStop - m_lapStart;
    }

    /**
     * Member function <b>getStart</b> returns the start value of the most recent lap.
     *
     * @return      a <b>T</b> value of the start of the most recent lap.
     */
    typename MetricTraits<T>::Type getStart() const {
      return m_lapStart;
    }

    /**
     * Member function <b>getStop</b> returns the stop value of the most recent lap.
     *
     * @return      a <b>T</b> value of the stop of the most recent lap.
     */
    typename MetricTraits<T>::Type getStop() const {
      return m_lapStop;
    }

    /**
     * Member function <b>getAccumulatedLap</b> returns the accumulated value of the metric.
     * If the <b>checkpoint</b> parameter if true, the value returned is the
     * difference between the accumulated value and the checkpointed value.
     *
     * @param checkpoint  a <b>bool</b> value of true of the checkpointed
     *        value is to be returned.
     *
     * @return      a <b>T</b> value of the accumulated or the
     *        checkpoint difference.
     */
    typename MetricTraits<T>::Type getAccumulatedLap(bool arg_checkpoint = false) const {
      if (arg_checkpoint)
        return m_accumulatedLap - m_checkpoint;
      else
        return m_accumulatedLap;
    }

    /**
     * Member function <b>dump</b> prints the value of the Metric to the
     * diagnostic writer.
     *
     * @param dout    a <b>Writer</b> reference to the diagnostic
     *        writer to write to.
     *
     * @return      a <b>Writer</b> reference to the diagnostic
     *        writer.
     */
    Writer &dump(Writer &dout) const;

    typename MetricTraits<T>::Type    m_lapStart;    ///< Most recent start time/count
    typename MetricTraits<T>::Type    m_lapStop;    ///< Most recent stop or lap time/count
    typename MetricTraits<T>::Type    m_accumulatedLap;  ///< Accumulated time/count
    mutable typename MetricTraits<T>::Type      m_checkpoint;    ///< Checkpointed time/count
  };

  /**
   * Creates a new <b>Timer</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the timer.
   *
   * @param parent    a <b>Timer</b> value of the parent timer.
   *
   */
  Timer(const std::string &name, const Timer parent);

  /**
   * Creates a new <b>Timer</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the timer.
   *
   * @param parent    a <b>Timer</b> value of the parent timer.
   *
   * @param timer_set    a <b>TimerSet</b> value of the timer set used to interpret the
   *                            TimerMask's of this and children timers.
   *
   */
  Timer(const std::string &name, const Timer parent, const TimerSet &timer_set);

  /**
   * Creates a new <b>Timer</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the timer.
   *
   * @param timer_mask    a <b>TimerMask</b> value which enables this timer.
   *
   * @param parent    a <b>Timer</b> value of the parent timer.
   *
   */
  Timer(const std::string &name, TimerMask timer_mask, const Timer parent);

  /**
   * Creates a new <b>Timer</b> instance.
   *
   * @param name    a <b>std::string</b> const reference to the name of
   *        the timer.
   *
   * @param timer_mask    a <b>TimerMask</b> value which enables this timer.
   *
   * @param parent    a <b>Timer</b> value of the parent timer.
   *
   * @param timer_set    a <b>TimerSet</b> value of the timer set used to interpret the
   *                            TimerMask's of this and children timers.
   *
   */
  Timer(const std::string &name, TimerMask timer_mask, const Timer parent, const TimerSet &timer_set);

  /**
   * Creates the root <b>Timer</b> timer instance.
   *
   */
  explicit Timer(TimerImpl &timer_impl)
    : m_timerImpl(&timer_impl)
  {}

  explicit Timer(TimerImpl *timer_impl)
    : m_timerImpl(timer_impl)
  {}

  Timer(const Timer &timer)
    : m_timerImpl(timer.m_timerImpl)
  {}

  Timer &operator=(const Timer &timer) {
    if (this != &timer)
      m_timerImpl = timer.m_timerImpl;

    return *this;
  }

  virtual ~Timer()
  {}

  const TimerList &getTimerList() const;

  TimerList::iterator begin();
  TimerList::const_iterator begin() const;
  TimerList::iterator end();
  TimerList::const_iterator end() const;

  /**
   * Member function <b>getName</b> returns the name of the timer.
   *
   * @return      a <b>std::string</b> const reference to the timer's
   *        name.
   */
  const std::string &getName() const;


  /**
   * Member function <b>getTimerMask</b> returns the timer mask of the timer.
   *
   * @return      a <b>TimerMask</b> value to the timer mask.
   */
  const TimerSet &getTimerSet() const;

  /**
   * Member function <b>getTimerMask</b> returns the timer mask of the timer.
   *
   * @return      a <b>TimerMask</b> value to the timer mask.
   */
  TimerMask getTimerMask() const;

  bool shouldRecord() const;

  /**
   * Member function <b>getSubtimerLapCount</b> returns the subtimer lap counter.
   *
   * @return      a <b>Counter</b> value of the subtimer lap
   *        counter.
   */
  double getSubtimerLapCount() const;

  /**
   * Member function <b>getLapCount</b> returns the lap counter metric.  The lap
   * count metric is the number of times the stop function has been executed.
   *
   * @return      a <b>CounterMetric</b> const reference of the lap counter
   *        metric.
   */
  template <class T>
  const Metric<T> &getMetric() const;

  /**
   * Member function <b>accumulateSubtimerLapCounts</b> accumulates the subtimer la
   * counts.
   *
   * @return      an <b>int</b> value of the count.
   */
  double accumulateSubtimerLapCounts() const;

  /**
   * Member function <b>start</b> starts the lap timer.
   *
   * @return      a <b>Timer</b> reference to this timer.
   */
  Timer &start();

  /**
   * Member function <b>lap</b> sets the lap stop time.
   *
   * @return      a <b>Timer</b> reference to the timer.
   */
  Timer &lap();

  /**
   * Member function <b>stop</b> sets the lap stop time and sums the just completed
   * lap time to the timer.
   *
   * @return      a <b>Timer</b> reference to the timer.
   */
  Timer &stop();

  /**
   * Member function <b>checkpoint</b> checkpoints the metrics by storing the
   * total time in the checkpoint value.
   *
   */
  void checkpoint() const;

  /**
   * Member function <b>dump</b> writes the timer to the specified
   * diagnostic writer.
   *
   * @param dout    a <b>Writer</b> variable reference to write the timer to.
   *
   * @return      a <b>Writer</b> reference to <i>dout</i>.
   */
  Writer &dump(Writer& dout) const;

private:
  TimerImpl *    m_timerImpl;      ///< Reference to the actual timer
};



/**
 * Class <b>TimeBlock</b> is a time sentry for timing a statement block.  The
 * timer is generally started upon construction. But, the start is delayed if the second
 * argument is false.  In this case, manually start the timer by calling the start()
 * function.  This gives the safety of using a sentry, but does not force to awkwardness
 * associated with local variables crossing the timed block.
 *
 */
class TimeBlock
{
public:
  /**
   * Creates a new <b>TimeBlock</b> instance.  The newly created instance will
   * start the timer if the <b>start</b> value is true, which is the default
   * case.  If the <b>start</b> value is false, the calling function is
   * responsible for starting the timer at the appropriate time.
   *
   * @param timer    a <b>Timer</b> reference to the timer accumulate
   *        block run times.
   *
   * @param start_timer  a <b>bool</b> value to have the timer started on
   *        construction.
   *
   */
  explicit TimeBlock(Timer &timer, bool start_timer = true)
    : m_timer(timer),
      m_started(start_timer)
  {
    if (start_timer)
      m_timer.start();
  }

private:
  TimeBlock(const TimeBlock &);
  TimeBlock &operator=(const TimeBlock &);

public:
  /**
   * Destroys a <b>TimeBlock</b> instance.  Stops the timer if is has been started.
   *
   */
  ~TimeBlock() {
    try {
      if (m_started)
        m_timer.stop();
    }
    catch (...) {
    }
  }

  /**
   * Member function <b>start</b> starts the timer associated with the time block.
   *
   */
  void start() {
    m_started = true;
    m_timer.start();
  }

  /**
   * Member function <b>lap</b> sets the stop time of the timer associated with
   * the time block.
   *
   */
  void lap() {
    m_timer.lap();
  }

  /**
   * Member function <b>stop</b> stops the timer associated with the time block.
   *
   */
  void stop() {
    m_started = false;
    m_timer.stop();
  }

private:
  Timer &               m_timer;  ///< Timer to accumulate block run times.
  bool      m_started;  ///< Timer has been started
};

/**
 * Class <b>TimeBlockSynchronized</b> is a time sentry for timing a statement
 * block.  The timer is generally started upon construction. But, the start is delayed
 * if the second argument is false.  In this case, manually start the timer by calling
 * the start() function.  This gives the safety of using a sentry, but does not force to
 * awkwardness associated with local variables crossing the timed block.
 *
 * Prior to starting the timer, an MPI synchronization barrier is set so that the
 * timing of routines which require MPI communication will all be at a known location
 * prior to executing.
 *
 */
class TimeBlockSynchronized
{
public:
  /**
   * Creates a new <b>TimeBlockSynchronized</b> instance.  If
   * <b>start_timer</b> is true, then the timer is started using the
   * <b>start()</b>.  An <b>MPI_Barrier</b> is called to synchronize the
   * start of the timer. The destructor will always stop a started timer.
   *
   * @param timer    a <b>Timer</b> reference to the timer to start.
   *
   * @param mpi_comm    a <b>MPI_Comm</b> value of the mpi communicator.

   * @param start_timer  a <b>bool</b> value to start the timer on construction.
   *
   */
  TimeBlockSynchronized(Timer &timer, ParallelMachine mpi_comm, bool start_timer = true);

  /**
   * Destroys a <b>TimeBlockSynchronized</b> instance.  Stops the timer if it has
   * been started.
   *
   */
  ~TimeBlockSynchronized();

  /**
   * Member function <b>start</b> starts the timer associated with the time block.
   * An <b>MPI_Barrier</b> is executed prior to starting the timer.
   *
   */
  void start();

  /**
   * Member function <b>stop</b> stops the timer associated with the time block.
   *
   */
  void stop();

private:
  Timer &      m_timer;  ////< Timer to accumulate block run times.
  ParallelMachine    m_mpiComm;  ////< MPI comm to synchronize across
  bool      m_started;  ////< Timer has been started
};


/**
 * @brief Function <b>operator<<</b> writes a timer to the diagnostic stream.
 *
 * @param dout      a <b>Writer</b> reference to the diagnostic writer to print
 *        to.
 *
 * @param timer      a <b>Timer::Metric</b> const reference to the timer
 *        to print.
 *
 * @return      a <b>Writer</b> reference to <b>dout</b>.
 */
template <class T>
inline Writer &operator<<(Writer &dout, const Timer::Metric<T> &timer) {
  return timer.dump(dout);
}

/**
 * Function <b>operator<<</b> writes a timer metric to the diagnostic stream.
 *
 * @param dout      a <b>Writer</b> reference to the diagnostic writer to print
 *        to.
 *
 * @param timer      a <b>Timer::Metric</b> const reference to the timer
 *        to print.
 *
 * @return      a <b>Writer</b> reference to <b>dout</b>.
 */
inline Writer &operator<<(Writer &dout, const Timer &timer) {
  return timer.dump(dout);
}

} // namespace diag
} // namespace stk


namespace sierra {
namespace Diag {

typedef stk::diag::Timer Timer;
typedef stk::diag::TimerSet TimerSet;
typedef stk::TimeFormat TimeFormat;
typedef stk::diag::TimeBlock TimeBlock;
typedef stk::diag::TimeBlockSynchronized TimeBlockSynchronized;

/**
 * @brief Enumeration <b><unnnamed></b> defines the bit mask values for the diagnostic
 * timer's in the <b>Diag</b> namespace.
 *
 */
enum TimerSetMask{
  TIMER_DOMAIN		= 0x00000001,		///< Enable domain timers
  TIMER_REGION		= 0x00000002,		///< Enable region timers
  TIMER_PROCEDURE	= 0x00000004,		///< Enable procedure timers
  TIMER_MECHANICS	= 0x00000008,		///< Enable mechanics timers
  TIMER_ALGORITHM	= 0x00000010,		///< Enable algorithm timers
  TIMER_SOLVER		= 0x00000020,		///< Enable solver timers
  TIMER_CONTACT		= 0x00000040,		///< Enable contact timers
  TIMER_MATERIAL	= 0x00000080,		///< Enable material timers
  TIMER_SEARCH		= 0x00000100,		///< Enable search timers
  TIMER_TRANSFER	= 0x00000200,		///< Enable transfer timers
  TIMER_ADAPTIVITY	= 0x00000400,		///< Enable adaptivity
  TIMER_RECOVERY	= 0x00000800,		///< Enable recovery
  TIMER_PROFILE_1	= 0x00001000,		///< Enable profile 1 timers
  TIMER_PROFILE_2	= 0x00002000,		///< Enable profile 2 timers
  TIMER_PROFILE_3	= 0x00004000,		///< Enable profile 3 timers
  TIMER_PROFILE_4	= 0x00008000,		///< Enable profile 4 timers
  TIMER_APP_1		= 0x00010000,		///< Enable application defined 1
  TIMER_APP_2		= 0x00020000,		///< Enable application defined 2
  TIMER_APP_3		= 0x00040000,		///< Enable application defined 3
  TIMER_APP_4		= 0x00080000,		///< Enable application defined 4
  TIMER_ALL		= 0x000FFFFF,		///< Enable all timers
  TIMER_NONE		= 0x00000000,		///< Enable no timers

  TIMER_FORCE		= 0x00000000		///< Force timer to be active
};


TimerSet &sierraTimerSet();

Timer &sierraTimer();

void sierraTimerDestroy();

class TimerParser;

/**
 * @brief Class <b>Timer</b> implements a diagnostic timer and timer container for the
 * collection and display of execution times.
 *
 */
typedef sierra::OptionMask TimerMask;			///< Bit mask for enabling timer

/**
 * @brief Enumeration <b><unnamed></b> defines some constants used for the
 * indented display of the diagnostic timers.
 *
 */
enum {
  DEFAULT_TIMER_NAME_MAX_WIDTH = 40			///< Width to truncate the name
};

/**
 * @brief Member function <b>theTimerParser</b> returns a reference to the timer
 * parser.
 *
 * @return			a <b>TimerParser</b> reference to the timer parser.
 */
TimerParser &theTimerParser();

/**
 * Function <b>setEnabledTimerMask</b> set the timer enable bit mask to
 * <b>timer_mask</b>.
 *
 * @param timer_mask		a <b>TimerMask</b> value to set the timer enable bit
 *				mask to.
 *
 */
void setEnabledTimerMask(TimerMask timer_mask);

/**
 * Function <b>getEnabledTimerMask</b> retruns the timer enable bit mask.
 *
 * @return			a <b>TimerMask</b> value of the timer enable bit
 *				mask.
 */
TimerMask getEnabledTimerMask();

void setTimeFormat(int time_format);

void setTimeFormatMillis();

int getTimeFormat();

/**
 * @brief Member function <b>setTimerNameMaxWidth</b> sets the maximum width for names
 * displayed in the timer output table.
 *
 * @param width		a <b>size_t</b> value to set for the maximum width for
 *				names displayed in the timer output table.
 *
 */
void setTimerNameMaxWidth(size_t width);

/**
 * @brief Member function <b>getTimeNameMaxWidth</b> returns the width to use for the
 * name cell of the table display.
 *
 * @return			a <b>size_t</b> value of the width to use for the name
 *				cell of the table display.
 */
size_t getTimerNameMaxWidth();

stk::diag::MetricTraits<stk::diag::CPUTime>::Type getCPULapTime(Timer timer);

stk::diag::MetricTraits<stk::diag::CPUTime>::Type getCPUAccumulatedLapTime(Timer timer);

stk::diag::MetricTraits<stk::diag::CPUTime>::Type getSierraCPUTime();
stk::diag::MetricTraits<stk::diag::CPUTime>::Type getSierraWallTime();


/**
 * @brief Class <b>TimerParser</b> implements the bit mask parser for the timer's bit masks.
 *
 */
class TimerParser : public OptionMaskParser
{
public:
  /**
   * Creates a new <b>TimerParser</b> instance.
   *
   */
  TimerParser();

  /**
   * @brief Member function <b>parse</b> parses the mask string and generates the
   * corresponding bit mask.
   *
   * @param mask_string		a <b>std::string</b> const reference to the mask string.
   *
   * @return			a <b>Mask</b> value of the bitmask corresponding to the
   *				mask string.
   */
  Mask parse(const char *mask_string) const;

  /**
   * Member function <b>parseArg</b> parses the argument and its argument values.
   *
   * @param name		a <b>std::string</b> const reference to the argument
   *				name.
   *
   * @param arg			a <b>std::string</b> const reference to the argument
   *				values.
   */
  virtual void parseArg(const std::string &name, const std::string &arg) const;  

  mutable stk::diag::MetricsMask        m_metricsSetMask;
  mutable stk::diag::MetricsMask        m_metricsMask;
};


class SierraRootTimer 
{
  public:
    SierraRootTimer();
    virtual ~SierraRootTimer();
    stk::diag::Timer & sierraTimer();

  private:
    stk::diag::Timer m_sierraTimer; 
};

} // namespace Diag
} // namespace sierra

///
/// @}
///

#endif // STK_UTIL_DIAG_Timer_hpp
