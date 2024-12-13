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

#include "stk_util/environment/Scheduler.hpp"

#include <float.h>  // for DBL_MAX, FLT_MAX

#include <cassert>  // for assert
#include <cmath>    // for ceil
#include <iomanip>
#include <iostream>  // for operator<<, basic_ostream::operator<<
#include <limits>    // for numeric_limits
#include <utility>   // for pair, swap

#include "stk_util/environment/EnvData.hpp"      // for EnvData
#include "stk_util/parallel/ParallelReduceBool.hpp"
#include "stk_util/util/Callback.hpp"            // for create_callback
#include "stk_util/util/SignalHandler.hpp"       // for SignalHandler

#ifndef TIME_MAX
#define TIME_MAX DBL_MAX
#endif

#ifndef TIME_EPSILON
#define TIME_EPSILON std::numeric_limits<Time>::epsilon()
#endif

namespace stk {
namespace util {

Scheduler::Scheduler() :
      tolerance_(1.0e-6),
      lastTime_(-TIME_MAX),
      firstTime_(-TIME_MAX),
      lastCalledTime_(-TIME_MAX),
      lastInterval_(-1),
      lookAhead_(0),
      startTime_(-TIME_MAX),
      terminationTime_(TIME_MAX),
      restartTime_(-TIME_MAX),
      forceSchedule_(false),
      synchronize_(false),
      initialized_(false) {}

Scheduler::~Scheduler() {}

Scheduler::Scheduler(const Scheduler &from) :
      timeIntervals_(from.timeIntervals_),
      stepIntervals_(from.stepIntervals_),
      times_(from.times_),
      steps_(from.steps_),
      tolerance_(from.tolerance_),
      lastTime_(from.lastTime_),
      firstTime_(from.firstTime_),
      lastCalledTime_(from.lastCalledTime_),
      lastInterval_(from.lastInterval_),
      lookAhead_(from.lookAhead_),
      startTime_(from.startTime_),
      terminationTime_(from.terminationTime_),
      restartTime_(from.restartTime_),
      forceSchedule_(from.forceSchedule_),
      synchronize_(from.synchronize_),
      initialized_(from.initialized_) {}

void Scheduler::reset_last_time()
{
  lastTime_ = -TIME_MAX;
  initialized_ = false;
}

bool Scheduler::internal_is_it_step(Step  step)
{
  assert(step >= 0);

  // See if this step is in the explicit steps_ list...
  if (steps_.find(step) != steps_.end())
    return true;

  // Find interval containing this step....
  // Intervals are sorted by step, so search for first
  // interval with step_num > step and then use previous (if any)
  if (stepIntervals_.empty()) {
    return false;
  }

  StepContainer::iterator begin = stepIntervals_.begin();
  StepContainer::iterator end   = stepIntervals_.end();

  // If output is not set to start until a ways into the simulation,
  // then the start step specified in 'begin' will be larger than the
  // current step.  In that case, return now.
  if ((*begin).first > step)
  {
    return false;
  }

  // Find 'next' interval.  If it is equal to end, then use begin interval,
  // If it is not equal to end, check the step it starts at and see if
  // it is now the active interval.  If it is, delete the first interval
  // and repeat the check on the next interval...
  StepContainer::iterator next = begin;
  ++next;

  while (next != end && (*next).first <= step) {
    stepIntervals_.erase(begin);
    begin = next++;
  }

  assert( (*begin).first <= step);
  assert(next == end || (*next).first > step);

  Step start = (*begin).first;
  Step inter = (*begin).second;
  assert(inter > 0);

  return (step - start) % inter == 0;
}

TolerancedTime Scheduler::get_toleranced_time_range(Time time) const
{
  TolerancedTime delta;
  if (firstTime_ > -TIME_MAX) {
    // In a transient analysis, the timesteps are typically very small
    // and if 'time' is large enough, multiple timesteps may fall
    // within the tolerance bounds about 'time' which will cause a
    // clustering of output steps around 'time'.  Typically, 'time' is
    // very small, so this won't happen, but in a coupled analysis
    // (e.g., quasistatic followed by transient), the times of the
    // transient analysis can be large and we will see this
    // problem. To solve this, we save the first time that this
    // Scheduler performed output and calculate the tolerance bounds
    // on the delta from this first time...
    delta.min = ((1.0 - tolerance_) * (time - firstTime_)) + firstTime_;
    delta.max = ((1.0 + tolerance_) * (time - firstTime_)) + firstTime_;
    if (delta.min == delta.max) {
      int i = 1;
      while (delta.min == delta.max) {
	delta.min = time - i*tolerance_;
	delta.max = time + i*tolerance_;
	i++;
      }
    }
  }
  else {
    delta.min = (1.0 - tolerance_) * time;
    delta.max = (1.0 + tolerance_) * time;
  }
  if (delta.max < delta.min) {
    // Dealing with negative times...
    std::swap(delta.min, delta.max);
  }
  return delta;
}

bool Scheduler::internal_is_it_time(Time time)
{
  // Returns true if:
  // 1. Parameter 'time' matches one of the times specified in the 'times_' list
  //    within a +/- tolerance.
  // 2. Parameter 'time' is greater than or equal to the next output
  //    time as specified in the timeIntervals_ list.
  // 3. Parameter 'time' is approximately equal to startTime_
  // 4. Parameter 'time' is approximately equal to terminationTime_
  // Also returns true if 'time' == 'lastTime_'.  This is to ensure that if
  // called multiple times with the same argument, it will return the
  // same response.

  STK_ThrowAssertMsg(time >= lastTime_, "time = " << time << ", lastTime_ = " << lastTime_);

  // If this is a restart, then calculate what the lastTime_ setting would
  // have been for this scheduler (based only on start time and deltas).
  if (!initialized_) {
    initialized_ = true;

    // If user has specified a start time, make sure none of the
    // "additional times" are previous to that time...
    if (startTime_ != -TIME_MAX && !times_.empty()) {
      std::set<Time>::iterator iter = times_.begin();
      while (iter != times_.end() && *iter < startTime_) {
	times_.erase(iter);
	iter = times_.begin();
      }
    }
      
    // This routine can be called from Region::initialize via a "will_output"
    // call prior to restart ocurring in which case restartTime_ will be > time.
    // In that case, return false and don't reset any "magic numbers"
    if (restartTime_ > time)
    {
      return false;
    }

    if (restartTime_ > -TIME_MAX) {
      assert(restartTime_ <= time);
      TimeContainer::const_iterator interval = get_time_interval(restartTime_, false);
      if (interval != timeIntervals_.end()) {
        Time start = (*interval).first;
        Time delta = (*interval).second;

        int inter = static_cast<int>((restartTime_ - start) / delta);
        if (inter >= 0) {
          firstTime_ = (static_cast<Time>(inter) * delta) + start;
          lastTime_  = restartTime_;
        }
      }
    }
  }

  // If called multiple times for same time, should return same
  // response each time.
  if (lastTime_ == time)
    return true;

  Time dt = time - lastCalledTime_;
  lastCalledTime_ = time;

  if (dt > 0 && dt < tolerance_) {
    tolerance_ = dt / 100.0L;
    if (tolerance_ < TIME_EPSILON) {
      tolerance_ = TIME_EPSILON;
    }
  }

  TolerancedTime delta = get_toleranced_time_range(time);

  // See if this time is in the explicit times_ list...
  // The 'times_' list is sorted and items are removed from
  // the list as they are used.  If the first item in the list
  // is less than or equal to the current time, then will do output.
  // Remove this item (and all others <= time) from the list.
  //
  // If list is empty or if first item is larger than time, don't
  // do output.
  {
    std::set<Time>::iterator iter = times_.begin();
    if (iter != times_.end()) {
      if (*iter <= delta.max) {
        while (iter != times_.end() && *iter <= delta.max) {
          times_.erase(iter);
          iter = times_.begin();
        }
        lastTime_ = time;
        if (firstTime_ == -TIME_MAX) firstTime_ = time;
        lastInterval_ = -1;
        return true;
      }
    }
  }

  if (delta.max < startTime_) {
    return false;
  }

  // Output once if time >= terminationTime_ and then turn off this
  // output stream.
  if (time >= terminationTime_ && terminationTime_ > -TIME_MAX) {
    if (lastTime_ < terminationTime_) {
      lastTime_ = time;
      if (firstTime_ == -TIME_MAX) firstTime_ = time;
      lastInterval_ = -1;
      return true;
    } else {
      return false;
    }
  }

  if (delta.min <= startTime_ && startTime_ <= delta.max) {
    lastTime_ = time;
    if (firstTime_ == -TIME_MAX) firstTime_ = time;
    lastInterval_ = -1;
    return true;
  }

  TimeContainer::const_iterator interval = get_time_interval(time, true);
  if (interval == timeIntervals_.end())
  {
    return false;
  }

  Time start  = (*interval).first;
  Time tdelta = (*interval).second;

  // If the time delta is equal to zero, then assume that
  // All times in this interval are to be printed...
  if (tdelta == 0.0 || (delta.min <= start && start <= delta.max)) {
    lastTime_ = time;
    if (firstTime_ == -TIME_MAX) firstTime_ = time;
    lastInterval_ = -1;
    return true;
  }

  // Calculate number of intervals needed to reach this time.
  int intervals = static_cast<int>((delta.max - start) / tdelta);

  if (lastInterval_ < 0 && lastTime_ > -TIME_MAX) {
    lastInterval_ = static_cast<int>((delta.min- start) / tdelta);
  }

  // If the last output time was in the same interval as the current time,
  // then don't output; otherwise, we want to output....
  if (intervals == lastInterval_) {
    // See if the time delta from the last write was > tdelta (Can
    // happen if time is large and tdelta small...)
    
    if ( time >= (lastTime_ + tdelta) ) {
      lastTime_ = time;
      lastInterval_ = intervals;
      if (firstTime_ == -TIME_MAX) firstTime_ = time;
      return true;
    } else {
      //std::cout << __FUNCTION__ << ", " << __LINE__ << " time= " << std::setprecision(16) << time << " >=  lastTime_ + tdelta= " << std::setprecision(16) << ( lastTime_ + tdelta) << std::endl; 
      return false;
    }
  } else {
    lastTime_ = time;
    lastInterval_ = intervals;
    if (firstTime_ == -TIME_MAX) firstTime_ = time;
    return true;
  }
}

void Scheduler::set_force_schedule()
{
  forceSchedule_ = true;
}

bool Scheduler::force_schedule()
{
  // It is possible that the forceSchedule_flag has been set on only one
  // processor so we need to see if it is true on any processor...
  bool result = forceSchedule_;
  if (EnvData::parallel_size() > 1) {
    result = stk::is_true_on_any_proc(EnvData::parallel_comm(), forceSchedule_);
  }
  forceSchedule_ = false;
  return result;
}

bool Scheduler::is_it_time(double t, Step step)
{
  Time time = t;
  // NOTE: It is possible that this routine is called multiple times
  // at the same time and step and it needs to return the same result
  // each time. To do this, it sets the lastTime_ variable whenever
  // the routine is run and returns with a "true" response.  Then, if
  // it is called again, it will compare time with lastTime_ and if
  // they match, return true again.

  // force_write always causes a write even if the time is outside
  // the bounds set by startTime and terminationTime...
  //    Needs to happen before the "time == lastime" check to make
  //    sure force flag is unset
  if (force_schedule()) {
    lastTime_ = time;
    return true;
  }

  // If called multiple times, return same response...
  if (time == lastTime_)
  {
    return true;
  }


  // If user specified a start time; that overrides
  // everything except for the force_write setting.
  // If user did not specify the start time, then startTime_ = -TIME_MAX
  // which doesn't force an output, but does enable the output
  // block...
  if (get_toleranced_time_range(time).max < startTime_)
  {
    return false;
  }

  // See if user specified an output time via time intervals
  // or explicit times.
  bool write = internal_is_it_time(time);

  if (!write) {
    // Check whether the 'step' intervals specify a write at this time...
    write = internal_is_it_step(step);

    if (write) {
      // Need to update lastTime and see if past terminationTime...
      // This is already done in the "internal_is_it_time(time)" routine.
      if (time >= terminationTime_ && lastTime_ >= terminationTime_) {
        write = false;
      }
    }
  }

  if (!write) {
    // Catch some strange cases here...
    // If the code has not specified a startTime, or a terminationTime, or
    // any step or time intervals..., Then assume that it is a non-transient
    // one-time call to this function and we need to output at this time...
    if (terminationTime_ == TIME_MAX &&
        startTime_ == -TIME_MAX &&
        timeIntervals_.empty() &&
        stepIntervals_.empty() &&
        times_.empty() &&
        steps_.empty()) {
      write = true;
    }
  }

  // This should be set already, but do it here again in case
  if (write)
    lastTime_ = time;

  return write;
}

/*!
 * Determine first output time >= 'time'
 * Steps to reach that time at current dt is: nstep = (time_out - time) / dt
 * [Account for overflow or extremely small dt....]
 * If 'nstep' > lookahead, return current dt,
 * else, new_dt = (time_out - time) / static_cast<int>(nstep)
 *
 * \pre dt > 0.0
 * \post 0.0 < returned_dt <= dt
 * ...
 */
double Scheduler::adjust_dt(double dt, double time)
{
  assert(dt > 0.0);

  if (lookAhead_ == 0)
  {
    return dt;
  }

  // See if this output block is still active...
  if (time >= terminationTime_)
  {
    return dt;
  }

  double next_imp = next_implicit_output_time(time);
  double next_exp = next_explicit_output_time(time);
  double next = next_imp < next_exp ? next_imp : next_exp;

  double delta = next - time;
  // Some codes call this routine prior to outputting the current step
  // In that case, we have already hit the desired output time and
  // don't need to adjust the dt
  if (delta < tolerance_ * time) {
    return dt;
  }

  TolerancedTime delta_range = get_toleranced_time_range(delta);
  const double steps = (dt >= delta_range.min) ? 1.0 : ceil(delta / dt);
  assert(steps > 0);

  // If 'steps' is less than 'lookAhead', then calculate
  // a new dt which will hit the time exactly...
  if (steps <= lookAhead_) {
    double new_dt = delta / steps;

    // Check for truncation errors....
    double proj_time = time + (new_dt * steps);
    if (proj_time < next) {
      new_dt *= (1.0 + TIME_EPSILON);
    }
    return new_dt;
  }
  return dt;
}

Time Scheduler::next_explicit_output_time(Time time) const
{
  // Return next output time greater than the passed in 'time'
  // from the list of 'also output at' times.
  // Return terminationTime_ if none.

  // Set up a tolerance range around the current time.
  TolerancedTime delta = get_toleranced_time_range(time);

  // See if this time is in the explicit times_ list...
  // List is sorted, so find first time > delta.max?
  std::set<Time>::const_iterator iter = times_.begin();
  std::set<Time>::const_iterator end  = times_.end();
  while (iter != end) {
    if (*iter > delta.max && *iter > startTime_ && *iter <= terminationTime_) {
      return *iter;
    }
    ++iter;
  }
  return terminationTime_;
}

Time Scheduler::next_implicit_output_time(Time time) const
{
  // [Very similar calculations to the 'is_it_time' function...Consolidate]

  if (!initialized_) {
    initialized_ = true;
    if (restartTime_ > -TIME_MAX) {
      assert(restartTime_ <= time);
      TimeContainer::const_iterator interval = get_time_interval(restartTime_, false);
      if (interval != timeIntervals_.end()) {
        Time start = (*interval).first;
        Time delta = (*interval).second;

        int inter = static_cast<int>((restartTime_ - start) / delta);
        if (inter >= 0) {
          firstTime_ = (static_cast<Time>(inter) * delta) + start;
          lastTime_  = restartTime_;
        }
      }
    }
  }

  // Return next output time after the passed in time.
  TolerancedTime delta = get_toleranced_time_range(time);
  if (delta.min < startTime_ && time > lastTime_)
    return startTime_;

  // Find interval containing this time....
  // Intervals are sorted by time, so search for first
  // interval with time_val > time and then use previous (if any)

  TimeContainer::const_iterator begin = timeIntervals_.begin();
  TimeContainer::const_iterator end   = timeIntervals_.end();
  if (begin == end)
    return FLT_MAX; // Empty

  // If output is not set to start until a ways into the simulation,
  // then the start time specified in 'begin' will be larger than the
  // current time.  Return that time.
  if ((*begin).first > time)
    return (*begin).first;

  // Find 'next' interval.  If it is equal to end, then use begin
  // interval, If it is not equal to end, check the time it starts at
  // and see if it is now the active interval.  If it is, repeat the
  // check on the next interval...
  TimeContainer::const_iterator next = begin;
  ++next;

  while (next != end && (*next).first <= delta.max) {
    begin = next++;
  }

#ifndef NDEBUG
  static constexpr double eps = std::numeric_limits<double>::epsilon();
  STK_ThrowAssertMsg( (*begin).first <= (time + eps), "begin.first="<<(*begin).first<<", time="<<time<<" diff "<<((*begin).first-time)<<" eps "<<std::numeric_limits<double>::epsilon());
#endif

  assert(next == end || (*next).first > delta.min);

  // At this point, have the correct time interval which specifies start and frequency
  Time start = (*begin).first;
  Time tdelta = (*begin).second;

  // If the time delta is equal to zero, then assume that
  // All times in this interval are to be printed...
  if (tdelta == 0.0) {
    return delta.max;
  }

  // Calculate number of intervals needed to reach this time.
  int intervals = static_cast<int>((delta.max - start) / tdelta);

  // Perform same calculation for the previous output time...
  int prev_inter = static_cast<int>((lastTime_ - start + tolerance_) / tdelta);

  // If the last output time was in the same interval as the current time,
  // then don't output; otherwise, we want to output....
  if (intervals == 0 || intervals == prev_inter) {
    intervals++;
  }

  // Calculate next output time....
  Time next_time = start + static_cast<Time>(intervals) * tdelta;
  Time next_interval_start = next != end ? (*next).first : next_time;
  if (next_interval_start < next_time) next_time = next_interval_start;
  if (get_toleranced_time_range(time).max > terminationTime_) {
    next_time = terminationTime_;
  }
  // Some codes call this routine prior to outputting the current step
  // In that case, we have already hit the desired output time and
  // don't need to adjust the dt. If this is the case, then next_time == time.
  //
  //assert(next_time >= time);
  //
  //if( next_time >= time )
  //{
  //  std::cout << std::setprecision(15) << "TRUE: next_time= " << next_time << ", " << std::setprecision(15) << time << std::endl;
  //}
  //else
  //{
  //  std::cout << std::setprecision(15) << "FALSE: next_time= " << next_time << ", " << std::setprecision(15) << time << std::endl;
  //}

  return next_time;
}

bool Scheduler::add_interval(Step step, Step interval)
{
  std::pair<StepContainer::iterator,bool> result;
  result=stepIntervals_.insert(StepContainer::value_type(step, interval));
  return result.second;
}

bool Scheduler::add_interval(double time, double delta)
{
  Time t = time;
  Time d = delta;
  return add_interval(t, d);
}

bool Scheduler::add_explicit(Step step)
{
  return steps_.insert(step).second;
}

bool Scheduler::add_explicit(double time)
{
  Time t = time;
  return add_explicit_internal(t);
}

bool Scheduler::add_interval(Time time, Time delta)
{
  // Adjust tolerance_ to be at least 3 orders of magnitude smaller
  // than delta... [3 orders of magnitude is arbitrary rule of thumb]
  // Note that this may cause problems if there are time intervals
  // with widely varying interval sizes...
  const Time scale_factor = 1000.0L;
  if (scale_factor * tolerance_ > delta) {
    tolerance_ = delta / scale_factor;
    if (tolerance_ < TIME_EPSILON) {
      tolerance_ = TIME_EPSILON;
    }
  }

  std::pair<TimeContainer::iterator,bool> result;
  result=timeIntervals_.insert(TimeContainer::value_type(time, delta));
  return result.second;
}

bool Scheduler::add_explicit_internal(Time time)
{
  return times_.insert(time).second;
}

bool Scheduler::set_lookahead(int lookahead)
{
  if (lookahead >= 0) {
    lookAhead_ = lookahead;
    return true;
  }
  return false;
}

bool Scheduler::set_start_time(double time)
{
  // This is a backwards compatibility function used when there was
  // only the capability to set the time interval and a default start
  // time of 0.0.  Now there are multiple time intervals each with a
  // start time. The function currently sets the variable startTime_
  // which is checked in a few routines.  The alternative which is to
  // modify the start time of the first time interval needs extra work
  // such as command ordering....
  startTime_ = time;
  return true;
}

bool Scheduler::set_termination_time(double time)
{
  if (terminationTime_ > time && time != -TIME_MAX)
  {
    terminationTime_ = time;
  }
  return true;
}

bool Scheduler::set_signal(const std::string& signal)
{
  // The signal that this scheduler should handle is passed in as a
  // string. This will be a valid signal name, but may not be handled
  // on this particular system.  The
  // SignalHandler::check_signal_name() function will determine
  // whether this is a valid signal.

  bool success = false;
  if (sierra::SignalHandler::instance().check_signal_name(signal.c_str())) {
    sierra::SignalHandler::instance().add_handler(signal.c_str(),
        *sierra::create_callback(*this, &Scheduler::set_force_schedule));
    success = true;
  } else {
  }
  return success;
}

// ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

void Scheduler::print(std::ostream &out) const
{
  out << "\nDump of Scheduler Object:\n";
  {
    out << "Step Intervals:\n";
    StepContainer::const_iterator iter = stepIntervals_.begin();
    StepContainer::const_iterator end  = stepIntervals_.end();
    if (iter != end) {
      while (iter != end) {
        out << "\tStarting at step " << (*iter).first << ", output every "
            << (*iter).second << " steps.\n";
        ++iter;
      }
    } else {
      out << "\t<EMPTY>\n";
    }
  }

  {
    out << "\nExplicit Steps:\n";
    std::set<Step>::const_iterator iter = steps_.begin();
    std::set<Step>::const_iterator end  = steps_.end();
    if (iter != end) {
      while (iter != end) {
        out << "\t" << (*iter++) << ", ";
      }
      out << "\n";
    } else {
      out << "\t<EMPTY>\n";
    }
  }

  {
    out << "\nTime Intervals:\n";
    TimeContainer::const_iterator iter = timeIntervals_.begin();
    TimeContainer::const_iterator end  = timeIntervals_.end();
    if (iter != end) {
      while (iter != end) {
        out << "\tStarting at time " << (*iter).first << ", output every "
            << (*iter).second << ".\n";
        ++iter;
      }
    } else {
      out << "\t<EMPTY>\n";
    }
  }

  {
    out << "\nExplicit Times:\n";
    std::set<Time>::const_iterator iter = times_.begin();
    std::set<Time>::const_iterator end  = times_.end();
    if (iter != end) {
      while (iter != end) {
        out << "\t" << (*iter++) << ", ";
      }
      out << "\n";
    } else {
      out << "\t<EMPTY>\n";
    }
  }

  {
    out << "\nOther Data:\n";
    out << "\tScheduler tolerance = " << tolerance_ << "\n";
    out << "\tLook Ahead = " << lookAhead_ << "\n";
    out << "\tStart Output at Time " << startTime_ << "\n";
    out << "\tStop  Output at Time " << terminationTime_ << "\n\n";
  }
}

TimeContainer::const_iterator Scheduler::get_time_interval(Time time, bool erase_old) const
{
  // Find interval containing this time....
  // Intervals are sorted by time, so search for first
  // interval with time_val > time and then use previous (if any)

  // Set up a tolerance range around the current time.
  TolerancedTime delta = get_toleranced_time_range(time);
  TimeContainer::iterator begin = timeIntervals_.begin();
  TimeContainer::iterator end   = timeIntervals_.end();
  if (begin == end)
    return end; // Empty

  // If output is not set to start until a ways into the simulation,
  // then the start time specified in 'begin' will be larger than the
  // current time.  In that case, return now.
  if ((*begin).first > delta.max)
    return end;

  // Find 'next' interval.  If it is equal to end, then use begin interval,
  // If it is not equal to end, check the time it starts at and see if
  // it is now the active interval.  If it is, delete the first interval
  // and repeat the check on the next interval...
  TimeContainer::iterator next = begin;
  ++next;

  while (next != end && (*next).first <= delta.max) {
    if (erase_old) {
      timeIntervals_.erase(begin);
      lastTime_ = -TIME_MAX; // Reset since in a new interval
    }
    begin = next++;
  }

  assert( (*begin).first <= delta.max);
  assert(next == end || (*next).first > delta.min);

  // At this point, have the correct time interval which specifies start and frequency
  return begin;
}

} // end namespace util
} // end namespace stk
