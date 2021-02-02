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

#ifndef STKUTIL_SCHEDULER_h
#define STKUTIL_SCHEDULER_h

#include <iosfwd>  // for ostream
#include <map>     // for map
#include <set>     // for set
#include <string>  // for string

namespace stk
{
namespace util
{

using Time = long double;
using Step = int;
using TimeContainer = std::map<double, double>;
using StepContainer = std::map<Step, Step>;

struct TolerancedTime {
    // A simple container to hold toleranced time values;
    double min;
    double max;
};

class Scheduler
  {
  public:
    Scheduler();

    Scheduler(const Scheduler &from);

    ~Scheduler();

    bool is_it_time(double time, Step step);

    bool add_interval(double time, double delta=0.0);
    bool add_interval(Time time, Time delta);
    bool add_interval(Step step, Step interval=1);

    bool set_termination_time(double time);

    bool add_explicit(double time);
    bool add_explicit(Step step);
    void set_force_schedule(); //!< Force true on next call to scheduler

    double adjust_dt(double dt, double time);
    bool set_lookahead(int lookahead);
    bool set_start_time(double time);
    void set_synchronize() {synchronize_ = true;}
    bool get_synchronize() {return synchronize_;}

    void set_restart_time(double time) { restartTime_ = time; }

    bool set_signal(const std::string& signal);

    /*! Function for use by encore.  They have a use case in which
     * they reset the time of a region back to zero multiple times
     * during an execution.  Normally, the Scheduler requires
     * that time be strictly increasing; this function resets the
     * lastTime_ value to -2.0 such that the Scheduler does not
     * throw an error about non-increasing time.
     */
    void reset_last_time();

    void print(std::ostream &out) const;

    const std::string & name() const { return name_; }
    void set_name(const std::string &n) { name_ = n; }

  private:

    bool internal_is_it_time(Time time);
    bool internal_is_it_step(Step step);
    bool force_schedule();
    bool add_explicit_internal(Time time);

    Time next_implicit_output_time(Time time) const;
    Time next_explicit_output_time(Time time) const;

    TimeContainer::const_iterator get_time_interval(Time time, bool erase_old) const;
    TolerancedTime get_toleranced_time_range(Time time) const;

    /*! timeIntervals_ is a container of time intervals where a
     * time interval is a pair <start_time, delta_time>.  A time
     * interval specifies that starting at time 'start_time',
     * output will be written every 'delta_time' intervals.  For
     * example, the interval <0.0, 0.1> would cause output to be
     * written at 0.0, 0.1, 0.2, 0.3, ....
     *
     * The stop or end time of an interval is the begin time of
     * the following interval.  Therefore, only a single time
     * interval is active at a single time and time intervals do
     * not overlap.
     */
    mutable TimeContainer timeIntervals_;

    //! See discussion for time intervals.
    StepContainer stepIntervals_;

    std::set<Time> times_;   ///<! List of addtional times at which output is desired.
    std::set<Step> steps_;   ///<! List of addtional steps at which output is desired.

    Time tolerance_;         ///<! tolerance for comparing times.
    mutable Time lastTime_;  ///<! Last time at which output was written.
    mutable Time firstTime_; ///<! First time at which output was written.
    Time lastCalledTime_;    ///<! Time routine called last; to estimate dt
    mutable Step lastInterval_;

    /*!
     * Number of simulation steps to adjust timestep
     * in order to hit output time exactly. If zero, don't adjust.
     * Used by adjust_dt to calculate timestep which will hit output time.
     */
    int lookAhead_;

    mutable Time startTime_; ///<! Used for backwards compatibility.

    /*! Time at which this scheduler ends.  A
     * is_it_time(time==terminationTime) call will return true, but
     * is_it_time(time >terminationTime) call will return false.
     */
    Time terminationTime_;

    /*! If this is a restarted run, the time at which it was
     *  restarted.  Used to make sure that output frequencies
     *  match between an original run and a restarted run.
     */
    Time restartTime_;

    bool forceSchedule_; ///<! Used to force scheduler to return true next time
    bool synchronize_;   ///<! Synchronize output with outputs from other regions in this procedure.

    mutable bool initialized_; ///<! True if this scheduler has been called.
    std::string name_;
  };
} // end namespace util
} // end namespace stk

#endif
