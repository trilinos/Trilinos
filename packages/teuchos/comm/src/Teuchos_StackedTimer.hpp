// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_STACKED_TIMER_HPP
#define TEUCHOS_STACKED_TIMER_HPP

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_PerformanceMonitorBase.hpp"
#include <string>
#include <vector>
#include <cassert>
#include <chrono>
#include <climits>
#include <cstdlib> // for std::getenv and atoi
#include <ctime> // for timestamp support
#include <iostream>

#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
namespace Kokkos {
namespace Tools {
extern void pushRegion (const std::string&);
extern void popRegion ();
} // namespace Profiling
} // namespace Kokkos
#endif


namespace Teuchos {

//! Error reporting function for stacked timer.
TEUCHOSCOMM_LIB_DLL_EXPORT void error_out(const std::string& msg, const bool fail_all = false);

/**
 * \brief The basic timer used internally, uses std::chrono::high_resolution_clock.
 *
 * This class hold  a time and number of times this timer is called, and a
 * count of how many "updates" this timer has services.  Example if you have a
 * mesh of 100 elements your count might be 7, but total iterations would be 700
 * this number is useful for dynamic systems where the total number of items services
 * might change in time
 */
class BaseTimer {

public:

  using Clock = std::chrono::high_resolution_clock;

  BaseTimer() : accumulation_(0.0), count_started_(0), count_updates_(0), running_(false) {}

  /// Start a currently stopped timer
  void start(){
    if (running_)
      error_out("Base_Timer:start Failed timer already running");
    start_time_ = Clock::now();

    count_started_++;
    running_ = true;
  }

  /// Stop a current running timer and accumulate time difference
  void stop(){
    if (!running_)
      error_out("Base_Timer:stop Failed timer not running");
    accumulation_ += std::chrono::duration_cast<std::chrono::duration<double>>(Clock::now() - start_time_).count();
    running_ = false;
  }

  /// Increment the total number of items updated between a start stop
  unsigned long long incrementUpdates(unsigned long long count=1) {count_updates_ += count; return count_updates_;}

  /// Get the total accumulated time since last reset or construction when the timer is running
  double accumulatedTime() const {return accumulation_;}

  /// Setter for accumulated time
  void setAccumulatedTime(double accum=0)  {accumulation_=accum;}

  /**
   * \brief return the average time per item updated
   *
   * This returns the time on average that the code spends updating an
   * iteration.  If it is running than it will not include the current time.  It
   * differs from accumulatedTimePerTimerCall in that it is meant to be timer per
   * event other that start/stop, like mesh update
   * @return average time per iteration pair
   */  double accumulatedTimePerUpdate() const {
    if (count_updates_ > 0) {
      return accumulation_/count_updates_;
    } else {
      return 0;
    }
  }


  /**
   * \brief return the average time per timer start/stop
   *
   * This returns the time on average that the code spends between
   * a call to start and stop.  If it is running than it will not include the current time
   * @return average time per start/stop pair
   */
  double accumulatedTimePerTimerCall() const {
    if (count_started_> 0) {
      return accumulation_/count_started_;
    } else {
      return 0;
    }
  }

  /**
   * \brief  Return the difference between two timers in seconds,
   * @param [in] from reference time you are computing difference from
   * @return my time - from times
   */

  double difference(const BaseTimer &from) const {
    return accumulation_ - from.accumulation_;
  }

  /// Reset all the timer stats, throws if it is already running
  void reset() {
    if (running_)
       error_out("BaseTimer, cannot reset a running timer");
    accumulation_ = 0.0;
    count_started_ = count_updates_ = 0;
  }

  /// Returns true if the timer is currently accumulating time.
  bool running() const { return running_;}

  /// Returns the number of calls to start().
  unsigned long numCalls() const { return count_started_; }

  /// Returns the number of updates added to this timer.
  unsigned long long numUpdates() const { return count_updates_; }

  /// Sets the number of calls to start() for this timer. This is only used for unit testing.
  void overrideNumCallsForUnitTesting(const unsigned long num_calls)
  { count_started_ = num_calls; }

  /// Sets the number of counts for this timer. This is only used for unit testing.
  void overrideNumUpdatesForUnitTesting(const unsigned long long num_updates)
  { count_updates_ = num_updates; }

  struct TimeInfo {
    TimeInfo():time(0.0), count(0), updates(0), running(false){}
    TimeInfo(BaseTimer* t): time(t->accumulation_), count(t->count_started_), updates(t->count_updates_), running(t->running()) {}
    double time;
    unsigned long count;
    unsigned long long updates;
    bool running;
  };

protected:
  double accumulation_;       // total time
  unsigned long count_started_; // Number of times this timer has been started
  unsigned long long count_updates_; // Total count of items updated during this timer
  Clock::time_point start_time_;
  bool running_;

  friend struct TimeInfo;
};

/**
 *\brief This class allows one to push and pop timers on and off a stack.
 *
 * To use this class one would do something like this
 *
 * StackedTimer timer(); // construct top level timer
 * timer.start("some comp");
 * do_work();
 * timer.start("sub comp");
 * more_work();
 * timer.stop(); // stopping sub comp
 * timer.stop(); // stopping comp
 *
 * timer.stop(); // stopping all timer stuff
 * timer.report(std::cout); // dump to screen
 *
 */
class TEUCHOSCOMM_LIB_DLL_EXPORT StackedTimer
{
protected:

  /**
   * \brief Timer info at a given level and all the children
   *
   * This holds the timer info for the leaf node of a timer stack.  It has
   * both its timer info and all the sub leaves bellow it.  You can start and
   * stop sub timers, get sub info and dump out things to a ostream
   *
   */
  class LevelTimer : public BaseTimer {
  protected:

    // TODO: implement operator=

    unsigned level_;
    std::string name_;
    LevelTimer *parent_;
    std::vector<LevelTimer> sub_timers_;
  public:
    /// Default constructor, shouldn't be used but needed for std::vector
    LevelTimer();

    /**
     * Standard constructor
     * @param [in] level Integer level of this timer
     * @param [in] name Name for this timer
     * @param [in] parent Parent of this timer
     * @param [in] start bool to start the timer on construction
     */
    LevelTimer(int level,
        const char* name = "RootTimer",
        LevelTimer *parent=nullptr,
        bool start_timer=true) :
          BaseTimer(),
          level_(level),
          name_(name),
          parent_(parent)
    {
      if ( start_timer )
        BaseTimer::start();

    }

    /// Copy constructor
    LevelTimer(const LevelTimer &src) :
      BaseTimer(src), level_(src.level_), name_(src.name_),parent_(src.parent_), sub_timers_(src.sub_timers_)
    {
      for (unsigned i=0;i<sub_timers_.size();++i)
        sub_timers_[i].parent_ = this;
    }

    /**
     * Start a sub timer of a given name, create if doesn't exist
     * @param [in] sub_name Name of subtimer
     * @return Pointer to the sub timer which was started
     */
    LevelTimer* start(const char* sub_name) {
      for (unsigned i=0;i<sub_timers_.size();i++ )
        if (sub_name == sub_timers_[i].name_ ) {
          sub_timers_[i].BaseTimer::start();
          return &sub_timers_[i];
        }
      sub_timers_.push_back(LevelTimer(level_+1,sub_name,this,true));
      return &sub_timers_[sub_timers_.size()-1];
    }

    /**
     * Stop the current running timer.  The timer name is required to verify that you are
     * stopping the current running timer and will error out if the names don't match
     *
     * @param [in] name the name of the timer you want to stop, used for matching start only
     * @return parent level timer
     */
    LevelTimer* stop(const std::string &name = "RootTimer") {
      if (name != name_)
        error_out("Stopping timer "+name+" But top level running timer is "+name_);
      BaseTimer::stop();
      return parent_;
    }


    /**
     * Return the full name of the timer with each level split by :
     * @return The full name of the timer
     */
    std::string get_full_name() const {
      std::string parent_name("");
      if ((parent_ != nullptr))
        parent_name = parent_->get_full_name() + "@";

      std::string my_name(name_);

      std::string full_name = parent_name + my_name;
      return full_name;
    }

    /**
     * Return the number of timers on this level
     * @return the number of timers and sub timers
     */

    int countTimers() {
      int count=1;
      for (unsigned i=0;i<sub_timers_.size(); ++i)
        count += sub_timers_[i].countTimers();
      return count;
    }

    void addTimerNames(Array<std::string> &names, unsigned &pos) {
      names[pos++] = get_full_name();
      for (unsigned i=0;i<sub_timers_.size(); ++i)
        sub_timers_[i].addTimerNames(names, pos);
    }

    /**
     * Return the time spent at a given level
     * @param [in] locate_name name of subtimer, if blank return current level time
     * @return time in seconds at provided level
     */
    double accumulatedTime(const std::string &locate_name="") {

      if (locate_name == "")
        return BaseTimer::accumulatedTime();

      std::string first_name,second_name;

      size_t i = locate_name.find_first_of('@');
      if ( i >= locate_name.size() ) {
        first_name = locate_name;
        second_name = "";
      } else {
        first_name.assign(locate_name,0,i);
        second_name.assign(locate_name,i+1,locate_name.size()-i-1);
      }
      for (unsigned j=0;j<sub_timers_.size();++j)
        if ( first_name == sub_timers_[j].name_)
          return sub_timers_[j].accumulatedTime(second_name);
      return 0;
    }

    /// Returns the level of the timer in the stack
    unsigned level() const
    {return level_;}

  protected:
    /**
     * \brief split a string into two parts split by a '@' if no '@' first gets the full string
     * @param [in] locate_name input string to split
     * @param [out] first_name Part of string before the first '@'
     * @param [out] second_name part of string after the first '@'
     */
    void splitString(const std::string &locate_name, std::string &first_name, std::string &second_name) {
      size_t i = locate_name.find_first_of('@');
      if ( i >= locate_name.size() ) {
        first_name = locate_name;
        second_name = "";
      } else {
        first_name.assign(locate_name,0,i);
        second_name.assign(locate_name,i+1,locate_name.size()-i-1);
      }
    }

  public:
    /**
     * Return the time spent per update at a given level
     * @param [in] locate_name name of subtimer, if blank return current level time
     * @return time in seconds per update at provided level
     */
    double accumulatedTimePerUpdate(const std::string &locate_name="") {

       if (locate_name == "")
         return BaseTimer::accumulatedTimePerUpdate();

       std::string first_name,second_name;
       splitString(locate_name, first_name, second_name);

       for (unsigned j=0;j<sub_timers_.size();j++)
         if ( first_name == sub_timers_[j].name_)
           return sub_timers_[j].accumulatedTimePerUpdate(second_name);
       return 0;
     }

     /**
      * Return the time spent per timer start/stop pair at a given level
      * @param [in] locate_name name of subtimer, if blank return current level time
      * @return time in seconds per timer start/stop pair at provided level
      */
     double accumulatedTimePerTimerCall(const std::string &locate_name="") {

       if (locate_name == "")
         return BaseTimer::accumulatedTimePerTimerCall();

       std::string first_name,second_name;
       splitString(locate_name, first_name, second_name);

       for (unsigned j=0;j<sub_timers_.size();j++)
         if ( first_name == sub_timers_[j].name_)
           return sub_timers_[j].accumulatedTimePerTimerCall(second_name);
       return 0;
     }

     /**
      * Pack up all the levels into a buffer and mpi send to rank 0
      */
     void pack();

     /**
      * Unpack the level timer stack from a mpi recv
      * @param [in] from rank you are sending from
      * @return pointer to level timer unpacked
      */
     LevelTimer* unpack(unsigned from);

     /**
      * Dump the timer stats in a pretty format to ostream
      * @param [in,out] os  Where are you dumping the stats, stdout??
      */
     void report(std::ostream &os);

    /**
     * Return pointer to the BaseTimer corresponding to a given string
     * @param name input string to search for
     * @return pointer to BaseTimer (nullptr if none found)
     */
    const BaseTimer* findBaseTimer(const std::string &name) const;
    
     /**
      * Return the time info for a given string
      * @param name input string to search for
      * @param set to true on exit if timer was found
      * @return Time data
      */
    BaseTimer::TimeInfo findTimer(const std::string &name,bool& found);

  }; // LevelTimer


public:
   /**
    * Construct a stacked timer
    * @param [in] name Top level name of the timer
    * @param [in] start_top_timer Automatically start the top level timer. If set to false, the user will have to start it manually.
    */
  explicit StackedTimer(const char *name, const bool start_base_timer = true)
    : timer_(0,name,nullptr,false),
      global_mpi_aggregation_called_(false),
      enable_verbose_(false),
      verbose_timestamp_levels_(0), // 0 disables
      verbose_ostream_(Teuchos::rcpFromRef(std::cout)),
      enable_timers_(true)
  {
    top_ = &timer_;
    if (start_base_timer)
      this->startBaseTimer();

    auto check_verbose = std::getenv("TEUCHOS_ENABLE_VERBOSE_TIMERS");
    if (check_verbose != nullptr)
      enable_verbose_ = true;

    auto check_timestamp = std::getenv("TEUCHOS_ENABLE_VERBOSE_TIMESTAMP_LEVELS");
    if (check_timestamp != nullptr) {
      verbose_timestamp_levels_ = std::atoi(check_timestamp);
    }
  }

  /**
   * Start the base level timer only
   */
  void startBaseTimer() {
    timer_.BaseTimer::start();
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
    ::Kokkos::Tools::pushRegion(timer_.get_full_name());
#endif
  }

  /**
   * Stop the base level timer only
   */
  void stopBaseTimer() {
    timer_.BaseTimer::stop();
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
    ::Kokkos::Tools::popRegion();
#endif
  }

  /**
   * Start a sublevel timer
   * @param [in] name Name of the timer you wish to start
   * @param [in] push_kokkos_profiling_region Optional parameter that if set to true, will pushRegion() in kokkos profiling for this timer. The TimeMonitor will always set this to false since it does its own pushRegion() in the Timer object (this prevents double registering with kokkos).
   */
  void start(const std::string name,
             const bool push_kokkos_profiling_region = true) {
    if (enable_timers_) {
      if (top_ == nullptr)
        top_ = timer_.start(name.c_str());
      else
        top_ = top_->start(name.c_str());
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
      if (push_kokkos_profiling_region) {
        ::Kokkos::Tools::pushRegion(name);
      }
#endif
    }
    if (enable_verbose_) {
      if (!verbose_timestamp_levels_) {
        *verbose_ostream_ << "STARTING: " << name << std::endl;
      }
      // gcc 4.X is incomplete in c++11 standard - missing
      // std::put_time. We'll disable this feature for gcc 4.
#if !defined(__GNUC__) || ( defined(__GNUC__) && (__GNUC__ > 4) )
      else if (top_ != nullptr) {
        if ( top_->level() <= verbose_timestamp_levels_) {
          auto now = std::chrono::system_clock::now();
          auto now_time = std::chrono::system_clock::to_time_t(now);
          auto gmt = gmtime(&now_time);
          auto timestamp = std::put_time(gmt, "%Y-%m-%d %H:%M:%S");
          auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
          *verbose_ostream_ << "STARTING: " << name << " LEVEL: " << top_->level() << " COUNT: " << timer_.numCalls() << " TIMESTAMP: " << timestamp << "." << ms.count() << std::endl;
        }
      }
#endif
    }
  }

  /**
   * Stop the current top running timer, should be called for the root timer prior to final output
   * @param [in] name Name of the timer you wish to stop
   * @param [in] pop_kokkos_profiling_region Optional parameter that if set to true, will popRegion() in kokkos profiling for this timer. The TimeMonitor will always set this to false since it does its own pushRegion() in the Timer object (this prevents double registering with kokkos).
   */
  void stop(const std::string &name,
            const bool pop_kokkos_profiling_region = true) {
    if (enable_timers_) {
      if (top_)
        top_ = top_->stop(name);
      else
        timer_.BaseTimer::stop();
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
      if (pop_kokkos_profiling_region) {
        ::Kokkos::Tools::popRegion();
      }
#endif
    }
    if (enable_verbose_) {
      if (!verbose_timestamp_levels_) {
        *verbose_ostream_ << "STOPPING: " << name << std::endl;
      }
      // gcc 4.X is incomplete in c++11 standard - missing
      // std::put_time. We'll disable this feature for gcc 4.
#if !defined(__GNUC__) || ( defined(__GNUC__) && (__GNUC__ > 4) )
      // The stop adjusts the top level, need to adjust by +1 for printing
      else if (top_ != nullptr) {
        if ( (top_->level()+1) <= verbose_timestamp_levels_) {
          auto now = std::chrono::system_clock::now();
          auto now_time = std::chrono::system_clock::to_time_t(now);
          auto gmt = gmtime(&now_time);
          auto timestamp = std::put_time(gmt, "%Y-%m-%d %H:%M:%S");
          auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
          *verbose_ostream_ << "STOPPING: " << name << " LEVEL: " << top_->level()+1 << " COUNT: " << timer_.numCalls() << " TIMESTAMP: " << timestamp << "." << ms.count() << std::endl;
        }
      }
#endif
    }
  }

  /**
   * Increment the iteration count for the running timer
   *   @param [in] i amount to increment the count
   */
  void incrementUpdates(const long long i = 1) {
    top_->incrementUpdates(i);
  }

  /**
   * Return the accumulated time the named timer has, if name is blank top level timer
   * @param [in] name Name of the timer to output
   * @return amount of time in seconds
   */
  double accumulatedTime(const std::string &name="") {
    if (top_) // Top is null for the head node when nothing is running
      return top_->accumulatedTime(name);
    else
      return timer_.accumulatedTime(name);
  }

  /**
   * Return the time spent per update at a given level
   * @param [in] name name of subtimer, if blank return current level time
   * @return time in seconds per update at provided level
   */
  double accumulatedTimePerUpdate(const std::string &name="") {
     if (top_) // Top is null for the head node when nothing is running
       return top_->accumulatedTimePerUpdate(name);
     else
       return timer_.accumulatedTimePerUpdate(name);
   }
  /**
   * Return the time spent per timer start/stop at a given level
   * @param [in] name name of subtimer, if blank return current level time
   * @return time in seconds per start/stop call at provided level
   */
  double accumulatedTimePerTimerCall(const std::string &name="") {
     if (top_) // Top is null for the head node when nothing is running
       return top_->accumulatedTimePerTimerCall(name);
     else
       return timer_.accumulatedTimePerTimerCall(name);
   }
  
  /**
   * Return pointer to the BaseTimer corresponding to a given string (full string name)
   * @param name input string to search for
   * @return BaseTimer
   */
  const BaseTimer* findBaseTimer(const std::string &name) const {
    const BaseTimer* baseTimer = timer_.findBaseTimer(name);
    TEUCHOS_TEST_FOR_EXCEPTION(baseTimer == nullptr, std::runtime_error,
                               "StackedTimer::findBaseTimer() failed to find a timer named \"" << name << "\"!\n");
    return baseTimer;
  }

  /**
   * Return the time info for a given string (full string name)
   * @param name input string to search for
   * @return Time data
   */
  BaseTimer::TimeInfo findTimer(const std::string &name) {
    bool foundTimer = false;
    const auto timeInfo = timer_.findTimer(name,foundTimer);
    TEUCHOS_TEST_FOR_EXCEPTION(!foundTimer, std::runtime_error,
      "StackedTimer::findTimer() failed to find a timer named \"" << name << "\"!\n");
    return timeInfo;
  }

  void report(std::ostream &os) {
    timer_.report(os);
  }

  /** Struct for controlling output options like histograms

      @param output_fraction Print the timer fractions within a level.
      @param output_total_updates Print the updates counter.
      @param output_historgram Print the histogram.
      @param output_minmax Print the min max and standard deviation across MPI processes.
      @param num_histogram The number of equally size bickets to use in the histogram.
      @param max_level The number of levels in the stacked timer to print (default prints all levels).
      @param print_warnings Print any relevant warnings on stacked timer use.
      @param align_columns Output will align the columsn of stacked timer data.
      @param print_names_before_values If set to true, writes the timer names before values.
      @param drop_time If a timer has a total time less that this value, the timer will not be printed and the total time of that timer will be added to the Remainder. Useful for ignoring negligible timers. Default is -1.0 to force printing of all timers even if they have zero accumulated time.
   */
  struct OutputOptions {
    OutputOptions() : output_fraction(false), output_total_updates(false), output_histogram(false),
                      output_minmax(false), output_proc_minmax(false), num_histogram(10), max_levels(INT_MAX),
                      print_warnings(true), align_columns(false), print_names_before_values(true),
                      drop_time(-1.0) {}
    bool output_fraction;
    bool output_total_updates;
    bool output_histogram;
    bool output_minmax;
    bool output_proc_minmax;
    int num_histogram;
    int max_levels;
    bool print_warnings;
    bool align_columns;
    bool print_names_before_values;
    double drop_time;
  };

  /**
   * Dump all the data from all the MPI ranks to an ostream
   * @param [in,out] os - Output stream
   * @param [in] comm - Teuchos comm pointer
   * @param [in] options - Options on what data to output and how to format the data
   */
  void report(std::ostream &os, Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options = OutputOptions());

  /**
   * Dump all the data from all the MPI ranks to the given output stream, in Watchr XML format.
   *  \c reportWatchrXML() is a wrapper for this function that creates \c os as a \c std::ofstream.
   * @param [in,out] os - stream where XML will be written
   * @param [in] datestamp - UTC datestamp in the format YYYY_MM_DD
   * @param [in] timestamp - UTC timestamp in ISO 8601 format: YYYY-MM-DDTHH:MM:SS (HH in 24-hour time, and T is just the character 'T')
   * @param [in] comm - Teuchos comm pointer
   * @return The complete filename, or empty string if no output was produced.
   */
  void reportXML(std::ostream &os, const std::string& datestamp, const std::string& timestamp, Teuchos::RCP<const Teuchos::Comm<int> > comm);

  /**
   * Dump all the data from all the MPI ranks to an XML file. This function
   * calls reportXML and is intended to be used directly by performance tests.
   *
   * In the XML report, the top-level timer name (representing the total) will be replaced with
   * \c name (the name of the test). This becomes the top-level chart name in Watchr.
   *
   * The output filename is controlled by two environment variables: \c $WATCHR_PERF_DIR (required for output)
   * and \c $WATCHR_BUILD_NAME (optional). \c $WATCHR_PERF_DIR is the directory where the output file
   * will be created. If it is not set or is empty, no output will be produced and the function returns an empty
   * string on all ranks.
   *
   * If \c $WATCHR_BUILD_NAME is set, the filename is
   * \c $WATCHR_BUILD_NAME-name_$DATESTAMP.xml. Additionally, the build name is prepended
   * (verbatim) to the the top-level chart name,
   * so that Watchr knows it is a different data series than runs of the same test from other builds.
   *
   * If \c $WATCHR_BUILD_NAME is not set or is empty, the filename is just \c name_$DATESTAMP.xml .
   * DATESTAMP is normally calculated from the current UTC time, in the format YYYY_MM_DD. DATESTAMP
   * can be overridden by setting $WATCHR_BUILD_DATE to a date (also in YYYY_MM_DD format).
   *
   * Optionally, the environment variable \c $TRILINOS_GIT_SHA can be set to the Git revision hash. This
   * is added to the performance report as metadata, and will appear in Watchr when mousing over a data point.
   * Only the first 10 characters are output.
   *
   * In the filename, all spaces in will be replaced by underscores.
   *
   * @param [in] name - Name of performance test
   * @param [in] comm - Teuchos comm pointer
   * @return If on rank 0 and output was produced, the complete output filename. Otherwise the empty string.
   */
  std::string reportWatchrXML(const std::string& name, Teuchos::RCP<const Teuchos::Comm<int> > comm);

  /// If set to true, print timer start/stop to verbose ostream.
  void enableVerbose(const bool enable_verbose);

  /// Enable timestamps in verbose mode for the number of levels specified.
  void enableVerboseTimestamps(const unsigned levels);

  /// Set the ostream for verbose mode(defaults to std::cout).
  void setVerboseOstream(const Teuchos::RCP<std::ostream>& os);

  /// Once called, the start and stop calls are no-ops. Used to stop
  /// timers during asynchronous execution.
  void disableTimers();

  /// Once called, the start and stop calls are reenabled. Used to
  /// restart timers after a call to disableTimers().
  void enableTimers();

  /**
   * Aggregate timer data from all MPI processes to rank 0 of the comm.
   * @param [in] comm - Teuchos comm pointer.
   * @param [in] options - Determines what data the user wants aggregated.
   */
  void aggregateMpiData(Teuchos::RCP<const Teuchos::Comm<int> > comm, OutputOptions options = OutputOptions());

  /**
   * Returns the total time for the timer averaged across all MPI
   * processes (inlcudes only the MPI processes where the timer
   * exists).
   *
   * @param[in] flat_timer_name Timer name that is flattened to include parent names in the string.
   * @return Total timer averaged across all MPI processes.
   */
  double getMpiAverageTime(const std::string& flat_timer_name);

  /**
   * Returns the number of times this timer was started, averaged
   * across all MPI processes (inlcudes only the MPI processes where
   * the timer exists).
   *
   * @param[in] flat_timer_name Timer name that is flattened to include parent names in the string.
   * @return Total number of times the timer was started, averaged across all MPI processes.
   */
  double getMpiAverageCount(const std::string& flat_timer_name);

  /**
   * Returns true if the timer name exists in the stacked timer. This
   * can only be called after aggregateMpiData() is called so that the
   * flattened names are generated.
   *
   * @param[in] flat_timer_name Timer name that is flattened to include parent names in the string.
   * @return true if the timer exists.
   */
  bool isTimer(const std::string& flat_timer_name);

protected:
  /// Current level running
  LevelTimer *top_;
  /// Base timer
  LevelTimer timer_;

  // Global MPI aggregation arrays
  Array<std::string> flat_names_;
  Array<double> min_;
  Array<double> max_;
  Array<int> procmin_;
  Array<int> procmax_;
  Array<double> sum_;
  Array<double> sum_sq_;
  Array<Array<int>> hist_;
  Array<unsigned long> count_;
  Array<unsigned long long> updates_;
  Array<int> active_;
  bool global_mpi_aggregation_called_;

  /// Stores the column widths for output alignment
  struct AlignmentWidths {
    std::string::size_type timer_names_;
    std::string::size_type average_time_;
    std::string::size_type fraction_;
    std::string::size_type count_;
    std::string::size_type total_updates_;
    std::string::size_type min_;
    std::string::size_type max_;
    std::string::size_type procmin_;
    std::string::size_type procmax_;
    std::string::size_type stddev_;
    std::string::size_type histogram_;
    AlignmentWidths() :
      timer_names_(0),
      average_time_(0),
      fraction_(0),
      count_(0),
      total_updates_(0),
      min_(0),
      max_(0),
      procmax_(0),
      stddev_(0),
      histogram_(0){}
  } alignments_;

  /// If set to true, prints to the debug ostream. At construction, default value is set from environment variable.
  bool enable_verbose_;

  /// If set to a value greater than 0, verbose mode will print that many levels of timers with timestamps. A value of zero disables timestamps.
  unsigned verbose_timestamp_levels_;

  /// For debugging, this is the ostream used for printing.
  Teuchos::RCP<std::ostream> verbose_ostream_;

  /// Used to disable timers for asynchronous work.
  bool enable_timers_;

  /**
    * Flatten the timers into a single array
    */
   void flatten();

   /**
    * Merge all the timers together into a single structure
    * @param [in] comm - Communicator to use
    */
   void merge(Teuchos::RCP<const Teuchos::Comm<int> > comm);

   /**
    * Migrate all the timer data to rank=0 if parallel
    */
   void collectRemoteData(Teuchos::RCP<const Teuchos::Comm<int> > comm, const OutputOptions &options );

   /**
    * Returns the index into the mpi aggregated timers for a timer name
    */
   int getFlatNameIndex(const std::string& flat_timer_name);

  /**
   * Compute the column widths to align the output from report() in
   * columns.
   *
   * \returns total time for this level
   */
  double computeColumnWidthsForAligment(std::string prefix, int print_level,
                                        std::vector<bool> &printed, double parent_time,
                                        const OutputOptions &options);

   /**
    * Recursive call to print a level of timer data.
    */
  double printLevel(std::string prefix, int level, std::ostream &os, std::vector<bool> &printed,
                    double parent_time, const OutputOptions &options);

   /**
    * Recursive call to print a level of timer data, in Watchr XML format.
    * If non-empty, rootName will replace the root level's timer name in the output. If empty or level > 0, it has no effect.
    */
  double printLevelXML(std::string prefix, int level, std::ostream &os, std::vector<bool> &printed, double parent_time, const std::string& rootName = "");

};  //StackedTimer


} //namespace Teuchos

#endif /* TEUCHOS_STACKED_TIMER_HPP */
