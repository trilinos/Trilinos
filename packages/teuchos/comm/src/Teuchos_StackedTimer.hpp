// @HEADER BEGIN
// @HEADER END

#ifndef TEUCHOS_STACKED_TIMER_HPP
#define TEUCHOS_STACKED_TIMER_HPP

#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include <string>
#include <vector>
#include <cassert>
#include <chrono>

namespace Teuchos {
  
//! Error reporting function for stacked timer.
void error_out(const std::string& msg, const bool fail_all = false);

/**
 * \brief the basic timer used elsewhere, uses MPI_Wtime for time
 *
 * This class hold  a time and number of times this timer is called, and a
 * count of how many "updates" this timer has services.  Example if you have a
 * mesh of 100 elements your count might be 7, but total iterations would be 700
 * this number is useful for dynamic systems where the total number of items services
 * might change in time
 */
class BaseTimer {

public:

  using Clock = std::chrono::steady_clock;
  
  BaseTimer() : accumulation_(0.0), count_started_(0), count_updates_(0), running_(false) {}

  /// Start a currently stopped timer
  void start(){
    if (running_)
      error_out("Base_Timer:start Failed timer already running");
    start_time_ = Clock::now(); // MPI_Wtime();
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
  unsigned long incrementCount(unsigned long count=1) {count_updates_ += count; return count_updates_;}

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
   */  double accumulatedTimePerUpdate()const {
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
  double accumulatedTimePerTimerCall()const {
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
    accumulation_=0.0;
    count_started_ = count_updates_ = 0;
  }

  unsigned long totalUpdates()const {return count_updates_;}

  bool running() const { return running_;}

protected:
  double accumulation_;       // total time
  unsigned long count_started_; // Number of times this timer has been started
  unsigned long count_updates_; // Total count of items updated during this timer
  Clock::time_point start_time_;
  bool running_;


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
class StackedTimer
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
        LevelTimer *parent=NULL,
        bool start=true) :
          BaseTimer(),
          level_(level),
          name_(name),
          parent_(parent),
          mpi_world_(Teuchos::DefaultComm<int>::getComm())
    {
      if ( start )
        BaseTimer::start();
    }

    /// Copy constructor
    LevelTimer(const LevelTimer &src) :
      BaseTimer(src), level_(src.level_), name_(src.name_),parent_(src.parent_), sub_timers_(src.sub_timers_), mpi_world_(src.mpi_world_)
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

#if 0
    // Thought this might be useful, not sure,
    std::string get_full_name() {
      std::string parent_name("");
      if ((parent_ != NULL) && (parent_->level_ > 0))
        parent_name = parent_->get_full_name() + ":";

      std::string my_name(name_);

      std::string full_name = parent_name + my_name;
      return full_name;
    }
#endif

    /**
     * Return the time spent at a given level
     * @param [in] locate_name name of subtimer, if blank return current level time
     * @return time in seconds at provided level
     */
    double accumulatedTime(const std::string &locate_name="") {

      if (locate_name == "")
        return BaseTimer::accumulatedTime();

      std::string first_name,second_name;

      size_t i = locate_name.find_first_of(':');
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

  protected:
    /**
     * \brief split a string into two parts split by a ':' if no ':' first gets the full string
     * @param [in] locate_name input string to split
     * @param [out] first_name Part of string before the first ':'
     * @param [out] second_name part of string after the first ':'
     */
    void splitString(const std::string &locate_name, std::string &first_name, std::string &second_name) {
      size_t i = locate_name.find_first_of(':');
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
      * @param [in/out] os  Where are you dumping the stats, stdout??
      */
     void report(std::ostream &os);


  protected:
     /**
      * Pack a single level and call recursively on sublevels
      * @param [inout] buff Buffer to pack to
      */
     void pack_level(std::vector<char> &buff);
     /**
      * Unpack a single level and call recursively on sublevels
      * @param [in] buff Buffer to unpack from
      * @param [inout] position Location in buffer to unpack from
      */
     void unpack_level(std::vector<char> &buff, int &position);

     constexpr static int BaseLevelTimer=13000;
     constexpr static int ready_to_recv_tag=BaseLevelTimer+1;
     constexpr static int send_size_tag=BaseLevelTimer+2;
     constexpr static int send_buffer_tag=BaseLevelTimer+3;

     Teuchos::RCP<const Teuchos::Comm<int>> mpi_world_;

  }; // LevelTimer


#if 0
  // Disabled for now, needed for parallel stuff
  /**
   * This is a flat version of a level_ timer
   */
  class TimerInstance{
   public:
     TimerInstance(std::string name_,
         unsigned level_,
         float time,
         float parent_,
         int count_started = 0,
         int count_updated = 0) :
           name_(name_),time_(time),parent_(parent_),level_(level_),count_started_(count_started), count_updated_(count_updated)
       { }
     std::string name_;
     float time_;
     float parent_;
     unsigned level_;
     int count_started_;
     int count_updated_;
   }; // TimerInstance
   std::vector< std::vector<TimerInstance> > flat_timers_;
#endif


public:
   /**
    * Construct a stacked timer
    * @param name Top level name of the timer
    */
  explicit StackedTimer(const char *name) :timer_(0,name,NULL,true) {top_ = &timer_;}
  //  ~StackedTimer();
  /**
   * Start the base level timer only, used really in testing only
   */
  void start() { timer_.BaseTimer::start();}
  /**
   * Start a sublevel timer
   * @param [in] name Name of the timer you wish to start
   */
  void start(const std::string name) {top_ = top_->start(name.c_str());}

  /**
   * Stop the current top running timer, should be called for the root timer prior to final output
   */
  void stop(const std::string &name = "RootTimer") {
    if (top_)
      top_ = top_->stop(name);
    else
      timer_.BaseTimer::stop( );
  }

  /**
   * Increment the iteration count for the running timer
   *   @param [in] i amount to increment the count
   */
  void incrementCount(int i) {
    top_->incrementCount(i);
  }

  /**
   * Return the accumulated time the named timer has, if name is blank top level timer
   * @param [in] name Name of the timer to output
   * @return amount of time in seconds
   */
  double accumulatedTime(const std::string &name=""){
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
  void report(std::ostream &os) {
    timer_.report(os);
  }

//  void compile_timers();
//  void compile_timers_collective();

  /// Base timer for all routines to access
  static Teuchos::RCP<StackedTimer> timer;

protected:
  /// Current level running
  LevelTimer *top_;
  /// Base timer
  LevelTimer timer_;

};  //StackedTimer


} //namespace Teuchos

#endif /* TEUCHOS_STACKED_TIMER_HPP */
