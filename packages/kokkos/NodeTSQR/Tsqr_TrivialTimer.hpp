#ifndef __TSQR_TrivialTimer_hpp
#define __TSQR_TrivialTimer_hpp

#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TrivialTimer
  /// \brief Satisfies TimerType concept trivially.
  ///
  /// This is a "prototype" for the TimerType concept; it satisfies
  /// the concept trivially.
  class TrivialTimer {
  public:
    /// Constructor.
    ///
    /// \param name [in] Timer label
    /// \param doStart [in] Whether to start timer on instantiation
    TrivialTimer (const std::string& theName, bool doStart = false);

    /// \brief Start the timer.
    ///
    /// This is a trivial timer, so this implementation does nothing.
    /// However, it satisfies our TimerType concept.
    void start (bool reset = false);
    
    /// \brief Stop the timer and return elapsed time.
    ///
    /// This is a trivial timer, so his implementation does nothing.
    /// This method should return 0.
    double stop ();

    //! Whether this timer is running
    bool isRunning () const { return isRunning_; }

    //! Name of this timer object
    const std::string& name() const { return name_; }

  private:
    //! Name of this timer
    std::string name_;
    
    //! Whether this timer is running
    bool isRunning_;

    //! Verify the TimerType concept
    static void verifyConcept();
  };

} // namespace TSQR

#endif // __TSQR_TrivialTimer_hpp
