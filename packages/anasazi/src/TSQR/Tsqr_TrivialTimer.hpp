#ifndef __TSQR_TrivialTimer_hpp
#define __TSQR_TrivialTimer_hpp

#include <string>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class TrivialTimer
  /// \brief Satisfies TimerType concept trivially
  class TrivialTimer {
  public:
    /// Constructor
    ///
    /// \param name [in] Timer label
    /// \param doStart [in] Whether to start timer on instantiation
    TrivialTimer (const std::string& theName, bool doStart = false);

    /// Start the timer (this implementation doesn't actually do
    /// anything).
    ///
    void start (bool reset = false);
    
    /// Stop the timer and return elapsed time (this implementation
    /// doesn't actually do anything).
    ///
    double stop ();

    /// Whether this timer is running
    ///
    bool isRunning () const { return isRunning_; }

    const std::string& name() const { return name_; }

  private:
    /// Name of this timer
    ///
    std::string name_;
    
    /// Whether this timer is running
    ///
    bool isRunning_;

    /// Verify the TimerType concept
    ///
    static void verifyConcept();
  };

} // namespace TSQR

#endif // __TSQR_TrivialTimer_hpp
