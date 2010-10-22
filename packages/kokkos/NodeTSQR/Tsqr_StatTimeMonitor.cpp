#include <Tsqr_StatTimeMonitor.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  StatTimeMonitor::StatTimeMonitor (Teuchos::Time& timer, TimeStats& stats)
    : timer_ (timer), stats_ (stats)
  {
    if (timer_.isRunning())
      timer_.stop(); // Just for sanity
    if (timer_.numCalls() == 0)
      stats_.init();
    timer_.start (true);
  }

  StatTimeMonitor::~StatTimeMonitor () 
  {
    const double curTime = timer_.stop();
    stats_.update (curTime);
  }

#if 0
  /// \brief Return total elapsed time of a particular timer
  ///
  /// Return the total elapsed time of a particular timer. 
  /// Ensures that the timer is not running (which would break
  /// totalElapsedTime()).
  static double
  fetchTime (const Teuchos::RCP< Teuchos::Time >& timer) 
  {
    if (timer->isRunning())
      timer->stop();
    return timer->totalElapsedTime();
  }
#endif // 0


} // namespace TSQR
