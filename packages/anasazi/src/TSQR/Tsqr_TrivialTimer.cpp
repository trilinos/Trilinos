#include <Tsqr_TrivialTimer.hpp>
#include <Tsqr_verifyTimerConcept.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  TrivialTimer::TrivialTimer (const std::string& name, bool doStart) :
    name_ (name), isRunning_ (false)
  {
    if (doStart)
      start();
  }

  void 
  TrivialTimer::verifyConcept()
  {
    TSQR::Test::verifyTimerConcept< TrivialTimer >();
  }

  void
  TrivialTimer::start (bool reset) 
  {
    isRunning_ = true;
  }

  double 
  TrivialTimer::stop () 
  { 
    isRunning_ = false;
    return double(0); 
  }

} // namespace TSQR
