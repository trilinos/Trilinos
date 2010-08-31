#ifndef __TSQR_StatTimeMonitor_hpp
#define __TSQR_StatTimeMonitor_hpp

#include <Teuchos_Time.hpp>
#include <Tsqr_TimeStats.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  class StatTimeMonitor {
  public:
    StatTimeMonitor (Teuchos::Time& timer, TimeStats& stats);
    ~StatTimeMonitor ();

  private:
    // Copying syntactically forbidden
    StatTimeMonitor (const StatTimeMonitor&);
    StatTimeMonitor& operator= (const StatTimeMonitor&);

    Teuchos::Time& timer_;
    TimeStats stats_;
  };

} // namespace TSQR

#endif // __TSQR_StatTimeMonitor_hpp
