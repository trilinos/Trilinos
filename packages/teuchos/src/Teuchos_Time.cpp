/*Paul
27-May-2002 Major Overhaul. Changed method names to fit naming Convention (already done).
28-May-2002 Heroux fixes things.
06-August-2002 Changed to images (nothing changed). Updated a few naming comnventions, and to use a templated Comm.
*/

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Time.hpp"

namespace Teuchos {
//=============================================================================
Time::Time(const Comm<double, int>& Comm) : Object("Teuchos::Time"), Comm_(&Comm) {
  startTime_ = wallTime();
}

//=============================================================================
Time::Time(const Time& Time) : Object(Time.label()), Comm_(Time.Comm_) {
  startTime_ = Time.startTime_; // We do want to do this in a copy constructor, after all, right?
}

//=============================================================================
double Time::wallTime() const {
  struct timeval tp;
  static long start = 0, startu;
  if (!start)
  {
    gettimeofday(&tp, NULL);
    start = tp.tv_sec;
    startu = tp.tv_usec;
    return(0.0);
  }
  gettimeofday(&tp, NULL);
  return( ((double) (tp.tv_sec - start)) + (tp.tv_usec-startu)/1000000.0 );
}

//=============================================================================
void Time::resetStartTime() {
  startTime_ = wallTime();
}

//=============================================================================
double Time::elapsedTime() const {
  return(wallTime() - startTime_);
}

}  // namespace Teuchos
