// 27-May-2002 Major Overhaul. Changed method names to fit namingConvention (already done).
// 28-May-2002 Heroux fixes things.

#ifndef _TPETRA_TIME_H_
#define _TPETRA_TIME_H_

#include "Tpetra_Object.h"
#include "Tpetra_Comm.h"
#include <sys/time.h>
#include <sys/resource.h>

/* To time a section of code, place it in between calls to resetStartTime and elapsedTime. 
   Note that resetStartTime is called by the constructor, so it is also possible to time a section of code 
   by placing it in between the creation of the Tpetra::Time object, and a call to elapsedTime. 
   A call to resetStartTime must be used for all subsequent timings though.

   A Tpetra::Comm object is also required to use Tpetra::Time, although I don't know why.
*/

namespace Tpetra {

class Time : public Tpetra::Object {

public:
  // Default constructor
  Time(const Tpetra::Comm &comm);
  
  // Copy Constructor
  Time(const Tpetra::Time &time);
  
  // Destructor
  virtual ~Time();
  
  // returns current wall-clock time in seconds
  double wallTime() const;
  
  // resets the timer to the current walltime
  void resetStartTime();
  
  // returns the elapsed time in between when the timer was set, and the current walltime.
  double elapsedTime() const;


private:

  double startTime_;
  const Tpetra::Comm *comm_;

};

} // namespace Tpetra

#include "Tpetra_Time.cpp"

#endif // TPETRA_TIME_H_
