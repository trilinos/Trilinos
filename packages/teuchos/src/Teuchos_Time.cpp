// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Time.hpp"
#include <sys/time.h>
#include <sys/resource.h>

using namespace Teuchos;

//=============================================================================
Time::Time(const string& name, bool start) 
  : startTime_(0), totalTime_(0), isRunning_(false), name_(name)
{
	if(start) this->start();
}

void Time::start()
{
  isRunning_ = true;
  startTime_ = wallTime();
}

double Time::stop()
{
  totalTime_ += wallTime() - startTime_;
  isRunning_ = false;
  startTime_ = 0;
  return totalTime_;
}



//=============================================================================
double Time::wallTime() 
{
  /* KL: warning: this code is probably not portable! */
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

