// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Time.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#elif ICL
#include <time.h>
#else
#include <sys/time.h>
#ifndef MINGW
#include <sys/resource.h>
#endif
#endif

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
  /* HT: have added some preprocessing to address problem compilers */
#ifdef HAVE_MPI

  return(MPI_Wtime());

#elif ICL

  clock_t start;

  start = clock();
  return( (double)( start ) / CLOCKS_PER_SEC );

#else

#ifndef MINGW
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
#else
  return( (double) clock() / CLOCKS_PER_SEC );
#endif

#endif

}

