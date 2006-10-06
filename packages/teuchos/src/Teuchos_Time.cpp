// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// Kris
// 07.08.03 -- Move into Teuchos package/namespace

#include "Teuchos_Time.hpp"

#if defined(__INTEL_COMPILER) && defined(_WIN32)

#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <assert.h>

namespace {

bool seconds_initialized = false;
LARGE_INTEGER start_count, count_freq;	// counts per sec.

inline void seconds_initialize() {
  if( seconds_initialized ) return;
  std::cout << "\nCalling Win32 version of Teuchos::seconds_initialize()!\n";
  // Figure out how often the performance counter increments
  ::QueryPerformanceFrequency( &count_freq );
  // Set this thread's priority as high as reasonably possible to prevent
  // timeslice interruptions
  ::SetThreadPriority( ::GetCurrentThread(), THREAD_PRIORITY_TIME_CRITICAL );
  // Get the first count.
  assert( QueryPerformanceCounter( &start_count ) );
  seconds_initialized = true;
}

}	// end namespace

#endif // defined(__INTEL_COMPILER) && defined(_WIN32)

namespace Teuchos {

//=============================================================================
Time::Time(const string& name, bool start) 
  : startTime_(0), totalTime_(0), isRunning_(false), name_(name), numCalls_(0)
{
  if(start) this->start();
}

void Time::start(bool reset)
{
  isRunning_ = true;
  if (reset) totalTime_ = 0;
  startTime_ = wallTime();
}

double Time::stop()
{
  if (isRunning_) {
    totalTime_ += ( wallTime() - startTime_ );
    isRunning_ = false;
    startTime_ = 0;
  }
  return totalTime_;
}

double Time::totalElapsedTime(bool readCurrentTime) const
{
  if(readCurrentTime)
    return wallTime() - startTime_ + totalTime_;
  return totalTime_;
}

//=============================================================================
double Time::wallTime() 
{
  /* KL: warning: this code is probably not portable! */
  /* HT: have added some preprocessing to address problem compilers */
	/* RAB: I modifed so that timer will work if MPI support is compiled in but not initialized */

#ifdef HAVE_MPI

	int mpiInitialized;
	MPI_Initialized(&mpiInitialized);

	if( mpiInitialized ) {

		return(MPI_Wtime());

	}
	else {

		clock_t start;

		start = clock();
		return( (double)( start ) / CLOCKS_PER_SEC );

	}

#elif defined(__INTEL_COMPILER) && defined(_WIN32)

  seconds_initialize();
  LARGE_INTEGER count;
  QueryPerformanceCounter( &count );
  // "QuadPart" is a 64 bit integer (__int64).  VC++ supports them!
  const double
    sec = (double)( count.QuadPart - start_count.QuadPart ) / count_freq.QuadPart;
  //std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
  return sec;

#elif ICL

  clock_t start;

  start = clock();
  return( (double)( start ) / CLOCKS_PER_SEC );

#else

#  ifndef MINGW
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
#  else // MINGW
  return( (double) clock() / CLOCKS_PER_SEC );
#  endif // MINGW

#endif

}

} // namespace Teuchos
