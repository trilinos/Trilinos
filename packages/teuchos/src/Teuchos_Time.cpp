// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
#include <cassert>

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
Time::Time(const std::string& name_in, bool start_in) 
  : startTime_(0), totalTime_(0), isRunning_(false), name_(name_in), numCalls_(0)
{
  if(start_in) this->start();
}

void Time::start(bool reset_in)
{
  isRunning_ = true;
  if (reset_in) totalTime_ = 0;
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

#elif ICL || defined(_WIN32)

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
