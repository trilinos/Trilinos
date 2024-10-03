// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
LARGE_INTEGER start_count, count_freq;  // counts per sec.

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

}       // end namespace

#endif // defined(__INTEL_COMPILER) && defined(_WIN32)

#ifdef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
#include <valgrind.h>
#include <algorithm>
#include <unistd.h>
#endif

#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
namespace Kokkos {
namespace Tools {
extern void pushRegion (const std::string&);
extern void popRegion ();
} // namespace Profiling
} // namespace Kokkos
#endif

#ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
#include "Kokkos_Core.hpp"
#endif

namespace Teuchos {

#ifdef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
  void removeIllegalChars(std::string& s){
    std::string illegalChars = "\\/:?\"<>|";
    for (auto it = s.begin() ; it < s.end() ; ++it){
      bool found = illegalChars.find(*it) != std::string::npos;
      if(found)
        *it = ' ';
    }
  }
#endif

//=============================================================================
Time::Time(const std::string& name_in, bool start_in)
  : startTime_(0), totalTime_(0), isRunning_(false), enabled_ (true), name_(name_in), numCalls_(0)
{
  if(start_in) this->start();
#ifdef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
  numCallsMassifSnapshots_ = 0;
#endif
}

void Time::start(bool reset_in)
{
  if (enabled_) {
    isRunning_ = true;
    if (reset_in) totalTime_ = 0;
#ifdef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
    if (numCallsMassifSnapshots_ < 100) {
      std::string filename = "massif.out." + std::to_string(::getpid()) + "." + name_ + "." + std::to_string(numCallsMassifSnapshots_) + ".start.out";
      removeIllegalChars(filename);
      std::replace(filename.begin(), filename.end(), ' ', '_');
      std::string cmd = "snapshot " + filename;
      VALGRIND_MONITOR_COMMAND(cmd.data());
    }
#endif
    startTime_ = wallTime();
#ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
    Kokkos::fence();
#endif
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
    ::Kokkos::Tools::pushRegion (name_);
#endif
  }
}

double Time::stop()
{
  if (enabled_) {
    if (isRunning_) {
      totalTime_ += ( wallTime() - startTime_ );
      isRunning_ = false;
      startTime_ = 0;
#ifdef HAVE_TEUCHOS_TIME_MASSIF_SNAPSHOTS
      if (numCallsMassifSnapshots_ < 100) {
        std::string filename = "massif.out." + std::to_string(::getpid()) + "." + name_ + "." + std::to_string(numCallsMassifSnapshots_) + ".stop.out";
        removeIllegalChars(filename);
        std::replace(filename.begin(), filename.end(), ' ', '_');
        std::string cmd = "snapshot " + filename;
        VALGRIND_MONITOR_COMMAND(cmd.data());
        numCallsMassifSnapshots_++;
      }
#endif
#ifdef HAVE_TEUCHOS_TIMER_KOKKOS_FENCE
      Kokkos::fence();
#endif
#if defined(HAVE_TEUCHOS_KOKKOS_PROFILING) && defined(HAVE_TEUCHOSCORE_KOKKOS)
      ::Kokkos::Tools::popRegion ();
#endif
    }
  }
  return totalTime_;
}

double Time::totalElapsedTime(bool readCurrentTime) const
{
  if(readCurrentTime)
    return wallTime() - startTime_ + totalTime_;
  return totalTime_;
}

void Time::reset () {
  totalTime_ = 0;
  numCalls_ = 0;
}

void Time::disable () {
  enabled_ = false;
}

void Time::enable () {
  enabled_ = true;
}

void Time::incrementNumCalls() {
  if (enabled_) {
    ++numCalls_;
  }
}

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
