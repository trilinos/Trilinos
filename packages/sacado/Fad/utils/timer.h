// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef _timer_h
#define _timer_h

#include <time.h>

//#include <sys/resource.h>
#include <iostream>
#include <iomanip>

namespace FAD {

class Timer {
public:
  Timer() : state_(uninitialized), t1_(0), t2_(0) {;}  
  Timer(double floating_ops1 ) : floating_ops(floating_ops1) {}
  ~Timer() {;}

  void start()
    { 
      state_ = running;
      t1_ = systemTime();
    }

  void stop()
    {
      t2_ = systemTime();
      state_ = stopped;
    }

  double elapsedSeconds()
    {
      return t2_ - t1_;
    }

  double sec()
    {
      return t2_ - t1_;
    }

  void report( char const* what ) const
    {
      double t = t2_ - t1_;
      std::cout << "Time for "
		<< what
		<< " = "
		<< setw(8)
		<< t
		<< " seconds [ "
		<< setw(6)
		<< floating_ops/t/1000000 
		<< " Mflops]"
		<< std::endl;
    }
 private:
  Timer(Timer&) { }
  void operator=(Timer&) { }

  double systemTime()
    {
    static const double secs_per_tick = ((double)1.0) / CLOCKS_PER_SEC;
	const clock_t ticks = clock();
	const double sec = ( (double) ticks ) * secs_per_tick;
	//std::cout << "ticks = " << ticks << ", sec = " << sec << std::endl;
    return sec;
/*
      getrusage(RUSAGE_SELF, &resourceUsage_);
      double seconds = resourceUsage_.ru_utime.tv_sec 
	+ resourceUsage_.ru_stime.tv_sec;
      double micros  = resourceUsage_.ru_utime.tv_usec 
	+ resourceUsage_.ru_stime.tv_usec;
      return seconds + micros/1.0e6;
*/
    }

  enum { uninitialized, running, stopped } state_;

  double floating_ops;

//  struct rusage resourceUsage_;

  double t1_, t2_;
};

} // namespace FAD

#endif
