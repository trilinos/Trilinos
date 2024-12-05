/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
//@HEADER
*/

#ifndef TIMER_DH_H
#define TIMER_DH_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include "euclid_common.h"

/*--------------------------------------------------------------*/
/* Stuff in this block isn't really needed for multi-processor
 * runs, since recording CPU time probably isn't useful.
 * if EUCLID_TIMING is defined in $PCPACK_DIR/bmake_XXX/common,
 * the times() function is used;
 * then MPI_Wtime() is used in preference to times().
 *
 * You may need to fiddle with some of these includes, depending
 * on your system.  Make sure and check the logFile to ensure
 * that CLK_TCK was properly defined.  See Timer_dhCreate()
 * for additional details. 
 *
 * if "JUNK_TIMING" is defined during compilation, timing functions
 * either do nothing, or return -1.0; this is primarily for debugging.
 */

#ifdef EUCLID_TIMING
#include <sys/times.h>
#include <sys/types.h>
#if HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#elif !defined(JUNK_TIMING)
/* #include <sys/types.h> 
#include <sys/sysconfig.h>
*/
#ifdef WIN32
#  include <time.h>
#endif
#if HAVE_UNISTD_H
#  include <unistd.h>		/* needed for sysconf(_SC_CLK_TCK) */
#endif /* HAVE_UNISTD_H */
#endif


/* 
   ??? may be needed for some compilers/platforms?
#include <limits.h>
#include <time.h>
#include <sys/resource.h>
*/

/*--------------------------------------------------------------*/

#ifdef __cplusplus
extern "C"
{
#endif

  struct _timer_dh
  {
    bool isRunning;
    long int sc_clk_tck;
    double begin_wall;
    double end_wall;

#ifdef EUCLID_TIMING
    struct tms begin_cpu;
    struct tms end_cpu;
#endif

  };

  extern void Timer_dhCreate (Timer_dh * t);
  extern void Timer_dhDestroy (Timer_dh t);
  extern void Timer_dhStart (Timer_dh t);
  extern void Timer_dhStop (Timer_dh t);
  extern double Timer_dhReadCPU (Timer_dh t);
  extern double Timer_dhReadWall (Timer_dh t);
  extern double Timer_dhReadUsage (Timer_dh t);

/* notes:
    (1)  unless compiled with EUCLID_TIMING defined, readCPU 
         and readUseage return -1.0.
    (2)  whenever start() is called, the timer is reset; you
         don't need to call stop() first.
    (3)  if stop() HAS been called, the readXX functions return
         timings between start() and stop(); , if start()
         was called but not stop(), they sample the time then
         return; if start() was never called, they return junk.
*/

#ifdef __cplusplus
}
#endif

#endif
