/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */

#ifndef __TIMER_H
#define __TIMER_H

#include "zoltan_timer.h"
#include <time.h> /* ANSI C; defines clock_t and clock() */

#ifndef CLOCKS_PER_SEC /* Should have been defined in time.h */
#define CLOCKS_PER_SEC 1000000 /* To prevent compile errors, not always the correct value. */
#endif

/*
 * POSIX compliant systems should use times() for user timing. 
 * This is the default in Zoltan. Make Zoltan with -DNO_TIMES if
 * your system does not have sys/times.h and times().
 * Note: BSD-like systems may use getrusage() instead for user timing,
 * but that has not been implemented here. 
 */

#if defined(__PUMAGON__) || defined(__LIBCATAMOUNT__) || defined(_MSC_VER)
/* Tflops with Cougar & Red Storm w/Catamount does not have sysconf() or times() */
/* Microsoft Visual Studio does not have times either */
#define NO_TIMES
#endif /* __PUMAGON__ */

#ifndef NO_TIMES
/* #include <sys/types.h> -- Included by sys/times.h on most systems. */
#include <sys/times.h>
#include <unistd.h> /* Needed for sysconf() and _SC_CLK_TCK */
#endif

#endif
