/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
/*-------------------------------------------------------------------------*/
/**  @file time.h
 *
 *   @brief API for calculating statistics.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 406 $.
 *   $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $.
 *
 */

#ifndef _TRIOS_TIMER_H_
#define _TRIOS_TIMER_H_

#include "Trios_logger.h"

#include <stdio.h>

#ifdef __cplusplus

#include <string>
#include <vector>
#include <ostream>


namespace Trios {

double GetTime();
long GetTimeMS();
long GetTimeNS();
long GetTimeUS();

int WriteTimings(const std::string &fname,
        const std::string &header,
        const std::vector<std::string> &timing_desc,
        const std::vector<double> &timings);

int WriteTimings(std::ostream &out,
        const std::string &header,
        const std::vector<std::string> &timing_desc,
        const std::vector<double> &timings,
        const bool write_header=true);

}
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)

extern log_level timer_debug_level;

    extern double trios_get_time();
    extern long trios_get_time_ns();
    extern long trios_get_time_ms();
    extern long trios_get_time_us();
    extern int trios_timer_test();
    extern const char *trios_timer_getimpl();


#if defined(TRIOS_USE_TIMERS)

/* always use this macro to declare timer variables */
#define trios_declare_timer(timer_var) double timer_var;

#define trios_start_timer(timer) { timer = trios_get_time(); }
#define trios_stop_timer(name, timer)  { timer = trios_get_time() - timer; log_debug(timer_debug_level, "%s Time = %10.8f", name, timer); }

#define trios_start_delay_timer(timer) { timer = trios_get_time(); }
#define trios_stop_delay_timer(timer)  { timer = trios_get_time() - timer; }
#define trios_log_delay_timer(name, timer)  { log_debug(timer_debug_level, "%s Time = %10.8f", name, timer); }

#else

#define trios_declare_timer(t)

#define trios_start_timer(timer)  {}
#define trios_stop_timer(name, timer)   {}

#define trios_start_delay_timer(timer)  {}
#define trios_stop_delay_timer(timer)   {}
#define trios_log_delay_timer(name, timer)   {}

#endif

#endif


#ifdef __cplusplus
}
#endif

#endif
