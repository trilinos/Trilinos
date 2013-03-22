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
/**  @file timer.c
 *
 *   @brief A portable timer for Trios.
 *
 *   @author Ron Oldfield (raoldfi@sandia.gov)
 *   @version $Revision: 406 $
 *   @date $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $
 *
 */

#include "Trios_config.h"

#include <math.h>      /* for fabsl() */
#include <stdlib.h>    /* for llabs() */
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include "Trios_logger.h"

log_level timer_debug_level = LOG_UNDEFINED;

/*-------------------- Global Variables ------------------ */

/* Even though configure scripts find clock_gettime, it is
 * not supported by Catamount and will cause problems if used.
 */
#ifdef __LIBCATAMOUNT__
#undef HAVE_TRIOS_CLOCK_GETTIME
#endif

#undef USING_PAPI
#undef USING_MACH_ABSOLUTE_TIME
#undef USING_CLOCK_GETTIME
#undef USING_GETTIMEOFDAY


#if defined(HAVE_TRIOS_PAPI)
#define USING_PAPI
#include <papi.h>
long_long init_usec, init_cycles;


#elif defined(HAVE_TRIOS_MACH_ABSOLUTE_TIME)
#define USING_MACH_ABSOLUTE_TIME
#include <mach/mach_time.h>

#elif defined(HAVE_TRIOS_GETTIMEOFDAY)
#define USING_GETTIMEOFDAY
#include <sys/time.h>
static struct timeval tv_init;


#elif defined(HAVE_TRIOS_CLOCK_GETTIME)
#define USING_CLOCK_GETTIME
#include <time.h>   /* values for various timers */

/* Options for clockid:
 *  - CLOCK_REALTIME,
 *  - CLOCK_MONOTONIC,
 *  - CLOCK_PROCESS_CPUTIME_ID,
 *  - CLOCK_THREAD_CPUTIME_ID
 */
static const clockid_t clockid = CLOCK_REALTIME;
static struct timespec ts_init;
static struct timespec ts_res;

#else
#error No supported timers.
#endif


/*-------------------- Local Functions ------------------ */

static void _timer_premain (void) __attribute__ ((constructor));
static void init_timer();

static void _timer_premain(void)
{
    init_timer();
}

static void init_timer() {
    static int initialized = 0;
    if (!initialized) {

#ifdef USING_PAPI
        PAPI_library_init(PAPI_VER_CURRENT);

        init_usec = PAPI_get_real_usec();
        init_cycles = PAPI_get_real_cyc();
#endif

#ifdef USING_CLOCK_GETTIME
        clock_getres(clockid, &ts_res);
        clock_gettime(clockid, &ts_init);
#endif

#ifdef USING_GETTIMEOFDAY
        gettimeofday(&tv_init, NULL);
#endif

        initialized = 1;
    }
}


#if defined(USING_CLOCK_GETTIME)
static uint64_t ts_diff_nsec(struct timespec start, struct timespec end)
{
    struct timespec diff;
    uint64_t result = 0;

    diff.tv_sec = end.tv_sec - start.tv_sec;
    diff.tv_nsec  = end.tv_nsec - start.tv_nsec;

    result = diff.tv_sec*1e9 + diff.tv_nsec;
    return result;
}
#endif

#ifdef USING_GETTIMEOFDAY

static uint64_t tv_diff_us(struct timeval start, struct timeval end)
{
    struct timeval diff;
    uint64_t result = 0;

    diff.tv_sec = end.tv_sec - start.tv_sec;
    diff.tv_usec  = end.tv_usec - start.tv_usec;

    result = diff.tv_sec*1000000 + diff.tv_usec;
    return result;
}
#endif


/*-------------------- Timer API ------------------ */

/**
 * @brief Return the time in nanoseconds from first call to trios_get_time_*.
 */
uint64_t trios_get_time_ns()
{
    uint64_t result = 0;

#ifdef USING_PAPI

    long_long cur_usec, diff_usec;

    init_timer();
    cur_usec = PAPI_get_real_usec();
    diff_usec = cur_usec - init_usec;

    /* to convert to ns, multiply by 1000 */
    result = diff_usec * 1000;
#endif

#ifdef USING_CLOCK_GETTIME
    struct timespec tp;

    init_timer();
    clock_gettime(clockid, &tp);
    result = ts_diff_nsec(ts_init, tp);
#endif

#ifdef USING_MACH_ABSOLUTE_TIME
    static mach_timebase_info_data_t info = {0,0};

    init_timer();
    if (info.denom == 0)
        mach_timebase_info(&info);

    result = mach_absolute_time() * (info.numer / info.denom);
#endif

#ifdef USING_GETTIMEOFDAY
    struct timeval tp;

    init_timer();
    gettimeofday( &tp, NULL );

    result = tv_diff_us(tv_init, tp);
    result = result*(uint64_t)1000;
#endif


    return result;
}

/**
 * @brief Return the time in seconds.
 *
 */
double trios_get_time()
{
    double result = 0.0;
    uint64_t ns = trios_get_time_ns();
    result = (double)ns*1.0e-9;
    return result;
}

long trios_get_time_us()
{
    uint64_t result;
    result = trios_get_time_ns() * 1e-3;
    return result;
}

long trios_get_time_ms()
{
    uint64_t result;

    result = trios_get_time_ns() * 1e-6;
    return result;
}

/**
 * @brief Test the timer.
 *
 * An easy way to test the timer is to time the sleep function
 * for a fixed amount of time.  If the timer results a time that
 * matches the sleep time within a reasonable error, the timer works.
 *
 * @return 0 on success, 1 on error.
 */
int trios_timer_test()
{
    int rc = 0;
/*    static const uint64_t sleep_us = 1000000; */  /* 1 second */
    static const unsigned int sleep_s = 1;  /* 1 second */

    /* TODO: what should be the acceptable error (1/2 sec?) */
    static const uint64_t error_ms = 500;

    uint64_t error_ns = error_ms * 1000000;
    uint64_t error_us = error_ms * 1000;
    double error_sec = (double)error_ms/1000.0;

    uint64_t start_ns=0, t_ns=0;
    uint64_t start_us, t_us;
    uint64_t start_ms, t_ms;
    double start_sec, t_sec;

    log_debug(timer_debug_level, "this is a timer test.  expect a 1 second delay.  starting...");

    /* run to initialize timer */
    start_ns = trios_get_time_ns();
    if (sleep(sleep_s)) {
        log_error(timer_debug_level, "Failed calling sleep()");
        return 1;
    }

    /* run to get measurements */
    start_ns = trios_get_time_ns();
    start_us = trios_get_time_us();
    start_ms = trios_get_time_ms();
    start_sec = trios_get_time();

    /* sleep one seconds */
    sleep(sleep_s);

    t_ns = trios_get_time_ns() - start_ns;
    t_us = trios_get_time_us() - start_us;
    t_ms = trios_get_time_ms() - start_ms;
    t_sec = trios_get_time() - start_sec;

    printf("slept for %u seconds:\n", sleep_s);
    printf("\tns = %lu\n", t_ns);
    printf("\tus = %lu\n", t_us);
    printf("\tms = %lu\n", t_ms);
    printf("\tsec = %f\n", t_sec);

    /* Make sure our values have a reasonable error */
    if (labs(t_ns - ((double)sleep_s * 1e9)) > error_ns) {
        uint64_t actual_err = labs(t_ns - ((double)sleep_s * 1e9));
        log_error(timer_debug_level,
                "trios_timer failed ns timer test: err = %llu ns"
                " > acceptable err = %llu",
                actual_err, error_ns);
        rc = 1;
    }

    if (labs(t_us - ((double)sleep_s * 1e6)) > error_us) {
        uint64_t actual_err = labs(t_us - ((double)sleep_s * 1e6));
        log_error(timer_debug_level,
                "trios_timer failed usec timer test: err = %llu usec"
                " > acceptable err = %llu",
                actual_err, error_us);
        rc |= 1;
    }

    if (labs(t_ms - ((double)sleep_s * 1e3)) > error_ms) {
        uint64_t actual_err = labs(t_ms - ((double)sleep_s * 1e3));
        log_error(timer_debug_level,
                "trios_timer failed ns timer test: err = %llu ms"
                " > acceptable err = %llu",
                actual_err, error_ms);
        rc |= 1;
    }

    if (fabs(t_sec - (double)sleep_s) > error_sec) {
        double actual_err = fabs(t_sec - (double)sleep_s);
        log_error(timer_debug_level,
                "trios_timer failed sec timer test: err = %g sec "
                " > acceptable err = %g",
                actual_err, error_sec);
        rc |= 1;
    }

    log_debug(timer_debug_level, "timer test complete");

    return rc;
}

const char *trios_timer_getimpl()
{
#ifdef USING_PAPI
    return "PAPI_get_real_usec()";
#elif defined(USING_CLOCK_GETTIME)
    return "clock_gettime()";
#elif defined (USING_MACH_ABSOLUTE_TIME)
    return "mach_absolute_time()";
#elif defined(USING_GETTIMEOFDAY)
    return "gettimeofday()";
#else
    return "no timer implementation";
#endif
}
