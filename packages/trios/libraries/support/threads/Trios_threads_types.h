/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/**
 * This API is modeled after pthreads.  pthreads is a POSIX API.  TRIOS runs in non-POSIX
 * environments such as Cray MTA.  These systems don't have pthreads, but may have other
 * threading libs.  This is lib creates a single API for TRIOS internal usage.
 */

#ifndef _TRIOS_THREADS_TYPES_H_
#define _TRIOS_THREADS_TYPES_H_

#include "Trios_config.h"


#include <stdint.h>

#ifdef HAVE_TRIOS_PTHREAD

#include <pthread.h>

typedef pthread_t nthread_id_t;

typedef struct {
    pthread_t thread;
} nthread_t;

typedef struct {
    pthread_attr_t attr;
} nthread_attr_t;

typedef struct {
    pthread_cond_t cond;
} nthread_cond_t;

typedef struct {
    pthread_condattr_t condattr;
} nthread_condattr_t;

typedef struct {
    pthread_mutex_t mutex;
} nthread_mutex_t;

typedef struct {
    nthread_mutex_t mutex;
    int64_t         value;
} nthread_counter_t;

#elif HAVE_TRIOS_MTA

#if __MTA__
#include <machine/runtime.h>
#define MTA_SYNC sync uint64_t
#else
#define MTA_SYNC uint64_t
#endif
#include <time.h>

typedef struct {
    uint64_t thread;
} nthread_t;

typedef struct {
    uint64_t attr;
} nthread_attr_t;

typedef struct {
    MTA_SYNC mutex;
    uint64_t      is_recursive;

//    MTA_SYNC locker_mutex;
    uint64_t      is_locked;
    nthread_id_t locker;
    int64_t      lock_count;
} nthread_mutex_t;

typedef struct {
    MTA_SYNC mutex;
    int64_t      value;
} nthread_counter_t;

typedef struct {
    uint64_t         cond;
    int64_t         waiter_count;
    nthread_mutex_t *bound_mutex;
} nthread_cond_t;

typedef struct {
    uint64_t condattr;
} nthread_condattr_t;

#else

#include <time.h>

typedef struct {
    uint32_t thread;
} nthread_t;

typedef struct {
    uint32_t attr;
} nthread_attr_t;

typedef struct {
    uint32_t cond;
} nthread_cond_t;

typedef struct {
    uint32_t condattr;
} nthread_condattr_t;

typedef struct {
    uint32_t mutex;
} nthread_mutex_t;

typedef struct {
    nthread_mutex_t mutex;
    int64_t         value;
} nthread_counter_t;

#endif



enum {
    NTHREAD_MUTEX_NORMAL,
    NTHREAD_MUTEX_RECURSIVE,
    NTHREAD_MUTEX_ERRORCHECK,
    NTHREAD_MUTEX_DEFAULT
};

#ifdef HAVE_TRIOS_PTHREAD

#define NTHREAD_MUTEX_INITIALIZER \
    { PTHREAD_MUTEX_INITIALIZER }
#ifdef __USE_GNU
#define NTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP \
    { PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP }
#define NTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP \
    { PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP }
#define NTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP \
    { PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP }
#endif


#define NTHREAD_COND_INITIALIZER \
    { PTHREAD_COND_INITIALIZER }

#define NTHREAD_COUNTER_INITIALIZER \
    { NTHREAD_MUTEX_INITIALIZER, 0 }

#elif HAVE_TRIOS_MTA

#define NTHREAD_MUTEX_INITIALIZER \
    { 0, 0, 0, 0, 0 }
#define NTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP \
    { 0, 1, 0, 0, 0 }
#define NTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP \
    { 0, 0, 0, 0, 0 }
#define NTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP \
    { 0, 0, 0, 0, 0 }


#define NTHREAD_COND_INITIALIZER \
    { 0, 0, NULL }

#define NTHREAD_COUNTER_INITIALIZER \
    { 0, 0 }

#else

#define NTHREAD_MUTEX_INITIALIZER \
    { 0 }
#ifdef __USE_GNU
#define NTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP \
    { 0 }
#define NTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP \
    { 0 }
#define NTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP \
    { 0 }
#endif


#define NTHREAD_COND_INITIALIZER \
    { 0 }

#define NTHREAD_COUNTER_INITIALIZER \
    { 0, 0 }

#endif

#endif /* _TRIOS_THREADS_TYPES_H_ */
