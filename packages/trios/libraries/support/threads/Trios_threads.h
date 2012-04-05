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

#ifndef _TRIOS_THREADS_H_
#define _TRIOS_THREADS_H_

#include "Trios_config.h"

#include "Trios_logger.h"
#include "Trios_threads_types.h"


#ifdef __cplusplus
extern "C" {
#endif

extern log_level thread_debug_level;



int nthread_create(
        nthread_t *thread,
        const nthread_attr_t *attr,
        void *(*start_routine)(void*),
        void *arg);

/**
 * Reclaim the space allocated for the thread when the thread exits.
 */
int nthread_detach(nthread_t thread);

void nthread_exit(void *retval);

int nthread_cond_init(
        nthread_cond_t *cond);
int nthread_cond_broadcast(nthread_cond_t *cond);
int nthread_cond_signal(nthread_cond_t *cond);
int nthread_cond_timedwait(
        nthread_cond_t *cond,
        nthread_mutex_t *mutex,
        const struct timespec *abstime);
int nthread_cond_wait(
        nthread_cond_t *cond,
        nthread_mutex_t *mutex);

int nthread_mutex_init(nthread_mutex_t *mutex, int mutex_type);
int nthread_lock(nthread_mutex_t *mutex);
int nthread_unlock(nthread_mutex_t *mutex);

nthread_id_t nthread_self(void);

void nthread_yield(void);
int nthread_join(
        nthread_t thread,
        void **value_ptr);

int64_t nthread_counter_increment(nthread_counter_t *c);
int64_t nthread_counter_decrement(nthread_counter_t *c);

#ifdef __cplusplus
}
#endif

#endif /* _TRIOS_THREADS_H_ */
