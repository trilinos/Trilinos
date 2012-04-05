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

#include "Trios_config.h"

#include <assert.h>
#include <unistd.h>

#if defined(HAVE_TRIOS_SCHED_YIELD)
#include <sched.h>
#endif


#include "Trios_logger.h"
#include "Trios_threads.h"

#if USING_MTGL
#include "mtgl/mtgl_config.h"
#include "mtgl/util.hpp"
#endif

#if HAVE_TRIOS_PTHREAD
#include <pthread.h>
#endif

#define _DEBUG_LOCKS_
#undef _DEBUG_LOCKS_
#define _DEBUG_COND_
#undef _DEBUG_COND_

/* set to LOG_UNDEFINED -- log commands will use default level */
log_level thread_debug_level = LOG_UNDEFINED;


int nthread_create(
        nthread_t *thread,
        const nthread_attr_t *attr,
        void *(*start_routine)(void*),
        void *arg)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_create(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_create(&thread->thread, &attr->attr, start_routine, arg);
#endif
    return rc;
}

int nthread_detach(
        nthread_t thread)
{
    int rc = 0;
    log_debug(thread_debug_level, "nthread_detach(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_detach(thread.thread);
#endif
    return rc;
}

void nthread_exit(void *retval)
{
    log_debug(thread_debug_level, "nthread_exit(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    pthread_exit(retval);
#endif
}

int nthread_cond_init(
    nthread_cond_t *cond)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_cond_init(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_cond_init(&cond->cond, NULL);
#elif HAVE_TRIOS_MTA
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_init(threadid=%lu): initializing cond(%p), bound_mutex(%p)\n",
            mta_get_threadid(), cond, cond->bound_mutex);
    fflush(stderr);
#endif
#if USING_MTGL
    mtgl::mtgl_qthread_init();
    mtgl::mt_purge(cond->cond);
#else
    purge(&cond->cond);
#endif
    cond->waiter_count=0;
    cond->bound_mutex=NULL;
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_init(threadid=%lu): initialized cond(%p), bound_mutex(%p)\n",
            mta_get_threadid(), cond, cond->bound_mutex);
    fflush(stderr);
#endif
#endif
    return(rc);
}

int nthread_cond_broadcast(nthread_cond_t *cond)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_cond_broadcast(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_cond_broadcast(&cond->cond);
#elif HAVE_TRIOS_MTA
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_broadcast(threadid=%lu): bcasting cond(%p), bound_mutex(%p) waiter_count(%d)\n",
            mta_get_threadid(), cond, cond->bound_mutex, cond->waiter_count);
    fflush(stderr);
#endif
    if (cond->waiter_count>0) {
        assert( cond->bound_mutex!=NULL );
        nthread_unlock(cond->bound_mutex);
#if USING_MTGL
        mtgl::mt_write(cond->cond, 0);
#else
        writeef(&cond->cond, 0);
#endif
    }
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_broadcast(threadid=%lu): bcasted cond(%p), bound_mutex(%p)\n",
            mta_get_threadid(), cond, cond->bound_mutex);
    fflush(stderr);
#endif
#endif
    return(rc);
}

int nthread_cond_signal(nthread_cond_t *cond)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_cond_signal(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_cond_signal(&cond->cond);
#elif HAVE_TRIOS_MTA
    rc = nthread_cond_broadcast(cond);
#endif
    return(rc);
}

int nthread_cond_timedwait(
        nthread_cond_t *cond,
        nthread_mutex_t *mutex,
        const struct timespec *abstime)
{
    int rc=0;

//    log_debug(thread_debug_level, "nthread_cond_timedwait(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_cond_timedwait(&cond->cond, &mutex->mutex, abstime);
#elif HAVE_TRIOS_MTA
    rc = nthread_cond_wait(cond, mutex);
#endif
    return(rc);
}

int nthread_cond_wait(
        nthread_cond_t *cond,
        nthread_mutex_t *mutex)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_cond_wait(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = pthread_cond_wait(&cond->cond, &mutex->mutex);
#elif HAVE_TRIOS_MTA
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_wait(threadid=%lu): waiting cond(%p), mutex(%p), bound_mutex(%p) waiter_count(%d)\n",
            mta_get_threadid(), cond, mutex, cond->bound_mutex, cond->waiter_count);
    fflush(stderr);
#endif
    assert( cond->bound_mutex==NULL || cond->bound_mutex==mutex );

    if (cond->bound_mutex == NULL) {
        cond->bound_mutex = mutex;
    }
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_wait(threadid=%lu): waiting cond(%p), mutex(%p), bound_mutex(%p) waiter_count(%d)\n",
            mta_get_threadid(), cond, mutex, cond->bound_mutex, cond->waiter_count);
    fflush(stderr);
#endif
    cond->waiter_count++;
    nthread_unlock(mutex);
#if USING_MTGL
    mtgl::mt_readfe(cond->cond);
#else
    readfe(&cond->cond);
#endif
    nthread_lock(mutex);
    cond->waiter_count--;
#ifdef _DEBUG_COND_
    fprintf(stderr, "nthread_cond_wait(threadid=%lu): waited cond(%p), mutex(%p), bound_mutex(%p) waiter_count(%d)\n",
            mta_get_threadid(), cond, mutex, cond->bound_mutex, cond->waiter_count);
    fflush(stderr);
#endif
#endif
    return(rc);
}

int nthread_mutex_init(nthread_mutex_t * mutex, int mutex_type)
{
#ifdef HAVE_TRIOS_PTHREAD
    pthread_mutexattr_t  mutex_attr;

    pthread_mutexattr_init(&mutex_attr);

    switch (mutex_type) {
    case NTHREAD_MUTEX_NORMAL:
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_NORMAL);
        break;
    case NTHREAD_MUTEX_RECURSIVE:
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);
        break;
    case NTHREAD_MUTEX_ERRORCHECK:
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
        break;
    case NTHREAD_MUTEX_DEFAULT:
    default:
        pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_DEFAULT);
        break;
    }

    pthread_mutex_init(&mutex->mutex, &mutex_attr);
#elif HAVE_TRIOS_MTA
    switch (mutex_type) {
        case NTHREAD_MUTEX_NORMAL:
            mutex->is_recursive = 0;
            break;
        case NTHREAD_MUTEX_RECURSIVE:
            mutex->is_recursive = 1;
            break;
        case NTHREAD_MUTEX_ERRORCHECK:
            mutex->is_recursive = 0;
            break;
        case NTHREAD_MUTEX_DEFAULT:
        default:
            mutex->is_recursive = 0;
            break;
    }

    mutex->mutex        = 0;
    mutex->is_locked    = 0;
    mutex->locker       = 0;
    mutex->lock_count   = 0;
#endif

    return(0);
}

#ifdef HAVE_TRIOS_MTA
static int lock_held(nthread_mutex_t *m){
  return (m->lock_count!=0 && m->locker==nthread_self());
}
static int lock_notheld(nthread_mutex_t *m){
  return (m->lock_count==0 || m->locker!=nthread_self());
}
#endif

int nthread_lock(nthread_mutex_t *mutex)
{
    int rc=0;

#ifdef HAVE_TRIOS_PTHREAD
#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: locking mutex(%p), pthread_mutex(%p)\n",
            mutex, &mutex->mutex);
    fflush(logger_get_file());
#endif
    rc = pthread_mutex_lock(&mutex->mutex);
#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: locked mutex(%p), pthread_mutex(%p)\n",
            mutex, &mutex->mutex);
    fflush(logger_get_file());
#endif
#elif HAVE_TRIOS_MTA
#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock(threadid=%lu): locking mutex(%p), lock_count(%d), is_recursive(%u)\n",
            mta_get_threadid(), mutex, mutex->lock_count, mutex->is_recursive);
    fflush(stderr);
#endif

    nthread_id_t me = nthread_self();
    if (mutex->is_recursive) {
        if ((mutex->lock_count>0) && (mutex->locker == me)) {
            mutex->lock_count++;
        } else {
            uint8_t m = mutex->mutex; // lock the mutex for the caller
            assert(mutex->lock_count==0);
            mutex->locker = me;
            mutex->lock_count=1;
        }
    } else {
        uint8_t m = mutex->mutex; // lock the mutex for the caller
        assert(mutex->lock_count==0);
        mutex->locker = me;
        mutex->lock_count=1;
    }
#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_lock(threadid=%lu): locked mutex(%p), lock_count(%d), is_recursive(%u)\n",
            mta_get_threadid(), mutex, mutex->lock_count, mutex->is_recursive);
    fflush(stderr);
#endif
//#elif _LIBCATAMOUNT_
//#warn No threading on Catamount.  Stubs only.
//#else
//#error No threading defined.  Must define either HAVE_TRIOS_PTHREAD or HAVE_TRIOS_MTA.
#endif
//    log_debug(thread_debug_level, "nthread_mutex_lock(STUB)");
    return(rc);
}

int nthread_unlock(nthread_mutex_t *mutex)
{
    int rc=0;

#ifdef HAVE_TRIOS_PTHREAD
#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: unlocking mutex(%p), pthread_mutex(%p)\n",
            mutex, &mutex->mutex);
    fflush(logger_get_file());
#endif
    rc = pthread_mutex_unlock(&mutex->mutex);
#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: unlocked mutex(%p), pthread_mutex(%p)\n",
            mutex, &mutex->mutex);
    fflush(logger_get_file());
#endif
#elif HAVE_TRIOS_MTA
#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_unlock(threadid=%lu): unlocking mutex(%p), lock_count(%d), is_recursive(%u)\n",
            mta_get_threadid(), mutex, mutex->lock_count, mutex->is_recursive);
    fflush(stderr);
#endif

    nthread_id_t me = nthread_self();
    if ((mutex->locker == me) && (lock_held(mutex))) {

        assert( lock_held(mutex) );
        mutex->lock_count--;
        assert( mutex->lock_count==0 || mutex->is_recursive==1 );

        if( mutex->lock_count==0 ){
            mutex->mutex = 0; // unlock the mutex for the caller
        }
    }
#ifdef _DEBUG_LOCKS_
    fprintf(stderr, "nthread_unlock(threadid=%lu): unlocked mutex(%p), lock_count(%d), is_recursive(%u)\n",
            mta_get_threadid(), mutex, mutex->lock_count, mutex->is_recursive);
    fflush(stderr);
#endif
#endif
//    log_debug(thread_debug_level, "nthread_mutex_unlock(STUB)");
    return(rc);
}

nthread_id_t nthread_self(void)
{
    nthread_id_t rc=0;

//    log_debug(thread_debug_level, "nthread_self(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    rc = (nthread_id_t)pthread_self();
#elif HAVE_TRIOS_MTA
#if USING_MTGL
    rc = 1;
#else
    rc = mta_get_threadid();
#endif
#endif
    return(rc);
}

void nthread_yield(void)
{
#if defined(HAVE_TRIOS_PTHREAD_YIELD)
    log_debug(thread_debug_level, "nthread_yield(STUB) - pthread_yield");
    pthread_yield();
#elif defined (HAVE_TRIOS_PTHREAD_YIELD_NP)
    log_debug(thread_debug_level, "nthread_yield(STUB) - pthread_yield_np");
    pthread_yield_np();
#elif defined(HAVE_TRIOS_SCHED_YIELD) && !defined(__LIBCATAMOUNT__)
    log_debug(thread_debug_level, "nthread_yield(STUB) - sched_yield");
    sched_yield();
#else
    log_debug(thread_debug_level, "nthread_yield(STUB) - usleep");
    usleep(0);
#endif
    return;
}

int nthread_join(
        nthread_t thread,
        void **value_ptr)
{
    int rc=0;

    log_debug(thread_debug_level, "nthread_join(STUB)");
#ifdef HAVE_TRIOS_PTHREAD
    if (thread.thread) {
        rc = pthread_join(thread.thread, value_ptr);
    }
#elif HAVE_TRIOS_MTA
    log_debug(thread_debug_level, "**********************");
    log_debug(thread_debug_level, "nthread_join");
    log_debug(thread_debug_level, "**********************");
#endif
    return(rc);
}

int64_t nthread_counter_increment(nthread_counter_t *c)
{
    int64_t t=0;

#ifdef HAVE_TRIOS_MTA
#if USING_MTGL
    t = mtgl::mt_incr(&c->value, 1);
#else
    t = int_fetch_add(&c->value, 1);
#endif
#else
    log_debug(thread_debug_level, "nthread_counter_increment(STUB)");
    nthread_lock(&c->mutex);
    t = c->value;
    c->value++;
    nthread_unlock(&c->mutex);
#endif
    return t;
}

int64_t nthread_counter_decrement(nthread_counter_t *c)
{
    int64_t t=0;

#ifdef HAVE_TRIOS_MTA
#if USING_MTGL
    t = mtgl::mt_incr(&c->value, -1);
#else
    t = int_fetch_add(&c->value, -1);
#endif
#else
    log_debug(thread_debug_level, "nthread_counter_decrement(STUB)");
    nthread_lock(&c->mutex);
    t = c->value;
    c->value--;
    nthread_unlock(&c->mutex);
#endif

    return t;
}
