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
/**
 * TRIOS does not have any internal threading, but it should run in a
 * multithreaded environment.  The nthread library provides a single API
 * for locks and atomic counters.
 */

#include "Trios_config.h"

#include <fcntl.h>
#include <sys/stat.h>

#ifdef HAVE_TRIOS_MALLOC_H
#include <malloc.h>
#endif

#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

#include <errno.h>
#include <string.h>

#if defined(HAVE_TRIOS_PTHREAD_H)
#include <pthread.h>
#endif

#if defined(HAVE_TRIOS_SEMAPHORE_H)
#include <semaphore.h>
#endif

#include "Trios_logger.h"
#include "Trios_threads.h"

#define _DEBUG_LOCKS_
#undef _DEBUG_LOCKS_

/* set to LOG_UNDEFINED -- log commands will use default level */
log_level thread_debug_level = LOG_UNDEFINED;


int nthread_lock_init(
        nthread_lock_t *lock)
{
    int  rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock_init: initializing lock(%p)\n", (void*)lock);
    fflush(logger_get_file());
#endif

#if defined(HAVE_TRIOS_PTHREAD_MUTEX_INIT)
    pthread_mutex_init(&lock->lock, NULL);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES)
    rc=sem_init(&lock->lock, 0, 1);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_lock_init: sem_init failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        lock->lock_ptr=NULL;
        return(-1);
    }
    lock->lock_ptr=&lock->lock;
#elif defined(HAVE_TRIOS_NAMED_SEMAPHORES)
    bool done=false;
    do {
        lock->name=tempnam("/tmp", "trios.");
        lock->lock_ptr=sem_open(lock->name+4, O_CREAT|O_EXCL, 0600, 1);
        if ((lock->lock_ptr == SEM_FAILED) && (errno == EEXIST)) {
            done=false;
        } else {
            done=true;
        }
    } while (!done);

    if (lock->lock_ptr == SEM_FAILED) {
        fprintf(logger_get_file(), "nthread_lock_init: sem_open failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        lock->lock_ptr=NULL;
        return(-1);
    }
#else
#warning No locking mechanism available on this system.
#endif

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock_init: initialized lock(%p), lock->name(%s)\n", (void*)lock, lock->name);
    fflush(logger_get_file());
#endif

    return(rc);
}

int nthread_lock(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: locking lock(%p)\n", (void*)lock);
    fflush(logger_get_file());
#endif

#if defined(HAVE_TRIOS_PTHREAD_MUTEX_LOCK)
    pthread_mutex_lock(&lock->lock);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
    if (lock->lock_ptr == NULL) {
        fprintf(logger_get_file(), "nthread_lock: lock not initialized\n");
        fflush(logger_get_file());
        return(-1);
    }

    rc=sem_wait(lock->lock_ptr);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_lock: sem_wait failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }
#else
#warning No locking mechanism available on this system.
#endif

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock: locked lock(%p)\n", (void*)lock);
    fflush(logger_get_file());
#endif

    return(rc);
}

int nthread_unlock(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_unlock: unlocking lock(%p)\n", (void*)lock);
    fflush(logger_get_file());
#endif

#if defined(HAVE_TRIOS_PTHREAD_MUTEX_UNLOCK)
    pthread_mutex_unlock(&lock->lock);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
    if (lock->lock_ptr == NULL) {
        fprintf(logger_get_file(), "nthread_unlock: lock not initialized\n");
        fflush(logger_get_file());
        return(-1);
    }

    rc=sem_post(lock->lock_ptr);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_unlock: sem_post failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }
#else
#warning No locking mechanism available on this system.
#endif

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_unlock: unlocked lock(%p)\n", (void*)lock);
    fflush(logger_get_file());
#endif

    return(rc);
}


int nthread_lock_fini(
        nthread_lock_t *lock)
{
    int rc=0;

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock_fini: finalizing lock(%p), lock->name(%s)", (void*)lock, lock->name);
    fflush(logger_get_file());
#endif

#if defined(HAVE_TRIOS_PTHREAD_MUTEX_DESTROY)
    // do nothing
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
    if (lock->lock_ptr == NULL) {
        fprintf(logger_get_file(), "nthread_lock_fini: lock not initialized\n");
        fflush(logger_get_file());
        return(-1);
    }
#endif

#if defined(HAVE_TRIOS_PTHREAD_MUTEX_DESTROY)
    pthread_mutex_destroy(&lock->lock);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES)
    rc=sem_destroy(lock->lock_ptr);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_lock_fini: sem_destroy failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }
#elif defined(HAVE_TRIOS_NAMED_SEMAPHORES)
    rc=sem_close(lock->lock_ptr);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_lock_fini: sem_close failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }
    rc=sem_unlink(lock->name+4);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_lock_fini: sem_unlink failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }

    free(lock->name);
    lock->name=NULL;
#else
#warning No locking mechanism available on this system.
#endif

#ifdef _DEBUG_LOCKS_
    fprintf(logger_get_file(), "nthread_lock_fini: finalized lock(%p), lock->name(%s)", (void*)lock, lock->name);
    fflush(logger_get_file());
#endif

    return(rc);
}


int nthread_cond_init(
        nthread_cond_t *condvar)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_INIT)
    pthread_cond_init(&condvar->condvar, NULL);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}

int nthread_wait(
        nthread_cond_t *condvar,
        nthread_lock_t *lock)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_WAIT)
    pthread_cond_wait(&condvar->condvar, &lock->lock);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}

int nthread_timedwait(
        nthread_cond_t *condvar,
        nthread_lock_t *lock,
        uint64_t        timeout)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_TIMEDWAIT)
    struct timespec abstime;

    int  sec=timeout/1000;
    long ms =timeout-sec*1000;

    clock_gettime(CLOCK_REALTIME, &abstime);

    // perform the addition
    abstime.tv_nsec+=ms*1000000;
    // adjust the time
    abstime.tv_sec+=abstime.tv_nsec/1000000000 + sec;
    abstime.tv_nsec=abstime.tv_nsec%1000000000;

    rc=pthread_cond_timedwait(&condvar->condvar, &lock->lock, &abstime);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}

int nthread_signal(
        nthread_cond_t *condvar)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_SIGNAL)
    pthread_cond_signal(&condvar->condvar);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}

int nthread_broadcast(
        nthread_cond_t *condvar)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_BROADCAST)
    pthread_cond_broadcast(&condvar->condvar);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}

int nthread_cond_fini(
        nthread_cond_t *condvar)
{
    int rc=0;

#if defined(HAVE_TRIOS_PTHREAD_COND_DESTROY)
    pthread_cond_destroy(&condvar->condvar);
#elif defined(HAVE_TRIOS_UNNAMED_SEMAPHORES) || defined(HAVE_TRIOS_NAMED_SEMAPHORES)
//#warning The semaphores implementation of NSSI threads does not have conditional variables.  Sorry.
#endif

    return(rc);
}


int nthread_counter_init(
        nthread_counter_t *c)
{
    int rc=0;

//    log_debug(thread_debug_level, "nthread_counter_init(STUB)");
    rc=nthread_lock_init(&c->lock);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_counter_init: nthread_lock_init failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }

    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_init: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        goto cleanup;
    }

    c->value = 0;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_init: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

cleanup:
    return rc;
}

int64_t nthread_counter_increment(
        nthread_counter_t *c)
{
    int64_t t=0;

//    log_debug(thread_debug_level, "nthread_counter_increment(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_increment: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        t = -1;
        goto cleanup;
    }

    t = c->value;
    c->value++;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_increment: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

cleanup:
    return t;
}

int64_t nthread_counter_decrement(
        nthread_counter_t *c)
{
    int64_t t=0;

//    log_debug(thread_debug_level, "nthread_counter_decrement(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_decrement: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        t = -1;
        goto cleanup;
    }

    t = c->value;
    c->value--;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_decrement: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

cleanup:
    return t;
}

int64_t nthread_counter_read(
        nthread_counter_t *c)
{
    int64_t t=0;

//    log_debug(thread_debug_level, "nthread_counter_read(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_read: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        t = -1;
        goto cleanup;
    }

    t = c->value;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_read: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

cleanup:
    return t;
}

int64_t nthread_counter_set(
        nthread_counter_t *c,
        int64_t            new_value)
{
    int64_t t=0;

//    log_debug(thread_debug_level, "nthread_counter_set(STUB)");

    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_set: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        t = -1;
        goto cleanup;
    }

    t = c->value;
    c->value=new_value;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_set: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

cleanup:
    return t;
}

int nthread_counter_fini(
        nthread_counter_t *c)
{
    int rc=0;

//    log_debug(thread_debug_level, "nthread_counter_fini(STUB)");
    if (nthread_lock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_fini: failed to lock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
        goto cleanup;
    }

    c->value = 0;

    if (nthread_unlock(&c->lock) != 0) {
        fprintf(logger_get_file(), "nthread_counter_fini: failed to unlock the counter lock: %s\n", strerror(errno));
        fflush(logger_get_file());
    }

    rc=nthread_lock_fini(&c->lock);
    if (rc == -1) {
        fprintf(logger_get_file(), "nthread_counter_fini: nthread_lock_fini failed: %s\n", strerror(errno));
        fflush(logger_get_file());
        return(-1);
    }

cleanup:
    return rc;
}
