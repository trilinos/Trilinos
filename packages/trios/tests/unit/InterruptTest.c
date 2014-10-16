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
/*
 * InterruptTest.c
 *
 *  Created on: March 7, 2014
 *      Author: thkorde
 */

#include "Trios_nnti.h"

#include <unistd.h>
#include <pthread.h>

#include "Trios_logger.h"
#include "Trios_timer.h"

NNTI_transport_t     trans_hdl;
NNTI_peer_t          server_hdl;
NNTI_buffer_t        mr1, mr2, mr3; /* registered memory regions */
NNTI_work_request_t  mr1_wr, mr2_wr, mr3_wr; /* work requests */
#define WR_COUNT 3
NNTI_work_request_t *mr_wr_list[WR_COUNT];
NNTI_result_t        err;
NNTI_status_t        status;

log_level interrupt_debug_level=LOG_UNDEFINED;

int expected_wait_time=0;
int actual_wait_time=0;

pthread_barrier_t barrier;

static int time_to_exit=FALSE;
static void *do_wait(void *args)
{
    NNTI_result_t rc;
    NNTI_status_t wait_status1;
    NNTI_status_t wait_status2;
    NNTI_status_t wait_status3;
    NNTI_status_t *wait_status_list[WR_COUNT];

    uint32_t which=0;

    double wait_time=0.0;

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_wait() #1 won't be interrupted so should timeout after 10000 ms.");

    wait_time=trios_get_time_ms();
    rc=NNTI_wait(&mr1_wr, 10000, &wait_status1);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_wait() #1: %lums", (uint64_t)wait_time);

    expected_wait_time+=10000;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);
    /* between these barriers the main thread is calling NNTI_interrupt() so that the next NNTI_wait() returns immediately. */
    pthread_barrier_wait(&barrier);

    log_debug(interrupt_debug_level, "NNTI_wait() #2 should return immediately.");

    wait_time=trios_get_time_ms();
    rc=NNTI_wait(&mr2_wr, 10000, &wait_status2);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_wait() #2: %lums", (uint64_t)wait_time);

    expected_wait_time+=0;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);

    log_debug(interrupt_debug_level, "NNTI_wait() #3 should return after 3000 ms.");

    wait_time=trios_get_time_ms();
    rc=NNTI_wait(&mr3_wr, 10000, &wait_status3);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_wait() #3: %lums", (uint64_t)wait_time);

    expected_wait_time+=3000;
    actual_wait_time+=wait_time;

    /* between these barriers the main thread is sleeping 3 seconds and then calling NNTI_interrupt(). */
    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitany() #1 won't be interrupted so should timeout after 10000 ms.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitany(mr_wr_list, WR_COUNT, 10000, &which, &wait_status1);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitany() #1: %lums", (uint64_t)wait_time);

    expected_wait_time+=10000;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);
    /* between these barriers the main thread is calling NNTI_interrupt() so that the next NNTI_waitany() returns immediately. */
    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitany() #2 should return immediately.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitany(mr_wr_list, WR_COUNT, 10000, &which, &wait_status2);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitany() #2: %lums", (uint64_t)wait_time);

    expected_wait_time+=0;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitany() #3 should return after 3000 ms.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitany(mr_wr_list, WR_COUNT, 10000, &which, &wait_status3);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitany() #3: %lums", (uint64_t)wait_time);

    expected_wait_time+=3000;
    actual_wait_time+=wait_time;

    /* between these barriers the main thread is sleeping 3 seconds and then calling NNTI_interrupt(). */
    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitall() #1 won't be interrupted so should timeout after 10000 ms.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_status_list[0]=&wait_status1;
    wait_status_list[1]=&wait_status2;
    wait_status_list[2]=&wait_status3;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitall(mr_wr_list, WR_COUNT, 10000, wait_status_list);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitall() #1: %lums", (uint64_t)wait_time);

    expected_wait_time+=10000;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);
    /* between these barriers the main thread is calling NNTI_interrupt() so that the next NNTI_waitall() returns immediately. */
    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitall() #2 should return immediately.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_status_list[0]=&wait_status1;
    wait_status_list[1]=&wait_status2;
    wait_status_list[2]=&wait_status3;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitall(mr_wr_list, WR_COUNT, 10000, wait_status_list);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitall() #2: %lums", (uint64_t)wait_time);

    expected_wait_time+=0;
    actual_wait_time+=wait_time;

    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    NNTI_create_work_request(&mr1, &mr1_wr);
    NNTI_create_work_request(&mr2, &mr2_wr);
    NNTI_create_work_request(&mr3, &mr3_wr);

    log_debug(interrupt_debug_level, "NNTI_waitall() #3 should return after 3000 ms.");

    mr_wr_list[0]=&mr1_wr;
    mr_wr_list[1]=&mr2_wr;
    mr_wr_list[2]=&mr3_wr;
    wait_status_list[0]=&wait_status1;
    wait_status_list[1]=&wait_status2;
    wait_status_list[2]=&wait_status3;
    wait_time=trios_get_time_ms();
    rc=NNTI_waitall(mr_wr_list, WR_COUNT, 10000, wait_status_list);
    wait_time=trios_get_time_ms()-wait_time;

    log_debug(interrupt_debug_level, "Time to complete NNTI_waitall() #3: %lums", (uint64_t)wait_time);

    expected_wait_time+=3000;
    actual_wait_time+=wait_time;

    /* between these barriers the main thread is sleeping 3 seconds and then calling NNTI_interrupt(). */
    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&mr1_wr);
    NNTI_destroy_work_request(&mr2_wr);
    NNTI_destroy_work_request(&mr3_wr);

    return(NULL);
}

static pthread_t wait_thread;
static void launch_wait_thread()
{
    /* Start up polling thread */
    fprintf (stdout, "Start wait thread.\n");
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    int rc = pthread_create(&wait_thread, &attr, do_wait, NULL);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot create polling thread, err code = %d\n", rc);
        return;
    }
    pthread_attr_destroy(&attr);
}
static void join_wait_thread()
{
    void *status;
    time_to_exit=TRUE;
    int rc = pthread_join(wait_thread, &status);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot join wait thread, err code = %d\n", rc);
    } else {
        fprintf (stdout, "Joined wait thread.\n");
    }
}


int main(int argc, char *argv[])
{
	int8_t success=TRUE;

    logger_init(LOG_ERROR, NULL/*"interrupt.log"*/);

    pthread_barrier_init(&barrier, NULL, 2);

    NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);

    NNTI_alloc(&trans_hdl, NNTI_RESULT_BUFFER_SIZE, 1, NNTI_RECV_DST, &mr1);
    NNTI_alloc(&trans_hdl, NNTI_RESULT_BUFFER_SIZE, 1, NNTI_RECV_DST, &mr2);
    NNTI_alloc(&trans_hdl, NNTI_RESULT_BUFFER_SIZE, 1, NNTI_RECV_DST, &mr3);

    launch_wait_thread();

    pthread_barrier_wait(&barrier);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_wait(). */
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_wait(). */
    sleep(3);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);

    pthread_barrier_wait(&barrier);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_waitany(). */
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_waitany(). */
    sleep(3);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);

    pthread_barrier_wait(&barrier);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_waitall(). */
    pthread_barrier_wait(&barrier);
    /* between these barriers the wait thread is calling NNTI_waitall(). */
    sleep(3);
    NNTI_interrupt(&trans_hdl);
    pthread_barrier_wait(&barrier);

    join_wait_thread();

    NNTI_free(&mr1);
    NNTI_free(&mr2);
    NNTI_free(&mr3);

    log_debug(interrupt_debug_level, "Time to complete NNTI_wait: expected=%lums ; actual=%lums", (uint64_t)expected_wait_time, (uint64_t)actual_wait_time);
    if ((actual_wait_time < (expected_wait_time-100)) ||
    	(actual_wait_time > (expected_wait_time+100))) {
        success=FALSE;
    }

    NNTI_fini(&trans_hdl);

    logger_fini();

    if (success)
        fprintf(stdout, "\nEnd Result: TEST PASSED\n");
    else
    	fprintf(stdout, "\nEnd Result: TEST FAILED\n");

    return (success ? 0 : 1 );
}
