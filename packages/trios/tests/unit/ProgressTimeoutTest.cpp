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
 * SelfSendTest.c
 *
 *  Created on: March 26, 2014
 *      Author: thkorde
 */

#include "Trios_nnti.h"

#include "Trios_logger.h"
#include "Trios_timer.h"

#include <unistd.h>
#include <string.h>
#include <pthread.h>

#include <iostream>

NNTI_transport_t     trans_hdl;
NNTI_peer_t          server_hdl;
NNTI_buffer_t        queue_mr, recv_mr, send_mr; /* registered memory regions */
NNTI_work_request_t  queue_wr, recv_wr, send_wr; /* registered memory regions */
NNTI_result_t        err;
NNTI_status_t        queue_status, recv_status, send_status;

pthread_barrier_t barrier2;
pthread_barrier_t barrier3;

bool success=true;

typedef struct {
    /** An integer value. */
    uint32_t int_val;
    /** A floating point value. */
    float float_val;
    /** A double value. */
    double double_val;
} data_t;

typedef struct {
    /* test data */
    data_t data;
    /** 32-bit checksum of the test data */
    uint64_t chksum;
} selfsend_args;

static uint64_t calc_checksum (const char * buf, const uint64_t size)
{
    unsigned long hash = 5381;
    const char* p = buf;
    uint64_t i = 0;

    for (i=0; i<size; i++) {
        hash = ((hash << 5) + hash) + *p; /* hash * 33 + c */
        p++;
    }

    i = (uint64_t) hash;
    return i;
}

static void *do_request_wait(void *args)
{
    NNTI_result_t rc;
    selfsend_args *ssa;

    double wait_time=0.0;

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 10, NNTI_RECV_QUEUE, &queue_mr);

    NNTI_create_work_request(&queue_mr, &queue_wr);

    /* client is waiting for us to initialize */
    pthread_barrier_wait(&barrier3);

    /* the client sends a message here */
    NNTI_wait(&queue_wr, -1, &queue_status);

    ssa=(selfsend_args *)queue_status.start;
    if (ssa->chksum != calc_checksum((const char *)&ssa->data, sizeof(data_t))) {
        fprintf(stdout, "checksum failure in the request thread");
        success=false;
    }

    pthread_barrier_wait(&barrier3);

    NNTI_destroy_work_request(&queue_wr);

    NNTI_free(&queue_mr);

    return(NULL);
}

#define RECV_TIMEOUT 3000
static void *do_recv_wait(void *args)
{
    NNTI_result_t rc;
    selfsend_args *ssa;

    double wait_time=0.0;

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_RECV_DST, &recv_mr);

    NNTI_create_work_request(&recv_mr, &recv_wr);

    /* client is waiting for us to initialize */
    pthread_barrier_wait(&barrier3);

    wait_time=trios_get_time_ms();
    /* the client sends a message here */
    NNTI_wait(&recv_wr, RECV_TIMEOUT, &recv_status);
    wait_time=trios_get_time_ms()-wait_time;

    // if wait time varies from the timeout by more than 100ms, then fail.
    if ((wait_time < (RECV_TIMEOUT-100)) || (wait_time > (RECV_TIMEOUT+100))) {
        fprintf(stdout, "Time to complete NNTI_wait: expected=%lums ; actual=%lums\n", (uint64_t)RECV_TIMEOUT, (uint64_t)wait_time);
        success=false;
    }

    pthread_barrier_wait(&barrier2);

    wait_time=trios_get_time_ms();
    /* the client sends a message here */
    NNTI_wait(&recv_wr, -1, &recv_status);
    wait_time=trios_get_time_ms()-wait_time;

    ssa=(selfsend_args *)recv_status.start;
    if (ssa->chksum != calc_checksum((const char *)&ssa->data, sizeof(data_t))) {
        fprintf(stdout, "checksum failure in the receive thread");
        success=false;
    }

    pthread_barrier_wait(&barrier3);

    NNTI_destroy_work_request(&recv_wr);

    NNTI_free(&recv_mr);

    return(NULL);
}

static pthread_t request_thread;
static pthread_t recv_thread;
static void launch_wait_threads()
{
    int rc=0;

    /* Start up polling threads */
    fprintf (stdout, "Start wait threads.\n");
    pthread_attr_t attr;

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    rc = pthread_create(&request_thread, &attr, do_request_wait, NULL);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot create request thread, err code = %d\n", rc);
        return;
    }
    rc = pthread_create(&recv_thread, &attr, do_recv_wait, NULL);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot create receive thread, err code = %d\n", rc);
        return;
    }
    pthread_attr_destroy(&attr);
}
static void join_wait_threads()
{
    int rc=0;

    void *status;
    rc = pthread_join(request_thread, &status);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot join request thread, err code = %d\n", rc);
    } else {
        fprintf (stdout, "Joined request thread.\n");
    }
    rc = pthread_join(recv_thread, &status);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot join receive thread, err code = %d\n", rc);
    } else {
        fprintf (stdout, "Joined receive thread.\n");
    }
}


int main(int argc, char *argv[])
{
    NNTI_result_t rc;
    selfsend_args *ssa;
    char server_url[NNTI_URL_LEN];

    logger_init(LOG_ERROR, NULL);

    pthread_barrier_init(&barrier2, NULL, 2);
    pthread_barrier_init(&barrier3, NULL, 3);

    rc=NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);
    rc=NNTI_get_url(&trans_hdl, server_url, NNTI_URL_LEN);

    launch_wait_threads();

    pthread_barrier_wait(&barrier3);

    rc=NNTI_connect(
            &trans_hdl,
            server_url,
            5000,
            &server_hdl);

    pthread_barrier_wait(&barrier2);

    rc=NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);

    ssa=(selfsend_args *)NNTI_BUFFER_C_POINTER(&send_mr);
    ssa->data.int_val   =10;
    ssa->data.float_val =10.0;
    ssa->data.double_val=10.0;
    ssa->chksum=calc_checksum((const char *)&ssa->data, sizeof(data_t));

    rc=NNTI_send(&server_hdl, &send_mr, &recv_mr, &send_wr);
    rc=NNTI_wait(&send_wr, 5000, &send_status);

    rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
    rc=NNTI_wait(&send_wr, 5000, &send_status);

    pthread_barrier_wait(&barrier3);

    NNTI_free(&send_mr);

    join_wait_threads();

    if (success)
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    else
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;

    return (success ? 0 : 1 );
}
