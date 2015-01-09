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
 *  Created on: March 20, 2014
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
NNTI_buffer_t        queue_mr, send_mr; /* registered memory regions */
NNTI_work_request_t  queue_wr, send_wr; /* registered memory regions */
NNTI_result_t        err;
NNTI_status_t        queue_status, send_status;

pthread_barrier_t barrier;

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

static void *do_wait(void *args)
{
    NNTI_result_t rc;
    selfsend_args *ssa;

    double wait_time=0.0;

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 10, NNTI_RECV_QUEUE, &queue_mr);

    NNTI_create_work_request(&queue_mr, &queue_wr);

    /* client is waiting for us to initialize */
    pthread_barrier_wait(&barrier);

    /* the client sends a message here */
    NNTI_wait(&queue_wr, -1, &queue_status);

    ssa=(selfsend_args *)queue_status.start;
    if (ssa->chksum != calc_checksum((const char *)&ssa->data, sizeof(data_t))) {
        success=false;
    }

    pthread_barrier_wait(&barrier);

    NNTI_destroy_work_request(&queue_wr);

    NNTI_free(&queue_mr);

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
    int rc = pthread_join(wait_thread, &status);
    if (rc) {
        fprintf (stdout, "ERROR: Cannot join wait thread, err code = %d\n", rc);
    } else {
        fprintf (stdout, "Joined wait thread.\n");
    }
}


int main(int argc, char *argv[])
{
    NNTI_result_t rc;
    selfsend_args *ssa;
    char server_url[NNTI_URL_LEN];

    logger_init(LOG_ERROR, NULL);

    pthread_barrier_init(&barrier, NULL, 2);

    rc=NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);
    rc=NNTI_get_url(&trans_hdl, server_url, NNTI_URL_LEN);

    launch_wait_thread();

    pthread_barrier_wait(&barrier);

    rc=NNTI_connect(
            &trans_hdl,
            server_url,
            5000,
            &server_hdl);

    rc=NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);

    ssa=(selfsend_args *)NNTI_BUFFER_C_POINTER(&send_mr);
    ssa->data.int_val   =10;
    ssa->data.float_val =10.0;
    ssa->data.double_val=10.0;
    ssa->chksum=calc_checksum((const char *)&ssa->data, sizeof(data_t));

    rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
    rc=NNTI_wait(&send_wr, 5000, &send_status);

    pthread_barrier_wait(&barrier);

    NNTI_free(&send_mr);

    join_wait_thread();

    NNTI_fini(&trans_hdl);

    if (success)
        std::cout << "\nEnd Result: TEST PASSED" << std::endl;
    else
        std::cout << "\nEnd Result: TEST FAILED" << std::endl;

    return (success ? 0 : 1 );
}
