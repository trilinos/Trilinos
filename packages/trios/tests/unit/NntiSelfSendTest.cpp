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

NNTI_transport_t     trans_hdl;
NNTI_peer_t          server_hdl;
NNTI_buffer_t        queue_mr, send_mr; /* registered memory regions */
NNTI_result_t        err;
NNTI_status_t        queue_status, send_status;

pthread_barrier_t barrier;

static void *do_wait(void *args)
{
    NNTI_result_t rc;
    NNTI_status_t wait_status1;

    double wait_time=0.0;

    char *queue_buf=(char *)malloc(10*NNTI_REQUEST_BUFFER_SIZE);
    memset(queue_buf, 0, 10*NNTI_REQUEST_BUFFER_SIZE);
    NNTI_register_memory(&trans_hdl, queue_buf, NNTI_REQUEST_BUFFER_SIZE, 10, NNTI_RECV_QUEUE, NULL, &queue_mr);

    /* client is waiting for us to initialize */
    pthread_barrier_wait(&barrier);

    /* the client sends a message here */
    NNTI_wait(&queue_mr, NNTI_RECV_QUEUE, -1, &queue_status);

    pthread_barrier_wait(&barrier);

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
    char server_url[NNTI_URL_LEN];

    logger_init(LOG_ERROR, "nntiselfsend.log");

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

    char *send_buf=(char *)malloc(NNTI_REQUEST_BUFFER_SIZE);
    memset(send_buf, 0, NNTI_REQUEST_BUFFER_SIZE);
    rc=NNTI_register_memory(&trans_hdl, send_buf, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, NULL, &send_mr);

    rc=NNTI_send(&server_hdl, &send_mr, NULL);
    rc=NNTI_wait(&send_mr, NNTI_SEND_SRC, 5000, &send_status);

    pthread_barrier_wait(&barrier);

    NNTI_unregister_memory(&send_mr);
    free(send_buf);

    join_wait_thread();

    NNTI_fini(&trans_hdl);

    return 0;
}
