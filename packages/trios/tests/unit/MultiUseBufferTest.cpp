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
 * NonContigTest.c
 *
 *  Created on: March 14, 2014
 *      Author: thkorde
 */

#include "Trios_nnti.h"

#include <unistd.h>
#include <memory.h>

#include <mpi.h>

#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nnti_fprint_types.h"

#include <sstream>
#include <ostream>


log_level multiuse_debug_level = LOG_UNDEFINED;


NNTI_transport_t     trans_hdl;
NNTI_peer_t          server_hdl;

NNTI_buffer_t        queue_mr;
NNTI_work_request_t  queue_wr;
NNTI_status_t        queue_status;

#define WR_COUNT 2
NNTI_work_request_t *mr_wr_list[WR_COUNT];
NNTI_result_t        err;
NNTI_status_t wait_status1;
NNTI_status_t wait_status2;
NNTI_status_t wait_status3;
NNTI_status_t *wait_status_list[WR_COUNT];

int one_kb=1024;
int one_mb=1024*1024;

#define NUM_SEGMENTS 10


/*
 * Can't use NNTI_REQUEST_BUFFER_SIZE or NNTI_RESULT_BUFFER_SIZE here,
 * because the encoded NNTI buffers are too big.  Define something
 * bigger here.
 */
#define NNTI_MULTIUSE_REQUEST_SIZE 2048
#define NNTI_MULTIUSE_RESULT_SIZE  2048



int client_test1(NNTI_buffer_t *local_multiuse_mr);
int client_test2(NNTI_buffer_t *local_multiuse_mr);
void server_test1(NNTI_buffer_t *local_multiuse_mr);
void server_test2(NNTI_buffer_t *local_multiuse_mr);



static inline uint64_t calc_checksum (char * buf, uint64_t size)
{
    unsigned long hash = 5381;
    char* p = buf;
    uint64_t i = 0;

    for (i=0; i<size; i++) {
        hash = ((hash << 5) + hash) + *p; /* hash * 33 + c */
        p++;
    }

    i = (uint64_t) hash;
    return i;
}

/*
 *  GET target
 */
int client_test1(NNTI_buffer_t *local_multiuse_mr)
{
    int success=TRUE;

    uint64_t *server_checksums=NULL;

    char *buf=NULL;

    int multiuse_mr_size;
    int recv_mr_size;

    NNTI_buffer_t        send_mr, recv_mr; /* contiguous send/recv registered memory regions */
    NNTI_work_request_t  send_wr, recv_wr; /* work requests */

    XDR send_xdrs;


    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_RESULT_SIZE, 1, NNTI_RECV_DST, &recv_mr);

    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_mb), 'A'+i, one_mb);
    }

    /* XDR encode the get, put and recv buffers */
    multiuse_mr_size  = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, local_multiuse_mr);
    recv_mr_size = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, &recv_mr);

    log_debug(multiuse_debug_level, "multiuse_mr_size=%d ; recv_mr_size=%d", multiuse_mr_size, recv_mr_size);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&send_mr),
            NNTI_BUFFER_SIZE(&send_mr), XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, local_multiuse_mr);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&send_mr)+multiuse_mr_size,
            NNTI_BUFFER_SIZE(&send_mr)-multiuse_mr_size, XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, &recv_mr);

    NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
    NNTI_wait(&send_wr, 1000, &wait_status1);

    NNTI_create_work_request(&recv_mr, &recv_wr);
    NNTI_wait(&recv_wr, 1000, &wait_status1);

    server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&recv_mr);

    /* checksums GET */
    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
        for (int j=0;j<5;j++) {
            out << (buf + (i*one_mb))[j];
        }
        log_debug(multiuse_debug_level, "client GET segments[%d]=%s", i, out.str().c_str());

        uint64_t my_checksum=calc_checksum(buf + (i*one_mb), one_mb);
        log_debug(multiuse_debug_level, "server_checksums[%d]=%X ; my_checksum=%X", i, server_checksums[i], my_checksum);
        if (server_checksums[i] != my_checksum) {
            success=FALSE;
        }
    }

    NNTI_free(&send_mr);
    NNTI_free(&recv_mr);

    log_debug(multiuse_debug_level, "exit");

    return success;
}

/*
 *  PUT target
 */
int client_test2(NNTI_buffer_t *local_multiuse_mr)
{
    int success=TRUE;

    uint64_t *server_checksums=NULL;

    char *buf=NULL;

    int multiuse_mr_size;
    int recv_mr_size;

    NNTI_buffer_t        send_mr, recv_mr; /* contiguous send/recv registered memory regions */
    NNTI_work_request_t  send_wr, recv_wr; /* work requests */

    XDR send_xdrs;


    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_RESULT_SIZE, 1, NNTI_RECV_DST, &recv_mr);

    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_mb), 'B'+i, one_mb);
    }

    /* XDR encode the get, put and recv buffers */
    multiuse_mr_size  = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, local_multiuse_mr);
    recv_mr_size = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, &recv_mr);

    log_debug(multiuse_debug_level, "multiuse_mr_size=%d ; recv_mr_size=%d", multiuse_mr_size, recv_mr_size);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&send_mr),
            NNTI_BUFFER_SIZE(&send_mr), XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, local_multiuse_mr);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&send_mr)+multiuse_mr_size,
            NNTI_BUFFER_SIZE(&send_mr)-multiuse_mr_size, XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, &recv_mr);

    NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
    NNTI_wait(&send_wr, 1000, &wait_status1);

    NNTI_create_work_request(&recv_mr, &recv_wr);
    NNTI_wait(&recv_wr, 1000, &wait_status1);

    server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&recv_mr);

    /* checksums PUT */
    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
        for (int j=0;j<5;j++) {
            out << (buf + (i*one_mb))[j];
        }
        log_debug(multiuse_debug_level, "client PUT segments[%d]=%s", i, out.str().c_str());

        uint64_t my_checksum=calc_checksum(buf + (i*one_mb), one_mb);
        log_debug(multiuse_debug_level, "server_checksums[%d]=%X ; my_checksum=%X", i, server_checksums[i], my_checksum);
        if (server_checksums[i] != my_checksum) {
            success=FALSE;
        }
    }

    NNTI_free(&send_mr);
    NNTI_free(&recv_mr);

    log_debug(multiuse_debug_level, "exit");

    return success;
}

int client(void)
{
    int success=TRUE;
    char url[NNTI_URL_LEN];

    NNTI_buffer_t local_multiuse_mr;

    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 10, NNTI_RECV_QUEUE, &queue_mr);

    NNTI_alloc(&trans_hdl, NUM_SEGMENTS*one_mb, 1, (NNTI_buf_ops_t)(NNTI_BOP_LOCAL_READ|NNTI_BOP_LOCAL_WRITE|NNTI_BOP_REMOTE_READ|NNTI_BOP_REMOTE_WRITE), &local_multiuse_mr);

    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    NNTI_connect(&trans_hdl, url, 5000, &server_hdl);

    /* I play the target role here */
    success = client_test1(&local_multiuse_mr);
    if (success == FALSE) {
        fprintf(stdout, "TEST #1 failed.  Aborting...\n");
        goto out;
    }
    fprintf(stdout, "TEST #1 passed\n");

    /* I play the target role here */
    success = client_test2(&local_multiuse_mr);
    if (success == FALSE) {
        fprintf(stdout, "TEST #2 failed.  Aborting...\n");
        goto out;
    }
    fprintf(stdout, "TEST #2 passed\n");

    NNTI_get_url(&trans_hdl, url, NNTI_URL_LEN);
    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 1, MPI_COMM_WORLD);
    log_debug(multiuse_debug_level, "multiuse client url is %s", url);

    /* I play the initiator role here */
    server_test1(&local_multiuse_mr);

    /* I play the initiator role here */
    server_test2(&local_multiuse_mr);

out:
    NNTI_free(&local_multiuse_mr);

    NNTI_free(&queue_mr);

    log_debug(multiuse_debug_level, "exit");

    return(success);
}

/*
 *  GET initiator
 */
void server_test1(NNTI_buffer_t *local_multiuse_mr)
{
    char *buf;

    NNTI_buffer_t        remote_multiuse_mr; /* contiguous registered memory regions */
    NNTI_work_request_t  multiuse_wr; /* work requests */

    NNTI_buffer_t        send_mr, recv_mr; /* contiguous send/recv registered memory regions */
    NNTI_work_request_t  send_wr;          /* work requests */

    XDR recv_xdrs;


    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 1, NNTI_SEND_SRC, &send_mr);

    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_mb), 'K'+i, one_mb);
    }

    NNTI_create_work_request(&queue_mr, &queue_wr);
    NNTI_wait(&queue_wr, -1, &queue_status);


    /* XDR decode the get, put and recv buffers */
    xdrmem_create(&recv_xdrs, (char *)queue_status.start+queue_status.offset,
            NNTI_MULTIUSE_REQUEST_SIZE, XDR_DECODE);

    memset(&remote_multiuse_mr, 0, sizeof(NNTI_buffer_t));
    memset(&recv_mr, 0, sizeof(NNTI_buffer_t));
    xdr_NNTI_buffer_t(&recv_xdrs, &remote_multiuse_mr);
    xdr_NNTI_buffer_t(&recv_xdrs, &recv_mr);

    if (logging_debug(multiuse_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "remote_multiuse_mr",
                "after XDR decode", &remote_multiuse_mr);
        fprint_NNTI_buffer(logger_get_file(), "recv_mr",
                "after XDR decode", &recv_mr);

        fprint_NNTI_buffer(logger_get_file(), "local_multiuse_mr",
                "before PUT", local_multiuse_mr);
    }

    NNTI_get(&remote_multiuse_mr, 0, NUM_SEGMENTS*one_mb, local_multiuse_mr, 0, &multiuse_wr);
    NNTI_wait(&multiuse_wr, 1000, &wait_status1);

    uint64_t *server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&send_mr);

    /* checksum GET */
    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
        for (int j=0;j<5;j++) {
            out << (buf + (i*one_mb))[j];
        }
        log_debug(multiuse_debug_level, "server GET segments[%d]=%s", i, out.str().c_str());

        server_checksums[i]=calc_checksum((buf + (i*one_mb)), one_mb);
        log_debug(multiuse_debug_level, "server GET checksums[%d]=%X", i, server_checksums[i]);
    }

    NNTI_send(&queue_status.src, &send_mr, &recv_mr, &send_wr);
    NNTI_wait(&send_wr, 1000, &wait_status1);

    NNTI_free(&send_mr);

    log_debug(multiuse_debug_level, "exit");
}

/*
 *  PUT initiator
 */
void server_test2(NNTI_buffer_t *local_multiuse_mr)
{
    char *buf;

    NNTI_buffer_t        remote_multiuse_mr; /* contiguous registered memory regions */
    NNTI_work_request_t  multiuse_wr; /* work requests */

    NNTI_buffer_t        send_mr, recv_mr; /* contiguous send/recv registered memory regions */
    NNTI_work_request_t  send_wr;          /* work requests */

    XDR recv_xdrs;


    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 1, NNTI_SEND_SRC, &send_mr);

    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_mb), 'L'+i, one_mb);
    }

    NNTI_create_work_request(&queue_mr, &queue_wr);
    NNTI_wait(&queue_wr, -1, &queue_status);


    /* XDR decode the get, put and recv buffers */
    xdrmem_create(&recv_xdrs, (char *)queue_status.start+queue_status.offset,
            NNTI_MULTIUSE_REQUEST_SIZE, XDR_DECODE);

    memset(&remote_multiuse_mr, 0, sizeof(NNTI_buffer_t));
    memset(&recv_mr, 0, sizeof(NNTI_buffer_t));
    xdr_NNTI_buffer_t(&recv_xdrs, &remote_multiuse_mr);
    xdr_NNTI_buffer_t(&recv_xdrs, &recv_mr);

    if (logging_debug(multiuse_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "remote_multiuse_mr",
                "after XDR decode", &remote_multiuse_mr);
        fprint_NNTI_buffer(logger_get_file(), "recv_mr",
                "after XDR decode", &recv_mr);

        fprint_NNTI_buffer(logger_get_file(), "local_multiuse_mr",
                "before PUT", local_multiuse_mr);
    }

    NNTI_put(local_multiuse_mr, 0, NUM_SEGMENTS*one_mb, &remote_multiuse_mr, 0, &multiuse_wr);
    NNTI_wait(&multiuse_wr, 1000, &wait_status1);

    uint64_t *server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&send_mr);

    /* checksum PUT */
    buf=NNTI_BUFFER_C_POINTER(local_multiuse_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
        for (int j=0;j<5;j++) {
            out << (buf + (i*one_mb))[j];
        }
        log_debug(multiuse_debug_level, "server PUT segments[%d]=%s", i, out.str().c_str());

        server_checksums[i]=calc_checksum((buf + (i*one_mb)), one_mb);
        log_debug(multiuse_debug_level, "server PUT checksums[%d]=%X", i, server_checksums[i]);
    }

    NNTI_send(&queue_status.src, &send_mr, &recv_mr, &send_wr);
    NNTI_wait(&send_wr, 1000, &wait_status1);

    NNTI_free(&send_mr);

    log_debug(multiuse_debug_level, "exit");
}

void server(void)
{
    int success=TRUE;
    char url[NNTI_URL_LEN];

    NNTI_buffer_t local_multiuse_mr;

    log_debug(multiuse_debug_level, "enter");

    NNTI_alloc(&trans_hdl, NNTI_MULTIUSE_REQUEST_SIZE, 10, NNTI_RECV_QUEUE, &queue_mr);

    NNTI_alloc(&trans_hdl, NUM_SEGMENTS*one_mb, 1, (NNTI_buf_ops_t)(NNTI_BOP_LOCAL_READ|NNTI_BOP_LOCAL_WRITE|NNTI_BOP_REMOTE_READ|NNTI_BOP_REMOTE_WRITE), &local_multiuse_mr);

    NNTI_get_url(&trans_hdl, url, NNTI_URL_LEN);
    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);
    log_debug(multiuse_debug_level, "multiuse server url is %s", url);

    /* I play the initiator role here */
    server_test1(&local_multiuse_mr);

    /* I play the initiator role here */
    server_test2(&local_multiuse_mr);

    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 1, MPI_COMM_WORLD);
    NNTI_connect(&trans_hdl, url, 5000, &server_hdl);

    /* I play the target role here */
    success = client_test1(&local_multiuse_mr);
    if (success == FALSE) {
        fprintf(stdout, "TEST #3 failed.  Aborting...\n");
        goto out;
    }
    fprintf(stdout, "TEST #3 passed\n");

    /* I play the target role here */
    success = client_test2(&local_multiuse_mr);
    if (success == FALSE) {
        fprintf(stdout, "TEST #4 failed.  Aborting...\n");
        goto out;
    }
    fprintf(stdout, "TEST #4 passed\n");

out:
    NNTI_free(&local_multiuse_mr);

    NNTI_free(&queue_mr);

    log_debug(multiuse_debug_level, "exit");

    return;
}

NNTI_transport_id_t get_transport_from_env()
{
    char                *transport=NULL;
    NNTI_transport_id_t  trans_id=NNTI_DEFAULT_TRANSPORT;

    transport=getenv("NNTI_TRANSPORT");
    if (transport != NULL) {
        if (!strcmp(transport, "GEMINI")) {
            trans_id = NNTI_TRANSPORT_GEMINI;
        } else if (!strcmp(transport, "IB")) {
            trans_id = NNTI_TRANSPORT_IB;
        } else if (!strcmp(transport, "MPI")) {
            trans_id = NNTI_TRANSPORT_MPI;
        }
    }

    return(trans_id);
}

int main(int argc, char *argv[])
{
	int success=0;
    int nprocs, rank;

    NNTI_transport_id_t trans_id=NNTI_DEFAULT_TRANSPORT;

    char logname[1024];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sprintf(logname, "multiuse.%03d.log", rank);
    logger_init(LOG_ERROR, NULL /*logname*/);

    trans_id = get_transport_from_env();

    NNTI_init(trans_id, NULL, &trans_hdl);

    if (rank==0) {
        server();
    } else {
    	success=client();
    }

    MPI_Bcast(&success, 1, MPI_INT, 1, MPI_COMM_WORLD);

    NNTI_fini(&trans_hdl);

    MPI_Finalize();

    logger_fini();

    if (rank == 1) {
		if (success==TRUE)
			fprintf(stdout, "\nEnd Result: TEST PASSED\n");
		else
			fprintf(stdout, "\nEnd Result: TEST FAILED\n");
    }

    return (success==TRUE ? 0 : 1 );
}
