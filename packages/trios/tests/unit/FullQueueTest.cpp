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
 *  Created on: March 24, 2014
 *      Author: thkorde
 */

#include "Trios_nnti.h"
#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nnti_fprint_types.h"

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_oblackholestream.hpp"

#include <unistd.h>

#include <mpi.h>



#define SYNC_WITH_BARRIER



NNTI_transport_t trans_hdl;
NNTI_peer_t      server_hdl;
char             url[NNTI_URL_LEN];

NNTI_buffer_t       queue_mr;
NNTI_work_request_t queue_wr;

NNTI_buffer_t       client_ack_mr;
NNTI_work_request_t client_ack_wr;

NNTI_buffer_t       server_ack_mr;
NNTI_work_request_t server_ack_wr;

NNTI_buffer_t       send_mr; /* contiguous registered memory regions */
NNTI_work_request_t send_wr; /* work requests */

#define WR_COUNT 2
NNTI_work_request_t *mr_wr_list[WR_COUNT];
NNTI_result_t        err;

int one_mb=1024*1024;

log_level fullqueue_debug_level = LOG_UNDEFINED;

#define QUEUE_SIZE 10

int sent_count   =0;
int success_count=0;
int dropped_count=0;
int fill_count   =0;

static int buffer_pack(void *input, void **output, int32_t *output_size, xdrproc_t pack_func)
{
	XDR pack_xdrs;

	*output_size = xdr_sizeof(pack_func, input);
	*output=malloc(*output_size);
    xdrmem_create(&pack_xdrs, (caddr_t)*output, *output_size, XDR_ENCODE);
    pack_func(&pack_xdrs, input);

    return(0);
}

static int buffer_pack_free(void *input, int32_t input_size, xdrproc_t free_func)
{
	XDR free_xdrs;

    xdrmem_create(&free_xdrs, (caddr_t)input, input_size, XDR_FREE);
    free_func(&free_xdrs, input);

    return(0);
}

static int buffer_unpack(void *input, int32_t input_size, void *output, xdrproc_t unpack_func)
{
	XDR unpack_xdrs;

    xdrmem_create(&unpack_xdrs, (caddr_t)input, input_size, XDR_DECODE);
    unpack_func(&unpack_xdrs, output);

    return(0);
}

void client(void) {
    NNTI_result_t rc=NNTI_OK;
    NNTI_status_t send_status;
    NNTI_status_t client_ack_status;
    char *c_ptr;
    void    *packed=NULL;
    int32_t  packed_size=0;


    sent_count   =0;
    success_count=0;
    dropped_count=0;
    fill_count   =0;

    NNTI_connect(&trans_hdl, url, 5000, &server_hdl);

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_RECV_DST, &client_ack_mr);

    log_debug(LOG_ALL, "Client Begin Phase #1 - exchange ACK buffers");
    /*
     * Phase 1 - exchange ACK buffers
     */
    buffer_pack(&client_ack_mr, &packed, &packed_size, (xdrproc_t)&xdr_NNTI_buffer_t);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(fullqueue_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    // send the server the recv_mr so it can send back it's ack_mr
    c_ptr=NNTI_BUFFER_C_POINTER(&send_mr);
    memcpy(c_ptr, packed, packed_size);

    buffer_pack_free(packed, packed_size, (xdrproc_t)&xdr_NNTI_buffer_t);

    rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
    if (rc == NNTI_OK) {
        sent_count++;
    } else {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, 1000, &send_status);
    if (rc == NNTI_OK) {
        success_count++;
    } else {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // wait for the server to send back it's recv_mr
    NNTI_create_work_request(&client_ack_mr, &client_ack_wr);
    NNTI_wait(&client_ack_wr, -1, &client_ack_status);

    c_ptr=(char*)client_ack_status.start+client_ack_status.offset;
    buffer_unpack(c_ptr, client_ack_status.length, &server_ack_mr, (xdrproc_t)&xdr_NNTI_buffer_t);

    NNTI_destroy_work_request(&client_ack_wr);

//    fprint_NNTI_buffer(logger_get_file(), "server_ack_mr",
//            "received server ack hdl", &server_ack_mr);

#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    // send an ACK so the server knows to proceed
    rc=NNTI_send(&server_hdl, &send_mr, &server_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, 1000, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif


    log_debug(LOG_ALL, "Client Begin Phase #2 - the client fills the queue");
    /*
     * Phase 2 - fill the server's queue
     */
    // the server's queue has QUEUE_SIZE slots available, so these NNTI_send() operations should all succeed.
    for (int i=0;i<QUEUE_SIZE;i++) {
        rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
        if (rc == NNTI_OK) {
            sent_count++;
        } else {
            log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&send_wr, 1000, &send_status);
        if (rc == NNTI_OK) {
            success_count++;
            fill_count++;
        } else if (rc == NNTI_EDROPPED) {
        	dropped_count++;
        	break;
        } else {
            log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    c_ptr=NNTI_BUFFER_C_POINTER(&send_mr);
    *(int*)c_ptr=fill_count;

    // send an ACK so the server knows to proceed
    rc=NNTI_send(&server_hdl, &send_mr, &server_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, 1000, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }


    log_debug(LOG_ALL, "Client Begin Phase #3 - the client sends to a full queue");
    /*
     * Phase 3 - send to a full queue
     */
    // the server's queue should be full, so these NNTI_send() operations should fail with NNTI_EDROPPED.
    for (int i=0;i<QUEUE_SIZE;i++) {
        rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
        if (rc == NNTI_OK) {
            sent_count++;
        } else {
            log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&send_wr, 10000, &send_status);
        if (rc == NNTI_EDROPPED) {
            dropped_count++;
        } else {
            log_error(fullqueue_debug_level, "NNTI_wait() did not return NNTI_EDROPPED: %d", rc);
//            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    // send an ACK so the server knows to proceed
    rc=NNTI_send(&server_hdl, &send_mr, &server_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, 1000, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif


    log_debug(LOG_ALL, "Client Begin Phase #4 - the server drains the queue (fill_count=%d)", fill_count);
    /*
     * Phase 4 - the server drains the queue
     */
#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    // wait for the server to drain the queue
    NNTI_create_work_request(&client_ack_mr, &client_ack_wr);
    NNTI_wait(&client_ack_wr, -1, &client_ack_status);
    NNTI_destroy_work_request(&client_ack_wr);
#endif


    log_debug(LOG_ALL, "Client Begin Phase #5 - the client fills the queue again");
    /*
     * Phase 5 - fill the server's queue again
     */
    fill_count=0;
    // the server's queue is empty again, so these NNTI_send() operations should all succeed.
    for (int i=0;i<QUEUE_SIZE;i++) {
        rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr);
        if (rc == NNTI_OK) {
            sent_count++;
        } else {
            log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&send_wr, 10000, &send_status);
        if (rc == NNTI_OK) {
            success_count++;
            fill_count++;
        } else if (rc == NNTI_EDROPPED) {
        	break;
        } else {
            log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    c_ptr=NNTI_BUFFER_C_POINTER(&send_mr);
    *(int*)c_ptr=fill_count;

    // send an ACK so the server knows to proceed
    rc=NNTI_send(&server_hdl, &send_mr, &server_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, 1000, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    log_debug(LOG_ALL, "Client Begin Phase #6 - the server drains the queue again (fill_count=%d)", fill_count);
    /*
     * Phase 6 - the server drains the queue again
     */
#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    // wait for the server to drain the queue
    NNTI_create_work_request(&client_ack_mr, &client_ack_wr);
    NNTI_wait(&client_ack_wr, -1, &client_ack_status);
    NNTI_destroy_work_request(&client_ack_wr);
#endif

    NNTI_free(&send_mr);
    NNTI_free(&client_ack_mr);

    return;
}

void server(void)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_status_t queue_status;
    NNTI_status_t send_status;
    NNTI_status_t server_ack_status;
    char *c_ptr;
    void    *packed=NULL;
    int32_t  packed_size=0;


    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, QUEUE_SIZE, NNTI_RECV_QUEUE, &queue_mr);
    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_RECV_DST, &server_ack_mr);

    log_debug(LOG_ALL, "Server Begin Phase #1 - exchange ACK buffers");
    /*
     * Phase 1 - exchange ACK buffers
     */
    // wait for the client to send it's recv_mr
    NNTI_create_work_request(&queue_mr, &queue_wr);
    NNTI_wait(&queue_wr, -1, &queue_status);

    c_ptr=(char*)queue_status.start+queue_status.offset;
    buffer_unpack(c_ptr, queue_status.length, &client_ack_mr, (xdrproc_t)&xdr_NNTI_buffer_t);

    // we're done with the queue element
    NNTI_destroy_work_request(&queue_wr);

//    fprint_NNTI_buffer(logger_get_file(), "client_ack_mr",
//            "received client ack hdl", &client_ack_mr);

    buffer_pack(&server_ack_mr, &packed, &packed_size, (xdrproc_t)&xdr_NNTI_buffer_t);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(fullqueue_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    // send our recv_mr back to the client
    c_ptr=NNTI_BUFFER_C_POINTER(&send_mr);
    memcpy(c_ptr, packed, packed_size);

    buffer_pack_free(packed, packed_size, (xdrproc_t)&xdr_NNTI_buffer_t);

    rc=NNTI_send(&queue_status.src, &send_mr, &client_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, -1, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    NNTI_create_work_request(&server_ack_mr, &server_ack_wr);
    NNTI_wait(&server_ack_wr, -1, &server_ack_status);
    NNTI_destroy_work_request(&server_ack_wr);
#endif


    log_debug(LOG_ALL, "Server Begin Phase #2 - the client fills the queue");
    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 2 - the client fills the queue
     */
    NNTI_create_work_request(&server_ack_mr, &server_ack_wr);
    NNTI_wait(&server_ack_wr, -1, &server_ack_status);
    NNTI_destroy_work_request(&server_ack_wr);

    fill_count=*(int*)((char*)server_ack_status.start+server_ack_status.offset);


    log_debug(LOG_ALL, "Server Begin Phase #3 - the client sends to a full queue");
    /*
     * Phase 3 - the client sends to the full queue
     */
#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    NNTI_create_work_request(&server_ack_mr, &server_ack_wr);
    NNTI_wait(&server_ack_wr, -1, &server_ack_status);
    NNTI_destroy_work_request(&server_ack_wr);
#endif


    log_debug(LOG_ALL, "Server Begin Phase #4 - the server drains the queue (fill_count=%d)", fill_count);
    /*
     * Phase 4 - drain the queue
     */
    for (int i=0;i<fill_count;i++) {
        NNTI_create_work_request(&queue_mr, &queue_wr);
        rc=NNTI_wait(&queue_wr, 1000, &queue_status);
        NNTI_destroy_work_request(&queue_wr);
        if (rc != NNTI_OK) {
            log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    rc=NNTI_send(&queue_status.src, &send_mr, &client_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, -1, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif


    log_debug(LOG_ALL, "Server Begin Phase #5 - the client fills the queue again");
    /*
     * Phase 5 - the client fills the queue again
     */
    NNTI_create_work_request(&server_ack_mr, &server_ack_wr);
    NNTI_wait(&server_ack_wr, -1, &server_ack_status);
    NNTI_destroy_work_request(&server_ack_wr);

    fill_count=*(int*)((char*)server_ack_status.start+server_ack_status.offset);


    log_debug(LOG_ALL, "Server Begin Phase #6 - the server drains the queue again (fill_count=%d)", fill_count);
    /*
     * Phase 6 - drain the queue again
     */
    for (int i=0;i<fill_count;i++) {
        NNTI_create_work_request(&queue_mr, &queue_wr);
        rc=NNTI_wait(&queue_wr, 1000, &queue_status);
        NNTI_destroy_work_request(&queue_wr);
        if (rc != NNTI_OK) {
            log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

#ifdef SYNC_WITH_BARRIER
    MPI_Barrier(MPI_COMM_WORLD);
#else
    rc=NNTI_send(&queue_status.src, &send_mr, &client_ack_mr, &send_wr);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr, -1, &send_status);
    if (rc != NNTI_OK) {
        log_error(fullqueue_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
#endif


    NNTI_free(&queue_mr);
    NNTI_free(&send_mr);
    NNTI_free(&server_ack_mr);

    return;
}

int main(int argc, char *argv[])
{
    int nprocs, rank;

    char logname[1024];

    bool success=true;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank==0) {
    	sprintf(logname, "fullqueue.%03d.log", rank);
        logger_init(LOG_ERROR, NULL /*logname*/);
    } else {
    	sprintf(logname, "fullqueue.%03d.log", rank);
        logger_init(LOG_ERROR, NULL /*logname*/);
    }

    NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);

    if (rank==0) {
        NNTI_get_url(&trans_hdl, url, NNTI_URL_LEN);
    }

    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

    log_debug(fullqueue_debug_level, "full queue server url is %s", url);

    if (rank==0) {
        server();
    } else {
        client();

        log_debug(fullqueue_debug_level, "sent_count=%d, success_count=%d, fill_count=%d, dropped_count=%d", sent_count, success_count, fill_count, dropped_count);

        if (sent_count != (success_count+dropped_count)) {
            success=false;
        }

//        if (sent_count != (1 + (3*QUEUE_SIZE))) {
//            success=false;
//        }
//        if (success_count != (1 + (2*QUEUE_SIZE))) {
//            success=false;
//        }
//        if (dropped_count != (1*QUEUE_SIZE)) {
//            success=false;
//        }
    }

    NNTI_fini(&trans_hdl);

    MPI_Finalize();

    logger_fini();

    Teuchos::oblackholestream blackhole;
    std::ostream &out = ( rank == 1 ? std::cout : blackhole );
    if (success)
        out << "\nEnd Result: TEST PASSED" << std::endl;
    else
        out << "\nEnd Result: TEST FAILED" << std::endl;

    return (success ? 0 : 1 );
}
