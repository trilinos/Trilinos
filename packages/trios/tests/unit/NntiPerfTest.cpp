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
 * NntiPerfTest.c
 *
 *  Created on: May 16, 2014
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


int nprocs, nclients, nservers;
int rank, client_rank;


NNTI_transport_t trans_hdl;
NNTI_peer_t      server_hdl;
char             url[NNTI_URL_LEN];

NNTI_buffer_t queue_mr;
NNTI_buffer_t client_ack_mr;
NNTI_buffer_t server_ack_mr;
NNTI_buffer_t send_mr;

NNTI_work_request_t  queue_wr;
NNTI_work_request_t  client_ack_wr;
NNTI_work_request_t  server_ack_wr;
NNTI_work_request_t *send_wr;  // need an array of work requests because we will be issuing multiple async sends

NNTI_buffer_t get_src_mr;
NNTI_buffer_t get_dst_mr;
NNTI_buffer_t put_src_mr;
NNTI_buffer_t put_dst_mr;

NNTI_work_request_t *get_wr;  // need an array of work requests because we will be issuing multiple async gets
NNTI_work_request_t *put_wr;  // need an array of work requests because we will be issuing multiple async puts

NNTI_result_t        err;

int one_mb=1024*1024;

log_level nntiperf_debug_level = LOG_UNDEFINED;

uint32_t num_sends;

uint32_t num_gets;
uint32_t get_size;

uint32_t num_puts;
uint32_t put_size;


static int buffer_pack(void *input, char **output, uint64_t *output_size)
{
    NNTI_dt_sizeof(&trans_hdl, input, output_size);
    *output=(char*)malloc(*output_size);
    NNTI_dt_pack(&trans_hdl, input, *output, *output_size);

    return(0);
}

static int buffer_free(void *input)
{
    NNTI_dt_free(&trans_hdl, input);

    return(0);
}

static int buffer_unpack(char *input, uint64_t input_size, void *output)
{
    NNTI_dt_unpack(&trans_hdl, output, input, input_size);

    return(0);
}

void client(void) {
    NNTI_result_t rc=NNTI_OK;
    NNTI_status_t rdma_status;
    NNTI_status_t send_status;
    NNTI_status_t client_ack_status;
    char    *c_ptr;
    char     *packed=NULL;
    uint64_t  packed_size=0;

    double op_timer;

    Teuchos::oblackholestream blackhole;
    std::ostream &out = ( rank == 1 ? std::cout : blackhole );

    NNTI_connect(&trans_hdl, url, 5000, &server_hdl);

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    char *send_buf=NNTI_BUFFER_C_POINTER(&send_mr);
    memset(send_buf, 0, NNTI_REQUEST_BUFFER_SIZE);

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_RECV_DST, &client_ack_mr);
    char *client_ack_buf=NNTI_BUFFER_C_POINTER(&client_ack_mr);
    memset(client_ack_buf, 0, NNTI_REQUEST_BUFFER_SIZE);

    NNTI_alloc(&trans_hdl, get_size, 1, NNTI_GET_DST, &get_dst_mr);
    char *get_dst_buf=NNTI_BUFFER_C_POINTER(&get_dst_mr);
    memset(get_dst_buf, 0, get_size);

    NNTI_alloc(&trans_hdl, put_size, 1, NNTI_PUT_SRC, &put_src_mr);
    char *put_src_buf=NNTI_BUFFER_C_POINTER(&put_src_mr);
    memset(put_src_buf, 0, put_size);

    send_wr=(NNTI_work_request_t*)malloc(num_sends*sizeof(NNTI_work_request_t));
    get_wr =(NNTI_work_request_t*)malloc(num_gets*sizeof(NNTI_work_request_t));
    put_wr =(NNTI_work_request_t*)malloc(num_puts*sizeof(NNTI_work_request_t));

    /*
     * Phase 1 - exchange buffer handles
     */
    buffer_pack(&client_ack_mr, &packed, &packed_size);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(nntiperf_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    // send the server the recv_mr so it can send back it's ack_mr
    memcpy(send_buf, packed, packed_size);

    NNTI_create_work_request(&client_ack_mr, &client_ack_wr);

    rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr[0]);
    if (rc != NNTI_OK) {
        log_error(nntiperf_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr[0], 5000, &send_status);
    if (rc != NNTI_OK) {
        log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    // wait for the server to send back it's recv_mr
    rc=NNTI_wait(&client_ack_wr, -1, &client_ack_status);

    char *ptr=(char*)client_ack_status.start+client_ack_status.offset;

    memcpy(&packed_size, ptr, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(packed, ptr, packed_size);
    ptr += packed_size;

    buffer_unpack(packed, packed_size, &server_ack_mr);

    memcpy(&packed_size, ptr, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(packed, ptr, packed_size);
    ptr += packed_size;

    buffer_unpack(packed, packed_size, &get_src_mr);

    memcpy(&packed_size, ptr, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(packed, ptr, packed_size);
    ptr += packed_size;

    buffer_unpack(packed, packed_size, &put_dst_mr);

    NNTI_destroy_work_request(&client_ack_wr);

//    fprint_NNTI_buffer(logger_get_file(), "server_ack_mr",
//            "received server ack hdl", &server_ack_mr);
//    fprint_NNTI_buffer(logger_get_file(), "get_src_mr",
//            "received get src hdl", &get_src_mr);
//    fprint_NNTI_buffer(logger_get_file(), "put_dst_mr",
//            "received put dst hdl", &put_dst_mr);

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 2 - test sync request performance
     */
    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_sends;i++) {
        rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr[0]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_send() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&send_wr[0], 1000, &send_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_sends > 0) {
    	out << " sync requests per second == " << num_sends/op_timer << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 3 - test async request performance
     */
    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_sends;i++) {
        rc=NNTI_send(&server_hdl, &send_mr, NULL, &send_wr[i]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_send() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    for (uint32_t i=0;i<num_sends;i++) {
        rc=NNTI_wait(&send_wr[i], 1000, &send_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_sends > 0) {
    	out << "async requests per second == " << num_sends/op_timer << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 4 - test sync get performance
     */
    // warm up the pipes
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_get(&get_src_mr, client_rank*get_size, get_size, &get_dst_mr, 0, &get_wr[0]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_get() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&get_wr[0], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_get(&get_src_mr, client_rank*get_size, get_size, &get_dst_mr, 0, &get_wr[0]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_get() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&get_wr[0], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_gets > 0) {
    	out << " sync get (" << get_size << " byte transfer) == " << (double)(num_gets*get_size)/one_mb/op_timer << " MBps" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 5 - test async get performance
     */
    // warm up the pipes
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_get(&get_src_mr, client_rank*get_size, get_size, &get_dst_mr, 0, &get_wr[i]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_get() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_wait(&get_wr[i], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_get(&get_src_mr, client_rank*get_size, get_size, &get_dst_mr, 0, &get_wr[i]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_get() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    for (uint32_t i=0;i<num_gets;i++) {
        rc=NNTI_wait(&get_wr[i], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_gets > 0) {
    	out << "async get (" << get_size << " byte transfer) == " << (double)(num_gets*get_size)/one_mb/op_timer << " MBps" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 6 - test sync put performance
     */
    // warm up the pipes
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_put(&put_src_mr, 0, put_size, &put_dst_mr, client_rank*put_size, &put_wr[0]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_put() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&put_wr[0], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_put(&put_src_mr, 0, put_size, &put_dst_mr, client_rank*put_size, &put_wr[0]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_put() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
        rc=NNTI_wait(&put_wr[0], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_puts > 0) {
    	out << " sync put (" << put_size << " byte transfer) == " << (double)(num_puts*put_size)/one_mb/op_timer << " MBps" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 7 - test async put performance
     */
    // warm up the pipes
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_put(&put_src_mr, 0, put_size, &put_dst_mr, client_rank*put_size, &put_wr[i]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_put() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_wait(&put_wr[i], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }

    op_timer=trios_get_time();
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_put(&put_src_mr, 0, put_size, &put_dst_mr, client_rank*put_size, &put_wr[i]);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_put() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    for (uint32_t i=0;i<num_puts;i++) {
        rc=NNTI_wait(&put_wr[i], 1000, &rdma_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() did not return NNTI_OK: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }
    }
    op_timer=trios_get_time()-op_timer;
    if (num_puts > 0) {
    	out << "async put (" << put_size << " byte transfer) == " << (double)(num_puts*put_size)/one_mb/op_timer << " MBps" << std::endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    free(send_wr);
    free(get_wr);
    free(put_wr);

    buffer_free(&server_ack_mr);
    buffer_free(&get_src_mr);
    buffer_free(&put_dst_mr);

    NNTI_free(&send_mr);
    NNTI_free(&client_ack_mr);
    NNTI_free(&get_dst_mr);
    NNTI_free(&put_src_mr);

    return;
}

void server(void)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_status_t queue_status;
    NNTI_status_t send_status;
    NNTI_status_t server_ack_status;
    char *c_ptr;
    char     *packed=NULL;
    uint64_t  packed_size=0;


    int num_elements=nclients+(4*nclients*num_sends);
    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, num_elements, NNTI_RECV_QUEUE, &queue_mr);

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &send_mr);
    char *send_buf=NNTI_BUFFER_C_POINTER(&send_mr);
    memset(send_buf, 0, NNTI_REQUEST_BUFFER_SIZE);

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_RECV_DST, &server_ack_mr);
    char *server_ack_buf=NNTI_BUFFER_C_POINTER(&server_ack_mr);
    memset(server_ack_buf, 0, NNTI_REQUEST_BUFFER_SIZE);

    char *get_src_buf=(char *)malloc(nclients*get_size);
    memset(get_src_buf, 0, nclients*get_size);
    NNTI_register_memory(&trans_hdl, get_src_buf, nclients*get_size, 1, NNTI_GET_SRC, &get_src_mr);

    char *put_dst_buf=(char *)malloc(nclients*put_size);
    memset(put_dst_buf, 0, nclients*put_size);
    NNTI_register_memory(&trans_hdl, put_dst_buf, nclients*put_size, 1, NNTI_PUT_DST, &put_dst_mr);

    send_wr=(NNTI_work_request_t*)malloc(num_sends*sizeof(NNTI_work_request_t));

    /*
     * Phase 1 - exchange buffers handles
     */
    // wait for the client to send it's recv_mr
    NNTI_create_work_request(&queue_mr, &queue_wr);

    NNTI_wait(&queue_wr, -1, &queue_status);

    c_ptr=(char*)queue_status.start+queue_status.offset;
    buffer_unpack(c_ptr, queue_status.length, &client_ack_mr);

    NNTI_destroy_work_request(&queue_wr);

//    fprint_NNTI_buffer(logger_get_file(), "client_ack_mr",
//            "received client ack hdl", &client_ack_mr);

    // send our server_ack_mr, get_src_mr and put_dst_mr back to the client
    buffer_pack(&server_ack_mr, &packed, &packed_size);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(nntiperf_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    char *ptr=send_buf;
    memcpy(ptr, &packed_size, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(ptr, packed, packed_size);
    ptr += packed_size;

    free(packed);

    buffer_pack(&get_src_mr, &packed, &packed_size);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(nntiperf_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    memcpy(ptr, &packed_size, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(ptr, packed, packed_size);
    ptr += packed_size;

    free(packed);

    buffer_pack(&put_dst_mr, &packed, &packed_size);
    if (packed_size > NNTI_REQUEST_BUFFER_SIZE) {
        log_error(nntiperf_debug_level, "buffer_pack() says encoded NNTI_buffer_t is larger than NNTI_REQUEST_BUFFER_SIZE");
    	MPI_Abort(MPI_COMM_WORLD, -10);
    }

    memcpy(ptr, &packed_size, sizeof(packed_size));
    ptr += sizeof(packed_size);
    memcpy(ptr, packed, packed_size);
    ptr += packed_size;

    free(packed);

    rc=NNTI_send(&queue_status.src, &send_mr, &client_ack_mr, &send_wr[0]);
    if (rc != NNTI_OK) {
        log_error(nntiperf_debug_level, "NNTI_send() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }
    rc=NNTI_wait(&send_wr[0], 5000, &send_status);
    if (rc != NNTI_OK) {
        log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
        MPI_Abort(MPI_COMM_WORLD, rc);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 2 - client sends sync requests
     */
    for (uint32_t i=0;i<nclients*num_sends;i++) {
        NNTI_create_work_request(&queue_mr, &queue_wr);

        rc=NNTI_wait(&queue_wr, 1000, &queue_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }

        NNTI_destroy_work_request(&queue_wr);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 3 - client sends async requests
     */
    for (uint32_t i=0;i<nclients*num_sends;i++) {
        NNTI_create_work_request(&queue_mr, &queue_wr);

        rc=NNTI_wait(&queue_wr, 1000, &queue_status);
        if (rc != NNTI_OK) {
            log_error(nntiperf_debug_level, "NNTI_wait() returned an error: %d", rc);
            MPI_Abort(MPI_COMM_WORLD, rc);
        }

        NNTI_destroy_work_request(&queue_wr);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 4 - client does sync gets
     */
    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 5 - client does async gets
     */
    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 6 - client does sync puts
     */
    MPI_Barrier(MPI_COMM_WORLD);

    /*
     * Phase 7 - client does async puts
     */
    MPI_Barrier(MPI_COMM_WORLD);

    free(send_wr);

    buffer_free(&client_ack_mr);

    NNTI_free(&queue_mr);
    NNTI_free(&send_mr);
    NNTI_free(&server_ack_mr);

    NNTI_unregister_memory(&get_src_mr);
    free(get_src_buf);

    NNTI_unregister_memory(&put_dst_mr);
    free(put_dst_buf);

    return;
}

int parse_args(int argc, char *argv[])
{
	num_sends=atol(argv[1]);

	num_gets=atol(argv[2]);
	get_size=atol(argv[3]);

	num_puts=atol(argv[4]);
	put_size=atol(argv[5]);

	return(0);
}

int main(int argc, char *argv[])
{
    char logname[1024];

    bool success=true;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    nservers=1;
    nclients=nprocs-nservers;
    client_rank=rank-1;

    if (nprocs != 2) {
    	if (rank == 0) {
        	fprintf(stderr, "%s only supports 2 ranks.\n", argv[0]);
            MPI_Barrier(MPI_COMM_WORLD);
    	} else {
    	    MPI_Barrier(MPI_COMM_WORLD);
    	}
    	MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (argc == 6) {
        parse_args(argc, argv);
    } else if (argc == 1) {
        // no args from user.  set some defaults.
        num_sends=10;
        num_gets=10;
        get_size=1024*1024;
        num_puts=10;
        put_size=1024*1024;
    } else {
        // partial args from user.  can't do that.  abort.
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <num sends> <num_gets> <get_size> <num_puts> <put_size>\n", argv[0]);
            MPI_Barrier(MPI_COMM_WORLD);
        } else {
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    sprintf(logname, "nntiperf.%03d.log", rank);
    logger_init(LOG_ERROR, NULL);

    NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);

    if (rank==0) {
        NNTI_get_url(&trans_hdl, url, NNTI_URL_LEN);
    }

    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

    log_debug(nntiperf_debug_level, "NNTI perfermance server url is %s", url);

    if (rank==0) {
        server();
    } else {
        client();
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
//    return 0;
}
