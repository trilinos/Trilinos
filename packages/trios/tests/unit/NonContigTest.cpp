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


log_level noncontig_debug_level = LOG_UNDEFINED;


NNTI_transport_t     trans_hdl;
NNTI_peer_t          server_hdl;

NNTI_buffer_t        contig_queue_mr;
NNTI_work_request_t  contig_queue_wr;
NNTI_status_t        queue_status;

NNTI_buffer_t        contig_send_mr, contig_recv_mr, contig_get_src_mr, contig_put_dst_mr; /* contiguous registered memory regions */
NNTI_work_request_t  contig_send_wr, contig_recv_wr, contig_get_dst_wr, contig_put_src_wr; /* work requests */

NNTI_buffer_t        noncontig_get_dst_mr, noncontig_put_src_mr; /* non-contiguous registered memory regions */
NNTI_work_request_t  noncontig_get_dst_wr, noncontig_put_src_wr; /* work requests */

#define WR_COUNT 2
NNTI_work_request_t *mr_wr_list[WR_COUNT];
NNTI_result_t        err;
NNTI_status_t wait_status1;
NNTI_status_t wait_status2;
NNTI_status_t wait_status3;
NNTI_status_t *wait_status_list[WR_COUNT];

int one_kb=1024;
int one_mb=1024*1024;

#define NUM_SEGMENTS 5

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


int client(void) {
	int success=TRUE;

	uint64_t *server_checksums=NULL;

	char *buf=NULL;

    int get_mr_size;
    int put_mr_size;
    int recv_mr_size;

    XDR send_xdrs;

    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 1, NNTI_SEND_SRC, &contig_send_mr);

    NNTI_alloc(&trans_hdl, NUM_SEGMENTS*one_kb, 1, NNTI_GET_SRC,  &contig_get_src_mr);
    NNTI_alloc(&trans_hdl, NUM_SEGMENTS*one_kb, 1, NNTI_PUT_DST,  &contig_put_dst_mr);

    NNTI_alloc(&trans_hdl, NNTI_RESULT_BUFFER_SIZE, 1, NNTI_RECV_DST, &contig_recv_mr);

    buf=NNTI_BUFFER_C_POINTER(&contig_get_src_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_kb), 'M'+i, one_kb);
    }
    buf=NNTI_BUFFER_C_POINTER(&contig_put_dst_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        memset(buf + (i*one_kb), 'S'+i, one_kb);
    }

    /* XDR encode the get, put and recv buffers */
    get_mr_size  = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, &contig_get_src_mr);
    put_mr_size  = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, &contig_put_dst_mr);
    recv_mr_size = xdr_sizeof((xdrproc_t)&xdr_NNTI_buffer_t, &contig_recv_mr);

    log_debug(noncontig_debug_level, "get_mr_size=%d ; put_mr_size=%d ; recv_mr_size=%d", get_mr_size, put_mr_size, recv_mr_size);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&contig_send_mr),
            NNTI_BUFFER_SIZE(&contig_send_mr), XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, &contig_get_src_mr);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&contig_send_mr)+get_mr_size,
            NNTI_BUFFER_SIZE(&contig_send_mr)-get_mr_size, XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, &contig_put_dst_mr);

    xdrmem_create(&send_xdrs, NNTI_BUFFER_C_POINTER(&contig_send_mr)+get_mr_size+put_mr_size,
            NNTI_BUFFER_SIZE(&contig_send_mr)-get_mr_size-put_mr_size, XDR_ENCODE);
    xdr_NNTI_buffer_t(&send_xdrs, &contig_recv_mr);

    NNTI_send(&server_hdl, &contig_send_mr, NULL, &contig_send_wr);
    NNTI_wait(&contig_send_wr, 1000, &wait_status1);

    NNTI_create_work_request(&contig_recv_mr, &contig_recv_wr);
    NNTI_wait(&contig_recv_wr, 1000, &wait_status1);

    server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&contig_recv_mr);

	/* checksums GET */
	buf=NNTI_BUFFER_C_POINTER(&contig_get_src_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
    	for (int j=0;j<5;j++) {
    		out << (buf + (i*one_kb))[j];
    	}
    	log_debug(noncontig_debug_level, "client GET src segments[%d]=%s", i, out.str().c_str());

    	uint64_t my_checksum=calc_checksum(buf + (i*one_kb), one_kb);
    	log_debug(noncontig_debug_level, "server_checksums[%d]=%X ; my_checksum=%X", i, server_checksums[i], my_checksum);
    	if (server_checksums[i] != my_checksum) {
    		success=FALSE;
    	}
    }

	/* checksums PUT */
    buf=NNTI_BUFFER_C_POINTER(&contig_put_dst_mr);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
    	for (int j=0;j<5;j++) {
    		out << (buf + (i*one_kb))[j];
    	}
    	log_debug(noncontig_debug_level, "client PUT dst segments[%d]=%s", i, out.str().c_str());

    	uint64_t my_checksum=calc_checksum(buf + (i*one_kb), one_kb);
    	log_debug(noncontig_debug_level, "server_checksums[NUM_SEGMENTS+%d]=%X ; my_checksum=%X", i, server_checksums[NUM_SEGMENTS+i], my_checksum);
    	if (server_checksums[NUM_SEGMENTS+i] != my_checksum) {
    		success=FALSE;
    	}
    }

    NNTI_free(&contig_send_mr);
    NNTI_free(&contig_get_src_mr);
    NNTI_free(&contig_put_dst_mr);
    NNTI_free(&contig_recv_mr);

    return(success);
}

void server(void) {

    char    *get_buf,*put_buf;
    char    *get_segments[NUM_SEGMENTS],*put_segments[NUM_SEGMENTS];
    uint64_t get_segment_lengths[NUM_SEGMENTS],put_segment_lengths[NUM_SEGMENTS];

    XDR recv_xdrs;


    NNTI_alloc(&trans_hdl, NNTI_REQUEST_BUFFER_SIZE, 10, NNTI_RECV_QUEUE, &contig_queue_mr);
    NNTI_alloc(&trans_hdl, NNTI_RESULT_BUFFER_SIZE, 1, NNTI_SEND_SRC, &contig_send_mr);

    get_buf=(char *)malloc(NUM_SEGMENTS*one_mb);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        get_segment_lengths[i]=one_kb;
        get_segments[i]=get_buf + (i*one_mb);

        memset(get_segments[i], 'A'+i, one_mb);
    }
    NNTI_register_segments(&trans_hdl, get_segments, get_segment_lengths, NUM_SEGMENTS, NNTI_GET_DST, &noncontig_get_dst_mr);

    put_buf=(char *)malloc(NUM_SEGMENTS*one_mb);
    for (int i=0;i<NUM_SEGMENTS;i++) {
        put_segment_lengths[i]=one_kb;
        put_segments[i]=put_buf + (i*one_mb);

        memset(put_segments[i], 'G'+i, one_mb);
    }
    NNTI_register_segments(&trans_hdl, put_segments, put_segment_lengths, NUM_SEGMENTS, NNTI_PUT_SRC, &noncontig_put_src_mr);


    NNTI_create_work_request(&contig_queue_mr, &contig_queue_wr);
    NNTI_wait(&contig_queue_wr, -1, &queue_status);


    /* XDR decode the get, put and recv buffers */
    xdrmem_create(&recv_xdrs, (char *)queue_status.start,
            NNTI_REQUEST_BUFFER_SIZE, XDR_DECODE);

    xdr_NNTI_buffer_t(&recv_xdrs, &contig_get_src_mr);
    xdr_NNTI_buffer_t(&recv_xdrs, &contig_put_dst_mr);
    xdr_NNTI_buffer_t(&recv_xdrs, &contig_recv_mr);

    if (logging_debug(noncontig_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "contig_get_src_mr",
                "after XDR decode", &contig_get_src_mr);
        fprint_NNTI_buffer(logger_get_file(), "contig_put_dst_mr",
                "after XDR decode", &contig_put_dst_mr);
        fprint_NNTI_buffer(logger_get_file(), "contig_recv_mr",
                "after XDR decode", &contig_recv_mr);

        fprint_NNTI_buffer(logger_get_file(), "noncontig_get_dst_mr",
                "before GET", &noncontig_get_dst_mr);
        fprint_NNTI_buffer(logger_get_file(), "noncontig_put_src_mr",
                "before PUT", &noncontig_put_src_mr);
    }

    NNTI_get(&contig_get_src_mr, 0, NUM_SEGMENTS*one_kb, &noncontig_get_dst_mr, 0, &noncontig_get_dst_wr);
    NNTI_wait(&noncontig_get_dst_wr, 1000, &wait_status1);

    NNTI_put(&noncontig_put_src_mr, 0, NUM_SEGMENTS*one_kb, &contig_put_dst_mr, 0, &noncontig_put_src_wr);
    NNTI_wait(&noncontig_put_src_wr, 1000, &wait_status1);

	uint64_t *server_checksums=(uint64_t*)NNTI_BUFFER_C_POINTER(&contig_send_mr);

	/* checksum GET */
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
    	for (int j=0;j<5;j++) {
    		out << get_segments[i][j];
    	}
    	log_debug(noncontig_debug_level, "server GET dst segments[%d]=%s", i, out.str().c_str());

    	server_checksums[i]=calc_checksum(get_segments[i], get_segment_lengths[i]);
    	log_debug(noncontig_debug_level, "server GET checksums[%d]=%X", i, server_checksums[i]);
    }
	/* checksum PUT */
    for (int i=0;i<NUM_SEGMENTS;i++) {
        std::stringstream out(std::stringstream::out);
    	for (int j=0;j<5;j++) {
    		out << put_segments[i][j];
    	}
    	log_debug(noncontig_debug_level, "server PUT src segments[%d]=%s", i, out.str().c_str());

    	server_checksums[NUM_SEGMENTS+i]=calc_checksum(put_segments[i], put_segment_lengths[i]);
    	log_debug(noncontig_debug_level, "server PUT checksums[NUM_SEGMENTS+%d]=%X", i, server_checksums[NUM_SEGMENTS+i]);
    }

    NNTI_send(&queue_status.src, &contig_send_mr, &contig_recv_mr, &contig_send_wr);
    NNTI_wait(&contig_send_wr, 1000, &wait_status1);


    NNTI_unregister_memory(&noncontig_get_dst_mr);
    free(get_buf);
    NNTI_unregister_memory(&noncontig_put_src_mr);
    free(put_buf);

    NNTI_free(&contig_send_mr);
    NNTI_free(&contig_queue_mr);

    return;
}

int main(int argc, char *argv[])
{
	int success=0;
    int nprocs, rank;

    char logname[1024];
    char url[NNTI_URL_LEN];

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    sprintf(logname, "noncontig.%03d.log", rank);
    logger_init(LOG_ERROR, NULL /*logname*/);

    NNTI_init(NNTI_DEFAULT_TRANSPORT, NULL, &trans_hdl);

    if (rank==0) {
        NNTI_get_url(&trans_hdl, url, NNTI_URL_LEN);
    }

    MPI_Bcast(&url[0], NNTI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

    log_debug(noncontig_debug_level, "noncontig server url is %s", url);

    if (rank==1) {
        NNTI_connect(&trans_hdl, url, 5000, &server_hdl);
    }

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
