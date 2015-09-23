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
/*-------------------------------------------------------------------------*/
/**  @file rpc_client.cc
 *
 *   @brief  Implementation of the \ref rpc_client_api "RPC Client API".
 *           for the NSSI.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 1640 $
 *   $Date: 2007-11-28 11:59:53 -0700 (Wed, 28 Nov 2007) $
 *
 */
#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nssi_request.h"
#include "Trios_nssi_types.h"
#include "Trios_nssi_fprint_types.h"
#include "Trios_nnti_fprint_types.h"
#include "Trios_logger.h"
#include "Trios_threads.h"
#include "Trios_timer.h"
#include "buffer_queue.h"
#include "Trios_nssi_rpc.h"
#include "Trios_nssi_xdr.h"
#include "Trios_nssi_debug.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>

#include <vector>

#ifdef HAVE_TRIOS_MALLOC_H
#include <malloc.h>
#endif

#include <inttypes.h>
#include <unistd.h>
#include <stdlib.h>
#include <mpi.h>


#include "nssi_opcodes.h"
#include "nssi_service_args.h"


#include "Trios_nnti.h"

#define MIN_TIMEOUT 1000

#define MAX_RETRIES 50

extern NNTI_transport_t transports[NSSI_RPC_COUNT];
extern nssi_config_t nssi_config;

static bool client_initialized = false;

static nthread_counter_t request_count;

extern trios_buffer_queue_t send_bq;
extern trios_buffer_queue_t recv_bq;
extern trios_buffer_queue_t rdma_target_bq;
extern trios_buffer_queue_t rdma_get_bq;
extern trios_buffer_queue_t rdma_put_bq;

/**
 *   @addtogroup rpc_ptl_impl
 *
 *   This section describes a portals implementation of the
 *   NSSI asynchrounous \ref rpc_api "RPC" mechanism.
 */

/*--------------------------- Private methods ------------------- */

#ifdef ANSI_FUNC

static int client_init(void)
#else

static int
client_init ()
#endif
{
    int rc = NSSI_OK;

    if (client_initialized) {
        return rc;
    }

    // Initialize the thread counter
    nthread_counter_init(&request_count);

    NSSI_REGISTER_CLIENT_STUB(NSSI_OP_GET_SERVICE, void, void, nssi_service);
    NSSI_REGISTER_CLIENT_STUB(NSSI_OP_KILL_SERVICE, nssi_kill_service_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(NSSI_OP_TRACE_RESET, nssi_trace_reset_args, void, void);

    client_initialized = true;

    return rc;
}

/**
 * @brief Clean up portals data structures associated with
 * operation requests.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method frees memory and data structures created in
 * the prepare_request and prepare_result methods.
 *
 *
 * If the arguments were too large for the short request, the
 * process_request method allocated memory for the encoded
 * arguments and all portals data structures for the arguments.
 *
 * For results, the process_request method allocates space for
 * the encoded short result and portals data structures associated
 * with the short result.
 */
static int cleanup(nssi_request *request)
{
    return NSSI_OK;
}


/**
 * @brief Process the result of an operation request.
 *
 * @ingroup rpc_ptl_impl
 *
 * After processing an operation request, a NSSI server
 * sends a short result back to the client. If the
 * contents of the result do not fit into the short result
 * the server stores the result until the client can fetch
 * it.
 *
 * This client-side method executes after the short
 * result arrives at the client. It decodes the
 * result header, possibly fetches the results from
 * the server, and decodes the result.
 *
 * @param encoded_short_res @input The result buffer received from the server.
 * @param request           @input The request data structure.
 */
static int process_result(char *encoded_short_res_buf, nssi_request *request)
{

    int rc = NSSI_OK;  /* return code */
    uint32_t result_size;
    nssi_result_header header;
    void *decoded_result = NULL;
    XDR hdr_xdrs;
    XDR res_xdrs;

    trios_declare_timer(call_time);

    NNTI_buffer_t       encoded_long_res_hdl;
    NNTI_work_request_t encoded_long_res_wr;
    NNTI_status_t       status;

    NNTI_buffer_t       long_res_ack_hdl;
    NNTI_work_request_t long_res_ack_wr;
    NNTI_status_t       long_res_ack_status;

    log_level debug_level = rpc_debug_level;

    /* assign the result decoding function */
    xdrproc_t xdr_decode_result = request->xdr_decode_result;

    log_debug(debug_level,"Start processing result, status=%d", request->status);

    /** set the status to processing */
    request->status = NSSI_PROCESSING_RESULT;

    /* initialize the header */
    memset(&header, 0, sizeof(nssi_result_header));

    xdrmem_create(&hdr_xdrs,
            encoded_short_res_buf,
            NSSI_SHORT_RESULT_SIZE,
            XDR_DECODE);

    /* decode the header */
    log_debug(debug_level,"decoding result header...");

    trios_start_timer(call_time);
    if (! xdr_nssi_result_header(&hdr_xdrs, &header)) {
        log_fatal(rpc_debug_level,"failed to decode the result header");
        rc = NSSI_EDECODE;
        goto cleanup;
    }
    trios_stop_timer("xdr_nssi_result_header - decode", call_time);

    if (logging_debug(debug_level)) {
        fprint_nssi_result_header(logger_get_file(), "header", "DEBUG", &header);
    }

    /* what to do if the remote code had an error */
    if (header.rc != NSSI_OK) {
        request->status = NSSI_REQUEST_ERROR;
        request->error_code = header.rc;
        rc = NSSI_OK;
        goto cleanup;
    }

    /* get result size from the header */
    result_size = header.result_size;

    if (result_size > 0) {

        /* decode the result */
        log_debug(debug_level,"getting the result (size == %d)...", result_size);
        decoded_result = request->result;

        /* extract result from header */
        if (!header.fetch_result) {
            log_debug(debug_level,"decoding the result...");
            trios_start_timer(call_time);
            if (!xdr_decode_result(&hdr_xdrs, decoded_result))   {
                log_fatal(debug_level,"failed to decode the result");
                rc = NSSI_EDECODE;
                goto cleanup;
            }
            trios_stop_timer("xdr_decode_result - decode", call_time);
        }

        /* else fetch the result from the server */
        else {
            log_debug(debug_level,"fetching result (%d bytes) "
                    "from server...", result_size);

            /* allocate space for result */
            trios_start_timer(call_time);
            rc=NNTI_alloc(
                    &transports[request->svc->transport_id],
                    result_size,
                    1,
                    NNTI_GET_DST,
                    &encoded_long_res_hdl);
            trios_stop_timer("NNTI_register_memory - long result", call_time);
            if (rc != NNTI_OK) {
                log_error(debug_level, "failed registering long result: %s",
                        nnti_err_str(rc));
            }

            /* fetch the data from the server */
            trios_start_timer(call_time);
            rc=NNTI_get(&header.result_addr,
                    0,
                    result_size,
                    &encoded_long_res_hdl,
                    0,
                    &encoded_long_res_wr);
            trios_stop_timer("NNTI_get - long result", call_time);
            if (rc != NNTI_OK) {
                log_error(debug_level, "failed getting long result: %s",
                        nnti_err_str(rc));
            }
            trios_start_timer(call_time);
            rc=NNTI_wait(&encoded_long_res_wr,
                    -1,
                    &status);
            trios_stop_timer("NNTI_wait - long result", call_time);
            if (rc != NNTI_OK) {
                log_error(debug_level, "failed waiting for the long result: %s",
                        nnti_err_str(rc));
            }
            /* create a memory stream for XDR-decoding the result */
            xdrmem_create(&res_xdrs, NNTI_BUFFER_C_POINTER(&encoded_long_res_hdl),
                    NNTI_BUFFER_SIZE(&encoded_long_res_hdl), XDR_DECODE);

            log_debug(debug_level,"decoding the fetched result...");
            trios_start_timer(call_time);
            if (!xdr_decode_result(&res_xdrs, decoded_result))   {
                log_fatal(debug_level,"failed to decode the result");
                rc = NSSI_EDECODE;
                goto cleanup;
            }
            trios_stop_timer("xdr_decode_result - decode", call_time);

            trios_start_timer(call_time);
            rc=NNTI_alloc(
                    &transports[request->svc->transport_id],
                    sizeof(int8_t),
                    1,
                    NNTI_SEND_SRC,
                    &long_res_ack_hdl);
            trios_stop_timer("NNTI_register_memory - long result ack", call_time);
            if (rc != NNTI_OK) {
                log_error(debug_level, "failed registering long result ack: %s",
                        nnti_err_str(rc));
            }

            trios_start_timer(call_time);
            rc=NNTI_send(
                    &request->svc->svc_host,
                    &long_res_ack_hdl,
                    &header.result_ack_addr,
                    &long_res_ack_wr);
            trios_stop_timer("NNTI_send - long result ack", call_time);
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed sending short result: %s",
                        nnti_err_str(rc));
            }
            trios_start_timer(call_time);
            rc=NNTI_wait(
                    &long_res_ack_wr,
                    -1,
                    &long_res_ack_status);
            trios_stop_timer("NNTI_wait - long result ack", call_time);
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed waiting for short result: %s",
                        nnti_err_str(rc));
            }
        }
    }

    request->status = NSSI_REQUEST_COMPLETE;

cleanup:
    if (header.fetch_result) {
        NNTI_free(&encoded_long_res_hdl);
        if (rc != NNTI_OK) {
            log_error(debug_level, "failed unregistering long result: %s",
                    nnti_err_str(rc));
        }

        NNTI_free(&long_res_ack_hdl);
        if (rc != NNTI_OK) {
            log_error(debug_level, "failed unregistering long result ack: %s",
                    nnti_err_str(rc));
        }
    }

    /* clean up portals data structures */
    log_debug(debug_level,"clean up data structures");
    cleanup(request);

    log_debug(debug_level,"end");

    return rc;
}


/**
 * @brief Initialize an RPC client and return the
 * service provided by the server. */
int nssi_rpc_clnt_init(
        const nssi_rpc_transport rpc_transport,
        const char              *my_url,
        const char              *svc_url,
        nssi_service            *result)
{
    int rc = NSSI_OK;
    int default_timeout = 1000;

    /* initialize RPC */
    rc = nssi_rpc_init(rpc_transport, NSSI_DEFAULT_ENCODE, my_url);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "could not initialize RPC");
        return rc;
    }

    nssi_init(rpc_transport);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "could not initialize nssi client");
        return rc;
    }

    /* ping the server */
    return nssi_get_service(rpc_transport, svc_url, default_timeout, result);
}

/**
 * @brief Send a ping request to the server.
 */
int nssi_get_services(
        const nssi_rpc_transport rpc_transport,
        const char             **url,
        const int                num_servers,
        const int                timeout,
        nssi_service            *result)
{
    int rc = NSSI_OK;
    int i;

    for (i=0; i<num_servers; i++) {
        rc = nssi_get_service(rpc_transport, url[i], timeout, &result[i]);
        if (rc != NSSI_OK) {
            log_error(rpc_debug_level, "could not get service %d",i);
            return rc;
        }
    }

    return rc;
}

/**
 * @brief Send a ping request to the server.
 */
int nssi_get_service(
        const nssi_rpc_transport rpc_transport,
        const char              *url,
        const int                timeout,
        nssi_service            *result)
{
    int rc = NSSI_OK;
    int cleanup_rc = NSSI_OK;

    /* local count does not need protection */
    static int64_t local_count;

    /* xdrs for the header and the args. */
    XDR hdr_xdrs, res_xdrs;

    NNTI_peer_t   peer_hdl;

    nssi_request_header req_header;   /* the request header */
    nssi_result_header  res_header;   /* the request header */
    char *buf;
    int short_req_len = 0;
    NNTI_buffer_t       short_req;
    NNTI_buffer_t       short_res;
    NNTI_buffer_t      *short_req_hdl=&short_req;
    NNTI_work_request_t short_req_wr;
    NNTI_buffer_t      *short_res_hdl=&short_res;
    NNTI_work_request_t short_res_wr;
    NNTI_status_t       wait_status;

    unsigned long len=0;

    log_debug(rpc_debug_level, "entered nssi_get_service");

    client_init();

    rc=NNTI_connect(
            &transports[rpc_transport],
            url,
            timeout,
            &peer_hdl);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed to connect to service at %s: %s",
                url, nnti_err_str(rc));
        return(rc);
    }


    /*------ Initialize variables and buffers ------*/
    memset(&req_header, 0, sizeof(nssi_request_header));

    if (nssi_config.use_buffer_queue) {
        short_req_hdl=trios_buffer_queue_pop(&send_bq);
        assert(short_req_hdl);

        short_res_hdl=trios_buffer_queue_pop(&recv_bq);
        assert(short_res_hdl);
    } else {
        /* allocate memory for the short request buffer */
        short_req_len = NSSI_SHORT_REQUEST_SIZE;
        rc=NNTI_alloc(
                &transports[rpc_transport],
                short_req_len,
                1,
                NNTI_SEND_SRC,
                short_req_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering short request: %s",
                    nnti_err_str(rc));
        }
        buf=NNTI_BUFFER_C_POINTER(short_req_hdl);
        memset(buf, 0, short_req_len);  // fixes uninitialized bytes error from valgrind

        rc=NNTI_alloc(
                &transports[rpc_transport],
                NSSI_SHORT_RESULT_SIZE,
                1,
                NNTI_RECV_DST,
                short_res_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering short result: %s",
                    nnti_err_str(rc));
        }
        buf=NNTI_BUFFER_C_POINTER(short_res_hdl);
        memset(buf, 0, NSSI_SHORT_RESULT_SIZE);
    }

    req_header.res_addr=*short_res_hdl;

    /* increment global counter */
    local_count = nthread_counter_increment(&request_count);
    if (local_count == -1) {
        log_error(rpc_debug_level, "Unable to increment counter");
        return NNTI_EIO;
    }

    /* set the opcode for the request header */
    req_header.opcode = NSSI_OP_GET_SERVICE;
    /* set the request ID for the header */
    req_header.id = local_count;
    req_header.fetch_args = FALSE;

    /* create an xdr memory stream for the short request buffer */
    xdrmem_create(&hdr_xdrs, NNTI_BUFFER_C_POINTER(short_req_hdl),
            NNTI_BUFFER_SIZE(short_req_hdl), XDR_ENCODE);
    /* encode the header  */
    log_debug(rpc_debug_level,"encoding request header");
    if (! xdr_nssi_request_header(&hdr_xdrs, &req_header)) {
        log_fatal(rpc_debug_level,"failed to encode the request header");
        return NSSI_EENCODE;
    }
    /* get the number of valid bytes in the request */
    len = xdr_sizeof((xdrproc_t)&xdr_nssi_request_header, &req_header);

    /* send the encoded short request buffer to the server */
    log_debug(rpc_debug_level,"sending short request, id=%lu, len=%d", req_header.id, len);

    /* Times out after DEFAULT_RPC_TIMEOUT */
    rc=NNTI_send(
            &peer_hdl,
            short_req_hdl,
            NULL,
            &short_req_wr);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed sending get service request: %s",
                nnti_err_str(rc));
    }
    rc=NNTI_wait(
            &short_req_wr,
            timeout,
            &wait_status);
    if (rc == NNTI_OK) {
        NNTI_create_work_request(
                short_res_hdl,
                &short_res_wr);
        rc=NNTI_wait(
                &short_res_wr,
                timeout,
                &wait_status);
        NNTI_destroy_work_request(
                &short_res_wr);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed waiting for short result: %s",
                    nnti_err_str(rc));
        }
    }
    else if (rc == NNTI_ETIMEDOUT) {
        log_info(rpc_debug_level, "put request timed out");
    }
    else  {
        log_error(rpc_debug_level, "failed waiting for get_service: %s",
                nnti_err_str(rc));
    }


    if (rc == NNTI_OK) {
        log_debug(rpc_debug_level,"message sent");

        xdrmem_create(
                &res_xdrs,
                NNTI_BUFFER_C_POINTER(short_res_hdl),
                NNTI_BUFFER_SIZE(short_res_hdl),
                XDR_DECODE);
        log_debug(rpc_debug_level,"decoding result header");
        memset(&res_header, 0, sizeof(nssi_result_header));
        if (! xdr_nssi_result_header(&res_xdrs, &res_header)) {
            log_fatal(rpc_debug_level,"failed to decode the result header");
            return NSSI_EDECODE;
        }
        memset(result, 0, sizeof(nssi_service));
        if (! xdr_nssi_service(&res_xdrs, result)) {
            log_fatal(rpc_debug_level,"failed to decode the nssi_service");
            return NSSI_EDECODE;
        }

        if (logging_debug(rpc_debug_level)) {
            fprint_nssi_service(logger_get_file(), "result", "end of nssi_get_service", result);
        }
    }

// Cleanup data structures

    if (nssi_config.use_buffer_queue) {
        trios_buffer_queue_push(&send_bq, short_req_hdl);
        short_req_hdl=NULL;

        trios_buffer_queue_push(&recv_bq, short_res_hdl);
        short_res_hdl=NULL;
    } else {
        cleanup_rc=NNTI_free(short_req_hdl);
        if (cleanup_rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering short request: %s",
                    nnti_err_str(cleanup_rc));
        }

        cleanup_rc=NNTI_free(short_res_hdl);
        if (cleanup_rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering short result: %s",
                    nnti_err_str(cleanup_rc));
        }
    }

    log_debug(rpc_debug_level, "finished nssi_get_service (rc=%d)", rc);

    return rc;
}

/**
 * @brief Send a ping request to the server.
 */
int nssi_free_service(
        const nssi_rpc_transport rpc_transport,
        nssi_service            *svc)
{
    int rc=0;

    rc=NNTI_disconnect(
            &transports[rpc_transport],
            &svc->svc_host);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed disconnecting from service: %s",
                nnti_err_str(rc));
    }
    return(rc);
}

/* --------------- CLIENT INTERFACE --------- */

/**
 * @brief Initialize the client.
 */
int nssi_init(const nssi_rpc_transport transport_id)
{
    int rc = NSSI_OK;

    client_init();

    client_initialized=true;

    return rc;
}

/**
 * @brief Complete a Nessie client.
 */
int nssi_fini(const nssi_rpc_transport transport_id)
{
    int rc = NSSI_OK;

    nthread_counter_fini(&request_count);

    client_initialized=false;

    return(rc);
}

/**
 * @brief Send a kill request to the server.
 */
int nssi_kill(
        const nssi_service *svc,
        const int sig,
        const long timeout)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;
    log_level debug_level = rpc_debug_level;
    debug_level = rpc_debug_level;

    nssi_kill_service_args args;

    memset(&args, 0, sizeof(nssi_kill_service_args));
    args.sig = sig;

    client_init();

    /* send the request */
    rc = nssi_call_rpc(svc, NSSI_OP_KILL_SERVICE, &args, NULL, 0, NULL, &req);
    if (rc != NSSI_OK) {
        log_error(debug_level, "unable to call remote method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* If we're forcing a kill, we don't want to wait for the response
     * from the server.  This will likely cause some memory leaks. */
    if (sig != 0) {
        req.status = NSSI_REQUEST_COMPLETE;
    }


    /* wait for completion */
    log_debug(debug_level, "calling nssi_timedwait");
    //rc2 = nssi_timedwait(&req, timeout, &rc);
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(rpc_debug_level, "failed waiting for request %lu: %s",
                req.id,
                nssi_err_str(rc2));

        /* In this case, we failed waiting for the request because
           the service exited before it had a chance to send
           the result. */
        /*return rc2;*/
    }



    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nssi_err_str(rc));

        /* this case is a legitimate error */
        return rc;
    }

    log_debug(debug_level, "Finished with nssi_kill");

    return rc;
}


/**
 * @brief Reset tracing on a server.
 */
int nssi_trace_reset(
        const nssi_service *svc,
        const char *fname,
        const int ftype,
        const char *enable,
        const long timeout)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;
    static const char *empty="";

    nssi_trace_reset_args args;

    client_init();

    memset(&args, 0, sizeof(nssi_trace_reset_args));

    /* initialize the args */
    args.fname = (char *)fname;
    args.ftype = ftype;
    args.enable = (char *)enable;

    log_debug(rpc_debug_level, "Calling nssi_trace_reset(%s, %d, %s)\n",
            fname, ftype, enable);

    if (!args.fname)
        args.fname = (char *)empty;

    if (!args.enable)
        args.enable = (char *)empty;

    /* send the request */
    rc = nssi_call_rpc(svc, NSSI_OP_TRACE_RESET, &args, NULL, 0, NULL, &req);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "unable to call remote method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    log_debug(rpc_debug_level, "calling nssi_timedwait");
    rc2 = nssi_timedwait(&req, timeout, &rc);
    if (rc2 != NSSI_OK) {
        log_error(rpc_debug_level, "failed waiting for request %lu: %s",
                req.id,
                nssi_err_str(rc2));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;
}

/** @brief Test a request for completion
 *
 * @param req the handle of the pending request.
 * @return status of the operation.
 */
int nssi_test(nssi_request *req, int *rc) {
    log_error(rpc_debug_level,"not supported");
    return (int)NSSI_ENOTSUP;
}

static int setup_short_request(nssi_request *request)
{
    int rc = NSSI_OK;
    int short_request_len = 0;

    trios_declare_timer(call_time);

    if (nssi_config.use_buffer_queue) {
        request->short_request_hdl=trios_buffer_queue_pop(&send_bq);
        log_debug(rpc_debug_level, "Popped short request buffer");
    } else {
        /* allocate memory for the short request buffer */
        short_request_len = NNTI_BUFFER_SIZE(&request->svc->req_addr);

        trios_start_timer(call_time);
        rc=NNTI_alloc(
                &transports[request->svc->transport_id],
                short_request_len,
                1,
                NNTI_SEND_SRC,
                request->short_request_hdl);
        trios_stop_timer("NNTI_register_memory - short request", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering short request: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
        log_debug(rpc_debug_level, "Allocated short request buffer");
    }

cleanup:
    return rc;
}

static int teardown_short_request(nssi_request *request)
{
    int rc = NSSI_OK;

    trios_declare_timer(call_time);

    if (nssi_config.use_buffer_queue) {
        trios_buffer_queue_push(&send_bq, request->short_request_hdl);
    } else {
        trios_start_timer(call_time);
        rc=NNTI_free(request->short_request_hdl);
        trios_stop_timer("NNTI_unregister_memory - send req", call_time);
    }
    request->short_request_hdl=NULL;

    return rc;
}

static int setup_bulk_data(nssi_request *request)
{
    int rc = NSSI_OK;

    trios_declare_timer(call_time);

    if ((nssi_config.use_buffer_queue) &&
            (nssi_config.rdma_buffer_queue_buffer_size >= request->data_size)) {
        log_debug(rpc_debug_level, "using buffer queue for TARGET buffer");
        request->bulk_data_hdl=trios_buffer_queue_pop(&rdma_target_bq);
        assert(request->bulk_data_hdl);
        NNTI_BUFFER_SIZE(request->bulk_data_hdl)=request->data_size;
        /* copy the user buffer contents into RDMA buffer */
        trios_start_timer(call_time);
        memcpy(NNTI_BUFFER_C_POINTER(request->bulk_data_hdl), (char *)request->data, request->data_size);
        trios_stop_timer("memcpy target buf to bq", call_time);
    } else {
        log_debug(rpc_debug_level, "using user buffer for TARGET buffer");
        log_debug (rpc_debug_level, "Registering data buffer (size=%d)", request->data_size);

        trios_start_timer(call_time);
        rc=NNTI_register_memory(
                &transports[request->svc->transport_id],
                (char *)request->data,
                request->data_size,
                1,
                (NNTI_buf_ops_t)(NNTI_GET_SRC|NNTI_PUT_DST),
                request->bulk_data_hdl);
        trios_stop_timer("NNTI_register_memory - data", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering data: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
        log_debug(rpc_debug_level, "Allocated bulk data buffer");
    }

    if (logging_debug(rpc_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "request->bulk_data_hdl",
                "NNTI_buffer_t", request->bulk_data_hdl);
    }

cleanup:
    return rc;
}

static int teardown_bulk_data(nssi_request *request)
{
    int rc = NSSI_OK;

    if ((nssi_config.use_buffer_queue) &&
            (nssi_config.rdma_buffer_queue_buffer_size >= request->data_size)) {
        trios_buffer_queue_push(&rdma_target_bq, request->bulk_data_hdl);
    } else {
        rc=NNTI_unregister_memory(request->bulk_data_hdl);
    }
    request->bulk_data_hdl=NULL;

    return rc;
}

static int setup_short_result(nssi_request *request)
{
    int rc = NSSI_OK;

    trios_declare_timer(call_time);

    if (nssi_config.use_buffer_queue) {
        request->short_result_hdl=trios_buffer_queue_pop(&recv_bq);
    } else {
        log_debug (rpc_debug_level, "allocating short result (size=%d)", NSSI_SHORT_RESULT_SIZE);

        trios_start_timer(call_time);
        rc=NNTI_alloc(
                &transports[request->svc->transport_id],
                NSSI_SHORT_RESULT_SIZE,
                1,
                NNTI_RECV_DST,
                request->short_result_hdl);
        trios_stop_timer("NNTI_register_memory - short result", call_time);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed registering short result: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
        log_debug(rpc_debug_level, "Allocated short result buffer");
    }

cleanup:
    return rc;
}

static int teardown_short_result(nssi_request *request)
{
    int rc = NSSI_OK;

    if (nssi_config.use_buffer_queue) {
        trios_buffer_queue_push(&recv_bq, request->short_result_hdl);
    } else {
        rc=NNTI_free(request->short_result_hdl);
    }
    request->short_result_hdl=NULL;

    return rc;
}

static int setup_long_args(nssi_request *request, nssi_size args_size)
{
    int rc = NSSI_OK;

    log_debug(rpc_debug_level,"allocating space for args");
    /* allocate memory for the encoded arguments. The request
     * structure keeps track of the buffer so it can free
     * the memory later. */
    rc=NNTI_alloc(
            &transports[request->svc->transport_id],
            args_size,
            1,
            NNTI_GET_SRC,
            request->long_args_hdl);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed registering long args: %s",
                nnti_err_str(rc));
        goto cleanup;
    }

cleanup:
    return rc;
}

static int teardown_long_args(nssi_request *req)
{
    int rc = NSSI_OK;
    char *buf;

    trios_declare_timer(call_time);

    if (!req->use_long_args) {
        return(rc);
    }

    log_debug(rpc_debug_level, "waiting for long args");

    if (req->app_pinned_long_args == FALSE) {
        rc=NNTI_free(req->long_args_hdl);
        if (rc != NNTI_OK) {
            log_error(rpc_debug_level, "failed unregistering long args: %s",
                    nnti_err_str(rc));
        }
    }

    return rc;
}

/**
 * @brief Wait for a request to complete.
 *
 * @param req the handle of the pending request.
 */
int nssi_timedwait(nssi_request *req, int timeout, int *remote_rc)
{
    int rc=NSSI_OK;
    log_level debug_level = rpc_debug_level;
    NNTI_status_t status;

    int retries=0;

    trios_declare_timer(total_time);
    trios_declare_timer(call_time);

    trios_start_timer(total_time);

    switch (req->status) {

        case NSSI_REQUEST_NULL:

            goto request_null;

        case NSSI_REQUEST_ERROR:
            log_debug(rpc_debug_level,"timedwait finished");
            *remote_rc = req->error_code;

            goto request_error;

        case NSSI_REQUEST_COMPLETE:
            log_debug(rpc_debug_level,"timedwait finished");
            *remote_rc = NSSI_OK;

            goto request_complete;

        default:
            log_debug(debug_level, "calling NNTI_wait for request");
            /* Times out after DEFAULT_RPC_TIMEOUT */
            trios_start_timer(call_time);
            rc=NNTI_wait(
                    &req->short_request_wr,
                    timeout,
                    &status);
            trios_stop_timer("NNTI_wait - send req", call_time);
            if (status.result == NNTI_EDROPPED) {
                log_debug(LOG_ALL, "request dropped");
            }
            while (status.result == NNTI_EDROPPED) {
                if (retries++ < MAX_RETRIES) {
                    trios_start_timer(call_time);
                    rc=NNTI_send(
                            &req->svc->svc_host,
                            req->short_request_hdl,
                            &req->svc->req_addr,
                            &req->short_request_wr);
                    trios_stop_timer("NNTI_send - resend req", call_time);
                    if (rc != NNTI_OK) {
                        log_error(rpc_debug_level, "failed resending short request: %s",
                                nnti_err_str(rc));
                        break;
                    }
                    trios_start_timer(call_time);
                    rc=NNTI_wait(
                            &req->short_request_wr,
                            timeout,
                            &status);
                    trios_stop_timer("NNTI_wait - resend req", call_time);
                    if (status.result == NNTI_EDROPPED) {
                        log_debug(LOG_ALL, "request dropped");
                    }
                } else {
                    log_error(rpc_debug_level, "retries exhausted: %s",
                            nnti_err_str(rc));
                    break;
                }
            }
            if (rc == NNTI_ETIMEDOUT) {
                log_info(rpc_debug_level, "send request timed out");

                goto wait_timedout;
            }
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed waiting for send: %s",
                        nnti_err_str(rc));

                goto wait_failed;
            }

            /* the service is processing the request. */
            req->status = NSSI_PROCESSING_REQUEST;

            log_debug(debug_level, "calling NNTI_wait for result");
            trios_start_timer(call_time);
            NNTI_create_work_request(
                    req->short_result_hdl,
                    &req->short_result_wr);
            rc=NNTI_wait(
                    &req->short_result_wr,
                    timeout,
                    &status);
            trios_stop_timer("NNTI_wait - short result", call_time);
            if (status.result == NNTI_ETIMEDOUT) {
                log_info(debug_level, "NNTI_wait for result timed out");
                rc = status.result;

                goto wait_timedout;
            }
            if (status.result != NNTI_OK) {
                log_info(debug_level, "NNTI_wait for result failed");
                rc = status.result;
                req->status = NSSI_REQUEST_ERROR;

                goto wait_failed;
            }

            log_debug(debug_level, "Processing result status.offset=%d", status.offset);

            /* we are now ready to process the result */
            req->status = NSSI_PROCESSING_RESULT;

            rc = process_result((char *)status.start+status.offset, req);
            if (rc != NSSI_OK) {
                log_fatal(debug_level,"unable to process result");
                fprint_NNTI_status(logger_get_file(), "status", "FATAL", &status);

                goto wait_failed;
            }

            NNTI_destroy_work_request(
                    &req->short_result_wr);

            goto wait_success;
    }




request_null:
    return rc;

request_error:
    return rc;

request_bad_status:
    /* this should only execute if something went wrong */
    log_fatal(debug_level,"invalid request status %d",req->status);
    return NSSI_EBADRPC;

request_complete:
    if (req->is_responseless == TRUE) {
        return rc;
    }
    goto buffer_cleanup;

wait_timedout:
    return rc;

wait_failed:
    /* check for an error */
    *remote_rc = req->error_code;
    return rc;

wait_success:
    /* call the callback function associated with the request */
    log_debug(debug_level, "calling callback for wait(): op=%d", req->opcode);
    if (req->callback) {
        req->callback(req);
    }

    /* check for completion */
    log_debug(debug_level,"timedwait finished");
    *remote_rc = NSSI_OK;
    goto buffer_cleanup;

buffer_cleanup:
    if (req->app_pinned_short_request == FALSE) {
        log_debug(debug_level, "Cleanup short request");
        rc = teardown_short_request(req);
        if (rc != NSSI_OK) {
            log_error(debug_level, "failed to cleanup long args");
        }
    }

    if (req->app_pinned_short_result == FALSE) {
        log_debug(debug_level, "Cleanup short result");
        rc = teardown_short_result(req);
        if (rc != NSSI_OK) {
            log_error(debug_level, "failed to cleanup long args");
        }
    }

    if (req->use_long_args) {
        log_debug(debug_level, "Cleanup long args");
        /* Now we need to clean up the long arguments (if they were used) */
        rc = teardown_long_args(req);
        if (rc != NSSI_OK) {
            log_error(debug_level, "failed to cleanup long args");
        }
    }

    /* If the request has data associated with it, the data should
     * be transferred by now (server would not have sent result).
     * We need to unlink the MD and free the event queue.
     */
    if (req->data != NULL) {
        if (req->app_pinned_bulk_data == FALSE) {
            if ((nssi_config.use_buffer_queue) &&
                (nssi_config.rdma_buffer_queue_buffer_size >= req->data_size)) {
                /* copy the RDMA buffer contents into the user buffer.
                 * we can't tell if the server op was get or put.
                 * if it was get, then this is a waste of time.
                 * if it was put, then this is required.
                 */
                trios_start_timer(call_time);
                memcpy(req->data, NNTI_BUFFER_C_POINTER(req->bulk_data_hdl), req->data_size);
                trios_stop_timer("memcpy bq to target buf", call_time);
            }
            log_debug(debug_level, "Cleanup bulk data");
            teardown_bulk_data(req);
        }
    }

    trios_stop_timer("nssi_timedwait - total", total_time);

    return rc;
}

int nssi_wait(nssi_request *req, int *rc)
{
    log_debug(rpc_debug_level, "calling nssi_timedwait(timeout=-1)");
    return nssi_timedwait(req, -1, rc);
}

/**
 * @brief Wait for any request to complete.
 *
 * A request is not complete unless we receive the short
 * result.
 *
 * @param req_array  @input_type  The array of pending requests.
 * @param size       @input_type  The number of pending requests.
 * @param timeout    @input_type  The time to wait for any one request.
 *
 */
int nssi_waitany(
        nssi_request *req_array,
        nssi_size     req_count,
        int           timeout,
        int          *which,
        int          *remote_rc)
{
	NNTI_result_t nnti_rc=NNTI_OK; /* result of NNTI functions */

    int rc = NSSI_OK;  /* return code */
    int rc2;
    int i;

	int errors  =0;
	int timeouts=0;

    NNTI_work_request_t **work_requests=new NNTI_work_request_t*[req_count];
    uint32_t              nnti_which=0;
    NNTI_status_t         status;

    log_level debug_level = rpc_debug_level;

    int retries=0;

    trios_declare_timer(call_time);


    for (i=0;i<req_count;i++) {
        switch (req_array[i].status) {
			case NSSI_REQUEST_NULL:
			case NSSI_REQUEST_ERROR:
			case NSSI_REQUEST_COMPLETE:
				work_requests[i]=NULL;
				break;
            case NSSI_SENDING_REQUEST:
                work_requests[i]=&req_array[i].short_request_wr;
                break;
            case NSSI_PROCESSING_REQUEST:
    	        nnti_rc = NNTI_create_work_request(
                                req_array[i].short_result_hdl,
                                &req_array[i].short_result_wr);
				work_requests[i]=&req_array[i].short_result_wr;
				break;
            default:
                break;
        }
    }

wait_again:
    nnti_rc=NNTI_waitany(
    		work_requests,
    		req_count,
    		timeout,
    		&nnti_which,
    		&status);

    *which=nnti_which;

    if (status.result == NNTI_EDROPPED) {
        if (retries++ < MAX_RETRIES) {
            trios_start_timer(call_time);
            rc=NNTI_send(
                    &req_array[*which].svc->svc_host,
                    req_array[*which].short_request_hdl,
                    &req_array[*which].svc->req_addr,
                    &req_array[*which].short_request_wr);
            trios_stop_timer("NNTI_send - resend req", call_time);
            if (rc != NNTI_OK) {
                log_error(rpc_debug_level, "failed resending short request: %s",
                        nnti_err_str(rc));
                req_array[*which].status = NSSI_REQUEST_ERROR;
            } else {
                goto wait_again;
            }
        } else {
            log_error(rpc_debug_level, "retries exhausted: %s",
                    nnti_err_str(rc));
        }
    }
    if (status.result == NNTI_ETIMEDOUT) {
        log_info(debug_level, "NNTI_wait for result timed out");
        rc = status.result;
        return rc;
    }
    if (status.result != NNTI_OK) {
        log_info(debug_level, "NNTI_wait for result failed");
        rc = status.result;
        req_array[*which].status = NSSI_REQUEST_ERROR;
        work_requests[*which]=NULL;
    }

    if (req_array[*which].status == NSSI_SENDING_REQUEST) {
        /* the send is complete.  now wait for the result. */
        req_array[*which].status = NSSI_PROCESSING_REQUEST;

        nnti_rc = NNTI_create_work_request(
                        req_array[*which].short_result_hdl,
                        &req_array[*which].short_result_wr);
        work_requests[*which]=&req_array[*which].short_result_wr;

        goto wait_again;
    }

    log_debug(debug_level, "Processing result status.offset=%d", status.offset);

    /* we are now ready to process the result */
    req_array[*which].status = NSSI_PROCESSING_RESULT;
    rc = process_result((char *)status.start+status.offset, &req_array[*which]);
    if (rc != NSSI_OK) {
        log_fatal(debug_level,"unable to process result");
        fprint_NNTI_status(logger_get_file(), "status", "FATAL", &status);
        return rc;
    }

    if (req_array[*which].use_long_args) {
    	log_debug(debug_level, "Cleanup long args");
    	/* Now we need to clean up the long arguments (if they were used) */
    	rc = teardown_long_args(&req_array[*which]);
    	if (rc != NSSI_OK) {
    		log_error(debug_level, "failed to cleanup long args");
    		return NSSI_EBADRPC;
    	}
    }

cleanup:

	for (i=0;i<req_count;i++) {
		if ((work_requests[i]) && (work_requests[i] == &req_array[i].short_result_wr)) {
			NNTI_destroy_work_request(work_requests[i]);
		}
	}

	delete(work_requests);

	/* Did we get here because of an error in this code? */
	if (rc != NSSI_OK) {
		return rc;
	}

	/* call the callback function associated with the request */
	log_debug(debug_level, "calling callback for wait(): op=%d", req_array[*which].opcode);
	if (req_array[*which].callback) {
	    req_array[*which].callback(&req_array[*which]);
	}

    if (req_array[*which].app_pinned_short_result == FALSE) {
        teardown_short_result(&req_array[*which]);
    }

	/* If the request has data associated with it, the data should
	 * be transferred by now (server would not have sent result).
	 * We need to unlink the MD and free the event queue.
	 */
	if (req_array[*which].data != NULL) {
        if (req_array[*which].app_pinned_bulk_data == FALSE) {
            if ((nssi_config.use_buffer_queue) &&
                (nssi_config.rdma_buffer_queue_buffer_size >= req_array[*which].data_size)) {
                /* copy the RDMA buffer contents into the user buffer.
                 * we can't tell if the server op was get or put.
                 * if it was get, then this is a waste of time.
                 * if it was put, then this is required.
                 */
                trios_start_timer(call_time);
                memcpy(req_array[*which].data, NNTI_BUFFER_C_POINTER(req_array[*which].bulk_data_hdl), req_array[*which].data_size);
                trios_stop_timer("memcpy bq to target buf", call_time);
            }
            teardown_bulk_data(&req_array[*which]);
        }
	}


	/* check for an error */
	if (req_array[*which].status == NSSI_REQUEST_ERROR) {
		*remote_rc = req_array[*which].error_code;
		return rc;
	}

	/* check for completion */
	if (req_array[*which].status == NSSI_REQUEST_COMPLETE) {
		log_debug(debug_level,"timedwait finished");
		*remote_rc = NSSI_OK;
		return rc;
	}

	/* this should only execute if something went wrong */
	log_fatal(debug_level,"invalid request status %d",req_array[*which].status);
	return NSSI_EBADRPC;

}

/**
 * @brief Wait for all requests to complete.
 *
 * A request is not complete unless we receive the short
 * result.
 *
 * @param req_array  @input_type  The array of pending requests.
 * @param size       @input_type  The number of pending requests.
 * @param timeout    @input_type  The time to wait for any one request.
 *
 */
int nssi_waitall(
        nssi_request *req_array,
        nssi_size     req_count,
        int           timeout)
{
	NNTI_result_t nnti_rc=NNTI_OK; /* result of NNTI functions */

    int rc = NSSI_OK;  /* return code */
    int rc2;
    int i;

	int errors  =0;
	int timeouts=0;

    NNTI_work_request_t **work_requests=new NNTI_work_request_t*[req_count];
    NNTI_status_t         status;

    uint32_t which=0;

    log_level debug_level = rpc_debug_level;

    int retries=0;

    trios_declare_timer(call_time);


    for (i=0;i<req_count;i++) {
        switch (req_array[i].status) {
			case NSSI_REQUEST_NULL:
			case NSSI_REQUEST_ERROR:
			case NSSI_REQUEST_COMPLETE:
				work_requests[i]=NULL;
				break;
            case NSSI_SENDING_REQUEST:
                work_requests[i]=&req_array[i].short_request_wr;
                break;
            case NSSI_PROCESSING_REQUEST:
                nnti_rc = NNTI_create_work_request(
                                req_array[i].short_result_hdl,
                                &req_array[i].short_result_wr);
                work_requests[i]=&req_array[i].short_result_wr;
                break;
            default:
                break;
        }
    }

    for (i=0;i<req_count;i++) {
wait_again:
        nnti_rc=NNTI_waitany(
    			work_requests,
    			req_count,
    			timeout,
    			&which,
    			&status);
        if (status.result == NNTI_EDROPPED) {
            if (retries++ < MAX_RETRIES) {
                trios_start_timer(call_time);
                rc=NNTI_send(
                        &req_array[which].svc->svc_host,
                        req_array[which].short_request_hdl,
                        &req_array[which].svc->req_addr,
                        &req_array[which].short_request_wr);
                trios_stop_timer("NNTI_send - resend req", call_time);
                if (rc != NNTI_OK) {
                    log_error(rpc_debug_level, "failed resending short request: %s",
                            nnti_err_str(rc));
                    break;
                }
                trios_start_timer(call_time);
                goto wait_again;
            } else {
                log_error(rpc_debug_level, "retries exhausted: %s",
                        nnti_err_str(rc));
            }
        }
        if (status.result == NNTI_ETIMEDOUT) {
            log_info(debug_level, "work_request[%d] timed out", which);
            req_array[which].status = NSSI_REQUEST_TIMEDOUT;
            timeouts++;
            continue;
        }
        if (status.result != NNTI_OK) {
            log_info(debug_level, "work_request[%d] failed", which);
            req_array[which].status = NSSI_REQUEST_ERROR;
            errors++;
            goto wr_error;
        }

        if (req_array[which].status == NSSI_SENDING_REQUEST) {
            /* the send is complete.  now wait for the result. */
            req_array[which].status = NSSI_PROCESSING_REQUEST;

            nnti_rc = NNTI_create_work_request(
                            req_array[which].short_result_hdl,
                            &req_array[which].short_result_wr);
            work_requests[which]=&req_array[which].short_result_wr;

            goto wait_again;
        }

        log_debug(debug_level, "Processing result status.offset=%d", status.offset);

        /* we are now ready to process the result */
        req_array[which].status = NSSI_PROCESSING_RESULT;
        rc = process_result((char *)status.start+status.offset, &req_array[which]);
        if (rc != NSSI_OK) {
            log_fatal(debug_level,"unable to process result");
            fprint_NNTI_status(logger_get_file(), "status", "FATAL", &status);
            req_array[which].status = NSSI_REQUEST_ERROR;
            errors++;
        }

        if (req_array[which].use_long_args) {
            log_debug(debug_level, "Cleanup long args");
            /* Now we need to clean up the long arguments (if they were used) */
            rc = teardown_long_args(&req_array[which]);
            if (rc != NSSI_OK) {
                log_error(debug_level, "failed to cleanup long args");
                req_array[which].status = NSSI_REQUEST_ERROR;
                errors++;
            }
        }

        NNTI_destroy_work_request(work_requests[which]);

        /* call the callback function associated with the request */
        log_debug(debug_level, "calling req_array[%d] callback: op=%d", which, req_array[which].opcode);
        if (req_array[which].callback) {
            req_array[which].callback(&req_array[which]);
        }

        if (req_array[which].app_pinned_short_result == FALSE) {
            teardown_short_result(&req_array[which]);
        }

        /* If the request has data associated with it, the data should
         * be transferred by now (server would not have sent result).
         * We need to unlink the MD and free the event queue.
         */
        if (req_array[which].data != NULL) {
            if (req_array[which].app_pinned_bulk_data == FALSE) {
                if ((nssi_config.use_buffer_queue) &&
                    (nssi_config.rdma_buffer_queue_buffer_size >= req_array[which].data_size)) {
                    /* copy the RDMA buffer contents into the user buffer.
                     * we can't tell if the server op was get or put.
                     * if it was get, then this is a waste of time.
                     * if it was put, then this is required.
                     */
                    trios_start_timer(call_time);
                    memcpy(req_array[which].data, NNTI_BUFFER_C_POINTER(req_array[which].bulk_data_hdl), req_array[which].data_size);
                    trios_stop_timer("memcpy bq to target buf", call_time);
                }
                teardown_bulk_data(&req_array[which]);
            }
        }

wr_error:
        work_requests[which] = NULL;
    }

	delete(work_requests);

    if (errors > 0) {
        return NSSI_EBADRPC;
    }

    if (timeouts > 0) {
        return NSSI_ETIMEDOUT;
    }

    return rc;
}

/**
 * @brief Initialize portals data structures associated with an
 * RPC request.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method checks to see if the operation arguments for a
 * request arguments fit into the short request buffer (along
 * with the request header).  If not, this method creates the necessary
 * portals data structures on the client that enable the server to
 * fetch the arguments explicitely.
 *
 */
static int encode_args(
        const nssi_service *svc,
        void *args,
        NNTI_buffer_t *short_req_hdl,
        nssi_request_header *header,
        nssi_request *request)
{

    int rc = NSSI_OK;  /* return code */

    trios_declare_timer(call_time);

    /* xdrs for the header and the args. */
    XDR hdr_xdrs, args_xdrs;

    nssi_size hdr_size    = 0;
    nssi_size args_size   = 0;
    nssi_size data_size   = 0;
    nssi_size data_offset = 0;
    nssi_size remaining   = 0;

    char *buf=NULL;


    if (logging_debug(rpc_debug_level)) {
        fprint_nssi_service(logger_get_file(), "service", "DEBUG", svc);
    }
    /* get the encoding functions for the request */
    switch (svc->rpc_encode) {
        case NSSI_RPC_XDR:
            rc = nssi_lookup_xdr_encoding(request->opcode,
                    &request->xdr_encode_args,
                    &request->xdr_encode_data,
                    &request->xdr_decode_result);
            if (rc != NSSI_OK) {
                log_error(rpc_debug_level,
                "could not find xdr encoding functions");
                goto cleanup;
            }
            break;

        default:
            log_error(rpc_debug_level, "invalid encoding scheme");
            rc = NSSI_EENCODE;
            goto cleanup;
    }

    /* set the request ID for the header */
    header->id = request->id;


    /* the size of the encoded header */
    hdr_size = xdr_sizeof((xdrproc_t)&xdr_nssi_request_header, header);

    if (args == NULL) {
        /* there are no args */
        header->fetch_args = FALSE;
        args_size = 0;
    } else {
        /* the size of the encoded args */
        args_size = xdr_sizeof(request->xdr_encode_args, args);
        /* what's left after encoding the header */
        remaining = NNTI_BUFFER_SIZE(short_req_hdl) - hdr_size;
        if (remaining > args_size) {
            /* args will fit in the request. fetch not required. */
            header->fetch_args = FALSE;
        } else {
            /* args will not fit in the request. server must fetch. */
            header->fetch_args = TRUE;

            if (request->is_responseless == TRUE) {
                log_debug(rpc_debug_level, "illegal combination: responseless + long args");
                return NSSI_EINVAL;
            }
        }
    }

    log_debug(rpc_debug_level, "initial hdr_size=%d", hdr_size);

    /**
     *   If the arguments are too large to fit in an
     *   \ref short request buffer, the client stores
     *   excess arguments and waits for the server to
     *   "fetch" the arguments using the \ref nssi_get_data method.
     */
    if (header->fetch_args == TRUE) {

        log_debug(rpc_debug_level,"putting args (len=%d) "
                "in long request", args_size);

        request->use_long_args=1;

        /* the app may provide its own registered long args buffer.  check here. */
        if (request->app_pinned_long_args == FALSE) {
            setup_long_args(request, args_size);
        }
        header->args_addr=*request->long_args_hdl;

        /* create an xdr memory stream for the encoded args */
        xdrmem_create(&args_xdrs, NNTI_BUFFER_C_POINTER(request->long_args_hdl),
                NNTI_BUFFER_SIZE(request->long_args_hdl), XDR_ENCODE);

        /* if we get here, args should not be NULL */
        assert(args);

        /* encode the args  */
        log_debug(rpc_debug_level,"encoding args");
        trios_start_timer(call_time);
        if (! request->xdr_encode_args(&args_xdrs, args)) {
            log_fatal(rpc_debug_level,"failed to encode the args");
            return NSSI_EENCODE;
        }
        trios_stop_timer("xdr_encode_args - encode", call_time);

    } else if ((header->fetch_args == FALSE) && (args != NULL)) {

        log_debug(rpc_debug_level,"putting args (len=%d) "
                "in short request", args_size);

        request->use_long_args=0;

        /* create an xdr memory stream inside the short request for the encoded args */
        xdrmem_create(&args_xdrs, NNTI_BUFFER_C_POINTER(short_req_hdl) + hdr_size,
                NNTI_BUFFER_SIZE(short_req_hdl) - hdr_size, XDR_ENCODE);

        /* if we get here, args should not be NULL */
        assert(args);

        /* encode the args  */
        log_debug(rpc_debug_level,"encoding args");
        trios_start_timer(call_time);
        if (! request->xdr_encode_args(&args_xdrs, args)) {
            log_fatal(rpc_debug_level,"failed to encode the args");
            return NSSI_EENCODE;
        }
        trios_stop_timer("xdr_encode_args - encode", call_time);

    }

    /* if the args are being fetched, then the header has changed.
     * recalculate the size of the encoded header. */
    if (header->fetch_args == TRUE) {
        hdr_size = xdr_sizeof((xdrproc_t)&xdr_nssi_request_header, header);
        log_debug(rpc_debug_level, "using long args. recalculated hdr_size=%d", hdr_size);
    }

    if (nssi_config.put_data_in_request == FALSE) {
        header->fetch_data = TRUE;
    } else {
        /* the size of the bulk data buffer */
        data_size=NNTI_BUFFER_SIZE(&header->data_addr);
        if (data_size > 0) {
            /* what's left after encoding the header into the short request */
            remaining = NNTI_BUFFER_SIZE(short_req_hdl) - hdr_size;
            data_offset = hdr_size;
            if (header->fetch_args == FALSE) {
                /* the args are in the short request.  reduce the remaining space by the size of the encoded args. */
                remaining -= args_size;
                data_offset += args_size;
            }
            if (remaining > data_size) {
                /* data will fit in the request. fetch not required. */
                header->fetch_data = FALSE;
            } else {
                /* data will not fit in the request.  server must fetch. */
                header->fetch_data = TRUE;
            }
        } else {
            /* there is no data */
            header->fetch_data = FALSE;
        }

        if ((data_size > 0) && (header->fetch_data == FALSE)) {
            log_debug(rpc_debug_level,"putting data (size=%d, offset=%d) in short request", data_size, data_offset);
            /* copy small data into the short request */
            memcpy(NNTI_BUFFER_C_POINTER(short_req_hdl)+data_offset, NNTI_BUFFER_C_POINTER(&header->data_addr), data_size);
        }
    }


    /*
     * encode the header
     */
    log_debug(rpc_debug_level,"encoding request header");

    /* create an xdr memory stream for the short request buffer */
    xdrmem_create(&hdr_xdrs, NNTI_BUFFER_C_POINTER(short_req_hdl),
            NNTI_BUFFER_SIZE(short_req_hdl), XDR_ENCODE);

    trios_start_timer(call_time);
    if (! xdr_nssi_request_header(&hdr_xdrs, header)) {
        log_fatal(rpc_debug_level,"failed to encode the request header");
        return NSSI_EENCODE;
    }
    trios_stop_timer("xdr_nssi_request_header - encode", call_time);


    log_debug(rpc_debug_level, "short_req_size=%d, "
            "hdr_size(final)=%d, args_size=%d, data_size=%d",
            NNTI_BUFFER_SIZE(short_req_hdl), hdr_size, args_size, data_size);

    /* print the header for debugging */
    if (logging_debug(rpc_debug_level)) {
        fprint_nssi_request_header(logger_get_file(), "header",
                "nssi_request_header", header);
    }

cleanup:

    /* done! */
    return rc;
}


/**
 * @brief Setup an RPC request.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method populates an RPC request.  All input parameters are copied
 * into the request.  The request can be used to send the RPC to the service
 * and wait for the result
 *
 * @param svc           @input descriptor for the NSSI service.
 * @param opcode        @input descriptor for the remote method.
 * @param args          @input pointer to the arguments.
 * @param data          @input pointer to data (for bulk data transfers).
 * @param data_size     @input length of data buffer
 * @param result        @input where to put results.
 * @param request       @output The request handle (used to test for
 *                              completion).
 */
int nssi_create_request(
        const nssi_service *svc,
        const int opcode,
        void *args,
        void *data,
        uint32_t data_size,
        void *result,
        nssi_request *request)
{
    /* local count does not need protection */
    static int64_t local_count;

    /* local variables */
    int rc=NSSI_OK;  /* return code */

    log_level debug_level = rpc_debug_level;

    char *buf;

    int retries=0;

    trios_declare_timer(total_time);

    unsigned long len=0;

    trios_start_timer(total_time);

    log_debug(rpc_debug_level, "enter");

    /* increment global counter */
    local_count = nthread_counter_increment(&request_count);
    if (local_count == -1) {
        log_error(debug_level, "Unable to increment counter");
        goto cleanup;
    }

    /*------ Initialize variables and buffers ------*/
    memset(request, 0, sizeof(nssi_request));

    /* set request fields */
    request->svc                      = svc;
    request->id                       = local_count;           /* id of the request (used for debugging) */
    request->opcode                   = opcode;                /* operation ID */
    request->is_responseless          = FALSE;
    request->is_sync                  = FALSE;
    request->sync_timeout             = -1;
    request->error_code               = NSSI_OK;               /* return code of remote method */
    request->status                   = NSSI_SENDING_REQUEST;  /* status of this request */
    request->app_pinned_short_request = FALSE;
    request->short_request_hdl        = &request->short_request;
    request->app_pinned_long_args     = FALSE;
    request->long_args_hdl            = &request->long_args;
    request->args                     = args;                  /* where to find the RPC args */
    request->app_pinned_bulk_data     = FALSE;
    request->bulk_data_hdl            = &request->bulk_data;
    request->data_size                = data_size;
    request->data                     = (data_size > 0) ? data : NULL;
    request->app_pinned_short_result  = FALSE;
    request->short_result_hdl         = &request->short_result;
    request->result                   = result;                /* where to put the result */

cleanup:
    trios_stop_timer("nssi_create_request - total", total_time);
    log_debug(rpc_debug_level, "exit");

    return rc;
}

/**
 * @brief Send an RPC request to an NSSI server.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method encodes and transfers an RPC request header and
 * operation arguments to an NSSI server using NNTI. If the
 * arguments are sufficiently small, \b nssi_send_request sends
 * the request header and the arguments in a single message.
 * If the arguments are large (i.e., too large for the request buffer),
 * the server fetches the arguments from a client-side buffer.
 *
 * @param request   @inout The request handle (used to test for completion).
 */
int nssi_send_request(
        nssi_request *request)
{
    /* local variables */
    int rc=NSSI_OK;  /* return code */

    NNTI_status_t status;

    log_level debug_level = rpc_debug_level;

    nssi_request_header header;   /* the request header */
    char *buf;

    trios_declare_timer(total_time);
    trios_declare_timer(call_time);

    unsigned long len=0;

    trios_start_timer(total_time);

    log_debug(rpc_debug_level, "enter");

    if (request->is_responseless == TRUE) {
        /* quick check for illegal request setup */
        if (request->app_pinned_short_result == TRUE) {
            log_debug(rpc_debug_level, "illegal combination: responseless + app pinned short result");
            return NSSI_EINVAL;
        }
        if ((request->app_pinned_short_result == FALSE) && (request->result)) {
            log_debug(rpc_debug_level, "illegal combination: responseless + short result NOT NULL");
            return NSSI_EINVAL;
        }
        if (request->data_size > 0) {
            log_debug(rpc_debug_level, "illegal combination: responseless + bulk data");
            return NSSI_EINVAL;
        }
    }


    /*------ Initialize variables and buffers ------*/
    memset(&header, 0, sizeof(nssi_request_header));

    /* set the opcode for the request header */
    header.opcode = request->opcode;

    /* the app may provide its own registered request buffer.  check here. */
    if (request->app_pinned_short_request == FALSE) {
        /* the app didn't provide a registered request buffer.  get one here. */
        setup_short_request(request);
    }

    if (request->data_size > 0) {

        /* the app may provide its own registered bulk data buffer.  check here. */
        if (request->app_pinned_bulk_data == FALSE) {
            /* the app didn't provide a registered bulk data buffer.  get one here. */
            setup_bulk_data(request);
        }
        header.data_addr=*request->bulk_data_hdl;
    }

    header.is_responseless = request->is_responseless;
    if (request->is_responseless == FALSE) {
        /* the app may provide its own registered result buffer.  check here. */
        if (request->app_pinned_short_result == FALSE) {
            /* the app didn't provide a registered result buffer.  get one here. */
            setup_short_result(request);
        }
        header.res_addr=*request->short_result_hdl;
    }

    /* --- encode the arguments (might place args in the short request) --- */
    trios_start_timer(call_time);
    rc = encode_args(request->svc, request->args, request->short_request_hdl, &header, request);
    trios_stop_timer("encode_args", call_time);
    if (rc != NSSI_OK) {
        log_fatal(rpc_debug_level, "unable to encode arguments");
        goto cleanup;
    }

    log_debug(LOG_OFF, "header.id=%d", header.id);

    /* get the number of valid bytes in the request */
    trios_start_timer(call_time);
    len = xdr_sizeof((xdrproc_t)&xdr_nssi_request_header, &header);
    trios_stop_timer("xdr_sizeof - request header", call_time);

    /* if the args in the short request, add args len */
    if (!header.fetch_args) {
        len += NNTI_BUFFER_SIZE(&header.args_addr);
    }

    /* send the encoded short request buffer to the server */
    log_debug(rpc_debug_level,"sending short request, id=%lu, len=%d", header.id, len);

    trios_start_timer(call_time);
    rc=NNTI_send(
            &request->svc->svc_host,
            request->short_request_hdl,
            &request->svc->req_addr,
            &request->short_request_wr);
    trios_stop_timer("NNTI_send - send req", call_time);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "failed sending short request: %s",
                nnti_err_str(rc));
        goto cleanup;
    }
    log_debug(rpc_debug_level,"message sent");

    if (request->is_responseless == TRUE) {
        trios_start_timer(call_time);
        rc=NNTI_wait(
                &request->short_request_wr,
                -1,
                &status);
        trios_stop_timer("NNTI_wait - send req", call_time);
        if (status.result == NNTI_EDROPPED) {
            log_debug(rpc_debug_level, "request dropped");
        }
        if (rc != NNTI_OK) {
            log_debug(rpc_debug_level, "failed waiting for send: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }

    if (request->is_sync == TRUE) {
        int remote_rc=NSSI_OK;
        rc=nssi_timedwait(request, request->sync_timeout, &remote_rc);
    }

cleanup:
    log_debug(rpc_debug_level, "finished nssi_send_request (req.opcode=%d, req.id=%d)",
            request->opcode, request->id);

    if (rc == NSSI_OK) {
        if (request->is_responseless == TRUE) {
            request->status = NSSI_REQUEST_COMPLETE;
        }
    } else {
        /* change the state of the pending request */
        request->status = NSSI_REQUEST_ERROR;

        if (request->app_pinned_short_request == FALSE) {
            teardown_short_request(request);
        }

        if (request->app_pinned_bulk_data == FALSE) {
            if (request->data_size > 0) {
                teardown_bulk_data(request);
            }
        }
        if (request->app_pinned_short_result == FALSE) {
            teardown_short_result(request);
        }
    }

    trios_stop_timer("nssi_send_request - total", total_time);

    return rc;
}


/**
 * @brief Setup and send an RPC request to an NSSI server.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method populates, encodes and transfers an RPC request header and
 * operation arguments to an NSSI server using NNTI. If the
 * arguments are sufficiently small, \b nssi_call_rpc sends
 * the request header and the arguments in a single message.
 * If the arguments are large (i.e., too large for the request buffer),
 * the server fetches the arguments from a client-side buffer.
 *
 * @param svc           @input descriptor for the NSSI service.
 * @param opcode        @input descriptor for the remote method.
 * @param args          @input pointer to the arguments.
 * @param data          @input pointer to data (for bulk data transfers).
 * @param data_size     @input length of data buffer
 * @param result        @input where to put results.
 * @param request       @output The request handle (used to test for
 *                              completion).
 */
int nssi_call_rpc(
        const nssi_service *svc,
        const int opcode,
        void *args,
        void *data,
        uint32_t data_size,
        void *result,
        nssi_request *request)
{
    /* local variables */
    int rc=NSSI_OK;  /* return code */

    log_level debug_level = rpc_debug_level;

    nssi_request_header header;   /* the request header */

    int retries=0;

    trios_declare_timer(total_time);
    trios_declare_timer(call_time);

    unsigned long len=0;

    trios_start_timer(total_time);

    log_debug(rpc_debug_level, "enter");

    rc=nssi_create_request(svc, opcode, args, data, data_size, result, request);
    if (rc != NSSI_OK) {
        goto cleanup;
    }
    rc=nssi_send_request(request);
    if (rc != NSSI_OK) {
        goto cleanup;
    }

cleanup:

    return rc;

}


/**
 * @brief Send an RPC request to an NSSI server.
 *
 * @ingroup rpc_ptl_impl
 *
 * This method encodes and transfers an RPC request header and
 * operation arguments to an NSSI server using Portals. If the
 * arguments are sufficiently small, \b nssi_call_rpc sends
 * the request header and the arguments in a single message.
 * If the arguments are large (i.e., too large for the request buffer),
 * the server to fetch the arguments from a client-side portal.
 *
 * @param rpc           @input descriptor for the remote method.
 * @param args          @input pointer to the arguments.
 * @param data          @input pointer to data (for bulk data transfers).
 * @param data_size     @input length of data buffer
 * @param result        @input where to put results.
 * @param req           @output The request handle (used to test for
 *                              completion).
 */

int nssi_call_rpc_sync(
        const nssi_service *svc,
        const int opcode,
        void *args,
        void *data,
        uint32_t data_size,
        void *result)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    rc = nssi_call_rpc(svc, opcode, args, data, data_size, result, &req);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "failed calling nssi_call_rpc: %s", nssi_err_str(rc));
        goto cleanup;
    }

    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(rpc_debug_level, "failed waiting for result from server: %s", nssi_err_str(rc));
        rc = rc2;
        goto cleanup;
    }

    if (rc != NSSI_OK) {
        log_warn(rpc_debug_level, "error in remote method: %s",
                nssi_err_str(rc));
        return rc;
    }

cleanup:
    return rc;
}

int nssi_multicast_rpc(
        const nssi_service *svcs,
        const nssi_size num_svcs,
        const int opcode,
        void *args,
        void *data,
        uint32_t data_size,
        void *results,
        nssi_size result_size, /* size in bytes of the result */
        nssi_request *requests)
{
    nssi_size i=0;
    int rc=0;

    for (i=0;i<num_svcs;i++) {
        if (results == NULL) {
            rc = nssi_call_rpc(&svcs[i], opcode, args, data, data_size, NULL, &requests[i]);
        } else {
            rc = nssi_call_rpc(&svcs[i], opcode, args, data, data_size, (char *)results+(i*result_size), &requests[i]);
        }
        if (rc != NSSI_OK) {
            log_error(rpc_debug_level, "failed calling nssi_call_rpc: %s", nssi_err_str(rc));
            goto cleanup;
        }
    }

cleanup:
    return rc;
}

int nssi_multicast_rpc_sync(
        const nssi_service *svcs,
        const nssi_size num_svcs,
        const int opcode,
        void *args,
        void *data,
        uint32_t data_size,
        void *results,
        nssi_size result_size) /* size in bytes of the result */
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    std::vector <nssi_request> reqs(num_svcs);

    rc = nssi_multicast_rpc(svcs, num_svcs, opcode, args, data, data_size, results, result_size, &reqs[0]);
    if (rc != NSSI_OK) {
        log_error(rpc_debug_level, "failed calling nssi_call_rpc: %s", nssi_err_str(rc));
        return rc;
    }

    rc2 = nssi_waitall(&reqs[0], num_svcs, -1);
    if (rc2 != NSSI_OK) {
        log_error(rpc_debug_level, "failed waiting for result from server: %s", nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_warn(rpc_debug_level, "error in remote method: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;
}

int nssi_atomic_increment(
        const nssi_service *svc,
        const uint64_t      remote_atomic,
        const uint64_t      local_atomic)
{
    int rc=NSSI_OK;  /* return code */

    NNTI_work_request_t wr;
    NNTI_status_t       status;

    rc=NNTI_atomic_fop(&transports[svc->transport_id],
                       &svc->svc_host,
                       remote_atomic,
                       local_atomic,
                       1,
                       NNTI_ATOMIC_FADD,
                       &wr);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nnti_err_str(rc));
        return rc;
    }
    rc=NNTI_wait(&wr, -1, &status);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nnti_err_str(rc));
        return rc;
    }

    return(rc);
}

int nssi_atomic_decrement(
        const nssi_service *svc,
        const uint64_t      remote_atomic,
        const uint64_t      local_atomic)
{
    int rc=NSSI_OK;  /* return code */

    NNTI_work_request_t wr;
    NNTI_status_t       status;

    rc=NNTI_atomic_fop(&transports[svc->transport_id],
                       &svc->svc_host,
                       remote_atomic,
                       local_atomic,
                       -1,
                       NNTI_ATOMIC_FADD,
                       &wr);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nnti_err_str(rc));
        return rc;
    }
    rc=NNTI_wait(&wr, -1, &status);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "remote method failed: %s",
                nnti_err_str(rc));
        return rc;
    }

    return(rc);
}

int nssi_atomic_read(
        const nssi_service *svc,
        const uint64_t      local_atomic,
        int64_t            *value)
{
    int rc=NSSI_OK;  /* return code */

    NNTI_work_request_t wr;
    NNTI_status_t       status;

    rc=NNTI_atomic_read(
            &transports[svc->transport_id],
            local_atomic,
            value);
    if (rc != NNTI_OK) {
        log_error(rpc_debug_level, "atomic read failed: %s",
                nnti_err_str(rc));
        return rc;
    }

    return(rc);
}
