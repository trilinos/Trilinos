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
/**  @file xfer-server.cpp
 *
 *   @brief Example data transfer server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 */

/**
 * @defgroup xfer_server Data Transfer Server
 *
 * @ingroup xfer_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
//#include "nssi_types.h"
#include "Trios_nssi_server.h"
#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nssi_debug.h"

#include <iostream>
#include <string>

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

#include "xfer_util.h"
#include "xfer_threads.h"

//#include <TPI.hpp>   // Thread Pool Interface (Trilinos Package)



#include <xfer_service_args.h>

log_level xfer_debug_level = LOG_UNDEFINED;


/**
 * @brief Emulate a write operation where all the data is sent
 *        through the function arguments.
 *
 * Transfer an array of data structures to the server through the
 * procedure arguments, forcing the client to encode the array before
 * sending and the server to decode the array when receiving.  We use
 * this method to evaluate the performance of the encoding/decoding the arguments.
 * For large arrays, this method also tests our two-phase transfer protocol in
 * which the client pushes a small header of arguments and lets the server pull the
 * remaining arguments on demand.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_write_encode_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_write_encode_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    const int len = args->len;
    const long int seed = args->seed;
    const bool validate = args->validate;

    int nbytes = len*sizeof(data_t);

    log_level debug_level = xfer_debug_level;
    int rc;

    /* process array (nothing to do) */

    log_debug(debug_level, "starting xfer_write_encode_srvr");

    /* Validate the array that was sent through the args */
    if (validate) {
        data_array_t tmp_array;
        tmp_array.data_array_t_len = len;

        tmp_array.data_array_t_val = (data_t *)malloc(nbytes);

        xfer_init_data_array(seed, &tmp_array);
        rc = xfer_compare_data_arrays(&tmp_array, const_cast<data_array_t *>(&args->array));

        if (rc != 0) {
            log_warn(debug_level, "Unable to validate array");
        }

        free(tmp_array.data_array_t_val);
    }



    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    return rc;
}

/**
 * @brief Transfer an array of data structures to the server using the data channel.
 *
 * This procedure passes the length of the array in the arguments. The server
 * then ``pulls'' the unencoded data from the client using the \ref nssi_get function.
 * This method evaluates the RDMA transfer performance for the \ref nssi_get_data function.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_write_rdma_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_write_rdma_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = xfer_debug_level;

    const int len = args->len;
    const long int seed = args->seed;
    const bool validate = args->validate;

    int nbytes = len*sizeof(data_t);

    log_debug(debug_level, "starting xfer_write_rdma_srvr");

    /* allocate space for the incoming buffer */
    data_array_t array;
    array.data_array_t_len = len;
    array.data_array_t_val = (data_t *)malloc(nbytes);
    memset(array.data_array_t_val, 0, nbytes);

    log_debug(debug_level, "getting data from client (%s)", caller->url);

    /* now we need to fetch the data from the client */
    rc = nssi_get_data(caller, array.data_array_t_val, nbytes, data_addr);
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not fetch data from client");
        return rc;
    }

    /* Validate the array */
    if (validate && (rc == 0)) {
        data_array_t tmp_array;
        tmp_array.data_array_t_len = len;
        tmp_array.data_array_t_val = (data_t *)malloc(nbytes);

        xfer_init_data_array(seed, &tmp_array);
        rc = xfer_compare_data_arrays(&tmp_array, &array);

        if (rc != 0) {
            log_warn(debug_level, "Unable to validate array");
        }

        free(tmp_array.data_array_t_val);
    }

    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    free(array.data_array_t_val);

    return rc;
}



/**
 * @brief Transfer an array of data structures to the client
 *        using the control channel.
 *
 * This method tells the server to send the data array to the
 * client through the result data structure, forcing the server
 * to encode the array before sending and the client to decode the array
 * when receiving. This procedure evaluates the performance of the encoding/decoding
 * the arguments.  For large arrays, this method also tests our two-phase transfer
 * protocol for the result structure in which the server pushes a small header of
 * the result and lets the client pull the remaining result on demand (at the
 * \ref nssi_wait function).
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_read_encode_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_read_encode_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = xfer_debug_level;

    const int len = args->len;
    const int seed = args->seed;
    const bool validate = args->validate;
    int nbytes = len*sizeof(data_t);

    xfer_read_encode_res res;

    log_debug(debug_level, "starting xfer_read_encode_srvr");

    memset(&res, 0, sizeof(res));

    /* allocate space for the outgoing buffer */
    res.array.data_array_t_val = (data_t *)malloc(nbytes);
    res.array.data_array_t_len = len;

    if (validate) {
        xfer_init_data_array(seed, &res.array);
    }
    log_debug(debug_level, "getting data from client (%s)", caller->url);

    rc = nssi_send_result(caller, request_id, NSSI_OK, &res, res_addr);

    free(res.array.data_array_t_val);

    return rc;
}

/**
 * @brief Emulate a read operation where the bulk data is sent using
 *        the \ref nssi_put function.
 *
 * Transfer an array of data structures to the client using the data
 * channel.  This procedure passes the length of the array in the arguments.
 * The server then ``puts'' the unencoded data into the client memory using
 * the \ref nssi_put_data function.  This method evaluates the RDMA
 * transfer performance for \ref nssi_put_data.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_read_rdma_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_read_rdma_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc;
    log_level debug_level = xfer_debug_level;

    const bool validate = args->validate;
    const int seed = args->seed;

    int nbytes = args->len*sizeof(data_t);

    /* allocate space for the outgoing buffer */
    data_array_t array;
    array.data_array_t_len = args->len;
    array.data_array_t_val = (data_t *)malloc(nbytes);

    /* process array (nothing to do) */
    log_debug(debug_level, "starting xfer_read_rdma_srvr");

    /* if we need to validate the array, it needs to be initialized */
    if (validate) {
        xfer_init_data_array(seed, &array);
    }

    log_debug(debug_level, "putting data in client buffer");

    /* now we need to put the data to the client */
    rc = nssi_put_data(caller, array.data_array_t_val, nbytes, data_addr, -1);
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not put data to client");
        return rc;
    }

//    for (int idx=0;idx<len;idx++) {
//        log_debug(debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                idx, array[idx].int_val,
//                idx, array[idx].float_val,
//                idx, array[idx].double_val);
//    }

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    free(array.data_array_t_val);

    return rc;
}



void make_progress(bool is_idle)
{
    log_debug(xfer_debug_level, "current_time(%llu) is_idle(%llu)", (uint64_t)trios_get_time_ms(), (uint64_t)is_idle);

    return;
}


/**
 * @brief The NSSI xfer-server.
 *
 * NSSI has already been initialized and the client already knows the URL of the
 * server.  This function simply registers the server methods and starts the
 * service loop.   The client will send a request to kill the service upon completion.
 *
 */
int xfer_server_main(nssi_rpc_transport transport, int num_threads, MPI_Comm server_comm)
{
    int rc = NSSI_OK;

    nssi_service xfer_svc;
    log_level debug_level;
    int server_rank;

    MPI_Comm_rank(server_comm, &server_rank);

    /* options that can be overriden by the command-line */
    int verbose = 3;  /* default debug_level */
    std::string server_url(NSSI_URL_LEN, '\0');          /* NNTI-style url of the server */
    std::string logfile("");
    const char *log_str=NULL;


    memset(&xfer_svc, 0, sizeof(nssi_service));


//    log_str=logfile.c_str();
//    if (logfile.c_str()[0]=='\0') {
//        log_str=NULL;
//    }
//    /* initialize and enable logging */
//    logger_init((log_level)verbose, NULL);
//    debug_level = (log_level)verbose;


    /* initialize the nssi service */
    rc = nssi_service_init(transport, NSSI_SHORT_REQUEST_SIZE, &xfer_svc);
    if (rc != NSSI_OK) {
        log_error(xfer_debug_level, "could not init xfer_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    // register callbacks for the service methods
    NSSI_REGISTER_SERVER_STUB(XFER_WRITE_ENCODE_OP, xfer_write_encode_srvr, xfer_write_encode_args, void);
    NSSI_REGISTER_SERVER_STUB(XFER_WRITE_RDMA_OP, xfer_write_rdma_srvr, xfer_write_rdma_args, void);
    NSSI_REGISTER_SERVER_STUB(XFER_READ_ENCODE_OP, xfer_read_encode_srvr, xfer_read_encode_args, xfer_read_encode_res);
    NSSI_REGISTER_SERVER_STUB(XFER_READ_RDMA_OP, xfer_read_rdma_srvr, xfer_read_rdma_args, void);


    // Get the Server URL
    std::string url(NSSI_URL_LEN, '\0');
    nssi_get_url(transport, &url[0], NSSI_URL_LEN);


    // Set the maxumum number of requests to handle (-1 == infinite)
    xfer_svc.max_reqs = -1;
    //        xfer_svc.progress_callback=(uint64_t)make_progress;
    //        xfer_svc.progress_callback_timeout=100;

    log_debug(xfer_debug_level, "Starting Server: url = %s", url.c_str());

    // Tell the NSSI server to output log data
    //rpc_debug_level = xfer_debug_level;

    // Start processing requests, the client will send a request to exit when done.
    // If we're running multithreaded, we need to replace the process_request function with the
    // enqueue_reqs function and start the process_queue_reqs thread.

    if (num_threads > 0) {
        log_debug(xfer_debug_level, "Starting server threads");

        rc = xfer_start_server_threads(num_threads, 1000);

        // The main thread will execute unit it receives a "kill" request
        rc = nssi_service_start_wfn(&xfer_svc, &xfer_enqueue_rpc_request);

        xfer_cancel_server_threads();
    }
    else {
        rc = nssi_service_start(&xfer_svc);
        if (rc != NSSI_OK) {
            log_info(xfer_debug_level, "exited xfer_svc: %s",
                    nssi_err_str(rc));
        }
    }



    sleep(5);

    /* shutdown the xfer_svc */
    log_debug(xfer_debug_level, "shutting down service library");
    nssi_service_fini(&xfer_svc);


    return rc;
}

/**
 * @}
 */
