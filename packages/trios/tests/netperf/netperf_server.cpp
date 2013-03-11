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
/**  @file netperf-server.cpp
 *
 *   @brief Example data transfer server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 */

/**
 * @defgroup netperf_server Data Transfer Server
 *
 * @ingroup netperf_example
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

#include "netperf_util.h"
#include "netperf_threads.h"



#include <netperf_service_args.h>

log_level netperf_debug_level = LOG_UNDEFINED;

#undef GNI_PERF

#ifdef GNI_PERF
#include <gemini.h>
extern gemini_state_t gni_state;
#endif

struct netperf_experiment_data {
    uint32_t      num_reqs;
    int32_t       array_len;
    bool          validate;
    data_array_t *transfer_array;
    data_array_t *validation_array;
};
typedef struct netperf_experiment_data netperf_experiment_data_t;

static netperf_experiment_data_t experiment_data;

/**
 * @brief
 *
 * @param request_id  ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int netperf_experiment_setup_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const netperf_experiment_setup_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = netperf_debug_level;

    experiment_data.num_reqs  = args->num_reqs;
    experiment_data.array_len = args->array_len;
    experiment_data.validate  = args->validate;

    int nbytes = experiment_data.array_len*sizeof(data_t);

    log_debug(debug_level, "starting netperf_experiment_setup_srvr");

    /* allocate and initialize space for the transfer buffers */
    experiment_data.transfer_array=(data_array_t *)malloc(experiment_data.num_reqs*sizeof(data_array_t));
    for (uint32_t i=0; i<experiment_data.num_reqs; i++) {
        experiment_data.transfer_array[i].data_array_t_len = experiment_data.array_len;
        experiment_data.transfer_array[i].data_array_t_val = (data_t *)malloc(nbytes);
        netperf_init_data_array(i, &experiment_data.transfer_array[i]);
    }

    /* allocate and initialize space for the validation buffers */
    experiment_data.validation_array=(data_array_t *)malloc(experiment_data.num_reqs*sizeof(data_array_t));
    for (uint32_t i=0; i<experiment_data.num_reqs; i++) {
        experiment_data.validation_array[i].data_array_t_len = experiment_data.array_len;
        experiment_data.validation_array[i].data_array_t_val = (data_t *)malloc(nbytes);
        netperf_init_data_array(i, &experiment_data.validation_array[i]);
    }

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    return rc;
}


/**
 * @brief
 *
 * @param request_id  ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int netperf_experiment_teardown_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const void *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = netperf_debug_level;

    log_debug(debug_level, "starting netperf_experiment_teardown_srvr");

    /* free space for the transfer buffers */
    for (uint32_t i=0; i<experiment_data.num_reqs; i++) {
        free(experiment_data.transfer_array[i].data_array_t_val);
    }
    free(experiment_data.transfer_array);

    /* free space for the validation buffers */
    for (uint32_t i=0; i<experiment_data.num_reqs; i++) {
        free(experiment_data.validation_array[i].data_array_t_val);
    }
    free(experiment_data.validation_array);

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
 * @param request_id  ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int netperf_write_rdma_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const netperf_write_rdma_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = netperf_debug_level;

    const uint32_t buffer_index = args->buffer_index;
    int nbytes = experiment_data.array_len*sizeof(data_t);

    double get_time=0.0;

    log_debug(debug_level, "starting netperf_write_rdma_srvr");

    log_debug(debug_level, "getting data from client (%s)", caller->url);

    /* now we need to fetch the data from the client */
    get_time = Trios::GetTime();
    rc = nssi_get_data(caller, experiment_data.transfer_array[buffer_index].data_array_t_val, nbytes, data_addr);
    get_time = Trios::GetTime() - get_time;
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not fetch data from client");
        return rc;
    }

    log_debug(LOG_ALL, "nssi_get_data() bandwidth = %e", (double)nbytes/1024/1024/get_time);

    /* Validate the array */
    if (experiment_data.validate && (rc == 0)) {
        rc = netperf_compare_data_arrays(&experiment_data.validation_array[buffer_index], &experiment_data.transfer_array[buffer_index]);

        if (rc != 0) {
            log_warn(debug_level, "Unable to validate array");
        }
    }

    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

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
 * @param request_id  ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int netperf_read_rdma_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const netperf_read_rdma_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc;
    log_level debug_level = netperf_debug_level;

    const uint32_t buffer_index = args->buffer_index;
    int nbytes = experiment_data.array_len*sizeof(data_t);

    double put_time=0.0;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting netperf_read_rdma_srvr");

    log_debug(debug_level, "putting data in client buffer");

    /* now we need to put the data to the client */
    put_time = Trios::GetTime();
    rc = nssi_put_data(caller, experiment_data.transfer_array[buffer_index].data_array_t_val, nbytes, data_addr, -1);
    put_time = Trios::GetTime() - put_time;
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not put data to client");
        return rc;
    }

    log_debug(LOG_ALL, "nssi_put_data() bandwidth = %e", (double)nbytes/1024/1024/put_time);

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    return rc;
}



void make_progress(bool is_idle)
{
    log_debug(netperf_debug_level, "current_time(%llu) is_idle(%llu)", (uint64_t)trios_get_time_ms(), (uint64_t)is_idle);

    return;
}


/**
 * @brief The NSSI netperf-server.
 *
 * NSSI has already been initialized and the client already knows the URL of the
 * server.  This function simply registers the server methods and starts the
 * service loop.   The client will send a request to kill the service upon completion.
 *
 */
int netperf_server_main(nssi_rpc_transport transport, int num_threads, MPI_Comm server_comm)
{
    int rc = NSSI_OK;

    nssi_service netperf_svc;
    int server_rank;

    MPI_Comm_rank(server_comm, &server_rank);

    /* options that can be overriden by the command-line */
    std::string server_url(NSSI_URL_LEN, '\0');          /* NNTI-style url of the server */
    std::string logfile("");


    memset(&netperf_svc, 0, sizeof(nssi_service));


//    log_str=logfile.c_str();
//    if (logfile.c_str()[0]=='\0') {
//        log_str=NULL;
//    }
//    /* initialize and enable logging */
//    logger_init((log_level)verbose, NULL);
//    debug_level = (log_level)verbose;


    /* initialize the nssi service */
    rc = nssi_service_init(transport, NSSI_SHORT_REQUEST_SIZE, &netperf_svc);
    if (rc != NSSI_OK) {
        log_error(netperf_debug_level, "could not init netperf_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    // register callbacks for the service methods
    NSSI_REGISTER_SERVER_STUB(NETPERF_EXPERIMENT_SETUP_OP, netperf_experiment_setup_srvr, netperf_experiment_setup_args, void);
    NSSI_REGISTER_SERVER_STUB(NETPERF_EXPERIMENT_TEARDOWN_OP, netperf_experiment_teardown_srvr, void, void);
    NSSI_REGISTER_SERVER_STUB(NETPERF_WRITE_RDMA_OP, netperf_write_rdma_srvr, netperf_write_rdma_args, void);
    NSSI_REGISTER_SERVER_STUB(NETPERF_READ_RDMA_OP, netperf_read_rdma_srvr, netperf_read_rdma_args, void);

#ifdef GNI_PERF
    gemini_init_state(server_comm, &gni_state);
#endif

    // Get the Server URL
    std::string url(NSSI_URL_LEN, '\0');
    nssi_get_url(transport, &url[0], NSSI_URL_LEN);


    // Set the maxumum number of requests to handle (-1 == infinite)
    netperf_svc.max_reqs = -1;
    //        netperf_svc.progress_callback=(uint64_t)make_progress;
    //        netperf_svc.progress_callback_timeout=100;

    log_debug(netperf_debug_level, "Starting Server: url = %s", url.c_str());

    // Tell the NSSI server to output log data
    //rpc_debug_level = netperf_debug_level;

    // Start processing requests, the client will send a request to exit when done.
    // If we're running multithreaded, we need to replace the process_request function with the
    // enqueue_reqs function and start the process_queue_reqs thread.

    if (num_threads > 0) {
        log_debug(netperf_debug_level, "Starting server threads");

        netperf_start_server_threads(num_threads, 1000);

        // The main thread will execute unit it receives a "kill" request
        rc = nssi_service_start_wfn(&netperf_svc, &netperf_enqueue_rpc_request);

        netperf_cancel_server_threads();
    }
    else {
        rc = nssi_service_start(&netperf_svc);
        if (rc != NSSI_OK) {
            log_info(netperf_debug_level, "exited netperf_svc: %s",
                    nssi_err_str(rc));
        }
    }



    sleep(5);

    /* shutdown the netperf_svc */
    log_debug(netperf_debug_level, "shutting down service library");
    nssi_service_fini(&netperf_svc);


    return rc;
}

/**
 * @}
 */
