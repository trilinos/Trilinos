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
/**  @file injection-server.cpp
 *
 *   @brief Example data transfer server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 */

/**
 * @defgroup injection_server Data Transfer Server
 *
 * @ingroup injection_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nssi_server.h"
#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nssi_debug.h"

#include "injection_test.h"

#include <iostream>
#include <string>

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

#include <injection_service_args.h>


log_level injection_debug_level = LOG_UNDEFINED;



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
int injection_empty_request_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const void *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc;
    log_level debug_level = injection_debug_level;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting injection_empty_request_srvr");

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    return rc;
}



void make_progress(bool is_idle)
{
    log_debug(injection_debug_level, "current_time(%llu) is_idle(%llu)", (uint64_t)trios_get_time_ms(), (uint64_t)is_idle);

    return;
}


/**
 * @brief The NSSI injection-server.
 *
 * NSSI has already been initialized and the client already knows the URL of the
 * server.  This function simply registers the server methods and starts the
 * service loop.   The client will send a request to kill the service upon completion.
 *
 */
int injection_server_main(struct injection_args &args, MPI_Comm server_comm)
{
    int rc = NSSI_OK;

    nssi_service injection_svc;
    int server_rank;

    MPI_Comm_rank(server_comm, &server_rank);

    /* options that can be overriden by the command-line */
    std::string server_url(NSSI_URL_LEN, '\0');          /* NNTI-style url of the server */
    std::string logfile("");


    memset(&injection_svc, 0, sizeof(nssi_service));


    /* initialize the nssi service */
    rc = nssi_service_init((nssi_rpc_transport)args.transport, NSSI_SHORT_REQUEST_SIZE, &injection_svc);
    if (rc != NSSI_OK) {
        log_error(injection_debug_level, "could not init injection_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    // register callbacks for the service methods
    NSSI_REGISTER_SERVER_STUB(INJECTION_EMPTY_REQUEST_OP, injection_empty_request_srvr, void, void);


    // Get the Server URL
    std::string url(NSSI_URL_LEN, '\0');
    nssi_get_url((nssi_rpc_transport)args.transport, &url[0], NSSI_URL_LEN);


    // Set the maxumum number of requests to handle (-1 == infinite)
    injection_svc.max_reqs = -1;
    //        injection_svc.progress_callback=(uint64_t)make_progress;
    //        injection_svc.progress_callback_timeout=100;

    log_debug(injection_debug_level, "Starting Server: url = %s", url.c_str());

    // Tell the NSSI server to output log data
    //rpc_debug_level = injection_debug_level;

    // start processing requests, the client will send a request to exit when done
    rc = nssi_service_start(&injection_svc);
    if (rc != NSSI_OK) {
        log_info(injection_debug_level, "exited injection_svc: %s",
                nssi_err_str(rc));
    }

    sleep(5);

    /* shutdown the injection_svc */
    log_debug(injection_debug_level, "shutting down service library");
    nssi_service_fini(&injection_svc);


    return rc;
}

/**
 * @}
 */
