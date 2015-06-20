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
/**  @file opreg-server.cpp
 *
 *   @brief Example data transfer server.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 */

/**
 * @defgroup opreg_server Data Transfer Server
 *
 * @ingroup opreg_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nssi_server.h"
#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nssi_debug.h"
#include "Trios_nssi_fprint_types.h"

#include "opreg_service_args.h"

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>



log_level opreg_debug_level = LOG_UNDEFINED;

uint64_t calc_checksum (const char * buf, const uint64_t size);


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
int opreg_request_processor(
        const unsigned long        request_id,
        const NNTI_peer_t         *caller,
        const opreg_args *args,
        const NNTI_buffer_t       *data_addr,
        const NNTI_buffer_t       *res_addr)
{
    int rc = NNTI_OK;
    log_level debug_level = opreg_debug_level;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting opreg_request_processor()");

    log_debug(debug_level, "args->data.int_val(%d) args->data.float_val(%f) args->data.double_val(%f) args->chksum(%lu)",
            args->data.int_val, args->data.float_val, args->data.double_val, args->chksum);

    /* checksum the data */
    if (args->chksum != calc_checksum((const char*)&args->data, sizeof(args->data))) {
        log_error(debug_level, "client checksum (%lu) != calculated checksum (%lu)",
                args->chksum, calc_checksum((const char*)&args->data, sizeof(args->data)));
        rc=NNTI_EBADRPC;
    }

    rc = nssi_send_result(caller, request_id, rc, NULL, res_addr);

    return rc;
}



/**
 * @brief The NSSI opreg-server.
 *
 * NSSI has already been initialized and the client already knows the URL of the
 * server.  This function simply registers the server methods and starts the
 * service loop.   The client will send a request to kill the service upon completion.
 *
 */
int opreg_c_server_main(nssi_rpc_transport transport, MPI_Comm server_comm)
{
    int rc = NSSI_OK;

    nssi_service opreg_svc;
    int server_rank;

    char url[NSSI_URL_LEN];

    MPI_Comm_rank(server_comm, &server_rank);


    memset(&opreg_svc, 0, sizeof(nssi_service));


    /* initialize the nssi service */
    rc = nssi_service_init(transport, NSSI_SHORT_REQUEST_SIZE, &opreg_svc);
    if (rc != NSSI_OK) {
        log_error(opreg_debug_level, "could not init opreg_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    /* register callbacks for the service methods */
    NSSI_REGISTER_SERVER_STUB(OPREG_REQUEST_OP, opreg_request_processor, opreg_args, void);


    /* Get the Server URL */
    nssi_get_url(transport, &url[0], NSSI_URL_LEN);


    /* Set the maxumum number of requests to handle (-1 == infinite) */
    opreg_svc.max_reqs = -1;

    log_debug(opreg_debug_level, "Starting Server: url = %s", url);

    /* start processing requests, the client will send a request to exit when done */
    rc = nssi_service_start(&opreg_svc);
    if (rc != NSSI_OK) {
        log_info(opreg_debug_level, "exited opreg_svc: %s",
                nssi_err_str(rc));
    }

    sleep(5);

    /* shutdown the opreg_svc */
    log_debug(opreg_debug_level, "shutting down service library");
    nssi_service_fini(&opreg_svc);


    return rc;
}

/**
 * @}
 */
