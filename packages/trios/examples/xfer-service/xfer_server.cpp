/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
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

#include <iostream>
#include <string>

#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>



#include <xfer_service_args.h>

log_level xfer_debug_level = LOG_UNDEFINED;


/**
 * @brief Transfer an array of \ref data_t structures through
 *        the process arguments.
 *
 * This function simply returns a result that tells the client
 * the request was received.  Before this function is called,
 * the request and the data array passed in the args structure
 * is decoded.  This function helps test the overhead of the
 * encoding/decoding process.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_push_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_push_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
//    const int len = args->array.data_array_t_len;
//    const data_t *array = args->array.data_array_t_val;
    log_level debug_level = xfer_debug_level;
    int rc;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting xfer_push_srvr");

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    return rc;
}

/**
 * @brief Transfer an array of \ref data_t structures through
 *        the data portal.
 *
 * This server-side function is a callback that executes when
 * the server receives an \ref XFER_PULL request.  This function
 * pulls data from the client using the nssi_get_data() function
 * and returns a result telling the calling process that this
 * function is complete.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_pull_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_pull_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = xfer_debug_level;

    const int len = args->len;
    int nbytes = len*sizeof(data_t);

    log_debug(debug_level, "starting xfer_pull_srvr");


    /* allocate space for the incoming buffer */
    data_t *buf = (data_t *)malloc(nbytes);

    log_debug(debug_level, "getting data from client (%s)", caller->url);

    /* now we need to fetch the data from the client */
    rc = nssi_get_data(caller, buf, nbytes, data_addr);
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not fetch data from client");
        return rc;
    }

    rc = nssi_send_result(caller, request_id, NSSI_OK, NULL, res_addr);

    free(buf);

    return rc;
}

/**
 * @brief Transfer an array of \ref data_t structures through
 *        the process arguments.  Send it back through the result.
 *
 * This function returns in the result the same data that the
 * client sent in the request.  Before this function is called,
 * the request and the data array passed in the args structure
 * is decoded.  When this function calls nssi_send_result(),
 * the data array passed in the results structure
 * is encoded.  This function helps test that both small and
 * large data can be sent via args and results.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_roundtrip_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_roundtrip_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc;
    log_level debug_level = xfer_debug_level;

    int len = 0;
    data_t *array = NULL;
    xfer_roundtrip_res res;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting xfer_roundtrip_srvr (args=%p)", args);

    len   = args->array.data_array_t_len;
    array = args->array.data_array_t_val;

    memset(&res, 0, sizeof(res));

    res.array.data_array_t_val = (data_t *)malloc(len * sizeof(data_t));
    res.array.data_array_t_len = len;

    log_debug(debug_level, "args->array.data_array_t_val=%p, array=%p, res.array.data_array_t_val=%p", args->array.data_array_t_val, array, res.array.data_array_t_val);

    log_debug(debug_level, "calling memcpy");
    memcpy(res.array.data_array_t_val, array, len * sizeof(data_t));
    log_debug(debug_level, "done memcpy");

//    log_debug(debug_level, "mta_get_num_teams ==%d", mta_get_num_teams());
//    log_debug(debug_level, "mta_get_team_index==%d", mta_get_team_index(0));
//    log_debug(debug_level, "mta_get_rt_teamid ==%d", mta_get_rt_teamid());
//    log_debug(debug_level, "yielding!!!!");
//    mta_yield();

//    for (i=0; i<len; i++) {
//        log_debug(LOG_ALL, "args array int(%d) float(%f) double(%f)",
//                args->array.data_array_t_val[i].int_val, args->array.data_array_t_val[i].float_val, args->array.data_array_t_val[i].double_val);
//    }
//    for (i=0; i<len; i++) {
//        log_debug(LOG_ALL, "result array int(%d) float(%f) double(%f)",
//                res.array.data_array_t_val[i].int_val, res.array.data_array_t_val[i].float_val, res.array.data_array_t_val[i].double_val);
//    }

    log_debug(debug_level, "calling nssi_send_result");
    rc = nssi_send_result(caller, request_id, NSSI_OK, &res, res_addr);

    free(res.array.data_array_t_val);

    return rc;
}

/**
 * @brief Transfer an array of \ref data_t structures through
 *        the data portal.  Send it back through the result.
 *
 * This server-side function is a callback that executes when
 * the server receives an \ref XFER_GET request.  This function
 * pulls data from the client using the nssi_get_data() function
 * and returns a result containing the same data.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data.
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_get_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_get_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc = NSSI_OK;
    log_level debug_level = xfer_debug_level;

    const int len = args->len;
    int nbytes = len*sizeof(data_t);

    xfer_get_res res;

    log_debug(debug_level, "starting xfer_get_srvr");

    memset(&res, 0, sizeof(res));

    res.array.data_array_t_val = (data_t *)malloc(nbytes);
    res.array.data_array_t_len = len;

    log_debug(debug_level, "getting data from client (%s)", caller->url);

    /* now we need to fetch the data from the client */
    rc = nssi_get_data(caller, res.array.data_array_t_val, nbytes, data_addr);
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not fetch data from client");
        return rc;
    }

//    for (int idx=0;idx<len;idx++) {
//        log_debug(LOG_ALL, "res.array[%d].int_val=%u ; res.array[%d].float_val=%f ; res.array[%d].double_val=%f",
//                idx, res.array.data_array_t_val[idx].int_val,
//                idx, res.array.data_array_t_val[idx].float_val,
//                idx, res.array.data_array_t_val[idx].double_val);
//    }

    rc = nssi_send_result(caller, request_id, NSSI_OK, &res, res_addr);

    free(res.array.data_array_t_val);

    return rc;
}

/**
 * @brief Transfer an array of \ref data_t structures through
 *        the process arguments.  Send it back through a data portal.
 *
 * This server-side function is a callback that executes when
 * the server receives an \ref XFER_PUT request.  Before this
 * function is called, the request and the data array passed
 * in the args structure is decoded.  This function
 * puts data to the client using the nssi_put_data() function
 * and returns a result telling the calling process that this
 * function is complete.
 *
 * @param request_id   ID of the request.
 * @param caller      The process ID of the calling process.
 * @param args        Arguments passed with the request.
 * @param data_addr   The remote memory descriptor for the data (not used).
 * @param res_addr    The remote memory descriptor for the result.
 */
int xfer_put_srvr(
        const unsigned long request_id,
        const NNTI_peer_t *caller,
        const xfer_put_args *args,
        const NNTI_buffer_t *data_addr,
        const NNTI_buffer_t *res_addr)
{
    int rc;
    log_level debug_level = xfer_debug_level;

    const int len = args->array.data_array_t_len;
    const data_t *array = args->array.data_array_t_val;
    int nbytes = len*sizeof(data_t);
    xfer_put_res res;

    memset(&res, 0, sizeof(res));

    res.len = len;

    /* process array (nothing to do) */
    log_debug(debug_level, "starting xfer_put_srvr");

    /* now we need to put the data to the client */
    rc = nssi_put_data(caller, array, nbytes, data_addr, -1);
    if (rc != NSSI_OK) {
        log_warn(debug_level, "could not put data to client");
        return rc;
    }

//    for (int idx=0;idx<len;idx++) {
//        log_debug(LOG_ALL, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                idx, array[idx].int_val,
//                idx, array[idx].float_val,
//                idx, array[idx].double_val);
//    }

    rc = nssi_send_result(caller, request_id, NSSI_OK, &res, res_addr);

    return rc;
}



void make_progress(bool is_idle)
{
    log_debug(LOG_ALL, "current_time(%llu) is_idle(%llu)", (uint64_t)trios_get_time_ms(), (uint64_t)is_idle);

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
int xfer_server_main(MPI_Comm server_comm)
{
    int rc = NSSI_OK;

    nssi_service xfer_svc;
    log_level debug_level;
    int server_rank;

    MPI_Comm_rank(server_comm, &server_rank);

    /* options that can be overriden by the command-line */
    int verbose = 3;  /* default debug_level */
    int num_threads = 0;
    std::string server_url(NSSI_URL_LEN, '\0');          /* NNTI-style url of the server */
    std::string logfile("");
    const char *log_str=NULL;


    memset(&xfer_svc, 0, sizeof(nssi_service));


    log_str=logfile.c_str();
    if (logfile.c_str()[0]=='\0') {
        log_str=NULL;
    }
    /* initialize and enable logging */
    logger_init((log_level)verbose, NULL);
    debug_level = (log_level)verbose;


    /* initialize the nssi service */
    rc = nssi_service_init(NSSI_DEFAULT_TRANSPORT, NSSI_SHORT_REQUEST_SIZE, &xfer_svc);
    if (rc != NSSI_OK) {
        log_error(xfer_debug_level, "could not init xfer_svc: %s",
                nssi_err_str(rc));
        return -1;
    }

    // register callbacks for the service methods
    NSSI_REGISTER_SERVER_STUB(XFER_PUSH, xfer_push_srvr, xfer_push_args, void);
    NSSI_REGISTER_SERVER_STUB(XFER_PULL, xfer_pull_srvr, xfer_pull_args, void);
    NSSI_REGISTER_SERVER_STUB(XFER_ROUNDTRIP, xfer_roundtrip_srvr, xfer_roundtrip_args, xfer_roundtrip_res);
    NSSI_REGISTER_SERVER_STUB(XFER_GET, xfer_get_srvr, xfer_get_args, xfer_get_res);
    NSSI_REGISTER_SERVER_STUB(XFER_PUT, xfer_put_srvr, xfer_put_args, xfer_put_res);


    // Get the Server URL
    std::string url(NSSI_URL_LEN, '\0');
    nssi_get_url(NSSI_DEFAULT_TRANSPORT, &url[0], NSSI_URL_LEN);


    // Set the maxumum number of requests to handle (-1 == infinite)
    xfer_svc.max_reqs = -1;
    //        xfer_svc.progress_callback=(uint64_t)make_progress;
    //        xfer_svc.progress_callback_timeout=100;

    log_debug(LOG_ALL, "Starting Server: url = %s", url.c_str());

    // start processing requests, the client will send a request to exit when done
    rc = nssi_service_start(&xfer_svc, num_threads);
    if (rc != NSSI_OK) {
        log_info(xfer_debug_level, "exited xfer_svc: %s",
                nssi_err_str(rc));
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
