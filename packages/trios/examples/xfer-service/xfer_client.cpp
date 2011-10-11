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
/**
 * @file xfer-client.cpp
 *
 */

/**
 * @defgroup xfer_client Data Transfer Client
 *
 * @ingroup xfer_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_logger.h"
#include "Trios_timer.h"

#include "xfer_client.h"


#include <iostream>
#include <vector>
#include <ostream>
#include <iomanip>
#include <string>


#include <mpi.h>

#include <math.h>
#include <limits.h>
#include <assert.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "Teuchos_CommandLineProcessor.hpp"





#include <xfer_service_args.h>


log_level client_debug_level = LOG_UNDEFINED;






int verify_result(int len, data_t *result)
{
    if ((int)result->int_val != (len-1)) {
        log_error(client_debug_level, "wrong int result: "
                "expected %d, got %d",
                (int)(len-1), result->int_val);
        return NSSI_EBADRPC;
    };
    if (result->float_val != (float)(len-1)) {
        log_error(client_debug_level, "wrong float result: "
                "expected %f, got %f",
                (float)(len-1), result->float_val);
        return NSSI_EBADRPC;
    };
    if (result->double_val != (double)(len-1)) {
        log_error(client_debug_level, "wrong double result: "
                "expected %g, got %g",
                (double)(len-1), result->double_val);
        return NSSI_EBADRPC;
    };

    return NSSI_OK;
}


/**
 * @brief Asynchronous transfer of \ref data_t array with the request structure.
 *
 * This function marshals the \ref data_array_t array of 16-byte
 * structures with the request and send it to the server using the
 * nssi_call_rpc() function.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_push_clnt(
    const nssi_service *svc,
    const data_array_t *array,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_push_args args;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_push_args));

    args.array.data_array_t_len = array->data_array_t_len;
    args.array.data_array_t_val = array->data_array_t_val;

    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_PUSH, &args, NULL, 0, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_push: %s",
            nssi_err_str(rc));
    }

    return rc;
}

/**
 * @brief Synchronous transfer of data with the request structure.
 *
 * This method attaches the entire request buffer to a request
 * and sends it to the server.  This forces the client to encode
 * the buffer before sending, and it makes the size of the request
 * large, depending on the size of the data array.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 */
int xfer_push_clnt_blk(
    const nssi_service *svc,
    const data_array_t *array)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_push_clnt(svc, array, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call async method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(client_debug_level, "failed waiting for request: %s",
                nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(client_debug_level, "remote method failed: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;

}

/**
 * @brief Asynchronous transfer of data via RDMA.
 *
 * This method sends the length of the data array to the server.  After
 * receiving the request, the server pulls the data from the memory
 * buffer using the nssi_get_data() function.  This function is asynchronous --
 * after sending the request to the server, this function returns the
 * request data structure.  The caller must at some point call the
 * nssi_wait() or nssi_test() function to find out if the request is
 * complete.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The service request data structure.
 */
int xfer_pull_clnt(
    const nssi_service *svc,
    const data_array_t *array,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_pull_args args;
    int nbytes;

    /* the buffer to send to the server */
    const data_t *buf = array->data_array_t_val;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_pull_args));

    /* the only argument to send is the size of the array */
    args.len = array->data_array_t_len;

    /* calculate the number of bytes in the buffer */
    nbytes = args.len*sizeof(data_t);

    /* call the remote methods (send buffer using the data portal) */
    rc = nssi_call_rpc(svc, XFER_PULL,  &args, (char *)buf, nbytes, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_push: %s",
                nssi_err_str(rc));
    }

    return rc;
}

/**
 * @brief Synchronous transfer of data via RDMA.
 *
 * This method sends the length of the data array to the server.  After
 * receiving the request, the server pulls the data from the memory
 * buffer using the nssi_get_data() function.  This function is blocking --
 * it calls the nssi_wait() function to make sure the server completed
 * the request before returning.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 */
int xfer_pull_clnt_blk(
    const nssi_service *svc,
    const data_array_t *array)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_pull_clnt(svc, array, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call async method: %s",
            nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(client_debug_level, "failed waiting for request: %s",
            nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(client_debug_level, "remote method failed: %s",
            nssi_err_str(rc));
        return rc;
    }

    return rc;

}

/**
 * @brief Asynchronous transfer of \ref data_t array with the request structure.  The server sends the data back in the result.
 *
 * This function marshals the \ref data_array_t array of 16-byte
 * structures with the request and send it to the server using the
 * nssi_call_rpc() function.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_roundtrip_clnt(
    const nssi_service *svc,
    const data_array_t *array,
    xfer_roundtrip_res *res,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_roundtrip_args args;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_roundtrip_args));
    memset(res, 0, sizeof(xfer_roundtrip_res));

    args.array.data_array_t_len = array->data_array_t_len;
    args.array.data_array_t_val = array->data_array_t_val;

    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_ROUNDTRIP, &args, NULL, 0, res, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_3: %s",
            nssi_err_str(rc));
    }

    return rc;
}

/**
 * @brief Synchronous transfer of \ref data_t array with the request structure.  The server sends the data back in the result.
 *
 * This method attaches the entire request buffer to a request
 * and sends it to the server.  This forces the client to encode
 * the buffer before sending, and it makes the size of the request
 * large, depending on the size of the data array.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 */
int xfer_roundtrip_clnt_blk(
    const nssi_service *svc,
    const data_array_t *array,
    xfer_roundtrip_res *res)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_roundtrip_clnt(svc, array, res, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call async method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(client_debug_level, "failed waiting for request: %s",
                nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(client_debug_level, "remote method failed: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;

}

/**
 * @brief Asynchronous transfer of \ref data_t array via RDMA.  The server sends the data back in the result.
 *
 * This method sends the length of the data array to the server.  After
 * receiving the request, the server pulls the data from the memory
 * buffer using the nssi_get_data() function.  This function is asynchronous --
 * after sending the request to the server, this function returns the
 * request data structure.  The caller must at some point call the
 * nssi_wait() or nssi_test() function to find out if the request is
 * complete.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_get_clnt(
    const nssi_service *svc,
    const data_array_t *array,
    xfer_get_res *res,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_get_args args;

    int nbytes;

    /* the buffer to send to the server */
    data_t *buf = array->data_array_t_val;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_get_args));
    memset(res, 0, sizeof(xfer_get_res));

    /* the only argument to send is the size of the array */
    args.len = array->data_array_t_len;

    /* calculate the number of bytes in the buffer */
    nbytes = args.len*sizeof(data_t);


    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_GET, &args, buf, nbytes, res, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_4: %s",
            nssi_err_str(rc));
    }

    return rc;
}

/**
 * @brief Synchronous transfer of \ref data_t array via RDMA.  The server sends the data back in the result.
 *
 * This method sends the length of the data array to the server.  After
 * receiving the request, the server pulls the data from the memory
 * buffer using the nssi_get_data() function.  This function is ssynchronous --
 * after sending the request to the server, this function calls nssi_wait()
 * to find out the result of the request.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 */
int xfer_get_clnt_blk(
    const nssi_service *svc,
    const data_array_t *array,
    xfer_get_res *res)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_get_clnt(svc, array, res, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call async method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(client_debug_level, "failed waiting for request: %s",
                nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(client_debug_level, "remote method failed: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;

}

/**
 * @brief Asynchronous transfer of \ref data_t array with the request structure.  The server sends the data back via RDMA.
 *
 * This function marshals the \ref data_array_t array of 16-byte
 * structures with the request and send it to the server using the
 * nssi_call_rpc() function.  The server is expected to send the data
 * back using nssi_put_data().
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_put_clnt(
    const nssi_service *svc,
    const data_array_t *array,
    data_array_t *put_buf,
    xfer_put_res *res,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_put_args args;

    int nbytes;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_put_args));
    memset(res, 0, sizeof(xfer_put_res));

    args.array.data_array_t_len = array->data_array_t_len;
    args.array.data_array_t_val = array->data_array_t_val;

    /* calculate the number of bytes in the buffer */
    nbytes = args.array.data_array_t_len*sizeof(data_t);

    put_buf->data_array_t_len = array->data_array_t_len;
    put_buf->data_array_t_val = (data_t *)malloc(nbytes);

    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_PUT, &args, put_buf->data_array_t_val, nbytes, res, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_5: %s",
            nssi_err_str(rc));
    }

    return rc;
}

/**
 * @brief Synchronous transfer of \ref data_t array with the request structure.  The server sends the data back in the result.
 *
 * This method attaches the entire request buffer to a request
 * and sends it to the server.  This forces the client to encode
 * the buffer before sending, and it makes the size of the request
 * large, depending on the size of the data array.
 *
 * @param svc  The service description of the remote service.
 * @param array The data array to transfer.
 */
int xfer_put_clnt_blk(
    const nssi_service *svc,
    const data_array_t *array,
    data_array_t *put_buf,
    xfer_put_res *res)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_put_clnt(svc, array, put_buf, res, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call async method: %s",
                nssi_err_str(rc));
        return rc;
    }

    /* wait for completion */
    rc2 = nssi_wait(&req, &rc);
    if (rc2 != NSSI_OK) {
        log_error(client_debug_level, "failed waiting for request: %s",
                nssi_err_str(rc));
        return rc2;
    }

    if (rc != NSSI_OK) {
        log_error(client_debug_level, "remote method failed: %s",
                nssi_err_str(rc));
        return rc;
    }

    return rc;

}



int read_contact_info(const char *fname, char *url, int maxlen)
{
    const char *contact_filename=NULL;
    FILE *cf=NULL;

    if ((fname==NULL) || (fname[0]=='\0')) {
        contact_filename=getenv("NNTI_CONTACT_FILENAME");
    } else {
        contact_filename=fname;
    }
    if (contact_filename==NULL) {
        url[0]='\0';
        return(-1);
    }
    cf=fopen(contact_filename, "r");
    if (cf == NULL) {
        url[0]='\0';
        return(1);
    }
    fgets(url, maxlen, cf);
    fclose(cf);

    return(0);
}




/**
 * @brief Main code for data transfer client.
 *
 * @param args      The options for the experiment, set at the command-line
 * @param xfer_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
xfer_client_main (struct xfer_args &args, nssi_service &xfer_svc, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int i,j;
    double time;
    data_array_t array;
    double start_time;
    FILE *result_fp = stdout;
    log_level debug_level = LOG_ALL;
    int client_rank, client_size;

    /* the array of results (for async experiments) */
    std::vector<data_t> results;

    xfer_roundtrip_res *res3=NULL;
    xfer_get_res *res4=NULL;
    xfer_put_res *res5=NULL;
    data_array_t *put_buf=NULL;

    std::vector<double> timings;
    std::vector<string> timings_desc;

    /* unique to each process */
    int num_reqs;


    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);

    /* initialize logger */
//    if (args.logfile.empty()) {
//        logger_init(args.debug_level, NULL);
//        debug_level = args.debug_level;
//    } else {
//        char fn[1024];
//        sprintf(fn, "%s.%03d.log", args.logfile.c_str(), client_rank);
//        logger_init(args.debug_level, fn);
//        debug_level = args.debug_level;
//    }


    /* divide the requests among the processors */
    num_reqs = args.num_reqs;
    /* the array of requests (for async experiments) */
    std::vector < nssi_request > reqs(num_reqs);

    /* open the result file */
    result_fp = stdout;
    if (client_rank == 0) {
        if (!args.result_file.empty()) {
            result_fp = fopen(args.result_file.c_str(), args.result_file_mode.c_str());
            if (result_fp == NULL) {
                log_warn(client_debug_level,
                        "invalid result file:"
                        "defaults to stdout");
                result_fp = stdout;
            }
        }
    }

    /* register the XDR encoding functions */
    NSSI_REGISTER_CLIENT_STUB(1, xfer_push_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(2, xfer_pull_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(3, xfer_roundtrip_args, void, xfer_roundtrip_res);
    NSSI_REGISTER_CLIENT_STUB(4, xfer_get_args, void, xfer_get_res);
    NSSI_REGISTER_CLIENT_STUB(5, xfer_put_args, void, xfer_put_res);


    /* initialize the data array */
    array.data_array_t_len = args.len;
    array.data_array_t_val = (data_t *)malloc(args.len*sizeof(data_t));
    if (array.data_array_t_val == NULL) {
        log_error(client_debug_level, "out of space");
        goto abort;
    }

    /* initialize the data */
    for (i=0; i<args.len; i++) {
        array.data_array_t_val[i].int_val = (int)i;
        array.data_array_t_val[i].float_val = (float)i + 0.1;
        array.data_array_t_val[i].double_val = (double)i + 0.2;
    }

    log_debug(debug_level, "%d: Starting experiment loop", client_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                client_rank, i, args.num_reqs, args.len);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();
        switch (args.io_method) {

            case PUSH_SYNC:
                log_debug(debug_level, "%d: push_sync %d reqs", client_rank, args.num_reqs);
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_push_clnt_blk(&xfer_svc, &array);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }
                break;

            case PUSH_ASYNC:
                log_debug(debug_level, "%d: push_async %d reqs", client_rank, args.num_reqs);
                /* submit requests */
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_push_clnt(&xfer_svc, &array, &reqs[j]);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                log_debug(debug_level, "%d: push_async wait", client_rank, args.num_reqs);
                /* wait for results */
                rc = nssi_waitall(&reqs[0], args.num_reqs, -1);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "error transferring data");
                    goto abort;
                }
                break;

            case PULL_SYNC:
                /* submit requests */
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_pull_clnt_blk(&xfer_svc, &array);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }
                break;

            case PULL_ASYNC:
                /* submit requests */
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_pull_clnt(&xfer_svc, &array, &reqs[j]);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                /* wait for results */
                rc = nssi_waitall(&reqs[0], args.num_reqs, -1);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "error transferring data");
                    goto abort;
                }

                /* store the average time for an async req */
                break;

            case ROUNDTRIP_SYNC:
                /* submit requests */
                res3 = (xfer_roundtrip_res *)malloc(1 * sizeof(xfer_roundtrip_res));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_roundtrip_clnt_blk(&xfer_svc, &array, res3);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                    if (0 != memcmp(array.data_array_t_val, res3->array.data_array_t_val, args.len*sizeof(data_t))) {
                        log_error(client_debug_level, "result array does NOT match original");
//                        for (i=0; i<args.len; i++) {
//                            log_debug(debug_level, "orig array int(%d) float(%f) double(%f)",
//                                    array.data_array_t_val[i].int_val, array.data_array_t_val[i].float_val, array.data_array_t_val[i].double_val);
//                        }
//                        for (i=0; i<args.len; i++) {
//                            log_debug(debug_level, "result array int(%d) float(%f) double(%f)",
//                                    res3->array.data_array_t_val[i].int_val, res3->array.data_array_t_val[i].float_val, res3->array.data_array_t_val[i].double_val);
//                        }
                    }
                    if (res3->array.data_array_t_val) free(res3->array.data_array_t_val);
                }

                break;

            case ROUNDTRIP_ASYNC:
                /* submit requests */
                res3 = (xfer_roundtrip_res *)malloc(args.num_reqs * sizeof(xfer_roundtrip_res));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_roundtrip_clnt(&xfer_svc, &array, &(res3[j]), &reqs[j]);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                /* wait for results */
                rc = nssi_waitall(&reqs[0], args.num_reqs, -1);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "error transferring data");
                    goto abort;
                }

                for (j=0; j<args.num_reqs; j++) {
                    if (res3[j].array.data_array_t_val) free(res3[j].array.data_array_t_val);
                }

                /* store the average time for an async req */
                break;

            case GET_SYNC:
                /* submit requests */
                res4 = (xfer_get_res *)malloc(1 * sizeof(xfer_get_res));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_get_clnt_blk(&xfer_svc, &array, res4);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }

                    if (0 != memcmp(array.data_array_t_val, res4->array.data_array_t_val, args.len*sizeof(data_t))) {
                        log_error(client_debug_level, "result array does NOT match original");
//                        for (int idx=0;idx<args.len;idx++) {
//                            log_debug(debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                                    idx, array.data_array_t_val[idx].int_val,
//                                    idx, array.data_array_t_val[idx].float_val,
//                                    idx, array.data_array_t_val[idx].double_val);
//                            log_debug(debug_level, "res4.array[%d].int_val=%u ; res4.array[%d].float_val=%f ; res4.array[%d].double_val=%f",
//                                    idx, res4->array.data_array_t_val[idx].int_val,
//                                    idx, res4->array.data_array_t_val[idx].float_val,
//                                    idx, res4->array.data_array_t_val[idx].double_val);
//                        }
                    }
                    if (res4->array.data_array_t_val) free(res4->array.data_array_t_val);
                }

                break;

            case GET_ASYNC:
                /* submit requests */
                res4 = (xfer_get_res *)malloc(args.num_reqs * sizeof(xfer_get_res));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_get_clnt(&xfer_svc, &array, &(res4[j]), &reqs[j]);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                /* wait for results */
                rc = nssi_waitall(&reqs[0], args.num_reqs, -1);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "error transferring data");
                    goto abort;
                }

                for (j=0; j<args.num_reqs; j++) {
//                    for (int idx=0;idx<args.len;idx++) {
//                        log_debug(debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                                idx, array.data_array_t_val[idx].int_val,
//                                idx, array.data_array_t_val[idx].float_val,
//                                idx, array.data_array_t_val[idx].double_val);
//                        log_debug(debug_level, "res4.array[%d].int_val=%u ; res4.array[%d].float_val=%f ; res4.array[%d].double_val=%f",
//                                idx, res4[j].array.data_array_t_val[idx].int_val,
//                                idx, res4[j].array.data_array_t_val[idx].float_val,
//                                idx, res4[j].array.data_array_t_val[idx].double_val);
//                    }

                    if (res4[j].array.data_array_t_val) free(res4[j].array.data_array_t_val);
                }

                /* store the average time for an async req */
                break;

            case PUT_SYNC:
                /* submit requests */
                res5 = (xfer_put_res *)malloc(1 * sizeof(xfer_put_res));
                put_buf = (data_array_t *)malloc(1 * sizeof(data_array_t));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_put_clnt_blk(&xfer_svc, &array, put_buf, res5);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                    if (0 != memcmp(array.data_array_t_val, put_buf->data_array_t_val, args.len*sizeof(data_t))) {
                        log_error(client_debug_level, "result array does NOT match original");
//                        for (int idx=0;idx<args.len;idx++) {
//                            log_debug(debug_level, "array[%d].int_val        =%u ; array[%d].float_val        =%f ; array[%d].double_val        =%f",
//                                    idx, array.data_array_t_val[idx].int_val,
//                                    idx, array.data_array_t_val[idx].float_val,
//                                    idx, array.data_array_t_val[idx].double_val);
//                            log_debug(debug_level, "put_buf.array[%d].int_val=%u ; put_buf.array[%d].float_val=%f ; put_buf.array[%d].double_val=%f",
//                                    idx, put_buf->data_array_t_val[idx].int_val,
//                                    idx, put_buf->data_array_t_val[idx].float_val,
//                                    idx, put_buf->data_array_t_val[idx].double_val);
//                        }
                    }
                    if (put_buf->data_array_t_val) free(put_buf->data_array_t_val);
                }
                free(put_buf);

                break;

            case PUT_ASYNC:
                /* submit requests */
                res5 = (xfer_put_res *)malloc(args.num_reqs * sizeof(xfer_put_res));
                put_buf = (data_array_t *)malloc(args.num_reqs * sizeof(data_array_t));
                for (j=0; j<args.num_reqs; j++) {
                    rc = xfer_put_clnt(&xfer_svc, &array, &(put_buf[j]), &(res5[j]), &reqs[j]);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                /* wait for results */
                rc = nssi_waitall(&reqs[0], args.num_reqs, -1);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "error transferring data");
                    goto abort;
                }

                for (j=0; j<args.num_reqs; j++) {
                    if (put_buf[j].data_array_t_val) free(put_buf[j].data_array_t_val);
                }
                free(put_buf);

                /* store the average time for an async req */
                break;

            default:
                log_error(client_debug_level, "unrecognized experiment type");
                return -1;
        }

        log_debug(debug_level, "%d: Finished inner loop", client_rank);
        MPI_Barrier(client_comm);
        time = Trios::GetTime() - start_time;

        // Put the number of clients as the first column in the timings info
        if (i == 0) timings_desc.push_back("Number of clients");
        timings.push_back(client_size);

        // Total reqs as the second column of the timings info
        if (i == 0) timings_desc.push_back("Aggregate number of operations per experiment");
        timings.push_back(args.num_reqs*client_size);

        // Number of structures/op
        if (i == 0) timings_desc.push_back("Data structures per operation per client");
        timings.push_back(args.len);

        // Aggregate bytes per experiment
        if (i == 0) timings_desc.push_back("Bytes per experiment");
        double nbytes = args.num_reqs*client_size* (
                sizeof(struct nssi_request) +
                args.len * sizeof(struct data_t) +
                sizeof(int));
        timings.push_back(nbytes);

        // Time
        if (i == 0) timings_desc.push_back("Total time for the experiment");
        timings.push_back(time);

        // Throughput
        if (i == 0) timings_desc.push_back("Throughput MB/sec");
        double tput = nbytes/time/(1024*1024);
        timings.push_back(tput);

        if (client_rank == 0) {
            bool write_header = (i == 0);
            Trios::WriteTimings(cout, "", timings_desc, timings, write_header);
        }
        timings_desc.clear();
        timings.clear();
    }

    log_debug(debug_level, "%d: Finished outer loop", client_rank);

    if (res3 != NULL) free(res3);
    if (res4 != NULL) free(res4);
    if (res5 != NULL) free(res5);

    free(array.data_array_t_val);

    MPI_Barrier(client_comm);
    if (client_rank==0)
        cout << "experiment complete" << endl;


    fclose(result_fp);

    /* sizeof request header */
    log_debug(debug_level, "sizeof(nssi_request_header)=%lu\n",
            (unsigned long)sizeof(nssi_request_header));
    log_debug(debug_level, "sizeof(nssi_result_header)=%lu\n",
            (unsigned long)sizeof(nssi_result_header));


    return 0;

abort:
    exit(2);
}

/**
 * @}
 */
