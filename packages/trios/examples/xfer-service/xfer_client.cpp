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
#include <cassert>


#include <mpi.h>

#include <math.h>
#include <limits.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include "Teuchos_CommandLineProcessor.hpp"

#include "xfer_service_args.h"
#include "xfer_util.h"



log_level client_debug_level = LOG_UNDEFINED;

/* prototype for a function to initialize buffers */
extern int print_args(
        std::ostream &out,
        const struct xfer_args &args,
        const char *prefix);



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
int xfer_write_encode(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_write_encode_args args;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_write_encode_args));
    args.len = len;
    args.seed = seed;
    args.validate = validate;

    args.array.data_array_t_len = array->data_array_t_len;
    args.array.data_array_t_val = array->data_array_t_val;

    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_WRITE_ENCODE_OP, &args, NULL, 0, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_write_encode: %s",
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
int xfer_write_encode_blk(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_write_encode(svc, len, seed, validate, array, &req);
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
int xfer_write_rdma(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_write_rdma_args args;
    int nbytes;

    /* the buffer to send to the server */
    const data_t *buf = array->data_array_t_val;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_write_rdma_args));

    /* the only argument to send is the size of the array */
    args.len = len;
    args.seed = seed;
    args.validate = validate;

    /* calculate the number of bytes in the buffer */
    nbytes = args.len*sizeof(data_t);

    /* call the remote methods (send buffer using the data portal) */
    rc = nssi_call_rpc(svc, XFER_WRITE_RDMA_OP,  &args, (char *)buf, nbytes, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call xfer_write_encode: %s",
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
int xfer_write_rdma_blk(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_write_rdma(svc, len, seed, validate, array, &req);
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
int xfer_read_encode(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    xfer_read_encode_res *res,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_read_encode_args args;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_read_encode_args));
    memset(res, 0, sizeof(xfer_read_encode_res));

    /* initialize the arguments */
    args.len = len;
    args.seed = seed;
    args.validate = validate;

    /* call the remote methods, array comes back in res */
    rc = nssi_call_rpc(svc, XFER_READ_ENCODE_OP, &args, NULL, 0, res, req);
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
int xfer_read_encode_blk(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    xfer_read_encode_res *res)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_read_encode(svc, len, seed, validate, res, &req);
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
int xfer_read_rdma(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    data_array_t *put_buf,
    nssi_request *req)
{
    int rc = NSSI_OK;
    xfer_read_rdma_args args;
    log_level debug_level = client_debug_level;

    int nbytes;

    /* initialize the arguments */
    memset(&args, 0, sizeof(xfer_read_rdma_args));

    assert(len == put_buf->data_array_t_len);


    /* tell the server how about the buffer */
    args.len = len;
    args.seed = seed;
    args.validate = validate;

    /* calculate the number of bytes in the buffer */
    nbytes = args.len*sizeof(data_t);

    data_t *buf;
    //buf = (data_t *)malloc(nbytes);
    buf = put_buf->data_array_t_val;

    log_debug(debug_level, "Calling RPC for XFER_READ_RDMA_OP, len=%d, val=%p, nbytes=%d",
            put_buf->data_array_t_len, buf, nbytes);

    /* call the remote methods */
    rc = nssi_call_rpc(svc, XFER_READ_RDMA_OP, &args, buf, nbytes, NULL, req);
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
int xfer_read_rdma_blk(
    const nssi_service *svc,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    data_array_t *put_buf)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = xfer_read_rdma(svc, len, seed, validate, put_buf, &req);
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
    if (fgets(url, maxlen, cf) == NULL) {
        log_error(client_debug_level, "failed to read URL from %s", fname);
    }
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
    double start_time;
    std::ofstream result_stream;
    log_level debug_level = LOG_WARN;
    int client_rank, client_size;

    /* the array of results (for async experiments) */
    std::vector<data_t> results;

    xfer_read_encode_res *res4=NULL;

    std::vector<double> timings;
    std::vector<string> timings_desc;

    /* unique to each process */
    int num_reqs;


    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);


    /* divide the requests among the processors */
    num_reqs = args.num_reqs;
    /* the array of requests (for async experiments) */
    std::vector < nssi_request > reqs(num_reqs);

    /* open the result file */
    if (client_rank == 0) {

        if (!args.result_file.empty()) {

            if (args.result_file_mode.compare("a") == 0)
                result_stream.open(args.result_file.c_str(), fstream::out | fstream::app);
            else
                result_stream.open(args.result_file.c_str(), fstream::out);

            if (!result_stream.is_open()) {
                log_warn(client_debug_level,
                        "invalid result file:"
                        "defaults to stdout");
            }
        }
    }

    /* register the XDR encoding functions */
    NSSI_REGISTER_CLIENT_STUB(XFER_WRITE_ENCODE_OP, xfer_write_encode_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(XFER_WRITE_RDMA_OP, xfer_write_rdma_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(XFER_READ_ENCODE_OP, xfer_read_encode_args, void, xfer_read_encode_res);
    NSSI_REGISTER_CLIENT_STUB(XFER_READ_RDMA_OP, xfer_read_rdma_args, void, void);


    /* allocate space for the data arrays used for testing */
    std::vector <data_array_t> array_vec;
    array_vec.reserve(args.num_reqs);   // async jobs need num_reqs arrays

    for (i=0; i<args.num_reqs; i++) {
        array_vec[i].data_array_t_len = args.len;
        array_vec[i].data_array_t_val = (data_t *)malloc(args.len*sizeof(data_t));
        if (array_vec[i].data_array_t_val == NULL) {
            log_error(client_debug_level, "out of space");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        memset(array_vec[i].data_array_t_val, 0, args.len*sizeof(data_t));
    }


    log_debug(debug_level, "%d: Starting experiment loop", client_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                client_rank, i, args.num_reqs, args.len);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();
        switch (args.io_method) {

        case XFER_WRITE_ENCODE_SYNC:
        {
            log_debug(debug_level, "%d: push_sync %d reqs", client_rank, args.num_reqs);
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_write_encode_blk(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }
            break;
        }

        case XFER_WRITE_ENCODE_ASYNC:
        {
            log_debug(debug_level, "%d: push_async %d reqs", client_rank, args.num_reqs);
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_write_encode(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j], &reqs[j]);
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
        }

        case XFER_WRITE_RDMA_SYNC:
        {
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_write_rdma_blk(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }
            break;
        }

        case XFER_WRITE_RDMA_ASYNC:
        {
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_write_rdma(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j], &reqs[j]);
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
        }

        case XFER_READ_ENCODE_SYNC:
        {
            /* submit requests */
            res4 = (xfer_read_encode_res *)malloc(1 * sizeof(xfer_read_encode_res));
            for (j=0; j<args.num_reqs; j++) {
                rc = xfer_read_encode_blk(&xfer_svc, args.len, j, args.validate_flag, res4);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }

                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                    rc = xfer_compare_data_arrays(&res4->array, &array_vec[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }

                if (res4->array.data_array_t_val) free(res4->array.data_array_t_val);
            }

            break;
        }

        case XFER_READ_ENCODE_ASYNC:
        {
            /* submit requests */
            res4 = (xfer_read_encode_res *)malloc(args.num_reqs * sizeof(xfer_read_encode_res));
            for (j=0; j<args.num_reqs; j++) {
                rc = xfer_read_encode(&xfer_svc, args.len, j, args.validate_flag, &(res4[j]), &reqs[j]);
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


            // Validate each of the results
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                    rc = xfer_compare_data_arrays(&res4[j].array, &array_vec[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }

                if (res4[j].array.data_array_t_val) free(res4[j].array.data_array_t_val);
            }

            /* store the average time for an async req */
            break;
        }

        case XFER_READ_RDMA_SYNC:
        {
            log_debug(debug_level, "Starting XFER_READ_RDMA_SYNC");

            // initialize the validation array, if necessary
            data_array_t tmp_array;
            if (args.validate_flag) {
                tmp_array.data_array_t_len = args.len;
                tmp_array.data_array_t_val = (data_t *)malloc(args.len * sizeof(data_t));
            }

            for (j=0; j<args.num_reqs; j++) {
                rc = xfer_read_rdma_blk(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }

                if (args.validate_flag) {

                    xfer_init_data_array(j, &tmp_array);

                    rc = xfer_compare_data_arrays(&tmp_array, &array_vec[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }
            }

            if (args.validate_flag) {
                free(tmp_array.data_array_t_val);
            }

            break;
        }

        case XFER_READ_RDMA_ASYNC:
        {
            /* submit requests */

            for (j=0; j<args.num_reqs; j++) {

                rc = xfer_read_rdma(&xfer_svc, args.len, j, args.validate_flag, &array_vec[j], &reqs[j]);
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

            if (args.validate_flag) {
                // initialize the validation array, if necessary
                data_array_t tmp_array;
                tmp_array.data_array_t_len = args.len;
                tmp_array.data_array_t_val = (data_t *)malloc(args.len * sizeof(data_t));

                for (j=0; j<args.num_reqs; j++) {
                    // initialize the temporary array
                    xfer_init_data_array(j, &tmp_array);

                    // compare the arrays
                    rc = xfer_compare_data_arrays(&tmp_array, &array_vec[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }

                free(tmp_array.data_array_t_val);
            }

            /* store the average time for an async req */
            break;
        }

        default:
        {
            log_error(client_debug_level, "unrecognized experiment type");
            return -1;
        }
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
                sizeof(struct nssi_request_header) +
                sizeof(struct xfer_write_rdma_args) +
                args.len * sizeof(struct data_t) +
                sizeof(struct nssi_result_header));
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
            if (result_stream.is_open() && write_header) print_args(result_stream, args, "%");
            Trios::WriteTimings(result_stream.is_open() ? result_stream : cout, "", timings_desc, timings, write_header);
        }
        timings_desc.clear();
        timings.clear();
    }

    log_debug(debug_level, "%d: Finished outer loop", client_rank);

    // Clean up the data arrays
    for (i=0; i<args.num_reqs; i++) {
        free (array_vec[i].data_array_t_val);
    }

    if (res4 != NULL) free(res4);


    MPI_Barrier(client_comm);
    if (client_rank==0)
        cout << "experiment complete" << endl;


    result_stream.close();

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
