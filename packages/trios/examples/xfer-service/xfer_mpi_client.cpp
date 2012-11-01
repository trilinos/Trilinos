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
#include "xfer_debug.h"
#include "xfer_util.h"


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

#include "xfer_service_args.h"



log_level client_debug_level = LOG_UNDEFINED;


extern int print_args(
        std::ostream &out,
        const struct xfer_args &args,
        const char *prefix);



/**
 * @brief Transfer an array of data structures using the MPI send/receive functions.
 *
 *
 * @param len  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_mpi_send(
    const int server_rank,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    data_array_t *array)
{
    int rc = 0;
    MPI_Status status;
    log_level debug_level = xfer_debug_level;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int req_buf[3];
    req_buf[0] = len;
    req_buf[1] = seed;
    req_buf[2] = (int)validate;

    int nbytes = len * sizeof(struct data_t);

    log_debug(debug_level, "%d: Sending send req (len=%d, seed=%d, validate=%d)", rank, len, seed, validate);

    // The request is just an MPI send with the number of elements
    rc = MPI_Send(req_buf, 3, MPI_INT, server_rank, XFER_MPI_SEND_REQ_TAG, MPI_COMM_WORLD);

    log_debug(debug_level, "%d: Waiting for ack from %d to before sending data", rank, server_rank);
    // Receive an ack that the server is ready
    MPI_Recv(&rc, 0, MPI_INT, server_rank, XFER_MPI_SEND_ACK_TAG, MPI_COMM_WORLD, &status);


    log_debug(debug_level, "%d: Send data: %d structs, %d bytes to %d", rank, len, nbytes, server_rank);

    // Now we send the data
    rc = MPI_Send(&array->data_array_t_val[0], nbytes,
            MPI_BYTE, server_rank, XFER_MPI_SEND_DATA_TAG, MPI_COMM_WORLD);

    return rc;
}


/**
 * @brief Transfer an array of data structures using the MPI send/receive functions.
 *
 *
 * @param len  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_mpi_isend(
    const int server_rank,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array,
    MPI_Request *req)
{
    int rc = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    log_level debug_level = xfer_debug_level;

    int req_buf[3];
    req_buf[0] = (int)len;
    req_buf[1] = (int)seed;
    req_buf[2] = (int)validate;

    int nbytes = len * sizeof(data_t);

    log_debug(debug_level, "%d: Sending isend req (len=%d, seed=%d, validate=%d)", rank, len, seed, validate);
    rc = MPI_Send(req_buf, 3, MPI_INT, server_rank, XFER_MPI_ISEND_REQ_TAG, MPI_COMM_WORLD);

    // Receive an ack that the server is ready
    MPI_Recv(&rc, 0, MPI_INT, server_rank, XFER_MPI_ISEND_ACK_TAG, MPI_COMM_WORLD, &status);

    log_debug(debug_level, "%d: Isend data: %d structs, %d bytes", rank, len, nbytes);
    rc = MPI_Isend(array->data_array_t_val, nbytes, MPI_BYTE, server_rank, XFER_MPI_ISEND_DATA_TAG, MPI_COMM_WORLD, req);

    return rc;
}


/**
 * @brief Transfer an array of data structures using the MPI put functions.
 *
 *
 * @param len  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_mpi_put(
    const int server_rank,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array)
{
    int rc = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Win win;

    int req_buf[3];
    req_buf[0] = (int)len;
    req_buf[1] = (int)seed;
    req_buf[2] = (int)validate;
    int nbytes = len * sizeof(data_t);

    log_level debug_level = xfer_debug_level;


    log_debug(debug_level, "%d: Sending put req (len=%d, seed=%d, validate=%d)", rank, len, seed, validate);
    rc = MPI_Send(req_buf, 3, MPI_INT, server_rank, XFER_MPI_PUT_REQ_TAG, MPI_COMM_WORLD);

    // The window on this side is not used to receive data
    MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MPI_Win_fence(0, win);


    log_debug(debug_level, "%d: Put data: %d structs, %d bytes, val[0]=%d", rank, len, nbytes,
            array->data_array_t_val[0].int_val);

    // Put the data from this window into a buffer on the server
    MPI_Put(array->data_array_t_val, nbytes, MPI_BYTE, server_rank, 0, nbytes, MPI_BYTE, win);

    // This makes sure the operations are complete
    MPI_Win_fence(0, win);

    MPI_Win_free(&win);

    return rc;
}


/**
 * @brief Transfer an array of data structures using the MPI put functions.
 *
 *
 * @param len  The service description of the remote service.
 * @param array The data array to transfer.
 * @param req  The \ref nssi_request returned to caller.
 */
int xfer_mpi_get(
    const int server_rank,
    const uint32_t len,
    const uint32_t seed,
    const bool validate,
    const data_array_t *array)
{
    int rc = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Win win;
    MPI_Group comm_group, group;

    int req_buf[3];
    req_buf[0] = (int)len;
    req_buf[1] = (int)seed;
    req_buf[2] = (int)validate;
    int nbytes = len * sizeof(data_t);

    log_level debug_level = xfer_debug_level;


    log_debug(debug_level, "%d: Sending put req (len=%d, seed=%d, validate=%d)", rank, len, seed, validate);
    rc = MPI_Send(req_buf, 3, MPI_INT, server_rank, XFER_MPI_GET_REQ_TAG, MPI_COMM_WORLD);

    // The window on this side is not used to receive data
    MPI_Win_create(NULL, 0, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);

    MPI_Win_fence(0, win);


    log_debug(debug_level, "%d: Get data: %d structs, %d bytes, val[0]=%d", rank, len, nbytes,
            array->data_array_t_val[0].int_val);

    // Put the data from this window into a buffer on the server
    MPI_Get(array->data_array_t_val, nbytes, MPI_BYTE, server_rank, 0, nbytes, MPI_BYTE, win);

    // This makes sure the operations are complete
    MPI_Win_fence(0, win);

    MPI_Win_free(&win);

    return rc;
}

/**
 * @brief Main code for data transfer client.
 *
 * @param args      The options for the experiment, set at the command-line
 * @param xfer_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
xfer_mpi_client_main (struct xfer_args &args, const int server_rank, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int i,j;
    double time;
    double start_time;
    std::ofstream result_stream;
    log_level debug_level = xfer_debug_level;
    int client_rank, client_size;
    int global_rank, global_size;

    /* the array of results (for async experiments) */
    std::vector<data_t> results;

    xfer_read_encode_res *res4=NULL;

    std::vector<double> timings;
    std::vector<string> timings_desc;

    /* unique to each process */
    int num_reqs;

    MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);

    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);

    log_debug(debug_level, "%d: Starting client", global_rank);

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
    std::vector < MPI_Request > reqs(num_reqs);
    std::vector < MPI_Status > status(num_reqs);

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


    log_debug(debug_level, "%d: Starting experiment loop", global_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                global_rank, i, args.num_reqs, args.len);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();
        switch (args.io_method) {

        case XFER_MPI_SEND:
        {
            log_debug(debug_level, "%d: push_sync %d reqs", global_rank, args.num_reqs);
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_mpi_send(server_rank, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != 0) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }
            break;
        }

        case XFER_MPI_ISEND:
        {
            log_debug(debug_level, "%d: push_async %d reqs", global_rank, args.num_reqs);
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_mpi_isend(server_rank, args.len, j, args.validate_flag, &array_vec[j], &reqs[j]);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }

            log_debug(debug_level, "%d: waiting for %d reqs to complete", global_rank, args.num_reqs);
            /* wait for results */
            rc = MPI_Waitall(args.num_reqs, &reqs[0], &status[0]);
            if (rc != MPI_SUCCESS) {
                log_error(client_debug_level, "error transferring data");
                MPI_Abort(MPI_COMM_WORLD, rc);
            }
            break;
        }

        case XFER_MPI_PUT:
        {
            log_debug(debug_level, "%d: put %d reqs", global_rank, args.num_reqs);
            for (j=0; j<args.num_reqs; j++) {
                if (args.validate_flag) {
                    xfer_init_data_array(j, &array_vec[j]);
                }
                rc = xfer_mpi_put(server_rank, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != 0) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }
            break;
        }

        case XFER_MPI_GET:
        {
            log_debug(debug_level, "%d: get %d reqs", global_rank, args.num_reqs);
            for (j=0; j<args.num_reqs; j++) {

                rc = xfer_mpi_get(server_rank, args.len, j, args.validate_flag, &array_vec[j]);
                if (rc != 0) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
                if (args.validate_flag) {
                    rc = xfer_validate_array(j, &array_vec[j]);
                    if (rc == 0) {
                        log_info(debug_level, "Validate Passed");
                    }
                    else {
                        log_error(debug_level, "Validate Failed: rc=%d, array[0].int_val=%d",
                                rc, array_vec[j].data_array_t_val[0].int_val);
                    }
                }
            }
            break;
        }


        default:
        {
            log_error(client_debug_level, "unrecognized experiment type");
            return -1;
        }
        }

        log_debug(debug_level, "%d: Finished inner loop", global_rank);
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

    log_debug(debug_level, "%d: Finished outer loop", global_rank);

    // Clean up the data arrays
    for (i=0; i<args.num_reqs; i++) {
        free (array_vec[i].data_array_t_val);
    }

    if (res4 != NULL) free(res4);

    MPI_Barrier(client_comm);


    // Send a request to exit the server
    int req_buf[3];
    if (client_rank == 0) {
        log_debug(debug_level, "%d: Sending exit request to %d", global_rank, server_rank);
        req_buf[0] = req_buf[1] = req_buf[2] = 0;
        rc = MPI_Send(req_buf, 3, MPI_INT, server_rank, XFER_MPI_FINI_REQ, MPI_COMM_WORLD);
    }


    result_stream.close();


    return 0;

abort:
    exit(2);
}

/**
 * @}
 */
