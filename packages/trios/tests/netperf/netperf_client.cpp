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
 * @file netperf-client.cpp
 *
 */

/**
 * @defgroup netperf_client Data Transfer Client
 *
 * @ingroup netperf_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_logger.h"
#include "Trios_timer.h"

#include "netperf_client.h"


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

#include "netperf_service_args.h"
#include "netperf_util.h"



log_level client_debug_level = LOG_UNDEFINED;

/* prototype for a function to initialize buffers */
extern int print_args(
        std::ostream &out,
        const struct netperf_args &args,
        const char *prefix);

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
int netperf_write_rdma(
    const nssi_service *svc,
    const uint32_t buffer_index,
    nssi_request *req)
{
    int rc = NSSI_OK;
    netperf_write_rdma_args args;
    int nbytes=experiment_data.array_len*sizeof(data_t);

    /* the buffer to send to the server */
    const data_t *buf = experiment_data.transfer_array[buffer_index].data_array_t_val;

    /* initialize the arguments */
    memset(&args, 0, sizeof(netperf_write_rdma_args));

    /* the only argument to send is the size of the array */
    args.buffer_index=buffer_index;

    /* call the remote methods (send buffer using the data portal) */
    rc = nssi_call_rpc(svc, NETPERF_WRITE_RDMA_OP,  &args, (char *)buf, nbytes, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call netperf_write_encode: %s",
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
int netperf_write_rdma_blk(
    const nssi_service *svc,
    const uint32_t buffer_index)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = netperf_write_rdma(svc, buffer_index, &req);
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
int netperf_read_rdma(
    const nssi_service *svc,
    const uint32_t buffer_index,
    nssi_request *req)
{
    int rc = NSSI_OK;
    netperf_read_rdma_args args;
    log_level debug_level = client_debug_level;

    int nbytes=experiment_data.array_len*sizeof(data_t);

    /* initialize the arguments */
    memset(&args, 0, sizeof(netperf_read_rdma_args));

    /* tell the server how about the buffer */
    args.buffer_index=buffer_index;

    data_t *buf=experiment_data.transfer_array[buffer_index].data_array_t_val;

    log_debug(debug_level, "Calling RPC for NETPERF_READ_RDMA_OP, len=%d, val=%p, nbytes=%d",
            experiment_data.transfer_array[buffer_index].data_array_t_len, buf, nbytes);

    /* call the remote methods */
    rc = nssi_call_rpc(svc, NETPERF_READ_RDMA_OP, &args, buf, nbytes, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call netperf_5: %s",
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
int netperf_read_rdma_blk(
    const nssi_service *svc,
    const uint32_t buffer_index)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = netperf_read_rdma(svc, buffer_index, &req);
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


int netperf_send_setup(
    const nssi_service *svc)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    netperf_experiment_setup_args args;
    log_level debug_level = client_debug_level;

    /* initialize the arguments */
    memset(&args, 0, sizeof(netperf_read_rdma_args));
    args.num_reqs  = experiment_data.num_reqs;
    args.array_len = experiment_data.array_len;
    args.validate  = experiment_data.validate;

    log_debug(debug_level, "Calling RPC for NETPERF_EXPERIMENT_SETUP_OP, num_reqs=%d, array_len=%p, validate=%d",
            experiment_data.num_reqs, experiment_data.array_len, experiment_data.validate);

    /* call the remote methods */
    rc = nssi_call_rpc(svc, NETPERF_EXPERIMENT_SETUP_OP, &args, NULL, 0, NULL, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call netperf_5: %s",
            nssi_err_str(rc));
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

int netperf_send_teardown(
    const nssi_service *svc)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    log_level debug_level = client_debug_level;

    log_debug(debug_level, "Calling RPC for NETPERF_EXPERIMENT_TEARDOWN_OP");

    /* call the remote methods */
    rc = nssi_call_rpc(svc, NETPERF_EXPERIMENT_TEARDOWN_OP, NULL, NULL, 0, NULL, &req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call netperf_5: %s",
            nssi_err_str(rc));
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

int netperf_experiment_setup(
        const uint32_t num_reqs,
        const int32_t  array_len,
        const bool     validate)
{
    int rc = NSSI_OK;
    log_level debug_level = client_debug_level;

    experiment_data.num_reqs  = num_reqs;
    experiment_data.array_len = array_len;
    experiment_data.validate  = validate;

    int nbytes = experiment_data.array_len*sizeof(data_t);

    log_debug(debug_level, "starting netperf_experiment_setup");

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

    return rc;
}


int netperf_experiment_teardown()
{
    int rc = NSSI_OK;
    log_level debug_level = client_debug_level;

    log_debug(debug_level, "starting netperf_experiment_teardown");

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
 * @param netperf_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
netperf_client_main (struct netperf_args &args, nssi_service &netperf_svc, MPI_Comm client_comm)
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
    NSSI_REGISTER_CLIENT_STUB(NETPERF_EXPERIMENT_SETUP_OP, netperf_experiment_setup_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(NETPERF_EXPERIMENT_TEARDOWN_OP, void, void, void);
    NSSI_REGISTER_CLIENT_STUB(NETPERF_WRITE_RDMA_OP, netperf_write_rdma_args, void, void);
    NSSI_REGISTER_CLIENT_STUB(NETPERF_READ_RDMA_OP, netperf_read_rdma_args, void, void);

    netperf_experiment_setup(args.num_reqs, args.len, args.validate_flag);
    netperf_send_setup(&netperf_svc);

    log_debug(debug_level, "%d: Starting experiment loop", client_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                client_rank, i, args.num_reqs, args.len);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();
        switch (args.io_method) {

        case NETPERF_WRITE_RDMA_SYNC:
        {
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                rc = netperf_write_rdma_blk(&netperf_svc, j);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }
            }
            break;
        }

        case NETPERF_WRITE_RDMA_ASYNC:
        {
            /* submit requests */
            for (j=0; j<args.num_reqs; j++) {
                rc = netperf_write_rdma(&netperf_svc, j, &reqs[j]);
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

        case NETPERF_READ_RDMA_SYNC:
        {
            log_debug(debug_level, "Starting NETPERF_READ_RDMA_SYNC");

            // initialize the validation array, if necessary
            for (j=0; j<args.num_reqs; j++) {
                rc = netperf_read_rdma_blk(&netperf_svc, j);
                if (rc != NSSI_OK) {
                    log_error(client_debug_level, "could not transfer data: %s",
                            nssi_err_str(rc));
                    goto abort;
                }

                if (args.validate_flag) {
                    rc = netperf_compare_data_arrays(&experiment_data.validation_array[j], &experiment_data.transfer_array[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }
            }

            break;
        }

        case NETPERF_READ_RDMA_ASYNC:
        {
            /* submit requests */

            for (j=0; j<args.num_reqs; j++) {

                rc = netperf_read_rdma(&netperf_svc, j, &reqs[j]);
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
                for (j=0; j<args.num_reqs; j++) {
                    rc = netperf_compare_data_arrays(&experiment_data.validation_array[j], &experiment_data.transfer_array[j]);
                    if (rc != 0) {
                        log_error(client_debug_level, "Validation failed");
                        MPI_Abort(MPI_COMM_WORLD, rc);
                    }
                }
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
                sizeof(struct netperf_write_rdma_args) +
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

    netperf_experiment_teardown();
    netperf_send_teardown(&netperf_svc);

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
