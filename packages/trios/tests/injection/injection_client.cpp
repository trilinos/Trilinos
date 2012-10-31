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
/*
 * injection_client.cpp
 *
 *  Created on: Nov 14, 2011
 *      Author: thkorde
 */

/**
 * @defgroup injection_client Data Transfer Client
 *
 * @ingroup injection_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_logger.h"
#include "Trios_timer.h"

#include "injection_test.h"


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

#include <injection_service_args.h>


log_level client_debug_level = LOG_UNDEFINED;

extern int print_args(
        std::ostream &out,
        const struct injection_args &args,
        const char *prefix);



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
int injection_empty_request(
    const nssi_service *svc,
    nssi_request *req)
{
    int rc = NSSI_OK;
    log_level debug_level = client_debug_level;

    log_debug(debug_level, "Calling RPC for INJECTION_EMPTY_REQUEST_OP");

    /* call the remote methods */
    rc = nssi_call_rpc(svc, INJECTION_EMPTY_REQUEST_OP, NULL, NULL, 0, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call INJECTION_EMPTY_REQUEST_OP: %s",
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
int injection_empty_request_blk(
    const nssi_service *svc)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = injection_empty_request(svc, &req);
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
 * @param injection_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
injection_client_main (struct injection_args &args, nssi_service &injection_svc, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int i,j;
    double time;
    double start_time;
    std::ofstream result_stream;
    log_level debug_level = LOG_WARN;
    int client_rank, client_size;

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
    NSSI_REGISTER_CLIENT_STUB(INJECTION_EMPTY_REQUEST_OP, void, void, void);


    log_debug(debug_level, "%d: Starting experiment loop", client_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                client_rank, i, args.num_reqs);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();

        switch (args.io_method) {
            case INJECTION_EMPTY_REQUEST_SYNC:
            {
                log_debug(debug_level, "Starting INJECTION_EMPTY_REQUEST_SYNC");

                for (j=0; j<args.num_reqs; j++) {
                    rc = injection_empty_request_blk(&injection_svc);
                    if (rc != NSSI_OK) {
                        log_error(client_debug_level, "could not transfer data: %s",
                                nssi_err_str(rc));
                        goto abort;
                    }
                }

                break;
            }

            case INJECTION_EMPTY_REQUEST_ASYNC:
            {
                /* submit requests */
                log_debug(debug_level, "Starting INJECTION_EMPTY_REQUEST_ASYNC");

                for (j=0; j<args.num_reqs; j++) {
                    rc = injection_empty_request(&injection_svc, &reqs[j]);
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

        // Time
        if (i == 0) timings_desc.push_back("Total time for the experiment");
        timings.push_back(time);

        // Throughput
        if (i == 0) timings_desc.push_back("Requests/sec");
        double tput = args.num_reqs/time;
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
