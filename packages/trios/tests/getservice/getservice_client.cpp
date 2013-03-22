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
 * getservice_client.cpp
 *
 *  Created on: Jul 20, 2012
 *      Author: thkorde
 */

/**
 * @defgroup getservice_client Data Transfer Client
 *
 * @ingroup getservice_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_nssi_debug.h"

#include "getservice_test.h"


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


log_level client_debug_level = LOG_UNDEFINED;

extern int print_args(
        std::ostream &out,
        const struct getservice_args &args,
        const char *prefix);


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


int
compare_services(nssi_service &service1, nssi_service &service2)
{
    int rc=NSSI_OK;

    if (service1.rpc_encode != service2.rpc_encode) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }

    if (strcmp(service1.svc_host.url, service2.svc_host.url)) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    if (service1.svc_host.peer.transport_id != service2.svc_host.peer.transport_id) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    switch (service1.svc_host.peer.transport_id) {
        case NNTI_TRANSPORT_PORTALS:
            if (service1.svc_host.peer.NNTI_remote_process_t_u.portals.nid != service2.svc_host.peer.NNTI_remote_process_t_u.portals.nid) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.svc_host.peer.NNTI_remote_process_t_u.portals.pid != service2.svc_host.peer.NNTI_remote_process_t_u.portals.pid) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_IB:
            if (service1.svc_host.peer.NNTI_remote_process_t_u.ib.addr != service2.svc_host.peer.NNTI_remote_process_t_u.ib.addr) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.svc_host.peer.NNTI_remote_process_t_u.ib.port != service2.svc_host.peer.NNTI_remote_process_t_u.ib.port) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.svc_host.peer.NNTI_remote_process_t_u.ib.qpn  != service2.svc_host.peer.NNTI_remote_process_t_u.ib.qpn) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_GEMINI:
            if (service1.svc_host.peer.NNTI_remote_process_t_u.gni.addr    != service2.svc_host.peer.NNTI_remote_process_t_u.gni.addr) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.svc_host.peer.NNTI_remote_process_t_u.gni.port    != service2.svc_host.peer.NNTI_remote_process_t_u.gni.port) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.svc_host.peer.NNTI_remote_process_t_u.gni.inst_id != service2.svc_host.peer.NNTI_remote_process_t_u.gni.inst_id) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_MPI:
            if (service1.svc_host.peer.NNTI_remote_process_t_u.mpi.rank != service2.svc_host.peer.NNTI_remote_process_t_u.mpi.rank) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_LUC:
        case NNTI_TRANSPORT_LOCAL:
        case NNTI_TRANSPORT_NULL:
            break;
    }

    if (service1.req_addr.transport_id != service2.req_addr.transport_id) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    switch (service1.req_addr.buffer_owner.peer.transport_id) {
        case NNTI_TRANSPORT_PORTALS:
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.portals.nid != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.portals.nid) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.portals.pid != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.portals.pid) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_IB:
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.addr != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.addr) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.port != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.port) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.qpn  != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.ib.qpn) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_GEMINI:
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.addr    != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.addr) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.port    != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.port) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.gni.inst_id) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_MPI:
            if (service1.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank != service2.req_addr.buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_LUC:
        case NNTI_TRANSPORT_LOCAL:
        case NNTI_TRANSPORT_NULL:
            break;
    }

    if (service1.req_addr.buffer_addr.transport_id != service2.req_addr.buffer_addr.transport_id) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    switch (service1.req_addr.buffer_addr.transport_id) {
        case NNTI_TRANSPORT_PORTALS:
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id  != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.buffer_id) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.match_bits != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.match_bits) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.size       != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.portals.size) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_IB:
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.buf      != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.buf) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.key      != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.key) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.size     != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.size) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.ack_buf  != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.ack_buf) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.ack_size != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.ib.ack_size) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_GEMINI:
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.type              != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.type) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.buf               != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.buf) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1    != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword1) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2    != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.mem_hdl.qword2) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.size              != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.size) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr           != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_addr) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1 != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword1) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2 != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.gni.wc_mem_hdl.qword2) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_MPI:
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.rtr_tag  != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.rtr_tag) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.rts_tag  != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.rts_tag) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            if (service1.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.size     != service2.req_addr.buffer_addr.NNTI_remote_addr_t_u.mpi.size) {
                log_error(getservice_debug_level, "service compare FAILED");
                rc=NSSI_EINVAL;
                goto out;
            }
            break;
        case NNTI_TRANSPORT_LUC:
        case NNTI_TRANSPORT_LOCAL:
        case NNTI_TRANSPORT_NULL:
            break;
    }

    if (service1.req_addr.ops          != service2.req_addr.ops) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    if (service1.req_addr.payload_size != service2.req_addr.payload_size) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }
    if (service1.req_addr.payload      != service2.req_addr.payload) {
        log_error(getservice_debug_level, "service compare FAILED");
        rc=NSSI_EINVAL;
        goto out;
    }

out:
    return(rc);
}


/**
 * @brief Main code for data transfer client.
 *
 * @param args      The options for the experiment, set at the command-line
 * @param getservice_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
getservice_client_main (struct getservice_args &args, nssi_service &getservice_svc, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int i,j;
    double time;
    double start_time;
    std::ofstream result_stream;
    log_level debug_level = getservice_debug_level;
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

    log_debug(debug_level, "%d: Starting experiment loop", client_rank);

    /* loop over the experiments (count == num_trials, num_reqs == ops_per_trial) */
    for (i=0; i<args.num_trials; i++) {

        log_debug(debug_level, "%d: trial=%d, reqs=%d, len=%d\n",
                client_rank, i, args.num_reqs);

        MPI_Barrier(client_comm);
        start_time = Trios::GetTime();

        log_debug(debug_level, "Starting ");

        for (j=0; j<args.num_reqs; j++) {
            nssi_service duplicate_svc;

            log_debug(debug_level, "Getting duplicate service");
            rc=nssi_get_service((nssi_rpc_transport)args.transport, args.server_url.c_str(), args.timeout, &duplicate_svc);
            if ((rc != NSSI_OK) && (rc != NSSI_ETIMEDOUT)) {
                log_error(debug_level, "could not get svc description: %s", nssi_err_str(rc));
                break;
            }

            if (compare_services(getservice_svc, duplicate_svc) != NSSI_OK) {
                log_error(debug_level, "duplicate service does not match original");
                break;
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
}

/**
 * @}
 */
