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
 * NssiWaitTest.c
 *
 *  Created on: August 14, 2014
 *      Author: thkorde
 */

#include "Trios_nssi_client.h"
#include "Trios_nssi_server.h"
#include "Trios_nssi_fprint_types.h"

#include <unistd.h>
#include <memory.h>

#include <signal.h>

#include <mpi.h>

#include "Trios_logger.h"
#include "Trios_timer.h"
#include "Trios_trace.h"

#include "nssi_opcodes.h"
#include "nssi_service_args.h"

#include <sstream>
#include <ostream>


log_level nssiwait_debug_level = LOG_UNDEFINED;


#define REQ_COUNT 5


int wait_test(const nssi_service *svc)
{
    static const char *empty="";

    int rc  = NSSI_OK;
    int rc2 = NSSI_OK;

    int timeout=5000;
    nssi_request req[REQ_COUNT];

    nssi_trace_reset_args args;

    memset(&args, 0, sizeof(nssi_trace_reset_args));

    /* initialize the args */
    args.fname  = (char *)empty;
    args.ftype  = 0;
    args.enable = (char *)empty;

    log_debug(nssiwait_debug_level, "Calling nssi_trace_reset(%s, %d, %s)\n",
    		args.fname, args.ftype, args.enable);

    /* send the request */
    for (int i=0;i<REQ_COUNT;i++) {
		rc = nssi_call_rpc(svc, NSSI_OP_TRACE_RESET, &args, NULL, 0, NULL, &req[i]);
		if (rc != NSSI_OK) {
			log_error(nssiwait_debug_level, "unable to call remote method: %s",
					nssi_err_str(rc));
			return rc;
		}
    }

    /* wait for completion */
    for (int i=0;i<REQ_COUNT;i++) {
		log_debug(nssiwait_debug_level, "calling nssi_timedwait(req[%d])", i);
		rc2 = nssi_timedwait(&req[i], timeout, &rc);
		if (rc2 != NSSI_OK) {
			log_error(nssiwait_debug_level, "failed waiting for request %lu: %s",
					req[i].id,
					nssi_err_str(rc2));
			return rc2;
		}

		if (rc != NSSI_OK) {
			log_error(nssiwait_debug_level, "remote method failed: %s",
					nssi_err_str(rc));
			return rc;
		}
    }

    return rc;
}


int waitany_test(const nssi_service *svc)
{
    static const char *empty="";

    int rc  = NSSI_OK;
    int rc2 = NSSI_OK;

    int timeout=5000;
    nssi_request req[REQ_COUNT];
    int which=-1;

    nssi_trace_reset_args args;

    memset(&args, 0, sizeof(nssi_trace_reset_args));

    /* initialize the args */
    args.fname  = (char *)empty;
    args.ftype  = 0;
    args.enable = (char *)empty;

    log_debug(nssiwait_debug_level, "Calling nssi_trace_reset(%s, %d, %s)\n",
    		args.fname, args.ftype, args.enable);

    /* send the request */
    for (int i=0;i<REQ_COUNT;i++) {
		rc = nssi_call_rpc(svc, NSSI_OP_TRACE_RESET, &args, NULL, 0, NULL, &req[i]);
		if (rc != NSSI_OK) {
			log_error(nssiwait_debug_level, "unable to call remote method: %s",
					nssi_err_str(rc));
			return rc;
		}
    }

    /* wait for completion */
    for (int i=0;i<REQ_COUNT;) {
		log_debug(nssiwait_debug_level, "calling nssi_waitany()");
		rc2 = nssi_waitany(&req[0], REQ_COUNT, timeout, &which, &rc);
		if (rc2 != NSSI_OK) {
			log_error(nssiwait_debug_level, "failed waiting for request %lu: %s",
					req[which].id,
					nssi_err_str(rc2));
			return rc2;
		}
		if (rc == NSSI_ETIMEDOUT) {
			continue;
		}
		if (rc != NSSI_OK) {
			log_error(nssiwait_debug_level, "remote method failed: %s",
					nssi_err_str(rc));
			return rc;
		}

		i++;

		log_debug(nssiwait_debug_level, "nssi_waitany() which=%d", which);
    }

    return rc;
}

int waitall_test(const nssi_service *svc)
{
    static const char *empty="";

    int rc  = NSSI_OK;

    int timeout=5000;
    nssi_request req[REQ_COUNT];

    nssi_trace_reset_args args;

    memset(&args, 0, sizeof(nssi_trace_reset_args));

    /* initialize the args */
    args.fname  = (char *)empty;
    args.ftype  = 0;
    args.enable = (char *)empty;

    log_debug(nssiwait_debug_level, "Calling nssi_trace_reset(%s, %d, %s)\n",
    		args.fname, args.ftype, args.enable);

    /* send the request */
    for (int i=0;i<REQ_COUNT;i++) {
		rc = nssi_call_rpc(svc, NSSI_OP_TRACE_RESET, &args, NULL, 0, NULL, &req[i]);
		if (rc != NSSI_OK) {
			log_error(nssiwait_debug_level, "unable to call remote method: %s",
					nssi_err_str(rc));
			return rc;
		}
    }

    /* wait for completion */
	log_debug(nssiwait_debug_level, "calling nssi_waitall()");
	rc = nssi_waitall(&req[0], REQ_COUNT, timeout);
	if (rc != NSSI_OK) {
		log_error(nssiwait_debug_level, "failed waiting for all requests: %s",
				nssi_err_str(rc));
		return rc;
	}
	log_debug(nssiwait_debug_level, "nssi_waitall() complete");

    return rc;
}



int main(int argc, char *argv[])
{
	int rc=0;
	int success=0;
    int nprocs, rank;

    char logname[1024];
    char url[NSSI_URL_LEN];

    nssi_service svc;


    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//    struct sigaction sigact;
//    sigact.sa_handler=SIG_DFL;
//    sigaction ( 6, &sigact, NULL);
//    sigaction (11, &sigact, NULL);

    sprintf(logname, "nssiwait.%03d.log", rank);
    logger_init(LOG_ERROR, NULL);

    // common init
    nssi_rpc_init(NSSI_DEFAULT_TRANSPORT, NSSI_DEFAULT_ENCODE, NULL);

    if (rank==0) {
        // service init
        nssi_service_init(NSSI_DEFAULT_TRANSPORT, NSSI_SHORT_REQUEST_SIZE, &svc);
        nssi_get_url(NSSI_DEFAULT_TRANSPORT, &url[0], NSSI_URL_LEN);

        MPI_Bcast(&url[0], NSSI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

        svc.max_reqs = -1;
        rc = nssi_service_start(&svc);
        if (rc != NSSI_OK) {
            log_info(nssiwait_debug_level, "exited selfsend_svc: %s",
                    nssi_err_str(rc));
        }

        /* finalize the service */
        nssi_service_fini(&svc);

    } else {
        // client init
        nssi_init(NSSI_DEFAULT_TRANSPORT);

        MPI_Bcast(&url[0], NSSI_URL_LEN, MPI_CHAR, 0, MPI_COMM_WORLD);

        for (int i=0; i < 3; i++) {
            log_debug(nssiwait_debug_level, "Try to connect to server: attempt #%d", i);
            rc=nssi_get_service(NSSI_DEFAULT_TRANSPORT, url, 5000, &svc);
            if (rc == NSSI_OK)
                break;
            else if (rc != NSSI_ETIMEDOUT) {
                log_error(nssiwait_debug_level, "could not get svc description: %s",
                        nssi_err_str(rc));
                break;
            }
        }

        wait_test(&svc);
        waitany_test(&svc);
        waitall_test(&svc);

        // shutdown the service
        rc = nssi_kill(&svc, 0, 5000);
        if (rc != NSSI_OK) {
            log_error(nssiwait_debug_level, "Error in nssi_kill");
        }

        // finalize the client
        nssi_fini(NSSI_DEFAULT_TRANSPORT);

        success=1;
    }

    // finalize the NSSI library
    rc = nssi_rpc_fini(NSSI_DEFAULT_TRANSPORT);
    if (rc != NSSI_OK) {
        log_error(nssiwait_debug_level, "Error in nssi_rpc_fini");
    }

    MPI_Bcast(&success, 1, MPI_INT, 1, MPI_COMM_WORLD);

    MPI_Finalize();

    logger_fini();

    if (rank == 1) {
		if (success==TRUE)
			fprintf(stdout, "\nEnd Result: TEST PASSED\n");
		else
			fprintf(stdout, "\nEnd Result: TEST FAILED\n");
    }

    return (success==TRUE ? 0 : 1 );
}
