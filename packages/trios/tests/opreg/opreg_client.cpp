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
 * opreg_client.cpp
 *
 *  Created on: Nov 14, 2011
 *      Author: thkorde
 */

/**
 * @defgroup opreg_client Data Transfer Client
 *
 * @ingroup opreg_example
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_logger.h"
#include "Trios_timer.h"

#include "Teuchos_CommandLineProcessor.hpp"

#include "opreg_test.h"
#include "opreg_service_args.h"


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



log_level client_debug_level = LOG_UNDEFINED;

extern int print_args(
        std::ostream &out,
        const struct opreg_args &args,
        const char *prefix);

extern "C" {
uint64_t calc_checksum (const char * buf, const uint64_t size);
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
int opreg_request(
    const nssi_service *svc,
    nssi_request *req)
{
    int rc = NSSI_OK;
    log_level debug_level = client_debug_level;
    opreg_args args;

    log_debug(debug_level, "Calling RPC for OPREG_REQUEST_OP");

    /* initialize the arguments */
    memset(&args, 0, sizeof(opreg_args));
    args.data.int_val    = 10;
    args.data.float_val  = 10.10;
    args.data.double_val = 10.10;
    args.chksum = calc_checksum((const char*)&args.data, sizeof(args.data));

    log_debug(debug_level, "args.data.int_val(%d) args.data.float_val(%f) args.data.double_val(%f) args.chksum(%lu)",
            args.data.int_val, args.data.float_val, args.data.double_val, args.chksum);

    /* call the remote methods */
    rc = nssi_call_rpc(svc, OPREG_REQUEST_OP, &args, NULL, 0, NULL, req);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "unable to call OPREG_EMPTY_REQUEST_OP: %s",
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
int opreg_request_blk(
    const nssi_service *svc)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;
    nssi_request req;

    /* call the async function */
    rc = opreg_request(svc, &req);
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
 * @param opreg_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
opreg_client_main (struct opreg_cmdline_args &args, nssi_service &opreg_svc, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int client_rank, client_size;


    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);


    /* register the XDR encoding functions */
    NSSI_REGISTER_CLIENT_STUB(OPREG_REQUEST_OP, opreg_args, void, void);

    rc = opreg_request_blk(&opreg_svc);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "could not transfer data: %s",
                nssi_err_str(rc));
        goto abort;
    }

    if (client_rank==0)
        cout << "experiment complete" << endl;


    return 0;

abort:
    exit(2);
}

/**
 * @}
 */
