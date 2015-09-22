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
 * send_client.cpp
 *
 *  Created on: Sep 25, 2014
 *      Author: thkorde
 */

/**
 * @defgroup send_client Data Transfer Client
 *
 * @ingroup send_test
 *
 * @{
 */

#include "Trios_config.h"
#include "Trios_nssi_client.h"
#include "Trios_nnti.h"
#include "Trios_logger.h"
#include "Trios_timer.h"

#include "Teuchos_CommandLineProcessor.hpp"

#include "nssi_send_test.h"
#include "nssi_send_service_args.h"


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
        const struct send_args &args,
        const char *prefix);
uint64_t calc_checksum (const char * buf, const uint64_t size);


#define MAX_COUNT 100


#define REQUEST_FLAG_APP_PINNED_SHORT_REQUEST     (1 << 0)
#define REQUEST_FLAG_APP_PINNED_LONG_ARGS         (1 << 1)
#define REQUEST_FLAG_APP_PINNED_BULK_DATA         (1 << 2)
#define REQUEST_FLAG_APP_PINNED_SHORT_RESULT      (1 << 3)

#define REQUEST_FLAG_SYNC                         (1 << 4)
#define REQUEST_FLAG_RESPONSELESS                 (1 << 5)

#define REQUEST_APP_PINNED_FLAG_COUNT 4


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
int send_no_rdma(
    const nssi_service *svc,
    uint32_t            flags,
    NNTI_buffer_t      *short_request_hdl,
    NNTI_buffer_t      *long_args_hdl,
    NNTI_buffer_t      *bulk_data_hdl,
    NNTI_buffer_t      *short_result_hdl,
    int32_t             count,
    send_test_args     *args,
    void               *data,
    uint32_t            data_size,
    send_test_res      *res,
    nssi_request       *req)
{
    int rc = NSSI_OK;
    log_level debug_level = client_debug_level;

    log_debug(debug_level, "Calling RPC for SEND_TEST_NO_RDMA_OP");

    rc = nssi_create_request(
            svc,
            SEND_TEST_NO_RDMA_OP,
            args,
            data,
            data_size,
            res,
            req);

    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_REQUEST) {
        if (short_request_hdl == NULL) {
            MPI_Abort(MPI_COMM_WORLD, -5);
        }
        req->app_pinned_short_request=TRUE;
        req->short_request_hdl = short_request_hdl;
    }
    if (flags & REQUEST_FLAG_APP_PINNED_LONG_ARGS) {
        if (long_args_hdl == NULL) {
            MPI_Abort(MPI_COMM_WORLD, -5);
        }
        req->app_pinned_long_args=TRUE;
        req->long_args_hdl = long_args_hdl;
    }
    if (flags & REQUEST_FLAG_APP_PINNED_BULK_DATA) {
        if (bulk_data_hdl == NULL) {
            MPI_Abort(MPI_COMM_WORLD, -5);
        }
        req->app_pinned_bulk_data=TRUE;
        req->bulk_data_hdl = bulk_data_hdl;
    }
    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_RESULT) {
        if (short_result_hdl == NULL) {
            MPI_Abort(MPI_COMM_WORLD, -5);
        }
        req->app_pinned_short_result=TRUE;
        req->short_result_hdl = short_result_hdl;
    }

    if (flags & REQUEST_FLAG_SYNC) {
        req->is_sync=TRUE;
    }
    if (flags & REQUEST_FLAG_RESPONSELESS) {
        req->is_responseless=TRUE;
        req->result          =NULL;
        req->short_result_hdl=NULL;
    }

    /* initialize the arguments */
    memset(args, 0, sizeof(send_test_args));

    args->count                  = count;
    args->array.data_array_t_len = count;
    args->array.data_array_t_val = (data_t *)calloc(count, sizeof(data_t));

    for (int i=0;i<count;i++) {
        args->array.data_array_t_val->int_val    = i + 10;
        args->array.data_array_t_val->float_val  = i + 10.10;
        args->array.data_array_t_val->double_val = i + 10.10;
    }
    args->args_chksum = calc_checksum((const char*)args->array.data_array_t_val, sizeof(count*sizeof(data_t)));

    args->bulk_chksum = 0;

    /* call the remote methods */
    rc = nssi_send_request(req);
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
int send_no_rdma_blk(
        const nssi_service *svc,
        uint32_t            flags,
        NNTI_buffer_t      *short_request_hdl,
        NNTI_buffer_t      *long_args_hdl,
        NNTI_buffer_t      *bulk_data_hdl,
        NNTI_buffer_t      *short_result_hdl,
        int32_t             count,
        send_test_args     *args,
        void               *data,
        uint32_t            data_size,
        send_test_res      *res)
{
    int rc = NSSI_OK;
    int rc2 = NSSI_OK;

    nssi_request req;

    if ((flags & REQUEST_FLAG_SYNC) || (flags & REQUEST_FLAG_RESPONSELESS)) {
        log_error(client_debug_level, "This is the nssi_send_request() + nssi_wait() test - cannot be combined with REQUEST_FLAG_SYNC or REQUEST_FLAG_RESPONSELESS");
        rc = NSSI_EINVAL;
        goto out;
    }

    /* call the async function */
    rc = send_no_rdma(
            svc,
            flags,
            short_request_hdl,
            long_args_hdl,
            bulk_data_hdl,
            short_result_hdl,
            count,
            args,
            data,
            data_size,
            res,
            &req);
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

out:
    return rc;

}

int execute_send_no_rdma(
        const nssi_service *svc,
        uint32_t            flags,
        NNTI_buffer_t      *short_request_hdl,
        NNTI_buffer_t      *long_args_hdl,
        NNTI_buffer_t      *bulk_data_hdl,
        NNTI_buffer_t      *short_result_hdl,
        int32_t             count,
        void               *data,
        uint32_t            data_size)
{
    int rc = NSSI_OK;

    send_test_args  args;
    send_test_res   res;

    if ((flags & REQUEST_FLAG_SYNC) || (flags & REQUEST_FLAG_RESPONSELESS)) {
        nssi_request req;

        rc = send_no_rdma(
                svc,
                flags,
                short_request_hdl,
                long_args_hdl,
                bulk_data_hdl,
                short_result_hdl,
                count,
                &args,
                data,
                data_size,
                &res,
                &req);
        if (rc != NSSI_OK) {
            log_error(client_debug_level, "unable to call async method: %s",
                    nssi_err_str(rc));
            return rc;
        }
    } else {
        rc = send_no_rdma_blk(
                svc,
                flags,
                short_request_hdl,
                long_args_hdl,
                bulk_data_hdl,
                short_result_hdl,
                count,
                &args,
                data,
                data_size,
                &res);
        if (rc != NSSI_OK) {
            log_error(client_debug_level, "unable to call async method: %s",
                    nssi_err_str(rc));
            return rc;
        }
    }

    if (!(flags & REQUEST_FLAG_RESPONSELESS)) {
        if (res.args_chksum != args.args_chksum) {
            log_error(client_debug_level, "server returned a different checksum: args.args_chksum(%llu) != res.args_chksum(%llu)", args.args_chksum, res.args_chksum);
            rc=NSSI_EINVAL;
        }
        if (res.bulk_chksum != args.bulk_chksum) {
            log_error(client_debug_level, "server returned a different checksum: args.bulk_chksum(%llu) != res.bulk_chksum(%llu)", args.bulk_chksum, res.bulk_chksum);
            rc=NSSI_EINVAL;
        }
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


extern NNTI_transport_t transports[NSSI_RPC_COUNT];

int setup_handles(
        const nssi_service *svc,
        uint32_t            flags,
        NNTI_buffer_t      *short_request_hdl,
        NNTI_buffer_t      *long_args_hdl,
        NNTI_buffer_t      *bulk_data_hdl,
        void               *data,
        uint32_t            data_size,
        NNTI_buffer_t      *short_result_hdl)
{
    int rc=NNTI_OK;

    int short_request_len = NNTI_BUFFER_SIZE(&svc->req_addr);
    int long_args_len     = 5120;


    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_REQUEST) {
        rc=NNTI_alloc(
                &transports[svc->transport_id],
                short_request_len,
                1,
                NNTI_SEND_SRC,
                short_request_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed registering short request: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_LONG_ARGS) {
        rc=NNTI_alloc(
                &transports[svc->transport_id],
                long_args_len,
                1,
                NNTI_GET_SRC,
                long_args_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed registering long args: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_BULK_DATA) {
        rc=NNTI_register_memory(
                &transports[svc->transport_id],
                (char *)data,
                data_size,
                1,
                (NNTI_buf_ops_t)(NNTI_GET_SRC|NNTI_PUT_DST),
                bulk_data_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed registering data: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_RESULT) {
        rc=NNTI_alloc(
                &transports[svc->transport_id],
                NSSI_SHORT_RESULT_SIZE,
                1,
                NNTI_RECV_DST,
                short_result_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed registering short result: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }

cleanup:
    return rc;
}


int teardown_handles(
        uint32_t       flags,
        NNTI_buffer_t *short_request_hdl,
        NNTI_buffer_t *long_args_hdl,
        NNTI_buffer_t *bulk_data_hdl,
        NNTI_buffer_t *short_result_hdl)
{
    int rc=NNTI_OK;

    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_REQUEST) {
        rc=NNTI_free(
                short_request_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed freeing short request: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_LONG_ARGS) {
        rc=NNTI_free(
                long_args_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed freeing long args: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_BULK_DATA) {
        rc=NNTI_unregister_memory(
                bulk_data_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed unregistering data: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }
    if (flags & REQUEST_FLAG_APP_PINNED_SHORT_RESULT) {
        rc=NNTI_free(
                short_result_hdl);
        if (rc != NNTI_OK) {
            log_error(client_debug_level, "failed freeing short result: %s",
                    nnti_err_str(rc));
            goto cleanup;
        }
    }

cleanup:
    return rc;
}


/**
 * @brief Main code for data transfer client.
 *
 * @param args      The options for the experiment, set at the command-line
 * @param send_svc  The nssi_service descriptor for the remote service (already connected)
 * @param comm      The communicator for the client application
 */
int
send_client_main (struct send_cmdline_args &args, nssi_service &send_svc, MPI_Comm client_comm)
{
    using namespace std;

    int rc;
    int client_rank, client_size;

    uint32_t      flags=0;

    NNTI_buffer_t short_request_hdl;
    NNTI_buffer_t long_args_hdl;
    NNTI_buffer_t bulk_data_hdl;
    NNTI_buffer_t short_result_hdl;

    int32_t       count=1;

    void         *data=NULL;
    uint32_t      data_size=0;


    MPI_Comm_rank(client_comm, &client_rank);
    MPI_Comm_size(client_comm, &client_size);


    /* register the XDR encoding functions */
    NSSI_REGISTER_CLIENT_STUB(SEND_TEST_NO_RDMA_OP, send_test_args, void, send_test_res);
    NSSI_REGISTER_CLIENT_STUB(SEND_TEST_SRV_GET_OP, send_test_args, void, send_test_res);
    NSSI_REGISTER_CLIENT_STUB(SEND_TEST_SRV_PUT_OP, send_test_args, void, send_test_res);

    data_size=MAX_COUNT*sizeof(data_t);
    data=malloc(data_size);



    count=1;


    /* First send a traditional RPC request.
     * NSSI is responsible for allocating and pinning all memory regions.
     * nssi_wait() is explicitly called to get the result.
     */
    flags=0;
    rc = execute_send_no_rdma(
            &send_svc,
            flags,
            &short_request_hdl,
            &long_args_hdl,
            &bulk_data_hdl,
            &short_result_hdl,
            count,
            data,
            data_size);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "could not transfer data: %s",
                nssi_err_str(rc));
        goto abort;
    }

    /* Next send a traditional RPC request but skip the nssi_wait().
     * NSSI is responsible for allocating and pinning all memory regions.
     * Set req->is_sync which causes an implicit nssi_wait() at the end of nssi_send_request().
     */
    flags = REQUEST_FLAG_SYNC;
    rc = execute_send_no_rdma(
            &send_svc,
            flags,
            &short_request_hdl,
            &long_args_hdl,
            &bulk_data_hdl,
            &short_result_hdl,
            count,
            data,
            data_size);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "could not transfer data: %s",
                nssi_err_str(rc));
        goto abort;
    }

    /* Next send a traditional RPC request but skip the nssi_wait().
     * NSSI is responsible for allocating and pinning all memory regions.
     * Set req->is_responseless which causes the server to skip sending the result.
     */
    flags = REQUEST_FLAG_RESPONSELESS;
    rc = execute_send_no_rdma(
            &send_svc,
            flags,
            &short_request_hdl,
            &long_args_hdl,
            &bulk_data_hdl,
            &short_result_hdl,
            count,
            NULL,
            0);
    if (rc != NSSI_OK) {
        log_error(client_debug_level, "could not transfer data: %s",
                nssi_err_str(rc));
        goto abort;
    }


    for (int i=0;i<REQUEST_APP_PINNED_FLAG_COUNT;i++) {
        flags=(1 << i);
        rc=setup_handles(
                &send_svc,
                flags,
                &short_request_hdl,
                &long_args_hdl,
                &bulk_data_hdl,
                data,
                data_size,
                &short_result_hdl);

        /* First send a traditional RPC request.
         * NSSI is responsible for allocating and pinning all memory regions.
         * nssi_wait() is explicitly called to get the result.
         */
        flags=(1 << i);
        rc = execute_send_no_rdma(
                &send_svc,
                flags,
                &short_request_hdl,
                &long_args_hdl,
                &bulk_data_hdl,
                &short_result_hdl,
                count,
                data,
                data_size);
        if (rc != NSSI_OK) {
            log_error(client_debug_level, "could not transfer data: %s",
                    nssi_err_str(rc));
            goto abort;
        }

        /* Next send a traditional RPC request but skip the nssi_wait().
         * NSSI is responsible for allocating and pinning all memory regions.
         * Set req->is_sync which causes an implicit nssi_wait() at the end of nssi_send_request().
         */
        flags=(1 << i) | REQUEST_FLAG_SYNC;
        rc = execute_send_no_rdma(
                &send_svc,
                flags,
                &short_request_hdl,
                &long_args_hdl,
                &bulk_data_hdl,
                &short_result_hdl,
                count,
                data,
                data_size);
        if (rc != NSSI_OK) {
            log_error(client_debug_level, "could not transfer data: %s",
                    nssi_err_str(rc));
            goto abort;
        }


        if (((1 << i)) == REQUEST_FLAG_APP_PINNED_SHORT_REQUEST) {
            /* Next send a traditional RPC request but skip the nssi_wait().
             * NSSI is responsible for allocating and pinning all memory regions.
             * Set req->is_responseless which causes the server to skip sending the result.
             */
            flags=(1 << i) | REQUEST_FLAG_RESPONSELESS;
            rc = execute_send_no_rdma(
                    &send_svc,
                    flags,
                    &short_request_hdl,
                    &long_args_hdl,
                    &bulk_data_hdl,
                    &short_result_hdl,
                    count,
                    NULL,
                    0);
            if (rc != NSSI_OK) {
                log_error(client_debug_level, "could not transfer data: %s",
                        nssi_err_str(rc));
                goto abort;
            }
        }

        flags=(1 << i);
        rc=teardown_handles(
                flags,
                &short_request_hdl,
                &long_args_hdl,
                &bulk_data_hdl,
                &short_result_hdl);
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
