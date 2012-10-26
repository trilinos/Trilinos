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
#include "Trios_nssi_debug.h"

#include "xfer_util.h"

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
 * @brief Send the data through the MPI_Send function
 *
 * @param len      Number of data structures sent.
 * @param source   MPI rank of the sending process.
 */
int process_mpi_send(
        const int len,
        const int seed,
        const int validate,
        const int source)
{
    int nbytes = len*sizeof(data_t);
    MPI_Request req;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log_level debug_level = xfer_debug_level;
    int rc = 0;
    MPI_Status status;

    log_debug(debug_level, "%d: starting process_mpi_send(len=%d, seed=%d, validate=%d)", rank, len, seed, validate);


    // Allocate space for the incoming buffer
    data_array_t array;
    array.data_array_t_len = len;
    array.data_array_t_val = new data_t[len];


    log_debug(debug_level, "%d: Server Posting IRECV for %d structs, %d bytes, tag=%d",
            rank, len, nbytes, XFER_MPI_SEND_DATA_TAG);

    // Post a receive for the buffer of ''len'' elements
    MPI_Irecv(array.data_array_t_val, nbytes, MPI_BYTE, source, XFER_MPI_SEND_DATA_TAG, MPI_COMM_WORLD, &req);

    log_debug(debug_level, "%d: Server Sending \"READY ACK\" to %d", rank, source);

    // Send an ack to the client that the server is ready
    rc = MPI_Send(&rc, 0, MPI_INT, source, XFER_MPI_SEND_ACK_TAG, MPI_COMM_WORLD);

    log_debug(debug_level, "%d: Server is waiting for data from %d, rc=%d", rank, source, rc);

    // Wait for the data to arrive
    rc = MPI_Wait(&req, &status);
    log_debug(debug_level, "%d: MPI_Wait complete: rc=%d, status.MPI_ERROR=%d, status.MPI_SOURCE=%d",
            rank, rc, status.MPI_ERROR, status.MPI_SOURCE);

    if (validate) {
        rc = xfer_validate_array(seed, &array);
        if (rc != 0)
            log_warn(debug_level, "Invalid array");
    }

    // free the memory
    delete [] array.data_array_t_val;

    return rc;
}


/**
 * @brief Send the data through the MPI_Send function
 *
 * @param len      Number of data structures sent.
 * @param source   MPI rank of the sending process.
 */
int process_mpi_isend(
        const int len,
        const int seed,
        const int validate,
        const int source)
{
    int nbytes = len*sizeof(data_t);

    log_level debug_level = xfer_debug_level;
    int rc=0;
    MPI_Status status;
    MPI_Request req;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log_debug(debug_level, "%d: starting process_mpi_isend(len=%d, seed=%d, validate=%d)", rank, len, seed, validate);


    // Allocate space for the incoming buffer
    data_array_t array;
    array.data_array_t_len = len;
    array.data_array_t_val = new data_t[len];

    log_debug(debug_level, "Posting IRECV");

    // Post a receive for the buffer of ''len'' elements
    MPI_Irecv(array.data_array_t_val, nbytes, MPI_BYTE, source, XFER_MPI_ISEND_DATA_TAG, MPI_COMM_WORLD, &req);

    log_debug(debug_level, "Sending \"READY ACK\"");

    // Send an ack to the client that the server is ready
    MPI_Send(&rc, 0, MPI_INT, source, XFER_MPI_ISEND_ACK_TAG, MPI_COMM_WORLD);

    log_debug(debug_level, "Waiting for data");

    // Wait for the data to arrive
    MPI_Wait(&req, &status);

    if (validate) {
        rc = xfer_validate_array(seed, &array);
        if (rc != 0)
            log_warn(debug_level, "Invalid array");
    }

    // Do we need to send a result?
    delete [] array.data_array_t_val;

    return rc;
}


/**
 * @brief Send the data through the MPI_Put function
 *
 * @param len      Number of data structures sent.
 * @param source   MPI rank of the sending process.
 */
int process_mpi_put(
        const int len,
        const int seed,
        const int validate,
        const int source)
{
    int nbytes = len*sizeof(data_t);

    log_level debug_level = xfer_debug_level;
    int rc = 0;
    MPI_Status status;
    MPI_Request req;
    MPI_Win win;
    MPI_Group comm_group, group;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log_debug(debug_level, "%d: starting process_mpi_put(len=%d, seed=%d, validate=%d)",
            rank, len, seed, validate);

    // Allocate space for the incoming buffer
    data_array_t array;
    array.data_array_t_len = len;
    array.data_array_t_val = new data_t[len];

    log_debug(debug_level, "%d: creating window", rank);


    // Collective call to create a window for the data being transferred (like registering the memory)
    MPI_Win_create(array.data_array_t_val, sizeof(data_t)*len,
            1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);


    MPI_Win_fence(0, win);
    log_debug(debug_level, "Waiting for data");
    MPI_Win_fence(0, win);


    log_debug(debug_level, "Data should be here... validate if necessary, val[0]=%d",array.data_array_t_val[0].int_val);

    if (validate) {
        rc = xfer_validate_array(seed, &array);
        if (rc == 0) {
            log_info(debug_level, "Validate Passed");
        }
        else {
            log_warn(debug_level, "Invalid array");
        }

        log_info(debug_level, "after validate, rc=%d", rc);
    }


    // clean up
    MPI_Win_free(&win);


    delete [] array.data_array_t_val;


    return rc;
}


/**
 * @brief Client reads the data using the MPI_Get function
 *
 * @param len      Number of data structures sent.
 * @param source   MPI rank of the sending process.
 */
int process_mpi_get(
        const int len,
        const int seed,
        const int validate,
        const int source)
{
    int nbytes = len*sizeof(data_t);

    log_level debug_level = xfer_debug_level;
    int rc = 0;
    MPI_Status status;
    MPI_Request req;
    MPI_Win win;
    MPI_Group comm_group, group;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    log_debug(debug_level, "%d: starting process_mpi_get(len=%d, seed=%d, validate=%d)",
            rank, len, seed, validate);

    // Allocate space for the incoming buffer
    data_array_t array;
    array.data_array_t_len = len;
    array.data_array_t_val = new data_t[len];

    log_debug(debug_level, "%d: creating window", rank);

    // initialize the array
    xfer_init_data_array(seed, &array);

    // Collective call to create a window for the data being transferred (like registering the memory)
    MPI_Win_create(array.data_array_t_val, sizeof(data_t)*len,
            1, MPI_INFO_NULL, MPI_COMM_WORLD, &win);


    MPI_Win_fence(0, win);
    log_debug(debug_level, "Waiting for client to get data");
    MPI_Win_fence(0, win);


    log_debug(debug_level, "Data should be transferred");


    // clean up
    MPI_Win_free(&win);


    delete [] array.data_array_t_val;


    return rc;
}


void make_progress(bool is_idle)
{
    log_debug(xfer_debug_level, "current_time(%llu) is_idle(%llu)", (uint64_t)trios_get_time_ms(), (uint64_t)is_idle);

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
int xfer_mpi_server_main(MPI_Comm server_comm)
{
    int rc = 0;

    log_level debug_level = xfer_debug_level;
    int server_rank;

    MPI_Comm_rank(server_comm, &server_rank);

    log_debug(debug_level, "%d: Starting server", server_rank);

    /* options that can be overriden by the command-line */
    std::string logfile("");

    bool done=false;
    MPI_Status status;

    // Server loop just waits for requests from any of the clients.  The
    // server isn't expecting any particular type of request.

    while (!done)  {
        int req_buf[3];

        log_debug(debug_level, "Waiting for request");

        // Receive the next request, the tag identifies the type of request
        MPI_Recv(&req_buf, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        int len = req_buf[0];
        int seed = req_buf[1];
        int validate = req_buf[2];

        log_debug(debug_level, "Server Received Request (source=%d, tag=%d, len=%d, seed=%d, validate=%d)",
                status.MPI_SOURCE, status.MPI_TAG, len, seed, validate);

        switch (status.MPI_TAG) {

        case XFER_MPI_SEND_REQ_TAG:
            rc = process_mpi_send(len, seed, validate, status.MPI_SOURCE);
            if (rc != 0)
                log_error(debug_level, "Error processing SEND request");
            break;

        case XFER_MPI_ISEND_REQ_TAG:
            log_debug(debug_level, "Server Calling ISEND");
            rc = process_mpi_isend(len, seed, validate, status.MPI_SOURCE);
            if (rc != 0)
                log_error(debug_level, "Error processing ISEND request");
            break;

        case XFER_MPI_PUT_REQ_TAG:
            log_debug(debug_level, "Server Calling PUT");
            rc = process_mpi_put(len, seed, validate, status.MPI_SOURCE);
            if (rc != 0)
                log_error(debug_level, "Error processing PUT request");
            break;

        case XFER_MPI_GET_REQ_TAG:
            log_debug(debug_level, "Server Calling PUT");
            rc = process_mpi_get(len, seed, validate, status.MPI_SOURCE);
            if (rc != 0)
                log_error(debug_level, "Error processing GET request");
            break;

        case XFER_MPI_FINI_REQ:
            log_debug(debug_level, "Server Received EXIT");
            done = true;
            break;

        default:
            log_error(debug_level, "tag=%d is not supported", status.MPI_TAG);
            return -1;
        }
    }

    log_debug(debug_level, "Exiting server");
    return rc;
}

/**
 * @}
 */
