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
/* -------------------------------------------------------------------------- */
/**
 *   @file netperf_xdr.x
 *
 *   @brief Type definitions for a simple rpc service.
 *
 *   @author Ron Oldfield (raoldfi\@cs.sandia.gov).
 *   $Revision: 342 $.
 *   $Date: 2005-05-01 23:30:57 -0600 (Sun, 01 May 2005) $.
 *
 */

/**
 * @defgroup netperf_example  Nessie Data Transfer Example
 *
 * The data-transfer example demonstrates a simple client
 * and server that transfer an array of 16-byte \ref data_t
 * data structures from a parallel application to a set of
 * servers.  We implemented two variations:
 * one that transfers the array of data structures
 * with the request, and a second method that has each server
 * pull the data using the \ref nssi_get_data() function.  Although
 * this example is fairly simple, it makes a decent benchmark code
 * to evaluate overheads of the Nessie transfer protocols and encoding
 * schemes.
 *
*/

/**
 * @defgroup netperf_types  Nessie Example Types
 * @ingroup netperf_example
 *
 * @{
 */

/* Extra stuff to put at the beginning of the header file */
#ifdef RPC_HDR
%#include "Trios_xdr.h"
#endif

/* Extra stuff to put at the beginning of the C file. */
#ifdef RPC_XDR
%#include "Trios_xdr.h"
#endif



/**
 * @brief Opcodes for the types of transfer operations.
 */
enum netperf_op {
    /** Opcode for client to tell server experiment parameters. */
    NETPERF_EXPERIMENT_SETUP_OP=1,
    /** Opcode for client to tell server to conclude experiment. */
    NETPERF_EXPERIMENT_TEARDOWN_OP,
    /** Opcode for the writing the data through the data channel. */
    NETPERF_WRITE_RDMA_OP,
    /** Opcode for reading the data throught the data channel. */
    NETPERF_READ_RDMA_OP
};

/**
 * @brief Opcodes for the MPI transfer tests
 */
enum netperf_mpi_tag {
    /** Opcode for blocking send request */
    NETPERF_MPI_SEND_REQ_TAG = 1001,
    NETPERF_MPI_SEND_ACK_TAG,
    NETPERF_MPI_SEND_DATA_TAG,
    /** Opcode for non-blocking send request */
    NETPERF_MPI_ISEND_REQ_TAG,
    NETPERF_MPI_ISEND_ACK_TAG,
    NETPERF_MPI_ISEND_DATA_TAG,
    /** Opcode for blocking recv request */
    NETPERF_MPI_RECV_REQ_TAG,
    NETPERF_MPI_RECV_ACK_TAG,
    NETPERF_MPI_RECV_DATA_TAG,
    /** Opcode for non-blocking recv request */
    NETPERF_MPI_IRECV_REQ_TAG,
    NETPERF_MPI_IRECV_ACK_TAG,
    NETPERF_MPI_IRECV_DATA_TAG,
    /** Opcode for one-sided put request */
    NETPERF_MPI_PUT_REQ_TAG,
    NETPERF_MPI_PUT_ACK_TAG,
    NETPERF_MPI_PUT_DATA_TAG,
    /** Opcode for one-sided get request */
    NETPERF_MPI_GET_REQ_TAG,
    NETPERF_MPI_GET_ACK_TAG,
    NETPERF_MPI_GET_DATA_TAG,
    /** Opcode to quit the server */
    NETPERF_MPI_FINI_REQ
};

/**
 * @brief A 16-byte structure that contains an int, float, and double.
 *
 * This structure contains an int, float, and double as an example
 * of a complex structure with multiple types.  This will exercise the
 * encoding/decoding features of Nessie.
 */
struct data_t {
    /** An integer value. */
    uint32_t int_val;
    /** A floating point value. */
    float float_val;
    /** A double value. */
    double double_val;
};

/**
 * @brief Array of 16-byte structures that we can send with the request.
 *
 * Rpcgen will use this definition to define encoding functions to
 * encode and decode an array of \ref data_t structures.  We will
 * use these functions when sending the array with the request.
 */
typedef data_t data_array_t<>;

/**
 * @brief Arguments for NETPERF_EXPERIMENT_SETUP.
 *
 * Communicate experiment parameters to server.
 */
struct netperf_experiment_setup_args {
    /** The number of requests per trial. The server will create and
     * initialize this many data arrays. */
    uint32_t num_reqs;
    /** The length of each data array. */
    int32_t array_len;
    /** A flag to perform validation.  Determines if server will
     * create validation arrays. */
    bool validate;
};

/**
 * @brief Arguments for the first transfer operation, NETPERF_WRITE_RDMA.
 *
 * The first transfer operation only needs to send the index of the
 * data array to transfer.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to get.
 */
struct netperf_write_rdma_args {
        /** Index the array */
        uint32_t buffer_index;
};

/**
 * @brief Arguments for the second transfer operation, NETPERF_READ_RDMA.
 *
 * The second transfer operation only needs to send the index of the
 * data array to transfer.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to put.
 */
struct netperf_read_rdma_args {
        /** Index the array */
        uint32_t buffer_index;
};

/**
 * @}
 */
