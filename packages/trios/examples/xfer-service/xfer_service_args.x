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
 *   @file xfer_xdr.x
 *
 *   @brief Type definitions for a simple rpc service.
 *
 *   @author Ron Oldfield (raoldfi\@cs.sandia.gov).
 *   $Revision: 342 $.
 *   $Date: 2005-05-01 23:30:57 -0600 (Sun, 01 May 2005) $.
 *
 */

/**
 * @defgroup xfer_example  Nessie Data Transfer Example
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
 * @defgroup xfer_types  Nessie Example Types
 * @ingroup xfer_example
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
enum xfer_op {
    /** Opcode for writing the data through the function arguments. */
    XFER_WRITE_ENCODE_OP = 1,
    /** Opcode for the writing the data through the data channel. */
    XFER_WRITE_RDMA_OP,
    /**  Opcode for reading the data through the result structure. */
    XFER_READ_ENCODE_OP,
    /** Opcode for reading the data throught the data channel. */
    XFER_READ_RDMA_OP
};

/**
 * @brief Opcodes for the MPI transfer tests
 */
enum xfer_mpi_tag {
    /** Opcode for blocking send request */
    XFER_MPI_SEND_REQ_TAG = 1001,
    XFER_MPI_SEND_ACK_TAG,
    XFER_MPI_SEND_DATA_TAG,
    /** Opcode for non-blocking send request */
    XFER_MPI_ISEND_REQ_TAG,
    XFER_MPI_ISEND_ACK_TAG,
    XFER_MPI_ISEND_DATA_TAG,
    /** Opcode for blocking recv request */
    XFER_MPI_RECV_REQ_TAG,
    XFER_MPI_RECV_ACK_TAG,
    XFER_MPI_RECV_DATA_TAG,
    /** Opcode for non-blocking recv request */
    XFER_MPI_IRECV_REQ_TAG,
    XFER_MPI_IRECV_ACK_TAG,
    XFER_MPI_IRECV_DATA_TAG,
    /** Opcode for one-sided put request */
    XFER_MPI_PUT_REQ_TAG,
    XFER_MPI_PUT_ACK_TAG,
    XFER_MPI_PUT_DATA_TAG,
    /** Opcode for one-sided get request */
    XFER_MPI_GET_REQ_TAG,
    XFER_MPI_GET_ACK_TAG,
    XFER_MPI_GET_DATA_TAG,
    /** Opcode to quit the server */
    XFER_MPI_FINI_REQ
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
 * @brief Arguments for the first transfer operation, XFER_WRITE_ENCODE.
 *
 * The first transfer operation includes the array of \ref data_t
 * structures as an argument of the remote operation.  This will
 * cause the array to be sent to the server as part of the request.
 * The client encodes the data before sending it to the server, the server
 * decodes the data structure when it arrives.
 */
struct xfer_write_encode_args {
        /** The length of the data array. */
        int32_t len;
        /** A seed for initializing the array */
        uint32_t seed;
        /** A flag to perform validation */
        bool validate;
        /** The array of \ref data_t structures, including length. */
        data_array_t array;

};

/**
 * @brief Arguments for the second transfer operation, XFER_WRITE_RDMA.
 *
 * The second transfer operation only needs to send the length
 * of the data array.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to fetch.
 */
struct xfer_write_rdma_args {
        /** The length of the data array. */
        int32_t len;
        /** A seed for initializing the array */
        uint32_t seed;
        /** A flag to perform validation */
        bool validate;
};


/**
 * @brief Arguments for the third transfer operation, XFER_READ_ENCODE.
 *
 * The third transfer operation only needs to send the length
 * of the data array.  It uses the data argument in the nssi_call_rpc()
 * to identify the raw data for the server to fetch.
 */
struct xfer_read_encode_args {
        /** The length of the data array. */
        int32_t len;
        /** A value used to initialize the data */
        uint32_t seed;
        /** A flag to perform validation */
        bool validate;
};

/**
 * @brief Results for the third transfer operation, XFER_READ_ENCODE.
 *
 * The result of the xfer_read_encode operation includes the
 * array of \ref data_t structures.  If the size of array is large,
 * the address of array will be sent as part of the result and the
 * client will fetch array via RDMA.  In either case, the structures
 * are encoded by the server and decoded by the client.
 */
struct xfer_read_encode_res {
        /** The array of \ref data_t structures, including length. */
        data_array_t array;
};

/**
 * @brief Arguments for the fourth transfer operation, XFER_READ_RDMA.
 *
 * The xfer_read_rdma operation includes the array of \ref data_t
 * structures as an argument of the remote operation.  This will
 * cause the array to be sent to the server as part of the request.
 */
struct xfer_read_rdma_args {
        /** The length of the data array. */
        int32_t len;
        /** A seed for initializing the array */
        uint32_t seed;
        /** A flag to perform validation */
        bool validate;
};



/**
 * @}
 */
