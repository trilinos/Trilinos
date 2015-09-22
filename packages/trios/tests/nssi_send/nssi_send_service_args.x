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
 *   @file send_service_args.x
 *
 *   @brief Type definitions for the send test.
 *
 *   @author Todd kordenbrock (thkorde\@sandia.gov).
 *
 *  Created on: Sep 25, 2014
 *
 */

/**
 * @defgroup send_test  Exercise NSSI Send Options
 *
 */

/**
 * @defgroup send_test_types  NSSI Send Types
 * @ingroup send_test
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



enum send_test_type {
    SEND_TEST_NO_RDMA = 0,
    SEND_TEST_SRV_GET,
    SEND_TEST_SRV_PUT
};

/**
 * @brief Opcodes for the types of transfer operations.
 */
enum send_test_op {
    /** Opcode for a request with no server-side RDMA of bulk data. */
    SEND_TEST_NO_RDMA_OP = 1,
    
    /** Opcode for a request with a server-side GET of bulk data. */
    SEND_TEST_SRV_GET_OP = 2,

    /** Opcode for a request with a server-side PUT of bulk data. */
    SEND_TEST_SRV_PUT_OP = 3
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



struct send_test_args {
    /** Number of elements in the data array */
    int32_t      count;
    /** The array of \ref data_t structures, including length. */
    data_array_t array;
    /** 64-bit checksum of the data array sent in the args */
    uint64_t     args_chksum;
    /** 64-bit checksum of the data array transfered by RDMA */
    uint64_t     bulk_chksum;
};

struct send_test_res {
    /** 64-bit checksum of the data array sent in the args */
    uint64_t     args_chksum;
    /** 64-bit checksum of the data array transfered by RDMA */
    uint64_t     bulk_chksum;
};

/**
 * @}
 */
