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
 *   @file opreg_service_args.x
 *
 *   @brief Type definitions for an request opreg test.
 *
 *   @author Todd kordenbrock (thkorde\@sandia.gov).
 *
 *  Created on: Aug 22, 2011
 *
 */

/**
 * @defgroup opreg_example  Nessie Data Transfer Example
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
 * @defgroup opreg_types  Nessie Example Types
 * @ingroup opreg_example
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



enum opreg_test_type {
    OPREG_C_TEST = 0,
    OPREG_CPP_TEST
};

/**
 * @brief Opcodes for the types of transfer operations.
 */
enum opreg_op {
    /** Opcode for sending a request. */
    OPREG_REQUEST_OP = 1
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

struct opreg_args {
        /* test data */
        data_t data;
        /** 32-bit checksum of the test data */
        uint64_t chksum;
};

/**
 * @}
 */
