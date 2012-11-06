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
/*-------------------------------------------------------------------------*/
/**  @file rpc_xdr.cc
 *
 *   @brief  XDR-encoding support for the
 *           \ref rpc_client_api "RPC Client API".
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 481 $
 *   $Date: 2005-11-03 16:21:36 -0700 (Thu, 03 Nov 2005) $
 *
 */

#include "Trios_config.h"
#include "Trios_nssi_types.h"

#include <map>

#include <stdio.h>
#ifdef HAVE_TRIOS_MALLOC_H
#include <malloc.h>
#endif
#include <signal.h>
#include <stdlib.h>  /* find malloc()/free() */
#include <string.h> /* find memset() */

#include "Trios_nssi_xdr.h"
#include "Trios_logger.h"
#include "Trios_threads.h"
#include "Trios_nssi_debug.h"

#include "nssi_opcodes.h"
#include "nssi_service_args.h"


/**
 * @brief XDR encodings for NSSI remote operations.
 *
 * The <tt>\ref xdr_encodings</tt> structure contains the
 * xdr-encoding function pointers required for a remote
 * NSSI function.
 */
struct xdr_encodings {
    /** @brief A function pointer to encode/decode the arguments. */
    xdrproc_t encode_args;
    /** @brief A function pointer to encode/decode the data. */
    xdrproc_t encode_data;
    /** @brief A function pointer to encode/decode the result. */
    xdrproc_t encode_result;

    /** Default constructor */
    xdr_encodings() {}

    /** Constructor */
    xdr_encodings(xdrproc_t args, xdrproc_t data, xdrproc_t result) :
        encode_args(args), encode_data(data), encode_result(result) {}

    /** Copy constructor */
    xdr_encodings(const struct xdr_encodings &encodings)
    {
        this->encode_args = encodings.encode_args;
        this->encode_data = encodings.encode_data;
        this->encode_result = encodings.encode_result;
    }


};


/**
 * @brief A hashtable for XDR encodings of NSSI remote functions.
 *
 * The <tt>\ref encodings_ht</tt> hashtable is a thread-safe global
 * data structure that holds encodings of registered NSSI functions.
 * It is used by both clients and servers.
 */
static std::map<int, struct xdr_encodings> encodings_map;
typedef std::map<int, struct xdr_encodings>::iterator encodings_map_iterator_t;
static nthread_lock_t encodings_map_mutex;
static bool xdr_initialized = false;



/* --------------------- Public functions ------------------- */

/**
 * @brief Initialize the XDR encoding mechanism.
 */
int nssi_xdr_init()
{
    if (xdr_initialized) return NSSI_OK;

    /* check to see if logging is enabled */
    if (logger_not_initialized()) {
        logger_set_file(stderr);
    }

    nthread_lock_init(&encodings_map_mutex);

    xdr_initialized = true;

    return NSSI_OK;
}

int register_service_encodings(void)
{
    int rc = NSSI_OK;

    /* get service */
    nssi_register_xdr_encoding(NSSI_OP_GET_SERVICE,
        (xdrproc_t)NULL,
        (xdrproc_t)NULL,
        (xdrproc_t)&xdr_nssi_service);

    /* kill service */
    nssi_register_xdr_encoding(NSSI_OP_KILL_SERVICE,
        (xdrproc_t)&xdr_nssi_kill_service_args,
        (xdrproc_t)NULL,
        (xdrproc_t)xdr_void);

    /* trace reset */
    nssi_register_xdr_encoding(NSSI_OP_TRACE_RESET,
        (xdrproc_t)&xdr_nssi_trace_reset_args,
        (xdrproc_t)NULL,
        (xdrproc_t)xdr_void);

    return rc;
}

/**
 * @brief Finalize the XDR encoding mechanism.
 */
int nssi_xdr_fini()
{
    int rc = NSSI_OK;

    nssi_clear_xdr_encodings();

    nthread_lock_fini(&encodings_map_mutex);

    xdr_initialized=false;

    return rc;
}


/**
 * @brief Register xdr encodings for an opcode.
 *
 * The <tt>\ref nssi_register_xdr_encodings</tt> function registers
 * the xdr encodings functions associated with a particular remote
 * NSSI operation.
 *
 * @param opcode        @input_type   The opcode to lookup.
 * @param encode_args   @input_type   The xdr function used to
 *                                    encode/decode the arguments.
 * @param encode_data   @input_type   The xdr function used to
 *                                    encode/decode the data.
 * @param encode_result @input_type   The xdr function used to
 *                                    encode/decode the result.
 */
int nssi_register_xdr_encoding(int opcode,
    xdrproc_t encode_args,
    xdrproc_t encode_data,
    xdrproc_t encode_result)
{
    int rc = NSSI_OK;
    struct xdr_encodings encodings(encode_args, encode_data, encode_result);

    nssi_xdr_init();

    log_debug(rpc_debug_level,"REGISTERING OPCODE=%d",opcode);

    /* insert structure into map */
    nthread_lock(&encodings_map_mutex);
    encodings_map[(int)opcode] = encodings;  /* should use copy constructor */
    nthread_unlock(&encodings_map_mutex);

    return rc;
}


/**
 * @brief Lookup the encodings for an opcode.
 *
 * The <tt>\ref nssi_lookup_xdr_encodings</tt> function looks up
 * the registered xdr encodings for an opcode.
 *
 * @param opcode        @input_type   The opcode to lookup.
 * @param encode_args   @output_type  Points to the xdr function used to
 *                                    encode/decode the arguments.
 * @param encode_data   @output_type  Points to the xdr function used to
 *                                    encode/decode the data.
 * @param encode_result @output_type  Points to the xdr function used to
 *                                    encode/decode the result.
 *
 * @return <b>\ref NSSI_OK</b> If successful.
 * @return <b>\ref NSSI_ENOENT</b> If the encodings for the specified
 *                                    opcode do not exist.
 */
int nssi_lookup_xdr_encoding(int opcode,
    xdrproc_t *encode_args,
    xdrproc_t *encode_data,
    xdrproc_t *encode_result)
{
    int rc = NSSI_OK;
    encodings_map_iterator_t iter;

    nssi_xdr_init();

    iter = encodings_map.find((int)opcode);
    if ( iter == encodings_map.end()) {
        log_warn(rpc_debug_level, "could not find encodings for opcode=%d",opcode);
        return NSSI_ENOENT;
    }

    *encode_args = iter->second.encode_args;
    *encode_data = iter->second.encode_data;
    *encode_result = iter->second.encode_result;

    return rc;
}


/**
 * @brief Remove all encodings.
 *
 * The <tt>\ref nssi_clear_xdr_encodings</tt> function removes
 * the xdr encodings functions associated with all remote
 * NSSI operations.
 *
 */
int nssi_clear_xdr_encodings()
{
    int rc = NSSI_OK;

    log_debug(rpc_debug_level,"DEREGISTERING ALL OPCODES");

    /* insert structure into map */
    nthread_lock(&encodings_map_mutex);
    encodings_map.clear();
    nthread_unlock(&encodings_map_mutex);

    return rc;
}
