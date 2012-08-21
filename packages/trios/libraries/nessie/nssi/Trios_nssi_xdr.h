/* ************************************************************************

                   Trios: Trilinos I/O Support
                 Copyright 2011 Sandia Corporation

 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 the U.S. Government retains certain rights in this software.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are
 met:

 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 3. Neither the name of the Corporation nor the names of the
 contributors may be used to endorse or promote products derived from
 this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)

*************************************************************************/
/*-------------------------------------------------------------------------*/
/**
 *   @file rpc_common.h
 *
 *   @brief Prototypes for the RPC methods used by both clients and servers.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   @author Rolf Riesen (rolf\@cs.sandia.gov)
 *   $Revision: 301 $
 *   $Date: 2005-03-16 00:21:00 -0700 (Wed, 16 Mar 2005) $
 */

#ifndef _NSSI_RPC_XDR_H_
#define _NSSI_RPC_XDR_H_

#include "Trios_config.h"
#include "Trios_xdr.h"



#ifdef __cplusplus
extern "C" {
#endif

#if defined(__STDC__) || defined(__cplusplus)


/**
 * @brief Initialize the XDR encoding mechanism.
 */
extern int nssi_xdr_init();

/**
 * @brief Finalize the XDR encoding mechanism.
 */
extern int nssi_xdr_fini();

/**
  * @brief Register encodings for standard service ops.
  */
extern int register_service_encodings(void);

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
extern int nssi_register_xdr_encoding(int opcode,
    xdrproc_t encode_args,
    xdrproc_t encode_data,
    xdrproc_t encode_result);


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
 * @return <b>\ref NSSI_ERR_NOENT</b> If the encodings for the specified
 *                                    opcode do not exist.
 */
int nssi_lookup_xdr_encoding(int opcode,
    xdrproc_t *encode_args,
    xdrproc_t *encode_data,
    xdrproc_t *encode_result);


/**
 * @brief Remove all encodings.
 *
 * The <tt>\ref nssi_clear_xdr_encodings</tt> function removes
 * the xdr encodings functions associated with all remote
 * NSSI operations.
 *
 */
int nssi_clear_xdr_encodings();

#else /* K&R C */
#endif



#ifdef __cplusplus
}
#endif

#endif

