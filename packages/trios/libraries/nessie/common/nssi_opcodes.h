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
/**
 *   @file rpc_opcodes.h
 *
 *   @brief Opcodes used by the all servers and clients.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 540 $
 *   $Date: 2006-01-11 16:57:56 -0700 (Wed, 11 Jan 2006) $
 */

#ifndef _RPC_OPCODES_H_
#define _RPC_OPCODES_H_

#include "Trios_nssi_types.h"

#ifdef __cplusplus
extern "C" {
#endif

	enum nssi_rpc_opcodes {
	    /** @brief Null NSSI operation. */
	    NSSI_OP_NULL = 1000,

		/** @brief Create a new object. */
		NSSI_OP_GET_SERVICE,

		/** @brief Kill an active service. */
		NSSI_OP_KILL_SERVICE,

		/** @brief Kill an active service. */
		NSSI_OP_TRACE_RESET
	};



#ifdef __cplusplus
}
#endif

#endif

