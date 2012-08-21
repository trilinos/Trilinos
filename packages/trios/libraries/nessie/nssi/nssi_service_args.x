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
/**  @file nssi_service_args.x
 *
 *   @brief XDR definitions for the argument structures for
 *   \ref generic service operations.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov).
 *   $Revision: 969 $.
 *   $Date: 2006-08-28 15:40:44 -0600 (Mon, 28 Aug 2006) $.
 *
 */

#ifdef RPC_HDR
%#include "Trios_xdr.h"
#endif

#ifdef RPC_XDR
%#include "Trios_xdr.h"
#endif



/**
 * @brief Arguments for the \ref nssi_create_container method that
 * have to be passed to the authorization server.
 */
struct nssi_get_service_args {

    /** @brief The ID of the service (unused). */
    int id;
};

/**
 * @brief Arguments for the \ref nssi_kill_service method that
 * have to be passed to the authorization server.
 */
struct nssi_kill_service_args {

    /** @brief The signal to use when killing the service. */
    int sig;
};


/**
 * @brief Arguments for the \ref nssi_reset_tracing method.
 */
struct nssi_trace_reset_args {

    /** @brief The file type for the new tracing file. */
    int ftype;

    /** @brief The name of the new tracing file. */
    string fname<PATH_MAX>;

    /** @brief Comma-separated list of traces to enable. */
    string enable<512>;
};

