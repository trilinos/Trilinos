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
 *   $Revision: 1640 $
 *   $Date: 2007-11-28 11:59:53 -0700 (Wed, 28 Nov 2007) $
 */

#ifndef _NSSI_RPC_COMMON_H_
#define _NSSI_RPC_COMMON_H_

#include "Trios_config.h"
#include "Trios_nssi_types.h"

#ifdef __cplusplus
extern "C" {
#endif

#if defined(HAVE_TRIOS_PORTALS) || defined(HAVE_TRIOS_CRAYPORTALS)
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_PTL
#elif defined(HAVE_TRIOS_INFINIBAND)
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_IB
#elif defined(HAVE_TRIOS_LUC)
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_LUC
#elif defined(HAVE_TRIOS_GEMINI)
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_GEMINI
#elif defined(HAVE_TRIOS_MPI)
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_MPI
#else
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_LOCAL
#endif

#define NSSI_DEFAULT_ENCODE NSSI_RPC_XDR

#if defined(__STDC__) || defined(__cplusplus)

    extern void *memdup(void *src, int size);

    /**
     * @brief Initialize the RPC mechanism.
     *
     * The <tt>\ref nssi_rpc_init</tt> method
     * initializes the rpc library.
     *
     * @param rpc_transport @input_type  Identifies the transport mechanism
     *                                   to use for communication.
     * @param rpc_encode    @input_type  Identifies the mechanism used to
     *                                   encode the rpc control messages.
         *
     * @return \ref NSSI_OK Indicates sucess.
     * @return \ref NSSI_ERR_RPC Indicates failure in the NSSI RPC library.
     */
    extern int nssi_rpc_init(
            const nssi_rpc_transport rpc_transport,
            const nssi_rpc_encode    rpc_encode,
            const char              *url);

    /**
     * @brief Finalize the RPC mechanism.
     *
     * The \b nssi_rpc_fini method performs any cleanup
     * required by the underlying communication protocol
     * before exiting the client application.
         *
     * @return \ref NSSI_OK Indicates sucess.
     * @return \ref NSSI_ERR_RPC Indicates failure in the NSSI RPC library.
     */
    extern int nssi_rpc_fini(const nssi_rpc_transport rpc_transport);

    /**
     * @brief Get the process ID of this process.
     *
     * @ingroup rpc_funcs
     *
     * The <tt>\ref nssi_get_my_pid</tt> method gets the
         * \ref nssi_remote_pid "process id" of the calling process.
     *
     * @param id @output_type If successful, points to the \ref nssi_remote_pid "process ID"
         *                   of the calling process. Undefined otherwise.
     *
     * @return \ref NSSI_OK Indicates sucess.
     * @return \ref NSSI_ERR_RPC Indicates failure in the NSSI RPC library.
         *
     */
    extern int nssi_get_url(
            const nssi_rpc_transport rpc_transport,
            char                    *url,
            uint32_t                 maxlen);

#else /* K&R C */
#endif



#ifdef __cplusplus
}
#endif

#endif
