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
/**
 *   @file nssi_client.h
 *
 *   @brief Prototypes for the client-side methods for RPC.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   @author Rolf Riesen (rolf\@cs.sandia.gov)
 *   $Revision: 1640 $
 *   $Date: 2007-11-28 11:59:53 -0700 (Wed, 28 Nov 2007) $
 */

#ifndef _NSSI_RPC_CLNT_H_
#define _NSSI_RPC_CLNT_H_

#include "Trios_nssi_request.h"
#include "Trios_nssi_rpc.h"
#include "Trios_nssi_xdr.h"

/**
 * @defgroup rpc_client_api  Nessie Client API
 *
 * An NSSI client communicates with remote services using an
 * asynchronous, remote procedure call (RPC) interface.
 * As described in @latexonly Section~\ref{sec:Data-movement},@endlatexonly
 * this interface uses the Portals message passing API to
 * take advantage of special features of the network such
 * as RDMA and OS bypass.
 *
 * The typical usage scenerio includes the following steps:
 *
 * -# Obtain the service description by calling the
 *   \ref nssi_lookup_service "\c nssi_lookup_service()" or
 *   \ref nssi_get_service "\c nssi_get_service()" functions.  The service
 *   description, represented by the ref nssi_service data structure,
 *   contains details about how to communicate with the
 *   server, including the remote memory address reserved for incoming
 *   requests and the encoding scheme to use for control messages.
 * -# Call the asynchronous \ref nssi_call_rpc "\c nssi_call_rpc()" function
 *   specifying the operation and providing buffer space for
 *   arguments, bulk data, and results.  The buffers must remain
 *   unmodified by the client until the remote operation completes.
 * -# To check for completion, call the functions \ref nssi_wait "\c nssi_wait()",
 *   \ref nssi_timedwait "\c nssi_timedwait()", or \ref nssi_test "\c nssi_test"
 *   with the <tt>\ref nssi_request</tt> data structure that was created
 *   by the \ref nssi_call_rpc "\c nssi_call_rpc()" function.
 *
 * This remainder of this section contains descriptions of the
 * data structures, special data types, and functions that support
 * client/service communication.
 *
 */

/**
 * @ingroup rpc_client_api
 * @brief Register a client stub.
 *
 * RPC calls require client and server stubs.
 * This function registers the client operation and
 * assigs the appropriate xdr-encoding functions for
 * the arguments, data, and result data structures.
 */
#define NSSI_REGISTER_CLIENT_STUB(opcode, argtype, datatype, rettype) \
    nssi_register_xdr_encoding(opcode, (xdrproc_t)&xdr_ ## argtype, \
            (xdrproc_t)&xdr_ ## datatype, \
            (xdrproc_t)&xdr_ ## rettype)


#ifdef __cplusplus
extern "C" {
#endif



#if defined(__STDC__) || defined(__cplusplus)


    /*
     * @brief Initialize an RPC client.
     *
     * The <em>\ref nssi_rpc_clnt_init</em> function is a blocking call
     * that initializes all internal data structures needed by an
     * NSSI RPC client. This method is called at least once for each
     * remote server that the client expects to use, and the implementation
     * of this function may require a communication with that server.
     *
     * @param server_id  The process ID of the remote server.
     * @param result  Points to a data structure that holds information
     *                       about how to send RPC requests to the remote server.
     *
     * @return <b>\ref NSSI_OK</b> Indicates success.
     * @return <b>\ref NSSI_EBADRPC</b> Indicates an failure in the
     *                                  communication library.
     * @return <b>\ref NSSI_ETIMEDOUT</b> Indicates that the client timed out
     *                                       trying to communicate with the server.
     *
     */
    extern int nssi_rpc_clnt_init(
            NNTI_transport_id_t transport,
            const char *url,
            nssi_service *result);

    /**
     * @brief Initialize an NSSI client.
     *
     * This function must be called before any other calls to
     * the Nessie API.
     */
    extern int nssi_init(const nssi_rpc_transport transport_id);

    /**
     * @brief Finalize an NSSI client.
     *
     */
    extern int nssi_fini(const nssi_rpc_transport transport_id);

    /**
     * @ingroup rpc_client_api
     *
     * @brief Get the service descripton from a known host.
     *
     * The <tt>\ref nssi_get_service</tt> function contacts a
     * remote process to fetch the service description that
     * describes how to send and receive RPC requests to the
     * service.
     *
     * @param server_id  Identifies the host and process ID of the
     *                              remote service.
     * @param timeout  The max amount of time (ms) for the operation. Use -1 for infinite.
     * @param result    Points to the service description of the
     *                              remote service.
     *
     * @return <b>\ref NSSI_OK</b>           Indicates success.
     * @return <b>\ref NSSI_EBADRPC</b>      Indicates a failure in the communication
     *                                library.
     * @return <b>\ref NSSI_ETIMEDOUT</b> Indicates that the client timed out
     *                                trying to communicate with the server.
     */
    extern int nssi_get_service(
            const nssi_rpc_transport rpc_transport,
            const char              *url,
            const int                timeout,
            nssi_service            *result);


    /**
     * @ingroup rpc_client_api
     *
     * @brief Get an array of service descriptons.
     *
     * The <tt>\ref nssi_get_services</tt> function contacts
     * remote processes to fetch service descriptions.
     *
     * @param server_id  Identifies the host and process ID of the
     *                              remote service.
     * @param count    The number of service descriptors to get.
     * @param timeout  The max amount of time (ms) for the operation. Use -1 for infinite.
     * @param result    Points to the service description of the
     *                              remote service.
     *
     * @return <b>\ref NSSI_OK</b>           Indicates success.
     * @return <b>\ref NSSI_EBADRPC</b>      Indicates a failure in the communication
     *                                library.
     * @return <b>\ref NSSI_ETIMEDOUT</b> Indicates that the client timed out
     *                                trying to communicate with the server.
     */
    extern int nssi_get_services(
            const nssi_rpc_transport rpc_transport,
            const char             **url,
            const int                num_servers,
            const int                timeout,
            nssi_service            *result);

    extern int nssi_free_service(
            const nssi_rpc_transport rpc_transport,
            nssi_service            *svc);

    /**
     * @brief Reset tracing on a designated service.
     */
    extern int nssi_trace_reset(
            const nssi_service *svc,
            const char *trace_fname,
            const int trace_ftype,
            const char *enable,
            const long timeout);

    extern int nssi_kill(
            const nssi_service *svc,
            const int sig,
            const long timeout);



    /**
     * @brief Call a remote procedure.
     *
     * @ingroup rpc_client_api
     *
     * The nssi_call_rpc() function encodes and
     * transfers an RPC request to an NSSI server.  It returns
     * an integer value that corresponds to a return code
     * defined in the enumerator \ref nssi_return_codes;
     * however, since the function does not wait for completion
     * of the remote method, the return value only indicates
     * success or failure of marshaling and tranferring
     * the request to the server, not the success of the
     * remote method.
     *
     * The arguments include the function-specific set of arguments
     * required by the remote function, optional data and result
     * arguments, and the <tt>\ref nssi_request</tt> data
     * structure that holds all required information about the
     * pending request.
     *
     * The ``data'' argument is reserved for functions that perform
     * bulk data transfers between the client and the server.
     * This argument  points to client-side memory buffer
     * that the server pulls from (e.g., for reads) or puts
     * into (e.g., for writes).  This buffer must remain
     * unmodified by the client until the remote function completes.
     *
     * The ``result'' argument is reserved for functions that
     * expect a non-integer ``control'' response from the server.
     * The type of the result argument is function specific and points to the
     * client-side buffer reserved for that data. At first glance,
     * this may seem similar to the data argument; however, there
     * are several important distinctions:
     *
     *    -# Results are always directed from the server to the client,
     *       but data could flow in either direction.
     *    -# Data is typically large in comparison to the ``control'' structures
     *       (arguments and results). A clever implementation may transfer
     *       args and results over a network channel optimized for small
     *       messages, but transfer data over a network channel optimized
     *       for bulk transfers.
     *    -# Results are encoded in a portable binary format before transfer,
     *       data is transferred as ``raw'' binary data.
     *
     * The final argument is a pointer to the \ref nssi_request structure that
     * contains all required information about the pending request, including
     * the \ref nssi_request_status "status" of the pending
     * request and the return code of the remote method (if complete).
     * The client calls the functions nssi_get_status(),
     * nssi_test(), or nssi_wait()
     * to asynchronously get the status of a pending request, test for
     * completion, or wait for completion of a request.
     *
     * @param svc            Points to the data structure that describes how
     *                             to communicate with the remote service.
     * @param opcode         Identifies the remote operation to execute.
     * @param args           Points to the unencoded arguments for the remote request.
     * @param data           In not null, points to memory reserved for bulk data
     *                             this buffer must remain available while the request
     *                             is pending.
     * @param data_len       The length of the data buffer (if not null).
     * @param result         Points to a memory reserved for result. This buffer must
     *                             remain available until the request completes.
     * @param req            Points to a data structure that holds information about
     *                              the pending request.
     *
     * @return <b>\ref NSSI_OK</b> Indicates success of encoding and transferring the request.
     * @return <b>\ref NSSI_EENCODE</b> Indicates an error encoding the request.
     * @return <b>\ref NSSI_EBADRPC</b> Indicates an error in the underlying transport library
     *                                  while transferring the request.
     * @return <b>\ref NSSI_ENOTSUP</b> Indicates an unsupported transport mechanism.
     */
    extern int nssi_call_rpc(
            const nssi_service *svc,
            const int opcode,
            void *args,
            void *data,
            uint32_t data_len,
            void *result,
            nssi_request *req);


    /**
     * @brief Blocking remote rpc call.
     *
     *
     */
    extern int nssi_call_rpc_sync(
            const nssi_service *svc,
            const int opcode,
            void *args,
            void *data,
            uint32_t data_len,
            void *result);

    /**
     * @brief Test for completion of an RPC request.
     *
     * @ingroup rpc_client_api
     *
     * The <tt>\ref nssi_test</tt> function checks the status of the specified
     * request and returns <b>TRUE</b> if the status is equal to \ref NSSI_REQUEST_COMPLETE
     * or \ref NSSI_REQUEST_ERROR .
     *
     * @param req    Points to the data structure that holds information
     *                     about the request.
     * @param rc    If the request is complete, points to the return code
     *                     of the completed request.  Undefined otherwise.
     *
     * @return <b>TRUE</b> Indicates that the request is complete.
     * @return <b>FALSE</b> Indicates that the request is not complete.
     */
    extern int nssi_test(
            nssi_request *req,
            int *rc);

    /**
     * @brief Wait for a fixed amount of time for an RPC request to complete.
     *
     * @ingroup rpc_client_api
     *
     * The <tt>\ref nssi_timedwait</tt> function blocks for no more than
     * \em timeout milliseconds waiting for the remote procedure to complete.
     *
     * @param req   Points to the request structure associated with the
     *                    remote function.
     * @param timeout  The maximum amount of time (milliseconds) that the
     *                       function will block waiting for the remote function
     *                       to complete.
     * @param remote_rc    If the remote function completes, this parameter is
     *                     the return code of the remote method. Undefined otherwise.
     *
     * @return <b>\ref NSSI_OK</b> Indicates that the remote function completed
     *                             (possibly with an error).
     * @return <b>\ref NSSI_EBADRPC</b> Indicates failure in the low-level transport mechanism.
     * @return <b>\ref NSSI_ETIMEDOUT</b> Indicates that the remote procedure failed to complete
     *                                within the alloted time.
     */
    extern int nssi_timedwait(
            nssi_request *req,
            int timeout,
            int *remote_rc);

    /**
     * @brief Wait for an RPC request to complete.
     *
     * @ingroup rpc_client_api
     *
     * The <tt>\ref nssi_wait</tt> function blocks until the remote procedure
     * associated with an RPC request completes.
     *
     * @param req   Points to the request structure associated with the
     *                    remote function.
     * @param rc    If the remote function completes, this parameter is
     *                     the return code of the remote method. Undefined otherwise.
     *
     * @return <b>\ref NSSI_OK</b> Indicates that the remote function is complete
     *                             (possibly with an error).
     * @return <b>\ref NSSI_EBADRPC</b> Indicates failure in the low-level transport mechanism.
     */
    extern int nssi_wait(
            nssi_request *req,
            int *rc);


    /**
     * @brief Wait for all requests to complete.
     *
     * A request is not complete unless we receive the short
     * result.
     *
     * @param req_array    The array of pending requests.
     * @param size         The number of pending requests.
     * @param timeout      The time to wait for any one request.
     *
     */
    extern int nssi_waitall(
            nssi_request *req_array,
            nssi_size size,
            int timeout);


    /**
     * @brief Wait for any request to complete.
     *
     * @ingroup rpc_client_api
     *
     * The <tt>\ref nssi_waitany</tt> function blocks for no more than
     * \em timeout milliseconds waiting for any one of an array of requests
     * to complete.
     *
     * @param req_array   Points to an array of requests.
     * @param size        The size of the request array.
     * @param timeout  The maximum amount of time (milliseconds) that the
     *                       function will block waiting for a request to complete.
     * @param which        The index of the complete request.
     * @param remote_rc    The return code of the completed request.
     *
     * @return <b>\ref NSSI_OK</b> Indicates that a request completed
     *                             (possibly with an error).
     * @return <b>\ref NSSI_EBADRPC</b> Indicates failure in the low-level transport mechanism.
     * @return <b>\ref NSSI_ETIMEDOUT</b> Indicates that no request completed
     *                                within the alloted time.
     */
    extern int nssi_waitany(
            nssi_request **req_array,
            nssi_size size,
            int timeout,
            int *which,
            int *remote_rc);

    /**
     * @brief Return the status of an RPC request.
     *
     * @ingroup rpc_client_api
     *
     * The <tt>\ref nssi_get_status</tt> function returns the status (if known)
     * of a remote request.  This function is primarily used for
     * debugging.
     *
     * @param req   Points to the request structure associated with the
     *                    remote function.
     * @param status  The state of the pending request.
     * @param rc      If the status of the remote function is complete,
     *                     this parameter holds the return code of the remote
     *                     method. Undefined otherwise.
     *
     * @return \ref NSSI_OK Indicates success.
     * @return \ref NSSI_EBADRPC Indicates failure to get the status.
     */
    extern int nssi_get_status(
            nssi_request *req,
            int *status,
            int *rc);

    extern int nssi_multicast_rpc(
            const nssi_service *svcs,
            const nssi_size num_svcs,
            const int opcode,
            void *args,
            void *data,
            uint32_t data_size,
            void *results,
            nssi_size result_size, /* size in bytes of the result */
            nssi_request *requests);

    extern int nssi_multicast_rpc_sync(
            const nssi_service *svcs,
            const nssi_size num_svcs,
            const int opcode,
            void *args,
            void *data,
            uint32_t data_size,
            void *results,
            nssi_size result_size);

#else /* K&R C */
#endif



#ifdef __cplusplus
}
#endif

#endif

