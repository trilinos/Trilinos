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
 *   @file nssi_server.h
 *
 *   @brief Prototypes for the server-side methods for RPC.
 *
 *   @author Ron Oldfield (raoldfi\@sandia.gov)
 *   $Revision: 1640 $
 *   $Date: 2007-11-28 11:59:53 -0700 (Wed, 28 Nov 2007) $
 */

#ifndef _NSSI_RPC_SRVR_H_
#define _NSSI_RPC_SRVR_H_

#include "Trios_nssi_rpc.h"



/**
 * @defgroup rpc_server_api  Nessie Server API
 *
 *  The server-side APIs for Nessie provide methods to identify
 *  supported call-back functions, register a service, transport
 *  data between client and server, and control the operations of
 *  an active service.
 */

/* forward declaration */
struct nssi_request;

/* this is only for C++ code */
#ifdef __cplusplus

/*#include "nssi_TypeManip.h"*/
#include "Trios_nssi_TypeManip.h"

/**
 * @ingroup rpc_server_api
 * @brief A macro for registering server callback functions.
 *
 * This macro constructs the service op structure with appropriate
 * XDR-encoding functions, then adds the function to the list of
 * supported ops for this service.
 *
 * @param id  The opcode of the function to register.
 * @param fptr  The function pointer.
 * @param arg_type Type for the function input arguments.
 * @param res_type Type of the result. The stub will use nssi_send_result.
 *
 * @note We need a special version for C++ because sizeof(void)
 * is not allowed in C++.  See nssi_tests.h for
 * the implementation of SizeOf<T>.
 */
#define NSSI_REGISTER_SERVER_STUB(id, fptr, arg_type, res_type) \
    { \
        nssi_svc_op op; \
        op.opcode = id; \
        op.func = (nssi_rpc_proc)&fptr; \
        op.sizeof_args = Nessie::SizeOf<arg_type>::value; \
        op.decode_args = (xdrproc_t)&xdr_ ## arg_type; \
        op.sizeof_res = Nessie::SizeOf<res_type>::value; \
        op.encode_res = (xdrproc_t)&xdr_ ## res_type; \
        nssi_service_add_op( NULL, &op); \
    }
#else

/**
 * @ingroup rpc_server_api
 * @brief A macro for registering server callback functions.
 *
 * This macro constructs the service op structure with appropriate
 * XDR-encoding functions, then adds the function to the list of
 * supported ops for this service.
 *
 * @param id  The opcode of the function to register.
 * @param fptr  The function pointer.
 * @param arg_type Types for the function input arguments.
 * @param res_type Type of the result. The stub will use nssi_send_result.
 *
 */
#define NSSI_REGISTER_SERVER_STUB(id, fptr, arg_type, res_type) \
    { \
        nssi_svc_op op; \
        op.opcode = id; \
        op.func = (nssi_rpc_proc)&fptr; \
        op.sizeof_args = sizeof(arg_type); \
        op.decode_args = (xdrproc_t)&xdr_ ## arg_type; \
        op.sizeof_res = sizeof(res_type); \
        op.encode_res = (xdrproc_t)&xdr_ ## res_type; \
        nssi_service_add_op( NULL, &op); \
    }
#endif




#ifdef __cplusplus
extern "C" {
#endif


        /**
         * @ingroup rpc_server_api
         * @brief Definition of a function pointer for RPC services.
         *
         * Each server-side stub uses the following parameters:
         *    - \ref nssi_remote_pid *caller : the remote PID of the calling process
         *    - void *args : arguments from the client
         *    - const nssi_rma *data_addr : remote address for get/put of raw data on the client
         *    - const nssi_rma *result_addr: remote address for result buffer (for \ref nssi_send_result)
         */
        typedef int (*nssi_rpc_proc) (
                const unsigned long request_id,
                const NNTI_peer_t *,
                const void *,
                const NNTI_buffer_t *,
                const NNTI_buffer_t *);


        /**
         * @brief A structure associated with an operation made available
         * by a registered service.
         *
         * The \b nssi_svc_op structure includes all fields required to
         * spawn call a local method when an RPC request arrives.
         */
        typedef struct nssi_svc_op {
                /** @brief The operation code.  */
                int opcode;

                /** @brief Function pointer for the server-side function.  */
                nssi_rpc_proc func;

                /** @brief Size of the arguments structure. */
                uint32_t sizeof_args;

                /** @brief A function to decode the arguments when a request arrives. */
                xdrproc_t decode_args;

                /** @brief Size of the result structure. */
                uint32_t sizeof_res;

                /** @brief A function to encode the result after servicing the request. */
                xdrproc_t encode_res;
        } nssi_svc_op;


        /**
         * @brief Structure for a list of RPC services.
         */
        typedef struct nssi_svc_op_list {

                /** @brief The service for this entry. */
                nssi_svc_op svc_op;

                /** @brief The remaining list of services. */
                struct nssi_svc_op_list *next;
        } nssi_svc_op_list;

        typedef struct rpc_request {
            nssi_service *svc;
            NNTI_peer_t caller;
            char *req_buf;
            nssi_size short_req_len;
            int id;
            double arrival_time;
        } nssi_svc_rpc_request;

#if defined(__STDC__) || defined(__cplusplus)

        /**
         * @ingroup rpc_server_api
         * @brief Create a daemon process.
         */
        extern void nssi_daemonize();

        /**
         * @ingroup rpc_server_api
         * @brief An abstract method to get data from a remote memory descriptor.
         *
         * The server stub uses this function to get or put data to a
         * client memory descriptor.
         *
         * @param caller    Remote process ID of the destination node.
         * @param buf          the buffer for the data.
         * @param len          the maximum length of the buffer.
         * @param data_addr    the remote memory descriptor.
         *
         * @todo Get rid of the caller parameter.
         */
        extern int nssi_get_data(
                const NNTI_peer_t *caller,
                void *buf,
                const int len,
                const NNTI_buffer_t *data_addr);

        /**
         * @ingroup rpc_server_api
         * @brief An abstract method to put data into a remote memory descriptor.
         *
         * The server stub uses this function to put data to a
         * client memory descriptor.
         *
         * @param caller      Remote process ID of the destination.
         * @param buf         Buffer for the data.
         * @param len         The amount of data to send.
         * @param data_addr   The remote memory descriptor.
         * @param timeout     Time to wait for response from the destination.
         *
         * @todo Get rid of the caller parameter.
         */
        extern int nssi_put_data(
                const NNTI_peer_t *caller,
                const void *buf,
                const int len,
                const NNTI_buffer_t *data_addr,
                const long timeout);

        /**
         * @ingroup rpc_server_api
         * @brief Send result back to the client.
         *
         * This function allows the service library to send a result to
         * a remote memory descriptor.  The result will be packaged as
         * an \ref nssi_result_header, then PUT on the memory descriptor of
         * the remote address.
         *
         * @param caller  The remote process ID of the caller.
         * @param request_id ID of the request associated with this result.
         * @param rc      The return code for the server function.
         * @param result  Pointer to the result data structure (NULL if not used)
         * @param result_addr The remote memory address of the client (where
         * to PUT the result)
         *
         * @returns \ref NSSI_OK if successful.
         */
        extern int nssi_send_result(
                const NNTI_peer_t *caller,
                unsigned long request_id,
                const int rc,
                void *result,
                const NNTI_buffer_t *result_addr);


        /**
         * @ingroup rpc_server_api
         * @brief Initialize an RPC server.
         *
         * @ingroup rpc_server_api
         *
         * This method initializes the portals library and allocates space
         * for incoming requests.
         *
         * @param match_bits  the portals match bits
         * @param short_req_len  the length of a portals short request queue
         * @param service     the service descriptor to register (to register for clients).
         */
        extern int nssi_service_init(
                const nssi_rpc_transport rpc_transport,
                const int short_req_len,
                nssi_service *service);


        /**
         * @ingroup rpc_server_api
         * @brief Cleanly abort the RPC server.
         *
         * The abort function kills a running service by sending a
         * SIGINT signal to the running process.  If the service
         * has already started,  the signal is caught by the
         * signal handler.
         */
        extern void nssi_service_abort();

        /**
         * @ingroup rpc_server_api
         * @brief Returns true if the service needs to shut down.
         */
        extern int nssi_exit_now();


        /**
         * @ingroup rpc_server_api
         * @brief Process an encoded RPC request.
         *
         * This function processes an encoded RPC request by decoding the
         * header, then calling a registered callback based on the opcode
         * sent in the header.
         */
        extern int nssi_process_rpc_request(
                nssi_svc_rpc_request *req);


        /**
         * @ingroup rpc_server_api
         * @brief Start an RPC service.
         *
         * The \b nssi_service_start implements a loop that waits for incoming
         * RPC requests, then calls the nssi_process_rpc_request function to
         * process those requests.
         *
         * @param service  The service descriptor.
         */
        extern int nssi_service_start(
                nssi_service *service);

        /**
         * This function is essentially the same as nssi_start, but it allows
         * the caller to pass in a different function to process encoded requests.
         * This is useful, for example to implement a multithreaded service. The
         * code could call an "enqueue_rpc_request" function and have a separate
         * thread that pops requests off the queue and processes them.
         */
        extern int nssi_service_start_wfn(
                nssi_service *svc,
                int (*process_req)(nssi_svc_rpc_request *req));


        /*
         * Private function for registering new callbacks for this service.
         */
        extern int nssi_service_add_op(
                const nssi_service *svc,
                const nssi_svc_op *ops);


        /**
         * @ingroup rpc_server_api
         * @brief Close down an active service.
         *
         * Shuts down the service and releases any associated resources.
         *
         * @param service  The service descriptor.
         */
        extern int nssi_service_fini(
                const nssi_service *service);

        /**
         * @ingroup rpc_server_api
         * @brief Register an RPC service.
         *
         * This method creates a named RPC service on the specified
         * registry server.  Along with the name of the service, the
         * server has to specify where (in the form of an \ref nssi_rma) the
         * client should "put" requests.
         *
         * @param registry_id  Process ID of the registry server.
         * @param name  Name of the service.
         * @param service  The service description to associate with the name.
         * @param req     The request handle (used to test for completion).
         *
         * @note This function is not yet implemented.
         */
        extern int nssi_register_service(
                const NNTI_peer_t registry_id,
                const char *name,
                const nssi_service *service,
                nssi_request *req);

        /**
         * @brief Return the age, in seconds, of the specified request.
         *
         * For active requests, the age is the current time minus the arrival time
         * of the request.  For expired request (if we support those), it is the
         * completion-arrival time.
         */
        extern double nssi_get_request_age(const NNTI_peer_t *caller, const int req_id);


#else /* K&R C */
#endif



#ifdef __cplusplus
}
#endif

#endif

