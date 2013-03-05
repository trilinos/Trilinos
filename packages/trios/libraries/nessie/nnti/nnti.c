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
 * @file nnti.c
 *
 * @brief nnti.c
 *
 * @author Todd Kordenbrock (thkorde\@sandia.gov).
 * Created on: Jan 13, 2011
 */

#include "Trios_config.h"

#include <string.h>

#include "Trios_nnti.h"
#include "nnti_internal.h"

#if defined(HAVE_TRIOS_PORTALS) || defined(HAVE_TRIOS_CRAYPORTALS)
#include "nnti_ptls.h"
#endif
#if defined(HAVE_TRIOS_INFINIBAND)
#include "nnti_ib.h"
#endif
#if defined(HAVE_TRIOS_LUC)
#include "nnti_luc.h"
#endif
#if defined(HAVE_TRIOS_GEMINI)
#include "nnti_gni.h"
#endif
#if defined(HAVE_TRIOS_MPI)
#include "nnti_mpi.h"
#endif

#include "Trios_logger.h"



/* set to LOG_UNDEFINED -- log commands will use default level */
log_level nnti_debug_level = LOG_UNDEFINED;


static NNTI_internal_transport_t available_transports[NNTI_TRANSPORT_COUNT];


/**
 * @brief Initialize NNTI to use a specific transport.
 *
 * Enable the use of a particular transport by this process.  <tt>my_url</tt>
 * allows the process to have some control (if possible) over the
 * URL assigned for the transport.  For example, a Portals URL to put
 * might be "ptl://-1,128".  This would tell Portals to use the default
 * network ID, but use PID=128.  If the transport
 * can be initialized without this info (eg. a Portals client), <tt>my_url</tt> can
 * be NULL or empty.
 *
 */
NNTI_result_t NNTI_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    NNTI_result_t rc=NNTI_OK;
    static int8_t first_init=TRUE;

    if (first_init==TRUE) {
        memset(&available_transports[0], 0, NNTI_TRANSPORT_COUNT*sizeof(NNTI_internal_transport_t));
        first_init=FALSE;
    }

#if defined(HAVE_TRIOS_PORTALS) || defined(HAVE_TRIOS_CRAYPORTALS)
    if (trans_id == NNTI_TRANSPORT_PORTALS) {
        available_transports[trans_id].initialized                   = 1;
        available_transports[trans_id].ops.nnti_init_fn              = NNTI_ptl_init;
        available_transports[trans_id].ops.nnti_get_url_fn           = NNTI_ptl_get_url;
        available_transports[trans_id].ops.nnti_connect_fn           = NNTI_ptl_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn        = NNTI_ptl_disconnect;
        available_transports[trans_id].ops.nnti_register_memory_fn   = NNTI_ptl_register_memory;
        available_transports[trans_id].ops.nnti_unregister_memory_fn = NNTI_ptl_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn              = NNTI_ptl_send;
        available_transports[trans_id].ops.nnti_put_fn               = NNTI_ptl_put;
        available_transports[trans_id].ops.nnti_get_fn               = NNTI_ptl_get;
        available_transports[trans_id].ops.nnti_wait_fn              = NNTI_ptl_wait;
        available_transports[trans_id].ops.nnti_waitany_fn           = NNTI_ptl_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn           = NNTI_ptl_waitall;
        available_transports[trans_id].ops.nnti_fini_fn              = NNTI_ptl_fini;
    }
#endif
#if defined(HAVE_TRIOS_INFINIBAND)
    if (trans_id == NNTI_TRANSPORT_IB) {
        available_transports[trans_id].initialized                   = 1;
        available_transports[trans_id].ops.nnti_init_fn              = NNTI_ib_init;
        available_transports[trans_id].ops.nnti_get_url_fn           = NNTI_ib_get_url;
        available_transports[trans_id].ops.nnti_connect_fn           = NNTI_ib_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn        = NNTI_ib_disconnect;
        available_transports[trans_id].ops.nnti_register_memory_fn   = NNTI_ib_register_memory;
        available_transports[trans_id].ops.nnti_unregister_memory_fn = NNTI_ib_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn              = NNTI_ib_send;
        available_transports[trans_id].ops.nnti_put_fn               = NNTI_ib_put;
        available_transports[trans_id].ops.nnti_get_fn               = NNTI_ib_get;
        available_transports[trans_id].ops.nnti_wait_fn              = NNTI_ib_wait;
        available_transports[trans_id].ops.nnti_waitany_fn           = NNTI_ib_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn           = NNTI_ib_waitall;
        available_transports[trans_id].ops.nnti_fini_fn              = NNTI_ib_fini;
    }
#endif
#if defined(HAVE_TRIOS_LUC)
    if (trans_id == NNTI_TRANSPORT_LUC) {
        available_transports[trans_id].initialized                   = 1;
        available_transports[trans_id].ops.nnti_init_fn              = NNTI_luc_init;
        available_transports[trans_id].ops.nnti_get_url_fn           = NNTI_luc_get_url;
        available_transports[trans_id].ops.nnti_connect_fn           = NNTI_luc_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn        = NNTI_luc_disconnect;
        available_transports[trans_id].ops.nnti_register_memory_fn   = NNTI_luc_register_memory;
        available_transports[trans_id].ops.nnti_unregister_memory_fn = NNTI_luc_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn              = NNTI_luc_send;
        available_transports[trans_id].ops.nnti_put_fn               = NNTI_luc_put;
        available_transports[trans_id].ops.nnti_get_fn               = NNTI_luc_get;
        available_transports[trans_id].ops.nnti_wait_fn              = NNTI_luc_wait;
        available_transports[trans_id].ops.nnti_waitany_fn           = NNTI_luc_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn           = NNTI_luc_waitall;
        available_transports[trans_id].ops.nnti_fini_fn              = NNTI_luc_fini;
    }
#endif
#if defined(HAVE_TRIOS_GEMINI)
    if (trans_id == NNTI_TRANSPORT_GEMINI) {
        available_transports[trans_id].initialized                   = 1;
        available_transports[trans_id].ops.nnti_init_fn              = NNTI_gni_init;
        available_transports[trans_id].ops.nnti_get_url_fn           = NNTI_gni_get_url;
        available_transports[trans_id].ops.nnti_connect_fn           = NNTI_gni_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn        = NNTI_gni_disconnect;
        available_transports[trans_id].ops.nnti_register_memory_fn   = NNTI_gni_register_memory;
        available_transports[trans_id].ops.nnti_unregister_memory_fn = NNTI_gni_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn              = NNTI_gni_send;
        available_transports[trans_id].ops.nnti_put_fn               = NNTI_gni_put;
        available_transports[trans_id].ops.nnti_get_fn               = NNTI_gni_get;
        available_transports[trans_id].ops.nnti_wait_fn              = NNTI_gni_wait;
        available_transports[trans_id].ops.nnti_waitany_fn           = NNTI_gni_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn           = NNTI_gni_waitall;
        available_transports[trans_id].ops.nnti_fini_fn              = NNTI_gni_fini;
    }
#endif
#if defined(HAVE_TRIOS_MPI)
    if (trans_id == NNTI_TRANSPORT_MPI) {
        available_transports[trans_id].initialized                   = 1;
        available_transports[trans_id].ops.nnti_init_fn              = NNTI_mpi_init;
        available_transports[trans_id].ops.nnti_get_url_fn           = NNTI_mpi_get_url;
        available_transports[trans_id].ops.nnti_connect_fn           = NNTI_mpi_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn        = NNTI_mpi_disconnect;
        available_transports[trans_id].ops.nnti_register_memory_fn   = NNTI_mpi_register_memory;
        available_transports[trans_id].ops.nnti_unregister_memory_fn = NNTI_mpi_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn              = NNTI_mpi_send;
        available_transports[trans_id].ops.nnti_put_fn               = NNTI_mpi_put;
        available_transports[trans_id].ops.nnti_get_fn               = NNTI_mpi_get;
        available_transports[trans_id].ops.nnti_wait_fn              = NNTI_mpi_wait;
        available_transports[trans_id].ops.nnti_waitany_fn           = NNTI_mpi_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn           = NNTI_mpi_waitall;
        available_transports[trans_id].ops.nnti_fini_fn              = NNTI_mpi_fini;
    }
#endif

    trans_hdl->id = trans_id;

    rc = available_transports[trans_id].ops.nnti_init_fn(
            trans_id,
            my_url,
            trans_hdl);

    return(rc);
}


/**
 * @brief Return the URL field of this transport.
 *
 * Return the URL field of this transport.  After initialization, the transport will
 * have a specific location on the network where peers can contact it.  The
 * transport will convert this location to a string that other instances of the
 * transport will recognize.
 *
 * URL format: "transport://address/memory_descriptor"
 *    - transport - (required) identifies how the URL should parsed
 *    - address   - (required) uniquely identifies a location on the network
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port", "luc://endpoint_id/"
 *    - memory_descriptor - (optional) transport-specific representation of RMA params
 *
 */
NNTI_result_t NNTI_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_get_url_fn(
                trans_hdl,
                url,
                maxlen);
    }

    return(rc);
}

/**
 * @brief Prepare for communication with the peer identified by <tt>url</tt>.
 *
 * Parse <tt>url</tt> in a transport specific way.  Perform any transport specific
 * actions necessary to begin communication with this peer.
 *
// * If the peer is found and responds
// * to a ping, a handle will be allocated and assigned to the pointer.  This
// * handle should be used to move data to/from the peer.
 *
 * Connectionless transport: parse and populate
 * Connected transport: parse, connection and populate
 *
 */
NNTI_result_t NNTI_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_connect_fn(
                trans_hdl,
                url,
                timeout,
                peer_hdl);
    }

    return(rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_disconnect_fn(
                trans_hdl,
                peer_hdl);
    }

    return(rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 *
 */
NNTI_result_t NNTI_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_register_memory_fn(
                trans_hdl,
                buffer,
                element_size,
                num_elements,
                ops,
                peer,
                reg_buf);
    }

    return(rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 *
 */
NNTI_result_t NNTI_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[reg_buf->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[reg_buf->transport_id].ops.nnti_unregister_memory_fn(
                reg_buf);
    }

    return(rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 *
 */
NNTI_result_t NNTI_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[msg_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[msg_hdl->transport_id].ops.nnti_send_fn(
                peer_hdl,
                msg_hdl,
                dest_hdl);
    }

    return(rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[src_buffer_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[src_buffer_hdl->transport_id].ops.nnti_put_fn(
                src_buffer_hdl,
                src_offset,
                src_length,
                dest_buffer_hdl,
                dest_offset);
    }

    return(rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[dest_buffer_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[dest_buffer_hdl->transport_id].ops.nnti_get_fn(
                src_buffer_hdl,
                src_offset,
                src_length,
                dest_buffer_hdl,
                dest_offset);
    }

    return(rc);
}


/**
 * @brief Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on <tt>reg_buf</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 */
NNTI_result_t NNTI_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[reg_buf->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[reg_buf->transport_id].ops.nnti_wait_fn(
                reg_buf,
                remote_op,
                timeout,
                status);
    }

    return(rc);
}


/**
 * @brief Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on any buffer in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the receive queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_waitany (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_transport_id_t id=NNTI_TRANSPORT_NULL;
    uint32_t i=0;

    for (i=0;i<buf_count;i++) {
        if (buf_list[i]) {
            id=buf_list[i]->transport_id;
        }
    }

    if (id == NNTI_TRANSPORT_NULL) {
        rc=NNTI_EINVAL;
    } else {
        if (available_transports[id].initialized==0) {
            rc=NNTI_ENOTINIT;
        } else {
            rc = available_transports[id].ops.nnti_waitany_fn(
                    buf_list,
                    buf_count,
                    remote_op,
                    timeout,
                    which,
                    status);
        }
    }

    return(rc);
}


/**
 * @brief Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on all buffers in <tt>buf_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in buf_list must be registered with the same transport.
 *   2) You can't wait on the receive queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_waitall (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_transport_id_t id=NNTI_TRANSPORT_NULL;
    uint32_t i=0;

    for (i=0;i<buf_count;i++) {
        if (buf_list[i]) {
            id=buf_list[i]->transport_id;
        }
    }

    if (id == NNTI_TRANSPORT_NULL) {
        rc=NNTI_EINVAL;
    } else {
        if (available_transports[id].initialized==0) {
            rc=NNTI_ENOTINIT;
        } else {
            rc = available_transports[id].ops.nnti_waitall_fn(
                    buf_list,
                    buf_count,
                    remote_op,
                    timeout,
                    status);
        }
    }

    return(rc);
}


/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_fini (
        const NNTI_transport_t *trans_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_fini_fn(
                trans_hdl);
        memset(&available_transports[trans_hdl->id], 0, sizeof(NNTI_internal_transport_t));
    }

    return(rc);
}
