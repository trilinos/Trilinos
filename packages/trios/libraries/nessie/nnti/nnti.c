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
#if defined(HAVE_TRIOS_GEMINI)
#include "nnti_gni.h"
#endif
#if defined(HAVE_TRIOS_BGPDCMF)
#include "nnti_dcmf.h"
#endif
#if defined(HAVE_TRIOS_BGQPAMI)
#include "nnti_pami.h"
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
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_ptl_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_ptl_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_ptl_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_ptl_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_ptl_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_ptl_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_ptl_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_ptl_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_ptl_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_ptl_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_ptl_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_ptl_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_ptl_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_ptl_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_ptl_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_ptl_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_ptl_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_ptl_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_ptl_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_ptl_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_ptl_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_ptl_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_ptl_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_ptl_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_ptl_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_ptl_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_ptl_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_ptl_fini;
    }
#endif
#if defined(HAVE_TRIOS_INFINIBAND)
    if (trans_id == NNTI_TRANSPORT_IB) {
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_ib_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_ib_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_ib_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_ib_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_ib_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_ib_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_ib_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_ib_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_ib_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_ib_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_ib_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_ib_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_ib_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_ib_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_ib_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_ib_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_ib_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_ib_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_ib_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_ib_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_ib_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_ib_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_ib_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_ib_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_ib_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_ib_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_ib_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_ib_fini;
    }
#endif
#if defined(HAVE_TRIOS_GEMINI)
    if (trans_id == NNTI_TRANSPORT_GEMINI) {
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_gni_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_gni_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_gni_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_gni_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_gni_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_gni_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_gni_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_gni_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_gni_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_gni_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_gni_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_gni_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_gni_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_gni_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_gni_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_gni_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_gni_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_gni_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_gni_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_gni_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_gni_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_gni_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_gni_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_gni_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_gni_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_gni_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_gni_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_gni_fini;
    }
#endif
#if defined(HAVE_TRIOS_BGPDCMF)
    if (trans_id == NNTI_TRANSPORT_DCMF) {
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_bgpdcmf_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_bgpdcmf_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_bgpdcmf_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_bgpdcmf_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_bgpdcmf_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_bgpdcmf_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_bgpdcmf_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_bgpdcmf_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_bgpdcmf_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_bgpdcmf_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_bgpdcmf_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_bgpdcmf_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_bgpdcmf_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_bgpdcmf_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_bgpdcmf_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_bgpdcmf_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_bgpdcmf_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_bgpdcmf_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_bgpdcmf_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_bgpdcmf_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_bgpdcmf_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_bgpdcmf_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_bgpdcmf_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_bgpdcmf_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_bgpdcmf_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_bgpdcmf_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_bgpdcmf_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_bgpdcmf_fini;
    }
#endif
#if defined(HAVE_TRIOS_BGQPAMI)
    if (trans_id == NNTI_TRANSPORT_PAMI) {
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_bgqpami_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_bgqpami_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_bgqpami_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_bgqpami_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_bgqpami_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_bgqpami_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_bgqpami_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_bgqpami_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_bgqpami_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_bgqpami_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_bgqpami_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_bgqpami_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_bgqpami_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_bgqpami_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_bgqpami_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_bgqpami_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_bgqpami_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_bgqpami_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_bgqpami_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_bgqpami_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_bgqpami_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_bgqpami_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_bgqpami_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_bgqpami_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_bgqpami_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_bgqpami_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_bgqpami_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_bgqpami_fini;
    }
#endif
#if defined(HAVE_TRIOS_MPI)
    if (trans_id == NNTI_TRANSPORT_MPI) {
        available_transports[trans_id].initialized                      = 1;
        available_transports[trans_id].ops.nnti_init_fn                 = NNTI_mpi_init;
        available_transports[trans_id].ops.nnti_get_url_fn              = NNTI_mpi_get_url;
        available_transports[trans_id].ops.nnti_connect_fn              = NNTI_mpi_connect;
        available_transports[trans_id].ops.nnti_disconnect_fn           = NNTI_mpi_disconnect;
        available_transports[trans_id].ops.nnti_alloc_fn                = NNTI_mpi_alloc;
        available_transports[trans_id].ops.nnti_free_fn                 = NNTI_mpi_free;
        available_transports[trans_id].ops.nnti_register_memory_fn      = NNTI_mpi_register_memory;
        available_transports[trans_id].ops.nnti_register_segments_fn    = NNTI_mpi_register_segments;
        available_transports[trans_id].ops.nnti_unregister_memory_fn    = NNTI_mpi_unregister_memory;
        available_transports[trans_id].ops.nnti_send_fn                 = NNTI_mpi_send;
        available_transports[trans_id].ops.nnti_put_fn                  = NNTI_mpi_put;
        available_transports[trans_id].ops.nnti_get_fn                  = NNTI_mpi_get;
        available_transports[trans_id].ops.nnti_scatter_fn              = NNTI_mpi_scatter;
        available_transports[trans_id].ops.nnti_gather_fn               = NNTI_mpi_gather;
        available_transports[trans_id].ops.nnti_atomic_set_callback_fn  = NNTI_mpi_atomic_set_callback;
        available_transports[trans_id].ops.nnti_atomic_read_fn          = NNTI_mpi_atomic_read;
        available_transports[trans_id].ops.nnti_atomic_fop_fn           = NNTI_mpi_atomic_fop;
        available_transports[trans_id].ops.nnti_atomic_cswap_fn         = NNTI_mpi_atomic_cswap;
        available_transports[trans_id].ops.nnti_create_work_request_fn  = NNTI_mpi_create_work_request;
        available_transports[trans_id].ops.nnti_clear_work_request_fn   = NNTI_mpi_clear_work_request;
        available_transports[trans_id].ops.nnti_destroy_work_request_fn = NNTI_mpi_destroy_work_request;
        available_transports[trans_id].ops.nnti_cancel_fn               = NNTI_mpi_cancel;
        available_transports[trans_id].ops.nnti_cancelall_fn            = NNTI_mpi_cancelall;
        available_transports[trans_id].ops.nnti_interrupt_fn            = NNTI_mpi_interrupt;
        available_transports[trans_id].ops.nnti_wait_fn                 = NNTI_mpi_wait;
        available_transports[trans_id].ops.nnti_waitany_fn              = NNTI_mpi_waitany;
        available_transports[trans_id].ops.nnti_waitall_fn              = NNTI_mpi_waitall;
        available_transports[trans_id].ops.nnti_fini_fn                 = NNTI_mpi_fini;
    }
#endif

    trans_hdl->datatype = NNTI_dt_transport;
    trans_hdl->id       = trans_id;

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
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port"
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

    peer_hdl->datatype = NNTI_dt_peer;

    return(rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 *
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
 * @brief Allocate a block of memory and prepare it for network operations.
 *
 * Allocate a block of memory in a transport specific way and wrap it in an
 * NNTI_buffer_t.  The transport may take additional actions to prepare the
 * memory for network send/receive.
 *
 */
NNTI_result_t NNTI_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_alloc_fn(
                trans_hdl,
                element_size,
                num_elements,
                ops,
                reg_buf);
    }

    reg_buf->datatype = NNTI_dt_buffer;

    return(rc);
}


/**
 * @brief Cleanup after network operations are complete and release the block
 * of memory.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_alloc().
 *
 */
NNTI_result_t NNTI_free (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[reg_buf->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[reg_buf->transport_id].ops.nnti_free_fn(
                reg_buf);
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
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        if ((ops == NNTI_BOP_SEND_SRC) || (ops == NNTI_BOP_RECV_DST) || (ops == NNTI_BOP_RECV_QUEUE)) {
            log_error(nnti_debug_level, "NNTI_BOP_SEND_SRC, NNTI_BOP_RECV_DST and NNTI_BOP_RECV_QUEUE types require the use of NNTI_alloc().");
            rc=NNTI_EINVAL;
        } else {
            rc = available_transports[trans_hdl->id].ops.nnti_register_memory_fn(
                    trans_hdl,
                    buffer,
                    element_size,
                    num_elements,
                    ops,
                    reg_buf);
        }
    }

    reg_buf->datatype = NNTI_dt_buffer;

    return(rc);
}


/**
 * @brief Prepare a list of memory segments for network operations.
 *
 * Wrap a list of user allocated memory segments in an NNTI_buffer_t.  The
 * transport may take additional actions to prepare the memory segments for
 * network send/receive.  If the memory segments don't meet the transport's
 * requirements for memory regions, then errors or poor performance may
 * result.
 *
 */
NNTI_result_t NNTI_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        if ((ops == NNTI_BOP_SEND_SRC) || (ops == NNTI_BOP_RECV_DST) || (ops == NNTI_BOP_RECV_QUEUE)) {
            log_error(nnti_debug_level, "NNTI_BOP_SEND_SRC, NNTI_BOP_RECV_DST and NNTI_BOP_RECV_QUEUE types cannot be segmented.");
            rc=NNTI_EINVAL;
        } else {
            rc = available_transports[trans_hdl->id].ops.nnti_register_segments_fn(
                    trans_hdl,
                    segments,
                    segment_lengths,
                    num_segments,
                    ops,
                    reg_buf);
        }
    }

    reg_buf->datatype = NNTI_dt_buffer;

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
 * @brief Calculate the number of bytes required to store an encoded NNTI datatype.
 *
 */
NNTI_result_t NNTI_dt_sizeof (
        const NNTI_transport_t *trans_hdl,
        void                   *nnti_dt,
        uint64_t               *packed_len)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        xdrproc_t sizeof_fn;
        NNTI_datatype_t *dt=(NNTI_datatype_t*)nnti_dt;

        switch (*dt) {
            case NNTI_dt_transport:
                sizeof_fn=(xdrproc_t)&xdr_NNTI_transport_t;
                break;
            case NNTI_dt_peer:
                sizeof_fn=(xdrproc_t)&xdr_NNTI_peer_t;
                break;
            case NNTI_dt_buffer:
                sizeof_fn=(xdrproc_t)&xdr_NNTI_buffer_t;
                break;
            case NNTI_dt_work_request:
                sizeof_fn=(xdrproc_t)&xdr_NNTI_work_request_t;
                break;
            case NNTI_dt_status:
                sizeof_fn=(xdrproc_t)&xdr_NNTI_status_t;
                break;
        }
        *packed_len = xdr_sizeof(sizeof_fn, nnti_dt);
        *packed_len += sizeof(NNTI_datatype_t);
    }

    return(rc);
}


/**
 * @brief Encode an NNTI datatype into an array of bytes.
 *
 */
NNTI_result_t NNTI_dt_pack (
        const NNTI_transport_t *trans_hdl,
        void                   *nnti_dt,
        char                   *packed_buf,
        uint64_t                packed_buflen)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        XDR       encode_xdrs;
        xdrproc_t encode_fn;
        NNTI_datatype_t *dt=(NNTI_datatype_t*)nnti_dt;

        switch (*dt) {
            case NNTI_dt_transport:
                encode_fn=(xdrproc_t)&xdr_NNTI_transport_t;
                break;
            case NNTI_dt_peer:
                encode_fn=(xdrproc_t)&xdr_NNTI_peer_t;
                break;
            case NNTI_dt_buffer:
                encode_fn=(xdrproc_t)&xdr_NNTI_buffer_t;
                break;
            case NNTI_dt_work_request:
                encode_fn=(xdrproc_t)&xdr_NNTI_work_request_t;
                break;
            case NNTI_dt_status:
                encode_fn=(xdrproc_t)&xdr_NNTI_status_t;
                break;
        }

        *(NNTI_datatype_t*)packed_buf = *dt;

        packed_buf += sizeof(NNTI_datatype_t);
        packed_buflen -= sizeof(NNTI_datatype_t);

        xdrmem_create(
                &encode_xdrs,
                packed_buf,
                packed_buflen,
                XDR_ENCODE);

        if (!encode_fn(&encode_xdrs, nnti_dt)) {
            log_fatal(nnti_debug_level,"packing failed");
            return NNTI_EENCODE;
        }
    }

    return(rc);
}


/**
 * @brief Decode an array of bytes into an NNTI datatype.
 *
 */
NNTI_result_t NNTI_dt_unpack (
        const NNTI_transport_t *trans_hdl,
        void                   *nnti_dt,
        char                   *packed_buf,
        uint64_t                packed_buflen)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        XDR       decode_xdrs;
        xdrproc_t decode_fn;
        NNTI_datatype_t *dt=(NNTI_datatype_t*)packed_buf;
        uint64_t dt_size;

        switch (*dt) {
            case NNTI_dt_transport:
                decode_fn=(xdrproc_t)&xdr_NNTI_transport_t;
                dt_size=sizeof(NNTI_transport_t);
                break;
            case NNTI_dt_peer:
                decode_fn=(xdrproc_t)&xdr_NNTI_peer_t;
                dt_size=sizeof(NNTI_peer_t);
                break;
            case NNTI_dt_buffer:
                decode_fn=(xdrproc_t)&xdr_NNTI_buffer_t;
                dt_size=sizeof(NNTI_buffer_t);
                break;
            case NNTI_dt_work_request:
                decode_fn=(xdrproc_t)&xdr_NNTI_work_request_t;
                dt_size=sizeof(NNTI_work_request_t);
                break;
            case NNTI_dt_status:
                decode_fn=(xdrproc_t)&xdr_NNTI_status_t;
                dt_size=sizeof(NNTI_status_t);
                break;
        }

        packed_buf += sizeof(NNTI_datatype_t);
        packed_buflen -= sizeof(NNTI_datatype_t);

        memset(nnti_dt, 0, dt_size);
        xdrmem_create(
                &decode_xdrs,
                packed_buf,
                packed_buflen,
                XDR_DECODE);
        if (!decode_fn(&decode_xdrs, nnti_dt)) {
            log_fatal(nnti_debug_level,"unpacking failed");
            return NNTI_EDECODE;
        }
    }

    return(rc);
}


/**
 * @brief Free a variable size NNTI datatype that was unpacked with NNTI_dt_unpack().
 *
 */
NNTI_result_t NNTI_dt_free (
        const NNTI_transport_t *trans_hdl,
        void                   *nnti_dt)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        xdrproc_t free_fn;
        NNTI_datatype_t *dt=(NNTI_datatype_t*)nnti_dt;

        switch (*dt) {
            case NNTI_dt_transport:
                free_fn=(xdrproc_t)&xdr_NNTI_transport_t;
                break;
            case NNTI_dt_peer:
                free_fn=(xdrproc_t)&xdr_NNTI_peer_t;
                break;
            case NNTI_dt_buffer:
                free_fn=(xdrproc_t)&xdr_NNTI_buffer_t;
                break;
            case NNTI_dt_work_request:
                free_fn=(xdrproc_t)&xdr_NNTI_work_request_t;
                break;
            case NNTI_dt_status:
                free_fn=(xdrproc_t)&xdr_NNTI_status_t;
                break;
        }

        xdr_free(free_fn, nnti_dt);
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
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[msg_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[msg_hdl->transport_id].ops.nnti_send_fn(
                peer_hdl,
                msg_hdl,
                dest_hdl,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

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
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
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
                dest_offset,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

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
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
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
                dest_offset,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Transfer data to a peer.
 *
 */
NNTI_result_t NNTI_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[src_buffer_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[src_buffer_hdl->transport_id].ops.nnti_scatter_fn(
                src_buffer_hdl,
                src_length,
                dest_buffer_list,
                dest_count,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Transfer data from a peer.
 *
 */
NNTI_result_t NNTI_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[dest_buffer_hdl->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[dest_buffer_hdl->transport_id].ops.nnti_gather_fn(
                src_buffer_list,
                src_length,
                src_count,
                dest_buffer_hdl,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


NNTI_result_t NNTI_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_atomic_set_callback_fn(
                trans_hdl,
                local_atomic,
                cbfunc,
                context);
    }

    return(rc);
}


NNTI_result_t NNTI_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_atomic_read_fn(
                trans_hdl,
                local_atomic,
                value);
    }

    return(rc);
}


NNTI_result_t NNTI_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_atomic_fop_fn(
                trans_hdl,
                peer_hdl,
                target_atomic,
                result_atomic,
                operand,
                op,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


NNTI_result_t NNTI_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_atomic_cswap_fn(
                trans_hdl,
                peer_hdl,
                target_atomic,
                result_atomic,
                compare_operand,
                swap_operand,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[reg_buf->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[reg_buf->transport_id].ops.nnti_create_work_request_fn(
                reg_buf,
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_clear_work_request (
        NNTI_work_request_t  *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[wr->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[wr->transport_id].ops.nnti_clear_work_request_fn(
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_destroy_work_request (
        NNTI_work_request_t  *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[wr->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[wr->transport_id].ops.nnti_destroy_work_request_fn(
                wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_cancel (
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[wr->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[wr->transport_id].ops.nnti_cancel_fn(wr);
    }

    wr->datatype = NNTI_dt_work_request;

    return(rc);
}


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_transport_id_t id=NNTI_TRANSPORT_NULL;
    uint32_t i=0;

    for (i=0;i<wr_count;i++) {
        if (wr_list[i]) {
            id=wr_list[i]->transport_id;
        }
    }

    if (id == NNTI_TRANSPORT_NULL) {
        rc=NNTI_EINVAL;
    } else {
        if (available_transports[id].initialized==0) {
            rc=NNTI_ENOTINIT;
        } else {
            rc = available_transports[id].ops.nnti_cancelall_fn(
                    wr_list,
                    wr_count);
        }
    }

    for (i=0;i<wr_count;i++) {
        if (wr_list[i]) {
            wr_list[i]->datatype = NNTI_dt_work_request;
        }
    }

    return(rc);
}


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_interrupt (
        const NNTI_transport_t *trans_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[trans_hdl->id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[trans_hdl->id].ops.nnti_interrupt_fn(
                trans_hdl);
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
        NNTI_work_request_t *wr,
        const int            timeout,
        NNTI_status_t       *status)
{
    NNTI_result_t rc=NNTI_OK;

    if (available_transports[wr->transport_id].initialized==0) {
        rc=NNTI_ENOTINIT;
    } else {
        rc = available_transports[wr->transport_id].ops.nnti_wait_fn(
                wr,
                timeout,
                status);
    }

    status->datatype = NNTI_dt_status;

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
 *
 */
NNTI_result_t NNTI_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_transport_id_t id=NNTI_TRANSPORT_NULL;
    uint32_t i=0;

    for (i=0;i<wr_count;i++) {
        if (wr_list[i]) {
            id=wr_list[i]->transport_id;
        }
    }

    if (id == NNTI_TRANSPORT_NULL) {
        rc=NNTI_EINVAL;
    } else {
        if (available_transports[id].initialized==0) {
            rc=NNTI_ENOTINIT;
        } else {
            rc = available_transports[id].ops.nnti_waitany_fn(
                    wr_list,
                    wr_count,
                    timeout,
                    which,
                    status);
        }
    }

    status->datatype = NNTI_dt_status;

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
 *
 */
NNTI_result_t NNTI_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status)
{
    NNTI_result_t rc=NNTI_OK;
    NNTI_transport_id_t id=NNTI_TRANSPORT_NULL;
    uint32_t i=0;

    for (i=0;i<wr_count;i++) {
        if (wr_list[i]) {
            id=wr_list[i]->transport_id;
        }
    }

    if (id == NNTI_TRANSPORT_NULL) {
        rc=NNTI_EINVAL;
    } else {
        if (available_transports[id].initialized==0) {
            rc=NNTI_ENOTINIT;
        } else {
            rc = available_transports[id].ops.nnti_waitall_fn(
                    wr_list,
                    wr_count,
                    timeout,
                    status);
        }
    }

    for (i=0;i<wr_count;i++) {
        if (status[i]) {
            status[i]->datatype = NNTI_dt_status;
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
