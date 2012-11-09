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
/*
 * nnti_internal.h
 *
 *  Created on: Feb 3, 2011
 *      Author: thkorde
 */

#ifndef NNTI_INTERNAL_H_
#define NNTI_INTERNAL_H_

#include "Trios_logger.h"

#include <Trios_nnti_xdr.h>

typedef NNTI_result_t (*NNTI_INIT_FN) (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl);

typedef NNTI_result_t (*NNTI_GET_URL_FN)(
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen);

typedef NNTI_result_t (*NNTI_CONNECT_FN) (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl);

typedef NNTI_result_t (*NNTI_DISCONNECT_FN) (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl);

typedef NNTI_result_t (*NNTI_REGISTER_MEMORY_FN) (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf);

typedef NNTI_result_t (*NNTI_UNREGISTER_MEMORY_FN) (
        NNTI_buffer_t    *reg_buf);

typedef NNTI_result_t (*NNTI_SEND_FN) (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl);

typedef NNTI_result_t (*NNTI_PUT_FN) (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

typedef NNTI_result_t (*NNTI_GET_FN) (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

typedef NNTI_result_t (*NNTI_WAIT_FN) (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status);

typedef NNTI_result_t (*NNTI_WAITANY_FN) (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status);

typedef NNTI_result_t (*NNTI_WAITALL_FN) (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status);

typedef NNTI_result_t (*NNTI_FINI_FN) (
        const NNTI_transport_t *trans_hdl);

typedef struct NNTI_transport_ops_t
{
    NNTI_INIT_FN              nnti_init_fn;
    NNTI_GET_URL_FN           nnti_get_url_fn;
    NNTI_CONNECT_FN           nnti_connect_fn;
    NNTI_DISCONNECT_FN        nnti_disconnect_fn;
    NNTI_REGISTER_MEMORY_FN   nnti_register_memory_fn;
    NNTI_UNREGISTER_MEMORY_FN nnti_unregister_memory_fn;
    NNTI_SEND_FN              nnti_send_fn;
    NNTI_PUT_FN               nnti_put_fn;
    NNTI_GET_FN               nnti_get_fn;
    NNTI_WAIT_FN              nnti_wait_fn;
    NNTI_WAITANY_FN           nnti_waitany_fn;
    NNTI_WAITALL_FN           nnti_waitall_fn;
    NNTI_FINI_FN              nnti_fini_fn;
} NNTI_transport_ops_t;


/**
 * @brief The internal representation of a configured transport.
 */
typedef struct {
    /** @brief The transport id. */
    NNTI_transport_id_t  id;

    /** @brief A reference to my process that can be sent to a peer so the peer can contact me. */
    NNTI_peer_t          me;

    /** @brief The transport ops. */
    NNTI_transport_ops_t ops;

    /** @brief Is the transport initialized? */
    uint8_t              initialized;
} NNTI_internal_transport_t;


extern log_level nnti_debug_level;


#endif /* NNTI_INTERNAL_H_ */
