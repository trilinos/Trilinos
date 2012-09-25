/*
 * nnti_gni.h
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#ifndef NNTI_GNI_H_
#define NNTI_GNI_H_


#include "Trios_nnti.h"
#include "nnti_internal.h"



#ifdef __cplusplus
extern "C" {
#endif

#define GNI_SET_MATCH_ANY(p) \
(p)->peer.transport_id = NNTI_TRANSPORT_GEMINI; \
(p)->peer.NNTI_remote_process_t_u.gni.addr   = 0; \
(p)->peer.NNTI_remote_process_t_u.gni.port   = 0;


NNTI_result_t NNTI_gni_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl);

NNTI_result_t NNTI_gni_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen);

NNTI_result_t NNTI_gni_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl);

NNTI_result_t NNTI_gni_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl);

NNTI_result_t NNTI_gni_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf);

NNTI_result_t NNTI_gni_unregister_memory (
        NNTI_buffer_t    *reg_buf);

NNTI_result_t NNTI_gni_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl);

NNTI_result_t NNTI_gni_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

NNTI_result_t NNTI_gni_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

NNTI_result_t NNTI_gni_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status);

NNTI_result_t NNTI_gni_waitany (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status);

NNTI_result_t NNTI_gni_waitall (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status);

NNTI_result_t NNTI_gni_fini (
        const NNTI_transport_t *trans_hdl);

#ifdef __cplusplus
}
#endif

#endif /* NNTI_GNI_H_*/
