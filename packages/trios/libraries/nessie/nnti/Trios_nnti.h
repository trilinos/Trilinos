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
 * @file Trios_nnti.h
 *
 * @brief Trios_nnti.h
 *
 * @author Todd Kordenbrock (thkorde\@sandia.gov).
 * Created on: Jan 13, 2011
 */

#ifndef NNTI_H_
#define NNTI_H_


#include <Trios_nnti_xdr.h>

#ifdef __cplusplus
extern "C" {
#endif

/** @addtogroup nnti_impl
 *  @{
 */

/**
 * @brief Initialize NNTI to use a specific transport.
 *
 * \param[in]  trans_id  The ID of the transport the client wants to use.
 * \param[in]  my_url    A string that describes the transport parameters.
 * \param[out] trans_hdl A handle to the configured transport.
 * \return A result code (NNTI_OK or an error)
 *
 */
NNTI_result_t NNTI_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl);

/**
 * @brief Return the URL field of this transport.
 *
 * \param[in]  trans_hdl A handle to the configured transport.
 * \param[out] url       A string that describes this process in a transport specific way.
 * \param[in]  maxlen    The length of the 'url' string parameter.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen);

/**
 * @brief Prepare for communication with the peer identified by url.
 *
 * \param[in]  trans_hdl A handle to the configured transport.
 * \param[in]  url       A string that describes a peer in a transport specific way.
 * \param[in]  timeout   The amount of time (in milliseconds) to wait before giving up.
 * \param[out] peer_hdl  A handle to a peer that can be used for network operations.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl);

/**
 * @brief Terminate communication with this peer.
 *
 * \param[in] trans_hdl A handle to the configured transport.
 * \param[in] peer_hdl  A handle to a peer.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl);

/**
 * @brief Prepare a block of memory for network operations.
 *
 * \param[in]  trans_hdl A handle to the configured transport.
 * \param[in]  buffer    Pointer to a memory block.
 * \param[in]  size      The size (in bytes) of 'buffer'.
 * \param[in]  ops       Allowable transport ops for 'buffer'.
 * \param[in]  peer      Only exchange data with this peer (anyone if NULL).
 * \param[out] reg_buf   A pointer to a NNTI_buffer_t.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf);

/**
 * @brief Cleanup after network operations are complete.
 *
 * \param[in]  reg_buf The buffer to cleanup.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_unregister_memory (
        NNTI_buffer_t    *reg_buf);

/**
 * @brief Send a message to a peer.
 *
 * \param[in] peer_hdl  The peer to send the message to.
 * \param[in] msg_hdl   A buffer containing the message to send.
 * \param[in] dest_hdl  A buffer to put the data into.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl);

/**
 * @brief Transfer data to a peer.
 *
 * \param[in] src_buffer_hdl   A buffer containing the data to put.
 * \param[in] src_offset       The offset (in bytes) into the src_buffer from which to put.
 * \param[in] src_length       The number of bytes to put.
 * \param[in] dest_buffer_hdl  A buffer to put the data into.
 * \param[in] dest_offset      The offset (in bytes) into the dest at which to put.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

/**
 * @brief Transfer data from a peer.
 *
 * \param[in] src_buffer_hdl   A buffer containing the data to get.
 * \param[in] src_offset       The offset (in bytes) into the src_buffer from which to get.
 * \param[in] src_length       The number of bytes to get.
 * \param[in] dest_buffer_hdl  A buffer to get the data into.
 * \param[in] dest_offset      The offset (in bytes) into the dest at which to get.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset);

/**
 * @brief Wait for a specific buffer to be acted on by a peer.
 *
 * \param[in]  reg_buf    The buffer that will be acted on by a peer.
 * \param[in]  remote_op  The remote operation that we expect to occur.
 * \param[in]  timeout    The amount of time to wait before giving up.
 * \param[out] status     The details of the completed (or timed out) operation.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status);

/**
 * @brief Wait for any buffer in the list to be acted on by a peer.
 *
 * \param[in]  buf_list   The array of buffers that will be acted on by a peer.
 * \param[in]  buf_count  The number of buffers in the array.
 * \param[in]  remote_op  The remote operation that we expect to occur.
 * \param[in]  timeout    The amount of time to wait before giving up.
 * \param[out] which      The index of the buffer that completed.
 * \param[out] status     The details of the completed (or timed out) operation.
 * \return A result code (NNTI_OK or an error)
 *
 * Caveats: All buffers in buf_list must be registered with the same transport.
 */
NNTI_result_t NNTI_waitany (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status);

/**
 * @brief Wait for all buffers in the list to be acted on by a peer.
 *
 * \param[in]  buf_list   The array of buffers that will be acted on by a peer.
 * \param[in]  buf_count  The number of buffers in the array.
 * \param[in]  remote_op  The remote operation that we expect to occur.
 * \param[in]  timeout    The amount of time to wait before giving up.
 * \param[out] status     The array of statuses that contain the details of the completed (or timed out) operations.
 * \return A result code (NNTI_OK or an error)
 *
 * Caveats: All buffers in buf_list must be registered with the same transport.
 */
NNTI_result_t NNTI_waitall (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status);

/**
 * @brief Disable this transport.
 *
 * \param[in] trans_hdl A handle to the configured transport.
 * \return A result code (NNTI_OK or an error)
*/
NNTI_result_t NNTI_fini (
        const NNTI_transport_t *trans_hdl);

/** @} */




#if defined(HAVE_TRIOS_PORTALS) || defined(HAVE_TRIOS_CRAYPORTALS)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_PORTALS
#elif defined(HAVE_TRIOS_INFINIBAND)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_IB
#elif defined(HAVE_TRIOS_GEMINI)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_GEMINI
#elif defined(HAVE_TRIOS_MPI)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_MPI
#else
#define NSSI_DEFAULT_TRANSPORT NSSI_RPC_LOCAL
#endif




/**
 * @brief Return TRUE if transport 't' has been successfully initialized.
 */
#define NNTI_TRANSPORT_INITIALIZED(t)
/**
 * @brief Return TRUE if buffer 'b' has been filled by the transport.
 */
#define NNTI_BUFFER_FULL(b)
/**
 * @brief Return TRUE if buffer 'b' is empty.
 */
#define NNTI_BUFFER_EMPTY(b)
/**
 * @brief Return a 'char *' pointer to the payload of buffer 'b'.
 */
#define NNTI_BUFFER_C_POINTER(b) ((char *)(b)->payload)
/**
 * @brief Return the size of buffer 'b' in bytes.
 */
#define NNTI_BUFFER_SIZE(b) (b)->payload_size

/**
 * @brief
 *
 * #define NNTI_INIT_BUFFER(b) \
    (b)->buffer_owner.url=""; \
    (b)->peer.url="";
 */
#define NNTI_INIT_BUFFER(b) \
    memset((b)->buffer_owner.url, 0, sizeof((b)->buffer_owner.url)); \
    memset((b)->peer.url, 0, sizeof((b)->peer.url));

/**
 * @brief
 * #define NNTI_INIT_PEER(p) \
    (p)->url="";
 */
#define NNTI_INIT_PEER(p) \
    memset((p)->url, 0, sizeof((p)->url));

#ifdef __cplusplus
}
#endif

#endif /* NNTI_H_*/
