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

typedef NNTI_result_t (*NNTI_callback_fn_t) (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          local_atomic,
        void                   *context);


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
 * @brief Allocate a block of memory and prepare it for network operations.
 *
 * \param[in]  trans_hdl A handle to the configured transport.
 * \param[in]  size      The size (in bytes) of 'buffer'.
 * \param[in]  ops       Allowable transport ops for 'buffer'.
 * \param[out] reg_buf   A pointer to a NNTI_buffer_t.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

/**
 * @brief Cleanup after network operations are complete and release the block
 * of memory.
 *
 * \param[in]  reg_buf The buffer to cleanup.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_free (
        NNTI_buffer_t *reg_buf);

/**
 * @brief Prepare a block of memory for network operations.
 *
 * \param[in]  trans_hdl A handle to the configured transport.
 * \param[in]  buffer    Pointer to a memory block.
 * \param[in]  size      The size (in bytes) of 'buffer'.
 * \param[in]  ops       Allowable transport ops for 'buffer'.
 * \param[out] reg_buf   A pointer to a NNTI_buffer_t.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

/**
 * @brief Prepare a list of memory segments for network operations.
 *
 * \param[in]  trans_hdl       A handle to the configured transport.
 * \param[in]  segments        List of memory segments.
 * \param[in]  segment_lengths List of segment lengths (in bytes).
 * \param[in]  num_segments    The number of segments in the list.
 * \param[in]  ops             Allowable transport ops.
 * \param[out] reg_buf         A pointer to a NNTI_buffer_t.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf);

/**
 * @brief Cleanup after network operations are complete.
 *
 * \param[in]  reg_buf The buffer to cleanup.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_unregister_memory (
        NNTI_buffer_t *reg_buf);

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
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr);

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
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr);

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
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr);

/**
 * @brief Transfer data to a peer.
 *
 * \param[in] src_buffer_hdl    A buffer containing the data to put.
 * \param[in] src_length        The number of bytes to put.
 * \param[in] dest_buffer_list  A list of buffers to put the data into.
 * \param[in] dest_count        The number of destination buffers.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr);


/**
 * @brief Transfer data from a peer.
 *
 * \param[in] src_buffer_list  A list of buffers containing the data to get.
 * \param[in] src_length       The number of bytes to get.
 * \param[in] src_count        The number of source buffers.
 * \param[in] dest_buffer_hdl  A buffer to get the data into.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr);


/**
 * assign a callback to an atomic variable
 *
 * \param[in]  trans_hdl     A handle to the configured transport.
 * \param[in]  local_atomic  index of the local result variable
 * \param[in]  cbfunc        callback function invoked when the atomic variable is modified
 * \param[in]  context       data passed to cbfunc() at every invocation
 *
 * This function assigns a callback to the atomic variable {local_atomic}.
 * {cbfunc} will be invoked after every successful atomic operation on
 * {local_atomic}.  Reads are not considered an atomic operation.
 */
NNTI_result_t NNTI_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context);

/**
 * read a 64-bit value from an local atomic variable
 *
 * \param[in]  trans_hdl     A handle to the configured transport.
 * \param[in]  local_atomic  index of the local atomic variable
 * \param[out] value         current value of the atomic variable
 *
 *
 */
NNTI_result_t NNTI_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value);


/**
 * perform a 64-bit atomic operation with get semantics
 *
 * \param[in]  trans_hdl     A handle to the configured transport.
 * \param[in]  peer          NNTI process that hosts the target atomic variable
 * \param[in]  target_atomic index of the target atomic variable
 * \param[in]  result_atomic index of the local result atomic variable
 * \param[in]  operand       64-bit operand to the atomic operation
 * \param[in]  op            atomic operation to execute
 * \param[out] wr            work request to wait on
 *
 * This function executes an atomic operation with get semantics. When the operation
 * is complete the result of the operation is visible in {target_atomic} and the
 * previous value is visible in {result_atomic}.
 */
NNTI_result_t NNTI_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr);

/**
 * perform a 64-bit compare-and-swap operation
 *
 * \param[in]  trans_hdl       A handle to the configured transport.
 * \param[in]  peer            NNTI process that hosts the target atomic variable
 * \param[in]  target_atomic   index of the target atomic variable
 * \param[in]  result_atomic   index of the local result atomic variable
 * \param[in]  compare_operand 64-bit operand to compare with
 * \param[in]  swap_operand    64-bit operand to swap in
 * \param[out] wr              work request to wait on
 *
 * This function executes an atomic operation with get semantics. When the operation
 * is complete the result of the operation is visible in (target_atomic} and the
 * previous value is visible in {result_atomic}.
 */
NNTI_result_t NNTI_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr);


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr);


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_clear_work_request (
        NNTI_work_request_t  *wr);


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_destroy_work_request (
        NNTI_work_request_t  *wr);


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_cancel (
        NNTI_work_request_t *wr);


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_interrupt (
        const NNTI_transport_t *trans_hdl);


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
        NNTI_work_request_t *wr,
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
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
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
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
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




#if defined(HAVE_TRIOS_INFINIBAND)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_IB
#elif defined(HAVE_TRIOS_GEMINI)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_GEMINI
#elif defined(HAVE_TRIOS_BGPDCMF)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_DCMF
#elif defined(HAVE_TRIOS_BGQPAMI)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_PAMI
#elif defined(HAVE_TRIOS_MPI)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_MPI
#elif defined(HAVE_TRIOS_PORTALS) || defined(HAVE_TRIOS_CRAYPORTALS)
#define NNTI_DEFAULT_TRANSPORT NNTI_TRANSPORT_PORTALS
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
