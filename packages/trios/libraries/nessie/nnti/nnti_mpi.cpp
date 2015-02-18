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
 * nnti_ptls.c
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#include "Trios_config.h"
#include "Trios_threads.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_nnti_fprint_types.h"

// MPI is only used to increment the PID
#ifdef HAVE_TRIOS_MPI
#include <mpi.h>
#endif

#include <assert.h>
#include <string.h>

#include <map>
#include <deque>
#include <algorithm>

#include "nnti_mpi.h"
#include "nnti_utils.h"



/* if defined, the RDMA initiator will send an ACK message to the RDMA
 * target when the RDMA op is complete.  the target process must wait
 * on the target buffer in order to get the ACK.  this creates two-sided
 * semantics for RDMA ops.   in this mode, when the wait returns the
 * the RDMA op is complete and status indicates what data was addressed.
 */
#define USE_RDMA_TARGET_ACK
/* if undefined, the ACK message is NOT sent to the RDMA target when
 * the RDMA op is complete.  this creates one-sided semantics for RDMA
 * ops.  in this mode, the target has no idea when the RDMA op is
 * complete and what data was addressed.  NNTI_wait() returns NNTI_EINVAL
 * if passed a target buffer.
 */
#undef USE_RDMA_TARGET_ACK


typedef struct {

	uint32_t min_atomics_vars;

} nnti_mpi_config;


#define NNTI_MPI_REQUEST_TAG          0x01
#define NNTI_MPI_ATOMICS_REQUEST_TAG  0x02
#define NNTI_MPI_ATOMICS_RESULT_TAG   0x03


#define MPI_OP_PUT_INITIATOR  1
#define MPI_OP_GET_INITIATOR  2
#define MPI_OP_PUT_TARGET     3
#define MPI_OP_GET_TARGET     4
#define MPI_OP_SEND_REQUEST   5
#define MPI_OP_SEND_BUFFER    6
#define MPI_OP_NEW_REQUEST    7
#define MPI_OP_RECEIVE        8
#define MPI_OP_FETCH_ADD      9
#define MPI_OP_COMPARE_SWAP  10


typedef enum {
    BUFFER_INIT=0,
    SEND_COMPLETE=1,
    RECV_COMPLETE,
    RDMA_WRITE_INIT,
    RDMA_RTS_COMPLETE,
    RDMA_WRITE_COMPLETE,
    RDMA_READ_INIT,
    RDMA_RTR_COMPLETE,
    RDMA_READ_COMPLETE,
    RDMA_TARGET_INIT,
    RDMA_TARGET_COMPLETE,
    RDMA_COMPLETE
} mpi_op_state_t;

typedef struct {
    int32_t  tag;
    uint32_t op;
    uint64_t offset;
    uint64_t length;
} mpi_ready_msg;

typedef enum {
	MPI_ATOMIC_FETCH_ADD  =1,
	MPI_ATOMIC_CMP_AND_SWP=2
} mpi_atomic_op_t;

typedef struct {
	uint8_t  op;
    uint32_t index;
    int64_t  compare_add;
    int64_t  swap;
} mpi_atomic_request_msg;

typedef struct {
    int64_t result;
} mpi_atomic_result_msg;

#define RTR_REQ_INDEX       0
#define RTS_REQ_INDEX       1
#define SEND_INDEX          2
#define GET_RECV_INDEX      3
#define GET_SEND_INDEX      4
#define PUT_RECV_INDEX      5
#define PUT_SEND_INDEX      6
#define RECV_INDEX          7
#define ATOMICS_SEND_INDEX  8
#define ATOMICS_RECV_INDEX  9

#define MAX_INDEX           10

#define RTR_REQUEST_ACTIVE           0x001
#define RTS_REQUEST_ACTIVE           0x002
#define SEND_REQUEST_ACTIVE          0x004
#define GET_RECV_REQUEST_ACTIVE      0x008
#define GET_SEND_REQUEST_ACTIVE      0x010
#define PUT_RECV_REQUEST_ACTIVE      0x020
#define PUT_SEND_REQUEST_ACTIVE      0x040
#define RECV_REQUEST_ACTIVE          0x080
#define ATOMICS_SEND_REQUEST_ACTIVE  0x100
#define ATOMICS_RECV_REQUEST_ACTIVE  0x200

typedef struct mpi_work_request {
    NNTI_work_request_t *nnti_wr;

    NNTI_buffer_t  *reg_buf;
    NNTI_peer_t     peer;
    uint64_t        src_offset;
    uint64_t        dst_offset;
    uint64_t        length;
    int32_t         tag;
    MPI_Request     request[MAX_INDEX];
    MPI_Request    *request_ptr;
    uint32_t        request_count;
    int             request_index;
    uint8_t         active_requests;

    mpi_ready_msg   rtr_msg;
    mpi_ready_msg   rts_msg;

    int             atomics_result_index;

    mpi_op_state_t  op_state;

    MPI_Status      last_event;
    uint8_t         last_op;
    uint8_t         is_last_op_complete;
} mpi_work_request;

typedef std::deque<mpi_work_request *>           wr_queue_t;
typedef std::deque<mpi_work_request *>::iterator wr_queue_iter_t;

typedef struct mpi_memory_handle {
    int64_t rtr_tag;
    int64_t rts_tag;
    int64_t data_tag;

    wr_queue_t     wr_queue;
    nthread_lock_t wr_queue_lock;
} mpi_memory_handle;


#define NUM_REQ_QUEUES 2
typedef struct mpi_request_queue_handle {
    NNTI_buffer_t *reg_buf;

    /* incoming queue */
    char *req_queue;

    /* each message is no larger than req_size */
    int req_size;

    /* number of requests slots in the queue */
    int req_count;

    MPI_Request *mpi_request_list;
    uint64_t     last_request_index;

} mpi_request_queue_handle;

typedef struct {
	nthread_lock_t lock;
	int64_t        value;
} mpi_atomic_t;

typedef struct mpi_transport_global {

    int  rank;
    int  size;
    char proc_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm nnti_comm;

    nthread_counter_t mbits;

    mpi_request_queue_handle req_queue;

    mpi_atomic_t           *atomics;
    mpi_atomic_request_msg  atomics_request_msg;
    mpi_atomic_result_msg   atomics_result_msg;
    MPI_Request             atomics_recv_request;
    MPI_Request             atomics_send_request;

    bool init_called_mpi_init;

} mpi_transport_global;



static nthread_lock_t nnti_mpi_lock;


static int process_event(
        mpi_work_request *mpi_wr,
        const MPI_Status *event);
static NNTI_result_t setup_atomics(void);
static int check_atomic_operation(void);
static int check_target_buffer_progress(void);
static NNTI_result_t post_atomics_recv_request(void);
static NNTI_result_t post_recv_dst_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         tag,
        uint64_t        length);
static NNTI_result_t post_recv_queue_work_request(
        NNTI_buffer_t    *reg_buf,
        int64_t           tag,
        uint64_t          length,
        uint64_t          count,
        MPI_Request      *request_list);
static NNTI_result_t post_RTR_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t post_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t post_RTR_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t repost_recv_work_request(
        mpi_work_request *mpi_wr);
static NNTI_result_t repost_RTR_recv_work_request(
        mpi_work_request *mpi_wr);
static NNTI_result_t repost_RTS_recv_work_request(
        mpi_work_request *mpi_wr);
static NNTI_result_t repost_RTR_RTS_recv_work_request(
        mpi_work_request *mpi_wr);
static int is_wr_complete(
        mpi_work_request *mpi_wr);
static int8_t is_any_wr_complete(
        mpi_work_request **wr_list,
        const uint32_t     wr_count,
        uint32_t          *which);
static int8_t is_all_wr_complete(
        mpi_work_request **wr_list,
        const uint32_t     wr_count);
static int8_t is_wr_complete(
        NNTI_work_request_t *wr);
static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which);
static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);

static NNTI_result_t insert_target_buffer(NNTI_buffer_t *buf);
static NNTI_buffer_t *del_target_buffer(NNTI_buffer_t *buf);
//static void print_target_buffer_deque();

static void create_status(
        NNTI_work_request_t *wr,
        mpi_work_request    *mpi_wr,
        int                  nnti_rc,
        NNTI_status_t       *status);
static void create_peer(
        NNTI_peer_t *peer,
        int          rank);

static void config_init(
        nnti_mpi_config *c);
static void config_get_from_env(
        nnti_mpi_config *c);


#define MPI_MEM_HDL(b) ((mpi_memory_handle *)((b)->transport_private))
#define MPI_WORK_REQUEST(wr) ((mpi_work_request *)((wr)->transport_private))


/* Thomas Wang's 64 bit to 32 bit Hash Function (http://www.concentric.net/~ttwang/tech/inthash.htm) */
static uint32_t hash6432shift(uint64_t key)
{
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21;             // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t)key;
}

static std::map<uint32_t, NNTI_buffer_t *> buffers_by_bufhash;
typedef std::map<uint32_t, NNTI_buffer_t *>::iterator buf_by_bufhash_iter_t;
typedef std::pair<uint32_t, NNTI_buffer_t *> buf_by_bufhash_t;
static nthread_lock_t nnti_buf_bufhash_lock;

static std::map<uint32_t, mpi_work_request *> wr_by_wrhash;
typedef std::map<uint32_t, mpi_work_request *>::iterator wr_by_wrhash_iter_t;
typedef std::pair<uint32_t, mpi_work_request *> wr_by_wrhash_t;
static nthread_lock_t nnti_wr_wrhash_lock;

typedef std::deque<NNTI_buffer_t *>           target_buffer_queue_t;
typedef std::deque<NNTI_buffer_t *>::iterator target_buffer_queue_iter_t;
static nthread_lock_t                        nnti_target_buffer_queue_lock;

target_buffer_queue_t target_buffers;


static nnti_mpi_config config;


static mpi_transport_global transport_global_data;
static const int MAX_SLEEP = 10;  /* in milliseconds */

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
 */
NNTI_result_t NNTI_mpi_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    static uint8_t initialized=FALSE;

    int mpi_initialized=FALSE;
    int name_len=0;


    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);


    if (!initialized) {

        nthread_lock_init(&nnti_mpi_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);
        nthread_lock_init(&nnti_wr_wrhash_lock);
        nthread_lock_init(&nnti_target_buffer_queue_lock);

        config_init(&config);
        config_get_from_env(&config);

        if (my_url != NULL) {
            log_error(nnti_debug_level,"The MPI transport does not accept a URL at init.  Ignoring URL.");
        }

        log_debug(nnti_debug_level, "initializing MPI transport");

        memset(&transport_global_data, 0, sizeof(mpi_transport_global));

        rc=MPI_Initialized(&mpi_initialized);
        if (rc) {
            log_fatal(nnti_debug_level,"MPI_Initialized() failed, %d", rc);
            abort();
        }

        if (mpi_initialized==FALSE) {
            int argc=0;
            char **argv=NULL;
            int provided=-1;

            log_debug(nnti_debug_level, "initializing MPI library");

            rc=MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
            if (rc) {
                log_fatal(nnti_debug_level,"MPI_Init_thread() failed, %d", rc);
                abort();
            }

            transport_global_data.init_called_mpi_init=true;
        }

//        MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
        MPI_Comm_size(MPI_COMM_WORLD, &transport_global_data.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &transport_global_data.rank);
        MPI_Get_processor_name(transport_global_data.proc_name, &name_len);

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "MPI Initialized: rank=%llu, size=%llu, proc_name=%s\n",
                    (unsigned long long)transport_global_data.rank,
                    (unsigned long long)transport_global_data.size,
                    transport_global_data.proc_name);
        }

        nthread_counter_init(&transport_global_data.mbits);
        for (int i=0;i<0x111;i++) {
        	// mbits 0x000-0x110 are reserved.  increment mbits to 0x111.
        	uint64_t v=nthread_counter_increment(&transport_global_data.mbits);
        }

        setup_atomics();

        create_peer(&trans_hdl->me, transport_global_data.rank);

        initialized = TRUE;
    }


    log_debug(nnti_debug_level, "exit");


    return(nnti_rc);
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
 */
NNTI_result_t NNTI_mpi_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    return(nnti_rc);
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
NNTI_result_t NNTI_mpi_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];

    int peer_rank;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(peer_hdl);

    if (url != NULL) {
        if ((nnti_rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            return(nnti_rc);
        }
        if (0!=strcmp(transport, "mpi")) {
            /* the peer described by 'url' is not an MPI peer */
            return(NNTI_EINVAL);
        }

        if ((nnti_rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            return(nnti_rc);
        }

        peer_rank=strtol(address, NULL, 0);
    } else {
        /*  */
        return(NNTI_EINVAL);
    }

    create_peer(
            peer_hdl,
            peer_rank);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_mpi_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    assert(trans_hdl);
    assert(peer_hdl);

    return(nnti_rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_mpi_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    char *buf=(char *)malloc(element_size*num_elements);
    assert(buf);

    nnti_rc=NNTI_mpi_register_memory(
            trans_hdl,
            buf,
            element_size,
            num_elements,
            ops,
            reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_mpi_alloc", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_mpi_free (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    char *buf=NNTI_BUFFER_C_POINTER(reg_buf);
    assert(buf);

    nnti_rc=NNTI_mpi_unregister_memory(reg_buf);

    free(buf);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_mpi_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
//    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    mpi_mem_hdl=new mpi_memory_handle();
    assert(mpi_mem_hdl);
    nthread_lock_init(&mpi_mem_hdl->wr_queue_lock);

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)mpi_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (ops == NNTI_BOP_RECV_QUEUE) {
        mpi_mem_hdl->rtr_tag  = 0;
        mpi_mem_hdl->rts_tag  = 0;
        mpi_mem_hdl->data_tag = NNTI_MPI_REQUEST_TAG;

    } else if (ops == NNTI_BOP_RECV_DST) {
        mpi_mem_hdl->rtr_tag  = 0;
        mpi_mem_hdl->rts_tag  = 0;
        mpi_mem_hdl->data_tag = nthread_counter_increment(&transport_global_data.mbits);

    } else if (ops == NNTI_BOP_SEND_SRC) {
        mpi_mem_hdl->rtr_tag  = 0;
        mpi_mem_hdl->rts_tag  = 0;
        mpi_mem_hdl->data_tag = nthread_counter_increment(&transport_global_data.mbits);
    } else {
        mpi_mem_hdl->rtr_tag  = nthread_counter_increment(&transport_global_data.mbits);
        mpi_mem_hdl->rts_tag  = nthread_counter_increment(&transport_global_data.mbits);
        mpi_mem_hdl->data_tag = nthread_counter_increment(&transport_global_data.mbits);
    }

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(1, sizeof(NNTI_remote_addr_t));
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=1;

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].transport_id                      = NNTI_TRANSPORT_MPI;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.size     = element_size;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.rtr_tag  = mpi_mem_hdl->rtr_tag;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.rts_tag  = mpi_mem_hdl->rts_tag;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag = mpi_mem_hdl->data_tag;

    if (ops == NNTI_BOP_RECV_QUEUE) {
        mpi_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        q_hdl->reg_buf=reg_buf;

        q_hdl->req_queue         =buffer;
        q_hdl->req_size          =element_size;
        q_hdl->req_count         =num_elements;
        q_hdl->last_request_index=0;
        q_hdl->mpi_request_list  =(MPI_Request*)calloc(q_hdl->req_count, sizeof(MPI_Request));

        /* initialize the buffer */
        memset(q_hdl->req_queue, 0, q_hdl->req_count*q_hdl->req_size);

        post_recv_queue_work_request(
                reg_buf,
                NNTI_MPI_REQUEST_TAG,
                q_hdl->req_size,
                q_hdl->req_count,
                q_hdl->mpi_request_list);

    } else if (ops == NNTI_BOP_RECV_DST) {
        post_recv_dst_work_request(
                reg_buf,
                mpi_mem_hdl->data_tag,
                element_size);

    } else if (ops == NNTI_BOP_REMOTE_WRITE) {
        post_RTS_recv_work_request(reg_buf);
        insert_target_buffer(reg_buf);

    } else if (ops == NNTI_BOP_REMOTE_READ) {
        post_RTR_recv_work_request(reg_buf);
        insert_target_buffer(reg_buf);

    } else if (ops == (NNTI_BOP_REMOTE_READ|NNTI_BOP_REMOTE_WRITE)) {
        post_RTR_RTS_recv_work_request(reg_buf);
        insert_target_buffer(reg_buf);

    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_mpi_register_memory", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
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
NNTI_result_t NNTI_mpi_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    return NNTI_ENOTSUP;
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_mpi_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    log_debug(nnti_debug_level, "unregistering reg_buf(%p) buf(%p)", reg_buf, reg_buf->payload);

    del_target_buffer(reg_buf);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    while (!mpi_mem_hdl->wr_queue.empty()) {
        mpi_work_request *mpi_wr=NULL;

        mpi_wr=mpi_mem_hdl->wr_queue.front();
        mpi_mem_hdl->wr_queue.pop_front();

        log_debug(nnti_debug_level, "removing pending mpi_wr=%p mpi_wr->active_requests=%X", mpi_wr, (uint32_t)mpi_wr->active_requests);

//        for (int i=0;i<MAX_INDEX;i++) {
//        	if (mpi_wr->active_requests & (1<<i)) {
//        		// this request is active.  cancel it.
//        		log_debug(nnti_debug_level, "canceling mpi_wr->request[%d]=%p", i, &mpi_wr->request[i]);
//        		nthread_lock(&nnti_mpi_lock);
//        		MPI_Cancel(&mpi_wr->request[i]);
//        		nthread_unlock(&nnti_mpi_lock);
//        	}
//        }

//        del_wr_wrhash(mpi_wr);
//        if ((mpi_mem_hdl->type!=REQUEST_BUFFER) &&
//            (mpi_mem_hdl->type!=RECEIVE_BUFFER)) {
//            free(mpi_wr);
//        }
    }
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    if (mpi_mem_hdl)
        nthread_lock_fini(&mpi_mem_hdl->wr_queue_lock);
        delete mpi_mem_hdl;
    if (reg_buf->buffer_segments.NNTI_remote_addr_array_t_val)
        free(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val);

    reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
    MPI_SET_MATCH_ANY(&reg_buf->buffer_owner);
    reg_buf->ops               = (NNTI_buf_ops_t)0;
    reg_buf->payload_size      = 0;
    reg_buf->payload           = 0;
    reg_buf->transport_private = 0;

    log_debug(debug_level, "Finished unregistering, rc=%d",rc);

    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t NNTI_mpi_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;
    int                dest_rank;
    uint32_t           tag;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "msg_hdl",
                "NNTI_mpi_send", msg_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_hdl",
                "NNTI_mpi_send", dest_hdl);
    }

    mpi_mem_hdl=MPI_MEM_HDL(msg_hdl);
    assert(mpi_mem_hdl);
    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);

    mpi_wr->nnti_wr   =wr;
    mpi_wr->reg_buf   =(NNTI_buffer_t *)msg_hdl;
    mpi_wr->src_offset=0;
    mpi_wr->dst_offset=0;
    mpi_wr->length    =msg_hdl->payload_size;
    mpi_wr->op_state  =BUFFER_INIT;

    if (dest_hdl == NULL) {
        mpi_wr->peer   =*peer_hdl;
        mpi_wr->last_op=MPI_OP_SEND_REQUEST;

        dest_rank=peer_hdl->peer.NNTI_remote_process_t_u.mpi.rank;
        tag      =NNTI_MPI_REQUEST_TAG;
    } else {
        mpi_wr->peer   =dest_hdl->buffer_owner;
        mpi_wr->last_op=MPI_OP_SEND_BUFFER;

        dest_rank=dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;
        tag      =dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag;
    }

    log_debug(nnti_debug_level, "sending to (rank=%d, tag=%d)", dest_rank, tag);

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Issend(
            (char*)msg_hdl->payload,
            msg_hdl->payload_size,
            MPI_BYTE,
            dest_rank,
            tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[SEND_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to send with Isend");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    mpi_wr->request_ptr  =&mpi_wr->request[SEND_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->active_requests |= SEND_REQUEST_ACTIVE;

    wr->transport_id     =msg_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)msg_hdl;
    wr->ops              =NNTI_BOP_SEND_SRC;
    wr->transport_private=(uint64_t)mpi_wr;

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_mpi_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;
    int                dest_rank;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_mpi_put", src_buffer_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_mpi_put", dest_buffer_hdl);
    }

    mpi_mem_hdl=MPI_MEM_HDL(src_buffer_hdl);
    assert(mpi_mem_hdl);
    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);

    mpi_wr->nnti_wr   =wr;
    mpi_wr->reg_buf   =(NNTI_buffer_t *)src_buffer_hdl;
    mpi_wr->peer      =dest_buffer_hdl->buffer_owner;
    mpi_wr->src_offset=src_offset;
    mpi_wr->dst_offset=dest_offset;
    mpi_wr->length    =src_length;
    mpi_wr->op_state  =RDMA_WRITE_INIT;
    mpi_wr->last_op   =MPI_OP_PUT_INITIATOR;

    dest_rank=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;

    mpi_wr->rts_msg.length=src_length;
    mpi_wr->rts_msg.offset=dest_offset;
    mpi_wr->rts_msg.tag   =dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag;

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Issend(
            &mpi_wr->rts_msg,
            sizeof(mpi_wr->rts_msg),
            MPI_BYTE,
            dest_rank,
            dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.rts_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTS_REQ_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend RTS msg");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Issend(
            (char*)src_buffer_hdl->payload+src_offset,
            src_length,
            MPI_BYTE,
            dest_rank,
            dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[PUT_SEND_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend region");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    mpi_wr->request_ptr=&mpi_wr->request[RTS_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->active_requests |= RTS_REQUEST_ACTIVE;
    mpi_wr->active_requests |= PUT_SEND_REQUEST_ACTIVE;

    log_debug(nnti_debug_level, "putting to (%s, dest_rank=%d)", dest_buffer_hdl->buffer_owner.url, dest_rank);

    wr->transport_id     =src_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)src_buffer_hdl;
    wr->ops              =NNTI_BOP_LOCAL_READ;
    wr->transport_private=(uint64_t)mpi_wr;

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_mpi_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;
    int                src_rank;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_mpi_get", src_buffer_hdl);
    }
    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_mpi_get", dest_buffer_hdl);
    }

    log_debug(nnti_debug_level, "getting from (%s, src_offset=%llu, src_length=%llu, dest_offset=%llu)",
            src_buffer_hdl->buffer_owner.url, src_offset, src_length, dest_offset);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "src_buffer_hdl",
                "NNTI_mpi_get", src_buffer_hdl);
        fprint_NNTI_buffer(logger_get_file(), "dest_buffer_hdl",
                "NNTI_mpi_get", dest_buffer_hdl);
    }

    mpi_mem_hdl=MPI_MEM_HDL(dest_buffer_hdl);
    assert(mpi_mem_hdl);
    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);

    mpi_wr->nnti_wr   =wr;
    mpi_wr->reg_buf   =(NNTI_buffer_t *)dest_buffer_hdl;
    mpi_wr->peer      =src_buffer_hdl->buffer_owner;
    mpi_wr->src_offset=src_offset;
    mpi_wr->dst_offset=dest_offset;
    mpi_wr->length    =src_length;
    mpi_wr->op_state  =RDMA_READ_INIT;
    mpi_wr->last_op   =MPI_OP_GET_INITIATOR;

    src_rank=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;

    mpi_wr->rtr_msg.length=src_length;
    mpi_wr->rtr_msg.offset=src_offset;
    mpi_wr->rtr_msg.tag=dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag;

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Irecv(
            (char*)dest_buffer_hdl->payload+dest_offset,
            src_length,
            MPI_BYTE,
            src_rank,
            dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.data_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[GET_RECV_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Irecv region");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Issend(
            &mpi_wr->rtr_msg,
            sizeof(mpi_wr->rtr_msg),
            MPI_BYTE,
            src_rank,
            src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.mpi.rtr_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTR_REQ_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend RTR msg");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    mpi_wr->request_ptr=&mpi_wr->request[RTR_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->active_requests |= GET_RECV_REQUEST_ACTIVE;
    mpi_wr->active_requests |= RTR_REQUEST_ACTIVE;

    log_debug(nnti_debug_level, "getting from (%s, src_rank=%d)", src_buffer_hdl->buffer_owner.url, src_rank);

    wr->transport_id     =dest_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)dest_buffer_hdl;
    wr->ops              =NNTI_BOP_LOCAL_WRITE;
    wr->transport_private=(uint64_t)mpi_wr;

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * \param[in] src_buffer_hdl    A buffer containing the data to put.
 * \param[in] src_length        The number of bytes to put.
 * \param[in] dest_buffer_list  A list of buffers to put the data into.
 * \param[in] dest_count        The number of destination buffers.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_mpi_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *wr)
{
    return NNTI_ENOTSUP;
}


/**
 * @brief Transfer data from a peer.
 *
 * \param[in] src_buffer_list  A list of buffers containing the data to get.
 * \param[in] src_length       The number of bytes to get.
 * \param[in] src_count        The number of source buffers.
 * \param[in] dest_buffer_hdl  A buffer to get the data into.
 * \return A result code (NNTI_OK or an error)
 */
NNTI_result_t NNTI_mpi_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *wr)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_mpi_atomic_set_callback (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		NNTI_callback_fn_t      cbfunc,
		void                   *context)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_mpi_atomic_read (
		const NNTI_transport_t *trans_hdl,
		const uint64_t          local_atomic,
		int64_t                *value)
{
	nthread_lock(&transport_global_data.atomics[local_atomic].lock);
	*value = transport_global_data.atomics[local_atomic].value;
	nthread_unlock(&transport_global_data.atomics[local_atomic].lock);

    return NNTI_OK;
}


NNTI_result_t NNTI_mpi_atomic_fop (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           operand,
		const NNTI_atomic_op_t  op,
		NNTI_work_request_t    *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_work_request  *mpi_wr=NULL;
    int                dest_rank;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);

    mpi_wr->nnti_wr   =wr;
    mpi_wr->op_state  =BUFFER_INIT;

    mpi_wr->atomics_result_index=result_atomic;

    mpi_wr->peer   =*peer_hdl;
    mpi_wr->last_op=MPI_OP_FETCH_ADD;
    dest_rank      =peer_hdl->peer.NNTI_remote_process_t_u.mpi.rank;

    transport_global_data.atomics_request_msg.op         =MPI_ATOMIC_FETCH_ADD;
    transport_global_data.atomics_request_msg.index      =target_atomic;
    transport_global_data.atomics_request_msg.compare_add=operand;

    log_debug(nnti_debug_level, "sending fetch-add to (rank=%d)", dest_rank);

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Isend(
            (char*)&transport_global_data.atomics_request_msg,
            sizeof(transport_global_data.atomics_request_msg),
            MPI_BYTE,
            dest_rank,
            NNTI_MPI_ATOMICS_REQUEST_TAG,
            MPI_COMM_WORLD,
            &mpi_wr->request[ATOMICS_SEND_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to send with Isend");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Irecv(
            (char*)&transport_global_data.atomics_result_msg,
            sizeof(transport_global_data.atomics_result_msg),
            MPI_BYTE,
            dest_rank,
            NNTI_MPI_ATOMICS_RESULT_TAG,
            MPI_COMM_WORLD,
            &mpi_wr->request[ATOMICS_RECV_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to post recv with Irecv");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    mpi_wr->request_ptr  =&mpi_wr->request[ATOMICS_SEND_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->active_requests |= ATOMICS_SEND_REQUEST_ACTIVE;
    mpi_wr->active_requests |= ATOMICS_RECV_REQUEST_ACTIVE;

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_BOP_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)mpi_wr;

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


NNTI_result_t NNTI_mpi_atomic_cswap (
		const NNTI_transport_t *trans_hdl,
		const NNTI_peer_t      *peer_hdl,
		const uint64_t          target_atomic,
		const uint64_t          result_atomic,
		const int64_t           compare_operand,
		const int64_t           swap_operand,
		NNTI_work_request_t    *wr)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_work_request  *mpi_wr=NULL;
    int                dest_rank;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);

    mpi_wr->nnti_wr   =wr;
    mpi_wr->op_state  =BUFFER_INIT;

    mpi_wr->atomics_result_index=result_atomic;

    mpi_wr->peer   =*peer_hdl;
    mpi_wr->last_op=MPI_OP_FETCH_ADD;
    dest_rank      =peer_hdl->peer.NNTI_remote_process_t_u.mpi.rank;

    transport_global_data.atomics_request_msg.op         =MPI_ATOMIC_CMP_AND_SWP;
    transport_global_data.atomics_request_msg.index      =target_atomic;
    transport_global_data.atomics_request_msg.compare_add=compare_operand;
    transport_global_data.atomics_request_msg.swap       =swap_operand;

    log_debug(nnti_debug_level, "sending compare-swap to (rank=%d)", dest_rank);

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Isend(
            (char*)&transport_global_data.atomics_request_msg,
            sizeof(transport_global_data.atomics_request_msg),
            MPI_BYTE,
            dest_rank,
            NNTI_MPI_ATOMICS_REQUEST_TAG,
            MPI_COMM_WORLD,
            &mpi_wr->request[ATOMICS_SEND_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to send with Isend");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Irecv(
            (char*)&transport_global_data.atomics_result_msg,
            sizeof(transport_global_data.atomics_result_msg),
            MPI_BYTE,
            dest_rank,
            NNTI_MPI_ATOMICS_RESULT_TAG,
            MPI_COMM_WORLD,
            &mpi_wr->request[ATOMICS_RECV_INDEX]);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to post recv with Irecv");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    mpi_wr->request_ptr  =&mpi_wr->request[ATOMICS_SEND_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->active_requests |= ATOMICS_SEND_REQUEST_ACTIVE;
    mpi_wr->active_requests |= ATOMICS_RECV_REQUEST_ACTIVE;

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_BOP_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)mpi_wr;

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_mpi_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p ; wr=%p)", reg_buf, wr);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    wr->transport_id     =reg_buf->transport_id;
    wr->reg_buf          =reg_buf;
    wr->ops              =reg_buf->ops;
    wr->transport_private=(uint64_t)NULL;

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr=%p)", reg_buf, wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_mpi_clear_work_request (
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    wr->transport_private=(uint64_t)NULL;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_mpi_destroy_work_request (
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    wr->transport_id     =NNTI_TRANSPORT_NULL;
    wr->reg_buf          =NULL;
    wr->ops              =(NNTI_buf_ops_t)0;
    wr->transport_private=(uint64_t)NULL;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_mpi_cancel (
        NNTI_work_request_t *wr)
{
    return NNTI_ENOTSUP;
}


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_mpi_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    return NNTI_ENOTSUP;
}


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_mpi_interrupt (
        const NNTI_transport_t *trans_hdl)
{
    char dummy=0xAA;

    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
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
NNTI_result_t NNTI_mpi_wait (
        NNTI_work_request_t  *wr,
        const int             timeout,
        NNTI_status_t        *status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;

    long elapsed_time=0;
//    long timeout_per_call;
    MPI_Status event;
    int done=FALSE;

    log_level debug_level=nnti_debug_level;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(wr);
    assert(status);

	mpi_wr=MPI_WORK_REQUEST(wr);
    if (wr->ops != NNTI_BOP_ATOMICS) {
    	mpi_mem_hdl=MPI_MEM_HDL(wr->reg_buf);
    	assert(mpi_mem_hdl);
    	if (mpi_wr==NULL) {
    		nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    		mpi_wr=mpi_mem_hdl->wr_queue.front();
    		nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
    		wr->transport_private=(uint64_t)mpi_wr;
    	}
    }
	assert(mpi_wr);

    if (is_wr_complete(mpi_wr) == TRUE) {
        log_debug(debug_level, "work request already complete");
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "work request NOT complete");

//        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_atomic_operation();
            check_target_buffer_progress();

            log_debug(debug_level, "waiting on wr(%p) request(%p)", wr , mpi_wr->request_ptr);

            memset(&event, 0, sizeof(MPI_Status));
            done=FALSE;
            trios_start_timer(call_time);
            nthread_lock(&nnti_mpi_lock);
            rc = MPI_Testany(mpi_wr->request_count, mpi_wr->request_ptr, &mpi_wr->request_index, &done, &event);
            nthread_unlock(&nnti_mpi_lock);
            trios_stop_timer("NNTI_mpi_wait - MPI_Test", call_time);
            if ((rc==MPI_SUCCESS) && (mpi_wr->request_index==MPI_UNDEFINED) && (done==TRUE)) {
                log_debug(debug_level, "MPI_Testany() says there a no active requests (rc=%d, which_req=%d, done=%d)", rc, mpi_wr->request_index, done);
            }
            if (rc != MPI_SUCCESS) {
            	log_debug(debug_level, "MPI_Testany() return %d", rc);
            	abort();
            }
            log_debug(debug_level, "polling status is %d, which_req=%d, done=%d", rc, mpi_wr->request_index, done);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\tsource  = %d", event.MPI_SOURCE);
            log_debug(debug_level, "\ttag     = %d", event.MPI_TAG);
            log_debug(debug_level, "\terror   = %d", event.MPI_ERROR);
            log_debug(debug_level, "}");


            if (rc == MPI_SUCCESS) {
                /* case 1: success */
                if (done == TRUE) {
                    nnti_rc = NNTI_OK;

                }
                /* case 2: timed out */
                else {
                    elapsed_time = (trios_get_time_ms() - entry_time);

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout))) {
                        log_debug(debug_level, "MPI_Test() timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }

                    int timeout_remaining=timeout-elapsed_time;
                    if ((timeout < 0) || (timeout_remaining > MAX_SLEEP)) {
                    	nnti_sleep(MAX_SLEEP);
                    } else {
                    	if (timeout_remaining > 0) {
                    		nnti_sleep(timeout_remaining);
                    	}
                    }

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Test() failed (request=%p): rc=%d",
                        mpi_wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

            process_event(mpi_wr, &event);

            if (is_wr_complete(mpi_wr) == TRUE) {
                break;
            }
        }
    }

    create_status(wr, mpi_wr, nnti_rc, status);


    if (nnti_rc==NNTI_OK) {
        switch (mpi_wr->last_op) {
            case MPI_OP_NEW_REQUEST:
            case MPI_OP_RECEIVE:
                mpi_mem_hdl=MPI_MEM_HDL(mpi_wr->reg_buf);
                assert(mpi_mem_hdl);
                nthread_lock(&mpi_mem_hdl->wr_queue_lock);
                mpi_mem_hdl->wr_queue.pop_front();
                nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
                repost_recv_work_request(mpi_wr);
                break;
            case MPI_OP_PUT_TARGET:
            case MPI_OP_GET_TARGET:
                if ((mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) && (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ)) {
                    repost_RTR_RTS_recv_work_request(mpi_wr);
                } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) {
                    repost_RTS_recv_work_request(mpi_wr);
                } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ) {
                    repost_RTR_recv_work_request(mpi_wr);
                }
                break;
            case MPI_OP_FETCH_ADD:
            case MPI_OP_COMPARE_SWAP:
                free(mpi_wr);
                break;
            default:
                log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu",
                        (uint64_t)status->offset, (uint64_t)status->length);
                uint32_t index=status->offset/status->length;
                log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu, index=%llu",
                        (uint64_t)status->offset, (uint64_t)status->length, (uint64_t)index);

                free(mpi_wr);
        }
    }

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_mpi_wait", status);
    }

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_mpi_wait", total_time);

    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on any buffer in <tt>wr_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on any buffer in <tt>wr_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in wr_list must be registered with the same transport.
 *   2) You can't wait on the request queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_mpi_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;

    uint32_t     mpi_request_count=0;
    MPI_Request *mpi_requests=NULL;
    uint32_t    *request_to_wr_index=NULL;
    MPI_Status   event;
    int done=FALSE;

    int which_req=0;

    long elapsed_time=0;
//    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_mpi_wait(wr_list[0], timeout, status);
        *which=0;
        goto cleanup;
    }

    for (uint32_t i=0;i<wr_count;i++) {
        if (wr_list[i] != NULL) {
        	mpi_mem_hdl=MPI_MEM_HDL(wr_list[i]->reg_buf);
        	assert(mpi_mem_hdl);
        	mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
        	if (mpi_wr==NULL) {
        		nthread_lock(&mpi_mem_hdl->wr_queue_lock);
        		mpi_wr=mpi_mem_hdl->wr_queue.front();
        		nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
        		wr_list[i]->transport_private=(uint64_t)mpi_wr;
        	}
        	assert(mpi_wr);
        }
    }

    if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
        log_debug(debug_level, "work request already complete (which=%u, wr_list[%d]=%p)", *which, *which, wr_list[*which]);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "work request NOT complete (wr_list=%p)", wr_list);

//        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_atomic_operation();
            check_target_buffer_progress();

            log_debug(debug_level, "waiting on wr_list(%p)", wr_list);

            /*
             * The list of MPI_Requests is recreated each time through this loop.  There is probably a better way.
             */
            /* find the number of MPI_Requests in wr_list */
            mpi_request_count=0;
            for (uint32_t i=0;i<wr_count;i++) {
                if (wr_list[i] != NULL) {
                    mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
                    if (mpi_wr != NULL) {
                        mpi_request_count += mpi_wr->request_count;
                    }
                }
            }
            log_debug(debug_level, "wr_list contains %lu MPI_Requests", mpi_request_count);

            mpi_requests=(MPI_Request *)realloc(mpi_requests, mpi_request_count * sizeof(MPI_Request));
            assert(mpi_requests);
            request_to_wr_index=(uint32_t *)realloc(request_to_wr_index, mpi_request_count * sizeof(uint32_t));
            assert(request_to_wr_index);

            uint32_t request_index=0;
            /* for each NNTI work request */
            for (uint32_t i=0;i<wr_count;i++) {
                if (wr_list[i] != NULL) {
                    mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
                    if (mpi_wr != NULL) {
                        for (uint32_t j=0;j<mpi_wr->request_count;j++) {
                            mpi_requests[request_index]=mpi_wr->request_ptr[j];
                            request_to_wr_index[request_index]=i;
                            request_index++;
                        }
                    }
                }
            }

            memset(&event, 0, sizeof(MPI_Status));
            done=FALSE;
            trios_start_timer(call_time);
            nthread_lock(&nnti_mpi_lock);
            rc = MPI_Testany(mpi_request_count, mpi_requests, &which_req, &done, &event);
            nthread_unlock(&nnti_mpi_lock);
            trios_stop_timer("NNTI_mpi_waitany - MPI_Testany", call_time);
            if ((rc==MPI_SUCCESS) && (which_req==MPI_UNDEFINED) && (done==TRUE)) {
                log_debug(debug_level, "MPI_Testany() says there a no active requests (rc=%d, which_req=%d, done=%d)", rc, which_req, done);
            }
            log_debug(debug_level, "polling status is %d, which_req=%d, done=%d", rc, which_req, done);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\tsource  = %d", event.MPI_SOURCE);
            log_debug(debug_level, "\ttag     = %d", event.MPI_TAG);
            log_debug(debug_level, "\terror   = %d", event.MPI_ERROR);
            log_debug(debug_level, "}");


            if (rc == MPI_SUCCESS) {
                /* case 1: success */
                if (done == TRUE) {
                    *which=request_to_wr_index[which_req];
                    nnti_rc = NNTI_OK;
                    log_debug(debug_level, "*which == %d", *which);
                }
                /* case 2: timed out */
                else {
                    elapsed_time = (trios_get_time_ms() - entry_time);

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout))) {
                        log_debug(debug_level, "MPI_Testany() timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }

                    int timeout_remaining=timeout-elapsed_time;
                    if ((timeout < 0) || (timeout_remaining > MAX_SLEEP)) {
                    	nnti_sleep(MAX_SLEEP);
                    } else {
                    	if (timeout_remaining > 0) {
                    		nnti_sleep(timeout_remaining);
                    	}
                    }

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Testany() failed (request=%p): rc=%d",
                        mpi_wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

            process_event(MPI_WORK_REQUEST(wr_list[*which]), &event);

            if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
                break;
            }
        }
    }


    create_status(wr_list[*which], MPI_WORK_REQUEST(wr_list[*which]), nnti_rc, status);


    if (nnti_rc==NNTI_OK) {
        mpi_wr=MPI_WORK_REQUEST(wr_list[*which]);
        assert(mpi_wr);
        mpi_mem_hdl=MPI_MEM_HDL(mpi_wr->reg_buf);
        assert(mpi_mem_hdl);

        switch (mpi_wr->last_op) {
            case MPI_OP_NEW_REQUEST:
                nthread_lock(&mpi_mem_hdl->wr_queue_lock);
                mpi_mem_hdl->wr_queue.pop_front();
                nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
                repost_recv_work_request(mpi_wr);
                break;
            case MPI_OP_RECEIVE:
                nthread_lock(&mpi_mem_hdl->wr_queue_lock);
                mpi_mem_hdl->wr_queue.pop_front();
                nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
                repost_recv_work_request(mpi_wr);
                break;
            case MPI_OP_PUT_TARGET:
            case MPI_OP_GET_TARGET:
                if ((mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) && (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ)) {
                    repost_RTR_RTS_recv_work_request(mpi_wr);
                } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) {
                    repost_RTS_recv_work_request(mpi_wr);
                } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ) {
                    repost_RTR_recv_work_request(mpi_wr);
                }
                break;
            case MPI_OP_FETCH_ADD:
            case MPI_OP_COMPARE_SWAP:
                free(mpi_wr);
                break;
            default:
                log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu",
                        (uint64_t)status->offset, (uint64_t)status->length);
                uint32_t index=status->offset/status->length;
                log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu, index=%llu",
                        (uint64_t)status->offset, (uint64_t)status->length, (uint64_t)index);

                free(mpi_wr);
        }
    }

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_mpi_wait", status);
    }

cleanup:
    if (mpi_requests != NULL)        free(mpi_requests);
    if (request_to_wr_index != NULL) free(request_to_wr_index);

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_mpi_waitany", total_time);

    return(nnti_rc);
}

/**
 * @brief Wait for <tt>remote_op</tt> on all buffers in <tt>wr_list</tt> to complete.
 *
 * Wait for <tt>remote_op</tt> on all buffers in <tt>wr_list</tt> to complete or timeout
 * waiting.  This is typically used to wait for a result or a bulk data
 * transfer.  The timeout is specified in milliseconds.  A timeout of <tt>-1</tt>
 * means wait forever.  A timeout of <tt>0</tt> means do not wait.
 *
 * Caveats:
 *   1) All buffers in wr_list must be registered with the same transport.
 *   2) You can't wait on the receive queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_mpi_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr=NULL;

    uint32_t     mpi_request_count=0;
    MPI_Request *mpi_requests=NULL;
    MPI_Status  *events=NULL;
    int done=FALSE;

    long elapsed_time=0;
//    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_mpi_wait(wr_list[0], timeout, status[0]);
        goto cleanup;
    }

    for (uint32_t i=0;i<wr_count;i++) {
        if (wr_list[i] != NULL) {
        	mpi_mem_hdl=MPI_MEM_HDL(wr_list[i]->reg_buf);
        	assert(mpi_mem_hdl);
        	mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
        	if (mpi_wr==NULL) {
        		nthread_lock(&mpi_mem_hdl->wr_queue_lock);
        		mpi_wr=mpi_mem_hdl->wr_queue.front();
        		nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
        		wr_list[i]->transport_private=(uint64_t)mpi_wr;
        	}
        	assert(mpi_wr);
        }
    }

    if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
        log_debug(debug_level, "all work requests already complete (wr_list=%p)", wr_list);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "all work_requests NOT complete (wr_list=%p)", wr_list);

//        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_atomic_operation();
            check_target_buffer_progress();

            log_debug(debug_level, "waiting on wr_list(%p)", wr_list);

            /*
             * The list of MPI_Requests is recreated each time through this loop.  There is probably a better way.
             */
            /* find the number of MPI_Requests in wr_list */
            mpi_request_count=0;
            for (uint32_t i=0;i<wr_count;i++) {
                if (wr_list[i] != NULL) {
                    mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
                    if (mpi_wr != NULL) {
                        mpi_request_count += mpi_wr->request_count;
                    }
                }
            }
            log_debug(debug_level, "wr_list (wr_count=%lu) contains %lu MPI_Requests", wr_count, mpi_request_count);

            mpi_requests=(MPI_Request *)realloc(mpi_requests, mpi_request_count * sizeof(MPI_Request));
            events      =(MPI_Status *)realloc(events, mpi_request_count * sizeof(MPI_Status));
            assert(mpi_requests);
            assert(events);
            uint32_t request_index=0;
            /* for each NNTI work request */
            for (uint32_t i=0;i<wr_count;i++) {
                if (wr_list[i] != NULL) {
                    mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
                    if (mpi_wr != NULL) {
                        for (uint32_t j=0;j<mpi_wr->request_count;j++) {
                            mpi_requests[request_index++]=mpi_wr->request_ptr[j];
                        }
                    }
                }
            }

            memset(events, 0, mpi_request_count*sizeof(MPI_Status));
            done=FALSE;
            trios_start_timer(call_time);
            nthread_lock(&nnti_mpi_lock);
            rc = MPI_Testall(mpi_request_count, mpi_requests, &done, events);
            nthread_unlock(&nnti_mpi_lock);
            trios_stop_timer("NNTI_mpi_waitall - MPI_Testall", call_time);
            log_debug(debug_level, "polling status is %d (rc=%d, done=%d)", rc, rc, done);

            if (done == TRUE) {
                for (uint32_t i=0;i<mpi_request_count;i++) {
                    log_debug(debug_level, "Poll Event[%d]= {", i);
                    log_debug(debug_level, "\tsource  = %d", events[i].MPI_SOURCE);
                    log_debug(debug_level, "\ttag     = %d", events[i].MPI_TAG);
                    log_debug(debug_level, "\terror   = %d", events[i].MPI_ERROR);
                    log_debug(debug_level, "}");
                }
            }

            if (rc == MPI_SUCCESS) {
                /* case 1: success */
                if (done == TRUE) {
                    nnti_rc = NNTI_OK;
                }
                /* case 2: timed out */
                else {
                    elapsed_time = (trios_get_time_ms() - entry_time);

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout))) {
                        log_debug(debug_level, "MPI_Testall() timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }

                    int timeout_remaining=timeout-elapsed_time;
                    if ((timeout < 0) || (timeout_remaining > MAX_SLEEP)) {
                    	nnti_sleep(MAX_SLEEP);
                    } else {
                    	if (timeout_remaining > 0) {
                    		nnti_sleep(timeout_remaining);
                    	}
                    }

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Testall() failed (request=%p): rc=%d",
                        mpi_wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

            for (uint32_t i=0;i<wr_count;i++) {
                log_debug(debug_level, "processing event #%lu of %lu", i, wr_count);
                process_event(MPI_WORK_REQUEST(wr_list[i]), &(events[i]));
            }

            if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
                break;
            }
        }
    }


    for (uint32_t i=0;i<wr_count;i++) {
        mpi_wr=MPI_WORK_REQUEST(wr_list[i]);
        assert(mpi_wr);
        mpi_mem_hdl=MPI_MEM_HDL(mpi_wr->reg_buf);
        assert(mpi_mem_hdl);

        create_status(wr_list[i], mpi_wr, nnti_rc, status[i]);

        if (nnti_rc == NNTI_OK) {
            switch (mpi_wr->last_op) {
                case MPI_OP_NEW_REQUEST:
                    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
                    mpi_mem_hdl->wr_queue.pop_front();
                    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
                    repost_recv_work_request(mpi_wr);
                    break;
                case MPI_OP_RECEIVE:
                    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
                    mpi_mem_hdl->wr_queue.pop_front();
                    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
                    repost_recv_work_request(mpi_wr);
                    break;
                case MPI_OP_PUT_TARGET:
                case MPI_OP_GET_TARGET:
                    if ((mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) && (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ)) {
                        repost_RTR_RTS_recv_work_request(mpi_wr);
                    } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_WRITE) {
                        repost_RTS_recv_work_request(mpi_wr);
                    } else if (mpi_wr->nnti_wr->ops & NNTI_BOP_REMOTE_READ) {
                        repost_RTR_recv_work_request(mpi_wr);
                    }
                    break;
                case MPI_OP_FETCH_ADD:
                case MPI_OP_COMPARE_SWAP:
                    free(mpi_wr);
                    break;
                default:
                    log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu",
                            (uint64_t)status[i]->offset, (uint64_t)status[i]->length);
                    uint32_t index=status[i]->offset/status[i]->length;
                    log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu, index=%llu",
                            (uint64_t)status[i]->offset, (uint64_t)status[i]->length, (uint64_t)index);

                    free(mpi_wr);
            }
        }

        if (logging_debug(debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status[i]",
                    "end of NNTI_mpi_wait", status[i]);
        }
    }

cleanup:
    if (mpi_requests != NULL)        free(mpi_requests);
    if (events != NULL)              free(events);

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_mpi_waitall", total_time);

    return(nnti_rc);
}

/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_mpi_fini (
        const NNTI_transport_t *trans_hdl)
{
    nthread_counter_fini(&transport_global_data.mbits);

    nthread_lock_fini(&nnti_mpi_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&nnti_wr_wrhash_lock);
    nthread_lock_fini(&nnti_target_buffer_queue_lock);

    if (transport_global_data.init_called_mpi_init) {
    	MPI_Finalize();
    }

    return(NNTI_OK);
}



static NNTI_result_t setup_atomics(void)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    uint32_t atomics_bytes;

    log_debug(nnti_debug_level, "enter");

    atomics_bytes=config.min_atomics_vars * sizeof(mpi_atomic_t);
    trios_start_timer(callTime);
    transport_global_data.atomics=(mpi_atomic_t*)malloc(atomics_bytes);
    if (transport_global_data.atomics == NULL) {
    	return(NNTI_ENOMEM);
    }
    memset(transport_global_data.atomics, 0, atomics_bytes);
    trios_stop_timer("malloc and memset", callTime);

    trios_start_timer(callTime);
    for (int i=0;i<config.min_atomics_vars;i++) {
    	nthread_lock_init(&transport_global_data.atomics[i].lock);
    	transport_global_data.atomics[i].value=0;
    }
    trios_stop_timer("init locks", callTime);

    post_atomics_recv_request();

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}




static int check_target_buffer_progress()
{
    int nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *mpi_wr   =NULL;

    int rc=MPI_SUCCESS;

    MPI_Status event;
    int        which_req=0;
    int        done=FALSE;


    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    NNTI_buffer_t *reg_buf=NULL;
    target_buffer_queue_iter_t i;

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    nthread_lock(&nnti_target_buffer_queue_lock);
    if (target_buffers.size() == 0) {
        log_debug(debug_level, "there are no registered target buffers.  we're done here.");
    }
    for (i=target_buffers.begin(); i != target_buffers.end(); i++) {
        reg_buf = *i;

        mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
        assert(mpi_mem_hdl);
        nthread_lock(&mpi_mem_hdl->wr_queue_lock);
        mpi_wr=mpi_mem_hdl->wr_queue.front();
        nthread_unlock(&mpi_mem_hdl->wr_queue_lock);
        log_debug(debug_level, "checking reg_buf(%p) mpi_wr(%p)", reg_buf, mpi_wr);
        if (mpi_wr==NULL) {
            log_debug(debug_level, "there are no work requests posted.  continuing.");
            continue;
        }
        if (is_wr_complete(mpi_wr) == TRUE) {
            log_debug(debug_level, "this work request is complete.  continuing.");
            continue;
        }

        log_debug(debug_level, "testing mpi_wr->reg_buf(%p) mpi_wr->request_ptr(%p)", mpi_wr->reg_buf , mpi_wr->request_ptr);

        memset(&event, 0, sizeof(MPI_Status));
        done=FALSE;
        trios_start_timer(call_time);
        nthread_lock(&nnti_mpi_lock);
        rc = MPI_Testany(mpi_wr->request_count, mpi_wr->request_ptr, &which_req, &done, &event);
        nthread_unlock(&nnti_mpi_lock);
        trios_stop_timer("check_target_buffer_progress - MPI_Test", call_time);
        if ((rc==MPI_SUCCESS) && (which_req==MPI_UNDEFINED) && (done==TRUE)) {
            log_debug(debug_level, "MPI_Testany() says there a no active requests (rc=%d, which_req=%d, done=%d)", rc, which_req, done);
        }
        log_debug(debug_level, "polling status is %d, which_req=%d, done=%d", rc, which_req, done);

        log_debug(debug_level, "Poll Event= {");
        log_debug(debug_level, "\tsource  = %d", event.MPI_SOURCE);
        log_debug(debug_level, "\ttag     = %d", event.MPI_TAG);
        log_debug(debug_level, "\terror   = %d", event.MPI_ERROR);
        log_debug(debug_level, "}");

        if (rc == MPI_SUCCESS) {
            /* case 1: success */
            if (done == TRUE) {
                nnti_rc = NNTI_OK;

                if (event.MPI_TAG == mpi_mem_hdl->rts_tag) {
                    mpi_wr->last_op = MPI_OP_PUT_TARGET;
                } else if (event.MPI_TAG == mpi_mem_hdl->rtr_tag) {
                    mpi_wr->last_op = MPI_OP_GET_TARGET;
                } else {
                }
            } else {
                continue;
            }
        }
        /* MPI_Test failure */
        else {
            log_error(debug_level, "MPI_Test() failed (request=%p): rc=%d",
                    mpi_wr->request_ptr, rc);
            nnti_rc = NNTI_EIO;
            break;
        }

        process_event(mpi_wr, &event);

        nthread_lock(&mpi_mem_hdl->wr_queue_lock);
        wr_queue_iter_t victim=find(mpi_mem_hdl->wr_queue.begin(), mpi_mem_hdl->wr_queue.end(), mpi_wr);
        if (victim != mpi_mem_hdl->wr_queue.end()) {
            log_debug(debug_level, "erasing mpi_wr(%p) from wr_queue", mpi_wr);
        	mpi_mem_hdl->wr_queue.erase(victim);
        }
        nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

        if ((reg_buf->ops & NNTI_BOP_REMOTE_WRITE) && (reg_buf->ops & NNTI_BOP_REMOTE_READ)) {
            repost_RTR_RTS_recv_work_request(mpi_wr);
        } else if (reg_buf->ops & NNTI_BOP_REMOTE_WRITE) {
            repost_RTS_recv_work_request(mpi_wr);
        } else if (reg_buf->ops & NNTI_BOP_REMOTE_READ) {
            repost_RTR_recv_work_request(mpi_wr);
        }
    }
    nthread_unlock(&nnti_target_buffer_queue_lock);

    trios_stop_timer("check_target_buffer_progress", total_time);

    log_debug(debug_level, "exit");

    return(nnti_rc);
}

static int check_atomic_operation(void)
{
    int nnti_rc=NNTI_OK;

    int rc=MPI_SUCCESS;

    MPI_Status event;
    int        done=FALSE;

    int atomics_index   =-1;
    mpi_atomic_t *atomic=NULL;

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    memset(&event, 0, sizeof(MPI_Status));
    done=FALSE;
    trios_start_timer(call_time);
    nthread_lock(&nnti_mpi_lock);
    rc = MPI_Test(&transport_global_data.atomics_recv_request, &done, &event);
    nthread_unlock(&nnti_mpi_lock);
    trios_stop_timer("check_atomic_operation - MPI_Test", call_time);
    log_debug(debug_level, "polling status is %d", rc);

    log_debug(debug_level, "Poll Event= {");
    log_debug(debug_level, "\tsource  = %d", event.MPI_SOURCE);
    log_debug(debug_level, "\ttag     = %d", event.MPI_TAG);
    log_debug(debug_level, "\terror   = %d", event.MPI_ERROR);
    log_debug(debug_level, "}");
    log_debug(debug_level, "done = %d", done);

    if (rc == MPI_SUCCESS) {
    	/* case 1: success */
    	if (done) {
    		nnti_rc = NNTI_OK;
    	} else {
    		goto cleanup;
    	}
    }
    /* MPI_Test failure */
    else {
    	log_error(debug_level, "MPI_Test(atomics_recv_request) failed: rc=%d", rc);
    	nnti_rc = NNTI_EIO;
    	goto cleanup;
    }

    atomics_index=transport_global_data.atomics_request_msg.index;
    atomic       =&transport_global_data.atomics[atomics_index];

    nthread_lock(&atomic->lock);

    switch (transport_global_data.atomics_request_msg.op) {
		case MPI_ATOMIC_FETCH_ADD:
			transport_global_data.atomics_result_msg.result = atomic->value;
			atomic->value += transport_global_data.atomics_request_msg.compare_add;
			break;
		case MPI_ATOMIC_CMP_AND_SWP:
			transport_global_data.atomics_result_msg.result=atomic->value;
			if (atomic->value == transport_global_data.atomics_request_msg.compare_add) {
				atomic->value = transport_global_data.atomics_request_msg.swap;
			}
			break;
		default:
	    	log_error(debug_level, "unknown atomic op: rc=%d", transport_global_data.atomics_request_msg.op);
			break;
    }

    nthread_lock(&nnti_mpi_lock);
    rc=MPI_Send(
            (char*)&transport_global_data.atomics_result_msg,
            sizeof(transport_global_data.atomics_result_msg),
            MPI_BYTE,
            event.MPI_SOURCE,
            NNTI_MPI_ATOMICS_RESULT_TAG,
            MPI_COMM_WORLD);
    nthread_unlock(&nnti_mpi_lock);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to send with Isend");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    nthread_unlock(&atomic->lock);

    post_atomics_recv_request();

cleanup:
    trios_stop_timer("check_atomic_operation", total_time);

    log_debug(debug_level, "exit");

    return(nnti_rc);
}


static int process_event(
        mpi_work_request *mpi_wr,
        const MPI_Status *event)
{
    int rc=NNTI_OK;
    NNTI_buffer_t     *reg_buf    =NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_level debug_level = nnti_debug_level;

    assert(mpi_wr);

    mpi_wr->last_event=*event;

    if ((mpi_wr->nnti_wr) && (mpi_wr->nnti_wr->ops == NNTI_BOP_ATOMICS)) {
    	if (mpi_wr->op_state == BUFFER_INIT) {
            log_debug(debug_level, "got NNTI_BOP_ATOMICS send completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = SEND_COMPLETE;
            mpi_wr->active_requests &= ~ATOMICS_SEND_REQUEST_ACTIVE;

            mpi_wr->request_ptr  =&mpi_wr->request[ATOMICS_RECV_INDEX];
            mpi_wr->request_count=1;

    	} else if (mpi_wr->op_state == SEND_COMPLETE) {
            log_debug(debug_level, "got NNTI_BOP_ATOMICS recv completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RECV_COMPLETE;
            mpi_wr->active_requests &= ~ATOMICS_RECV_REQUEST_ACTIVE;
    	}

    	nthread_lock(&transport_global_data.atomics[mpi_wr->atomics_result_index].lock);
    	transport_global_data.atomics[mpi_wr->atomics_result_index].value = transport_global_data.atomics_result_msg.result;
    	nthread_unlock(&transport_global_data.atomics[mpi_wr->atomics_result_index].lock);

        mpi_wr->nnti_wr->result=NNTI_OK;
        return NNTI_OK;
    }

    reg_buf=mpi_wr->reg_buf;
    assert(reg_buf);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    log_debug(debug_level, "mpi_wr=%p; mpi_wr->last_op=%d", mpi_wr, mpi_wr->last_op);
    if (mpi_wr->last_op == MPI_OP_SEND_REQUEST) {
        if (mpi_wr->op_state == BUFFER_INIT) {
            log_debug(debug_level, "got SEND_REQUEST completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = SEND_COMPLETE;
            mpi_wr->active_requests &= ~SEND_REQUEST_ACTIVE;
        }
    } else if (mpi_wr->last_op == MPI_OP_SEND_BUFFER) {
        if (mpi_wr->op_state == BUFFER_INIT) {
            log_debug(debug_level, "got SEND_BUFFER completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = SEND_COMPLETE;
            mpi_wr->active_requests &= ~SEND_REQUEST_ACTIVE;
        }
    } else if (mpi_wr->last_op == MPI_OP_NEW_REQUEST) {
        if (mpi_wr->op_state == BUFFER_INIT) {
            log_debug(debug_level, "got NEW REQUEST completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RECV_COMPLETE;
            mpi_wr->dst_offset = mpi_wr->request_index*mpi_wr->length;
        }
    } else if (mpi_wr->last_op == MPI_OP_RECEIVE) {
        if (mpi_wr->op_state == BUFFER_INIT) {
            log_debug(debug_level, "got RECEIVE completion - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RECV_COMPLETE;
        }
    } else if (mpi_wr->last_op == MPI_OP_PUT_INITIATOR) {
        if (mpi_wr->op_state == RDMA_WRITE_INIT) {
            log_debug(debug_level, "got put_src RTS completion (initiator) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_RTS_COMPLETE;
            mpi_wr->active_requests &= ~RTS_REQUEST_ACTIVE;

            mpi_wr->request_ptr  =&mpi_wr->request[PUT_SEND_INDEX];
            mpi_wr->request_count=1;

        } else if (mpi_wr->op_state == RDMA_RTS_COMPLETE) {
            log_debug(debug_level, "got put_src WRITE completion (initiator) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_WRITE_COMPLETE;
            mpi_wr->active_requests &= ~PUT_SEND_REQUEST_ACTIVE;
        }
    } else if (mpi_wr->last_op == MPI_OP_PUT_TARGET) {
        if ((mpi_wr->op_state == RDMA_WRITE_INIT) || (mpi_wr->op_state == RDMA_TARGET_INIT)) {
            log_debug(debug_level, "got put_dst RTS completion (target) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            if (mpi_wr->op_state == RDMA_TARGET_INIT) {
                nthread_lock(&nnti_mpi_lock);
                MPI_Cancel(&mpi_wr->request[RTR_REQ_INDEX]);
                mpi_wr->active_requests &= ~RTR_REQUEST_ACTIVE;
                nthread_unlock(&nnti_mpi_lock);
            }
            mpi_wr->op_state = RDMA_RTS_COMPLETE;
            mpi_wr->active_requests &= ~RTS_REQUEST_ACTIVE;

            log_debug(debug_level, "receiving data from PUT initiator - rank(%d) tag(%d) dst_offset(%llu) dst_length(%llu)",
                    event->MPI_SOURCE, mpi_wr->rts_msg.tag,
                    mpi_wr->rts_msg.offset, mpi_wr->rts_msg.length);
            nthread_lock(&nnti_mpi_lock);
            MPI_Irecv(
                    (char*)reg_buf->payload+mpi_wr->rts_msg.offset,
                    mpi_wr->rts_msg.length,
                    MPI_BYTE,
                    event->MPI_SOURCE,
                    mpi_mem_hdl->data_tag,
                    MPI_COMM_WORLD,
                    &mpi_wr->request[PUT_RECV_INDEX]);
            nthread_unlock(&nnti_mpi_lock);
            mpi_wr->request_ptr  =&mpi_wr->request[PUT_RECV_INDEX];
            mpi_wr->request_count=1;
            mpi_wr->active_requests |= PUT_RECV_REQUEST_ACTIVE;

            mpi_wr->dst_offset=mpi_wr->rts_msg.offset;
            mpi_wr->length    =mpi_wr->rts_msg.length;

        } else if (mpi_wr->op_state == RDMA_RTS_COMPLETE) {
            log_debug(debug_level, "got put_dst WRITE completion (target) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_WRITE_COMPLETE;
            mpi_wr->active_requests &= ~PUT_RECV_REQUEST_ACTIVE;
        }
    } else if (mpi_wr->last_op == MPI_OP_GET_INITIATOR) {
        if (mpi_wr->op_state == RDMA_READ_INIT) {
            log_debug(debug_level, "got get_dst RTR completion (initiator) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_RTR_COMPLETE;
            mpi_wr->active_requests &= ~RTR_REQUEST_ACTIVE;

            mpi_wr->request_ptr  =&mpi_wr->request[GET_RECV_INDEX];
            mpi_wr->request_count=1;

        } else if (mpi_wr->op_state == RDMA_RTR_COMPLETE) {
            log_debug(debug_level, "got get_dst READ completion (initiator) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_READ_COMPLETE;
            mpi_wr->active_requests &= ~GET_RECV_REQUEST_ACTIVE;
        }
    } else if (mpi_wr->last_op == MPI_OP_GET_TARGET) {
        if ((mpi_wr->op_state == RDMA_READ_INIT) || (mpi_wr->op_state == RDMA_TARGET_INIT)) {
            log_debug(debug_level, "got get_src RTR completion (target) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            if (mpi_wr->op_state == RDMA_TARGET_INIT) {
                nthread_lock(&nnti_mpi_lock);
                MPI_Cancel(&mpi_wr->request[RTS_REQ_INDEX]);
                mpi_wr->active_requests &= ~RTS_REQUEST_ACTIVE;
                nthread_unlock(&nnti_mpi_lock);
            }
            mpi_wr->op_state = RDMA_RTR_COMPLETE;
            mpi_wr->active_requests &= ~RTR_REQUEST_ACTIVE;

            log_debug(debug_level, "sending data to GET initiator - rank(%d) tag(%d) src_offset(%llu) src_length(%llu)",
                    event->MPI_SOURCE, mpi_wr->rtr_msg.tag,
                    mpi_wr->rtr_msg.offset, mpi_wr->rtr_msg.length);
            nthread_lock(&nnti_mpi_lock);
            MPI_Issend(
                    (char*)reg_buf->payload+mpi_wr->rtr_msg.offset,
                    mpi_wr->rtr_msg.length,
                    MPI_BYTE,
                    event->MPI_SOURCE,
                    mpi_wr->rtr_msg.tag,
                    MPI_COMM_WORLD,
                    &mpi_wr->request[GET_SEND_INDEX]);
            nthread_unlock(&nnti_mpi_lock);
            mpi_wr->request_ptr  =&mpi_wr->request[GET_SEND_INDEX];
            mpi_wr->request_count=1;
            mpi_wr->active_requests |= GET_SEND_REQUEST_ACTIVE;

            mpi_wr->src_offset=mpi_wr->rtr_msg.offset;
            mpi_wr->length    =mpi_wr->rtr_msg.length;
        } else if (mpi_wr->op_state == RDMA_RTR_COMPLETE) {
            log_debug(debug_level, "got get_src READ completion (target) - event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);

            mpi_wr->op_state = RDMA_READ_COMPLETE;
            mpi_wr->active_requests &= ~GET_SEND_REQUEST_ACTIVE;
        }
    } else {
        log_error(debug_level, "unrecognized event arrived from %d - tag %4d",
                event->MPI_SOURCE, event->MPI_TAG);
    }

    return (rc);
}

static NNTI_result_t post_recv_dst_work_request(
        NNTI_buffer_t    *reg_buf,
        int64_t           tag,
        uint64_t          length)
{
    mpi_work_request    *mpi_wr     =NULL;
    mpi_memory_handle   *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);
    mpi_wr->reg_buf          = reg_buf;
    mpi_wr->dst_offset       = 0;
    mpi_wr->length           = length;
    mpi_wr->tag              = tag;
    mpi_wr->request_ptr      = &mpi_wr->request[RECV_INDEX];
    mpi_wr->request_count    = 1;
    mpi_wr->op_state         = BUFFER_INIT;
    mpi_wr->active_requests |= RECV_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_RECV_DST) {
        mpi_wr->last_op=MPI_OP_RECEIVE;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_RECV_DST(32) : actual %u", reg_buf->ops);

    }

    log_debug(nnti_debug_level, "posting irecv (reg_buf=%p ; mpi_wr=%p ; request_ptr=%p, tag=%lld)", reg_buf, mpi_wr, mpi_wr->request_ptr, mpi_wr->tag);
    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            (char*)reg_buf->payload,
            length,
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_wr->tag,
            MPI_COMM_WORLD,
            mpi_wr->request_ptr);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t post_recv_queue_work_request(
        NNTI_buffer_t    *reg_buf,
        int64_t           tag,
        uint64_t          length,
        uint64_t          count,
        MPI_Request      *request_list)
{
    uint64_t i=0;
    mpi_work_request    *mpi_wr=NULL;
    mpi_memory_handle   *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr->reg_buf      =reg_buf;
    mpi_wr->dst_offset   =0;
    mpi_wr->length       =length;
    mpi_wr->tag          =tag;
    mpi_wr->request_ptr  =request_list;
    mpi_wr->request_count=count;
    mpi_wr->op_state     =BUFFER_INIT;

    for (i=0;i<count;i++) {
        uint32_t offset =i*length;

        if (reg_buf->ops==NNTI_BOP_RECV_QUEUE) {
            mpi_wr->last_op=MPI_OP_NEW_REQUEST;

        } else {
            log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_RECV_QUEUE(64) : actual %u", reg_buf->ops);

        }

        log_debug(nnti_debug_level, "posting irecv (reg_buf=%p ; mpi_wr=%p ; request_ptr=%p, tag=%lld)", reg_buf, mpi_wr, &mpi_wr->request_ptr[i], mpi_wr->tag);
        nthread_lock(&nnti_mpi_lock);
        MPI_Irecv(
                (char*)reg_buf->payload + offset,
                length,
                MPI_BYTE,
                MPI_ANY_SOURCE,
                mpi_wr->tag,
                MPI_COMM_WORLD,
                &(request_list[i]));
        nthread_unlock(&nnti_mpi_lock);
    }

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t post_RTR_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request  *mpi_wr     =NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);
    mpi_wr->reg_buf      =reg_buf;
    mpi_wr->request_ptr  =&mpi_wr->request[RTR_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTR_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_REMOTE_READ) {
        mpi_wr->last_op =MPI_OP_GET_TARGET;
        mpi_wr->op_state=RDMA_READ_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_READ(2) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rtr_msg,
            sizeof(mpi_wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            mpi_wr->request_ptr);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t post_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request *mpi_wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);
    mpi_wr->reg_buf      =reg_buf;
    mpi_wr->request_ptr  =&mpi_wr->request[RTS_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTS_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_REMOTE_WRITE) {
        mpi_wr->last_op =MPI_OP_PUT_TARGET;
        mpi_wr->op_state=RDMA_WRITE_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_WRITE(8) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rts_msg,
            sizeof(mpi_wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            mpi_wr->request_ptr);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t post_RTR_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request *mpi_wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(mpi_wr);
    mpi_wr->reg_buf      =reg_buf;
    mpi_wr->request_ptr  =&mpi_wr->request[RTR_REQ_INDEX];
    mpi_wr->request_count=2;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTR_REQUEST_ACTIVE;
    mpi_wr->active_requests |= RTS_REQUEST_ACTIVE;

    if ((reg_buf->ops & NNTI_BOP_REMOTE_WRITE) && (reg_buf->ops & NNTI_BOP_REMOTE_READ)) {
        mpi_wr->op_state=RDMA_TARGET_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_WRITE|NNTI_BOP_REMOTE_READ(10) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rtr_msg,
            sizeof(mpi_wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTR_REQ_INDEX]);

    MPI_Irecv(
            &mpi_wr->rts_msg,
            sizeof(mpi_wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTS_REQ_INDEX]);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t post_atomics_recv_request(void)
{
    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &transport_global_data.atomics_request_msg,
            sizeof(transport_global_data.atomics_request_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            NNTI_MPI_ATOMICS_REQUEST_TAG,
            MPI_COMM_WORLD,
            &transport_global_data.atomics_recv_request);
    nthread_unlock(&nnti_mpi_lock);

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}


static NNTI_result_t repost_recv_work_request(
        mpi_work_request *mpi_wr)
{
    NNTI_buffer_t     *reg_buf    =mpi_wr->reg_buf;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    assert(reg_buf);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr->op_state=BUFFER_INIT;
    mpi_wr->active_requests |= RECV_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_RECV_QUEUE) {
        mpi_wr->last_op=MPI_OP_NEW_REQUEST;

    } else if (reg_buf->ops==NNTI_BOP_RECV_DST) {
        mpi_wr->last_op=MPI_OP_RECEIVE;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_RECV_DST(32) or NNTI_BOP_RECV_QUEUE(64) : actual %u", reg_buf->ops);

    }

    log_debug(nnti_debug_level, "posting irecv (reg_buf=%p ; mpi_wr=%p ; request_ptr=%p, tag=%lld)", reg_buf, mpi_wr, mpi_wr->request_ptr, mpi_wr->tag);
    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            (char*)reg_buf->payload + (mpi_wr->request_index * mpi_wr->length),
            mpi_wr->length,
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_wr->tag,
            MPI_COMM_WORLD,
            &mpi_wr->request_ptr[mpi_wr->request_index]);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTR_recv_work_request(
        mpi_work_request *mpi_wr)
{
    NNTI_buffer_t     *reg_buf    =mpi_wr->reg_buf;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    assert(reg_buf);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr->request_ptr  =&mpi_wr->request[RTR_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTR_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_REMOTE_READ) {
        mpi_wr->last_op =MPI_OP_GET_TARGET;
        mpi_wr->op_state=RDMA_READ_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_READ(2) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rtr_msg,
            sizeof(mpi_wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            mpi_wr->request_ptr);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTS_recv_work_request(
        mpi_work_request *mpi_wr)
{
    NNTI_buffer_t     *reg_buf    =mpi_wr->reg_buf;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    assert(reg_buf);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr->request_ptr  =&mpi_wr->request[RTS_REQ_INDEX];
    mpi_wr->request_count=1;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTS_REQUEST_ACTIVE;

    if (reg_buf->ops==NNTI_BOP_REMOTE_WRITE) {
        mpi_wr->last_op =MPI_OP_PUT_TARGET;
        mpi_wr->op_state=RDMA_WRITE_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_WRITE(8) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rts_msg,
            sizeof(mpi_wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            mpi_wr->request_ptr);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTR_RTS_recv_work_request(
        mpi_work_request *mpi_wr)
{
    NNTI_buffer_t     *reg_buf    =mpi_wr->reg_buf;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    assert(reg_buf);
    mpi_mem_hdl=MPI_MEM_HDL(reg_buf);
    assert(mpi_mem_hdl);

    mpi_wr->request_ptr  =&mpi_wr->request[RTR_REQ_INDEX];
    mpi_wr->request_count=2;
    mpi_wr->op_state     =BUFFER_INIT;
    mpi_wr->active_requests |= RTR_REQUEST_ACTIVE;
    mpi_wr->active_requests |= RTS_REQUEST_ACTIVE;

    if ((reg_buf->ops & NNTI_BOP_REMOTE_WRITE) && (reg_buf->ops & NNTI_BOP_REMOTE_READ)) {
        mpi_wr->op_state=RDMA_TARGET_INIT;

    } else {
        log_warn(nnti_debug_level, "incompatible buffer ops (expected NNTI_BOP_REMOTE_WRITE|NNTI_BOP_REMOTE_READ(10) : actual %u", reg_buf->ops);

    }

    nthread_lock(&nnti_mpi_lock);
    MPI_Irecv(
            &mpi_wr->rtr_msg,
            sizeof(mpi_wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTR_REQ_INDEX]);

    MPI_Irecv(
            &mpi_wr->rts_msg,
            sizeof(mpi_wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            &mpi_wr->request[RTS_REQ_INDEX]);
    nthread_unlock(&nnti_mpi_lock);

    nthread_lock(&mpi_mem_hdl->wr_queue_lock);
    mpi_mem_hdl->wr_queue.push_back(mpi_wr);
    nthread_unlock(&mpi_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; mpi_wr->request_ptr=%p)", reg_buf, mpi_wr->request_ptr);

    return(NNTI_OK);
}


static int is_wr_complete(
        mpi_work_request *mpi_wr)
{
    int8_t rc=FALSE;

    assert(mpi_wr);

    switch (mpi_wr->last_op) {
        case  MPI_OP_PUT_INITIATOR:
            if (mpi_wr->op_state == RDMA_WRITE_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_GET_INITIATOR:
            if (mpi_wr->op_state == RDMA_READ_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_PUT_TARGET:
            if (mpi_wr->op_state == RDMA_WRITE_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_GET_TARGET:
            if (mpi_wr->op_state == RDMA_READ_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_SEND_REQUEST:
        case  MPI_OP_SEND_BUFFER:
            if (mpi_wr->op_state == SEND_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_NEW_REQUEST:
        case  MPI_OP_RECEIVE:
            if (mpi_wr->op_state == RECV_COMPLETE) {
                rc=TRUE;
            }
            break;
        case  MPI_OP_FETCH_ADD:
        case  MPI_OP_COMPARE_SWAP:
            if (mpi_wr->op_state == RECV_COMPLETE) {
                rc=TRUE;
            }
            break;
        default:
            break;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static int8_t is_any_wr_complete(
        mpi_work_request **wr_list,
        const uint32_t     wr_count,
        uint32_t          *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_wr_complete(
        mpi_work_request **wr_list,
        const uint32_t     wr_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_wr_complete(
        NNTI_work_request_t *wr)
{
    mpi_work_request *mpi_wr=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    mpi_wr=MPI_WORK_REQUEST(wr);
    assert(mpi_wr);

    return(is_wr_complete(mpi_wr));
}

static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<wr_count;i++) {
        if ((wr_list[i] != NULL) &&
            (is_wr_complete(wr_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}



static NNTI_result_t insert_target_buffer(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter (buf=%p)", buf);

    nthread_lock(&nnti_target_buffer_queue_lock);
    target_buffers.push_back(buf);
    nthread_unlock(&nnti_target_buffer_queue_lock);

    log_debug(nnti_debug_level, "bufhash buffer added (buf=%p)", buf);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}
static NNTI_buffer_t *del_target_buffer(NNTI_buffer_t *buf)
{
    NNTI_buffer_t             *found=NULL;
    target_buffer_queue_iter_t victim;

    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter (buf=%p)", buf);

    nthread_lock(&nnti_target_buffer_queue_lock);
    victim=find(target_buffers.begin(), target_buffers.end(), buf);
    if (victim != target_buffers.end()) {
        log_debug(debug_level, "buffer found.  erasing from target_buffers.");
        target_buffers.erase(victim);
    } else {
        log_debug(debug_level, "buffer NOT found");
    }
    nthread_unlock(&nnti_target_buffer_queue_lock);

    log_debug(nnti_debug_level, "exit");

    return(buf);
}
//static void print_target_buffer_deque()
//{
//    if (!logging_debug(nnti_debug_level)) {
//        return;
//    }
//
//    if (target_buffers.empty()) {
//        log_debug(nnti_debug_level, "target_buffers is empty");
//        return;
//    }
//
//    target_buffer_queue_iter_t i;
//    for (i=target_buffers.begin(); i != target_buffers.end(); i++) {
//        log_debug(nnti_debug_level, "target buffer (%p)", *i);
//    }
//}

static void create_status(
        NNTI_work_request_t  *wr,
        mpi_work_request     *mpi_wr,
        int                   nnti_rc,
        NNTI_status_t        *status)
{
    log_debug(nnti_debug_level, "enter");

    status->op    =wr->ops;
    status->result=(NNTI_result_t)nnti_rc;
    if (nnti_rc==NNTI_OK) {
        if (mpi_wr->reg_buf) {
        	status->start =mpi_wr->reg_buf->payload;
            status->length=mpi_wr->length;
        }
        status->result=(NNTI_result_t)nnti_rc;
        switch (mpi_wr->last_op) {
            case MPI_OP_PUT_INITIATOR:
            case MPI_OP_GET_TARGET:
            case MPI_OP_SEND_REQUEST:
            case MPI_OP_SEND_BUFFER:
                status->offset=mpi_wr->src_offset;
                create_peer(&status->src, transport_global_data.rank); // allocates url
                create_peer(&status->dest, mpi_wr->last_event.MPI_SOURCE); // allocates url
                break;
            case MPI_OP_GET_INITIATOR:
            case MPI_OP_PUT_TARGET:
            case MPI_OP_NEW_REQUEST:
            case MPI_OP_RECEIVE:
                status->offset=mpi_wr->dst_offset;
                create_peer(&status->src, mpi_wr->last_event.MPI_SOURCE); // allocates url
                create_peer(&status->dest, transport_global_data.rank); // allocates url
                break;
        }
    }

    log_debug(nnti_debug_level, "exit");
}

static void create_peer(NNTI_peer_t *peer, int rank)
{
    log_debug(nnti_debug_level, "enter");

    sprintf(peer->url, "mpi://%u/", rank);

    peer->peer.transport_id                     = NNTI_TRANSPORT_MPI;
    peer->peer.NNTI_remote_process_t_u.mpi.rank = rank;

    log_debug(nnti_debug_level, "exit");
}

static void config_init(nnti_mpi_config *c)
{
    c->min_atomics_vars    = 512;
}

static void config_get_from_env(nnti_mpi_config *c)
{
    char *env_str=NULL;

    // defaults
    c->min_atomics_vars    = 512;

    if ((env_str=getenv("TRIOS_NNTI_MIN_ATOMIC_VARS")) != NULL) {
        errno=0;
        uint32_t min_vars=strtoul(env_str, NULL, 0);
        if (errno == 0) {
            log_debug(nnti_debug_level, "setting c->min_atomics_vars to %lu", min_vars);
            c->min_atomics_vars=min_vars;
        } else {
            log_debug(nnti_debug_level, "TRIOS_NNTI_MIN_ATOMIC_VARS value conversion failed (%s).  using c->min_atomics_vars default.", strerror(errno));
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_MIN_ATOMIC_VARS is undefined.  using c->min_atomics_vars default");
    }
}
