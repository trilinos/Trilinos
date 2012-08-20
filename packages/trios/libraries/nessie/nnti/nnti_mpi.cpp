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

#include "nnti_mpi.h"
#include "nnti_utils.h"



/* if undefined, the ACK message is NOT sent to the RDMA target when
 * the RDMA op is complete.  this creates one-sided semantics for RDMA
 * ops.  in this mode, the target has no idea when the RDMA op is
 * complete and what data was addressed.  NNTI_wait() returns NNTI_EINVAL
 * if passed a target buffer.
 */
#undef USE_RDMA_TARGET_ACK
/* if defined, the RDMA initiator will send an ACK message to the RDMA
 * target when the RDMA op is complete.  the target process must wait
 * on the target buffer in order to get the ACK.  this creates two-sided
 * semantics for RDMA ops.   in this mode, when the wait returns the
 * the RDMA op is complete and status indicates what data was addressed.
 */
#define USE_RDMA_TARGET_ACK


#define NNTI_MPI_REQUEST_TAG  0x01


#define MPI_OP_PUT_INITIATOR  1
#define MPI_OP_GET_INITIATOR  2
#define MPI_OP_PUT_TARGET     3
#define MPI_OP_GET_TARGET     4
#define MPI_OP_SEND_REQUEST   5
#define MPI_OP_SEND_BUFFER    6
#define MPI_OP_NEW_REQUEST    7
#define MPI_OP_RECEIVE        8

typedef enum {
    REQUEST_BUFFER,
    RECEIVE_BUFFER,
    SEND_BUFFER,
    GET_SRC_BUFFER,
    GET_DST_BUFFER,
    PUT_SRC_BUFFER,
    PUT_DST_BUFFER,
    RDMA_TARGET_BUFFER,
    UNKNOWN_BUFFER
} mpi_buffer_type;

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

#define RTR_REQ_INDEX  0
#define RTS_REQ_INDEX  1
#define SEND_INDEX     2
#define GET_RECV_INDEX 3
#define GET_SEND_INDEX 4
#define PUT_RECV_INDEX 5
#define PUT_SEND_INDEX 6
typedef struct mpi_work_request {
    NNTI_buffer_t  *reg_buf;
    NNTI_peer_t     peer;
    uint64_t        src_offset;
    uint64_t        dst_offset;
    uint64_t        length;
    int32_t         tag;
    MPI_Request     request[7];
    MPI_Request    *request_ptr;
    uint32_t        request_count;

    mpi_ready_msg   rtr_msg;
    mpi_ready_msg   rts_msg;

    mpi_op_state_t  op_state;

    MPI_Status      last_event;
    uint8_t         last_op;
    uint8_t         is_last_op_complete;
} mpi_work_request;

typedef std::deque<mpi_work_request *>           wr_queue_t;
typedef std::deque<mpi_work_request *>::iterator wr_queue_iter_t;

typedef struct mpi_memory_handle {
    mpi_buffer_type type;

    int32_t rtr_tag;
    int32_t rts_tag;
    int32_t data_tag;

    wr_queue_t wr_queue;
} mpi_memory_handle;


#define NUM_REQ_QUEUES 2
typedef struct mpi_request_queue_handle {
    NNTI_buffer_t *reg_buf;

    /* incoming queues */
    char *req_queue;

    /* each message is no larger than req_size */
    int req_size;

    /* each md can recv reqs_per_queue messages */
    int req_count;

    MPI_Request *requests;

} mpi_request_queue_handle;

typedef struct mpi_transport_global {

    int  rank;
    int  size;
    char proc_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Comm nnti_comm;

    mpi_request_queue_handle req_queue;

} mpi_transport_global;



static nthread_lock_t nnti_mpi_lock;


//static const NNTI_buffer_t *decode_event_buffer(
//        const NNTI_buffer_t *wait_buf,
//        const MPI_Status    *event);
static int process_event(
        const NNTI_buffer_t *reg_buf,
        const MPI_Status    *event);
static int check_target_buffer_progress(void);
static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         tag,
        uint64_t        offset,
        uint64_t        length,
        MPI_Request    *request);
static NNTI_result_t post_RTR_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t post_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t post_RTR_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr);
static NNTI_result_t repost_RTR_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr);
static NNTI_result_t repost_RTS_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr);
static NNTI_result_t repost_RTR_RTS_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr);
static int is_wr_complete(
        mpi_work_request *wr);
static mpi_work_request *first_incomplete_wr(
        mpi_memory_handle *mpi_mem_hdl);
static int8_t is_wr_queue_empty(
        const NNTI_buffer_t *reg_buf);
static int is_buf_op_complete(
        const NNTI_buffer_t *reg_buf);
static int8_t is_any_buf_op_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        uint32_t             *which);
static int8_t is_all_buf_ops_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count);

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf);
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash);
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *buf);
static void print_bufhash_map(void);

static NNTI_result_t insert_wr_wrhash(mpi_work_request *);
static mpi_work_request *get_wr_wrhash(const uint32_t bufhash);
static mpi_work_request *del_wr_wrhash(mpi_work_request *);
static void print_wrhash_map(void);

static NNTI_result_t insert_target_buffer(NNTI_buffer_t *buf);
static NNTI_buffer_t *del_target_buffer(NNTI_buffer_t *buf);
static void print_target_buffer_deque();

static void create_status(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        int                   nnti_rc,
        NNTI_status_t        *status);
static void create_peer(
        NNTI_peer_t *peer,
        int          rank);
static void copy_peer(
        NNTI_peer_t *src,
        NNTI_peer_t *dest);


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


static mpi_transport_global transport_global_data;
static const int MIN_TIMEOUT = 0;  /* in milliseconds */

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


//    char transport[NNTI_URL_LEN];
//    char address[NNTI_URL_LEN];
//    char memdesc[NNTI_URL_LEN];
//    char *sep;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);


    if (!initialized) {

        nthread_lock_init(&nnti_mpi_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);
        nthread_lock_init(&nnti_wr_wrhash_lock);
        nthread_lock_init(&nnti_target_buffer_queue_lock);

        if (my_url != NULL) {
            log_error(nnti_debug_level,"The MPI transport does not accept a URL at init.  Ignoring URL.");
//            if ((rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
//                return(rc);
//            }
//            if (0!=strcmp(transport, "mpi")) {
//                return(NNTI_EINVAL);
//            }
//
//            if ((rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
//                return(rc);
//            }
//
//            if ((rc=nnti_url_get_memdesc(my_url, memdesc, NNTI_URL_LEN)) != NNTI_OK) {
//                return(rc);
//            }
//
//            nid=strtol(address, NULL, 0);
        }


        log_debug(nnti_debug_level, "initializing MPI");

        memset(&transport_global_data, 0, sizeof(mpi_transport_global));

        rc=MPI_Initialized(&mpi_initialized);
        if (rc) {
            log_fatal(nnti_debug_level,"MPI_Initialized() failed, %d", rc);
            abort();
        }

        if (mpi_initialized==FALSE) {
            int argc=0;
            char **argv=NULL;

            log_debug(nnti_debug_level, "initializing MPI library");
            rc=MPI_Init(&argc, &argv);
            if (rc) {
                log_fatal(nnti_debug_level,"MPI_Init() failed, %d", rc);
                abort();
            }
        }

        MPI_Comm_size(MPI_COMM_WORLD, &transport_global_data.size);
        MPI_Comm_rank(MPI_COMM_WORLD, &transport_global_data.rank);
        MPI_Get_processor_name(transport_global_data.proc_name, &name_len);

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "MPI Initialized: rank=%llu, size=%llu, proc_name=%s\n",
                    (unsigned long long)transport_global_data.rank,
                    (unsigned long long)transport_global_data.size,
                    transport_global_data.proc_name);
        }

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
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port", "luc://endpoint_id/"
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
//    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
//    char memdesc[NNTI_URL_LEN];
//    char *sep;

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

//        if ((rc=nnti_url_get_memdesc(url, memdesc, NNTI_URL_LEN)) != NNTI_OK) {
//            return(rc);
//        }

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
NNTI_result_t NNTI_mpi_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf)
{
//    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    static uint64_t mbits=0x111;

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

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)mpi_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (ops == NNTI_RECV_QUEUE) {
        mpi_mem_hdl->type=REQUEST_BUFFER;
        mpi_mem_hdl->data_tag = NNTI_MPI_REQUEST_TAG;

    } else if (ops == NNTI_RECV_DST) {
        mpi_mem_hdl->type=RECEIVE_BUFFER;
        mpi_mem_hdl->rtr_tag  = 0;
        mpi_mem_hdl->rts_tag  = 0;
        mpi_mem_hdl->data_tag = mbits++;

    } else if (ops == NNTI_SEND_SRC) {
        mpi_mem_hdl->type=SEND_BUFFER;
        mpi_mem_hdl->rtr_tag  = 0;
        mpi_mem_hdl->rts_tag  = 0;
        mpi_mem_hdl->data_tag = mbits++;
    } else {
        if (ops == NNTI_GET_DST) {
            mpi_mem_hdl->type=GET_DST_BUFFER;
        } else if (ops == NNTI_GET_SRC) {
            mpi_mem_hdl->type=GET_SRC_BUFFER;
        } else if (ops == NNTI_PUT_SRC) {
            mpi_mem_hdl->type=PUT_SRC_BUFFER;
        } else if (ops == NNTI_PUT_DST) {
            mpi_mem_hdl->type=PUT_DST_BUFFER;
        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            mpi_mem_hdl->type=RDMA_TARGET_BUFFER;
        } else {
            mpi_mem_hdl->type=UNKNOWN_BUFFER;
        }
        mpi_mem_hdl->rtr_tag  = mbits++;
        mpi_mem_hdl->rts_tag  = mbits++;
        mpi_mem_hdl->data_tag = mbits++;
    }

    reg_buf->buffer_addr.transport_id                  = NNTI_TRANSPORT_MPI;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.mpi.size = element_size;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.mpi.rtr_tag  = mpi_mem_hdl->rtr_tag;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.mpi.rts_tag  = mpi_mem_hdl->rts_tag;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag = mpi_mem_hdl->data_tag;

    if (ops == NNTI_RECV_QUEUE) {
        int32_t index=0;
        mpi_request_queue_handle *q_hdl=&transport_global_data.req_queue;

        q_hdl->reg_buf=reg_buf;

        q_hdl->req_queue=buffer;
        q_hdl->req_size =element_size;
        q_hdl->req_count=num_elements;
        q_hdl->requests =(MPI_Request*)calloc(q_hdl->req_count, sizeof(MPI_Request));

        /* initialize the buffer */
        memset(q_hdl->req_queue, 0, q_hdl->req_count*q_hdl->req_size);

        for (index=0; index<q_hdl->req_count; index++) {
            post_recv_work_request(
                    reg_buf,
                    NNTI_MPI_REQUEST_TAG,
                    (index*q_hdl->req_size),
                    q_hdl->req_size,
                    &q_hdl->requests[index]);
        }

    } else if (ops == NNTI_RECV_DST) {
        post_recv_work_request(
                reg_buf,
                mpi_mem_hdl->data_tag,
                0,
                element_size,
                NULL);

    } else if (ops == NNTI_PUT_DST) {
        post_RTS_recv_work_request(reg_buf);
        insert_target_buffer(reg_buf);

    } else if (ops == NNTI_GET_SRC) {
        post_RTR_recv_work_request(reg_buf);
        insert_target_buffer(reg_buf);

    } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
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

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;

    assert(mpi_mem_hdl);

    log_debug(nnti_debug_level, "unregistering reg_buf(%p) buf(%p)", reg_buf, reg_buf->payload);

    wr_queue_iter_t i;
    for (i=mpi_mem_hdl->wr_queue.begin(); i != mpi_mem_hdl->wr_queue.end(); i++) {
        assert(*i);
        if (is_wr_complete(*i) == FALSE) {
        }
    }

    while (!mpi_mem_hdl->wr_queue.empty()) {
        mpi_work_request *wr=mpi_mem_hdl->wr_queue.front();
        log_debug(nnti_debug_level, "removing pending wr=%p", wr);

        MPI_Cancel(wr->request_ptr);

        mpi_mem_hdl->wr_queue.pop_front();
        del_wr_wrhash(wr);
        free(wr);
    }

    del_target_buffer(reg_buf);

    if (mpi_mem_hdl) delete mpi_mem_hdl;

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
        const NNTI_buffer_t *dest_hdl)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;
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

    mpi_mem_hdl=(mpi_memory_handle *)msg_hdl->transport_private;
    assert(mpi_mem_hdl);
    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);

    wr->reg_buf=(NNTI_buffer_t *)msg_hdl;

    wr->op_state = BUFFER_INIT;

    if (dest_hdl == NULL) {
        wr->peer  =*peer_hdl;
        dest_rank = peer_hdl->peer.NNTI_remote_process_t_u.mpi.rank;
        tag       = NNTI_MPI_REQUEST_TAG;

        wr->last_op=MPI_OP_SEND_REQUEST;

    } else {
        wr->peer  = dest_hdl->buffer_owner;
        dest_rank = dest_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;
        tag       = dest_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag;

        wr->last_op=MPI_OP_SEND_BUFFER;

    }

    wr->src_offset =0;
    wr->dst_offset =0;
    wr->length     =msg_hdl->payload_size;

    log_debug(nnti_debug_level, "sending to (rank=%d, tag=%d)", dest_rank, tag);

    rc=MPI_Isend(
            (char*)msg_hdl->payload,
            msg_hdl->payload_size,
            MPI_BYTE,
            dest_rank,
            tag,
            MPI_COMM_WORLD,
            &wr->request[SEND_INDEX]);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to send with Isend");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    wr->request_ptr=&wr->request[SEND_INDEX];
    wr->request_count=1;

    mpi_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

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
        const uint64_t       dest_offset)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;
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

    mpi_mem_hdl=(mpi_memory_handle *)src_buffer_hdl->transport_private;
    assert(mpi_mem_hdl);
    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);

    wr->reg_buf   =(NNTI_buffer_t *)src_buffer_hdl;
    wr->peer      =dest_buffer_hdl->buffer_owner;
    wr->src_offset=src_offset;
    wr->dst_offset=dest_offset;
    wr->length    =src_length;

    wr->op_state = RDMA_WRITE_INIT;

    dest_rank=dest_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;

    wr->rts_msg.length=src_length;
    wr->rts_msg.offset=dest_offset;
    wr->rts_msg.tag   =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag;

    rc=MPI_Isend(
            &wr->rts_msg,
            sizeof(wr->rts_msg),
            MPI_BYTE,
            dest_rank,
            dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.rts_tag,
            MPI_COMM_WORLD,
            &wr->request[RTS_REQ_INDEX]);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend RTS msg");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    wr->request_ptr=&wr->request[RTS_REQ_INDEX];
    wr->request_count=1;

    rc=MPI_Isend(
            (char*)src_buffer_hdl->payload+src_offset,
            src_length,
            MPI_BYTE,
            dest_rank,
            dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag,
            MPI_COMM_WORLD,
            &wr->request[PUT_SEND_INDEX]);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend region");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    log_debug(nnti_debug_level, "putting to (%s, dest_rank=%d)", dest_buffer_hdl->buffer_owner.url, dest_rank);

    wr->last_op=MPI_OP_PUT_INITIATOR;

    mpi_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

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
        const uint64_t       dest_offset)
{
    int rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;
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

    mpi_mem_hdl=(mpi_memory_handle *)dest_buffer_hdl->transport_private;
    assert(mpi_mem_hdl);
    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);

    wr->reg_buf   =(NNTI_buffer_t *)dest_buffer_hdl;
    wr->peer      =src_buffer_hdl->buffer_owner;
    wr->src_offset=src_offset;
    wr->dst_offset=dest_offset;
    wr->length    =src_length;

    wr->op_state = RDMA_READ_INIT;

    src_rank=src_buffer_hdl->buffer_owner.peer.NNTI_remote_process_t_u.mpi.rank;

    wr->rtr_msg.length=src_length;
    wr->rtr_msg.offset=src_offset;
    wr->rtr_msg.tag=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag;

    rc=MPI_Irecv(
            (char*)dest_buffer_hdl->payload+dest_offset,
            src_length,
            MPI_BYTE,
            src_rank,
            dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.data_tag,
            MPI_COMM_WORLD,
            &wr->request[GET_RECV_INDEX]);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Irecv region");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    rc=MPI_Isend(
            &wr->rtr_msg,
            sizeof(wr->rtr_msg),
            MPI_BYTE,
            src_rank,
            src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.mpi.rtr_tag,
            MPI_COMM_WORLD,
            &wr->request[RTR_REQ_INDEX]);
    if (rc != MPI_SUCCESS) {
        log_error(nnti_debug_level, "failed to Isend RTR msg");
        nnti_rc = NNTI_EBADRPC;
        goto cleanup;
    }

    wr->request_ptr=&wr->request[RTR_REQ_INDEX];
    wr->request_count=1;

    log_debug(nnti_debug_level, "getting from (%s, src_rank=%d)", src_buffer_hdl->buffer_owner.url, src_rank);

    wr->last_op=MPI_OP_GET_INITIATOR;

    mpi_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(nnti_rc);
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
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;

//    const NNTI_buffer_t  *wait_buf=NULL;

    long elapsed_time=0;
    long timeout_per_call;
    MPI_Status event;
    int done=FALSE;

//    MPI_Request *reqs     =NULL;
//    uint32_t     req_count=0;
    int          which_req=0;

    log_level debug_level=nnti_debug_level;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(reg_buf);
    assert(status);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);
    wr=first_incomplete_wr(mpi_mem_hdl);
    assert(wr);

    if (is_buf_op_complete(reg_buf) == TRUE) {
        log_debug(debug_level, "buffer op already complete");
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "buffer op NOT complete");

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_target_buffer_progress();
            if (is_buf_op_complete(reg_buf) == TRUE) {
                log_debug(debug_level, "buffer op completed during check_target_buffer_progress()");
                nnti_rc = NNTI_OK;
                break;
            }

            log_debug(debug_level, "waiting on reg_buf(%p) request(%p)", reg_buf , wr->request_ptr);

            memset(&event, 0, sizeof(MPI_Status));
            done=FALSE;
            trios_start_timer(call_time);
            rc = MPI_Testany(wr->request_count, wr->request_ptr, &which_req, &done, &event);
            trios_stop_timer("NNTI_mpi_wait - MPI_Test", call_time);
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

//                    trios_start_timer(call_time);
//                    nnti_sleep(timeout_per_call);
//                    trios_stop_timer("NNTI_mpi_wait - nnti_sleep", call_time);
//
//                    elapsed_time += timeout_per_call;

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Test() failed (request=%p): rc=%d",
                        wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

            process_event(reg_buf, &event);

            if (is_buf_op_complete(reg_buf) == TRUE) {
                break;
            }
        }
    }

    create_status(reg_buf, remote_op, nnti_rc, status);


    if (nnti_rc==NNTI_OK) {
        mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
        if (mpi_mem_hdl->type == REQUEST_BUFFER) {
            wr=mpi_mem_hdl->wr_queue.front();
            mpi_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);

        } else if (mpi_mem_hdl->type == RECEIVE_BUFFER) {
            wr=mpi_mem_hdl->wr_queue.front();
            mpi_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);

        } else if (mpi_mem_hdl->type == PUT_DST_BUFFER) {
            wr=mpi_mem_hdl->wr_queue.front();
            mpi_mem_hdl->wr_queue.pop_front();
            repost_RTS_recv_work_request((NNTI_buffer_t *)reg_buf, wr);

        } else if (mpi_mem_hdl->type == GET_SRC_BUFFER) {
            wr=mpi_mem_hdl->wr_queue.front();
            mpi_mem_hdl->wr_queue.pop_front();
            repost_RTR_recv_work_request((NNTI_buffer_t *)reg_buf, wr);

        }
#if defined(USE_RDMA_TARGET_ACK)
        else if (mpi_mem_hdl->type == RDMA_TARGET_BUFFER) {
            wr=mpi_mem_hdl->wr_queue.front();
            mpi_mem_hdl->wr_queue.pop_front();
            repost_RTR_RTS_recv_work_request((NNTI_buffer_t *)reg_buf, wr);

        }
#endif
        else {
            assert(mpi_mem_hdl);
            wr=mpi_mem_hdl->wr_queue.front();
            assert(wr);

            log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu",
                    (uint64_t)status->offset, (uint64_t)status->length);
            uint32_t index=status->offset/status->length;
            log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu, index=%llu",
                    (uint64_t)status->offset, (uint64_t)status->length, (uint64_t)index);

            mpi_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            free(wr);
        }
    }


//    if (nnti_rc == NNTI_OK) {
//        mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
//        assert(mpi_mem_hdl);
//        wr=mpi_mem_hdl->wr_queue.front();
//        assert(wr);
//
//        log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu",
//                (uint64_t)status->offset, (uint64_t)status->length);
//        uint32_t index=status->offset/status->length;
//        log_debug(nnti_debug_level, "status->offset=%llu, status->length=%llu, index=%llu",
//                (uint64_t)status->offset, (uint64_t)status->length, (uint64_t)index);
//
//        mpi_mem_hdl->wr_queue.pop_front();
//        del_wr_wrhash(wr);
//        free(wr);
//
//        if (mpi_mem_hdl->type == REQUEST_BUFFER) {
//            mpi_request_queue_handle *q_hdl=&transport_global_data.req_queue;
//            post_recv_work_request(
//                    (NNTI_buffer_t *)reg_buf,
//                    NNTI_MPI_REQUEST_TAG,
//                    (index*q_hdl->req_size),
//                    q_hdl->req_size,
//                    &q_hdl->requests[index]);
//        }
//    }

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_mpi_wait", status);
    }

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_mpi_wait", total_time);

    return(nnti_rc);
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
 *   2) You can't wait on the request queue and RDMA buffers in the same call.  Will probably be fixed in the future.
 */
NNTI_result_t NNTI_mpi_waitany (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;

    MPI_Status   event;
    MPI_Request *requests=NULL;
    int which_req=-1;
    int done=FALSE;

//    int i=0;

//    const NNTI_buffer_t  *wait_buf=NULL;

    long elapsed_time=0;
    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(buf_list);
    assert(buf_count > 0);
    if (buf_count > 1) {
        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
        for (uint32_t i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                assert(((mpi_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_mpi_wait(buf_list[0], remote_op, timeout, status);
        *which=0;
        goto cleanup;
    }

    if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
        log_debug(debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "buffer op NOT complete (buf_list=%p)", buf_list);

        requests=(MPI_Request *)calloc(buf_count, sizeof(MPI_Request*));
        assert(requests);
        for (uint32_t i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                mpi_mem_hdl=(mpi_memory_handle *)buf_list[i]->transport_private;
                assert(mpi_mem_hdl);
                wr=mpi_mem_hdl->wr_queue.front();
                if (wr != NULL) {
                    requests[i]=*wr->request_ptr;
                }
            }
        }

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_target_buffer_progress();
            if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                log_debug(debug_level, "buffer op completed during check_target_buffer_progress() (which=%u, buf_list[%d]=%p)",
                        *which, *which, buf_list[*which]);
                nnti_rc = NNTI_OK;
                break;
            }

            log_debug(debug_level, "waiting on buf_list(%p)", buf_list);

            trios_start_timer(call_time);
            rc = MPI_Testany(buf_count, requests, &which_req, &done, &event);
            trios_stop_timer("NNTI_mpi_waitany - MPI_Testany", call_time);
            log_debug(debug_level, "polling status is %d", rc);

            log_debug(debug_level, "Poll Event= {");
            log_debug(debug_level, "\tsource  = %d", event.MPI_SOURCE);
            log_debug(debug_level, "\ttag     = %d", event.MPI_TAG);
            log_debug(debug_level, "}");


            if (rc == MPI_SUCCESS) {
                /* case 1: success */
                if (done == TRUE) {
                    *which=which_req;
                    nnti_rc = NNTI_OK;
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

//                    nnti_sleep(timeout_per_call);
//
//                    elapsed_time += timeout_per_call;

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Testany() failed (request=%p): rc=%d",
                        wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

//            wait_buf=decode_event_buffer(buf_list[0], &event);
            process_event(buf_list[*which], &event);

//            if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
//                break;
//            }
        }
    }


//    if ((rc!=PTL_OK) && (wr->last_event.ni_fail_type != PTL_NI_OK)) {
//        log_error(debug_level, "NI reported error: ni_fail_type=%s",
//                PtlNIFailStr(transport_global_data.ni_h, wr->last_event.ni_fail_type));
//        nnti_rc = NNTI_EIO;
//    }

    create_status(buf_list[*which], remote_op, nnti_rc, status);

    if (nnti_rc == NNTI_OK) {
        mpi_mem_hdl=(mpi_memory_handle *)buf_list[*which]->transport_private;
        assert(mpi_mem_hdl);
        wr=mpi_mem_hdl->wr_queue.front();
        assert(wr);
        mpi_mem_hdl->wr_queue.pop_front();
        del_wr_wrhash(wr);
        free(wr);
    }

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_mpi_wait", status);
    }

cleanup:
    if (requests != NULL) free(requests);

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_mpi_waitany", total_time);

    return(nnti_rc);
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
NNTI_result_t NNTI_mpi_waitall (
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status)
{
    int rc=MPI_SUCCESS;
    NNTI_result_t nnti_rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;

    MPI_Status  *events=NULL;
    MPI_Request *requests=NULL;
    int done=FALSE;

//    int i=0;

    long elapsed_time=0;
    long timeout_per_call;

    long entry_time=trios_get_time_ms();

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    assert(buf_list);
    assert(buf_count > 0);
    if (buf_count > 1) {
        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
        for (uint32_t i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                assert(((mpi_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_mpi_wait(buf_list[0], remote_op, timeout, status[0]);
        goto cleanup;
    }

    if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
        log_debug(debug_level, "all buffer ops already complete (buf_list=%p)", buf_list);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "all buffer ops NOT complete (buf_list=%p)", buf_list);

        events  =(MPI_Status *)calloc(buf_count, sizeof(MPI_Status*));
        requests=(MPI_Request *)calloc(buf_count, sizeof(MPI_Request*));
        assert(requests);
        for (uint32_t i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                mpi_mem_hdl=(mpi_memory_handle *)buf_list[i]->transport_private;
                assert(mpi_mem_hdl);
                wr=mpi_mem_hdl->wr_queue.front();
                if (wr != NULL) {
                    requests[i]=*wr->request_ptr;
                }
            }
        }

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                return NNTI_ECANCELED;
            }

            check_target_buffer_progress();
            if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                log_debug(debug_level, "all buffer ops completed during check_target_buffer_progress()");
                nnti_rc = NNTI_OK;
                break;
            }

            log_debug(debug_level, "waiting on buf_list(%p)", buf_list);

            trios_start_timer(call_time);
            rc = MPI_Testall(buf_count, requests, &done, events);
            trios_stop_timer("NNTI_mpi_waitall - MPI_Testall", call_time);
            log_debug(debug_level, "polling status is %d", rc);

            if (done == TRUE) {
                for (uint32_t i=0;i<buf_count;i++) {
                    log_debug(debug_level, "Poll Event[%d]= {", i);
                    log_debug(debug_level, "\tsource  = %d", events[i].MPI_SOURCE);
                    log_debug(debug_level, "\ttag     = %d", events[i].MPI_TAG);
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

//                    nnti_sleep(timeout_per_call);
//
//                    elapsed_time += timeout_per_call;

                    /* continue if the timeout has not expired */
                    /* log_debug(debug_level, "timedout... continuing"); */

                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Testall() failed (request=%p): rc=%d",
                        wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

//            wait_buf=decode_event_buffer(buf_list[0], &event);
            for (uint32_t i=0;i<buf_count;i++) {
                process_event(buf_list[i], &events[i]);
            }

            if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                break;
            }
        }
    }

//    mpi_mem_hdl=(mpi_memory_handle *)buf_list[0]->transport_private;
//    if ((rc!=PTL_OK) && (wr->last_event.ni_fail_type != PTL_NI_OK)) {
//        log_error(debug_level, "NI reported error: ni_fail_type=%s",
//                PtlNIFailStr(transport_global_data.ni_h, wr->last_event.ni_fail_type));
//        nnti_rc = NNTI_EIO;
//    }



    for (uint32_t i=0;i<buf_count;i++) {
        create_status(buf_list[i], remote_op, nnti_rc, status[i]);

        if (nnti_rc == NNTI_OK) {
            mpi_mem_hdl=(mpi_memory_handle *)buf_list[i]->transport_private;
            assert(mpi_mem_hdl);
            wr=mpi_mem_hdl->wr_queue.front();
            assert(wr);
            mpi_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            free(wr);
        }

        if (logging_debug(debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status[i]",
                    "end of NNTI_mpi_wait", status[i]);
        }
    }

cleanup:
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
//    PtlFini();
    nthread_lock_fini(&nnti_mpi_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&nnti_wr_wrhash_lock);
    nthread_lock_fini(&nnti_target_buffer_queue_lock);

    return(NNTI_OK);
}







//static mpi_work_request *decode_work_request(
//        const MPI_Status *event)
//{
//    log_level debug_level = nnti_debug_level;
//
//    const NNTI_buffer_t *event_buf=NULL;
//    mpi_memory_handle *mpi_mem_hdl=NULL;
//    mpi_work_request  *wr=NULL;
//    mpi_work_request  *debug_wr=NULL;
//
//    log_debug(debug_level, "enter");
//
//    event_buf=(NNTI_buffer_t *)event->md.user_ptr;
//    assert(event_buf);
//    mpi_mem_hdl=(mpi_memory_handle *)event_buf->transport_private;
//    assert(mpi_mem_hdl);
//
//    wr_queue_iter_t i;
//    for (i=mpi_mem_hdl->wr_queue.begin(); i != mpi_mem_hdl->wr_queue.end(); i++) {
//        assert(*i);
//        if (is_wr_complete(*i) == FALSE) {
//            // work request is incomplete, check if it matches this event
//            switch(mpi_mem_hdl->type) {
//                case REQUEST_BUFFER:
//                    if (((*i)->src_offset == event->offset) &&
//                        ((*i)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", *i);
//                    }
//                    break;
//                case SEND_BUFFER:
//                case PUT_SRC_BUFFER:
//                    if (((*i)->src_offset == event->offset) &&
//                        ((*i)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", *i);
//                    }
//                    break;
//                case GET_DST_BUFFER:
//                    if (((*i)->dst_offset == event->offset) &&
//                        ((*i)->length == event->mlength)) {
//
//                        wr=*i;
//                    } else {
//                        log_debug(debug_level, "work request doesn't match (wr=%p)", *i);
//                    }
//                    break;
//                case RECEIVE_BUFFER:
//                case GET_SRC_BUFFER:
//                case PUT_DST_BUFFER:
//                case RDMA_TARGET_BUFFER:
//                    wr=*i;
//                    break;
//                default:
//                    log_debug(debug_level, "unknown event type %d (event_buf==%p)", mpi_mem_hdl->type, event_buf);
//                    break;
//            }
//            if (wr) {
//                break;
//            }
//        } else {
//            log_debug(debug_level, "work request is already complete (wr=%p)", *i);
//        }
//    }
//
//    if (!wr) {
//        for (i=mpi_mem_hdl->wr_queue.begin(); i != mpi_mem_hdl->wr_queue.end(); i++) {
//            debug_wr=*i;
//            log_debug(LOG_ALL, "e.offset=%llu, e.rlength=%llu, e.mlength=%llu, wr=%p, wr.length=%llu, wr.src_offset=%llu, wr.dst_offset=%llu, wr.is_complete=%d",
//                    (uint64_t)event->offset, (uint64_t)event->rlength, (uint64_t)event->mlength,
//                    debug_wr,
//                    (uint64_t)debug_wr->length, (uint64_t)debug_wr->src_offset, (uint64_t)debug_wr->dst_offset,
//                    (is_wr_complete(debug_wr)==TRUE));
//        }
//    }
//    assert(wr);
//
//    log_debug(debug_level, "exit (wr==%p)", wr);
//
//    return(wr);
//}

//static const NNTI_buffer_t *decode_event_buffer(
//        const NNTI_buffer_t *wait_buf,
//        const MPI_Status    *event)
//{
//    const NNTI_buffer_t *event_buf=NULL;
//    mpi_memory_handle *mpi_mem_hdl=NULL;
//
//    log_debug(nnti_debug_level, "enter");
//
//    if ((wait_buf != NULL) && (wait_buf->transport_private != NULL) && (((mpi_memory_handle *)wait_buf->transport_private)->type == REQUEST_BUFFER)) {
//        event_buf=wait_buf;
//        mpi_mem_hdl=(mpi_memory_handle *)event_buf->transport_private;
//        assert(mpi_mem_hdl);
//
//        log_debug(nnti_debug_level, "the wait buffer is a REQUEST BUFFER, so event.md.user_ptr is the index of the request buffer.");
//    } else {
//        event_buf=(NNTI_buffer_t *)event->md.user_ptr;
//        mpi_mem_hdl=(mpi_memory_handle *)event_buf->transport_private;
//        assert(mpi_mem_hdl);
//
//        if (event_buf == wait_buf) {
//            log_debug(nnti_debug_level, "the wc matches the wait buffer (eq=%d, user_ptr=%p, wait_buf=%p)",
//                    mpi_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_buf);
//        } else {
//            log_debug(nnti_debug_level, "the wc does NOT match the wait buffer (eq=%d, user_ptr=%p, wait_buf=%p)",
//                    mpi_mem_hdl->eq_h, (void *)event->md.user_ptr, wait_buf);
//        }
//    }
//
//    log_debug(nnti_debug_level, "exit (event_buf==%p)", event_buf);
//
//    return(event_buf);
//}


static int check_target_buffer_progress()
{
    int nnti_rc=NNTI_OK;

    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;

    int rc=MPI_SUCCESS;

    MPI_Status event;
    int        which_req=0;
    int        done=FALSE;


    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    NNTI_buffer_t             *reg_buf=NULL;
    target_buffer_queue_iter_t i;

    nthread_lock(&nnti_target_buffer_queue_lock);
    for (i=target_buffers.begin(); i != target_buffers.end(); i++) {
        reg_buf = *i;

        if (is_buf_op_complete(reg_buf) == TRUE) {
            continue;
        } else {
            mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
            assert(mpi_mem_hdl);
            wr=first_incomplete_wr(mpi_mem_hdl);
            if (wr==NULL) {
                continue;
            }

            log_debug(debug_level, "testing reg_buf(%p) request(%p)", reg_buf , wr->request_ptr);

            memset(&event, 0, sizeof(MPI_Status));
            done=FALSE;
            trios_start_timer(call_time);
            rc = MPI_Testany(wr->request_count, wr->request_ptr, &which_req, &done, &event);
            trios_stop_timer("check_target_buffer_progress - MPI_Test", call_time);
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
                } else {
                    continue;
                }
            }
            /* MPI_Test failure */
            else {
                log_error(debug_level, "MPI_Test() failed (request=%p): rc=%d",
                        wr->request_ptr, rc);
                nnti_rc = NNTI_EIO;
                break;
            }

            process_event(reg_buf, &event);
        }
    }
    nthread_unlock(&nnti_target_buffer_queue_lock);

    trios_stop_timer("check_target_buffer_progress", total_time);

    return(nnti_rc);
}


static int process_event(
        const NNTI_buffer_t  *reg_buf,
        const MPI_Status     *event)
{
    int rc=NNTI_OK;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;
    log_level debug_level = nnti_debug_level;

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

//    if (mpi_mem_hdl->type != REQUEST_BUFFER) {
//        wr = decode_work_request(event);
//    } else {
//        wr=mpi_mem_hdl->wr_queue.front();
//    }

//    wr=mpi_mem_hdl->wr_queue.front();
    wr=first_incomplete_wr(mpi_mem_hdl);
    assert(wr);

    wr->last_event=*event;

    log_debug(debug_level, "reg_buf=%p; wr->last_op=%d", reg_buf, wr->last_op);
    switch (mpi_mem_hdl->type) {
        case SEND_BUFFER:
            if (wr->op_state == BUFFER_INIT) {
                log_debug(debug_level, "got SEND_BUFFER completion - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = SEND_COMPLETE;
            }
            break;

        case REQUEST_BUFFER:
            if (wr->op_state == BUFFER_INIT) {
                log_debug(debug_level, "got NEW REQUEST completion - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RECV_COMPLETE;
            }
            break;

        case RECEIVE_BUFFER:
            if (wr->op_state == BUFFER_INIT) {
                log_debug(debug_level, "got RECEIVE completion - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RECV_COMPLETE;
            }
            break;

        case PUT_SRC_BUFFER: // initiator
            if (wr->op_state == RDMA_WRITE_INIT) {
                log_debug(debug_level, "got put_src RTS completion (initiator) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_RTS_COMPLETE;

                wr->request_ptr=&wr->request[PUT_SEND_INDEX];
                wr->request_count=1;

            } else if (wr->op_state == RDMA_RTS_COMPLETE) {
                log_debug(debug_level, "got put_src WRITE completion (initiator) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_WRITE_COMPLETE;
            }
            break;

        case PUT_DST_BUFFER: // target
            if (wr->op_state == RDMA_WRITE_INIT) {
                log_debug(debug_level, "got put_dst RTS completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_RTS_COMPLETE;

                log_debug(debug_level, "receiving data from PUT initiator - rank(%d) tag(%d) dst_offset(%llu) dst_length(%llu)",
                        event->MPI_SOURCE, wr->rts_msg.tag,
                        wr->rts_msg.offset, wr->rts_msg.length);
                MPI_Irecv(
                        (char*)reg_buf->payload+wr->rts_msg.offset,
                        wr->rts_msg.length,
                        MPI_BYTE,
                        event->MPI_SOURCE,
                        mpi_mem_hdl->data_tag,
                        MPI_COMM_WORLD,
                        &wr->request[PUT_RECV_INDEX]);
                wr->request_ptr=&wr->request[PUT_RECV_INDEX];
                wr->request_count=1;

                wr->dst_offset=wr->rts_msg.offset;
                wr->length    =wr->rts_msg.length;

            } else if (wr->op_state == RDMA_RTS_COMPLETE) {
                log_debug(debug_level, "got put_dst WRITE completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_WRITE_COMPLETE;
            }
            break;

        case GET_DST_BUFFER: // initiator
            if (wr->op_state == RDMA_READ_INIT) {
                log_debug(debug_level, "got get_dst RTR completion (initiator) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_RTR_COMPLETE;

                wr->request_ptr = &wr->request[GET_RECV_INDEX];
                wr->request_count=1;

            } else if (wr->op_state == RDMA_RTR_COMPLETE) {
                log_debug(debug_level, "got get_dst READ completion (initiator) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_READ_COMPLETE;
            }
            break;

        case GET_SRC_BUFFER: // target
            if (wr->op_state == RDMA_READ_INIT) {
                log_debug(debug_level, "got get_src RTR completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_RTR_COMPLETE;

                log_debug(debug_level, "sending data to GET initiator - rank(%d) tag(%d) src_offset(%llu) src_length(%llu)",
                        event->MPI_SOURCE, wr->rtr_msg.tag,
                        wr->rtr_msg.offset, wr->rtr_msg.length);
                MPI_Isend(
                        (char*)reg_buf->payload+wr->rtr_msg.offset,
                        wr->rtr_msg.length,
                        MPI_BYTE,
                        event->MPI_SOURCE,
                        wr->rtr_msg.tag,
                        MPI_COMM_WORLD,
                        &wr->request[GET_SEND_INDEX]);
                wr->request_ptr=&wr->request[GET_SEND_INDEX];
                wr->request_count=1;

                wr->src_offset=wr->rtr_msg.offset;
                wr->length    =wr->rtr_msg.length;
            } else if (wr->op_state == RDMA_RTR_COMPLETE) {
                log_debug(debug_level, "got get_src READ completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_READ_COMPLETE;
            }
            break;

        case RDMA_TARGET_BUFFER:
            if (wr->op_state == RDMA_TARGET_INIT) {
                if (event->MPI_TAG == mpi_mem_hdl->rts_tag) {
                    log_debug(debug_level, "got target_buffer RTS completion (target) - event arrived from %d - tag %4d",
                            event->MPI_SOURCE, event->MPI_TAG);
                    wr->op_state = RDMA_RTS_COMPLETE;

                    MPI_Cancel(&wr->request[RTR_REQ_INDEX]);

                    log_debug(debug_level, "receiving data from PUT initiator - rank(%d) tag(%d) dst_offset(%llu) dst_length(%llu)",
                            event->MPI_SOURCE, wr->rts_msg.tag,
                            wr->rts_msg.offset, wr->rts_msg.length);
                    MPI_Irecv(
                            (char*)reg_buf->payload+wr->rts_msg.offset,
                            wr->rts_msg.length,
                            MPI_BYTE,
                            event->MPI_SOURCE,
                            mpi_mem_hdl->data_tag,
                            MPI_COMM_WORLD,
                            &wr->request[PUT_RECV_INDEX]);
                    wr->request_ptr=&wr->request[PUT_RECV_INDEX];
                    wr->request_count=1;

                    wr->dst_offset=wr->rts_msg.offset;
                    wr->length    =wr->rts_msg.length;

                } else if (event->MPI_TAG == mpi_mem_hdl->rtr_tag) {
                    log_debug(debug_level, "got target_buffer RTR completion (target) - event arrived from %d - tag %4d",
                            event->MPI_SOURCE, event->MPI_TAG);
                    wr->op_state = RDMA_RTR_COMPLETE;

                    MPI_Cancel(&wr->request[RTS_REQ_INDEX]);

                    log_debug(debug_level, "sending data to GET initiator - rank(%d) tag(%d) src_offset(%llu) src_length(%llu)",
                            event->MPI_SOURCE, wr->rtr_msg.tag,
                            wr->rtr_msg.offset, wr->rtr_msg.length);
                    MPI_Isend(
                            (char*)reg_buf->payload+wr->rtr_msg.offset,
                            wr->rtr_msg.length,
                            MPI_BYTE,
                            event->MPI_SOURCE,
                            wr->rtr_msg.tag,
                            MPI_COMM_WORLD,
                            &wr->request[GET_SEND_INDEX]);
                    wr->request_ptr=&wr->request[GET_SEND_INDEX];
                    wr->request_count=1;

                    wr->src_offset=wr->rtr_msg.offset;
                    wr->length    =wr->rtr_msg.length;
                } else {

                }

            } else if (wr->op_state == RDMA_RTS_COMPLETE) {
                log_debug(debug_level, "got target_buffer WRITE completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_TARGET_COMPLETE;

            } else if (wr->op_state == RDMA_RTR_COMPLETE) {
                log_debug(debug_level, "got target_buffer READ completion (target) - event arrived from %d - tag %4d",
                        event->MPI_SOURCE, event->MPI_TAG);
                wr->op_state = RDMA_TARGET_COMPLETE;
            }
            break;

        case UNKNOWN_BUFFER:
        default:
            log_error(debug_level, "unrecognized event arrived from %d - tag %4d",
                    event->MPI_SOURCE, event->MPI_TAG);
            break;
    }

    return (rc);
}

static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         tag,
        uint64_t        offset,
        uint64_t        length,
        MPI_Request    *request)
{
    mpi_work_request *wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);
    wr->reg_buf    = reg_buf;
    wr->dst_offset = offset;
    wr->length     = length;
    wr->tag        = tag;
    if (request != NULL) {
        wr->request_ptr = request;
    } else {
        wr->request_ptr = &wr->request[0];
    }
    wr->request_count=1;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==REQUEST_BUFFER) {
        wr->last_op=MPI_OP_NEW_REQUEST;

    } else if (mpi_mem_hdl->type==RECEIVE_BUFFER) {
        wr->last_op=MPI_OP_RECEIVE;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    log_debug(nnti_debug_level, "posting irecv (reg_buf=%p ; wr=%p ; request_ptr=%p, tag=%lld)", reg_buf, wr, wr->request_ptr, wr->tag);
    MPI_Irecv(
            (char*)reg_buf->payload + offset,
            length,
            MPI_BYTE,
            MPI_ANY_SOURCE,
            wr->tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t post_RTR_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request *wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);
    wr->reg_buf      =reg_buf;
    wr->request_ptr  =&wr->request[RTR_REQ_INDEX];
    wr->request_count=1;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rtr_msg,
            sizeof(wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t post_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request *wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);
    wr->reg_buf      =reg_buf;
    wr->request_ptr  =&wr->request[RTS_REQ_INDEX];
    wr->request_count=1;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rts_msg,
            sizeof(wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t post_RTR_RTS_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    mpi_work_request *wr=NULL;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr=(mpi_work_request *)calloc(1, sizeof(mpi_work_request));
    assert(wr);
    wr->reg_buf      =reg_buf;
    wr->request_ptr  =&wr->request[RTR_REQ_INDEX];
    wr->request_count=2;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rtr_msg,
            sizeof(wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            &wr->request[RTR_REQ_INDEX]);

    MPI_Irecv(
            &wr->rts_msg,
            sizeof(wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            &wr->request[RTS_REQ_INDEX]);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}


static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==REQUEST_BUFFER) {
        wr->last_op=MPI_OP_NEW_REQUEST;

    } else if (mpi_mem_hdl->type==RECEIVE_BUFFER) {
        wr->last_op=MPI_OP_RECEIVE;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    log_debug(nnti_debug_level, "posting irecv (reg_buf=%p ; wr=%p ; request_ptr=%p, tag=%lld)", reg_buf, wr, wr->request_ptr, wr->tag);
    MPI_Irecv(
            (char*)reg_buf->payload + wr->dst_offset,
            wr->length,
            MPI_BYTE,
            MPI_ANY_SOURCE,
            wr->tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTR_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr->request_ptr  =&wr->request[RTR_REQ_INDEX];
    wr->request_count=1;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rtr_msg,
            sizeof(wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTS_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr->request_ptr  =&wr->request[RTS_REQ_INDEX];
    wr->request_count=1;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rts_msg,
            sizeof(wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            wr->request_ptr);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}

static NNTI_result_t repost_RTR_RTS_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        mpi_work_request *wr)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    wr->request_ptr  =&wr->request[RTR_REQ_INDEX];
    wr->request_count=2;

    wr->op_state = BUFFER_INIT;

    if (mpi_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    } else {
        log_warn(nnti_debug_level, "bad buffer type: %llu", (uint64_t)mpi_mem_hdl->type);

    }

    MPI_Irecv(
            &wr->rtr_msg,
            sizeof(wr->rtr_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rtr_tag,
            MPI_COMM_WORLD,
            &wr->request[RTR_REQ_INDEX]);

    MPI_Irecv(
            &wr->rts_msg,
            sizeof(wr->rts_msg),
            MPI_BYTE,
            MPI_ANY_SOURCE,
            mpi_mem_hdl->rts_tag,
            MPI_COMM_WORLD,
            &wr->request[RTS_REQ_INDEX]);

    mpi_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr->request_ptr=%p)", reg_buf, wr->request_ptr);

    return(NNTI_OK);
}


static int is_wr_complete(
        mpi_work_request *wr)
{
    int8_t rc=FALSE;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    mpi_mem_hdl=(mpi_memory_handle *)wr->reg_buf->transport_private;
    assert(mpi_mem_hdl);

    switch (mpi_mem_hdl->type) {
        case SEND_BUFFER:
            if (wr->op_state == SEND_COMPLETE) {
                rc=TRUE;
            }
            break;
        case PUT_SRC_BUFFER:
            if (wr->op_state == RDMA_WRITE_COMPLETE) {
                rc=TRUE;
            }
            break;
        case GET_DST_BUFFER:
            if (wr->op_state == RDMA_READ_COMPLETE) {
                rc=TRUE;
            }
            break;
        case REQUEST_BUFFER:
        case RECEIVE_BUFFER:
            if (wr->op_state == RECV_COMPLETE) {
                rc=TRUE;
            }
            break;
        case PUT_DST_BUFFER:
            if (wr->op_state == RDMA_WRITE_COMPLETE) {
                rc=TRUE;
            }
            break;
        case GET_SRC_BUFFER:
            if (wr->op_state == RDMA_READ_COMPLETE) {
                rc=TRUE;
            }
            break;
        case RDMA_TARGET_BUFFER:
            if (wr->op_state == RDMA_TARGET_COMPLETE) {
                rc=TRUE;
            }
            break;
        case UNKNOWN_BUFFER:
            break;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static mpi_work_request *first_incomplete_wr(
        mpi_memory_handle *mpi_mem_hdl)
{
    mpi_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(mpi_mem_hdl);

    if (mpi_mem_hdl->wr_queue.empty()) {
        log_debug(nnti_debug_level, "work request queue is empty");
    } else {
        wr_queue_iter_t i;
        for (i=mpi_mem_hdl->wr_queue.begin(); i != mpi_mem_hdl->wr_queue.end(); i++) {
            wr=*i;
            assert(wr);
            if (is_wr_complete(wr) == FALSE) {
                break;
            }
        }
    }

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);
    return(wr);
}

static int8_t is_wr_queue_empty(
        const NNTI_buffer_t *reg_buf)
{
    int8_t rc=FALSE;
    mpi_memory_handle *mpi_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    if (mpi_mem_hdl->wr_queue.empty()) {
        rc=TRUE;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}


static int is_buf_op_complete(
        const NNTI_buffer_t *reg_buf)
{
    int8_t rc=FALSE;
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr=NULL;
//    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
    assert(mpi_mem_hdl);

    if (is_wr_queue_empty(reg_buf) == TRUE) {
        log_debug(nnti_debug_level, "work request queue is empty - return FALSE");
        rc=FALSE;
    } else {
        wr=mpi_mem_hdl->wr_queue.front();
        assert(wr);

        rc = is_wr_complete(wr);
    }

    if (rc==TRUE) {
        log_debug(nnti_debug_level, "op is complete");
    }
    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(rc);
}

static int8_t is_any_buf_op_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        uint32_t             *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<buf_count;i++) {
        if ((buf_list[i] != NULL) &&
            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
            (is_buf_op_complete(buf_list[i]) == TRUE)) {

            *which=i;
            rc = TRUE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_all_buf_ops_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count)
{
    int8_t rc=TRUE;

    log_debug(nnti_debug_level, "enter");

    for (uint32_t i=0;i<buf_count;i++) {
        if ((buf_list[i] != NULL) &&
            (is_wr_queue_empty(buf_list[i]) == FALSE) &&
            (is_buf_op_complete(buf_list[i]) == FALSE)) {

            rc = FALSE;
            break;
        }
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)buf->payload);

    nthread_lock(&nnti_buf_bufhash_lock);
    assert(buffers_by_bufhash.find(h) == buffers_by_bufhash.end());
    buffers_by_bufhash[h] = buf;
    nthread_unlock(&nnti_buf_bufhash_lock);

    log_debug(nnti_debug_level, "bufhash buffer added (buf=%p bufhash=%lx)", buf, h);

    return(rc);
}
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash)
{
    NNTI_buffer_t *buf=NULL;

    log_debug(nnti_debug_level, "looking for bufhash=%x", (uint64_t)bufhash);
    nthread_lock(&nnti_buf_bufhash_lock);
    if (buffers_by_bufhash.find(bufhash) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[bufhash];
    }
    nthread_unlock(&nnti_buf_bufhash_lock);

    if (buf != NULL) {
        log_debug(nnti_debug_level, "buffer found (buf=%p)", buf);
        return buf;
    }

    log_debug(nnti_debug_level, "buffer NOT found");
//    print_bufhash_map();

    return(NULL);
}
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *buf)
{
    uint32_t h=hash6432shift((uint64_t)buf->payload);
    log_level debug_level = nnti_debug_level;

    nthread_lock(&nnti_buf_bufhash_lock);
    if (buffers_by_bufhash.find(h) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[h];
    }
    nthread_unlock(&nnti_buf_bufhash_lock);

    if (buf != NULL) {
        log_debug(debug_level, "buffer found");
        buffers_by_bufhash.erase(h);
    } else {
        log_debug(debug_level, "buffer NOT found");
    }

    return(buf);
}
static void print_bufhash_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    if (buffers_by_bufhash.empty()) {
        log_debug(nnti_debug_level, "bufhash_map is empty");
        return;
    }

    buf_by_bufhash_iter_t i;
    for (i=buffers_by_bufhash.begin(); i != buffers_by_bufhash.end(); i++) {
        log_debug(nnti_debug_level, "bufhash_map key=%x buf=%p", i->first, i->second);
    }
}

static NNTI_result_t insert_wr_wrhash(mpi_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)wr);

    nthread_lock(&nnti_wr_wrhash_lock);
    assert(wr_by_wrhash.find(h) == wr_by_wrhash.end());
    wr_by_wrhash[h] = wr;
    nthread_unlock(&nnti_wr_wrhash_lock);

    log_debug(nnti_debug_level, "wrhash work request added (wr=%p hash=%x)", wr, h);

    return(rc);
}
static mpi_work_request *get_wr_wrhash(const uint32_t wrhash)
{
    mpi_work_request *wr=NULL;

    log_debug(nnti_debug_level, "looking for wrhash=%x", (uint64_t)wrhash);
    nthread_lock(&nnti_wr_wrhash_lock);
    if (wr_by_wrhash.find(wrhash) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[wrhash];
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

    if (wr != NULL) {
        log_debug(nnti_debug_level, "work request found (wr=%p)", wr);
        return wr;
    }

    log_debug(nnti_debug_level, "work request NOT found");
//    print_wrhash_map();

    return(NULL);
}
static mpi_work_request *del_wr_wrhash(mpi_work_request *wr)
{
    uint32_t h=hash6432shift((uint64_t)wr);
    log_level debug_level = nnti_debug_level;

    nthread_lock(&nnti_wr_wrhash_lock);
    if (wr_by_wrhash.find(h) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[h];
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

    if (wr != NULL) {
        log_debug(debug_level, "work request found");
        wr_by_wrhash.erase(h);
    } else {
        log_debug(debug_level, "work request NOT found");
    }

    return(wr);
}
static void print_wrhash_map()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    if (wr_by_wrhash.empty()) {
        log_debug(nnti_debug_level, "wrhash_map is empty");
        return;
    }

    wr_by_wrhash_iter_t i;
    for (i=wr_by_wrhash.begin(); i != wr_by_wrhash.end(); i++) {
        log_debug(nnti_debug_level, "wrhash_map key=%lx wr=%p", i->first, i->second);
    }
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
    target_buffer_queue_iter_t i;

    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter (buf=%p)", buf);

    nthread_lock(&nnti_target_buffer_queue_lock);
    for (i=target_buffers.begin(); i != target_buffers.end(); i++) {
        if (*i == buf) {
            found = *i;
            break;
        }
    }
    if (found != NULL) {
        log_debug(debug_level, "buffer found");
        target_buffers.erase(i);
    } else {
        log_debug(debug_level, "buffer NOT found");
    }
    nthread_unlock(&nnti_target_buffer_queue_lock);

    log_debug(nnti_debug_level, "exit");

    return(buf);
}
static void print_target_buffer_deque()
{
    if (!logging_debug(nnti_debug_level)) {
        return;
    }

    if (target_buffers.empty()) {
        log_debug(nnti_debug_level, "target_buffers is empty");
        return;
    }

    target_buffer_queue_iter_t i;
    for (i=target_buffers.begin(); i != target_buffers.end(); i++) {
        log_debug(nnti_debug_level, "target buffer (%p)", *i);
    }
}

static void create_status(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        int                   nnti_rc,
        NNTI_status_t        *status)
{
    mpi_memory_handle *mpi_mem_hdl=NULL;
    mpi_work_request  *wr        =NULL;

    log_debug(nnti_debug_level, "enter");

    memset(status, 0, sizeof(NNTI_status_t));
    status->op     = remote_op;
    status->result = (NNTI_result_t)nnti_rc;
    if (nnti_rc==NNTI_OK) {
        mpi_mem_hdl=(mpi_memory_handle *)reg_buf->transport_private;
        assert(mpi_mem_hdl);
        wr=mpi_mem_hdl->wr_queue.front();
        assert(wr);

        status->start  = reg_buf->payload;
        status->length = wr->length;
        status->result = (NNTI_result_t)nnti_rc;
        switch (wr->last_op) {
            case MPI_OP_PUT_INITIATOR:
            case MPI_OP_GET_TARGET:
            case MPI_OP_SEND_REQUEST:
            case MPI_OP_SEND_BUFFER:
                status->offset = wr->src_offset;
                create_peer(&status->src, transport_global_data.rank); // allocates url
                create_peer(&status->dest, wr->last_event.MPI_SOURCE); // allocates url
                break;
            case MPI_OP_GET_INITIATOR:
            case MPI_OP_PUT_TARGET:
            case MPI_OP_NEW_REQUEST:
            case MPI_OP_RECEIVE:
                status->offset = wr->dst_offset;
                create_peer(&status->src, wr->last_event.MPI_SOURCE); // allocates url
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

static void copy_peer(NNTI_peer_t *src, NNTI_peer_t *dest)
{
    log_debug(nnti_debug_level, "enter");

    strcpy(dest->url, src->url);

    dest->peer.transport_id                     = NNTI_TRANSPORT_MPI;
    dest->peer.NNTI_remote_process_t_u.mpi.rank = src->peer.NNTI_remote_process_t_u.mpi.rank;

    log_debug(nnti_debug_level, "exit");
}
