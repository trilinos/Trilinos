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
 * nnti_ib.c
 *
 *  Created on: Jan 13, 2011
 *      Author: thkorde
 */

#include "Trios_config.h"

#include "Trios_threads.h"
#include "Trios_timer.h"
#include "Trios_signal.h"
#include "Trios_nnti_fprint_types.h"

#include <assert.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/poll.h>
#include <sys/mman.h>
#include <netinet/in.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <verbs.h>

#include <map>
#include <deque>
#include <algorithm>

#ifdef HAVE_TRIOS_HPCTOOLKIT
#include <hpctoolkit.h>
#define SAMPLING_IS_ACTIVE() hpctoolkit_sampling_is_active()
#define SAMPLING_STOP() hpctoolkit_sampling_stop()
#define SAMPLING_START() hpctoolkit_sampling_start()
#else
#define SAMPLING_IS_ACTIVE() 0
#define SAMPLING_STOP()
#define SAMPLING_START()
#endif

#include "nnti_ib.h"
#include "nnti_utils.h"



//#define USE_WR_POOL
//#undef USE_WR_POOL
//
//
///* if undefined, the ACK message is NOT sent to the RDMA target when
// * the RDMA op is complete.  this creates one-sided semantics for RDMA
// * ops.  in this mode, the target has no idea when the RDMA op is
// * complete and what data was addressed.  NNTI_wait() returns NNTI_EINVAL
// * if passed a target buffer.
// */
//#undef USE_RDMA_TARGET_ACK
///* if defined, the RDMA initiator will send an ACK message to the RDMA
// * target when the RDMA op is complete.  the target process must wait
// * on the target buffer in order to get the ACK.  this creates two-sided
// * semantics for RDMA ops.   in this mode, when the wait returns the
// * the RDMA op is complete and status indicates what data was addressed.
// */
//#define USE_RDMA_TARGET_ACK



typedef struct {

    uint32_t min_atomics_vars;

    bool     use_wr_pool;

    bool     use_rdma_target_ack;

    bool     use_mlock;
    bool     use_memset;

    bool     drop_if_full_queue;

} nnti_ib_config;


/**
 * These states are used to signal events between the completion handler
 * and the main client or server thread.
 *
 * Once CONNECTED, they cycle through RDMA_READ_ADV, RDMA_WRITE_ADV,
 * and RDMA_WRITE_COMPLETE for each ping.
 */
typedef enum {
    IDLE = 1,
    CONNECT_REQUEST,
    ADDR_RESOLVED,
    ROUTE_RESOLVED,
    CONNECTED,
    DISCONNECTED,
    ERROR
} ib_connection_state;

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
} ib_buffer_type;

#define IB_OP_PUT_INITIATOR  1
#define IB_OP_GET_INITIATOR  2
#define IB_OP_PUT_TARGET     3
#define IB_OP_GET_TARGET     4
#define IB_OP_SEND_REQUEST   5
#define IB_OP_SEND_BUFFER    6
#define IB_OP_NEW_REQUEST    7
#define IB_OP_RECEIVE        8
#define IB_OP_FETCH_ADD      9
#define IB_OP_COMPARE_SWAP  10

typedef enum {
    NNTI_IB_WR_STATE_RESET,         // this work request is idle
    NNTI_IB_WR_STATE_POSTED,        // target is ready to receive data
    NNTI_IB_WR_STATE_STARTED,       // initiator has posted the sq_wr descriptor
    NNTI_IB_WR_STATE_ATTACHED,      // target has attached this work request to an NNTI_work_request_t
    NNTI_IB_WR_STATE_CANCELING,     // NNTI_cancel() called, but not NNTI_wait()
    NNTI_IB_WR_STATE_RDMA_COMPLETE, // RDMA op complete
    NNTI_IB_WR_STATE_ACK_COMPLETE,   // work completion transfer complete
    NNTI_IB_WR_STATE_WAIT_COMPLETE, // NNTI_wait() called and completed successfully
    NNTI_IB_WR_STATE_ERROR          // something went wrong.  check wr->result
} nnti_ib_wr_state_t;

#define SRQ_DEPTH 20480
#define CQ_DEPTH 20480

typedef struct {
    struct ibv_qp           *qp;
    uint32_t                 qpn;
    uint32_t                 peer_qpn;
} conn_qp;
typedef struct {
    NNTI_peer_t   peer;
    char         *peer_name;
    NNTI_ip_addr  peer_addr;
    NNTI_tcp_port peer_port;
    uint16_t      peer_lid;
    uint32_t      peer_req_qpn;

    conn_qp       req_qp;
    conn_qp       data_qp;

    uint32_t      atomics_rkey;
    uint64_t      atomics_addr;

    ib_connection_state state;

    int8_t disconnect_requested;
} ib_connection;

typedef struct {
    uint32_t op;
    uint64_t offset;
    uint64_t length;
} ib_rdma_ack;

typedef struct {
    nthread_lock_t lock;

    NNTI_work_request_t *nnti_wr;

    nnti_ib_wr_state_t state;

    uint32_t key;

    NNTI_buffer_t *reg_buf;
    ib_connection *conn;

    struct ibv_send_wr sq_wr;
    struct ibv_recv_wr rq_wr;
    struct ibv_sge     sge;
    struct ibv_sge    *sge_list;
    uint32_t           sge_count;

    struct ibv_comp_channel *comp_channel;
    struct ibv_cq           *cq;
    struct ibv_qp           *qp;
    uint32_t                 qpn;
    uint32_t                 peer_qpn;

    struct ibv_send_wr ack_sq_wr;
    struct ibv_recv_wr ack_rq_wr;
    struct ibv_sge     ack_sge;
    struct ibv_mr     *ack_mr;
    ib_rdma_ack        ack;

    /* this is a copy of the last work completion that arrived for this buffer */
    struct ibv_wc    last_wc;

    uint8_t       last_op;
    uint64_t      offset;
    uint64_t      length;

} ib_work_request;

typedef std::deque<ib_work_request *>           wr_queue_t;
typedef std::deque<ib_work_request *>::iterator wr_queue_iter_t;

typedef struct {
    ib_buffer_type  type;
    struct ibv_mr  *mr;
    struct ibv_mr **mr_list;
    uint32_t        mr_count;
    wr_queue_t      wr_queue;
    nthread_lock_t  wr_queue_lock;
    uint32_t        ref_count;
} ib_memory_handle;

typedef struct {
    NNTI_buffer_t *reg_buf;

    char *req_buffer;  /* incoming queue */
    int   req_count;   /* incoming queue */
    int   req_size;    /* each message is no larger than req_size */

    uint32_t req_received;
} ib_request_queue_handle;

typedef struct {
    struct ibv_device       *dev;
    struct ibv_context      *ctx;
    struct ibv_pd           *pd;

    uint32_t                 cqe_count;
    uint32_t                 srq_count;
    uint32_t                 sge_count;
    uint32_t                 qp_count;

    struct ibv_comp_channel *req_comp_channel;
    struct ibv_cq           *req_cq;
    struct ibv_srq          *req_srq;

    struct ibv_comp_channel *data_comp_channel;
    struct ibv_cq           *data_cq;
    struct ibv_srq          *data_srq;

    uint32_t req_srq_count;
    uint32_t data_srq_count;

    uint16_t nic_lid;
    int      nic_port;

    int      listen_sock;
    char     listen_name[NNTI_HOSTNAME_LEN];
    uint32_t listen_addr;  /* in NBO */
    uint16_t listen_port;  /* in NBO */

    int interrupt_pipe[2];

    int64_t        *atomics;
    struct ibv_mr  *atomics_mr;
    nthread_lock_t  atomics_lock;

    ib_request_queue_handle req_queue;
} ib_transport_global;




static nthread_lock_t nnti_ib_lock;
static nthread_lock_t nnti_progress_lock;
static nthread_cond_t nnti_progress_cond;

static void *aligned_malloc(
        size_t size);
static struct ibv_mr *register_memory_segment(
        void *buf,
        uint64_t len,
        enum ibv_access_flags access);
static int unregister_memory_segment(
        struct ibv_mr *mr);
static int register_ack(
        ib_work_request *ib_wr);
static int unregister_ack(
        ib_work_request *ib_wr);
static void send_ack (
        ib_work_request *ib_wr);
static NNTI_result_t setup_data_channel(void);
static NNTI_result_t setup_request_channel(void);
static NNTI_result_t setup_interrupt_pipe(void);
static NNTI_result_t setup_atomics(void);
static ib_work_request *decode_work_request(
        const struct ibv_wc *wc);
static int cancel_wr(
        ib_work_request *ib_wr);
static int process_event(
        ib_work_request     *ib_wr,
        const struct ibv_wc *wc);
static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         wr_id,
        uint64_t        offset);
static NNTI_result_t post_ack_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t repost_recv_work_request(
        NNTI_work_request_t  *wr,
        ib_work_request      *ib_wr);
static NNTI_result_t repost_ack_recv_work_request(
        NNTI_work_request_t *wr,
        ib_work_request     *ib_wr);
static int8_t is_wr_canceling(
        ib_work_request *ib_wr);
static int8_t is_wr_complete(
        ib_work_request *ib_wr);
static int8_t is_any_wr_complete(
        ib_work_request **wr_list,
        const uint32_t    wr_count,
        uint32_t         *which);
static int8_t is_all_wr_complete(
        ib_work_request **wr_list,
        const uint32_t    wr_count);
static int8_t is_wr_complete(
        NNTI_work_request_t *wr);
static int8_t is_any_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        uint32_t             *which);
static int8_t is_all_wr_complete(
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);
static void create_status(
        NNTI_work_request_t *wr,
        ib_work_request     *ib_wr,
        NNTI_result_t        nnti_rc,
        NNTI_status_t       *status);
static void create_peer(
        NNTI_peer_t *peer,
        char *name,
        NNTI_ip_addr addr,
        NNTI_tcp_port port);
//static void copy_peer(
//        NNTI_peer_t *src,
//        NNTI_peer_t *dest);
static int init_server_listen_socket(void);
static NNTI_result_t check_listen_socket_for_new_connections();
static struct ibv_device *get_ib_device(void);
static void get_qp_state(struct ibv_qp *qp);
static void print_all_qp_state(void);
static int tcp_read(int sock, void *incoming, size_t len);
static int tcp_write(int sock, const void *outgoing, size_t len);
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len);
static void transition_connection_to_ready(
        int sock,
        ib_connection *conn);
static void transition_qp_from_reset_to_ready(
        struct ibv_qp *qp,
        uint32_t       peer_qpn,
        int            peer_lid,
        int            min_rnr_timer,
        int            ack_timeout,
        int            retry_count);
static void transition_qp_from_error_to_ready(
        struct ibv_qp *qp,
        uint32_t       peer_qpn,
        int            peer_lid,
        int            min_rnr_timer,
        int            ack_timeout,
        int            retry_count);
static void transition_connection_to_error(
        ib_connection *conn);
static void transition_qp_to_error(
        struct ibv_qp *qp,
        uint32_t peer_qpn,
        int peer_lid);
static NNTI_result_t init_connection(
        ib_connection **conn,
        const int sock,
        const int is_server);
static void close_connection(ib_connection *c);
static void print_wc(
        const struct ibv_wc *wc,
        bool force);
static void print_ib_conn(ib_connection *c);
//static void print_qpn_map(void);
//static void print_peer_map(void);
static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, ib_connection *conn);
static NNTI_result_t insert_conn_qpn(const NNTI_qp_num qpn, ib_connection *conn);
static ib_connection *get_conn_peer(const NNTI_peer_t *peer);
static ib_connection *get_conn_qpn(const NNTI_qp_num qpn);
static NNTI_peer_t *get_peer_by_url(const char *url);
static ib_connection *del_conn_peer(const NNTI_peer_t *peer);
static ib_connection *del_conn_qpn(const NNTI_qp_num qpn);

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf);
static NNTI_buffer_t *get_buf_bufhash(const uint32_t bufhash);
static NNTI_buffer_t *del_buf_bufhash(NNTI_buffer_t *buf);
//static void print_bufhash_map(void);

static NNTI_result_t wr_pool_init(uint32_t pool_size);
static ib_work_request *wr_pool_rdma_pop(void);
static ib_work_request *wr_pool_sendrecv_pop(void);
static void wr_pool_rdma_push(ib_work_request *ib_wr);
static void wr_pool_sendrecv_push(ib_work_request *ib_wr);
static NNTI_result_t wr_pool_fini(void);

static void close_all_conn(void);

static NNTI_result_t poll_all(
        int timeout);
static NNTI_result_t progress(
        int                   timeout,
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count);

static void config_init(
        nnti_ib_config *c);
static void config_get_from_env(
        nnti_ib_config *c);

//static void print_wr(ib_work_request *ib_wr);
//static void print_xfer_buf(void *buf, uint32_t size);
//static void print_ack_buf(ib_rdma_ack *ack);

#define IB_MEM_HDL(b) ((ib_memory_handle *)((b)->transport_private))
#define IB_WORK_REQUEST(ib_wr) ((ib_work_request *)((ib_wr)->transport_private))

static bool ib_initialized=false;

static ib_transport_global transport_global_data;
static const int MIN_TIMEOUT = 0;  /* in milliseconds */


/**
 * This custom key is used to look up existing connections.
 */
struct addrport_key {
    NNTI_ip_addr    addr;       /* part1 of a compound key */
    NNTI_tcp_port   port;       /* part2 of a compound key */

    // Need this operators for the hash map
    bool operator<(const addrport_key &key1) const {
        if (addr < key1.addr) {
            return true;
        } else if (addr == key1.addr) {
            if (port < key1.port) {
                return true;
            }
        }
        return false;
    }
    bool operator>(const addrport_key &key1) const {
        if (addr > key1.addr) {
            return true;
        } else if (addr == key1.addr) {
            if (port > key1.port) {
                return true;
            }
        }
        return false;
    }
};


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

/*
 * We need a couple of maps to keep track of connections.  Servers need to find
 * connections by QP number when requests arrive.  Clients need to find connections
 * by peer address and port.  Setup those maps here.
 */
static std::map<addrport_key, ib_connection *> connections_by_peer;
typedef std::map<addrport_key, ib_connection *>::iterator conn_by_peer_iter_t;
typedef std::pair<addrport_key, ib_connection *> conn_by_peer_t;
static nthread_lock_t nnti_conn_peer_lock;

static std::map<NNTI_qp_num, ib_connection *> connections_by_qpn;
typedef std::map<NNTI_qp_num, ib_connection *>::iterator conn_by_qpn_iter_t;
typedef std::pair<NNTI_qp_num, ib_connection *> conn_by_qpn_t;
static nthread_lock_t nnti_conn_qpn_lock;

static std::map<uint32_t, NNTI_buffer_t *> buffers_by_bufhash;
typedef std::map<uint32_t, NNTI_buffer_t *>::iterator buf_by_bufhash_iter_t;
typedef std::pair<uint32_t, NNTI_buffer_t *> buf_by_bufhash_t;
static nthread_lock_t nnti_buf_bufhash_lock;

typedef std::deque<ib_work_request *>           wr_pool_t;
typedef std::deque<ib_work_request *>::iterator wr_pool_iter_t;
static nthread_lock_t nnti_wr_pool_lock;

typedef uint32_t wr_key_t;
static std::map<wr_key_t, ib_work_request *> wrmap;
typedef std::map<wr_key_t, ib_work_request *>::iterator wrmap_iter_t;
//typedef std::pair<wr_key_t, ib_work_request *> wrmap_pair_t;
static nthread_lock_t nnti_wrmap_lock;
static nthread_counter_t nnti_wrmap_counter;

static wr_pool_t rdma_wr_pool;
static wr_pool_t sendrecv_wr_pool;


static nnti_ib_config config;


/* ---------------- Wrappers to protect HPCToolkit issues ---------- */
static
void ibv_ack_cq_events_wrapper(struct ibv_cq *cq, unsigned int nevents)
{
    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    ibv_ack_cq_events(cq, nevents);
    if (sampling) SAMPLING_START();

    return;
}
static
struct ibv_pd *ibv_alloc_pd_wrapper(struct ibv_context *context)
{
    struct ibv_pd *pd=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    pd=ibv_alloc_pd(context);
    if (sampling) SAMPLING_START();

    return pd;
}
static
int ibv_close_device_wrapper(struct ibv_context *context)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_close_device(context);
    if (sampling) SAMPLING_START();

    return rc;
}
static
struct ibv_comp_channel *ibv_create_comp_channel_wrapper(struct ibv_context *context)
{
    struct ibv_comp_channel *channel=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    channel=ibv_create_comp_channel(context);
    if (sampling) SAMPLING_START();

    return channel;
}
static
struct ibv_cq *ibv_create_cq_wrapper(struct ibv_context *context, int cqe,
                             void *cq_context,
                             struct ibv_comp_channel *channel,
                             int comp_vector)
{
    struct ibv_cq *cq=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    cq=ibv_create_cq(context, cqe, cq_context, channel, comp_vector);
    if (sampling) SAMPLING_START();

    return cq;
}
static
struct ibv_qp *ibv_create_qp_wrapper(struct ibv_pd *pd,
                             struct ibv_qp_init_attr *qp_init_attr)
{
    struct ibv_qp *qp=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    qp=ibv_create_qp(pd, qp_init_attr);
    if (sampling) SAMPLING_START();

    return qp;
}
static
struct ibv_srq *ibv_create_srq_wrapper(struct ibv_pd *pd,
                               struct ibv_srq_init_attr *srq_init_attr)
{
    struct ibv_srq *srq=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    srq=ibv_create_srq(pd, srq_init_attr);
    if (sampling) SAMPLING_START();

    return srq;
}
static
int ibv_dealloc_pd_wrapper(struct ibv_pd *pd)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_dealloc_pd(pd);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_dereg_mr_wrapper(struct ibv_mr *mr)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_dereg_mr(mr);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_destroy_comp_channel_wrapper(
        struct ibv_comp_channel *channel)
{
    int rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_destroy_comp_channel(channel);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_destroy_cq_wrapper(
        struct ibv_cq *cq)
{
    int rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_destroy_cq(cq);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_destroy_qp_wrapper(struct ibv_qp *qp)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_destroy_qp(qp);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_destroy_srq_wrapper(
        struct ibv_srq *srq)
{
    int rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_destroy_srq(srq);
    if (sampling) SAMPLING_START();

    return rc;
}
static
void ibv_free_device_list_wrapper(struct ibv_device **list)
{
    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    ibv_free_device_list(list);
    if (sampling) SAMPLING_START();

    return;
}
static
int ibv_get_cq_event_wrapper(struct ibv_comp_channel *channel,
                             struct ibv_cq **cq,
                             void **cq_context)
{
    int rc;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_get_cq_event(channel, cq, cq_context);
    if (sampling) SAMPLING_START();

    return rc;
}
static
struct ibv_device **ibv_get_device_list_wrapper(int *num_devices)
{
    struct ibv_device **dev=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    dev=ibv_get_device_list(num_devices);
    if (sampling) SAMPLING_START();

    return dev;
}
static
int ibv_modify_qp_wrapper(struct ibv_qp *qp,
                  struct ibv_qp_attr *attr,
                  int attr_mask)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_modify_qp(qp, attr, attr_mask);
    if (sampling) SAMPLING_START();

    return rc;
}
static
struct ibv_context *ibv_open_device_wrapper(struct ibv_device *device)
{
    struct ibv_context *context=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    context=ibv_open_device(device);
    if (sampling) SAMPLING_START();

    return context;
}
static
int ibv_poll_cq_wrapper(struct ibv_cq *cq, int num_entries, struct ibv_wc *wc)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_poll_cq(cq, num_entries, wc);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_post_send_wrapper(struct ibv_qp *qp,
                  struct ibv_send_wr *ib_wr,
                  struct ibv_send_wr **bad_wr)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_post_send(qp, ib_wr, bad_wr);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_post_srq_recv_wrapper(struct ibv_srq *srq,
                      struct ibv_recv_wr *recv_wr,
                      struct ibv_recv_wr **bad_recv_wr)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_post_srq_recv(srq,  recv_wr, bad_recv_wr);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_query_device_wrapper(struct ibv_context *context,
                     struct ibv_device_attr *device_attr)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_query_device(context, device_attr);
    if (sampling) SAMPLING_START();

    return rc;
}
static
int ibv_query_port_wrapper(struct ibv_context *context, uint8_t port_num,
                   struct ibv_port_attr *port_attr)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_query_port(context, port_num, port_attr);
    if (sampling) SAMPLING_START();

    return rc;
}
static
struct ibv_mr *ibv_reg_mr_wrapper(struct ibv_pd *pd, void *addr,
                          size_t length, int access)
{
    struct ibv_mr *mr=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    mr=ibv_reg_mr(pd, addr, length, access);
    if (sampling) SAMPLING_START();

    return mr;
}
static
int ibv_req_notify_cq_wrapper(struct ibv_cq *cq, int solicited_only)
{
    int rc=0;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    rc=ibv_req_notify_cq(cq, solicited_only);
    if (sampling) SAMPLING_START();

    return rc;
}
static
const char *ibv_wc_status_str_wrapper(enum ibv_wc_status status)
{
    const char *str=NULL;

    bool sampling = SAMPLING_IS_ACTIVE();
    if (sampling) SAMPLING_STOP();
    str=ibv_wc_status_str(status);
    if (sampling) SAMPLING_START();

    return str;
}



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
NNTI_result_t NNTI_ib_init (
        const NNTI_transport_id_t  trans_id,
        const char                *my_url,
        NNTI_transport_t          *trans_hdl)
{
    int rc=0;

    struct ibv_device_attr dev_attr;
    struct ibv_port_attr   dev_port_attr;
//    uint32_t cqe_count;
//    uint32_t srq_count;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char *sep;
//    char *endptr;

    char hostname[NNTI_HOSTNAME_LEN];
//    NNTI_ip_addr  addr;
//    NNTI_tcp_port port;

    assert(trans_hdl);


    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "my_url=%s", my_url);
    log_debug(nnti_debug_level, "initialized=%d, FALSE==%d", (int)ib_initialized, (int)FALSE);

    if (!ib_initialized) {

        nthread_lock_init(&nnti_ib_lock);

        nthread_lock_init(&nnti_progress_lock);
        nthread_cond_init(&nnti_progress_cond);

        nthread_lock_init(&nnti_conn_peer_lock);
        nthread_lock_init(&nnti_conn_qpn_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);

        nthread_lock_init(&nnti_wrmap_lock);
        nthread_counter_init(&nnti_wrmap_counter);

        nthread_lock_init(&transport_global_data.atomics_lock);

        nthread_lock_init(&nnti_wr_pool_lock);

        config_init(&config);
        config_get_from_env(&config);

        log_debug(nnti_debug_level, "my_url=%s", my_url);

        if (my_url != NULL) {
            if ((rc=nnti_url_get_transport(my_url, transport, NNTI_URL_LEN)) != NNTI_OK) {
                return(NNTI_EINVAL);
            }
            if (0!=strcmp(transport, "ib")) {
                return(NNTI_EINVAL);
            }

            if ((rc=nnti_url_get_address(my_url, address, NNTI_URL_LEN)) != NNTI_OK) {
                return(NNTI_EINVAL);
            }

            sep=strchr(address, ':');
            if (sep == address) {
                /* no hostname given; try gethostname */
                gethostname(hostname, NNTI_HOSTNAME_LEN);
            } else {
                strncpy(hostname, address, sep-address);
                hostname[sep-address]='\0';
            }
//            sep++;
//            port=strtol(sep, &endptr, 0);
//            if (endptr == sep) {
//                /* no port given; use -1 */
//                port=-1;
//            }
        } else {
            gethostname(hostname, NNTI_HOSTNAME_LEN);
//            port=-1;
        }
        strcpy(transport_global_data.listen_name, hostname);


        log_debug(nnti_debug_level, "initializing InfiniBand");

//        /* register trace groups (let someone else enable) */
//        trace_counter_gid = trace_get_gid(TRACE_RPC_COUNTER_GNAME);
//        trace_interval_gid = trace_get_gid(TRACE_RPC_INTERVAL_GNAME);


        memset(&transport_global_data, 0, sizeof(ib_transport_global));

        struct ibv_device *dev=get_ib_device();

        /* open the device */
        transport_global_data.ctx = ibv_open_device_wrapper(dev);
        if (!transport_global_data.ctx) {
            log_error(nnti_debug_level, "ibv_open_device failed");
            return NNTI_EIO;
        }

        transport_global_data.nic_port = 1;

        /* get the lid and verify port state */
        rc = ibv_query_port_wrapper(transport_global_data.ctx, transport_global_data.nic_port, &dev_port_attr);
        if (rc) {
            log_error(nnti_debug_level, "ibv_query_port failed");
            return NNTI_EIO;
        }

        transport_global_data.nic_lid = dev_port_attr.lid;

        if (dev_port_attr.state != IBV_PORT_ACTIVE) {
            log_error(nnti_debug_level, "port is not active.  cannot continue.");
            return NNTI_EIO;
        }

        /* Query the device for the max_ requests and such */
        rc = ibv_query_device_wrapper(transport_global_data.ctx, &dev_attr);
        if (rc) {
            log_error(nnti_debug_level, "ibv_query_device failed");
            return NNTI_EIO;
        }

        log_debug(nnti_debug_level, "max %d completion queue entries", dev_attr.max_cqe);
        transport_global_data.cqe_count = dev_attr.max_cqe;

        log_debug(nnti_debug_level, "max %d shared receive queue work requests", dev_attr.max_srq_wr);
        transport_global_data.srq_count = dev_attr.max_srq_wr/2;

        log_debug(nnti_debug_level, "max %d shared receive queue scatter gather elements", dev_attr.max_srq_sge);
        transport_global_data.sge_count = dev_attr.max_srq_sge;

        log_debug(nnti_debug_level, "max %d queue pair work requests", dev_attr.max_qp_wr);
        transport_global_data.qp_count = dev_attr.max_qp_wr/2;

        /* Allocate a Protection Domain (global) */
        transport_global_data.pd = ibv_alloc_pd_wrapper(transport_global_data.ctx);
        if (!transport_global_data.pd) {
            log_error(nnti_debug_level, "ibv_alloc_pd failed");
            return NNTI_EIO;
        }

        setup_request_channel();
        setup_data_channel();

        if (config.use_wr_pool) {
            rc=wr_pool_init(101);
            if (rc!=NNTI_OK) {
                log_error(nnti_debug_level, "wr_pool_init(): %d", rc);
                goto cleanup;
            }
        }

        setup_interrupt_pipe();

        setup_atomics();

        init_server_listen_socket();

        if (logging_info(nnti_debug_level)) {
            fprintf(logger_get_file(), "InfiniBand Initialized: host(%s) port(%u)\n",
                    transport_global_data.listen_name,
                    ntohs(transport_global_data.listen_port));
        }

        create_peer(
                &trans_hdl->me,
                transport_global_data.listen_name,
                transport_global_data.listen_addr,
                transport_global_data.listen_port);

        ib_initialized = true;
    }

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
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
NNTI_result_t NNTI_ib_get_url (
        const NNTI_transport_t *trans_hdl,
        char                   *url,
        const uint64_t          maxlen)
{
    NNTI_result_t rc=NNTI_OK;

    assert(trans_hdl);
    assert(url);
    assert(maxlen>0);

    strncpy(url, trans_hdl->me.url, maxlen);
    url[maxlen-1]='\0';

    return(rc);
}


/**
 * @brief Prepare for communication with the peer identified by <tt>url</tt>.
 *
 * Parse <tt>url</tt> in a transport specific way.  Perform any transport specific
 * actions necessary to begin communication with this peer.
 *
 *
 * Connectionless transport: parse and populate
 * Connected transport: parse, connection and populate
 *
 */
NNTI_result_t NNTI_ib_connect (
        const NNTI_transport_t *trans_hdl,
        const char             *url,
        const int               timeout,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(callTime);

    NNTI_peer_t *existing_peer=NULL;

    char transport[NNTI_URL_LEN];
    char address[NNTI_URL_LEN];
    char *sep;

    char          hostname[NNTI_HOSTNAME_LEN];
    NNTI_tcp_port port;
    int s;
    struct hostent *host_entry;
    struct sockaddr_in skin;

    ib_connection *conn=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(peer_hdl);

    existing_peer=get_peer_by_url(url);
    if (existing_peer!=NULL) {
        *peer_hdl=*existing_peer;
        return(NNTI_OK);
    }

    if (url != NULL) {
        if ((rc=nnti_url_get_transport(url, transport, NNTI_URL_LEN)) != NNTI_OK) {
            return(rc);
        }
        if (0!=strcmp(transport, "ib")) {
            /* the peer described by 'url' is not an IB peer */
            return(NNTI_EINVAL);
        }

        if ((rc=nnti_url_get_address(url, address, NNTI_URL_LEN)) != NNTI_OK) {
            return(rc);
        }

        sep=strchr(address, ':');
        strncpy(hostname, address, sep-address);
        hostname[sep-address]='\0';
        port=strtol(sep+1, NULL, 0);
    } else {
        /*  */
        return(NNTI_EINVAL);
    }



    s = socket(AF_INET, SOCK_STREAM, 0);
    if (s < 0) {
        log_warn(nnti_debug_level, "failed to create tcp socket: errno=%d (%s)", errno, strerror(errno));
        return NNTI_EIO;
    }

    host_entry = gethostbyname(hostname);
    if (!host_entry) {
        log_warn(nnti_debug_level, "failed to resolve server name (%s): %s", hostname, strerror(errno));
        return NNTI_ENOENT;
    }
    memset(&skin, 0, sizeof(skin));
    skin.sin_family = host_entry->h_addrtype;
    memcpy(&skin.sin_addr, host_entry->h_addr_list[0], (size_t) host_entry->h_length);
    skin.sin_port = htons(port);

    trios_start_timer(callTime);
retry:
    if (connect(s, (struct sockaddr *) &skin, sizeof(skin)) < 0) {
        if (errno == EINTR) {
            goto retry;
        } else {
            log_warn(nnti_debug_level, "failed to connect to server (%s:%u): errno=%d (%s)", hostname, port, errno, strerror(errno));
            return NNTI_EIO;
        }
    }
    trios_stop_timer("socket connect", callTime);

    conn = (ib_connection *)calloc(1, sizeof(ib_connection));
    log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
    if (conn == NULL) {
        log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
        rc=NNTI_ENOMEM;
        goto cleanup;
    }

    trios_start_timer(callTime);
    init_connection(&conn, s, 0);
    trios_stop_timer("ib init connection", callTime);

    create_peer(
            peer_hdl,
            conn->peer_name,
            conn->peer_addr,
            conn->peer_port);

    conn->peer=*peer_hdl;

    insert_conn_qpn(conn->req_qp.qpn, conn);
    insert_conn_qpn(conn->data_qp.qpn, conn);
    insert_conn_peer(peer_hdl, conn);

    transition_connection_to_ready(s, conn);

    if (close(s) < 0) {
        log_warn(nnti_debug_level, "failed to close tcp socket: errno=%d (%s)", errno, strerror(errno));
        return NNTI_EIO;
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer_hdl",
                "end of NNTI_ib_connect", peer_hdl);
    }

cleanup:
    log_debug(nnti_debug_level, "enter");

    return(rc);
}


/**
 * @brief Terminate communication with this peer.
 *
 * Perform any transport specific actions necessary to end communication with
 * this peer.
 */
NNTI_result_t NNTI_ib_disconnect (
        const NNTI_transport_t *trans_hdl,
        NNTI_peer_t            *peer_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "Disconnecting from IB");

    assert(trans_hdl);
    assert(peer_hdl);

    ib_connection *conn=get_conn_peer(peer_hdl);
    close_connection(conn);
    del_conn_peer(peer_hdl);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_ib_alloc (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    char *buf=(char *)aligned_malloc(element_size*num_elements);
    assert(buf);

    nnti_rc=NNTI_ib_register_memory(
            trans_hdl,
            buf,
            element_size,
            num_elements,
            ops,
            reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ib_alloc", reg_buf);
    }

    if (config.use_memset) {
        trios_start_timer(callTime);
        memset(buf, 0, element_size*num_elements);
        trios_stop_timer("memset", callTime);
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
NNTI_result_t NNTI_ib_free (
        NNTI_buffer_t *reg_buf)
{
    NNTI_result_t nnti_rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    char *buf=NNTI_BUFFER_C_POINTER(reg_buf);
    assert(buf);

    nnti_rc=NNTI_ib_unregister_memory(reg_buf);

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
NNTI_result_t NNTI_ib_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    uint32_t cqe_num;

    NNTI_buffer_t    *old_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

    old_buf=get_buf_bufhash(hash6432shift((uint64_t)buffer));
    if (old_buf==NULL) {
        ib_mem_hdl=new ib_memory_handle();
        assert(ib_mem_hdl);
        ib_mem_hdl->ref_count=1;
        nthread_lock_init(&ib_mem_hdl->wr_queue_lock);
    } else {
        ib_mem_hdl=IB_MEM_HDL(old_buf);
        ib_mem_hdl->ref_count++;
    }

    log_debug(nnti_debug_level, "ib_mem_hdl->ref_count==%lu", ib_mem_hdl->ref_count);

    memset(reg_buf, 0, sizeof(NNTI_buffer_t));

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size      = element_size;
    reg_buf->payload           = (uint64_t)buffer;
    reg_buf->transport_private = (uint64_t)ib_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (ib_mem_hdl->ref_count==1) {

        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(1, sizeof(NNTI_remote_addr_t));
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=1;

        ib_mem_hdl->mr_list=&ib_mem_hdl->mr;
        ib_mem_hdl->mr_count=1;

        if (ops == NNTI_RECV_QUEUE) {
            ib_request_queue_handle *q_hdl=&transport_global_data.req_queue;

            ib_mem_hdl->type=REQUEST_BUFFER;

            q_hdl->reg_buf=reg_buf;

            q_hdl->req_buffer  =buffer;
            q_hdl->req_size    =element_size;
            q_hdl->req_count   =num_elements;
            q_hdl->req_received=0;

            reg_buf->payload_size=q_hdl->req_size;

            cqe_num=q_hdl->req_count;
            if (cqe_num >= transport_global_data.srq_count) {
                cqe_num = transport_global_data.srq_count;
            }

            ib_mem_hdl->mr=register_memory_segment(
                                buffer,
                                num_elements*element_size,
                                IBV_ACCESS_LOCAL_WRITE);

            log_debug(nnti_debug_level, "mr=%p mr_list[0]=%p", ib_mem_hdl->mr, ib_mem_hdl->mr_list[0]);

            for (uint32_t i=0;i<cqe_num;i++) {
                post_recv_work_request(
                        reg_buf,
                        -1,
                        (i*q_hdl->req_size));
            }

        } else {
            ib_mem_hdl->mr=register_memory_segment(
                                buffer,
                                element_size,
                                (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

            log_debug(nnti_debug_level, "mr=%p mr_list[0]=%p", ib_mem_hdl->mr, ib_mem_hdl->mr_list[0]);

            if (ops == NNTI_RECV_DST) {
                ib_mem_hdl->type=RECEIVE_BUFFER;

                post_recv_work_request(
                        reg_buf,
                        -1,
                        0);

            } else if (ops == NNTI_SEND_SRC) {
                ib_mem_hdl->type=SEND_BUFFER;

            } else if (ops == NNTI_GET_DST) {
                ib_mem_hdl->type=GET_DST_BUFFER;

            } else if (ops == NNTI_GET_SRC) {
                ib_mem_hdl->type=GET_SRC_BUFFER;

            } else if (ops == NNTI_PUT_SRC) {
                ib_mem_hdl->type=PUT_SRC_BUFFER;

            } else if (ops == NNTI_PUT_DST) {
                ib_mem_hdl->type=PUT_DST_BUFFER;

            } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
                ib_mem_hdl->type=RDMA_TARGET_BUFFER;

            } else {
                ib_mem_hdl->type=UNKNOWN_BUFFER;
            }
        }
    } else {
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(1, sizeof(NNTI_remote_addr_t));
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=1;
    }

    if (config.use_rdma_target_ack) {
        if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
            (ib_mem_hdl->type == GET_SRC_BUFFER) ||
            (ib_mem_hdl->type == PUT_DST_BUFFER)) {
            post_ack_recv_work_request(reg_buf);
        }
    }

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].transport_id                 = NNTI_TRANSPORT_IB;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.size = ib_mem_hdl->mr_list[0]->length;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf  = (uint64_t)ib_mem_hdl->mr_list[0]->addr;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.key  = ib_mem_hdl->mr_list[0]->rkey;

    if (ib_mem_hdl->ref_count==1) {
        insert_buf_bufhash(reg_buf);
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ib_register_memory", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Prepare a block of memory for network operations.
 *
 * Wrap a user allocated block of memory in an NNTI_buffer_t.  The transport
 * may take additional actions to prepare the memory for network send/receive.
 * If the memory block doesn't meet the transport's requirements for memory
 * regions, then errors or poor performance may result.
 */
NNTI_result_t NNTI_ib_register_segments (
        const NNTI_transport_t *trans_hdl,
        char                  **segments,
        const uint64_t         *segment_lengths,
        const uint64_t          num_segments,
        const NNTI_buf_ops_t    ops,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    uint32_t cqe_num;

    NNTI_buffer_t     *old_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(trans_hdl);
    assert(segments);
    assert(segment_lengths>0);
    assert(num_segments>0);
    assert(ops>0);
    assert(reg_buf);

    if ((ops == NNTI_SEND_SRC) || (ops == NNTI_RECV_DST) || (ops == NNTI_RECV_QUEUE)) {
        log_debug(nnti_debug_level, "NNTI_SEND_SRC, NNTI_RECV_DST and NNTI_RECV_QUEUE types cannot be segmented.");
        return(NNTI_EINVAL);
    }

    old_buf=get_buf_bufhash(hash6432shift((uint64_t)segments[0]));
    if (old_buf==NULL) {
        ib_mem_hdl=new ib_memory_handle();
        assert(ib_mem_hdl);
        ib_mem_hdl->ref_count=1;
        nthread_lock_init(&ib_mem_hdl->wr_queue_lock);
    } else {
        // check that the number of old_buf segments equals num_segments.
        if (old_buf->buffer_segments.NNTI_remote_addr_array_t_len != num_segments) {
            log_debug(nnti_debug_level, "Segment count mismatch (old=%d new=%d).  Aborting...",
                    old_buf->buffer_segments.NNTI_remote_addr_array_t_len, num_segments);
            return(NNTI_EINVAL);
        }
        // check that all the old_buf segments match the current list of segments
        for (int i=0;i<old_buf->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
            if (old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.buf != (uint64_t)segments[i]) {
                log_debug(nnti_debug_level, "Segment address mismatch (old[%d]=%p new[%d]=%p).  Aborting...",
                        i, old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.buf, i, segments[i]);
                return(NNTI_EINVAL);
            }
            if (old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size != segment_lengths[i]) {
                log_debug(nnti_debug_level, "Segment length mismatch (old[%d]=%d new[%d]=%d).  Aborting...",
                        i, old_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size, i, segment_lengths[i]);
                return(NNTI_EINVAL);
            }
        }
        ib_mem_hdl=IB_MEM_HDL(old_buf);
        ib_mem_hdl->ref_count++;
    }

    log_debug(nnti_debug_level, "ib_mem_hdl->ref_count==%lu", ib_mem_hdl->ref_count);

    memset(reg_buf, 0, sizeof(NNTI_buffer_t));

    reg_buf->transport_id      = trans_hdl->id;
    reg_buf->buffer_owner      = trans_hdl->me;
    reg_buf->ops               = ops;
    reg_buf->payload_size=0;
    for (int i=0;i<num_segments;i++) {
        reg_buf->payload_size += segment_lengths[i];
    }
    reg_buf->payload           = (uint64_t)segments[0];
    reg_buf->transport_private = (uint64_t)ib_mem_hdl;

    log_debug(nnti_debug_level, "rpc_buffer->payload_size=%ld",
            reg_buf->payload_size);

    if (ib_mem_hdl->ref_count==1) {

        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val=(NNTI_remote_addr_t *)calloc(num_segments, sizeof(NNTI_remote_addr_t));
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_len=num_segments;

        ib_mem_hdl->mr_list=(struct ibv_mr **)calloc(num_segments, sizeof(struct ibv_mr *));
        ib_mem_hdl->mr_count=num_segments;

        for (int i=0;i<num_segments;i++) {
            ib_mem_hdl->mr_list[i]=register_memory_segment(
                                        segments[i],
                                        segment_lengths[i],
                                        (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

            log_debug(nnti_debug_level, "mr_list[%d]=%p", i, ib_mem_hdl->mr_list[i]);
        }

        if (ops == NNTI_GET_DST) {
            ib_mem_hdl->type=GET_DST_BUFFER;

        } else if (ops == NNTI_GET_SRC) {
            ib_mem_hdl->type=GET_SRC_BUFFER;

        } else if (ops == NNTI_PUT_SRC) {
            ib_mem_hdl->type=PUT_SRC_BUFFER;

        } else if (ops == NNTI_PUT_DST) {
            ib_mem_hdl->type=PUT_DST_BUFFER;

        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            ib_mem_hdl->type=RDMA_TARGET_BUFFER;

        } else {
            ib_mem_hdl->type=UNKNOWN_BUFFER;
        }
    }

    if (config.use_rdma_target_ack) {
        if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
            (ib_mem_hdl->type == GET_SRC_BUFFER) ||
            (ib_mem_hdl->type == PUT_DST_BUFFER)) {
            post_ack_recv_work_request(reg_buf);
        }
    }

    for (int i=0;i<num_segments;i++) {
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].transport_id                 = NNTI_TRANSPORT_IB;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size = ib_mem_hdl->mr_list[i]->length;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.buf  = (uint64_t)ib_mem_hdl->mr_list[i]->addr;
        reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.key  = ib_mem_hdl->mr_list[i]->rkey;
    }

    if (ib_mem_hdl->ref_count==1) {
        insert_buf_bufhash(reg_buf);
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ib_register_memory", reg_buf);
    }

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_ib_unregister_memory (
        NNTI_buffer_t *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;
    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "start of NNTI_ib_unregister_memory", reg_buf);
    }

    ib_mem_hdl=IB_MEM_HDL(reg_buf);
    assert(ib_mem_hdl);
    ib_mem_hdl->ref_count--;

    log_debug(nnti_debug_level, "ib_mem_hdl->ref_count==%lu", ib_mem_hdl->ref_count);

    if (ib_mem_hdl->ref_count==0) {
        log_debug(nnti_debug_level, "ib_mem_hdl->ref_count is 0.  release all resources.");
        log_debug(nnti_debug_level, "This buffer has %d segments.", ib_mem_hdl->mr_count);
        for (int i=0;i<ib_mem_hdl->mr_count;i++) {
            log_debug(nnti_debug_level, "Unregistering segment #%d.", i);
            unregister_memory_segment(ib_mem_hdl->mr_list[i]);
        }

        del_buf_bufhash(reg_buf);

        nthread_lock(&ib_mem_hdl->wr_queue_lock);
        while (!ib_mem_hdl->wr_queue.empty()) {
            ib_work_request *ib_wr=ib_mem_hdl->wr_queue.front();
            log_debug(nnti_debug_level, "popping pending (reg_buf=%p, ib_wr=%p)", reg_buf, ib_wr);
            if (ib_wr==NULL)
                break;
            ib_mem_hdl->wr_queue.pop_front();
            if (config.use_wr_pool) {
                if (ib_wr->ack_mr!=NULL) {
                    wr_pool_rdma_push(ib_wr);
                } else {
                    wr_pool_sendrecv_push(ib_wr);
                }
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(ib_wr);
                }
                log_debug(nnti_debug_level, "freeing ib_wr=%p", ib_wr);
                free(ib_wr);
            }
        }
        nthread_unlock(&ib_mem_hdl->wr_queue_lock);

        nthread_lock_fini(&ib_mem_hdl->wr_queue_lock);

        if (ib_mem_hdl) delete ib_mem_hdl;

        if (reg_buf->buffer_segments.NNTI_remote_addr_array_t_val != NULL) {
            free(reg_buf->buffer_segments.NNTI_remote_addr_array_t_val);
        }
        if ((ib_mem_hdl->mr_count > 1) && (ib_mem_hdl->mr_list != NULL)) {
            free(ib_mem_hdl->mr_list);
        }

        reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
        IB_SET_MATCH_ANY(&reg_buf->buffer_owner);
        reg_buf->ops               = (NNTI_buf_ops_t)0;
        reg_buf->payload_size      = 0;
        reg_buf->payload           = 0;
        reg_buf->transport_private = 0;
    }

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Send a message to a peer.
 *
 * Send a message (<tt>msg_hdl</tt>) to a peer (<tt>peer_hdl</tt>).  It is expected that the
 * message is small, but the exact maximum size is transport dependent.
 */
NNTI_result_t NNTI_ib_send (
        const NNTI_peer_t   *peer_hdl,
        const NNTI_buffer_t *msg_hdl,
        const NNTI_buffer_t *dest_hdl,
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    log_level debug_level=nnti_debug_level;

    ib_mem_hdl=IB_MEM_HDL(msg_hdl);
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        ib_wr=wr_pool_sendrecv_pop();
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        log_debug(nnti_debug_level, "allocated ib_wr (wr=%p ; ib_wr=%p)", wr, ib_wr);
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);

    ib_wr->conn = get_conn_peer(peer_hdl);
    assert(ib_wr->conn);

    ib_wr->nnti_wr = wr;
    ib_wr->reg_buf = (NNTI_buffer_t *)msg_hdl;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_STARTED;

    if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE)) {
        ib_wr->comp_channel=transport_global_data.req_comp_channel;
        ib_wr->cq          =transport_global_data.req_cq;
        ib_wr->qp          =ib_wr->conn->req_qp.qp;
        ib_wr->qpn         =(uint64_t)ib_wr->conn->req_qp.qpn;
        ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->peer_req_qpn;

        ib_wr->sq_wr.opcode    =IBV_WR_SEND;
        ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

        ib_wr->last_op=IB_OP_SEND_REQUEST;
    } else {
        ib_wr->comp_channel=transport_global_data.data_comp_channel;
        ib_wr->cq          =transport_global_data.data_cq;
        ib_wr->qp          =ib_wr->conn->data_qp.qp;
        ib_wr->qpn         =(uint64_t)ib_wr->conn->data_qp.qpn;
        ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->data_qp.peer_qpn;

        ib_wr->sq_wr.wr.rdma.rkey       =dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.key;
        ib_wr->sq_wr.wr.rdma.remote_addr=dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf;

        ib_wr->sq_wr.imm_data  =hash6432shift((uint64_t)dest_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
        ib_wr->sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

        ib_wr->last_op=IB_OP_SEND_BUFFER;
    }

    // setup the scatter-gather list for this work request
    ib_wr->sge_count=msg_hdl->buffer_segments.NNTI_remote_addr_array_t_len;
    if (ib_wr->sge_count == 1) {
        ib_wr->sge_list=&ib_wr->sge;
        ib_wr->sge_list[0].addr  =(uint64_t)ib_mem_hdl->mr_list[0]->addr;
        ib_wr->sge_list[0].length=msg_hdl->payload_size;
        ib_wr->sge_list[0].lkey  =ib_mem_hdl->mr_list[0]->lkey;
    } else {
        ib_wr->sge_list=(struct ibv_sge *)calloc(ib_wr->sge_count, sizeof(struct ibv_sge));
        for (int i=0;i<ib_wr->sge_count;i++) {
            ib_wr->sge_list[i].addr  =(uint64_t)ib_mem_hdl->mr_list[i]->addr;
            ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[i]->length;
            ib_wr->sge_list[i].lkey  =ib_mem_hdl->mr_list[i]->lkey;
        }
    }
    ib_wr->sq_wr.wr_id   = (uint64_t)ib_wr->key;
    ib_wr->sq_wr.next    = NULL;  // RAOLDFI ADDED
    ib_wr->sq_wr.sg_list = ib_wr->sge_list;
    ib_wr->sq_wr.num_sge = ib_wr->sge_count;


    log_debug(nnti_debug_level, "sending to (%s, qp=%p, qpn=%d, sge.addr=%p, sge.length=%llu, sq_wr.imm_data=0x%llx, sq_wr.ib_wr.rdma.rkey=0x%x, sq_wr.ib_wr.rdma.remote_addr=%p)",
            peer_hdl->url,
            ib_wr->qp,
            ib_wr->qpn,
            (void *)  ib_wr->sge.addr,
            (uint64_t)ib_wr->sge.length,
            (uint64_t)ib_wr->sq_wr.imm_data,
                      ib_wr->sq_wr.wr.rdma.rkey,
            (void *)  ib_wr->sq_wr.wr.rdma.remote_addr);

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    wr->transport_id     =msg_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)msg_hdl;
    wr->ops              =NNTI_SEND_SRC;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)ib_wr;

    log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
    trios_start_timer(call_time);
    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_send - ibv_post_send", call_time);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Transfer data to a peer.
 *
 * Put the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_ib_put (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ib_mem_hdl=IB_MEM_HDL(src_buffer_hdl);
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        if (config.use_rdma_target_ack) {
            ib_wr=wr_pool_rdma_pop();
        } else {
            ib_wr=wr_pool_sendrecv_pop();
        }
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        log_debug(nnti_debug_level, "allocated ib_wr (wr=%p ; ib_wr=%p)", wr, ib_wr);
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);

    ib_wr->conn = get_conn_peer(&dest_buffer_hdl->buffer_owner);
    assert(ib_wr->conn);

    ib_wr->nnti_wr = wr;
    ib_wr->reg_buf = (NNTI_buffer_t *)src_buffer_hdl;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_STARTED;

    ib_wr->comp_channel=transport_global_data.data_comp_channel;
    ib_wr->cq          =transport_global_data.data_cq;
    ib_wr->qp          =ib_wr->conn->data_qp.qp;
    ib_wr->qpn         =(uint64_t)ib_wr->conn->data_qp.qpn;
    ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->data_qp.peer_qpn;

    if (dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len == 1) {
        // this is the easy case.  the destination (remote) buffer is contiguous so we can complete this PUT with one ibv_send_wr.

        // the src_offset could exclude some local segments from this operation.  find the first segment.
        uint64_t segment_offset=src_offset;
        uint32_t first_segment =0;
        for (int i=0;i<src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
            if (segment_offset > src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size) {
                segment_offset -= src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size;
                continue;
            } else if (segment_offset == src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size) {
                segment_offset = 0;
                first_segment  = i+1;
                break;
            } else {
                first_segment  = i;
                break;
            }
        }

        log_debug(nnti_debug_level, "first_segment=%u, segment_offset=%lu", first_segment, segment_offset);

        // setup the scatter-gather list for this work request
        ib_wr->sge_count=src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len - first_segment;
        if (ib_wr->sge_count == 1) {
            ib_wr->sge_list=&ib_wr->sge;
            ib_wr->sge_list[0].addr  =(uint64_t)ib_mem_hdl->mr_list[first_segment]->addr + segment_offset;
            ib_wr->sge_list[0].length=src_length;
            ib_wr->sge_list[0].lkey  =ib_mem_hdl->mr_list[first_segment]->lkey;
        } else {
            ib_wr->sge_list=(struct ibv_sge *)calloc(ib_wr->sge_count, sizeof(struct ibv_sge));
            for (int i=0,j=first_segment;i<ib_wr->sge_count;i++,j++) {
                ib_wr->sge_list[i].addr=(uint64_t)ib_mem_hdl->mr_list[j]->addr + segment_offset;
                if (i == 0) {
                    ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[j]->length - segment_offset;
                } else {
                    ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[j]->length;
                }
                ib_wr->sge_list[i].lkey=ib_mem_hdl->mr_list[j]->lkey;
            }
        }
        ib_wr->sq_wr.sg_list=ib_wr->sge_list;
        ib_wr->sq_wr.num_sge=ib_wr->sge_count;

        ib_wr->sq_wr.wr.rdma.rkey        = dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.key;
        ib_wr->sq_wr.wr.rdma.remote_addr = dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf+dest_offset;

        ib_wr->sq_wr.opcode    =IBV_WR_RDMA_WRITE;
        ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;
        ib_wr->sq_wr.wr_id     =(uint64_t)ib_wr->key;
        ib_wr->sq_wr.imm_data  =hash6432shift((uint64_t)dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
        ib_wr->sq_wr.next      =NULL;  // RAOLDFI ADDED

    } else {
        // this is the hard case.  the destination (remote) buffer is non-contiguous so we need multiple ibv_send_wr to complete this PUT.

    }

    if (config.use_rdma_target_ack) {
        ib_wr->ack.op    =IB_OP_PUT_TARGET;
        ib_wr->ack.offset=dest_offset;
        ib_wr->ack.length=src_length;

        if (!config.use_wr_pool) {
            register_ack(ib_wr);
        }
        ib_wr->ack_sge.addr  =(uint64_t)ib_wr->ack_mr->addr;
        ib_wr->ack_sge.length=ib_wr->ack_mr->length;
        ib_wr->ack_sge.lkey  =ib_wr->ack_mr->lkey;

        ib_wr->ack_sq_wr.sg_list=&ib_wr->ack_sge;
        ib_wr->ack_sq_wr.num_sge=1;

        ib_wr->ack_sq_wr.wr.rdma.rkey       =dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_key;
        ib_wr->ack_sq_wr.wr.rdma.remote_addr=dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_buf;

        ib_wr->ack_sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        ib_wr->ack_sq_wr.send_flags=IBV_SEND_SIGNALED|IBV_SEND_FENCE;
        ib_wr->sq_wr.wr_id         =(uint64_t)ib_wr->key;
        ib_wr->ack_sq_wr.imm_data  =hash6432shift((uint64_t)dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
        ib_wr->ack_sq_wr.next      =NULL;  // RAOLDFI ADDED
    }

    ib_wr->last_op=IB_OP_PUT_INITIATOR;
    ib_wr->length=src_length;
    ib_wr->offset=src_offset;

    log_debug(nnti_debug_level, "putting to (%s, qp=%p, qpn=%lu)",
            dest_buffer_hdl->buffer_owner.url,
            ib_wr->qp,
            ib_wr->qpn);

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    wr->transport_id     =src_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)src_buffer_hdl;
    wr->ops              =NNTI_PUT_SRC;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)ib_wr;

    log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
    trios_start_timer(call_time);
    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_put - ibv_post_send", call_time);

    if (config.use_rdma_target_ack) {
        send_ack(ib_wr);
    }

    if (ib_wr->sge_list != &ib_wr->sge) {
        free(ib_wr->sge_list);
    }

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Transfer data from a peer.
 *
 * Get the contents of <tt>src_buffer_hdl</tt> into <tt>dest_buffer_hdl</tt>.  It is
 * assumed that the destination is at least <tt>src_length</tt> bytes in size.
 *
 */
NNTI_result_t NNTI_ib_get (
        const NNTI_buffer_t *src_buffer_hdl,
        const uint64_t       src_offset,
        const uint64_t       src_length,
        const NNTI_buffer_t *dest_buffer_hdl,
        const uint64_t       dest_offset,
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;
    log_level debug_level = nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    struct ibv_send_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    trios_start_timer(total_time);

    log_debug(debug_level, "enter (wr=%p)", wr);

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ib_mem_hdl=IB_MEM_HDL(dest_buffer_hdl);
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        if (config.use_rdma_target_ack) {
            ib_wr=wr_pool_rdma_pop();
        } else {
            ib_wr=wr_pool_sendrecv_pop();
        }
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        log_debug(nnti_debug_level, "allocated ib_wr (wr=%p ; ib_wr=%p)", wr, ib_wr);
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);

    ib_wr->conn = get_conn_peer(&src_buffer_hdl->buffer_owner);
    assert(ib_wr->conn);

    ib_wr->nnti_wr = wr;
    ib_wr->reg_buf = (NNTI_buffer_t *)dest_buffer_hdl;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_STARTED;

    ib_wr->comp_channel=transport_global_data.data_comp_channel;
    ib_wr->cq          =transport_global_data.data_cq;
    ib_wr->qp          =ib_wr->conn->data_qp.qp;
    ib_wr->qpn         =(uint64_t)ib_wr->conn->data_qp.qpn;
    ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->data_qp.peer_qpn;

    if (src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len == 1) {
        // this is the easy case.  the source (remote) buffer is contiguous so we can complete this GET with one ibv_send_wr.

        // the src_offset could exclude some local segments from this operation.  find the first segment.
        uint64_t segment_offset=dest_offset;
        uint32_t first_segment =0;
        for (int i=0;i<dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len;i++) {
            if (segment_offset > dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size) {
                segment_offset -= dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size;
                continue;
            } else if (segment_offset == dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[i].NNTI_remote_addr_t_u.ib.size) {
                segment_offset = 0;
                first_segment  = i+1;
                break;
            } else {
                first_segment  = i;
                break;
            }
        }

        log_debug(nnti_debug_level, "first_segment=%u, segment_offset=%lu", first_segment, segment_offset);

        // setup the scatter-gather list for this work request
        ib_wr->sge_count=dest_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_len - first_segment;
        if (ib_wr->sge_count == 1) {
            ib_wr->sge_list=&ib_wr->sge;
            ib_wr->sge_list[0].addr  =(uint64_t)ib_mem_hdl->mr_list[first_segment]->addr + segment_offset;
            ib_wr->sge_list[0].length=src_length;
            ib_wr->sge_list[0].lkey  =ib_mem_hdl->mr_list[first_segment]->lkey;
        } else {
            ib_wr->sge_list=(struct ibv_sge *)calloc(ib_wr->sge_count, sizeof(struct ibv_sge));
            for (int i=0,j=first_segment;i<ib_wr->sge_count;i++,j++) {
                ib_wr->sge_list[i].addr=(uint64_t)ib_mem_hdl->mr_list[j]->addr + segment_offset;
                if (i == 0) {
                    ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[j]->length - segment_offset;
                } else {
                    ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[j]->length;
                }
                ib_wr->sge_list[i].lkey=ib_mem_hdl->mr_list[j]->lkey;
            }
        }
        ib_wr->sq_wr.sg_list=ib_wr->sge_list;
        ib_wr->sq_wr.num_sge=ib_wr->sge_count;

        for (int i=0;i<ib_wr->sge_count;i++) {
            log_debug(nnti_debug_level, "ib_wr->sge_list[%d].addr=0x%lX ; ib_wr->sge_list[%d].length=%u ; ib_wr->sge_list[%d].lkey=0x%X",
                    i, ib_wr->sge_list[i].addr, i, ib_wr->sge_list[i].length, i, ib_wr->sge_list[i].lkey);
        }
        for (int i=0;i<ib_wr->sq_wr.num_sge;i++) {
            log_debug(nnti_debug_level, "ib_wr->sq_wr.sg_list[%d].addr=0x%lX ; ib_wr->sq_wr.sg_list[%d].length=%u ; ib_wr->sq_wr.sg_list[%d].lkey=0x%X",
                    i, ib_wr->sq_wr.sg_list[i].addr, i, ib_wr->sq_wr.sg_list[i].length, i, ib_wr->sq_wr.sg_list[i].lkey);
        }


        ib_wr->sq_wr.wr.rdma.rkey        = src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.key;
        ib_wr->sq_wr.wr.rdma.remote_addr = src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf+src_offset;

        ib_wr->sq_wr.opcode    =IBV_WR_RDMA_READ;
        ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;
        ib_wr->sq_wr.wr_id     =(uint64_t)ib_wr->key;
        ib_wr->sq_wr.imm_data  =hash6432shift((uint64_t)src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
        ib_wr->sq_wr.next      =NULL;  // RAOLDFI ADDED

    } else {
        // this is the hard case.  the source (remote) buffer is non-contiguous so we need multiple ibv_send_wr to complete this GET.

    }

    if (config.use_rdma_target_ack) {
        ib_wr->ack.op    =IB_OP_GET_TARGET;
        ib_wr->ack.offset=src_offset;
        ib_wr->ack.length=src_length;

        if (!config.use_wr_pool) {
            register_ack(ib_wr);
        }
        ib_wr->ack_sge.addr  =(uint64_t)ib_wr->ack_mr->addr;
        ib_wr->ack_sge.length=ib_wr->ack_mr->length;
        ib_wr->ack_sge.lkey  =ib_wr->ack_mr->lkey;

        ib_wr->ack_sq_wr.sg_list=&ib_wr->ack_sge;
        ib_wr->ack_sq_wr.num_sge=1;

        ib_wr->ack_sq_wr.wr.rdma.rkey       =src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_key;
        ib_wr->ack_sq_wr.wr.rdma.remote_addr=src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_buf;

        ib_wr->ack_sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        ib_wr->ack_sq_wr.send_flags=IBV_SEND_SIGNALED|IBV_SEND_FENCE;
        ib_wr->sq_wr.wr_id         =(uint64_t)ib_wr->key;
        ib_wr->ack_sq_wr.imm_data  =hash6432shift((uint64_t)src_buffer_hdl->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
        ib_wr->ack_sq_wr.next      =NULL;  // RAOLDFI ADDED
    }

    ib_wr->last_op=IB_OP_GET_INITIATOR;
    ib_wr->length=src_length;
    ib_wr->offset=dest_offset;

    log_debug(debug_level, "getting from (%s, qp=%p, qpn=%lu)",
            src_buffer_hdl->buffer_owner.url,
            ib_wr->qp,
            ib_wr->qpn);

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    wr->transport_id     =dest_buffer_hdl->transport_id;
    wr->reg_buf          =(NNTI_buffer_t*)dest_buffer_hdl;
    wr->ops              =NNTI_GET_DST;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)ib_wr;

    log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
    trios_start_timer(call_time);
    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_get - ibv_post_send", call_time);

    if (config.use_rdma_target_ack) {
        send_ack(ib_wr);
    }

//    print_wr(ib_wr);

    log_debug(nnti_debug_level, "exit");

    trios_stop_timer("NNTI_ib_get - total", total_time);

    return(rc);
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
NNTI_result_t NNTI_ib_scatter (
        const NNTI_buffer_t  *src_buffer_hdl,
        const uint64_t        src_length,
        const NNTI_buffer_t **dest_buffer_list,
        const uint64_t        dest_count,
        NNTI_work_request_t  *ib_wr)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

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
NNTI_result_t NNTI_ib_gather (
        const NNTI_buffer_t **src_buffer_list,
        const uint64_t        src_length,
        const uint64_t        src_count,
        const NNTI_buffer_t  *dest_buffer_hdl,
        NNTI_work_request_t  *ib_wr)
{
    log_debug(nnti_debug_level, "enter");

    log_debug(nnti_debug_level, "exit");

    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_ib_atomic_set_callback (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          local_atomic,
        NNTI_callback_fn_t      cbfunc,
        void                   *context)
{
    return NNTI_ENOTSUP;
}


NNTI_result_t NNTI_ib_atomic_read (
        const NNTI_transport_t *trans_hdl,
        const uint64_t          local_atomic,
        int64_t                *value)
{
    nthread_lock(&transport_global_data.atomics_lock);
    *value = transport_global_data.atomics[local_atomic];
    nthread_unlock(&transport_global_data.atomics_lock);

    return NNTI_OK;
}


NNTI_result_t NNTI_ib_atomic_fop (
        const NNTI_transport_t *trans_hdl,
        const NNTI_peer_t      *peer_hdl,
        const uint64_t          target_atomic,
        const uint64_t          result_atomic,
        const int64_t           operand,
        const NNTI_atomic_op_t  op,
        NNTI_work_request_t    *wr)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr;

    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);

    log_level debug_level=nnti_debug_level;

    if (config.use_wr_pool) {
        ib_wr=wr_pool_sendrecv_pop();
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        log_debug(nnti_debug_level, "allocated ib_wr (wr=%p ; ib_wr=%p)", wr, ib_wr);
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);

    ib_wr->conn = get_conn_peer(peer_hdl);
    assert(ib_wr->conn);

    ib_wr->nnti_wr = wr;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_STARTED;

    ib_wr->comp_channel=transport_global_data.data_comp_channel;
    ib_wr->cq          =transport_global_data.data_cq;
    ib_wr->qp          =ib_wr->conn->data_qp.qp;
    ib_wr->qpn         =(uint64_t)ib_wr->conn->data_qp.qpn;
    ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->data_qp.peer_qpn;

    ib_wr->last_op=IB_OP_FETCH_ADD;

    ib_wr->sq_wr.wr.atomic.rkey       =ib_wr->conn->atomics_rkey;
    ib_wr->sq_wr.wr.atomic.remote_addr=ib_wr->conn->atomics_addr+(target_atomic*sizeof(int64_t));
    ib_wr->sq_wr.wr.atomic.compare_add=operand;

    ib_wr->sq_wr.opcode    =IBV_WR_ATOMIC_FETCH_AND_ADD;
    ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

    ib_wr->sge_count=1;
    ib_wr->sge_list=&ib_wr->sge;
    ib_wr->sge_list[0].addr  =(uint64_t)&transport_global_data.atomics[result_atomic];
    ib_wr->sge_list[0].length=sizeof(int64_t);
    ib_wr->sge_list[0].lkey  =transport_global_data.atomics_mr->lkey;

    ib_wr->sq_wr.wr_id   = (uint64_t)ib_wr->key;
    ib_wr->sq_wr.next    = NULL;
    ib_wr->sq_wr.sg_list = ib_wr->sge_list;
    ib_wr->sq_wr.num_sge = ib_wr->sge_count;


    log_debug(nnti_debug_level, "sending to (%s, qp=%p, qpn=%d, sge.addr=%p, sge.length=%llu, sq_wr.ib_wr.rdma.rkey=0x%x, sq_wr.ib_wr.rdma.remote_addr=%p)",
            peer_hdl->url,
            ib_wr->qp,
            ib_wr->qpn,
            (void *)  ib_wr->sge.addr,
            (uint64_t)ib_wr->sge.length,
                      ib_wr->sq_wr.wr.atomic.rkey,
            (void *)  ib_wr->sq_wr.wr.atomic.remote_addr);

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)ib_wr;

    log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
    trios_start_timer(call_time);
    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_send - ibv_post_send", call_time);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


NNTI_result_t NNTI_ib_atomic_cswap (
        const NNTI_transport_t *trans_hdl,
        const NNTI_peer_t      *peer_hdl,
        const uint64_t          target_atomic,
        const uint64_t          result_atomic,
        const int64_t           compare_operand,
        const int64_t           swap_operand,
        NNTI_work_request_t    *wr)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr;

    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);

    log_level debug_level=nnti_debug_level;

    if (config.use_wr_pool) {
        ib_wr=wr_pool_sendrecv_pop();
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        log_debug(nnti_debug_level, "allocated ib_wr (wr=%p ; ib_wr=%p)", wr, ib_wr);
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);

    ib_wr->conn = get_conn_peer(peer_hdl);
    assert(ib_wr->conn);

    ib_wr->nnti_wr = wr;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_STARTED;

    ib_wr->comp_channel=transport_global_data.data_comp_channel;
    ib_wr->cq          =transport_global_data.data_cq;
    ib_wr->qp          =ib_wr->conn->data_qp.qp;
    ib_wr->qpn         =(uint64_t)ib_wr->conn->data_qp.qpn;
    ib_wr->peer_qpn    =(uint64_t)ib_wr->conn->data_qp.peer_qpn;

    ib_wr->last_op=IB_OP_COMPARE_SWAP;

    ib_wr->sq_wr.wr.atomic.rkey       =ib_wr->conn->atomics_rkey;
    ib_wr->sq_wr.wr.atomic.remote_addr=ib_wr->conn->atomics_addr+(target_atomic*sizeof(int64_t));
    ib_wr->sq_wr.wr.atomic.compare_add=compare_operand;
    ib_wr->sq_wr.wr.atomic.swap       =swap_operand;

    ib_wr->sq_wr.opcode    =IBV_WR_ATOMIC_CMP_AND_SWP;
    ib_wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

    ib_wr->sge_count=1;
    ib_wr->sge_list=&ib_wr->sge;
    ib_wr->sge_list[0].addr  =(uint64_t)&transport_global_data.atomics[result_atomic];
    ib_wr->sge_list[0].length=sizeof(int64_t);
    ib_wr->sge_list[0].lkey  =transport_global_data.atomics_mr->lkey;

    ib_wr->sq_wr.wr_id   = (uint64_t)ib_wr->key;
    ib_wr->sq_wr.next    = NULL;
    ib_wr->sq_wr.sg_list = ib_wr->sge_list;
    ib_wr->sq_wr.num_sge = ib_wr->sge_count;


    log_debug(nnti_debug_level, "sending to (%s, qp=%p, qpn=%d, sge.addr=%p, sge.length=%llu, sq_wr.ib_wr.rdma.rkey=0x%x, sq_wr.ib_wr.rdma.remote_addr=%p)",
            peer_hdl->url,
            ib_wr->qp,
            ib_wr->qpn,
            (void *)  ib_wr->sge.addr,
            (uint64_t)ib_wr->sge.length,
                      ib_wr->sq_wr.wr.atomic.rkey,
            (void *)  ib_wr->sq_wr.wr.atomic.remote_addr);

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    wr->transport_id     =trans_hdl->id;
    wr->reg_buf          =(NNTI_buffer_t*)NULL;
    wr->ops              =NNTI_ATOMICS;
    wr->result           =NNTI_OK;
    wr->transport_private=(uint64_t)ib_wr;

    log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
    trios_start_timer(call_time);
    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_send - ibv_post_send", call_time);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}


/**
 * @brief Create a receive work request that can be used to wait for buffer
 * operations to complete.
 *
 */
NNTI_result_t NNTI_ib_create_work_request (
        NNTI_buffer_t        *reg_buf,
        NNTI_work_request_t  *wr)
{
    wr_queue_iter_t   iter;
    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p ; wr=%p)", reg_buf, wr);

    ib_mem_hdl=IB_MEM_HDL(reg_buf);
    assert(ib_mem_hdl);

    nthread_lock(&ib_mem_hdl->wr_queue_lock);

    for ( iter=ib_mem_hdl->wr_queue.begin() ; iter != ib_mem_hdl->wr_queue.end() ; ++iter ) {
        if ((*iter)->nnti_wr == NULL) {
            ib_wr=*iter;
            break;
        }
    }
    assert(ib_wr);

    wr->transport_private=(uint64_t)ib_wr;
    ib_wr->nnti_wr=wr;

    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    wr->transport_id     =reg_buf->transport_id;
    wr->reg_buf          =reg_buf;
    wr->ops              =reg_buf->ops;
    wr->result           =NNTI_OK;

    log_debug(nnti_debug_level, "assigned ib_work_request(%p) to NNTI_work_request(%p)", ib_wr, wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p ; wr=%p)", reg_buf, wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from a previous receive
 * and prepares it for reuse.
 *
 */
NNTI_result_t NNTI_ib_clear_work_request (
        NNTI_work_request_t  *wr)
{
    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    wr->result           =NNTI_OK;
    wr->transport_private=NULL;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Disassociates a receive work request from reg_buf.
 *
 */
NNTI_result_t NNTI_ib_destroy_work_request (
        NNTI_work_request_t  *wr)
{
    ib_memory_handle *ib_mem_hdl;
    ib_work_request *ib_wr;

    log_debug(nnti_debug_level, "enter (wr=%p ; ib_wr=%p)", wr, ib_wr);

    ib_mem_hdl=IB_MEM_HDL(wr->reg_buf);
    assert(ib_mem_hdl);
    if ((ib_mem_hdl->type == SEND_BUFFER) || (ib_mem_hdl->type == GET_DST_BUFFER) || (ib_mem_hdl->type == PUT_SRC_BUFFER)) {
        // work requrests from initiator buffers are cleaned up in NNTI_ib_wait*().  nothing to do here.
        return(NNTI_OK);
    }

    ib_wr=IB_WORK_REQUEST(wr);
    assert(ib_wr);

    if (ib_mem_hdl->type == REQUEST_BUFFER) {
        if (ib_wr->state == NNTI_IB_WR_STATE_WAIT_COMPLETE) {
            repost_recv_work_request(wr, ib_wr);

            nthread_lock(&ib_mem_hdl->wr_queue_lock);
            wr_queue_iter_t q_victim=find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), ib_wr);
            if (q_victim != ib_mem_hdl->wr_queue.end()) {
                log_debug(nnti_debug_level, "erasing ib_wr=%p from the wr_queue", ib_wr);
                ib_mem_hdl->wr_queue.erase(q_victim);
            }
            assert(find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), ib_wr) == ib_mem_hdl->wr_queue.end());

            ib_mem_hdl->wr_queue.push_back(ib_wr);

            nthread_unlock(&ib_mem_hdl->wr_queue_lock);
        }
        ib_wr->nnti_wr=NULL;

    } else if (ib_mem_hdl->type == RECEIVE_BUFFER) {
        if (ib_wr->state == NNTI_IB_WR_STATE_WAIT_COMPLETE) {
            repost_recv_work_request(wr, ib_wr);

            nthread_lock(&ib_mem_hdl->wr_queue_lock);
            wr_queue_iter_t q_victim=find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), ib_wr);
            if (q_victim != ib_mem_hdl->wr_queue.end()) {
                log_debug(nnti_debug_level, "erasing ib_wr=%p from the wr_queue", ib_wr);
                ib_mem_hdl->wr_queue.erase(q_victim);
            }
            assert(find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), ib_wr) == ib_mem_hdl->wr_queue.end());

            ib_mem_hdl->wr_queue.push_back(ib_wr);

            nthread_unlock(&ib_mem_hdl->wr_queue_lock);
        }
        ib_wr->nnti_wr=NULL;

    }

    wr->transport_id     =NNTI_TRANSPORT_NULL;
    wr->reg_buf          =NULL;
    wr->ops              =(NNTI_buf_ops_t)0;
    wr->result           =NNTI_OK;
    wr->transport_private=NULL;

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}


/**
 * @brief Attempts to cancel an NNTI opertion.
 *
 */
NNTI_result_t NNTI_ib_cancel (
        NNTI_work_request_t *wr)
{
    NNTI_result_t rc=NNTI_OK;
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ib_wr=IB_WORK_REQUEST(wr);
    assert(ib_wr);

    nthread_lock(&ib_wr->lock);
    if ((ib_wr->state == NNTI_IB_WR_STATE_STARTED) || (ib_wr->state == NNTI_IB_WR_STATE_POSTED)) {
        ib_wr->state = NNTI_IB_WR_STATE_CANCELING;
    } else {
        log_warn(nnti_debug_level, "wr=%p is not in progress (current state is %d).  Cannot cancel.", wr, ib_wr->state);
        rc=NNTI_EINVAL;
    }
    nthread_unlock(&ib_wr->lock);

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return rc;
}


/**
 * @brief Attempts to cancel a list of NNTI opertions.
 *
 */
NNTI_result_t NNTI_ib_cancelall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    NNTI_result_t rc=NNTI_OK;
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    for (int i=0;i<wr_count;i++) {
        ib_wr=IB_WORK_REQUEST(wr_list[i]);
        assert(ib_wr);

        nthread_lock(&ib_wr->lock);
        if ((ib_wr->state == NNTI_IB_WR_STATE_STARTED) || (ib_wr->state == NNTI_IB_WR_STATE_POSTED)) {
            ib_wr->state = NNTI_IB_WR_STATE_CANCELING;
        } else {
            log_warn(nnti_debug_level, "wr_list[%d]=%p is not in progress (current state is %d).  Cannot cancel.", i, wr_list[i], ib_wr->state);
            rc=NNTI_EINVAL;
        }
        nthread_unlock(&ib_wr->lock);
    }

    log_debug(nnti_debug_level, "exit");

    return rc;
}


/**
 * @brief Interrupts NNTI_wait*()
 *
 */
NNTI_result_t NNTI_ib_interrupt (
        const NNTI_transport_t *trans_hdl)
{
    uint32_t dummy=0xAAAAAAAA;

    log_debug(nnti_debug_level, "enter");

    write(transport_global_data.interrupt_pipe[1], &dummy, 4);

    log_debug(nnti_debug_level, "exit");

    return NNTI_OK;
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
NNTI_result_t NNTI_ib_wait (
        NNTI_work_request_t  *wr,
        const int             timeout,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *ib_wr=NULL;

    const NNTI_buffer_t  *wait_buf=NULL;

    int ibv_rc=0;
    NNTI_result_t rc;
    long elapsed_time = 0;
    long timeout_per_call;

    log_level debug_level=nnti_debug_level;

    log_level old_log_level=logger_get_default_level();

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    struct ibv_wc wc;

    log_debug(debug_level, "enter (wr=%p ; ib_wr=%p ; timeout=%d)", wr, IB_WORK_REQUEST(wr), timeout);

    trios_start_timer(total_time);

    assert(wr);
    assert(status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "wr->reg_buf",
                "start of NNTI_ib_wait", wr->reg_buf);
    }

    ib_wr=IB_WORK_REQUEST(wr);
    assert(ib_wr);

    if (is_wr_canceling(ib_wr)) {
        cancel_wr(ib_wr);
    } else if (is_wr_complete(ib_wr) == TRUE) {
        log_debug(debug_level, "wr already complete (wr=%p ; ib_wr=%p)", wr, IB_WORK_REQUEST(wr));
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "wr NOT complete (wr=%p ; ib_wr=%p)", wr, IB_WORK_REQUEST(wr));

        trios_start_timer(call_time);

        while (1) {
            rc=progress(timeout-elapsed_time, &wr, 1);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc == NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(debug_level, "progress() successful...");
                nnti_rc = rc;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc = rc;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_debug(debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = rc;
                break;
            }

            if (is_wr_complete(ib_wr) == TRUE) {
                log_debug(debug_level, "wr completed (wr=%p ; ib_wr=%p)", wr, IB_WORK_REQUEST(wr));
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }
        }

        trios_stop_timer("NNTI_ib_wait - wr complete", call_time);
    }

    trios_start_timer(call_time);
    create_status(wr, ib_wr, nnti_rc, status);
    trios_stop_timer("NNTI_ib_wait - create_status", call_time);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ib_wait", status);
    }

    if (is_wr_complete(ib_wr)) {

        trios_start_timer(call_time);

        ib_wr->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;

        if (ib_wr->nnti_wr->ops == NNTI_ATOMICS) {
            nthread_lock(&nnti_wrmap_lock);
            wrmap_iter_t m_victim=wrmap.find(ib_wr->key);
            if (m_victim != wrmap.end()) {
                log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", ib_wr, ib_wr->key);
                wrmap.erase(m_victim);
            }
            nthread_unlock(&nnti_wrmap_lock);

            if (config.use_wr_pool) {
                wr_pool_sendrecv_push(ib_wr);
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(ib_wr);
                }
                log_debug(nnti_debug_level, "freeing ib_wr=%p", ib_wr);
                free(ib_wr);
            }
            /*
             * This work request (wr) has reached a completed (final) state.  wr is reset here.
             */
            wr->transport_private=(uint64_t)NULL;

        } else {
            ib_mem_hdl=IB_MEM_HDL(wr->reg_buf);
            assert(ib_mem_hdl);

            if ((ib_mem_hdl->type == REQUEST_BUFFER) || (ib_mem_hdl->type == RECEIVE_BUFFER)) {
                // defer cleanup to NNTI_ib_destroy_work_request()
            }
            else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                    (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                    (ib_mem_hdl->type == PUT_DST_BUFFER))     {
                if (config.use_rdma_target_ack) {
                    repost_ack_recv_work_request(wr, ib_wr);
                }
            }
            else {
                nthread_lock(&ib_mem_hdl->wr_queue_lock);
                wr_queue_iter_t q_victim=find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), ib_wr);
                if (q_victim != ib_mem_hdl->wr_queue.end()) {
                    log_debug(nnti_debug_level, "erasing ib_wr=%p from the wr_queue", ib_wr);
                    ib_mem_hdl->wr_queue.erase(q_victim);
                }
                nthread_unlock(&ib_mem_hdl->wr_queue_lock);

                nthread_lock(&nnti_wrmap_lock);
                wrmap_iter_t m_victim=wrmap.find(ib_wr->key);
                if (m_victim != wrmap.end()) {
                    log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", ib_wr, ib_wr->key);
                    wrmap.erase(m_victim);
                }
                nthread_unlock(&nnti_wrmap_lock);

                if (config.use_wr_pool) {
                    wr_pool_sendrecv_push(ib_wr);
                } else {
                    if (config.use_rdma_target_ack) {
                        unregister_ack(ib_wr);
                    }
                    log_debug(nnti_debug_level, "freeing ib_wr=%p", ib_wr);
                    free(ib_wr);
                }
                /*
                 * This work request (wr) has reached a completed (final) state.  wr is reset here.
                 */
                wr->transport_private=(uint64_t)NULL;
            }
        }

        trios_stop_timer("NNTI_ib_wait - wr cleanup", call_time);
    }

    trios_stop_timer("NNTI_ib_wait", total_time);

    log_debug(debug_level, "exit (wr=%p ; ib_wr=%p)", wr, IB_WORK_REQUEST(wr));

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
NNTI_result_t NNTI_ib_waitany (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *ib_wr=NULL;
    const NNTI_buffer_t     *wait_buf=NULL;

    int ibv_rc=0;
    NNTI_result_t rc;
    long elapsed_time = 0;
    long timeout_per_call;

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    struct ibv_wc wc;

    long entry_time=trios_get_time_ms();

    log_debug(debug_level, "enter");

    trios_start_timer(total_time);

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_ib_wait(wr_list[0], timeout, status);
        *which=0;
        goto cleanup;
    }

    if (wr_count > 1) {
        for (uint32_t i=0;i<wr_count;i++) {
            if (wr_list[i] != NULL) {
                ib_wr=IB_WORK_REQUEST(wr_list[i]);
                assert(ib_wr);
                if (is_wr_canceling(ib_wr)) {
                    cancel_wr(ib_wr);
                }
            }
        }
    }

    if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
        log_debug(debug_level, "wr already complete (which=%u, wr_list[%d]=%p, ib_wr=%p)", *which, *which, wr_list[*which], IB_WORK_REQUEST(wr_list[*which]));
        IB_WORK_REQUEST(wr_list[*which])->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "wr NOT complete (wr_list=%p)", wr_list);

        while (1) {
            rc=progress(timeout-elapsed_time, wr_list, wr_count);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc == NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(debug_level, "progress() successful...");
                nnti_rc = NNTI_OK;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc=NNTI_EINTR;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_error(debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = NNTI_EIO;
                break;
            }

            if (is_any_wr_complete(wr_list, wr_count, which) == TRUE) {
                log_debug(debug_level, "wr completed (which=%u, wr_list[%d]=%p, ib_wr=%p)", *which, *which, wr_list[*which], IB_WORK_REQUEST(wr_list[*which]));
                IB_WORK_REQUEST(wr_list[*which])->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }
        }
    }

    create_status(wr_list[*which], IB_WORK_REQUEST(wr_list[*which]), nnti_rc, status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ib_waitany", status);
    }

    if (is_wr_complete(IB_WORK_REQUEST(wr_list[*which]))) {

        IB_WORK_REQUEST(wr_list[*which])->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;

        ib_mem_hdl=IB_MEM_HDL(wr_list[*which]->reg_buf);
        assert(ib_mem_hdl);

        if ((ib_mem_hdl->type == REQUEST_BUFFER) || (ib_mem_hdl->type == RECEIVE_BUFFER)) {
            // defer cleanup to NNTI_ib_destroy_work_request()
        }
        else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                (ib_mem_hdl->type == PUT_DST_BUFFER))     {
            if (config.use_rdma_target_ack) {
                repost_ack_recv_work_request(wr_list[*which], IB_WORK_REQUEST(wr_list[*which]));
            }
        }
        else {
            nthread_lock(&ib_mem_hdl->wr_queue_lock);
            wr_queue_iter_t q_victim=find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), IB_WORK_REQUEST(wr_list[*which]));
            if (q_victim != ib_mem_hdl->wr_queue.end()) {
                log_debug(nnti_debug_level, "erasing ib_wr=%p from the wr_queue", IB_WORK_REQUEST(wr_list[*which]));
                ib_mem_hdl->wr_queue.erase(q_victim);
            }
            nthread_unlock(&ib_mem_hdl->wr_queue_lock);

            nthread_lock(&nnti_wrmap_lock);
            wrmap_iter_t m_victim=wrmap.find(IB_WORK_REQUEST(wr_list[*which])->key);
            if (m_victim != wrmap.end()) {
                log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", IB_WORK_REQUEST(wr_list[*which]), IB_WORK_REQUEST(wr_list[*which])->key);
                wrmap.erase(m_victim);
            }
            nthread_unlock(&nnti_wrmap_lock);

            if (config.use_wr_pool) {
                wr_pool_sendrecv_push(IB_WORK_REQUEST(wr_list[*which]));
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(IB_WORK_REQUEST(wr_list[*which]));
                }
                log_debug(nnti_debug_level, "freeing ib_wr=%p", IB_WORK_REQUEST(wr_list[*which]));
                free(IB_WORK_REQUEST(wr_list[*which]));
            }
            /*
             * This work request (wr) has reached a completed (final) state.  wr is reset here.
             */
            wr_list[*which]->transport_private=(uint64_t)NULL;
        }
    }

cleanup:
    trios_stop_timer("NNTI_ib_waitany", total_time);

    log_debug(debug_level, "exit");

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
NNTI_result_t NNTI_ib_waitall (
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count,
        const int             timeout,
        NNTI_status_t       **status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *ib_wr=NULL;
    const NNTI_buffer_t     *wait_buf=NULL;

    int ibv_rc=0;
    NNTI_result_t rc;
    long elapsed_time = 0;
    long timeout_per_call;

    log_level debug_level=nnti_debug_level;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    struct ibv_wc wc;

    long entry_time=trios_get_time_ms();

    log_debug(debug_level, "enter");

    trios_start_timer(total_time);

    assert(wr_list);
    assert(wr_count > 0);
    assert(status);

    if (wr_count == 1) {
        nnti_rc=NNTI_ib_wait(wr_list[0], timeout, status[0]);
        goto cleanup;
    }

    if (wr_count > 1) {
        for (uint32_t i=0;i<wr_count;i++) {
            if (wr_list[i] != NULL) {
                ib_wr=IB_WORK_REQUEST(wr_list[i]);
                assert(ib_wr);
                if (is_wr_canceling(ib_wr)) {
                    cancel_wr(ib_wr);
                }
            }
        }
    }

    if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
        log_debug(debug_level, "all wr already complete (wr_list=%p)", wr_list);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "all wr NOT complete (wr_list=%p)", wr_list);

        while (1) {
            rc=progress(timeout-elapsed_time, wr_list, wr_count);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc == NNTI_OK) {
//                logger_set_default_level(old_log_level);
                log_debug(debug_level, "progress() successful...");
                nnti_rc = NNTI_OK;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                log_debug(debug_level, "progress() timed out...");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: interrupted */
            else if (rc==NNTI_EINTR) {
                log_debug(debug_level, "progress() interrupted...");
//                logger_set_default_level(LOG_OFF);
                nnti_rc=NNTI_EINTR;
                break;
            }
            /* case 4: failure */
            else {
//                logger_set_default_level(old_log_level);
                log_error(debug_level, "progress() failed: %s", strerror(errno));
                nnti_rc = NNTI_EIO;
                break;
            }

            if (is_all_wr_complete(wr_list, wr_count) == TRUE) {
                log_debug(debug_level, "all wr completed (wr_list=%p)", wr_list);
                nnti_rc = NNTI_OK;
                break;
            }

            /* if the caller asked for a legitimate timeout, we need to exit */
            if ((timeout >= 0) && (elapsed_time >= timeout)) {
                nnti_rc = NNTI_ETIMEDOUT;
                break;
            }

            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }
        }
    }

    for (uint32_t i=0;i<wr_count;i++) {

        if (wr_list[i] == NULL) {
            continue;
        }

        create_status(wr_list[i], IB_WORK_REQUEST(wr_list[i]), nnti_rc, status[i]);

        if (logging_debug(debug_level)) {
            fprint_NNTI_status(logger_get_file(), "status",
                    "end of NNTI_ib_waitall", status[i]);
        }

        if (is_wr_complete(IB_WORK_REQUEST(wr_list[i]))) {

            IB_WORK_REQUEST(wr_list[i])->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;

            ib_mem_hdl=IB_MEM_HDL(wr_list[i]->reg_buf);
            assert(ib_mem_hdl);

            if ((ib_mem_hdl->type == REQUEST_BUFFER) || (ib_mem_hdl->type == RECEIVE_BUFFER)) {
                // defer cleanup to NNTI_ib_destroy_work_request()
            }
            else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                    (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                    (ib_mem_hdl->type == PUT_DST_BUFFER))     {
                if (config.use_rdma_target_ack) {
                    repost_ack_recv_work_request(wr_list[i], IB_WORK_REQUEST(wr_list[i]));
                }
            }
            else {
                nthread_lock(&ib_mem_hdl->wr_queue_lock);
                wr_queue_iter_t q_victim=find(ib_mem_hdl->wr_queue.begin(), ib_mem_hdl->wr_queue.end(), IB_WORK_REQUEST(wr_list[i]));
                if (q_victim != ib_mem_hdl->wr_queue.end()) {
                    log_debug(nnti_debug_level, "erasing ib_wr=%p from the wr_queue", IB_WORK_REQUEST(wr_list[i]));
                    ib_mem_hdl->wr_queue.erase(q_victim);
                }
                nthread_unlock(&ib_mem_hdl->wr_queue_lock);

                nthread_lock(&nnti_wrmap_lock);
                wrmap_iter_t m_victim=wrmap.find(IB_WORK_REQUEST(wr_list[i])->key);
                if (m_victim != wrmap.end()) {
                    log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", IB_WORK_REQUEST(wr_list[i]), IB_WORK_REQUEST(wr_list[i])->key);
                    wrmap.erase(m_victim);
                }
                nthread_unlock(&nnti_wrmap_lock);

                if (config.use_wr_pool) {
                    wr_pool_sendrecv_push(IB_WORK_REQUEST(wr_list[i]));
                } else {
                    if (config.use_rdma_target_ack) {
                        unregister_ack(IB_WORK_REQUEST(wr_list[i]));
                    }
                    log_debug(nnti_debug_level, "freeing ib_wr=%p", IB_WORK_REQUEST(wr_list[i]));
                    free(IB_WORK_REQUEST(wr_list[i]));
                }
                /*
                 * This work request (wr) has reached a completed (final) state.  wr is reset here.
                 */
                wr_list[i]->transport_private=(uint64_t)NULL;
            }
        }
    }

cleanup:
    trios_stop_timer("NNTI_ib_waitall", total_time);

    log_debug(debug_level, "exit");

    return(nnti_rc);
}


/**
 * @brief Disable this transport.
 *
 * Shutdown the transport.  Any outstanding sends, gets and puts will be
 * canceled.  Any new transport requests will fail.
 *
 */
NNTI_result_t NNTI_ib_fini (
        const NNTI_transport_t *trans_hdl)
{
    NNTI_result_t rc=NNTI_OK;;

    log_debug(nnti_debug_level, "enter");

    close_all_conn();

    if (config.use_wr_pool) {
        rc=wr_pool_fini();
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "wr_pool_fini() failed: %d", rc);
            rc=NNTI_EINVAL;
        }
    }

    close(transport_global_data.listen_sock);
    transport_global_data.listen_name[0]='\0';
    transport_global_data.listen_addr=0;
    transport_global_data.listen_port=0;

    ibv_destroy_comp_channel_wrapper(transport_global_data.data_comp_channel);
    ibv_destroy_cq_wrapper(transport_global_data.data_cq);
    ibv_destroy_srq_wrapper(transport_global_data.data_srq);

    ibv_destroy_comp_channel_wrapper(transport_global_data.req_comp_channel);
    ibv_destroy_cq_wrapper(transport_global_data.req_cq);
    ibv_destroy_srq_wrapper(transport_global_data.req_srq);

    ibv_dealloc_pd_wrapper(transport_global_data.pd);

    ibv_close_device_wrapper(transport_global_data.ctx);

    nthread_lock_fini(&nnti_ib_lock);

    nthread_lock_fini(&nnti_progress_lock);
    nthread_cond_fini(&nnti_progress_cond);

    nthread_lock_fini(&nnti_conn_peer_lock);
    nthread_lock_fini(&nnti_conn_qpn_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&transport_global_data.atomics_lock);
    nthread_lock_fini(&nnti_wr_pool_lock);

    ib_initialized=false;

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}





static NNTI_result_t setup_request_channel(void)
{
    int flags;

    transport_global_data.req_comp_channel = ibv_create_comp_channel_wrapper(transport_global_data.ctx);
    if (!transport_global_data.req_comp_channel) {
        log_error(nnti_debug_level, "ibv_create_comp_channel failed");
        return NNTI_EIO;
    }
    transport_global_data.req_cq = ibv_create_cq_wrapper(
            transport_global_data.ctx,
            transport_global_data.cqe_count,
            NULL,
            transport_global_data.req_comp_channel,
            0);
    if (!transport_global_data.req_cq) {
        log_error(nnti_debug_level, "ibv_create_cq failed");
        return NNTI_EIO;
    }

    transport_global_data.req_srq_count=0;

    struct ibv_srq_init_attr srq_attr;
    memset(&srq_attr, 0, sizeof(srq_attr)); // initialize to avoid valgrind uninitialized warning
    srq_attr.attr.max_wr = transport_global_data.srq_count;
    srq_attr.attr.max_sge = transport_global_data.sge_count;

    transport_global_data.req_srq = ibv_create_srq_wrapper(transport_global_data.pd, &srq_attr);
    if (!transport_global_data.req_srq)  {
        log_error(nnti_debug_level, "ibv_create_srq failed");
        return NNTI_EIO;
    }

    if (ibv_req_notify_cq_wrapper(transport_global_data.req_cq, 0)) {
        log_error(nnti_debug_level, "ibv_req_notify_cq failed");
        return NNTI_EIO;
    }

    /* use non-blocking IO on the async fd and completion fd */
    flags = fcntl(transport_global_data.ctx->async_fd, F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get async_fd flags");
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.ctx->async_fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set async_fd to nonblocking");
        return NNTI_EIO;
    }

    flags = fcntl(transport_global_data.req_comp_channel->fd, F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get completion fd flags");
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.req_comp_channel->fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set completion fd to nonblocking");
        return NNTI_EIO;
    }


    return(NNTI_OK);
}

static NNTI_result_t setup_data_channel(void)
{
    int flags;

    transport_global_data.data_comp_channel = ibv_create_comp_channel_wrapper(transport_global_data.ctx);
    if (!transport_global_data.data_comp_channel) {
        log_error(nnti_debug_level, "ibv_create_comp_channel failed");
        return NNTI_EIO;
    }
    transport_global_data.data_cq = ibv_create_cq_wrapper(
            transport_global_data.ctx,
            transport_global_data.cqe_count,
            NULL,
            transport_global_data.data_comp_channel,
            0);
    if (!transport_global_data.data_cq) {
        log_error(nnti_debug_level, "ibv_create_cq failed");
        return NNTI_EIO;
    }

    transport_global_data.data_srq_count=0;

    struct ibv_srq_init_attr srq_attr;
    memset(&srq_attr, 0, sizeof(srq_attr));   // initialize to avoid valgrind warning
    srq_attr.attr.max_wr  = transport_global_data.srq_count;
    srq_attr.attr.max_sge = transport_global_data.sge_count;

    transport_global_data.data_srq = ibv_create_srq_wrapper(transport_global_data.pd, &srq_attr);
    if (!transport_global_data.data_srq)  {
        log_error(nnti_debug_level, "ibv_create_srq failed");
        return NNTI_EIO;
    }

    if (ibv_req_notify_cq_wrapper(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "ibv_req_notify_cq failed");
        return NNTI_EIO;
    }

    /* use non-blocking IO on the async fd and completion fd */
    flags = fcntl(transport_global_data.ctx->async_fd, F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get async_fd flags");
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.ctx->async_fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set async_fd to nonblocking");
        return NNTI_EIO;
    }

    flags = fcntl(transport_global_data.data_comp_channel->fd, F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get completion fd flags");
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.data_comp_channel->fd, F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set completion fd to nonblocking");
        return NNTI_EIO;
    }


    return(NNTI_OK);
}

static NNTI_result_t setup_interrupt_pipe(void)
{
    int rc=0;
    int flags;

    rc=pipe(transport_global_data.interrupt_pipe);
    if (rc < 0) {
        log_error(nnti_debug_level, "pipe() failed: %s", strerror(errno));
        return NNTI_EIO;
    }

    /* use non-blocking IO on the read side of the interrupt pipe */
    flags = fcntl(transport_global_data.interrupt_pipe[0], F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get interrupt_pipe flags: %s", strerror(errno));
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.interrupt_pipe[0], F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set interrupt_pipe to nonblocking: %s", strerror(errno));
        return NNTI_EIO;
    }

    /* use non-blocking IO on the write side of the interrupt pipe */
    flags = fcntl(transport_global_data.interrupt_pipe[1], F_GETFL);
    if (flags < 0) {
        log_error(nnti_debug_level, "failed to get interrupt_pipe flags: %s", strerror(errno));
        return NNTI_EIO;
    }
    if (fcntl(transport_global_data.interrupt_pipe[1], F_SETFL, flags | O_NONBLOCK) < 0) {
        log_error(nnti_debug_level, "failed to set interrupt_pipe to nonblocking: %s", strerror(errno));
        return NNTI_EIO;
    }

    return(NNTI_OK);
}

static NNTI_result_t setup_atomics(void)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    uint32_t atomics_bytes;

    int flags = IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_ATOMIC;

    log_debug(nnti_debug_level, "enter");

    atomics_bytes=config.min_atomics_vars * sizeof(int64_t);
    trios_start_timer(callTime);
    transport_global_data.atomics=(int64_t*)aligned_malloc(atomics_bytes);
    if (transport_global_data.atomics == NULL) {
        return(NNTI_ENOMEM);
    }
    memset(transport_global_data.atomics, 0, atomics_bytes);
    trios_stop_timer("malloc and memset", callTime);

    trios_start_timer(callTime);
    mr = ibv_reg_mr_wrapper(transport_global_data.pd, transport_global_data.atomics, atomics_bytes, flags);
    if (!mr) {
        log_error(
                nnti_debug_level,
                "failed to register memory region (IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ "
                "| IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_ATOMIC): %s",
                strerror(errno));
        return(NNTI_ENOMEM);
    }
    trios_stop_timer("register", callTime);

    transport_global_data.atomics_mr=mr;

    log_debug(nnti_debug_level, "exit (buf==%p, mr==%p, lkey %x, rkey %x)...", transport_global_data.atomics, mr, mr->lkey, mr->rkey);

    return(NNTI_OK);
}

static void *aligned_malloc(
        size_t size)
{
    int rc = 0;
    void * ptr = NULL;

#ifdef NO_MEMALIGN
      ptr = malloc( size );
#else
      rc = posix_memalign( &ptr , 128 , size );
#endif

    if (rc!=0) {
        log_error(nnti_debug_level, "posix_memalign() failed: rc=%d (%s)", rc, strerror(errno));
        ptr=NULL;
    }

    return ptr;
}

static struct ibv_mr *register_memory_segment(
        void *buf,
        uint64_t len,
        enum ibv_access_flags access)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    log_debug(nnti_debug_level, "enter buffer(%p) len(%d)", buf, len);

    if (config.use_mlock) {
        trios_start_timer(callTime);
        mlock(buf, len);
        trios_stop_timer("mlock", callTime);
    }

    trios_start_timer(callTime);
    mr = ibv_reg_mr_wrapper(transport_global_data.pd, buf, len, access);
    if (!mr) {
        if (errno == EFAULT) {
            log_debug(nnti_debug_level, "ibv_reg_mr failed with EFAULT.  trying to register with IBV_ACCESS_REMOTE_READ.");
            mr = ibv_reg_mr_wrapper(transport_global_data.pd, buf, len, IBV_ACCESS_REMOTE_READ);
            if (!mr) {
                log_error(nnti_debug_level, "failed to register memory region with IBV_ACCESS_REMOTE_READ: %s", strerror(errno));
                return(NULL);
            }
        } else {
            log_error(nnti_debug_level, "failed to register memory region (IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE): %s", strerror(errno));
            return(NULL);
        }
    }
    trios_stop_timer("register", callTime);

    log_debug(nnti_debug_level, "exit (buf==%p, mr==%p, lkey %x, rkey %x)...", buf, mr, mr->lkey, mr->rkey);

    return(mr);
}

static int register_ack(ib_work_request *ib_wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    uint32_t len;

    log_debug(nnti_debug_level, "enter");

    len = sizeof(ib_wr->ack);

    if (config.use_memset) {
        trios_start_timer(callTime);
        memset(&ib_wr->ack, 0, len);
        trios_stop_timer("memset", callTime);
    }
    if (config.use_mlock) {
        trios_start_timer(callTime);
        mlock(&ib_wr->ack, len);
        trios_stop_timer("mlock", callTime);
    }

    trios_start_timer(callTime);
    mr = ibv_reg_mr_wrapper(transport_global_data.pd, &ib_wr->ack, len, (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));
    if (!mr) {
        log_error(nnti_debug_level, "failed to register memory region");
        perror("errno");
        return errno;
    }
    trios_stop_timer("register", callTime);

    ib_wr->ack_mr=mr;

    log_debug(nnti_debug_level, "exit (mr==%p, addr %p, length %lu, lkey %x, rkey %x)...", mr, mr->addr, mr->length, mr->lkey, mr->rkey);

    return (rc);
}

static int unregister_memory_segment(struct ibv_mr *mr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    trios_start_timer(callTime);
    ibv_rc=ibv_dereg_mr_wrapper(mr);
    if (ibv_rc != 0) {
        log_error(nnti_debug_level, "deregistering the memory buffer failed");
    }
    trios_stop_timer("deregister", callTime);

    if (config.use_mlock) {
        trios_start_timer(callTime);
        munlock(mr->addr, mr->length);
        trios_stop_timer("munlock", callTime);
    }

    log_debug(nnti_debug_level, "exit");

    return (rc);
}

static int unregister_ack(ib_work_request *ib_wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    if (ib_wr->ack_mr!=NULL) {
        trios_start_timer(callTime);
        ibv_rc=ibv_dereg_mr_wrapper(ib_wr->ack_mr);
        if (ibv_rc != 0) {
            log_error(nnti_debug_level, "deregistering the ACK buffer failed");
        }
        trios_stop_timer("deregister", callTime);
    }

    if (config.use_mlock) {
        trios_start_timer(callTime);
        munlock(ib_wr->ack_mr->addr, ib_wr->ack_mr->length);
        trios_stop_timer("munlock", callTime);
    }

    log_debug(nnti_debug_level, "exit");

    return (rc);
}

static void send_ack (
        ib_work_request *ib_wr)
{
    struct ibv_send_wr *bad_wr;

    NNTI_buffer_t    *reg_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    reg_buf=ib_wr->reg_buf;
    assert(reg_buf);
    ib_mem_hdl=IB_MEM_HDL(reg_buf);
    assert(ib_mem_hdl);

    log_debug(nnti_debug_level, "sending ACK to (ack_sge.addr=%p, ack_sge.length=%llu, ack_sq_wr.ib_wr.rdma.rkey=%x, ack_sq_wr.ib_wr.rdma.remote_addr=%p)",
            (void *)  ib_wr->ack_sge.addr,
            (uint64_t)ib_wr->ack_sge.length,
                      ib_wr->ack_sq_wr.wr.rdma.rkey,
            (void *)  ib_wr->ack_sq_wr.wr.rdma.remote_addr);

    if (ibv_post_send_wrapper(ib_wr->qp, &ib_wr->ack_sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
    }

    return;
}

static ib_work_request *decode_work_request(
        const struct ibv_wc *wc)
{
    ib_connection   *conn=NULL;
    ib_work_request *ib_wr=NULL;
    NNTI_buffer_t   *event_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    if (wc->opcode==IBV_WC_RECV) {
        log_debug(nnti_debug_level, "wc->opcode is IBV_WC_RECV, so this is a new request.  wc.wr_id is a key to the work request map.");

//        ib_mem_hdl=IB_MEM_HDL(transport_global_data.req_queue.reg_buf);
//        assert(ib_mem_hdl);

        nthread_lock(&nnti_wrmap_lock);
        if (wrmap.find(wc->wr_id) != wrmap.end()) {
            ib_wr = wrmap[wc->wr_id];
        }
        nthread_unlock(&nnti_wrmap_lock);
        assert(ib_wr);

    } else if (wc->opcode==IBV_WC_SEND) {
        log_debug(nnti_debug_level, "wc->opcode is IBV_WC_SEND, so this is a send request (IB_OP_SEND_REQUEST).  wc.wr_id is a key to the work request map.");

        nthread_lock(&nnti_wrmap_lock);
        if (wrmap.find(wc->wr_id) != wrmap.end()) {
            ib_wr = wrmap[wc->wr_id];
        }
        nthread_unlock(&nnti_wrmap_lock);
        assert(ib_wr);

    } else {
        if (wc->imm_data == 0) {
            log_debug(nnti_debug_level, "wc->imm_data == 0, so I am the initiator.  wc.wr_id is a key to the work request map.");

            // This is not a request buffer and I am the initiator, so wc.wr_id is the hash of the work request
            nthread_lock(&nnti_wrmap_lock);
            if (wrmap.find(wc->wr_id) != wrmap.end()) {
                ib_wr = wrmap[wc->wr_id];
            }
            nthread_unlock(&nnti_wrmap_lock);
            assert(ib_wr);

        } else {
            log_debug(nnti_debug_level, "wc->imm_data != 0, so I am the target.  wc->imm_data is either the hash of a buffer or the hash of an ib_wr.");

            // This is not a request buffer and I am the target, so wc.imm_data is the hash of either the buffer or the work request
            event_buf = get_buf_bufhash(wc->imm_data);
            assert(event_buf);
            ib_mem_hdl=IB_MEM_HDL(event_buf);
            assert(ib_mem_hdl);
            ib_wr = ib_mem_hdl->wr_queue.front();
            assert(ib_wr);
        }
    }

    log_debug(nnti_debug_level, "exit (ib_wr==%p)", ib_wr);

    return(ib_wr);
}

static int cancel_wr(
        ib_work_request *ib_wr)
{
    NNTI_result_t rc=NNTI_OK;

    log_debug(nnti_debug_level, "enter (ib_wr=%p)", ib_wr);

    ib_wr->state=NNTI_IB_WR_STATE_WAIT_COMPLETE;
    ib_wr->nnti_wr->result=NNTI_ECANCELED;

    log_debug(nnti_debug_level, "exit (ib_wr==%p)", ib_wr);

    return(rc);
}

static int process_event(
        ib_work_request     *ib_wr,
        const struct ibv_wc *wc)
{
    NNTI_result_t rc=NNTI_OK;

    ib_memory_handle *ib_mem_hdl=NULL;

    log_level debug_level=nnti_debug_level;

    log_debug(nnti_debug_level, "enter (ib_wr=%p)", ib_wr);

    if ((ib_wr->nnti_wr) && (ib_wr->nnti_wr->ops == NNTI_ATOMICS)) {
        ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;
        ib_wr->nnti_wr->result=NNTI_OK;
        return NNTI_OK;
    }

    if (wc->status == IBV_WC_RNR_RETRY_EXC_ERR) {
        ib_wr->state =NNTI_IB_WR_STATE_WAIT_COMPLETE;
        ib_wr->nnti_wr->result=NNTI_EDROPPED;
        return NNTI_EDROPPED;
    } else if (wc->status != IBV_WC_SUCCESS) {
        ib_wr->state =NNTI_IB_WR_STATE_WAIT_COMPLETE;
        ib_wr->nnti_wr->result=NNTI_EIO;
        return NNTI_EIO;
    }

    ib_mem_hdl=IB_MEM_HDL(ib_wr->reg_buf);
    assert(ib_mem_hdl);

    log_debug(nnti_debug_level, "ib_wr=%p; ib_wr->last_op=%d", ib_wr, ib_wr->last_op);
    log_debug(nnti_debug_level, "ib_mem_hdl->type==%llu", (uint64_t)ib_mem_hdl->type);

    ib_wr->last_wc = *wc;

    switch (ib_mem_hdl->type) {
        case SEND_BUFFER:
            if (wc->opcode==IBV_WC_SEND) {
                log_debug(debug_level, "send completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;
                ib_wr->nnti_wr->result=NNTI_OK;
            }
            if (wc->opcode==IBV_WC_RDMA_WRITE) {
                log_debug(debug_level, "send completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;
                ib_wr->nnti_wr->result=NNTI_OK;
            }
            break;
        case PUT_SRC_BUFFER:
            log_debug(debug_level, "RDMA write event - wc==%p, ib_wr==%p, state==%d", wc, ib_wr, ib_wr->state);
            if (wc->opcode==IBV_WC_RDMA_WRITE) {
                if (wc->wc_flags==0) {
                    if (ib_wr->state==NNTI_IB_WR_STATE_STARTED) {
                        log_debug(debug_level, "RDMA write (initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                        ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;
                        if (!config.use_rdma_target_ack) {
                            ib_wr->nnti_wr->result=NNTI_OK;
                        }
                    }
                }
                if (wc->wc_flags==IBV_WC_WITH_IMM) {
                    if ((config.use_rdma_target_ack) &&
                        (ib_wr->state == NNTI_IB_WR_STATE_RDMA_COMPLETE)) {
                        log_debug(debug_level, "RDMA write ACK (initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                        ib_wr->last_op=IB_OP_PUT_INITIATOR;
                        ib_wr->state=NNTI_IB_WR_STATE_ACK_COMPLETE;
                        ib_wr->nnti_wr->result=NNTI_OK;
                    }
                }
            }
//            if (ib_wr->state == NNTI_IB_WR_STATE_RDMA_COMPLETE) {
//                print_xfer_buf((void *)ib_wr->reg_buf->payload, ib_wr->reg_buf->payload_size);
//                print_ack_buf(&ib_wr->ack);
//            }
            break;
        case GET_DST_BUFFER:
            log_debug(debug_level, "RDMA read event - wc==%p, ib_wr==%p, state==%d", wc, ib_wr, ib_wr->state);
            if ((wc->opcode==IBV_WC_RDMA_READ) &&
                (wc->wc_flags==0)) {
                if (ib_wr->state==NNTI_IB_WR_STATE_STARTED) {
                    log_debug(debug_level, "RDMA read (initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                    ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;
                    if (!config.use_rdma_target_ack) {
                        ib_wr->nnti_wr->result=NNTI_OK;
                    }
                }
            }
            else if ((config.use_rdma_target_ack) &&
                     ((wc->opcode==IBV_WC_RDMA_WRITE) &&
                      (wc->wc_flags==IBV_WC_WITH_IMM))) {
                if (ib_wr->state == NNTI_IB_WR_STATE_RDMA_COMPLETE) {
                    log_debug(debug_level, "RDMA read ACK (initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                    ib_wr->last_op=IB_OP_GET_INITIATOR;
                    ib_wr->state=NNTI_IB_WR_STATE_ACK_COMPLETE;
                    ib_wr->nnti_wr->result=NNTI_OK;
                }
            }
//            if (ib_wr->state == NNTI_IB_WR_STATE_RDMA_COMPLETE) {
//                print_xfer_buf((void *)ib_wr->reg_buf->payload, ib_wr->reg_buf->payload_size);
//                print_ack_buf(&ib_wr->ack);
//            }
            break;
        case REQUEST_BUFFER:
            if (wc->opcode==IBV_WC_RECV) {
                log_debug(debug_level, "recv completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->last_op=IB_OP_NEW_REQUEST;
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;

//                if (transport_global_data.req_queue.req_received == transport_global_data.srq_count) {
//                    log_debug(debug_level, "resetting req_queue.req_received to 0");
//                    transport_global_data.req_queue.req_received=0;
//                }
//                if (transport_global_data.req_queue.req_received != (ib_wr->offset/ib_wr->length)) {
//                    log_warn(debug_level, "req_queue.req_received(%llu) != (ib_wr->offset(%llu)/ib_wr->length(%llu))",
//                            transport_global_data.req_queue.req_received, ib_wr->offset, ib_wr->length);
//                }
                transport_global_data.req_queue.req_received++;

                if (ib_wr->cq == transport_global_data.req_cq) {
                    transport_global_data.req_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
                }
            }
            break;
        case RECEIVE_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "recv completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->last_op=IB_OP_RECEIVE;
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;

                if (ib_wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
            break;
        case PUT_DST_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "RDMA write (target) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->last_op=IB_OP_PUT_TARGET;
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;

                if (ib_wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
//            if (ib_wr->op_state == RDMA_WRITE_COMPLETE) {
//                print_xfer_buf((void *)ib_wr->reg_buf->payload, ib_wr->reg_buf->payload_size);
//                print_ack_buf(&ib_wr->ack);
//            }
            break;
        case GET_SRC_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "RDMA read (target) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                ib_wr->last_op=IB_OP_GET_TARGET;
                ib_wr->state=NNTI_IB_WR_STATE_RDMA_COMPLETE;

                if (ib_wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
//            if (ib_wr->op_state == RDMA_READ_COMPLETE) {
//                print_xfer_buf((void *)ib_wr->reg_buf->payload, ib_wr->reg_buf->payload_size);
//                print_ack_buf(&ib_wr->ack);
//            }
            break;
        case RDMA_TARGET_BUFFER:
            if (ib_wr->last_op==IB_OP_GET_INITIATOR) {
                if ((wc->opcode==IBV_WC_RDMA_READ) &&
                    (wc->wc_flags==0)) {
                    if (ib_wr->state==NNTI_IB_WR_STATE_STARTED) {
                        log_debug(debug_level, "RDMA target (read initiator) completion - ib_wr==%p", ib_wr);
                        ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE;
                        if (!config.use_rdma_target_ack) {
                            ib_wr->nnti_wr->result=NNTI_OK;
                        }
                    }
                }
                else if ((config.use_rdma_target_ack) &&
                         ((wc->opcode==IBV_WC_RDMA_WRITE) &&
                          (wc->wc_flags==IBV_WC_WITH_IMM))) {
                    if (ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE) {
                        log_debug(debug_level, "RDMA target ACK (read initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                        ib_wr->last_op=IB_OP_GET_INITIATOR;
                        ib_wr->state=NNTI_IB_WR_STATE_ACK_COMPLETE;
                        ib_wr->nnti_wr->result=NNTI_OK;
                    }
                }
            } else if (ib_wr->last_op==IB_OP_PUT_INITIATOR) {
                if (wc->opcode==IBV_WC_RDMA_WRITE) {
                    if (wc->wc_flags==0) {
                        if (ib_wr->state==NNTI_IB_WR_STATE_STARTED) {
                            log_debug(debug_level, "RDMA target (write initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                            ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE;
                            if (config.use_rdma_target_ack) {
                                ib_wr->nnti_wr->result=NNTI_OK;
                            }
                        }
                    }
                    if (wc->wc_flags==IBV_WC_WITH_IMM) {
                        if ((config.use_rdma_target_ack) &&
                            (ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE)) {
                            log_debug(debug_level, "RDMA target ACK (write initiator) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                            ib_wr->last_op=IB_OP_PUT_INITIATOR;
                            ib_wr->state=NNTI_IB_WR_STATE_ACK_COMPLETE;
                            ib_wr->nnti_wr->result=NNTI_OK;
                        }
                    }
                }
            } else {
                if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                    log_debug(debug_level, "RDMA target (target) completion - wc==%p, ib_wr==%p", wc, ib_wr);
                    ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE;

                    if (config.use_rdma_target_ack) {
                        ib_wr->last_op=ib_wr->ack.op;
                    }
                    if (ib_wr->cq == transport_global_data.data_cq) {
                        transport_global_data.data_srq_count--;
                        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                    }
                }
            }
//            if (ib_wr->op_state == RDMA_TARGET_COMPLETE) {
//                print_xfer_buf((void *)ib_wr->reg_buf->payload, ib_wr->reg_buf->payload_size);
//                print_ack_buf(&ib_wr->ack);
//            }
        case UNKNOWN_BUFFER:
        default:

            break;
    }

    return (rc);
}

static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         wr_id,
        uint64_t        offset)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    struct ibv_srq *srq=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=IB_MEM_HDL(reg_buf);
    assert(ib_mem_hdl);

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        nthread_lock_init(&ib_wr->lock);
    } else {
        if (config.use_wr_pool) {
            ib_wr=wr_pool_sendrecv_pop();
        } else {
            ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
            nthread_lock_init(&ib_wr->lock);
        }
    }
    assert(ib_wr);
    log_debug(nnti_debug_level, "allocated ib_wr (ib_wr=%p)", ib_wr);
    ib_wr->reg_buf = reg_buf;
    ib_wr->offset  = offset;
    ib_wr->length  = reg_buf->payload_size;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_POSTED;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        ib_wr->last_op=IB_OP_NEW_REQUEST;
    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq                =transport_global_data.req_srq;
        ib_wr->comp_channel=transport_global_data.req_comp_channel;
        ib_wr->cq          =transport_global_data.req_cq;
    } else {
        srq                =transport_global_data.data_srq;
        ib_wr->comp_channel=transport_global_data.data_comp_channel;
        ib_wr->cq          =transport_global_data.data_cq;
    }

    ib_wr->sge_count=reg_buf->buffer_segments.NNTI_remote_addr_array_t_len;

    if (ib_wr->sge_count == 1) {

        ib_wr->sge_list=&ib_wr->sge;

        ib_wr->sge_list[0].addr  =(uint64_t)ib_mem_hdl->mr_list[0]->addr+offset;
        ib_wr->sge_list[0].length=reg_buf->payload_size;
        ib_wr->sge_list[0].lkey  =ib_mem_hdl->mr_list[0]->lkey;

    } else {

        ib_wr->sge_list=(struct ibv_sge *)calloc(ib_wr->sge_count, sizeof(struct ibv_sge));

        for (int i=0;i<ib_wr->sge_count;i++) {
            ib_wr->sge_list[i].addr  =(uint64_t)ib_mem_hdl->mr_list[i]->addr;
            ib_wr->sge_list[i].length=ib_mem_hdl->mr_list[i]->length;
            ib_wr->sge_list[i].lkey  =ib_mem_hdl->mr_list[i]->lkey;
        }

    }

    if (wr_id == -1) {
        ib_wr->rq_wr.wr_id=(uint64_t)ib_wr->key;
    } else {
        ib_wr->rq_wr.wr_id=wr_id;
    }
    ib_wr->rq_wr.sg_list=ib_wr->sge_list;
    ib_wr->rq_wr.num_sge=ib_wr->sge_count;

    log_debug(nnti_debug_level, "sge_count=%d, sge_list=%p", ib_wr->sge_count, ib_wr->sge_list);

    if ((ib_wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((ib_wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t post_ack_recv_work_request(
        NNTI_buffer_t  *reg_buf)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;
    struct ibv_srq   *srq;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=IB_MEM_HDL(reg_buf);
    assert(ib_mem_hdl);

    if (config.use_wr_pool) {
        ib_wr=wr_pool_rdma_pop();
    } else {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        nthread_lock_init(&ib_wr->lock);
    }
    assert(ib_wr);
    log_debug(nnti_debug_level, "allocated ib_wr (ib_wr=%p)", ib_wr);
    ib_wr->reg_buf = reg_buf;

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);

    ib_wr->state=NNTI_IB_WR_STATE_POSTED;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        ib_wr->last_op=IB_OP_NEW_REQUEST;
    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        ib_wr->comp_channel=transport_global_data.req_comp_channel;
        ib_wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        ib_wr->comp_channel=transport_global_data.data_comp_channel;
        ib_wr->cq          =transport_global_data.data_cq;
    }

    if (!config.use_wr_pool) {
        register_ack(ib_wr);
    }

    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_size = ib_wr->ack_mr->length;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_buf  = (uint64_t)ib_wr->ack_mr->addr;
    reg_buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.ack_key  = ib_wr->ack_mr->rkey;

    ib_wr->sge.addr  =(uint64_t)ib_wr->ack_mr->addr;
    ib_wr->sge.length=ib_wr->ack_mr->length;
    ib_wr->sge.lkey  =ib_wr->ack_mr->lkey;

    ib_wr->rq_wr.wr_id  =(uint64_t)ib_wr->key;
    ib_wr->rq_wr.sg_list=&ib_wr->sge;
    ib_wr->rq_wr.num_sge=1;

    if ((ib_wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((ib_wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    nthread_lock(&nnti_wrmap_lock);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_recv_work_request(
        NNTI_work_request_t  *wr,
        ib_work_request      *ib_wr)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;

    struct ibv_srq *srq=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p ; ib_wr=%p)", wr, ib_wr);

    assert(wr);
    assert(ib_wr);
    ib_mem_hdl=IB_MEM_HDL(wr->reg_buf);
    assert(ib_mem_hdl);

    nthread_lock(&nnti_wrmap_lock);
    wrmap_iter_t m_victim=wrmap.find(ib_wr->key);
    if (m_victim != wrmap.end()) {
        log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", ib_wr, ib_wr->key);
        wrmap.erase(m_victim);
    }

    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);
    ib_wr->rq_wr.wr_id=(uint64_t)ib_wr->key;

    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    ib_wr->state=NNTI_IB_WR_STATE_POSTED;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        ib_wr->last_op=IB_OP_NEW_REQUEST;
    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        ib_wr->comp_channel=transport_global_data.req_comp_channel;
        ib_wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        ib_wr->comp_channel=transport_global_data.data_comp_channel;
        ib_wr->cq          =transport_global_data.data_cq;
    }

    if ((ib_wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((ib_wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}

static NNTI_result_t repost_ack_recv_work_request(
        NNTI_work_request_t *wr,
        ib_work_request     *ib_wr)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;
    struct ibv_srq   *srq;

    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ib_wr=IB_WORK_REQUEST(wr);
    assert(ib_wr);
    ib_mem_hdl=IB_MEM_HDL(wr->reg_buf);
    assert(ib_mem_hdl);

    ib_wr->state=NNTI_IB_WR_STATE_POSTED;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        ib_wr->last_op=IB_OP_NEW_REQUEST;
    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        ib_wr->comp_channel=transport_global_data.req_comp_channel;
        ib_wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        ib_wr->comp_channel=transport_global_data.data_comp_channel;
        ib_wr->cq          =transport_global_data.data_cq;
    }

    if ((ib_wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((ib_wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        log_debug(nnti_debug_level, "posting ib_wr=%p (key=%lx)", ib_wr, ib_wr->key);
        ibv_rc=ibv_post_srq_recv_wrapper(srq, &ib_wr->rq_wr, &bad_wr);
        if (ibv_rc) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &ib_wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    log_debug(nnti_debug_level, "pushing ib_wr=%p", ib_wr);
    nthread_lock(&ib_mem_hdl->wr_queue_lock);
    ib_mem_hdl->wr_queue.push_back(ib_wr);
    nthread_unlock(&ib_mem_hdl->wr_queue_lock);

    nthread_lock(&nnti_wrmap_lock);
    wrmap_iter_t m_victim=wrmap.find(ib_wr->key);
    if (m_victim != wrmap.end()) {
        log_debug(nnti_debug_level, "erasing ib_wr=%p (key=%lx) from the wrmap", ib_wr, ib_wr->key);
        wrmap.erase(m_victim);
    }
    ib_wr->key = nthread_counter_increment(&nnti_wrmap_counter);
    log_debug(nnti_debug_level, "wrmap[key(%lx)]=ib_wr(%p)", ib_wr->key, ib_wr);
    assert(wrmap.find(ib_wr->key) == wrmap.end());
    wrmap[ib_wr->key] = ib_wr;
    nthread_unlock(&nnti_wrmap_lock);

    log_debug(nnti_debug_level, "exit (wr=%p)", wr);

    return(NNTI_OK);
}

static int8_t is_wr_canceling(
        ib_work_request *ib_wr)
{
    int8_t rc=FALSE;

    nthread_lock(&ib_wr->lock);
    if (ib_wr->state == NNTI_IB_WR_STATE_CANCELING) {
        rc=TRUE;
    }
    nthread_unlock(&ib_wr->lock);

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_wr_complete(
        ib_work_request *ib_wr)
{
    int8_t rc=FALSE;

    nthread_lock(&ib_wr->lock);
    if ((ib_wr->state==NNTI_IB_WR_STATE_RDMA_COMPLETE) ||
        (ib_wr->state==NNTI_IB_WR_STATE_WAIT_COMPLETE)) {
        rc=TRUE;
    }
    nthread_unlock(&ib_wr->lock);

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);

    return(rc);
}

static int8_t is_any_wr_complete(
        ib_work_request **wr_list,
        const uint32_t    wr_count,
        uint32_t         *which)
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
        ib_work_request **wr_list,
        const uint32_t    wr_count)
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
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter (wr=%p)", wr);

    ib_wr=IB_WORK_REQUEST(wr);
    assert(ib_wr);

    return(is_wr_complete(ib_wr));
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

static void create_status(
        NNTI_work_request_t *wr,
        ib_work_request     *ib_wr,
        NNTI_result_t        nnti_rc,
        NNTI_status_t       *status)
{
    ib_connection    *conn      =NULL;

    NNTI_buffer_t *reg_buf;

    trios_declare_timer(total_time);

    assert(wr);
    assert(ib_wr);
    assert(status);

    trios_start_timer(total_time);

    memset(status, 0, sizeof(NNTI_status_t));
    status->op     = wr->ops;
    if (is_wr_complete(ib_wr)) {
        status->result = wr->result;
    } else {
        status->result = nnti_rc;
    }
    if (status->result==NNTI_OK) {

//        print_ack_buf(&ib_wr->ack);
//        print_wr(ib_wr);

        conn = get_conn_qpn(ib_wr->last_wc.qp_num);

        if (ib_wr->nnti_wr->ops != NNTI_ATOMICS) {
            status->start  = (uint64_t)ib_wr->reg_buf->payload;
        }
        switch (ib_wr->last_op) {
            case IB_OP_PUT_INITIATOR:
            case IB_OP_SEND_REQUEST:
            case IB_OP_SEND_BUFFER:
                create_peer(&status->src, transport_global_data.listen_name, transport_global_data.listen_addr, transport_global_data.listen_port);
                create_peer(&status->dest, conn->peer_name, conn->peer_addr, conn->peer_port);
                break;
            case IB_OP_GET_TARGET:
                if (config.use_rdma_target_ack) {
                    create_peer(&status->src, transport_global_data.listen_name, transport_global_data.listen_addr, transport_global_data.listen_port);
                    create_peer(&status->dest, conn->peer_name, conn->peer_addr, conn->peer_port);
                }
                break;
            case IB_OP_GET_INITIATOR:
            case IB_OP_NEW_REQUEST:
            case IB_OP_RECEIVE:
                create_peer(&status->src, conn->peer_name, conn->peer_addr, conn->peer_port);
                create_peer(&status->dest, transport_global_data.listen_name, transport_global_data.listen_addr, transport_global_data.listen_port);
                break;
            case IB_OP_PUT_TARGET:
                if (config.use_rdma_target_ack) {
                    create_peer(&status->src, conn->peer_name, conn->peer_addr, conn->peer_port);
                    create_peer(&status->dest, transport_global_data.listen_name, transport_global_data.listen_addr, transport_global_data.listen_port);
                }
                break;
        }
        switch (ib_wr->last_op) {
            case IB_OP_NEW_REQUEST:
                status->offset = ib_wr->offset;
                status->length = ib_wr->last_wc.byte_len;
                break;
            case IB_OP_SEND_REQUEST:
            case IB_OP_SEND_BUFFER:
            case IB_OP_RECEIVE:
                status->offset = 0;
                status->length = ib_wr->last_wc.byte_len;
                break;
            case IB_OP_GET_INITIATOR:
            case IB_OP_PUT_INITIATOR:
                status->offset = ib_wr->offset;
                status->length = ib_wr->length;
                break;
            case IB_OP_GET_TARGET:
            case IB_OP_PUT_TARGET:
                if (config.use_rdma_target_ack) {
                    status->offset = ib_wr->ack.offset;
                    status->length = ib_wr->ack.length;
                }
                break;
        }
    }

    trios_stop_timer("create_status - total", total_time);
}

static void create_peer(NNTI_peer_t *peer, char *name, NNTI_ip_addr addr, NNTI_tcp_port port)
{
    log_debug(nnti_debug_level, "enter");

    sprintf(peer->url, "ib://%s:%u/", name, ntohs(port));

    peer->peer.transport_id                   =NNTI_TRANSPORT_IB;
    peer->peer.NNTI_remote_process_t_u.ib.addr=addr;
    peer->peer.NNTI_remote_process_t_u.ib.port=port;

    log_debug(nnti_debug_level, "exit");
}

//static void copy_peer(NNTI_peer_t *src, NNTI_peer_t *dest)
//{
//    log_debug(nnti_debug_level, "enter");
//
//    strncpy(dest->url, src->url, NNTI_URL_LEN);
//    dest->url[NNTI_URL_LEN-1]='\0';
//
//    src->peer.transport_id                    =NNTI_TRANSPORT_IB;
//    dest->peer.NNTI_remote_process_t_u.ib.addr=src->peer.NNTI_remote_process_t_u.ib.addr;
//    dest->peer.NNTI_remote_process_t_u.ib.port=src->peer.NNTI_remote_process_t_u.ib.port;
//
//    log_debug(nnti_debug_level, "exit");
//}

static int init_server_listen_socket()
{
    NNTI_result_t rc=NNTI_OK;
    int flags;
    struct hostent *host_entry;
    struct sockaddr_in skin;
    socklen_t skin_size=sizeof(struct sockaddr_in);

    transport_global_data.listen_sock = socket(AF_INET, SOCK_STREAM, 0);
    if (transport_global_data.listen_sock < 0)
        log_error(nnti_debug_level, "failed to create tcp socket: %s", strerror(errno));

    flags = 1;
    if (setsockopt(transport_global_data.listen_sock, SOL_SOCKET, SO_REUSEADDR, &flags, sizeof(flags)) < 0)
        log_error(nnti_debug_level, "failed to set tcp socket REUSEADDR flag: %s", strerror(errno));

    flags=fcntl(transport_global_data.listen_sock, F_GETFL, 0);
    fcntl(transport_global_data.listen_sock, F_SETFL, flags | O_NONBLOCK);

    if (transport_global_data.listen_name[0]!='\0') {
        log_debug(nnti_debug_level, "using hostname from command-line (%s).", transport_global_data.listen_name);
    } else {
        gethostname(transport_global_data.listen_name, NNTI_HOSTNAME_LEN);
        log_debug(nnti_debug_level, "hostname not given on command-line.  using gethostname() result (%s).", transport_global_data.listen_name);
    }
    /* lookup the host provided on the command line */
    host_entry = gethostbyname(transport_global_data.listen_name);
    if (!host_entry) {
        log_warn(nnti_debug_level, "failed to resolve server name (%s): %s", transport_global_data.listen_name, strerror(errno));
        return NNTI_ENOENT;
    }

    memset(&skin, 0, sizeof(skin));
    skin.sin_family = AF_INET;
    memcpy(&skin.sin_addr, host_entry->h_addr_list[0], (size_t) host_entry->h_length);
    /* 0 here means to bind to a random port assigned by the kernel */
    skin.sin_port = 0;

retry:
    if (bind(transport_global_data.listen_sock, (struct sockaddr *) &skin, sizeof(skin)) < 0) {
        if (errno == EINTR) {
            goto retry;
        } else {
            log_error(nnti_debug_level, "failed to bind tcp socket: %s", strerror(errno));
        }
    }
    /* after the bind, get the "name" for the socket.  the "name" contains the port assigned by the kernel. */
    getsockname(transport_global_data.listen_sock, (struct sockaddr *)&skin, &skin_size);
    transport_global_data.listen_addr = (uint32_t)skin.sin_addr.s_addr;
    transport_global_data.listen_port = (uint16_t)skin.sin_port;
    log_debug(nnti_debug_level, "listening on ip(%s) addr(%u) port(%u)",
            transport_global_data.listen_name,
            (unsigned int)ntohl(skin.sin_addr.s_addr),
            (unsigned int)ntohs(skin.sin_port));
    if (listen(transport_global_data.listen_sock, 1024) < 0)
        log_error(nnti_debug_level, "failed to listen on tcp socket: %s", strerror(errno));

    return rc;
}


static void transition_connection_to_ready(
        int sock,
        ib_connection *conn)
{
    int rc=NNTI_OK;
    int min_rnr_timer;
    int ack_timeout;
    int retry_count;
    trios_declare_timer(callTime);

    /* bring the two QPs up to RTR */
    trios_start_timer(callTime);

    min_rnr_timer=1;  /* means 0.01ms delay before sending RNR NAK */
    ack_timeout  =17; /* time to wait for ACK/NAK before retransmitting.  4.096us * 2^17 == 0.536s */
    if (config.drop_if_full_queue) {
        retry_count=1;  /* number of retries if no answer on primary path or if remote sends RNR NAK */
    } else {
        retry_count=7;   /* number of retries if no answer on primary path or if remote sends RNR NAK.  7 has special meaning of infinite retries. */
    }
    transition_qp_from_reset_to_ready(conn->req_qp.qp, conn->peer_req_qpn, conn->peer_lid, min_rnr_timer, ack_timeout, retry_count);


    min_rnr_timer=31;  /* means 491.52ms delay before sending RNR NAK */
    ack_timeout  =17;  /* time to wait for ACK/NAK before retransmitting.  4.096us * 2^17 == 0.536s */
    retry_count  =7;   /* number of retries if no answer on primary path or if remote sends RNR NAK.  7 has special meaning of infinite retries. */
    transition_qp_from_reset_to_ready(conn->data_qp.qp, conn->data_qp.peer_qpn, conn->peer_lid, min_rnr_timer, ack_timeout, retry_count);

    trios_stop_timer("transition_qp_from_reset_to_ready", callTime);

    trios_start_timer(callTime);
    /* final sychronization to ensure both sides have posted RTRs */
    rc = tcp_exchange(sock, 0, &rc, &rc, sizeof(rc));
    trios_stop_timer("exch data", callTime);
}

static void transition_qp_from_reset_to_ready(
        struct ibv_qp *qp,
        uint32_t       peer_qpn,
        int            peer_lid,
        int            min_rnr_timer,
        int            ack_timeout,
        int            retry_count)
{
    enum ibv_qp_attr_mask mask;
    struct ibv_qp_attr attr;

    log_debug(nnti_debug_level, "enter (qp=%p ; qp->qp_num=%u ; peer_qpn=%u ; peer_lid=%u)", qp, qp->qp_num, peer_qpn, peer_lid);

    /* Transition QP to Init */
    mask = (enum ibv_qp_attr_mask)
       ( IBV_QP_STATE
     | IBV_QP_ACCESS_FLAGS
     | IBV_QP_PKEY_INDEX
     | IBV_QP_PORT );
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_INIT;
    attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE |
            IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ |
            IBV_ACCESS_REMOTE_ATOMIC;
    attr.pkey_index = 0;
    attr.port_num = transport_global_data.nic_port;
    if (ibv_modify_qp_wrapper(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp from RESET to INIT state");
    }

    /* Transition QP to Ready-to-Receive (RTR) */
    mask = (enum ibv_qp_attr_mask)
       ( IBV_QP_STATE
     | IBV_QP_MAX_DEST_RD_ATOMIC
     | IBV_QP_AV
     | IBV_QP_PATH_MTU
     | IBV_QP_RQ_PSN
     | IBV_QP_DEST_QPN
     | IBV_QP_MIN_RNR_TIMER );
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_RTR;
    attr.max_dest_rd_atomic = 1;
    attr.ah_attr.dlid = peer_lid;
    attr.ah_attr.port_num = transport_global_data.nic_port;
    attr.path_mtu = IBV_MTU_1024;
    attr.rq_psn = 0;
    attr.dest_qp_num = peer_qpn;
    attr.min_rnr_timer = min_rnr_timer; /* delay before sending RNR NAK */
    if (ibv_modify_qp_wrapper(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp from INIT to RTR state");
    }

    /* Transition QP to Ready-to-Send (RTS) */
    mask = (enum ibv_qp_attr_mask)
       ( IBV_QP_STATE
     | IBV_QP_SQ_PSN
     | IBV_QP_MAX_QP_RD_ATOMIC
     | IBV_QP_TIMEOUT
     | IBV_QP_RETRY_CNT
     | IBV_QP_RNR_RETRY );
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_RTS;
    attr.sq_psn = 0;
    attr.max_rd_atomic = 1;
    attr.timeout = ack_timeout;   /* time to wait for ACK/NAK before retransmitting.  4.096us * 2^ack_timeout */
    attr.retry_cnt = retry_count; /* number of retries if no answer on primary path */
    attr.rnr_retry = retry_count; /* number of retries if remote sends RNR NAK */
    if (ibv_modify_qp_wrapper(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp from RTR to RTS state");
    }

    log_debug(nnti_debug_level, "exit");
}

static void transition_qp_from_error_to_ready(
        struct ibv_qp *qp,
        uint32_t       peer_qpn,
        int            peer_lid,
        int            min_rnr_timer,
        int            ack_timeout,
        int            retry_count)
{
    enum ibv_qp_attr_mask mask;
    struct ibv_qp_attr attr;

    log_debug(nnti_debug_level, "enter (qp=%p ; qp->qp_num=%u ; peer_qpn=%u ; peer_lid=%u)", qp, qp->qp_num, peer_qpn, peer_lid);

    /* Transition QP to Reset */
    mask = (enum ibv_qp_attr_mask)
       ( IBV_QP_STATE );
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_RESET;
    if (ibv_modify_qp_wrapper(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp from ERROR to RESET state");
    }

    transition_qp_from_reset_to_ready(qp, peer_qpn, peer_lid, min_rnr_timer, ack_timeout, retry_count);

    log_debug(nnti_debug_level, "exit");
}

static void transition_connection_to_error(
        ib_connection *conn)
{
//    int i;
    trios_declare_timer(callTime);

    /* bring the two QPs up to RTR */
    trios_start_timer(callTime);
    transition_qp_to_error(conn->req_qp.qp, conn->peer_req_qpn, conn->peer_lid);
    transition_qp_to_error(conn->data_qp.qp, conn->data_qp.peer_qpn, conn->peer_lid);
    trios_stop_timer("transition_qp_to_error", callTime);
}

static void transition_qp_to_error(
        struct ibv_qp *qp,
        uint32_t peer_qpn,
        int peer_lid)
{
    struct ibv_qp_attr attr;

    /* Transition QP to Error */
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_ERR;
    if (ibv_modify_qp_wrapper(qp, &attr, IBV_QP_STATE)) {
        log_error(nnti_debug_level, "failed to modify qp to ERROR state");
    }
}


/*
 * Try hard to read the whole buffer.  Abort on read error.
 */
static int tcp_read(int sock, void *incoming, size_t len)
{
    int bytes_this_read=0;
    int bytes_left=len;
    int bytes_read=0;

    while (bytes_left > 0) {
        bytes_this_read = read(sock, (char *)incoming + bytes_read, bytes_left);
        if (bytes_this_read < 0) {
            return bytes_this_read;
        }
        if (bytes_this_read == 0) {
            break;
        }
        bytes_left -= bytes_this_read;
        bytes_read += bytes_this_read;
    }
    return bytes_read;
}

/*
 * Try hard to write the whole buffer.  Abort on write error.
 */
static int tcp_write(int sock, const void *outgoing, size_t len)
{
    int bytes_this_write=0;
    int bytes_left=len;
    int bytes_written=0;

    while (bytes_left > 0) {
        bytes_this_write = write(sock, (const char *)outgoing + bytes_written, bytes_left);
        if (bytes_this_write < 0) {
            return bytes_this_write;
        }
        bytes_left    -= bytes_this_write;
        bytes_written += bytes_this_write;
    }
    return bytes_written;
}

/*
 * Two processes exchange data over a TCP socket.  Both sides send and receive the
 * same amount of data.  Only one process can declare itself the server (is_server!=0),
 * otherwise this will hang because both will wait for the read to complete.
 *
 * Server receives, then sends.
 * Client sends, then receives.
 */
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len)
{
    int rc;

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to read IB connection info: errno=%d", errno);
            goto out;
        }
        if (rc != (int) len) {
            log_error(nnti_debug_level, "partial read, %d/%d bytes", rc, (int) len);
            rc = 1;
            goto out;
        }
    } else {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_write(sock, outgoing, len);
        trios_stop_timer("tcp_write", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "client failed to write IB connection info: errno=%d", errno);
            goto out;
        }
    }

    if (is_server) {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_write(sock, outgoing, len);
        trios_stop_timer("tcp_write", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "server failed to write IB connection info: errno=%d", errno);
            goto out;
        }
    } else {
        trios_declare_timer(callTime);
        trios_start_timer(callTime);
        rc = tcp_read(sock, incoming, len);
        trios_stop_timer("tcp_read", callTime);
        if (rc < 0) {
            log_warn(nnti_debug_level, "client failed to read IB connection info: errno=%d", errno);
            goto out;
        }
        if (rc != (int) len) {
            log_error(nnti_debug_level, "partial read, %d/%d bytes", rc, (int) len);
            rc = 1;
            goto out;
        }
    }

    rc = 0;

out:
    return rc;
}

static int new_client_connection(
        ib_connection *c,
        int sock)
{
    int rc;
//    int flags=0;
    struct ibv_qp_init_attr att;

    /*
     * IB parameters passed through TCP to establish the IB connection.
     */
    struct {
        char     name[NNTI_HOSTNAME_LEN];
        uint32_t addr;
        uint32_t port;
        uint32_t lid;
        uint32_t req_qpn;
        uint32_t my_qpn;
        uint32_t peer_qpn;
        uint32_t atomics_rkey;
        uint64_t atomics_addr;
    } param_in, param_out;

    trios_declare_timer(callTime);

    // Initialize to avoid an "uninitialized bytes messages from valgrind
    memset(&param_out, 0, sizeof(param_out));
    memset(&param_in, 0, sizeof(param_in));

    strcpy(param_out.name, transport_global_data.listen_name);
    param_out.addr = htonl((uint32_t)transport_global_data.listen_addr);
    param_out.port = htonl((uint32_t)transport_global_data.listen_port);

    /* create the main queue pair */
    memset(&att, 0, sizeof(att));
    att.qp_context       = c;
    att.send_cq          = transport_global_data.req_cq;
    att.recv_cq          = transport_global_data.req_cq;
    att.srq              = transport_global_data.req_srq;
    att.cap.max_recv_wr  = transport_global_data.qp_count;
    att.cap.max_send_wr  = transport_global_data.qp_count;
    att.cap.max_recv_sge = 32;
    att.cap.max_send_sge = 32;
    att.qp_type          = IBV_QPT_RC;

    trios_start_timer(callTime);
    c->req_qp.qp = ibv_create_qp_wrapper(transport_global_data.pd, &att);
    if (!c->req_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    c->req_qp.qpn = c->req_qp.qp->qp_num;
    trios_stop_timer("create_qp", callTime);

    /* exchange data, converting info to network order and back */
    param_out.lid     = htonl(transport_global_data.nic_lid);
    param_out.req_qpn = htonl(c->req_qp.qpn);

    memset(&att, 0, sizeof(att));
    att.qp_context       = c;
    att.send_cq          = transport_global_data.data_cq;
    att.recv_cq          = transport_global_data.data_cq;
    att.srq              = transport_global_data.data_srq;
    att.cap.max_recv_wr  = transport_global_data.qp_count;
    att.cap.max_send_wr  = transport_global_data.qp_count;
    att.cap.max_recv_sge = 32;
    att.cap.max_send_sge = 32;
    att.qp_type          = IBV_QPT_RC;
    c->data_qp.qp = ibv_create_qp_wrapper(transport_global_data.pd, &att);
    if (!c->data_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    if (ibv_req_notify_cq_wrapper(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
    }

    c->data_qp.qpn     = c->data_qp.qp->qp_num;
    param_out.my_qpn = htonl(c->data_qp.qpn);

    param_out.atomics_rkey=transport_global_data.atomics_mr->rkey;
    param_out.atomics_addr=(uint64_t)transport_global_data.atomics_mr->addr;

    trios_start_timer(callTime);
    rc = tcp_exchange(sock, 0, &param_in, &param_out, sizeof(param_in));
    trios_stop_timer("exch data", callTime);
    if (rc)
        goto out;

    c->peer_name    = strdup(param_in.name);
    c->peer_addr    = ntohl(param_in.addr);
    c->peer_port    = ntohl(param_in.port);
    c->peer_lid     = ntohl(param_in.lid);
    c->peer_req_qpn = ntohl(param_in.req_qpn);
    c->data_qp.peer_qpn = ntohl(param_in.my_qpn);
    c->atomics_rkey = param_in.atomics_rkey;
    c->atomics_addr = param_in.atomics_addr;

out:
    return rc;
}

static int new_server_connection(
        ib_connection *c,
        int sock)
{
    int rc;
//    int flags=0;
    struct ibv_qp_init_attr att;

    /*
     * IB parameters passed through TCP to establish the IB connection.
     */
    struct {
        char     name[NNTI_HOSTNAME_LEN];
        uint32_t addr;
        uint32_t port;
        uint32_t lid;
        uint32_t req_qpn;
        uint32_t my_qpn;
        uint32_t peer_qpn;
        uint32_t atomics_rkey;
        uint64_t atomics_addr;
    } param_in, param_out;

    // initialize structs to avoid valgrind warnings
    memset(&param_out, 0, sizeof(param_out));
    memset(&param_in, 0, sizeof(param_in));

    trios_declare_timer(callTime);

    strcpy(param_out.name, transport_global_data.listen_name);
    param_out.addr = htonl((uint32_t)transport_global_data.listen_addr);
    param_out.port = htonl((uint32_t)transport_global_data.listen_port);

    /* create the main queue pair */
    memset(&att, 0, sizeof(att));
    att.qp_context       = c;
    att.send_cq          = transport_global_data.req_cq;
    att.recv_cq          = transport_global_data.req_cq;
    att.srq              = transport_global_data.req_srq;
    att.cap.max_recv_wr  = transport_global_data.qp_count;
    att.cap.max_send_wr  = transport_global_data.qp_count;
    att.cap.max_recv_sge = 32;
    att.cap.max_send_sge = 32;
    att.qp_type          = IBV_QPT_RC;

    trios_start_timer(callTime);
    c->req_qp.qp = ibv_create_qp_wrapper(transport_global_data.pd, &att);
    if (!c->req_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    c->req_qp.qpn = c->req_qp.qp->qp_num;
    trios_stop_timer("create_qp", callTime);

    /* exchange data, converting info to network order and back */
    param_out.lid     = htonl(transport_global_data.nic_lid);
    param_out.req_qpn = htonl(c->req_qp.qpn);

    memset(&att, 0, sizeof(att));
    att.qp_context       = c;
    att.send_cq          = transport_global_data.data_cq;
    att.recv_cq          = transport_global_data.data_cq;
    att.srq              = transport_global_data.data_srq;
    att.cap.max_recv_wr  = transport_global_data.qp_count;
    att.cap.max_send_wr  = transport_global_data.qp_count;
    att.cap.max_recv_sge = 32;
    att.cap.max_send_sge = 32;
    att.qp_type          = IBV_QPT_RC;
    c->data_qp.qp = ibv_create_qp_wrapper(transport_global_data.pd, &att);
    if (!c->data_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    if (ibv_req_notify_cq_wrapper(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
    }

    c->data_qp.qpn   = c->data_qp.qp->qp_num;
    param_out.my_qpn = htonl(c->data_qp.qpn);

    param_out.atomics_rkey=transport_global_data.atomics_mr->rkey;
    param_out.atomics_addr=(uint64_t)transport_global_data.atomics_mr->addr;

    trios_start_timer(callTime);
    rc = tcp_exchange(sock, 1, &param_in, &param_out, sizeof(param_in));
    trios_stop_timer("exch data", callTime);
    if (rc)
        goto out;

    c->peer_name    = strdup(param_in.name);
    c->peer_addr    = ntohl(param_in.addr);
    c->peer_port    = ntohl(param_in.port);
    c->peer_lid     = ntohl(param_in.lid);
    c->peer_req_qpn = ntohl(param_in.req_qpn);
    c->data_qp.peer_qpn = ntohl(param_in.my_qpn);
    c->atomics_rkey = param_in.atomics_rkey;
    c->atomics_addr = param_in.atomics_addr;

out:
    return rc;
}

static NNTI_result_t insert_conn_peer(const NNTI_peer_t *peer, ib_connection *conn)
{
    NNTI_result_t  rc=NNTI_OK;
    addrport_key key;

    key.addr = peer->peer.NNTI_remote_process_t_u.ib.addr;
    key.port = peer->peer.NNTI_remote_process_t_u.ib.port;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "insert_conn_peer", peer);
    }

    nthread_lock(&nnti_conn_peer_lock);
//    assert(connections_by_peer.find(key) == connections_by_peer.end());
    if (connections_by_peer.find(key) == connections_by_peer.end()) {
        connections_by_peer[key] = conn;   // add to connection map
    }
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(nnti_debug_level, "peer connection added (conn=%p)", conn);

//    print_peer_map();

    return(rc);
}
static NNTI_result_t insert_conn_qpn(const NNTI_qp_num qpn, ib_connection *conn)
{
    NNTI_result_t  rc=NNTI_OK;

    nthread_lock(&nnti_conn_qpn_lock);
    assert(connections_by_qpn.find(qpn) == connections_by_qpn.end());
    connections_by_qpn[qpn] = conn;
    nthread_unlock(&nnti_conn_qpn_lock);

    log_debug(nnti_debug_level, "qpn connection added (conn=%p)", conn);

//    print_qpn_map();

    return(rc);
}
static ib_connection *get_conn_peer(const NNTI_peer_t *peer)
{
    ib_connection *conn = NULL;

    addrport_key   key;

    assert(peer);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "get_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.ib.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.ib.port;

    nthread_lock(&nnti_conn_peer_lock);
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }
    nthread_unlock(&nnti_conn_peer_lock);

//    print_peer_map();

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found");
//    print_peer_map();

    return(NULL);
}
static ib_connection *get_conn_qpn(const NNTI_qp_num qpn)
{
    ib_connection *conn=NULL;

    log_debug(nnti_debug_level, "looking for qpn=%llu", (unsigned long long)qpn);
    nthread_lock(&nnti_conn_qpn_lock);
    if (connections_by_qpn.find(qpn) != connections_by_qpn.end()) {
        conn = connections_by_qpn[qpn];
    }
    nthread_unlock(&nnti_conn_qpn_lock);

//    print_qpn_map();

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        return conn;
    }

    log_debug(nnti_debug_level, "connection NOT found");
//    print_qpn_map();

    return(NULL);
}
static NNTI_peer_t *get_peer_by_url(const char *url)
{
    ib_connection *conn = NULL;

    log_debug(nnti_debug_level, "looking for url=%s", url);

    conn_by_peer_iter_t i;
    nthread_lock(&nnti_conn_peer_lock);
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p",
                (uint64_t)i->first.addr, (uint64_t)i->first.port, i->second);
        if (strcmp(i->second->peer.url, url) == 0) {
            conn=i->second;
            break;
        }
    }
    nthread_unlock(&nnti_conn_peer_lock);

    if (conn != NULL) {
        log_debug(nnti_debug_level, "peer found");
        NNTI_peer_t *peer=(NNTI_peer_t*)malloc(sizeof(NNTI_peer_t));
        *peer=conn->peer;
        return peer;
    }

    return(NULL);
}
static ib_connection *del_conn_peer(const NNTI_peer_t *peer)
{
    ib_connection  *conn=NULL;
    addrport_key    key;

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_peer(logger_get_file(), "peer",
                "del_conn_peer", peer);
    }

    memset(&key, 0, sizeof(addrport_key));
    key.addr=peer->peer.NNTI_remote_process_t_u.ib.addr;
    key.port=peer->peer.NNTI_remote_process_t_u.ib.port;

    nthread_lock(&nnti_conn_peer_lock);
    if (connections_by_peer.find(key) != connections_by_peer.end()) {
        conn = connections_by_peer[key];
    }

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        connections_by_peer.erase(key);
        del_conn_qpn(conn->req_qp.qpn);
        del_conn_qpn(conn->data_qp.qpn);
    } else {
        log_debug(nnti_debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_peer_lock);

    return(conn);
}
static ib_connection *del_conn_qpn(const NNTI_qp_num qpn)
{
    ib_connection  *conn=NULL;

    nthread_lock(&nnti_conn_qpn_lock);
    if (connections_by_qpn.find(qpn) != connections_by_qpn.end()) {
        conn = connections_by_qpn[qpn];
    }

    if (conn != NULL) {
        log_debug(nnti_debug_level, "connection found");
        connections_by_qpn.erase(qpn);
    } else {
        log_debug(nnti_debug_level, "connection NOT found");
    }
    nthread_unlock(&nnti_conn_qpn_lock);

    return(conn);
}
//static void print_qpn_map()
//{
//    ib_connection     *conn=NULL;
//    conn_by_qpn_iter_t i;
//    for (i=connections_by_qpn.begin(); i != connections_by_qpn.end(); i++) {
//        conn=i->second;
//        log_debug(nnti_debug_level, "qpn_map key=%llu conn=%p (name=%s, addr=%llu, port=%llu)",
//                i->first, conn, conn->peer_name, (uint64_t)conn->peer_addr, (uint64_t)conn->peer_port);
//    }
//}
//static void print_peer_map()
//{
//    ib_connection      *conn=NULL;
//    conn_by_peer_iter_t i;
//    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
//        addrport_key key=i->first;
//        conn=i->second;
//        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p (name=%s, addr=%llu, port=%llu)",
//                (uint64_t)key.addr, (uint64_t)key.port, conn, conn->peer_name, (uint64_t)conn->peer_addr, (uint64_t)conn->peer_port);
//    }
//}

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);

    nthread_lock(&nnti_buf_bufhash_lock);
    assert(buffers_by_bufhash.find(h) == buffers_by_bufhash.end());
    buffers_by_bufhash[h] = buf;
    nthread_unlock(&nnti_buf_bufhash_lock);

    log_debug(nnti_debug_level, "bufhash buffer added (buf=%p bufhash=%x)", buf, (uint64_t)h);

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
    uint32_t h=hash6432shift((uint64_t)buf->buffer_segments.NNTI_remote_addr_array_t_val[0].NNTI_remote_addr_t_u.ib.buf);
    log_level debug_level = nnti_debug_level;

    nthread_lock(&nnti_buf_bufhash_lock);
    if (buffers_by_bufhash.find(h) != buffers_by_bufhash.end()) {
        buf = buffers_by_bufhash[h];
    }

    if (buf != NULL) {
        log_debug(debug_level, "buffer found");
        buffers_by_bufhash.erase(h);
    } else {
        log_debug(debug_level, "buffer NOT found");
    }
    nthread_unlock(&nnti_buf_bufhash_lock);

    return(buf);
}
//static void print_bufhash_map()
//{
//    if (!logging_debug(nnti_debug_level)) {
//        return;
//    }
//
//    if (buffers_by_bufhash.empty()) {
//        log_debug(nnti_debug_level, "bufhash_map is empty");
//        return;
//    }
//
//    buf_by_bufhash_iter_t i;
//    for (i=buffers_by_bufhash.begin(); i != buffers_by_bufhash.end(); i++) {
//        log_debug(nnti_debug_level, "bufhash_map key=%x buf=%p", i->first, i->second);
//    }
//}

static NNTI_result_t wr_pool_register(
        ib_work_request *ib_wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    uint32_t len;

    log_debug(nnti_debug_level, "enter");

    len = sizeof(ib_wr->ack);

    if (config.use_memset) {
        trios_start_timer(callTime);
        memset(&ib_wr->ack, 0, len);
        trios_stop_timer("memset", callTime);
    }
    if (config.use_mlock) {
        trios_start_timer(callTime);
        mlock(&ib_wr->ack, len);
        trios_stop_timer("mlock", callTime);
    }

    trios_start_timer(callTime);
    mr = ibv_reg_mr_wrapper(
            transport_global_data.pd,
            &ib_wr->ack,
            len,
            (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));
    if (!mr) {
        log_error(nnti_debug_level, "failed to register memory region");
        perror("errno");
        return NNTI_EIO;
    }
    trios_stop_timer("register", callTime);

    ib_wr->ack_mr=mr;

    log_debug(nnti_debug_level, "exit (mr==%p, addr %p, length %lu, lkey %x, rkey %x)...", mr, mr->addr, mr->length, mr->lkey, mr->rkey);

    return (rc);
}
static NNTI_result_t wr_pool_deregister(
        ib_work_request *ib_wr)
{
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    if (ib_wr->ack_mr!=NULL) {
        trios_start_timer(callTime);
        ibv_rc=ibv_dereg_mr_wrapper(ib_wr->ack_mr);
        if (ibv_rc != 0) {
            log_error(nnti_debug_level, "deregistering the ACK buffer failed");
        }
        trios_stop_timer("deregister", callTime);
        ib_wr->ack_mr=NULL;
    } else {
        log_debug(nnti_debug_level, "exit ib_wr(%p) - not registered", ib_wr);
        return(NNTI_OK);
    }

    if (config.use_mlock) {
        trios_start_timer(callTime);
        munlock(ib_wr->ack_mr->addr, ib_wr->ack_mr->length);
        trios_stop_timer("munlock", callTime);
    }

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_init(uint32_t pool_size)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    for (i=0;i<pool_size;i++) {
        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        assert(ib_wr);
        log_debug(nnti_debug_level, "allocated ib_wr (ib_wr=%p)", ib_wr);
        nthread_lock_init(&ib_wr->lock);
        rc=wr_pool_register(ib_wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register target work request: rc=%d", rc);
            goto cleanup;
        }
        wr_pool_rdma_push(ib_wr);

        ib_wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        assert(ib_wr);
        log_debug(nnti_debug_level, "allocated ib_wr (ib_wr=%p)", ib_wr);
        nthread_lock_init(&ib_wr->lock);
        wr_pool_sendrecv_push(ib_wr);
    }

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(rc);
}
static ib_work_request *wr_pool_rdma_pop(void)
{
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    if (!rdma_wr_pool.empty()) {
        ib_wr=rdma_wr_pool.front();
        rdma_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(ib_wr);
}
static ib_work_request *wr_pool_sendrecv_pop(void)
{
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    if (!sendrecv_wr_pool.empty()) {
        ib_wr=sendrecv_wr_pool.front();
        sendrecv_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(ib_wr);
}
static void wr_pool_rdma_push(ib_work_request *ib_wr)
{
    log_debug(nnti_debug_level, "enter");

    ib_wr->last_op=0;
    ib_wr->state  =NNTI_IB_WR_STATE_RESET;

    nthread_lock(&nnti_wr_pool_lock);
    rdma_wr_pool.push_front(ib_wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static void wr_pool_sendrecv_push(ib_work_request *ib_wr)
{
    log_debug(nnti_debug_level, "enter");

    ib_wr->last_op=0;
    ib_wr->state  =NNTI_IB_WR_STATE_RESET;

    nthread_lock(&nnti_wr_pool_lock);
    sendrecv_wr_pool.push_front(ib_wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static NNTI_result_t wr_pool_fini(void)
{
    NNTI_result_t  rc=NNTI_OK;
    ib_work_request *ib_wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    while (!rdma_wr_pool.empty()) {
        ib_wr=rdma_wr_pool.front();
        rdma_wr_pool.pop_front();
        assert(ib_wr);
        rc=wr_pool_deregister(ib_wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister target work request: rc=%d", rc);
            goto cleanup;
        }
        log_debug(nnti_debug_level, "freeing ib_wr=%p", ib_wr);
        free(ib_wr);
    }
    while (!sendrecv_wr_pool.empty()) {
        ib_wr=sendrecv_wr_pool.front();
        sendrecv_wr_pool.pop_front();
        assert(ib_wr);
        log_debug(nnti_debug_level, "freeing ib_wr=%p", ib_wr);
        free(ib_wr);
    }

cleanup:
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}

static void close_all_conn(void)
{
    log_level debug_level = nnti_debug_level;

    log_debug(nnti_debug_level, "enter (%d qpn connections, %d peer connections)",
            connections_by_qpn.size(), connections_by_peer.size());

    nthread_lock(&nnti_conn_qpn_lock);
    conn_by_qpn_iter_t qpn_iter = connections_by_qpn.begin();
    while (qpn_iter != connections_by_qpn.end()) {
        log_debug(debug_level, "close connection (qpn=%llu)", qpn_iter->first);
        close_connection(qpn_iter->second);

        connections_by_qpn.erase(qpn_iter++);
    }

    nthread_unlock(&nnti_conn_qpn_lock);

    nthread_lock(&nnti_conn_peer_lock);
    conn_by_peer_iter_t peer_iter = connections_by_peer.begin();
    while (peer_iter != connections_by_peer.end()) {
        log_debug(debug_level, "close connection (peer.addr=%llu)", peer_iter->first.addr);
//        close_connection(peer_iter->second);

        connections_by_peer.erase(peer_iter++);
    }
    nthread_unlock(&nnti_conn_peer_lock);

    log_debug(debug_level, "exit (%d qpn connections, %d peer connections)",
            connections_by_qpn.size(), connections_by_peer.size());

    return;
}

/**
 * @brief initialize
 */
static NNTI_result_t init_connection(
        ib_connection **conn,
        const int sock,
        const int is_server)
{
    int rc=0; /* return code */

    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "initializing ib connection");

    (*conn)->disconnect_requested = FALSE;

    trios_start_timer(callTime);
    if (is_server) {
        rc = new_server_connection(*conn, sock);
    } else {
        rc = new_client_connection(*conn, sock);
    }
    if (rc) {
        close_connection(*conn);
        goto out;
    }
    trios_stop_timer("new connection", callTime);

    print_ib_conn(*conn);

out:
    return((NNTI_result_t)rc);
}

/*
 * At an explicit BYE message, or at finalize time, shut down a connection.
 * If descriptors are posted, defer and clean up the connection structures
 * later.
 */
static void close_connection(ib_connection *c)
{
    int rc;
//    int i;

    if (c==NULL) return;

    log_debug(nnti_debug_level, "enter");

    if (c->state==DISCONNECTED) {
        log_debug(nnti_debug_level, "conn(%p) is already closed", c);
        return;
    }

    print_ib_conn(c);

    transition_connection_to_error(c);

    if (c->peer_name) free(c->peer_name);
    if (c->req_qp.qp) {
        rc=ibv_destroy_qp_wrapper(c->req_qp.qp);
        if (rc < 0)
            log_error(nnti_debug_level, "failed to destroy QP");
    }
    if (c->data_qp.qp) {
        rc=ibv_destroy_qp_wrapper(c->data_qp.qp);
        if (rc < 0)
            log_error(nnti_debug_level, "failed to destroy QP");
    }
    c->state=DISCONNECTED;

    log_debug(nnti_debug_level, "exit");
}

static NNTI_result_t check_for_waiting_connection()
{
    NNTI_result_t rc = NNTI_OK;

    struct sockaddr_in ssin;
    socklen_t len;
    int s;
    NNTI_peer_t peer;   // this doesn't need to be a pointer
    ib_connection *conn = NULL;

    len = sizeof(ssin);
    s = accept(transport_global_data.listen_sock, (struct sockaddr *) &ssin, &len);
    if (s < 0) {
        if ((errno == EAGAIN) || (errno == EWOULDBLOCK)) {
            log_debug(nnti_debug_level, "no connections waiting to be accepted: %s", strerror(errno));
            rc = NNTI_EWOULDBLOCK;
            goto cleanup;
        } else {
            log_error(nnti_debug_level, "failed to accept tcp socket connection: %s", strerror(errno));
            rc = NNTI_EIO;
            goto cleanup;
        }
    } else {
        static int connection_count = 0;
        char         *peer_hostname = strdup(inet_ntoa(ssin.sin_addr));
//        NNTI_ip_addr  peer_addr  = ssin.sin_addr.s_addr;
        NNTI_tcp_port peer_port  = ntohs(ssin.sin_port);

        conn = (ib_connection *)calloc(1, sizeof(ib_connection));    // TODO: FIND OUT WHERE THIS IS FREED
        log_debug(nnti_debug_level, "calloc returned conn=%p.", conn);
        if (conn == NULL) {
            log_error(nnti_debug_level, "calloc returned NULL.  out of memory?: %s", strerror(errno));
            rc=NNTI_ENOMEM;
            goto cleanup;
        }

//        nthread_lock(&nnti_ib_lock);
        rc=init_connection(&conn, s, 1);
        if (rc!=NNTI_OK) {
            goto cleanup;
        }
        create_peer(
                &peer,
                conn->peer_name,
                conn->peer_addr,
                conn->peer_port);

        conn->peer=peer;

        insert_conn_qpn(conn->req_qp.qpn, conn);
        insert_conn_qpn(conn->data_qp.qpn, conn);
        insert_conn_peer(&peer, conn);

        log_debug(nnti_debug_level, "Allocating new connection count=%d", ++connection_count);

        transition_connection_to_ready(s, conn);
//        nthread_unlock(&nnti_ib_lock);

        log_debug(nnti_debug_level, "accepted new connection from %s:%u", peer_hostname, peer_port);
        if (peer_hostname) free(peer_hostname);   // fixed to avoid lost bytes reported by valgrind

        if (close(s) < 0) {
            log_error(nnti_debug_level, "failed to close new tcp socket");
        }

        if (logging_debug(nnti_debug_level)) {
            fprint_NNTI_peer(logger_get_file(), "peer",
                    "end of check_listen_socket_for_new_connections", &peer);
        }
    }

cleanup:
    return rc;
}

/**
 * Check for new connections.  The listening socket is left nonblocking
 * so this test can be quick; but accept is not really that quick compared
 * to polling an IB interface, for instance.  Returns >0 if an accept worked.
 */
static NNTI_result_t check_listen_socket_for_new_connections()
{
    bool done=false;
    trios_declare_timer(callTime);

    trios_start_timer(callTime);
    while(!done) {
        if (check_for_waiting_connection() != NNTI_OK) {
            done=true;
        }
    }
    trios_stop_timer("check_for_waiting_connection", callTime);

    return(NNTI_OK);
}


static struct ibv_device *get_ib_device(void)
{
    struct ibv_device *dev;
    struct ibv_device **dev_list;
    int dev_count=0;

    dev_list = ibv_get_device_list_wrapper(&dev_count);
    if (dev_count == 0)
        return NULL;
    if (dev_count > 1) {
        log_debug(nnti_debug_level, "found %d devices, defaulting the dev_list[0] (%p)", dev_count, dev_list[0]);
    }
    dev = dev_list[0];
    ibv_free_device_list_wrapper(dev_list);

    return dev;
}

static void get_qp_state(struct ibv_qp *qp)
{
    struct ibv_qp_attr attr;
    struct ibv_qp_init_attr init_attr;

    if (ibv_query_qp(qp, &attr,
              IBV_QP_STATE, &init_attr)) {
        log_error(nnti_debug_level, "Failed to query QP state\n");
        return;
    }
    log_debug(nnti_debug_level, "qp(%p)->state=%d", qp, attr.qp_state);
}

static void print_all_qp_state(void)
{
    if (logging_debug(nnti_debug_level)) {
        ib_connection     *conn=NULL;
        conn_by_qpn_iter_t i=connections_by_qpn.begin();
        if (i != connections_by_qpn.end()) {
            conn=i->second;
            get_qp_state(conn->req_qp.qp);
        }
    }
}

static void print_wc(const struct ibv_wc *wc, bool force)
{
    if (force) {
        log_debug(LOG_ALL, "wc=%p, wc.opcode=%d, wc.flags=%d, wc.status=%d (%s), wc.wr_id=%lx, wc.vendor_err=%u, wc.byte_len=%u, wc.qp_num=%u, wc.imm_data=%x, wc.src_qp=%u",
            wc,
            wc->opcode,
            wc->wc_flags,
            wc->status,
            ibv_wc_status_str_wrapper(wc->status),
            wc->wr_id,
            wc->vendor_err,
            wc->byte_len,
            wc->qp_num,
            wc->imm_data,
            wc->src_qp);
    } else if ((wc->status != IBV_WC_SUCCESS) && (wc->status != IBV_WC_RNR_RETRY_EXC_ERR)) {
        log_error(nnti_debug_level, "wc=%p, wc.opcode=%d, wc.flags=%d, wc.status=%d (%s), wc.wr_id=%lx, wc.vendor_err=%u, wc.byte_len=%u, wc.qp_num=%u, wc.imm_data=%x, wc.src_qp=%u",
            wc,
            wc->opcode,
            wc->wc_flags,
            wc->status,
            ibv_wc_status_str_wrapper(wc->status),
            wc->wr_id,
            wc->vendor_err,
            wc->byte_len,
            wc->qp_num,
            wc->imm_data,
            wc->src_qp);
    } else {
        log_debug(nnti_debug_level, "wc=%p, wc.opcode=%d, wc.flags=%d, wc.status=%d (%s), wc.wr_id=%lx, wc.vendor_err=%u, wc.byte_len=%u, wc.qp_num=%u, wc.imm_data=%x, wc.src_qp=%u",
            wc,
            wc->opcode,
            wc->wc_flags,
            wc->status,
            ibv_wc_status_str_wrapper(wc->status),
            wc->wr_id,
            wc->vendor_err,
            wc->byte_len,
            wc->qp_num,
            wc->imm_data,
            wc->src_qp);
    }
}

static NNTI_result_t process_comp_channel_event(
        struct ibv_comp_channel *comp_channel,
        struct ibv_cq           *cq)
{
    NNTI_result_t rc=NNTI_OK;

    int retries_left=3;

    struct ibv_cq *ev_cq;
    void          *ev_ctx;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);


    log_debug(nnti_debug_level, "enter");

    trios_start_timer(call_time);
    trios_start_timer(total_time);

try_again:
    if (ibv_get_cq_event_wrapper(comp_channel, &ev_cq, &ev_ctx) == 0) {
        trios_stop_timer("ibv_get_cq_event", call_time);
        log_debug(nnti_debug_level, "got event from comp_channel=%p for cq=%p", comp_channel, ev_cq);
        trios_start_timer(call_time);
        ibv_ack_cq_events_wrapper(ev_cq, 1);
        trios_stop_timer("ibv_ack_cq_events", call_time);
        log_debug(nnti_debug_level, "ACKed event on cq=%p", ev_cq);
        rc = NNTI_OK;
    } else {
        if (errno == EAGAIN) {
            if (retries_left > 0) {
                trios_stop_timer("ibv_get_cq_event...EAGAIN...retry", call_time);
                retries_left--;
                trios_start_timer(call_time);
                goto try_again;
            } else {
                trios_stop_timer("ibv_get_cq_event...EAGAIN...retries exhausted", call_time);
                rc = NNTI_EAGAIN;
                goto cleanup;
            }
        }
        trios_stop_timer("ibv_get_cq_event...failed", call_time);
        log_error(nnti_debug_level, "ibv_get_cq_event failed (ev_cq==%p): %s",
                ev_cq, strerror(errno));
        rc = NNTI_EIO;
        goto cleanup;
    }

cleanup:
    trios_stop_timer("process_comp_channel_event", total_time);

    log_debug(nnti_debug_level, "exit");

    return(rc);
}

#define CONNECTION_SOCKET_INDEX 0
#define DATA_CQ_SOCKET_INDEX    1
#define REQ_CQ_SOCKET_INDEX     2
#define INTERRUPT_PIPE_INDEX    3
#define FD_COUNT 4
static NNTI_result_t poll_all(
        int timeout)
{
    NNTI_result_t rc=NNTI_OK;
    int poll_rc=0;
    struct pollfd my_pollfd[FD_COUNT];

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);


    log_debug(nnti_debug_level, "enter (timeout=%d)", timeout);

    trios_start_timer(total_time);

    /*
     * poll the channel until it has an event and sleep ms_timeout
     * milliseconds between any iteration
     */
    my_pollfd[CONNECTION_SOCKET_INDEX].fd      = transport_global_data.listen_sock;
    my_pollfd[CONNECTION_SOCKET_INDEX].events  = POLLIN;
    my_pollfd[CONNECTION_SOCKET_INDEX].revents = 0;
    my_pollfd[DATA_CQ_SOCKET_INDEX].fd         = transport_global_data.data_comp_channel->fd;
    my_pollfd[DATA_CQ_SOCKET_INDEX].events     = POLLIN;
    my_pollfd[DATA_CQ_SOCKET_INDEX].revents    = 0;
    my_pollfd[REQ_CQ_SOCKET_INDEX].fd          = transport_global_data.req_comp_channel->fd;
    my_pollfd[REQ_CQ_SOCKET_INDEX].events      = POLLIN;
    my_pollfd[REQ_CQ_SOCKET_INDEX].revents     = 0;
    my_pollfd[INTERRUPT_PIPE_INDEX].fd         = transport_global_data.interrupt_pipe[0];
    my_pollfd[INTERRUPT_PIPE_INDEX].events     = POLLIN;
    my_pollfd[INTERRUPT_PIPE_INDEX].revents    = 0;
    log_debug(nnti_debug_level, "polling with timeout==%d", timeout);

    // Test for errno==EINTR to deal with timing interrupts from HPCToolkit
    do {
        trios_start_timer(call_time);
        poll_rc = poll(&my_pollfd[0], FD_COUNT, timeout);
        trios_stop_timer("poll", call_time);
    } while ((poll_rc < 0) && (errno == EINTR));

    if (poll_rc == 0) {
        log_debug(nnti_debug_level, "poll() timed out: poll_rc=%d", poll_rc);
        rc = NNTI_ETIMEDOUT;
        goto cleanup;
    } else if (poll_rc < 0) {
        if (errno == EINTR) {
            log_error(nnti_debug_level, "poll() interrupted by signal: poll_rc=%d (%s)", poll_rc, strerror(errno));
            rc = NNTI_EINTR;
        } else if (errno == ENOMEM) {
            log_error(nnti_debug_level, "poll() out of memory: poll_rc=%d (%s)", poll_rc, strerror(errno));
            rc = NNTI_ENOMEM;
        } else {
            log_error(nnti_debug_level, "poll() invalid args: poll_rc=%d (%s)", poll_rc, strerror(errno));
            rc = NNTI_EINVAL;
        }
        goto cleanup;
    } else {
        log_debug(nnti_debug_level, "polled on %d file descriptor(s).  events occurred on %d file descriptor(s).", FD_COUNT, poll_rc);
        log_debug(nnti_debug_level, "poll success: poll_rc=%d ; my_pollfd[CONNECTION_SOCKET_INDEX].revents=%d",
                poll_rc, my_pollfd[CONNECTION_SOCKET_INDEX].revents);
        log_debug(nnti_debug_level, "poll success: poll_rc=%d ; my_pollfd[DATA_CQ_SOCKET_INDEX].revents=%d",
                poll_rc, my_pollfd[DATA_CQ_SOCKET_INDEX].revents);
        log_debug(nnti_debug_level, "poll success: poll_rc=%d ; my_pollfd[REQ_CQ_SOCKET_INDEX].revents=%d",
                poll_rc, my_pollfd[REQ_CQ_SOCKET_INDEX].revents);
        log_debug(nnti_debug_level, "poll success: poll_rc=%d ; my_pollfd[INTERRUPT_PIPE_INDEX].revents=%d",
                poll_rc, my_pollfd[INTERRUPT_PIPE_INDEX].revents);
    }

    if (my_pollfd[CONNECTION_SOCKET_INDEX].revents == POLLIN) {
        check_listen_socket_for_new_connections();
    }
    if (my_pollfd[DATA_CQ_SOCKET_INDEX].revents == POLLIN) {
        process_comp_channel_event(transport_global_data.data_comp_channel, transport_global_data.data_cq);
    }
    if (my_pollfd[REQ_CQ_SOCKET_INDEX].revents == POLLIN) {
        process_comp_channel_event(transport_global_data.req_comp_channel, transport_global_data.req_cq);
    }
    if (my_pollfd[INTERRUPT_PIPE_INDEX].revents == POLLIN) {
        log_debug(nnti_debug_level, "poll() interrupted by NNTI_ib_interrupt");
        // read all bytes from the pipe
        ssize_t bytes_read=0;
        uint32_t dummy=0;
        do {
            bytes_read=read(transport_global_data.interrupt_pipe[0], &dummy, 4);
            if (dummy != 0xAAAAAAAA) {
                log_warn(nnti_debug_level, "interrupt byte is %X, should be 0xAAAAAAAA", dummy);
            }
            log_debug(nnti_debug_level, "bytes_read==%lu", (uint64_t)bytes_read);
        } while(bytes_read > 0);
        rc = NNTI_EINTR;
    }
    if (ibv_req_notify_cq_wrapper(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
        rc = NNTI_EIO;
    }
    if (ibv_req_notify_cq_wrapper(transport_global_data.req_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
        rc = NNTI_EIO;
    }

cleanup:
    trios_stop_timer("poll_all", total_time);

    log_debug(nnti_debug_level, "exit");
    return(rc);
}


#define CQ_COUNT 2
static NNTI_result_t progress(
        int                   timeout,
        NNTI_work_request_t **wr_list,
        const uint32_t        wr_count)
{
    int           rc=0;
    NNTI_result_t nnti_rc=NNTI_OK;
    int           ibv_rc=0;

    uint32_t which=0;

    struct ibv_wc            wc;
    struct ibv_comp_channel *comp_channel;
    struct ibv_cq           *cq_list[CQ_COUNT];
    static int               cq_index=0;

    long entry_time  =trios_get_time_ms();
    long elapsed_time=0;

    log_level debug_level  =nnti_debug_level;
    log_level old_log_level=logger_get_default_level();

    static bool in_progress  =false;   // if true, another thread is already making progress.
    bool        made_progress=false;
    int8_t      wr_complete  =FALSE;

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);


    log_debug(debug_level, "enter (timeout=%d)", timeout);

    trios_start_timer(total_time);

    cq_list[0]=transport_global_data.data_cq;
    cq_list[1]=transport_global_data.req_cq;

    /*
     * Only one thread is allowed to make progress at a time.  All others
     * wait for the progress maker to finish, then everyone returns at once.
     */
    nthread_lock(&nnti_progress_lock);

    wr_complete = is_any_wr_complete(wr_list, wr_count, &which);

    if (!in_progress) {
        log_debug(debug_level, "making progress");
        // no other thread is making progress.  we'll do it.
        in_progress=true;
        log_debug(debug_level, "set in_progress=true");
        nthread_unlock(&nnti_progress_lock);
    } else {
        if (wr_complete == TRUE) {
            nthread_unlock(&nnti_progress_lock);
            goto cleanup;
        }

        // another thread is making progress.  we'll wait until they are done.
        rc=0;
        elapsed_time=0;
        if (in_progress) {
            if (timeout > 0) {
                log_debug(debug_level, "waiting for progress with timeout=%d", timeout-elapsed_time);
                // wait for progress or until timeout
                rc=nthread_timedwait(&nnti_progress_cond, &nnti_progress_lock, timeout-elapsed_time);
            } else if (timeout < 0) {
                log_debug(debug_level, "waiting infinitely for progress");
                // infinite wait for progress
                rc=nthread_wait(&nnti_progress_cond, &nnti_progress_lock);
            } else {
                log_debug(debug_level, "waiting for progress with timeout=0.  immediate timeout.");
                // timeout == 0 and we are not the progress maker.  report a timeout.
                rc=ETIMEDOUT;
            }
            elapsed_time = (trios_get_time_ms() - entry_time);
            log_debug(debug_level, "rc=%d, elapsed_time=%d", rc, elapsed_time);
        }
        nthread_unlock(&nnti_progress_lock);
        if (rc == ETIMEDOUT) {
            log_debug(debug_level, "timed out waiting for progress");
            nnti_rc = NNTI_ETIMEDOUT;
        } else if (rc == 0) {
            log_debug(debug_level, "someone made progress");
            nnti_rc=NNTI_OK;
        }
        goto cleanup;
    }

    while (!made_progress)   {

        if (trios_exit_now()) {
            log_debug(debug_level, "caught abort signal");
            nnti_rc=NNTI_ECANCELED;
            break;
        }

        check_listen_socket_for_new_connections();

        for (int i=0;i<CQ_COUNT;i++) {

            struct ibv_cq *cq=cq_list[i];

            memset(&wc, 0, sizeof(struct ibv_wc));

            log_debug(debug_level, "polling for 1 work completion on cq=%p", cq);
            trios_start_timer(call_time);
            ibv_rc = ibv_poll_cq_wrapper(cq, 1, &wc);
            trios_stop_timer("progress - ibv_poll_cq", call_time);

            log_debug(debug_level, "ibv_poll_cq(cq=%p) rc==%d", cq, ibv_rc);

            if (ibv_rc < 0) {
                // an error occurred
                log_debug(debug_level, "ibv_poll_cq failed: %d", ibv_rc);
            } else if (ibv_rc == 0) {
                // no work completions in the queue
            } else  if (ibv_rc > 0) {
                // got a work completion
                log_debug(debug_level, "got wc from cq=%p", cq);
                log_debug(debug_level, "polling status is %s", ibv_wc_status_str(wc.status));

                print_wc(&wc, false);
                if ((wc.status == IBV_WC_RNR_RETRY_EXC_ERR) ||
                    (wc.status == IBV_WC_RETRY_EXC_ERR)) {
                    nnti_rc=NNTI_EDROPPED;

                    ib_work_request *ib_wr=decode_work_request(&wc);
                    int min_rnr_timer=1;  /* means 0.01ms delay before sending RNR NAK */
                    int ack_timeout  =17; /* time to wait for ACK/NAK before retransmitting.  4.096us * 2^17 == 0.536ss */
                    int retry_count;
                    if (config.drop_if_full_queue) {
                        retry_count=1; /* number of retries if no answer on primary path or if remote sends RNR NAK */
                    } else {
                        retry_count=7; /* number of retries if no answer on primary path or if remote sends RNR NAK.  7 has special meaning of infinite retries. */
                    }
                    transition_qp_from_error_to_ready(
                            ib_wr->conn->req_qp.qp,
                            ib_wr->conn->peer_req_qpn,
                            ib_wr->conn->peer_lid,
                            min_rnr_timer,
                            ack_timeout,
                            retry_count);

                } else if (wc.status != IBV_WC_SUCCESS) {
                    log_error(debug_level, "Failed status %s (%d) for wr_id %lx",
                            ibv_wc_status_str(wc.status),
                            wc.status, wc.wr_id);
                    nnti_rc=NNTI_EIO;
                }

                trios_start_timer(call_time);
                ib_work_request *ib_wr=decode_work_request(&wc);
                trios_stop_timer("progress - decode_work_request", call_time);
                trios_start_timer(call_time);
                nthread_lock(&ib_wr->lock);
                process_event(ib_wr, &wc);
                nthread_unlock(&ib_wr->lock);
                trios_stop_timer("progress - process_event", call_time);

                made_progress=true;
            }
        }

        if ((!made_progress) && (wr_complete == FALSE)) {
            trios_start_timer(call_time);
            rc = poll_all(/*100*/ timeout-elapsed_time);
            trios_stop_timer("progress - poll_all", call_time);

            elapsed_time = (trios_get_time_ms() - entry_time);

            /* case 1: success */
            if (rc == NNTI_OK) {
                logger_set_default_level(old_log_level);
                nnti_rc = NNTI_OK;
                continue;
            }
            /* case 2: timed out */
            else if (rc==NNTI_ETIMEDOUT) {
                /* if the caller asked for a legitimate timeout, we need to exit */
                if (((timeout >= 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                    logger_set_default_level(old_log_level);
                    log_debug(debug_level, "poll_all timed out");
                    nnti_rc = NNTI_ETIMEDOUT;
                    break;
                }
                /* continue if the timeout has not expired */
                log_debug(debug_level, "poll_all timedout... retrying");

//                log_debug(debug_level, "***** disable debug logging.  will enable after polling success. *****");
//                logger_set_default_level(LOG_OFF);
            }
            /* case 3: poll was interrupted.  could be a signal or NNTI_interrupt(). */
            else if (rc==NNTI_EINTR) {
                logger_set_default_level(old_log_level);
                nnti_rc = NNTI_EINTR;
                break;
            }
            /* case 4: poll out of memory */
            else if (rc==NNTI_ENOMEM) {
                logger_set_default_level(old_log_level);
                nnti_rc = NNTI_ENOMEM;
                break;
            }
            /* case 5: poll got invalid arguments */
            else if (rc==NNTI_EINVAL) {
                logger_set_default_level(old_log_level);
                nnti_rc = NNTI_EINVAL;
                break;
            }
            /* case 6: failure */
            else {
                logger_set_default_level(old_log_level);
                log_error(debug_level, "poll_all failed: %s", strerror(errno));
                nnti_rc = NNTI_EIO;
                break;
            }
        }
    }

unlock:
    nthread_lock(&nnti_progress_lock);
    in_progress=false;
    log_debug(debug_level, "set in_progress=false");
    nthread_broadcast(&nnti_progress_cond);
    log_debug(debug_level, "broadcasted on nnti_progress_cond");
    nthread_unlock(&nnti_progress_lock);

cleanup:
    trios_stop_timer("progress", total_time);

    log_debug(debug_level, "exit");

    return(nnti_rc);
}


static void print_ib_conn(ib_connection *c)
{
//    int i=0;
    log_level debug_level=nnti_debug_level;

    log_debug(debug_level, "c->peer_name       =%s", c->peer_name);
    log_debug(debug_level, "c->peer_addr       =%u", c->peer_addr);
    log_debug(debug_level, "c->peer_port       =%u", (uint32_t)c->peer_port);
    log_debug(debug_level, "c->peer_lid        =%llu", (uint64_t)c->peer_lid);
    log_debug(debug_level, "c->peer_req_qpn    =%llu", (uint64_t)c->peer_req_qpn);

    log_debug(debug_level, "c->req_comp_channel=%p", transport_global_data.req_comp_channel);
    log_debug(debug_level, "c->req_cq          =%p", transport_global_data.req_cq);
    log_debug(debug_level, "c->req_srq         =%p", transport_global_data.req_srq);

    log_debug(debug_level, "c->req_qp.qp          =%p", c->req_qp.qp);
    log_debug(debug_level, "c->req_qp.qpn         =%llu", (uint64_t)c->req_qp.qpn);

    log_debug(debug_level, "c->data_qp.qp          =%p",     c->data_qp.qp);
    log_debug(debug_level, "c->data_qp.qpn         =%llu",   (uint64_t)c->data_qp.qpn);
    log_debug(debug_level, "c->data_qp.peer_qpn    =%llu",   (uint64_t)c->data_qp.peer_qpn);

    log_debug(debug_level, "c->state           =%d", c->state);

    log_debug(debug_level, "c->disconnect_requested=%d", c->disconnect_requested);

}

static void config_init(nnti_ib_config *c)
{
    c->min_atomics_vars    = 512;
    c->use_wr_pool         = false;
    c->use_rdma_target_ack = false;
    c->use_mlock           = true;
    c->use_memset          = true;
    c->drop_if_full_queue  = false;
}

static void config_get_from_env(nnti_ib_config *c)
{
    char *env_str=NULL;

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
    if ((env_str=getenv("TRIOS_NNTI_USE_WR_POOL")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_wr_pool to TRUE");
            c->use_wr_pool=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_wr_pool to FALSE");
            c->use_wr_pool=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_WR_POOL is undefined.  using c->use_wr_pool default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_RDMA_TARGET_ACK")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(nnti_debug_level, "setting c->use_rdma_target_ack to TRUE");
            c->use_rdma_target_ack=true;
        } else {
            log_debug(nnti_debug_level, "setting c->use_rdma_target_ack to FALSE");
            c->use_rdma_target_ack=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_RDMA_TARGET_ACK is undefined.  using c->use_rdma_target_ack default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_MLOCK")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(LOG_ALL, "setting c->use_mlock to TRUE");
            c->use_mlock=true;
        } else {
            log_debug(LOG_ALL, "setting c->use_mlock to FALSE");
            c->use_mlock=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_MLOCK is undefined.  using c->mlock default");
    }
    if ((env_str=getenv("TRIOS_NNTI_USE_MEMSET")) != NULL) {
        if ((!strcasecmp(env_str, "TRUE")) ||
            (!strcmp(env_str, "1"))) {
            log_debug(LOG_ALL, "setting c->use_memset to TRUE");
            c->use_memset=true;
        } else {
            log_debug(LOG_ALL, "setting c->use_memset to FALSE");
            c->use_memset=false;
        }
    } else {
        log_debug(nnti_debug_level, "TRIOS_NNTI_USE_MEMSET is undefined.  using c->memset default");
    }
}

//static void print_wr(ib_work_request *ib_wr)
//{
//    log_debug(nnti_debug_level, "ib_wr (op=%llu ; offset=%llu ; length=%llu)",
//            (uint64_t)ib_wr->last_op,
//            ib_wr->offset,
//            ib_wr->length);
//}
//
//static void print_ack_buf(ib_rdma_ack *ack)
//{
//    log_debug(nnti_debug_level, "ack (op=%llu ; offset=%llu ; length=%llu)",
//            (uint64_t)ack->op,
//            ack->offset,
//            ack->length);
//}
//
//static void print_xfer_buf(void *buf, uint32_t size)
//{
//    struct data_t {
//            uint32_t int_val;
//            float float_val;
//            double double_val;
//    };
//
//    struct  data_array_t {
//            u_int data_array_t_len;
//            struct data_t *data_array_t_val;
//    };
//
////    struct data_array_t *da=(struct data_array_t *)buf;
//    const struct data_t *array = (struct data_t *)buf;
//    const int len = size/sizeof(struct data_t);
//    int idx=0;
//
//    for (idx=0;idx<len;idx++) {
//        log_debug(nnti_debug_level, "array[%d].int_val=%u ; array[%d].float_val=%f ; array[%d].double_val=%f",
//                idx, array[idx].int_val,
//                idx, array[idx].float_val,
//                idx, array[idx].double_val);
//    }
//
//}
