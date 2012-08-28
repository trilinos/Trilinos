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

    bool use_wr_pool;

    bool use_rdma_target_ack;

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

typedef enum {
    BUFFER_INIT=0,
    SEND_COMPLETE=1,
    RECV_COMPLETE,
    RDMA_WRITE_INIT,
    RDMA_WRITE_NEED_ACK,
    RDMA_WRITE_COMPLETE,
    RDMA_READ_INIT,
    RDMA_READ_NEED_ACK,
    RDMA_READ_COMPLETE,
    RDMA_TARGET_INIT,
    RDMA_TARGET_NEED_ACK,
    RDMA_TARGET_COMPLETE,
    RDMA_COMPLETE
} ib_op_state_t;

#define SRQ_DEPTH 20480
#define CQ_DEPTH 20480

typedef struct {
    struct ibv_qp           *qp;
    uint32_t                 qpn;
    uint32_t                 peer_qpn;
} conn_qp;
typedef struct {
    NNTI_peer_t       peer;
    char         *peer_name;
    NNTI_ip_addr  peer_addr;
    NNTI_tcp_port peer_port;
    uint16_t      peer_lid;
    uint32_t      peer_req_qpn;

    conn_qp       req_qp;
    conn_qp       data_qp;

    ib_connection_state state;

    int8_t disconnect_requested;
} ib_connection;

typedef struct {
    uint32_t op;
    uint64_t offset;
    uint64_t length;
} ib_rdma_ack;

typedef struct {
    NNTI_buffer_t *reg_buf;
    ib_connection *conn;

    struct ibv_send_wr sq_wr;
    struct ibv_recv_wr rq_wr;
    struct ibv_sge     sge;

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
    ib_op_state_t op_state;
    uint8_t       is_last_op_complete;

} ib_work_request;

typedef std::deque<ib_work_request *>           wr_queue_t;
typedef std::deque<ib_work_request *>::iterator wr_queue_iter_t;

typedef struct {
    ib_buffer_type type;
    struct ibv_mr *mr;
    wr_queue_t     wr_queue;
    uint32_t       ref_count;
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
    uint32_t                 qp_count;

    struct ibv_comp_channel *req_comp_channel;
    struct ibv_cq           *req_cq;
    struct ibv_srq          *req_srq;

    struct ibv_comp_channel *data_comp_channel;
    struct ibv_cq           *data_cq;
    struct ibv_srq          *data_srq;

    int32_t req_srq_count;
    int32_t data_srq_count;

    uint16_t nic_lid;
    int      nic_port;

    int      listen_sock;
    char     listen_name[NNTI_HOSTNAME_LEN];
    uint32_t listen_addr;  /* in NBO */
    uint16_t listen_port;  /* in NBO */

    ib_request_queue_handle req_queue;
} ib_transport_global;




static nthread_lock_t nnti_ib_lock;


static int register_memory(
        ib_memory_handle *hdl,
        void *buf,
        uint64_t len,
        enum ibv_access_flags access);
static int unregister_memory(
        ib_memory_handle *hdl);
static int register_ack(
        ib_work_request *wr);
static int unregister_ack(
        ib_work_request *wr);
static void send_ack (
        ib_work_request *wr);
static NNTI_result_t setup_data_channel(void);
static NNTI_result_t setup_request_channel(void);
static ib_work_request *decode_work_request(
        const struct ibv_wc *wc);
static const NNTI_buffer_t *decode_event_buffer(
        const NNTI_buffer_t *wait_buf,
        const struct ibv_wc *wc);
static int process_event(
        const NNTI_buffer_t  *reg_buf,
        const struct ibv_wc  *wc);
static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         wr_id,
        uint64_t        offset,
        uint64_t        length);
static NNTI_result_t post_ack_recv_work_request(
        NNTI_buffer_t  *reg_buf);
static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        ib_work_request  *wr);
static NNTI_result_t repost_ack_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        ib_work_request  *wr);
static ib_work_request *first_incomplete_wr(
        ib_memory_handle *ib_mem_hdl);
static int8_t is_buf_op_complete(
        const NNTI_buffer_t *reg_buf);
static int8_t is_any_buf_op_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        uint32_t             *which);
static int8_t is_all_buf_ops_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count);
static void create_status(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        NNTI_result_t         nnti_rc,
        NNTI_status_t        *status);
static void create_peer(
        NNTI_peer_t *peer,
        char *name,
        NNTI_ip_addr addr,
        NNTI_tcp_port port);
static void copy_peer(
        NNTI_peer_t *src,
        NNTI_peer_t *dest);
static int init_server_listen_socket(void);
static int start_connection_listener_thread(void);
static NNTI_result_t check_listen_socket_for_new_connections();
static struct ibv_device *get_ib_device(void);
static int tcp_read(int sock, void *incoming, size_t len);
static int tcp_write(int sock, const void *outgoing, size_t len);
static int tcp_exchange(int sock, int is_server, void *incoming, void *outgoing, size_t len);
static void transition_connection_to_ready(
        int sock,
        ib_connection *conn);
static void transition_qp_to_ready(
        struct ibv_qp *qp,
        uint32_t peer_qpn,
        int peer_lid);
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
        const struct ibv_wc *wc);
static NNTI_result_t poll_comp_channel(
        struct ibv_comp_channel *comp_channel,
        struct ibv_cq           *cq,
        int timeout);
static void print_ib_conn(ib_connection *c);
static void print_qpn_map(void);
static void print_peer_map(void);
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
static void print_bufhash_map(void);

static NNTI_result_t insert_wr_wrhash(ib_work_request *);
static ib_work_request *get_wr_wrhash(const uint32_t bufhash);
static ib_work_request *del_wr_wrhash(ib_work_request *);
static void print_wrhash_map(void);

static NNTI_result_t wr_pool_init(uint32_t pool_size);
static ib_work_request *wr_pool_rdma_pop(void);
static ib_work_request *wr_pool_sendrecv_pop(void);
static void wr_pool_rdma_push(ib_work_request *wr);
static void wr_pool_sendrecv_push(ib_work_request *wr);
static NNTI_result_t wr_pool_fini(void);

static void close_all_conn(void);

static void config_init(
        nnti_ib_config *c);
static void config_get_from_env(
        nnti_ib_config *c);

//static void print_wr(ib_work_request *wr);
//static void print_xfer_buf(void *buf, uint32_t size);
//static void print_ack_buf(ib_rdma_ack *ack);

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

static std::map<uint32_t, ib_work_request *> wr_by_wrhash;
typedef std::map<uint32_t, ib_work_request *>::iterator wr_by_wrhash_iter_t;
typedef std::pair<uint32_t, ib_work_request *> wr_by_wrhash_t;
static nthread_lock_t nnti_wr_wrhash_lock;

typedef std::deque<ib_work_request *>           wr_pool_t;
typedef std::deque<ib_work_request *>::iterator wr_pool_iter_t;
static nthread_lock_t nnti_wr_pool_lock;

static wr_pool_t rdma_wr_pool;
static wr_pool_t sendrecv_wr_pool;


static nnti_ib_config config;




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
    uint32_t cqe_count;
    uint32_t srq_count;

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

        nthread_lock_init(&nnti_conn_peer_lock);
        nthread_lock_init(&nnti_conn_qpn_lock);
        nthread_lock_init(&nnti_wr_wrhash_lock);
        nthread_lock_init(&nnti_buf_bufhash_lock);

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
        transport_global_data.ctx = ibv_open_device(dev);
        if (!transport_global_data.ctx) {
            log_error(nnti_debug_level, "ibv_open_device failed");
            return NNTI_EIO;
        }

        transport_global_data.nic_port = 1;

        /* get the lid and verify port state */
        rc = ibv_query_port(transport_global_data.ctx, transport_global_data.nic_port, &dev_port_attr);
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
        rc = ibv_query_device(transport_global_data.ctx, &dev_attr);
        if (rc) {
            log_error(nnti_debug_level, "ibv_query_device failed");
            return NNTI_EIO;
        }

        log_debug(nnti_debug_level, "max %d completion queue entries", dev_attr.max_cqe);
        transport_global_data.cqe_count = dev_attr.max_cqe;

        log_debug(nnti_debug_level, "max %d shared receive queue work requests", dev_attr.max_srq_wr);
        transport_global_data.srq_count = dev_attr.max_srq_wr/2;

        log_debug(nnti_debug_level, "max %d queue pair work requests", dev_attr.max_qp_wr);
        transport_global_data.qp_count = dev_attr.max_qp_wr/2;

        /* Allocate a Protection Domain (global) */
        transport_global_data.pd = ibv_alloc_pd(transport_global_data.ctx);
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

        init_server_listen_socket();
//        start_connection_listener_thread();

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
 *                - ex. "ptl://nid:pid/", "ib://ip_addr:port", "luc://endpoint_id/"
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

    assert(trans_hdl);
    assert(peer_hdl);

    ib_connection *conn=get_conn_peer(peer_hdl);
    close_connection(conn);
    del_conn_peer(peer_hdl);

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
NNTI_result_t NNTI_ib_register_memory (
        const NNTI_transport_t *trans_hdl,
        char                   *buffer,
        const uint64_t          element_size,
        const uint64_t          num_elements,
        const NNTI_buf_ops_t    ops,
        const NNTI_peer_t      *peer,
        NNTI_buffer_t          *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;

    uint32_t cqe_num;

    struct ibv_recv_wr *bad_wr=NULL;

    NNTI_buffer_t     *old_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    assert(trans_hdl);
    assert(buffer);
    assert(element_size>0);
    assert(num_elements>0);
    assert(ops>0);
    assert(reg_buf);

//    if (ops==NNTI_PUT_SRC) print_xfer_buf(buffer, element_size);

    old_buf=get_buf_bufhash(hash6432shift((uint64_t)buffer));
    if (old_buf==NULL) {
        ib_mem_hdl=new ib_memory_handle();
        assert(ib_mem_hdl);
        ib_mem_hdl->ref_count=1;
    } else {
        ib_mem_hdl=(ib_memory_handle*)old_buf->transport_private;
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

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    num_elements*element_size,
                    IBV_ACCESS_LOCAL_WRITE);

            for (int i=0;i<cqe_num;i++) {
                post_recv_work_request(
                        reg_buf,
                        i,
                        (i*q_hdl->req_size),
                        q_hdl->req_size);
            }

        } else if (ops == NNTI_RECV_DST) {
            ib_mem_hdl->type=RECEIVE_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else if (ops == NNTI_SEND_SRC) {
            ib_mem_hdl->type=SEND_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else if (ops == NNTI_GET_DST) {
            ib_mem_hdl->type=GET_DST_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else if (ops == NNTI_GET_SRC) {
            ib_mem_hdl->type=GET_SRC_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else if (ops == NNTI_PUT_SRC) {
            //        print_xfer_buf(buffer, element_size);

            ib_mem_hdl->type=PUT_SRC_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

            //        print_xfer_buf(buffer, element_size);

        } else if (ops == NNTI_PUT_DST) {
            ib_mem_hdl->type=PUT_DST_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else if (ops == (NNTI_GET_SRC|NNTI_PUT_DST)) {
            ib_mem_hdl->type=RDMA_TARGET_BUFFER;

            register_memory(
                    ib_mem_hdl,
                    buffer,
                    element_size,
                    (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));

        } else {
            ib_mem_hdl->type=UNKNOWN_BUFFER;
        }
    }

    if (ops == NNTI_RECV_DST) {
        post_recv_work_request(
                reg_buf,
                -1,
                0,
                element_size);
    }

    if (config.use_rdma_target_ack) {
        if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
            (ib_mem_hdl->type == GET_SRC_BUFFER) ||
            (ib_mem_hdl->type == PUT_DST_BUFFER)) {
            post_ack_recv_work_request(reg_buf);
        }
    }

    reg_buf->buffer_addr.transport_id                     = NNTI_TRANSPORT_IB;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.size     = ib_mem_hdl->mr->length;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.buf      = (uint64_t)ib_mem_hdl->mr->addr;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.key      = ib_mem_hdl->mr->rkey;

    if (ib_mem_hdl->ref_count==1) {
        insert_buf_bufhash(reg_buf);
    }

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "end of NNTI_ib_register_memory", reg_buf);
    }

    return(rc);
}


/**
 * @brief Cleanup after network operations are complete.
 *
 * Destroy an NNTI_buffer_t that was previously created by NNTI_regsiter_buffer().
 * It is the user's responsibility to release the the memory region.
 */
NNTI_result_t NNTI_ib_unregister_memory (
        NNTI_buffer_t    *reg_buf)
{
    NNTI_result_t rc=NNTI_OK;
    ib_memory_handle *ib_mem_hdl=NULL;

    assert(reg_buf);

    if (logging_debug(nnti_debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "start of NNTI_ib_unregister_memory", reg_buf);
    }

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);
    ib_mem_hdl->ref_count--;

    log_debug(nnti_debug_level, "ib_mem_hdl->ref_count==%lu", ib_mem_hdl->ref_count);

    if (ib_mem_hdl->ref_count==0) {
        log_debug(nnti_debug_level, "ib_mem_hdl->ref_count is 0.  release all resources.");
        unregister_memory(ib_mem_hdl);

        del_buf_bufhash(reg_buf);

        while (!ib_mem_hdl->wr_queue.empty()) {
            ib_work_request *wr=ib_mem_hdl->wr_queue.front();
            log_debug(nnti_debug_level, "removing pending wr=%p", wr);
            ib_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            if (config.use_wr_pool) {
                if (wr->ack_mr!=NULL) {
                    wr_pool_rdma_push(wr);
                } else {
                    wr_pool_sendrecv_push(wr);
                }
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(wr);
                }
                free(wr);
            }
        }

        if (ib_mem_hdl) delete ib_mem_hdl;

        reg_buf->transport_id      = NNTI_TRANSPORT_NULL;
        IB_SET_MATCH_ANY(&reg_buf->buffer_owner);
        reg_buf->ops               = (NNTI_buf_ops_t)0;
        //    IB_SET_MATCH_ANY(&reg_buf->peer);
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
        const NNTI_buffer_t *dest_hdl)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(peer_hdl);
    assert(msg_hdl);

    log_level debug_level=nnti_debug_level; //LOG_ALL;

    ib_mem_hdl=(ib_memory_handle *)msg_hdl->transport_private;
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        wr=wr_pool_sendrecv_pop();
    } else {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
    }
    assert(wr);

    wr->conn = get_conn_peer(peer_hdl);
    assert(wr->conn);

    wr->reg_buf = (NNTI_buffer_t *)msg_hdl;

    wr->op_state = BUFFER_INIT;

    if ((dest_hdl == NULL) || (dest_hdl->ops == NNTI_RECV_QUEUE)) {
        wr->comp_channel=transport_global_data.req_comp_channel;
        wr->cq          =transport_global_data.req_cq;
        wr->qp          =wr->conn->req_qp.qp;
        wr->qpn         =(uint64_t)wr->conn->req_qp.qpn;
        wr->peer_qpn    =(uint64_t)wr->conn->peer_req_qpn;

        wr->sge.addr  =(uint64_t)ib_mem_hdl->mr->addr;
        wr->sge.length=ib_mem_hdl->mr->length;
        wr->sge.lkey  =ib_mem_hdl->mr->lkey;

        wr->sq_wr.wr_id  =hash6432shift((uint64_t)wr);
        wr->sq_wr.sg_list=&wr->sge;
        wr->sq_wr.num_sge=1;

        wr->sq_wr.opcode    =IBV_WR_SEND;
        wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

        wr->last_op=IB_OP_SEND_REQUEST;

    } else {
        wr->comp_channel=transport_global_data.data_comp_channel;
        wr->cq          =transport_global_data.data_cq;
        wr->qp          =wr->conn->data_qp.qp;
        wr->qpn         =(uint64_t)wr->conn->data_qp.qpn;
        wr->peer_qpn    =(uint64_t)wr->conn->data_qp.peer_qpn;

        wr->sge.addr  =(uint64_t)ib_mem_hdl->mr->addr;
        wr->sge.length=ib_mem_hdl->mr->length;
        wr->sge.lkey  =ib_mem_hdl->mr->lkey;

        wr->sq_wr.wr_id  =hash6432shift((uint64_t)wr);
        wr->sq_wr.sg_list=&wr->sge;
        wr->sq_wr.num_sge=1;

        wr->sq_wr.wr.rdma.rkey       =dest_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.key;
        wr->sq_wr.wr.rdma.remote_addr=dest_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf;

        wr->sq_wr.imm_data  =hash6432shift((uint64_t)dest_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf);
        wr->sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        wr->sq_wr.send_flags=IBV_SEND_SIGNALED;

        wr->last_op=IB_OP_SEND_BUFFER;
    }

    log_debug(nnti_debug_level, "sending to (%s, qp=%p, qpn=%lu, sge.addr=%p, sge.length=%llu, sq_wr.imm_data=%llx, sq_wr.wr.rdma.rkey=%x, sq_wr.wr.rdma.remote_addr=%p)",
            peer_hdl->url,
            wr->qp,
            wr->qpn,
            (void *)  wr->sge.addr,
            (uint64_t)wr->sge.length,
            (uint64_t)wr->sq_wr.imm_data,
                      wr->sq_wr.wr.rdma.rkey,
            (void *)  wr->sq_wr.wr.rdma.remote_addr);

    nthread_lock(&nnti_wr_wrhash_lock);

    trios_start_timer(call_time);
    if (ibv_post_send(wr->qp, &wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_send - ibv_post_send", call_time);

    ib_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&nnti_wr_wrhash_lock);

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
        const uint64_t       dest_offset)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ib_mem_hdl=(ib_memory_handle *)src_buffer_hdl->transport_private;
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        if (config.use_rdma_target_ack) {
            wr=wr_pool_rdma_pop();
        } else {
            wr=wr_pool_sendrecv_pop();
        }
    } else {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
    }
    assert(wr);

    wr->conn = get_conn_peer(&dest_buffer_hdl->buffer_owner);
    assert(wr->conn);

    wr->reg_buf = (NNTI_buffer_t *)src_buffer_hdl;

    wr->op_state=RDMA_WRITE_INIT;

    wr->comp_channel=transport_global_data.data_comp_channel;
    wr->cq          =transport_global_data.data_cq;
    wr->qp          =wr->conn->data_qp.qp;
    wr->qpn         =(uint64_t)wr->conn->data_qp.qpn;
    wr->peer_qpn    =(uint64_t)wr->conn->data_qp.peer_qpn;

    wr->sge.addr  =(uint64_t)ib_mem_hdl->mr->addr+src_offset;
    wr->sge.length=src_length;
    wr->sge.lkey  =ib_mem_hdl->mr->lkey;

    wr->sq_wr.sg_list=&wr->sge;
    wr->sq_wr.num_sge=1;

    wr->sq_wr.wr.rdma.rkey        = dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.key;
    wr->sq_wr.wr.rdma.remote_addr = dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf+dest_offset;

    wr->sq_wr.opcode    =IBV_WR_RDMA_WRITE;
    wr->sq_wr.send_flags=IBV_SEND_SIGNALED;
    wr->sq_wr.wr_id     =hash6432shift((uint64_t)wr);
    wr->sq_wr.imm_data  =hash6432shift((uint64_t)dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf);

    if (config.use_rdma_target_ack) {
        wr->ack.op    =IB_OP_PUT_TARGET;
        wr->ack.offset=dest_offset;
        wr->ack.length=src_length;

        if (!config.use_wr_pool) {
            register_ack(wr);
        }
        wr->ack_sge.addr  =(uint64_t)wr->ack_mr->addr;
        wr->ack_sge.length=wr->ack_mr->length;
        wr->ack_sge.lkey  =wr->ack_mr->lkey;

        wr->ack_sq_wr.sg_list=&wr->ack_sge;
        wr->ack_sq_wr.num_sge=1;

        wr->ack_sq_wr.wr.rdma.rkey       =dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.ack_key;
        wr->ack_sq_wr.wr.rdma.remote_addr=dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.ack_buf;

        wr->ack_sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        wr->ack_sq_wr.send_flags=IBV_SEND_SIGNALED|IBV_SEND_FENCE;
        wr->ack_sq_wr.wr_id     =hash6432shift((uint64_t)wr);
        wr->ack_sq_wr.imm_data  =hash6432shift((uint64_t)dest_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf);
    }

    wr->last_op=IB_OP_PUT_INITIATOR;
    wr->length=src_length;
    wr->offset=src_offset;

    log_debug(nnti_debug_level, "putting to (%s, qp=%p, qpn=%lu)",
            dest_buffer_hdl->buffer_owner.url,
            wr->qp,
            wr->qpn);

    nthread_lock(&nnti_wr_wrhash_lock);

    trios_start_timer(call_time);
    if (ibv_post_send(wr->qp, &wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_put - ibv_post_send", call_time);

    if (config.use_rdma_target_ack) {
        send_ack(wr);
    }

    ib_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&nnti_wr_wrhash_lock);

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
        const uint64_t       dest_offset)
{
    NNTI_result_t rc=NNTI_OK;

    trios_declare_timer(call_time);

    struct ibv_send_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(src_buffer_hdl);
    assert(dest_buffer_hdl);

    ib_mem_hdl=(ib_memory_handle *)dest_buffer_hdl->transport_private;
    assert(ib_mem_hdl);
    if (config.use_wr_pool) {
        if (config.use_rdma_target_ack) {
            wr=wr_pool_rdma_pop();
        } else {
            wr=wr_pool_sendrecv_pop();
        }
    } else {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
    }
    assert(wr);

    wr->conn = get_conn_peer(&src_buffer_hdl->buffer_owner);
    assert(wr->conn);

    wr->reg_buf = (NNTI_buffer_t *)dest_buffer_hdl;

    wr->op_state=RDMA_READ_INIT;

    wr->comp_channel=transport_global_data.data_comp_channel;
    wr->cq          =transport_global_data.data_cq;
    wr->qp          =wr->conn->data_qp.qp;
    wr->qpn         =(uint64_t)wr->conn->data_qp.qpn;
    wr->peer_qpn    =(uint64_t)wr->conn->data_qp.peer_qpn;

    wr->sge.addr  =(uint64_t)ib_mem_hdl->mr->addr+dest_offset;
    wr->sge.length=src_length;
    wr->sge.lkey  =ib_mem_hdl->mr->lkey;

    wr->sq_wr.sg_list=&wr->sge;
    wr->sq_wr.num_sge=1;

    wr->sq_wr.wr.rdma.rkey       =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.key;
    wr->sq_wr.wr.rdma.remote_addr=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf+src_offset;

    wr->sq_wr.opcode    =IBV_WR_RDMA_READ;
    wr->sq_wr.send_flags=IBV_SEND_SIGNALED;
    wr->sq_wr.wr_id     =hash6432shift((uint64_t)wr);
    wr->sq_wr.imm_data  =hash6432shift((uint64_t)src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf);

    if (config.use_rdma_target_ack) {
        wr->ack.op    =IB_OP_GET_TARGET;
        wr->ack.offset=src_offset;
        wr->ack.length=src_length;

        if (!config.use_wr_pool) {
            register_ack(wr);
        }
        wr->ack_sge.addr  =(uint64_t)wr->ack_mr->addr;
        wr->ack_sge.length=wr->ack_mr->length;
        wr->ack_sge.lkey  =wr->ack_mr->lkey;

        wr->ack_sq_wr.sg_list=&wr->ack_sge;
        wr->ack_sq_wr.num_sge=1;

        wr->ack_sq_wr.wr.rdma.rkey       =src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.ack_key;
        wr->ack_sq_wr.wr.rdma.remote_addr=src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.ack_buf;

        wr->ack_sq_wr.opcode    =IBV_WR_RDMA_WRITE_WITH_IMM;
        wr->ack_sq_wr.send_flags=IBV_SEND_SIGNALED|IBV_SEND_FENCE;
        wr->ack_sq_wr.wr_id     =hash6432shift((uint64_t)wr);
        wr->ack_sq_wr.imm_data  =hash6432shift((uint64_t)src_buffer_hdl->buffer_addr.NNTI_remote_addr_t_u.ib.buf);
    }

    wr->last_op=IB_OP_GET_INITIATOR;
    wr->length=src_length;
    wr->offset=dest_offset;

    log_debug(nnti_debug_level, "getting from (%s, qp=%p, qpn=%lu)",
            src_buffer_hdl->buffer_owner.url,
            wr->qp,
            wr->qpn);

    nthread_lock(&nnti_wr_wrhash_lock);

    trios_start_timer(call_time);
    if (ibv_post_send(wr->qp, &wr->sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
        rc=NNTI_EIO;
    }
    trios_stop_timer("NNTI_ib_get - ibv_post_send", call_time);

    if (config.use_rdma_target_ack) {
        send_ack(wr);
    }

//    print_wr(wr);

    ib_mem_hdl->wr_queue.push_back(wr);
    insert_wr_wrhash(wr);

    nthread_unlock(&nnti_wr_wrhash_lock);

//    print_wr(wr);

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
NNTI_result_t NNTI_ib_wait (
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *wr=NULL;
    ib_request_queue_handle *q_hdl=NULL;
    ib_connection           *conn=NULL;

    const NNTI_buffer_t  *wait_buf=NULL;

    int ibv_rc=0;
    NNTI_result_t rc;
    long elapsed_time = 0;
    long timeout_per_call;

    struct ibv_comp_channel *comp_channel;
    struct ibv_cq           *cq;

    log_level debug_level=nnti_debug_level;

    long entry_time=trios_get_time_ms();

    trios_declare_timer(call_time);
    trios_declare_timer(total_time);

    struct ibv_wc wc;

    trios_start_timer(total_time);

    log_debug(debug_level, "enter (reg_buf=%p)", reg_buf);

    assert(reg_buf);
    assert(status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_buffer(logger_get_file(), "reg_buf",
                "start of NNTI_ib_wait", reg_buf);
    }

    q_hdl     =&transport_global_data.req_queue;
    assert(q_hdl);

    check_listen_socket_for_new_connections();

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);
    wr=first_incomplete_wr(ib_mem_hdl);
    assert(wr);

    comp_channel=wr->comp_channel;
    cq          =wr->cq;

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            memset(status, 0, sizeof(NNTI_status_t));
            status->op     = remote_op;
            status->result = NNTI_EINVAL;
            return(NNTI_EINVAL);
        }
    }

    if (is_buf_op_complete(reg_buf) == TRUE) {
        log_debug(debug_level, "buffer op already complete (reg_buf=%p)", reg_buf);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "buffer op NOT complete (reg_buf=%p)", reg_buf);

        timeout_per_call = MIN_TIMEOUT;

        log_level old_log_level=logger_get_default_level();
        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

            if (is_buf_op_complete(reg_buf) == TRUE) {
                log_debug(debug_level, "buffer op completed outside this thread (reg_buf=%p)", reg_buf);
                nnti_rc = NNTI_OK;
                break;
            }

            memset(&wc, 0, sizeof(struct ibv_wc));
            trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
            ibv_rc = ibv_poll_cq(cq, 1, &wc);
            trios_stop_timer("NNTI_ib_wait - ibv_poll_cq", call_time);
            if (ibv_rc < 0) {
                log_debug(debug_level, "ibv_poll_cq failed: %d", ibv_rc);
                break;
            }
            log_debug(debug_level, "ibv_poll_cq(cq=%p) rc==%d", cq, ibv_rc);

            if (ibv_rc > 0) {
                log_debug(debug_level, "got wc from cq=%p", cq);
                log_debug(debug_level, "polling status is %s", ibv_wc_status_str(wc.status));

                print_wc(&wc);

                if (wc.status != IBV_WC_SUCCESS) {
                    log_error(debug_level, "Failed status %s (%d) for wr_id %lx",
                            ibv_wc_status_str(wc.status),
                            wc.status, wc.wr_id);
                    nnti_rc=NNTI_EIO;
                    break;
                }

                wait_buf=decode_event_buffer(reg_buf, &wc);
                process_event(wait_buf, &wc);
nthread_unlock(&nnti_ib_lock);

                if (is_buf_op_complete(reg_buf) == TRUE) {
                    nnti_rc = NNTI_OK;
                    break;
                }
            } else {
nthread_unlock(&nnti_ib_lock);
retry:
                check_listen_socket_for_new_connections();

                if (is_buf_op_complete(reg_buf) == TRUE) {
                    log_debug(debug_level, "buffer op completed outside this thread (reg_buf=%p)", reg_buf);
                    nnti_rc = NNTI_OK;
                    break;
                }

                trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
                rc = poll_comp_channel(comp_channel, cq, timeout_per_call);
nthread_unlock(&nnti_ib_lock);
                trios_stop_timer("NNTI_ib_wait - poll_comp_channel", call_time);
                /* case 1: success */
                if (rc == NNTI_OK) {
                    logger_set_default_level(old_log_level);
                    nnti_rc = NNTI_OK;
                    continue;
                }
                /* case 2: timed out */
                else if (rc==NNTI_ETIMEDOUT) {
                    elapsed_time = (trios_get_time_ms() - entry_time);

//                    elapsed_time += timeout_per_call;

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                        logger_set_default_level(old_log_level);
                        log_debug(debug_level, "poll_comp_channel timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }
                    /* continue if the timeout has not expired */
                    log_debug(debug_level, "poll_comp_channel timedout... retrying");

                    log_debug(debug_level, "***** disable debug logging.  will enable after polling success. *****");
                    logger_set_default_level(LOG_OFF);

                    goto retry;
                }
                /* case 3: poll was interupted */
                else if (rc==NNTI_EAGAIN) {
                    logger_set_default_level(old_log_level);
                    nnti_rc = NNTI_EAGAIN;
                    break;
                }
                /* case 4: failure */
                else {
                    logger_set_default_level(old_log_level);
                    log_error(debug_level, "poll_comp_channel failed (cq==%p): %s",
                            cq, strerror(errno));
                    nnti_rc = NNTI_EIO;
                    break;
                }
            }
        }
    }

    create_status(reg_buf, remote_op, nnti_rc, status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ib_wait", status);
    }

    if (nnti_rc==NNTI_OK) {
        ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;

        if (ib_mem_hdl->type == REQUEST_BUFFER) {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
        } else if (ib_mem_hdl->type == RECEIVE_BUFFER) {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
        }
        else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                (ib_mem_hdl->type == PUT_DST_BUFFER))     {
            if (config.use_rdma_target_ack) {
                wr=ib_mem_hdl->wr_queue.front();
                ib_mem_hdl->wr_queue.pop_front();
                repost_ack_recv_work_request((NNTI_buffer_t *)reg_buf, wr);
            }
        }
        else {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            if (config.use_wr_pool) {
                wr_pool_sendrecv_push(wr);
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(wr);
                }
                free(wr);
            }
        }

    }

    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_ib_wait", total_time);

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
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        uint32_t             *which,
        NNTI_status_t        *status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *wr=NULL;
    ib_request_queue_handle *q_hdl=NULL;
    ib_connection           *conn=NULL;
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

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    check_listen_socket_for_new_connections();

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            memset(status, 0, sizeof(NNTI_status_t));
            status->op     = remote_op;
            status->result = NNTI_EINVAL;
            return(NNTI_EINVAL);
        }
    }

    assert(buf_list);
    assert(buf_count > 0);
    if (buf_count > 1) {
        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
        for (int i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                assert(((ib_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_ib_wait(buf_list[0], remote_op, timeout, status);
        *which=0;
        goto cleanup;
    }

    if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
        log_debug(debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "buffer op NOT complete (buf_list=%p)", buf_list);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

            if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                log_debug(debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
                nnti_rc = NNTI_OK;
                break;
            }

            memset(&wc, 0, sizeof(struct ibv_wc));
            trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
            ibv_rc = ibv_poll_cq(transport_global_data.data_cq, 1, &wc);
            trios_stop_timer("NNTI_ib_waitany - ibv_poll_cq", call_time);
            if (ibv_rc < 0) {
                log_debug(debug_level, "ibv_poll_cq failed: %d", ibv_rc);
                break;
            }
            log_debug(debug_level, "ibv_poll_cq(cq=%p) rc==%d", transport_global_data.data_cq, ibv_rc);

            if (ibv_rc > 0) {
                log_debug(debug_level, "got wc from cq=%p", transport_global_data.data_cq);
                log_debug(debug_level, "polling status is %s", ibv_wc_status_str(wc.status));

                print_wc(&wc);

                if (wc.status != IBV_WC_SUCCESS) {
                    log_error(debug_level, "Failed status %s (%d) for wr_id %lu",
                            ibv_wc_status_str(wc.status),
                            wc.status, wc.wr_id);
                    nnti_rc=NNTI_EIO;
                    break;
                }

                wait_buf=decode_event_buffer(buf_list[0], &wc);
                process_event(wait_buf, &wc);
nthread_unlock(&nnti_ib_lock);

                if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                    nnti_rc = NNTI_OK;
                    break;
                }
            } else {
nthread_unlock(&nnti_ib_lock);
retry:
                check_listen_socket_for_new_connections();

                if (is_any_buf_op_complete(buf_list, buf_count, which) == TRUE) {
                    log_debug(debug_level, "buffer op already complete (which=%u, buf_list[%d]=%p)", *which, *which, buf_list[*which]);
                    nnti_rc = NNTI_OK;
                    break;
                }

                trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
                rc = poll_comp_channel(transport_global_data.data_comp_channel, transport_global_data.data_cq, timeout_per_call);
nthread_unlock(&nnti_ib_lock);
                trios_stop_timer("NNTI_ib_waitany - poll_comp_channel", call_time);
                /* case 1: success */
                if (rc == NNTI_OK) {
                    nnti_rc = NNTI_OK;
                    continue;
                }
                /* case 2: timed out */
                else if (rc==NNTI_ETIMEDOUT) {
                    elapsed_time = (trios_get_time_ms() - entry_time);

//                    elapsed_time += timeout_per_call;

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                        log_debug(debug_level, "poll_comp_channel timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }
                    /* continue if the timeout has not expired */
                    log_debug(debug_level, "poll_comp_channel timedout... retrying");

                    goto retry;
                }
                /* case 3: poll was interupted */
                else if (rc==NNTI_EAGAIN) {
                    nnti_rc = NNTI_EAGAIN;
                    break;
                }
                /* case 4: failure */
                else {
                    log_error(debug_level, "poll_comp_channel failed (cq==%p): %s",
                            transport_global_data.data_cq, strerror(errno));
                    nnti_rc = NNTI_EIO;
                    break;
                }
            }
        }
    }

    create_status(buf_list[*which], remote_op, nnti_rc, status);

    if (logging_debug(debug_level)) {
        fprint_NNTI_status(logger_get_file(), "status",
                "end of NNTI_ib_waitany", status);
    }

    if (nnti_rc==NNTI_OK) {
        ib_mem_hdl=(ib_memory_handle *)buf_list[*which]->transport_private;

        if (ib_mem_hdl->type == REQUEST_BUFFER) {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)buf_list[*which], wr);
        } else if (ib_mem_hdl->type == RECEIVE_BUFFER) {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            repost_recv_work_request((NNTI_buffer_t *)buf_list[*which], wr);
        }
        else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                (ib_mem_hdl->type == PUT_DST_BUFFER))     {
            if (config.use_rdma_target_ack) {
                wr=ib_mem_hdl->wr_queue.front();
                ib_mem_hdl->wr_queue.pop_front();
                repost_ack_recv_work_request((NNTI_buffer_t *)buf_list[*which], wr);
            }
        }
        else {
            wr=ib_mem_hdl->wr_queue.front();
            ib_mem_hdl->wr_queue.pop_front();
            del_wr_wrhash(wr);
            if (config.use_wr_pool) {
                wr_pool_sendrecv_push(wr);
            } else {
                if (config.use_rdma_target_ack) {
                    unregister_ack(wr);
                }
                free(wr);
            }
        }
    }

cleanup:
    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_ib_waitany", total_time);

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
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        const NNTI_buf_ops_t  remote_op,
        const int             timeout,
        NNTI_status_t       **status)
{
    NNTI_result_t nnti_rc=NNTI_OK;
    ib_memory_handle        *ib_mem_hdl=NULL;
    ib_work_request         *wr=NULL;
    ib_request_queue_handle *q_hdl=NULL;
    ib_connection           *conn=NULL;
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

    trios_start_timer(total_time);

    log_debug(debug_level, "enter");

    check_listen_socket_for_new_connections();

    if (!config.use_rdma_target_ack) {
        if ((remote_op==NNTI_GET_SRC) || (remote_op==NNTI_PUT_DST) || (remote_op==(NNTI_GET_SRC|NNTI_PUT_DST))) {
            for (int i=0;i<buf_count;i++) {
                memset(status[i], 0, sizeof(NNTI_status_t));
                status[i]->op     = remote_op;
                status[i]->result = NNTI_EINVAL;
            }
            return(NNTI_EINVAL);
        }
    }

    assert(buf_list);
    assert(buf_count > 0);
    if (buf_count > 1) {
        /* if there is more than 1 buffer in the list, none of them can be a REQUEST_BUFFER */
        for (int i=0;i<buf_count;i++) {
            if (buf_list[i] != NULL) {
                assert(((ib_memory_handle *)buf_list[i]->transport_private)->type != REQUEST_BUFFER);
            }
        }
    }
    assert(status);

    if (buf_count == 1) {
        nnti_rc=NNTI_ib_wait(buf_list[0], remote_op, timeout, status[0]);
        goto cleanup;
    }

    if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
        log_debug(debug_level, "all buffer ops already complete (buf_list=%p)", buf_list);
        nnti_rc = NNTI_OK;
    } else {
        log_debug(debug_level, "all buffer ops NOT complete (buf_list=%p)", buf_list);

        timeout_per_call = MIN_TIMEOUT;

        while (1)   {
            if (trios_exit_now()) {
                log_debug(debug_level, "caught abort signal");
                nnti_rc=NNTI_ECANCELED;
                break;
            }

            check_listen_socket_for_new_connections();

            if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                log_debug(debug_level, "all buffer ops already complete (buf_list=%p)", buf_list);
                nnti_rc = NNTI_OK;
                break;
            }

            memset(&wc, 0, sizeof(struct ibv_wc));
            trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
            ibv_rc = ibv_poll_cq(transport_global_data.data_cq, 1, &wc);
            trios_stop_timer("NNTI_ib_waitany - ibv_poll_cq", call_time);
            if (ibv_rc < 0) {
                log_debug(debug_level, "ibv_poll_cq failed: %d", ibv_rc);
                break;
            }
            log_debug(debug_level, "ibv_poll_cq(cq=%p) rc==%d", transport_global_data.data_cq, ibv_rc);

            if (ibv_rc > 0) {
                log_debug(debug_level, "got wc from cq=%p", transport_global_data.data_cq);
                log_debug(debug_level, "polling status is %s", ibv_wc_status_str(wc.status));

                print_wc(&wc);

                if (wc.status != IBV_WC_SUCCESS) {
                    log_error(debug_level, "Failed status %s (%d) for wr_id %lu",
                            ibv_wc_status_str(wc.status),
                            wc.status, wc.wr_id);
                    nnti_rc=NNTI_EIO;
                    break;
                }

                wait_buf=decode_event_buffer(buf_list[0], &wc);
                process_event(wait_buf, &wc);
nthread_unlock(&nnti_ib_lock);

                if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                    nnti_rc = NNTI_OK;
                    break;
                }
            } else {
nthread_unlock(&nnti_ib_lock);
retry:
                check_listen_socket_for_new_connections();

                if (is_all_buf_ops_complete(buf_list, buf_count) == TRUE) {
                    log_debug(debug_level, "all buffer ops already complete (buf_list=%p)", buf_list);
                    nnti_rc = NNTI_OK;
                    break;
                }

                trios_start_timer(call_time);
nthread_lock(&nnti_ib_lock);
                rc = poll_comp_channel(transport_global_data.data_comp_channel, transport_global_data.data_cq, timeout_per_call);
nthread_unlock(&nnti_ib_lock);
                trios_stop_timer("NNTI_ib_waitany - poll_comp_channel", call_time);
                /* case 1: success */
                if (rc == NNTI_OK) {
                    nnti_rc = NNTI_OK;
                    continue;
                }
                /* case 2: timed out */
                else if (rc==NNTI_ETIMEDOUT) {
                    elapsed_time = (trios_get_time_ms() - entry_time);

//                  elapsed_time += timeout_per_call;

                    /* if the caller asked for a legitimate timeout, we need to exit */
                    if (((timeout > 0) && (elapsed_time >= timeout)) || trios_exit_now()) {
                        log_debug(debug_level, "poll_comp_channel timed out");
                        nnti_rc = NNTI_ETIMEDOUT;
                        break;
                    }
                    /* continue if the timeout has not expired */
                    log_debug(debug_level, "poll_comp_channel timedout... retrying");

                    goto retry;
                }
                /* case 3: poll was interupted */
                else if (rc==NNTI_EAGAIN) {
                    nnti_rc = NNTI_EAGAIN;
                    break;
                }
                /* case 4: failure */
                else {
                    log_error(debug_level, "poll_comp_channel failed (cq==%p): %s",
                            transport_global_data.data_cq, strerror(errno));
                    nnti_rc = NNTI_EIO;
                    break;
                }
            }
        }
    }


    for (int i=0;i<buf_count;i++) {
        create_status(buf_list[i], remote_op, nnti_rc, status[i]);

        if (nnti_rc==NNTI_OK) {
            ib_mem_hdl=(ib_memory_handle *)buf_list[i]->transport_private;

            if (ib_mem_hdl->type == REQUEST_BUFFER) {
                wr=ib_mem_hdl->wr_queue.front();
                ib_mem_hdl->wr_queue.pop_front();
                repost_recv_work_request((NNTI_buffer_t *)buf_list[i], wr);
            } else if (ib_mem_hdl->type == RECEIVE_BUFFER) {
                wr=ib_mem_hdl->wr_queue.front();
                ib_mem_hdl->wr_queue.pop_front();
                repost_recv_work_request((NNTI_buffer_t *)buf_list[i], wr);
            }
            else if ((ib_mem_hdl->type == RDMA_TARGET_BUFFER) ||
                    (ib_mem_hdl->type == GET_SRC_BUFFER)      ||
                    (ib_mem_hdl->type == PUT_DST_BUFFER))     {
                if (config.use_rdma_target_ack) {
                    wr=ib_mem_hdl->wr_queue.front();
                    ib_mem_hdl->wr_queue.pop_front();
                    repost_ack_recv_work_request((NNTI_buffer_t *)buf_list[i], wr);
                }
            }
            else {
                wr=ib_mem_hdl->wr_queue.front();
                ib_mem_hdl->wr_queue.pop_front();
                del_wr_wrhash(wr);
                if (config.use_wr_pool) {
                    wr_pool_sendrecv_push(wr);
                } else {
                    if (config.use_rdma_target_ack) {
                        unregister_ack(wr);
                    }
                    free(wr);
                }
            }
        }
    }

cleanup:
    log_debug(debug_level, "exit");

    trios_stop_timer("NNTI_ib_waitall", total_time);

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

    ibv_destroy_comp_channel(transport_global_data.data_comp_channel);
    ibv_destroy_cq(transport_global_data.data_cq);
    ibv_destroy_srq(transport_global_data.data_srq);

    ibv_destroy_comp_channel(transport_global_data.req_comp_channel);
    ibv_destroy_cq(transport_global_data.req_cq);
    ibv_destroy_srq(transport_global_data.req_srq);

    ibv_dealloc_pd(transport_global_data.pd);

    ibv_close_device(transport_global_data.ctx);

    nthread_lock_fini(&nnti_ib_lock);
    nthread_lock_fini(&nnti_conn_peer_lock);
    nthread_lock_fini(&nnti_conn_qpn_lock);
    nthread_lock_fini(&nnti_wr_wrhash_lock);
    nthread_lock_fini(&nnti_buf_bufhash_lock);
    nthread_lock_fini(&nnti_wr_pool_lock);

    ib_initialized=false;

    return(NNTI_OK);
}





static NNTI_result_t setup_request_channel(void)
{
    int flags;

    transport_global_data.req_comp_channel = ibv_create_comp_channel(transport_global_data.ctx);
    if (!transport_global_data.req_comp_channel) {
        log_error(nnti_debug_level, "ibv_create_comp_channel failed");
        return NNTI_EIO;
    }
    transport_global_data.req_cq = ibv_create_cq(
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

    struct ibv_srq_init_attr attr;
    memset(&attr, 0, sizeof(attr)); // initialize to avoid valgrind uninitialized warning
    attr.attr.max_wr = transport_global_data.srq_count;
    attr.attr.max_sge = 1;

    transport_global_data.req_srq = ibv_create_srq(transport_global_data.pd, &attr);
    if (!transport_global_data.req_srq)  {
        log_error(nnti_debug_level, "ibv_create_srq failed");
        return NNTI_EIO;
    }

    if (ibv_req_notify_cq(transport_global_data.req_cq, 0)) {
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

    transport_global_data.data_comp_channel = ibv_create_comp_channel(transport_global_data.ctx);
    if (!transport_global_data.data_comp_channel) {
        log_error(nnti_debug_level, "ibv_create_comp_channel failed");
        return NNTI_EIO;
    }
    transport_global_data.data_cq = ibv_create_cq(
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

    struct ibv_srq_init_attr attr;
    memset(&attr, 0, sizeof(attr));   // initialize to avoid valgrind warning
    attr.attr.max_wr  = transport_global_data.srq_count;
    attr.attr.max_sge = 1;

    transport_global_data.data_srq = ibv_create_srq(transport_global_data.pd, &attr);
    if (!transport_global_data.data_srq)  {
        log_error(nnti_debug_level, "ibv_create_srq failed");
        return NNTI_EIO;
    }

    if (ibv_req_notify_cq(transport_global_data.data_cq, 0)) {
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

static int register_memory(
        ib_memory_handle *hdl,
        void *buf,
        uint64_t len,
        enum ibv_access_flags access)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    if (hdl->type != REQUEST_BUFFER) log_debug(nnti_debug_level, "enter buffer(%p) len(%d)", buf, len);

    trios_start_timer(callTime);
    mlock(buf, len);
    munlock(buf, len);
    trios_stop_timer("mlock", callTime);

    trios_start_timer(callTime);
    mr = ibv_reg_mr(transport_global_data.pd, buf, len, access);
    if (!mr) {
        if (errno == EFAULT) {
            log_debug(nnti_debug_level, "ibv_reg_mr failed with EFAULT.  trying to register with IBV_ACCESS_REMOTE_READ.");
            mr = ibv_reg_mr(transport_global_data.pd, buf, len, IBV_ACCESS_REMOTE_READ);
            if (!mr) {
                log_error(nnti_debug_level, "failed to register memory region with IBV_ACCESS_REMOTE_READ: %s", strerror(errno));
                return(errno);
            }
        } else {
            log_error(nnti_debug_level, "failed to register memory region (IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE): %s", strerror(errno));
            return errno;
        }
    }
    trios_stop_timer("register", callTime);

    hdl->mr=mr;

    if (hdl->type != REQUEST_BUFFER) log_debug(nnti_debug_level, "exit (buf==%p, mr==%p, lkey %x, rkey %x)...", buf, mr, mr->lkey, mr->rkey);

    return (rc);
}

static int register_ack(ib_work_request *wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    uint32_t len;

    log_debug(nnti_debug_level, "enter");

    len = sizeof(wr->ack);

    trios_start_timer(callTime);
    mlock(&wr->ack, len);
    munlock(&wr->ack, len);
    trios_stop_timer("mlock", callTime);

    trios_start_timer(callTime);
    mr = ibv_reg_mr(transport_global_data.pd, &wr->ack, len, (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));
    if (!mr) {
        log_error(nnti_debug_level, "failed to register memory region");
        perror("errno");
        return errno;
    }
    trios_stop_timer("register", callTime);

    wr->ack_mr=mr;

    log_debug(nnti_debug_level, "exit (mr==%p, addr %p, length %lu, lkey %x, rkey %x)...", mr, mr->addr, mr->length, mr->lkey, mr->rkey);

    return (rc);
}

static int unregister_memory(ib_memory_handle *hdl)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int i=0;
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    trios_start_timer(callTime);
    ibv_rc=ibv_dereg_mr(hdl->mr);
    if (ibv_rc != 0) {
        log_error(nnti_debug_level, "deregistering the memory buffer failed");
    }
    trios_stop_timer("deregister", callTime);

    log_debug(nnti_debug_level, "exit");

    return (rc);
}

static int unregister_ack(ib_work_request *wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    if (wr->ack_mr!=NULL) {
        trios_start_timer(callTime);
        ibv_rc=ibv_dereg_mr(wr->ack_mr);
        if (ibv_rc != 0) {
            log_error(nnti_debug_level, "deregistering the ACK buffer failed");
        }
        trios_stop_timer("deregister", callTime);
    }

    log_debug(nnti_debug_level, "exit");

    return (rc);
}

static void send_ack (
        ib_work_request *wr)
{
    struct ibv_send_wr *bad_wr;

    NNTI_buffer_t    *reg_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;

    reg_buf=wr->reg_buf;
    assert(reg_buf);
    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    log_debug(nnti_debug_level, "sending ACK to (ack_sge.addr=%p, ack_sge.length=%llu, ack_sq_wr.wr.rdma.rkey=%x, ack_sq_wr.wr.rdma.remote_addr=%p)",
            (void *)  wr->ack_sge.addr,
            (uint64_t)wr->ack_sge.length,
                      wr->ack_sq_wr.wr.rdma.rkey,
            (void *)  wr->ack_sq_wr.wr.rdma.remote_addr);

    if (ibv_post_send(wr->qp, &wr->ack_sq_wr, &bad_wr)) {
        log_error(nnti_debug_level, "failed to post send: %s", strerror(errno));
    }

    return;
}

static ib_work_request *decode_work_request(
        const struct ibv_wc *wc)
{
    ib_work_request  *wr=NULL;

    if (wc->imm_data != 0) {
        wr = get_wr_wrhash(wc->imm_data);
    } else {
        wr = get_wr_wrhash(wc->wr_id);
    }

    return(wr);
}

static const NNTI_buffer_t *decode_event_buffer(
        const NNTI_buffer_t *wait_buf,
        const struct ibv_wc *wc)
{
    const NNTI_buffer_t *event_buf=NULL;
    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    if (wc->opcode==IBV_WC_RECV) {
        log_debug(nnti_debug_level, "wc->opcode is IBV_WC_RECV, so wc.wr_id is the request buffer index.");

        event_buf=transport_global_data.req_queue.reg_buf;
        ib_mem_hdl=(ib_memory_handle *)event_buf->transport_private;
        assert(ib_mem_hdl);

    } else if (wc->opcode==IBV_WC_SEND) {
        log_debug(nnti_debug_level, "wc->opcode is IBV_WC_SEND, so wc.wr_id is wr hash (IB_OP_SEND_REQUEST).");

        wr = get_wr_wrhash(wc->wr_id);
        assert(wr);
        event_buf=wr->reg_buf;

    } else {
        if (wc->imm_data == 0) {
            // This is not a request buffer and I am the initiator, so wc.wr_id is the hash of the work request
            wr = get_wr_wrhash(wc->wr_id);
            assert(wr);
            event_buf=wr->reg_buf;
        } else {
            // This is not a request buffer and I am the target, so wc.imm_data is the hash of either the buffer or the work request
            event_buf = get_buf_bufhash(wc->imm_data);
            if (event_buf == NULL) {
                wr = get_wr_wrhash(wc->imm_data);
                assert(wr);
                event_buf=wr->reg_buf;
            }
            assert(event_buf);
        }
        assert(event_buf);

        if (event_buf == wait_buf) {
            log_debug(nnti_debug_level, "the wc matches the wait buffer (cq=%p, wr_id=%p, wait_buf=%p)",
                    (void *)transport_global_data.data_cq, wc->wr_id, wait_buf);
        } else {
            log_debug(nnti_debug_level, "the wc does NOT match the wait buffer (cq=%p, wr_id=%p, wait_buf=%p)",
                    (void *)transport_global_data.data_cq, wc->wr_id, wait_buf);
        }
    }

    log_debug(nnti_debug_level, "exit (event_buf==%p)", event_buf);

    return(event_buf);
}

int process_event(
        const NNTI_buffer_t  *wait_buf,
        const struct ibv_wc  *wc)
{
    NNTI_result_t rc=NNTI_OK;

    const NNTI_buffer_t *event_buf=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_level debug_level=nnti_debug_level;

    if (wc->status != IBV_WC_SUCCESS) {
        return NNTI_EIO;
    }


    event_buf=wait_buf;
    ib_mem_hdl=(ib_memory_handle *)event_buf->transport_private;
    assert(ib_mem_hdl);

    wr = decode_work_request(wc);
    if (wr == NULL) {
        wr=first_incomplete_wr(ib_mem_hdl);
    }
    assert(wr);

    log_debug(nnti_debug_level, "event_buf=%p; wr=%p; wr->last_op=%d", event_buf, wr, wr->last_op);
    log_debug(nnti_debug_level, "ib_mem_hdl->type==%llu", (uint64_t)ib_mem_hdl->type);

    wr->last_wc = *wc;

    debug_level=nnti_debug_level;
    switch (ib_mem_hdl->type) {
        case SEND_BUFFER:
            if (wc->opcode==IBV_WC_SEND) {
                log_debug(debug_level, "send completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->op_state = SEND_COMPLETE;
            }
            if (wc->opcode==IBV_WC_RDMA_WRITE) {
                log_debug(debug_level, "send completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->op_state = SEND_COMPLETE;
            }
            break;
        case PUT_SRC_BUFFER:
            if (wc->opcode==IBV_WC_RDMA_WRITE) {
                log_debug(debug_level, "RDMA write event - wc==%p, event_buf==%p, op_state==%d", wc, event_buf, wr->op_state);
                if (wr->op_state==RDMA_WRITE_INIT) {
                    log_debug(debug_level, "RDMA write (initiator) completion - wc==%p, event_buf==%p", wc, event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_WRITE_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_WRITE_COMPLETE;
                    }
                }
                else if ((config.use_rdma_target_ack) &&
                        (wr->op_state==RDMA_WRITE_NEED_ACK)) {
                    log_debug(debug_level, "RDMA write ACK (initiator) completion - wc==%p, event_buf==%p", wc, event_buf);
                    wr->last_op=IB_OP_PUT_INITIATOR;
                    wr->op_state = RDMA_WRITE_COMPLETE;
                }
            }
//            if (wr->op_state == RDMA_WRITE_COMPLETE) {
//                print_xfer_buf((void *)wr->reg_buf->payload, wr->reg_buf->payload_size);
//                print_ack_buf(&wr->ack);
//            }
            break;
        case GET_DST_BUFFER:
            if (wc->opcode==IBV_WC_RDMA_READ) {
                log_debug(debug_level, "RDMA read event - wc==%p, event_buf==%p, op_state==%d", wc, event_buf, wr->op_state);
                if (wr->op_state==RDMA_READ_INIT) {
                    log_debug(debug_level, "RDMA read (initiator) completion - wc==%p, event_buf==%p", wc, event_buf);
                    if (config.use_rdma_target_ack) {
                        wr->op_state=RDMA_READ_NEED_ACK;
                    } else {
                        wr->op_state = RDMA_READ_COMPLETE;
                    }
                }
            }
            else if ((config.use_rdma_target_ack) &&
                    (wc->opcode==IBV_WC_RDMA_WRITE)) {
                if (wr->op_state==RDMA_READ_NEED_ACK) {
                    log_debug(debug_level, "RDMA read ACK (initiator) completion - wc==%p, event_buf==%p", wc, event_buf);
                    wr->last_op=IB_OP_GET_INITIATOR;
                    wr->op_state = RDMA_READ_COMPLETE;
                }
            }
//            if (wr->op_state == RDMA_READ_COMPLETE) {
//                print_xfer_buf((void *)wr->reg_buf->payload, wr->reg_buf->payload_size);
//                print_ack_buf(&wr->ack);
//            }
            break;
        case REQUEST_BUFFER:
            if (wc->opcode==IBV_WC_RECV) {
                log_debug(debug_level, "recv completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->last_op=IB_OP_NEW_REQUEST;
                wr->op_state = RECV_COMPLETE;
                if (transport_global_data.req_queue.req_received == transport_global_data.srq_count) {
                    log_debug(debug_level, "resetting req_queue.req_received to 0");
                    transport_global_data.req_queue.req_received=0;
                }
                if (transport_global_data.req_queue.req_received != wc->wr_id) {
                    log_debug(debug_level, "req_queue.req_received(%llu) != wc->wr_id(%llu)", transport_global_data.req_queue.req_received, wc->wr_id);
                }
                transport_global_data.req_queue.req_received++;

                if (wr->cq == transport_global_data.req_cq) {
                    transport_global_data.req_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
                }
            }
            break;
        case RECEIVE_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "recv completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->last_op=IB_OP_RECEIVE;
                wr->op_state = RECV_COMPLETE;

                if (wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
            break;
        case PUT_DST_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "RDMA write (target) completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->last_op=IB_OP_PUT_TARGET;
                wr->op_state = RDMA_WRITE_COMPLETE;

                if (wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
//            if (wr->op_state == RDMA_WRITE_COMPLETE) {
//                print_xfer_buf((void *)wr->reg_buf->payload, wr->reg_buf->payload_size);
//                print_ack_buf(&wr->ack);
//            }
            break;
        case GET_SRC_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "RDMA read (target) completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->last_op=IB_OP_GET_TARGET;
                wr->op_state = RDMA_READ_COMPLETE;

                if (wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
//            if (wr->op_state == RDMA_READ_COMPLETE) {
//                print_xfer_buf((void *)wr->reg_buf->payload, wr->reg_buf->payload_size);
//                print_ack_buf(&wr->ack);
//            }
            break;
        case RDMA_TARGET_BUFFER:
            if (wc->opcode==IBV_WC_RECV_RDMA_WITH_IMM) {
                log_debug(debug_level, "RDMA target completion - wc==%p, event_buf==%p", wc, event_buf);
                wr->op_state = RDMA_TARGET_COMPLETE;
                if (config.use_rdma_target_ack) {
                    wr->last_op=wr->ack.op;
                }

                if (wr->cq == transport_global_data.data_cq) {
                    transport_global_data.data_srq_count--;
                    log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
                }
            }
//            if (wr->op_state == RDMA_TARGET_COMPLETE) {
//                print_xfer_buf((void *)wr->reg_buf->payload, wr->reg_buf->payload_size);
//                print_ack_buf(&wr->ack);
//            }
            break;
    }

    return (rc);
}

static NNTI_result_t post_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        int64_t         wr_id,
        uint64_t        offset,
        uint64_t        length)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    struct ibv_srq *srq=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
    } else {
        if (config.use_wr_pool) {
            wr=wr_pool_sendrecv_pop();
        } else {
            wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        }
    }
    assert(wr);
    wr->reg_buf = reg_buf;
    wr->offset  = offset;
    wr->length  = length;

    wr->op_state = BUFFER_INIT;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        wr->op_state=BUFFER_INIT;
        wr->last_op=IB_OP_NEW_REQUEST;

    } else if (ib_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (ib_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        wr->comp_channel=transport_global_data.req_comp_channel;
        wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        wr->comp_channel=transport_global_data.data_comp_channel;
        wr->cq          =transport_global_data.data_cq;
    }

    wr->sge.addr  =(uint64_t)ib_mem_hdl->mr->addr+offset;
    wr->sge.length=length;
    wr->sge.lkey  =ib_mem_hdl->mr->lkey;

    if (wr_id == -1) {
        wr->rq_wr.wr_id  =(uint64_t)wr;
    } else {
        wr->rq_wr.wr_id  =wr_id;
    }
    wr->rq_wr.sg_list=&wr->sge;
    wr->rq_wr.num_sge=1;

    if ((wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    ib_mem_hdl->wr_queue.push_back(wr);

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
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    if (config.use_wr_pool) {
        wr=wr_pool_rdma_pop();
    } else {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
    }
    assert(wr);
    wr->reg_buf = reg_buf;

    wr->op_state = BUFFER_INIT;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        wr->op_state=BUFFER_INIT;
        wr->last_op=IB_OP_NEW_REQUEST;

    } else if (ib_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (ib_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        wr->comp_channel=transport_global_data.req_comp_channel;
        wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        wr->comp_channel=transport_global_data.data_comp_channel;
        wr->cq          =transport_global_data.data_cq;
    }

    if (!config.use_wr_pool) {
        register_ack(wr);
    }

    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.ack_size = wr->ack_mr->length;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.ack_buf  = (uint64_t)wr->ack_mr->addr;
    reg_buf->buffer_addr.NNTI_remote_addr_t_u.ib.ack_key  = wr->ack_mr->rkey;

    wr->sge.addr  =(uint64_t)wr->ack_mr->addr;
    wr->sge.length=wr->ack_mr->length;
    wr->sge.lkey  =wr->ack_mr->lkey;

    wr->rq_wr.wr_id  =(uint64_t)wr;
    wr->rq_wr.sg_list=&wr->sge;
    wr->rq_wr.num_sge=1;

    if ((wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    ib_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_recv_work_request(
        NNTI_buffer_t  *reg_buf,
        ib_work_request  *wr)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;

    ib_memory_handle *ib_mem_hdl=NULL;

    struct ibv_srq *srq=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    assert(wr);

    wr->op_state = BUFFER_INIT;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        wr->op_state=BUFFER_INIT;
        wr->last_op=IB_OP_NEW_REQUEST;

    } else if (ib_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (ib_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        wr->comp_channel=transport_global_data.req_comp_channel;
        wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        wr->comp_channel=transport_global_data.data_comp_channel;
        wr->cq          =transport_global_data.data_cq;
    }

    if ((wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    ib_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static NNTI_result_t repost_ack_recv_work_request(
        NNTI_buffer_t    *reg_buf,
        ib_work_request  *wr)
{
    int ibv_rc=0;

    struct ibv_recv_wr *bad_wr=NULL;
    struct ibv_srq   *srq;

    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter (reg_buf=%p)", reg_buf);

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    wr->op_state = BUFFER_INIT;

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        wr->op_state=BUFFER_INIT;
        wr->last_op=IB_OP_NEW_REQUEST;

    } else if (ib_mem_hdl->type==RECEIVE_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==GET_SRC_BUFFER) {
        wr->op_state=RDMA_READ_INIT;

    } else if (ib_mem_hdl->type==PUT_DST_BUFFER) {
        wr->op_state=RDMA_WRITE_INIT;

    } else if (ib_mem_hdl->type==RDMA_TARGET_BUFFER) {
        wr->op_state=RDMA_TARGET_INIT;

    }

    if (ib_mem_hdl->type==REQUEST_BUFFER) {
        srq             =transport_global_data.req_srq;
        wr->comp_channel=transport_global_data.req_comp_channel;
        wr->cq          =transport_global_data.req_cq;
    } else {
        srq             =transport_global_data.data_srq;
        wr->comp_channel=transport_global_data.data_comp_channel;
        wr->cq          =transport_global_data.data_cq;
    }

    if ((wr->cq == transport_global_data.data_cq) && (transport_global_data.data_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.data_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.data_srq_count==%ld", transport_global_data.data_srq_count);
    }
    if ((wr->cq == transport_global_data.req_cq) && (transport_global_data.req_srq_count < (transport_global_data.srq_count/2))) {
        if (ibv_rc=ibv_post_srq_recv(srq, &wr->rq_wr, &bad_wr)) {
            log_error(nnti_debug_level, "failed to post SRQ recv (rq_wr=%p ; bad_wr=%p): %s",
                    &wr->rq_wr, bad_wr, strerror(ibv_rc));
            return (NNTI_result_t)errno;
        }

        transport_global_data.req_srq_count++;
        log_debug(nnti_debug_level, "transport_global_data.req_srq_count==%ld", transport_global_data.req_srq_count);
    }

    ib_mem_hdl->wr_queue.push_back(wr);

    log_debug(nnti_debug_level, "exit (reg_buf=%p)", reg_buf);

    return(NNTI_OK);
}

static int8_t is_wr_complete(
        const ib_work_request *wr)
{
    int8_t rc=FALSE;
    ib_memory_handle *ib_mem_hdl=NULL;

    ib_mem_hdl=(ib_memory_handle *)wr->reg_buf->transport_private;
    assert(ib_mem_hdl);

    switch (ib_mem_hdl->type) {
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
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static ib_work_request *first_incomplete_wr(
        ib_memory_handle *ib_mem_hdl)
{
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    assert(ib_mem_hdl);

    if (ib_mem_hdl->wr_queue.empty()) {
        log_debug(nnti_debug_level, "work request queue is empty");
    } else {
        wr_queue_iter_t i;
        for (i=ib_mem_hdl->wr_queue.begin(); i != ib_mem_hdl->wr_queue.end(); i++) {
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
    ib_memory_handle *ib_mem_hdl=NULL;

    log_debug(nnti_debug_level, "enter");

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    if (ib_mem_hdl->wr_queue.empty()) {
        rc=TRUE;
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static int8_t is_buf_op_complete(
        const NNTI_buffer_t *reg_buf)
{
    int8_t rc=FALSE;
    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
    assert(ib_mem_hdl);

    if (is_wr_queue_empty(reg_buf) == TRUE) {
        log_debug(nnti_debug_level, "work request queue is empty - return FALSE");
        rc=FALSE;
    } else {
        wr=ib_mem_hdl->wr_queue.front();
        assert(wr);

        rc = is_wr_complete(wr);
    }

    log_debug(nnti_debug_level, "exit (rc=%d)", rc);
    return(rc);
}

static int8_t is_any_buf_op_complete(
        const NNTI_buffer_t **buf_list,
        const uint32_t        buf_count,
        uint32_t             *which)
{
    int8_t rc=FALSE;

    log_debug(nnti_debug_level, "enter");

    for (int i=0;i<buf_count;i++) {
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

    for (int i=0;i<buf_count;i++) {
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

static void create_status(
        const NNTI_buffer_t  *reg_buf,
        const NNTI_buf_ops_t  remote_op,
        NNTI_result_t         nnti_rc,
        NNTI_status_t        *status)
{
    ib_connection    *conn      =NULL;
    ib_memory_handle *ib_mem_hdl=NULL;
    ib_work_request  *wr        =NULL;

    memset(status, 0, sizeof(NNTI_status_t));
    status->op     = remote_op;
    status->result = nnti_rc;
    if (nnti_rc==NNTI_OK) {
        ib_mem_hdl=(ib_memory_handle *)reg_buf->transport_private;
        assert(ib_mem_hdl);
        wr=ib_mem_hdl->wr_queue.front();
        assert(wr);

//        print_ack_buf(&wr->ack);
//        print_wr(wr);

        conn = get_conn_qpn(wr->last_wc.qp_num);

        status->start  = (uint64_t)reg_buf->payload;
        switch (wr->last_op) {
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
        switch (wr->last_op) {
            case IB_OP_NEW_REQUEST:
//                status->offset = wr->offset;
                status->offset = wr->last_wc.wr_id*transport_global_data.req_queue.req_size;
                status->length = wr->last_wc.byte_len;
                break;
            case IB_OP_SEND_REQUEST:
            case IB_OP_SEND_BUFFER:
            case IB_OP_RECEIVE:
                status->offset = 0;
                status->length = wr->last_wc.byte_len;
                break;
            case IB_OP_GET_INITIATOR:
            case IB_OP_PUT_INITIATOR:
                status->offset = wr->offset;
                status->length = wr->length;
                break;
            case IB_OP_GET_TARGET:
            case IB_OP_PUT_TARGET:
                if (config.use_rdma_target_ack) {
                    status->offset = wr->ack.offset;
                    status->length = wr->ack.length;
                }
                break;
        }
    }
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

static void copy_peer(NNTI_peer_t *src, NNTI_peer_t *dest)
{
    log_debug(nnti_debug_level, "enter");

    strncpy(dest->url, src->url, NNTI_URL_LEN);
    dest->url[NNTI_URL_LEN-1]='\0';

    src->peer.transport_id                    =NNTI_TRANSPORT_IB;
    dest->peer.NNTI_remote_process_t_u.ib.addr=src->peer.NNTI_remote_process_t_u.ib.addr;
    dest->peer.NNTI_remote_process_t_u.ib.port=src->peer.NNTI_remote_process_t_u.ib.port;

    log_debug(nnti_debug_level, "exit");
}

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
    int i;
    int rc=NNTI_OK;
    trios_declare_timer(callTime);

    /* bring the two QPs up to RTR */
    trios_start_timer(callTime);
    transition_qp_to_ready(conn->req_qp.qp, conn->peer_req_qpn, conn->peer_lid);
    transition_qp_to_ready(conn->data_qp.qp, conn->data_qp.peer_qpn, conn->peer_lid);
    trios_stop_timer("transition_qp_to_ready", callTime);

    trios_start_timer(callTime);
    /* final sychronization to ensure both sides have posted RTRs */
    rc = tcp_exchange(sock, 0, &rc, &rc, sizeof(rc));
    trios_stop_timer("exch data", callTime);
}

static void transition_qp_to_ready(
        struct ibv_qp *qp,
        uint32_t peer_qpn,
        int peer_lid)
{
    enum ibv_qp_attr_mask mask;
    struct ibv_qp_attr attr;

    /* Transition QP to Init */
    mask = (enum ibv_qp_attr_mask)
       ( IBV_QP_STATE
     | IBV_QP_ACCESS_FLAGS
     | IBV_QP_PKEY_INDEX
     | IBV_QP_PORT );
    memset(&attr, 0, sizeof(attr));
    attr.qp_state = IBV_QPS_INIT;
    attr.qp_access_flags = IBV_ACCESS_LOCAL_WRITE |
            IBV_ACCESS_REMOTE_WRITE | IBV_ACCESS_REMOTE_READ;
    attr.pkey_index = 0;
    attr.port_num = transport_global_data.nic_port;
    if (ibv_modify_qp(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp to INIT state");
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
    attr.min_rnr_timer = 31;
    if (ibv_modify_qp(qp, &attr, mask)) {
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
    attr.timeout = 24;  /* 4.096us * 2^24 */
    attr.retry_cnt = 7;
    attr.rnr_retry = 7;
    if (ibv_modify_qp(qp, &attr, mask)) {
        log_error(nnti_debug_level, "failed to modify qp from RTR to RTS state");
    }

}

static void transition_connection_to_error(
        ib_connection *conn)
{
    int i;
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
    if (ibv_modify_qp(qp, &attr, IBV_QP_STATE)) {
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
    int i, rc;
    int flags=0;
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
    att.cap.max_recv_sge = 1;
    att.cap.max_send_sge = 1;
    att.qp_type          = IBV_QPT_RC;

    trios_start_timer(callTime);
    c->req_qp.qp = ibv_create_qp(transport_global_data.pd, &att);
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
    att.cap.max_recv_sge = 1;
    att.cap.max_send_sge = 1;
    att.qp_type          = IBV_QPT_RC;
    c->data_qp.qp = ibv_create_qp(transport_global_data.pd, &att);
    if (!c->data_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    if (ibv_req_notify_cq(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
    }
    c->data_qp.qpn     = c->data_qp.qp->qp_num;
    param_out.my_qpn = htonl(c->data_qp.qpn);

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

out:
    return rc;
}

static int new_server_connection(
        ib_connection *c,
        int sock)
{
    int i, rc;
    int flags=0;
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
    att.cap.max_recv_sge = 1;
    att.cap.max_send_sge = 1;
    att.qp_type          = IBV_QPT_RC;

    trios_start_timer(callTime);
    c->req_qp.qp = ibv_create_qp(transport_global_data.pd, &att);
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
    att.cap.max_recv_sge = 1;
    att.cap.max_send_sge = 1;
    att.qp_type          = IBV_QPT_RC;
    c->data_qp.qp = ibv_create_qp(transport_global_data.pd, &att);
    if (!c->data_qp.qp) {
        log_error(nnti_debug_level, "failed to create QP: %s", strerror(errno));
    }
    if (ibv_req_notify_cq(transport_global_data.data_cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
    }

    c->data_qp.qpn   = c->data_qp.qp->qp_num;
    param_out.my_qpn = htonl(c->data_qp.qpn);

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
    assert(connections_by_peer.find(key) == connections_by_peer.end());
    connections_by_peer[key] = conn;   // add to connection map
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
static void print_qpn_map()
{
    ib_connection     *conn=NULL;
    conn_by_qpn_iter_t i;
    for (i=connections_by_qpn.begin(); i != connections_by_qpn.end(); i++) {
        conn=i->second;
        log_debug(nnti_debug_level, "qpn_map key=%llu conn=%p (name=%s, addr=%llu, port=%llu)",
                i->first, conn, conn->peer_name, (uint64_t)conn->peer_addr, (uint64_t)conn->peer_port);
    }
}
static void print_peer_map()
{
    ib_connection      *conn=NULL;
    conn_by_peer_iter_t i;
    for (i=connections_by_peer.begin(); i != connections_by_peer.end(); i++) {
        addrport_key key=i->first;
        conn=i->second;
        log_debug(nnti_debug_level, "peer_map key=(%llu,%llu) conn=%p (name=%s, addr=%llu, port=%llu)",
                (uint64_t)key.addr, (uint64_t)key.port, conn, conn->peer_name, (uint64_t)conn->peer_addr, (uint64_t)conn->peer_port);
    }
}

static NNTI_result_t insert_buf_bufhash(NNTI_buffer_t *buf)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)buf->buffer_addr.NNTI_remote_addr_t_u.ib.buf);

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
    uint32_t h=hash6432shift((uint64_t)buf->buffer_addr.NNTI_remote_addr_t_u.ib.buf);
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

static NNTI_result_t insert_wr_wrhash(ib_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t h=hash6432shift((uint64_t)wr);

//    nthread_lock(&nnti_wr_wrhash_lock);
    assert(wr_by_wrhash.find(h) == wr_by_wrhash.end());
    wr_by_wrhash[h] = wr;
//    nthread_unlock(&nnti_wr_wrhash_lock);

    log_debug(nnti_debug_level, "wrhash work request added (wrhash=%x, wr=%p)", h, wr);

    return(rc);
}
static ib_work_request *get_wr_wrhash(const uint32_t wrhash)
{
    ib_work_request *wr=NULL;

    log_debug(nnti_debug_level, "looking for wrhash=%x", (uint64_t)wrhash);
    nthread_lock(&nnti_wr_wrhash_lock);
    if (wr_by_wrhash.find(wrhash) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[wrhash];
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

    if (wr != NULL) {
        log_debug(nnti_debug_level, "work request found (wrhash=%x, wr=%p)", (uint64_t)wrhash, wr);
        return wr;
    }

    log_debug(nnti_debug_level, "work request NOT found (wrhash=%x)", (uint64_t)wrhash);
//    print_wrhash_map();

    return(NULL);
}
static ib_work_request *del_wr_wrhash(ib_work_request *wr)
{
    uint32_t h=hash6432shift((uint64_t)wr);
    log_level debug_level = nnti_debug_level;

    nthread_lock(&nnti_wr_wrhash_lock);
    if (wr_by_wrhash.find(h) != wr_by_wrhash.end()) {
        wr = wr_by_wrhash[h];
    }

    if (wr != NULL) {
        log_debug(debug_level, "work request found (wrhash=%x, wr=%p)", h, wr);
        wr_by_wrhash.erase(h);
    } else {
        log_debug(debug_level, "work request NOT found (wrhash=%x)", h);
    }
    nthread_unlock(&nnti_wr_wrhash_lock);

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

static NNTI_result_t wr_pool_register(
        ib_work_request *wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */

    trios_declare_timer(callTime);

    struct ibv_mr *mr=NULL;

    uint32_t len;

    log_debug(nnti_debug_level, "enter");

    len = sizeof(wr->ack);

    trios_start_timer(callTime);
    mlock(&wr->ack, len);
    munlock(&wr->ack, len);
    trios_stop_timer("mlock", callTime);

    trios_start_timer(callTime);
    mr = ibv_reg_mr(
            transport_global_data.pd,
            &wr->ack,
            len,
            (ibv_access_flags)(IBV_ACCESS_LOCAL_WRITE | IBV_ACCESS_REMOTE_READ | IBV_ACCESS_REMOTE_WRITE));
    if (!mr) {
        log_error(nnti_debug_level, "failed to register memory region");
        perror("errno");
        return NNTI_EIO;
    }
    trios_stop_timer("register", callTime);

    wr->ack_mr=mr;

    log_debug(nnti_debug_level, "exit (mr==%p, addr %p, length %lu, lkey %x, rkey %x)...", mr, mr->addr, mr->length, mr->lkey, mr->rkey);

    return (rc);
}
static NNTI_result_t wr_pool_deregister(
        ib_work_request *wr)
{
    NNTI_result_t rc=NNTI_OK; /* return code */
    int i=0;
    int ibv_rc=0;
    trios_declare_timer(callTime);

    log_debug(nnti_debug_level, "enter");

    if (wr->ack_mr!=NULL) {
        trios_start_timer(callTime);
        ibv_rc=ibv_dereg_mr(wr->ack_mr);
        if (ibv_rc != 0) {
            log_error(nnti_debug_level, "deregistering the ACK buffer failed");
        }
        trios_stop_timer("deregister", callTime);
        wr->ack_mr=NULL;
    } else {
        log_debug(nnti_debug_level, "exit wr(%p) - not registered", wr);
        return(NNTI_OK);
    }

    log_debug(nnti_debug_level, "exit");

    return(NNTI_OK);
}
static NNTI_result_t wr_pool_init(uint32_t pool_size)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    ib_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    for (i=0;i<pool_size;i++) {
        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        assert(wr);
        rc=wr_pool_register(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to register target work request: rc=%d", rc);
            goto cleanup;
        }
        wr_pool_rdma_push(wr);

        wr=(ib_work_request *)calloc(1, sizeof(ib_work_request));
        assert(wr);
        wr_pool_sendrecv_push(wr);
    }

cleanup:
    log_debug(nnti_debug_level, "exit");

    return(rc);
}
static ib_work_request *wr_pool_rdma_pop(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    ib_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    if (!rdma_wr_pool.empty()) {
        wr=rdma_wr_pool.front();
        rdma_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static ib_work_request *wr_pool_sendrecv_pop(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    ib_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    if (!sendrecv_wr_pool.empty()) {
        wr=sendrecv_wr_pool.front();
        sendrecv_wr_pool.pop_front();
    }
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return(wr);
}
static void wr_pool_rdma_push(ib_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;

    log_debug(nnti_debug_level, "enter");

    wr->last_op            =0;
    wr->is_last_op_complete=FALSE;
    wr->op_state           =BUFFER_INIT;

    nthread_lock(&nnti_wr_pool_lock);
    rdma_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static void wr_pool_sendrecv_push(ib_work_request *wr)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;

    log_debug(nnti_debug_level, "enter");

    wr->last_op            =0;
    wr->is_last_op_complete=FALSE;
    wr->op_state           =BUFFER_INIT;

    nthread_lock(&nnti_wr_pool_lock);
    sendrecv_wr_pool.push_front(wr);
    nthread_unlock(&nnti_wr_pool_lock);

    log_debug(nnti_debug_level, "exit");

    return;
}
static NNTI_result_t wr_pool_fini(void)
{
    NNTI_result_t  rc=NNTI_OK;
    uint32_t i;
    ib_work_request *wr=NULL;

    log_debug(nnti_debug_level, "enter");

    nthread_lock(&nnti_wr_pool_lock);
    while (!rdma_wr_pool.empty()) {
        wr=rdma_wr_pool.front();
        rdma_wr_pool.pop_front();
        assert(wr);
        rc=wr_pool_deregister(wr);
        if (rc!=NNTI_OK) {
            log_error(nnti_debug_level, "failed to deregister target work request: rc=%d", rc);
            goto cleanup;
        }
        free(wr);
    }
    while (!sendrecv_wr_pool.empty()) {
        wr=sendrecv_wr_pool.front();
        sendrecv_wr_pool.pop_front();
        assert(wr);
        free(wr);
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
    conn_by_qpn_iter_t qpn_iter;
    for (qpn_iter = connections_by_qpn.begin(); qpn_iter != connections_by_qpn.end(); qpn_iter++) {
        log_debug(debug_level, "close connection (qpn=%llu)", qpn_iter->first);
        close_connection(qpn_iter->second);
        connections_by_qpn.erase(qpn_iter);
    }
    nthread_unlock(&nnti_conn_qpn_lock);

    nthread_lock(&nnti_conn_peer_lock);
    conn_by_peer_iter_t peer_iter;
    for (peer_iter = connections_by_peer.begin(); peer_iter != connections_by_peer.end(); peer_iter++) {
        log_debug(debug_level, "close connection (peer.addr=%llu)", peer_iter->first.addr);
//        close_connection(peer_iter->second);
        connections_by_peer.erase(peer_iter);
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
    int i;

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
        rc=ibv_destroy_qp(c->req_qp.qp);
        if (rc < 0)
            log_error(nnti_debug_level, "failed to destroy QP");
    }
    if (c->data_qp.qp) {
        rc=ibv_destroy_qp(c->data_qp.qp);
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
        char         *peer_hostname = strdup(inet_ntoa(ssin.sin_addr));
        NNTI_ip_addr  peer_addr  = ssin.sin_addr.s_addr;
        NNTI_tcp_port peer_port  = ntohs(ssin.sin_port);

        conn = (ib_connection *)calloc(1, sizeof(ib_connection));
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
        insert_conn_qpn(conn->req_qp.qpn, conn);
        insert_conn_qpn(conn->data_qp.qpn, conn);
        insert_conn_peer(&peer, conn);

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
    while(!done) {
        if (check_for_waiting_connection() != NNTI_OK) {
            done=true;
        }
    }

    return(NNTI_OK);
}


static struct ibv_device *get_ib_device(void)
{
    struct ibv_device *dev;
    struct ibv_device **dev_list;
    int dev_count=0;

    dev_list = ibv_get_device_list(&dev_count);
    if (dev_count == 0)
        return NULL;
    if (dev_count > 1) {
                log_debug(nnti_debug_level, "found %d devices, defaulting the dev_list[0] (%p)", dev_count, dev_list[0]);
    }
    dev = dev_list[0];
    ibv_free_device_list(dev_list);

    return dev;
}

static void print_wc(const struct ibv_wc *wc)
{
    if (wc->status != 0) {
        log_error(nnti_debug_level, "wc=%p, wc.opcode=%d, wc.status=%d (%s), wc.wr_id=%lx, wc.vendor_err=%u, wc.byte_len=%u, wc.qp_num=%u, wc.imm_data=%x, wc.src_qp=%u",
            wc,
            wc->opcode,
            wc->status,
            ibv_wc_status_str(wc->status),
            wc->wr_id,
            wc->vendor_err,
            wc->byte_len,
            wc->qp_num,
            wc->imm_data,
            wc->src_qp);
    } else {
        log_debug(nnti_debug_level, "wc=%p, wc.opcode=%d, wc.status=%d (%s), wc.wr_id=%lx, wc.vendor_err=%u, wc.byte_len=%u, wc.qp_num=%u, wc.imm_data=%x, wc.src_qp=%u",
            wc,
            wc->opcode,
            wc->status,
            ibv_wc_status_str(wc->status),
            wc->wr_id,
            wc->vendor_err,
            wc->byte_len,
            wc->qp_num,
            wc->imm_data,
            wc->src_qp);
    }
}

static NNTI_result_t poll_comp_channel(
        struct ibv_comp_channel *comp_channel,
        struct ibv_cq           *cq,
        int timeout)
{
    NNTI_result_t rc=NNTI_OK;
    int poll_rc=0;
    struct pollfd my_pollfd;

    int retries_left=3;

    struct ibv_cq *ev_cq;
    void          *ev_ctx;

    log_debug(nnti_debug_level, "enter comp_channel(%p)", comp_channel);

    /*
     * poll the channel until it has an event and sleep ms_timeout
     * milliseconds between any iteration
     */
    my_pollfd.fd      = comp_channel->fd;
    my_pollfd.events  = POLLIN;
    my_pollfd.revents = 0;
    log_debug(nnti_debug_level, "polling with timeout==%d", timeout);

    // TODO Put a check for errno==EINTR to deal with timing interrupts from HPCToolkit
    do {
        poll_rc = poll(&my_pollfd, 1, timeout);
    } while ((poll_rc < 0) && (errno == EINTR));

    if (poll_rc == 0) {
        log_debug(nnti_debug_level, "poll timed out: %d", my_pollfd.revents);
        rc = NNTI_ETIMEDOUT;
        goto cleanup;
    } else if (poll_rc < 0) {
        log_error(nnti_debug_level, "poll error: poll_rc=%d (%s)", poll_rc, strerror(errno));
        rc = NNTI_EIO;
        goto cleanup;
    } else {
        log_debug(nnti_debug_level, "poll success: poll_rc=%d ; revents=%d", poll_rc, my_pollfd.revents);
    }

    log_debug(nnti_debug_level, "completion channel poll complete - %d event(s) waiting", my_pollfd.revents);

try_again:
    if (ibv_get_cq_event(comp_channel, &ev_cq, &ev_ctx) == 0) {
        log_debug(nnti_debug_level, "got event from comp_channel for cq=%p", ev_cq);
        ibv_ack_cq_events(ev_cq, 1);
        log_debug(nnti_debug_level, "ACKed event on cq=%p", ev_cq);
        rc = NNTI_OK;
    } else {
        if (errno == EAGAIN) {
            if (retries_left > 0) {
                retries_left--;
                goto try_again;
            } else {
                rc = NNTI_EAGAIN;
                goto cleanup;
            }
        }
        log_error(nnti_debug_level, "ibv_get_cq_event failed (ev_cq==%p): %s",
                ev_cq, strerror(errno));
        rc = NNTI_EIO;
        goto cleanup;
    }

    if (ev_cq != cq) {
        log_error(nnti_debug_level, "message received for unknown CQ (ev_cq==%p): %s",
                ev_cq, strerror(errno));
        rc = NNTI_EIO;
        goto cleanup;
    }

cleanup:
    if (ibv_req_notify_cq(cq, 0)) {
        log_error(nnti_debug_level, "Couldn't request CQ notification: %s", strerror(errno));
        rc = NNTI_EIO;
    }

    log_debug(nnti_debug_level, "exit");
    return(rc);
}

static void print_ib_conn(ib_connection *c)
{
    int i=0;
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
    c->use_wr_pool        =false;
    c->use_rdma_target_ack=true;
}

static void config_get_from_env(nnti_ib_config *c)
{
    char *env_str=NULL;

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
}

//static void print_wr(ib_work_request *wr)
//{
//    log_debug(nnti_debug_level, "wr (op=%llu ; offset=%llu ; length=%llu)",
//            (uint64_t)wr->last_op,
//            wr->offset,
//            wr->length);
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
